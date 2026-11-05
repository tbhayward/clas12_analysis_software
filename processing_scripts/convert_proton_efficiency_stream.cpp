#include <TFile.h>
#include <TTree.h>

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <memory>
#include <poll.h>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

std::vector<std::string> split(const std::string &line) {
    std::istringstream input(line);
    std::vector<std::string> fields;
    std::string field;
    while (input >> field) fields.push_back(field);
    return fields;
}

struct StreamTree {
    std::string tree_name;
    std::string buffer;
    std::vector<std::string> columns;
    std::vector<double> values;
    std::unique_ptr<TTree> tree;
    long long rows = 0;
    bool eof = false;

    explicit StreamTree(std::string name) : tree_name(std::move(name)) {}

    void process_line(const std::string &raw_line) {
        std::string line = raw_line;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) return;
        const auto fields = split(line);
        if (columns.empty()) {
            if (fields.empty()) throw std::runtime_error(tree_name + " stream has an empty header");
            std::set<std::string> unique(fields.begin(), fields.end());
            if (unique.size() != fields.size()) throw std::runtime_error(tree_name + " header has duplicate columns");
            columns = fields;
            values.assign(columns.size(), 0.0);
            tree = std::make_unique<TTree>(tree_name.c_str(), tree_name.c_str());
            // The StreamTree owns this object. Detach it from the current TFile
            // so closing the file cannot delete it a second time.
            tree->SetDirectory(nullptr);
            for (std::size_t i = 0; i < columns.size(); ++i) {
                const std::string leaf = columns[i] + "/D";
                tree->Branch(columns[i].c_str(), &values[i], leaf.c_str());
            }
            return;
        }
        if (fields.size() != columns.size()) {
            throw std::runtime_error(tree_name + " row has " + std::to_string(fields.size()) +
                                     " fields; expected " + std::to_string(columns.size()));
        }
        for (std::size_t i = 0; i < fields.size(); ++i) {
            char *end = nullptr;
            errno = 0;
            values[i] = std::strtod(fields[i].c_str(), &end);
            if (errno != 0 || end == fields[i].c_str() || *end != '\0') {
                throw std::runtime_error(tree_name + " contains a nonnumeric value in " + columns[i] +
                                         ": " + fields[i]);
            }
        }
        tree->Fill();
        ++rows;
    }

    void consume_lines() {
        std::size_t newline = 0;
        while ((newline = buffer.find('\n')) != std::string::npos) {
            process_line(buffer.substr(0, newline));
            buffer.erase(0, newline + 1);
        }
    }

    void finish() {
        if (!buffer.empty()) {
            process_line(buffer);
            buffer.clear();
        }
        if (!tree) throw std::runtime_error(tree_name + " stream ended before its header");
        eof = true;
    }
};

void set_nonblocking(int fd) {
    const int flags = fcntl(fd, F_GETFL, 0);
    if (flags < 0 || fcntl(fd, F_SETFL, flags | O_NONBLOCK) < 0) {
        throw std::runtime_error("cannot make FIFO nonblocking: " + std::string(std::strerror(errno)));
    }
}

void drain(int fd, StreamTree &stream) {
    char chunk[65536];
    while (true) {
        const ssize_t count = read(fd, chunk, sizeof(chunk));
        if (count > 0) {
            stream.buffer.append(chunk, static_cast<std::size_t>(count));
            stream.consume_lines();
            continue;
        }
        if (count == 0) {
            stream.finish();
            return;
        }
        if (errno == EAGAIN || errno == EWOULDBLOCK) return;
        if (errno == EINTR) continue;
        throw std::runtime_error(stream.tree_name + " FIFO read failed: " + std::string(std::strerror(errno)));
    }
}

}  // namespace

int validate_root(const char *path) {
    TFile input(path, "READ");
    if (input.IsZombie()) throw std::runtime_error("cannot open ROOT file for validation");
    auto *tags = dynamic_cast<TTree *>(input.Get("Tags"));
    auto *protons = dynamic_cast<TTree *>(input.Get("Protons"));
    if (!tags || !protons) throw std::runtime_error("missing Tags or Protons tree");
    for (const char *branch : {"tag_sequence", "rec_e_p_gt_2", "truth_e_p_gt_2"}) {
        if (!tags->GetBranch(branch)) throw std::runtime_error(std::string("Tags missing branch ") + branch);
    }
    for (const char *branch : {"tag_sequence", "selected", "match_angle"}) {
        if (!protons->GetBranch(branch)) throw std::runtime_error(std::string("Protons missing branch ") + branch);
    }
    std::cout << "[PROTON_EFFICIENCY_ROOT_VALID] tags=" << tags->GetEntries()
              << " protons=" << protons->GetEntries() << "\n";
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 3 && std::string(argv[1]) == "--validate") {
        try {
            return validate_root(argv[2]);
        } catch (const std::exception &error) {
            std::cerr << "ERROR: " << error.what() << "\n";
            return 1;
        }
    }
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " TAG_FIFO PROTON_FIFO OUTPUT.root\n"
                  << "   or: " << argv[0] << " --validate OUTPUT.root\n";
        return 2;
    }
    int tag_fd = -1, proton_fd = -1;
    try {
        // Blocking opens deliberately synchronize with FileWriter's tag-then-proton
        // open order. Once both ends exist, poll drains both streams concurrently.
        tag_fd = open(argv[1], O_RDONLY);
        if (tag_fd < 0) throw std::runtime_error("cannot open tag FIFO: " + std::string(std::strerror(errno)));
        proton_fd = open(argv[2], O_RDONLY);
        if (proton_fd < 0) throw std::runtime_error("cannot open proton FIFO: " + std::string(std::strerror(errno)));
        set_nonblocking(tag_fd);
        set_nonblocking(proton_fd);

        TFile output(argv[3], "RECREATE");
        if (output.IsZombie()) throw std::runtime_error("cannot create ROOT output");
        StreamTree tags("Tags"), protons("Protons");
        while (!tags.eof || !protons.eof) {
            pollfd descriptors[2] = {
                {tags.eof ? -1 : tag_fd, POLLIN, 0},
                {protons.eof ? -1 : proton_fd, POLLIN, 0},
            };
            const int ready = poll(descriptors, 2, -1);
            if (ready < 0) {
                if (errno == EINTR) continue;
                throw std::runtime_error("poll failed: " + std::string(std::strerror(errno)));
            }
            if (!tags.eof && (descriptors[0].revents & (POLLIN | POLLHUP | POLLERR))) drain(tag_fd, tags);
            if (!protons.eof && (descriptors[1].revents & (POLLIN | POLLHUP | POLLERR))) drain(proton_fd, protons);
        }
        output.cd();
        tags.tree->Write();
        protons.tree->Write();
        output.Close();
        close(tag_fd);
        close(proton_fd);
        std::cout << "[PROTON_EFFICIENCY_ROOT_ROWS] tags=" << tags.rows
                  << " protons=" << protons.rows << "\n";
        return 0;
    } catch (const std::exception &error) {
        if (tag_fd >= 0) close(tag_fd);
        if (proton_fd >= 0) close(proton_fd);
        std::cerr << "ERROR: " << error.what() << "\n";
        return 1;
    }
}
