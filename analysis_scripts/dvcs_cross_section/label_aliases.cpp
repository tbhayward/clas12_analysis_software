// Example usage inside update_total_counts_csv():
// We need to find columns like:
// "raw yield, ep->epg, (FD, FD), exp, fa18_inb, unpol" (and pos/neg)
// but we also accept "Fa18 Inb" etc.

auto find_col = [&](const std::vector<std::string>& header,
                    const std::string& topology,
                    const std::string& period,
                    const std::string& helicity) -> int {
    for (const auto& topo_try : topology_aliases(topology)) {
        for (const auto& per_try : period_aliases(period)) {
            for (const auto& hel_try : helicity_aliases(helicity)) {
                std::string candidate = "raw yield, ep->epg, " + topo_try + ", exp, " + per_try + ", " + hel_try;
                auto it = std::find(header.begin(), header.end(), candidate);
                if (it != header.end()) return int(it - header.begin());
            }
        }
    }
    return -1;
};

// And **before** any numeric read:
int c_un = find_col(header, "(FD, FD)", "fa18_inb", "unpol");
// ... similarly c_pos, c_neg
if (c_un < 0 /*|| c_pos < 0 || c_neg < 0*/) {
    std::cerr << "[total_counts] WARN: missing raw-yield column(s) for (topo=" << "(FD, FD)"
              << ", period=" << "fa18_inb" << "). Trying aliases failed; skipping this combination.\n";
    // skip this (period,topology) cleanly instead of reading invalid columns
    continue; // or return false if you want to fail fast
}

// When reading a row value, guard bounds and parse defensively:
auto safe_as_double = [&](const std::vector<std::string>& row, int col) -> double {
    if (col < 0 || col >= int(row.size())) return 0.0;
    const char* s = row[col].c_str();
    char* endp = nullptr;
    double v = std::strtod(s, &endp);
    if (endp == s) { // no conversion
        // treat empty / NaN-like strings as zero (or log once)
        return 0.0;
    }
    return v;
};