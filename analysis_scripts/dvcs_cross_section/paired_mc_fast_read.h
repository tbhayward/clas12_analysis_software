#ifndef DVCS_PAIRED_MC_FAST_READ_H
#define DVCS_PAIRED_MC_FAST_READ_H

// Read the validity branch first; decompress full reconstructed kinematics only
// for valid entries. Works with TTree and TChain, including file transitions.
// The caller MUST already SetBranchAddress(validity_name, &validity_value),
// enable that branch and bind its other branches before constructing this.
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <stdexcept>
#include <string>

class PairedMcFastRead {
public:
    PairedMcFastRead(TTree* tree, const char* validity_name, double* validity_value)
        : tree_(tree), name_(validity_name), value_(validity_value) {}

    bool read(Long64_t global_entry) {
        const Long64_t local_entry = tree_->LoadTree(global_entry);
        if (local_entry < 0) return false;
        TTree* current = tree_;
        if (auto* chain = dynamic_cast<TChain*>(tree_)) current = chain->GetTree();
        if (!current) return false;
        if (current != last_tree_) {
            branch_ = current->GetBranch(name_.c_str());
            if (!branch_) throw std::runtime_error("Missing paired MC validity branch: " + name_);
            last_tree_ = current;
        }
        if (branch_->GetEntry(local_entry) <= 0) return false;
        if (*value_ < 0.5) return false;
        return tree_->GetEntry(global_entry) > 0;
    }
private:
    TTree* tree_;
    std::string name_;
    double* value_;
    TTree* last_tree_ = nullptr;
    TBranch* branch_ = nullptr;
};
#endif
