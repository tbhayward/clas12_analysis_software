// periods.h
#pragma once
#include <vector>
#include <string>

struct PeriodDef {
    const char* label;    // used in filenames, logs, CLI (e.g. "fa18_inb")
    const char* tree_key; // exact TTree map key (e.g. "DVCS_Fa18_inb")
};

inline const std::vector<PeriodDef>& CANONICAL_PERIODS() {
    static const std::vector<PeriodDef> v = {
        {"fa18_inb",      "DVCS_Fa18_inb"},
        {"fa18_inb_supp", "DVCS_Fa18_inb_supp"},
        {"fa18_out",      "DVCS_Fa18_out"},
        {"sp18_inb",      "DVCS_Sp18_inb"},
        {"sp18_out",      "DVCS_Sp18_out"},
        {"sp19_inb",      "DVCS_Sp19_inb"}
    };
    return v;
}

inline const PeriodDef* find_period(const std::string& label) {
    for (const auto& p : CANONICAL_PERIODS()) {
        if (label == p.label) return &p;
    }
    return nullptr;
}

static constexpr const char* SUF_EPPI0  = "_eppi0";
static constexpr const char* SUF_REC_MC = "_rec_mc";
static constexpr const char* SUF_BKG    = "_bkg";