// label_aliases.h
#ifndef LABEL_ALIASES_H
#define LABEL_ALIASES_H

#include <string>
#include <vector>
#include <algorithm>

inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
    return s;
}

inline std::string squash_ws(std::string s) {
    std::string out; out.reserve(s.size());
    bool was_ws = false;
    for (char c : s) {
        bool is_ws = (c==' ' || c=='\t' || c=='\n' || c=='\r');
        if (is_ws) { if (!was_ws) out.push_back(' '); was_ws = true; }
        else { out.push_back(c); was_ws = false; }
    }
    return out;
}

inline std::string canon_period(std::string s) {
    s = squash_ws(s);
    s = to_lower(s);
    // normalize separators
    for (char& c : s) if (c==' ') c = '_';
    // common aliases
    if (s == "fa18_inb_supplemental" || s == "fa18_inb_supp." || s == "fa18_inb_suppl") s = "fa18_inb_supp";
    return s;
}

inline std::vector<std::string> period_aliases(std::string s) {
    std::vector<std::string> out;
    std::string c = canon_period(s);
    out.push_back(c);
    // also offer space form and title-case for robustness
    std::string spacey = c; for (char& ch : spacey) if (ch=='_') ch=' ';
    out.push_back(spacey);                    // "Fa18 Inb"
    // simple title-case
    if (!spacey.empty()) {
        std::string t = spacey;
        bool neww = true;
        for (char& ch : t) {
            if (neww && std::isalpha((unsigned char)ch)) ch = std::toupper((unsigned char)ch);
            else ch = std::tolower((unsigned char)ch);
            neww = (ch==' ');
        }
        out.push_back(t);                     // "Fa18 Inb"
    }
    // also accept exact user known variants
    if (c=="fa18_inb") { out.push_back("Fa18 Inb"); }
    if (c=="fa18_out") { out.push_back("Fa18 Out"); }
    if (c=="sp18_inb") { out.push_back("Sp18 Inb"); }
    if (c=="sp18_out") { out.push_back("Sp18 Out"); }
    if (c=="fa18_inb_supp") { out.push_back("Fa18 Inb supp"); out.push_back("Fa18 Inb Supp"); }
    return out;
}

inline std::string canon_topology(std::string s) {
    // Keep the "(X, Y)" formatting but tolerate extra spaces
    // Collapse internal spaces around comma
    std::string t;
    for (char c : s) {
        if (c=='\t' || c=='\n' || c=='\r') continue;
        t.push_back(c);
    }
    // remove spaces after '(' and before ')'
    auto trim_inside = [](std::string& z){
        // remove multiple spaces
        std::string y; y.reserve(z.size());
        bool prev_space = false;
        for (char c : z) {
            if (c==' ') {
                if (!prev_space) y.push_back(' ');
                prev_space = true;
            } else { y.push_back(c); prev_space = false; }
        }
        z.swap(y);
    };
    trim_inside(t);
    // Normalize to "(FD, FD)" exactly if tokens match ignoring spaces.
    auto ts = to_lower(t);
    for (char& c : ts) if (c==' ') c = '\0';
    if (ts.find("(fd,fd)") != std::string::npos) return "(FD, FD)";
    if (ts.find("(cd,fd)") != std::string::npos) return "(CD, FD)";
    if (ts.find("(cd,ft)") != std::string::npos) return "(CD, FT)";
    return s; // fallback
}

inline std::vector<std::string> topology_aliases(std::string s) {
    std::vector<std::string> out;
    std::string c = canon_topology(s);
    out.push_back(c);
    // accept a few spacing variants explicitly
    if (c=="(FD, FD)") { out.push_back("(FD,FD)"); out.push_back("( FD , FD )"); }
    if (c=="(CD, FD)") { out.push_back("(CD,FD)"); out.push_back("( CD , FD )"); }
    if (c=="(CD, FT)") { out.push_back("(CD,FT)"); out.push_back("( CD , FT )"); }
    return out;
}

inline std::vector<std::string> helicity_aliases(std::string h) {
    std::string c = to_lower(squash_ws(h));
    if (c=="unpol" || c=="0") return {"unpol","Unpol","0"};
    if (c=="pos"   || c=="+" ) return {"pos","Pos","+","plus","Plus"};
    if (c=="neg"   || c=="-" ) return {"neg","Neg","-","minus","Minus"};
    return {h};
}

#endif