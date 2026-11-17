#ifndef DVCS_CURRENT_RESPONSE_PROVENANCE_H
#define DVCS_CURRENT_RESPONSE_PROVENANCE_H
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
namespace current_response_provenance {
inline const std::array<const char*,8>& bin_columns() {
    static const std::array<const char*,8> columns={"xBmin","xBmax","Q2min","Q2max","t_abs_min","t_abs_max","phimin","phimax"};
    return columns;
}
inline bool same_value(double a,double b) {
    return std::isfinite(a) && std::isfinite(b) &&
        std::fabs(a-b)<=1e-8*std::max({1.0,std::fabs(a),std::fabs(b)});
}
// Deterministic identity checksum, not a security/authentication primitive.
inline std::string file_hash(const std::string& path) {
    std::ifstream input(path,std::ios::binary);
    if(!input) throw std::runtime_error("Cannot read current-response provenance input: "+path);
    std::uint64_t hash=14695981039346656037ULL;
    char buffer[65536];
    while(input.read(buffer,sizeof buffer) || input.gcount())
        for(std::streamsize i=0;i<input.gcount();++i) {
            hash^=static_cast<unsigned char>(buffer[i]); hash*=1099511628211ULL;
        }
    if(!input.eof()) throw std::runtime_error("Cannot finish hashing: "+path);
    std::ostringstream out; out<<std::hex<<std::setw(16)<<std::setfill('0')<<hash;
    return out.str();
}
}
#endif
