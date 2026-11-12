#ifndef EXCLUSIVITY_RESOLUTION_CORE_H
#define EXCLUSIVITY_RESOLUTION_CORE_H
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace resolution_study {
constexpr double DEG=180./3.14159265358979323846;
constexpr size_t NV=7;
inline const std::array<std::string,NV> VARIABLES={"Delta_phi","theta","theta_gamma_gamma","pTmiss","Emiss2","Mx2","Mx2_2"};
struct Window {
    bool active=false,upper=false;
    double low=0,high=0;
    bool pass(double value) const {
        return !active || (std::isfinite(value) && value<=high && (upper || value>=low));
    }
};
inline unsigned failures(const std::array<double,NV>& values,const std::array<Window,NV>& cuts) {
    unsigned mask=0;
    for(size_t v=0;v<NV;++v) if(!cuts[v].pass(values[v])) mask|=1u<<v;
    return mask;
}
inline bool leave_one_out(unsigned mask,size_t variable) {return (mask & ~(1u<<variable))==0;}
struct Histogram {
    int bins; double low,high;
    // Includes under/overflow; nonfinite values have a separate counter.
    std::vector<uint64_t> counts;
    std::vector<double> edges;
    uint64_t nonfinite=0,pass_own_cut=0;
    Histogram(int n,double l,double h):bins(n),low(l),high(h),counts(n+2) {
        if(n<2 || !std::isfinite(l) || !std::isfinite(h) || h<=l) throw std::runtime_error("Bad histogram range");
    }
    explicit Histogram(const std::vector<double>& e):Histogram(int(e.size())-1,e.front(),e.back()) {
        if(!std::is_sorted(e.begin(),e.end()) || std::adjacent_find(e.begin(),e.end())!=e.end()) throw std::runtime_error("Invalid nonuniform histogram edges");
        edges=e;
    }
    void add(double value,bool own_pass) {
        if(!std::isfinite(value)) {++nonfinite;return;}
        int bin=value<low?0:value>=high?bins+1:1+std::min(bins-1,int((value-low)/(high-low)*bins));
        if(!edges.empty() && value>=low && value<high) bin=int(std::upper_bound(edges.begin(),edges.end(),value)-edges.begin());
        ++counts[bin]; if(own_pass) ++pass_own_cut;
    }
};
struct JointGrid {
    int offset;
    std::vector<int> indices;
    std::vector<std::vector<double>> edges;
    int find(const std::array<double,8>& values) const {
        int cell=0;
        for(size_t i=0;i<indices.size();++i) {
            double value=values[indices[i]];const auto& e=edges[i];
            if(!std::isfinite(value) || value<e.front() || value>=e.back()) return -1;
            cell=cell*int(e.size()-1)+int(std::upper_bound(e.begin(),e.end(),value)-e.begin())-1;
        }
        return offset+cell;
    }
};
struct Stratum {
    std::string axis;
    double low=0,high=0;
    bool contains(double value) const {return std::isfinite(value) && value>=low && value<high;}
};
}
#endif
