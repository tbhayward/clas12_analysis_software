#pragma once
// Generated-domain integration. Detector cuts are corrected through acceptance.
// Phi factors out of the purely kinematic volume.
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <utility>

namespace dvcs_phase_space {
inline constexpr double proton_mass = 0.9382720813;
inline constexpr double generated_y_min = 0.19;
inline constexpr double generated_y_max = 0.85;
inline constexpr double energy_10p6 = 10.60;
inline constexpr double energy_10p2 = 10.1998;
struct Box { double x0,x1,q0,q1,t0,t1; };
inline double minus_t_min(double x,double q) {
    const double eps=4*proton_mass*proton_mass*x*x/q;
    return q*(2*(1-x)*(1-std::sqrt(1+eps))+eps)/(4*x*(1-x)+eps);
}
inline std::pair<double,double> q_range(double x,const Box& b,double energy,
                                      bool apply_w,double w_min) {
    if (!(x>0 && x<1)) return {0,0};
    double lo=std::max(b.q0,2*proton_mass*energy*x*generated_y_min);
    const double hi=std::min(b.q1,2*proton_mass*energy*x*generated_y_max);
    if(apply_w)lo=std::max(lo,(w_min*w_min-proton_mass*proton_mass)*x/(1-x));
    return {lo,hi};
}
inline std::vector<double> x_breaks(const Box& b,double energy,bool apply_w,double w_min) {
    std::vector<double> v{b.x0,b.x1};
    auto add=[&](double x){if(std::isfinite(x)&&x>b.x0&&x<b.x1)v.push_back(x);};
    for(double y:{generated_y_min,generated_y_max}) {
        const double slope=2*proton_mass*energy*y;
        for(double q:{b.q0,b.q1})add(q/slope);
        if(apply_w)add(1-(w_min*w_min-proton_mass*proton_mass)/slope);
    }
    if(apply_w)for(double q:{b.q0,b.q1})add(q/(q+w_min*w_min-proton_mass*proton_mass));
    std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());return v;
}
template<class F> double refine(const F& f,double a,double b,double fa,double fm,double fb,
                               double whole,double eps,int depth) {
    const double m=(a+b)/2,l=(a+m)/2,r=(m+b)/2;
    const double fl=f(l),fr=f(r);
    const double left=(m-a)*(fa+4*fl+fm)/6,right=(b-m)*(fm+4*fr+fb)/6;
    const double delta=left+right-whole;
    if(std::abs(delta)<=15*eps)return left+right+delta/15;
    if(depth==0)throw std::runtime_error("Kinematic volume integration did not converge");
    return refine(f,a,m,fa,fl,fm,left,eps/2,depth-1)+refine(f,m,b,fm,fr,fb,right,eps/2,depth-1);
}
template<class F> double integrate(const F& f,double a,double b,double eps) {
    if(!(b>a))return 0;
    const double fa=f(a),fm=f((a+b)/2),fb=f(b);
    return refine(f,a,b,fa,fm,fb,(b-a)*(fa+4*fm+fb)/6,eps,26);
}
inline double allowed_fraction(const Box& b,double energy,bool apply_w=true,double w_min=2) {
    if(!(b.x1>b.x0&&b.q1>b.q0&&b.t1>b.t0&&energy>0))return 0;
    auto at_x=[&](double x){
        const auto bounds=q_range(x,b,energy,apply_w,w_min);
        if(!(bounds.second>bounds.first))return 0.;
        auto at_q=[&](double q){
            return std::max(0.,b.t1-std::max(b.t0,minus_t_min(x,q)))/(b.t1-b.t0);
        };
        return integrate(at_q,bounds.first,bounds.second,1e-11*(bounds.second-bounds.first))/(b.q1-b.q0);
    };
    const auto breaks=x_breaks(b,energy,apply_w,w_min);double integral=0;
    for(size_t i=1;i<breaks.size();++i)
        integral+=integrate(at_x,breaks[i-1],breaks[i],1e-10*(breaks[i]-breaks[i-1]));
    return std::clamp(integral/(b.x1-b.x0),0.,1.);
}
} // namespace dvcs_phase_space
