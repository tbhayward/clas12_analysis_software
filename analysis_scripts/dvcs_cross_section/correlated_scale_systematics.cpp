// correlated_scale_systematics.cpp
// -----------------------------------------------------------------------------
// Final correlated, kinematically dependent scale systematic for the pass-2
// unpolarized DVCS cross section.
//
// The result combines two conceptually distinct studies:
//
//   1) the fitted current-dependent efficiency calibration uncertainty;
//   2) the residual period-to-period spread that remains after subtracting
//      the expected statistical and current-calibration contributions.
//
// The second term is learned from the four 10.6-GeV periods and parameterized
// versus the average proton polar angle.  The same residual detector-performance
// envelope is applied to Sp19 at its own average proton angle.
//
// IMPORTANT: the 2.16% charge/target-thickness normalization is written as a
// separate, uncorrelated normalization category.  It is NOT folded into the
// correlated-scale nuisance.
//
// The old combination-systematic columns are retained as diagnostics only.
// The authoritative final columns created here are:
//
//   run period residual sys frac, <sample>
//   correlated scale sys frac, <sample>
//   uncorrelated normalization sys frac, <sample>
//
// together with corresponding absolute cross-section errors.
// -----------------------------------------------------------------------------

#include "correlated_scale_systematics.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Csv {
    std::vector<std::string> header;
    std::unordered_map<std::string,int> index;
    std::vector<std::vector<std::string>> rows;
};

struct TupleValue {
    bool ok=false;
    double value=0.0;
    double stat=0.0;
};

static std::string trim(const std::string& s) {
    size_t a=0,b=s.size();
    while(a<b && std::isspace((unsigned char)s[a])) ++a;
    while(b>a && std::isspace((unsigned char)s[b-1])) --b;
    return s.substr(a,b-a);
}

static std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool q=false;
    for(size_t i=0;i<line.size();++i){
        const char c=line[i];
        if(c=='"'){
            if(q && i+1<line.size() && line[i+1]=='"'){
                cur.push_back('"'); ++i;
            } else q=!q;
        } else if(c==',' && !q){
            out.push_back(cur); cur.clear();
        } else cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

static std::string escape_csv(const std::string& s) {
    if(s.find_first_of(",\"\n\r")==std::string::npos) return s;
    std::string o="\"";
    for(char c:s){ if(c=='"')o+="\"\""; else o.push_back(c); }
    o+="\"";
    return o;
}

static Csv read_csv(const std::string& path) {
    std::ifstream in(path);
    if(!in) throw std::runtime_error("cannot open CSV: "+path);
    Csv c;
    std::string line;
    if(!std::getline(in,line)) throw std::runtime_error("empty CSV: "+path);
    c.header=split_csv(line);
    for(int i=0;i<(int)c.header.size();++i) c.index[c.header[i]]=i;
    while(std::getline(in,line)){
        auto r=split_csv(line);
        r.resize(c.header.size());
        c.rows.push_back(std::move(r));
    }
    return c;
}

static void write_csv(const std::string& path,const Csv& c) {
    const std::string tmp=path+".tmp";
    std::ofstream out(tmp);
    if(!out) throw std::runtime_error("cannot write CSV: "+tmp);
    for(size_t i=0;i<c.header.size();++i){
        if(i)out<<',';
        out<<escape_csv(c.header[i]);
    }
    out<<'\n';
    for(const auto&r:c.rows){
        for(size_t i=0;i<c.header.size();++i){
            if(i)out<<',';
            out<<escape_csv(i<r.size()?r[i]:std::string());
        }
        out<<'\n';
    }
    out.close();
    fs::rename(tmp,path);
}

static double num(const std::string& s) {
    const std::string t=trim(s);
    if(t.empty()) return std::numeric_limits<double>::quiet_NaN();
    char* e=nullptr;
    const double v=std::strtod(t.c_str(),&e);
    return e==t.c_str()?std::numeric_limits<double>::quiet_NaN():v;
}

static TupleValue tuple_value(const std::string& raw) {
    TupleValue out;
    std::string s=trim(raw);
    if(s.size()<5 || s.front()!='(' || s.back()!=')') return out;
    s=s.substr(1,s.size()-2);

    std::vector<std::string> fields;
    std::stringstream ss(s);
    std::string field;
    while(std::getline(ss,field,',')) fields.push_back(trim(field));
    if(fields.size()<2) return out;

    const double v=num(fields[0]);
    const double e=num(fields[1]);
    if(!std::isfinite(v)||!std::isfinite(e)||e<=0.0) return out;

    out.ok=true;
    out.value=v;
    out.stat=e;
    return out;
}

static int ensure_col(Csv& c,const std::string& name) {
    auto it=c.index.find(name);
    if(it!=c.index.end()) return it->second;
    const int idx=(int)c.header.size();
    c.header.push_back(name);
    c.index[name]=idx;
    for(auto&r:c.rows) r.push_back("");
    return idx;
}

static std::string fmt(double v) {
    if(!std::isfinite(v)) return "";
    std::ostringstream o;
    o<<std::setprecision(12)<<v;
    return o.str();
}

static double quantile(std::vector<double> v,double q) {
    if(v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(),v.end());
    const double p=q*(v.size()-1);
    const size_t i=(size_t)std::floor(p);
    const size_t j=(size_t)std::ceil(p);
    const double f=p-i;
    return v[i]*(1.0-f)+v[j]*f;
}

static const std::array<std::string,4>& ten6_periods() {
    static const std::array<std::string,4> p={{
        "Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"
    }};
    return p;
}

static std::string xs_col(const std::string& period) {
    return "normed cross sections, ep->epg, exp, "+period+", unpol";
}

static std::string current_col(const std::string& period) {
    return "current dependence sys frac, "+period;
}

struct RatioPoint {
    std::string period;
    double theta=0.0;
    double xB=std::numeric_limits<double>::quiet_NaN();
    double Q2=std::numeric_limits<double>::quiet_NaN();
    double t_abs=std::numeric_limits<double>::quiet_NaN();
    double phi=std::numeric_limits<double>::quiet_NaN();
    double e_theta=std::numeric_limits<double>::quiet_NaN();
    double g_theta=std::numeric_limits<double>::quiet_NaN();
    double ratio=0.0;
    double ratio_stat=0.0;
    double current_frac=0.0;
};

struct PeriodMean {
    bool ok=false;
    int n=0;
    double ratio=0.0;
    double stat=0.0;
    double current_frac=0.0;
};

struct ThetaResult {
    double lo=0.0,hi=0.0,center=0.0;
    std::array<PeriodMean,4> period;
    double mean_scale=0.0;
    double s_obs=0.0;
    double s_stat=0.0;
    double s_current=0.0;
    double s_resid=0.0;
};

static bool row_reference(const Csv& c,
                          const std::vector<std::string>& row,
                          double& ref) {
    double sw=0.0,swx=0.0;
    for(const auto&p:ten6_periods()){
        const auto it=c.index.find(xs_col(p));
        if(it==c.index.end()) continue;
        const TupleValue v=tuple_value(row[(size_t)it->second]);
        if(!v.ok) continue;
        const double w=1.0/(v.stat*v.stat);
        sw+=w;
        swx+=w*v.value;
    }
    if(!(sw>0.0)) return false;
    ref=swx/sw;
    return std::isfinite(ref)&&std::fabs(ref)>0.0;
}

static std::vector<RatioPoint> build_ratio_points(const Csv& c) {
    const auto it_theta=c.index.find("p_theta, 10.6 GeV");
    if(it_theta==c.index.end())
        throw std::runtime_error("missing p_theta, 10.6 GeV");

    const auto ixB=c.index.find("xBavg, 10.6 GeV");
    const auto iQ2=c.index.find("Q2avg, 10.6 GeV");
    const auto itabs=c.index.find("t_abs_avg, 10.6 GeV");
    const auto iphi=c.index.find("phiavg, 10.6 GeV");
    const auto ieth=c.index.find("e_theta, 10.6 GeV");
    const auto igth=c.index.find("g_theta, 10.6 GeV");

    std::vector<RatioPoint> out;

    for(const auto&row:c.rows){
        const double theta=num(row[(size_t)it_theta->second]);
        if(!std::isfinite(theta)) continue;

        double ref=0.0;
        if(!row_reference(c,row,ref)) continue;

        for(const auto&p:ten6_periods()){
            const auto ix=c.index.find(xs_col(p));
            const auto ic=c.index.find(current_col(p));
            if(ix==c.index.end()||ic==c.index.end()) continue;

            const TupleValue v=tuple_value(row[(size_t)ix->second]);
            const double cf=num(row[(size_t)ic->second]);
            if(!v.ok||!std::isfinite(cf)||cf<0.0) continue;

            RatioPoint rp;
            rp.period=p;
            rp.theta=theta;
            rp.xB=(ixB==c.index.end()) ? std::numeric_limits<double>::quiet_NaN()
                                       : num(row[(size_t)ixB->second]);
            rp.Q2=(iQ2==c.index.end()) ? std::numeric_limits<double>::quiet_NaN()
                                       : num(row[(size_t)iQ2->second]);
            rp.t_abs=(itabs==c.index.end()) ? std::numeric_limits<double>::quiet_NaN()
                                            : num(row[(size_t)itabs->second]);
            rp.phi=(iphi==c.index.end()) ? std::numeric_limits<double>::quiet_NaN()
                                         : num(row[(size_t)iphi->second]);
            rp.e_theta=(ieth==c.index.end()) ? std::numeric_limits<double>::quiet_NaN()
                                             : num(row[(size_t)ieth->second]);
            rp.g_theta=(igth==c.index.end()) ? std::numeric_limits<double>::quiet_NaN()
                                             : num(row[(size_t)igth->second]);
            rp.ratio=v.value/ref;
            rp.ratio_stat=std::fabs(v.stat/ref);
            rp.current_frac=cf;

            if(std::isfinite(rp.ratio)&&
               std::isfinite(rp.ratio_stat)&&
               rp.ratio_stat>0.0) out.push_back(rp);
        }
    }
    return out;
}

static PeriodMean weighted_period_mean(
    const std::vector<RatioPoint>& points,
    const std::string& period,
    double lo,double hi) {

    double sw=0.0,swr=0.0,swf2=0.0;
    int n=0;

    for(const auto&p:points){
        if(p.period!=period || !(p.theta>=lo && p.theta<hi)) continue;
        const double w=1.0/(p.ratio_stat*p.ratio_stat);
        sw+=w;
        swr+=w*p.ratio;
        // Current calibration is correlated across bins, so do not average it
        // down as 1/sqrt(N).  A weighted RMS retains a conservative typical
        // magnitude for this theta interval.
        swf2+=w*p.current_frac*p.current_frac;
        ++n;
    }

    PeriodMean out;
    out.n=n;
    if(!(sw>0.0)||n==0) return out;

    out.ratio=swr/sw;
    out.stat=1.0/std::sqrt(sw);
    out.current_frac=std::sqrt(std::max(0.0,swf2/sw));
    out.ok=std::isfinite(out.ratio)&&std::isfinite(out.stat)&&
           std::isfinite(out.current_frac);
    return out;
}

static std::vector<ThetaResult> build_theta_reference(
    const std::vector<RatioPoint>& points,
    double width,
    int min_points) {

    if(points.empty()) return {};

    double xmin=std::numeric_limits<double>::infinity();
    double xmax=-std::numeric_limits<double>::infinity();
    for(const auto&p:points){
        xmin=std::min(xmin,p.theta);
        xmax=std::max(xmax,p.theta);
    }

    const double start=std::floor(xmin/width)*width;
    const double stop=std::ceil(xmax/width)*width;

    std::vector<ThetaResult> out;

    for(double lo=start;lo<stop;lo+=width){
        ThetaResult r;
        r.lo=lo;
        r.hi=lo+width;
        r.center=0.5*(r.lo+r.hi);

        bool usable=true;
        for(size_t ip=0;ip<ten6_periods().size();++ip){
            r.period[ip]=weighted_period_mean(
                points,ten6_periods()[ip],r.lo,r.hi);
            if(!r.period[ip].ok || r.period[ip].n<min_points) usable=false;
        }
        if(!usable) continue;

        for(const auto&p:r.period) r.mean_scale+=p.ratio;
        r.mean_scale/=4.0;
        if(!std::isfinite(r.mean_scale)||std::fabs(r.mean_scale)<=0.0) continue;

        double obs2=0.0;
        double stat_var_sum=0.0;
        double current_var_sum=0.0;

        for(const auto&p:r.period){
            const double rel=p.ratio/r.mean_scale-1.0;
            obs2+=rel*rel;
            stat_var_sum+=(p.stat/r.mean_scale)*(p.stat/r.mean_scale);
            current_var_sum+=p.current_frac*p.current_frac;
        }

        r.s_obs=std::sqrt(obs2/4.0);

        // Expected RMS spread across N independent measurements about their
        // sample mean: (N-1)/N^2 * sum(sigma_p^2), with N=4 here.
        constexpr double kScatterFactor=3.0/16.0;
        r.s_stat=std::sqrt(std::max(0.0,kScatterFactor*stat_var_sum));
        r.s_current=std::sqrt(std::max(0.0,kScatterFactor*current_var_sum));

        const double residual2=
            r.s_obs*r.s_obs-
            r.s_stat*r.s_stat-
            r.s_current*r.s_current;

        r.s_resid=std::sqrt(std::max(0.0,residual2));
        out.push_back(r);
    }

    return out;
}


enum class DiagnosticVariable {
    XB,
    Q2,
    TABS,
    PHI,
    ETHETA,
    PTHETA,
    GTHETA
};

static double diagnostic_value(
    const RatioPoint& p,
    DiagnosticVariable variable) {

    switch(variable) {
        case DiagnosticVariable::XB:     return p.xB;
        case DiagnosticVariable::Q2:     return p.Q2;
        case DiagnosticVariable::TABS:   return p.t_abs;
        case DiagnosticVariable::PHI:    return p.phi;
        case DiagnosticVariable::ETHETA: return p.e_theta;
        case DiagnosticVariable::PTHETA: return p.theta;
        case DiagnosticVariable::GTHETA: return p.g_theta;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

static PeriodMean weighted_period_mean_variable(
    const std::vector<RatioPoint>& points,
    const std::string& period,
    DiagnosticVariable variable,
    double lo,
    double hi) {

    double sw=0.0,swr=0.0,swf2=0.0;
    int n=0;

    for(const auto& p:points) {
        if(p.period!=period) continue;
        const double x=diagnostic_value(p,variable);
        if(!std::isfinite(x) || !(x>=lo && x<hi)) continue;

        const double w=1.0/(p.ratio_stat*p.ratio_stat);
        sw+=w;
        swr+=w*p.ratio;
        // Current-response uncertainties are correlated calibration effects.
        // Use the weighted RMS magnitude rather than allowing them to shrink as
        // 1/sqrt(N) when forming this diagnostic interval.
        swf2+=w*p.current_frac*p.current_frac;
        ++n;
    }

    PeriodMean out;
    out.n=n;
    if(!(sw>0.0)||n==0) return out;

    out.ratio=swr/sw;
    out.stat=1.0/std::sqrt(sw);
    out.current_frac=std::sqrt(std::max(0.0,swf2/sw));
    out.ok=std::isfinite(out.ratio)&&std::isfinite(out.stat)&&
           std::isfinite(out.current_frac);
    return out;
}

static std::vector<ThetaResult> build_variable_reference(
    const std::vector<RatioPoint>& points,
    DiagnosticVariable variable,
    double start,
    double stop,
    double width,
    int min_points) {

    std::vector<ThetaResult> out;

    for(double lo=start;lo<stop;lo+=width) {
        ThetaResult r;
        r.lo=lo;
        r.hi=std::min(stop,lo+width);
        r.center=0.5*(r.lo+r.hi);

        bool usable=true;
        for(size_t ip=0;ip<ten6_periods().size();++ip) {
            r.period[ip]=weighted_period_mean_variable(
                points,ten6_periods()[ip],variable,r.lo,r.hi);
            if(!r.period[ip].ok || r.period[ip].n<min_points) usable=false;
        }
        if(!usable) continue;

        for(const auto& p:r.period) r.mean_scale+=p.ratio;
        r.mean_scale/=4.0;
        if(!std::isfinite(r.mean_scale)||std::fabs(r.mean_scale)<=0.0) continue;

        double obs2=0.0;
        double stat_var_sum=0.0;
        double current_var_sum=0.0;

        for(const auto& p:r.period) {
            const double rel=p.ratio/r.mean_scale-1.0;
            obs2+=rel*rel;
            stat_var_sum+=(p.stat/r.mean_scale)*(p.stat/r.mean_scale);
            current_var_sum+=p.current_frac*p.current_frac;
        }

        r.s_obs=std::sqrt(obs2/4.0);

        constexpr double kScatterFactor=3.0/16.0;
        r.s_stat=std::sqrt(std::max(0.0,kScatterFactor*stat_var_sum));
        r.s_current=std::sqrt(std::max(0.0,kScatterFactor*current_var_sum));

        const double residual2=
            r.s_obs*r.s_obs-
            r.s_stat*r.s_stat-
            r.s_current*r.s_current;

        r.s_resid=std::sqrt(std::max(0.0,residual2));
        out.push_back(r);
    }

    return out;
}

static void draw_residual_all_variables(
    const fs::path& path,
    const std::vector<RatioPoint>& points,
    int min_points) {

    struct Spec {
        DiagnosticVariable variable;
        const char* title;
        double start;
        double stop;
        double width;
    };

    // These are diagnostic display intervals only.  They do not alter the
    // production theta_p parameterization or any analysis binning.
    const std::array<Spec,7> specs={{
        {DiagnosticVariable::XB,     "x_{B}",                 0.06, 0.60, 0.06},
        {DiagnosticVariable::Q2,     "Q^{2} (GeV^{2})",       1.0,  8.2,  0.8},
        {DiagnosticVariable::TABS,   "|t| (GeV^{2})",         0.0,  1.20, 0.15},
        {DiagnosticVariable::PHI,    "#phi (deg)",            0.0,  360., 30.0},
        {DiagnosticVariable::ETHETA, "#theta_{e} (deg)",      8.0,  26.0, 2.0},
        {DiagnosticVariable::PTHETA, "#theta_{p} (deg)",      20.0, 66.0, 4.0},
        {DiagnosticVariable::GTHETA, "#theta_{#gamma} (deg)", 2.0,  40.0, 4.0}
    }};

    TCanvas cv("c_residual_all_variables","",1500,900);
    cv.Divide(4,2,0.002,0.002);

    for(int is=0;is<7;++is) {
        const auto ref=build_variable_reference(
            points,
            specs[is].variable,
            specs[is].start,
            specs[is].stop,
            specs[is].width,
            min_points);

        cv.cd(is+1);
        gPad->SetLeftMargin((is%4==0)?0.16:0.13);
        gPad->SetRightMargin(0.035);
        gPad->SetBottomMargin((is>=4)?0.17:0.13);
        gPad->SetTopMargin(0.13);
        gPad->SetTicks(1,1);

        double ymax=0.0;
        for(const auto& r:ref) ymax=std::max(ymax,r.s_resid);
        ymax=std::max(0.03,1.25*ymax);

        TH1D frame(
            ("h_residual_diag_"+std::to_string(is)).c_str(),
            "",
            100,
            specs[is].start,
            specs[is].stop
        );
        frame.SetMinimum(0.0);
        frame.SetMaximum(100.0*ymax);
        frame.GetXaxis()->SetTitle(specs[is].title);
        frame.GetYaxis()->SetTitle("Residual period systematic (%)");
        frame.GetXaxis()->SetTitleSize(0.050);
        frame.GetYaxis()->SetTitleSize(0.046);
        frame.GetXaxis()->SetLabelSize(0.042);
        frame.GetYaxis()->SetLabelSize(0.039);
        frame.GetYaxis()->SetTitleOffset((is%4==0)?1.45:1.22);
        frame.DrawCopy();

        TGraph g((int)ref.size());
        for(int i=0;i<(int)ref.size();++i) {
            g.SetPoint(i,ref[(size_t)i].center,100.0*ref[(size_t)i].s_resid);
        }
        g.SetMarkerStyle(20);
        g.SetMarkerSize(0.9);
        g.SetMarkerColor(kBlue+1);
        g.SetLineColor(kBlue+1);
        g.SetLineWidth(2);
        if(!ref.empty()) g.DrawClone("PL SAME");

        TLatex lab;
        lab.SetNDC();
        lab.SetTextFont(42);
        lab.SetTextSize(0.035);
        const std::string panel=std::string("(")+char('a'+is)+")";
        lab.DrawLatex(0.18,0.84,panel.c_str());
    }

    // Use the final empty pad for a compact statement of what the plot tests.
    cv.cd(8);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.10);

    TLatex info;
    info.SetNDC();
    info.SetTextFont(42);
    info.SetTextSize(0.050);
    info.DrawLatex(0.08,0.82,"Diagnostic only");
    info.SetTextSize(0.037);
    info.DrawLatex(0.08,0.67,"s_{resid} = #sqrt{max(0,");
    info.DrawLatex(0.11,0.59,"s_{obs}^{2}-s_{stat}^{2}-s_{current}^{2})}");
    info.SetTextSize(0.033);
    info.DrawLatex(0.08,0.42,"Production assignment:");
    info.DrawLatex(0.08,0.34,"4^{#circ} #theta_{p} intervals");
    info.DrawLatex(0.08,0.22,"Other panels test whether");
    info.DrawLatex(0.08,0.14,"a qualitatively different");
    info.DrawLatex(0.08,0.06,"dependence is being missed.");

    cv.cd(0);
    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextAlign(22);
    title.SetTextSize(0.021);
    title.DrawLatex(
        0.50,0.975,
        "Residual run-period consistency systematic versus analysis variables"
    );
    title.SetTextSize(0.0155);
    title.DrawLatex(
        0.50,0.949,
        "Four 10.6 GeV periods; statistical and current-calibration contributions removed"
    );

    cv.SaveAs(path.string().c_str());
}

static double interpolate_residual(
    const std::vector<ThetaResult>& ref,double theta) {

    if(ref.empty()||!std::isfinite(theta))
        return std::numeric_limits<double>::quiet_NaN();

    if(theta<=ref.front().center) return ref.front().s_resid;
    if(theta>=ref.back().center) return ref.back().s_resid;

    for(size_t i=1;i<ref.size();++i){
        if(theta<=ref[i].center){
            const double x0=ref[i-1].center,x1=ref[i].center;
            const double f=(theta-x0)/(x1-x0);
            return ref[i-1].s_resid*(1.0-f)+ref[i].s_resid*f;
        }
    }
    return ref.back().s_resid;
}

static void write_theta_reference(
    const fs::path& path,
    const std::vector<ThetaResult>& ref) {

    std::ofstream out(path);
    out<<"theta_min,theta_max,theta_center";
    for(const auto&p:ten6_periods())
        out<<","<<escape_csv(p+" n")
           <<","<<escape_csv(p+" mean ratio")
           <<","<<escape_csv(p+" mean-ratio stat")
           <<","<<escape_csv(p+" current frac");
    out<<",mean_scale,s_obs,s_stat,s_current,s_residual,"
          "s_obs_percent,s_stat_percent,s_current_percent,s_residual_percent\n";

    for(const auto&r:ref){
        out<<r.lo<<','<<r.hi<<','<<r.center;
        for(const auto&p:r.period)
            out<<','<<p.n<<','<<p.ratio<<','<<p.stat<<','<<p.current_frac;
        out<<','<<r.mean_scale<<','<<r.s_obs<<','<<r.s_stat<<','
           <<r.s_current<<','<<r.s_resid<<','
           <<100*r.s_obs<<','<<100*r.s_stat<<','
           <<100*r.s_current<<','<<100*r.s_resid<<'\n';
    }
}

static void draw_theta_decomposition(
    const fs::path& path,
    const std::vector<ThetaResult>& ref) {

    if(ref.empty()) return;

    TCanvas c("c_correlated_scale_theta","",1100,780);
    c.SetLeftMargin(.13);c.SetRightMargin(.035);
    c.SetBottomMargin(.13);c.SetTopMargin(.15);c.SetTicks(1,1);

    const double xmin=ref.front().lo;
    const double xmax=ref.back().hi;
    double ymax=0.0;
    for(const auto&r:ref) ymax=std::max(ymax,r.s_obs);
    ymax=std::max(.05,1.25*ymax);

    TH1D frame("h_correlated_scale_theta","",100,xmin,xmax);
    frame.SetMinimum(0.0);frame.SetMaximum(100*ymax);
    frame.GetXaxis()->SetTitle("#theta_{p} (deg)");
    frame.GetYaxis()->SetTitle("Relative period spread (%)");
    frame.GetXaxis()->SetTitleSize(.045);frame.GetYaxis()->SetTitleSize(.043);
    frame.GetXaxis()->SetLabelSize(.039);frame.GetYaxis()->SetLabelSize(.038);
    frame.GetYaxis()->SetTitleOffset(1.25);
    frame.Draw();

    TGraph gobs((int)ref.size()),gstat((int)ref.size()),
           gcur((int)ref.size()),gres((int)ref.size());

    for(int i=0;i<(int)ref.size();++i){
        const double x=ref[(size_t)i].center;
        gobs.SetPoint(i,x,100*ref[(size_t)i].s_obs);
        gstat.SetPoint(i,x,100*ref[(size_t)i].s_stat);
        gcur.SetPoint(i,x,100*ref[(size_t)i].s_current);
        gres.SetPoint(i,x,100*ref[(size_t)i].s_resid);
    }

    gobs.SetMarkerStyle(20);gobs.SetMarkerColor(kBlack);gobs.SetLineColor(kBlack);
    gstat.SetMarkerStyle(24);gstat.SetMarkerColor(kGray+2);gstat.SetLineColor(kGray+2);
    gcur.SetMarkerStyle(22);gcur.SetMarkerColor(kBlue+1);gcur.SetLineColor(kBlue+1);
    gres.SetMarkerStyle(21);gres.SetMarkerColor(kRed+1);gres.SetLineColor(kRed+1);
    for(TGraph*g:{&gobs,&gstat,&gcur,&gres}){g->SetLineWidth(2);g->SetMarkerSize(.95);}

    gobs.Draw("PL SAME");gstat.Draw("PL SAME");gcur.Draw("PL SAME");gres.Draw("PL SAME");

    TLegend leg(.55,.62,.93,.82);
    leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextFont(42);leg.SetTextSize(.029);
    leg.AddEntry(&gobs,"Observed period spread","pl");
    leg.AddEntry(&gstat,"Expected statistical contribution","pl");
    leg.AddEntry(&gcur,"Expected current-calibration contribution","pl");
    leg.AddEntry(&gres,"Residual run-period uncertainty","pl");
    leg.Draw();

    TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.031);
    t.DrawLatex(.13,.955,"Residual run-period consistency systematic");
    t.SetTextSize(.022);
    t.DrawLatex(.13,.915,
                "Four 10.6 GeV periods; 4^{#circ} #theta_{p} intervals");

    c.SaveAs(path.string().c_str());
}

struct AxisSpec {
    std::string lo,hi,title;
};

static void draw_final_kinematic_summary(
    const Csv& c,
    const std::vector<double>& corr10,
    const std::vector<double>& corrsp,
    const fs::path& path) {

    const std::array<AxisSpec,4> specs={{
        {"xBmin","xBmax","x_{B}"},
        {"Q2min","Q2max","Q^{2} (GeV^{2})"},
        {"t_abs_min","t_abs_max","|t| (GeV^{2})"},
        {"phimin","phimax","#phi (deg)"}
    }};

    TCanvas cv("c_final_corr_scale_kin","",1450,1050);
    cv.Divide(2,2,.002,.002);

    for(int ia=0;ia<4;++ia){
        cv.cd(ia+1);
        gPad->SetLeftMargin((ia%2==0)?.14:.12);
        gPad->SetRightMargin(.035);
        gPad->SetBottomMargin((ia>=2)?.15:.12);
        gPad->SetTopMargin(.11);
        gPad->SetTicks(1,1);

        auto ilo=c.index.find(specs[ia].lo);
        auto ihi=c.index.find(specs[ia].hi);
        if(ilo==c.index.end()||ihi==c.index.end()) continue;

        struct Bucket{
            double lo=0,hi=0;
            std::vector<double> a,b;
        };
        std::map<std::pair<double,double>,Bucket> buckets;

        for(size_t r=0;r<c.rows.size();++r){
            const double lo=num(c.rows[r][ilo->second]);
            const double hi=num(c.rows[r][ihi->second]);
            if(!std::isfinite(lo)||!std::isfinite(hi)||!(hi>lo))continue;

            // For phi only, use 30-degree display bins to avoid the clutter
            // caused by unequal native phi-bin widths.
            double dlo=lo,dhi=hi;
            if(ia==3){
                const double center=.5*(lo+hi);
                dlo=30.0*std::floor(center/30.0);
                dhi=dlo+30.0;
            }

            auto&b=buckets[{dlo,dhi}];
            b.lo=dlo;b.hi=dhi;
            if(r<corr10.size()&&std::isfinite(corr10[r]))b.a.push_back(100*corr10[r]);
            if(r<corrsp.size()&&std::isfinite(corrsp[r]))b.b.push_back(100*corrsp[r]);
        }

        TGraphErrors ga,gb;
        int na=0,nb=0;
        double xmin=INFINITY,xmax=-INFINITY,ymax=0.0;

        for(const auto&kv:buckets){
            const auto&b=kv.second;
            const double x=.5*(b.lo+b.hi),ex=.5*(b.hi-b.lo);
            xmin=std::min(xmin,b.lo);xmax=std::max(xmax,b.hi);

            if(!b.a.empty()){
                const double med=quantile(b.a,.5);
                const double q16=quantile(b.a,.16),q84=quantile(b.a,.84);
                ga.SetPoint(na,x,med);
                ga.SetPointError(na,ex,.5*(q84-q16));
                ymax=std::max(ymax,q84);++na;
            }
            if(!b.b.empty()){
                const double med=quantile(b.b,.5);
                const double q16=quantile(b.b,.16),q84=quantile(b.b,.84);
                gb.SetPoint(nb,x,med);
                gb.SetPointError(nb,ex,.5*(q84-q16));
                ymax=std::max(ymax,q84);++nb;
            }
        }
        if(!(xmax>xmin)){xmin=0;xmax=1;}

        TH1D frame(("h_final_corr_scale_"+std::to_string(ia)).c_str(),"",100,xmin,xmax);
        frame.SetMinimum(0.0);frame.SetMaximum(std::max(5.0,1.20*ymax));
        frame.GetXaxis()->SetTitle(specs[ia].title.c_str());
        frame.GetYaxis()->SetTitle("Correlated scale systematic (%)");
        frame.GetXaxis()->SetTitleSize(.050);frame.GetYaxis()->SetTitleSize(.047);
        frame.GetXaxis()->SetLabelSize(.043);frame.GetYaxis()->SetLabelSize(.041);
        frame.GetYaxis()->SetTitleOffset(1.20);
        frame.DrawCopy();

        ga.SetMarkerStyle(20);ga.SetMarkerColor(kBlue+1);ga.SetLineColor(kBlue+1);
        gb.SetMarkerStyle(24);gb.SetMarkerColor(kRed+1);gb.SetLineColor(kRed+1);
        ga.SetMarkerSize(.95);gb.SetMarkerSize(.95);
        if(na)ga.DrawClone("P SAME");
        if(nb)gb.DrawClone("P SAME");

        TLatex p;p.SetNDC();p.SetTextFont(42);p.SetTextSize(.036);
        const std::string lab=std::string("(")+char('a'+ia)+")";
        p.DrawLatex(.18,.84,lab.c_str());

        if(ia==0){
            TLegend l(.56,.68,.92,.82);
            l.SetBorderSize(0);l.SetFillStyle(0);l.SetTextFont(42);l.SetTextSize(.026);
            l.AddEntry(&ga,"Combined 10.6 GeV","p");
            l.AddEntry(&gb,"Sp19 Inb","p");
            l.DrawClone();
        }
    }

    cv.cd(0);
    TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextAlign(22);t.SetTextSize(.022);
    t.DrawLatex(.50,.975,"Final correlated kinematic-dependent scale systematic");
    t.SetTextSize(.016);
    t.DrawLatex(.50,.949,
                "Current-dependent efficiency #oplus residual run-period consistency");

    cv.SaveAs(path.string().c_str());
}


static void draw_systematic_category_summary(
    const fs::path& path,
    const Csv& c,
    const std::vector<double>& corr10,
    const std::vector<double>& corrsp,
    double norm_frac) {

    struct Summary {
        double ptp=0.0;
        double current=0.0;
        double residual=0.0;
        double correlated=0.0;
    };

    auto summarize=[&](
        const std::string& sample,
        const std::vector<double>& corr,
        const std::string& xscol,
        const std::string& currentcol)->Summary {

        std::vector<double> ptp,current,residual,corrv;
        const auto ix=c.index.find(xscol);
        const auto iptp=c.index.find("Syst. err (point-to-point total)");
        const auto icur=c.index.find(currentcol);
        const auto ires=c.index.find("run period residual sys frac, "+sample);

        for(size_t i=0;i<c.rows.size();++i) {
            if(i<corr.size()&&std::isfinite(corr[i]))corrv.push_back(100.0*corr[i]);

            if(ix!=c.index.end()&&iptp!=c.index.end()) {
                const TupleValue x=tuple_value(c.rows[i][ix->second]);
                const double e=num(c.rows[i][iptp->second]);
                if(x.ok&&std::isfinite(e)&&std::fabs(x.value)>0.0)
                    ptp.push_back(100.0*std::fabs(e/x.value));
            }

            if(icur!=c.index.end()) {
                const double v=num(c.rows[i][icur->second]);
                if(std::isfinite(v))current.push_back(100.0*v);
            }

            if(ires!=c.index.end()) {
                const double v=num(c.rows[i][ires->second]);
                if(std::isfinite(v))residual.push_back(100.0*v);
            }
        }

        Summary out;
        out.ptp=quantile(ptp,.5);
        out.current=quantile(current,.5);
        out.residual=quantile(residual,.5);
        out.correlated=quantile(corrv,.5);
        return out;
    };

    const Summary a=summarize(
        "10.6 GeV",
        corr10,
        "normed cross sections, ep->epg, exp, 10.6 GeV, unpol",
        "current dependence sys frac, 10.6 GeV"
    );

    const Summary b=summarize(
        "Sp19 Inb",
        corrsp,
        "normed cross sections, ep->epg, exp, Sp19 Inb, unpol",
        "current dependence sys frac, Sp19 Inb"
    );

    const std::array<std::string,5> labels={{
        "Point-to-point",
        "Overall norm.",
        "Current",
        "Period residual",
        "Final corr. scale"
    }};

    const std::array<double,5> va={{
        a.ptp,100.0*norm_frac,a.current,a.residual,a.correlated
    }};
    const std::array<double,5> vb={{
        b.ptp,100.0*norm_frac,b.current,b.residual,b.correlated
    }};

    TCanvas cv("c_systematic_category_summary","",1100,760);
    cv.SetLeftMargin(0.13);
    cv.SetRightMargin(0.04);
    cv.SetBottomMargin(0.18);
    cv.SetTopMargin(0.15);
    cv.SetTicks(1,1);

    const double ymax=1.30*std::max(
        *std::max_element(va.begin(),va.end()),
        *std::max_element(vb.begin(),vb.end())
    );

    TH1D frame("h_systematic_category_summary","",5,0.5,5.5);
    frame.SetMinimum(0.0);
    frame.SetMaximum(std::max(5.5,ymax));
    frame.GetYaxis()->SetTitle("Median relative uncertainty (%)");
    frame.GetYaxis()->SetTitleSize(0.045);
    frame.GetYaxis()->SetLabelSize(0.039);
    frame.GetYaxis()->SetTitleOffset(1.25);
    frame.GetXaxis()->SetLabelSize(0.038);
    for(int i=1;i<=5;++i)
        frame.GetXaxis()->SetBinLabel(i,labels[(size_t)i-1].c_str());
    frame.Draw();

    TGraph ga(5),gb(5);
    for(int i=0;i<5;++i) {
        ga.SetPoint(i,i+1-0.08,va[(size_t)i]);
        gb.SetPoint(i,i+1+0.08,vb[(size_t)i]);
    }

    ga.SetMarkerStyle(20);
    ga.SetMarkerSize(1.25);
    ga.SetMarkerColor(kBlue+1);
    gb.SetMarkerStyle(24);
    gb.SetMarkerSize(1.25);
    gb.SetMarkerColor(kRed+1);
    ga.Draw("P SAME");
    gb.Draw("P SAME");

    TLegend leg(0.62,0.73,0.92,0.84);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.029);
    leg.AddEntry(&ga,"Combined 10.6 GeV","p");
    leg.AddEntry(&gb,"Sp19 Inb","p");
    leg.Draw();

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(0.031);
    t.DrawLatex(0.13,0.955,"Summary of systematic-uncertainty categories");
    t.SetTextSize(0.020);
    t.DrawLatex(
        0.13,0.915,
        "Overall normalization (2.16%) is shown separately and is not included in the correlated scale"
    );

    cv.SaveAs(path.string().c_str());
}

static void write_high_level_summary(
    const fs::path& path,
    const Csv& c,
    const std::vector<double>& corr10,
    const std::vector<double>& corrsp,
    double norm_frac) {

    auto summarize=[&](const std::string& sample,
                       const std::vector<double>& corr,
                       const std::string& xscol,
                       const std::string& currentcol) {
        std::vector<double> cvals,ptp,current,resid;
        const auto ix=c.index.find(xscol);
        const auto iptp=c.index.find("Syst. err (point-to-point total)");
        const auto icur=c.index.find(currentcol);
        const auto ires=c.index.find("run period residual sys frac, "+sample);

        for(size_t i=0;i<c.rows.size();++i){
            if(i<corr.size()&&std::isfinite(corr[i]))cvals.push_back(100*corr[i]);

            if(ix!=c.index.end()&&iptp!=c.index.end()){
                const TupleValue x=tuple_value(c.rows[i][ix->second]);
                const double e=num(c.rows[i][iptp->second]);
                if(x.ok&&std::isfinite(e)&&std::fabs(x.value)>0)
                    ptp.push_back(100*std::fabs(e/x.value));
            }

            if(icur!=c.index.end()){
                const double v=num(c.rows[i][icur->second]);
                if(std::isfinite(v))current.push_back(100*v);
            }
            if(ires!=c.index.end()){
                const double v=num(c.rows[i][ires->second]);
                if(std::isfinite(v))resid.push_back(100*v);
            }
        }

        std::array<double,5> ans={{
            quantile(ptp,.5),
            100*norm_frac,
            quantile(current,.5),
            quantile(resid,.5),
            quantile(cvals,.5)
        }};
        return ans;
    };

    const auto a=summarize(
        "10.6 GeV",corr10,
        "normed cross sections, ep->epg, exp, 10.6 GeV, unpol",
        "current dependence sys frac, 10.6 GeV");
    const auto b=summarize(
        "Sp19 Inb",corrsp,
        "normed cross sections, ep->epg, exp, Sp19 Inb, unpol",
        "current dependence sys frac, Sp19 Inb");

    std::ofstream out(path);
    out<<"sample,median_point_to_point_percent,uncorrelated_normalization_percent,"
          "median_current_percent,median_run_period_residual_percent,"
          "median_correlated_scale_percent\n";
    out<<"10.6 GeV";
    for(double x:a)out<<','<<x;
    out<<'\n';
    out<<"Sp19 Inb";
    for(double x:b)out<<','<<x;
    out<<'\n';
}

} // namespace

bool correlated_scale_systematics(
    const std::string& csv_path,
    const CorrelatedScaleSystematicsOptions& options) {

    try {
        fs::create_directories(options.output_dir);

        Csv c=read_csv(csv_path);

        for(const auto&p:ten6_periods()){
            if(c.index.find(xs_col(p))==c.index.end())
                throw std::runtime_error("missing cross section column for "+p);
            if(c.index.find(current_col(p))==c.index.end())
                throw std::runtime_error(
                    "missing current-systematic column for "+p+
                    "; rerun the updated current systematic first");
        }

        const std::vector<RatioPoint> ratio_points=build_ratio_points(c);
        const std::vector<ThetaResult> reference=build_theta_reference(
            ratio_points,
            options.theta_bin_width_deg,
            options.min_ratio_points_per_period);

        if(reference.empty())
            throw std::runtime_error("no usable theta_p intervals for residual period systematic");

        write_theta_reference(
            fs::path(options.output_dir)/"run_period_residual_theta_p_reference.csv",
            reference);

        draw_theta_decomposition(
            fs::path(options.output_dir)/"run_period_residual_theta_p_decomposition.png",
            reference);

        const int r10=ensure_col(c,"run period residual sys frac, 10.6 GeV");
        const int ar10=ensure_col(c,"Syst.err (run period residual), 10.6 GeV");
        const int cs10=ensure_col(c,"correlated scale sys frac, 10.6 GeV");
        const int acs10=ensure_col(c,"Syst.err (correlated scale), 10.6 GeV");
        const int no10=ensure_col(c,"uncorrelated normalization sys frac, 10.6 GeV");
        const int ano10=ensure_col(c,"Syst.err (uncorrelated normalization), 10.6 GeV");

        const int rsp=ensure_col(c,"run period residual sys frac, Sp19 Inb");
        const int arsp=ensure_col(c,"Syst.err (run period residual), Sp19 Inb");
        const int cssp=ensure_col(c,"correlated scale sys frac, Sp19 Inb");
        const int acssp=ensure_col(c,"Syst.err (correlated scale), Sp19 Inb");
        const int nosp=ensure_col(c,"uncorrelated normalization sys frac, Sp19 Inb");
        const int anosp=ensure_col(c,"Syst.err (uncorrelated normalization), Sp19 Inb");

        const auto ith10=c.index.find("p_theta, 10.6 GeV");
        const auto ithsp=c.index.find("p_theta, Sp19 Inb");
        const auto ic10=c.index.find("current dependence sys frac, 10.6 GeV");
        const auto icsp=c.index.find("current dependence sys frac, Sp19 Inb");
        const auto ix10=c.index.find("normed cross sections, ep->epg, exp, 10.6 GeV, unpol");
        const auto ixsp=c.index.find("normed cross sections, ep->epg, exp, Sp19 Inb, unpol");

        if(ith10==c.index.end()||ithsp==c.index.end()||
           ic10==c.index.end()||icsp==c.index.end())
            throw std::runtime_error("missing theta/current columns required for final scale systematic");

        std::vector<double> corr10(c.rows.size(),std::numeric_limits<double>::quiet_NaN());
        std::vector<double> corrsp(c.rows.size(),std::numeric_limits<double>::quiet_NaN());

        for(size_t i=0;i<c.rows.size();++i){
            const double th10=num(c.rows[i][ith10->second]);
            const double thsp=num(c.rows[i][ithsp->second]);
            const double cur10=num(c.rows[i][ic10->second]);
            const double cursp=num(c.rows[i][icsp->second]);

            if(std::isfinite(th10)&&std::isfinite(cur10)&&cur10>=0.0){
                const double res=interpolate_residual(reference,th10);
                const double corr=std::hypot(cur10,res);
                corr10[i]=corr;
                c.rows[i][r10]=fmt(res);
                c.rows[i][cs10]=fmt(corr);
                c.rows[i][no10]=fmt(options.uncorrelated_normalization_fraction);

                if(ix10!=c.index.end()){
                    const TupleValue x=tuple_value(c.rows[i][ix10->second]);
                    if(x.ok){
                        c.rows[i][ar10]=fmt(std::fabs(x.value)*res);
                        c.rows[i][acs10]=fmt(std::fabs(x.value)*corr);
                        c.rows[i][ano10]=fmt(
                            std::fabs(x.value)*options.uncorrelated_normalization_fraction);
                    }
                }
            }

            if(std::isfinite(thsp)&&std::isfinite(cursp)&&cursp>=0.0){
                const double res=interpolate_residual(reference,thsp);
                const double corr=std::hypot(cursp,res);
                corrsp[i]=corr;
                c.rows[i][rsp]=fmt(res);
                c.rows[i][cssp]=fmt(corr);
                c.rows[i][nosp]=fmt(options.uncorrelated_normalization_fraction);

                if(ixsp!=c.index.end()){
                    const TupleValue x=tuple_value(c.rows[i][ixsp->second]);
                    if(x.ok){
                        c.rows[i][arsp]=fmt(std::fabs(x.value)*res);
                        c.rows[i][acssp]=fmt(std::fabs(x.value)*corr);
                        c.rows[i][anosp]=fmt(
                            std::fabs(x.value)*options.uncorrelated_normalization_fraction);
                    }
                }
            }
        }

        write_csv(csv_path,c);

        // Re-read after writing so the new columns are visible to the summary helper.
        const Csv cfinal=read_csv(csv_path);

        draw_final_kinematic_summary(
            cfinal,corr10,corrsp,
            fs::path(options.output_dir)/"correlated_scale_kinematic_summary.png");

        draw_residual_all_variables(
            fs::path(options.output_dir)/
                "run_period_residual_all_variables_diagnostic.png",
            ratio_points,
            options.min_ratio_points_per_period);

        draw_systematic_category_summary(
            fs::path(options.output_dir)/"systematic_category_summary.png",
            cfinal,corr10,corrsp,
            options.uncorrelated_normalization_fraction);

        write_high_level_summary(
            fs::path(options.output_dir)/"systematic_category_summary.csv",
            cfinal,corr10,corrsp,
            options.uncorrelated_normalization_fraction);

        std::cout
            << "[correlated-scale] Wrote residual period reference and final "
            << "correlated-scale columns to " << csv_path << '\n'
            << "[correlated-scale] Supported theta_p range: "
            << reference.front().center << " -- "
            << reference.back().center << " deg; outside this range the "
            << "nearest supported residual is used.\n";

        return true;
    } catch(const std::exception&e){
        std::cerr<<"[correlated-scale] ERROR: "<<e.what()<<std::endl;
        return false;
    }
}
