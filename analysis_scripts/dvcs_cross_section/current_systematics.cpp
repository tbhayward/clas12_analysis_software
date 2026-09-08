// current_systematics.cpp
// -----------------------------------------------------------------------------
// Correlated current-dependent efficiency systematic.
//
// total_counts.cpp writes the signed one-standard-deviation response of every
// current-corrected DATA/MC yield to each fitted calibration parameter.
// This module propagates those shared nuisance responses through
//
//   S = N(epgamma data) - N(misID MC) * N(ep pi0 data) / N(ep pi0 rec MC)
//   U = S * N(epgamma gen MC) / N(epgamma rec MC)
//
// so U is the acceptance-corrected DVCS yield. Frad, Fbin, luminosity and bin
// volume do not depend on the current-response calibration; therefore the
// fractional response of U is also the fractional current systematic on sigma.
//
// The scalar uncertainty in a bin is the quadrature sum of independent
// one-sigma nuisance responses. Correlation is NOT encoded by that scalar.
// It is retained explicitly in current_systematic_nuisance_responses.csv,
// which stores the signed response vector for every named nuisance.
// -----------------------------------------------------------------------------

#include "current_systematics.h"

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
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
#include <set>
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
    for(char c:s){ if(c=='"') o+="\"\""; else o.push_back(c); }
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
        if(i) out<<',';
        out<<escape_csv(c.header[i]);
    }
    out<<'\n';
    for(const auto&r:c.rows){
        for(size_t i=0;i<c.header.size();++i){
            if(i) out<<',';
            out<<escape_csv(i<r.size()?r[i]:std::string());
        }
        out<<'\n';
    }
    out.close();
    fs::rename(tmp,path);
}

static double num(const std::string& s) {
    if(trim(s).empty()) return std::numeric_limits<double>::quiet_NaN();
    char* e=nullptr;
    const double v=std::strtod(s.c_str(),&e);
    return e==s.c_str()?std::numeric_limits<double>::quiet_NaN():v;
}

static bool tuple_first(const std::string& raw,double& v) {
    std::string s=trim(raw);
    if(s.empty()) return false;
    if(s.front()=='(' && s.back()==')'){
        s=s.substr(1,s.size()-2);
        const auto comma=s.find(',');
        if(comma!=std::string::npos) s=s.substr(0,comma);
    }
    v=num(s);
    return std::isfinite(v);
}

static int ensure_col(Csv& c,const std::string& name) {
    auto it=c.index.find(name);
    if(it!=c.index.end()) return it->second;
    const int i=(int)c.header.size();
    c.header.push_back(name);
    c.index[name]=i;
    for(auto&r:c.rows) r.push_back("");
    return i;
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
    const size_t a=(size_t)std::floor(p),b=(size_t)std::ceil(p);
    const double f=p-a;
    return v[a]*(1-f)+v[b]*f;
}

struct YieldState {
    double nominal=std::numeric_limits<double>::quiet_NaN();
    std::map<std::string,double> delta;
};

using SampleMap=std::map<std::string,YieldState>;
using PeriodMap=std::map<std::string,SampleMap>;
using ResponseByRow=std::map<int,PeriodMap>;

static ResponseByRow read_responses(const std::string& path) {
    Csv c=read_csv(path);
    auto req=[&](const std::string& n)->int{
        auto it=c.index.find(n);
        if(it==c.index.end()) throw std::runtime_error("missing response column: "+n);
        return it->second;
    };
    const int cr=req("row"),cp=req("period"),cs=req("sample"),
              cn=req("nuisance"),cv=req("nominal_yield"),cd=req("delta_yield_1sigma");

    ResponseByRow out;
    for(const auto&r:c.rows){
        const int row=(int)std::llround(num(r[cr]));
        const std::string period=r[cp],sample=r[cs],nuisance=r[cn];
        const double nominal=num(r[cv]),delta=num(r[cd]);
        if(row<0 || period.empty() || sample.empty()) continue;
        YieldState& y=out[row][period][sample];
        if(std::isfinite(nominal)) y.nominal=nominal;
        if(nuisance!="NOMINAL" && !nuisance.empty() && std::isfinite(delta)){
            y.delta[nuisance]+=delta;
        }
    }
    return out;
}

static double read_dmax(const std::string& path) {
    if(!fs::exists(path)) return 0.0;
    Csv c=read_csv(path);
    if(c.header.size()<2) return 0.0;
    for(const auto&r:c.rows){
        if(r.size()>=2 && trim(r[0])=="Dmax"){
            const double d=num(r[1]);
            return std::isfinite(d) && d>=0.0 ? d : 0.0;
        }
    }
    return 0.0;
}

static bool sample_nom(const PeriodMap& p,const std::string& period,
                       const std::string& sample,double& v) {
    auto ip=p.find(period); if(ip==p.end()) return false;
    auto is=ip->second.find(sample); if(is==ip->second.end()) return false;
    v=is->second.nominal;
    return std::isfinite(v);
}

static double sample_delta(const PeriodMap& p,const std::string& period,
                           const std::string& sample,const std::string& nuisance) {
    auto ip=p.find(period); if(ip==p.end()) return 0.0;
    auto is=ip->second.find(sample); if(is==ip->second.end()) return 0.0;
    auto in=is->second.delta.find(nuisance);
    return in==is->second.delta.end()?0.0:in->second;
}

static std::set<std::string> nuisances_for_period(
    const PeriodMap& p,const std::string& period) {
    std::set<std::string> n;
    auto ip=p.find(period); if(ip==p.end()) return n;
    for(const auto&s:ip->second){
        for(const auto&kv:s.second.delta) n.insert(kv.first);
    }
    return n;
}

static bool unfolded_for_period(const PeriodMap& p,
                                const std::string& period,
                                const std::string* nuisance,
                                double sign,
                                double& U) {
    double D,P,Rpi,M,Rg,G;
    if(!sample_nom(p,period,"epg_data",D) ||
       !sample_nom(p,period,"eppi0_data",P) ||
       !sample_nom(p,period,"eppi0_rec",Rpi) ||
       !sample_nom(p,period,"misid_rec",M) ||
       !sample_nom(p,period,"epg_rec",Rg) ||
       !sample_nom(p,period,"epg_gen",G)) return false;

    if(nuisance){
        D   += sign*sample_delta(p,period,"epg_data",*nuisance);
        P   += sign*sample_delta(p,period,"eppi0_data",*nuisance);
        Rpi += sign*sample_delta(p,period,"eppi0_rec",*nuisance);
        M   += sign*sample_delta(p,period,"misid_rec",*nuisance);
        Rg  += sign*sample_delta(p,period,"epg_rec",*nuisance);
    }

    if(!(Rpi>0.0) || !(Rg>0.0) || !(G>0.0)) return false;

    const double S=D-M*P/Rpi;
    if(!(S>0.0) || !std::isfinite(S)) return false;

    U=S*G/Rg;
    return std::isfinite(U) && U>0.0;
}

struct BinResult {
    bool valid=false;
    double scalar=0.0;       // relative
    double data_scalar=0.0;
    double mc_scalar=0.0;
    double transfer_scalar=0.0;
    std::map<std::string,double> signed_response; // (U+ - U-)/(2U0)
};

static BinResult period_result(
    const PeriodMap& p,
    const std::string& period,
    double sp19_dmax) {

    BinResult out;
    double U0=0.0;
    if(!unfolded_for_period(p,period,nullptr,0.0,U0)) return out;

    const auto nuis=nuisances_for_period(p,period);
    double data2=0.0,mc2=0.0;

    for(const auto& n:nuis){
        double up=0.0,dn=0.0;
        if(!unfolded_for_period(p,period,&n,+1.0,up) ||
           !unfolded_for_period(p,period,&n,-1.0,dn)) continue;

        const double fsym=
            (std::fabs(up-U0)+std::fabs(dn-U0))/(2.0*std::fabs(U0));
        const double signed_f=(up-dn)/(2.0*U0);

        out.signed_response[n]=signed_f;

        if(n.rfind("data:",0)==0) data2+=fsym*fsym;
        else if(n.rfind("mc:",0)==0) mc2+=fsym*fsym;
    }

    out.data_scalar=std::sqrt(std::max(0.0,data2));
    out.mc_scalar=std::sqrt(std::max(0.0,mc2));

    if(period=="Sp19 Inb" && sp19_dmax>0.0){
        // The Fa18->Sp19 transfer has no trustworthy low-current Sp19 scan
        // from which to derive an independent calibration.  Use the largest
        // normalized 50-nA shape disparity directly as a conservative Sp19-
        // only correlated scale allowance.  This is deliberately simple:
        // Dmax=0.025 corresponds to a 2.5% transfer uncertainty on every Sp19
        // cross-section bin.  The shared Fa18 calibration nuisance vectors
        // remain separate and keep their genuine Fa18<->Sp19 correlation.
        out.transfer_scalar=sp19_dmax;
        out.signed_response["data:Sp19 transfer allowance:Dmax"]=
            out.transfer_scalar;
    }

    out.scalar=std::sqrt(
        out.data_scalar*out.data_scalar +
        out.mc_scalar*out.mc_scalar +
        out.transfer_scalar*out.transfer_scalar);

    out.valid=std::isfinite(out.scalar);
    return out;
}

static BinResult combined_10p6_result(const PeriodMap& p) {
    static const std::array<std::string,4> periods={
        "Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"
    };

    std::map<std::string,double> U0;
    double total0=0.0;
    std::set<std::string> alln;

    for(const auto& per:periods){
        double u=0.0;
        if(!unfolded_for_period(p,per,nullptr,0.0,u)) continue;
        U0[per]=u;
        total0+=u;
        const auto n=nuisances_for_period(p,per);
        alln.insert(n.begin(),n.end());
    }

    BinResult out;
    if(!(total0>0.0)) return out;

    double data2=0.0,mc2=0.0;

    for(const auto& n:alln){
        double plus=0.0,minus=0.0;
        bool ok=true;

        for(const auto& kv:U0){
            double up=0.0,dn=0.0;
            if(!unfolded_for_period(p,kv.first,&n,+1.0,up) ||
               !unfolded_for_period(p,kv.first,&n,-1.0,dn)){
                // If this nuisance does not occur in this period, nominal is
                // the correct contribution.
                const double dtest =
                    sample_delta(p,kv.first,"epg_data",n)+
                    sample_delta(p,kv.first,"eppi0_data",n)+
                    sample_delta(p,kv.first,"eppi0_rec",n)+
                    sample_delta(p,kv.first,"misid_rec",n)+
                    sample_delta(p,kv.first,"epg_rec",n);
                if(dtest==0.0){ up=dn=kv.second; }
                else { ok=false; break; }
            }
            plus+=up;
            minus+=dn;
        }

        if(!ok) continue;

        const double fsym=
            (std::fabs(plus-total0)+std::fabs(minus-total0))/
            (2.0*std::fabs(total0));
        const double signed_f=(plus-minus)/(2.0*total0);
        out.signed_response[n]=signed_f;

        if(n.rfind("data:",0)==0) data2+=fsym*fsym;
        else if(n.rfind("mc:",0)==0) mc2+=fsym*fsym;
    }

    out.data_scalar=std::sqrt(std::max(0.0,data2));
    out.mc_scalar=std::sqrt(std::max(0.0,mc2));
    out.scalar=std::sqrt(data2+mc2);
    out.valid=true;
    return out;
}

struct AxisSpec {
    std::string lo,hi,title;
};

static void make_kinematic_summary(
    const Csv& c,
    const std::vector<BinResult>& results,
    const std::string& label,
    const std::string& token,
    const std::string& outdir) {

    const std::array<AxisSpec,4> specs={{
        {"xBmin","xBmax","x_{B}"},
        {"Q2min","Q2max","Q^{2} (GeV^{2})"},
        {"t_abs_min","t_abs_max","|t| (GeV^{2})"},
        {"phimin","phimax","#phi (deg)"}
    }};

    TCanvas cv(("c_current_"+token).c_str(),"",1450,1050);
    cv.Divide(2,2,0.002,0.002);

    for(int ia=0;ia<4;++ia){
        cv.cd(ia+1);
        gPad->SetLeftMargin((ia%2==0)?0.14:0.12);
        gPad->SetRightMargin(0.035);
        gPad->SetBottomMargin((ia>=2)?0.15:0.12);
        gPad->SetTopMargin((ia<2)?0.16:0.10);
        gPad->SetTicks(1,1);

        auto ilo=c.index.find(specs[ia].lo),ihi=c.index.find(specs[ia].hi);
        if(ilo==c.index.end()||ihi==c.index.end()) continue;

        struct Bucket{double lo=0,hi=0;std::vector<double> v;};
        std::map<std::pair<double,double>,Bucket> buckets;

        for(size_t r=0;r<c.rows.size() && r<results.size();++r){
            if(!results[r].valid) continue;
            double lo=num(c.rows[r][ilo->second]);
            double hi=num(c.rows[r][ihi->second]);
            if(!std::isfinite(lo)||!std::isfinite(hi)||!(hi>lo)) continue;
            if(ia==3){
                double phi=0.5*(lo+hi);
                while(phi<0.0)phi+=360.0;
                while(phi>=360.0)phi-=360.0;
                const int ib=std::min(11,std::max(0,int(phi/30.0)));
                lo=30.0*ib; hi=lo+30.0;
            }
            auto&b=buckets[{lo,hi}];b.lo=lo;b.hi=hi;
            b.v.push_back(100.0*results[r].scalar);
        }

        TGraphAsymmErrors g;
        int ip=0;
        double xmin=INFINITY,xmax=-INFINITY,ymax=0.0;
        for(const auto&kv:buckets){
            const auto&b=kv.second;
            if(b.v.empty()) continue;
            const double x=.5*(b.lo+b.hi),ex=.5*(b.hi-b.lo);
            const double md=quantile(b.v,.5),q16=quantile(b.v,.16),q84=quantile(b.v,.84);
            g.SetPoint(ip,x,md);
            g.SetPointError(ip,ex,ex,md-q16,q84-md);
            xmin=std::min(xmin,b.lo);xmax=std::max(xmax,b.hi);ymax=std::max(ymax,q84);
            ++ip;
        }
        if(!(xmax>xmin)){xmin=0;xmax=1;}

        TH1D frame(("h_current_"+token+"_"+std::to_string(ia)).c_str(),"",100,xmin,xmax);
        frame.SetMinimum(0.0);
        frame.SetMaximum(std::max(2.0,1.20*ymax));
        frame.GetXaxis()->SetTitle(specs[ia].title.c_str());
        frame.GetYaxis()->SetTitle("Current-response systematic (%)");
        frame.GetXaxis()->SetTitleSize(.050);
        frame.GetYaxis()->SetTitleSize(.047);
        frame.GetXaxis()->SetLabelSize(.043);
        frame.GetYaxis()->SetLabelSize(.041);
        frame.GetYaxis()->SetTitleOffset(1.20);
        frame.DrawCopy();

        g.SetMarkerStyle(20);
        g.SetMarkerSize(.95);
        g.SetMarkerColor(kBlue+1);
        g.SetLineColor(kBlue+1);
        g.SetLineWidth(2);
        g.DrawClone("PZ SAME");

        TLatex p;
        p.SetNDC();p.SetTextFont(42);p.SetTextSize(.036);
        const std::string lab=std::string("(")+char('a'+ia)+")";
        p.DrawLatex(.18,.80,lab.c_str());
    }

    cv.cd(0);
    TLatex t;
    t.SetNDC();t.SetTextFont(42);t.SetTextAlign(22);
    t.SetTextSize(.021);
    const std::string title="Current-dependent efficiency systematic: "+label;
    t.DrawLatex(.50,.975,title.c_str());
    t.SetTextSize(.015);
    t.DrawLatex(.50,.946,
                "Median and central 68% range; #phi projection grouped in 30^{#circ} intervals");

    cv.SaveAs(
        (fs::path(outdir)/("current_systematic_"+token+"_kinematic_summary.png"))
            .string().c_str());
}

} // namespace

bool evaluate_current_dependence_systematics(
    const std::string& analysis_csv,
    const CurrentSystematicsOptions& options) {

    try {
        fs::create_directories(options.output_dir);

        Csv csv=read_csv(analysis_csv);
        const ResponseByRow resp=read_responses(options.nuisance_response_csv);
        const double dmax=read_dmax(options.transfer_shape_csv);

        std::vector<BinResult> r10(csv.rows.size());
        std::vector<BinResult> rfa(csv.rows.size());
        std::vector<BinResult> rfaout(csv.rows.size());
        std::vector<BinResult> rsp18in(csv.rows.size());
        std::vector<BinResult> rsp18out(csv.rows.size());
        std::vector<BinResult> rsp(csv.rows.size());

        std::ofstream response_out(
            fs::path(options.output_dir)/
            "current_systematic_nuisance_responses.csv");
        response_out
            <<"row,target,nuisance,signed_fractional_response\n";

        for(size_t i=0;i<csv.rows.size();++i){
            auto it=resp.find((int)i);
            if(it==resp.end()) continue;

            r10[i]=combined_10p6_result(it->second);
            rfa[i]=period_result(it->second,"Fa18 Inb",0.0);
            rfaout[i]=period_result(it->second,"Fa18 Out",0.0);
            rsp18in[i]=period_result(it->second,"Sp18 Inb",0.0);
            rsp18out[i]=period_result(it->second,"Sp18 Out",0.0);
            rsp[i]=period_result(it->second,"Sp19 Inb",dmax);

            const std::array<std::pair<std::string,const BinResult*>,6> targets={{
                {"10.6 GeV",&r10[i]},
                {"Fa18 Inb",&rfa[i]},
                {"Fa18 Out",&rfaout[i]},
                {"Sp18 Inb",&rsp18in[i]},
                {"Sp18 Out",&rsp18out[i]},
                {"Sp19 Inb",&rsp[i]}
            }};
            for(const auto&t:targets){
                if(!t.second->valid) continue;
                for(const auto&n:t.second->signed_response){
                    response_out<<i<<','
                                <<escape_csv(t.first)<<','
                                <<escape_csv(n.first)<<','
                                <<std::setprecision(12)<<n.second<<'\n';
                }
            }
        }

        if(options.write_csv_columns){
            struct TargetCols{
                std::string label;
                std::string token;
                std::string xs_col;
                std::vector<BinResult>* result;
            };
            std::array<TargetCols,6> targets={{
                {"10.6 GeV","10.6 GeV",
                 "cross sections, ep->epg, exp, 10.6 GeV, unpol",&r10},
                {"Fa18 Inb","Fa18 Inb",
                 "cross sections, ep->epg, exp, Fa18 Inb, unpol",&rfa},
                {"Fa18 Out","Fa18 Out",
                 "cross sections, ep->epg, exp, Fa18 Out, unpol",&rfaout},
                {"Sp18 Inb","Sp18 Inb",
                 "cross sections, ep->epg, exp, Sp18 Inb, unpol",&rsp18in},
                {"Sp18 Out","Sp18 Out",
                 "cross sections, ep->epg, exp, Sp18 Out, unpol",&rsp18out},
                {"Sp19 Inb","Sp19 Inb",
                 "cross sections, ep->epg, exp, Sp19 Inb, unpol",&rsp}
            }};

            for(auto&t:targets){
                const int cf=ensure_col(
                    csv,"current dependence sys frac, "+t.token);
                const int ca=ensure_col(
                    csv,"Syst.err (current dependence), "+t.token);
                const int cd=ensure_col(
                    csv,"current dependence DATA component frac, "+t.token);
                const int cm=ensure_col(
                    csv,"current dependence MC component frac, "+t.token);
                const int ct=ensure_col(
                    csv,"current dependence transfer component frac, "+t.token);

                auto ix=csv.index.find(t.xs_col);
                for(size_t r=0;r<csv.rows.size();++r){
                    const BinResult& br=(*t.result)[r];
                    if(!br.valid) continue;
                    csv.rows[r][cf]=fmt(br.scalar);
                    csv.rows[r][cd]=fmt(br.data_scalar);
                    csv.rows[r][cm]=fmt(br.mc_scalar);
                    csv.rows[r][ct]=fmt(br.transfer_scalar);

                    if(ix!=csv.index.end()){
                        double xs=0.0;
                        if(tuple_first(csv.rows[r][ix->second],xs)){
                            csv.rows[r][ca]=fmt(std::fabs(xs)*br.scalar);
                        }
                    }
                }
            }

            write_csv(analysis_csv,csv);
        }

        // Numerical summary.
        {
            std::ofstream out(
                fs::path(options.output_dir)/"current_systematic_summary.csv");
            out<<"sample,n,mean_percent,median_percent,p16_percent,p84_percent,p95_percent,"
                  "median_data_component_percent,median_mc_component_percent,"
                  "median_transfer_component_percent,Sp19_Dmax\n";

            auto emit=[&](const std::string& name,const std::vector<BinResult>& v){
                std::vector<double> all,data,mc,tr;
                for(const auto&r:v){
                    if(!r.valid) continue;
                    all.push_back(100*r.scalar);
                    data.push_back(100*r.data_scalar);
                    mc.push_back(100*r.mc_scalar);
                    tr.push_back(100*r.transfer_scalar);
                }
                if(all.empty()) return;
                const double mean=std::accumulate(all.begin(),all.end(),0.0)/all.size();
                out<<name<<','<<all.size()<<','<<mean<<','
                   <<quantile(all,.5)<<','<<quantile(all,.16)<<','
                   <<quantile(all,.84)<<','<<quantile(all,.95)<<','
                   <<quantile(data,.5)<<','<<quantile(mc,.5)<<','
                   <<quantile(tr,.5)<<','<<dmax<<'\n';
            };
            emit("10.6 GeV",r10);
            emit("Fa18 Inb",rfa);
            emit("Fa18 Out",rfaout);
            emit("Sp18 Inb",rsp18in);
            emit("Sp18 Out",rsp18out);
            emit("Sp19 Inb",rsp);
        }

        // Distribution overlay.
        {
            std::vector<double>a,b,c;
            for(size_t i=0;i<csv.rows.size();++i){
                if(r10[i].valid)a.push_back(100*r10[i].scalar);
                if(rfa[i].valid)b.push_back(100*rfa[i].scalar);
                if(rsp[i].valid)c.push_back(100*rsp[i].scalar);
            }
            std::vector<double> all=a;all.insert(all.end(),b.begin(),b.end());all.insert(all.end(),c.begin(),c.end());
            const double xmax=std::max(3.0,std::min(30.0,std::ceil(quantile(all,.99)+1.0)));

            TH1D h10("h_current_dist_10","",40,0,xmax);
            TH1D hfa("h_current_dist_fa","",40,0,xmax);
            TH1D hsp("h_current_dist_sp","",40,0,xmax);
            for(double x:a)if(x<xmax)h10.Fill(x);
            for(double x:b)if(x<xmax)hfa.Fill(x);
            for(double x:c)if(x<xmax)hsp.Fill(x);
            for(TH1D*h:{&h10,&hfa,&hsp})if(h->Integral()>0)h->Scale(1.0/h->Integral());

            TCanvas cv("c_current_dist","",1080,760);
            cv.SetLeftMargin(.13);cv.SetRightMargin(.035);cv.SetBottomMargin(.13);cv.SetTopMargin(.13);cv.SetTicks(1,1);
            h10.SetLineColor(kBlack);h10.SetLineWidth(3);
            hfa.SetLineColor(kBlue+1);hfa.SetLineWidth(3);
            hsp.SetLineColor(kRed+1);hsp.SetLineWidth(3);
            hsp.SetLineStyle(2);
            h10.SetMaximum(1.25*std::max({h10.GetMaximum(),hfa.GetMaximum(),hsp.GetMaximum()}));
            h10.GetXaxis()->SetTitle("Relative current-response systematic (%)");
            h10.GetYaxis()->SetTitle("Fraction of populated bins");
            h10.GetXaxis()->SetTitleSize(.045);h10.GetYaxis()->SetTitleSize(.043);
            h10.GetXaxis()->SetLabelSize(.039);h10.GetYaxis()->SetLabelSize(.038);
            h10.GetYaxis()->SetTitleOffset(1.25);
            h10.Draw("HIST");hfa.Draw("HIST SAME");hsp.Draw("HIST SAME");

            TLegend l(.58,.68,.92,.84);
            l.SetBorderSize(0);l.SetFillStyle(0);l.SetTextFont(42);l.SetTextSize(.030);
            l.AddEntry(&h10,"Combined 10.6 GeV","l");
            l.AddEntry(&hfa,"Fa18 Inb","l");
            l.AddEntry(&hsp,"Sp19 Inb","l");
            l.Draw();

            TLatex tt;tt.SetNDC();tt.SetTextFont(42);tt.SetTextSize(.031);
            tt.DrawLatex(.13,.955,"Current-dependent efficiency systematic");

            cv.SaveAs(
                (fs::path(options.output_dir)/"current_systematic_distribution.png")
                    .string().c_str());
        }

        // Diagnose whether each named correlated nuisance acts mostly as a
        // scale shift or changes the kinematic shape.  This does not replace
        // the full signed response vectors used to reconstruct the covariance.
        {
            std::ofstream out(fs::path(options.output_dir)/"current_nuisance_scale_shape_summary.csv");
            out<<"target,nuisance,n_bins,mean_signed_percent,rms_shape_percent,rms_total_percent\n";
            const std::array<std::pair<std::string,const std::vector<BinResult>*>,3> targets={{
                {"10.6 GeV",&r10},{"Fa18 Inb",&rfa},{"Sp19 Inb",&rsp}
            }};
            for(const auto&t:targets){
                std::map<std::string,std::vector<double>> by;
                for(const auto&r:*t.second){
                    if(!r.valid)continue;
                    for(const auto&kv:r.signed_response)by[kv.first].push_back(100.0*kv.second);
                }
                for(const auto&kv:by){
                    if(kv.second.empty())continue;
                    const double mu=std::accumulate(kv.second.begin(),kv.second.end(),0.0)/kv.second.size();
                    double ss=0.0,tt=0.0;
                    for(double x:kv.second){ss+=(x-mu)*(x-mu);tt+=x*x;}
                    out<<escape_csv(t.first)<<','<<escape_csv(kv.first)<<','<<kv.second.size()<<','
                       <<mu<<','<<std::sqrt(ss/kv.second.size())<<','<<std::sqrt(tt/kv.second.size())<<'\n';
                }
            }
        }

        // Visualize how each named correlated nuisance decomposes into a\n        // scale-like mean response and a residual kinematic-shape response.\n        // This is a diagnostic representation only; the production covariance\n        // continues to use the complete signed response vector bin-by-bin.\n        {\n            struct Point {\n                std::string nuisance;\n                double scale=0.0;\n                double shape=0.0;\n                double total=0.0;\n            };\n\n            const std::array<std::pair<std::string,const std::vector<BinResult>*>,3> targets={{\n                {"Combined 10.6 GeV",&r10},{"Fa18 Inb",&rfa},{"Sp19 Inb",&rsp}\n            }};\n\n            TCanvas cv("c_current_scale_shape","",1450,500);\n            cv.Divide(3,1,0.002,0.002);\n\n            for(int ip=0; ip<3; ++ip){\n                std::map<std::string,std::vector<double>> by;\n                for(const auto&r:*targets[ip].second){\n                    if(!r.valid) continue;\n                    for(const auto&kv:r.signed_response)\n                        by[kv.first].push_back(100.0*kv.second);\n                }\n\n                std::vector<Point> points;\n                for(const auto&kv:by){\n                    if(kv.second.empty()) continue;\n                    const double mu=std::accumulate(kv.second.begin(),kv.second.end(),0.0)/kv.second.size();\n                    double ss=0.0,tt=0.0;\n                    for(double x:kv.second){\n                        ss+=(x-mu)*(x-mu);\n                        tt+=x*x;\n                    }\n                    Point p;\n                    p.nuisance=kv.first;\n                    p.scale=std::fabs(mu);\n                    p.shape=std::sqrt(ss/kv.second.size());\n                    p.total=std::sqrt(tt/kv.second.size());\n                    points.push_back(p);\n                }\n\n                cv.cd(ip+1);\n                gPad->SetLeftMargin(.15);gPad->SetRightMargin(.045);\n                // Reserve an explicit title band above the frame.  Keep
                // period labels in that band rather than sitting on the
                // upper frame line.
                gPad->SetBottomMargin(.16);
                gPad->SetTopMargin(.19);
                gPad->SetTicks(1,1);\n\n                double xmax=0.0,ymax=0.0;\n                for(const auto&p:points){xmax=std::max(xmax,p.scale);ymax=std::max(ymax,p.shape);}\n                xmax=std::max(0.5,1.25*xmax);\n                ymax=std::max(0.5,1.25*ymax);\n\n                TH1D frame(("h_current_scale_shape_"+std::to_string(ip)).c_str(),"",100,0,xmax);\n                frame.SetMinimum(0.0);frame.SetMaximum(ymax);\n                frame.GetXaxis()->SetTitle("Scale-like response |mean| (%)");\n                frame.GetYaxis()->SetTitle("Residual shape RMS (%)");\n                frame.GetXaxis()->SetTitleSize(.046);frame.GetYaxis()->SetTitleSize(.044);\n                frame.GetXaxis()->SetLabelSize(.040);frame.GetYaxis()->SetLabelSize(.039);\n                frame.GetYaxis()->SetTitleOffset(1.30);\n                frame.DrawCopy();\n\n                TGraph gd,gm,gt;\n                int nd=0,nm=0,nt=0;\n                for(const auto&p:points){\n                    if(p.nuisance.rfind("data:Sp19 transfer",0)==0){\n                        gt.SetPoint(nt++,p.scale,p.shape);\n                    } else if(p.nuisance.rfind("data:",0)==0){\n                        gd.SetPoint(nd++,p.scale,p.shape);\n                    } else if(p.nuisance.rfind("mc:",0)==0){\n                        gm.SetPoint(nm++,p.scale,p.shape);\n                    }\n                }\n                gd.SetMarkerStyle(20);gd.SetMarkerSize(1.0);gd.SetMarkerColor(kBlue+1);\n                gm.SetMarkerStyle(22);gm.SetMarkerSize(1.0);gm.SetMarkerColor(kGreen+2);\n                gt.SetMarkerStyle(24);gt.SetMarkerSize(1.2);gt.SetMarkerColor(kRed+1);\n                if(nd)gd.DrawClone("P SAME");\n                if(nm)gm.DrawClone("P SAME");\n                if(nt)gt.DrawClone("P SAME");\n\n                TLatex lab;lab.SetNDC();lab.SetTextFont(42);lab.SetTextSize(.037);\n                lab.DrawLatex(.18,.85,targets[ip].first.c_str());\n\n                if(ip==0){\n                    TLegend l(.48,.60,.91,.79);\n                    l.SetBorderSize(0);l.SetFillStyle(0);l.SetTextFont(42);l.SetTextSize(.026);\n                    l.AddEntry(&gd,"DATA calibration nuisances","p");\n                    l.AddEntry(&gm,"MC calibration nuisances","p");\n                    l.DrawClone();\n                }\n                if(ip==2 && nt){\n                    TLegend l(.43,.54,.91,.79);\n                    l.SetBorderSize(0);l.SetFillStyle(0);l.SetTextFont(42);l.SetTextSize(.025);\n                    l.AddEntry(&gd,"DATA calibration nuisances","p");\n                    l.AddEntry(&gm,"MC calibration nuisances","p");\n                    l.AddEntry(&gt,"Sp19 transfer nuisance","p");\n                    l.DrawClone();\n                }\n\n                // Label only the three largest nuisance directions to keep the\n                // plot readable.  Use a compact detector-region suffix.\n                std::sort(points.begin(),points.end(),[](const Point&a,const Point&b){return a.total>b.total;});\n                const int nlabel=std::min<int>(3,points.size());\n                for(int il=0;il<nlabel;++il){\n                    std::string name=points[il].nuisance;\n                    const auto pos=name.find_last_of(':');\n                    if(pos!=std::string::npos) name=name.substr(pos+1);\n                    TLatex tx;tx.SetTextFont(42);tx.SetTextSize(.025);\n                    tx.DrawLatex(points[il].scale+0.015*xmax,\n                                 points[il].shape+0.020*ymax,\n                                 name.c_str());\n                }\n            }\n\n            cv.cd(0);\n            TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextAlign(22);t.SetTextSize(.021);\n            t.DrawLatex(.50,.976,"Scale-like and kinematic-shape content of current-response nuisances");\n            cv.SaveAs((fs::path(options.output_dir)/"current_nuisance_scale_shape.png").string().c_str());\n        }\n\n        // Show the DATA-calibration, MC-calibration, and Sp19 transfer pieces
        // separately.  The scalar total shown elsewhere is their quadrature.
        {
            TCanvas cv("c_current_components","",1450,470);
            cv.Divide(3,1,0.002,0.002);
            const std::array<std::pair<std::string,const std::vector<BinResult>*>,3> targets={{
                {"Combined 10.6 GeV",&r10},{"Fa18 Inb",&rfa},{"Sp19 Inb",&rsp}
            }};
            for(int ip=0;ip<3;++ip){
                cv.cd(ip+1);
                gPad->SetLeftMargin(.16);
                gPad->SetRightMargin(.035);
                gPad->SetBottomMargin(.16);
                // Leave a dedicated band above the plotting frame for the
                // period label so glyphs cannot be clipped by the frame.
                gPad->SetTopMargin(.19);
                gPad->SetTicks(1,1);
                std::vector<double>d,m,tr;
                for(const auto&r:*targets[ip].second){
                    if(!r.valid)continue;
                    d.push_back(100*r.data_scalar);m.push_back(100*r.mc_scalar);
                    if(r.transfer_scalar>0)tr.push_back(100*r.transfer_scalar);
                }
                std::vector<double> all=d;all.insert(all.end(),m.begin(),m.end());all.insert(all.end(),tr.begin(),tr.end());
                const double xmax=all.empty()?5.0:std::max(2.0,std::min(20.0,std::ceil(quantile(all,.99)+.5)));
                TH1D hd(("h_cur_comp_d_"+std::to_string(ip)).c_str(),"",36,0,xmax);
                TH1D hm(("h_cur_comp_m_"+std::to_string(ip)).c_str(),"",36,0,xmax);
                TH1D ht(("h_cur_comp_t_"+std::to_string(ip)).c_str(),"",36,0,xmax);
                for(double x:d)if(x<xmax)hd.Fill(x);for(double x:m)if(x<xmax)hm.Fill(x);for(double x:tr)if(x<xmax)ht.Fill(x);
                for(TH1D*h:{&hd,&hm,&ht})if(h->Integral()>0)h->Scale(1.0/h->Integral());
                hd.SetLineColor(kBlue+1);hd.SetLineWidth(3);hm.SetLineColor(kGreen+2);hm.SetLineWidth(3);ht.SetLineColor(kRed+1);ht.SetLineWidth(3);ht.SetLineStyle(2);
                hd.SetMaximum(1.25*std::max({hd.GetMaximum(),hm.GetMaximum(),ht.GetMaximum(),0.01}));
                hd.GetXaxis()->SetTitle("Relative component uncertainty (%)");hd.GetYaxis()->SetTitle("Fraction of bins");
                hd.GetXaxis()->SetTitleSize(.050);hd.GetYaxis()->SetTitleSize(.047);hd.GetXaxis()->SetLabelSize(.043);hd.GetYaxis()->SetLabelSize(.041);hd.GetYaxis()->SetTitleOffset(1.35);
                // These histograms are local to this loop iteration while the
                // canvas is saved after all three pads are filled.  Draw pad-
                // owned copies so the histograms survive until SaveAs().
                hd.DrawCopy("HIST");
                hm.DrawCopy("HIST SAME");
                if(ht.Integral()>0) ht.DrawCopy("HIST SAME");

                TLatex lab;
                lab.SetNDC();
                lab.SetTextFont(42);
                lab.SetTextSize(.035);
                lab.DrawLatex(.18,.865,targets[ip].first.c_str());

                if(ip==0){
                    TLegend l(.48,.60,.90,.78);
                    l.SetBorderSize(0);l.SetFillStyle(0);l.SetTextFont(42);l.SetTextSize(.026);
                    l.AddEntry(&hd,"DATA response calibration","l");
                    l.AddEntry(&hm,"MC response calibration","l");
                    l.DrawClone();
                }
                if(ip==2){
                    TLegend l(.43,.54,.90,.78);
                    l.SetBorderSize(0);l.SetFillStyle(0);l.SetTextFont(42);l.SetTextSize(.025);
                    l.AddEntry(&hd,"DATA response calibration","l");
                    l.AddEntry(&hm,"MC response calibration","l");
                    l.AddEntry(&ht,"Fa18 #rightarrow Sp19 transfer","l");
                    l.DrawClone();
                }
            }
            cv.cd(0);
            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextAlign(22);
            t.SetTextSize(.020);
            t.DrawLatex(.50,.965,
                        "Components of the current-dependent efficiency systematic");
            const std::string component_path =
                (fs::path(options.output_dir) /
                 "current_systematic_component_distributions.png").string();
            cv.SaveAs(component_path.c_str());
            std::cout
                << "[current-systematics] Wrote component-distribution diagnostic: "
                << component_path << std::endl;
        }

        make_kinematic_summary(csv,r10,"combined 10.6 GeV","10p6",options.output_dir);
        make_kinematic_summary(csv,rfa,"Fa18 Inb","fa18_inb",options.output_dir);
        make_kinematic_summary(csv,rsp,"Sp19 Inb","sp19_inb",options.output_dir);

        std::cout
            << "[current-systematics] Completed. Sp19 transfer Dmax="
            << dmax << ".\n"
            << "[current-systematics] Correlations are retained in "
            << (fs::path(options.output_dir)/
                "current_systematic_nuisance_responses.csv").string()
            << std::endl;

        return true;
    } catch(const std::exception&e){
        std::cerr<<"[current-systematics] ERROR: "<<e.what()<<std::endl;
        return false;
    }
}
