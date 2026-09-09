
#include "pi0_systematics.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <cstdio>
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

namespace {

namespace fs = std::filesystem;

constexpr double kPi0BackgroundRelativeUncertainty = 0.072;

struct Csv {
    std::vector<std::string> header;
    std::unordered_map<std::string,int> index;
    std::vector<std::vector<std::string>> rows;
};

struct Triple {
    double value = std::numeric_limits<double>::quiet_NaN();
    double stat  = std::numeric_limits<double>::quiet_NaN();
    double sys   = std::numeric_limits<double>::quiet_NaN();
    bool ok = false;
};

struct Summary {
    double median = std::numeric_limits<double>::quiet_NaN();
    double mean = std::numeric_limits<double>::quiet_NaN();
    double q16 = std::numeric_limits<double>::quiet_NaN();
    double q84 = std::numeric_limits<double>::quiet_NaN();
    double q95 = std::numeric_limits<double>::quiet_NaN();
    size_t n = 0;
};

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted=false;

    for(size_t i=0;i<line.size();++i) {
        const char c=line[i];
        if(c=='"') {
            if(quoted && i+1<line.size() && line[i+1]=='"') {
                cur.push_back('"');
                ++i;
            } else {
                quoted=!quoted;
            }
        } else if(c==',' && !quoted) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

static std::string csv_escape(const std::string& s) {
    if(s.find_first_of(",\"\n\r")==std::string::npos) return s;
    std::string out="\"";
    for(char c:s) {
        if(c=='"') out+="\"\"";
        else out.push_back(c);
    }
    out.push_back('"');
    return out;
}

static Csv read_csv(const std::string& path) {
    std::ifstream in(path);
    if(!in) throw std::runtime_error("cannot open CSV: "+path);

    Csv c;
    std::string line;
    if(!std::getline(in,line)) throw std::runtime_error("empty CSV: "+path);
    c.header=split_csv_line(line);

    for(int i=0;i<(int)c.header.size();++i) c.index[c.header[(size_t)i]]=i;

    while(std::getline(in,line)) {
        auto row=split_csv_line(line);
        if(row.size()<c.header.size()) row.resize(c.header.size());
        if(row.size()!=c.header.size())
            throw std::runtime_error("CSV row-width mismatch in "+path);
        c.rows.push_back(std::move(row));
    }
    return c;
}

static void write_csv_atomic(const std::string& path,const Csv& c) {
    const std::string tmp=path+".tmp_pi0_systematics";
    {
        std::ofstream out(tmp);
        if(!out) throw std::runtime_error("cannot write temporary CSV: "+tmp);
        for(size_t i=0;i<c.header.size();++i) {
            if(i) out<<',';
            out<<csv_escape(c.header[i]);
        }
        out<<'\n';
        for(const auto& row:c.rows) {
            for(size_t i=0;i<c.header.size();++i) {
                if(i) out<<',';
                out<<csv_escape(row[i]);
            }
            out<<'\n';
        }
    }
    std::remove(path.c_str());
    if(std::rename(tmp.c_str(),path.c_str())!=0)
        throw std::runtime_error("atomic rename failed for "+path);
}

static int ensure_column(Csv& c,const std::string& name) {
    auto it=c.index.find(name);
    if(it!=c.index.end()) return it->second;
    const int idx=(int)c.header.size();
    c.index[name]=idx;
    c.header.push_back(name);
    for(auto& row:c.rows) row.push_back("");
    return idx;
}

static int require_column(const Csv& c,const std::string& name) {
    auto it=c.index.find(name);
    if(it==c.index.end()) throw std::runtime_error("missing column: "+name);
    return it->second;
}

static double number(const std::string& s) {
    if(s.empty()) return std::numeric_limits<double>::quiet_NaN();
    char* end=nullptr;
    const double v=std::strtod(s.c_str(),&end);
    if(end==s.c_str()) return std::numeric_limits<double>::quiet_NaN();
    return v;
}

static Triple parse_triple(const std::string& s) {
    Triple out;
    std::string x=s;
    x.erase(std::remove_if(x.begin(),x.end(),[](unsigned char c){return std::isspace(c);}),x.end());
    if(x.empty()) return out;
    if(x.front()=='(' && x.back()==')') x=x.substr(1,x.size()-2);

    std::vector<double> vals;
    std::stringstream ss(x);
    std::string tok;
    while(std::getline(ss,tok,',')) {
        char* end=nullptr;
        const double v=std::strtod(tok.c_str(),&end);
        if(end==tok.c_str()) return out;
        vals.push_back(v);
    }
    if(vals.empty()) return out;
    out.value=vals[0];
    out.stat=(vals.size()>1)?vals[1]:0.0;
    out.sys=(vals.size()>2)?vals[2]:0.0;
    out.ok=std::isfinite(out.value);
    return out;
}

static std::string fmt(double x) {
    if(!std::isfinite(x)) return "";
    std::ostringstream os;
    os<<std::setprecision(12)<<x;
    return os.str();
}

static double quantile(std::vector<double> v,double q) {
    v.erase(std::remove_if(v.begin(),v.end(),[](double x){return !std::isfinite(x);}),v.end());
    if(v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(),v.end());
    if(q<=0.0) return v.front();
    if(q>=1.0) return v.back();
    const double pos=q*(v.size()-1);
    const size_t lo=(size_t)std::floor(pos);
    const size_t hi=(size_t)std::ceil(pos);
    const double f=pos-lo;
    return v[lo]*(1.0-f)+v[hi]*f;
}

static Summary summarize(const std::vector<double>& vin) {
    std::vector<double> v;
    for(double x:vin) if(std::isfinite(x)) v.push_back(x);
    Summary s;
    s.n=v.size();
    if(v.empty()) return s;
    s.mean=std::accumulate(v.begin(),v.end(),0.0)/v.size();
    s.median=quantile(v,.50);
    s.q16=quantile(v,.16);
    s.q84=quantile(v,.84);
    s.q95=quantile(v,.95);
    return s;
}

static std::vector<std::string> ten6_periods() {
    return {"Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"};
}

static std::string contam_col(const std::string& period) {
    return "contamination ratio, "+period;
}

static std::string unfolded_col(const std::string& period) {
    return "acceptance corrected yield, ep->epg, exp, "+period+", unpol";
}

static std::string xs_col(const std::string& sample) {
    return "normed cross sections, ep->epg, exp, "+sample+", unpol";
}

static bool period_background_to_signal(
    const Csv& c,
    const std::vector<std::string>& row,
    const std::string& period,
    double& signal_corr,
    double& background_corr) {

    const Triple u=parse_triple(row[(size_t)require_column(c,unfolded_col(period))]);
    const Triple ct=parse_triple(row[(size_t)require_column(c,contam_col(period))]);

    if(!u.ok || !ct.ok || !(u.value>0.0) ||
       !std::isfinite(ct.value) || ct.value<0.0 || ct.value>=0.95) {
        return false;
    }

    signal_corr=u.value;
    background_corr=u.value*ct.value/(1.0-ct.value);
    return std::isfinite(background_corr) && background_corr>=0.0;
}

static double recompute_point_to_point_total(const Csv& c,const std::vector<std::string>& row) {
    static const std::array<const char*,6> cols={{
        "Syst. err (pi0 subtraction)",
        "Syst. err (Acceptance)",
        "Syst.err (Frad)",
        "Syst.err (Fbin)",
        "Syst. err (exclusivity cuts)",
        "Syst. err (fiducial cuts)"
    }};
    double sum2=0.0;
    for(const char* name:cols) {
        const int j=require_column(c,name);
        const double v=number(row[(size_t)j]);
        if(!std::isfinite(v) || v<0.0) return std::numeric_limits<double>::quiet_NaN();
        sum2+=v*v;
    }
    return std::sqrt(sum2);
}

static void style_hist(TH1D& h,int color,int style) {
    h.SetLineColor(color);
    h.SetLineWidth(3);
    h.SetLineStyle(style);
}

static void normalize_hist(TH1D& h) {
    const double a=h.Integral();
    if(a>0.0) h.Scale(1.0/a);
}

static void draw_distribution_comparison(
    const fs::path& path,
    const std::vector<double>& pass1_frac,
    const std::vector<double>& pass2_frac,
    const Summary& s1,
    const Summary& s2) {

    double xmax=0.0;
    for(double x:pass1_frac) if(std::isfinite(x)) xmax=std::max(xmax,100.0*x);
    for(double x:pass2_frac) if(std::isfinite(x)) xmax=std::max(xmax,100.0*x);
    xmax=std::min(20.0,std::max(4.0,1.10*xmax));

    TCanvas cv("c_pi0_pass1_comparison","",1050,720);
    cv.SetLeftMargin(.14);
    cv.SetRightMargin(.04);
    cv.SetBottomMargin(.14);
    cv.SetTopMargin(.16);
    cv.SetTicks(1,1);

    TH1D h1("h_pi0_pass1","",50,0,xmax);
    TH1D h2("h_pi0_pass2","",50,0,xmax);
    for(double x:pass1_frac) if(std::isfinite(x) && 100.0*x<xmax) h1.Fill(100.0*x);
    for(double x:pass2_frac) if(std::isfinite(x) && 100.0*x<xmax) h2.Fill(100.0*x);
    normalize_hist(h1);
    normalize_hist(h2);
    style_hist(h1,kGray+2,2);
    style_hist(h2,kBlue+1,1);

    const double ymax=1.25*std::max(h1.GetMaximum(),h2.GetMaximum());
    TH1D frame("h_pi0_dist_frame","",50,0,xmax);
    frame.SetMinimum(0.0);
    frame.SetMaximum(ymax);
    frame.GetXaxis()->SetTitle("Relative #pi^{0}-subtraction systematic (%)");
    frame.GetYaxis()->SetTitle("Fraction of populated bins");
    frame.GetXaxis()->SetTitleSize(.046);
    frame.GetYaxis()->SetTitleSize(.046);
    frame.GetXaxis()->SetLabelSize(.040);
    frame.GetYaxis()->SetLabelSize(.040);
    frame.GetYaxis()->SetTitleOffset(1.35);
    frame.Draw();
    h1.Draw("HIST SAME");
    h2.Draw("HIST SAME");

    TLegend leg(.59,.67,.93,.82);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(.031);
    leg.AddEntry(&h1,"Pass-1 assigned systematic","l");
    leg.AddEntry(&h2,"Pass-2 from 7.2% background variation","l");
    leg.Draw();

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(.037);
    t.DrawLatex(.14,.945,"#pi^{0}-background systematic: pass-1 / pass-2 comparison");
    t.SetTextSize(.026);
    std::ostringstream ss;
    ss<<std::fixed<<std::setprecision(2)
      <<"Median: pass-1 "<<100.0*s1.median<<"%, pass-2 "<<100.0*s2.median<<"%";
    t.DrawLatex(.14,.895,ss.str().c_str());

    cv.SaveAs(path.string().c_str());
}

struct KPoint {
    double x=0.0;
    double f=0.0;
};

static std::vector<std::pair<double,double>> binned_medians(
    const std::vector<KPoint>& pts,
    double lo,double hi,double width) {

    std::vector<std::pair<double,double>> out;
    for(double a=lo;a<hi;a+=width) {
        const double b=std::min(hi,a+width);
        std::vector<double> vals;
        for(const auto& p:pts)
            if(std::isfinite(p.x)&&std::isfinite(p.f)&&p.x>=a&&p.x<b)
                vals.push_back(100.0*p.f);
        if(vals.size()<5) continue;
        out.push_back({0.5*(a+b),quantile(vals,.5)});
    }
    return out;
}

static void draw_kinematic_comparison(
    const fs::path& path,
    const std::array<std::vector<KPoint>,4>& p1,
    const std::array<std::vector<KPoint>,4>& p2) {

    struct Spec { const char* x; double lo,hi,w; };
    const std::array<Spec,4> specs={{
        {"x_{B}",0.05,0.60,0.055},
        {"Q^{2} (GeV^{2})",1.0,8.0,0.7},
        {"|t| (GeV^{2})",0.0,1.05,0.105},
        {"#phi (deg)",0.0,360.0,30.0}
    }};

    TCanvas cv("c_pi0_kinematic_comparison","",1250,850);
    cv.Divide(2,2,.003,.003);

    for(int ip=0;ip<4;++ip) {
        cv.cd(ip+1);
        gPad->SetLeftMargin(.14);
        gPad->SetRightMargin(.04);
        gPad->SetBottomMargin(.15);
        gPad->SetTopMargin(.13);
        gPad->SetTicks(1,1);

        const auto a=binned_medians(p1[(size_t)ip],specs[(size_t)ip].lo,specs[(size_t)ip].hi,specs[(size_t)ip].w);
        const auto b=binned_medians(p2[(size_t)ip],specs[(size_t)ip].lo,specs[(size_t)ip].hi,specs[(size_t)ip].w);

        double ymax=0.0;
        for(const auto& q:a) ymax=std::max(ymax,q.second);
        for(const auto& q:b) ymax=std::max(ymax,q.second);
        ymax=std::max(1.0,1.30*ymax);

        TH1D frame(("h_pi0_kframe_"+std::to_string(ip)).c_str(),"",100,specs[(size_t)ip].lo,specs[(size_t)ip].hi);
        frame.SetMinimum(0.0);
        frame.SetMaximum(ymax);
        frame.GetXaxis()->SetTitle(specs[(size_t)ip].x);
        frame.GetYaxis()->SetTitle("Median relative #pi^{0} systematic (%)");
        frame.GetXaxis()->SetTitleSize(.048);
        frame.GetYaxis()->SetTitleSize(.044);
        frame.GetXaxis()->SetLabelSize(.040);
        frame.GetYaxis()->SetLabelSize(.038);
        frame.GetYaxis()->SetTitleOffset(1.38);
        frame.Draw();

        TGraph g1((int)a.size()),g2((int)b.size());
        for(int i=0;i<(int)a.size();++i) g1.SetPoint(i,a[(size_t)i].first,a[(size_t)i].second);
        for(int i=0;i<(int)b.size();++i) g2.SetPoint(i,b[(size_t)i].first,b[(size_t)i].second);
        g1.SetMarkerStyle(24); g1.SetMarkerSize(.9); g1.SetMarkerColor(kGray+2); g1.SetLineColor(kGray+2); g1.SetLineStyle(2);
        g2.SetMarkerStyle(20); g2.SetMarkerSize(.9); g2.SetMarkerColor(kBlue+1); g2.SetLineColor(kBlue+1);
        if(!a.empty()) g1.Draw("PL SAME");
        if(!b.empty()) g2.Draw("PL SAME");

        TLatex lab;
        lab.SetNDC();
        lab.SetTextFont(42);
        lab.SetTextSize(.035);
        const std::string panel=std::string("(")+char('a'+ip)+")";
        lab.DrawLatex(.17,.84,panel.c_str());

        if(ip==0) {
            TLegend leg(.55,.66,.93,.81);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextFont(42);
            leg.SetTextSize(.027);
            leg.AddEntry(&g1,"Pass-1","pl");
            leg.AddEntry(&g2,"Pass-2","pl");
            leg.Draw();
        }
    }

    cv.cd(0);
    TLatex t;
    t.SetNDC();
    t.SetTextAlign(22);
    t.SetTextFont(42);
    t.SetTextSize(.022);
    t.DrawLatex(.50,.978,"Kinematic comparison of the #pi^{0}-background systematic");

    cv.SaveAs(path.string().c_str());
}

static void draw_background_fraction_comparison(
    const fs::path& path,
    const std::vector<double>& c1,
    const std::vector<double>& c2,
    const Summary& s1,
    const Summary& s2) {

    TCanvas cv("c_pi0_background_fraction","",1050,720);
    cv.SetLeftMargin(.14);
    cv.SetRightMargin(.04);
    cv.SetBottomMargin(.14);
    cv.SetTopMargin(.16);
    cv.SetTicks(1,1);

    TH1D h1("h_pi0_c1","",50,0,0.45);
    TH1D h2("h_pi0_c2","",50,0,0.45);
    for(double x:c1) if(std::isfinite(x)&&x>=0.0&&x<.45) h1.Fill(x);
    for(double x:c2) if(std::isfinite(x)&&x>=0.0&&x<.45) h2.Fill(x);
    normalize_hist(h1);
    normalize_hist(h2);
    style_hist(h1,kGray+2,2);
    style_hist(h2,kBlue+1,1);

    TH1D frame("h_pi0_cframe","",50,0,.45);
    frame.SetMinimum(0.0);
    frame.SetMaximum(1.25*std::max(h1.GetMaximum(),h2.GetMaximum()));
    frame.GetXaxis()->SetTitle("Effective acceptance-corrected #pi^{0} contamination");
    frame.GetYaxis()->SetTitle("Fraction of populated bins");
    frame.GetXaxis()->SetTitleSize(.046);
    frame.GetYaxis()->SetTitleSize(.046);
    frame.GetXaxis()->SetLabelSize(.040);
    frame.GetYaxis()->SetLabelSize(.040);
    frame.GetYaxis()->SetTitleOffset(1.35);
    frame.Draw();
    h1.Draw("HIST SAME");
    h2.Draw("HIST SAME");

    TLegend leg(.59,.67,.93,.82);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(.031);
    leg.AddEntry(&h1,"Pass-1 inferred from assigned systematic","l");
    leg.AddEntry(&h2,"Pass-2 directly calculated","l");
    leg.Draw();

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(.037);
    t.DrawLatex(.14,.945,"Effective #pi^{0} contamination: pass-1 / pass-2 consistency");
    t.SetTextSize(.025);
    std::ostringstream ss;
    ss<<std::fixed<<std::setprecision(3)
      <<"Median: pass-1 "<<s1.median<<", pass-2 "<<s2.median;
    t.DrawLatex(.14,.895,ss.str().c_str());

    cv.SaveAs(path.string().c_str());
}

} // namespace

bool pi0_systematics(
    const std::string& csv_path,
    const std::string& pass1_summary_path,
    const std::string& output_dir) {

    try {
        fs::create_directories(output_dir);

        Csv c=read_csv(csv_path);
        const Csv p1=read_csv(pass1_summary_path);

        const int c_pi0_abs=ensure_column(c,"Syst. err (pi0 subtraction)");
        const int c_pi0_frac=ensure_column(c,"pi0 subtraction sys frac, 10.6 GeV");
        const int c_pi0_ceff=ensure_column(c,"pi0 effective contamination, 10.6 GeV");
        const int c_pi0_abs_sp=ensure_column(c,"Syst. err (pi0 subtraction), Sp19 Inb (10.2 GeV)");
        const int c_pi0_frac_sp=ensure_column(c,"pi0 subtraction sys frac, Sp19 Inb (10.2 GeV)");
        const int c_pi0_ceff_sp=ensure_column(c,"pi0 effective contamination, Sp19 Inb (10.2 GeV)");
        const int c_total=require_column(c,"Syst. err (point-to-point total)");

        const int c_xs10=require_column(c,xs_col("10.6 GeV"));
        const int c_xssp=require_column(c,xs_col("Sp19 Inb"));

        std::vector<double> pass2_frac,pass2_ceff,sp_frac,sp_ceff;
        std::array<std::vector<KPoint>,4> p2k;

        const std::array<std::string,4> mean_cols={{
            "xBavg, 10.6 GeV",
            "Q2avg, 10.6 GeV",
            "t_abs_avg, 10.6 GeV",
            "phiavg, 10.6 GeV"
        }};

        std::array<int,4> mean_idx;
        for(int i=0;i<4;++i) mean_idx[(size_t)i]=require_column(c,mean_cols[(size_t)i]);

        size_t n10=0,nsp=0;
        for(auto& row:c.rows) {
            double sig_sum=0.0,bkg_sum=0.0;
            int nperiod=0;
            for(const auto& per:ten6_periods()) {
                double sig=0.0,bkg=0.0;
                if(period_background_to_signal(c,row,per,sig,bkg)) {
                    sig_sum+=sig;
                    bkg_sum+=bkg;
                    ++nperiod;
                }
            }

            const Triple xs10=parse_triple(row[(size_t)c_xs10]);
            if(nperiod>0 && sig_sum>0.0 && xs10.ok) {
                const double ratio=bkg_sum/sig_sum;
                const double frac=kPi0BackgroundRelativeUncertainty*ratio;
                const double ceff=ratio/(1.0+ratio);
                const double abs=std::fabs(xs10.value)*frac;

                row[(size_t)c_pi0_abs]=fmt(abs);
                row[(size_t)c_pi0_frac]=fmt(frac);
                row[(size_t)c_pi0_ceff]=fmt(ceff);
                pass2_frac.push_back(frac);
                pass2_ceff.push_back(ceff);
                for(int iv=0;iv<4;++iv) {
                    const double x=number(row[(size_t)mean_idx[(size_t)iv]]);
                    if(std::isfinite(x)) p2k[(size_t)iv].push_back({x,frac});
                }
                ++n10;
            } else {
                row[(size_t)c_pi0_abs].clear();
                row[(size_t)c_pi0_frac].clear();
                row[(size_t)c_pi0_ceff].clear();
            }

            double sigsp=0.0,bkgsp=0.0;
            const Triple xssp=parse_triple(row[(size_t)c_xssp]);
            if(period_background_to_signal(c,row,"Sp19 Inb",sigsp,bkgsp) &&
               sigsp>0.0 && xssp.ok) {
                const double ratio=bkgsp/sigsp;
                const double frac=kPi0BackgroundRelativeUncertainty*ratio;
                const double ceff=ratio/(1.0+ratio);
                row[(size_t)c_pi0_abs_sp]=fmt(std::fabs(xssp.value)*frac);
                row[(size_t)c_pi0_frac_sp]=fmt(frac);
                row[(size_t)c_pi0_ceff_sp]=fmt(ceff);
                sp_frac.push_back(frac);
                sp_ceff.push_back(ceff);
                ++nsp;
            } else {
                row[(size_t)c_pi0_abs_sp].clear();
                row[(size_t)c_pi0_frac_sp].clear();
                row[(size_t)c_pi0_ceff_sp].clear();
            }

            // The production point-to-point total is the combined 10.6-GeV
            // quantity.  Recompute it after replacing the inherited pass-1
            // pi0 term with the pass-2 value.
            const double total=recompute_point_to_point_total(c,row);
            row[(size_t)c_total]=fmt(total);
        }

        // Pass-1 validation inputs.  Since pass 1 used the same 7.2% relative
        // background uncertainty, its assigned cross-section systematic can be
        // inverted to an effective acceptance-corrected contamination:
        //
        //   f_pi0 = 0.072 B/S,
        //   c_eff = B/(S+B) = f_pi0/(0.072+f_pi0).
        const int p1_val=require_column(p1,"val");
        const int p1_pi0=require_column(p1,"Syst. err (pi0 subtraction)");
        const std::array<int,8> p1_edge_idx={{
            require_column(p1,"xBmin"),require_column(p1,"xBmax"),
            require_column(p1,"Q2min"),require_column(p1,"Q2max"),
            require_column(p1,"tmin"),require_column(p1,"tmax"),
            require_column(p1,"phimin"),require_column(p1,"phimax")
        }};

        std::vector<double> pass1_frac,pass1_ceff;
        std::array<std::vector<KPoint>,4> p1k;

        for(const auto& row:p1.rows) {
            const double xs=number(row[(size_t)p1_val]);
            const double e=number(row[(size_t)p1_pi0]);
            if(!std::isfinite(xs)||!std::isfinite(e)||std::fabs(xs)<=0.0||e<0.0) continue;
            const double frac=e/std::fabs(xs);
            const double ceff=frac/(kPi0BackgroundRelativeUncertainty+frac);
            pass1_frac.push_back(frac);
            pass1_ceff.push_back(ceff);

            const std::array<double,4> x={{
                0.5*(number(row[(size_t)p1_edge_idx[0]])+number(row[(size_t)p1_edge_idx[1]])),
                0.5*(number(row[(size_t)p1_edge_idx[2]])+number(row[(size_t)p1_edge_idx[3]])),
                0.5*(number(row[(size_t)p1_edge_idx[4]])+number(row[(size_t)p1_edge_idx[5]])),
                0.5*(number(row[(size_t)p1_edge_idx[6]])+number(row[(size_t)p1_edge_idx[7]]))
            }};
            for(int iv=0;iv<4;++iv)
                if(std::isfinite(x[(size_t)iv])) p1k[(size_t)iv].push_back({x[(size_t)iv],frac});
        }

        const Summary s1=summarize(pass1_frac);
        const Summary s2=summarize(pass2_frac);
        const Summary ssp=summarize(sp_frac);
        const Summary c1=summarize(pass1_ceff);
        const Summary c2=summarize(pass2_ceff);
        const Summary csp=summarize(sp_ceff);

        write_csv_atomic(csv_path,c);

        {
            std::ofstream out(fs::path(output_dir)/"pi0_systematic_summary.csv");
            out<<"sample,n,mean_frac,median_frac,q16_frac,q84_frac,q95_frac,median_effective_contamination\n";
            auto write=[&](const std::string& name,const Summary& s,const Summary& cs) {
                out<<name<<','<<s.n<<','<<s.mean<<','<<s.median<<','<<s.q16<<','<<s.q84<<','<<s.q95<<','<<cs.median<<'\n';
            };
            write("Pass-1",s1,c1);
            write("Pass-2 combined 10.6 GeV",s2,c2);
            write("Pass-2 Sp19 Inb",ssp,csp);
        }

        draw_distribution_comparison(
            fs::path(output_dir)/"pi0_systematic_pass1_comparison.png",
            pass1_frac,pass2_frac,s1,s2);

        draw_background_fraction_comparison(
            fs::path(output_dir)/"pi0_effective_contamination_pass1_comparison.png",
            pass1_ceff,pass2_ceff,c1,c2);

        draw_kinematic_comparison(
            fs::path(output_dir)/"pi0_systematic_kinematic_pass1_comparison.png",
            p1k,p2k);

        std::cout<<"[pi0-systematics] Assigned "<<100.0*kPi0BackgroundRelativeUncertainty
                 <<"% relative uncertainty to the acceptance-corrected pi0 background.\n";
        std::cout<<"[pi0-systematics] Combined 10.6 GeV bins: "<<n10
                 <<"; median pi0 cross-section systematic = "<<100.0*s2.median<<"%.\n";
        std::cout<<"[pi0-systematics] Sp19 Inb bins: "<<nsp
                 <<"; median pi0 cross-section systematic = "<<100.0*ssp.median<<"%.\n";
        std::cout<<"[pi0-systematics] Pass-1 comparison median = "<<100.0*s1.median
                 <<"% (relative cross-section effect).\n";
        std::cout<<"[pi0-systematics] Validation outputs: "<<output_dir<<"\n";
        return true;

    } catch(const std::exception& e) {
        std::cerr<<"[pi0-systematics] FATAL: "<<e.what()<<"\n";
        return false;
    }
}
