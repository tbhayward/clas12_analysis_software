
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
    // Reserve a dedicated title/subtitle band above the frame.
    cv.SetTopMargin(.18);
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
    frame.SetStats(0);
    h1.SetStats(0);
    h2.SetStats(0);
    frame.DrawCopy();
    h1.DrawClone("HIST SAME");
    h2.DrawClone("HIST SAME");

    TLegend leg(.60,.69,.92,.82);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(.029);
    leg.AddEntry(&h1,"Pass-1","l");
    leg.AddEntry(&h2,"Pass-2 (7.2% background uncertainty)","l");
    leg.DrawClone();

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
        frame.SetStats(0);
        frame.DrawCopy();

        TGraph g1((int)a.size()),g2((int)b.size());
        for(int i=0;i<(int)a.size();++i) g1.SetPoint(i,a[(size_t)i].first,a[(size_t)i].second);
        for(int i=0;i<(int)b.size();++i) g2.SetPoint(i,b[(size_t)i].first,b[(size_t)i].second);
        g1.SetMarkerStyle(24); g1.SetMarkerSize(.9); g1.SetMarkerColor(kGray+2); g1.SetLineColor(kGray+2); g1.SetLineStyle(2);
        g2.SetMarkerStyle(20); g2.SetMarkerSize(.9); g2.SetMarkerColor(kBlue+1); g2.SetLineColor(kBlue+1);
        if(!a.empty()) g1.DrawClone("PL SAME");
        if(!b.empty()) g2.DrawClone("PL SAME");

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
            leg.DrawClone();
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
    // Reserve a dedicated title/subtitle band above the frame.
    cv.SetTopMargin(.18);
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
    frame.SetStats(0);
    h1.SetStats(0);
    h2.SetStats(0);
    frame.DrawCopy();
    h1.DrawClone("HIST SAME");
    h2.DrawClone("HIST SAME");

    TLegend leg(.59,.67,.93,.82);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(.029);
    leg.AddEntry(&h1,"Pass-1 inferred","l");
    leg.AddEntry(&h2,"Pass-2 calculated","l");
    leg.DrawClone();

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

struct BinEdges {
    double xb0=0.0, xb1=0.0;
    double q20=0.0, q21=0.0;
    double t0=0.0, t1=0.0;
    double ph0=0.0, ph1=0.0;
};

static long long edge_code(double x) {
    return std::llround(1000000.0*x);
}

static std::string bin_key(const BinEdges& b) {
    std::ostringstream ss;
    ss<<edge_code(b.xb0)<<':'
      <<edge_code(b.xb1)<<':'
      <<edge_code(b.q20)<<':'
      <<edge_code(b.q21)<<':'
      <<edge_code(b.t0)<<':'
      <<edge_code(b.t1)<<':'
      <<edge_code(b.ph0)<<':'
      <<edge_code(b.ph1);
    return ss.str();
}

static BinEdges edges_from_row(
    const Csv& c,
    const std::vector<std::string>& row) {

    BinEdges b;
    b.xb0=number(row[(size_t)require_column(c,"xBmin")]);
    b.xb1=number(row[(size_t)require_column(c,"xBmax")]);
    b.q20=number(row[(size_t)require_column(c,"Q2min")]);
    b.q21=number(row[(size_t)require_column(c,"Q2max")]);
    b.t0=number(row[(size_t)require_column(c,"t_abs_min")]);
    b.t1=number(row[(size_t)require_column(c,"t_abs_max")]);
    b.ph0=number(row[(size_t)require_column(c,"phimin")]);
    b.ph1=number(row[(size_t)require_column(c,"phimax")]);
    return b;
}

struct MatchedPi0Point {
    BinEdges edges;
    std::string period;
    double pass1=std::numeric_limits<double>::quiet_NaN();
    double pass2=std::numeric_limits<double>::quiet_NaN();
    double ratio=std::numeric_limits<double>::quiet_NaN();
};

static void write_matched_pi0_csv(
    const fs::path& path,
    const std::vector<MatchedPi0Point>& pts) {

    std::ofstream out(path);
    out<<"period,xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phimin,phimax,"
          "pass1_contamination,pass2_contamination,pass2_over_pass1\n";
    for(const auto& p:pts) {
        out<<csv_escape(p.period)<<','
           <<p.edges.xb0<<','<<p.edges.xb1<<','
           <<p.edges.q20<<','<<p.edges.q21<<','
           <<p.edges.t0<<','<<p.edges.t1<<','
           <<p.edges.ph0<<','<<p.edges.ph1<<','
           <<p.pass1<<','<<p.pass2<<','<<p.ratio<<'\n';
    }
}

static void draw_matched_ratio_distributions(
    const fs::path& path,
    const std::vector<MatchedPi0Point>& pts) {

    const std::array<std::string,4> periods={{
        "Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"
    }};

    TCanvas cv("c_pi0_matched_pass1_ratio","",1300,900);
    cv.Divide(2,2,.003,.003);

    for(int ip=0;ip<4;++ip) {
        std::vector<double> ratios;
        for(const auto& p:pts) {
            if(p.period==periods[(size_t)ip] &&
               std::isfinite(p.ratio) && p.ratio>0.0 &&
               p.ratio<8.0) {
                ratios.push_back(p.ratio);
            }
        }

        cv.cd(ip+1);
        gPad->SetLeftMargin(.14);
        gPad->SetRightMargin(.04);
        gPad->SetBottomMargin(.15);
        gPad->SetTopMargin(.17);
        gPad->SetTicks(1,1);

        TH1D h(("h_pi0_matched_ratio_"+std::to_string(ip)).c_str(),"",60,0.0,4.0);
        h.SetStats(0);
        for(double r:ratios) h.Fill(r);
        normalize_hist(h);
        h.SetLineWidth(3);
        h.SetLineColor(kBlue+1);

        TH1D frame(("h_pi0_matched_ratio_frame_"+std::to_string(ip)).c_str(),"",60,0.0,4.0);
        frame.SetStats(0);
        frame.SetMinimum(0.0);
        frame.SetMaximum(std::max(0.05,1.25*h.GetMaximum()));
        frame.GetXaxis()->SetTitle("Pass-2 / pass-1 contamination");
        frame.GetYaxis()->SetTitle("Fraction of matched bins");
        frame.GetXaxis()->SetTitleSize(.047);
        frame.GetYaxis()->SetTitleSize(.044);
        frame.GetXaxis()->SetLabelSize(.040);
        frame.GetYaxis()->SetLabelSize(.038);
        frame.GetYaxis()->SetTitleOffset(1.35);
        frame.DrawCopy();
        h.DrawClone("HIST SAME");

        TLine one(1.0,0.0,1.0,frame.GetMaximum());
        one.SetLineStyle(2);
        one.SetLineWidth(2);
        one.SetLineColor(kGray+2);
        one.DrawClone();

        TLatex lab;
        lab.SetNDC();
        lab.SetTextFont(42);
        lab.SetTextSize(.034);
        lab.DrawLatex(.18,.87,periods[(size_t)ip].c_str());

        if(!ratios.empty()) {
            const Summary sr=summarize(ratios);
            std::ostringstream ss;
            ss<<std::fixed<<std::setprecision(2)
              <<"median ratio = "<<sr.median
              <<", N = "<<sr.n;
            lab.SetTextSize(.026);
            lab.DrawLatex(.18,.81,ss.str().c_str());
        }
    }

    cv.cd(0);
    TLatex title;
    title.SetNDC();
    title.SetTextAlign(22);
    title.SetTextFont(42);
    title.SetTextSize(.021);
    title.DrawLatex(.50,.975,
        "Matched-bin comparison of pass-2 and preliminary pass-1 #pi^{0} contamination");
    title.SetTextSize(.0155);
    title.DrawLatex(.50,.949,
        "Pass-1 inbending/outbending contamination is compared with the corresponding pass-2 polarity");

    cv.SaveAs(path.string().c_str());
}

static void write_high_contamination_bins(
    const fs::path& path,
    const Csv& c) {

    const int iceff=require_column(c,"pi0 effective contamination, 10.6 GeV");
    const int ispeff=require_column(c,"pi0 effective contamination, Sp19 Inb (10.2 GeV)");

    const std::array<std::string,8> edge_names={{
        "xBmin","xBmax","Q2min","Q2max",
        "t_abs_min","t_abs_max","phimin","phimax"
    }};
    std::array<int,8> eidx;
    for(int i=0;i<8;++i) eidx[(size_t)i]=require_column(c,edge_names[(size_t)i]);

    std::ofstream out(path);
    out<<"sample,xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phimin,phimax,"
          "effective_contamination\n";

    for(const auto& row:c.rows) {
        const std::array<std::pair<std::string,int>,2> samples={{
            {"10.6 GeV",iceff},
            {"Sp19 Inb",ispeff}
        }};
        for(const auto& s:samples) {
            const double ceff=number(row[(size_t)s.second]);
            if(!std::isfinite(ceff)||ceff<=0.50) continue;
            out<<csv_escape(s.first);
            for(int j=0;j<8;++j) out<<','<<row[(size_t)eidx[(size_t)j]];
            out<<','<<ceff<<'\n';
        }
    }
}

static void draw_high_contamination_kinematics(
    const fs::path& path,
    const Csv& c) {

    const int iceff=require_column(c,"pi0 effective contamination, 10.6 GeV");
    const std::array<std::string,4> mean_names={{
        "xBavg, 10.6 GeV",
        "Q2avg, 10.6 GeV",
        "t_abs_avg, 10.6 GeV",
        "phiavg, 10.6 GeV"
    }};
    std::array<int,4> midx;
    for(int i=0;i<4;++i) midx[(size_t)i]=require_column(c,mean_names[(size_t)i]);

    struct Spec { const char* title; double lo,hi; };
    const std::array<Spec,4> specs={{
        {"x_{B}",0.05,0.60},
        {"Q^{2} (GeV^{2})",1.0,8.0},
        {"|t| (GeV^{2})",0.0,1.05},
        {"#phi (deg)",0.0,360.0}
    }};

    TCanvas cv("c_pi0_high_contamination","",1250,850);
    cv.Divide(2,2,.003,.003);

    for(int iv=0;iv<4;++iv) {
        std::vector<std::pair<double,double>> pts;
        for(const auto& row:c.rows) {
            const double ceff=number(row[(size_t)iceff]);
            const double x=number(row[(size_t)midx[(size_t)iv]]);
            if(std::isfinite(ceff)&&std::isfinite(x)&&ceff>0.50)
                pts.push_back({x,ceff});
        }

        cv.cd(iv+1);
        gPad->SetLeftMargin(.14);
        gPad->SetRightMargin(.04);
        gPad->SetBottomMargin(.15);
        gPad->SetTopMargin(.15);
        gPad->SetTicks(1,1);

        TH1D frame(("h_pi0_high_frame_"+std::to_string(iv)).c_str(),"",100,
                   specs[(size_t)iv].lo,specs[(size_t)iv].hi);
        frame.SetStats(0);
        frame.SetMinimum(.48);
        frame.SetMaximum(1.02);
        frame.GetXaxis()->SetTitle(specs[(size_t)iv].title);
        frame.GetYaxis()->SetTitle("Effective #pi^{0} contamination");
        frame.GetXaxis()->SetTitleSize(.048);
        frame.GetYaxis()->SetTitleSize(.044);
        frame.GetXaxis()->SetLabelSize(.040);
        frame.GetYaxis()->SetLabelSize(.038);
        frame.GetYaxis()->SetTitleOffset(1.35);
        frame.DrawCopy();

        TGraph g((int)pts.size());
        for(int i=0;i<(int)pts.size();++i)
            g.SetPoint(i,pts[(size_t)i].first,pts[(size_t)i].second);
        g.SetMarkerStyle(20);
        g.SetMarkerSize(.8);
        g.SetMarkerColor(kRed+1);
        if(!pts.empty()) g.DrawClone("P SAME");

        TLatex lab;
        lab.SetNDC();
        lab.SetTextFont(42);
        lab.SetTextSize(.034);
        const std::string panel=std::string("(")+char('a'+iv)+")";
        lab.DrawLatex(.18,.85,panel.c_str());
    }

    cv.cd(0);
    TLatex title;
    title.SetNDC();
    title.SetTextAlign(22);
    title.SetTextFont(42);
    title.SetTextSize(.021);
    title.DrawLatex(.50,.975,
        "Kinematics of the rare high-#pi^{0}-contamination bins");
    title.SetTextSize(.0155);
    title.DrawLatex(.50,.949,
        "Combined 10.6 GeV bins with effective contamination greater than 0.50");

    cv.SaveAs(path.string().c_str());
}


} // namespace

bool pi0_systematics(
    const std::string& csv_path,
    const std::string& pass1_summary_path,
    const std::string& preliminary_pass1_csv_path,
    const std::string& output_dir) {

    try {
        fs::create_directories(output_dir);

        Csv c=read_csv(csv_path);
        const Csv p1=read_csv(pass1_summary_path);
        const Csv p1_prelim=read_csv(preliminary_pass1_csv_path);

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

        // Direct matched-bin comparison to the preliminary Sangbaek/pass-1 CSV.
        // That CSV stores one inbending and one outbending contamination value
        // per nominal kinematic bin.  We compare each pass-2 period to the
        // pass-1 value with the same torus polarity.
        const int ppre_valid=require_column(p1_prelim,"valid bin");
        const int ppre_cin=require_column(p1_prelim,"contamination ratio, inbending");
        const int ppre_cout=require_column(p1_prelim,"contamination ratio, outbending");

        std::unordered_map<std::string,std::pair<double,double>> preliminary_by_bin;
        for(const auto& row:p1_prelim.rows) {
            const std::string valid=row[(size_t)ppre_valid];
            if(!(valid=="True"||valid=="true"||valid=="1"||valid=="TRUE")) continue;

            BinEdges b;
            b.xb0=number(row[(size_t)require_column(p1_prelim,"xBmin")]);
            b.xb1=number(row[(size_t)require_column(p1_prelim,"xBmax")]);
            b.q20=number(row[(size_t)require_column(p1_prelim,"Q2min")]);
            b.q21=number(row[(size_t)require_column(p1_prelim,"Q2max")]);
            b.t0=number(row[(size_t)require_column(p1_prelim,"t_abs_min")]);
            b.t1=number(row[(size_t)require_column(p1_prelim,"t_abs_max")]);
            b.ph0=number(row[(size_t)require_column(p1_prelim,"phimin")]);
            b.ph1=number(row[(size_t)require_column(p1_prelim,"phimax")]);

            const double cin=number(row[(size_t)ppre_cin]);
            const double cout=number(row[(size_t)ppre_cout]);
            preliminary_by_bin[bin_key(b)]={cin,cout};
        }

        std::vector<MatchedPi0Point> matched_points;
        const std::array<std::pair<std::string,bool>,4> matched_periods={{
            {"Fa18 Inb",true},
            {"Fa18 Out",false},
            {"Sp18 Inb",true},
            {"Sp18 Out",false}
        }};

        for(const auto& row:c.rows) {
            const BinEdges b=edges_from_row(c,row);
            const auto it=preliminary_by_bin.find(bin_key(b));
            if(it==preliminary_by_bin.end()) continue;

            for(const auto& per:matched_periods) {
                const auto ic=c.index.find("contamination ratio, "+per.first);
                if(ic==c.index.end()) continue;

                const Triple tr=parse_triple(row[(size_t)ic->second]);
                if(!tr.ok||!std::isfinite(tr.value)||tr.value<0.0) continue;

                const double oldc=per.second ? it->second.first : it->second.second;
                if(!std::isfinite(oldc)||oldc<=0.0) continue;

                MatchedPi0Point mp;
                mp.edges=b;
                mp.period=per.first;
                mp.pass1=oldc;
                mp.pass2=tr.value;
                mp.ratio=tr.value/oldc;
                matched_points.push_back(mp);
            }
        }

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

        write_matched_pi0_csv(
            fs::path(output_dir)/"pi0_matched_pass1_pass2_contamination.csv",
            matched_points);

        draw_matched_ratio_distributions(
            fs::path(output_dir)/"pi0_matched_pass1_pass2_ratio_distributions.png",
            matched_points);

        write_high_contamination_bins(
            fs::path(output_dir)/"pi0_high_contamination_bins.csv",
            c);

        draw_high_contamination_kinematics(
            fs::path(output_dir)/"pi0_high_contamination_kinematics.png",
            c);

        std::cout<<"[pi0-systematics] Matched preliminary pass-1/pass-2 contamination points: "
                 <<matched_points.size()<<".\n";
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
