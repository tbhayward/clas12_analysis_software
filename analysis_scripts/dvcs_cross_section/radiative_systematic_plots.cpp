#include "radiative_systematic_plots.h"

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>

#include <algorithm>
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
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

static const std::string kColFrad10 =
    "Syst.err (Frad)";
static const std::string kColFrad19 =
    "Syst.err (Frad), Sp19 Inb (10.2 GeV)";
static const std::string kColXs10 =
    "normed cross sections, ep->epg, exp, 10.6 GeV, unpol";
static const std::string kColXs19 =
    "normed cross sections, ep->epg, exp, Sp19 Inb, unpol";

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    std::unordered_map<std::string,int> index;
};

struct RelPair {
    double ten6 = std::numeric_limits<double>::quiet_NaN();
    double sp19 = std::numeric_limits<double>::quiet_NaN();
};

static std::string trim(const std::string& s) {
    size_t b=0;
    while(b<s.size() && std::isspace((unsigned char)s[b])) ++b;
    size_t e=s.size();
    while(e>b && std::isspace((unsigned char)s[e-1])) --e;
    return s.substr(b,e-b);
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes=false;

    for(size_t i=0;i<line.size();++i){
        const char c=line[i];

        if(c=='"'){
            if(in_quotes && i+1<line.size() && line[i+1]=='"'){
                cur.push_back('"');
                ++i;
            } else {
                in_quotes=!in_quotes;
            }
        } else if(c==',' && !in_quotes){
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }

    out.push_back(cur);
    return out;
}

static CsvTable read_csv(const std::string& path) {
    std::ifstream f(path);
    if(!f) throw std::runtime_error("Could not open CSV: "+path);

    CsvTable t;
    std::string line;

    if(!std::getline(f,line)){
        throw std::runtime_error("Empty CSV: "+path);
    }

    t.header=split_csv_line(line);
    for(int i=0;i<(int)t.header.size();++i){
        t.index[t.header[(size_t)i]]=i;
    }

    while(std::getline(f,line)){
        if(line.empty()) continue;
        auto row=split_csv_line(line);
        row.resize(t.header.size());
        t.rows.push_back(std::move(row));
    }

    return t;
}

static int col(const CsvTable& t,const std::string& name) {
    const auto it=t.index.find(name);
    if(it==t.index.end()){
        throw std::runtime_error("Missing required column: "+name);
    }
    return it->second;
}

static double number(const std::string& raw) {
    const std::string s=trim(raw);
    if(s.empty()) return std::numeric_limits<double>::quiet_NaN();

    char* end=nullptr;
    const double v=std::strtod(s.c_str(),&end);

    return end==s.c_str()
        ? std::numeric_limits<double>::quiet_NaN()
        : v;
}

static bool tuple_first(const std::string& raw,double& value) {
    std::string s=trim(raw);

    while(s.size()>=2 && s.front()=='"' && s.back()=='"'){
        s=trim(s.substr(1,s.size()-2));
    }

    if(s.empty()) return false;

    if(s.front()!='('){
        value=number(s);
        return std::isfinite(value);
    }

    if(s.back()!=')') return false;

    s=s.substr(1,s.size()-2);
    const size_t comma=s.find(',');
    value=number(comma==std::string::npos ? s : s.substr(0,comma));
    return std::isfinite(value);
}

static double cell_number(const CsvTable& t,
                          const std::vector<std::string>& row,
                          const std::string& name) {
    return number(row[(size_t)col(t,name)]);
}

static double quantile(std::vector<double> v,double q) {
    if(v.empty()) return std::numeric_limits<double>::quiet_NaN();

    std::sort(v.begin(),v.end());

    if(q<=0.0) return v.front();
    if(q>=1.0) return v.back();

    const double p=q*(v.size()-1);
    const size_t i=(size_t)std::floor(p);
    const size_t j=(size_t)std::ceil(p);
    const double f=p-i;

    return v[i]*(1.0-f)+v[j]*f;
}

static double mean(const std::vector<double>& v) {
    if(v.empty()) return std::numeric_limits<double>::quiet_NaN();
    return std::accumulate(v.begin(),v.end(),0.0)/(double)v.size();
}

static bool approx(double a,double b,double eps=5.0e-6) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a-b)<eps;
}

static RelPair relative_uncertainties(
    const CsvTable& t,
    const std::vector<std::string>& row) {

    RelPair r;

    const double fr10=number(row[(size_t)col(t,kColFrad10)]);
    const double fr19=number(row[(size_t)col(t,kColFrad19)]);

    double xs10=std::numeric_limits<double>::quiet_NaN();
    double xs19=std::numeric_limits<double>::quiet_NaN();

    const bool ok10=tuple_first(row[(size_t)col(t,kColXs10)],xs10);
    const bool ok19=tuple_first(row[(size_t)col(t,kColXs19)],xs19);

    if(std::isfinite(fr10) && fr10>=0.0 &&
       ok10 && std::fabs(xs10)>0.0){
        r.ten6=100.0*fr10/std::fabs(xs10);
    }

    if(std::isfinite(fr19) && fr19>=0.0 &&
       ok19 && std::fabs(xs19)>0.0){
        r.sp19=100.0*fr19/std::fabs(xs19);
    }

    return r;
}

static void style_canvas(TCanvas& c) {
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.035);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.11);
    c.SetTicks(1,1);
}

static void note_label(const char* text) {
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.028);
    l.DrawLatex(0.13,0.925,text);
}

} // namespace

bool make_radiative_systematic_analysis_note_plots(
    const std::string& csv_path,
    const std::string& output_dir) {

    try {
        fs::create_directories(output_dir);

        const CsvTable t=read_csv(csv_path);

        for(const auto& name:std::vector<std::string>{
                kColFrad10,kColFrad19,kColXs10,kColXs19,
                "xBmin","xBmax","Q2min","Q2max",
                "t_abs_min","t_abs_max","phimin","phimax"}){
            (void)col(t,name);
        }

        std::vector<double> all10;
        std::vector<double> all19;

        for(const auto& row:t.rows){
            const RelPair r=relative_uncertainties(t,row);
            if(std::isfinite(r.ten6) && r.ten6>=0.0) all10.push_back(r.ten6);
            if(std::isfinite(r.sp19) && r.sp19>=0.0) all19.push_back(r.sp19);
        }

        if(all10.empty()){
            throw std::runtime_error("No valid 10.6-GeV relative Frad systematic values.");
        }

        {
            std::ofstream o(fs::path(output_dir)/"radiative_systematic_summary.csv");
            if(!o) throw std::runtime_error("Could not write summary CSV.");

            o<<"sample,n,mean_percent,median_percent,p16_percent,"
             <<"p84_percent,p95_percent,prescription\n";

            auto write_row=[&](const std::string& label,
                               const std::vector<double>& v,
                               const std::string& prescription){
                o<<label<<","<<v.size()<<","
                 <<std::setprecision(10)
                 <<mean(v)<<","
                 <<quantile(v,0.50)<<","
                 <<quantile(v,0.16)<<","
                 <<quantile(v,0.84)<<","
                 <<quantile(v,0.95)<<","
                 <<prescription<<"\n";
            };

            write_row("10.6 GeV",all10,"pass-1 bin-by-bin result");

            if(!all19.empty()){
                write_row(
                    "Sp19 Inb 10.2 GeV",
                    all19,
                    "1.05 x 10.6-GeV fractional uncertainty"
                );
            }
        }

        {
            const double p99_10=quantile(all10,0.99);
            const double p99_19=all19.empty() ? p99_10 : quantile(all19,0.99);

            double xmax=std::max(
                10.0,
                std::ceil(std::max(p99_10,p99_19)+1.0)
            );
            xmax=std::min(30.0,xmax);

            TH1D h10("h_rad_sys_10p6","",36,0.0,xmax);
            TH1D h19("h_rad_sys_10p2","",36,0.0,xmax);

            for(double v:all10) if(v>=0.0 && v<xmax) h10.Fill(v);
            for(double v:all19) if(v>=0.0 && v<xmax) h19.Fill(v);

            if(h10.Integral()>0.0) h10.Scale(1.0/h10.Integral());
            if(h19.Integral()>0.0) h19.Scale(1.0/h19.Integral());

            TCanvas c("c_rad_sys_distribution","",1050,760);
            style_canvas(c);

            h10.SetLineWidth(3);
            h10.SetLineColor(kBlue+1);
            h10.SetFillStyle(0);

            h19.SetLineWidth(3);
            h19.SetLineColor(kRed+1);
            h19.SetLineStyle(2);
            h19.SetFillStyle(0);

            const double ymax=1.25*std::max(h10.GetMaximum(),h19.GetMaximum());

            h10.SetMaximum(ymax>0.0 ? ymax : 1.0);
            h10.GetXaxis()->SetTitle("Relative radiative systematic (%)");
            h10.GetYaxis()->SetTitle("Fraction of populated analysis bins");
            h10.GetXaxis()->SetTitleSize(0.045);
            h10.GetYaxis()->SetTitleSize(0.045);
            h10.GetXaxis()->SetLabelSize(0.039);
            h10.GetYaxis()->SetLabelSize(0.039);
            h10.GetYaxis()->SetTitleOffset(1.25);

            h10.Draw("HIST");
            if(!all19.empty()) h19.Draw("HIST SAME");

            TLine med10(
                quantile(all10,0.50),0.0,
                quantile(all10,0.50),h10.GetMaximum()
            );
            med10.SetLineColor(kBlue+1);
            med10.SetLineStyle(3);
            med10.SetLineWidth(2);
            med10.Draw();

            if(!all19.empty()){
                TLine med19(
                    quantile(all19,0.50),0.0,
                    quantile(all19,0.50),h10.GetMaximum()
                );
                med19.SetLineColor(kRed+1);
                med19.SetLineStyle(3);
                med19.SetLineWidth(2);
                med19.Draw();
            }

            TLegend leg(0.57,0.72,0.93,0.86);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextSize(0.031);
            leg.AddEntry(&h10,"10.6 GeV: inherited pass-1","l");
            if(!all19.empty()){
                leg.AddEntry(
                    &h19,
                    "Sp19 10.2 GeV: 1.05 #times fractional",
                    "l"
                );
            }
            leg.Draw();

            note_label("Bin-by-bin radiative-correction systematic");

            c.SaveAs(
                (fs::path(output_dir)/"radiative_systematic_distribution.png")
                    .string().c_str()
            );
        }

        {
            struct Bucket {
                double lo=0.0;
                double hi=0.0;
                std::vector<double> ten6;
                std::vector<double> sp19;
            };

            std::map<std::pair<double,double>,Bucket> buckets;

            for(const auto& row:t.rows){
                const double lo=cell_number(t,row,"xBmin");
                const double hi=cell_number(t,row,"xBmax");
                if(!std::isfinite(lo) || !std::isfinite(hi)) continue;

                auto& b=buckets[{lo,hi}];
                b.lo=lo;
                b.hi=hi;

                const RelPair r=relative_uncertainties(t,row);
                if(std::isfinite(r.ten6)) b.ten6.push_back(r.ten6);
                if(std::isfinite(r.sp19)) b.sp19.push_back(r.sp19);
            }

            TGraphAsymmErrors g10;
            TGraphAsymmErrors g19;

            int i10=0;
            int i19=0;
            double ymax=0.0;

            for(const auto& kv:buckets){
                const Bucket& b=kv.second;
                const double x=0.5*(b.lo+b.hi);
                const double ex=0.5*(b.hi-b.lo);

                if(!b.ten6.empty()){
                    const double med=quantile(b.ten6,0.50);
                    const double p16=quantile(b.ten6,0.16);
                    const double p84=quantile(b.ten6,0.84);

                    g10.SetPoint(i10,x,med);
                    g10.SetPointError(i10,ex,ex,med-p16,p84-med);
                    ymax=std::max(ymax,p84);
                    ++i10;
                }

                if(!b.sp19.empty()){
                    const double med=quantile(b.sp19,0.50);
                    const double p16=quantile(b.sp19,0.16);
                    const double p84=quantile(b.sp19,0.84);

                    g19.SetPoint(i19,x,med);
                    g19.SetPointError(i19,ex,ex,med-p16,p84-med);
                    ymax=std::max(ymax,p84);
                    ++i19;
                }
            }

            TCanvas c("c_rad_sys_xb","",1050,760);
            style_canvas(c);

            TH1D frame("h_rad_sys_xb_frame","",100,0.04,0.62);
            frame.SetMinimum(0.0);
            frame.SetMaximum(std::max(5.0,1.18*ymax));
            frame.GetXaxis()->SetTitle("x_{B}");
            frame.GetYaxis()->SetTitle("Median radiative systematic (%)");
            frame.GetXaxis()->SetTitleSize(0.045);
            frame.GetYaxis()->SetTitleSize(0.045);
            frame.GetXaxis()->SetLabelSize(0.039);
            frame.GetYaxis()->SetLabelSize(0.039);
            frame.GetYaxis()->SetTitleOffset(1.22);
            frame.Draw();

            g10.SetMarkerStyle(20);
            g10.SetMarkerSize(1.0);
            g10.SetMarkerColor(kBlue+1);
            g10.SetLineColor(kBlue+1);
            g10.SetLineWidth(2);
            g10.Draw("PZ SAME");

            g19.SetMarkerStyle(24);
            g19.SetMarkerSize(1.0);
            g19.SetMarkerColor(kRed+1);
            g19.SetLineColor(kRed+1);
            g19.SetLineWidth(2);
            if(i19>0) g19.Draw("PZ SAME");

            TLegend leg(0.57,0.72,0.93,0.86);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextSize(0.031);
            leg.AddEntry(&g10,"10.6 GeV: inherited pass-1","pe");
            if(i19>0){
                leg.AddEntry(
                    &g19,
                    "Sp19 10.2 GeV: 1.05 #times fractional",
                    "pe"
                );
            }
            leg.Draw();

            note_label("Kinematic dependence of the radiative systematic");

            TLatex note;
            note.SetNDC();
            note.SetTextFont(42);
            note.SetTextSize(0.025);
            note.DrawLatex(
                0.13,0.88,
                "Points: median over populated (Q^{2}, |t|, #phi) bins; bars: central 68% interval"
            );

            c.SaveAs(
                (fs::path(output_dir)/"radiative_systematic_vs_xB.png")
                    .string().c_str()
            );
        }

        {
            struct P {
                double phi=0.0;
                double ten6=0.0;
                double sp19=0.0;
                bool has19=false;
            };

            std::vector<P> points;

            for(const auto& row:t.rows){
                const double xb0=cell_number(t,row,"xBmin");
                const double xb1=cell_number(t,row,"xBmax");
                const double q0=cell_number(t,row,"Q2min");
                const double q1=cell_number(t,row,"Q2max");
                const double t0=cell_number(t,row,"t_abs_min");
                const double t1=cell_number(t,row,"t_abs_max");

                if(!approx(xb0,0.204) || !approx(xb1,0.268) ||
                   !approx(q0,1.912) || !approx(q1,2.510) ||
                   !approx(t0,0.250) || !approx(t1,0.400)){
                    continue;
                }

                const double p0=cell_number(t,row,"phimin");
                const double p1=cell_number(t,row,"phimax");
                if(!std::isfinite(p0) || !std::isfinite(p1)) continue;

                const RelPair r=relative_uncertainties(t,row);
                if(!std::isfinite(r.ten6)) continue;

                P p;
                p.phi=0.5*(p0+p1);
                p.ten6=r.ten6;
                p.has19=std::isfinite(r.sp19);
                p.sp19=p.has19 ? r.sp19 : 0.0;
                points.push_back(p);
            }

            std::sort(
                points.begin(),points.end(),
                [](const P& a,const P& b){ return a.phi<b.phi; }
            );

            if(!points.empty()){
                TGraphErrors g10;
                TGraphErrors g19;

                double ymin=std::numeric_limits<double>::infinity();
                double ymax=0.0;
                int i10=0;
                int i19=0;

                for(const auto& p:points){
                    g10.SetPoint(i10,p.phi,p.ten6);
                    g10.SetPointError(i10,0.0,0.0);
                    ymin=std::min(ymin,p.ten6);
                    ymax=std::max(ymax,p.ten6);
                    ++i10;

                    if(p.has19){
                        g19.SetPoint(i19,p.phi,p.sp19);
                        g19.SetPointError(i19,0.0,0.0);
                        ymin=std::min(ymin,p.sp19);
                        ymax=std::max(ymax,p.sp19);
                        ++i19;
                    }
                }

                TCanvas c("c_rad_sys_phi","",1050,760);
                style_canvas(c);

                const double span=std::max(0.25,ymax-ymin);

                TH1D frame("h_rad_sys_phi_frame","",100,0.0,360.0);
                frame.SetMinimum(std::max(0.0,ymin-0.25*span));
                frame.SetMaximum(ymax+0.35*span);
                frame.GetXaxis()->SetTitle("#phi (deg)");
                frame.GetYaxis()->SetTitle("Radiative systematic (%)");
                frame.GetXaxis()->SetTitleSize(0.045);
                frame.GetYaxis()->SetTitleSize(0.045);
                frame.GetXaxis()->SetLabelSize(0.039);
                frame.GetYaxis()->SetLabelSize(0.039);
                frame.GetYaxis()->SetTitleOffset(1.20);
                frame.Draw();

                g10.SetMarkerStyle(20);
                g10.SetMarkerSize(1.0);
                g10.SetMarkerColor(kBlue+1);
                g10.SetLineColor(kBlue+1);
                g10.SetLineWidth(2);
                g10.Draw("PL SAME");

                g19.SetMarkerStyle(24);
                g19.SetMarkerSize(1.0);
                g19.SetMarkerColor(kRed+1);
                g19.SetLineColor(kRed+1);
                g19.SetLineWidth(2);
                if(i19>0) g19.Draw("PL SAME");

                TLegend leg(0.57,0.70,0.93,0.84);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextSize(0.031);
                leg.AddEntry(&g10,"10.6 GeV: inherited pass-1","pl");
                if(i19>0){
                    leg.AddEntry(
                        &g19,
                        "Sp19 10.2 GeV: 1.05 #times fractional",
                        "pl"
                    );
                }
                leg.Draw();

                note_label(
                    "Representative #phi dependence of the radiative systematic"
                );

                TLatex kin;
                kin.SetNDC();
                kin.SetTextFont(42);
                kin.SetTextSize(0.025);
                kin.DrawLatex(
                    0.13,0.875,
                    "0.204 < x_{B} < 0.268, 1.912 < Q^{2} < 2.510 GeV^{2}, 0.250 < |t| < 0.400 GeV^{2}"
                );

                c.SaveAs(
                    (fs::path(output_dir)/"radiative_systematic_phi_example.png")
                        .string().c_str()
                );
            }
        }

        std::cout
            << "[radiative-systematics] Wrote analysis-note outputs to "
            << output_dir << "\n";

        return true;

    } catch(const std::exception& e){
        std::cerr
            << "[radiative-systematics] ERROR: "
            << e.what() << "\n";
        return false;
    }
}
