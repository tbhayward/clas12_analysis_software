#include "systematic_projection_plots.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TPad.h>
#include <TH1D.h>

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
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string> > rows;
    std::unordered_map<std::string, int> index;
};

struct Component {
    std::string label;
    int marker = 20;
    int color = 1;
    int line_style = 1;
};

struct VariableSpec {
    std::string name;
    std::string x_title;
    std::string min_col;
    std::string max_col;
};

struct BinAccumulator {
    double x_sum = 0.0;
    int x_count = 0;
    std::vector<std::vector<double> > relative_values;
};

static std::string trim(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;
    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
    return s.substr(b, e - b);
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == ',' && !in_quotes) {
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
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Could not open CSV: " + path);

    CsvTable t;
    std::string line;
    if (!std::getline(fin, line)) throw std::runtime_error("Empty CSV: " + path);

    t.header = split_csv_line(line);
    for (int i = 0; i < (int)t.header.size(); ++i)
        t.index[t.header[(size_t)i]] = i;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        auto row = split_csv_line(line);
        row.resize(t.header.size());
        t.rows.push_back(std::move(row));
    }
    return t;
}

static double number(const std::string& s) {
    // Systematics modules rewrite the same CSV several times.  Be deliberately
    // tolerant of harmless wrapper characters that can remain around scalar
    // cells after those passes, while still requiring one unambiguous number.
    std::string x = trim(s);
    while (!x.empty() &&
           (x.front() == '"' || x.front() == '\'' ||
            x.front() == '(' || x.front() == '[' || x.front() == '{')) {
        x = trim(x.substr(1));
    }
    while (!x.empty() &&
           (x.back() == '"' || x.back() == '\'' ||
            x.back() == ')' || x.back() == ']' || x.back() == '}')) {
        x = trim(x.substr(0, x.size() - 1));
    }

    if (x.empty()) return std::numeric_limits<double>::quiet_NaN();

    char* end = nullptr;
    const double v = std::strtod(x.c_str(), &end);
    if (end == x.c_str()) return std::numeric_limits<double>::quiet_NaN();

    while (*end != '\0' && std::isspace((unsigned char)*end)) ++end;
    return *end == '\0' ? v : std::numeric_limits<double>::quiet_NaN();
}

static std::string cell(const CsvTable& t,
                        const std::vector<std::string>& row,
                        const std::string& col) {
    const auto it = t.index.find(col);
    if (it == t.index.end()) return "";
    return row[(size_t)it->second];
}

static bool parse_tuple_first(const std::string& raw, double& value) {
    // Accept the normal "(value,stat,sys)" representation as well as harmless
    // quote nesting left by repeated CSV read/write stages.  We only need the
    // first field here.
    std::string s = trim(raw);
    while (!s.empty() && (s.front() == '"' || s.front() == '\''))
        s = trim(s.substr(1));
    while (!s.empty() && (s.back() == '"' || s.back() == '\''))
        s = trim(s.substr(0, s.size() - 1));

    if (!s.empty() && s.front() == '(') s = trim(s.substr(1));
    if (!s.empty() && s.back() == ')')  s = trim(s.substr(0, s.size() - 1));

    const size_t comma = s.find(',');
    value = number(comma == std::string::npos ? s : s.substr(0, comma));
    return std::isfinite(value);
}

static double median(std::vector<double> values) {
    values.erase(
        std::remove_if(values.begin(), values.end(),
                       [](double x){ return !std::isfinite(x); }),
        values.end());
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();

    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    if (n % 2U) return values[n / 2U];
    return 0.5 * (values[n/2U - 1U] + values[n/2U]);
}

static std::vector<Component> components() {
    return {
        {"#pi^{0} subtraction", 20, kBlue + 1, 1},
        {"Acceptance",          21, kRed + 1, 1},
        {"F_{rad}",             22, kGreen + 2, 1},
        {"F_{bin}",             23, kMagenta + 1, 1},
        {"Exclusivity",         24, kOrange + 7, 1},
        {"Fiducial",            25, kCyan + 2, 1},
        {"Point-to-point total",29, kBlack, 1}
    };
}

static std::vector<VariableSpec> variables() {
    return {
        {"xB",  "x_{B}",           "xBmin",      "xBmax"},
        {"Q2",  "Q^{2} (GeV^{2})", "Q2min",      "Q2max"},
        {"t",   "|t| (GeV^{2})",   "t_abs_min",  "t_abs_max"},
        {"phi", "#phi (deg)",      "phimin",     "phimax"}
    };
}

static bool row_x(const CsvTable& table,
                  const std::vector<std::string>& row,
                  const VariableSpec& var,
                  double& x,
                  std::string& key) {
    const double lo = number(cell(table, row, var.min_col));
    const double hi = number(cell(table, row, var.max_col));
    if (!std::isfinite(lo) || !std::isfinite(hi)) return false;

    x = 0.5 * (lo + hi);
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(8) << lo << '|' << hi;
    key = ss.str();
    return true;
}

static bool finite_fraction(double x) {
    return std::isfinite(x) && x >= 0.0;
}

// Return the eight point-to-point fractional uncertainties for one CSV row.
// Values are returned as fractions, not percentages.
static std::vector<double> row_component_fractions(
    const CsvTable& table,
    const std::vector<std::string>& row,
    bool sp19) {

    double xs10 = std::numeric_limits<double>::quiet_NaN();
    double xs19 = std::numeric_limits<double>::quiet_NaN();

    parse_tuple_first(
        cell(table,row,"normed cross sections, ep->epg, exp, 10.6 GeV, unpol"),
        xs10);
    parse_tuple_first(
        cell(table,row,"normed cross sections, ep->epg, exp, Sp19 Inb, unpol"),
        xs19);

    const double xs = sp19 ? xs19 : xs10;
    if (!std::isfinite(xs) || xs == 0.0)
        return std::vector<double>(7,std::numeric_limits<double>::quiet_NaN());

    auto frac_from_abs_10p6 = [&](const std::string& col)->double {
        if (!std::isfinite(xs10) || xs10 == 0.0)
            return std::numeric_limits<double>::quiet_NaN();
        const double a = number(cell(table,row,col));
        return std::isfinite(a) ? std::fabs(a/xs10)
                                : std::numeric_limits<double>::quiet_NaN();
    };

    std::vector<double> f(7,std::numeric_limits<double>::quiet_NaN());

    if (!sp19) {
        f[0] = number(cell(table,row,"pi0 subtraction sys frac, 10.6 GeV"));
        f[1] = frac_from_abs_10p6("Syst. err (Acceptance)");
        f[2] = frac_from_abs_10p6("Syst.err (Frad)");
        f[3] = frac_from_abs_10p6("Syst.err (Fbin)");
        f[4] = frac_from_abs_10p6("Syst. err (exclusivity cuts)");
        f[5] = frac_from_abs_10p6("Syst. err (fiducial cuts)");
        // Proton-efficiency uncertainty is an overall normalization source,
        // so it is intentionally absent from the point-to-point projection.
        f[6] = frac_from_abs_10p6("Syst. err (point-to-point total)");

        // The displayed total is mathematically just the quadrature sum of
        // these six point-to-point components.  Reconstruct it if a downstream
        // CSV rewrite left the stored total temporarily unreadable.
        if (!finite_fraction(f[6])) {
            double sum2 = 0.0;
            bool complete = true;
            for (int j = 0; j < 6; ++j) {
                if (!finite_fraction(f[(size_t)j])) {
                    complete = false;
                    break;
                }
                sum2 += f[(size_t)j] * f[(size_t)j];
            }
            if (complete) f[6] = std::sqrt(sum2);
        }
        return f;
    }

    // Dedicated Sp19 quantities are used where available.
    f[0] = number(cell(
        table,row,"pi0 subtraction sys frac, Sp19 Inb (10.2 GeV)"));

    const double frad_abs = number(cell(
        table,row,"Syst.err (Frad), Sp19 Inb (10.2 GeV)"));
    if (std::isfinite(frad_abs))
        f[2] = std::fabs(frad_abs/xs19);

    // Acceptance now has a dedicated Sp19 result from the pass-2
    // DATA-reweighting + transfer-closure study.  Fbin, exclusivity and
    // fiducial retain the established bin-wise fractional transfer from the
    // corresponding 10.6-GeV bins.
    f[1] = number(cell(
        table,row,"acceptance reweighting conservative sys frac, Sp19 Inb"));
    const double fbin_sp_abs = number(cell(
        table,row,"Syst.err (Fbin), Sp19 Inb (10.2 GeV)"));
    if (std::isfinite(fbin_sp_abs))
        f[3] = std::fabs(fbin_sp_abs/xs19);

    // Exclusivity/fiducial retain the 10.6-GeV fractional transfer.  If the
    // source term is exactly zero, its Sp19 contribution is also exactly zero
    // even when that row has no combined 10.6-GeV cross section.
    const auto transferred_or_zero = [&](const std::string& col)->double {
        const double a = number(cell(table,row,col));
        if (!std::isfinite(a))
            return std::numeric_limits<double>::quiet_NaN();
        if (a == 0.0)
            return 0.0;
        return frac_from_abs_10p6(col);
    };
    f[4] = transferred_or_zero("Syst. err (exclusivity cuts)");
    f[5] = transferred_or_zero("Syst. err (fiducial cuts)");
    // Proton-efficiency uncertainty is an overall normalization source,
    // so it is intentionally absent from the point-to-point projection.

    // Prefer the dedicated production Sp19 total written by main_systematics.
    // Reconstruct it from the six displayed fractions only as a backward-safe
    // fallback for older CSVs.
    const double ptp_sp_abs = number(cell(
        table,row,"Syst. err (point-to-point total), Sp19 Inb (10.2 GeV)"));
    if (std::isfinite(ptp_sp_abs)) {
        f[6] = std::fabs(ptp_sp_abs/xs19);
    } else {
        double sum2 = 0.0;
        bool complete = true;
        for (int j=0;j<6;++j) {
            if (!finite_fraction(f[(size_t)j])) {
                complete = false;
                break;
            }
            sum2 += f[(size_t)j]*f[(size_t)j];
        }
        if (complete) f[6] = std::sqrt(sum2);
    }

    return f;
}

struct ProjectionPoint {
    double x = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> medians;
};

static std::vector<ProjectionPoint> build_projection(
    const CsvTable& table,
    const VariableSpec& var,
    bool sp19) {

    const auto comps = components();
    std::map<std::string,BinAccumulator> bins;
    size_t n_bad_total = 0;
    size_t n_bad_x = 0;
    size_t n_used = 0;

    for (const auto& row : table.rows) {
        const auto fractions = row_component_fractions(table,row,sp19);

        // Require a valid total for this energy before using the row.
        // fractions[6] is the point-to-point total; fractions[6] is the
        // proton-efficiency component.  The old index-6 check accidentally
        // rejected every row when the legacy proton-fraction column was absent.
        if (!finite_fraction(fractions[6])) {
            ++n_bad_total;
            continue;
        }

        double x = 0.0;
        std::string key;
        if (!row_x(table,row,var,x,key)) {
            ++n_bad_x;
            continue;
        }
        ++n_used;

        auto& b = bins[key];
        if (b.relative_values.empty())
            b.relative_values.resize(comps.size());

        b.x_sum += x;
        ++b.x_count;

        for (size_t j=0;j<comps.size();++j) {
            if (finite_fraction(fractions[j]))
                b.relative_values[j].push_back(100.0*fractions[j]);
        }
    }

    std::vector<ProjectionPoint> out;
    for (const auto& kv:bins) {
        const auto& b=kv.second;
        if (b.x_count<=0) continue;

        ProjectionPoint p;
        p.x=b.x_sum/b.x_count;
        p.medians.resize(comps.size(),
                         std::numeric_limits<double>::quiet_NaN());
        for (size_t j=0;j<comps.size();++j)
            p.medians[j]=median(b.relative_values[j]);

        out.push_back(std::move(p));
    }

    std::sort(out.begin(),out.end(),
              [](const ProjectionPoint& a,const ProjectionPoint& b){
                  return a.x<b.x;
              });

    std::cout << "[systematic-projections] "
              << (sp19 ? "10.2 GeV" : "10.6 GeV")
              << " " << var.name
              << ": rows used=" << n_used
              << ", invalid total/xs=" << n_bad_total
              << ", invalid x-bin=" << n_bad_x
              << ", projected bins=" << out.size() << "\n";

    return out;
}

static double panel_ymax(const std::vector<ProjectionPoint>& points) {
    double ymax=0.0;
    for (const auto& p:points)
        for (double y:p.medians)
            if (std::isfinite(y)) ymax=std::max(ymax,y);

    // Linear scale with enough room for data and a legend above the frame.
    if (!(ymax>0.0)) return 10.0;
    return std::max(5.0,1.18*ymax);
}

static void style_pad(TPad* pad,bool left) {
    pad->SetLeftMargin(left ? 0.13 : 0.105);
    pad->SetRightMargin(0.03);
    pad->SetBottomMargin(0.14);
    pad->SetTopMargin(0.20);
    pad->SetGridx(false);
    pad->SetGridy(true);
    pad->SetTicks(1,1);
    pad->SetLogy(false);
}

static void draw_energy_panel(
    TPad* pad,
    const std::vector<ProjectionPoint>& points,
    const VariableSpec& var,
    const std::string& energy_label,
    bool draw_y_title,
    bool draw_legend) {

    style_pad(pad,draw_y_title);

    const auto comps=components();
    if (points.empty()) return;

    double xmin=points.front().x;
    double xmax=points.back().x;
    if (!(xmax>xmin)) {
        xmin-=0.5;
        xmax+=0.5;
    } else {
        const double dx=0.04*(xmax-xmin);
        xmin-=dx;
        xmax+=dx;
    }

    TH1D frame(
        ("h_syst_projection_"+var.name+"_"+energy_label).c_str(),
        "",
        100,xmin,xmax);
    frame.SetStats(0);
    frame.SetMinimum(0.0);
    frame.SetMaximum(panel_ymax(points));
    frame.GetXaxis()->SetTitle(var.x_title.c_str());
    frame.GetYaxis()->SetTitle(
        draw_y_title ? "Median relative systematic uncertainty (%)" : "");
    frame.GetXaxis()->SetTitleFont(42);
    frame.GetYaxis()->SetTitleFont(42);
    frame.GetXaxis()->SetLabelFont(42);
    frame.GetYaxis()->SetLabelFont(42);
    frame.GetXaxis()->SetTitleSize(0.049);
    frame.GetYaxis()->SetTitleSize(0.045);
    frame.GetXaxis()->SetLabelSize(0.040);
    frame.GetYaxis()->SetLabelSize(0.038);
    frame.GetXaxis()->SetTitleOffset(1.04);
    frame.GetYaxis()->SetTitleOffset(draw_y_title ? 1.35 : 1.0);
    frame.DrawCopy();

    std::vector<std::unique_ptr<TGraph> > graphs;
    graphs.reserve(comps.size());

    for (size_t j=0;j<comps.size();++j) {
        std::vector<double> x,y;
        for (const auto& p:points) {
            const double yy=p.medians[j];
            if (!std::isfinite(yy)) continue;
            x.push_back(p.x);
            y.push_back(yy);
        }
        if (x.empty()) continue;

        auto g=std::make_unique<TGraph>(
            static_cast<int>(x.size()),x.data(),y.data());
        g->SetMarkerStyle(comps[j].marker);
        g->SetMarkerColor(comps[j].color);
        g->SetLineColor(comps[j].color);
        g->SetLineStyle(comps[j].line_style);
        g->SetMarkerSize(j==6 ? 1.15 : 0.85);
        g->SetLineWidth(j==6 ? 4 : 2);
        g->DrawClone("LP SAME");
        graphs.push_back(std::move(g));
    }

    TLatex lab;
    lab.SetNDC();
    lab.SetTextFont(42);
    lab.SetTextSize(0.035);
    lab.DrawLatex(0.16,0.845,energy_label.c_str());

    if (draw_legend) {
        // The legend sits in the reserved top band and therefore never covers data.
        TLegend leg(0.11,0.875,0.97,0.985);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextFont(42);
        leg.SetTextSize(0.025);
        leg.SetNColumns(4);
        leg.SetMargin(0.16);

        // Add entries using temporary style-holder graphs.
        std::vector<std::unique_ptr<TGraph> > holders;
        for (size_t j=0;j<comps.size();++j) {
            auto h=std::make_unique<TGraph>();
            h->SetMarkerStyle(comps[j].marker);
            h->SetMarkerColor(comps[j].color);
            h->SetLineColor(comps[j].color);
            h->SetLineWidth(j==6 ? 4 : 2);
            leg.AddEntry(h.get(),comps[j].label.c_str(),"lp");
            holders.push_back(std::move(h));
        }
        leg.DrawClone();
    }

    pad->RedrawAxis();
}

static bool make_one(const CsvTable& table,
                     const VariableSpec& var,
                     const std::string& output_dir) {

    const auto p10=build_projection(table,var,false);
    const auto p19=build_projection(table,var,true);
    if (p10.empty() && p19.empty()) return false;

    TCanvas c(("c_point_to_point_"+var.name).c_str(),"",1500,720);
    c.Divide(2,1,0.002,0.0);

    draw_energy_panel(
        static_cast<TPad*>(c.cd(1)),p10,var,"10.6 GeV combined",true,true);
    draw_energy_panel(
        static_cast<TPad*>(c.cd(2)),p19,var,"10.2 GeV Sp19 Inb",false,false);

    c.cd(0);
    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextAlign(22);
    title.SetTextSize(0.023);
    const std::string main_title=
        "Kinematic dependence of point-to-point systematic uncertainties";
    title.DrawLatex(0.50,0.985,main_title.c_str());

    const std::string png=
        output_dir+"/point_to_point_systematics_vs_"+var.name+".png";
    c.SaveAs(png.c_str());

    // Machine-readable medians for the note and later cross checks.
    const std::string csv=
        output_dir+"/point_to_point_systematics_vs_"+var.name+".csv";
    std::ofstream out(csv);
    out<<"energy,x";
    for (const auto& comp:components())
        out<<','<<comp.label<<"_relative_percent_median";
    out<<'\n';

    auto write_points=[&](const std::string& e,
                          const std::vector<ProjectionPoint>& points) {
        for (const auto& p:points) {
            out<<e<<','<<std::setprecision(12)<<p.x;
            for (double y:p.medians) out<<','<<y;
            out<<'\n';
        }
    };
    write_points("10.6 GeV",p10);
    write_points("10.2 GeV",p19);

    return true;
}

} // namespace

bool make_systematic_projection_plots(
    const std::string& csv_path,
    const std::string& output_dir) {

    try {
        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);
        fs::create_directories(output_dir);

        const CsvTable table=read_csv(csv_path);

        int made=0;
        for (const auto& var:variables()) {
            if (make_one(table,var,output_dir)) ++made;
        }

        std::cout<<"[systematic-projections] Made "<<made
                 <<" note-quality point-to-point projection canvases in "
                 <<output_dir<<"\n";

        return made==(int)variables().size();
    } catch (const std::exception& e) {
        std::cerr<<"[systematic-projections] ERROR: "<<e.what()<<"\n";
        return false;
    }
}
