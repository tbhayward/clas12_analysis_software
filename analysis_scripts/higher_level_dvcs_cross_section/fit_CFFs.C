/**
 * fit_CFFs.C
 * ──────────
 * program to fit DVCS Compton Form Factors (CFFs) imaginary then real parts
 *
 * Strategies:
 *   1) fit only Im–parts → BSA data
 *   2) simultaneous: fit Im–parts + renormReal + subtraction constants → BSA + xsec data
 *
 * Usage:
 *   ./fit_CFFs --strategy <1|2> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1>
 *             [--constraint <0|1|2>] [--scale <0|1>] [--input <BSA_file> [XSEC_file]] [--plot-fits]
 *
 *   If only BSA_file is given, xsec uses the built-in default.  If neither
 *   is given, both defaults are used.
 *
 * Compile:
 *   g++ -O2 fit_CFFs.C `root-config --cflags --libs` -lMinuit -o fit_CFFs
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <map>
#include <functional>
#include <limits>

// ROOT headers
#include "TMinuit.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TStyle.h"

// pull in full BMK_DVCS + CFF code, with globals
#include "DVCS_xsec.C"

// extern flags & renormalizations
extern bool   hasH, hasHt, hasE, hasEt;
extern double renormImag, renormReal;

// ----------------------------------------------------------------------------
/* GPD–H defaults (in DVCS_xsec.C) */
extern double r_H,   n_H,   alpha0_H,   alpha1_H,   b_H,   M2_H,   P_H;
/* GPD–Htilde */
extern double r_Ht,  n_Ht,  alpha0_Ht,  alpha1_Ht,  b_Ht,  M2_Ht,  P_Ht;
/* GPD–E */
extern double r_E,   n_E,   alpha0_E,   alpha1_E,   b_E,   M2_E,   P_E;
/* GPD–Etilde */
extern double r_Et,  n_Et,  alpha0_Et,  alpha1_Et,  b_Et,  M2_Et,  P_Et;

/* subtraction‐relation constants (fit in strategy 2) */
extern double C0_H,    MD2_H,    lambda_H;
extern double C0_Ht,   MD2_Ht,   lambda_Ht;
extern double C0_E,    MD2_E,    lambda_E;
extern double C0_Et,   MD2_Et,   lambda_Et;

// ----------------------------------------------------------------------------
// Control flags & data containers
static int   gStrategy    = 0;
static int   gStage       = 1;  // 1 = Im-only, 0 = global simultaneous
static int   gConstraint  = 0;  // 0,1,2
static bool  gPlotBinFits = true;
static int   gScale       = 0;  // 0 off, 1 on

// default data files:
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static std::string gXsFile  = "imports/rga_pass1_xsec_2018.txt";

struct DataPoint { double phi, Q2, xB, t, Eb, A, sigA; };
// NOTE: in the input files, 't' is the physical Mandelstam t and is NEGATIVE.
//       (So -t is positive.)
static std::vector<DataPoint> bsaData, xsData;
static std::vector<std::vector<DataPoint>> keptBins;

// per-bin averages & fits (averages use weights 1/sigA^2)
static std::vector<double> bin_xB, bin_Q2, bin_t, bin_Eb;
static std::vector<double> bin_A, bin_dA, bin_redChi2;
static std::vector<int>    bin_M;
static int     Nbins          = 0;
static double  reducedAmpChi2 = 0.0;

// ----------------------------------------------------------------------------
// Helpers: parse command line, load & bin data, plot per-bin fits
void PlotBinFit(int ibin, const std::string &ts);

void parse_args(int argc, char** argv){
    static struct option opts[] = {
      {"strategy",   required_argument,nullptr,'s'},
      {"H",          required_argument,nullptr,'h'},
      {"Ht",         required_argument,nullptr,'t'},
      {"E",          required_argument,nullptr,'e'},
      {"Et",         required_argument,nullptr,'x'},
      {"constraint", required_argument,nullptr,'C'},
      {"scale",      required_argument,nullptr,'S'},
      {"input",      required_argument,nullptr,'i'},
      {"plot-fits",  no_argument,      nullptr,'p'},
      {nullptr,0,nullptr,0}
    };
    int c;
    while((c=getopt_long(argc,argv,"s:h:t:e:x:C:S:i:p",opts,nullptr))!=-1){
        switch(c){
          case 's': gStrategy   = std::atoi(optarg); break;
          case 'h': hasH        = std::atoi(optarg); break;
          case 't': hasHt       = std::atoi(optarg); break;
          case 'e': hasE        = std::atoi(optarg); break;
          case 'x': hasEt       = std::atoi(optarg); break;
          case 'C': gConstraint = std::atoi(optarg); break;
          case 'S': gScale      = std::atoi(optarg); break;
          case 'p': gPlotBinFits= true;              break;
          case 'i':
            // first arg → BSA file
            gBsaFile = optarg;
            // if there's another non-option arg, treat it as the xsec file:
            if(optind<argc && argv[optind][0] != '-') {
              gXsFile = argv[optind++];
            }
            break;
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" --strategy<1|2> -H<0|1> -Ht<0|1>"
                     <<" -E<0|1> -Et<0|1>"
                     <<" [--constraint<0|1|2>] [--scale<0|1>]"
                     <<" [--input <BSA_file> [XSEC_file]]"
                     <<" [--plot-fits]\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>2
     ||gConstraint<0||gConstraint>2
     ||(gScale!=0 && gScale!=1)){
        std::cerr<<"Invalid strategy, constraint, or scale flag\n";
        std::exit(1);
    }
}

static bool pass_base_cuts(const DataPoint& d){
    // Require -t > 0  and  -t < 1.0  (inputs provide t < 0)
    if (d.t >= 0.0)   return false; // -t must be > 0
    if (-d.t >= 1.0)  return false; // -t must be < 1
    // Basic sanity for xB
    if (d.xB <= 0.0 || d.xB >= 1.0) return false;
    return true;
}

static bool pass_constraint(int mode, const DataPoint& d){
    // mode 0: no extra constraint
    // mode 1: require (-t)/Q2 > 0.2
    // mode 2: require (-t)/Q2 > 0.2 AND (-t) < 0.45
    if (mode == 0) return true;
    const double mt_over_Q2 = (-d.t) / d.Q2;
    if (mode >= 1){
        if (mt_over_Q2 <= 0.2) return false;
    }
    if (mode == 2){
        if ((-d.t) >= 0.45) return false;
    }
    return true;
}

void LoadData(){
    auto read=[&](const std::string &fn, auto &v){
        std::ifstream in(fn);
        if(!in){ std::cerr<<"ERROR: cannot open "<<fn<<"\n"; std::exit(1); }
        std::string line;
        while(std::getline(in,line)){
            if(line.empty()||line[0]=='#') continue;
            std::istringstream iss(line);
            DataPoint d; iss>>d.phi>>d.Q2>>d.xB>>d.t>>d.Eb>>d.A>>d.sigA;
            // Inputs provide t < 0 (so -t > 0). Apply base window and constraints on -t.
            if(!pass_base_cuts(d))      continue;
            if(!pass_constraint(gConstraint,d)) continue;
            v.push_back(d);
        }
    };
    bsaData.clear(); xsData.clear();
    read(gBsaFile, bsaData);
    read(gXsFile,  xsData);
}

void BinBsaData(){
    keptBins.clear();
    bin_xB.clear(); bin_Q2.clear(); bin_t.clear(); bin_Eb.clear();
    bin_A.clear();  bin_dA.clear();  bin_redChi2.clear(); bin_M.clear();

    if(bsaData.empty()) return;
    size_t start=0;
    for(size_t i=1;i<=bsaData.size();++i){
        bool newbin = (i==bsaData.size()|| bsaData[i].phi < bsaData[i-1].phi);
        if(!newbin) continue;

        auto pts = std::vector<DataPoint>(bsaData.begin()+start,
                                          bsaData.begin()+i);
        start = i;

        const int M = pts.size();
        if(M<=3) continue; // need NDF > 0 for a 3-parameter fit

        TGraphErrors gr(M);
        for(int j=0;j<M;++j){
            gr.SetPoint(j, pts[j].phi, pts[j].A);
            gr.SetPointError(j,0,pts[j].sigA);
        }
        TF1 ftmp("ftmp",
          "[0] + [1]*sin(x*TMath::Pi()/180.)"
               "/(1 + [2]*cos(x*TMath::Pi()/180.))",
          0,360);
        ftmp.SetParameters(0.,0.5,0.);
        gr.Fit(&ftmp,"Q0R");
        const double Afit    = ftmp.GetParameter(1);
        const double dA      = ftmp.GetParError(1);
        const double Bfit    = ftmp.GetParameter(2);
        const double chi2    = ftmp.GetChisquare();
        const double redchi2 = (M>3) ? (chi2 / (M - 3)) : std::numeric_limits<double>::infinity();

        if(gConstraint==1 && (redchi2>=2.0 || std::fabs(Bfit)>=0.9))
            continue;
        if(gConstraint==2 && (redchi2>=2.0 || std::fabs(Bfit)>=0.9))
            continue;

        keptBins.push_back(pts);
        bin_M.push_back(M);
        bin_A.push_back(Afit);
        bin_dA.push_back(dA);
        bin_redChi2.push_back(redchi2);

        double sumw=0, Sx=0, Sq=0, St=0, Se=0;
        for(const auto &d: pts){
            const double w = 1./(d.sigA*d.sigA);
            sumw+=w; Sx+=w*d.xB; Sq+=w*d.Q2; St+=w*d.t; Se+=w*d.Eb;
        }
        bin_xB .push_back(Sx/sumw);
        bin_Q2 .push_back(Sq/sumw);
        bin_t  .push_back(St/sumw);  // this is (negative) t
        bin_Eb .push_back(Se/sumw);
    }

    Nbins = (int)bin_A.size();

    // χ² of amplitude-only (A) fits, using per-bin dA:
    double totChi2=0; int totDof=0;
    for(int k=0;k<Nbins;++k){
        totChi2 += std::pow(bin_A[k]/bin_dA[k],2);
        totDof  += (bin_M[k] - 1);
    }
    reducedAmpChi2 = totDof>0 ? totChi2/totDof : 0.0;
}

// ----------------------------------------------------------------------------
// Build parameter-lists for MINUIT
static std::vector<std::string> parNamesIm, parNamesRe, parNamesAll;

void build_par_listIm(){
    parNamesIm.clear();
    if(hasH )
      parNamesIm.insert(parNamesIm.end(),
        {"r_H","alpha0_H","alpha1_H","b_H","M2_H","P_H"});
    if(hasHt)
      parNamesIm.insert(parNamesIm.end(),
        {"r_Ht","alpha0_Ht","alpha1_Ht","b_Ht","M2_Ht","P_Ht"});
    if(hasE )
      parNamesIm.insert(parNamesIm.end(),
        {"r_E","alpha0_E","alpha1_E","b_E","M2_E","P_E"});
    if(hasEt)
      parNamesIm.insert(parNamesIm.end(),
        {"r_Et","alpha0_Et","alpha1_Et","b_Et","M2_Et","P_Et"});
}

void build_par_listRe(){
    parNamesRe.clear();
    parNamesRe.push_back("renormReal");
    if(hasH )
      parNamesRe.insert(parNamesRe.end(),
        {"C0_H","MD2_H","lambda_H"});
    if(hasHt)
      parNamesRe.insert(parNamesRe.end(),
        {"C0_Ht","MD2_Ht","lambda_Ht"});
    if(hasE )
      parNamesRe.insert(parNamesRe.end(),
        {"C0_E","MD2_E","lambda_E"});
    if(hasEt)
      parNamesRe.insert(parNamesRe.end(),
        {"C0_Et","MD2_Et","lambda_Et"});
}

// ─────────────────────────────────────────────────────────────────────────────
// USER-TWEAKABLE BOUNDS for Stage-1 Im-CFF parameters.
// By default: loose bounds (±1e6) and P_* fixed to 1.0.
// Edit the lo/hi/step entries below to constrain as you wish.
// If you want to *free* a P_* later, give it a real range (e.g. {0.8,1.2,1e-3}).
struct Bounds { double lo, hi, step; };
// static std::map<std::string, Bounds> gImBounds = {
//   // 𝓗
//   {"r_H",       {-1e6, +1e6, 1e-3}}, {"alpha0_H",  {-1e6, +1e6, 1e-3}},
//   {"alpha1_H",  {-1e6, +1e6, 1e-3}}, {"b_H",       {-1e6, +1e6, 1e-3}},
//   {"M2_H",      {-1e6, +1e6, 1e-3}}, {"P_H",       { 1.0,  1.0,  0.0}}, // fixed

//   // 𝓗̃
//   {"r_Ht",      {-1e6, +1e6, 1e-3}}, {"alpha0_Ht", {-1e6, +1e6, 1e-3}},
//   {"alpha1_Ht", {-1e6, +1e6, 1e-3}}, {"b_Ht",      {-1e6, +1e6, 1e-3}},
//   {"M2_Ht",     {-1e6, +1e6, 1e-3}}, {"P_Ht",      { 1.0,  1.0,  0.0}}, // fixed

//   // 𝓔
//   {"r_E",       {-1e6, +1e6, 1e-3}}, {"alpha0_E",  {-1e6, +1e6, 1e-3}},
//   {"alpha1_E",  {-1e6, +1e6, 1e-3}}, {"b_E",       {-1e6, +1e6, 1e-3}},
//   {"M2_E",      {-1e6, +1e6, 1e-3}}, {"P_E",       { 1.0,  1.0,  0.0}}, // fixed

//   // 𝓔̃
//   {"r_Et",      {-1e6, +1e6, 1e-3}}, {"alpha0_Et", {-1e6, +1e6, 1e-3}},
//   {"alpha1_Et", {-1e6, +1e6, 1e-3}}, {"b_Et",      {-1e6, +1e6, 1e-3}},
//   {"M2_Et",     {-1e6, +1e6, 1e-3}}, {"P_Et",      { 1.0,  1.0,  0.0}}, // fixed
// };
static std::map<std::string, Bounds> gImBounds = {
  // 𝓗
  {"r_H",       {-1e6, +1e6, 1e-3}}, {"alpha0_H",  {-1e6, +1e6, 1e-3}},
  {"alpha1_H",  {-1e6, +1e6, 1e-3}}, {"b_H",       {-1e6, +1e6, 1e-3}},
  {"M2_H",      {-1e6, +1e6, 1e-3}}, {"P_H",       { 1.0,  1.0,  0.0}}, // fixed

  // 𝓗̃
  {"r_Ht",      {-1e6, +1e6, 1e-3}}, {"alpha0_Ht", {-1e6, +1e6, 1e-3}},
  {"alpha1_Ht", {-1e6, +1e6, 1e-3}}, {"b_Ht",      {-1e6, +1e6, 1e-3}},
  {"M2_Ht",     {-1e6, +1e6, 1e-3}}, {"P_Ht",      { 1.0,  1.0,  0.0}}, // fixed

  // 𝓔
  {"r_E",       {-1e6, +1e6, 1e-3}}, {"alpha0_E",  {-1e6, +1e6, 1e-3}},
  {"alpha1_E",  {-1e6, +1e6, 1e-3}}, {"b_E",       {-1e6, +1e6, 1e-3}},
  {"M2_E",      {-1e6, +1e6, 1e-3}}, {"P_E",       { 1.0,  1.0,  0.0}}, // fixed

  // 𝓔̃
  {"r_Et",      {-1e6, +1e6, 1e-3}}, {"alpha0_Et", {-1e6, +1e6, 1e-3}},
  {"alpha1_Et", {-1e6, +1e6, 1e-3}}, {"b_Et",      {-1e6, +1e6, 1e-3}},
  {"M2_Et",     {-1e6, +1e6, 1e-3}}, {"P_Et",      { 1.0,  1.0,  0.0}}, // fixed
};
// ─────────────────────────────────────────────────────────────────────────────

// ----------------------------------------------------------------------------
// χ² function (no alpha0/alpha1 penalty terms)
void fcn(int&, double*, double &f, double *par, int){
    int ip=0;
    if(gStage==1){
        // assign Im-part parameters
        if(hasH){
          r_H       = par[ip++];
          alpha0_H  = par[ip++];
          alpha1_H  = par[ip++];
          b_H       = par[ip++];
          M2_H      = par[ip++];
          P_H       = par[ip++];
        }
        if(hasHt){
          r_Ht      = par[ip++];
          alpha0_Ht = par[ip++];
          alpha1_Ht = par[ip++];
          b_Ht      = par[ip++];
          M2_Ht     = par[ip++];
          P_Ht      = par[ip++];
        }
        if(hasE){
          r_E       = par[ip++];
          alpha0_E  = par[ip++];
          alpha1_E  = par[ip++];
          b_E       = par[ip++];
          M2_E      = par[ip++];
          P_E       = par[ip++];
        }
        if(hasEt){
          r_Et      = par[ip++];
          alpha0_Et = par[ip++];
          alpha1_Et = par[ip++];
          b_Et      = par[ip++];
          M2_Et     = par[ip++];
          P_Et      = par[ip++];
        }

        // χ² vs. BSA-derived amplitudes (used bins)
        double chi2=0;
        for(int k=0;k<Nbins;++k){
            BMK_DVCS dvcs(-1,1,0, bin_Eb[k],bin_xB[k],bin_Q2[k],bin_t[k],0.0);
            double modelA = dvcs.s1_I() / dvcs.c0_BH();
            double resid  = (bin_A[k] - modelA)/bin_dA[k];
            chi2 += resid*resid;
        }
        f = chi2;
    }
    else { // global simultaneous (strategy 2)
        // assign Im-part parameters
        if(hasH){
          r_H       = par[ip++];
          alpha0_H  = par[ip++];
          alpha1_H  = par[ip++];
          b_H       = par[ip++];
          M2_H      = par[ip++];
          P_H       = par[ip++];
        }
        if(hasHt){
          r_Ht      = par[ip++];
          alpha0_Ht = par[ip++];
          alpha1_Ht = par[ip++];
          b_Ht      = par[ip++];
          M2_Ht     = par[ip++];
          P_Ht      = par[ip++];
        }
        if(hasE){
          r_E       = par[ip++];
          alpha0_E  = par[ip++];
          alpha1_E  = par[ip++];
          b_E       = par[ip++];
          M2_E      = par[ip++];
          P_E       = par[ip++];
        }
        if(hasEt){
          r_Et      = par[ip++];
          alpha0_Et = par[ip++];
          alpha1_Et = par[ip++];
          b_Et      = par[ip++];
          M2_Et     = par[ip++];
          P_Et      = par[ip++];
        }
        // assign Re-part parameters
        renormReal = par[ip++];
        if(hasH){
          C0_H      = par[ip++];
          MD2_H     = par[ip++];
          lambda_H  = par[ip++];
        }
        if(hasHt){
          C0_Ht     = par[ip++];
          MD2_Ht    = par[ip++];
          lambda_Ht = par[ip++];
        }
        if(hasE){
          C0_E      = par[ip++];
          MD2_E     = par[ip++];
          lambda_E  = par[ip++];
        }
        if(hasEt){
          C0_Et     = par[ip++];
          MD2_Et    = par[ip++];
          lambda_Et = par[ip++];
        }

        // χ² total = χ²_BSA + χ²_xsec
        double chi2 = 0;

        // BSA part (used bins)
        for(int k=0;k<Nbins;++k){
            BMK_DVCS dvcs(-1,1,0, bin_Eb[k],bin_xB[k],bin_Q2[k],bin_t[k],0.0);
            double modelA = dvcs.s1_I() / dvcs.c0_BH();
            double resid  = (bin_A[k] - modelA)/bin_dA[k];
            chi2 += resid*resid;
        }
        // xsec part (apply same base/constraint cuts when loading)
        for(const auto &d: xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double xs    = dvcs.CrossSection();
            double resid = (d.A - renormReal*xs)/d.sigA;
            chi2 += resid*resid;
        }
        f = chi2;
    }
}

// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
    parse_args(argc, argv);

    // timestamp
    time_t now = time(nullptr);
    char tb[32];
    strftime(tb,sizeof(tb),"%Y%m%d_%H%M%S",localtime(&now));

    std::cout<<"\n=== Strategy="<<gStrategy
             <<"  H="<<hasH<<" Ht="<<hasHt
             <<"  E="<<hasE<<" Et="<<hasEt
             <<"  constraint="<<gConstraint
             <<"  scale="<<gScale
             <<"  input="<<gBsaFile
             <<"  plot-fits="<<(gPlotBinFits?"ON":"OFF")
             <<"  ts="<<tb<<" ===\n\n";

    // load and initial binning
    LoadData();
    BinBsaData();

    // compute average per-bin reduced chi2 before any amplitude scaling
    double avgBinChi2_before = bin_redChi2.empty()?0.0:
            std::accumulate(bin_redChi2.begin(), bin_redChi2.end(), 0.0) / bin_redChi2.size();
    std::cout<<"Average χ²/bin-fit before amplitude scaling = "<<avgBinChi2_before<<"\n";

    if(gScale==1){
        // scale the original BSA point uncertainties by sqrt(mean bin reduced chi2)
        double scale_amp = std::sqrt(avgBinChi2_before > 0 ? avgBinChi2_before : 1.0);
        std::cout<<"Scaling BSA point uncertainties by sqrt("<<avgBinChi2_before<<") = "<<scale_amp<<"\n";
        for(auto &d: bsaData){
            d.sigA *= scale_amp;
        }
        // redo binning / amplitude fits
        BinBsaData();
        double avgBinChi2_after = bin_redChi2.empty()?0.0:
            std::accumulate(bin_redChi2.begin(), bin_redChi2.end(), 0.0) / bin_redChi2.size();
        std::cout<<"Average χ²/bin-fit after amplitude scaling = "<<avgBinChi2_after<<"\n";
    }

    // report binned fits
    std::cout<<"Data points entering Im-fit:\n";
    for(int k=0;k<Nbins;++k){
      std::cout<<" Bin "<<(k+1)
               <<", Eb="<<bin_Eb[k]
               <<", xB="<<bin_xB[k]
               <<", Q2="<<bin_Q2[k]
               <<", -t="<<-bin_t[k]
               <<", A="<<bin_A[k]
               <<" ±"<<bin_dA[k]
               <<", χ²_red="<<bin_redChi2[k]
               <<"\n";
    }
    double avgBinChi2 = bin_redChi2.empty()?0.0:
                        std::accumulate(bin_redChi2.begin(),
                                        bin_redChi2.end(),0.0)
                        / bin_redChi2.size();
    std::cout<<"\n BSA bins="<<Nbins
             <<"  Avg χ²/bin-fit="<<avgBinChi2
             <<"  χ²_per_amp-fit="<<reducedAmpChi2<<"\n\n";

    if(gPlotBinFits){
      gStyle->SetOptStat(0);
      for(int ib=0; ib<Nbins; ++ib) PlotBinFit(ib,tb);
      std::cout<<"Wrote per-bin fits to output/plots/binned_fits/\n\n";
    }

    std::map<std::string,double> finalValMap, finalErrMap;
    double chi2_total=0;
    int    ndf_total=0;

    // We'll store intermediate results if scaling CFF-fit is turned on
    std::map<std::string,double> interimValMap, interimErrMap;
    double chi2_total1=0, ndf_total1=0;
    double chi2_total2=0, ndf_total2=0;

    if(gStrategy==1){
      // ─── Stage 1: Im-fit only ────────────────────────────────────────────────
      gStage = 1;
      build_par_listIm();
      int nim = parNamesIm.size();
      std::vector<double> imVal(nim), imErr(nim);
      double chi2_im, edm, errdef;
      int nv,nx,ic, ndf_im;

      auto run_im_fit = [&](){
        TMinuit minu(nim);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);

        for(int i=0;i<nim;++i){
          const std::string &nm = parNamesIm[i];

          // default bounds (very loose)
          double lo   = -1e6;
          double hi   = +1e6;
          double step = 1e-3;

          // override from user-tweakable table if present
          if(auto it = gImBounds.find(nm); it != gImBounds.end()){
            lo   = it->second.lo;
            hi   = it->second.hi;
            step = it->second.step;
          }

          // initial values from current globals (unchanged from your code)
          double init=0.0;
               if(nm=="r_H")        init = r_H;
          else if(nm=="r_Ht")       init = r_Ht;
          else if(nm=="r_E")        init = r_E;
          else if(nm=="r_Et")       init = r_Et;

          else if(nm=="alpha0_H")   init = alpha0_H;
          else if(nm=="alpha0_Ht")  init = alpha0_Ht;
          else if(nm=="alpha0_E")   init = alpha0_E;
          else if(nm=="alpha0_Et")  init = alpha0_Et;

          else if(nm=="alpha1_H")   init = alpha1_H;
          else if(nm=="alpha1_Ht")  init = alpha1_Ht;
          else if(nm=="alpha1_E")   init = alpha1_E;
          else if(nm=="alpha1_Et")  init = alpha1_Et;

          else if(nm=="b_H")        init = b_H;
          else if(nm=="b_Ht")       init = b_Ht;
          else if(nm=="b_E")        init = b_E;
          else if(nm=="b_Et")       init = b_Et;

          else if(nm=="M2_H")       init = M2_H;
          else if(nm=="M2_Ht")      init = M2_Ht;
          else if(nm=="M2_E")       init = M2_E;
          else if(nm=="M2_Et")      init = M2_Et;

          else if(nm=="P_H")        init = P_H;
          else if(nm=="P_Ht")       init = P_Ht;
          else if(nm=="P_E")        init = P_E;
          else if(nm=="P_Et")       init = P_Et;

          // define parameter with bounds
          minu.DefineParameter(i, nm.c_str(), init, step, lo, hi);

          // if step==0 (e.g. P_* fixed), fix it
          if(step==0.0 || lo==hi) minu.FixParameter(i);
        }

        std::cout<<" Stage1: fitting Im-CFF parameters…\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_im,edm,errdef,nv,nx,ic);
        for(int i=0;i<nim;++i) minu.GetParameter(i,imVal[i],imErr[i]);
        ndf_im = Nbins - nim;

        // capture
        interimValMap.clear(); interimErrMap.clear();
        for(int i=0;i<nim;++i){
          interimValMap[parNamesIm[i]] = imVal[i];
          interimErrMap[parNamesIm[i]] = imErr[i];
        }
        chi2_total1 = chi2_im;
        ndf_total1  = ndf_im;
      };

      // first CFF fit (Im-only)
      run_im_fit();
      std::cout<<"CFF fit (Im-only) χ²/ndf before scaling amplitude uncertainties = "
               <<chi2_total1<<"/"<<ndf_total1
               <<" = "<<(ndf_total1>0?chi2_total1/ndf_total1:0)<<"\n";

      // if scale flagged, scale bin_dA by sqrt(chi2/ndf) and redo Im-fit
      if(gScale==1){
          double scale_cff = std::sqrt((ndf_total1>0?chi2_total1/ndf_total1:1.0));
          std::cout<<"Scaling amplitude uncertainties (bin_dA) by sqrt("<<chi2_total1<<"/"<<ndf_total1<<") = "<<scale_cff<<"\n";
          for(auto &dA : bin_dA) dA *= scale_cff;

          // recompute reducedAmpChi2 with new bin_dA
          double totChi2=0; int totDof=0;
          for(int k=0;k<Nbins;++k){
              totChi2 += std::pow(bin_A[k]/bin_dA[k],2);
              totDof  += (bin_M[k]-1);
          }
          reducedAmpChi2 = totDof>0 ? totChi2/totDof : 0.0;

          // rerun fit
          run_im_fit();
          chi2_total2 = chi2_total1;
          ndf_total2  = ndf_total1;
          std::cout<<"CFF fit (Im-only) χ²/ndf after scaling amplitude uncertainties = "
                   <<chi2_total2<<"/"<<ndf_total2
                   <<" = "<<(ndf_total2>0?chi2_total2/ndf_total2:0)<<"\n";
          // final results from second fit
          finalValMap = interimValMap;
          finalErrMap = interimErrMap;
          chi2_total = chi2_total2;
          ndf_total  = ndf_total2;
      } else {
          // without second scaling
          finalValMap = interimValMap;
          finalErrMap = interimErrMap;
          chi2_total = chi2_total1;
          ndf_total  = ndf_total1;
      }

    } else {
      // ─── Strategy 2: global simultaneous fit ────────────────────────────────
      gStage = 0;
      build_par_listIm();
      build_par_listRe();
      parNamesAll = parNamesIm;
      parNamesAll.insert(parNamesAll.end(),
                         parNamesRe.begin(), parNamesRe.end());
      int nAll = parNamesAll.size();
      std::vector<double> allVal(nAll), allErr(nAll);
      double chi2_glob, edm, errdef;
      int nv,nx,ic;

      auto run_global_fit = [&](){
        TMinuit minu(nAll);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);
        for(int i=0;i<nAll;++i){
          const auto &nm = parNamesAll[i];
          double init=0, lo=-1e6, hi=+1e6, step=0.01;

          // fix P = 1 (keep previous behavior in global stage, too)
          if(nm=="P_H"  || nm=="P_Ht"
          || nm=="P_E"  || nm=="P_Et"){
            init = 1.0; lo = 1.0; hi = 1.0; step = 0.0;
          }

          if     (nm=="r_H")        init=r_H;
          else if(nm=="alpha0_H")   init=alpha0_H;
          else if(nm=="alpha1_H")   init=alpha1_H;
          else if(nm=="b_H")        init=b_H;
          else if(nm=="M2_H")       init=M2_H;
          else if(nm=="P_H")        init=P_H;
          else if(nm=="r_Ht")       init=r_Ht;
          else if(nm=="alpha0_Ht")  init=alpha0_Ht;
          else if(nm=="alpha1_Ht")  init=alpha1_Ht;
          else if(nm=="b_Ht")       init=b_Ht;
          else if(nm=="M2_Ht")      init=M2_Ht;
          else if(nm=="P_Ht")       init=P_Ht;
          else if(nm=="r_E")        init=r_E;
          else if(nm=="alpha0_E")   init=alpha0_E;
          else if(nm=="alpha1_E")   init=alpha1_E;
          else if(nm=="b_E")        init=b_E;
          else if(nm=="M2_E")       init=M2_E;
          else if(nm=="P_E")        init=P_E;
          else if(nm=="r_Et")       init=r_Et;
          else if(nm=="alpha0_Et")  init=alpha0_Et;
          else if(nm=="alpha1_Et")  init=alpha1_Et;
          else if(nm=="b_Et")       init=b_Et;
          else if(nm=="M2_Et")      init=M2_Et;
          else if(nm=="P_Et")       init=P_Et;
          else if(nm=="renormReal") init=renormReal;
          else if(nm=="C0_H")       init=C0_H;
          else if(nm=="MD2_H")      init=MD2_H;
          else if(nm=="lambda_H")   init=lambda_H;
          else if(nm=="C0_Ht")      init=C0_Ht;
          else if(nm=="MD2_Ht")     init=MD2_Ht;
          else if(nm=="lambda_Ht")  init=lambda_Ht;
          else if(nm=="C0_E")       init=C0_E;
          else if(nm=="MD2_E")      init=MD2_E;
          else if(nm=="lambda_E")   init=lambda_E;
          else if(nm=="C0_Et")      init=C0_Et;
          else if(nm=="MD2_Et")     init=MD2_Et;
          else if(nm=="lambda_Et")  init=lambda_Et;

          minu.DefineParameter(i, nm.c_str(), init, step, lo, hi);
          if(step==0.0) minu.FixParameter(i);
        }
        std::cout<<" Stage global: fitting Im + renormReal + subtraction constants…\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_glob,edm,errdef,nv,nx,ic);
        for(int i=0;i<nAll;++i) minu.GetParameter(i,allVal[i],allErr[i]);
        // capture
        interimValMap.clear(); interimErrMap.clear();
        for(int i=0;i<nAll;++i){
          interimValMap[parNamesAll[i]] = allVal[i];
          interimErrMap[parNamesAll[i]] = allErr[i];
        }
        chi2_total1 = chi2_glob;
        ndf_total1  = Nbins + (int)xsData.size() - nAll;
      };

      // first global fit
      run_global_fit();
      std::cout<<"CFF global fit χ²/ndf before scaling amplitude uncertainties = "
               <<chi2_total1<<"/"<<ndf_total1
               <<" = "<<(ndf_total1>0?chi2_total1/ndf_total1:0)<<"\n";

      if(gScale==1){
          double scale_cff = std::sqrt((ndf_total1>0?chi2_total1/ndf_total1:1.0));
          std::cout<<"Scaling amplitude uncertainties (bin_dA) by sqrt("<<chi2_total1<<"/"<<ndf_total1<<") = "<<scale_cff<<"\n";
          for(auto &dA : bin_dA) dA *= scale_cff;
          // recompute reducedAmpChi2 after scaling
          double totChi2=0; int totDof=0;
          for(int k=0;k<Nbins;++k){
              totChi2 += std::pow(bin_A[k]/bin_dA[k],2);
              totDof  += (bin_M[k]-1);
          }
          reducedAmpChi2 = totDof>0 ? totChi2/totDof : 0.0;

          // rerun global fit with updated bin_dA
          run_global_fit();
          chi2_total2 = chi2_total1;
          ndf_total2  = Nbins + (int)xsData.size() - (int)parNamesAll.size();
          std::cout<<"CFF global fit χ²/ndf after scaling amplitude uncertainties = "
                   <<chi2_total2<<"/"<<ndf_total2
                   <<" = "<<(ndf_total2>0?chi2_total2/ndf_total2:0)<<"\n";
          finalValMap = interimValMap;
          finalErrMap = interimErrMap;
          chi2_total = chi2_total2;
          ndf_total  = ndf_total2;
      } else {
          finalValMap = interimValMap;
          finalErrMap = interimErrMap;
          chi2_total = chi2_total1;
          ndf_total  = ndf_total1;
      }
    }

    // ─── Kinematic ranges over USED BINS (those entering the fit) ────────────
    double xb_min =  std::numeric_limits<double>::infinity();
    double xb_max = -std::numeric_limits<double>::infinity();
    double xi_min =  std::numeric_limits<double>::infinity();
    double xi_max = -std::numeric_limits<double>::infinity();
    double mt_min =  std::numeric_limits<double>::infinity(); // -t minimum
    double mt_max = -std::numeric_limits<double>::infinity(); // -t maximum

    for (int k = 0; k < Nbins; ++k) {
        const double xb = bin_xB[k];
        if (xb > 0.0 && xb < 1.0) {
            xb_min = std::min(xb_min, xb);
            xb_max = std::max(xb_max, xb);
            const double xi = xb / (2.0 - xb);
            xi_min = std::min(xi_min, xi);
            xi_max = std::max(xi_max, xi);
        }
        const double mt = -bin_t[k]; // positive
        mt_min = std::min(mt_min, mt);
        mt_max = std::max(mt_max, mt);
    }
    const bool have_xb = std::isfinite(xb_min) && std::isfinite(xb_max);

    // ─── Output results ───────────────────────────────────────────────────────
    std::vector<std::string> outNames =
      (gStrategy==1 ? (std::vector<std::string>(parNamesIm))
                    : parNamesAll);

    system("mkdir -p output/fit_results");
    std::ofstream fout("output/fit_results/fit_results_"+std::string(tb)+".txt");
    fout<<"# fit_CFFs results\n"
        <<"timestamp   "<<tb<<"\n"
        <<"strategy    "<<gStrategy<<"\n"
        <<"constraint  "<<gConstraint<<"\n"
        <<"scale       "<<gScale<<"\n"
        <<"input       "<<gBsaFile<<"\n"
        <<"H "<<hasH<<" Ht "<<hasHt<<" E "<<hasE<<" Et "<<hasEt<<"\n"
        <<"# parameters:";
    for(auto &n: outNames) fout<<" "<<n;
    fout<<"\n# values:\n";
    for(auto &n: outNames) fout<<finalValMap[n]<<" ";
    fout<<"\n# errors:\n";
    for(auto &n: outNames) fout<<finalErrMap[n]<<" ";
    fout<<"\n# chi2 ndf chi2/ndf\n"
        <<chi2_total<<" "<<ndf_total<<" "<<(ndf_total>0?chi2_total/ndf_total:0)<<"\n";

    // NEW: write kinematic ranges for USED BINS
    fout<<"# kinematic ranges over USED BINS (after cuts)\n";
    if (have_xb){
        fout<<"xB_min "<<xb_min<<"  xB_max "<<xb_max<<"\n";
        fout<<"xi_min "<<xi_min<<"  xi_max "<<xi_max<<"\n";
    } else {
        fout<<"xB_min NA  xB_max NA\n";
        fout<<"xi_min NA  xi_max NA\n";
    }
    if (std::isfinite(mt_min) && std::isfinite(mt_max))
        fout<<"-t_min "<<mt_min<<"  -t_max "<<mt_max<<"\n";
    else
        fout<<"-t_min NA  -t_max NA\n";

    fout.close();

    // Console recap
    std::cout<<"\n--- Fit Results ---\n";
    for(auto &n: outNames){
      std::cout<<" "<<n<<" = "<<finalValMap[n]
               <<" ± "<<finalErrMap[n]<<"\n";
    }
    std::cout<<" χ²/ndf = "<<chi2_total<<"/"<<ndf_total
             <<" = "<<(ndf_total>0?chi2_total/ndf_total:0)<<"\n";
    std::cout<<" Average χ²/bin-fit = "<<avgBinChi2<<"\n";
    std::cout<<" χ²_per_amp-fit = "<<reducedAmpChi2<<"\n";

    if (have_xb){
        std::cout<<" xi/xB range over USED input bins:  min(xi) = "<<xi_min
                 <<" (from xB = "<<xb_min<<")"
                 <<", max(xi) = "<<xi_max
                 <<" (from xB = "<<xb_max<<")\n";
    } else {
        std::cout<<" xi/xB range over USED input bins:  (no valid xB to compute ξ)\n";
    }
    if (std::isfinite(mt_min) && std::isfinite(mt_max)){
        std::cout<<" -t range over USED input bins:     min(-t) = "<<mt_min
                 <<", max(-t) = "<<mt_max<<"\n\n";
    } else {
        std::cout<<" -t range over USED input bins:     (no valid t values)\n\n";
    }

    return 0;
}

// ──────────────────────────────────────────────────────────────────────────────
void PlotBinFit(int ibin, const std::string &ts) {
    if (!gPlotBinFits) return;
    auto &pts = keptBins[ibin];
    int n = pts.size(); if(n<=3) return;
    system("mkdir -p output/plots/binned_fits");
    gStyle->SetOptStat(0);

    TGraphErrors gr(n);
    for(int i=0;i<n;++i){
      gr.SetPoint(i, pts[i].phi, pts[i].A);
      gr.SetPointError(i,0,pts[i].sigA);
    }
    gr.SetMarkerStyle(20);

    TF1 f1(Form("f_bin%d",ibin),
      "[0] + [1]*sin(x*TMath::Pi()/180.)"
           "/(1 + [2]*cos(x*TMath::Pi()/180.))",
      0,360);
    f1.SetParameters(0., bin_A[ibin], 0.);
    f1.SetLineColor(kRed);
    f1.SetLineWidth(2);
    gr.Fit(&f1,"RQ0");

    double chi2 = f1.GetChisquare();
    double ndf  = f1.GetNDF();
    double Bfit = f1.GetParameter(2);
    if(gConstraint==1 && (chi2/ndf>=2.0 || std::fabs(Bfit)>=0.9)) return;
    if(gConstraint==2 && (chi2/ndf>=2.0 || std::fabs(Bfit)>=0.9)) return;

    TCanvas c(Form("c_bin%d",ibin),"",600,500);
    TH1F frame(Form("frame%d",ibin),"",360,0,360);
    frame.SetMinimum(-0.6); frame.SetMaximum(0.6);
    frame.GetXaxis()->SetTitle("#phi (deg)");
    frame.GetYaxis()->SetTitle("A_{LU}");
    frame.Draw("AXIS");
    gr.Draw("P same");
    f1.Draw("L same");

    TLegend leg(0.6,0.75,0.9,0.9);
    leg.SetBorderSize(1);
    leg.SetFillStyle(1001);
    leg.SetFillColor(0);
    leg.AddEntry(&gr, "data", "p");
    leg.AddEntry(&f1,
      Form("C = %.3f #pm %.3f",
           f1.GetParameter(0),f1.GetParError(0)), "l");
    leg.AddEntry(&f1,
      Form("A = %.3f #pm %.3f",
           f1.GetParameter(1),f1.GetParError(1)), "l");
    leg.AddEntry(&f1,
      Form("B = %.3f #pm %.3f",
           Bfit, f1.GetParError(2)), "l");
    leg.AddEntry((TObject*)0,
      Form("#chi^{2}/ndf = %.2f", chi2/ndf), "");
    leg.Draw();

    c.SaveAs(Form("output/plots/binned_fits/BinFit_%s_bin%d.pdf",
                  ts.c_str(), ibin));
}