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
 *             [--constraint <0|1>] [--input <BSA_file>] [--plot-fits]
 *
 *   constraint=0 : no cuts
 *   constraint=1 : per-point (-t/Q2>=0.2 or -t>=1) AND per-bin
 *                  (χ²_red>=2 or |B|>=0.9) cuts
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
// GPD–H defaults (in DVCS_xsec.C)
extern double r_H,   n_H,   alpha0_H,   alpha1_H,   b_H,   M2_H,   P_H;
// GPD–Htilde
extern double r_Ht,  n_Ht,  alpha0_Ht,  alpha1_Ht,  b_Ht,  M2_Ht,  P_Ht;
// GPD–E
extern double r_E,   n_E,   alpha0_E,   alpha1_E,   b_E,   M2_E,   P_E;
// GPD–Etilde
extern double r_Et,  n_Et,  alpha0_Et,  alpha1_Et,  b_Et,  M2_Et,  P_Et;

// subtraction‐relation constants (fit in strategy 2)
extern double C0_H,    MD2_H,    lambda_H;
extern double C0_Ht,   MD2_Ht,   lambda_Ht;
extern double C0_E,    MD2_E,    lambda_E;
extern double C0_Et,   MD2_Et,   lambda_Et;

// ----------------------------------------------------------------------------
// Control flags & data containers
static int   gStrategy    = 0;
static int   gStage       = 1;  // 1 = Im‐only, 0 = global simultaneous
static int   gConstraint  = 0;  // 0 or 1
static bool  gPlotBinFits = false;
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static const char* gXsFile  = "imports/rga_pass1_xsec_2018.txt";

struct DataPoint { double phi, Q2, xB, t, Eb, A, sigA; };
static std::vector<DataPoint> bsaData, xsData;
static std::vector<std::vector<DataPoint>> keptBins;

// per-bin averages & fits
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
      {"input",      required_argument,nullptr,'i'},
      {"plot-fits",  no_argument,      nullptr,'p'},
      {nullptr,0,nullptr,0}
    };
    int c;
    while((c=getopt_long(argc,argv,"s:h:t:e:x:C:i:p",opts,nullptr))!=-1){
        switch(c){
          case 's': gStrategy   = std::atoi(optarg); break;
          case 'h': hasH        = std::atoi(optarg); break;
          case 't': hasHt       = std::atoi(optarg); break;
          case 'e': hasE        = std::atoi(optarg); break;
          case 'x': hasEt       = std::atoi(optarg); break;
          case 'C': gConstraint = std::atoi(optarg); break;
          case 'i': gBsaFile    = optarg;            break;
          case 'p': gPlotBinFits= true;              break;
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" --strategy<1|2> -H<0|1> -Ht<0|1>"
                     <<" -E<0|1> -Et<0|1> [--constraint<0|1>]"
                     <<" [--input <BSA_file>] [--plot-fits]\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>2
     ||gConstraint<0||gConstraint>1){
        std::cerr<<"Invalid strategy or constraint\n";
        std::exit(1);
    }
}

void LoadData(){
    auto read=[&](const char* fn, auto &v){
        std::ifstream in(fn);
        if(!in){ std::cerr<<"ERROR: cannot open "<<fn<<"\n"; std::exit(1); }
        std::string line;
        while(std::getline(in,line)){
            if(line.empty()||line[0]=='#') continue;
            std::istringstream iss(line);
            DataPoint d; iss>>d.phi>>d.Q2>>d.xB>>d.t>>d.Eb>>d.A>>d.sigA;
            if(gConstraint==1){
              if((-d.t/d.Q2)>=0.2) continue;
              if((-d.t)    >=1.0) continue;
            }
            v.push_back(d);
        }
    };
    read(gBsaFile.c_str(), bsaData);
    read(gXsFile,          xsData);
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
        int M = pts.size();
        if(M<3) continue;
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
        double Afit    = ftmp.GetParameter(1);
        double dA      = ftmp.GetParError(1);
        double Bfit    = ftmp.GetParameter(2);
        double chi2    = ftmp.GetChisquare();
        double redchi2 = chi2 / (M - 3);
        if(gConstraint==1 && (redchi2>=2.0 || std::fabs(Bfit)>=0.9))
            continue;
        keptBins.push_back(pts);
        bin_M.push_back(M);
        bin_A.push_back(Afit);
        bin_dA.push_back(dA);
        bin_redChi2.push_back(redchi2);
        double sumw=0, Sx=0, Sq=0, St=0, Se=0;
        for(auto &d: pts){
            double w = 1./(d.sigA*d.sigA);
            sumw+=w; Sx+=w*d.xB; Sq+=w*d.Q2; St+=w*d.t; Se+=w*d.Eb;
        }
        bin_xB .push_back(Sx/sumw);
        bin_Q2 .push_back(Sq/sumw);
        bin_t  .push_back(St/sumw);
        bin_Eb .push_back(Se/sumw);
    }
    Nbins = bin_A.size();
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

// ----------------------------------------------------------------------------
// χ² function
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
        // χ² vs. BSA-derived amplitudes
        double chi2=0;
        for(int k=0;k<Nbins;++k){
            BMK_DVCS dvcs(-1,1,0,
                          bin_Eb[k],bin_xB[k],bin_Q2[k],bin_t[k],0.0);
            double modelA = dvcs.s1_I() / dvcs.c0_BH();
            double resid  = (bin_A[k] - modelA)/bin_dA[k];
            chi2 += resid*resid;
        }
        f = chi2;
    }
    else { // global simultaneous
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
        // BSA part
        for(int k=0;k<Nbins;++k){
            BMK_DVCS dvcs(-1,1,0,
                          bin_Eb[k],bin_xB[k],bin_Q2[k],bin_t[k],0.0);
            double modelA = dvcs.s1_I() / dvcs.c0_BH();
            double resid  = (bin_A[k] - modelA)/bin_dA[k];
            chi2 += resid*resid;
        }
        // xsec part
        for(auto &d: xsData){
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
             <<"  input="<<gBsaFile
             <<"  plot-fits="<<(gPlotBinFits?"ON":"OFF")
             <<"  ts="<<tb<<" ===\n\n";

    LoadData();
    BinBsaData();

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
      std::cout<<" Wrote per-bin fits to output/plots/binned_fits/\n\n";
    }

    std::map<std::string,double> valMap, errMap;
    double chi2_total=0;
    int    ndf_total=0;

    if(gStrategy==1){
      // ─── Stage 1: Im-fit only ────────────────────────────────────────────────
      gStage = 1;
      build_par_listIm();
      int nim = parNamesIm.size();
      std::vector<double> imVal(nim), imErr(nim);
      double chi2_im, edm, errdef;
      int nv,nx,ic, ndf_im;
      {
        TMinuit minu(nim);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);
        for(int i=0;i<nim;++i){
          const auto &nm = parNamesIm[i];
          double init=0, lo=-1e3, hi=+1e3, step=0.01;
          // impose positivity
          if(nm=="alpha0_H"  || nm=="alpha1_H"
          || nm=="alpha0_Ht" || nm=="alpha1_Ht"
          || nm=="alpha0_E"  || nm=="alpha1_E"
          || nm=="alpha0_Et" || nm=="alpha1_Et"
          || nm=="b_H"  || nm=="b_Ht"
          || nm=="b_E"  || nm=="b_Et"
          || nm=="M2_H" || nm=="M2_Ht"
          || nm=="M2_E" || nm=="M2_Et"){
            lo = 0.0;
          }
          // fix P parameters
          if(nm=="P_H"  || nm=="P_Ht"
          || nm=="P_E"  || nm=="P_Et"){
            init = 1.0; lo = 1.0; hi = 1.0; step = 0.0;
          }
          // defaults
          if     (nm=="r_H")       init=r_H;
          else if(nm=="r_Ht")      init=r_Ht;
          else if(nm=="r_E")       init=r_E;
          else if(nm=="r_Et")      init=r_Et;
          minu.DefineParameter(i, nm.c_str(), init, step, lo, hi);
          if(step==0.0) minu.FixParameter(i);
        }
        std::cout<<" Stage1: fitting Im-CFF parameters…\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_im,edm,errdef,nv,nx,ic);
        for(int i=0;i<nim;++i) minu.GetParameter(i,imVal[i],imErr[i]);
        ndf_im = Nbins - nim;
      }
      for(int i=0;i<nim;++i){
        valMap[parNamesIm[i]] = imVal[i];
        errMap[parNamesIm[i]] = imErr[i];
      }
      chi2_total = chi2_im;
      ndf_total  = ndf_im;
    }
    else {
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
      {
        TMinuit minu(nAll);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);
        for(int i=0;i<nAll;++i){
          const auto &nm = parNamesAll[i];
          double init=0, lo=-1e3, hi=+1e3, step=0.01;
          // positivity for all shape parameters:
          if(nm.find("alpha0")!=std::string::npos
          || nm.find("alpha1")!=std::string::npos
          || nm.find("b_")!=std::string::npos
          || nm.find("M2_")!=std::string::npos){
            lo = 0.0;
          }
          // fix P's:
          if(nm=="P_H"  || nm=="P_Ht"
          || nm=="P_E"  || nm=="P_Et"){
            init = 1.0; lo = 1.0; hi = 1.0; step = 0.0;
          }
          // initialize
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
          else if(nm=="MD2_H")      init=MD2_H, lo=0.0;
          else if(nm=="lambda_H")   init=lambda_H, lo=0.0;
          else if(nm=="C0_Ht")      init=C0_Ht;
          else if(nm=="MD2_Ht")     init=MD2_Ht, lo=0.0;
          else if(nm=="lambda_Ht")  init=lambda_Ht, lo=0.0;
          else if(nm=="C0_E")       init=C0_E;
          else if(nm=="MD2_E")      init=MD2_E, lo=0.0;
          else if(nm=="lambda_E")   init=lambda_E, lo=0.0;
          else if(nm=="C0_Et")      init=C0_Et;
          else if(nm=="MD2_Et")     init=MD2_Et, lo=0.0;
          else if(nm=="lambda_Et")  init=lambda_Et, lo=0.0;
          minu.DefineParameter(i, nm.c_str(), init, step, lo, hi);
          if(step==0.0) minu.FixParameter(i);
        }
        std::cout<<" Stage global: fitting Im + renormReal + subtraction constants…\n";
        minu.Migrad();
        minu.Command("HESSE");
        minu.mnstat(chi2_glob,edm,errdef,nv,nx,ic);
        for(int i=0;i<nAll;++i) minu.GetParameter(i,allVal[i],allErr[i]);
      }
      for(int i=0;i<(int)parNamesAll.size();++i){
        valMap[parNamesAll[i]] = allVal[i];
        errMap[parNamesAll[i]] = allErr[i];
      }
      chi2_total = chi2_glob;
      // degrees of freedom = N_bsa_bins + N_xs_points - N_pars
      ndf_total = Nbins + xsData.size() - parNamesAll.size();
    }

    // ─── Output results ─────────────────────────────────────────────────────────
    std::vector<std::string> outNames =
      (gStrategy==1 ? parNamesIm : parNamesAll);

    system("mkdir -p output/fit_results");
    std::ofstream fout("output/fit_results/fit_results_"+std::string(tb)+".txt");
    fout<<"# fit_CFFs results\n"
        <<"timestamp   "<<tb<<"\n"
        <<"strategy    "<<gStrategy<<"\n"
        <<"constraint  "<<gConstraint<<"\n"
        <<"input       "<<gBsaFile<<"\n"
        <<"H "<<hasH<<" Ht "<<hasHt<<" E "<<hasE<<" Et "<<hasEt<<"\n"
        <<"# parameters:";
    for(auto &n: outNames) fout<<" "<<n;
    fout<<"\n# values:\n";
    for(auto &n: outNames) fout<<valMap[n]<<" ";
    fout<<"\n# errors:\n";
    for(auto &n: outNames) fout<<errMap[n]<<" ";
    fout<<"\n# chi2 ndf chi2/ndf\n"
        <<chi2_total<<" "<<ndf_total<<" "<<(chi2_total/ndf_total)<<"\n";
    fout.close();

    std::cout<<"\n--- Fit Results ---\n";
    for(auto &n: outNames){
      std::cout<<" "<<n<<" = "<<valMap[n]
               <<" ± "<<errMap[n]<<"\n";
    }
    std::cout<<" χ²/ndf = "<<chi2_total<<"/"<<ndf_total
             <<" = "<<(chi2_total/ndf_total)<<"\n";
    std::cout<<" Average χ²/bin-fit = "<<avgBinChi2<<"\n";
    std::cout<<" χ²_per_amp-fit = "<<reducedAmpChi2<<"\n\n";
    return 0;
}

// ──────────────────────────────────────────────────────────────────────────────
void PlotBinFit(int ibin, const std::string &ts) {
    if (!gPlotBinFits) return;
    auto &pts = keptBins[ibin];
    int n = pts.size(); if(n<3) return;
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