/**
 * fit_CFFs.C
 * ──────────
 * program to fit DVCS Compton Form Factors (CFFs) imaginary then real parts
 *
 * Strategies:
 *   1) fit only Im–parts → BSA data
 *   2) two-step: (a) fit Im–parts → BSA, (b) fit renormReal → xsec
 *
 * Usage:
 *   ./fit_CFFs --strategy <1|2> -H <0|1> -Ht <0|1> -E <0|1> -Et <0|1>
 *             [--constraint <0|1>] [--input <BSA_file>] [--plot-fits]
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

// fixed integer–spin prefactors
extern double n_H, n_Ht, n_E, n_Et;

// **new** imaginary-part model parameters (six floats per CFF)
extern double r_H,      alpha0_H,  alpha1_H,  beta0_H,  beta1_H,  M2_H,  P_H;
extern double r_Ht,     alpha0_Ht, alpha1_Ht, beta0_Ht, beta1_Ht, M2_Ht, P_Ht;
extern double r_E,      alpha0_E,  alpha1_E,  beta0_E,  beta1_E,  M2_E,  P_E;
extern double r_Et,     alpha0_Et, alpha1_Et, beta0_Et, beta1_Et, M2_Et, P_Et;

// ──────────────────────────────────────────────────────────────────────────────
// control flags
static int   gStrategy     = 0;     // 1 or 2
static int   gStage        = 1;     // 1 = Im‐fit, 2 = Re‐fit
static int   gConstraint   = 0;     // 0 = no cut, 1 = apply -t/Q2<0.2
static bool  gPlotBinFits  = false; // toggle per‐bin φ‐fits plotting
static std::string gBsaFile = "imports/rga_prl_bsa.txt";
static const char* gXsFile  = "imports/rga_pass1_xsec_2018.txt";

// ──────────────────────────────────────────────────────────────────────────────
// raw data + binned observables
struct DataPoint { double phi, Q2, xB, t, Eb, A, sigA; };
static std::vector<DataPoint> bsaData, xsData;
static std::vector<std::vector<DataPoint>> binnedPoints;
static std::vector<double> bin_xB, bin_Q2, bin_t, bin_Eb;
static std::vector<double> bin_A, bin_dA, bin_chi2;
static std::vector<int>    bin_M;      // points per bin
static int     Nbins          = 0;
static double  reducedAmpChi2 = 0.0;

// ──────────────────────────────────────────────────────────────────────────────
// forward‐declare helper
void PlotBinFit(int ibin, const std::string &ts);

// ──────────────────────────────────────────────────────────────────────────────
// parse_args(): --strategy, -H, -Ht, -E, -Et, [--constraint], [--input], [--plot-fits]
void parse_args(int argc, char** argv){
    static struct option opts[] = {
        {"strategy",   required_argument, nullptr, 's'},
        {"H",          required_argument, nullptr, 'h'},
        {"Ht",         required_argument, nullptr, 't'},
        {"E",          required_argument, nullptr, 'e'},
        {"Et",         required_argument, nullptr, 'x'},
        {"constraint", required_argument, nullptr, 'C'},
        {"input",      required_argument, nullptr, 'i'},
        {"plot-fits",  no_argument,       nullptr, 'p'},
        {nullptr,0,nullptr,0}
    };
    int c;
    while((c = getopt_long(argc, argv, "s:h:t:e:x:C:i:p", opts, nullptr)) != -1){
        switch(c){
          case 's': gStrategy   = std::atoi(optarg); break;
          case 'h': hasH        = std::atoi(optarg); break;
          case 't': hasHt       = std::atoi(optarg); break;
          case 'e': hasE        = std::atoi(optarg); break;
          case 'x': hasEt       = std::atoi(optarg); break;
          case 'C': gConstraint = std::atoi(optarg); break;
          case 'i': gBsaFile    = std::string(optarg); break;
          case 'p': gPlotBinFits= true;                break;
          default:
            std::cerr<<"Usage: "<<argv[0]
                     <<" --strategy<1|2> -H<0|1> -Ht<0|1>"
                     <<" -E<0|1> -Et<0|1> [--constraint<0|1>]"
                     <<" [--input <BSA_file>] [--plot-fits]\n";
            std::exit(1);
        }
    }
    if(gStrategy<1||gStrategy>2){
        std::cerr<<"Invalid strategy\n"; std::exit(1);
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Load raw BSA + XSC, apply constraint if requested
void LoadData(){
    auto read = [&](const char* fn, auto &v){
        std::ifstream in(fn);
        if(!in){ std::cerr<<"ERROR: cannot open "<<fn<<"\n"; std::exit(1); }
        std::string line;
        while(std::getline(in,line)){
            if(line.empty()||line[0]=='#') continue;
            std::istringstream iss(line);
            DataPoint d;
            iss>>d.phi>>d.Q2>>d.xB>>d.t>>d.Eb>>d.A>>d.sigA;
            if(gConstraint==1 && (-d.t/d.Q2) >= 0.2) continue;
            v.push_back(d);
        }
    };
    read(gBsaFile.c_str(), bsaData);
    read(gXsFile,          xsData);
}

// ──────────────────────────────────────────────────────────────────────────────
// Bin BSA by φ-drop, fit A·sinφ/(1+B·cosφ) → extract A,dA,χ² and weighted kin.
void BinBsaData(){
    bin_xB.clear(); bin_Q2.clear(); bin_t.clear(); bin_Eb.clear();
    bin_A.clear();  bin_dA.clear();  bin_chi2.clear(); bin_M.clear();
    binnedPoints.clear();
    if(bsaData.empty()) return;

    size_t start = 0;
    for(size_t i=1; i<=bsaData.size(); ++i){
        bool newbin = (i==bsaData.size() || bsaData[i].phi < bsaData[i-1].phi);
        if(newbin){
            auto pts = std::vector<DataPoint>(
                bsaData.begin()+start, bsaData.begin()+i
            );
            binnedPoints.push_back(pts);
            int M = pts.size();

            // single‐parameter sinφ fit → A_bin ± dA_bin
            double SwA=0, Sw2=0;
            for(auto &d: pts){
                double s = std::sin(d.phi * TMath::Pi()/180.);
                double w = 1.0/(d.sigA*d.sigA);
                SwA += w*d.A*s;
                Sw2 += w*s*s;
            }
            double A_bin  = SwA / Sw2;
            double dA_bin = 1.0/std::sqrt(Sw2);

            // χ²
            double chi2=0;
            for(auto &d: pts){
                double s    = std::sin(d.phi * TMath::Pi()/180.);
                double diff = d.A - A_bin*s;
                chi2 += diff*diff/(d.sigA*d.sigA);
            }

            // weighted kinematics
            double sumw=0, Sx=0, Sq=0, St=0, Se=0;
            for(auto &d: pts){
                double w = 1.0/(d.sigA*d.sigA);
                sumw += w;
                Sx   += w*d.xB;
                Sq   += w*d.Q2;
                St   += w*d.t;
                Se   += w*d.Eb;
            }
            bin_xB .push_back(Sx/sumw);
            bin_Q2 .push_back(Sq/sumw);
            bin_t  .push_back(St/sumw);
            bin_Eb .push_back(Se/sumw);

            bin_A    .push_back(A_bin);
            bin_dA   .push_back(dA_bin);
            bin_chi2 .push_back(chi2);
            bin_M    .push_back(M);

            start = i;
        }
    }
    Nbins = bin_A.size();

    double totalChi2 = std::accumulate(bin_chi2.begin(), bin_chi2.end(), 0.0);
    int totalDof     = std::accumulate(bin_M.begin(), bin_M.end(), 0) - Nbins;
    reducedAmpChi2   = totalDof>0 ? totalChi2/totalDof : 0.0;
}

// ──────────────────────────────────────────────────────────────────────────────
// build which Im‐parameters to fit: **no** renormImag, **no** n_ fixed
static std::vector<std::string> parNamesIm;
void build_par_list(){
    parNamesIm.clear();
    if(hasH )  parNamesIm.insert(parNamesIm.end(),
       {"r_H","alpha0_H","alpha1_H","beta0_H","beta1_H","M2_H","P_H"});
    if(hasHt)  parNamesIm.insert(parNamesIm.end(),
       {"r_Ht","alpha0_Ht","alpha1_Ht","beta0_Ht","beta1_Ht","M2_Ht","P_Ht"});
    if(hasE )  parNamesIm.insert(parNamesIm.end(),
       {"r_E","alpha0_E","alpha1_E","beta0_E","beta1_E","M2_E","P_E"});
    if(hasEt)  parNamesIm.insert(parNamesIm.end(),
       {"r_Et","alpha0_Et","alpha1_Et","beta0_Et","beta1_Et","M2_Et","P_Et"});
}

// ──────────────────────────────────────────────────────────────────────────────
// χ² function: Im‐fit (gStage=1) or renormReal‐fit (gStage=2)
void fcn(int&, double*, double &f, double *par, int){
    int ip = 0;
    if(gStage==1){
        // load each floated parameter in the same order as parNamesIm
        if(hasH ){
          r_H      = par[ip++];
          alpha0_H = par[ip++];
          alpha1_H = par[ip++];
          beta0_H  = par[ip++];
          beta1_H  = par[ip++];
          M2_H     = par[ip++];
          P_H      = par[ip++];
        }
        if(hasHt){
          r_Ht      = par[ip++];
          alpha0_Ht = par[ip++];
          alpha1_Ht = par[ip++];
          beta0_Ht  = par[ip++];
          beta1_Ht  = par[ip++];
          M2_Ht     = par[ip++];
          P_Ht      = par[ip++];
        }
        if(hasE ){
          r_E      = par[ip++];
          alpha0_E = par[ip++];
          alpha1_E = par[ip++];
          beta0_E  = par[ip++];
          beta1_E  = par[ip++];
          M2_E     = par[ip++];
          P_E      = par[ip++];
        }
        if(hasEt){
          r_Et      = par[ip++];
          alpha0_Et = par[ip++];
          alpha1_Et = par[ip++];
          beta0_Et  = par[ip++];
          beta1_Et  = par[ip++];
          M2_Et     = par[ip++];
          P_Et      = par[ip++];
        }

        // compute χ²
        double chi2 = 0;
        for(int k=0; k<Nbins; ++k){
            BMK_DVCS dvcs(-1,1,0,
                          bin_Eb[k], bin_xB[k], bin_Q2[k], bin_t[k], 0.0);
            double modelA = dvcs.s1_I()/dvcs.c0_BH();
            double r      = (bin_A[k] - modelA)/bin_dA[k];
            chi2 += r*r;
        }
        f = chi2;
    }
    else {
        // renormReal‐fit unchanged
        renormReal = par[ip++];
        double chi2 = 0;
        for(auto &d: xsData){
            BMK_DVCS dvcs(-1,0,0,d.Eb,d.xB,d.Q2,d.t,d.phi);
            double m = dvcs.CrossSection();
            double r = (d.A - renormReal*m)/d.sigA;
            chi2 += r*r;
        }
        f = chi2;
    }
}

// ──────────────────────────────────────────────────────────────────────────────
int main(int argc, char** argv){
    parse_args(argc,argv);

    // one timestamp for both text + plots
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
    std::cout<<" BSA bins="<<Nbins
             <<"  (raw="<<bsaData.size()<<")\n";
    std::cout<<" Reduced χ² per amplitude‐fit = "<<reducedAmpChi2<<"\n\n";

    if(gPlotBinFits){
        gStyle->SetOptStat(0);
        for(int ib=0; ib<Nbins; ++ib) PlotBinFit(ib, tb);
        std::cout<<"Wrote φ‐fit plots to output/plots/binned_fits/\n\n";
    }

    // ─── Stage 1: Im‐fit ────────────────────────────────────────────────────────
    gStage = 1;
    build_par_list();
    int nim = parNamesIm.size();
    std::vector<double> imVal(nim), imErr(nim);
    double chi2_im, edm, errdef;
    int nv,nx,ic, ndf_im;

    {
        TMinuit minu(nim);
        minu.SetPrintLevel(1);
        minu.SetFCN(fcn);
        // define all parameters
        for(int i=0; i<nim; ++i){
            const auto &nm = parNamesIm[i];
            double init=0.0, step=0.01;
            #define GETINIT(NAME) if(nm==#NAME) init = NAME;
            GETINIT(r_H)       GETINIT(alpha0_H)
            GETINIT(alpha1_H)  GETINIT(beta0_H)
            GETINIT(beta1_H)   GETINIT(M2_H)
            GETINIT(P_H)
            GETINIT(r_Ht)      GETINIT(alpha0_Ht)
            GETINIT(alpha1_Ht) GETINIT(beta0_Ht)
            GETINIT(beta1_Ht)  GETINIT(M2_Ht)
            GETINIT(P_Ht)
            GETINIT(r_E)       GETINIT(alpha0_E)
            GETINIT(alpha1_E)  GETINIT(beta0_E)
            GETINIT(beta1_E)   GETINIT(M2_E)
            GETINIT(P_E)
            GETINIT(r_Et)      GETINIT(alpha0_Et)
            GETINIT(alpha1_Et) GETINIT(beta0_Et)
            GETINIT(beta1_Et)  GETINIT(M2_Et)
            GETINIT(P_Et)
            #undef GETINIT

            double lo = (nm.rfind("r_",0)==0 ? 0.0 : -1e3);
            minu.DefineParameter(i, nm.c_str(), init, step, lo, 1e3);
        }
        std::cout<<"Stage1: fitting Im‐CFF parameters...\n";
        minu.Migrad(); minu.Command("HESSE");
        minu.mnstat(chi2_im,edm,errdef,nv,nx,ic);
        for(int i=0; i<nim; ++i) minu.GetParameter(i, imVal[i], imErr[i]);
        ndf_im = Nbins - nim;
    }

    // collect results
    std::map<std::string,double> valMap, errMap;
    for(int i=0;i<nim;++i){
        valMap[parNamesIm[i]] = imVal[i];
        errMap[parNamesIm[i]] = imErr[i];
    }

    // ─── Stage 2: renormReal‐fit (if requested) ────────────────────────────────
    if(gStrategy==2){
        gStage = 2;
        double chi2_re, edm2, errdef2; int nv2,nx2,ic2, ndf_re;
        double reVal, reErr;
        {
            TMinuit m2(1);
            m2.SetPrintLevel(1);
            m2.SetFCN(fcn);
            m2.DefineParameter(0,"renormReal",renormReal,0.01,-1e3,1e3);
            std::cout<<"Stage2: fitting renormReal...\n";
            m2.Migrad(); m2.Command("HESSE");
            m2.mnstat(chi2_re,edm2,errdef2,nv2,nx2,ic2);
            m2.GetParameter(0, reVal, reErr);
            ndf_re = xsData.size() - 1;
        }
        valMap["renormReal"] = reVal;
        errMap["renormReal"] = reErr;
        chi2_im = chi2_re;
        ndf_im  = ndf_re;
    }

    // write results file
    std::vector<std::string> outNames = parNamesIm;
    if(gStrategy==2) outNames.push_back("renormReal");

    system("mkdir -p output/fit_results");
    std::ofstream fout("output/fit_results/fit_results_"+std::string(tb)+".txt");
    fout<<"# fit_CFFs results\n"
        <<"timestamp   "<<tb<<"\n"
        <<"strategy    "<<gStrategy<<"\n"
        <<"constraint  "<<gConstraint<<"\n"
        <<"input       "<<gBsaFile<<"\n"
        <<"H "<<hasH<<"  Ht "<<hasHt
        <<"  E "<<hasE<<"  Et "<<hasEt<<"\n"
        <<"# parameters:";
    for(auto &n: outNames) fout<<" "<<n;
    fout<<"\n# values:\n";
    for(auto &n: outNames) fout<<valMap[n]<<" ";
    fout<<"\n# errors:\n";
    for(auto &n: outNames) fout<<errMap[n]<<" ";
    fout<<"\n# chi2 ndf chi2/ndf\n"
        <<chi2_im<<" "<<ndf_im<<" "<<(chi2_im/ndf_im)<<"\n";
    fout.close();

    // echo to stdout
    std::cout<<"\n--- Fit Results ---\n";
    for(auto &n: outNames){
        std::cout<<" "<<n<<" = "
                 <<valMap[n]<<" ± "<<errMap[n]<<"\n";
    }
    std::cout<<" χ²/ndf = "<<chi2_im<<"/"<<ndf_im
             <<" = "<<(chi2_im/ndf_im)<<"\n";
    std::cout<<" Reduced χ² per amplitude‐fit = "<<reducedAmpChi2<<"\n";

    return 0;
}


// ──────────────────────────────────────────────────────────────────────────────
// ──────────────────────────────────────────────────────────────────────────────
// helper: plot the phi-distribution and fit it to C + A sinφ/(1+B cosφ)
void PlotBinFit(int ibin, const std::string &ts) {
    if (!gPlotBinFits) return;

    // regroup raw points into φ-bins
    std::vector<std::vector<DataPoint>> bins;
    size_t start = 0;
    for (size_t i = 1; i <= bsaData.size(); ++i) {
        bool newbin = (i == bsaData.size() || bsaData[i].phi < bsaData[i-1].phi);
        if (newbin) {
            bins.emplace_back(bsaData.begin() + start, bsaData.begin() + i);
            start = i;
        }
    }
    if (ibin < 0 || ibin >= int(bins.size())) return;
    auto &dp = bins[ibin];
    int n = dp.size();
    if (n < 2) return;

    // ensure directory
    system("mkdir -p output/plots/binned_fits");

    // disable default stat box & grid
    gStyle->SetOptStat(0);

    // graph
    TGraphErrors *gr = new TGraphErrors(n);
    for (int i = 0; i < n; ++i) {
        gr->SetPoint(i, dp[i].phi, dp[i].A);
        gr->SetPointError(i, 0.0, dp[i].sigA);
    }
    gr->SetMarkerStyle(20);

    // fit with offset C + A sinφ/(1+B cosφ)
    TF1 *f1 = new TF1(Form("f_bin%d", ibin),
        "[0] + [1]*sin(x*TMath::Pi()/180)/(1+[2]*cos(x*TMath::Pi()/180))",
        0, 360);
    f1->SetParameter(0, 0.0);               // C initial
    f1->SetParameter(1, bin_A[ibin]);       // A initial
    f1->SetParameter(2, 0.0);               // B initial
    // **constrain A and B**
    f1->SetParLimits(1, -1.0, 1.0);
    f1->SetParLimits(2, -1.0, 1.0);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    gr->Fit(f1, "RQN");

    // get chi2 and ndf
    double chi2 = f1->GetChisquare();
    double ndf  = f1->GetNDF();

    // canvas + frame
    TCanvas *c = new TCanvas(Form("c_bin%d", ibin), "", 600, 500);
    TH1F *frame = new TH1F(Form("frame%d", ibin), "", 360, 0, 360);
    frame->SetMinimum(-0.6);
    frame->SetMaximum(0.6);
    frame->GetXaxis()->SetTitle("#phi (deg)");
    frame->GetYaxis()->SetTitle("A_{LU}");
    frame->Draw("AXIS");
    frame->GetXaxis()->SetLimits(0, 360);
    frame->GetXaxis()->SetRangeUser(0, 360);
    gPad->Modified(); gPad->Update();

    // draw points + fit
    gr->Draw("P same");
    f1->Draw("L same");

    // legend
    TLegend *leg = new TLegend(0.60, 0.75, 0.90, 0.90);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(0);
    leg->AddEntry(gr, "data", "p");
    leg->AddEntry(f1,
        Form("C = %.3f +/- %.3f", f1->GetParameter(0), f1->GetParError(0)),
        "l");
    leg->AddEntry(f1,
        Form("A = %.3f +/- %.3f", f1->GetParameter(1), f1->GetParError(1)),
        "l");
    leg->AddEntry(f1,
        Form("B = %.3f +/- %.3f", f1->GetParameter(2), f1->GetParError(2)),
        "l");
    // **add chi2/ndf text**
    leg->AddEntry((TObject*)0,
        Form("#chi^{2}/ndf = %.2f", chi2/ndf),
        "");
    leg->Draw();

    // save with timestamp format
    c->SaveAs(Form("output/plots/binned_fits/BinFit_%s_bin%d.pdf", ts.c_str(), ibin));

    // cleanup
    delete leg;
    delete frame;
    delete c;
    delete gr;
    delete f1;
} 