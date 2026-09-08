#include "bin_centering_systematic_plots.h"
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
namespace fs=std::filesystem;
namespace {
struct Table{std::vector<std::string> h;std::vector<std::vector<std::string>> r;std::unordered_map<std::string,int> i;};
std::vector<std::string> split(const std::string& l){std::vector<std::string> o;std::string c;bool q=false;for(size_t k=0;k<l.size();++k){char x=l[k];if(x=='"'){if(q&&k+1<l.size()&&l[k+1]=='"'){c+='"';++k;}else q=!q;}else if(x==','&&!q){o.push_back(c);c.clear();}else c+=x;}o.push_back(c);return o;}
Table read(const std::string&p){std::ifstream f(p);if(!f)throw std::runtime_error("Could not open "+p);Table t;std::string l;if(!std::getline(f,l))throw std::runtime_error("Empty CSV: "+p);t.h=split(l);if(!t.h.empty()&&(t.h[0].empty()||t.h[0].find("Unnamed")!=std::string::npos))t.h[0]="bin index";for(int k=0;k<(int)t.h.size();++k)t.i[t.h[k]]=k;while(std::getline(f,l)){if(l.empty())continue;auto v=split(l);v.resize(t.h.size());t.r.push_back(std::move(v));}return t;}
int col(const Table&t,const std::string&n){auto it=t.i.find(n);if(it==t.i.end())throw std::runtime_error("Missing column: "+n);return it->second;}
double num(const std::string&s){if(s.empty())return NAN;char*e=nullptr;double v=strtod(s.c_str(),&e);return e==s.c_str()?NAN:v;}
double qtile(std::vector<double>v,double q){if(v.empty())return NAN;std::sort(v.begin(),v.end());double x=q*(v.size()-1);size_t a=(size_t)floor(x),b=std::min(a+1,v.size()-1);return v[a]+(x-a)*(v[b]-v[a]);}
double mean(const std::vector<double>&v){return v.empty()?NAN:std::accumulate(v.begin(),v.end(),0.0)/v.size();}
void style(TCanvas&c){
 c.SetLeftMargin(.15);
 c.SetRightMargin(.035);
 c.SetBottomMargin(.135);
 c.SetTopMargin(.165);
 c.SetTicks(1,1);
}
void axes(TAxis*x,TAxis*y){
 x->SetTitleSize(.040);y->SetTitleSize(.040);
 x->SetLabelSize(.035);y->SetLabelSize(.035);
 x->SetTitleOffset(1.08);y->SetTitleOffset(1.50);
}
}
bool make_bin_centering_systematic_analysis_note_plots(const std::string& sys_csv,const std::string& corr_csv,const std::string& out){
 try{
  fs::create_directories(out);Table s=read(sys_csv),c=read(corr_csv);
  int sval=col(s,"val"), sf=col(s,"Syst.err (Fbin)"), sx0=col(s,"xBmin"),sx1=col(s,"xBmax"),sq0=col(s,"Q2min"),sq1=col(s,"Q2max"),st0=col(s,"tmin"),st1=col(s,"tmax"),sp0=col(s,"phimin"),sp1=col(s,"phimax");
  std::vector<double> rel;for(auto&r:s.r){double v=num(r[sval]),e=num(r[sf]);if(std::isfinite(v)&&fabs(v)>0&&std::isfinite(e)&&e>=0)rel.push_back(100*e/fabs(v));}if(rel.empty())throw std::runtime_error("No valid Fbin systematic values");
  {std::ofstream o(fs::path(out)/"bin_centering_systematic_summary.csv");o<<"n,mean_percent,median_percent,p16_percent,p84_percent,p95_percent\n"<<rel.size()<<","<<std::setprecision(10)<<mean(rel)<<","<<qtile(rel,.5)<<","<<qtile(rel,.16)<<","<<qtile(rel,.84)<<","<<qtile(rel,.95)<<"\n";}
  {TCanvas cv("c_fbin_sys_dist","",1050,740);style(cv);TH1D h("h_fbin_sys_dist","",45,0,std::max(3.0,std::ceil(qtile(rel,.99)+.5)));for(double x:rel)h.Fill(x);if(h.Integral()>0)h.Scale(1/h.Integral());h.SetLineColor(kMagenta+2);h.SetLineWidth(3);h.SetFillColorAlpha(kMagenta+1,.18);h.GetXaxis()->SetTitle("Bin-centering systematic uncertainty (%)");h.GetYaxis()->SetTitle("Fraction of analysis bins");axes(h.GetXaxis(),h.GetYaxis());h.SetMaximum(1.2*h.GetMaximum());h.Draw("HIST");TLine m(qtile(rel,.5),0,qtile(rel,.5),h.GetMaximum());m.SetLineStyle(7);m.SetLineWidth(2);m.SetLineColor(kMagenta+3);m.Draw();TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.032);t.DrawLatex(.15,.952,"Bin-centering model uncertainty");TLegend l(.66,.78,.93,.88);l.SetBorderSize(0);l.SetFillStyle(0);std::ostringstream ss;ss<<"Median = "<<std::fixed<<std::setprecision(2)<<qtile(rel,.5)<<"%";l.AddEntry(&m,ss.str().c_str(),"l");l.Draw();cv.SaveAs((fs::path(out)/"bin_centering_systematic_distribution.png").string().c_str());}
  {std::vector<std::pair<double,double>> edges;for(auto&r:s.r){double a=num(r[sx0]),b=num(r[sx1]);if(std::isfinite(a)&&std::isfinite(b)&&std::find(edges.begin(),edges.end(),std::make_pair(a,b))==edges.end())edges.push_back({a,b});}std::sort(edges.begin(),edges.end());TCanvas cv("c_fbin_sys_xb","",1050,740);style(cv);cv.SetGridy();TH1F fr("h_fbin_sys_xb","",100,.05,.60);fr.SetMinimum(0);fr.SetMaximum(std::max(2.5,1.15*qtile(rel,.95)));fr.GetXaxis()->SetTitle("x_{B}");fr.GetYaxis()->SetTitle("Bin-centering systematic uncertainty (%)");axes(fr.GetXaxis(),fr.GetYaxis());fr.Draw();TGraphAsymmErrors g;g.SetMarkerStyle(20);g.SetMarkerColor(kMagenta+2);g.SetLineColor(kMagenta+2);int n=0;for(auto&e:edges){std::vector<double>v;for(auto&r:s.r)if(fabs(num(r[sx0])-e.first)<1e-8&&fabs(num(r[sx1])-e.second)<1e-8){double xs=num(r[sval]),u=num(r[sf]);if(std::isfinite(xs)&&fabs(xs)>0&&std::isfinite(u))v.push_back(100*u/fabs(xs));}if(v.empty())continue;double md=qtile(v,.5),lo=qtile(v,.16),hi=qtile(v,.84);g.SetPoint(n,.5*(e.first+e.second),md);g.SetPointError(n,0,0,md-lo,hi-md);++n;}g.Draw("PE SAME");TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.032);t.DrawLatex(.15,.952,"Kinematic dependence of the bin-centering model uncertainty");TLatex nte;nte.SetNDC();nte.SetTextFont(42);nte.SetTextSize(.023);nte.DrawLatex(.15,.875,"Points show the median; bars show the central 68% interval within each x_{B} bin");cv.SaveAs((fs::path(out)/"bin_centering_systematic_vs_xB.png").string().c_str());}
  // 2x2 kinematic summary; the xB-only plot is retained for the current note.
  {
   struct AxisSpec{int lo,hi;const char* title;};
   const std::array<AxisSpec,4> specs={{{sx0,sx1,"x_{B}"},{sq0,sq1,"Q^{2} (GeV^{2})"},{st0,st1,"|t| (GeV^{2})"},{sp0,sp1,"#phi (deg)"}}};
   TCanvas cv("c_fbin_sys_kinematic","",1450,1050);cv.Divide(2,2,.002,.002);
   for(int ia=0;ia<4;++ia){
    cv.cd(ia+1);gPad->SetLeftMargin((ia%2==0)?.14:.12);gPad->SetRightMargin(.035);gPad->SetBottomMargin((ia>=2)?.15:.12);gPad->SetTopMargin(.10);gPad->SetTicks(1,1);
    std::vector<std::pair<double,double>> edges;for(auto&r:s.r){double a=num(r[specs[ia].lo]),b=num(r[specs[ia].hi]);if(std::isfinite(a)&&std::isfinite(b)&&b>a&&std::find(edges.begin(),edges.end(),std::make_pair(a,b))==edges.end())edges.push_back({a,b});}std::sort(edges.begin(),edges.end());
    TGraphAsymmErrors g;g.SetMarkerStyle(20);g.SetMarkerSize(.95);g.SetMarkerColor(kMagenta+2);g.SetLineColor(kMagenta+2);g.SetLineWidth(2);double xmin=INFINITY,xmax=-INFINITY,ymax=0;int n=0;
    for(auto&e:edges){std::vector<double>v;for(auto&r:s.r)if(fabs(num(r[specs[ia].lo])-e.first)<1e-8&&fabs(num(r[specs[ia].hi])-e.second)<1e-8){double xs=num(r[sval]),u=num(r[sf]);if(std::isfinite(xs)&&fabs(xs)>0&&std::isfinite(u))v.push_back(100*u/fabs(xs));}if(v.empty())continue;double md=qtile(v,.5),lo=qtile(v,.16),hi=qtile(v,.84),x=.5*(e.first+e.second),ex=.5*(e.second-e.first);g.SetPoint(n,x,md);g.SetPointError(n,ex,ex,md-lo,hi-md);++n;xmin=std::min(xmin,e.first);xmax=std::max(xmax,e.second);ymax=std::max(ymax,hi);}
    if(!(xmax>xmin)){xmin=0;xmax=1;}TH1F fr(("h_fbin_sys_kinematic_"+std::to_string(ia)).c_str(),"",100,xmin,xmax);fr.SetMinimum(0);fr.SetMaximum(std::max(2.5,1.18*ymax));fr.GetXaxis()->SetTitle(specs[ia].title);fr.GetYaxis()->SetTitle("Bin-centering systematic (%)");fr.GetXaxis()->SetTitleSize(.050);fr.GetYaxis()->SetTitleSize(.048);fr.GetXaxis()->SetLabelSize(.043);fr.GetYaxis()->SetLabelSize(.042);fr.GetXaxis()->SetTitleOffset(1.08);fr.GetYaxis()->SetTitleOffset(1.16);fr.Draw();g.Draw("PZ SAME");TLatex p;p.SetNDC();p.SetTextFont(42);p.SetTextSize(.042);std::string lab=std::string("(")+char('a'+ia)+")";p.DrawLatex(.16,.92,lab.c_str());
   }
   cv.cd(0);TLatex tt;tt.SetNDC();tt.SetTextFont(42);tt.SetTextAlign(22);tt.SetTextSize(.024);tt.DrawLatex(.50,.992,"Kinematic dependence of the bin-centering model uncertainty");tt.SetTextSize(.018);tt.DrawLatex(.50,.965,"Points: median in each kinematic interval; bars: central 68% bin-to-bin range");cv.SaveAs((fs::path(out)/"bin_centering_systematic_kinematic_summary.png").string().c_str());
  }
  // Representative phi dependence: central Fbin from correction table with uncertainty inferred from the model spread.
  int cb=col(c,"Fbin"),cx0=col(c,"xBmin"),cx1=col(c,"xBmax"),cq0=col(c,"Q2min"),cq1=col(c,"Q2max"),ct0=col(c,"t_abs_min"),ct1=col(c,"t_abs_max"),cp=col(c,"phiavg");
  const double X0=.204,X1=.268,Q0=1.912,Q1=2.510,T0=.250,T1=.400;TGraphAsymmErrors g;g.SetMarkerStyle(20);g.SetMarkerColor(kMagenta+2);g.SetLineColor(kMagenta+2);g.SetLineWidth(2);int n=0;double ymin=2,ymax=0;
  for(auto&rr:c.r){if(fabs(num(rr[cx0])-X0)>1e-8||fabs(num(rr[cx1])-X1)>1e-8||fabs(num(rr[cq0])-Q0)>1e-8||fabs(num(rr[cq1])-Q1)>1e-8||fabs(num(rr[ct0])-T0)>1e-8||fabs(num(rr[ct1])-T1)>1e-8)continue;double ph=num(rr[cp]),fb=num(rr[cb]);if(!std::isfinite(ph)||!std::isfinite(fb))continue;double relp=NAN;for(auto&rs:s.r){if(fabs(num(rs[sx0])-X0)<1e-8&&fabs(num(rs[sx1])-X1)<1e-8&&fabs(num(rs[sq0])-Q0)<1e-8&&fabs(num(rs[sq1])-Q1)<1e-8&&fabs(num(rs[st0])-T0)<1e-8&&fabs(num(rs[st1])-T1)<1e-8&&ph>=num(rs[sp0])&&ph<=num(rs[sp1])){double xs=num(rs[sval]),u=num(rs[sf]);if(std::isfinite(xs)&&fabs(xs)>0&&std::isfinite(u))relp=u/fabs(xs);break;}}double ey=std::isfinite(relp)?fb*relp:0;g.SetPoint(n,ph,fb);g.SetPointError(n,0,0,ey,ey);ymin=std::min(ymin,fb-ey);ymax=std::max(ymax,fb+ey);++n;}
  if(n>0){TCanvas cv("c_fbin_sys_phi","",1050,740);style(cv);cv.SetGridy();TH1F fr("h_fbin_sys_phi","",100,0,360);fr.SetMinimum(std::min(.94,ymin-.015));fr.SetMaximum(std::max(1.04,ymax+.015));fr.GetXaxis()->SetTitle("#phi (deg)");fr.GetYaxis()->SetTitle("Bin-centering factor, F_{bin}");axes(fr.GetXaxis(),fr.GetYaxis());fr.Draw();g.Draw("PE SAME");TLine one(0,1,360,1);one.SetLineStyle(2);one.SetLineColor(kGray+2);one.Draw();TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.032);t.DrawLatex(.15,.952,"Representative bin-centering correction and model uncertainty");TLatex a;a.SetNDC();a.SetTextFont(42);a.SetTextSize(.022);a.DrawLatex(.15,.875,"0.204 < x_{B} < 0.268, 1.912 < Q^{2} < 2.510 GeV^{2}, 0.250 < |t| < 0.400 GeV^{2}");cv.SaveAs((fs::path(out)/"bin_centering_systematic_phi_example.png").string().c_str());}
  return true;
 }catch(const std::exception&e){std::cerr<<"[bin-centering-systematic-plots] ERROR: "<<e.what()<<"\n";return false;}
}
