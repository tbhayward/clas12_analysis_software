#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TROOT.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <filesystem>
namespace fs=std::filesystem;
int main(int argc,char**argv){
 if(argc!=4){std::cerr<<"Usage: validate_data_root manifest.tsv report.csv plots_dir\n";return 2;}
 gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
 std::ifstream in(argv[1]);if(!in){std::cerr<<"Cannot read manifest\n";return 2;}
 std::ofstream out(argv[2]);if(!out){std::cerr<<"Cannot write report\n";return 2;}
 fs::create_directories(argv[3]);
 out<<"mode,run,root_file,expected_rows,actual_rows,branches,nonfinite_kinematics,invalid_run,invalid_event,negative_q2,negative_w,negative_xb,negative_t1,negative_t,nonzero_t1,invalid_alias,status\n";
 struct H {std::unique_ptr<TH1D> q2,w,xb,t,t1,phi;};
 auto make=[&](const std::string&m){H h;h.q2=std::make_unique<TH1D>((m+"_q2").c_str(),"Q^{2};Q^{2} (GeV^{2});Candidates",100,0,12);h.w=std::make_unique<TH1D>((m+"_w").c_str(),"W;W (GeV);Candidates",100,0,6);h.xb=std::make_unique<TH1D>((m+"_xb").c_str(),"x_{B};x_{B};Candidates",100,0,1);h.t=std::make_unique<TH1D>((m+"_t").c_str(),"|t| from photon;|t| (GeV^{2});Candidates",100,0,2);h.t1=std::make_unique<TH1D>((m+"_t1").c_str(),"|t| from proton;|t_{p}| (GeV^{2});Candidates",100,0,2);h.phi=std::make_unique<TH1D>((m+"_phi").c_str(),"#phi_{h};#phi_{h} (rad);Candidates",100,-3.3,6.4);return h;};
 H dvcs=make("dvcs"),pi0=make("pi0");
 std::string mode,run,path,expectedstr;long long files=0,bad=0,entries=0;
 while(in>>mode>>run>>expectedstr>>path){
  ++files;long long expected=std::stoll(expectedstr);long long n=-1,nonfinite=0,badrun=0,badev=0,negq=0,negw=0,negxb=0,negt1=0,negt=0,nonzerot1=0,invalidalias=0;int nb=-1;std::string status="OK";
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(),"READ"));TTree*tree=nullptr;
  if(!f||f->IsZombie())status="OPEN_FAILED";
  else {tree=nullptr;f->GetObject("DataEvents",tree);if(!tree)status="MISSING_TREE";}
  const int required=mode=="dvcs_data"?70:77;
  if(tree){n=tree->GetEntries();nb=tree->GetListOfBranches()->GetEntries();if(n!=expected)status="ROW_MISMATCH";if(nb!=required)status="BRANCH_MISMATCH";
   std::vector<TLeaf*> leaf(required,nullptr);for(int i=0;i<required;i++){char name[20];std::snprintf(name,sizeof(name),"col%02d",i+1);leaf[i]=tree->GetLeaf(name);if(!leaf[i])status="MISSING_COLUMN";}
   const char* aliases[]={"t","t1","t2"};const char* names[]={"col32","col33","col34"};for(int i=0;i<3;i++){auto*alias=tree->GetAlias(aliases[i]);if(!alias||std::string(alias)!=names[i])++invalidalias;}if(invalidalias)status="BAD_ALIAS";
   if(status!="MISSING_COLUMN"){
    H&h=mode=="dvcs_data"?dvcs:pi0;
    for(long long j=0;j<n;j++){
     tree->GetEntry(j);auto get=[&](int c){return leaf[c-1]->GetValue();};
     double q=get(26),w=get(27),xb=get(31),t=get(32),t1=get(33),ph=get(60),rn=get(5),ev=get(6);
     if(!std::isfinite(q)||!std::isfinite(w)||!std::isfinite(xb)||!std::isfinite(t)||!std::isfinite(t1)||!std::isfinite(ph))++nonfinite;
     if(!std::isfinite(rn)||std::llround(rn)!=std::stoll(run))++badrun;
     if(!std::isfinite(ev)||ev<0)++badev;
     if(std::isfinite(q)&&q<0)++negq;if(std::isfinite(w)&&w<0)++negw;if(std::isfinite(xb)&&xb<0)++negxb;
     if(std::isfinite(t1)&&t1<0)++negt1;if(std::isfinite(t)&&t<0)++negt;if(std::isfinite(t1)&&std::abs(t1)>1e-10)++nonzerot1;
     if(std::isfinite(q))h.q2->Fill(q);if(std::isfinite(w))h.w->Fill(w);if(std::isfinite(xb))h.xb->Fill(xb);if(std::isfinite(t))h.t->Fill(std::abs(t));if(std::isfinite(t1))h.t1->Fill(std::abs(t1));if(std::isfinite(ph))h.phi->Fill(ph);
    }
    if(nonfinite||badrun||badev)status="INVALID_VALUES";
    if(n>0&&nonzerot1==0)status="ALL_T1_ZERO";
   }
  }
  if(status!="OK")++bad;if(n>0)entries+=n;
  out<<mode<<','<<run<<','<<path<<','<<expected<<','<<n<<','<<nb<<','<<nonfinite<<','<<badrun<<','<<badev<<','<<negq<<','<<negw<<','<<negxb<<','<<negt1<<','<<negt<<','<<nonzerot1<<','<<invalidalias<<','<<status<<'\n';
  if(status!="OK")std::cerr<<"[ISSUE] "<<mode<<" "<<run<<" "<<status<<"\n";
 }
 for(auto pair : {std::pair<const char*,H*>("dvcs_data",&dvcs),{"pi0_data",&pi0}}){
  auto save=[&](TH1D*hist,const char*key){TCanvas c("c","",850,550);hist->SetLineWidth(2);hist->Draw("HIST");c.SaveAs((fs::path(argv[3])/(std::string(pair.first)+"_"+key+".png")).string().c_str());};
  save(pair.second->q2.get(),"q2");save(pair.second->w.get(),"w");save(pair.second->xb.get(),"xb");save(pair.second->t.get(),"t");save(pair.second->t1.get(),"t1_proton");save(pair.second->phi.get(),"phi");
 }
 std::cout<<"[ROOT_AUDIT] files="<<files<<" entries="<<entries<<" issues="<<bad<<" report="<<argv[2]<<" plots="<<argv[3]<<"\n";
 return bad?1:0;
}
