#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TKey.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <stdexcept>
#include <memory>
#include "beam_pol.h"
#include "RunMetadata.h"
using namespace std;
const string schema="runnum evnum helicity beam_energy p_index detector sector_e sector fiducial_status e_pcal_fiducial_1 e_pcal_fiducial_2 e_dc_fiducial pion_dc_fiducial rec_proton_count e_p e_theta e_phi vz_e p_p p_theta p_phi vz_p chi2pid beta Q2 W x y nu t tmin tprime Mx Mx2 Mx2_deuteron Mx2_proton Emiss missing_px missing_py missing_pz missing_p pTmiss z xF pT phi open_angle DepA DepB DepC DepV DepW";
vector<string> split(const string& s){istringstream in(s);vector<string> v;string x;while(in>>x)v.push_back(x);return v;}
void validate(const char* path,long long expected){unique_ptr<TFile> f(TFile::Open(path));if(!f||f->IsZombie()||f->TestBit(TFile::kRecovered))throw runtime_error("Invalid ROOT file");auto t=dynamic_cast<TTree*>(f->Get("PhysicsEvents"));if(!t||t->GetEntries()!=expected)throw runtime_error("Tree entry mismatch");auto names=split(schema+" beam_pol target_pol");if(t->GetListOfBranches()->GetEntries()!=int(names.size()))throw runtime_error("Unexpected branch count");for(auto& n:names)if(!t->GetBranch(n.c_str()))throw runtime_error("Missing branch: "+n);}
int main(int argc,char** argv){try{
 if(argc==4&&string(argv[1])=="--validate"){validate(argv[2],stoll(argv[3]));return 0;}
 if(argc!=6)throw runtime_error("Usage: converter input.txt output.root run_info.csv target expected_run; or --validate file.root entries");
 int expectedRun=stoi(argv[5]);string target=argv[4];
 double target_pol=targetPolarization(argv[3],target,expectedRun),beam_pol=getPol(expectedRun);
 string line;ifstream in(argv[1]);if(!in||!getline(in,line)||line!="# "+schema)throw runtime_error("Text schema mismatch");auto names=split(schema);vector<double> d(names.size());vector<int> ints(names.size());set<string> integerNames;for(auto& n:split("runnum evnum helicity p_index detector sector_e sector fiducial_status e_pcal_fiducial_1 e_pcal_fiducial_2 e_dc_fiducial pion_dc_fiducial rec_proton_count"))integerNames.insert(n);
 TFile out(argv[2],"RECREATE");if(out.IsZombie())throw runtime_error("Cannot create ROOT file");TTree tree("PhysicsEvents","Inclusive e' pi- candidates; missing proton; stationary neutron reference");
 for(size_t i=0;i<names.size();i++){string leaf=names[i]+(integerNames.count(names[i])?"/I":"/D");if(integerNames.count(names[i]))tree.Branch(names[i].c_str(),&ints[i],leaf.c_str());else tree.Branch(names[i].c_str(),&d[i],leaf.c_str());}
 tree.Branch("beam_pol",&beam_pol);tree.Branch("target_pol",&target_pol);long long rows=0;
 while(getline(in,line)){if(line.empty())throw runtime_error("Empty data row");istringstream row(line);for(size_t i=0;i<names.size();i++){if(!(row>>d[i])||!isfinite(d[i]))throw runtime_error("Malformed row "+to_string(rows+1));if(integerNames.count(names[i])){if(d[i]!=trunc(d[i])||abs(d[i])>2147483647.)throw runtime_error("Invalid integer");ints[i]=int(d[i]);}}string extra;if(row>>extra)throw runtime_error("Extra text columns");if(ints[0]!=expectedRun||abs(ints[2])!=1||(ints[5]<0||ints[5]>2)||ints[4]<1||ints[13]<0)throw runtime_error("Invalid run/helicity/detector");for(const string name:{"e_pcal_fiducial_1","e_pcal_fiducial_2","e_dc_fiducial","pion_dc_fiducial"}){size_t index=find(names.begin(),names.end(),name)-names.begin();if(ints[index]!=0&&ints[index]!=1)throw runtime_error("Invalid fiducial flag");}int status=ints[9]&&ints[10]&&ints[11]?(ints[12]?3:1):(ints[12]?0:-1);if(ints[8]!=status)throw runtime_error("Inconsistent fiducial status");if(tree.Fill()<0)throw runtime_error("ROOT Fill failed");rows++;}
 if(!in.eof())throw runtime_error("Input read failed");if(tree.Write()<=0)throw runtime_error("ROOT Write failed");TNamed convention("kinematic_convention","stationary neutron; t=(q-pi)^2; exact gamma*n->pi-p forward tmin; tprime=t-tmin; Mx2 uses stationary neutron; Mx2_deuteron uses full deuteron; electron PCAL2/DC applied through shared pi+ test; pion fiducials tagged only; no proton required; radians; opening angles degrees");convention.Write();tree.SetDirectory(nullptr);out.Close();validate(argv[2],rows);cout<<rows<<endl;return 0;
 }catch(const exception& e){cerr<<e.what()<<endl;return 1;}}
