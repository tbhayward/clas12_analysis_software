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
using namespace std;
const string schema="runnum evnum helicity beam_energy p1_index p2_index detector1 detector2 sector_e sector1 sector2 fiducial_status e_p e_theta e_phi vz_e p1_p p1_theta p1_phi vz_p1 p2_p p2_theta p2_phi vz_p2 chi2pid1 chi2pid2 beta1 beta2 Q2 W x y nu t tmin tprime Mx2 Mx2_1 Mx2_2 Emiss pTmiss Mh z2 xF2 pT2 phi1 phi2 Delta_phi open_angle_ep1 open_angle_ep2 DepA DepB DepC DepV DepW";
vector<string> split(const string& s){istringstream in(s);vector<string> v;string x;while(in>>x)v.push_back(x);return v;}
void validate(const char* path,long long expected){unique_ptr<TFile> f(TFile::Open(path));if(!f||f->IsZombie()||f->TestBit(TFile::kRecovered))throw runtime_error("Invalid ROOT file");auto t=dynamic_cast<TTree*>(f->Get("PhysicsEvents"));if(!t||t->GetEntries()!=expected)throw runtime_error("Tree entry mismatch");auto names=split(schema+" beam_pol target_pol");if(t->GetListOfBranches()->GetEntries()!=int(names.size()))throw runtime_error("Unexpected branch count");for(auto& n:names)if(!t->GetBranch(n.c_str()))throw runtime_error("Missing branch: "+n);}
int main(int argc,char** argv){try{
 if(argc==4&&string(argv[1])=="--validate"){validate(argv[2],stoll(argv[3]));return 0;}
 if(argc!=6)throw runtime_error("Usage: converter input.txt output.root run_info.csv target expected_run; or --validate file.root entries");
 int expectedRun=stoi(argv[5]);string target=argv[4];
 map<int,double> pol;ifstream csv(argv[3]);if(!csv)throw runtime_error("Missing run metadata CSV");string line;
 while(getline(csv,line)){if(line.empty()||line[0]=='#')continue;istringstream row(line);vector<string> cells;string c;while(getline(row,c,','))cells.push_back(c);if(cells.size()!=6)throw runtime_error("Malformed run metadata");int run=stoi(cells[0]);double p=stod(cells[4]);if(!isfinite(p)||abs(p)>1)throw runtime_error("Invalid target polarization");if(pol.count(run)&&pol[run]!=p)throw runtime_error("Conflicting run metadata");pol[run]=p;}
 double target_pol=0,beam_pol=getPol(expectedRun);
 if(target=="ND3"){if(!pol.count(expectedRun))throw runtime_error("Missing ND3 target polarization for run "+to_string(expectedRun));target_pol=pol.at(expectedRun);}
 else if(target!="D2"&&target!="C"&&target!="CH2"&&target!="CD2"&&target!="ET"&&target!="He")throw runtime_error("Unknown target");
 ifstream in(argv[1]);if(!in||!getline(in,line)||line!="# "+schema)throw runtime_error("Text schema mismatch");auto names=split(schema);vector<double> d(names.size());vector<int> ints(names.size());set<string> integerNames;for(auto& n:split("runnum evnum helicity p1_index p2_index detector1 detector2 sector_e sector1 sector2 fiducial_status"))integerNames.insert(n);
 TFile out(argv[2],"RECREATE");if(out.IsZombie())throw runtime_error("Cannot create ROOT file");TTree tree("PhysicsEvents","Exclusive e' p pi-; stationary neutron reference");
 for(size_t i=0;i<names.size();i++){string leaf=names[i]+(integerNames.count(names[i])?"/I":"/D");if(integerNames.count(names[i]))tree.Branch(names[i].c_str(),&ints[i],leaf.c_str());else tree.Branch(names[i].c_str(),&d[i],leaf.c_str());}
 tree.Branch("beam_pol",&beam_pol);tree.Branch("target_pol",&target_pol);long long rows=0;
 while(getline(in,line)){if(line.empty())throw runtime_error("Empty data row");istringstream row(line);for(size_t i=0;i<names.size();i++){if(!(row>>d[i])||!isfinite(d[i]))throw runtime_error("Malformed row "+to_string(rows+1));if(integerNames.count(names[i])){if(d[i]!=trunc(d[i])||abs(d[i])>2147483647.)throw runtime_error("Invalid integer");ints[i]=int(d[i]);}}string extra;if(row>>extra)throw runtime_error("Extra text columns");if(ints[0]!=expectedRun||abs(ints[2])!=1||ints[7]!=1||(ints[6]!=1&&ints[6]!=2))throw runtime_error("Invalid run/helicity/detector");if(tree.Fill()<0)throw runtime_error("ROOT Fill failed");rows++;}
 if(!in.eof())throw runtime_error("Input read failed");if(tree.Write()<=0)throw runtime_error("ROOT Write failed");TNamed convention("kinematic_convention","stationary neutron; t=(q-pi)^2; exact gamma*n->pi-p forward tmin; tprime=t-tmin; radians; opening angles degrees");convention.Write();tree.SetDirectory(nullptr);out.Close();validate(argv[2],rows);cout<<rows<<endl;return 0;
 }catch(const exception& e){cerr<<e.what()<<endl;return 1;}}
