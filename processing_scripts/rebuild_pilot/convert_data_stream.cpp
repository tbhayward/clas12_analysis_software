#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cerrno>
#include <cstdio>
int main(int argc,char**argv){
 if(argc!=4){std::cerr<<"usage: convert_data_stream fifo output.root dvcs_data|pi0_data\n";return 2;}
 std::string mode=argv[3];size_t n=(mode=="dvcs_data"?70:mode=="pi0_data"?77:0);
 if(!n){std::cerr<<"Unknown mode\n";return 2;}
 std::ifstream in(argv[1]);if(!in){std::cerr<<"Cannot open FIFO\n";return 2;}
 TFile out(argv[2],"RECREATE");if(out.IsZombie())return 3;
 TTree tree("DataEvents","Original data processor selected-event columns");
 std::vector<double> v(n);
 const std::string common="fiducial_status num_pos num_neg num_neutrals runnum evnum helicity detector1 detector2 e_p e_theta e_phi vz_e p1_p p1_theta p1_phi vz_p1 p2_p p2_theta p2_phi vz_p2 open_angle_ep open_angle_ep1 open_angle_ep2 open_angle_p1p2 Q2 W Mx2 Mx2_1 Mx2_2 x t t1 t2 tmin y z z1 z2 Mh xF xF1 xF2 pT pT1 pT2 pTpT xi xi1 xi2 eta eta1 eta2 Delta_eta eta1_gN eta2_gN phi1 phi2 Delta_phi phih phiR theta DepA DepB DepC DepV DepW ";
 const std::string extra=(mode=="dvcs_data") ? "Emiss2 theta_gamma_gamma pTmiss" : "Mh_gammagamma detector_gamma1 detector_gamma2 open_angle_egamma1 open_angle_egamma2 gamma_phi1 gamma_phi2 Emiss2 theta_pi0_pi0 pTmiss";
 std::istringstream schema(common+extra);std::vector<std::string> names;std::string branch;
 while(schema>>branch)names.push_back(branch);
 if(names.size()!=n){std::cerr<<"Invalid schema size "<<names.size()<<" expected "<<n<<"\n";return 5;}
 for(size_t i=0;i<n;++i)tree.Branch(names[i].c_str(),&v[i],(names[i]+"/D").c_str());
 tree.SetAutoFlush(100000);
 std::string line;unsigned long long rows=0;bool bad=false;
 while(std::getline(in,line)){
  if(line.empty())continue;
  std::istringstream ss(line);std::string token;
  for(size_t i=0;i<n;++i){
   if(!(ss>>token)){bad=true;break;}
   char*end=nullptr;errno=0;v[i]=std::strtod(token.c_str(),&end);
   if(end==token.c_str()||*end!='\0'){bad=true;break;}
  }
  if(!bad && (ss>>token))bad=true;
  if(bad){std::cerr<<"Invalid column count or numeric token at selected row "<<(rows+1)<<" expected="<<n<<"\n";break;}
  tree.Fill();++rows;
 }
 if(in.bad())bad=true;
 if(!bad){tree.Write("",TObject::kOverwrite);out.Close();std::cout<<"[DATA_ROOT_ROWS] "<<rows<<"\n";return 0;}
 out.Close();std::remove(argv[2]);return 4;
}
