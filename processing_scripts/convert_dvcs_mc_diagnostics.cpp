#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <cstdint>
#include <cmath>
#include <iostream>
int main(int argc,char** argv){
 if(argc!=3){std::cerr<<"Usage: convert_dvcs_mc_diagnostics input.tsv output.root\n";return 2;}
 std::ifstream in(argv[1]); if(!in){std::cerr<<"Cannot read "<<argv[1]<<"\n";return 2;}
 std::string header;std::getline(in,header);std::istringstream hs(header);std::vector<std::string> names;std::string n;
 while(hs>>n) names.push_back(n);
 if(names.empty()){std::cerr<<"Empty header\n";return 2;}
 TFile out(argv[2],"RECREATE");TTree tree("MCEvents","Generated and truth matched MC diagnostic rows");
 std::vector<double> values(names.size(),0.);
 std::vector<Long64_t> integerValues(names.size(),0);
 const std::set<std::string> integerFields = {
   "file_index", "event_sequence", "run", "event", "n_mc", "n_rec", "n_match",
   "rec_final_state", "correction_run", "rec_kinematics_valid", "gen_kinematics_valid", "mcindex", "truth_pid", "rec_pindex",
   "rec_pid", "rec_status", "rec_detector", "fitted_selected"
 };
 for(size_t i=0;i<names.size();++i){
   if(integerFields.count(names[i])) tree.Branch(names[i].c_str(),&integerValues[i],(names[i]+"/L").c_str());
   else tree.Branch(names[i].c_str(),&values[i],(names[i]+"/D").c_str());
 }
 std::string line;size_t rows=0,bad=0;
 while(std::getline(in,line)){
  if(line.empty())continue;std::istringstream ls(line);bool valid=true;
  for(auto &v:values){std::string token;if(!(ls>>token)){valid=false;break;}try{size_t consumed=0;v=std::stod(token,&consumed);if(consumed!=token.size()){valid=false;break;}}catch(...){valid=false;break;}}
  double extra;if(ls>>extra)valid=false;
  if(!valid){++bad;continue;}
  for(size_t i=0;i<names.size();++i){
    if(integerFields.count(names[i])){
      if(!std::isfinite(values[i]) || values[i] != std::trunc(values[i]) ||
         std::fabs(values[i]) >= 9007199254740992.0){valid=false;break;}
      integerValues[i] = static_cast<Long64_t>(values[i]);
    }
  }
  if(!valid){++bad;continue;}
  tree.Fill();++rows;
 }
 tree.Write();out.Close();std::cout<<"Rows: "<<rows<<", malformed: "<<bad<<"\n";return bad?1:0;
}
