#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cerrno>
#include <cmath>
#include <cstring>
#include "reconstructed_data_schema.h"
// Stream one row per HIPO event from a named pipe, no intermediate TSV on disk.
int main(int argc, char** argv) {
 if(argc!=3){std::cerr<<"usage: convert_paired_stream input.fifo output.root\n";return 2;}
 std::ifstream in(argv[1]); if(!in){std::cerr<<"Cannot open stream\n";return 2;}
 std::string header; if(!std::getline(in,header)){std::cerr<<"No header\n";return 2;}
 std::istringstream hs(header);std::vector<std::string> names;std::string name;
 while(hs>>name)names.push_back(name);
 const auto* base = names.size()==146 ? &reconstructed_data_schema::dvcs() : names.size()==160 ? &reconstructed_data_schema::pi0() : nullptr;
 if(!base){std::cerr<<"Expected 146 or 160 paired columns\n";return 2;}
 std::vector<std::string> fields=*base;fields.push_back("weight");
 std::vector<std::string> expected={"event_sequence","weight","rec_kinematics_valid","gen_kinematics_valid"};
 for(const auto& prefix:{"gen_","rec_"})for(const auto& field:fields)expected.push_back(std::string(prefix)+field);
 if(names!=expected){std::cerr<<"Paired schema names/order mismatch\n";return 2;}
 TFile out(argv[2],"RECREATE");if(out.IsZombie())return 3;
 TTree tree("MCEvents","Event-aligned generated and reconstructed MC");
 std::vector<double> values(names.size());
 TTree reconstructed("PhysicsEvents","Valid reconstructed MC candidates; no additional analysis cuts");
 const size_t recOffset=4+fields.size();
 for(size_t i=0;i<fields.size();++i)reconstructed.Branch(fields[i].c_str(),&values[recOffset+i],(fields[i]+"/D").c_str());
 reconstructed.Branch("event_sequence",&values[0],"event_sequence/D");
 reconstructed.SetAutoFlush(100000);
 for(size_t i=0;i<names.size();++i)tree.Branch(names[i].c_str(),&values[i],(names[i]+"/D").c_str());
 tree.SetAutoFlush(100000);
 std::string line; unsigned long long count=0,recCount=0;bool bad=false;
 while(std::getline(in,line)) {
  if(line.empty())continue;
  const char* p=line.c_str();char* end=nullptr;
  for(size_t i=0;i<values.size();++i){
   while(*p==' '||*p=='\t')++p;
   if(!*p){bad=true;break;}
   errno=0; values[i]=std::strtod(p,&end);
   if(end==p||(*end && *end!=' ' && *end!='\t' && *end!='\r')){bad=true;break;}
   p=end;
  }
  while(*p==' '||*p=='\t'||*p=='\r')++p;
  if(*p)bad=true;
  if(bad){std::cerr<<"Malformed row at "<<(count+1)<<"\n";break;}
  if(values[0]!=static_cast<double>(count+1) || (values[2]!=0 && values[2]!=1) || (values[3]!=0 && values[3]!=1) || !std::isfinite(values[1])){
   std::cerr<<"Invalid sequence, validity flag or event weight at row "<<count+1<<"\n";bad=true;break;
  }
  if((values[3]==1 && values[4+fields.size()-1]!=values[1]) || (values[2]==1 && values[recOffset+fields.size()-1]!=values[1])){
   std::cerr<<"Inconsistent generator weights at row "<<count+1<<"\n";bad=true;break;
  }
  if(tree.Fill()<0){bad=true;break;}
  if(values[2]==1){if(reconstructed.Fill()<0){bad=true;break;}++recCount;}
  ++count;
 }
 if(in.bad()) {std::cerr<<"Input read error\n";bad=true;}
 if(!bad){if(tree.Write("",TObject::kOverwrite)<=0 || reconstructed.Write("",TObject::kOverwrite)<=0)bad=true;out.Close();}
 else {out.Close();std::remove(argv[2]);return 4;}
 if(bad){std::remove(argv[2]);return 4;}
 std::cout<<"[PAIRED_ROOT_ROWS] "<<count<<"\n";
 std::cout<<"[RECONSTRUCTED_ROOT_ROWS] "<<recCount<<"\n";
 return 0;
}
