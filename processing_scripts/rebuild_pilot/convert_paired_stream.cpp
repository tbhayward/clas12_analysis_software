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
// Stream one row per HIPO event from a named pipe, no intermediate TSV on disk.
int main(int argc, char** argv) {
 if(argc!=3){std::cerr<<"usage: convert_paired_stream input.fifo output.root\n";return 2;}
 std::ifstream in(argv[1]); if(!in){std::cerr<<"Cannot open stream\n";return 2;}
 std::string header; if(!std::getline(in,header)){std::cerr<<"No header\n";return 2;}
 std::istringstream hs(header);std::vector<std::string> names;std::string name;
 while(hs>>name)names.push_back(name);
 if(names.empty()){std::cerr<<"Empty schema\n";return 2;}
 TFile out(argv[2],"RECREATE");if(out.IsZombie())return 3;
 TTree tree("MCEvents","Event-aligned generated and reconstructed MC");
 std::vector<double> values(names.size());
 for(size_t i=0;i<names.size();++i)tree.Branch(names[i].c_str(),&values[i],(names[i]+"/D").c_str());
 tree.SetAutoFlush(100000);
 std::string line; unsigned long long count=0;bool bad=false;
 while(std::getline(in,line)) {
  if(line.empty())continue;
  const char* p=line.c_str();char* end=nullptr;
  for(size_t i=0;i<values.size();++i){
   while(*p==' '||*p=='\t')++p;
   if(!*p){bad=true;break;}
   errno=0; values[i]=std::strtod(p,&end);
   if(end==p){bad=true;break;}
   p=end;
  }
  while(*p==' '||*p=='\t'||*p=='\r')++p;
  if(*p)bad=true;
  if(bad){std::cerr<<"Malformed row at "<<(count+1)<<"\n";break;}
  tree.Fill();++count;
 }
 if(in.bad()) {std::cerr<<"Input read error\n";bad=true;}
 if(!bad){tree.Write("",TObject::kOverwrite);out.Close();}
 else {out.Close();std::remove(argv[2]);return 4;}
 std::cout<<"[PAIRED_ROOT_ROWS] "<<count<<"\n";
 return 0;
}
