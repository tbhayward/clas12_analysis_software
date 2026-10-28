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
 std::vector<double> v(n);std::vector<std::string> names(n);
 for(size_t i=0;i<n;++i){char buf[16];std::snprintf(buf,sizeof(buf),"col%02zu",i+1);names[i]=buf;tree.Branch(names[i].c_str(),&v[i],(names[i]+"/D").c_str());}
 // Explicit alias: original legacy data column 33 is proton-derived t1.
 tree.SetAlias("t1","col33");tree.SetAlias("t","col32");tree.SetAlias("t2","col34");
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
