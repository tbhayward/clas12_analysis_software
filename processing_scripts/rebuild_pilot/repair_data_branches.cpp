#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TObject.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdio>
int main(int argc,char**argv){
 if(argc!=3){std::cerr<<"usage: repair_data_branches file.root schema.txt\n";return 2;}
 std::ifstream schema(argv[2]);if(!schema){std::cerr<<"Cannot read schema\n";return 2;}
 std::vector<std::string> names;std::string n;while(schema>>n)names.push_back(n);
 TFile f(argv[1],"UPDATE");if(f.IsZombie()){std::cerr<<"Cannot open ROOT file\n";return 3;}
 TTree *t=nullptr;f.GetObject("DataEvents",t);if(!t)f.GetObject("PhysicsEvents",t);if(!t){std::cerr<<"Missing DataEvents/PhysicsEvents\n";return 4;}
 if((size_t)t->GetListOfBranches()->GetEntries()!=names.size()){std::cerr<<"Unexpected branch count\n";return 5;}
 bool old=true,newfmt=true;
 for(size_t i=0;i<names.size();++i){
  TBranch *b=(TBranch*)t->GetListOfBranches()->At(i);
  char buf[32];std::snprintf(buf,sizeof(buf,"col%02zu",i+1);
  if(std::string(b->GetName())!=buf)old=false;
  if(std::string(b->GetName())!=names[i])newfmt=false;
 }
 if(newfmt && std::string(t->GetName())=="PhysicsEvents"){std::cout<<"[ALREADY_COMPATIBLE] "<<argv[1]<<"\n";return 0;}
 if(!old && !newfmt){std::cerr<<"Unexpected or mixed branch names; refusing to change\n";return 6;}
 auto entries=t->GetEntries();
 if(old)for(size_t i=0;i<names.size();++i){
  TBranch *b=(TBranch*)t->GetListOfBranches()->At(i);
  b->SetName(names[i].c_str());
  if(b->GetListOfLeaves()->GetEntries()!=1){std::cerr<<"Unexpected leaf count\n";return 7;}
  TLeaf *leaf=(TLeaf*)b->GetListOfLeaves()->At(0);leaf->SetName(names[i].c_str());
 }
 if(t->GetListOfAliases())t->GetListOfAliases()->Clear();
 t->SetName("PhysicsEvents");
 f.cd();if(t->Write("",TObject::kOverwrite)<=0){std::cerr<<"Failed to write tree header\n";return 8;}
 f.Close();
 TFile verify(argv[1],"READ");TTree *check=nullptr;verify.GetObject("PhysicsEvents",check);
 if(!check||check->GetEntries()!=entries){std::cerr<<"Verification failed: entry count\n";return 9;}
 for(auto &name:names)if(!check->GetBranch(name.c_str())){std::cerr<<"Verification failed: "<<name<<"\n";return 10;}
 std::cout<<"[REPAIRED] "<<argv[1]<<" entries="<<entries<<"\n";
 return 0;
}
