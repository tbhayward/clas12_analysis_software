#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
// Read metadata for every file in one ROOT process. No full event scan.
int main(int argc,char** argv){
 if(argc!=2)return 2;
 std::ifstream in(argv[1]);std::string line;std::getline(in,line);
 std::vector<std::pair<std::string,std::string>> schema;std::istringstream header(line);std::string field;
 while(header>>field){auto p=field.find(':');schema.emplace_back(field.substr(0,p),field.substr(p+1));}
 int files=0,failures=0;
 while(std::getline(in,line)){
  std::istringstream row(line);std::string path,expectedText,label;std::getline(row,path,'\t');std::getline(row,expectedText,'\t');std::getline(row,label);
  try{
   long long expected=std::stoll(expectedText);
   std::unique_ptr<TFile> f(TFile::Open(path.c_str(),"READ"));
   if(!f||f->IsZombie()||f->TestBit(TFile::kRecovered))throw std::runtime_error("missing, damaged or recovered ROOT file");
   auto* tree=dynamic_cast<TTree*>(f->Get("PhysicsEvents"));
   if(!tree||tree->GetEntries()!=expected)throw std::runtime_error("tree/entry count mismatch");
   if(tree->GetListOfBranches()->GetEntries()!=int(schema.size()))throw std::runtime_error("branch count mismatch");
   for(const auto& s:schema){auto* branch=tree->GetBranch(s.first.c_str());auto* leaf=tree->GetLeaf(s.first.c_str());if(!branch||!leaf||s.second!=leaf->GetTypeName())throw std::runtime_error("branch name/type mismatch: "+s.first);}
   files++;
  }catch(const std::exception& error){failures++;std::cout<<"ERROR "<<label<<": "<<error.what()<<" ["<<path<<"]\n";}
 }
 std::cout<<"ROOT metadata: "<<files<<" files passed; "<<failures<<" failed\n";
 return failures?1:0;
}
