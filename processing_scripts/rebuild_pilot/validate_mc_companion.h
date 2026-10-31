#ifndef VALIDATE_MC_COMPANION_H
#define VALIDATE_MC_COMPANION_H
#include "reconstructed_data_schema.h"
#include <TLeaf.h>
#include <cmath>
#include <iostream>
inline bool validate_mc_companion(TFile& file,TTree* paired){
 auto* rec=dynamic_cast<TTree*>(file.Get("PhysicsEvents"));
 if(!rec){std::cerr<<"MISSING_PhysicsEvents\n";return false;}
 const auto* base=paired->GetListOfBranches()->GetEntries()==146 ? &reconstructed_data_schema::dvcs() : paired->GetListOfBranches()->GetEntries()==160 ? &reconstructed_data_schema::pi0() : nullptr;
 if(!base)return false;
 std::vector<std::string> fields=*base;fields.push_back("weight");
 std::vector<std::string> expected={"event_sequence","weight","rec_kinematics_valid","gen_kinematics_valid"};
 for(const auto& prefix:{"gen_","rec_"})for(const auto& field:fields)expected.push_back(std::string(prefix)+field);
 auto validSchema=[](TTree* tree,const std::vector<std::string>& names){
  if(tree->GetListOfBranches()->GetEntries()!=static_cast<int>(names.size()))return false;
  for(size_t i=0;i<names.size();++i){
   if(names[i]!=tree->GetListOfBranches()->At(i)->GetName())return false;
   auto* leaf=tree->GetLeaf(names[i].c_str());
   if(!leaf || std::string(leaf->GetTypeName())!="Double_t" || leaf->GetLeafCount() || leaf->GetLenStatic()!=1)return false;
  }return true;
 };
 auto compact=fields;compact.push_back("event_sequence");
 if(!validSchema(paired,expected)||!validSchema(rec,compact)){std::cerr<<"BAD_MC_SCHEMA\n";return false;}
 Long64_t selected=0;
 for(Long64_t i=0;i<paired->GetEntries();++i){
  if(paired->GetEntry(i)<=0)return false;
  double g=paired->GetLeaf("gen_kinematics_valid")->GetValue(),r=paired->GetLeaf("rec_kinematics_valid")->GetValue();
  if((g!=0 && g!=1)||(r!=0 && r!=1))return false;
  if(r!=1)continue;
  if(selected>=rec->GetEntries()||rec->GetEntry(selected++)<=0)return false;
  for(const auto& name:compact){
   std::string source=name=="event_sequence"?name:"rec_"+name;
   double a=paired->GetLeaf(source.c_str())->GetValue(),b=rec->GetLeaf(name.c_str())->GetValue();
   if(a!=b && !(std::isnan(a)&&std::isnan(b))){std::cerr<<"COMPANION_VALUE_MISMATCH "<<name<<"\n";return false;}
  }
 }
 if(selected!=rec->GetEntries())return false;
 std::cout<<"COMPANION_VALID entries="<<selected<<"\n";return true;
}
#endif
