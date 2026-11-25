#pragma once
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <set>
#include <vector>
inline void parseLooseRow(const std::string& line,const std::vector<std::string>& names,std::vector<double>& d,std::vector<int>& ints,const std::set<std::string>& integerNames,int expectedRun){
 std::istringstream row(line);
 const std::set<std::string> optional={"vz_e","vz_p1","vz_p2","chi2pid1","chi2pid2","beta1","beta2","chi2pid_e","beta_e"};
 for(size_t i=0;i<names.size();i++){
  std::string token;if(!(row>>token))throw std::runtime_error("Missing text column");size_t used=0;d[i]=std::stod(token,&used);
  bool nullable=(i>=28&&i<55)||optional.count(names[i]);
  if(used!=token.size()||std::isinf(d[i])||(std::isnan(d[i])&&!nullable))throw std::runtime_error("Nonfinite or malformed "+names[i]);
  if(integerNames.count(names[i])){if(!std::isfinite(d[i])||d[i]!=std::trunc(d[i])||std::abs(d[i])>2147483647.)throw std::runtime_error("Invalid integer "+names[i]);ints[i]=int(d[i]);}
  if((names[i].rfind("pass_",0)==0||names[i]=="electron_correction"||names[i]=="proton_correction")&&(ints[i]<-1||ints[i]>1))throw std::runtime_error("Invalid selection flag");
 }
 std::string extra;if(row>>extra)throw std::runtime_error("Extra text columns");
 if(ints[0]!=expectedRun||std::abs(ints[2])!=1||ints[6]<0||ints[6]>3||ints[7]<0||ints[7]>3)throw std::runtime_error("Invalid run/helicity/detector");
 double rawP2=0;for(const std::string name:{"e_raw_px","e_raw_py","e_raw_pz"}){auto k=std::find(names.begin(),names.end(),name)-names.begin();rawP2+=d[k]*d[k];}
 if(!(rawP2>4))throw std::runtime_error("Electron momentum must exceed 2 GeV");
 size_t valid=std::find(names.begin(),names.end(),"kinematics_valid")-names.begin();if(ints[valid]!=0&&ints[valid]!=1)throw std::runtime_error("Invalid kinematics flag");
 if(ints[valid])for(size_t i=28;i<55;i++)if(!std::isfinite(d[i]))throw std::runtime_error("Inconsistent kinematics_valid");
}
