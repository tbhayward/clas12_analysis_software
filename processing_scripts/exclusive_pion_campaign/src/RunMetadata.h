#pragma once
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cmath>
// The shared CSV mixes six-column RGC polarization records with five-column
// legacy charge records. Only a requested polarized-target row is relevant.
inline double targetPolarization(const std::string& filename,const std::string& target,int expectedRun){
 if(target!="ND3"){
  if(target!="D2"&&target!="C"&&target!="CH2"&&target!="CD2"&&target!="ET"&&target!="He")throw std::runtime_error("Unknown target: "+target);
  return 0.;
 }
 std::ifstream file(filename);if(!file)throw std::runtime_error("Missing run metadata CSV");
 bool found=false;double result=0.;std::string line;int lineNumber=0;
 while(std::getline(file,line)){
  ++lineNumber;if(line.find_first_not_of(" \t\r")==std::string::npos)continue;
  if(line[line.find_first_not_of(" \t\r")]=='#')continue;
  std::istringstream row(line);std::string cell;std::getline(row,cell,',');size_t used=0;int run;
  try{run=std::stoi(cell,&used);}catch(...){throw std::runtime_error("Invalid run number at CSV line "+std::to_string(lineNumber));}
  if(cell.find_first_not_of(" \t\r",used)!=std::string::npos)throw std::runtime_error("Invalid CSV run field");
  if(run!=expectedRun)continue;
  std::string cells[5];for(int i=0;i<5;i++)if(!std::getline(row,cells[i],','))throw std::runtime_error("Missing ND3 polarization columns for run "+std::to_string(run));
  if(std::getline(row,cell,','))throw std::runtime_error("Extra ND3 metadata columns");
  double p;try{p=std::stod(cells[3],&used);}catch(...){throw std::runtime_error("Invalid ND3 target polarization");}
  if(cells[3].find_first_not_of(" \t\r",used)!=std::string::npos||!std::isfinite(p)||std::abs(p)>1)throw std::runtime_error("Invalid ND3 target polarization");
  if(found&&p!=result)throw std::runtime_error("Conflicting ND3 polarization for run "+std::to_string(run));
  found=true;result=p;
 }
 if(!found)throw std::runtime_error("Missing ND3 target polarization for run "+std::to_string(expectedRun));return result;
}
