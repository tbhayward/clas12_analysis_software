#include "dilution_factor.h"
#include <cmath>
#include <iostream>

double dilution_factor(double currentVariable, const std::string& prefix) {

  // epi+X
  if (prefix == "epipPT") {
    return 0.184542-0.0499585*currentVariable+0.163844*std::pow(currentVariable,2)+
      0.157106*std::pow(currentVariable,3);
  }

  if (prefix == "x") { return 0.086655+0.0332876*currentVariable+0.268944*std::pow(currentVariable,2); }
  if (prefix == "z") { return -0.0728496+1.3191*currentVariable+-2.86053*std::pow(currentVariable,2)+1.96652*std::pow(currentVariable,3); }
  if (prefix == "zeta") { return 1.11574+-5.33344*currentVariable+9.55552*std::pow(currentVariable,2)+-5.80119*std::pow(currentVariable,3); }
  if (prefix == "PT") { return 0.106575+0.0128819*currentVariable+-0.014333*std::pow(currentVariable,2); }
  // if (prefix == "xF") { return 1; }
  if (prefix == "xF") { return 0.11673+-0.0406361*currentVariable+-0.106417*std::pow(currentVariable,2)+0.151004*std::pow(currentVariable,3); }
  if (prefix == "Mx") { return 0.0366071*std::exp(-0.5*std::pow((currentVariable-0.769794)/0.0562864,2)) + 0.0260418*std::exp(-0.5*std::pow((currentVariable-1.1057)/0.290501,2)) + 0.10761*std::exp(-0.5*std::pow((currentVariable-1.68863)/1.8129,2)); }
  
  if (prefix == "Q2y4z1") { return 0.0667422; }
  if (prefix == "Q2y4z2") { return 0.131855; }
  if (prefix == "Q2y4z3") { return 0.136825; }
  if (prefix == "Q2y4z4") { return 0.135157; }
  if (prefix == "Q2y4z5") { return 0.149875; }
  if (prefix == "Q2y8z1") { return 0.013595; }
  if (prefix == "Q2y8z2") { return 0.109565; }
  if (prefix == "Q2y8z3") { return 0.139691; }
  if (prefix == "Q2y8z4") { return 0.136057; }
  if (prefix == "Q2y8z5") { return 0.134013; }
  if (prefix == "Q2y12z1") { return 0; }
  if (prefix == "Q2y12z2") { return 0.174856; }
  if (prefix == "Q2y12z3") { return 0.165048; }
  if (prefix == "Q2y12z4") { return 0.140834; }
  if (prefix == "Q2y12z5") { return 0.150643; }
  if (prefix == "Q2y3z1") { return 0.0984995; }
  if (prefix == "Q2y3z2") { return 0.123269; }
  if (prefix == "Q2y3z3") { return 0.123317; }
  if (prefix == "Q2y3z4") { return 0.1172; }
  if (prefix == "Q2y3z5") { return 0.112746; }
  if (prefix == "Q2y7z1") { return 0.0845784; }
  if (prefix == "Q2y7z2") { return 0.115614; }
  if (prefix == "Q2y7z3") { return 0.120013; }
  if (prefix == "Q2y7z4") { return 0.112607; }
  if (prefix == "Q2y7z5") { return 0.0991583; }
  if (prefix == "Q2y11z1") { return 0.0950531; }
  if (prefix == "Q2y11z2") { return 0.131137; }
  if (prefix == "Q2y11z3") { return 0.14558; }
  if (prefix == "Q2y11z4") { return 0.13086; }
  if (prefix == "Q2y11z5") { return 0.102845; }
  if (prefix == "Q2y15z1") { return 0.0517377; }
  if (prefix == "Q2y15z2") { return 0.125339; }
  if (prefix == "Q2y15z3") { return 0.158101; }
  if (prefix == "Q2y15z4") { return 0.150249; }
  if (prefix == "Q2y15z5") { return 0.115582; }
  if (prefix == "Q2y2z1") { return 0.0995716; }
  if (prefix == "Q2y2z2") { return 0.117411; }
  if (prefix == "Q2y2z3") { return 0.102522; }
  if (prefix == "Q2y2z4") { return 0.0972157; }
  if (prefix == "Q2y2z5") { return 0.0732962; }
  if (prefix == "Q2y6z1") { return 0.0896362; }
  if (prefix == "Q2y6z2") { return 0.111472; }
  if (prefix == "Q2y6z3") { return 0.103157; }
  if (prefix == "Q2y6z4") { return 0.0850751; }
  if (prefix == "Q2y6z5") { return 0.0729225; }
  if (prefix == "Q2y10z1") { return 0.104736; }
  if (prefix == "Q2y10z2") { return 0.127091; }
  if (prefix == "Q2y10z3") { return 0.119057; }
  if (prefix == "Q2y10z4") { return 0.0845961; }
  if (prefix == "Q2y10z5") { return 0.0584329; }
  if (prefix == "Q2y14z1") { return 0.107384; }
  if (prefix == "Q2y14z2") { return 0.132264; }
  if (prefix == "Q2y14z3") { return 0.142461; }
  if (prefix == "Q2y14z4") { return 0.127831; }
  if (prefix == "Q2y14z5") { return 0.107375; }
  if (prefix == "Q2y17z1") { return 0.114428; }
  if (prefix == "Q2y17z2") { return 0.157584; }
  if (prefix == "Q2y17z3") { return 0.167516; }
  if (prefix == "Q2y17z4") { return 0.193292; }
  if (prefix == "Q2y17z5") { return 0.114817; }
  if (prefix == "Q2y1z1") { return 0.0890598; }
  if (prefix == "Q2y1z2") { return 0.101465; }
  if (prefix == "Q2y1z3") { return 0.0787156; }
  if (prefix == "Q2y1z4") { return 0.0692479; }
  if (prefix == "Q2y1z5") { return 0.0368427; }
  if (prefix == "Q2y5z1") { return 0.0880512; }
  if (prefix == "Q2y5z2") { return 0.101888; }
  if (prefix == "Q2y5z3") { return 0.0848115; }
  if (prefix == "Q2y5z4") { return 0.0736795; }
  if (prefix == "Q2y5z5") { return 0.0198957; }
  if (prefix == "Q2y9z1") { return 0.0979945; }
  if (prefix == "Q2y9z2") { return 0.112202; }
  if (prefix == "Q2y9z3") { return 0.0881914; }
  if (prefix == "Q2y9z4") { return 0.0896477; }
  if (prefix == "Q2y9z5") { return 0.0489943; }
  if (prefix == "Q2y13z1") { return 0.112754; }
  if (prefix == "Q2y13z2") { return 0.136534; }
  if (prefix == "Q2y13z3") { return 0.122508; }
  if (prefix == "Q2y13z4") { return 0.131878; }
  if (prefix == "Q2y13z5") { return 0.123804; }
  if (prefix == "Q2y16z1") { return 0.127336; }
  if (prefix == "Q2y16z2") { return 0.158212; }
  if (prefix == "Q2y16z3") { return 0.153284; }
  if (prefix == "Q2y16z4") { return 0.108524; }
  if (prefix == "Q2y16z5") { return 0.145373; }

  // epX
  if (prefix == "xF") { 
    return 0.186121-0.0263337*currentVariable-0.175587*std::pow(currentVariable,2)+
      0.0522814*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "xFpip") { 
    return 0.122453+0.189509*currentVariable-0.133621*std::pow(currentVariable,2)-
      0.0401427*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "xFpim") { 
    return 0.128348+0.195055*currentVariable-0.242617*std::pow(currentVariable,2)-
      0.204807*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "Mx") { 
    return 0.0847657+0.0762168*currentVariable-0.0128988*std::pow(currentVariable,2)+
      0.00274429*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "Q2TFR") {
    return 0.0884319+0.0414953*currentVariable-0.00584857*std::pow(currentVariable,2)+
      0.000500127*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "Q2bin") {
    return -0.341032+0.762811*currentVariable-0.399944*std::pow(currentVariable,2)+
      0.0686534*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "xTFR") {
    return 0.111702+0.0858432*currentVariable+0.880331*std::pow(currentVariable,2)-
      0.990298*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "xTFRpip") {
    return 0.117706-0.194421*currentVariable+0.977489*std::pow(currentVariable,2)-
      0.926193*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "xTFRpim") {
    return 0.0787795-0.263136*currentVariable+1.378*std::pow(currentVariable,2)-
      1.65335*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "PTTFR") {
    return 0.184491-0.161007*currentVariable+0.298733*std::pow(currentVariable,2)-
      0.187826*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "PTTFRpip") {
    return 0.176079-0.328598*currentVariable+0.475598*std::pow(currentVariable,2)-
      0.167004*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "PTTFRpim") {
    return 0.0275594+0.276354*currentVariable-0.471179*std::pow(currentVariable,2)-
      0.21497*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "zetaTFR") {
    return 1.52544-7.07359*currentVariable+12.5954*std::pow(currentVariable,2)-
      7.72548*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "zTFRpip") {
    return 0.0565765+0.882732*currentVariable-3.33409*std::pow(currentVariable,2)+
      5.51154*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "zTFRpim") {
    return -0.0253779+1.62183*currentVariable-6.76455*std::pow(currentVariable,2)+
      8.56005*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "Q2TFR") {
    return 0.093586+0.0370678*currentVariable-0.00373394*std::pow(currentVariable,2)+
      0.000215739*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "xCFR") {
    return 0.089331+0.429008*currentVariable-0.617364*std::pow(currentVariable,2)+
      0.7584*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "xCFRpip") {
    return 0.119971+0.416041*currentVariable-0.922544*std::pow(currentVariable,2)+
      1.01908*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "xCFRpim") {
    return 0.121553-0.12187*currentVariable+0.923064*std::pow(currentVariable,2)-
      0.949773*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "PTCFR") {
    return 0.151263+0.170759*currentVariable-0.439815*std::pow(currentVariable,2)+
      0.278509*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "PTCFRpip") {
    return 0.184542-0.0499585*currentVariable+0.163844*std::pow(currentVariable,2)+
      0.157106*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "PTCFRpim") {
    return 0.147254-0.134125*currentVariable+0.407317*std::pow(currentVariable,2)-
      0.339619*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "zetaCFR") {
    return 1.32783-6.22826*currentVariable+11.2985*std::pow(currentVariable,2)-
      7.01171*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "zCFRpip") {
    return 1.32783-6.22826*currentVariable+11.2985*std::pow(currentVariable,2)-
      7.01171*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "zCFRpim") {
    return 0.0997319+0.28069*currentVariable-0.547782*std::pow(currentVariable,2)+
      0.244802*std::pow(currentVariable,3);
  }
  // eX
  if (prefix == "xeX") {
    return 0.111702+0.0858432*currentVariable+0.880331*std::pow(currentVariable,2)-
      0.990298*std::pow(currentVariable,3);
  }

  // rho
  if (prefix == "rho") {
    return 0.913-0.192*currentVariable+0.104*std::pow(currentVariable,2);
  }
  // if (prefix == "rho") {
  //   return 0.882+0.033*currentVariable-0.060*std::pow(currentVariable,2);
  // }

  // epi+pX
  if (prefix == "b2banalysisx") {
    return 0.123+0.184*currentVariable-0.095*std::pow(currentVariable,2);
  }
  // epi+pX
  if (prefix == "b2banalysispTpT") {
    return 0.163-0.143*currentVariable+0.319*std::pow(currentVariable,2);
  }
  return 0.12;
}