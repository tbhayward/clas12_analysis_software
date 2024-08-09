#include "dilution_factor.h"
#include <cmath>
#include <iostream>
#include <TRandom3.h>

double dilution_factor(double currentVariable, const std::string& prefix) {

  // epi+X
  if (prefix == "epipPT") {
    return 0.184542-0.0499585*currentVariable+0.163844*std::pow(currentVariable,2)+
      0.157106*std::pow(currentVariable,3);
  }

  if (prefix == "x") { double sigma = 0.001; return 0.143882+0.0753706*currentVariable+0.654032*std::pow(currentVariable,2) + rand_gen.Gaus(0, sigma); }
  if (prefix == "PT") { double sigma = 0.001; return 0.199431-0.223368*currentVariable+0.644003*std::pow(currentVariable,2) + rand_gen.Gaus(0, sigma); }
  if (prefix == "xF") { double sigma = 0.001; return 0.195272-0.0135855*currentVariable-0.0736621*std::pow(currentVariable,2) + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z1") { double sigma = 0.0413275; return 0.242197 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z2") { double sigma = 0.00377823; return 0.135962 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z3") { double sigma = 0.00844362; return 0.158756 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z4") { double sigma = 0.00862382; return 0.14849 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z5") { double sigma = 0.00981565; return 0.156161 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z1") { double sigma = 0.0379279; return 0.294189 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z2") { double sigma = 0.00327553; return 0.18422 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z3") { double sigma = 0.00349915; return 0.213911 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z4") { double sigma = 0.00850302; return 0.22171 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z5") { double sigma = 0.00999568; return 0.228514 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z1") { double sigma = -nan; return 0 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z2") { double sigma = 0.0135898; return 0.246862 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z3") { double sigma = 0.00992952; return 0.251256 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z4") { double sigma = 0.00620339; return 0.232663 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z5") { double sigma = 0.0150385; return 0.210134 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z1") { double sigma = 0.00694702; return 0.13002 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z2") { double sigma = 0.00590863; return 0.159932 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z3") { double sigma = 0.00580286; return 0.161149 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z4") { double sigma = 0.00968058; return 0.146144 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z5") { double sigma = 0.0228311; return 0.175958 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z1") { double sigma = 0.00501468; return 0.155728 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z2") { double sigma = 0.00451093; return 0.201032 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z3") { double sigma = 0.0042392; return 0.206787 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z4") { double sigma = 0.00502108; return 0.205502 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z5") { double sigma = 0.00998189; return 0.190105 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z1") { double sigma = 0.0113599; return 0.163101 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z2") { double sigma = 0.00368236; return 0.216486 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z3") { double sigma = 0.00531377; return 0.227377 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z4") { double sigma = 0.00659144; return 0.226683 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z5") { double sigma = 0.00756681; return 0.242188 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z1") { double sigma = 0.0209033; return 0.185026 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z2") { double sigma = 0.0052656; return 0.232319 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z3") { double sigma = 0.012072; return 0.251985 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z4") { double sigma = 0.0117963; return 0.247976 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z5") { double sigma = 0.0443307; return 0.273208 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z1") { double sigma = 0.00736045; return 0.148595 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z2") { double sigma = 0.00387207; return 0.175807 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z3") { double sigma = 0.00398186; return 0.160963 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z4") { double sigma = 0.0140496; return 0.144791 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z5") { double sigma = 0.0292463; return 0.18984 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z1") { double sigma = 0.00718851; return 0.171444 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z2") { double sigma = 0.00281353; return 0.198539 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z3") { double sigma = 0.00359738; return 0.190327 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z4") { double sigma = 0.00860447; return 0.190697 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z5") { double sigma = 0.022291; return 0.170989 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z1") { double sigma = 0.00397289; return 0.193468 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z2") { double sigma = 0.0024578; return 0.214237 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z3") { double sigma = 0.00494081; return 0.214983 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z4") { double sigma = 0.0114092; return 0.184261 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z5") { double sigma = 0.0414947; return 0.216189 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z1") { double sigma = 0.00738746; return 0.184706 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z2") { double sigma = 0.00366386; return 0.208434 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z3") { double sigma = 0.00967299; return 0.230743 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z4") { double sigma = 0.0176664; return 0.252315 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z5") { double sigma = 0.0413178; return 0.246462 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z1") { double sigma = 0.0254881; return 0.215869 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z2") { double sigma = 0.00560271; return 0.238701 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z3") { double sigma = 0.00965464; return 0.262723 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z4") { double sigma = 0.0158283; return 0.268868 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z5") { double sigma = 0.0549176; return 0.37505 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z1") { double sigma = 0.00571837; return 0.165969 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z2") { double sigma = 0.00399325; return 0.180655 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z3") { double sigma = 0.009156; return 0.166592 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z4") { double sigma = 0.0255769; return 0.224066 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z5") { double sigma = 0.0419718; return 0.171837 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z1") { double sigma = 0.00662353; return 0.178628 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z2") { double sigma = 0.00454322; return 0.196224 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z3") { double sigma = 0.00739584; return 0.179131 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z4") { double sigma = 0.0282124; return 0.193116 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z5") { double sigma = 0.0344613; return 0.24764 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z1") { double sigma = 0.00563494; return 0.18516 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z2") { double sigma = 0.00740232; return 0.196564 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z3") { double sigma = 0.00915024; return 0.177837 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z4") { double sigma = 0.0290215; return 0.214407 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z5") { double sigma = 0.0353604; return 0.296948 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z1") { double sigma = 0.00720393; return 0.19919 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z2") { double sigma = 0.00606249; return 0.228949 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z3") { double sigma = 0.0146543; return 0.219792 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z4") { double sigma = 0.0235241; return 0.201572 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z5") { double sigma = 0.0417879; return 0.421139 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z1") { double sigma = 0.0102113; return 0.208276 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z2") { double sigma = 0.00648511; return 0.248289 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z3") { double sigma = 0.0116444; return 0.237183 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z4") { double sigma = 0.0392161; return 0.150532 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z5") { double sigma = 0.0566743; return 0.343884 + rand_gen.Gaus(0, sigma); }

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