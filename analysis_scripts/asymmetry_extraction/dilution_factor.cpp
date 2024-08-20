#include "dilution_factor.h"
#include <cmath>
#include <iostream>
#include <TRandom3.h>

double dilution_factor(double currentVariable, const std::string& prefix) {

  TRandom3 rand_gen;

  // epi+X
  if (prefix == "epipPT") {
    return 0.184542-0.0499585*currentVariable+0.163844*std::pow(currentVariable,2)+
      0.157106*std::pow(currentVariable,3);
  }

  if (prefix == "runnum") { return 0.1949;}
  if (prefix == "Q2") { return 0.0662175+0.0912281*currentVariable+-0.0188005*std::pow(currentVariable,2); }
  if (prefix == "x") { return 0.136311+0.283958*currentVariable+-0.125658*std::pow(currentVariable,2); }
  if (prefix == "y") { return 0.0541015+0.913993*currentVariable+-1.73469*std::pow(currentVariable,2); }
  if (prefix == "z") { return 0.0116813+1.35182*currentVariable+-3.00502*std::pow(currentVariable,2); }
  if (prefix == "zeta") { return 0.817679+-2.78097*currentVariable+4.26881*std::pow(currentVariable,2); }
  if (prefix == "PT") { return 0.208315+-0.221463*currentVariable+0.651064*std::pow(currentVariable,2); }
  if (prefix == "xF") { return 0.207138+-0.00265586*currentVariable+-0.0427432*std::pow(currentVariable,2); }
  if (prefix == "Mx") { return 0.135672*exp(-0.5*std::pow((currentVariable - 0.12) / 0.234466, 2)) + 0.0278955*exp(-0.5*std::pow((currentVariable - 0.785) / -0.0544874, 2)) + 0.0978352 + 0.186203*currentVariable + -0.104758*std::pow(currentVariable, 2); }
  
  if (prefix == "Q2all") { return 0.0515679+0.100782*currentVariable+-0.019643*std::pow(currentVariable,2); }
  if (prefix == "xall") { return 0.163076+-0.0361384*currentVariable+0.804483*std::pow(currentVariable,2); }
  if (prefix == "yall") { return 0.310092+-0.304066*currentVariable+0.135163*std::pow(currentVariable,2); }
  if (prefix == "zall") { return 0.0757085+0.687572*currentVariable+-1.05333*std::pow(currentVariable,2); }
  if (prefix == "zetaall") { return 1.26316+-5.40117*currentVariable+9.42118*std::pow(currentVariable,2); }
  if (prefix == "PTall") { return 0.216951+-0.274862*currentVariable+0.730387*std::pow(currentVariable,2); }
  if (prefix == "xFall") { return 0.219074+0.00890221*currentVariable+-0.0921533*std::pow(currentVariable,2); }

  if (prefix == "Q2y4z1") { double sigma = 0.0274445; return 0.184089 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z2") { double sigma = 0.00902867; return 0.163258 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z3") { double sigma = 0.0071842; return 0.184709 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z4") { double sigma = 0.0067977; return 0.182006 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z5") { double sigma = 0.00780683; return 0.185804 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z1") { double sigma = 0.0337462; return 0.177533 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z2") { double sigma = 0.0075944; return 0.196864 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z3") { double sigma = 0.00344734; return 0.225978 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z4") { double sigma = 0.00720619; return 0.234865 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z5") { double sigma = 0.00581591; return 0.247905 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z2") { double sigma = 0.00947351; return 0.246897 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z3") { double sigma = 0.00767194; return 0.256558 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z4") { double sigma = 0.00895384; return 0.246275 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z5") { double sigma = 0.0193475; return 0.250308 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z1") { double sigma = 0.00645845; return 0.14725 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z2") { double sigma = 0.00626156; return 0.177352 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z3") { double sigma = 0.00554662; return 0.179113 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z4") { double sigma = 0.00872289; return 0.164078 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z5") { double sigma = 0.0191617; return 0.192845 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z1") { double sigma = 0.00785352; return 0.161456 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z2") { double sigma = 0.00508948; return 0.21293 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z3") { double sigma = 0.00417176; return 0.215924 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z4") { double sigma = 0.00444297; return 0.215889 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z5") { double sigma = 0.0129412; return 0.199264 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z1") { double sigma = 0.00869239; return 0.18159 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z2") { double sigma = 0.00244425; return 0.232978 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z3") { double sigma = 0.00719626; return 0.235411 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z4") { double sigma = 0.00810567; return 0.239474 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z5") { double sigma = 0.00775354; return 0.252125 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z1") { double sigma = 0.00553508; return 0.256934 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z2") { double sigma = 0.00481419; return 0.243236 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z3") { double sigma = 0.0111252; return 0.25813 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z4") { double sigma = 0.00675023; return 0.256196 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z5") { double sigma = 0.0264374; return 0.275725 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z1") { double sigma = 0.00804756; return 0.158207 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z2") { double sigma = 0.00307527; return 0.187336 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z3") { double sigma = 0.00675439; return 0.169752 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z4") { double sigma = 0.00989174; return 0.168636 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z5") { double sigma = 0.0388246; return 0.218049 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z1") { double sigma = 0.00709041; return 0.182044 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z2") { double sigma = 0.00378046; return 0.208185 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z3") { double sigma = 0.0037177; return 0.20001 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z4") { double sigma = 0.00642133; return 0.206926 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z5") { double sigma = 0.0245985; return 0.172981 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z1") { double sigma = 0.0048396; return 0.200179 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z2") { double sigma = 0.00303875; return 0.224349 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z3") { double sigma = 0.00481792; return 0.224789 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z4") { double sigma = 0.0101067; return 0.198917 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z5") { double sigma = 0.0437512; return 0.240294 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z1") { double sigma = 0.0111535; return 0.192329 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z2") { double sigma = 0.0049102; return 0.219994 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z3") { double sigma = 0.00887562; return 0.239783 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z4") { double sigma = 0.0149301; return 0.260102 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z5") { double sigma = 0.038012; return 0.255697 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z1") { double sigma = 0.0187904; return 0.212049 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z2") { double sigma = 0.00773389; return 0.251115 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z3") { double sigma = 0.00782416; return 0.273362 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z4") { double sigma = 0.00915181; return 0.281845 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z5") { double sigma = 0.0346296; return 0.344907 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z1") { double sigma = 0.0062149; return 0.172189 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z2") { double sigma = 0.00459529; return 0.187632 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z3") { double sigma = 0.0101788; return 0.177895 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z4") { double sigma = 0.0296317; return 0.217736 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z5") { double sigma = 0.0432561; return 0.218876 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z1") { double sigma = 0.00719994; return 0.186163 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z2") { double sigma = 0.00449874; return 0.202302 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z3") { double sigma = 0.0071413; return 0.181966 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z4") { double sigma = 0.0281577; return 0.208254 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z5") { double sigma = 0.0353484; return 0.26305 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z1") { double sigma = 0.00715909; return 0.190122 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z2") { double sigma = 0.00801141; return 0.203101 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z3") { double sigma = 0.0085226; return 0.188999 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z4") { double sigma = 0.0296864; return 0.256718 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z5") { double sigma = 0.0342572; return 0.328478 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z1") { double sigma = 0.0074808; return 0.204487 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z2") { double sigma = 0.00651027; return 0.238415 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z3") { double sigma = 0.0147262; return 0.221938 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z4") { double sigma = 0.0225536; return 0.233242 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z5") { double sigma = 0.0363611; return 0.409126 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z1") { double sigma = 0.0108632; return 0.212234 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z2") { double sigma = 0.00515751; return 0.25681 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z3") { double sigma = 0.0132167; return 0.247501 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z4") { double sigma = 0.0376663; return 0.161071 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z5") { double sigma = 0.0520532; return 0.333858 + rand_gen.Gaus(0, sigma); }
  // // epX
  // if (prefix == "xF") { 
  //   return 0.186121-0.0263337*currentVariable-0.175587*std::pow(currentVariable,2)+
  //     0.0522814*std::pow(currentVariable,3);
  // }
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
  // // epX
  // if (prefix == "Mx") { 
  //   return 0.0847657+0.0762168*currentVariable-0.0128988*std::pow(currentVariable,2)+
  //     0.00274429*std::pow(currentVariable,3);
  // }
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