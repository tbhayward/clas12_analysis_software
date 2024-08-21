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

  if (prefix == "Q2") { return 0.0588289+0.0990336*currentVariable+-0.0213364*std::pow(currentVariable,2); }
  if (prefix == "x") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
  if (prefix == "y") { return 0.0534692+0.916586*currentVariable+-1.73992*std::pow(currentVariable,2); }
  if (prefix == "z") { return 0.0119851+1.34778*currentVariable+-2.99383*std::pow(currentVariable,2); }
  if (prefix == "zeta") { return 0.849162+-2.95586*currentVariable+4.58654*std::pow(currentVariable,2); }
  if (prefix == "PT") { return 0.208234+-0.222063*currentVariable+0.652242*std::pow(currentVariable,2); }
  if (prefix == "xF") { return 0.206961+-0.0028204*currentVariable+-0.0428585*std::pow(currentVariable,2); }
  if (prefix == "Mx") { return 0.136851*exp(-0.5*std::pow((currentVariable - 0.12) / 0.235624, 2)) + 0.0280814*exp(-0.5*std::pow((currentVariable - 0.785) / -0.0544159, 2)) + 0.0970699 + 0.185273*currentVariable + -0.103331*std::pow(currentVariable, 2); }

  // if (prefix == "Q2") { return 0.042234+0.11043*currentVariable+-0.0225887*std::pow(currentVariable,2); }
  // if (prefix == "x") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
  // if (prefix == "y") { return 0.320777+-0.365915*currentVariable+0.246253*std::pow(currentVariable,2); }
  // if (prefix == "z") { return 0.0751743+0.693025*currentVariable+-1.06889*std::pow(currentVariable,2); }
  // if (prefix == "zeta") { return 1.2463+-5.31976*currentVariable+9.29149*std::pow(currentVariable,2); }
  // if (prefix == "PT") { return 0.216343+-0.273528*currentVariable+0.729059*std::pow(currentVariable,2); }
  // if (prefix == "xF") { return 0.218748+0.0086161*currentVariable+-0.0906104*std::pow(currentVariable,2); }

  if (prefix == "Q2all") { return 0.042234+0.11043*currentVariable+-0.0225887*std::pow(currentVariable,2); }
  if (prefix == "xall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
  if (prefix == "yall") { return 0.320777+-0.365915*currentVariable+0.246253*std::pow(currentVariable,2); }
  if (prefix == "zall") { return 0.0751743+0.693025*currentVariable+-1.06889*std::pow(currentVariable,2); }
  if (prefix == "zetaall") { return 1.2463+-5.31976*currentVariable+9.29149*std::pow(currentVariable,2); }
  if (prefix == "PTall") { return 0.216343+-0.273528*currentVariable+0.729059*std::pow(currentVariable,2); }
  if (prefix == "xFall") { return 0.218748+0.0086161*currentVariable+-0.0906104*std::pow(currentVariable,2); }


  if (prefix == "Q2x1z1") { double sigma = 0.00941327; return 0.149075 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x1z2") { double sigma = 0.0132642; return 0.162577 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x1z3") { double sigma = 0.0229451; return 0.141929 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x1z4") { double sigma = 0.0665136; return 0.181142 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x2z1") { double sigma = 0.00902454; return 0.141377 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x2z2") { double sigma = 0.00594617; return 0.15723 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x2z3") { double sigma = 0.0121427; return 0.164111 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x2z4") { double sigma = 0.00819397; return 0.184566 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x3z1") { double sigma = 0.0111792; return 0.159587 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x3z2") { double sigma = 0.00459304; return 0.180858 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x3z3") { double sigma = 0.0104746; return 0.17691 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x3z4") { double sigma = 0.0159035; return 0.160347 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x4z1") { double sigma = 0.00977397; return 0.175046 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x4z2") { double sigma = 0.00518098; return 0.186986 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x4z3") { double sigma = 0.0110679; return 0.182736 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x4z4") { double sigma = 0.0221351; return 0.176834 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x5z1") { double sigma = 0.146381; return 0.348494 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x5z2") { double sigma = 0.0108159; return 0.154464 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x5z3") { double sigma = 0.0123293; return 0.168112 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x5z4") { double sigma = 0.0113772; return 0.164353 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x6z1") { double sigma = 0.0219882; return 0.132161 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x6z2") { double sigma = 0.0106923; return 0.160582 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x6z3") { double sigma = 0.00607725; return 0.167237 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x6z4") { double sigma = 0.0110221; return 0.145863 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x7z1") { double sigma = 0.00527885; return 0.141097 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x7z2") { double sigma = 0.00544389; return 0.181025 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x7z3") { double sigma = 0.00507036; return 0.189089 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x7z4") { double sigma = 0.00954266; return 0.186328 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x8z1") { double sigma = 0.00882734; return 0.168493 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x8z2") { double sigma = 0.00430375; return 0.199761 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x8z3") { double sigma = 0.00440225; return 0.192723 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x8z4") { double sigma = 0.00854078; return 0.198004 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x9z2") { double sigma = 0.00738543; return 0.184486 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x9z3") { double sigma = 0.00710463; return 0.195825 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x9z4") { double sigma = 0.00505832; return 0.197091 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x10z1") { double sigma = 0.042577; return 0.0624643 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x10z2") { double sigma = 0.00775736; return 0.198887 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x10z3") { double sigma = 0.00395318; return 0.212409 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x10z4") { double sigma = 0.0023216; return 0.22681 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x11z1") { double sigma = 0.00856085; return 0.165176 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x11z2") { double sigma = 0.00353569; return 0.202314 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x11z3") { double sigma = 0.00314671; return 0.214343 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x11z4") { double sigma = 0.00615722; return 0.208119 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x12z2") { double sigma = 0.009042; return 0.207689 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x12z3") { double sigma = 0.00594662; return 0.223018 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x12z4") { double sigma = 0.00583763; return 0.242526 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x13z1") { double sigma = 0.00990338; return 0.174887 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x13z2") { double sigma = 0.00306666; return 0.218621 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x13z3") { double sigma = 0.00275945; return 0.240082 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x13z4") { double sigma = 0.00645786; return 0.241681 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x14z1") { double sigma = 0.0253505; return 0.115816 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x14z2") { double sigma = 0.0123555; return 0.246474 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x14z3") { double sigma = 0.00707888; return 0.266938 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2x14z4") { double sigma = 0.00480178; return 0.261037 + rand_gen.Gaus(0, sigma); }
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