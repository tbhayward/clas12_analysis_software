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

  if (prefix == "runnum") { return 0.1886;}
  if (prefix == "Q2") { return 0.0586566+0.0924459*currentVariable+-0.0191629*std::pow(currentVariable,2); }
  if (prefix == "x") { return 0.131499+0.274516*currentVariable+-0.119065*std::pow(currentVariable,2); }
  if (prefix == "y") { return 0.0274489+1.00838*currentVariable+-1.88842*std::pow(currentVariable,2); }
  if (prefix == "z") { return 0.00670293+1.35226*currentVariable+-3.03035*std::pow(currentVariable,2); }
  if (prefix == "zeta") { return 0.778565+-2.58993*currentVariable+3.89568*std::pow(currentVariable,2); }
  if (prefix == "PT") { return 0.201133+-0.220436*currentVariable+0.656143*std::pow(currentVariable,2); }
  if (prefix == "xF") { return 0.199983+-0.0058117*currentVariable+-0.0425802*std::pow(currentVariable,2); }
  if (prefix == "Mx") { return 0.139383*exp(-0.5*std::pow((currentVariable - 0.120003) / 0.23574, 2)) + 0.0289296*exp(-0.5*std::pow((currentVariable - 0.785) / -0.0542062, 2)) + 0.0892425 + 0.189842*currentVariable + -0.106771*std::pow(currentVariable, 2); }

  if (prefix == "Q2") { return 0.069162+0.0852323*currentVariable+-0.016417*std::pow(currentVariable,2); }
  if (prefix == "xall") { return 0.199592+-0.378803*currentVariable+1.66622*std::pow(currentVariable,2); }
  if (prefix == "yall") { return 0.25066+-0.0182866*currentVariable+-0.42064*std::pow(currentVariable,2); }
  if (prefix == "zall") { return 0.113339+0.396386*currentVariable+-0.458942*std::pow(currentVariable,2); }
  if (prefix == "zetaall") { return 1.17109+-4.89527*currentVariable+8.4361*std::pow(currentVariable,2); }
  if (prefix == "PTall") { return 0.220685+-0.38361*currentVariable+0.96351*std::pow(currentVariable,2); }
  if (prefix == "xFall") { return 0.213126+0.0084475*currentVariable+-0.110573*std::pow(currentVariable,2); }

  if (prefix == "Q2exclusive") { return 0.0785674+0.082314*currentVariable+-0.0156957*std::pow(currentVariable,2); }
  if (prefix == "xexclusive") { return 0.353541+-1.88346*currentVariable+5.97796*std::pow(currentVariable,2); }
  if (prefix == "yexclusive") { return 0.0605892+1.32055*currentVariable+-3.43166*std::pow(currentVariable,2); }
  if (prefix == "zexclusive") { return 0.179395+-0.196159*currentVariable+1.13571*std::pow(currentVariable,2); }
  if (prefix == "zetaexclusive") { return 0.20643+0.755261*currentVariable+-1.76453*std::pow(currentVariable,2); }
  if (prefix == "PTexclusive") { return 0.25699+-0.69242*currentVariable+1.62295*std::pow(currentVariable,2); }
  if (prefix == "xFexclusive") { return 0.245947+-0.00360302*currentVariable+-0.259858*std::pow(currentVariable,2); }

  if (prefix == "Q2y4z1") { double sigma = 0.0300989; return 0.183729 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z2") { double sigma = 0.00959404; return 0.155779 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z3") { double sigma = 0.00739077; return 0.176862 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z4") { double sigma = 0.00682892; return 0.173158 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z5") { double sigma = 0.00813325; return 0.177581 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z1") { double sigma = 0.0323061; return 0.196968 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z2") { double sigma = 0.00759313; return 0.188365 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z3") { double sigma = 0.00346376; return 0.218019 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z4") { double sigma = 0.00726652; return 0.227776 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z5") { double sigma = 0.006145; return 0.239521 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z2") { double sigma = 0.0101037; return 0.242016 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z3") { double sigma = 0.00766527; return 0.249692 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z4") { double sigma = 0.00904429; return 0.238759 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z5") { double sigma = 0.0198533; return 0.242168 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z1") { double sigma = 0.00634032; return 0.140333 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z2") { double sigma = 0.00645917; return 0.171227 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z3") { double sigma = 0.00586263; return 0.172025 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z4") { double sigma = 0.00860682; return 0.155746 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z5") { double sigma = 0.018945; return 0.183599 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z1") { double sigma = 0.0079919; return 0.154717 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z2") { double sigma = 0.00535261; return 0.206773 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z3") { double sigma = 0.0044391; return 0.208906 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z4") { double sigma = 0.00422322; return 0.208627 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z5") { double sigma = 0.0130988; return 0.190384 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z1") { double sigma = 0.00889463; return 0.173019 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z2") { double sigma = 0.0029398; return 0.227325 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z3") { double sigma = 0.00701893; return 0.228699 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z4") { double sigma = 0.00767197; return 0.233023 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z5") { double sigma = 0.00791856; return 0.245471 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z1") { double sigma = 0.00915917; return 0.244723 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z2") { double sigma = 0.00501937; return 0.237058 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z3") { double sigma = 0.0116841; return 0.250107 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z4") { double sigma = 0.00650333; return 0.247137 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z5") { double sigma = 0.0269787; return 0.266435 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z1") { double sigma = 0.00808594; return 0.152169 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z2") { double sigma = 0.00316616; return 0.180812 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z3") { double sigma = 0.00703419; return 0.161469 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z4") { double sigma = 0.0102655; return 0.160687 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z5") { double sigma = 0.0388156; return 0.205713 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z1") { double sigma = 0.00716223; return 0.175822 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z2") { double sigma = 0.00393985; return 0.202306 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z3") { double sigma = 0.00365109; return 0.192781 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z4") { double sigma = 0.00678465; return 0.199473 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z5") { double sigma = 0.0247816; return 0.167823 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z1") { double sigma = 0.00471587; return 0.195018 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z2") { double sigma = 0.00329406; return 0.217805 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z3") { double sigma = 0.00518693; return 0.217901 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z4") { double sigma = 0.0106822; return 0.19078 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z5") { double sigma = 0.0441512; return 0.231607 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z1") { double sigma = 0.0118995; return 0.185545 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z2") { double sigma = 0.00471725; return 0.213421 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z3") { double sigma = 0.00895872; return 0.232042 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z4") { double sigma = 0.0154667; return 0.253428 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z5") { double sigma = 0.0396196; return 0.249491 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z1") { double sigma = 0.0180113; return 0.209211 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z2") { double sigma = 0.00786166; return 0.243704 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z3") { double sigma = 0.00879327; return 0.268846 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z4") { double sigma = 0.00873356; return 0.274875 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z5") { double sigma = 0.0363972; return 0.339988 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z1") { double sigma = 0.00630227; return 0.166736 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z2") { double sigma = 0.00477693; return 0.181605 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z3") { double sigma = 0.0102332; return 0.169801 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z4") { double sigma = 0.031231; return 0.210979 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z5") { double sigma = 0.0427374; return 0.215477 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z1") { double sigma = 0.00734848; return 0.180773 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z2") { double sigma = 0.00444057; return 0.196631 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z3") { double sigma = 0.00742509; return 0.175581 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z4") { double sigma = 0.0290777; return 0.204038 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z5") { double sigma = 0.0356724; return 0.263743 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z1") { double sigma = 0.00722627; return 0.184244 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z2") { double sigma = 0.00831526; return 0.196523 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z3") { double sigma = 0.00896459; return 0.182607 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z4") { double sigma = 0.0308985; return 0.248755 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z5") { double sigma = 0.0378749; return 0.326304 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z1") { double sigma = 0.00769362; return 0.197648 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z2") { double sigma = 0.0067651; return 0.232681 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z3") { double sigma = 0.0149076; return 0.214735 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z4") { double sigma = 0.0227755; return 0.227448 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z5") { double sigma = 0.0350961; return 0.407584 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z1") { double sigma = 0.0115227; return 0.207518 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z2") { double sigma = 0.00523784; return 0.251222 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z3") { double sigma = 0.0143091; return 0.239871 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z4") { double sigma = 0.0372947; return 0.151176 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z5") { double sigma = 0.0599256; return 0.335865 + rand_gen.Gaus(0, sigma); }
  
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