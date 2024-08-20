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

  if (prefix == "Q2all") { return 0.0756166+0.084845*currentVariable+-0.0163104*std::pow(currentVariable,2); }
  if (prefix == "xall") { return 0.203144+-0.362312*currentVariable+1.64139*std::pow(currentVariable,2); }
  if (prefix == "yall") { return 0.26201+-0.0383504*currentVariable+-0.390395*std::pow(currentVariable,2); }
  if (prefix == "zall") { return 0.116095+0.413281*currentVariable+-0.479925*std::pow(currentVariable,2); }
  if (prefix == "zetaall") { return 1.21156+-5.09334*currentVariable+8.81791*std::pow(currentVariable,2); }
  if (prefix == "PTall") { return 0.228164+-0.386966*currentVariable+0.963299*std::pow(currentVariable,2); }
  if (prefix == "xFall") { return 0.220114+0.0117653*currentVariable+-0.107717*std::pow(currentVariable,2); }

  if (prefix == "Q2exclusive") { return 0.0785674+0.082314*currentVariable+-0.0156957*std::pow(currentVariable,2); }
  if (prefix == "xexclusive") { return 0.353541+-1.88346*currentVariable+5.97796*std::pow(currentVariable,2); }
  if (prefix == "yexclusive") { return 0.0605892+1.32055*currentVariable+-3.43166*std::pow(currentVariable,2); }
  if (prefix == "zexclusive") { return 0.179395+-0.196159*currentVariable+1.13571*std::pow(currentVariable,2); }
  if (prefix == "zetaexclusive") { return 0.20643+0.755261*currentVariable+-1.76453*std::pow(currentVariable,2); }
  if (prefix == "PTexclusive") { return 0.25699+-0.69242*currentVariable+1.62295*std::pow(currentVariable,2); }
  if (prefix == "xFexclusive") { return 0.245947+-0.00360302*currentVariable+-0.259858*std::pow(currentVariable,2); }

  if (prefix == "Q2y4z1") { double sigma = 0.0300909; return 0.183912 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z2") { double sigma = 0.00959726; return 0.155972 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z3") { double sigma = 0.00738975; return 0.177047 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z4") { double sigma = 0.00682639; return 0.173348 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y4z5") { double sigma = 0.00813086; return 0.177757 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z1") { double sigma = 0.0323008; return 0.197157 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z2") { double sigma = 0.00758883; return 0.188543 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z3") { double sigma = 0.00346358; return 0.218191 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z4") { double sigma = 0.00726378; return 0.227942 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y8z5") { double sigma = 0.00614267; return 0.239683 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z1") { double sigma = -nan; return 0 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z2") { double sigma = 0.0100995; return 0.242177 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z3") { double sigma = 0.00766089; return 0.24985 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z4") { double sigma = 0.0090415; return 0.238917 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y12z5") { double sigma = 0.0198458; return 0.24235 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z1") { double sigma = 0.00633819; return 0.140529 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z2") { double sigma = 0.00645702; return 0.171414 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z3") { double sigma = 0.00586112; return 0.172212 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z4") { double sigma = 0.00860361; return 0.155939 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y3z5") { double sigma = 0.0189378; return 0.183772 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z1") { double sigma = 0.00798994; return 0.15491 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z2") { double sigma = 0.0053508; return 0.206948 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z3") { double sigma = 0.00443767; return 0.209081 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z4") { double sigma = 0.00422135; return 0.2088 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y7z5") { double sigma = 0.0130951; return 0.190561 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z1") { double sigma = 0.00888966; return 0.173208 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z2") { double sigma = 0.00293896; return 0.227492 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z3") { double sigma = 0.00701601; return 0.228866 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z4") { double sigma = 0.00767819; return 0.233188 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y11z5") { double sigma = 0.00792404; return 0.245633 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z1") { double sigma = 0.00915535; return 0.244837 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z2") { double sigma = 0.00501852; return 0.237227 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z3") { double sigma = 0.0116817; return 0.250268 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z4") { double sigma = 0.0064994; return 0.247303 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y15z5") { double sigma = 0.0269718; return 0.266595 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z1") { double sigma = 0.00808301; return 0.152362 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z2") { double sigma = 0.0031651; return 0.180996 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z3") { double sigma = 0.00703154; return 0.16166 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z4") { double sigma = 0.0102618; return 0.16088 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y2z5") { double sigma = 0.0388028; return 0.205882 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z1") { double sigma = 0.0071602; return 0.176008 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z2") { double sigma = 0.00393917; return 0.202484 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z3") { double sigma = 0.00365014; return 0.19296 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z4") { double sigma = 0.00678369; return 0.199642 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y6z5") { double sigma = 0.0247716; return 0.168021 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z1") { double sigma = 0.00471479; return 0.195198 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z2") { double sigma = 0.00329318; return 0.217976 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z3") { double sigma = 0.00518471; return 0.218075 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z4") { double sigma = 0.0106784; return 0.190958 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y10z5") { double sigma = 0.0441318; return 0.231748 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z1") { double sigma = 0.0118966; return 0.185726 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z2") { double sigma = 0.00471605; return 0.213598 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z3") { double sigma = 0.00895934; return 0.232214 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z4") { double sigma = 0.0154583; return 0.253573 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y14z5") { double sigma = 0.0396214; return 0.249673 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z1") { double sigma = 0.0180062; return 0.209384 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z2") { double sigma = 0.00785739; return 0.243871 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z3") { double sigma = 0.00879313; return 0.269003 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z4") { double sigma = 0.00872951; return 0.275013 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y17z5") { double sigma = 0.0363961; return 0.340081 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z1") { double sigma = 0.00630035; return 0.166926 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z2") { double sigma = 0.00477602; return 0.181787 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z3") { double sigma = 0.0102295; return 0.169987 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z4") { double sigma = 0.0312266; return 0.211149 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y1z5") { double sigma = 0.0427086; return 0.215615 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z1") { double sigma = 0.0073461; return 0.180957 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z2") { double sigma = 0.00443885; return 0.196809 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z3") { double sigma = 0.00742233; return 0.175765 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z4") { double sigma = 0.0290689; return 0.204221 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y5z5") { double sigma = 0.0356579; return 0.26392 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z1") { double sigma = 0.00722362; return 0.184427 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z2") { double sigma = 0.00831347; return 0.196703 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z3") { double sigma = 0.00896256; return 0.182785 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z4") { double sigma = 0.030899; return 0.248921 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y9z5") { double sigma = 0.0378476; return 0.326427 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z1") { double sigma = 0.00769142; return 0.197828 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z2") { double sigma = 0.00676175; return 0.232846 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z3") { double sigma = 0.0149048; return 0.21491 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z4") { double sigma = 0.0227644; return 0.227619 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y13z5") { double sigma = 0.0350796; return 0.407768 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z1") { double sigma = 0.0115209; return 0.207691 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z2") { double sigma = 0.00523744; return 0.251386 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z3") { double sigma = 0.0143049; return 0.240038 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z4") { double sigma = 0.037285; return 0.151392 + rand_gen.Gaus(0, sigma); }

  if (prefix == "Q2y16z5") { double sigma = 0.0599228; return 0.336022 + rand_gen.Gaus(0, sigma); }
  
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