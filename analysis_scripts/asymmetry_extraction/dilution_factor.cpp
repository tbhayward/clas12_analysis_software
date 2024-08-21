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

  // if (prefix == "Q2") { return 0.0588289+0.0990336*currentVariable+-0.0213364*std::pow(currentVariable,2); }
  // if (prefix == "x") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
  // if (prefix == "y") { return 0.0534692+0.916586*currentVariable+-1.73992*std::pow(currentVariable,2); }
  // if (prefix == "z") { return 0.0119851+1.34778*currentVariable+-2.99383*std::pow(currentVariable,2); }
  // if (prefix == "zeta") { return 0.849162+-2.95586*currentVariable+4.58654*std::pow(currentVariable,2); }
  // if (prefix == "PT") { return 0.208234+-0.222063*currentVariable+0.652242*std::pow(currentVariable,2); }
  // if (prefix == "xF") { return 0.206961+-0.0028204*currentVariable+-0.0428585*std::pow(currentVariable,2); }
  // if (prefix == "Mx") { return 0.136851*exp(-0.5*std::pow((currentVariable - 0.12) / 0.235624, 2)) + 0.0280814*exp(-0.5*std::pow((currentVariable - 0.785) / -0.0544159, 2)) + 0.0970699 + 0.185273*currentVariable + -0.103331*std::pow(currentVariable, 2); }

  // if (prefix == "Q2all") { return 0.042234+0.11043*currentVariable+-0.0225887*std::pow(currentVariable,2); }
  // if (prefix == "xall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
  // if (prefix == "yall") { return 0.320777+-0.365915*currentVariable+0.246253*std::pow(currentVariable,2); }
  // if (prefix == "zall") { return 0.0751743+0.693025*currentVariable+-1.06889*std::pow(currentVariable,2); }
  // if (prefix == "zetaall") { return 1.2463+-5.31976*currentVariable+9.29149*std::pow(currentVariable,2); }
  // if (prefix == "PTall") { return 0.216343+-0.273528*currentVariable+0.729059*std::pow(currentVariable,2); }
  // if (prefix == "xFall") { return 0.218748+0.0086161*currentVariable+-0.0906104*std::pow(currentVariable,2); }

  if (prefix == "Q2") {
    double p0 = 0.0588289;
    double p1 = 0.0990336;
    double p2 = -0.0213364;
    double p3 = 0.00167213;
    // double p0_err = 0.0183319;
    // double p1_err = 0.0172302;
    // double p2_err = 0.00492875;
    // double p3_err = 0.000430101;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.0183319* 0.0183319 + 0.0172302* 0.0172302 * std::pow(currentVariable, 2) + 0.00492875* 0.00492875 * std::pow(currentVariable, 4) + 0.000430101* 0.000430101 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "x") {
    double p0 = 0.136288;
    double p1 = 0.283353;
    double p2 = -0.125745;
    double p3 = 0.212206;
    // double p0_err = 0.0133012;
    // double p1_err = 0.162195;
    // double p2_err = 0.607058;
    // double p3_err = 0.701473;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.0133012* 0.0133012 + 0.162195* 0.162195 * std::pow(currentVariable, 2) + 0.607058* 0.607058 * std::pow(currentVariable, 4) + 0.701473* 0.701473 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "y") {
    double p0 = 0.0534692;
    double p1 = 0.916586;
    double p2 = -1.73992;
    double p3 = 1.00709;
    // double p0_err = 0.113619;
    // double p1_err = 0.624366;
    // double p2_err = 1.12157;
    // double p3_err = 0.659314;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.113619* 0.113619 + 0.624366* 0.624366 * std::pow(currentVariable, 2) + 1.12157* 1.12157 * std::pow(currentVariable, 4) + 0.659314* 0.659314 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "z") {
    double p0 = 0.0119851;
    double p1 = 1.34778;
    double p2 = -2.99383;
    double p3 = 2.176;
    // double p0_err = 0.0157019;
    // double p1_err = 0.144574;
    // double p2_err = 0.411366;
    // double p3_err = 0.364181;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.0157019* 0.0157019 + 0.144574* 0.144574 * std::pow(currentVariable, 2) + 0.411366* 0.411366 * std::pow(currentVariable, 4) + 0.364181* 0.364181 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "zeta") {
    double p0 = 0.849162;
    double p1 = -2.95586;
    double p2 = 4.58654;
    double p3 = -2.5085;
    // double p0_err = 0.294276;
    // double p1_err = 1.65841;
    // double p2_err = 3.08362;
    // double p3_err = 1.89114;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.294276* 0.294276 + 1.65841* 1.65841 * std::pow(currentVariable, 2) + 3.08362* 3.08362 * std::pow(currentVariable, 4) + 1.89114* 1.89114 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "PT") {
    double p0 = 0.208234;
    double p1 = -0.222063;
    double p2 = 0.652242;
    double p3 = -0.478728;
    // double p0_err = 0.00607506;
    // double p1_err = 0.0468951;
    // double p2_err = 0.107988;
    // double p3_err = 0.0745258;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.00607506* 0.00607506 + 0.0468951* 0.0468951 * std::pow(currentVariable, 2) + 0.107988* 0.107988 * std::pow(currentVariable, 4) + 0.0745258* 0.0745258 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "xF") {
    double p0 = 0.206961;
    double p1 = -0.0028204;
    double p2 = -0.0428585;
    double p3 = 0.123092;
    // double p0_err = 0.0014299;
    // double p1_err = 0.00711853;
    // double p2_err = 0.0309251;
    // double p3_err = 0.0425113;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.0014299* 0.0014299 + 0.00711853* 0.00711853 * std::pow(currentVariable, 2) + 0.0309251* 0.0309251 * std::pow(currentVariable, 4) + 0.0425113* 0.0425113 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  
  if (prefix == "Mx") {
    double amp1 = 0.136851 + rand_gen.Gaus(0, 0.0205586);
    double mean1 = 0.12 + rand_gen.Gaus(0, 0.0296479);
    double sigma1 = 0.235624 + rand_gen.Gaus(0, 0.0126034);
    double amp2 = 0.0280814 + rand_gen.Gaus(0, 0.00494047);
    double mean2 = 0.785 + rand_gen.Gaus(0, 0.0237513);
    double sigma2 = -0.0544159 + rand_gen.Gaus(0, 0.0101243);
    double constTerm = 0.0970699 + rand_gen.Gaus(0, 0.0287095);
    double linearTerm = 0.185273 + rand_gen.Gaus(0, 0.0545798);
    double quadTerm = -0.103331 + rand_gen.Gaus(0, 0.0334958);
    double gauss1 = amp1 * exp(-0.5 * std::pow((currentVariable - mean1) / sigma1, 2));
    double gauss2 = amp2 * exp(-0.5 * std::pow((currentVariable - mean2) / sigma2, 2));
    double poly = constTerm + linearTerm * currentVariable + quadTerm * std::pow(currentVariable, 2);
    return gauss1 + gauss2 + poly;
  }

  if (prefix == "Q2all") {
    double p0 = 0.042234;
    double p1 = 0.11043;
    double p2 = -0.0225887;
    double p3 = 0.00159672;
    // double p0_err = 0.0166364;
    // double p1_err = 0.0154501;
    // double p2_err = 0.00437267;
    // double p3_err = 0.000377922;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.0166364* 0.0166364 + 0.0154501* 0.0154501 * std::pow(currentVariable, 2) + 0.00437267* 0.00437267 * std::pow(currentVariable, 4) + 0.000377922* 0.000377922 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "xall") {
    double p0 = 0.163228;
    double p1 = -0.0367291;
    double p2 = 0.804498;
    double p3 = -0.621569;
    // double p0_err = 0.0102424;
    // double p1_err = 0.116196;
    // double p2_err = 0.401017;
    // double p3_err = 0.425727;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.0102424* 0.0102424 + 0.116196* 0.116196 * std::pow(currentVariable, 2) + 0.401017* 0.401017 * std::pow(currentVariable, 4) + 0.425727* 0.425727 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "yall") {
    double p0 = 0.320777;
    double p1 = -0.365915;
    double p2 = 0.246253;
    double p3 = 0.00109763;
    // double p0_err = 0.0410817;
    // double p1_err = 0.243585;
    // double p2_err = 0.466392;
    // double p3_err = 0.28929;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.0410817* 0.0410817 + 0.243585* 0.243585 * std::pow(currentVariable, 2) + 0.466392* 0.466392 * std::pow(currentVariable, 4) + 0.28929* 0.28929 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "zall") {
    double p0 = 0.0751743;
    double p1 = 0.693025;
    double p2 = -1.06889;
    double p3 = 0.590062;
    // double p0_err = 0.0128945;
    // double p1_err = 0.111959;
    // double p2_err = 0.293561;
    // double p3_err = 0.234875;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.0128945* 0.0128945 + 0.111959* 0.111959 * std::pow(currentVariable, 2) + 0.293561* 0.293561 * std::pow(currentVariable, 4) + 0.234875* 0.234875 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "zetaall") {
    double p0 = 1.2463;
    double p1 = -5.31976;
    double p2 = 9.29149;
    double p3 = -5.56437;
    // double p0_err = 0.219781;
    // double p1_err = 1.22643;
    // double p2_err = 2.25747;
    // double p3_err = 1.37042;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.219781* 0.219781 + 1.22643* 1.22643 * std::pow(currentVariable, 2) + 2.25747* 2.25747 * std::pow(currentVariable, 4) + 1.37042* 1.37042 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "PTall") {
    double p0 = 0.216343;
    double p1 = -0.273528;
    double p2 = 0.729059;
    double p3 = -0.48757;
    // double p0_err = 0.00567176;
    // double p1_err = 0.0421257;
    // double p2_err = 0.0940832;
    // double p3_err = 0.0631752;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.00567176* 0.00567176 + 0.0421257* 0.0421257 * std::pow(currentVariable, 2) + 0.0940832* 0.0940832 * std::pow(currentVariable, 4) + 0.0631752* 0.0631752 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }
  if (prefix == "xFall") {
    double p0 = 0.218748;
    double p1 = 0.0086161;
    double p2 = -0.0906104;
    double p3 = 0.00809141;
    // double p0_err = 0.0011874;
    // double p1_err = 0.00500779;
    // double p2_err = 0.0192647;
    // double p3_err = 0.0234072;
    double central_value = p0 + p1 * currentVariable + 
                           p2 * std::pow(currentVariable, 2) + 
                           p3 * std::pow(currentVariable, 3);
    double uncertainty = std::sqrt(0.0011874* 0.0011874 + 0.00500779* 0.00500779 * std::pow(currentVariable, 2) + 0.0192647* 0.0192647 * std::pow(currentVariable, 4) + 0.0234072* 0.0234072 * std::pow(currentVariable, 6));
    double result = rand_gen.Gaus(central_value, uncertainty);
    return result;
  }

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