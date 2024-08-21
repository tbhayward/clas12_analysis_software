#include "dilution_factor.h"
#include <cmath>
#include <iostream>
#include <TRandom3.h>

double dilution_factor(double currentVariable, const std::string& prefix) {

  TRandom3 rand_gen;

  if (prefix == "Q2") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
  if (prefix == "x") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
  if (prefix == "y") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
  if (prefix == "z") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
  if (prefix == "zeta") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
  if (prefix == "PT") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
  if (prefix == "xF") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
  if (prefix == "Mx") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }

  if (prefix == "Q2all") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
  if (prefix == "xall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
  if (prefix == "yall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
  if (prefix == "zall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
  if (prefix == "zetaall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
  if (prefix == "PTall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
  if (prefix == "xFall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }


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
  
  return 0.19;
}