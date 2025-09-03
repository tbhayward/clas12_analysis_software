#include <TSystem.h>
#include <iostream>

// InclusiveKinematics header
#include <iguana/algorithms/physics/InclusiveKinematics/Algorithm.h>

void iguana_ex_cpp_ROOT_macro() {

  // If your include path is not already set by the env script, uncomment next line:
  // gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("IGUANA")));

  // Load the Iguana algorithms library
  gSystem->Load("libIguanaAlgorithms");

  // Start the algorithm
  iguana::physics::InclusiveKinematics algo;
  algo.Start();

  // Prepare event (e.g. run number 5032), then compute from a toy scattered lepton 4-vector
  auto key = algo.PrepareEvent(5032);
  auto result = algo.ComputeFromLepton(0.3, 0.3, 5.0, key); // px, py, pz [GeV]

  std::cout << "kinematics:"
            << "\n Q2 = " << result.Q2
            << "\n  x = " << result.x
            << "\n  W = " << result.W
            << std::endl;

  // Stop the algorithm
  algo.Stop();
}