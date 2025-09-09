#include <iostream>
#include <string>

// HIPO
#include <hipo4/reader.h>
#include <hipo4/dictionary.h>
#include <hipo4/bank.h>

// Prefer the new RGAFiducialFilter; fall back to old FiducialFilter if needed
#if __has_include(<iguana/algorithms/clas12/RGAFiducialFilter/Algorithm.h>)
  #include <iguana/algorithms/clas12/RGAFiducialFilter/Algorithm.h>
  using FidAlg = iguana::clas12::RGAFiducialFilter;
#elif __has_include(<iguana/algorithms/clas12/FiducialFilter/Algorithm.h>)
  #include <iguana/algorithms/clas12/FiducialFilter/Algorithm.h>
  using FidAlg = iguana::clas12::FiducialFilter;
#else
  #error "No fiducial filter algorithm header found (RGAFiducialFilter or FiducialFilter)."
#endif

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " input.hipo [maxEvents]\n";
    return 1;
  }
  const std::string inFile = argv[1];
  long long maxEvents = (argc > 2) ? std::stoll(argv[2]) : -1;

  // Open HIPO & dictionary
  hipo::reader r;
  r.open(inFile.c_str());
  hipo::dictionary dict;
  r.readDictionary(dict);

  // Banks
  hipo::bank bPart   (dict.getSchema("REC::Particle"));
  hipo::bank bConfig (dict.getSchema("RUN::config"));
  hipo::bank bCal    (dict.hasSchema("REC::Calorimeter")   ? dict.getSchema("REC::Calorimeter")   : hipo::schema{});
  hipo::bank bFT     (dict.hasSchema("REC::ForwardTagger") ? dict.getSchema("REC::ForwardTagger") : hipo::schema{});
  hipo::bank bTraj   (dict.hasSchema("REC::Traj")          ? dict.getSchema("REC::Traj")          : hipo::schema{});

  // Algorithm
  FidAlg fid;

  // Start once with schemas
  hipo::banklist banks0;
  banks0.push_back(bPart);
  if (bCal.getSchema().getName()   != "") banks0.push_back(bCal);
  if (bFT.getSchema().getName()    != "") banks0.push_back(bFT);
  if (bTraj.getSchema().getName()  != "") banks0.push_back(bTraj);
  banks0.push_back(bConfig);
  fid.Start(banks0);

  long long evt = 0, keptTot = 0, cutTot = 0;
  hipo::event ev;

  while (r.next()) {
    if (maxEvents >= 0 && evt >= maxEvents) break;
    r.read(ev);

    // Fill banks for this event
    ev.getStructure(bPart);
    if (bCal.getSchema().getName()   != "") ev.getStructure(bCal);
    if (bFT.getSchema().getName()    != "") ev.getStructure(bFT);
    if (bTraj.getSchema().getName()  != "") ev.getStructure(bTraj);
    ev.getStructure(bConfig);

    // Event banklist
    hipo::banklist banks;
    banks.push_back(bPart);
    if (bCal.getSchema().getName()   != "") banks.push_back(bCal);
    if (bFT.getSchema().getName()    != "") banks.push_back(bFT);
    if (bTraj.getSchema().getName()  != "") banks.push_back(bTraj);
    banks.push_back(bConfig);

    const int before = (int)bPart.getRowList().size();
    fid.Run(banks);                    // must prune REC::Particle in-place
    const int after  = (int)bPart.getRowList().size();
    const int cut    = before - after;

    keptTot += after;
    cutTot  += cut;

    std::cout << "evt " << evt << ": kept " << after << " / " << before
              << " (cut " << cut << ")\n";
    ++evt;
  }

  fid.Stop();
  std::cout << "\nTOTAL over " << evt << " events: kept=" << keptTot
            << ", cut=" << cutTot << std::endl;
  return 0;
}