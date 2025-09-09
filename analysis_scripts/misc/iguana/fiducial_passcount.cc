// Run with:
//   root -l -b -q 'rga_fiducial_passcount.C("input.hipo", 100000)'

#include <TSystem.h>
#include <iostream>

// HIPO
#include <hipo4/reader.h>
#include <hipo4/dictionary.h>
#include <hipo4/bank.h>

// The *new* algorithm:
#include <iguana/algorithms/clas12/RGAFiducialFilter/Algorithm.h>

void rga_fiducial_passcount(const char* inFile = "", Long64_t maxEvents = -1)
{
  if (!inFile || !*inFile) {
    std::cerr << "usage: rga_fiducial_passcount(\"/path/to/input.hipo\", [maxEvents])\n";
    return;
  }

  // Load Iguana libs
  gSystem->Load("libIguanaAlgorithms");

  // Open HIPO & dictionary
  hipo::reader r; r.open(inFile);
  hipo::dictionary dict; r.readDictionary(dict);

  // Banks
  hipo::bank bPart   (dict.getSchema("REC::Particle"));
  hipo::bank bConfig (dict.getSchema("RUN::config"));
  hipo::bank bCal    (dict.hasSchema("REC::Calorimeter")   ? dict.getSchema("REC::Calorimeter")   : hipo::schema{});
  hipo::bank bFT     (dict.hasSchema("REC::ForwardTagger") ? dict.getSchema("REC::ForwardTagger") : hipo::schema{});
  hipo::bank bTraj   (dict.hasSchema("REC::Traj")          ? dict.getSchema("REC::Traj")          : hipo::schema{});

  // Algorithm instance; uses installed YAML defaults automatically
  iguana::clas12::RGAFiducialFilter algo;

  // Start once with schemas
  hipo::banklist banks_start;
  banks_start.push_back(bPart);
  if (bCal.getSchema().getName()   != "") banks_start.push_back(bCal);
  if (bFT.getSchema().getName()    != "") banks_start.push_back(bFT);
  if (bTraj.getSchema().getName()  != "") banks_start.push_back(bTraj);
  banks_start.push_back(bConfig);
  algo.Start(banks_start);

  Long64_t evt = 0, keptTot = 0, cutTot = 0;
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

    // Build event banklist and run algorithm
    hipo::banklist banks_evt;
    banks_evt.push_back(bPart);
    if (bCal.getSchema().getName()   != "") banks_evt.push_back(bCal);
    if (bFT.getSchema().getName()    != "") banks_evt.push_back(bFT);
    if (bTraj.getSchema().getName()  != "") banks_evt.push_back(bTraj);
    banks_evt.push_back(bConfig);

    const int before = (int)bPart.getRowList().size();
    algo.Run(banks_evt);  // should prune REC::Particle in-place
    const int after  = (int)bPart.getRowList().size();
    const int cut    = before - after;

    keptTot += after; cutTot += cut;
    std::cout << "evt " << evt << ": kept " << after << " / " << before
              << " (cut " << cut << ")\n";
    ++evt;
  }

  algo.Stop();
  std::cout << "\nTOTAL over " << evt << " events: kept=" << keptTot
            << ", cut=" << cutTot << std::endl;
}