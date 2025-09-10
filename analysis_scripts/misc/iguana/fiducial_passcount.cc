#include <iostream>
#include <string>

#include <hipo4/reader.h>
#include <hipo4/bank.h>
#include <hipo4/dictionary.h>

#include <iguana/algorithms/clas12/RGAFiducialFilter/Algorithm.h>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " input.hipo [maxEvents]\n";
    return 1;
  }

  const std::string inFile = argv[1];
  long long maxEvents = (argc > 2) ? std::stoll(argv[2]) : -1;

  // Open file and read dictionary
  hipo::reader r;
  r.open(inFile.c_str());
  hipo::dictionary dict;
  r.readDictionary(dict);

  // Banks (schemas)
  hipo::bank bPart   (dict.getSchema("REC::Particle"));
  hipo::bank bCal    (dict.hasSchema("REC::Calorimeter")   ? dict.getSchema("REC::Calorimeter")   : hipo::schema{});
  hipo::bank bFT     (dict.hasSchema("REC::ForwardTagger") ? dict.getSchema("REC::ForwardTagger") : hipo::schema{});
  hipo::bank bTraj   (dict.hasSchema("REC::Traj")          ? dict.getSchema("REC::Traj")          : hipo::schema{});
  hipo::bank bConfig (dict.getSchema("RUN::config"));

  // Algorithm instance (auto-loads installed Config.yaml via IGUANA_ETCDIR)
  iguana::clas12::RGAFiducialFilter fid;
  fid.SetLogLevel("debug"); // change to "info" or remove when you're done checking

  // Start once (with schema-only banks)
  {
    hipo::banklist banks0;
    banks0.push_back(bPart);
    if (bCal.getSchema().getName()   != "") banks0.push_back(bCal);
    if (bFT.getSchema().getName()    != "") banks0.push_back(bFT);
    if (bTraj.getSchema().getName()  != "") banks0.push_back(bTraj);
    banks0.push_back(bConfig);
    fid.Start(banks0);
  }

  long long evt = 0, keptTot = 0, cutTot = 0;

  // Event loop
  hipo::event ev;
  while (r.next()) {
    if (maxEvents >= 0 && evt >= maxEvents) break;

    // Read event data into the standalone bank objects
    r.read(ev);
    ev.getStructure(bPart);
    if (bCal.getSchema().getName()   != "") ev.getStructure(bCal);
    if (bFT.getSchema().getName()    != "") ev.getStructure(bFT);
    if (bTraj.getSchema().getName()  != "") ev.getStructure(bTraj);
    ev.getStructure(bConfig);

    // Build the banklist for THIS event (copies the banks as they are now)
    hipo::banklist banks;
    banks.push_back(bPart);
    if (bCal.getSchema().getName()   != "") banks.push_back(bCal);
    if (bFT.getSchema().getName()    != "") banks.push_back(bFT);
    if (bTraj.getSchema().getName()  != "") banks.push_back(bTraj);
    banks.push_back(bConfig);

    // IMPORTANT: measure before/after sizes from the bank **inside** `banks`
    const int before = static_cast<int>(banks.at(0).getRowList().size());
    fid.Run(banks); // prunes REC::Particle (banks.at(0)) in-place
    const int after  = static_cast<int>(banks.at(0).getRowList().size());
    const int cut    = before - after;

    keptTot += after;
    cutTot  += cut;

    std::cout << "evt " << evt << ": kept " << after << " / " << before
              << " (cut " << cut << ")\n";
    ++evt;
  }

  std::cout << "\nTOTAL: kept " << keptTot
            << " | cut " << cutTot
            << " | events " << evt << std::endl;

  return 0;
}