#include <hipo4/reader.h>
#include <hipo4/event.h>
#include <hipo4/bank.h>
#include <hipo4/dictionary.h>

#include <iguana/algorithms/clas12/RGAFiducialFilter/Algorithm.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>

#include <unordered_set>
#include <string>

int main(int argc, char** argv)
{
  const char* in_file         = (argc > 1 ? argv[1] : "data.hipo");
  const int   num_events_req  = (argc > 2 ? std::stoi(argv[2]) : 1000);
  const bool  interactive     = (argc > 3 ? std::string(argv[3]) == "true" : false);
  TApplication* app = interactive ? new TApplication("app", &argc, argv) : nullptr;

  hipo::reader     reader;
  hipo::event      event;
  hipo::dictionary dict;
  reader.open(in_file);
  reader.readDictionary(dict);

  // Build ONE banklist and keep using those same objects
  hipo::banklist banks;
  const int idx_particle = 0, idx_calor = 1, idx_config = 2;
  banks.emplace_back(dict.getSchema("REC::Particle"));
  banks.emplace_back(dict.getSchema("REC::Calorimeter"));
  banks.emplace_back(dict.getSchema("RUN::config"));

  // Start algorithm with that banklist
  iguana::clas12::RGAFiducialFilter fid;
  fid.Start(banks);

  // Histograms: PCAL lv before/after (0..400)
  TH1D* h_lv_before = new TH1D("h_lv_before", "PCAL lv;lv;counts", 400, 0, 400);
  TH1D* h_lv_after  = new TH1D("h_lv_after" , "PCAL lv;lv;counts", 400, 0, 400);
  h_lv_before->SetLineWidth(2);
  h_lv_after ->SetLineWidth(2);
  h_lv_after ->SetLineColor(kRed);

  int ev = 0;
  while (reader.next() && (num_events_req == 0 || ev < num_events_req)) {
    reader.read(event);

    // Read event data INTO the banks INSIDE the banklist
    event.read(banks.at(idx_particle));
    event.read(banks.at(idx_calor));
    event.read(banks.at(idx_config));

    auto& recPart = banks.at(idx_particle);
    auto& recCal  = banks.at(idx_calor);

    // BEFORE: all PCAL hits (layer==1)
    const int ncal = recCal.getRows();
    for (int i = 0; i < ncal; ++i) {
      if (recCal.getInt("layer", i) == 1) {
        h_lv_before->Fill(recCal.getFloat("lv", i));
      }
    }

    // Run filter (mutates REC::Particle inside banklist)
    fid.Run(banks);

    // Collect kept track indices after filtering
    std::unordered_set<int> kept;
    for (auto const& row : recPart.getRowList()) kept.insert(row);

    // AFTER: only PCAL hits whose pindex survived
    for (int i = 0; i < ncal; ++i) {
      if (recCal.getInt("layer", i) == 1) {
        if (kept.count(recCal.getInt("pindex", i)))
          h_lv_after->Fill(recCal.getFloat("lv", i));
      }
    }

    ++ev;
  }

  fid.Stop();

  TCanvas* c = new TCanvas("c", "RGAFiducialFilter PCAL lv (before vs after)", 900, 700);
  h_lv_before->Draw("hist");
  h_lv_after->Draw("hist same");
  auto leg = new TLegend(0.60, 0.75, 0.88, 0.88);
  leg->AddEntry(h_lv_before, "before", "l");
  leg->AddEntry(h_lv_after , "after" , "l");
  leg->Draw();

  if (interactive) {
    app->Run();
  } else {
    c->SaveAs("out-rga-fiducial-example.png");
  }
  return 0;
}