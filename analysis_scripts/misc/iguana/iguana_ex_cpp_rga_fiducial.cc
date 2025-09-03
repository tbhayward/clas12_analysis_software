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
#include <iostream>

int main(int argc, char** argv)
{
  // -------- args --------
  const char* in_file         = (argc > 1 ? argv[1] : "data.hipo");
  const int   num_events_req  = (argc > 2 ? std::stoi(argv[2]) : 1000);
  const bool  interactive     = (argc > 3 ? std::string(argv[3]) == "true" : false);

  // -------- ROOT app if interactive --------
  TApplication* app = interactive ? new TApplication("app", &argc, argv) : nullptr;

  // -------- open HIPO, dictionary, banks --------
  hipo::reader      reader;
  hipo::event       event;
  hipo::dictionary  dict;

  reader.open(in_file);
  reader.readDictionary(dict);

  hipo::bank recPart (dict.getSchema("REC::Particle"));
  hipo::bank recCal  (dict.getSchema("REC::Calorimeter"));
  hipo::bank runConf (dict.getSchema("RUN::config"));

  // -------- build a banklist for Iguana --------
  hipo::banklist banks;
  banks.add(recPart);
  banks.add(recCal);
  banks.add(runConf);

  // -------- start algorithm --------
  iguana::clas12::RGAFiducialFilter fid;
  fid.Start(banks);

  // -------- histograms (before/after PCAL lv, 0..400) --------
  TH1D* h_lv_before = new TH1D("h_lv_before", "PCAL lv;lv;counts", 400, 0, 400);
  TH1D* h_lv_after  = new TH1D("h_lv_after" , "PCAL lv;lv;counts", 400, 0, 400);
  h_lv_before->SetLineWidth(2);
  h_lv_after ->SetLineWidth(2);
  h_lv_after ->SetLineColor(kRed);

  // -------- event loop --------
  int ev = 0;
  while (reader.read(event) && (num_events_req == 0 || ev < num_events_req)) {
    // read banks for this event
    event.read(recPart);
    event.read(recCal);
    event.read(runConf);

    // fill "before" (all PCAL hits: layer == 1)
    const int ncal = recCal.getRows();
    for (int i = 0; i < ncal; ++i) {
      if (recCal.getInt("layer", i) == 1) {
        h_lv_before->Fill(recCal.getFloat("lv", i));
      }
    }

    // run the fiducial filter (filters REC::Particle rows in-place)
    fid.Run(banks);

    // collect kept track indices from filtered REC::Particle
    std::unordered_set<int> kept;
    for (auto const& row : recPart.getRowList())
      kept.insert(row);

    // fill "after": only PCAL hits whose pindex survived the filter
    for (int i = 0; i < ncal; ++i) {
      if (recCal.getInt("layer", i) == 1) {
        if (kept.count(recCal.getInt("pindex", i))) {
          h_lv_after->Fill(recCal.getFloat("lv", i));
        }
      }
    }

    ++ev;
  }

  // -------- stop algorithm --------
  fid.Stop();

  // -------- draw / save --------
  TCanvas* c = new TCanvas("c", "RGAFiducialFilter PCAL lv (before vs after)", 900, 700);
  h_lv_before->Draw("hist");
  h_lv_after->Draw("hist same");
  auto leg = new TLegend(0.60, 0.75, 0.88, 0.88);
  leg->AddEntry(h_lv_before, "before", "l");
  leg->AddEntry(h_lv_after , "after" , "l");
  leg->Draw();

  if (interactive) {
    std::cout << "\nShowing plot interactively; press Ctrl+C to exit.\n";
    app->Run();
  } else {
    c->SaveAs("out-rga-fiducial-example.png");
    std::cout << "Wrote out-rga-fiducial-example.png\n";
  }

  return 0;
}