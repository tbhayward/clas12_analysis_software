// run `calculateGapCharge.groovy` first, then this macro

std::optional<Double_t> GetChargeFromValerii(int run_num);

void analyzeGapCharge(TString datFileN = "charge_gaps.dat") {

  // read the input
  auto tr = new TTree("tr", "tr");
  tr->ReadFile(datFileN);

  // set tree branch addresses
  Int_t runnum, binnum, runset, golden, lastBin;
  Char_t comment[1000];
  Double_t chargeGap, chargeMin, chargeMax;
  Double_t chargeGapMin=1e6;
  Double_t chargeGapMax=-1e6;
  Double_t chargeRatMin=1e6;
  Double_t chargeRatMax=-1e6;
  tr->SetBranchAddress("runnum", &runnum);
  tr->SetBranchAddress("binnum", &binnum);
  tr->SetBranchAddress("runset", &runset);
  tr->SetBranchAddress("golden", &golden);
  tr->SetBranchAddress("lastBin", &lastBin);
  tr->SetBranchAddress("comment", comment);
  tr->SetBranchAddress("fcChargeMin", &chargeMin);
  tr->SetBranchAddress("fcChargeMax", &chargeMax);
  tr->SetBranchAddress("fcChargeGapToNextBin", &chargeGap);

  // get plot ranges and initialize charge accumulators
  std::map<Int_t, TString> runsets;
  runsets[-1] = "full_run_list";
  enum charge_calc { kAllFiles, kGoldenFiles, kFullRunMin, kFullRunMax };
  std::map<Int_t, std::map<charge_calc, Double_t>> chargeTot;
  std::map<Int_t, Int_t> numBinsPerRun;
  for(Long64_t e = 0; e < tr->GetEntries(); e++) {
    tr->GetEntry(e);
    runsets[runset] = comment;
    chargeTot[runnum] = {
      {kAllFiles,    0},
      {kGoldenFiles, 0},
      {kFullRunMin,  1e7},
      {kFullRunMax,  0}
    };
    if(!lastBin) {
      chargeGapMin = std::min(chargeGapMin, chargeGap);
      chargeGapMax = std::max(chargeGapMax, chargeGap);
      chargeRatMin = std::min(chargeRatMin, chargeGap / (chargeMax-chargeMin));
      chargeRatMax = std::max(chargeRatMax, chargeGap / (chargeMax-chargeMin));
    }
    else {
      int numFiles = binnum/5 + 1;
      // std::cout << runnum << " " << binnum << " " << numFiles << std::endl;
      numBinsPerRun.insert({runnum, numFiles});
    }
  }

  // define plots
  std::map<Int_t, TH1D*> gapDists;
  std::map<Int_t, TH1D*> ratDists;
  std::map<Int_t, TGraph*> gapVsBinNum;
  std::map<Int_t, TGraph*> ratVsBinNum;
  for(auto const& [r, c] : runsets) {
    gapDists.insert({r, new TH1D(
          Form("gapDist_%d", r),
          TString("Gap charge distribution: ") + c + ";gap [nC]",
          100000,
          chargeGapMin,
          chargeGapMax
          )});
    ratDists.insert({r, new TH1D(
          Form("ratDist_%d", r),
          TString("Gap charge / bin charge: ") + c + ";gap/bin",
          100000,
          chargeRatMin,
          chargeRatMax
          )});
    gapVsBinNum.insert({r, new TGraph()});
    gapVsBinNum[r]->SetName(Form("gapVsBinNum_%d", r));
    gapVsBinNum[r]->SetTitle(TString("Gap charge vs bin num: ") + c + TString(";runnum + binnum/1e5;gap [nC]"));
    gapVsBinNum[r]->SetMarkerStyle(kFullCircle);
    ratVsBinNum.insert({r, new TGraph()});
    ratVsBinNum[r]->SetName(Form("ratVsBinNum_%d", r));
    ratVsBinNum[r]->SetTitle(TString("Gap charge / bin charge vs bin num: ") + c + TString(";runnum + binnum/1e5;gap/bin"));
    ratVsBinNum[r]->SetMarkerStyle(kFullCircle);
  }

  // fill plots and accumulate charge
  for(Long64_t e = 0; e < tr->GetEntries(); e++) {
    tr->GetEntry(e);
    auto charge = chargeMax - chargeMin;
    if(golden==1 && lastBin==0) {
      for(auto const& [r, c] : runsets) {
        if(r >= 0 && runset != r) continue;
        auto chargeRat = chargeGap / charge;
        auto runbin = runnum + static_cast<double>(binnum)/1e5;
        gapDists.at(r)->Fill(chargeGap);
        ratDists.at(r)->Fill(chargeRat);
        gapVsBinNum.at(r)->AddPoint(runbin, chargeGap);
        ratVsBinNum.at(r)->AddPoint(runbin, chargeRat);
      }
    }
    chargeTot[runnum][kAllFiles] += charge;
    if(golden) chargeTot[runnum][kGoldenFiles] += charge;
    chargeTot[runnum][kFullRunMin] = std::min(chargeMin, chargeTot[runnum][kFullRunMin]);
    chargeTot[runnum][kFullRunMax] = std::max(chargeMax, chargeTot[runnum][kFullRunMax]);
  }

  // cross check the accumulated charge
  enum chargeGrEnum { gGoldenSum, gAllBinsSum, gRunMaxMinusMin, gFracDiff, nChargeGr };
  TGraph* chargeGr[nChargeGr];
  for(int g=0; g<nChargeGr; g++)
    chargeGr[g] = new TGraph();
  chargeGr[gGoldenSum]->SetName("chargeGrGoldenSum");
  chargeGr[gAllBinsSum]->SetName("chargeGrAllBinsSum");
  chargeGr[gRunMaxMinusMin]->SetName("chargeGrRunMaxMinusMin");
  chargeGr[gFracDiff]->SetName("chargeGrFracDiff");
  auto chargeFracDiffDist = new TH1D("chargeFracDiffDist", "Charge Fractional Difference per run, weighted by num. QA bins;frac. diff.", 50, -0.1, 0.1);
  std::cout << "CHARGE CROSS CHECK" << std::endl << "==================" << std::endl;
  std::cout << "    run num       (1) golden QA bins' sum       (2) all QA bins' sum       (3) run's max-min       100 * [1 - (2)/(3)]" << std::endl;
  for(auto const& [run_num, charge_hash] : chargeTot) {
    auto runCharge = charge_hash.at(kFullRunMax) - charge_hash.at(kFullRunMin);
    auto fracDiff  = (1 - charge_hash.at(kAllFiles) / runCharge);
    chargeGr[gGoldenSum]->AddPoint(run_num,charge_hash.at(kGoldenFiles));
    chargeGr[gAllBinsSum]->AddPoint(run_num,charge_hash.at(kAllFiles));
    chargeGr[gRunMaxMinusMin]->AddPoint(run_num,runCharge);
    chargeGr[gFracDiff]->AddPoint(run_num,fracDiff);
    chargeFracDiffDist->Fill(fracDiff, numBinsPerRun.at(run_num));
    std::cout << ">>> " <<
      run_num << "       " <<
      charge_hash.at(kGoldenFiles) << "       " <<
      charge_hash.at(kAllFiles) << "       " <<
      runCharge << "       " <<
      100.0 * fracDiff << "%" <<
      std::endl;
  }
  // cross check with Valerii
  /*
  std::cout << "CROSS CHECK WITH VALERII" << std::endl << "==================" << std::endl;
  std::cout << "QADB       Valerii       Diff" << std::endl;
  for(auto const& [run_num, charge_hash] : chargeTot) {
    auto qadbCharge    = charge_hash.at(kGoldenFiles);
    auto valeriiCharge = GetChargeFromValerii(run_num);
    if(valeriiCharge)
      std::cout << ">>> " << std::setprecision(10) <<
        qadbCharge << "       " <<
        valeriiCharge.value() << "       " <<
        std::round(qadbCharge)-std::round(valeriiCharge.value()) <<
        std::endl;
  }
  */

  // draw plots
  Int_t padnum;

  auto canv0 = new TCanvas();
  canv0->Divide(1,3);
  for(int i=1; i<=3; i++)
    canv0->GetPad(i)->SetGrid(1,1);
  chargeGr[gGoldenSum]->SetMarkerColor(kOrange);
  chargeGr[gAllBinsSum]->SetMarkerColor(kRed+1);
  chargeGr[gRunMaxMinusMin]->SetMarkerColor(kBlue-4);
  chargeGr[gFracDiff]->SetMarkerColor(kGreen+1);
  for(int g=0; g<nChargeGr; g++) {
    chargeGr[g]->SetMarkerStyle(kFullCircle);
  }
  auto chargeMgr = new TMultiGraph();
  chargeMgr->SetName("chargeMgr");
  chargeMgr->Add(chargeGr[gGoldenSum]);
  chargeMgr->Add(chargeGr[gAllBinsSum]);
  chargeMgr->Add(chargeGr[gRunMaxMinusMin]);
  chargeMgr->SetTitle("Charge per run;run number;charge [nC]");
  chargeGr[gFracDiff]->SetTitle("Charge Fractional Difference per run;run number;frac. diff.");
  canv0->cd(1);
  chargeMgr->Draw("ap");
  canv0->cd(2);
  chargeGr[gFracDiff]->Draw("ap");
  canv0->cd(3);
  chargeFracDiffDist->SetLineColor(kGreen+1);
  chargeFracDiffDist->SetLineWidth(3);
  chargeFracDiffDist->Draw();

  auto canv1 = new TCanvas();
  canv1->Divide(2,2);
  padnum = 1;
  for(auto const& [r, c] : runsets) {
    auto pad = canv1->GetPad(padnum++);
    // pad->SetLogy();
    pad->cd();
    gapDists.at(r)->Draw();
    // gapDists.at(r)->GetXaxis()->SetRangeUser(-1000, 1000);
  }
  auto canv2 = new TCanvas();
  canv2->Divide(2,2);
  padnum = 1;
  for(auto const& [r, c] : runsets) {
    auto pad = canv2->GetPad(padnum++);
    // pad->SetLogy();
    pad->cd();
    ratDists.at(r)->Draw();
    // ratDists.at(r)->GetXaxis()->SetRangeUser(-1000, 1000);
  }
  auto canv3 = new TCanvas();
  canv3->Divide(2,2);
  padnum = 1;
  for(auto const& [r, c] : runsets) {
    auto pad = canv3->GetPad(padnum++);
    pad->cd();
    gapVsBinNum.at(r)->Draw("ap");
  }
  auto canv4 = new TCanvas();
  canv4->Divide(2,2);
  padnum = 1;
  for(auto const& [r, c] : runsets) {
    auto pad = canv4->GetPad(padnum++);
    pad->cd();
    ratVsBinNum.at(r)->Draw("ap");
  }
}

// Valerii's charge per run, for golden QA cuts
std::optional<Double_t> GetChargeFromValerii(int run_num) {
  std::unordered_map<Int_t, Double_t> chg = {
    {5032, 34245.6},
    {5036, 114022},
    {5038, 378709},
    {5039, 96850.5},
    {5040, 147343},
    {5041, 175940},
    {5043, 140423},
    {5045, 137370},
    {5052, 34382.7},
    {5053, 52192.9},
    {5116, 38416.6},
    {5117, 393660},
    {5119, 12884.6},
    {5120, 78582.2},
    {5124, 395035},
    {5125, 397571},
    {5126, 397503},
    {5127, 57527.1},
    {5139, 139615},
    {5153, 69038.4},
    {5158, 12711.2},
    {5162, 80069.8},
    {5163, 38739.3},
    {5164, 5041.72},
    {5181, 21149},
    {5191, 86529.9},
    {5193, 154744},
    {5195, 75158.6},
    {5196, 411810},
    {5197, 414557},
    {5198, 411532},
    {5199, 412208},
    {5200, 154042},
    {5201, 80925.5},
    {5202, 410589},
    {5203, 414655},
    {5204, 398399},
    {5205, 90713.7},
    {5206, 289584},
    {5208, 333470},
    {5211, 87206.4},
    {5212, 391124},
    {5215, 386497},
    {5216, 165194},
    {5219, 385594},
    {5220, 350442},
    {5221, 389766},
    {5222, 388895},
    {5223, 135114},
    {5230, 371039},
    {5231, 390388},
    {5232, 391028},
    {5233, 389224},
    {5234, 391996},
    {5235, 145449},
    {5237, 156465},
    {5238, 310085},
    {5248, 401352},
    {5345, 382028},
    {5346, 379346},
    {5347, 178928},
    {5349, 242454},
    {5351, 52604.3},
    {5354, 377576},
    {5355, 170245},
    {5367, 373381},
    {5342, 115739},
    {5343, 60278.6},
    {5344, 152308},
    {5356, 377070},
    {5357, 377122},
    {5358, 382403},
    {5359, 382255},
    {5360, 378618},
    {5361, 372562},
    {5362, 28149.7},
    {5366, 393710},
    {5368, 327150},
    {5369, 94684.2},
    {5372, 369081},
    {5373, 361669},
    {5374, 382198},
    {5375, 387490},
    {5376, 290796},
    {5377, 12941.5},
    {5378, 36996.3},
    {5379, 387531},
    {5380, 357707},
    {5381, 371711},
    {5383, 281662},
    {5386, 187210},
    {5390, 26329.6},
    {5391, 365707},
    {5392, 198954},
    {5393, 363850},
    {5398, 42407.6},
    {5401, 31765.7},
    {5403, 110510},
    {5404, 14970.5},
    {5406, 179818},
    {5407, 416414},
    {5249, 434539},
    {5252, 272836},
    {5253, 98629.5},
    {5257, 430276},
    {5258, 429167},
    {5259, 161727},
    {5261, 430571},
    {5262, 202647},
    {5303, 420564},
    {5304, 406558},
    {5305, 157358},
    {5306, 401701},
    {5307, 47816.4},
    {5310, 105109},
    {5311, 138911},
    {5315, 97427.7},
    {5317, 401966},
    {5318, 400666},
    {5319, 405694},
    {5320, 200445},
    {5323, 27479.1},
    {5324, 33812.9},
    {5335, 13391.9},
    {5339, 95397.9},
    {5340, 35553.8}
  };
  try {
    return chg.at(run_num);
  } catch(std::exception const& ex) {
    std::cerr << "ERROR: Valerii did not give a charge for run " << run_num << endl;
    return {};
  }
}
