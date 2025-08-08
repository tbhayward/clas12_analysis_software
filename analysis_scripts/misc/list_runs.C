void list_runs() {
    TFile f("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_CH2_eX.root");
    TTree *t = (TTree*)f.Get("PhysicsEvents");
    if(!t) { std::cerr << "Tree not found\n"; return; }

    std::set<int> runs;
    int runnum;
    t->SetBranchAddress("runnum",&runnum);

    Long64_t n = t->GetEntries();
    for(Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        runs.insert(runnum);
    }

    for(auto r : runs) std::cout << r << "\n";
}