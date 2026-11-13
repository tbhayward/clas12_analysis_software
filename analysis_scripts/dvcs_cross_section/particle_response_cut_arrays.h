// Optional individual-photon export for frozen-response kinematic diagnostics.
static std::array<double,8> diagnostic_photon_coordinates(const particle_response::ControlRecord& c,const std::array<double,particle_response::NP>& pars,const std::array<double,12>& noise){
 using namespace particle_response;std::array<double,8> out;out.fill(NAN);out[3]=c.photons[6];out[7]=c.photons[7];
 if(c.photons[9]<.5)return out;
 for(int photon=0;photon<2;++photon){auto v=c.base.v;for(int j=0;j<3;++j)v[11+j]=c.photons[photon*3+j];perturb(v,11,int(c.photons[6+photon])==1?3:4,pars,noise,6+3*photon);if(v[11]<=0)continue;for(int j=0;j<3;++j)out[photon*4+j]=v[11+j];}
 return out;
}
// Export aligned pre-exclusivity event rows for the unchanged production fitter.
// Columns: topology, seven observables, full-global-pass, |t|, production-minimal-pass, MC half.
static void export_cut_arrays(const json& spec,const json& job,const std::string& out){
 using namespace particle_response;const std::string period=job.at("period"),sample=job.at("sample");
 bool control=sample=="pi0_mc"||sample=="pi0_data";double beam=beam_energy_for_period_label(period);auto cfg=config(spec.at("global_cuts"));
 std::array<double,NP> pars{};if(job.contains("parameters"))pars=job.at("parameters").get<std::array<double,NP>>();
 const bool photons=job.value("export_individual_photons",false);const bool extended=job.value("export_kinematics",false)||photons;const int ncols=photons?28:extended?20:12;
 std::ofstream binary(out,std::ios::binary);uint64_t count=0;
 auto emit=[&](const std::array<double,NC>& v,uint64_t key,const std::array<double,8>* photon_values=nullptr){
  if(!std::isfinite(v[3])||!std::isfinite(v[4])||!std::isfinite(v[23]))return;
  int d1=int(v[3]),d2=int(v[4]),top=d1==1&&d2==1?0:d1==2&&d2==1?1:d1==2&&d2==0?2:-1;if(top<0)return;
  if(job.contains("export_topology")&&top!=job.at("export_topology").get<int>())return;
  bool global=!is_excluded_run(int(v[23]),cfg)&&passes_analysis_phase_space_cuts(period,v[0],v[1],v[2],v[5],cfg)&&passes_proton_momentum_cut(d1,period,v[8],cfg)&&passes_global_cuts(v[2],v[14],v[15],d1,d2,period,v[5],v[6],v[7],v[9],v[10],v[11],v[12],v[13],cfg);
  if(job.value("export_global_only",false)&&!global)return;
  int es=fd_sector_from_phi_rad(v[7]),ps=fd_sector_from_phi_rad(v[10]),gs=fd_sector_from_phi_rad(v[13]);
  bool minimal=std::isfinite(v[2])&&std::isfinite(v[14])&&(-v[2])<1.&&v[14]>5.&&es>=1&&
   (d1!=1||(ps>=1&&ps!=es))&&(d2!=1||(gs>=1&&gs!=es))&&(d1!=1||d2!=1||ps!=gs);
  std::array<double,28> row{};row[0]=top;for(int j=0;j<7;++j)row[1+j]=v[16+j];row[8]=global;row[9]=std::abs(v[2]);row[10]=minimal;row[11]=(mix(key^0xb234fe3adbf120ddULL)>>40)&1;
  if(extended){int indices[8]={0,1,5,6,8,9,11,12};for(int j=0;j<8;++j)row[12+j]=v[indices[j]];}
  if(photons){for(int j=20;j<28;++j)row[j]=NAN;if(photon_values){for(int j=0;j<8;++j)row[20+j]=(*photon_values)[j];}else if(!control){row[20]=v[11];row[21]=v[12];row[22]=v[13];row[23]=v[4];}}
  binary.write(reinterpret_cast<const char*>(row.data()),sizeof(double)*ncols);++count;
 };
 if(job.contains("cache")){
  std::ifstream input(job.at("cache").get<std::string>(),std::ios::binary);if(!input)throw std::runtime_error("Missing cut-export cache");
  if(control){ControlRecord c;while(input.read(reinterpret_cast<char*>(&c),sizeof c)){auto base=derived(c.base.v,beam,c.base.key&1,MPI0);std::array<double,12> noise;for(int j=0;j<12;++j)noise[j]=gaussian(c.base.key,j);double mass;std::array<double,8> pp{};if(photons)pp=diagnostic_photon_coordinates(c,pars,noise);emit(vary_control(c,pars,beam,base,noise,mass),c.base.key,photons?&pp:nullptr);}}
  else{Record r;while(input.read(reinterpret_cast<char*>(&r),sizeof r)){if(job.contains("export_topology")){if(!std::isfinite(r.v[3])||!std::isfinite(r.v[4]))continue;int d1=int(r.v[3]),d2=int(r.v[4]);int t=d1==1&&d2==1?0:d1==2&&d2==1?1:d1==2&&d2==0?2:-1;if(t!=job.at("export_topology").get<int>())continue;}auto base=derived(r.v,beam,r.key&1);emit(vary(r,pars,beam,r.key&1,base),r.key);}}
  if(!input.eof()||input.gcount())throw std::runtime_error("Truncated cut-export cache");
 }else{
  // DATA are read once, without sampling, then shared by every response model.
  for(const auto& file:job.at("files")){
   std::string path=file;TFile f(path.c_str());auto* tree=dynamic_cast<TTree*>(f.Get("PhysicsEvents"));if(f.IsZombie()||!tree)throw std::runtime_error("Missing DATA PhysicsEvents");Reader b(tree);if(control&&photons)b.pi0_mc=true;
   tree->SetCacheSize(32*1024*1024);std::vector<std::string> names={"x","Q2","t1","detector1","detector2","e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi","open_angle_ep2","pTmiss"};
   for(const auto& n:names)b.bind(n,b.reconstructed);bool direct=b.bind("Delta_phi",b.reconstructed,true);
   for(const auto& n:VARIABLES)if(n!="Delta_phi"&&n!="pTmiss")b.bind(control&&n=="theta_gamma_gamma"?"theta_pi0_pi0":n,b.reconstructed);
   if(control&&photons)for(auto n:{"gamma_phi1","gamma_phi2","open_angle_egamma1","open_angle_egamma2","Mh_gammagamma","detector_gamma2"})b.bind(n,b.reconstructed);
   bool hasrun=b.bind("runnum",b.reconstructed,true);if(!hasrun&&!cfg.excluded_runs.empty())throw std::runtime_error("DATA blacklist requires runnum");tree->StopCacheLearningPhase();
   std::vector<const double*> ptr;for(const auto& n:names)ptr.push_back(&b.v.at(n));for(const auto& n:VARIABLES)ptr.push_back(n=="Delta_phi"&&!direct?nullptr:&b.v.at(control&&n=="theta_gamma_gamma"?"theta_pi0_pi0":n));ptr.push_back(hasrun?&b.v.at("runnum"):nullptr);
   for(Long64_t i=0;i<tree->GetEntries();++i){b.read(b.reconstructed,i);Record r{};for(int j=0;j<24;++j)r.v[j]=ptr[j]?*ptr[j]:0;if(!direct)r.v[16]=std::abs(std::remainder(r.v[10]-r.v[13],2*PI));if(control&&photons){ControlRecord c{};c.base=r;c.photons[6]=r.v[4];c.photons[7]=b.get("detector_gamma2");double error;recover_photons(r,beam,b.get("gamma_phi1"),b.get("gamma_phi2"),b.get("open_angle_egamma1"),b.get("open_angle_egamma2"),b.get("Mh_gammagamma"),c.photons,error);std::array<double,12> z{};std::array<double,NP> zero{};auto pp=diagnostic_photon_coordinates(c,zero,z);emit(r.v,mix(pathkey(path)^uint64_t(i)),&pp);}else emit(r.v,mix(pathkey(path)^uint64_t(i)));}
  }
 }
 if(!binary)throw std::runtime_error("Cut-export write failed");std::ofstream meta(out+".json");meta<<json{{"shape",{count,ncols}},{"columns",{"topology","Delta_phi","theta","theta_gamma_gamma","pTmiss","Emiss2","Mx2","Mx2_2","global_pass","t_abs","minimal_pass","half"}},{"individual_photon_columns",photons?json({"gamma1_p","gamma1_theta_rad","gamma1_phi_rad","gamma1_detector","gamma2_p","gamma2_theta_rad","gamma2_phi_rad","gamma2_detector"}):json::array()},{"kinematic_columns",extended?json({"xB","Q2","electron_p","electron_theta_rad","proton_p","proton_theta_rad","photon_p","photon_theta_rad"}):json::array()}};
}
