static std::mutex pi0_recovery_mutex;
static std::map<std::string,std::array<uint64_t,3>> pi0_recovery_counts;
static std::map<std::string,std::array<uint64_t,2>> pi0_recovery_bins;
static thread_local std::map<std::string,std::array<uint64_t,3>> pi0_recovery_thread_counts;
static thread_local std::map<std::string,std::array<uint64_t,2>> pi0_recovery_thread_bins;
// Included only in an isolated copy of total_counts.cpp, after BranchBinder.
static void apply_pi0_test_response(BranchBinder& b,const WorkConfig& work,const PeriodTags& tags,Long64_t entry){
 using namespace particle_response;
 b.response_unresolved=false;
 if(pi0_test_parameters.empty()||work.sample_kind!=SampleKind::MC_REC)return;
 if(work.channel_cfg.channel!=Channel::EPPI0&&work.channel_cfg.channel!=Channel::EPPI0_BKG_AS_DVCS)return;
 auto pars=pi0_test_parameters.at(tags.period_label).get<std::array<double,NP>>();
 bool active=false;for(double x:pars)active|=x!=0.;if(!active)return;
 if(!(b.has_e_p&&b.has_e_theta&&b.has_e_phi&&b.has_p1_p&&b.has_p1_theta&&b.has_p1_phi&&b.has_p2_p&&b.has_p2_theta&&b.has_p2_phi))throw std::runtime_error("Morph requires all measured four-vector coordinates");
 Record r{};r.v={b.x,b.Q2,b.t1,b.detector1,b.detector2,b.e_p,b.e_theta,b.e_phi,b.p1_p,b.p1_theta,b.p1_phi,b.p2_p,b.p2_theta,b.p2_phi,b.open_angle_ep2,b.pTmiss,b.Delta_phi,b.theta,work.channel_cfg.channel==Channel::EPPI0?b.theta_pi0_pi0:b.theta_gamma_gamma,b.pTmiss,b.Emiss2,b.Mx2,b.Mx2_2,b.runnum};
 uint64_t key=1469598103934665603ULL;for(unsigned char c:tags.tree_key){key^=c;key*=1099511628211ULL;}r.key=mix(key^uint64_t(entry)^pi0_test_seed);
 bool control=work.channel_cfg.channel==Channel::EPPI0;double beam=b.bin_beam_energy;
 auto original=derived(r.v,beam,!b.has_Delta_phi,control?MPI0:0.);std::array<double,NC> varied;
 if(control){
  for(bool ok:b.response_extra_ok)if(!ok)throw std::runtime_error("Direct pi0 morph requires individual-photon recovery branches");
  ControlRecord c{};c.base=r;c.photons[6]=b.detector2;c.photons[7]=b.response_extra[6];double error;
  bool recovered=recover_photons(r,beam,b.response_extra[0],b.response_extra[1],b.response_extra[2],b.response_extra[3],b.response_extra[4],c.photons,error);
  {auto& a=pi0_recovery_thread_counts[tags.tree_key];++a[0];if(!recovered)++a[1];}
  if(!recovered){b.response_unresolved=true;return;} // Match prior vary_control behavior: retain original event exactly.

  std::array<double,12> noise;for(int i=0;i<12;++i)noise[i]=gaussian(r.key,i);double mass;varied=vary_control(c,pars,beam,original,noise,mass);
 }else varied=vary(r,pars,beam,!b.has_Delta_phi,original);
 auto plane=[&](const std::array<double,NC>& v){auto k=particle(std::sqrt(beam*beam-ME*ME),0,0,ME);auto e=particle(v[5],v[6],v[7],ME);auto g=particle(v[11],v[12],v[13],control?MPI0:0.);return trento((k-e).p,e.p,g.p);};
 // Preserve the stored phi convention/rounding offset, applying only its physical change.
 b.phi2=wrap(b.phi2+std::remainder(plane(varied)-plane(r.v),2*particle_response::PI));
 b.x=varied[0];b.Q2=varied[1];b.t1=varied[2];b.e_p=varied[5];b.e_theta=varied[6];b.e_phi=varied[7];b.p1_p=varied[8];b.p1_theta=varied[9];b.p1_phi=varied[10];b.p2_p=varied[11];b.p2_theta=varied[12];b.p2_phi=varied[13];b.open_angle_ep2=varied[14];b.pTmiss=varied[15];b.Delta_phi=varied[16];b.theta=varied[17];if(control)b.theta_pi0_pi0=varied[18];else b.theta_gamma_gamma=varied[18];b.Emiss2=varied[20];b.Mx2=varied[21];b.Mx2_2=varied[22];
}

static void audit_pi0_recovery_row(const BranchBinder& b,const WorkConfig& work,const PeriodTags& tags,int row){
 if(pi0_test_parameters.empty()||work.channel_cfg.channel!=Channel::EPPI0||work.sample_kind!=SampleKind::MC_REC)return;
 auto& c=pi0_recovery_thread_bins[tags.tree_key+"|"+std::to_string(row)];++c[0];if(b.response_unresolved)++c[1];
}

static void flush_pi0_recovery_audit(){
 std::lock_guard<std::mutex> lock(pi0_recovery_mutex);
 for(const auto& p:pi0_recovery_thread_counts)for(int i=0;i<3;++i)pi0_recovery_counts[p.first][i]+=p.second[i];
 for(const auto& p:pi0_recovery_thread_bins)for(int i=0;i<2;++i)pi0_recovery_bins[p.first][i]+=p.second[i];
 pi0_recovery_thread_counts.clear();pi0_recovery_thread_bins.clear();
}
