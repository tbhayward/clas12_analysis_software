#include "particle_control_core.h"
static bool control_pass(const std::array<double,24>& v,const std::string& period,const GlobalCutConfig& cfg,const std::array<std::array<Window,NV>,3>& cuts,int& top){
 if(!std::isfinite(v[3])||!std::isfinite(v[4])||!std::isfinite(v[23]))return false;
 int d1=int(v[3]),d2=int(v[4]);top=d1==1&&d2==1?0:d1==2&&d2==1?1:d1==2&&d2==0?2:-1;
 if(top<0||is_excluded_run(int(v[23]),cfg))return false;
 if(!passes_analysis_phase_space_cuts(period,v[0],v[1],v[2],v[5],cfg)||!passes_proton_momentum_cut(d1,period,v[8],cfg)||
 !passes_global_cuts(v[2],v[14],v[15],d1,d2,period,v[5],v[6],v[7],v[9],v[10],v[11],v[12],v[13],cfg))return false;
 for(int j=0;j<int(NV);++j)if(!cuts[top][j].pass(v[16+j]))return false;return true;
}
static void cache_control(const json& spec,const json& job,const std::string& out){
 using namespace particle_response;uint64_t entries=0;for(const auto& p:job.at("files")){TFile f(p.get<std::string>().c_str());auto* t=dynamic_cast<TTree*>(f.Get("PhysicsEvents"));if(!t)throw std::runtime_error("Missing control PhysicsEvents");entries+=t->GetEntries();}
 if(!entries)throw std::runtime_error("Empty control input");double probability=std::min(1.,spec.at("expected_cache_events").get<double>()/entries);
 std::ofstream binary(out,std::ios::binary);uint64_t saved=0;json audits=json::array();std::string period=job.at("period");double beam=beam_energy_for_period_label(period);
 for(const auto& p:job.at("files")){
  std::string path=p;TFile f(path.c_str());auto* tree=dynamic_cast<TTree*>(f.Get("PhysicsEvents"));tree->SetCacheSize(32*1024*1024);Reader b(tree);b.pi0_mc=true;
  std::vector<std::string> names={"x","Q2","t1","detector1","detector2","e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi","open_angle_ep2","pTmiss"};
  for(const auto& n:names)b.bind(n,b.reconstructed);bool direct=b.bind("Delta_phi",b.reconstructed,true);
  for(const auto& n:VARIABLES)if(n!="Delta_phi"&&n!="pTmiss")b.bind(n=="theta_gamma_gamma"?"theta_pi0_pi0":n,b.reconstructed);
  bool hasrun=b.bind("runnum",b.reconstructed,true);if(!hasrun&&!spec.at("global_cuts").at("excluded_runs").empty())throw std::runtime_error("Control blacklist needs runnum");
  for(auto n:{"Mh_gammagamma","detector_gamma2","open_angle_egamma1","open_angle_egamma2","gamma_phi1","gamma_phi2"})b.bind(n,b.reconstructed);
  tree->StopCacheLearningPhase();uint64_t key=pathkey(path)^spec.at("seed").get<uint64_t>(),unresolved=0,number=0;
  std::vector<const double*> ptr;for(auto& n:names)ptr.push_back(&b.v.at(n));
  for(auto& n:VARIABLES)ptr.push_back(n=="Delta_phi"&&!direct?nullptr:&b.v.at(n=="theta_gamma_gamma"?"theta_pi0_pi0":n));ptr.push_back(hasrun?&b.v.at("runnum"):nullptr);
  for(Long64_t i=0;i<tree->GetEntries();++i){uint64_t samplekey=mix(key^uint64_t(i));if(uniform(samplekey^0x93c467e37db0c7a5ULL)>=probability)continue;
   b.read(b.reconstructed,i);ControlRecord c{};for(int j=0;j<24;++j)c.base.v[j]=ptr[j]?*ptr[j]:0;
   if(!direct)c.base.v[16]=std::abs(std::remainder(c.base.v[10]-c.base.v[13],2*PI));c.base.key=(samplekey&~uint64_t(1))|uint64_t(!direct);
   c.photons[6]=c.base.v[4];c.photons[7]=b.get("detector_gamma2");if((c.photons[6]!=0&&c.photons[6]!=1)||(c.photons[7]!=0&&c.photons[7]!=1))throw std::runtime_error("Photon detector tags must be 0=FT or 1=FD");double error;
   bool ok=recover_photons(c.base,beam,b.get("gamma_phi1"),b.get("gamma_phi2"),b.get("open_angle_egamma1"),b.get("open_angle_egamma2"),b.get("Mh_gammagamma"),c.photons,error);
   if(!ok)++unresolved;binary.write(reinterpret_cast<const char*>(&c),sizeof c);++saved;++number;
  }
  audits.push_back({{"path",path},{"saved",number},{"unresolved_photons",unresolved},{"audited",0},{"absolute_difference_sum",std::array<double,12>{}},{"absolute_difference_max",std::array<double,12>{}}});
 }
 if(!binary)throw std::runtime_error("Control cache write failed");std::ofstream meta(out+".json");meta<<json{{"kind","two-photon-control-v1"},{"records",saved},{"record_bytes",sizeof(ControlRecord)},{"probability",probability},{"entries",entries},{"audit",audits}};
}
static void evaluate_control(const json& spec,const json& job,const std::string& out){
 using namespace particle_response;std::string period=job.at("period"),code=spec.at("period_codes").at(period);auto cfg=config(spec.at("global_cuts"));
 std::array<std::array<Window,NV>,3> cuts;std::array<std::string,3> topnames={"FD_FD","CD_FD","CD_FT"};
 for(int t=0;t<3;++t){const auto& block=spec.at("exclusivity").at("eppi0_"+code+"_"+topnames[t]).at("mc");for(auto it=block.begin();it!=block.end();++it){auto v=std::find(VARIABLES.begin(),VARIABLES.end(),it.key());if(v==VARIABLES.end())throw std::runtime_error("Unknown control cut");cuts[t][v-VARIABLES.begin()]=window(it.value());}}
 std::vector<std::array<double,NP>> pars;for(auto& p:job.at("parameters"))pars.push_back(p.get<std::array<double,NP>>());
 if(pars.empty()||std::any_of(pars[0].begin(),pars[0].end(),[](double p){return p!=0;}))throw std::runtime_error("Control needs nominal first");
 bool joint=spec.value("joint_control_fit",false);int maxbins=0;std::array<std::vector<double>,NV> edges;
 std::array<std::array<Window,NV>,3> shape_cuts;
 if(joint)for(int t=0;t<3;++t){const auto& block=spec.at("exclusivity").at("eppi0_"+code+"_"+topnames[t]).at("data");for(auto it=block.begin();it!=block.end();++it){auto j=std::find(VARIABLES.begin(),VARIABLES.end(),it.key())-VARIABLES.begin();if(j==NV)throw std::runtime_error("Unknown control shape cut");shape_cuts[t][j]=window(it.value());}}
 if(joint)for(int v=0;v<int(NV);++v){const auto& j=spec.at("variables")[v];int n=j.at("bins");maxbins=std::max(maxbins,n+2);if(j.contains("edges"))edges[v]=j.at("edges").get<std::vector<double>>();else{double lo=j.at("low"),hi=j.at("high");for(int b=0;b<=n;++b)edges[v].push_back(lo+(hi-lo)*b/n);}}
 std::vector<int64_t> shapes(joint?pars.size()*2*3*NV*maxbins:0);
 auto add_shape=[&](const std::array<double,24>& v,int scenario,int half){
   int top=-1;std::array<std::array<Window,NV>,3> empty;
   if(!control_pass(v,period,cfg,empty,top))return;
   std::array<double,NV> obs;for(int j=0;j<int(NV);++j)obs[j]=v[16+j];unsigned failed=failures(obs,shape_cuts[top]);
   for(int j=0;j<int(NV);++j)if(leave_one_out(failed,j)&&std::isfinite(obs[j])){int bin=std::upper_bound(edges[j].begin(),edges[j].end(),obs[j])-edges[j].begin();size_t index=(((size_t(scenario)*2+half)*3+top)*NV+j)*maxbins+bin;++shapes[index];}
 };
 std::vector<int64_t> selected(pars.size()*6),mass_out(pars.size()*2);std::array<int64_t,2> unresolved{},half_records{};ControlRecord c;uint64_t count=0;double beam=beam_energy_for_period_label(period);
 std::ifstream input(job.at("cache").get<std::string>(),std::ios::binary);if(!input)throw std::runtime_error("Missing control cache");
 while(input.read(reinterpret_cast<char*>(&c),sizeof c)){
  ++count;int half=(mix(c.base.key^0xb234fe3adbf120ddULL)>>40)&1,t=-1;bool passed=control_pass(c.base.v,period,cfg,cuts,t);++half_records[half];
  if(joint)add_shape(c.base.v,0,half);
  if(passed){++selected[half*3+t];if(c.photons[9]<.5)++unresolved[half];}
  std::array<double,12> noise;for(int j=0;j<12;++j)noise[j]=gaussian(c.base.key,j);
  auto original=derived(c.base.v,beam,c.base.key&1,MPI0);
  for(size_t s=1;s<pars.size();++s){double mass;auto v=vary_control(c,pars[s],beam,original,noise,mass);if(joint)add_shape(v,s,half);if(control_pass(v,period,cfg,cuts,t))++selected[(s*2+half)*3+t];if(mass<.11||mass>.16)++mass_out[s*2+half];}
  if(c.photons[8]<.11||c.photons[8]>.16)++mass_out[half];
 }
 if(!input.eof()||input.gcount())throw std::runtime_error("Truncated control cache");
 std::ofstream binary(out,std::ios::binary);if(joint)binary.write(reinterpret_cast<char*>(shapes.data()),shapes.size()*sizeof(int64_t));else binary.write(reinterpret_cast<char*>(selected.data()),selected.size()*sizeof(int64_t));
 std::ofstream meta(out+".json");meta<<json{{"hist_shape",{int(pars.size()),2,3,1,int(NV),maxbins}},{"shape",{int(pars.size()),2,3}},{"records",count},{"half_records",half_records},{"selected",selected},{"unresolved_selected",unresolved},{"mass_outside_input_window",mass_out}};
}
