// Reuse the tested scalar ROOT reader, saved-config loader and cut parser.
#define main exclusivity_diagnostic_main_unused
#include "exclusivity_resolution_worker.cpp"
#undef main
#include "particle_response_core.h"
#include <cstring>
#include <numeric>
using particle_response::Record;
static uint64_t pathkey(const std::string& s){uint64_t h=1469598103934665603ULL;for(unsigned char c:s){h^=c;h*=1099511628211ULL;}return h;}
static void cache_events(const json& spec,const json& job,const std::string& out){
 uint64_t entries=0;for(const auto& p:job.at("files")){TFile f(p.get<std::string>().c_str());auto* t=dynamic_cast<TTree*>(f.Get("PhysicsEvents"));if(f.IsZombie()||!t)throw std::runtime_error("Missing reconstructed input");entries+=t->GetEntries();}
 if(!entries)throw std::runtime_error("Empty reconstructed sample");
 double probability=std::min(1.,spec.at("expected_cache_events").get<double>()/entries);
 std::ofstream binary(out,std::ios::binary);if(!binary)throw std::runtime_error("Cannot write cache");
 uint64_t saved=0;json audits=json::array();
 for(const auto& input:job.at("files")){
  std::string path=input;TFile f(path.c_str());auto* tree=dynamic_cast<TTree*>(f.Get("PhysicsEvents"));
  tree->SetCacheSize(32*1024*1024);Reader b(tree);
  std::vector<std::string> names={"x","Q2","t1","detector1","detector2","e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi","open_angle_ep2","pTmiss"};
  for(const auto& n:names)b.bind(n,b.reconstructed);
  bool direct=b.bind("Delta_phi",b.reconstructed,true);
  for(const auto& n:VARIABLES)if(n!="Delta_phi"&&n!="pTmiss")b.bind(n,b.reconstructed);
  bool hasrun=b.bind("runnum",b.reconstructed,true);
  if(!hasrun&&!spec.at("global_cuts").at("excluded_runs").empty())throw std::runtime_error("Run blacklist needs runnum");
  tree->StopCacheLearningPhase();uint64_t key=pathkey(path)^spec.at("seed").get<uint64_t>();
  std::array<double,12> abs_sum{},abs_max{};uint64_t n_audit=0;
  std::vector<const double*> ptr;for(const auto& n:names)ptr.push_back(&b.v.at(n));
  for(const auto& n:VARIABLES)ptr.push_back(n=="Delta_phi"&&!direct?nullptr:&b.v.at(n));
  ptr.push_back(hasrun?&b.v.at("runnum"):nullptr);
  uint64_t file_saved=0;
  for(Long64_t i=0;i<tree->GetEntries();++i){uint64_t samplekey=particle_response::mix(key^uint64_t(i));if(particle_response::uniform(samplekey^0x93c467e37db0c7a5ULL)>=probability)continue;
   b.read(b.reconstructed,i);Record r{};
   for(int j=0;j<24;++j)r.v[j]=ptr[j]?*ptr[j]:0;
   if(!direct)r.v[16]=std::abs(std::remainder(r.v[10]-r.v[13],2*particle_response::PI));
   // Flags occupy low bits; separate stable hash bits determine train/test.
   r.key=(samplekey&~uint64_t(1))|uint64_t(!direct);
   binary.write(reinterpret_cast<const char*>(&r),sizeof r);++saved;++file_saved;
   if(std::isfinite(r.v[5])&&r.v[5]>0&&r.v[8]>0&&r.v[11]>0){
    auto calc=particle_response::derived(r.v,beam_energy_for_period_label(job.at("period")),!direct);
    int dest[12]={0,1,2,14,15,16,17,18,19,20,21,22};
    bool finite=true;for(int j=0;j<12;++j)finite=finite&&std::isfinite(calc[j])&&std::isfinite(r.v[dest[j]]);
    if(finite){++n_audit;for(int j=0;j<12;++j){double d=calc[j]-r.v[dest[j]];if(j==5)d=std::remainder(d,2*particle_response::PI);d=std::abs(d);abs_sum[j]+=d;abs_max[j]=std::max(abs_max[j],d);}}
   }
  }
  audits.push_back({{"path",path},{"saved",file_saved},{"audited",n_audit},{"absolute_difference_sum",abs_sum},{"absolute_difference_max",abs_max},{"delta_phi_fallback",!direct}});
 }
 if(!binary)throw std::runtime_error("Cache write failed");
 json meta={{"records",saved},{"record_bytes",sizeof(Record)},{"probability",probability},{"entries",entries},{"audit",audits}};
 std::ofstream m(out+".json");m<<meta;
}
static void evaluate(const json& spec,const json& job,const std::string& out){
 using namespace particle_response;
 auto cfg=config(spec.at("global_cuts"));auto no_pt=cfg;no_pt.enable_pTmiss_cut=false;
 std::string period=job.at("period"),code=spec.at("period_codes").at(period);
 std::array<std::array<Window,NV>,3> cuts;
 std::array<std::string,3> topnames={"FD_FD","CD_FD","CD_FT"};
 for(int t=0;t<3;++t)for(auto it=spec.at("exclusivity").at("DVCS_"+code+"_"+topnames[t]).at("data").begin();it!=spec.at("exclusivity").at("DVCS_"+code+"_"+topnames[t]).at("data").end();++it){
  auto pos=std::find(VARIABLES.begin(),VARIABLES.end(),it.key());if(pos==VARIABLES.end())throw std::runtime_error("Unsupported production cut");cuts[t][pos-VARIABLES.begin()]=window(it.value());}
 std::vector<std::array<double,NP>> parameters;for(const auto& p:job.at("parameters"))parameters.push_back(p.get<std::array<double,NP>>());
 if(parameters.empty()||std::any_of(parameters.front().begin(),parameters.front().end(),[](double p){return p!=0;}))throw std::runtime_error("First response scenario must be nominal");
 int ns=parameters.size(),ng=job.value("detailed",false)?int(spec.at("strata").size()):1;
 int maxbins=0;std::array<std::vector<double>,NV> edges;
 for(int v=0;v<int(NV);++v){const auto& j=spec.at("variables")[v];int n=j.at("bins");maxbins=std::max(maxbins,n+2);if(j.contains("edges"))edges[v]=j.at("edges").get<std::vector<double>>();else{double lo=j.at("low"),hi=j.at("high");for(int b=0;b<=n;++b)edges[v].push_back(lo+(hi-lo)*b/n);}}
 size_t stride=size_t(2)*3*ng*NV*maxbins;
 std::vector<int64_t> counts(size_t(ns)*stride),selected(size_t(ns)*2*3);
 std::vector<JointGrid> grids;std::vector<Stratum> strata;std::vector<int> ai;
 if(ng>1){for(const auto& j:spec.at("joint_grids"))grids.push_back({j.at("offset"),j.at("indices").get<std::vector<int>>(),j.at("edges").get<std::vector<std::vector<double>>>()});for(const auto& j:spec.at("strata"))strata.push_back({j.at("axis"),j.value("low",0.),j.value("high",0.)});ai=spec.at("strata_axis_index").get<std::vector<int>>();}
 std::vector<size_t> baseline_bins;baseline_bins.reserve(70);int baseline_selected=-1;
 auto add=[&](const std::array<double,NC>& v,int scenario,int half,int sign){
  if(!std::isfinite(v[3])||!std::isfinite(v[4]))return;
  int d1=int(v[3]),d2=int(v[4]),t=d1==1&&d2==1?0:d1==2&&d2==1?1:d1==2&&d2==0?2:-1;
  if(t<0||is_excluded_run(int(v[23]),cfg))return;
  auto global=[&](const GlobalCutConfig& c){return passes_analysis_phase_space_cuts(period,v[0],v[1],v[2],v[5],c)&&passes_proton_momentum_cut(d1,period,v[8],c)&&passes_global_cuts(v[2],v[14],v[15],d1,d2,period,v[5],v[6],v[7],v[9],v[10],v[11],v[12],v[13],c);};
  bool gp=global(cfg);if(!gp&&!(cfg.enable_pTmiss_cut&&global(no_pt)))return;
  std::array<double,NV> values;for(int j=0;j<int(NV);++j)values[j]=v[16+j];unsigned failed=failures(values,cuts[t]);
  if(gp&&!failed){selected[(scenario*2+half)*3+t]+=sign;if(scenario==0)baseline_selected=half*3+t;}
  std::vector<int> groups={0};
  if(ng>1){std::array<double,8> a={v[5],v[6]*DEG,v[8],v[9]*DEG,v[11],v[12]*DEG,v[0],v[1]};for(int g=1;g<ng;++g)if(ai[g]>=0&&strata[g].contains(a[ai[g]]))groups.push_back(g);for(const auto& grid:grids){int g=grid.find(a);if(g>=0)groups.push_back(g);}}
  for(int j=0;j<int(NV);++j){if(!leave_one_out(failed,j)||(!gp&&j!=3)||!std::isfinite(values[j]))continue;const auto& e=edges[j];int bin=std::upper_bound(e.begin(),e.end(),values[j])-e.begin();for(int g:groups){size_t index=((((size_t(scenario)*2+half)*3+t)*ng+g)*NV+j)*maxbins+bin;counts[index]+=sign;if(scenario==0)baseline_bins.push_back(index);}}
 };
 std::ifstream input(job.at("cache").get<std::string>(),std::ios::binary);if(!input)throw std::runtime_error("Cannot read event cache");
 Record r;double beam=beam_energy_for_period_label(period);uint64_t nr=0;
 while(input.read(reinterpret_cast<char*>(&r),sizeof r)){
  ++nr;if(!std::isfinite(r.v[3])||!std::isfinite(r.v[4])||!std::isfinite(r.v[23]))continue;int half=(mix(r.key^0xb234fe3adbf120ddULL)>>40)&1;bool fallback=r.key&1;
  auto base=derived(r.v,beam,fallback);std::array<double,9> noise;for(int j=0;j<9;++j)noise[j]=gaussian(r.key,j);
  // Scenario zero must be nominal. Store differences only for affected tracks.
  baseline_bins.clear();baseline_selected=-1;add(r.v,0,half,1);
  for(int s=1;s<ns;++s){bool affected=false;int regions[3]={0,int(r.v[3])==1?1:2,int(r.v[4])==1?3:4};for(int reg:regions)for(int c=0;c<4;++c)affected=affected||parameters[s][4*reg+c]!=0;if(!affected)continue;
   for(size_t bin:baseline_bins)--counts[size_t(s)*stride+bin];if(baseline_selected>=0)--selected[s*6+baseline_selected];
   add(vary(r,parameters[s],beam,fallback,base,&noise),s,half,1);
  }
 }
 if(!input.eof()||input.gcount()!=0)throw std::runtime_error("Truncated event cache");
 for(int s=1;s<ns;++s){for(size_t i=0;i<stride;++i)counts[size_t(s)*stride+i]+=counts[i];for(int i=0;i<6;++i)selected[s*6+i]+=selected[i];}
 if(std::any_of(counts.begin(),counts.end(),[](auto n){return n<0;}))throw std::runtime_error("Negative response histogram");
 std::ofstream binary(out,std::ios::binary);binary.write(reinterpret_cast<const char*>(counts.data()),counts.size()*sizeof(int64_t));if(!binary)throw std::runtime_error("Response write failed");
 std::ofstream meta(out+".json");meta<<json{{"shape",{ns,2,3,ng,int(NV),maxbins}},{"records",nr},{"selected",selected}};
}
#include "particle_control_worker.h"

int main(int argc,char** argv){try{if(argc!=5)throw std::runtime_error("worker cache/evaluate SPEC JOB OUTPUT");auto spec=read_json(argv[2]),job=read_json(argv[3]);std::string mode=argv[1];if(mode=="cache")cache_events(spec,job,argv[4]);else if(mode=="cache-control")cache_control(spec,job,argv[4]);else if(mode=="evaluate"){if(job.at("sample")=="pi0_mc")evaluate_control(spec,job,argv[4]);else evaluate(spec,job,argv[4]);}else throw std::runtime_error("Unknown mode");}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
