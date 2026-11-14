// Isolated selection-first evaluator; all production response laws are reused.
#include "particle_response_worker.cpp"
#include "unordered_control_arrays.h"
using namespace particle_response;
static int pair_group(const ControlRecord& c){
 int p=int(c.base.v[3]),a=int(c.photons[6]),b=int(c.photons[7]);
 if((p!=1&&p!=2)||(a!=0&&a!=1)||(b!=0&&b!=1))return -1;
 int raw=(p==2?3:0)+2-a-b;return raw==0?0:raw>=3?raw-2:-1;
}
static bool selected_global(const std::array<double,NC>& v,const std::string& period,const GlobalCutConfig& cfg,int tag){
 return std::isfinite(v[0])&&std::isfinite(v[1])&&std::isfinite(v[2])&&v[5]>0&&v[8]>0&&v[11]>0&&std::isfinite(v[3])&&std::isfinite(v[23])&&!is_excluded_run(int(v[23]),cfg)&&passes_analysis_phase_space_cuts(period,v[0],v[1],v[2],v[5],cfg)&&passes_proton_momentum_cut(int(v[3]),period,v[8],cfg)&&passes_global_cuts(v[2],v[14],v[15],int(v[3]),tag,period,v[5],v[6],v[7],v[9],v[10],v[11],v[12],v[13],cfg);
}
static double photon_mass(const ControlRecord& c){
 auto a=particle(c.photons[0],c.photons[1],c.photons[2],0),b=particle(c.photons[3],c.photons[4],c.photons[5],0);return std::sqrt(std::max(0.,(a+b).m2()));
}
static void selected_evaluate(const json& spec,const json& job,const std::string& output){
 std::string period=job.at("period"),sample=job.at("sample");bool control=sample=="pi0_mc";int groups=control?4:3,nv=control?8:7,maxbins=0;
 auto cfg=config(spec.at("global_cuts"));double beam=beam_energy_for_period_label(period);
 std::vector<std::array<double,NP>> pars;for(auto& p:job.at("parameters"))pars.push_back(p.get<std::array<double,NP>>());
 if(pars.empty()||std::any_of(pars[0].begin(),pars[0].end(),[](double x){return x!=0;}))throw std::runtime_error("Nominal must be scenario zero");
 std::vector<bool> active(pars.size());for(size_t s=0;s<pars.size();++s)active[s]=std::any_of(pars[s].begin(),pars[s].end(),[](double x){return x!=0;});
 std::vector<std::vector<double>> edges;for(int j=0;j<nv;++j){const auto& v=j==7?spec.at("mass_variable"):spec.at("variables")[j];std::vector<double> e;if(v.contains("edges"))e=v.at("edges").get<std::vector<double>>();else{int n=v.at("bins");double lo=v.at("low"),hi=v.at("high");for(int b=0;b<=n;++b)e.push_back(lo+(hi-lo)*b/n);}maxbins=std::max(maxbins,int(e.size())+1);edges.push_back(e);}
 std::vector<std::array<Window,7>> windows(groups);const auto& ws=spec.at(control?"control_windows":"dvcs_windows").at(period);
 for(int g=0;g<groups;++g)for(auto it=ws[g].at("data").begin();it!=ws[g].at("data").end();++it){auto pos=std::find(VARIABLES.begin(),VARIABLES.end(),it.key());if(pos==VARIABLES.end())throw std::runtime_error("Unknown fixed window");windows[g][pos-VARIABLES.begin()]=window(it.value());}
 size_t stride=size_t(2)*groups*nv*maxbins;std::vector<int64_t> hist(pars.size()*stride),selected(pars.size()*2*groups),unresolved(pars.size()*2*groups);uint64_t records=0;
 auto add=[&](const std::array<double,NC>& v,int scenario,int half,int group,bool global,double mass,bool bad){
  if(group<0||!global)return;std::array<double,7> obs;for(int j=0;j<7;++j)obs[j]=v[16+j];unsigned failed=failures(obs,windows[group]);
  if(!failed){++selected[(scenario*2+half)*groups+group];if(bad)++unresolved[(scenario*2+half)*groups+group];}
  for(int j=0;j<nv;++j){if(j==7?bool(failed):!leave_one_out(failed,j))continue;double value=j==7?mass:obs[j];if(!std::isfinite(value))continue;int bin=std::upper_bound(edges[j].begin(),edges[j].end(),value)-edges[j].begin();++hist[(((size_t(scenario)*2+half)*groups+group)*nv+j)*maxbins+bin];}
 };
 std::ifstream input(job.at("cache").get<std::string>(),std::ios::binary);if(!input)throw std::runtime_error("Missing cache");
 if(control){ControlRecord c;while(input.read(reinterpret_cast<char*>(&c),sizeof c)){
  ++records;int group=pair_group(c);if(group<0)continue;int half=(mix(c.base.key^0xb234fe3adbf120ddULL)>>40)&1;std::array<double,12> noise;for(int j=0;j<12;++j)noise[j]=gaussian(c.base.key,j);auto original=derived(c.base.v,beam,c.base.key&1,MPI0);double zero_mass=c.photons[9]>.5?photon_mass(c):c.photons[8];
  for(size_t s=0;s<pars.size();++s){double mass=c.photons[8];auto v=s?vary_control(c,pars[s],beam,original,noise,mass):c.base.v;
   if(s&&active[s]&&c.photons[9]>.5)mass=c.photons[8]+mass-zero_mass;
   bool gp=selected_global(v,period,cfg,int(c.photons[6]))&&selected_global(v,period,cfg,int(c.photons[7]));add(v,s,half,group,gp,mass,c.photons[9]<.5);
  }
 }}else{Record r;while(input.read(reinterpret_cast<char*>(&r),sizeof r)){
  ++records;if(!std::isfinite(r.v[3])||!std::isfinite(r.v[4]))continue;int d1=int(r.v[3]),d2=int(r.v[4]),group=d1==1&&d2==1?0:d1==2&&d2==1?1:d1==2&&d2==0?2:-1;if(group<0)continue;int half=(mix(r.key^0xb234fe3adbf120ddULL)>>40)&1;auto original=derived(r.v,beam,r.key&1);std::array<double,9> noise;for(int j=0;j<9;++j)noise[j]=gaussian(r.key,j);
  for(size_t s=0;s<pars.size();++s){auto v=s?vary(r,pars[s],beam,r.key&1,original,&noise):r.v;add(v,s,half,group,selected_global(v,period,cfg,d2),NAN,false);}
 }}
 if(!input.eof()||input.gcount())throw std::runtime_error("Truncated cache");
 std::ofstream binary(output,std::ios::binary);binary.write(reinterpret_cast<char*>(hist.data()),hist.size()*sizeof(int64_t));if(!binary)throw std::runtime_error("Histogram write failed");std::ofstream meta(output+".json");meta<<json{{"shape",{int(pars.size()),2,groups,nv,maxbins}},{"selected",selected},{"unresolved",unresolved},{"records",records}};
}
int main(int argc,char**argv){try{if(argc!=5)throw std::runtime_error("selected worker MODE SPEC JOB OUT");auto spec=read_json(argv[2]),job=read_json(argv[3]);std::string mode=argv[1];if(mode=="evaluate-selected")selected_evaluate(spec,job,argv[4]);else if(mode=="export-selected"){std::string sample=job.at("sample");if(sample=="pi0_data"||sample=="pi0_mc")export_unordered_control_arrays(spec,job,argv[4]);else export_cut_arrays(spec,job,argv[4]);}else throw std::runtime_error("Unknown mode");}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
