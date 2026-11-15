// Persistent native batch evaluator. Python streams ROOT baskets and preserves
// metadata types; this process applies the shared response and all physics updates.
#include "production_response_kinematics.h"
#include "production_response_config.h"
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <cstring>
using namespace production_response;
static json read(const char* p){std::ifstream f(p);if(!f)throw std::runtime_error(std::string("Cannot read ")+p);json j;f>>j;return j;}
static uint64_t pathkey(const std::string& s){uint64_t h=1469598103934665603ULL;for(unsigned char c:s){h^=c;h*=1099511628211ULL;}return h;}
static bool global(const std::array<double,NC>& v,const std::string& period,const GlobalCutConfig& c){return std::isfinite(v[0])&&std::isfinite(v[1])&&std::isfinite(v[2])&&std::isfinite(v[23])&&v[5]>0&&v[8]>0&&v[11]>0&&!is_excluded_run(int(v[23]),c)&&passes_analysis_phase_space_cuts(period,v[0],v[1],v[2],v[5],c)&&passes_proton_momentum_cut(int(v[3]),period,v[8],c)&&passes_global_cuts(v[2],v[14],v[15],int(v[3]),int(v[4]),period,v[5],v[6],v[7],v[9],v[10],v[11],v[12],v[13],c);}
struct Evaluator {
 std::string source,period;bool control,direct,active;double beam,probability;uint64_t seed,source_key=0,offset=0;GlobalCutConfig cfg;
 std::array<double,NP> pars;std::vector<std::array<double,NP>> scenarios;std::vector<std::string> columns;std::map<std::string,size_t> lookup;std::vector<std::pair<int,int>> derived_columns,primary_columns;std::array<int,24> fields;
 Evaluator(const json& j):source(j.at("source")),period(j.at("period")),control(j.at("control")),beam(j.at("beam_energy")),probability(j.at("note_probability")),seed(j.at("seed")),cfg(config(j.at("global_cuts"))),pars(j.at("parameters").get<std::array<double,NP>>()),columns(j.at("columns").get<std::vector<std::string>>()){
  if(std::abs(beam-beam_energy_for_period_label(period))>1e-10)throw std::runtime_error("Response beam energy differs from production period definition");
  for(size_t i=0;i<columns.size();++i)if(!lookup.emplace(columns[i],i).second)throw std::runtime_error("Duplicate response column");
  source_key=pathkey(source)^seed;direct=lookup.count("Delta_phi");active=std::any_of(pars.begin(),pars.end(),[](double x){return x!=0;});
  const std::array<std::string,24> n={"x","Q2","t1","detector1","detector2","e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi","open_angle_ep2","pTmiss","Delta_phi","theta","theta_gamma_gamma","pTmiss","Emiss2","Mx2","Mx2_2","runnum"};
  fields.fill(-1);for(int k=0;k<24;++k){std::string name=control&&k==4?"detector_gamma1":control&&k==18?"theta_pi0_pi0":n[k];auto it=lookup.find(name);if(it!=lookup.end())fields[k]=it->second;else if(k!=16&&k!=23)throw std::runtime_error("Missing response input: "+name);}
  const std::set<std::string> extra={"Mh_gammagamma","gamma_phi1","gamma_phi2","open_angle_egamma1","open_angle_egamma2"};
  for(size_t i=0;i<columns.size();++i){int d=index(columns[i]),primary=-1;for(int k=5;k<=13;++k)if(columns[i]==n[k])primary=k;
   if(primary>=0)primary_columns.emplace_back(i,primary);else if(d>=0)derived_columns.emplace_back(i,d);else if(!metadata(columns[i])&&!extra.count(columns[i]))throw std::runtime_error("Unclassified physics branch: "+columns[i]);}
  if(control)for(auto& name:extra)if(!lookup.count(name))throw std::runtime_error("Missing two-photon metadata: "+name);
  scenarios={std::array<double,NP>{},pars};for(int k=0;k<NP;++k){auto a=pars;a[k]=0;scenarios.push_back(a);}
 }
 void row(double* values,int32_t* flags,std::vector<double>& notes){
  Record r{};for(int k=0;k<24;++k)r.v[k]=fields[k]>=0?values[fields[k]]:0;
  auto get=[&](const std::string& n){return values[lookup.at(n)];};auto set=[&](const std::string& n,double x){values[lookup.at(n)]=x;};
  // Detector assignments are immutable; unsupported tags are retained but excluded.
  bool tags=std::isfinite(r.v[3])&&std::isfinite(r.v[4])&&(r.v[3]==1||r.v[3]==2)&&(r.v[4]==0||r.v[4]==1);
  uint64_t key=mix(source_key^offset++);r.key=(key&~uint64_t(1))|uint64_t(!direct);
  if(!direct)r.v[16]=std::abs(std::remainder(r.v[10]-r.v[13],2*PI));
  ControlRecord c{};c.base=r;std::array<double,12> noise{};
  if(tags){int regions[4]={0,r.v[3]==1?1:2,r.v[4]==1?3:4,control&&get("detector_gamma2")==1?3:4};for(int track=0;track<(control?4:3);++track)for(int coordinate=0;coordinate<3;++coordinate)if(pars[4*regions[track]+coordinate]!=0)noise[3*track+coordinate]=gaussian(r.key,3*track+coordinate);}
  double mass=NAN,original_mass=NAN;bool resolved=true;
  if(control){c.photons[6]=r.v[4];c.photons[7]=get("detector_gamma2");if(!std::isfinite(c.photons[7])||(c.photons[7]!=0&&c.photons[7]!=1))tags=false;double error;resolved=tags&&recover_photons(r,beam,get("gamma_phi1"),get("gamma_phi2"),get("open_angle_egamma1"),get("open_angle_egamma2"),get("Mh_gammagamma"),c.photons,error);c.photons[8]=get("Mh_gammagamma");mass=c.photons[8];if(resolved){auto one=particle(c.photons[0],c.photons[1],c.photons[2],0),two=particle(c.photons[3],c.photons[4],c.photons[5],0);original_mass=std::sqrt(std::max(0.,(one+two).m2()));}}
  auto original=derived(r.v,beam,!direct,control?MPI0:0);
  auto vary_one=[&](const std::array<double,NP>& a,double& m){if(!tags)return r.v;if(control){auto v=vary_control(c,a,beam,original,noise,m);if(resolved&&std::any_of(a.begin(),a.end(),[](double x){return x!=0;}))m=c.photons[8]+m-original_mass;return v;}std::array<double,9> z;std::copy_n(noise.begin(),9,z.begin());return vary(r,a,beam,!direct,original,&z);};
  auto v=vary_one(pars,mass);bool valid=tags&&std::isfinite(v[0])&&std::isfinite(v[1])&&v[0]>0&&v[1]>0&&v[5]>0&&std::sqrt(v[5]*v[5]+ME*ME)<beam&&v[8]>0&&v[11]>0;bool unresolved=control&&!resolved;bool applied=tags&&active&&!unresolved;flags[0]=applied;flags[1]=valid;flags[2]=unresolved;
  if(uniform(key^0xeca7d931ULL)<probability){for(auto& a:scenarios){double mm=control?c.photons[8]:NAN;auto w=vary_one(a,mm);std::array<double,31> z{};int d1=tags?int(w[3]):-1,d2=tags?int(w[4]):-1;z[0]=d1==1&&d2==1?0:d1==2&&d2==1?1:d1==2&&d2==0?2:-1;for(int k=0;k<7;++k)z[k+1]=w[16+k];z[8]=tags&&global(w,period,cfg);z[9]=std::abs(w[2]);z[10]=(mix(r.key^0xb234fe3adbf120ddULL)>>40)&1;z[11]=control&&tags?(d1==2?3:0)+2-c.photons[6]-c.photons[7]:-1;z[12]=mm;z[13]=resolved;for(int j=0;j<16;++j)z[14+j]=w[j];z[30]=w[23];notes.insert(notes.end(),z.begin(),z.end());}}
  if(applied){auto before=calculate(r.v,beam,!direct,control,&original),after=calculate(v,beam,!direct,control);
   for(auto [i,k]:primary_columns)values[i]=v[k];
   for(auto [i,k]:derived_columns){double delta=after[k]-before[k];if(circular(columns[i]))delta=std::remainder(delta,2*PI);values[i]+=delta;if(circular(columns[i]))values[i]=columns[i]=="Delta_phi"&&!direct?std::abs(std::remainder(values[i],2*PI)):wrap(values[i]);}
   for(int k=0;k<24;++k)if((k<3||k>=14&&k<=22)&&fields[k]>=0)values[fields[k]]=v[k];
   if(control&&resolved){set("Mh_gammagamma",mass);for(int photon=0;photon<2;++photon){auto a=r.v;for(int k=0;k<3;++k)a[11+k]=c.photons[3*photon+k];auto oldg=particle(a[11],a[12],a[13],0),olde=particle(r.v[5],r.v[6],r.v[7],ME);F k={{0,0,std::sqrt(beam*beam-ME*ME)},beam};perturb(a,11,int(c.photons[6+photon])==1?3:4,pars,noise,6+3*photon);auto newg=particle(a[11],a[12],a[13],0),newe=particle(v[5],v[6],v[7],ME);std::string az="gamma_phi"+std::to_string(photon+1),op="open_angle_egamma"+std::to_string(photon+1);set(az,wrap(get(az)+std::remainder(trento((k-newe).p,newe.p,newg.p)-trento((k-olde).p,olde.p,oldg.p),2*PI)));set(op,get(op)+(angle(newe.p,newg.p)-angle(olde.p,oldg.p))*180/PI);}}
  }
  if(!valid&&active){set("x",NAN);set("Q2",NAN);}
 }
};
int main(int argc,char** argv){try{
 if(argc!=2)throw std::runtime_error("production_response_writer JOB.json");std::ios::sync_with_stdio(false);std::cin.tie(nullptr);Evaluator evaluator(read(argv[1]));
 while(true){uint64_t n=0;std::cin.read(reinterpret_cast<char*>(&n),sizeof n);if(std::cin.eof()&&!std::cin.gcount())break;if(!std::cin||n>65536)throw std::runtime_error("Invalid response batch");size_t size=n*evaluator.columns.size();std::vector<double> values(size),notes;std::vector<int32_t> flags(n*3);std::cin.read(reinterpret_cast<char*>(values.data()),size*sizeof(double));if(!std::cin)throw std::runtime_error("Truncated response batch");
 for(size_t i=0;i<n;++i)evaluator.row(values.data()+i*evaluator.columns.size(),flags.data()+3*i,notes);
 uint64_t note_rows=notes.size()/(22*31);std::cout.write(reinterpret_cast<char*>(values.data()),size*sizeof(double));std::cout.write(reinterpret_cast<char*>(flags.data()),flags.size()*sizeof(int32_t));std::cout.write(reinterpret_cast<char*>(&note_rows),sizeof note_rows);std::cout.write(reinterpret_cast<char*>(notes.data()),notes.size()*sizeof(double));std::cout.flush();if(!std::cout)throw std::runtime_error("Response pipe write failed");}
 return 0;
}catch(const std::exception& e){std::cerr<<"[production_response_writer] "<<e.what()<<'\n';return 1;}}
