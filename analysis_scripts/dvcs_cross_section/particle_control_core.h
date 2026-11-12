#ifndef PARTICLE_CONTROL_CORE_H
#define PARTICLE_CONTROL_CORE_H
#include "particle_response_core.h"
#include <vector>
namespace particle_response {
constexpr double MPI0=.1349766;
// Photon p/theta/phi x2, detector tags x2, measured mgg, reconstruction-valid flag.
struct ControlRecord {Record base;std::array<double,10> photons;};
inline std::vector<V> photon_directions(V q,V e,double azimuth,double opening_deg){
 V z=q*(1/q.norm()),eh=e*(1/e.norm()),x=e-z*e.dot(z);if(x.norm()<1e-12)return {};
 x=x*(1/x.norm());V y=z.cross(x),trans=x*std::cos(azimuth)+y*std::sin(azimuth);
 double a=eh.dot(z),b=eh.dot(trans),radius=std::hypot(a,b),c=std::cos(opening_deg*PI/180);
 if(radius<1e-12||std::abs(c)>radius+1e-5)return {};
 double phi=std::atan2(b,a),d=std::acos(std::clamp(c/radius,-1.,1.));std::vector<V> answer;
 for(double sign:{-1.,1.})for(int turn=-1;turn<=1;++turn){double theta=phi+sign*d+2*PI*turn;if(theta>=0&&theta<=PI)answer.push_back(z*std::cos(theta)+trans*std::sin(theta));}
 return answer;
}
inline bool recover_photons(const Record& r,double beam,double az1,double az2,double open1,double open2,double mass,std::array<double,10>& out,double& error){
 F e=particle(r.v[5],r.v[6],r.v[7],ME),k={{0,0,std::sqrt(beam*beam-ME*ME)},beam};
 V P=particle(r.v[11],r.v[12],r.v[13],MPI0).p;
 auto one=photon_directions((k-e).p,e.p,az1,open1),two=photon_directions((k-e).p,e.p,az2,open2);
 error=1e100;double best_p_error=1e100,best_mass_error=1e100;
 for(V u:one)for(V v:two){double c=u.dot(v),den=1-c*c;if(den<1e-10)continue;
  double E1=(P.dot(u)-c*P.dot(v))/den,E2=(P.dot(v)-c*P.dot(u))/den;if(E1<=0||E2<=0)continue;
  double pe=(u*E1+v*E2-P).norm(),me=std::abs(std::sqrt(std::max(0.,2*E1*E2*(1-c)))-mass);
  double score=pe/std::max(1.,P.norm())+me/.135;
  if(score<error){error=score;best_p_error=pe;best_mass_error=me;out[0]=E1;out[1]=angle(u,{0,0,1});out[2]=wrap(std::atan2(u.y,u.x));out[3]=E2;out[4]=angle(v,{0,0,1});out[5]=wrap(std::atan2(v.y,v.x));}
 }
 out[8]=mass;out[9]=best_p_error<.005&&best_mass_error<.003;return out[9]>0;
}
inline void perturb(std::array<double,NC>& v,int col,int region,const std::array<double,NP>& pars,const std::array<double,12>& z,int first){
 int j=4*region;v[col]*=1+pars[j+3]+pars[j]*z[first];v[col+1]+=pars[j+1]*z[first+1];v[col+2]+=pars[j+2]*z[first+2];
 if(v[col]>0){auto c=particle(v[col],v[col+1],v[col+2],0);v[col+1]=angle(c.p,{0,0,1});v[col+2]=wrap(std::atan2(c.p.y,c.p.x));}
}
inline std::array<double,NC> vary_control(const ControlRecord& c,const std::array<double,NP>& pars,double beam,const std::array<double,12>& original,const std::array<double,12>& noise,double& mass){
 auto v=c.base.v;mass=c.photons[8];bool active=false;for(double x:pars)active=active||x!=0;
 if(!active||c.photons[9]<.5)return v; // unresolved records retained, never dropped
 perturb(v,5,0,pars,noise,0);perturb(v,8,int(v[3])==1?1:2,pars,noise,3);
 auto g1=v,g2=v;for(int j=0;j<3;++j){g1[11+j]=c.photons[j];g2[11+j]=c.photons[3+j];}
 perturb(g1,11,int(c.photons[6])==1?3:4,pars,noise,6);perturb(g2,11,int(c.photons[7])==1?3:4,pars,noise,9);
 if(v[5]<=0||v[8]<=0||g1[11]<=0||g2[11]<=0){v[0]=NAN;return v;}
 F sum=particle(g1[11],g1[12],g1[13],0)+particle(g2[11],g2[12],g2[13],0);mass=std::sqrt(std::max(0.,sum.m2()));
 v[11]=sum.p.norm();v[12]=angle(sum.p,{0,0,1});v[13]=wrap(std::atan2(sum.p.y,sum.p.x));
 auto next=derived(v,beam,c.base.key&1,MPI0);int dest[12]={0,1,2,14,15,16,17,18,19,20,21,22};
 for(int i=0;i<12;++i){double delta=next[i]-original[i];if(i==5)delta=std::remainder(delta,2*PI);v[dest[i]]=c.base.v[dest[i]]+delta;}
 v[16]=(c.base.key&1)?std::abs(std::remainder(v[16],2*PI)):wrap(v[16]);return v;
}
}
#endif
