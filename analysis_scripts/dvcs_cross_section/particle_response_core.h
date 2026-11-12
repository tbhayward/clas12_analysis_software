#ifndef PARTICLE_RESPONSE_CORE_H
#define PARTICLE_RESPONSE_CORE_H
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
namespace particle_response {
constexpr double PI=3.14159265358979323846, MP=.9382720813, ME=.00051099895;
constexpr int NC=24, NP=20;
struct Record {std::array<double,NC> v;uint64_t key;};
inline uint64_t mix(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
inline double uniform(uint64_t x){return ((mix(x)>>11)+.5)*0x1.0p-53;}
inline double gaussian(uint64_t key,int coordinate){return std::sqrt(-2*std::log(uniform(key+2*coordinate)))*std::cos(2*PI*uniform(key+2*coordinate+1));}
inline double wrap(double x){x=std::fmod(x,2*PI);return x<0?x+2*PI:x;}
struct V {
 double x=0,y=0,z=0;
 V operator+(V b)const{return {x+b.x,y+b.y,z+b.z};}
 V operator-(V b)const{return {x-b.x,y-b.y,z-b.z};}
 V operator*(double s)const{return {x*s,y*s,z*s};}
 double dot(V b)const{return x*b.x+y*b.y+z*b.z;}
 V cross(V b)const{return {y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x};}
 double norm()const{return std::sqrt(dot(*this));}
};
struct F {
 V p;double e;
 F operator+(F b)const{return {p+b.p,e+b.e};}
 F operator-(F b)const{return {p-b.p,e-b.e};}
 double m2()const{return e*e-p.dot(p);}
 F boost(V b)const{double b2=b.dot(b);if(b2==0)return *this;if(b2>=1)return {{NAN,NAN,NAN},NAN};double g=1/std::sqrt(1-b2),bp=b.dot(p);return {p+b*((g-1)*bp/b2+g*e),g*(e+bp)};}
};
inline double angle(V a,V b){double n=a.norm()*b.norm();return n>0?std::acos(std::clamp(a.dot(b)/n,-1.,1.)):NAN;}
inline F particle(double p,double t,double f,double m){return {{p*std::sin(t)*std::cos(f),p*std::sin(t)*std::sin(f),p*std::cos(t)},std::sqrt(p*p+m*m)};}
inline double trento(V q,V e,V h){V a=q.cross(e),b=q.cross(h);return wrap(std::atan2(e.cross(h).dot(q)*q.norm(),a.dot(b)));}
// Derived columns: x,Q2,t1,open_angle,pT,then the seven production variables.
inline std::array<double,12> derived(const std::array<double,NC>& v,double beam,bool fallback=false,double neutral_mass=0){
 F k={{0,0,std::sqrt(beam*beam-ME*ME)},beam},target={{0,0,0},MP};
 F e=particle(v[5],v[6],v[7],ME),p=particle(v[8],v[9],v[10],MP),g=particle(v[11],v[12],v[13],neutral_mass);
 F q=k-e,missing=k+target-e-p-g,h=p+g;
 auto pb=p.boost(h.p*(-1/h.e));
 double delta=fallback?std::abs(std::remainder(v[10]-v[13],2*PI)):wrap(trento(q.p,e.p,g.p)-trento(q.p,e.p,p.p));
 double transverse=std::hypot(missing.p.x,missing.p.y);
 return {-q.m2()/(2*MP*q.e),-q.m2(),(target-p).m2(),angle(e.p,g.p)*180/PI,transverse,
         delta,angle(pb.p,h.p),angle((k+target-e-p).p,g.p)*180/PI,transverse,missing.e,missing.m2(),(k+target-e-g).m2()};
}
inline std::array<double,NC> vary(const Record& r,const std::array<double,NP>& pars,double beam,bool fallback,const std::array<double,12>& original,const std::array<double,9>* noise=nullptr){
 auto v=r.v;bool active=false;
 int regions[3]={0,int(v[3])==1?1:2,int(v[4])==1?3:4};
 int start[3]={5,8,11};
 auto z=[&](int i){return noise?(*noise)[i]:gaussian(r.key,i);};
 for(int i=0;i<3;++i){int off=4*regions[i],col=start[i];
  if(pars[off]==0&&pars[off+1]==0&&pars[off+2]==0&&pars[off+3]==0)continue;
  for(int j=0;j<4;++j)active=active||pars[off+j]!=0;
  v[col]*=1+pars[off+3]+pars[off]*z(3*i);
  v[col+1]+=pars[off+1]*z(3*i+1);
  v[col+2]=wrap(v[col+2]+pars[off+2]*z(3*i+2));
  if(v[col]<=0){v[0]=NAN;return v;}
  // Canonical spherical coordinates also handle a smear through the beam axis.
  F canonical=particle(v[col],v[col+1],v[col+2],0);
  v[col+1]=angle(canonical.p,{0,0,1});v[col+2]=wrap(std::atan2(canonical.p.y,canonical.p.x));
 }
 if(!active)return r.v; // exact nominal, including stored rounding/calibrations
 auto shifted=derived(v,beam,fallback);
 int dest[12]={0,1,2,14,15,16,17,18,19,20,21,22};
 for(int i=0;i<12;++i){double delta=shifted[i]-original[i];if(i==5)delta=std::remainder(delta,2*PI);v[dest[i]]=r.v[dest[i]]+delta;}
 v[16]=fallback?std::abs(std::remainder(v[16],2*PI)):wrap(v[16]);
 return v;
}
}
#endif
