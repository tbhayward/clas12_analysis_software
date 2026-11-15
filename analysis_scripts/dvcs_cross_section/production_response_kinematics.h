#ifndef PRODUCTION_RESPONSE_KINEMATICS_H
#define PRODUCTION_RESPONSE_KINEMATICS_H
// Production extension of the fitted response: preserve saved calibrations by
// adding the change calculated from the original/modified four-vectors.
#include "particle_control_core.h"
#include <string>
#include <vector>
#include <unordered_map>
namespace production_response {
using namespace particle_response;
inline const std::vector<std::string>& names(){
 static const std::vector<std::string> n={"x","Q2","W","y","nu","gamma","t","t1","t2","tmin","open_angle_ep","open_angle_ep1","open_angle_ep2","open_angle_p1p2","Delta_phi","theta","theta_gamma_gamma","theta_pi0_pi0","pTmiss","Emiss2","Mx2","Mx2_1","Mx2_2","Mh","z","z1","z2","xF","xF1","xF2","pT","pT1","pT2","pTpT","xi","xi1","xi2","zeta","zeta1","zeta2","eta","eta1","eta2","eta_gN","eta1_gN","eta2_gN","Delta_eta","phi1","phi2","phih","phiR","DepA","DepB","DepC","DepV","DepW","e_px","e_py","e_pz","e_e","p1_px","p1_py","p1_pz","p1_e","p2_px","p2_py","p2_pz","p2_e","p_px","p_py","p_pz","p_p","p_e","p_theta","p_phi"};return n;
}
inline int index(const std::string& name){static const auto lookup=[](){std::unordered_map<std::string,int> x;int i=0;for(auto& n:names())x[n]=i++;return x;}();auto it=lookup.find(name);return it==lookup.end()?-1:it->second;}
inline double signed_mass(double m2){return m2>=0?std::sqrt(m2):-std::sqrt(-m2);}
inline double rapidity(F a,V axis){double z=a.p.dot(axis)/axis.norm();return .5*std::log((a.e+z)/(a.e-z));}
inline bool circular(const std::string& name){return name=="phi1"||name=="phi2"||name=="phih"||name=="phiR"||name=="p_phi"||name=="Delta_phi";}
inline std::vector<double> calculate(const std::array<double,NC>& v,double beam,bool fallback,bool control,const std::array<double,12>* known=nullptr){
 F k={{0,0,std::sqrt(beam*beam-ME*ME)},beam},target={{0,0,0},MP};
 F e=particle(v[5],v[6],v[7],ME),p=particle(v[8],v[9],v[10],MP),g=particle(v[11],v[12],v[13],control?MPI0:0),h=p+g,q=k-e,miss=k+target-e-h;
 double Q2=-q.m2(),nu=q.e,x=Q2/(2*MP*nu),W=std::sqrt((q+target).m2()),y=nu/beam,gamma=2*MP*x/std::sqrt(Q2),eps=gamma*gamma;
 double tmin=-Q2*(2*(1-x)*(1-std::sqrt(1+eps))+eps)/(4*x*(1-x)+eps);
 auto d=known?*known:derived(v,beam,fallback,control?MPI0:0);V boost=(q+target).p*(-1/(q+target).e);F pg=p.boost(boost),gg=g.boost(boost),hg=h.boost(boost),qg=q.boost(boost),tg=target.boost(boost),eg=e.boost(boost);
 F breit=q+F{{0,0,0},2*x*MP};V bb=breit.p*(-1/breit.e);F qb=q.boost(bb),pb=p.boost(bb),gb=g.boost(bb),hb=h.boost(bb);
 auto pt=[&](F a){return qg.p.cross(a.p).norm()/qg.p.norm();};auto xf=[&](F a){return 2*a.p.dot(qg.p)/(qg.p.norm()*W);};auto inner=[](F a,F b){return a.e*b.e-a.p.dot(b.p);};double denominator=inner(target,q);
 V rel=pg.p*(nu/p.e)-gg.p*(nu/g.e);double phi1=trento(q.p,e.p,p.p),phi2=trento(q.p,e.p,g.p),phih=trento(q.p,e.p,h.p),phiR=trento(qg.p,eg.p,rel);
 double a=(1-y+y*y/2+y*y*eps/4)/(1+eps),b=(1-y-y*y*eps/4)/(1+eps),c=y*(1-y/2)/std::sqrt(1+eps),root=std::sqrt(1-y-y*y*eps/4),vv=(2-y)*root/(1+eps),ww=y*root/std::sqrt(1+eps);
 std::vector<double> z={x,Q2,W,y,nu,gamma,(target-h).m2(),(target-p).m2(),(target-g).m2(),tmin,angle(e.p,h.p)*180/PI,angle(e.p,p.p)*180/PI,d[3],angle(p.p,g.p)*180/PI,d[5],d[6],d[7],d[7],d[8],d[9],d[10],(k+target-e-p).m2(),d[11],signed_mass(h.m2()),h.e/nu,p.e/nu,g.e/nu,xf(hg),xf(pg),xf(gg),pt(hg),pt(pg),pt(gg),pt(pg)*pt(gg),inner(h,q)/denominator,inner(p,q)/denominator,inner(g,q)/denominator,hg.e/tg.e,pg.e/tg.e,gg.e/tg.e,-rapidity(hb,qb.p),-rapidity(pb,qb.p),-rapidity(gb,qb.p),rapidity(hg,qg.p),rapidity(pg,qg.p),rapidity(gg,qg.p),rapidity(pb,qb.p)-rapidity(gb,qb.p),phi1,phi2,phih,phiR,a,b,c,vv,ww,e.p.x,e.p.y,e.p.z,e.e,p.p.x,p.p.y,p.p.z,p.e,g.p.x,g.p.y,g.p.z,g.e,h.p.x,h.p.y,h.p.z,h.p.norm(),h.e,angle(h.p,{0,0,1}),wrap(std::atan2(h.p.y,h.p.x))};
 if(z.size()!=names().size())throw std::logic_error("Production response column inventory mismatch");return z;
}
inline bool metadata(const std::string& n){return n=="runnum"||n=="evnum"||n=="event_sequence"||n=="helicity"||n=="weight"||n=="fiducial_status"||n.rfind("detector",0)==0||n.rfind("num_",0)==0||n.rfind("vz_",0)==0||n.rfind("vx_",0)==0||n.rfind("vy_",0)==0||n.rfind("gen_",0)==0||n.rfind("mc_",0)==0||n=="rec_kinematics_valid";}
}
#endif
