#!/usr/bin/env python3
import tempfile,subprocess
from pathlib import Path
import numpy as np
import study_particle_response as base
import study_photon_order_audit as audit


def tables():
    a=np.zeros((5,28));a[:,0]=[1,2,0,1,1];a[:,8]=1;a[:,20]=[3,1,1,2,np.nan];a[:,24]=[1,3,3,2,np.nan];a[:,23]=[1,0,1,1,1];a[:,27]=[0,1,0,1,0]
    rr,mm=audit.summarize(a,'test','DATA','global',np.ones(5,dtype=bool))
    allrow=next(r for r in rr if r['photon_pair']=='all resolved');assert allrow['events']==4 and allrow['unresolved_events']==1 and allrow['equal_energy']==1
    assert next(r for r in mm if r['saved_topology']=='CD_FT' and r['harder_photon_topology']=='CD_FD')['events']==1
    assert next(r for r in mm if r['saved_topology']=='FD_FD' and r['harder_photon_topology'].startswith('FD_FT'))['events']==1
    assert sum(r['events'] for r in mm)==3 # unresolved and energy ties explicitly excluded
    print('PASS: mixed-detector grouping, unsupported FD–FT group, unresolved coverage and tie exclusion')


def native(tmp):
    cpp=tmp/'permutation.cpp';cpp.write_text(r'''
#include "particle_control_core.h"
#include <cassert>
using namespace particle_response;
int main(){Record r{};r.v[3]=2;r.v[4]=1;r.v[5]=5;r.v[6]=.2;r.v[7]=.3;r.v[8]=.7;r.v[9]=.8;r.v[11]=5;r.v[12]=.12;r.v[13]=1.3;
F pi=particle(5,.12,1.3,MPI0),e=particle(5,.2,.3,ME),k={{0,0,std::sqrt(10.604*10.604-ME*ME)},10.604};
F a=particle(MPI0/2,std::acos(.5),1.,0),b={a.p*(-1),a.e};a=a.boost(pi.p*(1/pi.e));b=b.boost(pi.p*(1/pi.e));double f1=trento((k-e).p,e.p,a.p),f2=trento((k-e).p,e.p,b.p),o1=angle(e.p,a.p)*180/PI,o2=angle(e.p,b.p)*180/PI;
ControlRecord c{},sw{};c.base=sw.base=r;double error;
assert(recover_photons(r,10.604,f1,f2,o1,o2,MPI0,c.photons,error));assert(recover_photons(r,10.604,f2,f1,o2,o1,MPI0,sw.photons,error));
for(int j=0;j<3;++j){assert(std::abs(c.photons[j]-sw.photons[3+j])<1e-7);assert(std::abs(c.photons[3+j]-sw.photons[j])<1e-7);}assert(std::abs(c.photons[0]-a.e)<1e-7);
c.photons[6]=1;c.photons[7]=0;sw.photons[6]=0;sw.photons[7]=1;std::array<double,NP> pars{};pars[15]=.01;pars[19]=-.02;std::array<double,12> z{};double m1,m2;auto original=derived(r.v,10.604,false,MPI0);auto v=vary_control(c,pars,10.604,original,z,m1),w=vary_control(sw,pars,10.604,original,z,m2);for(int j=11;j<14;++j)assert(std::abs(v[j]-w[j])<1e-7);assert(std::abs(m1-m2)<1e-7);
// Swapping only tags must be detectable for detector-dependent responses.
sw.photons[6]=1;sw.photons[7]=0;auto wrong=vary_control(sw,pars,10.604,original,z,m2);assert(std::abs(wrong[11]-v[11])>.001);
}
''');exe=tmp/'permutation';subprocess.run(['c++','-std=c++17','-O2','-I'+str(base.SOURCE),str(cpp),'-o',str(exe)],check=True);subprocess.run([str(exe)],check=True)
    print('PASS: recovered coordinates preserve input photon roles; swapping coordinates and tags preserves response; swapping tags alone does not')

if __name__=='__main__':
    tables()
    with tempfile.TemporaryDirectory() as folder:native(Path(folder))
