#!/usr/bin/env python3
"""Domain, volume, centering convention and plotting-reference regressions.
No ROOT data needed. Optional --real-km15 checks the installed physics model.
"""
import os
for n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[n]='1'
import argparse,csv,importlib,math,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
SOURCE=Path(__file__).resolve().parent
sys.path.insert(0,str(SOURCE/'external_scripts'))
from kinematic_domain import Domain,minus_t_min
import km15_bin_center_batch as bc

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--real-km15',action='store_true');args=ap.parse_args()
    domain=Domain();interior=(.155,.204,1.912,2.510,.11,.15);boundary=(.118,.155,1.912,2.510,.11,.15)
    assert abs(domain.fraction(interior,10.6)-1)<1e-9
    assert .51<domain.fraction(boundary,10.1998)<.52
    impossible=(.6,.7,1.,1.2,.11,.15);assert domain.fraction(impossible,10.6)==0
    for box in (interior,boundary,(.155,.204,2.510,3.295,.25,.4)):
        pts,w=domain.quadrature(box,(345.,360.),10.1998,5)
        assert all(domain.contains(p) for p in pts)
        assert all(w>0)
    # One standalone C++ implementation must agree with Python, including custom W cuts.
    with tempfile.TemporaryDirectory() as tmp:
        t=Path(tmp);cpp=t/'check.cpp';exe=t/'check'
        cpp.write_text('#include "kinematic_bin_integration.h"\n#include <iostream>\n#include <iomanip>\nint main(){dvcs_phase_space::Box b;double e,w;int a;std::cout<<std::setprecision(17);while(std::cin>>b.x0>>b.x1>>b.q0>>b.q1>>b.t0>>b.t1>>e>>a>>w)std::cout<<dvcs_phase_space::allowed_fraction(b,e,a,w)<<"\\n";}')
        subprocess.run(['c++','-O2','-std=c++17','-I'+str(SOURCE),str(cpp),'-o',str(exe)],check=True)
        cases=[(box,e,apply,w) for box in (interior,boundary,impossible,(.155,.204,2.510,3.295,.25,.4))
               for e in (10.6,10.1998) for apply,w in ((True,2.),(True,2.2),(False,2.))]
        inp='\n'.join(' '.join(map(str,(*b,e,int(a),w))) for b,e,a,w in cases)+'\n'
        vals=np.array([float(v) for v in subprocess.check_output([str(exe)],input=inp,text=True).splitlines()])
        expected=np.array([Domain(a,w).fraction(b,e) for b,e,a,w in cases]);assert np.allclose(vals,expected,atol=1e-8,rtol=1e-7)
    # Constant model: clipping the bin must not manufacture a centering correction.
    saved=bc._km15
    bc._km15=lambda p:2.
    center=(.14,2.,.13,352.5,10.1998)
    assert domain.contains(center)
    v=domain.volume(boundary,(345.,360.),10.1998)[0]
    task=(0,'10.2 GeV',center,boundary,(345.,360.),10.1998,3,11,.002,domain,v)
    result=bc.integrate_task(task);assert result['status']=='converged';assert abs(result['factor']-1)<1e-12
    bc._km15=saved
    # Lee denominator convention: 0.8/0.9 * 2/1.8 = 80/81.
    original={'Fbin, 10.6 GeV':'(.8,0,0)','bin_volume, 10.6 GeV':'(2,0,0)',
              'cross sections, ep->epg, exp, 10.6 GeV, unpol':'(10,2,0)'}
    updated=dict(original);updated['Fbin, 10.6 GeV']='(.9,0,0)';updated['bin_volume, 10.6 GeV']='(1.8,0,0)'
    bc.rescale_cross_sections([original],[updated]);value,stat,_=bc.triple(updated['cross sections, ep->epg, exp, 10.6 GeV, unpol'])
    assert abs(value-10*80/81)<1e-10 and abs(stat-2*80/81)<1e-10
    assert original['cross sections, ep->epg, exp, 10.6 GeV, unpol']=='(10,2,0)'
    missing=dict(updated);missing['Fbin, 10.6 GeV']=''
    bc.rescale_cross_sections([original],[missing])
    assert missing['cross sections, ep->epg, exp, 10.6 GeV, unpol']==''
    # Nonconstant model catches a reciprocal factor error; constant F=1 cannot.
    bc._km15=lambda p: 1.+p[3]/360.
    asymmetric=(.18,2.,.13,3.,10.6)
    box=(.155,.204,1.912,2.510,.11,.15)
    vol=domain.volume(box,(0.,30.),10.6)[0]
    result=bc.integrate_task((0,'10.6 GeV',asymmetric,box,(0.,30.),10.6,3,11,.002,domain,vol))
    assert result['status']=='converged',result
    assert abs(result['factor']-(1+15/360)/(1+3/360))<1e-10
    bc._km15=saved
    # Check the exact arithmetic used by the production C++ calculation.
    with tempfile.TemporaryDirectory() as tmp:
        t=Path(tmp);cpp=t/'factor.cpp';exe=t/'factor'
        cpp.write_text('#include "correction_convention.h"\n#include <cassert>\n#include <cmath>\nint main(){assert(std::abs(dvcs_corrections::cross_section(80.,10.,.8,.5)-20.)<1e-12);}')
        subprocess.run(['c++','-std=c++17','-I'+str(SOURCE),str(cpp),'-o',str(exe)],check=True)
        subprocess.run([str(exe)],check=True)
    legacy={'Frad, 10.6 GeV':'(.8,0,0)','Fbin, 10.6 GeV':'(.5,0,0)',
            'cross sections, ep->epg, exp, 10.6 GeV, unpol':'(8,2,0)'}
    try:bc.convert_input_convention([legacy])
    except ValueError:pass
    else:raise AssertionError('Unmarked CSV accepted without a declared convention')
    fixed,convention=bc.convert_input_convention([legacy],'legacy-multiply')
    assert abs(bc.triple(fixed[0]['cross sections, ep->epg, exp, 10.6 GeV, unpol'])[0]-50)<1e-12
    again,_=bc.convert_input_convention(fixed);assert again==fixed
    assert legacy['cross sections, ep->epg, exp, 10.6 GeV, unpol']=='(8,2,0)'
    # Import plotting dependencies only for the actual reference-coordinate check.
    plot=importlib.import_module('plot_pass2_final_results_by_t');base=plot.base
    assert plot.make_parser().parse_args(['input.csv']).reference_kinematics=='imported'
    assert plot.make_parser().parse_args(['input.csv','--reference-kinematics','midpoints']).reference_kinematics=='midpoints'
    config=base.ModelConfig(10.6,'unused','unused',sys.executable,24,False,False,False)
    refs={'xs':((.2,2.,.15,7.,10.6),(.25,2.2,.17,23.,10.6)),
          'bsa':((.2,2.,.15,7.,10.6),)}
    oldmany,oldbh=base.km15_xs_many,base.bh_xs
    base.km15_xs_many=lambda points,cfg:[x+q+t+phi/100+h/10 for x,q,t,phi,h in points]
    base.bh_xs=lambda x,q,t,phi,cfg:x+q+t+phi/50
    _,curves=plot.compute_cell_models(('key',config,('xs','bsa'),refs))
    assert curves['xs']['phi']==[7.,23.]
    assert abs(curves['xs']['km15'][1]-(.25+2.2+.17+.23))<1e-12
    assert plot.interpolate_model(curves['xs'],'km15',23.)==curves['xs']['km15'][1]
    assert math.isnan(plot.interpolate_model(curves['xs'],'km15',15.))
    assert math.isfinite(plot.interpolate_model(curves['bsa'],'km15',7.))
    base.km15_xs_many,base.bh_xs=oldmany,oldbh
    if args.real_km15:
        bc._init_worker();result=bc.integrate_task(task)
        assert result['status']=='converged',result
        print('Real KM15 boundary example:',result)
        cached=[bc._km15((.14,2.,.13,phi,10.1998)) for phi in (1.,7.,19.,177.,355.)]
        os.environ['KM15_DISABLE_CFF_CACHE']='1';bc._init_worker()
        reference=[bc._km15((.14,2.,.13,phi,10.1998)) for phi in (1.,7.,19.,177.,355.)]
        os.environ.pop('KM15_DISABLE_CFF_CACHE')
        assert cached==reference, 'Cached/scalar CFF predictions differ'
    # End-to-end failure policy: strict refuses a combined CSV; preview omits
    # only the failed energy while preserving the original and volume control.
    import run_numerical_corrections as runner
    import json
    saved_recalculate=bc.recalculate;saved_argv=sys.argv[:]
    with tempfile.TemporaryDirectory() as tmp:
        t=Path(tmp);input_path=t/'input.csv';config_path=t/'cuts.json'
        config_path.write_text(json.dumps(dict(enable_analysis_phase_space_cuts=True,analysis_W_min_GeV=2.)))
        record=dict(zip(('xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max'),map(str,boundary)))
        record.update({'bin index':'1','phimin':'0','phimax':'15'})
        for label in ('10.6 GeV','10.2 GeV'):
            record['Frad, '+label]='(.9,0,0)';record['Fbin, '+label]='(.8,0,0)';record['bin_volume, '+label]='(2,0,0)'
            record['cubic bin_volume, '+label]='(2,0,0)'
        for label in ('10.6 GeV','Sp19 Inb'):
            for prefix,value in zip(('xBavg','Q2avg','t_abs_avg','phiavg'),(.14,2.,.13,7.)):
                record[prefix+', '+label]=str(value)
            record[f'cross sections, ep->epg, exp, {label}, unpol']='(10,2,0)'
        with input_path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(record));w.writeheader();w.writerow(record)
        before=input_path.read_bytes()
        bc.recalculate=lambda rows,*a,**k:([dict(r) for r in rows],[dict(row=0,label='10.2 GeV',status='unconverged',order=11)])
        for policy in ('error','omit'):
            destination=t/policy
            sys.argv=['runner','--input-csv',str(input_path),'--output-dir',str(destination),
                      '--cuts-json',str(config_path),'--on-unconverged',policy,'--input-correction-convention','legacy-multiply']
            try:runner.main()
            except SystemExit:
                assert policy=='error'
            assert (destination/'volumes_only.csv').is_file()
            assert (destination/'original.csv').read_bytes()==before
            if policy=='error':assert not (destination/'new_fbin_and_volumes.csv').exists()
            else:
                row=next(csv.DictReader((destination/'new_fbin_and_volumes.csv').open()))
                assert row['Fbin, 10.2 GeV']==''
                assert row['cross sections, ep->epg, exp, Sp19 Inb, unpol']==''
                assert row['cross sections, ep->epg, exp, 10.6 GeV, unpol']!=''
        assert input_path.read_bytes()==before
        # Fast mode must not call the numerical integrator, even once.
        bc.recalculate=lambda *a,**k: (_ for _ in ()).throw(AssertionError('Unexpected integration'))
        sys.argv=['runner','--input-csv',str(input_path),'--output-dir',str(t/'fast'),
                  '--cuts-json',str(config_path),'--input-correction-convention','legacy-multiply','--convention-only']
        assert runner.main()==0
        assert (t/'fast'/'original.csv').read_bytes()==before
        repaired=next(csv.DictReader((t/'fast'/'convention_fixed.csv').open()))
        assert repaired['cross section correction convention']=='pass1-divide'
        assert abs(bc.triple(repaired['cross sections, ep->epg, exp, 10.6 GeV, unpol'])[0]-10/(.9*.8)**2)<1e-9
    bc.recalculate=saved_recalculate;sys.argv=saved_argv
    print('PASS: allowed nodes, constant-model normalization, C++/Python volumes, reciprocal convention, frozen-input rescaling, exact mean-reference predictions')

if __name__=='__main__':main()
