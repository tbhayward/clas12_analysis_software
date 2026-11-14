#!/usr/bin/env python3
import argparse,copy,json,subprocess,tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import study_selected_particle_response as study
base=study.base


def scoring(tmp):
    variables=[dict(name=n,label=n,bins=80,low=-1,high=1) for n in base.OBS];spec=dict(periods=['fa18_inb'],variables=variables,mass_variable=dict(bins=80,low=-1,high=1))
    x=np.linspace(-1,1,80);good=np.r_[0,100*np.exp(-x*x/.08),0];bad=np.roll(good,8);hist={s:np.tile(good,(2,2,g,nv,1)) for s,g,nv in [('mc',3,7),('background',3,7),('pi0_mc',4,8)]};selected={s:np.ones((2,2,g))*1000 for s,g in [('mc',3),('background',3),('pi0_mc',4)]};data={ch:np.tile(good,(g,nv,1)) for ch,g,nv in [('epgamma',3,7),('epgammagamma',4,8)]};counts={'epgamma':np.ones(3)*1000,'epgammagamma':np.ones(4)*100}
    engine=study.Engine(None,spec,{s:dict(probability=1.) for s in hist},data,counts,10,0);before,_=engine.score_one(hist,selected,1,1.)
    hist['pi0_mc'][1,0,2,7]=bad;after,_=engine.score_one(hist,selected,1,1.);assert after['train']>before['train'] and after['validation_mc']==before['validation_mc']
    hist['pi0_mc'][1,1,2,7]=bad;check,_=engine.score_one(hist,selected,1,1.);assert check['train']==after['train'] and check['validation_mc']>after['validation_mc']
    hist['pi0_mc'][1,0,2,6]=bad;check2,_=engine.score_one(hist,selected,1,1.);assert check2['train']==check['train'] and check2['reserved']>check['reserved']
    p=[[0.]*20,[.01]*20];assert study.select_training(p,[dict(train=2,validation_mc=0),dict(train=1,validation_mc=100)],p)==1
    assert study.select_training(p,[dict(train=2),dict(train=1)],[p[0]])==0
    print('PASS: joint control-mass training, MC-half isolation, reserved-variable exclusion and training-only/allowed-candidate selection')


def native(tmp):
    import test_frozen_response_exports as fixture
    fixture.native(tmp);out=tmp/'selected';out.mkdir();work=tmp/'scratch';work.mkdir();worker=study.compile_worker(out,tmp/'reference',work)
    spec=json.loads((tmp/'spec.json').read_text());spec.update(periods=['fa18_inb'],dvcs_windows={'fa18_inb':[dict(data={'Emiss2':dict(cut_low=-.05,cut_high=.05)}) for _ in range(3)]},control_windows={'fa18_inb':[dict(data={}) for _ in range(4)]},mass_variable=dict(bins=140,low=.10,high=.17));base.dump(tmp/'selected_spec.json',spec)
    control_dtype=np.dtype([('v','<f8',(24,)),('key','<u8'),('photons','<f8',(10,))]);controls=np.fromfile(tmp/'cache.bin',dtype=control_dtype);controls=controls[1:3].copy();controls['photons'][:,6:8]=1;controls['photons'][:,8]=.135;controls['photons'][:,9]=1;controls.tofile(tmp/'ctrl.bin')
    pars=[0.]*20;pars[15]=.05
    job=dict(period='fa18_inb',sample='pi0_mc',cache=str(tmp/'ctrl.bin'),parameters=[[0.]*20,pars]);base.dump(tmp/'selected_job.json',job)
    target=tmp/'control_hist.bin';subprocess.run([str(worker),'evaluate-selected',str(tmp/'selected_spec.json'),str(tmp/'selected_job.json'),str(target)],check=True);m=json.loads(Path(str(target)+'.json').read_text());assert m['shape'][:4]==[2,2,4,8] and sum(m['selected'])==4
    # Calibrated zero-mass preservation, and consistency of the export path.
    for label,vector in [('zero',[0.]*20),('shift',pars)]:
        job.update(parameters=vector,export_individual_photons=True,export_global_only=False);base.dump(tmp/'selected_job.json',job);target=tmp/(label+'.bin');subprocess.run([str(worker),'export-selected',str(tmp/'selected_spec.json'),str(tmp/'selected_job.json'),str(target)],check=True);a=np.fromfile(target,dtype='<f8').reshape(-1,29)
        if label=='zero':assert np.array_equal(a[:,28],np.full(2,.135))
        else:
            original=controls['photons'][0];p1=original[0];p2=original[3];t1,t2=original[1],original[4];df=original[2]-original[5];m0=np.sqrt(4*p1*p2*(np.sin((t1-t2)/2)**2+np.sin(t1)*np.sin(t2)*np.sin(df/2)**2));assert np.allclose(a[:,28],.135+.05*m0,atol=1e-12)
    dt=np.dtype([('v','<f8',(24,)),('key','<u8')]);a=np.zeros(2,dtype=dt);a['v']=controls['v'];a['v'][:,3]=2;a['v'][:,4]=1;a['v'][:,11]=3;a['v'][:,20]=[0,.17];a.tofile(tmp/'signal.bin');job.update(sample='mc',cache=str(tmp/'signal.bin'),parameters=[[0.]*20,pars],export_individual_photons=False);base.dump(tmp/'selected_job.json',job);target=tmp/'signal_hist.bin';subprocess.run([str(worker),'evaluate-selected',str(tmp/'selected_spec.json'),str(tmp/'selected_job.json'),str(target)],check=True);meta=json.loads(Path(str(target)+'.json').read_text());assert sum(meta['selected'][:6])==1 and sum(meta['selected'][6:])==1
    for label,vector in [('old',[0.]*20),('new',pars)]:
        job['parameters']=vector;base.dump(tmp/'selected_job.json',job);target=tmp/(label+'_signal.bin');subprocess.run([str(worker),'export-selected',str(tmp/'selected_spec.json'),str(tmp/'selected_job.json'),str(target)],check=True);z=np.fromfile(target,dtype='<f8').reshape(-1,12);mask=(z[:,5]>=-.05)&(z[:,5]<=.05)
        if label=='old':assert mask.tolist()==[True,False]
        else:assert mask.tolist()==[False,True]
    # Real Python Engine + native jobs: support, normalization, mass scoring and files.
    rng=np.random.default_rng(81);n=400
    ctrl=np.repeat(controls[:1],4*n);sig=np.repeat(a[:1],3*n)
    for g in range(4):
        part=slice(g*n,(g+1)*n);ctrl['v'][part,3]=1 if g==0 else 2;pair=[(1,1),(1,1),(1,0),(0,0)][g];ctrl['photons'][part,6:8]=pair;ctrl['v'][part,4]=pair[0]
    for g in range(3):
        part=slice(g*n,(g+1)*n);sig['v'][part,3]=1 if g==0 else 2;sig['v'][part,4]=0 if g==2 else 1
    for ar in (sig,ctrl):
        ar['v'][:,16:23]=rng.uniform(-.5,.5,(len(ar),7));ar['key']=np.arange(len(ar))*4
    ctrl['photons'][:,8]=rng.uniform(.12,.15,len(ctrl));ctrl.tofile(tmp/'engine_control.bin');sig.tofile(tmp/'engine_signal.bin')
    spec['variables']=[dict(name=name,label=name,bins=80,low=-1,high=1) for name in base.OBS];spec['dvcs_windows']={'fa18_inb':[dict(data={}) for _ in range(3)]};spec['control_windows']={'fa18_inb':[dict(data={}) for _ in range(4)]};base.dump(tmp/'engine_spec.json',spec)
    exported={}
    for sample,path in [('mc',tmp/'engine_signal.bin'),('pi0_mc',tmp/'engine_control.bin')]:
        jb=dict(period='fa18_inb',sample=sample,cache=str(path),parameters=[0.]*20,export_individual_photons=sample=='pi0_mc');base.dump(tmp/'engine_job.json',jb);target=tmp/(sample+'_engine.bin');subprocess.run([str(worker),'export-selected',str(tmp/'engine_spec.json'),str(tmp/'engine_job.json'),str(target)],check=True);exported[sample]=study.arrays_from_export(target)
    raw=dict(mc=exported['mc'],background=exported['mc'],data=np.concatenate([exported['mc'],exported['mc']]),pi0_mc=exported['pi0_mc'],pi0_data=exported['pi0_mc']);windows=dict(dvcs=spec['dvcs_windows']['fa18_inb'],control=spec['control_windows']['fa18_inb']);data,counts=study.data_histograms('fa18_inb',raw,windows,spec);cache={'mc':dict(path=str(tmp/'engine_signal.bin'),probability=1.),'background':dict(path=str(tmp/'engine_signal.bin'),probability=1.),'pi0_mc':dict(path=str(tmp/'engine_control.bin'),probability=1.)}
    pw=tmp/'engine_work';pw.mkdir();pool=base.Pool(worker,tmp/'engine_spec.json',pw,3);engine=study.Engine(pool,spec,cache,data,counts,10,.1);p=[0.]*20;p[15]=.001
    pars,scores,records,hist,selected=engine.evaluate([[0.]*20,p],'integration');assert len(pars)==2 and all(np.isfinite(s['train']) for s in scores) and ('epgammagamma',2,7,0) in engine.support;base.discard_histograms(hist)
    print('PASS: real native evaluator, full pre-cut migration, mass calibration, Python/native Engine integration, control normalization and fixed training support')

def orchestration(tmp):
    from unittest.mock import patch
    run=tmp/'joint';ref=tmp/'reference';run.mkdir();(ref/'backend').mkdir(parents=True);variables=[dict(name=n,label=n,bins=80,low=-1,high=1) for n in base.OBS]
    hashes={}
    for key,name in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]:
        path=ref/'backend'/name;path.write_text('fixture');hashes[key]=base.digest(path)
    spec=dict(periods=['fa18_inb'],variables=variables,input_hashes=hashes,control_to_leakage_exposure_ratio=1.)
    inventory=[]
    for sample in ('data','pi0_data'):
        path=tmp/(sample+'.root');path.write_text('fixture');st=path.stat();inventory.append(dict(period='fa18_inb',sample=sample,path=str(path),size=st.st_size,mtime_ns=st.st_mtime_ns))
    base.dump(ref/'manifest.json',dict(spec=spec,inventory=inventory));base.dump(run/'spec.json',spec);base.dump(run/'manifest.json',dict(arguments=dict(reference_study=str(ref),background_transfer_uncertainty=.1,scan_range='wide'),reference_manifest_sha256=base.digest(ref/'manifest.json'),source_hashes={n:base.digest(base.SOURCE/n) for n in ('particle_response_worker.cpp','particle_response_core.h','particle_control_core.h','particle_control_worker.h')}))
    samples={}
    for sample in ('mc','background','pi0_mc'):
        path=tmp/(sample+'.cache');path.write_bytes(b'fixture');samples['fa18_inb__'+sample]=dict(file=str(path),sha256=base.digest(path),probability=1.)
    base.dump(run/'cache_manifest.json',dict(samples=samples,reference_sha256=base.digest(ref/'manifest.json')));(run/'fa18_inb').mkdir();base.dump(run/'fa18_inb/candidate_parameters.json',dict(parameters=dict.fromkeys(base.PARAMS,0.),bias_mode='fit',fit_topologies=base.TOPS))
    events=[]
    def fake_compile(out,reference,work):
        (out/'backend').mkdir();path=out/'backend/selected_response_worker';path.write_text('fixture worker');return path
    class Pool:
        def __init__(self,worker,spec,output,workers,budget):self.spec=spec;self.output=output
        def run(self,mode,jobs,label):
            assert mode=='export-selected';events.append('exports');result=[];folder=self.output/label;folder.mkdir()
            for i,job in enumerate(jobs):
                control=job['sample'].startswith('pi0');groups=[0,3,4,5] if control else [0,1,2];a=np.zeros((100*len(groups),29 if control else 12));a[:,0]=np.repeat(groups,100);a[:,8]=1;a[:,9]=.2;a[:,10]=1;a[:,11]=np.arange(len(a))%2
                for j in range(1,8):a[:,j]=np.linspace(-.3,.3,len(a))
                if control:a[:,28]=.135
                path=folder/(str(i)+'.bin');a.tofile(path);base.dump(Path(str(path)+'.json'),dict(shape=list(a.shape)));result.append((i,job,path))
            return result
    def fake_cuts(period,a,folder):
        stage='initial' if 'initial_cuts' in str(folder) else 'final';events.append(stage);bound=.5 if stage=='initial' else .2;block=dict(data={'Emiss2':dict(cut_low=-bound,cut_high=bound)},mc={'Emiss2':dict(cut_low=-bound,cut_high=bound)});result=dict(dvcs=[copy.deepcopy(block) for _ in range(3)],control=[copy.deepcopy(block) for _ in range(4)],orders=[['Emiss2']]*3);folder.mkdir(parents=True,exist_ok=True);base.dump(folder/'cuts_95.json',result);return result
    def fake_fit(engine,start,folder,rounds,scan_range):
        events.append('fit');assert engine.spec['dvcs_windows']['fa18_inb'][0]['data']['Emiss2']['cut_high']==.5
        hist={'mc':engine.data['epgamma'][None,None].repeat(2,axis=1),'background':engine.data['epgamma'][None,None].repeat(2,axis=1),'pi0_mc':engine.data['epgammagamma'][None,None].repeat(2,axis=1)};selected={'mc':np.ones((1,2,3))*50,'background':np.ones((1,2,3))*10,'pi0_mc':np.ones((1,2,4))*50}
        base.dump(folder/'candidate_parameters.json',dict(parameters=dict(zip(base.PARAMS,start)),bias_mode='fit',fit_topologies=base.TOPS,background_multiplier=1.));score={ch+'_'+k:0. for ch in ('epgamma','epgammagamma') for k in ('train','validation_mc','reserved')};score['background_multiplier']=1.
        return start,[tuple(start)],[score],[],hist,selected
    out=tmp/'output';args=['--joint-run',str(run),'--output',str(out),'--workers','1','--rounds','1','--scratch-root',str(tmp/'scratch')]
    with patch.object(study,'compile_worker',fake_compile),patch.object(base,'Pool',Pool),patch.object(study,'production_cut_fit',fake_cuts),patch.object(study,'fit_response',fake_fit),patch.object(study,'render_fit',return_value=[]):
        study.main(args);assert events.index('initial')<events.index('fit')<events.index('final');assert (out/'relative_acceptance.csv').is_file() and (out/'coarse_pi0_transfer.csv').is_file() and (out/'window_changes.csv').is_file()
        n=len(events);study.main(args+['--resume']);assert len(events)==n
        base.dump(out/'response/fa18_inb/candidate_parameters.json',dict(parameters={}))
        try:study.main(args+['--resume'])
        except ValueError:pass
        else:raise AssertionError('Changed completed fit was accepted')
    print('PASS: full driver stage ordering, fixed initial cuts, final recount outputs, cache-only checked resume and changed-candidate rejection')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--native',action='store_true');args=p.parse_args()
    with tempfile.TemporaryDirectory() as folder:scoring(Path(folder))
    with tempfile.TemporaryDirectory() as folder:orchestration(Path(folder))
    if args.native:
        with tempfile.TemporaryDirectory() as folder:native(Path(folder))
