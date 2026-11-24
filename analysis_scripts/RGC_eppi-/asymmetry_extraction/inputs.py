"""Strict readers for nominal selected files and charge/polarization states."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
PERIODS=['RGB_spring2019_D2','RGB_fall2019_D2','RGB_spring2020_D2',
         'RGC_summer2022_ND3','RGC_fall2022_ND3','RGC_spring2023_ND3']


def beam_pol(run):
    # Identical to the campaign beam_pol.h, including its 0.83 fallback.
    intervals=[(6616,6783,.8453),(6142,6149,.81132),(6150,6188,.82137),
     (6189,6260,.83598),(6261,6339,.80770),(6340,6342,.85536),(6344,6399,.87038),
     (6420,6476,.88214),(6479,6532,.86580),(6533,6603,.87887),
     (11013,11309,.84983),(11323,11334,.87135),(11335,11387,.85048),(11389,11571,.84262),
     (16658,16695,.74052),(16043,16771,.8384),(16843,17407,.8372),(17477,17810,.8040)]
    return next((v for lo,hi,v in intervals if lo<=run<=hi),.83)


def read_rows(path):
    rows={}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):continue
        fields=next(csv.reader([line]))
        try:run=int(fields[0])
        except ValueError:continue
        if len(fields) not in [5,6,7]:raise ValueError(f'Run {run}: expected 5, 6 or 7 fields')
        values=list(map(float,fields[1:]))
        if not np.all(np.isfinite(values)):raise ValueError(f'Run {run}: nonfinite metadata')
        if run in rows and rows[run]!=values:raise ValueError(f'Conflicting metadata for run {run}')
        rows[run]=values
    return rows


def states_for(period,runs,rows):
    records={};states=[]
    for run in runs:
        if run not in rows:raise ValueError(f'{period}: missing charge metadata for run {run}')
        v=rows[run]
        if period.startswith('RGB'):
            if len(v)!=6 or any(x!=0 for x in v[1:]):
                raise ValueError(f'{period} run {run}: expected seven-column RGB total-charge row with zero measured helicity charges and target polarization')
            qp=qm=v[0]/2;pt=0.
        else:
            if len(v)!=5:raise ValueError(f'{period} run {run}: expected six-column RGC row')
            qp,qm,pt=v[1:4]
        if qp<0 or qm<0 or qp+qm<=0:raise ValueError(f'{period} run {run}: nonpositive source-run exposure')
        if abs(pt)>1 or (period.startswith('RGC') and pt==0):raise ValueError(f'{period} run {run}: invalid target polarization')
        pb=beam_pol(run);records[run]=dict(q_plus=qp,q_minus=qm,pt=pt,pb=pb)
        for h,q in [(1,qp),(-1,qm)]:
            if q>0:states.append([q,h*pb,pt,h*pb*pt])
    return np.array(states),records


def load(registry_path,metadata_path,dilution_path,selection_code):
    import sys
    sys.path.insert(0,str(selection_code))
    from selection import load_selected_inputs
    import uproot
    registry=load_selected_inputs(registry_path)
    dilution=json.loads(Path(dilution_path).read_text())
    if dilution['cutset_id']!=registry['cutset_id']:raise ValueError('Dilution and selected data have different cut sets')
    df={int(r['bin']):(float(r['df']),float(r['stat_std'])) for r in dilution['bins']}
    if set(df)!=set(range(1,25)) or any(not 0<f<1 or not 0<s<1 for f,s in df.values()):raise ValueError('Need 24 finite physical dilution factors with symmetric statistical uncertainties')
    rows=read_rows(metadata_path);groups={};audit={}
    for name in PERIODS:
        item=registry['files'][name]
        with uproot.open(item['path']) as f:
            report=json.loads(str(f['selection_manifest']))
            if report['partial'] or report['read_rows']!=report['source_entries']:raise ValueError(name+': incomplete export')
            arrays=f[registry['tree']].arrays(['runnum','evnum','analysis_bin','phi2','helicity','beam_pol','target_pol','DepA','DepB','DepC','DepV','DepW'],library='np')
        states,records=states_for(name,sorted(map(int,report['runs'])),rows)
        if not all(np.all(np.isfinite(a)) for a in arrays.values()):raise ValueError(name+': nonfinite event input')
        if np.any(arrays['DepA']<=0):raise ValueError(name+': nonpositive DepA')
        if np.any(~np.isin(arrays['helicity'],[-1,1])):raise ValueError(name+': invalid helicity')
        if np.any((arrays['phi2']<0)|(arrays['phi2']>2*np.pi+1e-8)):raise ValueError(name+': phi2 is not radians')
        for run in np.unique(arrays['runnum']):
            m=arrays['runnum']==run;r=records[int(run)]
            if not np.allclose(arrays['beam_pol'][m],r['pb'],atol=1e-8,rtol=0):raise ValueError(f'{name}: beam polarization mismatch for {run}')
            if not np.allclose(arrays['target_pol'][m],r['pt'],atol=5e-4,rtol=0):raise ValueError(f'{name}: target polarization mismatch for {run}')
            if (r['q_plus']==0 and np.any(arrays['helicity'][m]>0)) or (r['q_minus']==0 and np.any(arrays['helicity'][m]<0)):raise ValueError(f'{name}: event in zero-exposure state')
        events=dict(runnum=arrays['runnum'],evnum=arrays['evnum'],phi=arrays['phi2'],h=arrays['helicity'],pb=arrays['beam_pol'],pt=np.array([records[int(r)]['pt'] for r in arrays['runnum']]),bin=arrays['analysis_bin'])
        for suffix,branch in [('b','B'),('c','C'),('v','V'),('w','W')]:events['r'+suffix]=arrays['Dep'+branch]/arrays['DepA']
        groups[name]=dict(events=events,states=states)
        pairs=np.rec.fromarrays([arrays['runnum'],arrays['evnum']])
        audit[name]=dict(source_runs=len(records),source_run_numbers=sorted(records),rows=len(pairs),unique_events=len(np.unique(pairs)),
          charge=sum(r['q_plus']+r['q_minus'] for r in records.values()),
          charge_assumption='Qplus=Qminus=total/2 (user-authorized RGB assumption)' if name.startswith('RGB') else 'Measured helicity charges, RGC raw sign equals stored analysis sign',
          bin_counts=np.bincount(arrays['analysis_bin'],minlength=25)[1:25].tolist())
    identity=dict(cutset_id=registry['cutset_id'],metadata_sha256=hashlib.sha256(Path(metadata_path).read_bytes()).hexdigest(),
      dilution_sha256=hashlib.sha256(Path(dilution_path).read_bytes()).hexdigest(),
      files={n:dict(path=registry['files'][n]['path'],size=Path(registry['files'][n]['path']).stat().st_size,mtime=Path(registry['files'][n]['path']).stat().st_mtime_ns) for n in PERIODS})
    return groups,df,audit,identity
