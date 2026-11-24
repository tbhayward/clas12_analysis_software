"""Stream selected TTrees with all input branches and auditable cut provenance."""
import csv
import json
from pathlib import Path
import time
import numpy as np
from candidate_cache import save_json
from selection import ChannelSelection
from kinematics import bin_indices

DOWNSTREAM=['runnum','evnum','helicity','beam_pol','target_pol','phi2','DepA','DepB','DepC','DepV','DepW']

EXTRA_TYPES=dict(analysis_bin=np.dtype('int32'),exclusivity_group=np.dtype('int32'),Mx2_deuteron=np.dtype('float64'),p_miss=np.dtype('float64'),delta_theta_deg=np.dtype('float64'),exclusivity_reference_assisted=np.dtype('int32'))

def audit_file(path,tree_name,expected,original_types,cutset_id):
    import uproot
    with uproot.open(path) as f:
        tree=f[tree_name]
        if tree.num_entries!=expected:raise ValueError('Selected entry count mismatch')
        if str(f['selection_cutset_id'])!=cutset_id:raise ValueError('Selected file has wrong cut-set identifier')
        for name,dtype in {**original_types,**EXTRA_TYPES}.items():
            if name not in tree or tree[name].interpretation.numpy_dtype!=dtype:raise ValueError('Missing or changed branch '+name)

def export_dataset(job):
    import uproot
    dataset,source,tree_name,cutfile,output,step_size,max_entries=job
    selection=ChannelSelection(cutfile)
    root=Path(output)/'selected';root.mkdir(parents=True,exist_ok=True)
    reports=Path(output)/'datasets';reports.mkdir(exist_ok=True)
    source=Path(source);stat=source.stat()
    identity=dict(path=str(source.resolve()),size=stat.st_size,mtime_ns=stat.st_mtime_ns)
    final=root/('selected_'+dataset+'.root');report=reports/(dataset+'.json')
    if final.exists() and report.exists():
        old=json.loads(report.read_text())
        if old['source']==identity and old['cutset_id']==selection.cutset_id and old['max_entries']==max_entries:
            audit_file(final,tree_name,old['selected_rows'],{k:np.dtype(v) for k,v in old['original_types'].items()},selection.cutset_id)
            print('SELECTED CACHE '+dataset,flush=True);return old
    temporary=final.with_suffix('.partial.root')
    selected_rows=read=0;counts=np.zeros((24,2),dtype=np.int64)
    runs={}
    started=time.monotonic()
    save_json(reports/(dataset+'.progress.json'),dict(state='opening',updated=time.time()))
    with uproot.open(source) as f,uproot.recreate(temporary,compression=uproot.ZLIB(1)) as target:
        tree=f[tree_name]
        original_types={name:branch.interpretation.numpy_dtype for name,branch in tree.items()}
        if set(DOWNSTREAM)-set(original_types):raise ValueError('Missing downstream metadata '+str(sorted(set(DOWNSTREAM)-set(original_types))))
        if set(original_types)&set(EXTRA_TYPES):raise ValueError('Input already has selection annotations')
        if any(dtype.kind not in 'biuf' or dtype.shape for dtype in original_types.values()):raise ValueError('Expected scalar numeric input branches')
        target.mktree(tree_name,{**original_types,**EXTRA_TYPES},title='Nominal preliminary exclusive e prime p pi minus candidates')
        stop=min(tree.num_entries,max_entries or tree.num_entries)
        for arrays in tree.iterate(library='np',step_size=step_size,entry_stop=stop):
            read+=len(arrays['x'])
            bins=bin_indices(arrays['x'],-arrays['tprime'])
            mask,extra=selection.evaluate(arrays,dataset)
            selected_rows+=int(mask.sum())
            in_bin=bins>=0
            counts[:,0]+=np.bincount(bins[in_bin],minlength=24)
            counts[:,1]+=np.bincount(bins[mask],minlength=24)
            if np.any(mask):
                if any(not np.all(np.isfinite(arrays[k][mask])) for k in DOWNSTREAM):raise ValueError('Nonfinite downstream metadata')
                if not np.all(np.isin(arrays['helicity'][mask],[-1,1])):raise ValueError('Invalid selected beam helicity')
                if (dataset.startswith('RGB_') or not dataset.endswith('_ND3')) and np.any(arrays['target_pol'][mask]!=0):raise ValueError('Unpolarized target has nonzero target polarization')
            if np.any(mask):target[tree_name].extend({**{k:v[mask] for k,v in arrays.items()},**extra})
            if 'runnum' in arrays:
                for r in np.unique(arrays['runnum']):
                    key=str(int(r));n=runs.setdefault(key,dict(source_combinations=0,selected_combinations=0))
                    n['source_combinations']+=int(np.count_nonzero(arrays['runnum']==r))
                    n['selected_combinations']+=int(np.count_nonzero(mask&(arrays['runnum']==r)))
            save_json(reports/(dataset+'.progress.json'),dict(state='selecting',percent=100*read/max(1,stop),selected_rows=selected_rows,updated=time.time()))
            print(f'{dataset}: selected {selected_rows:,}, read {100*read/max(1,stop):.1f}%',flush=True)
        active=selection.active_bins(dataset)
        value=dict(dataset=dataset,fit_group=selection.group(dataset),source=identity,tree=tree_name,cutset_id=selection.cutset_id,version=selection.payload['version'],max_entries=max_entries,source_entries=tree.num_entries,read_rows=read,selected_rows=selected_rows,partial=read<tree.num_entries,active_bins=active.tolist(),reference_assisted_bins=(np.flatnonzero(selection.assisted_bins(dataset))+1).tolist(),unavailable_bins=(np.flatnonzero(~active)+1).tolist(),counts=counts.tolist(),runs=runs,original_types={k:v.str for k,v in original_types.items()},elapsed_seconds=time.monotonic()-started,selected_file=str(final.resolve()),normalization='Combination counts are not luminosity or charge. Retain original QA-selected run/charge normalization independently of event exclusivity cuts.')
        target['selection_cutset_id']=selection.cutset_id
        target['selection_manifest']=json.dumps(value,sort_keys=True)
        if 'kinematic_convention' in f:target['source_kinematic_convention']=f['kinematic_convention'].member('fTitle')
    audit_file(temporary,tree_name,selected_rows,original_types,selection.cutset_id)
    temporary.replace(final)
    save_json(report,value)
    with (reports/(dataset+'.csv')).open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['bin','status','source_combinations_in_bin','selected_combinations','cutset_id'])
        for b in range(24):w.writerow([b+1,'available' if active[b] else 'unavailable_no_nominal_cut',counts[b,0],counts[b,1] if active[b] else '',selection.cutset_id])
    save_json(reports/(dataset+'.progress.json'),dict(state='complete',selected_rows=selected_rows,updated=time.time()))
    return value
