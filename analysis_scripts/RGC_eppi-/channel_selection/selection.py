"""Shared, versioned cuts for all e' p pi- target samples and downstream analyses."""
import hashlib
import inspect
import json
from pathlib import Path
import numpy as np
from candidate_cache import select,derive
from kinematics import XB,TP,BRANCHES

VERSION='eppi_minus_preliminary_v1'
GROUPS={
 'RGB_spring2019':['RGB_spring2019_D2'],
 'RGB_fall2019':['RGB_fall2019_D2'],
 'RGB_spring2020':['RGB_spring2020_D2'],
 'RGC_summer2022':['RGC_summer2022_ND3'],
 'RGC_fall2022_spring2023':['RGC_fall2022_ND3','RGC_spring2023_ND3'],
}

def group_for_dataset(name):
    if name in [n for v in GROUPS.values() for n in v]:
        return next(g for g,names in GROUPS.items() if name in names)
    for period,group in [('summer2022','RGC_summer2022'),('fall2022','RGC_fall2022_spring2023'),('spring2023','RGC_fall2022_spring2023')]:
        if name in [f'RGC_{period}_{t}' for t in ['C','CH2','CD2','ET','He']]:return group
    raise ValueError('Unknown dataset: '+name)

def fingerprint():
    import candidate_cache,kinematics
    return hashlib.sha256(Path(candidate_cache.__file__).read_bytes()+Path(kinematics.__file__).read_bytes()+Path(__file__).read_bytes()).hexdigest()

def seal(payload):
    return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

class ChannelSelection:
    def __init__(self,path,allow_partial=False):
        self.path=Path(path);self.payload=json.loads(self.path.read_text())
        cutset_id=self.payload.pop('cutset_id')
        if seal(self.payload)!=cutset_id:raise ValueError('Cut-set checksum mismatch')
        self.cutset_id=cutset_id
        if self.payload['version']!=VERSION:raise ValueError('Unsupported selection version')
        if self.payload['base_selection_sha256']!=fingerprint():raise ValueError('Selection/kinematics code does not match cut set')
        if self.payload['partial'] and not allow_partial:raise ValueError('Partial calibration cannot be used as nominal selection')
        self.groups=self.payload['groups']
        self.group_ids=self.payload.get('group_ids',{name:i+1 for i,name in enumerate(GROUPS)})
    def group(self,dataset):
        g=group_for_dataset(dataset)
        if g not in self.groups:raise ValueError('Missing fitted group '+g)
        return g
    def active_bins(self,dataset):
        return np.array([f['active'] for f in self.groups[self.group(dataset)]['bins']],dtype=bool)
    def assisted_bins(self,dataset):
        return np.array([f.get('reference_assisted',False) for f in self.groups[self.group(dataset)]['bins']],dtype=bool)
    def candidate_mask(self,mass,momentum,bins,dataset):
        entries=self.groups[self.group(dataset)]['bins']
        active=np.array([v['active'] for v in entries])
        ml=np.array([v['mass_low'] if v['active'] else np.inf for v in entries]);mh=np.array([v['mass_high'] if v['active'] else -np.inf for v in entries])
        pl=np.array([v['momentum_low'] if v['active'] else np.inf for v in entries]);ph=np.array([v['momentum_high'] if v['active'] else -np.inf for v in entries])
        valid=(bins>=0)&(bins<24)&np.isfinite(mass)&np.isfinite(momentum)
        safe=np.clip(bins,0,23)
        return valid&active[safe]&(mass>=ml[safe])&(mass<=mh[safe])&(momentum>=pl[safe])&(momentum<=ph[safe])
    def evaluate(self,arrays,dataset):
        """Return full-length boolean mask and annotations. No weights/polarization cuts."""
        good,bins=select(arrays)
        ids=np.flatnonzero(good)
        a={k:arrays[k][good] for k in BRANCHES}
        if 'Emiss' in arrays:
            from kinematics import MD,MN
            emiss=arrays['Emiss'][good]
            if not np.all(np.isfinite(emiss)):raise ValueError('Nonfinite saved missing energy')
            p2=emiss*emiss-a['Mx2']
            if np.any(p2 < -1e-8*np.maximum(1,emiss*emiss+np.abs(a['Mx2']))):raise ValueError('Inconsistent saved missing energy and mass')
            momentum=np.sqrt(np.maximum(p2,0))
            mass=a['Mx2']+2*(MD-MN)*emiss+(MD-MN)**2
            # Verify the saved-energy shortcut independently for each input chunk.
            n=min(32,len(mass))
            mcheck,pcheck,_=derive({k:v[:n] for k,v in a.items()})
            if not np.allclose(mass[:n],mcheck,rtol=0,atol=1e-6) or not np.allclose(momentum[:n],pcheck,rtol=0,atol=1e-6):raise ValueError('Saved missing energy does not match four-vector reconstruction')
        else:
            mass,momentum,_=derive(a)
        accepted=self.candidate_mask(mass,momentum,bins[good],dataset)
        mask=np.zeros(len(good),dtype=bool);mask[ids[accepted]]=True
        # Angular residual is needed only for retained rows during full-branch exports.
        _,_,angle=derive({k:v[accepted] for k,v in a.items()})
        # Annotation arrays correspond exactly to mask-selected original rows.
        extra=dict(analysis_bin=(bins[good][accepted]+1).astype(np.int32),exclusivity_group=np.full(accepted.sum(),self.group_ids[self.group(dataset)],dtype=np.int32),Mx2_deuteron=mass[accepted],p_miss=momentum[accepted],delta_theta_deg=angle,exclusivity_reference_assisted=self.assisted_bins(dataset)[bins[good][accepted]].astype(np.int32))
        return mask,extra


def load_selected_inputs(path,verify_files=True):
    """Resolve a portable output registry and reject mixed cut-set provenance."""
    path=Path(path)
    registry=json.loads(path.read_text())
    cut_path=Path(registry['cut_file'])
    if not cut_path.is_absolute():cut_path=path.parent/cut_path
    selection=ChannelSelection(cut_path)
    if registry['cutset_id']!=selection.cutset_id:raise ValueError('Registry has a different cut set')
    registry['cut_file']=str(cut_path.resolve())
    for name,item in registry['files'].items():
        filename=Path(item['path'])
        if not filename.is_absolute():filename=path.parent/filename
        item['path']=str(filename.resolve())
        if verify_files:
            import uproot
            with uproot.open(filename) as f:
                if str(f['selection_cutset_id'])!=selection.cutset_id:raise ValueError('Mixed selected-file cut sets')
                if f[registry['tree']].num_entries!=item['selected_rows']:raise ValueError('Selected entry count mismatch')
                stored=json.loads(str(f['selection_manifest']))
                if stored['dataset']!=name or stored['fit_group']!=selection.group(name):raise ValueError('Selected dataset identity mismatch')
    return registry
