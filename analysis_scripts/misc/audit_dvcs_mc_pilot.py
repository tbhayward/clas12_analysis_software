#!/usr/bin/env python3
"""Read-only pilot audit for DVCS/π0 reconstructed MC ROOT triplets. Requires PyROOT."""
import argparse, collections, csv, math
from pathlib import Path
import ROOT
ROOT.gROOT.SetBatch(True)

def tree_open(path, name):
    f=ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie(): raise RuntimeError(f'Cannot open {path}')
    t=f.Get(name)
    if not t: f.Close(); raise RuntimeError(f'Missing {name}: {path}')
    return f,t

def branches(t): return {b.GetName() for b in t.GetListOfBranches()}
def val(t,k): return getattr(t,k)
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',default='/volatile/clas12/thayward/dvcs_reprocessing_20261022/pilot_10k/fa18_inb')
    ap.add_argument('--output',default='dvcs_mc_pilot_audit.csv')
    ap.add_argument('--max-files',type=int,default=0,help='0 audits all files')
    args=ap.parse_args()
    base=Path(args.root); files=sorted(p for p in base.rglob('*.mc_events.root') if p.parent.parent.name in ('dvcs_rec','pi0_rec','pi0_as_dvcs'))
    if args.max_files: files=files[:args.max_files]
    if not files: raise SystemExit(f'No MC event sidecars under {base}')
    results=[]
    for ep in files:
        prefix=str(ep)[:-len('.mc_events.root')]
        pp=Path(prefix+'.mc_particles.root'); rp=Path(prefix+'.root')
        r={'mode':ep.parent.parent.name,'current':ep.parent.name,'file':Path(prefix).name,'status':'PASS','issues':''}
        issues=[]; opened=[]
        try:
            for path,name in ((ep,'MCEvents'),(pp,'MCEvents'),(rp,'PhysicsEvents')):
                f,t=tree_open(path,name); opened.append(f)
                if path==ep: et=t
                elif path==pp: pt=t
                else: rt=t
            eb=branches(et); pb=branches(pt); rb=branches(rt)
            for k in ('file_index','event_sequence','run','event','n_mc','n_rec','n_match','gen_Q2','gen_xB','gen_y','gen_t_proton'):
                if k not in eb: issues.append('missing event '+k)
            for k in ('file_index','event_sequence','run','event','mcindex','truth_pid','rec_pindex','rec_pid','fitted_selected','gen_px','gen_py','gen_pz'):
                if k not in pb: issues.append('missing particle '+k)
            for k in ('runnum','evnum','Q2','x','t1'):
                if k not in rb: issues.append('missing reco '+k)
            r.update(events=int(et.GetEntries()),particles=int(pt.GetEntries()),reco_candidates=int(rt.GetEntries()),event_branches=len(eb),particle_branches=len(pb),reco_branches=len(rb))
            if issues: raise RuntimeError('; '.join(issues))
            ids=set(); expected_particles=0; bad_kin=0; no_truth=0; event_evnums=collections.Counter()
            for e in et:
                key=(int(e.file_index),int(e.event_sequence))
                if key in ids: issues.append('duplicate MCEvents sequence'); break
                ids.add(key); event_evnums[(int(e.run),int(e.event))]+=1
                expected_particles+=int(e.n_mc)
                if int(e.n_mc)<3: no_truth+=1
                if not all(math.isfinite(float(getattr(e,k))) for k in ('gen_Q2','gen_xB','gen_y','gen_t_proton')): bad_kin+=1
            particle_ids=set(); n_part=0; bad_indices=0; matched=0; misid=0; selected=0; proton=0; proton_matched=0
            for p in pt:
                key=(int(p.file_index),int(p.event_sequence))
                if key not in ids: bad_indices+=1
                particle_ids.add(key); n_part+=1
                if int(p.rec_pindex)>=0:
                    matched+=1
                    if int(p.rec_pid)!=int(p.truth_pid): misid+=1
                if int(p.fitted_selected)==1: selected+=1
                if int(p.truth_pid)==2212:
                    proton+=1
                    if int(p.rec_pindex)>=0: proton_matched+=1
            r.update(matched_particles=matched,misidentified_matched=misid,fitted_selected_particles=selected,truth_protons=proton,matched_protons=proton_matched,nonfinite_gen_kinematics=bad_kin,missing_mc_truth_events=no_truth)
            if bad_indices: issues.append(f'{bad_indices} particle rows with missing event sequence')
            if n_part>expected_particles: issues.append('more diagnostic particle rows than MC::Particle rows')
            if not particle_ids.issubset(ids): issues.append('particle event keys not subset of MCEvents')
            unmatched_reco=0; ambiguous_reco=0
            for e in rt:
                key=(int(e.runnum),int(e.evnum))
                if key not in event_evnums: unmatched_reco+=1
                elif event_evnums[key]>1: ambiguous_reco+=1
            r.update(reco_not_in_mc_events=unmatched_reco,reco_ambiguous_event_ids=ambiguous_reco)
            if unmatched_reco: issues.append(f'{unmatched_reco} reco entries have no MCEvents run/event key')
            if ambiguous_reco: issues.append(f'{ambiguous_reco} reco entries have ambiguous run/event key')
            if bad_kin: issues.append(f'{bad_kin} nonfinite generated kinematics')
            if no_truth: issues.append(f'{no_truth} MCEvents with fewer than three MC particles')
            if not ('gen_phi' in eb): r['gen_phi_present']='NO'
            else: r['gen_phi_present']='YES'
        except Exception as exc: issues.append('ERROR '+str(exc))
        finally:
            for f in opened: f.Close()
        r['status']='REVIEW' if issues else 'PASS';r['issues']=' | '.join(dict.fromkeys(issues))
        results.append(r)
        print(f"{r['status']:6} {r['mode']:12} {r['current']:7} {r['file']} events={r.get('events','?')} reco={r.get('reco_candidates','?')} {r['issues']}",flush=True)
    cols=['status','mode','current','file','events','particles','reco_candidates','event_branches','particle_branches','reco_branches','matched_particles','misidentified_matched','fitted_selected_particles','truth_protons','matched_protons','nonfinite_gen_kinematics','missing_mc_truth_events','reco_not_in_mc_events','reco_ambiguous_event_ids','gen_phi_present','issues']
    with open(args.output,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(results)
    print(f'Wrote {args.output}; {len(results)} files; {sum(r["status"]!="PASS" for r in results)} require review.')
if __name__=='__main__':main()
