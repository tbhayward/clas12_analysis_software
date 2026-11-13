#!/usr/bin/env python3
"""Isolated full-statistics production-counter replay, then fixed-response pi0 test."""
import argparse,csv,json,os,shlex,shutil,subprocess
from pathlib import Path
import numpy as np
import study_particle_response as response
import study_response_cut_redetermination as candidates
import study_production_pi0_audit as audit
SOURCE=Path(__file__).resolve().parent


def replace_once(text,old,new):
    if text.count(old)!=1:raise ValueError('Production source layout changed; refusing ambiguous instrumentation: '+old[:60])
    return text.replace(old,new,1)


def instrument(text):
    text=replace_once(text,'#include "total_counts.h"','#include "total_counts.h"\n#include "particle_control_core.h"')
    text=replace_once(text,'struct BranchBinder {','static nlohmann::json pi0_test_parameters;\nstatic uint64_t pi0_test_seed=1031;\nstruct BranchBinder {\nbool response_unresolved=false;\nstd::array<double,7> response_extra{};\nstd::array<bool,7> response_extra_ok{};')
    # Force measured coordinates only in the isolated response scenario.
    text=replace_once(text,'        t->StopCacheLearningPhase();','''        if(!pi0_test_parameters.empty()&&work_cfg.sample_kind==SampleKind::MC_REC){
            for(auto n:{"e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi"})ena(n);
            if(work_cfg.channel_cfg.channel==Channel::EPPI0)for(auto n:{"gamma_phi1","gamma_phi2","open_angle_egamma1","open_angle_egamma2","Mh_gammagamma","detector_gamma1","detector_gamma2"})ena(n);
        }
        t->StopCacheLearningPhase();''')
    text=replace_once(text,'        if (analysis_uses_photon_t() &&', '''        if(!pi0_test_parameters.empty()&&work_cfg.sample_kind==SampleKind::MC_REC&&work_cfg.channel_cfg.channel==Channel::EPPI0){
            const char* names[7]={"gamma_phi1","gamma_phi2","open_angle_egamma1","open_angle_egamma2","Mh_gammagamma","detector_gamma1","detector_gamma2"};
            for(int i=0;i<7;++i)bD(names[i],&response_extra[i],response_extra_ok[i]);
        }
        if (analysis_uses_photon_t() &&''')
    text=replace_once(text,'struct BranchBinder {','struct BranchBinder {') # Also confirms uniqueness after instrumentation.
    start=text.index('struct BranchBinder {');end=text.index('\n};',start)+3
    text=text[:end]+'\n#include "production_pi0_response_adapter.h"\n'+text[end:]
    anchor='''        TopologyIndex topo = TopologyIndex::INVALID;
        const std::string* topo_dir_ptr = nullptr;'''
    start=text.index('static WorkCounts accumulate_counts_for_tree_legacy(');end=text.index('\nstatic WorkCounts ',start+1)
    part=text[start:end]
    if anchor not in part:raise ValueError('Legacy counter event-loop anchor missing')
    part=part.replace(anchor,'        apply_pi0_test_response(b,work_cfg,tags,i);\n'+anchor,1)
    row_anchor='            matched_any = true;'
    if row_anchor not in part:raise ValueError('Production row-match anchor missing')
    part=part.replace(row_anchor,'            audit_pi0_recovery_row(b,work_cfg,tags,r);\n'+row_anchor,1)
    part=part.replace('    return out;','    flush_pi0_recovery_audit();\n    return out;')
    text=text[:start]+part+text[end:]
    return text


def build(folder):
    folder.mkdir();text=(SOURCE/'total_counts.cpp').read_text();(folder/'counter.cpp').write_text(instrument(text))
    for path in SOURCE.glob('*.h'):shutil.copy2(path,folder/path.name)
    worker=(SOURCE/'exclusivity_resolution_worker.cpp').read_text();cfg=worker[worker.index('static GlobalCutConfig config('):worker.index('struct Reader {')]
    main=r'''
#include "load_trees.h"
using json=nlohmann::json;
static nlohmann::json readj(const std::string& p){std::ifstream f(p);if(!f)throw std::runtime_error("Cannot read "+p);nlohmann::json j;f>>j;return j;}
'''+cfg+r'''
int main(int argc,char**argv){try{
 if(argc!=7)throw std::runtime_error("counter CSV CUTS GLOBALS CURRENT RESPONSE_JSON_OR_NONE WORKERS");
 set_default_global_cuts(config(readj(argv[3])));
 if(std::string(argv[5])!="none"){auto j=readj(argv[5]);pi0_test_parameters=j.at("periods");pi0_test_seed=j.at("seed");}
 std::map<std::string,TTree*> d,g,r,pd,pg,pr,b,rg,rr,empty;
 if(!loadTrees(d,g,r,pd,pg,pr,b,rg,rr))throw std::runtime_error("Production tree loader failed");
 TotalCountsOptions o;o.make_plots=false;o.make_note_outputs=false;o.write_current_nuisance_responses=false;o.current_response_model_json=argv[4];o.require_sp18_out_epg_e_theta_current_model=false;
 if(!update_total_counts_csv(argv[1],empty,empty,empty,empty,empty,pr,b,argv[2],"counter_diagnostics",std::stoi(argv[6]),o))return 2;
 nlohmann::json audit;for(const auto& a:pi0_recovery_counts)audit["trees"][a.first]={{"attempted",a.second[0]},{"unresolved_retained",a.second[1]}};
 for(const auto& a:pi0_recovery_bins)audit["production_rows"][a.first]={{"selected",a.second[0]},{"unresolved_retained",a.second[1]}};
 std::ofstream("photon_recovery_audit.json")<<audit.dump(2);
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
'''
    (folder/'counter.cpp').write_text((folder/'counter.cpp').read_text()+main)
    flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True));libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True));exe=folder/'production_pi0_counter'
    cmd=[os.environ.get('CXX','c++'),*flags,'-O2','-std=c++17','-pthread','-I'+str(folder),str(folder/'counter.cpp'),str(SOURCE/'global_cuts.cpp'),str(SOURCE/'load_trees.cpp'),'-o',str(exe),*libs]
    # ROOT's worker budget is capped by the production counter. Enable its existing OpenMP path on Linux.
    if os.uname().sysname=='Linux':cmd.insert(1,'-fopenmp')
    subprocess.run(cmd,check=True);return exe


def yield_columns(header):
    return [c for c in header if c.startswith(('reconstructed yield, ep->eppi0,','reconstructed yield, ep->eppi0->epg,','reconstructed current corrected yield, ep->eppi0,','reconstructed current corrected yield, ep->eppi0->epg,'))]


def gate(baseline,replay,out):
    a=list(csv.DictReader(baseline.open()));b=list(csv.DictReader(replay.open()))
    if len(a)!=len(b):raise ValueError('Row count changed')
    mismatches=[];checks=0
    for i,(x,y) in enumerate(zip(a,b)):
        if any(x[k]!=y[k] for k in audit.EDGES):raise ValueError('Production bin edges changed')
        if str(x.get('valid bin','1')).strip()=='0':continue
        for c in yield_columns(list(x)):
            try:v=audit.triple(x[c]);w=audit.triple(y[c])
            except ValueError as exc:raise ValueError(f'Yield parse failure at row {i}, column {c!r}; saved={x[c]!r}; replay={y[c]!r}') from exc
            checks+=1
            if not np.allclose(v,w,rtol=2e-8,atol=2e-7):mismatches.append(dict(row=i,column=c,saved=str(v.tolist()),replay=str(w.tolist())))
    audit.table(out/'nominal_mismatches.csv',mismatches)
    response.dump(out/'nominal_gate.json',dict(scope='Saved valid production rows; terminal invalid rows excluded',checked_cells=checks,mismatches=len(mismatches),passed=not mismatches))
    if mismatches:raise ValueError(f'Nominal replay differs in {len(mismatches)} cells: inspect nominal_mismatches.csv. Morph not run.')


def prepare_csv(baseline,target):
    with baseline.open() as f:
        reader=csv.DictReader(f);header=reader.fieldnames;rows=list(reader)
    columns=yield_columns(header)
    for row in rows:
        for col in columns:row[col]='(0.00000000,0.00000000,0.00000000)'
    with target.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=header);writer.writeheader();writer.writerows(rows)


def topology_comparison(nominal,morphed,out):
    with nominal.open() as f:a=list(csv.DictReader(f))
    with morphed.open() as f:b=list(csv.DictReader(f))
    rows=[]
    for old,new in zip(a,b):
        for period in audit.PERIODS:
            for top in audit.TOPS:
                row=dict(period=period,topology=top,**{k:old[k] for k in audit.EDGES})
                for channel,label in [('ep->eppi0','control'),('ep->eppi0->epg','leakage')]:
                    col=f'reconstructed current corrected yield, {channel}, {top}, mc, {period}'
                    v=audit.triple(old[col])[0];w=audit.triple(new[col])[0];row['nominal_'+label]=v;row['morphed_'+label]=w;row[label+'_ratio']=w/v if v>0 else float('nan')
                rows.append(row)
    audit.table(out/'topology_mc_comparison.csv',rows)


def make_plots(rows,out):
    import matplotlib.pyplot as plt
    summary=[];(out/'plots').mkdir(exist_ok=True)
    for period in audit.PERIODS:
        a=[r for r in rows if r['period']==period];bins=sorted({(float(r['t_abs_min']),float(r['t_abs_max'])) for r in a});fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
        for tag,color in [('nominal','#777777'),('morphed','#22946b')]:
            yy=[]
            for lo,hi in bins:
                rr=[r for r in a if float(r['t_abs_min'])==lo and float(r['t_abs_max'])==hi];den=sum(r['data_dvcs'] for r in rr);fraction=sum(r[tag+'_contamination']*r['data_dvcs'] for r in rr)/den if den>0 else float('nan');yy.append(100*fraction)
                summary.append(dict(period=period,scenario=tag,t_low=lo,t_high=hi,data_dvcs=den,contamination_fraction=fraction))
            ax.plot([(lo+hi)/2 for lo,hi in bins],yy,'o-',color=color,label=tag)
        ax.set(xlabel='|t| (GeV²)',ylabel='Contamination (%)',title=period+' · production-bin calculation, DATA-yield-weighted aggregation');ax.legend();ax.grid(alpha=.2);fig.savefig(out/'plots'/(period.replace(' ','_')+'.png'),dpi=120);plt.close(fig)
    audit.table(out/'summary.csv',summary)


def completed_count(out,scenario):
    target=out/(scenario+'.csv');marker=out/(scenario+'_count_complete.json')
    if not target.is_file():return False
    if marker.is_file():
        if json.loads(marker.read_text())['csv_sha256']!=response.digest(target):raise ValueError('Completed counts changed: '+scenario)
        return True
    # Recovery of the original parser crash: the C++ process returned successfully
    # before the Python gate. Require its post-write success marker for this CSV.
    log=out/(scenario+'.log')
    return log.is_file() and ('[total_counts] Updated data and MC count columns in: '+str(target)) in log.read_text()


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--resume',action='store_true');p.add_argument('--rebuild-backend',action='store_true');p.add_argument('--workers',type=int,default=10);p.add_argument('--production-run',type=Path,required=True);p.add_argument('--response-followup',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args(argv)
    run=args.production_run.resolve();out=args.output.resolve();follow=args.response_followup.resolve()/'cut_redetermination'
    if out.exists() and not args.resume:p.error('Output exists; use --resume for checked recovery')
    if args.resume and not (out/'manifest.json').is_file():p.error('Resume requires the original manifest')
    if not 1<=args.workers<=10:p.error('Workers must be 1–10 total')
    baseline=run/'output/csvs/dvcs_pass2_analysis.csv';cuts=run/'output/jsons/combined_cuts.json';globals=run/'output/jsons/global_cuts_config.json';current=run/'output/dvcs_current_dependence/calibration/current_response_model.json'
    for path in (baseline,cuts,globals,current):
        if not path.is_file():p.error('Missing saved production input: '+str(path))
    meta=json.loads((follow/'manifest.json').read_text());training=Path(meta['arguments']['training_study']);spec=json.loads((follow/'spec.json').read_text());pars={}
    for period in spec['periods']:
        path=training/period/'add_pt/candidate_parameters.json'
        if response.digest(path)!=meta['candidate_hashes']['add_pt'][period]:p.error('Response parameters changed')
        pars[period]=candidates.load_candidate_parameters(path)
    current_manifest=dict(arguments={k:str(v) for k,v in vars(args).items()},production_inputs={str(x):response.digest(x) for x in (baseline,cuts,globals,current)},counter_sha256=response.digest(SOURCE/'total_counts.cpp'),response_parameters=pars)
    if args.resume:
        previous=json.loads((out/'manifest.json').read_text())
        for key in ('production_inputs','counter_sha256','response_parameters'):
            if previous[key]!=current_manifest[key]:p.error('Cannot resume: changed '+key)
        saved_response=json.loads((out/'response.json').read_text())
        if saved_response!=dict(seed=spec.get('seed',1031),periods=pars):p.error('Saved response configuration changed')
    else:
        out.mkdir(parents=True);response.dump(out/'response.json',dict(seed=spec.get('seed',1031),periods=pars));response.dump(out/'manifest.json',current_manifest)
    exe=out/'backend/production_pi0_counter'
    if args.rebuild_backend:
        if completed_count(out,'morphed'):p.error('Morphed counts already complete; use a new output for a changed backend')
        if (out/'backend').exists():
            backup=out/'backend_before_rebuild'
            if backup.exists():p.error('Backend backup already exists; preserve it before another rebuild')
            (out/'backend').rename(backup)
    if not exe.is_file():
        if (out/'backend').exists():raise ValueError('Incomplete backend; preserve it and choose a fresh output')
        exe=build(out/'backend')
    for scenario in ('nominal','morphed'):
        target=out/(scenario+'.csv')
        if args.resume and completed_count(out,scenario):
            print('[reuse] '+scenario+' completed production counts; no ROOT reread',flush=True)
        else:
            prepare_csv(baseline,target)
            print('[count] '+scenario+' full-statistics pi0 control/leakage only; production binning, cuts and current weights; ten-worker cap',flush=True)
            try:
                with (out/(scenario+'.log')).open('w') as log:subprocess.run([str(exe),str(target),str(cuts),str(globals),str(current),'none' if scenario=='nominal' else str(out/'response.json'),str(args.workers)],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
            except subprocess.CalledProcessError as exc:
                tail='\n'.join((out/(scenario+'.log')).read_text().splitlines()[-30:])
                raise RuntimeError(scenario+' counter failed; completed nominal counts remain reusable. Log tail:\n'+tail) from exc
        response.dump(out/(scenario+'_count_complete.json'),dict(csv_sha256=response.digest(target)))
        if scenario=='nominal':gate(baseline,target,out)
    recovery_rows=[]
    recovery_path=out/'photon_recovery_audit.json'
    if recovery_path.is_file():
        recovery=json.loads(recovery_path.read_text()) or {}
        for key,value in recovery.get('production_rows',{}).items():
            recovery_rows.append(dict(tree_row=key,selected=value['selected'],unresolved_retained=value['unresolved_retained'],unresolved_fraction=value['unresolved_retained']/value['selected'] if value['selected'] else 0.))
        audit.table(out/'photon_recovery_rows.csv',recovery_rows)
        print('[recovery] Selected unresolved control events retained unchanged: '+str(sum(r['unresolved_retained'] for r in recovery_rows)),flush=True)
    topology_comparison(out/'nominal.csv',out/'morphed.csv',out)
    rows0=audit.replay(out/'nominal.csv');rows1=audit.replay(out/'morphed.csv');result=[]
    for a,b in zip(rows0,rows1):
        result.append(dict(period=a['period'],**{k:a[k] for k in audit.EDGES},nominal_contamination=a['replayed_fraction'],morphed_contamination=b['replayed_fraction'],change_percentage_points=100*(b['replayed_fraction']-a['replayed_fraction']),nominal_mc_leakage=a['mc_leakage'],morphed_mc_leakage=b['mc_leakage'],nominal_mc_pi0=a['mc_pi0'],morphed_mc_pi0=b['mc_pi0'],data_dvcs=a['data_dvcs'],signal_yield_ratio=(1-b['replayed_fraction'])/(1-a['replayed_fraction']) if a['replayed_fraction']<1 else float('nan')))
    audit.table(out/'contamination_comparison.csv',result)
    make_plots(result,out)
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><h1>Production pi0 response comparison</h1><p>Full MC statistics, all production rows/periods/topologies. Original DATA yields, original separate DATA/MC production windows, exact production counter and current weighting retained. No control reweighting or fit.</p><p>The nominal gate must pass before morphing. The morph changes reconstructed coordinates before global/exclusivity cuts, phi matching and event current weights. Individual control photons are varied and the pi0 is rebuilt. Unresolved recovery retains the original event unchanged, as in previous response studies; photon_recovery_audit.json reports attempts and unresolved events selected in every production row. This is an incomplete morph for those events, not validated by the nominal gate. Missing required branches remain fatal. Shared trees are read only for the two relevant MC channels, once per scenario. DATA/other channels are not counted.</p><p><a href="photon_recovery_rows.csv">Unresolved recovery by selected production row</a> · <a href="photon_recovery_audit.json">Recovery audit</a> · <a href="topology_mc_comparison.csv">Per-topology MC changes</a> · <a href="summary.csv">Yield-weighted summary</a> · <a href="nominal_gate.json">Nominal gate</a> · <a href="contamination_comparison.csv">Per-bin contamination and MC factors</a> · <a href="nominal.csv">Nominal CSV</a> · <a href="morphed.csv">Morphed-count CSV</a></p><p>Stored contamination columns in morphed.csv are inherited baseline values: use contamination_comparison.csv for recalculated fractions. Signal-yield ratios isolate subtraction only, not acceptance or final cross sections. Plots are under plots/PERIOD.png. The nominal gate covers saved valid production bins; terminal invalid bins are not used as physics validation. See logs for actual production loader inputs. A gate failure can indicate changed processing/calibration/code provenance and must not be bypassed.</p>')
    print('[done] '+str(out/'REPORT.html'),flush=True)

if __name__=='__main__':main()
