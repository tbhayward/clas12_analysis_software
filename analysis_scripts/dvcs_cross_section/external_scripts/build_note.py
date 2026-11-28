"""Build the v2 analysis note: python3 build_note.py review, full or both."""
import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
import sys
import importlib.util
from difflib import SequenceMatcher

sys.dont_write_bytecode = True
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np
import pandas as pd


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def number(value, index=0):
    if not value.strip():
        return np.nan
    return float(value.strip("()").split(",")[index])


def json_safe(value):
    """Use JSON null for diagnostics undefined by a fit scenario."""
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    return value


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected_rows(rows, profile):
    # Membership follows nominal edges, never shifted reporting coordinates.
    return [r for r in rows if profile == "full" or float(r["t_abs_min"]) >= .25]


def summary(rows, label):
    base = f"normed cross sections, ep->epg, exp, {label}, unpol"
    rows = [r for r in rows if np.isfinite(number(r[base])) and number(r[base]) > 0]
    sigma = np.array([number(r[base]) for r in rows])

    def column(key, index=0):
        return np.array([number(r[key], index) for r in rows])

    values = [column(base, 1) / sigma]
    for key in ["Syst. err (pi0 subtraction)", "Syst. err (Acceptance)",
                "Syst.err (Frad)", "Syst.err (Fbin)", "Syst. err (exclusivity cuts)",
                "Syst. err (fiducial cuts)", "Syst. err (point-to-point total)"]:
        denominator = sigma
        if label == "Sp19 Inb":
            if "cuts)" in key:
                denominator = column(base.replace(label, "10.6 GeV"))
            else:
                key += ", Sp19 Inb (10.2 GeV)"
        with np.errstate(invalid="ignore", divide="ignore"):
            values.append(column(key) / denominator)
    for key in ["current dependence", "run period residual", "correlated scale",
                "overall normalization"]:
        values.append(column(f"{key} sys frac, {label}"))
    values.append(np.sqrt(values[7]**2 + values[10]**2 + values[11]**2))
    stats = []
    for v in values:
        finite = 100 * v[np.isfinite(v)]
        if not len(finite):
            raise ValueError(f"No finite support in {label} uncertainty summary")
        stats.append([np.mean(finite), *np.percentile(finite, [50, 16, 84])])
    return rows, np.array(stats), np.array(values)


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"Expected one v2 source passage: {old[:100]}")
    return text.replace(old, new, 1)


def edit_paragraph(text, prefix, transform):
    paragraphs = text.split("\n\n")
    indices = [i for i, p in enumerate(paragraphs) if p.startswith(prefix)]
    if len(indices) != 1:
        raise ValueError(f"Expected one paragraph beginning {prefix}")
    i = indices[0]
    paragraphs[i] = transform(paragraphs[i])
    return "\n\n".join(paragraphs)


def replace_table_rows(text, label, rows):
    for match in re.finditer(r"\\begin\{table\}.*?\\end\{table\}", text, re.S):
        if r"\label{" + label + "}" not in match[0]:
            continue
        table = re.sub(r"(?m)^\d[^\n]* & [^\n]*\n", "", match[0])
        marker = list(re.finditer(r"\\hline", table))[1].end()
        table = table[:marker] + "\n" + "\n".join(" & ".join(row) + r" \\" for row in rows) + table[marker:]
        return text[:match.start()] + table + text[match.end():]
    raise ValueError(f"Missing table {label}")


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def projection_points(sample, values, key):
    # Equal centers do not imply equal intervals, especially for phi.
    edges = [(float(r[key]), float(r[key.replace("min", "max")])) for r in sample]
    points = []
    for lo, hi in sorted(set(edges)):
        mask = np.array([edge == (lo, hi) for edge in edges])
        medians = [100*np.median(v[mask & np.isfinite(v)]) for v in values[1:8]]
        points.append(((lo+hi)/2, lo, hi, medians))
    return sorted(points, key=lambda point: point[0])


def update_summary(root, base, rows):
    selected = selected_rows(rows, "review")
    folder = root/"chapters/systematics/figures/summary"
    samples=[summary(selected,label) for label in ['10.6 GeV','Sp19 Inb']]
    original=(root/'chapters/systematics/summary.tex').read_text()
    s=original
    for label,oldvalues,newvalues in [
     ('10.6 GeV',['9.60','9.05','7.65','7.52'],[samples[0][1][i,j] for i,j in [(0,0),(7,0),(0,1),(7,1)]]),
     ('medians',['2.27','5.20','8.87','9.83'],[samples[e][1][i,1] for e,i in [(0,10),(1,10),(0,12),(1,12)]])]:
     for old,new in zip(oldvalues,newvalues):
      assert s.count('$'+old+r'\%$')==1,old
      s=s.replace('$'+old+r'\%$',f'${new:.2f}\\%$',1)
    tables=iter(samples)
    def update_table(match):
     _,stats,_=next(tables);lines=match[0].splitlines();indices=[i for i,l in enumerate(lines) if re.search(r' & [\d.]+ &',l)]
     assert len(indices)==len(stats)
     for i,(mean,median,lo,hi) in zip(indices,stats):
      cells=lines[i].split(' & ');cells[1:4]=[f'{mean:.2f}',f'{median:.2f}',f'{lo:.2f}--{hi:.2f}'];lines[i]=' & '.join(cells)
     return '\n'.join(lines)
    s=re.sub(r'\\begin\{table\}.*?\\end\{table\}',update_table,s,flags=re.S)
    for command in ['label','caption','includegraphics']:
     pattern=r'\\'+command+r'(?:\[[^\]]*\])?\{([^\n]+)'
     assert re.findall(pattern,s)==re.findall(pattern,original),command
    (root/'chapters/systematics/summary.tex').write_text(s)
    
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    styles=[("#0000cc","o",r"$\pi^0$ subtraction"),("#cc0000","s","Acceptance"),
            ("#009900","^",r"$F_{\rm rad}$"),("#cc00cc","v",r"$F_{\rm bin}$"),
            ("#ff6600","o","Exclusivity"),("#009999","s","Fiducial"),("black","*","Point-to-point total")]
    projections={}
    for axis,key in [("xB","xBmin"),("Q2","Q2min"),("t","t_abs_min"),("phi","phimin")]:
        # First reproduce the full-sample medians stored by the production plotting code.
        production=read_rows(base/("output/systematics/point_to_point_projections/point_to_point_systematics_vs_"+axis+".csv"))
        for energy,label in [("10.6 GeV","10.6 GeV"),("10.2 GeV","Sp19 Inb")]:
            sample,_,values=summary(rows,label)
            expected=projection_points(sample,values,key)
            actual=[point for point in production if point["energy"]==energy]
            assert len(actual)==len(expected),(axis,energy,len(actual),len(expected))
            for x in sorted(set(p[0] for p in expected)):
                # Different phi intervals can share a center. Preserve their edges.
                want=sorted(tuple(p[3]) for p in expected if p[0]==x)
                got=sorted(tuple(float(p[c]) for c in list(p)[2:]) for p in actual if np.isclose(float(p["x"]),x,atol=1e-10,rtol=0))
                assert np.allclose(got,want,rtol=1e-7,atol=1e-8),(axis,energy,x,got,want)
        fig,axes=plt.subplots(1,2,figsize=(15,7.2))
        fig.subplots_adjust(left=.075,right=.97,bottom=.14,top=.73,wspace=.16)
        entries=[]
        for ax,(sample,_,values),energy in zip(axes,samples,["10.6 GeV combined","10.2 GeV Sp19 Inb"]):
            projected=projection_points(sample,values,key)
            x=np.array([p[0] for p in projected]);medians=np.array([p[3] for p in projected]).T
            for j,((color,marker,label),y) in enumerate(zip(styles,medians)):
                ax.plot(x,y,color=color,marker=marker,markersize=4 if j<6 else 6,lw=.8 if j<6 else 2.5,
                        markerfacecolor="none" if j in [4,5] else color,markeredgewidth=.8,label=label)
            dx=.04*(x[-1]-x[0]) if len(x)>1 else .5
            ax.set_xlim(x[0]-dx,x[-1]+dx);ax.set_ylim(0,max(5,1.18*np.max(medians)))
            ax.set_xlabel({"xB":r"$x_B$","Q2":r"$Q^2$ (GeV$^2$)","t":r"$|t|$ (GeV$^2$)","phi":r"$\phi$ (deg)"}[axis],fontsize=16)
            ax.tick_params(which="both",direction="in",top=True,right=True,labelsize=13)
            ax.xaxis.set_minor_locator(AutoMinorLocator());ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.grid(axis="y",linestyle=":",color="black",linewidth=.55)
            ax.set_title(energy,loc="left",fontsize=14,pad=18)
            for i,xx in enumerate(x):entries.append({"energy":energy,"x":float(xx),"low":projected[i][1],"high":projected[i][2],"medians_percent":medians[:,i].tolist()})
        axes[0].set_ylabel("Median relative systematic uncertainty (%)",fontsize=16)
        # Retain the production legend's two-row, four-column arrangement.
        handles,labels=axes[0].get_legend_handles_labels();order=[0,4,1,5,2,6,3]
        fig.legend([handles[i] for i in order],[labels[i] for i in order],ncol=4,loc="upper left",bbox_to_anchor=(.065,.9),frameon=False,fontsize=10,columnspacing=3.6,labelspacing=1.3)
        fig.suptitle("Kinematic dependence of point-to-point systematic uncertainties",fontsize=11,y=.985)
        fig.savefig(folder/("point_to_point_systematics_vs_"+axis+".png"),dpi=180);plt.close(fig)
        projections[axis]=entries
    return {"projections": projections}


def period_consistency(rows):
    periods = ["Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"]
    acc = {p: [0., 0.] for p in periods}
    for row in rows:
        values = {}
        for period in periods:
            cell = row[f"normed cross sections, ep->epg, exp, {period}, unpol"]
            if not cell.strip():
                continue
            value, stat = number(cell), number(cell, 1)
            if np.isfinite(value) and np.isfinite(stat) and stat > 0:
                values[period] = value, stat
        if not values:
            continue
        reference = sum(v/e**2 for v, e in values.values())/sum(1/e**2 for v, e in values.values())
        if not np.isfinite(reference) or reference == 0:
            continue
        for period, (value, stat) in values.items():
            weight = (reference/stat)**2
            acc[period][0] += weight
            acc[period][1] += weight*value/reference
    means = [total/weight for weight, total in acc.values()]
    variances = [1/weight for weight, _ in acc.values()]
    return float(100*np.sqrt(max(0, np.mean((np.array(means)-1)**2)-np.mean(variances))))


def world_comparisons(root, base, nominal, scripts):
    out = base/"output"
    tables = out/"global_comparison_current_001/tables"
    proposal = root/"chapters/systematics/figures/global_comparison"
    world = load_module(scripts/"compare_world_dvcs_cross_sections.py", "note_world")
    ratios = load_module(scripts/"plot_cross_section_ratio_pull_summary.py", "note_ratios")
    retained={key:r for key,r in nominal.items() if float(r['t_abs_min'])>=.25}
    ids={'pass2:'+key for key in retained}
    matches=pd.read_csv(tables/'pairwise_matched_points.csv',low_memory=False)
    anchors=matches[matches.point_id_b.isin(ids)].copy()
    assert set(anchors.dataset_b)=={'pass2'}
    for r in anchors.itertuples():
        row=retained[r.point_id_b.split(':')[-1]]
        value=float(row['normed cross sections, ep->epg, exp, 10.6 GeV, unpol'].strip('()').split(',')[0])
        assert np.isclose(value,r.xs_b,rtol=1e-8,atol=0)
        assert np.isclose(float(row['t_abs_min']),r.t_abs_min_b,atol=0)
    
    metrics={}
    # Refit the earlier-experiment cross-checks without their low-|t| points.
    for key in ['defurne2015','defurne2017']:
        sample=matches[(matches.dataset_a=='jo2015')&(matches.dataset_b==key)&(matches.t_abs_a>=.25)&(matches.t_abs_b>=.25)].copy()
        fit=world.fit_pair_relative_normalization(*[sample[c].to_numpy() for c in ['xs_a_to_b_km15','xs_b','point_unc_a_to_b_km15','point_unc_b','transport_model_unc_abs']],float(sample.norm_frac_a.iloc[0]),float(sample.norm_frac_b.iloc[0]))
        frac=np.sqrt((sample.point_unc_a_to_b_km15/sample.xs_a_to_b_km15)**2+(sample.point_unc_b/sample.xs_b)**2+(sample.transport_model_unc_abs/sample.xs_a_to_b_km15)**2)
        pull=(np.log(sample.xs_a_to_b_km15/sample.xs_b)+np.log(fit['scale_a_to_b']))/frac
        metrics[key]={'N':len(sample),'scale':fit['scale_a_to_b'],'profiled_pull_rms':float(np.sqrt(np.mean(pull**2)))}
    print('Earlier-experiment cross-checks',metrics,flush=True)
    _,_,_,fitmetrics=world.run_pass2_anchor_nuisance_scenarios(anchors)
    metrics['world_fits']=json.loads(fitmetrics.to_json(orient='records'))
    model=pd.read_csv(tables/'pass2_model_publication_point_residuals.csv')
    model=model[(model.model_name=='KM15')&(model.scenario=='raw')&model.point_id.isin(ids)]
    positive_ids = {"pass2:"+key for key, row in retained.items() if number(row["normed cross sections, ep->epg, exp, 10.6 GeV, unpol"]) > 0}
    assert len(model) == model.point_id.nunique() and set(model.point_id) == positive_ids
    metrics["pass2_count"] = len(model)
    inputs=[model[c].to_numpy() for c in ['data','model','point_unc','norm_frac','corr_scale_frac']]
    fits={}
    for name,n,c in [('raw',False,False),('normalization',True,False),('correlated',True,True)]:
        result=world.fit_pass2_model_publication_nuisances(*inputs,include_norm=n,include_corr=c)
        fits[name]={k:v for k,v in result.items() if k!='point_table'}
    metrics['km15']=fits
    canonical=pd.read_csv(tables/'canonical_world_data_with_models.csv',low_memory=False)
    lee=canonical[(canonical.dataset=='lee2026')&(canonical.t_abs_min>=.25)]
    metrics['pass1_native']=world.fit_one_normalization_nuisance(lee["xs"].to_numpy(),lee.km15_native.to_numpy(),lee.point_unc_abs.to_numpy(),float(lee.norm_frac.iloc[0]))
    metrics['current_median_pct']=float(100*np.median([float(r['current dependence sys frac, 10.6 GeV']) for r in retained.values()]))
    metrics['correlated_median_pct']=float(100*np.median([float(r['correlated scale sys frac, 10.6 GeV']) for r in retained.values()]))
    for reference in ['BH','KM15']:
        points=[p for p in ratios.source_points(out/'global_comparison_current_001',reference) if p['point_id'] in ids]
        assert {p["point_id"] for p in points} == positive_ids
        samples=[(points,None,None)]
        if reference=='BH':samples.append(([p for p in points if p['phi_deg']%360<30 or p['phi_deg']%360>330],'bh_phi0',r'Pass-2 / BH for $\phi<30^\circ$ or $\phi>330^\circ$'))
        for sample,slug,title in samples:
            with tempfile.TemporaryDirectory() as tmp:
                render=ratios.render(sample,reference,Path(tmp),.31,220,slug,title)
                file=Path(render['figure']);(proposal/file.name).write_bytes(file.read_bytes())
            if slug:
                summary=ratios.summarize(sample);ts=[s for s in summary if s['variable']=='t']
                metrics['bh_phi0']={'N':len(sample),'t_ratio_medians':[float(s['ratio_median']) for s in ts],'point_ids':[p['point_id'] for p in sample]}
    print('Updated benchmarks',json.dumps({k:v for k,v in metrics.items() if k not in ['world_fits','km15','bh_phi0']},default=lambda x:int(x) if isinstance(x,np.integer) else x),flush=True)
    
    # Retain the same physical xB canvases as the original note, irrespective of page numbering.
    fullsummary=pd.read_csv(tables/'pass2_anchor_panel_summary.csv')
    oldspecs=world._world_canvas_specs_by_xb(fullsummary)
    xbkeys=[oldspecs[i]['xb_key'] for i in [0,2,4]]
    summary=world.make_pass2_anchor_panel_summary(anchors)
    chosen=summary[summary.apply(lambda r:any(np.isclose(r.xBmin_ref,lo) and np.isclose(r.xBmax_ref,hi) for lo,hi in xbkeys),axis=1)].copy()
    newspecs=world._world_canvas_specs_by_xb(chosen)
    assert [x['xb_key'] for x in newspecs]==xbkeys
    assert (chosen.t_abs_min_ref>=.25).all()
    metrics['world_figures']=[{'filename':f'pass2_world_raw_page{n:02d}.png','xB_edges':list(spec['xb_key']),'panel_count':len(spec['placements']),'t_edges':[list(k) for k in spec['t_keys']]} for n,spec in zip([1,3,5],newspecs)]
    emff=load_module(scripts/'extract_emff_from_dvcs_bh.py','note_emff')
    # Check the installed model against saved predictions before drawing new curves.
    check=anchors.drop_duplicates('published_bin_b').iloc[[0,-1]]
    for r in check.itertuples():
        result=emff.evaluate_km15_point((0,r.xB_b,r.Q2_b,r.t_abs_b,r.phi_b,r.ebeam_b,'identity'))
        assert np.isclose(result['km15_ep'],r.km15_b,rtol=1e-7,atol=0),(result['km15_ep'],r.km15_b)
    with tempfile.TemporaryDirectory() as tmp:
        world.plot_pass2_anchor_world_panels(anchors,chosen,Path(tmp),emff,{})
        for index,target in enumerate([1,3,5],1):
            (proposal/f'pass2_world_raw_page{target:02d}.png').write_bytes((Path(tmp)/f'pass2_world_raw_page{index:02d}.png').read_bytes())
    
    return metrics


def update_world_text(text, metrics):
    fits = metrics["km15"]
    def substitute(old,new):
        nonlocal text
        assert text.count(old)==1,(old,text.count(old))
        text=text.replace(old,new,1)
    for key,oldscale,oldpull in [('defurne2015','0.96','0.89'),('defurne2017','0.97','0.86')]:
        values=metrics[key]
        substitute(f'scale near ${oldscale}$',f"scale near ${values['scale']:.2f}$")
        substitute(f'pull RMS near ${oldpull}$',f"pull RMS near ${values['profiled_pull_rms']:.2f}$")
    for r,oldshift,oldchi in zip(metrics['world_fits'],['16.5','18.8','16.1','14.4'],['2.18','2.11','1.97','1.82']):
        substitute('$'+oldshift+r'\%$',f"${r['pass2_total_shift_median_pct']:.1f}\\%$")
        substitute(r'\chi^2/\mathrm{ndf}='+oldchi,rf"\chi^2/\mathrm{{ndf}}={r['chi2_per_ndf']:.2f}")
    substitute(r'$\chi^2/N=15.45$ over 1,445 points',f"$\\chi^2/N={fits['raw']['chi2_per_point']:.2f}$ over {metrics['pass2_count']:,} points")
    # Two identical chi2 values and shifts appear in separate sentences.
    text=text.replace(r'$20.1\%$',f"${fits['normalization']['normalization_shift_pct']:.1f}\\%$",1)
    text=text.replace(r'$\chi^2/N=9.78$',f"$\\chi^2/N={fits['normalization']['chi2_per_point']:.2f}$",1)
    text=text.replace(r'$\chi^2/N=9.78$',f"$\\chi^2/N={fits['correlated']['chi2_per_point']:.2f}$",1)
    text=text.replace(r'$20.1\%$',f"${fits['correlated']['total_shift_median_pct']:.1f}\\%$",1)
    substitute(r'$\chi^2/N=1.24$',f"$\\chi^2/N={metrics['pass1_native']['chi2_per_point']:.2f}$")
    substitute(r'$1.60\%$',f"${metrics['current_median_pct']:.2f}\\%$")
    substitute(r'$2.27\%$',f"${metrics['correlated_median_pct']:.2f}\\%$")
    substitute(r'approximately $20\%$',f"approximately ${fits['normalization']['normalization_shift_pct']:.0f}\\%$")
    substitute('all 1,445 supported production points',f"all {metrics['pass2_count']:,} supported points")
    substitute('the 1,289 matched points',f"the {metrics['pass1_count']:,} matched points")
    text=text.replace('277 points',str(metrics['bh_phi0']['N'])+' points')
    lo=min(metrics['bh_phi0']['t_ratio_medians']);hi=max(metrics['bh_phi0']['t_ratio_medians'])
    substitute(r'$0.72$ to $0.85$ across the six',f'${lo:.2f}$ to ${hi:.2f}$ across the four')
    return text


def update_direct_pass1(section, frame):
    tstats = [float(frame.loc[(frame.hayward_t_abs >= lo) & (frame.hayward_t_abs < hi), "hayward_over_lee_direct"].median())
              for lo, hi in zip([.25, .4, .6, .8], [.4, .6, .8, 1.])]
    ratios=frame.hayward_over_lee_direct.to_numpy();pulls=frame.direct_pull_pointwise.to_numpy()
    text=section
    old=re.search(r'The same-bin comparison identifies .*?(?=\n\n)',text,re.S)[0]
    new=(f'The same-bin comparison identifies {len(frame):,} positive matched points. The median ratio of the present cross section to pass~1 is ${np.median(ratios):.3f}$, '
         f'with a central $68\\%$ interval of ${np.percentile(ratios,16):.3f}$--${np.percentile(ratios,84):.3f}$. The arithmetic mean is ${np.mean(ratios):.3f}$. '
         f'The mean signed pull is ${np.mean(pulls):.3f}$, its median is ${np.median(pulls):.3f}$ and its RMS is ${np.sqrt(np.mean(pulls**2)):.2f}$. '
         'These summarize the matched points rather than an uncertainty on their median or a universal normalization factor. The analyses contain some overlapping experimental events and pass~1 carries a separate $31\\%$ correlated normalization uncertainty.')
    text=text.replace(old,new)
    old=re.search(r'Figure~\\ref\{fig:global_pass1_kinematics\} and Table.*?(?=\n\n)',text,re.S)[0]
    new=(r'Figure~\ref{fig:global_pass1_kinematics} and Table~\ref{tab:global_pass1_t_ratios} show the main differential trend. '
         f'The median pass-2/pass-1 ratios are ${tstats[0]:.2f}$, ${tstats[1]:.2f}$ and ${tstats[2]:.2f}$ over $0.25<|t|<0.80~\\mathrm{{GeV}}^2$ and ${tstats[3]:.2f}$ in the highest interval. '
         r'Most of the intermediate region is within roughly $10\%$ of pass~1. The highest-$|t|$ sample has fewer points and a broad ratio distribution. '
         r'These patterns should be considered alongside the pass-1 normalization allowance and the incomplete proton and photon DATA/MC efficiency ratios.')
    text=text.replace(old,new)
    old=re.search(r'In the first two \$\|t\|\$ intervals,.*?(?=\n\n)',text,re.S)[0]
    new=(r'The median pulls over $0.25<|t|<0.80~\mathrm{GeV}^2$ are between $-0.23$ and $-0.38$, while the highest interval gives $+0.29$. '
         r'Thus the central values in this region are generally below pass~1, although a typical point is less than half a combined statistical and point-to-point uncertainty below pass~1. '
         r'This is not a significance for their combined difference because the normalization and shared-event correlations have not been included in that way.')
    text=text.replace(old,new)
    text=text.replace(r'from approximately $-0.88$ below $Q^2=2~\mathrm{GeV}^2$',r'from approximately $-0.72$ below $Q^2=2~\mathrm{GeV}^2$').replace(r'to $+0.15$--$+0.35$ above',r'to $+0.16$--$+0.35$ above')
    
    for label,col,edges,fmt in [
     ('tab:global_pass1_t_ratios','hayward_t_abs',[.25,.4,.6,.8,1.],lambda x:f'{x:.2f}'),
     ('tab:global_pass1_xB_pulls','hayward_xB',[.062,.09,.118,.155,.204,.268,.357,.446,.581],lambda x:f'{x:g}'),
     ('tab:global_pass1_Q2_pulls','hayward_Q2',[1,1.5,2,3,4,6],lambda x:f'{x:g}')]:
        replacement=[]
        total=0
        for lo,hi in zip(edges[:-1],edges[1:]):
            sample=frame[(frame[col]>=lo)&(frame[col]<hi)]
            total+=len(sample)
            ratio=sample.hayward_over_lee_direct;pull=sample.direct_pull_pointwise
            cells=[fmt(lo)+'--'+fmt(hi),str(len(sample)),f'{np.median(ratio):.3f}']
            if col=='hayward_t_abs':cells.append(f'{np.percentile(ratio,16):.3f}--{np.percentile(ratio,84):.3f}')
            cells += [f'{np.mean(pull):+.2f}',f'{np.median(pull):+.2f}']
            replacement.append(' & '.join(cells)+r' \\')
        assert total == len(frame)
        for match in re.finditer(r'\\begin\{table\}.*?\\end\{table\}',text,re.S):
            if r'\label{'+label+'}' not in match[0]:continue
            table=match[0]
            table=re.sub(r'(?m)^\d[^\n]* & [^\n]*\n','',table)
            marker=list(re.finditer(r'\\hline',table))[1].end()
            table=table[:marker]+'\n'+'\n'.join(replacement)+table[marker:]
            text=text[:match.start()]+table+text[match.end():]
            break
        else:raise ValueError(label)
    
    return text


def direct_pass1_figures(frame, output, scripts):
    import matplotlib.pyplot as plt
    ratios = frame.hayward_over_lee_direct.to_numpy()
    # Keep the original scatter-plot style. Only the data population changes.
    for col,xlabel in [('hayward_t_abs',r'$|t|$ (GeV$^2$)'),('hayward_xB',r'$x_B$'),('hayward_Q2',r'$Q^2$ (GeV$^2$)')]:
        fig,ax=plt.subplots(figsize=(8.4,5.4))
        ax.scatter(frame[col],frame.hayward_over_lee_direct,s=15,alpha=.45)
        ax.axhline(1,color='black',lw=1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\sigma_{\rm Hayward}/\sigma_{\rm Lee}$ (same analysis bin)')
        ax.set_ylim(0,min(3.,max(2.,float(np.percentile(ratios,99))*1.1)))
        ax.grid(alpha=.18)
        ax.set_title('CLAS12 pass-2 Hayward / Lee pass-1: exact-bin central-value ratio')
        fig.tight_layout()
        fig.savefig(output/f'exact_bin_hayward_over_lee_vs_{col}.png',dpi=220)
        plt.close(fig)
    script=load_module(scripts/"plot_cross_section_ratio_pull_summary.py", "direct_ratios")
    points=[dict(point_id=r.hayward_point_id,phi_deg=r.hayward_phi,xB=r.hayward_xB,Q2=r.hayward_Q2,t_abs=r.hayward_t_abs,ratio=r.hayward_over_lee_direct,pull=r.direct_pull_pointwise) for r in frame.itertuples()]
    with tempfile.TemporaryDirectory() as tmp:
        result=script.render(points,'Pass-1',Path(tmp),.31,220)
        (output/'pass1_ratio_pull_kinematic_summary.png').write_bytes(Path(result['figure']).read_bytes())
    


def update_global(root, base, rows, scripts):
    import pandas as pd
    import matplotlib.pyplot as plt
    folder = root/"chapters/systematics/figures/global_comparison"
    path = root/"chapters/systematics/global_comparison.tex"
    original = path.read_text()
    text = original
    retained = {r["bin index"]: r for r in selected_rows(rows, "review")}
    exact = pd.read_csv(base/"output/global_comparison_current_001/tables/lee_hayward_exact_same_bin_points.csv")
    frame = exact[exact.hayward_point_id.str.removeprefix("pass2:").isin(retained)].copy()
    validate_comparison(frame, retained, "hayward_point_id", "hayward_xs")
    if frame.hayward_point_id.duplicated().any():
        raise ValueError("Duplicate pass-1 matches")
    a = text.index(r"\subsubsection{Direct comparison")
    b = text.index(r"\subsubsection{Laboratory-angle dependence}")
    text = text[:a] + update_direct_pass1(text[a:b], frame) + text[b:]
    direct_pass1_figures(frame, folder, scripts)
    group_masks = [frame.p_theta >= 40, frame.p_theta < 35, (frame.p_theta >= 35) & (frame.p_theta < 40)]
    groups = [frame.loc[mask, "hayward_over_lee_direct"] for mask in group_masks]
    def update_angles(paragraph):
        for old, new in [
            ("$0.809$ (1,022 points)", f"${groups[0].median():.3f}$ ({len(groups[0]):,} points)"),
            ("$1.502$ (195 points,", f"${groups[1].median():.3f}$ ({len(groups[1])} points,"),
            ("contains 72 points and has a median of $1.177$", f"contains {len(groups[2])} points and has a median of ${groups[2].median():.3f}$")]:
            paragraph = replace_once(paragraph, old, new)
        return paragraph
    text = edit_paragraph(text, "The laboratory-angle study", update_angles)
    fig, ax = plt.subplots(figsize=(8.4,5.4))
    ax.scatter(frame.p_theta, frame.hayward_over_lee_direct, s=15, alpha=.45)
    ax.axhline(1,color="black",lw=1)
    ax.axvline(35,color="black",lw=.8,ls="--")
    ax.axvline(40,color="black",lw=.8,ls=":")
    ax.set_xlabel(r"$\theta_p$ (deg)")
    ax.set_ylabel(r"$\sigma_{\rm Hayward}/\sigma_{\rm Lee}$ (same analysis bin)")
    ax.set_ylim(0,min(3.,max(2.,float(np.percentile(frame.hayward_over_lee_direct,99))*1.1)))
    ax.grid(alpha=.18)
    ax.set_title("CLAS12 pass-2 Hayward / Lee pass-1: exact-bin central-value ratio")
    fig.tight_layout()
    fig.savefig(folder/"exact_bin_hayward_over_lee_vs_p_theta.png",dpi=220)
    plt.close(fig)
    metrics = world_comparisons(root, base, {r["bin index"]:r for r in rows}, scripts)
    metrics["pass1_count"] = len(frame)
    a = text.index(r"\subsubsection{Comparison with earlier CLAS")
    b = text.index(r"\subsubsection{Normalization strategy", a)
    text = text[:a] + update_world_text(text[a:b], metrics) + text[b:]
    canonical = pd.read_csv(base/"output/global_comparison_current_001/tables/canonical_world_data_with_models.csv",low_memory=False).set_index("point_id")
    model = pd.read_csv(base/"output/global_comparison_current_001/tables/pass2_model_publication_point_residuals.csv")
    bh = model[(model.model_name=="BH") & (model.scenario=="raw") & model.point_id.isin(metrics["bh_phi0"]["point_ids"])].copy()
    bh = bh.join(canonical[["p_theta", "g_theta"]], on="point_id")
    ft, other = bh[bh.g_theta < 5], bh[bh.g_theta >= 5]
    consistency = period_consistency(selected_rows(rows, "review"))
    replacements = {
        "The strong $|t|$ and proton-angle trends": "The observed $|t|$ and proton-angle trends",
        "In that 277-point sample, 271 points": f"In that {len(bh)}-point sample, {int((bh.p_theta>40).sum())} points",
        "$0.776$ for 133 points": f"${np.median(ft.data/ft.model):.3f}$ for {len(ft)} points",
        "The remaining 144 points have a median of $0.803$": f"The remaining {len(other)} points have a median of ${np.median(other.data/other.model):.3f}$",
        r"consistency diagnostic of approximately $1.8\%$": f"consistency diagnostic of approximately ${consistency:.1f}\\%$",
        " The first two intervals remain lower, with ratios near $0.72$--$0.73$.": ""}
    for old, new in replacements.items():
        text = replace_once(text, old, new)
    # Preserve the original structure, equation, labels and all eleven graphics.
    for command in ["subsubsection", "label", "includegraphics"]:
        pattern = r"\\" + command + r"(?:\[[^\]]*\])?\{([^}]+)\}"
        if re.findall(pattern, text) != re.findall(pattern, original):
            raise ValueError(f"Global-comparison structure changed: {command}")
    path.write_text(text)
    metrics["pass1_point_ids"] = frame.hayward_point_id.tolist()
    metrics["pass1_ratios_by_id"] = dict(zip(frame.hayward_point_id, frame.hayward_over_lee_direct))
    metrics["pass1_median_ratio"] = float(frame.hayward_over_lee_direct.median())
    metrics["period_consistency_percent"] = consistency
    metrics["efficiency_subsets"] = {"N":len(bh), "high_proton_angle_N":int((bh.p_theta>40).sum()),
        "FT_like_N":len(ft), "FT_like_median":float(np.median(ft.data/ft.model)),
        "other_N":len(other), "other_median":float(np.median(other.data/other.model))}
    return metrics


def update_results(root, profile):
    xs_path = root/"chapters/results/cross_section.tex"
    xs = xs_path.read_text()
    if profile == "review":
        xs = re.sub(r"\\begin\{figure\}.*?\\end\{figure\}\s*\\clearpage",
                    lambda m: "" if re.search(r"xs_tbin_0[12]_", m[0]) else m[0], xs, flags=re.S)
        xs = replace_once(xs, r"\ref{fig:results_xs_10p6_t01}", r"\ref{fig:results_xs_10p6_t03}")
        xs = replace_once(xs,
            "In particular, the anomalous behavior in some low-$|t|$ kinematic bins requires a final kinematic correction and corresponding revalidation before the extraction can be considered complete.",
            "The lower-$|t|$ cross-section bins remain under study.")
        xs = xs.replace(" These v0 cross sections remain subject to the outstanding kinematic correction, especially in the low-$|t|$ region.", "")
        xs = replace_once(xs, "across the full acceptance; the remaining tensions", "across the full acceptance. The remaining tensions")
        for p in (root/"chapters/results/figures/v0").glob("*_xs_tbin_0[12]_*"):
            p.unlink()
        extraction = root/"chapters/bsa_extraction/bsa_determination.tex"
        extraction.write_text(replace_once(extraction.read_text(),
            r"The BSA is determined in the same four-dimensional $(x_B,Q^2,|t|,\phi)$ bins used for the cross-section measurement.",
            r"The BSA is determined on the same four-dimensional $(x_B,Q^2,|t|,\phi)$ extraction grid as the cross section and retains all six $|t|$ intervals."))
    xs_path.write_text(xs)
    # The input-run BSA figures are common to both profiles.
    bsa = root/"chapters/results/bsa.tex"
    text = replace_once(bsa.read_text(), "The measurements are shown together",
        "The BSA results cover all six $|t|$ intervals. The vertical bars show statistical uncertainties. The measurements are shown together")
    text = replace_once(text, "pending the final analysis review and associated cross-section corrections", "pending the final analysis review")
    text = replace_once(text, "is likewise striking; importantly, KM15", "is likewise striking. Importantly, KM15")
    bsa.write_text(text)


def update_abstract(root, rows, comparisons):
    samples = [summary(selected_rows(rows, "review"), label) for label in ["10.6 GeV", "Sp19 Inb"]]
    t_ratios = []
    # Values are already checked against the nominal CSV by update_global.
    for lo, hi in zip([.25, .4, .6, .8], [.4, .6, .8, 1.]):
        ids = {"pass2:"+r["bin index"] for r in rows if float(r["t_abs_min"]) == lo and float(r["t_abs_max"]) == hi}
        values = [v for key, v in comparisons["pass1_ratios_by_id"].items() if key in ids]
        t_ratios.append(float(np.median(values)))
    path = root/"main.tex"
    text = path.read_text()
    a = text.index(r"\begin{abstract}")
    b = text.index(r"\end{abstract}") + len(r"\end{abstract}")
    abstract = text[a:b]
    replacements = {
        "exxclusive channel": "exclusive channel",
        r"The four-fold differential measurements cover over 1300 kinematic points in the valence quark region, $0.06 < x_B < 0.58$, $1.00 < Q^2 < 5.76~\mathrm{GeV}^2$ and $0.11 < |t| < 1.00~\mathrm{GeV}^2$.":
            f"The four-fold differential cross sections cover {len(samples[0][0]):,} kinematic points at 10.6 GeV and {len(samples[1][0]):,} at 10.2 GeV in the valence quark region, "
            + r"$0.06 < x_B < 0.58$, $1.00 < Q^2 < 5.76~\mathrm{GeV}^2$ and $0.25 < |t| < 1.00~\mathrm{GeV}^2$. The beam-spin asymmetries retain coverage down to $|t|=0.11~\mathrm{GeV}^2$.",
        r"$9.6\%$ and $9.1\%$": f"${samples[0][1][0,0]:.1f}\\%$ and ${samples[0][1][7,0]:.1f}\\%$",
        r"are analyzed separately and median point-to-point differences are at the sub-2\% level.":
            f"are analyzed separately and give a statistical-subtracted consistency diagnostic of approximately ${comparisons['period_consistency_percent']:.1f}\\%$.",
        r"The cross sections exhibit a kinematic dependence in their agreement with pass-1, with median pass-2-to-pass-1 ratios of approximately $0.7$ in the two lowest $|t|$ intervals, approaching unity at intermediate-to-larger $|t|$ where statistics dominate.":
            f"The cross sections exhibit a kinematic dependence in their agreement with pass-1, with median pass-2-to-pass-1 ratios ranging from approximately ${min(t_ratios):.2f}$ to ${max(t_ratios):.2f}$ across the reported $|t|$ intervals. "
            + "The lower-$|t|$ cross-section bins remain under study."}
    for old, new in replacements.items():
        abstract = replace_once(abstract, old, new)
    path.write_text(text[:a]+abstract+text[b:])


def validate_comparison(frame, retained, id_column, value_column):
    base = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol"
    for point_id, value in zip(frame[id_column], frame[value_column]):
        row = retained.get(point_id.removeprefix("pass2:"))
        if row is None or not np.isclose(number(row[base]), float(value), rtol=1e-8, atol=0):
            raise ValueError(f"Comparison input differs from nominal extraction: {point_id}")


def validate_sources(root, out, rows, profile):
    if len(rows) != len({r["bin index"] for r in rows}):
        raise ValueError("Duplicate nominal bin IDs")
    samples = {}
    for label in ["10.6 GeV", "Sp19 Inb"]:
        sample, stats, values = summary(selected_rows(rows, profile), label)
        samples[label] = {"point_ids": [r["bin index"] for r in sample], "count": len(sample),
                          "uncertainty_percent": stats.tolist(), "finite_support_counts": np.isfinite(values).sum(axis=1).tolist()}
    tables = re.findall(r"\\begin\{table\}.*?\\end\{table\}", (root/"chapters/systematics/summary.tex").read_text(), re.S)
    if len(tables) != 2:
        raise ValueError("Expected two uncertainty summary tables")
    for label, table in zip(samples, tables):
        actual = [line.split(" & ")[1:4] for line in table.splitlines() if re.search(r" & [\d.]+ &", line)]
        expected = [[f"{mean:.2f}", f"{median:.2f}", f"{lo:.2f}--{hi:.2f}"] for mean, median, lo, hi in samples[label]["uncertainty_percent"]]
        if actual != expected:
            raise ValueError(f"Uncertainty table disagrees with selected input: {label}")
    figures = {}
    for image in (root/"chapters/results/figures/v0").glob("*.png"):
        _, energy, kind, _, interval, *_ = image.name.split("_")
        source = out/"final_results_by_t"/energy.replace("GeV", "_GeV")/{"bsa":"BSA", "xs":"sigma"}[kind]/"pass2_vs_KM15/values"/image.name
        if sha256(source) != sha256(image):
            raise ValueError(f"Final figure changed: {image.name}")
        if profile == "review" and kind == "xs" and int(interval) < 3:
            raise ValueError("Unreported cross-section canvas in review output")
        figures[str(image.relative_to(root))] = {"source": str(source), "sha256": sha256(image)}
    for kind, expected in [("bsa",12), ("xs",8 if profile == "review" else 12)]:
        if sum(f"_{kind}_" in name for name in figures) != expected:
            raise ValueError(f"Unexpected {kind} figure coverage")
    for tex in root.rglob("*.tex"):
        for name in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex.read_text()):
            if not (root/name).is_file():
                raise ValueError(f"Missing graphic {name} in {tex.relative_to(root)}")
    return samples, figures


def compile_note(root):
    for command in [["pdflatex","-draftmode"], ["bibtex","main"], ["pdflatex","-draftmode"], ["pdflatex"]]:
        if command[0] == "pdflatex":
            command += ["-interaction=nonstopmode", "-halt-on-error", "main.tex"]
        result = subprocess.run(command, cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
        if result.returncode:
            raise RuntimeError(result.stdout[-5000:])
    log = (root/"main.log").read_text(errors="replace")
    bad = [line for line in log.splitlines() if re.search(r"Citation .* undefined|Reference .* undefined|multiply defined|There were undefined references", line)]
    if bad:
        raise ValueError("\n".join(bad))
    return {"pages": int(re.search(r"Output written on main.pdf \((\d+) pages", log)[1]),
            "unresolved_references": 0, "unresolved_citations": 0}


def build(args):
    import matplotlib
    matplotlib.use("Agg")
    args.run = args.run.expanduser().resolve()
    args.note = args.note.expanduser().resolve()
    args.dest = args.dest.expanduser().resolve()
    args.dest.mkdir(parents=True, exist_ok=True)
    own_scripts = Path(__file__).resolve().parent
    scripts = own_scripts if (own_scripts/"compare_world_dvcs_cross_sections.py").is_file() else args.run/"external_scripts"
    with tempfile.TemporaryDirectory(dir=args.dest) as staging:
        staging = Path(staging)
        root = staging/f"overleaf_{args.profile}"
        with ZipFile(args.note) as archive:
            archive.extractall(root)
        out = args.run/"output"
        csv_path = out/"csvs/dvcs_pass2_analysis.csv"
        rows = read_rows(csv_path)
        for image in (root/"chapters/results/figures/v0").glob("*.png"):
            _, energy, kind, *_ = image.name.split("_")
            source = out/"final_results_by_t"/energy.replace("GeV","_GeV")/{"bsa":"BSA", "xs":"sigma"}[kind]/"pass2_vs_KM15/values"/image.name
            shutil.copyfile(source, image)
        audit = {}
        if args.profile == "review":
            audit["systematics"] = update_summary(root, args.run, rows)
            audit["comparisons"] = update_global(root, args.run, rows, scripts)
            update_abstract(root, rows, audit["comparisons"])
        update_results(root, args.profile)
        samples, figures = validate_sources(root, out, rows, args.profile)
        # Only technical compilation repairs accompany the approved source edits.
        main = root/"main.tex"
        main.write_text(main.read_text().replace(r"\usepackage{braket}", r"\usepackage{braket,booktabs}", 1))
        all_tex = "\n".join(p.read_text() for p in root.rglob("*.tex"))
        if r"\ref{fig:subfig_a}" not in all_tex:
            for tex in root.rglob("*.tex"):
                text = tex.read_text()
                if r"\label{fig:subfig_a}" in text:
                    tex.write_text(text.replace(r"\label{fig:subfig_a}", ""))
        # Retain graphics for inactive chapters as well as the active note.
        all_tex = "\n".join(p.read_text() for p in root.rglob("*.tex"))
        for asset in list(root.rglob("*")):
            if asset.is_file() and asset.suffix.lower() in {".png", ".pdf", ".jpg", ".jpeg"} and str(asset.relative_to(root)) not in all_tex:
                asset.unlink()
        inputs = [csv_path, args.note]
        if args.profile == "review":
            tables = out/"global_comparison_current_001/tables"
            inputs += [tables/name for name in ["lee_hayward_exact_same_bin_points.csv", "pairwise_matched_points.csv", "pass2_model_publication_point_residuals.csv", "canonical_world_data_with_models.csv", "pass2_anchor_panel_summary.csv"]]
            inputs += [scripts/name for name in ["plot_cross_section_ratio_pull_summary.py", "compare_world_dvcs_cross_sections.py", "extract_emff_from_dvcs_bh.py"]]
            inputs += [out/f"systematics/point_to_point_projections/point_to_point_systematics_vs_{axis}.csv" for axis in ["xB", "Q2", "t", "phi"]]
        original_files = {p.relative_to(root) for p in root.rglob("*") if p.is_file()}
        compilation = compile_note(root)
        pdf = staging/f"DVCS-{args.profile}.pdf"
        shutil.copyfile(root/"main.pdf", pdf)
        # Keep the deliverable source directory free of compilation byproducts.
        for file in list(root.rglob("*")):
            if file.is_file() and file.relative_to(root) not in original_files:
                file.unlink()
        meta = {"profile":args.profile, "run":str(args.run), "source_note":str(args.note),
                "selection":"nominal t_abs_min >= 0.25" if args.profile == "review" else "all nominal bins",
                "bsa_selection":"all nominal bins, input-run figures in both profiles",
                "uncertainty_order":["statistical","pi0","acceptance","radiative","bin_centering","exclusivity","fiducial","point_to_point","current","run_period","correlated","normalization","total_systematic"],
                "summary_columns":["mean","median","p16","p84"], "samples":samples,
                "final_figures":figures, "audit":audit, "compilation":compilation,
                "input_sha256":{str(p):sha256(p) for p in inputs}, "builder_sha256":sha256(Path(__file__)),
                "asset_sha256":{str(p.relative_to(root)):sha256(p) for p in root.rglob("*") if p.is_file()}}
        (root/"build-info.json").write_text(json.dumps(json_safe(meta),indent=2,allow_nan=False))
        archive_path = staging/f"DVCS-{args.profile}.zip"
        with ZipFile(archive_path,"w",ZIP_DEFLATED,compresslevel=1) as archive:
            for file in sorted(root.rglob("*")):
                if file.is_file():
                    archive.write(file,file.relative_to(root))
        destination = args.dest/root.name
        if destination.exists():
            destination.rename(staging/"previous_source")
        root.rename(destination)
        pdf.replace(args.dest/pdf.name)
        archive_path.replace(args.dest/archive_path.name)
        print(f"Built {args.profile}: {compilation['pages']} pages, clean Overleaf directory, ZIP and PDF.")


def package_single_project(args):
    """Package a review/full switch for one existing Overleaf project."""
    build=args.dest.expanduser().resolve();full=build/'overleaf_full';review=build/'overleaf_review'
    addon=build/'single_project_additions';folder=addon/'note_profiles'
    if addon.exists():shutil.rmtree(addon)
    folder.mkdir(parents=True)
    
    def balanced(text,index,opening='{',closing='}'):
        depth=0;i=index
        while i<len(text):
            if text[i]=='\\':i+=2;continue
            if text[i]==opening:depth+=1
            if text[i]==closing:
                depth-=1
                if depth==0:return i+1
            i+=1
        raise ValueError('Unbalanced TeX group')
    
    def atoms(text):
        result=[];i=0
        while i<len(text):
            start=i
            # Whole figures and tables keep conditional branches structurally complete.
            env=re.match(r'\\begin\{(figure|table|equation)\}',text[i:])
            if env:
                end=r'\end{'+env[1]+'}'
                j=text.index(end,i)+len(end);result.append(text[i:j]);i=j;continue
            if text[i]=='%':
                j=text.find('\n',i);i=len(text) if j<0 else j+1
            elif text[i]=='$':
                marker='$$' if text.startswith('$$',i) else '$';i+=len(marker)
                while i<len(text):
                    if text[i]=='\\':i+=2;continue
                    if text.startswith(marker,i):i+=len(marker);break
                    i+=1
            elif text[i]=='\\':
                m=re.match(r'\\(?:[A-Za-z@]+\*?|.)',text[i:]);i+=len(m[0])
                while i<len(text) and text[i] in '{[':
                    i=balanced(text,i,text[i],'}' if text[i]=='{' else ']')
            elif text[i]=='{':i=balanced(text,i)
            elif text[i].isspace():
                while i<len(text) and text[i].isspace():i+=1
            else:
                while i<len(text) and not text[i].isspace() and text[i] not in '\\${%':i+=1
            result.append(text[start:i])
        assert ''.join(result)==text
        return result
    
    def merge(a,b):
        aa,bb=atoms(a),atoms(b);result=[]
        for tag,i,j,k,l in SequenceMatcher(None,aa,bb,autojunk=False).get_opcodes():
            if tag=='equal':result.extend(aa[i:j])
            else:result.append(r'\NoteChoice{'+''.join(aa[i:j])+'}{'+''.join(bb[k:l])+'}')
        return ''.join(result)
    
    files=['chapters/systematics/summary.tex','chapters/systematics/global_comparison.tex','chapters/results/cross_section.tex','chapters/results/bsa.tex','chapters/bsa_extraction/bsa_determination.tex']
    for name in files:
        a=(full/name).read_text();b=(review/name).read_text()
        shared=merge(a,b).replace(r'\includegraphics',r'\NoteIncludeGraphics')
        target=folder/'shared'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_text(shared)
    for parent,children in [('chapters/systematics.tex',['chapters/systematics/summary','chapters/systematics/global_comparison']),('chapters/results.tex',['chapters/results/cross_section','chapters/results/bsa']),('chapters/bsa_extraction.tex',['chapters/bsa_extraction/bsa_determination'])]:
        text=(full/parent).read_text()
        for name in children:text=text.replace(r'\input{'+name+'}',r'\input{note_profiles/shared/'+name+'}')
        target=folder/'shared'/parent;target.parent.mkdir(parents=True,exist_ok=True);target.write_text(text)
    
    def abstract(text):
        a=text.index(r'\begin{abstract}');b=text.index(r'\end{abstract}')+len(r'\end{abstract}')
        return text[a:b]
    (folder/'shared/abstract.tex').write_text(merge(abstract((full/'main.tex').read_text()),abstract((review/'main.tex').read_text()))+'\n')
    
    # New assets live entirely under note_profiles, leaving existing graphics intact.
    assets={}
    with ZipFile(args.note.expanduser().resolve()) as original:
        for file in full.rglob('*.png'):
            name=str(file.relative_to(full))
            if original.read(name)!=file.read_bytes():
                target=folder/'common'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(file,target);assets[name]='common'
        for file in review.rglob('*.png'):
            name=str(file.relative_to(review))
            if file.read_bytes()!=(full/name).read_bytes():
                target=folder/'review'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(file,target);assets[name]='review'
    
    (folder/'selection.tex').write_text('% Change only this line to switch the entire note.\n\\ReviewNotetrue  % review\n% \\ReviewNotefalse % full\n')
    (folder/'setup.tex').write_text(r'''% Load once at the end of the preamble, before \begin{document}.
\newif\ifReviewNote
\input{note_profiles/selection}
\usepackage{booktabs}
\newcommand{\NoteChoice}[2]{\ifReviewNote#2\else#1\fi}
\newcommand{\NoteCommonGraphics}[2]{%
  \IfFileExists{note_profiles/common/#2}{%
    \includegraphics[#1]{note_profiles/common/#2}%
  }{\includegraphics[#1]{#2}}%
}
\newcommand{\NoteIncludeGraphics}[2][]{%
  \ifReviewNote
    \IfFileExists{note_profiles/review/#2}{%
      \includegraphics[#1]{note_profiles/review/#2}%
    }{\NoteCommonGraphics{#1}{#2}}%
  \else
    \NoteCommonGraphics{#1}{#2}%
  \fi
}
''')
    # A tested example main is supplied outside the additions-only folder for comparison.
    with ZipFile(args.note.expanduser().resolve()) as original:main=original.read('main.tex').decode()
    a=main.index(r'\begin{abstract}');b=main.index(r'\end{abstract}')+len(r'\end{abstract}')
    main=main[:a]+r'\input{note_profiles/shared/abstract}'+main[b:]
    for name in ['systematics','results','bsa_extraction']:
        main=main.replace(r'\input{chapters/'+name+'}',r'\input{note_profiles/shared/chapters/'+name+'}')
    assert main.count(r'\begin{document}')==1
    main=main.replace(r'\begin{document}',r'\input{note_profiles/setup}'+'\n\n'+r'\begin{document}',1)
    (build/'main_single_project_example.tex').write_text(main)
    (folder/'README.txt').write_text('''ONE EXISTING OVERLEAF PROJECT

Upload this note_profiles folder into your existing project. It adds files under a new folder and does not replace your existing chapter or figure files.

In your existing main.tex, make these five edits:

1. Immediately before \\begin{document}, add:
   \\input{note_profiles/setup}

2. Replace the complete abstract block, including \\begin{abstract} and \\end{abstract}, with:
   \\input{note_profiles/shared/abstract}

3. Replace \\input{chapters/systematics} with:
   \\input{note_profiles/shared/chapters/systematics}

4. Replace \\input{chapters/bsa_extraction} with:
   \\input{note_profiles/shared/chapters/bsa_extraction}

5. Replace \\input{chapters/results} with:
   \\input{note_profiles/shared/chapters/results}

Keep main.tex as the Overleaf main document.

SWITCHING

Open note_profiles/selection.tex. Use \\ReviewNotetrue for review or \\ReviewNotefalse for full, then recompile. You do not need Python or the analysis repository to switch in Overleaf.

EDITING

The five profile-sensitive chapter files are now edited under note_profiles/shared/chapters. Edit the abstract in note_profiles/shared/abstract.tex. Ordinary shared text appears once. \\NoteChoice{full text}{review text} contains differences between the two versions. Generated table and figure blocks can have separate branches. Edit common wording outside those branches to update both versions together. Other chapters and the bibliography remain in their existing locations.

The new figures are stored under note_profiles/common or note_profiles/review. The original figures are kept in place. Both profiles use the updated Desktop BSA figures.

This overlay was prepared from the supplied v2 archive and the approved revisions. If you have edited a profile-sensitive chapter in Overleaf since that archive, incorporate those edits into its shared replacement before switching to it. The original copies remain available for comparison.

To stop using the overlay, undo the five edits in main.tex. The original chapter and figure files have not been replaced.
''')
    archive=build/'DVCS-single-project-additions.zip'
    with ZipFile(archive,'w',ZIP_DEFLATED) as z:
        for file in sorted(addon.rglob('*')):
            if file.is_file():z.write(file,file.relative_to(addon))
    print('Additions ZIP:',archive,'MB',round(archive.stat().st_size/1024**2,2))
    print('Shared sources',len(list((folder/'shared').rglob('*.tex'))),'new assets',len(assets))
    print('Example main:',build/'main_single_project_example.tex')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", choices=["review", "full", "both"])
    desktop = Path.home()/"Desktop"
    parser.add_argument("--note", type=Path, default=desktop/"DVCS_RGA_pass2_analysis_note_v2.zip")
    parser.add_argument("--run", type=Path, default=desktop/"dvcs_cross_section")
    parser.add_argument("--dest", type=Path, default=Path(__file__).resolve().parent/"note_build")
    args = parser.parse_args()
    if args.profile == "both":
        for profile in ["full", "review"]:
            build(argparse.Namespace(**{**vars(args), "profile": profile}))
        package_single_project(args)
    else:
        build(args)
