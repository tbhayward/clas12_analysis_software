"""Render the production note solely from small saved histograms and tables."""
import csv,html,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import study_particle_response as base
CORE=(0,3,4,5)
COLORS=dict(DATA='#205cab',nominal='#888888',modified='#22946b')
LABELS={'Delta_phi':r'$\Delta\phi$ (rad)','pTmiss':r'$p_T^{\mathrm{miss}}$ (GeV)','Emiss2':r'$E_{\mathrm{miss}}$ (GeV)','Mx2':r'$M_X^2(ep\gamma)$ (GeV$^2$)','Mx2_2':r'$M_X^2(e\gamma)$ (GeV$^2$)','m_gamma_gamma':r'$M_{\gamma\gamma}$ (GeV)'}
REGIONS=['Electron FD','Proton FD','Proton CD','Photon FD','Photon FT']
COORDS=('p_sigma_fraction','theta_sigma_rad','phi_sigma_rad','p_bias_fraction')

def render(out):
    out=Path(out);params=list(csv.DictReader((out/'parameters.csv').open()));scores=list(csv.DictReader((out/'shape_scores.csv').open()));links=[]
    plt.rcParams.update({'font.family':'serif','font.size':11,'axes.labelsize':12,'axes.titlesize':12,'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','pdf.fonttype':42})
    def save(fig,folder,name,pdf):
        folder.mkdir(parents=True,exist_ok=True);path=folder/(name+'.pdf');fig.savefig(path,bbox_inches='tight');fig.savefig(path.with_suffix('.png'),dpi=220,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig);links.append(path.relative_to(out).as_posix())
    with np.load(out/'plot_histograms.npz',allow_pickle=False) as h,PdfPages(out/'particle_response_figures.pdf') as pdf:
      for period in dict.fromkeys(r['period'] for r in params):
        label=period.replace('_',' ').title();folder=out/'plots'/period;records=[r for r in scores if r['period']==period];par=np.zeros((5,4))
        for r in params:
            if r['period']==period:par[base.REGIONS.index(r['particle']),COORDS.index(r['coordinate'])]=float(r['value'])
        fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
        for j,ax in enumerate(axes.flat):
            ax.bar(range(5),par[:,j]*(100 if j in (0,3) else 1000),color=['#537ba6','#b87a3f','#8b5c99','#48937e','#c45c64']);ax.axhline(0,color='black',lw=.7);ax.set(xticks=range(5),xticklabels=REGIONS,ylabel=['Added momentum width (%)','Added polar-angle width (mrad)','Added azimuth width (mrad)','Momentum bias (%)'][j]);ax.tick_params(axis='x',rotation=25);ax.grid(axis='y',alpha=.2)
        fig.suptitle(label+' — frozen reconstructed-MC response');save(fig,folder,'parameters',pdf)
        support={(r['channel'],r['topology'],r['variable'],int(r['half'])) for r in records if r['scenario']=='nominal' and r['status']=='supported'}
        def average(scenario,channel,half,reserved=False):
            names=['Mx2_2'] if reserved else [base.OBS[j] for j in CORE]+(['m_gamma_gamma'] if channel=='epgammagamma' else []);values=[float(r['score']) for r in records if r['scenario']==scenario and r['channel']==channel and int(r['half'])==half and r['variable'] in names and (r['channel'],r['topology'],r['variable'],int(r['half'])) in support];return float(np.mean(values)) if values else np.nan
        fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
        for ax,ch in zip(axes,('epgamma','epgammagamma')):
            for i,(scenario,color) in enumerate((('nominal',COLORS['nominal']),('modified',COLORS['modified']))):ax.bar(np.arange(3)+(i-.5)*.35,[average(scenario,ch,0),average(scenario,ch,1),average(scenario,ch,1,True)],.35,label=scenario,color=color)
            ax.set(xticks=range(3),xticklabels=['Training variables\nMC half 0','Training variables\nMC half 1','Reserved mass\nMC half 1'],ylabel='Shape discrepancy score',title=r'$ep\gamma$' if ch=='epgamma' else r'$ep\gamma\gamma$');ax.legend();ax.grid(axis='y',alpha=.2)
        fig.suptitle(label+' — final common DATA windows; frozen response');save(fig,folder,'validation',pdf)
        fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained');baseline=.5*(average('modified','epgamma',0)+average('modified','epgammagamma',0))
        for j,ax in enumerate(axes.flat):
            effect=[.5*(average('without_'+base.PARAMS[4*r+j],'epgamma',0)+average('without_'+base.PARAMS[4*r+j],'epgammagamma',0))-baseline for r in range(5)];ax.bar(range(5),effect,color='#537ba6');ax.axhline(0,color='black',lw=.7);ax.set(xticks=range(5),xticklabels=REGIONS,title=['Momentum width','Polar-angle width','Azimuth width','Momentum bias'][j],ylabel='Score increase when omitted');ax.tick_params(axis='x',rotation=25);ax.grid(axis='y',alpha=.2)
        fig.suptitle(label+' — conditional influence of each fitted term\nOther parameters fixed; positive means this term improves the joint score');save(fig,folder,'parameter_influence',pdf)
        for ch in ('epgamma','epgammagamma'):
          for g,topology in enumerate(base.TOPS):
            fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
            for panel,(ax,j) in enumerate(zip(axes.flat,CORE)):
                key=f'{period}_{ch}_{g}_{j}_';edges=h[f'edges_{j}'];d=h[key+'DATA'];centers=(edges[1:]+edges[:-1])/2;width=np.diff(edges);ds=d.sum();curves=[h[key+n] for n in ('nominal','modified')]
                if ds>0:ax.errorbar(centers,d[1:-1]/ds/width,yerr=np.sqrt(d[1:-1])/ds/width,fmt='o',ms=2,color=COLORS['DATA'],elinewidth=.6,label='DATA')
                for name,m in zip(('nominal','modified'),curves):
                    total=m.sum()
                    if total>0 and np.all(np.isfinite(m)):
                        ax.stairs(m[1:-1]/total/width,edges,color=COLORS[name],lw=1.6,label=('Nominal' if name=='nominal' else 'Modified')+(' signal + background' if ch=='epgamma' else ' MC'))
                        if name=='modified' and ch=='epgamma':ax.stairs(h[key+name+'_signal'][1:-1]/total/width,edges,color=COLORS[name],lw=1,ls='--',label='Modified signal');ax.stairs(h[key+name+'_background'][1:-1]/total/width,edges,color=COLORS[name],lw=1,ls=':',label=r'Modified $\pi^0$ background')
                    else:ax.text(.5,.5,'Insufficient/invalid MC normalization',transform=ax.transAxes,ha='center',fontsize=8)
                if j==5:ax.set_xlim(*base.mass_plot_limits(edges,[d]+curves))
                else:
                    q=[base.diag.quantile(v[1:-1],edges,z) for v in [d]+curves if v.sum()>0 and np.all(np.isfinite(v)) for z in (.01,.99)];q=[x for x in q if np.isfinite(x)]
                    if q:lo,hi=min(q),max(q);pad=max(.005,(hi-lo)*.1);ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                xlabel=r'$M_X^2(ep\pi^0)$ (GeV$^2$)' if ch=='epgammagamma' and j==5 else LABELS[base.OBS[j]];ax.set(xlabel=xlabel,ylabel='Unit-area density');ax.text(.97,.96,'('+chr(97+panel)+')',transform=ax.transAxes,va='top',ha='right');ax.grid(alpha=.15)
                if panel==0:ax.legend(fontsize=8)
            fig.suptitle(label+' — '+topology.replace('_','–')+' — '+(r'$ep\gamma$' if ch=='epgamma' else r'$ep\gamma\gamma$')+'\nCommon final DATA selection; displayed exclusivity cut omitted');save(fig,folder/('dvcs' if ch=='epgamma' else 'two_photon'),'exclusivity_'+topology,pdf)
        fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout='constrained');edges=h['edges_mass'];width=np.diff(edges)
        for g,ax in enumerate(axes):
            for name,color in [('DATA',COLORS['DATA']),('nominal',COLORS['nominal']),('modified',COLORS['modified'])]:
                values=h[f'{period}_mass_{g}_{name}']
                if values.sum()>0:ax.stairs(values[1:-1]/values.sum()/width,edges,label=name,color=color,lw=1.4)
            ax.set(title=base.TOPS[g].replace('_','–'),xlabel=LABELS['m_gamma_gamma'],ylabel='Unit-area density');ax.legend();ax.grid(alpha=.15)
        fig.suptitle(label+' — two-photon mass; common final selection');save(fig,folder/'two_photon','mass',pdf)
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><style>body{font:17px system-ui;max-width:1200px;margin:30px auto}img{max-width:100%}</style><h1>Production particle-response note figures</h1><p><a href="particle_response_figures.pdf">Combined vector PDF</a> · <a href="parameters.csv">Parameters</a> · <a href="shape_scores.csv">Scores and support flags</a></p><p>DATA coordinates and generated truth are unchanged. Frozen response modifies full-statistics reconstructed inputs before selection. Note figures use a bounded paired MC sample with saved probabilities. Common final DATA windows define N−1 shapes; production yields use their actual separate DATA/MC windows. Both photons share detector response laws.</p><p>Control labels follow the existing production photon-1 convention, identically in the cut finder and counting. They are not the diagnostic unordered groups. Unresolved controls retain their original vectors and are tabulated. Background transfer here uses raw selected control counts and the saved exposure ratio, not production current/proton weights. The fitted exploratory background nuisance is not promoted to a production correction. MC halves share the same DATA; scores are not formal chi-square probabilities.</p><p>Parameter-influence bars omit one term at a time without refitting; correlations mean effects are not additive. Parameters have no fitted uncertainty estimate. Normalization follows selected populations per scenario, and projections are shown as unit-area shapes (including histogram under/overflow). The remaining two-photon mass mismatch is displayed explicitly. Original processing thresholds and the control mass window restrict available events; this stage cannot restore events absent from the trees.</p>'+''.join('<p><a href="'+x+'">'+html.escape(x)+'</a><br><img loading="lazy" src="'+x.replace('.pdf','.png')+'"></p>' for x in links));base.dump(out/'render_provenance.json',dict(renderer_sha256=base.digest(Path(__file__)),histograms_sha256=base.digest(out/'plot_histograms.npz')))
    print('[response note] '+str(out/'REPORT.html'),flush=True)
