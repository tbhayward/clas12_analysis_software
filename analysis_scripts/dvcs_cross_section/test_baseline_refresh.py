#!/usr/bin/env python3
"""Safe reset/cache preservation and bounded pull-display regressions; no ROOT."""
import importlib.util,json,sys,tempfile
from pathlib import Path
import numpy as np
SOURCE=Path(__file__).resolve().parent
sys.path.insert(0,str(SOURCE/'external_scripts'))
import plot_pass2_final_results_by_t as plot
import rerun_nominal_baseline as runner

def main():
    values=np.array([-10000.,-10.,-2.,0.,10.,100000.])
    display,outside=plot.pull_display_points(values,10.)
    assert outside.tolist()==[True,False,False,False,False,True]
    assert max(abs(display))<=10 and np.array_equal(values,[-10000.,-10.,-2.,0.,10.,100000.])
    assert plot.make_parser().parse_args(['input.csv']).pull_limit==10
    for bad in ['0','-1','nan','inf']:
        try:plot.positive_finite(bad)
        except Exception:pass
        else:raise AssertionError(bad)
    with tempfile.TemporaryDirectory() as tmp:
        args=plot.make_parser().parse_args(['input.csv','--reference-kinematics','midpoints'])
        args._active_energy=10.6
        k=plot.base.BinKey(.118,.155,1.912,2.510,.11,.15)
        panels={k:{'xs':[(30.,-1000.,1.),(90.,2.,1.),(150.,1000.,1.)]}}
        cache={plot.model_cache_key(k,args,'xs'):{'phi':[30.,90.,150.],'km15':[0.,0.,0.]}}
        close=plot.plt.close;held=[];plot.plt.close=lambda fig:held.append(fig)
        try:
            pulls=plot.draw_pull_canvas([k],(.11,.15),'xs',10.6,panels,cache,args,Path(tmp)/'pulls.png','km15')
            assert pulls==[-1000.,2.,1000.]
            assert held[-1].axes[0].get_ylim()==(-10.,10.)
        finally:
            plot.plt.close=close
            for fig in held:close(fig)
    saved_source=runner.SOURCE
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);src=root/'source';src.mkdir();runs=root/'test_runs'
        for name in ('external_scripts','external_pi0_normalization','imports'):(src/name).mkdir()
        (src/'dvcs_analysis').write_text('#!/usr/bin/env python3\nimport json,sys\nfrom pathlib import Path\nPath("invoked.json").write_text(json.dumps(sys.argv[1:]))\n')
        (src/'dvcs_analysis').chmod(0o755);(src/'main_systematics').write_text('unused')
        (src/'external_scripts/plot_pass2_final_results_by_t.py').write_text('from pathlib import Path\nimport json,sys\nPath("plotted").write_text(json.dumps(sys.argv[1:]))\n')
        for name in runner.RUN_NAMES:
            p=runs/name/'output/final_results_by_t/10.6_GeV';p.mkdir(parents=True)
            (p/'model_curves_cache_10.6GeV.json').write_text(json.dumps({name:{'cached':True}}))
            (runs/name/'obsolete.csv').write_text('old')
        (runs/'unrelated').mkdir();(runs/'unrelated/keep').write_text('untouched')
        runner.SOURCE=src
        assert runner.main(['--test-root',str(runs),'--reset-test-runs','--plot'])==0
        nominal=runs/'proton_t_default_001'
        assert json.loads((nominal/'invoked.json').read_text())==['--skip-systematics']
        assert (nominal/'plotted').exists()
        assert '--no-pass1' not in json.loads((nominal/'plotted').read_text())
        assert not (nominal/'obsolete.csv').exists()
        assert not (runs/'merged_t_001').exists() and not (runs/'photon_t_001').exists()
        cache=json.loads((nominal/'output/final_results_by_t/10.6_GeV/model_curves_cache_10.6GeV.json').read_text())
        assert set(cache)==set(runner.RUN_NAMES)
        assert (runs/'theory_cache_preserved/model_curves_cache_10.6GeV.json').is_file()
        assert (runs/'unrelated/keep').read_text()=='untouched'
        assert not (runs/'.nominal_baseline_reset_lock').exists()
        # An unrequested overwrite must be refused, keeping the result intact.
        try:runner.main(['--test-root',str(runs)])
        except SystemExit:pass
        else:raise AssertionError('Existing run overwritten')
        assert (nominal/'plotted').exists()
        assert '--no-pass1' not in json.loads((nominal/'plotted').read_text())
    runner.SOURCE=saved_source
    print('PASS: bounded display, original pull values, nominal flags, cache preservation, exact reset scope, overwrite refusal')

if __name__=='__main__':main()
