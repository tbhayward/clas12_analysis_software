#!/usr/bin/env python3
"""Imported Fbin coordinates must control XS points, models and cache identity."""
import csv
import hashlib
import sys
import tempfile
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE/'external_scripts'))
import plot_pass2_final_results_by_t as plot


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    assert plot.make_parser().parse_args(['input.csv']).reference_kinematics == 'imported'
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary); imports = root/'imports'; imports.mkdir()
        output = root/'output/csvs'; output.mkdir(parents=True)
        reference = imports/'all_bin_v3.csv'; data = output/'analysis.csv'
        donors = []; rows = []
        for ident,phi0,phi1,phi_ref,phi_measured in [(11,0,15,6,8),(12,15,30,21,23)]:
            box = dict(xBmin=.204,xBmax=.268,Q2min=1.912,Q2max=2.510,t_abs_min=.25,t_abs_max=.4,phimin=phi0,phimax=phi1)
            donors.append({'':ident,**box,'xBavg':.23,'Q2avg':2.2,'t_abs_avg':.32,'phiavg':phi_ref,'Fbin':1.05})
            rows.append({'bin index':ident,'valid bin':1,**box,
                'Fbin, 10.6 GeV':1.05,'Fbin, 10.2 GeV':1.05,
                'normed cross sections, ep->epg, exp, 10.6 GeV, unpol':'(10, 1, 0)',
                'normed cross sections, ep->epg, exp, Sp19 Inb, unpol':'(8, 1, 0)',
                plot.base.PASS2_BSA_TARGET_COLUMNS['10.6 GeV']:'(.1, .02, 0)',
                plot.base.PASS2_BSA_TARGET_COLUMNS['Sp19 Inb']:'(.1, .03, 0)',
                **{f'{name}, {label}':value for label in ['10.6 GeV','Sp19 Inb']
                   for name,value in [('xBavg',.24),('Q2avg',2.3),('t_abs_avg',.33),('phiavg',phi_measured)]}})
        write_csv(reference,donors); write_csv(data,rows)
        before = hashlib.sha256(data.read_bytes()).hexdigest()
        args = plot.make_parser().parse_args([str(data)])
        panels = plot.parse_pass2(args)
        for energy,cells in panels.items():
            key,item = next(iter(cells.items()))
            assert [p[0] for p in item['xs']] == [6,21]
            assert [p[1:] for p in item['xs']] == ([(10,1)]*2 if energy==10.6 else [(8,1)]*2)
            assert [p[0] for p in item['bsa']] == [8,23]
            assert plot.curve_reference(args,energy,key,'xs') == (.23,2.2,.32)
            assert np.allclose(plot.curve_reference(args,energy,key,'bsa'),(.24,2.3,.33))
            refs = args._model_references[(energy,key,'xs')]
            assert all(p[:3] == (.23,2.2,.32) for p in refs)
            assert all(p[4] == (10.6 if energy==10.6 else 10.1998) for p in refs)
            args._active_energy = energy
            new_key = plot.model_cache_key(key,args,'xs')
            old_key = plot.model_cache_key(key,args,'bsa')
            args.phi_dense = 48
            assert plot.model_cache_key(key,args,'xs') == new_key
            args.reference_kinematics = 'means'
            assert plot.model_cache_key(key,args,'bsa') == old_key
            assert plot.model_cache_key(key,args,'xs') != new_key
            args.reference_kinematics = 'imported'
        assert hashlib.sha256(data.read_bytes()).hexdigest() == before
        # Evaluate both observables in one cell: display curves must use their
        # own anchors, exact pulls must never interpolate a different reference.
        energy = 10.6; key,item = next(iter(panels[energy].items()))
        cfg = plot.base.ModelConfig(energy,'unused','unused',sys.executable,24,False,False,False)
        requests = []
        old_many,old_bh = plot.base.km15_xs_many,plot.base.bh_xs
        def model(points, config):
            requests.extend(points)
            return [1+x+q+t+phi/100+h/10 for x,q,t,phi,h in points]
        plot.base.km15_xs_many = model
        plot.base.bh_xs = lambda x,q,t,phi,cfg: x+q+t+phi/100
        references = {kind:args._model_references[(energy,key,kind)] for kind in ('xs','bsa')}
        references['__display__'] = {kind:plot.curve_reference(args,energy,key,kind) for kind in ('xs','bsa')}
        references['__modes__'] = dict(xs='imported_centering',bsa='measured_means')
        try:
            _,curves = plot.compute_cell_models((key,cfg,('xs','bsa'),references))
            assert curves['xs']['reference'] == 'imported_centering'
            assert len(curves['xs']['display']['phi']) == 24
            for x,q,t,phi,h in requests:
                assert np.allclose((x,q,t),(.23,2.2,.32) if h==0 else (.24,2.3,.33))
            expected = 1+.23+2.2+.32+.06
            assert abs(plot.interpolate_model(curves['xs'],'km15',6)-expected)<1e-12
            assert np.isnan(plot.interpolate_model(curves['xs'],'km15',8))
            pulls = plot.km15_pulls(item['xs'],curves['xs'])
            assert abs(pulls[0][1]-(10-expected))<1e-12
            args._active_energy = energy
            plot.draw_value_canvas([key],(.25,.4),'xs',energy,{key:item},None,
                                   {plot.model_cache_key(key,args,'xs'):curves['xs']},
                                   args,root/'smoke.png')
            assert (root/'smoke.png').is_file()
        finally:
            plot.base.km15_xs_many,plot.base.bh_xs = old_many,old_bh
        # Wrong/missing donors and independently changed Fbin must fail clearly.
        for bad in ['factor','id','geometry','coordinate']:
            changed = [dict(r) for r in rows]
            if bad=='factor':changed[0]['Fbin, 10.6 GeV']=.8
            if bad=='id':changed[0]['bin index']=999
            if bad=='geometry':changed[0]['Q2min']=1.8
            write_csv(data,changed)
            ds = [dict(r) for r in donors]
            if bad=='coordinate':ds[0]['Q2avg']=''
            write_csv(reference,ds)
            try:plot.parse_pass2(plot.make_parser().parse_args([str(data)]))
            except ValueError:pass
            else:raise AssertionError(f'Accepted inconsistent {bad}')
        write_csv(data,rows);write_csv(reference,donors)
        reference.unlink()
        try:plot.parse_pass2(plot.make_parser().parse_args([str(data),'--centering-reference-csv',str(reference)]))
        except FileNotFoundError:pass
        else:raise AssertionError('Missing explicit reference table silently replaced')
        means = plot.parse_pass2(plot.make_parser().parse_args([str(data),'--reference-kinematics','means']))
        assert next(iter(means[10.6].values()))['xs'][0][0]==8
    print('PASS: imported XS references/phi, unchanged yields/errors, separate BSA reference, exact pulls, dense curves, cache isolation/reuse, donor validation, read-only inputs')


if __name__ == '__main__':
    main()
