#!/usr/bin/env python3
"""Small compiled checks for runtime bin overrides; no ROOT files required."""
import csv
from pathlib import Path
import subprocess
import tempfile

SOURCE = Path(__file__).resolve().parent

def main():
    with tempfile.TemporaryDirectory() as temp:
        d = Path(temp)
        harness = d / 'check.cpp'
        harness.write_text(r'''
#include "sebastian_tgamma_binning.h"
#include "initialize_pass2_csv.h"
#include <cassert>
#include <cmath>
int main(int argc,char** argv) {
 assert(analysis_bin_t(-.2,1,.2,.3,2,.4,.5,10.2)==.2);
 analysis_photon_t=true;
 double a=analysis_bin_t(-.2,3,.2,.3,4,.1,.5,10.2);
 double b=analysis_bin_t(-.2,3,.2,.3,4,.1,.5,10.604);
 assert(std::isfinite(a) && std::abs(a-b)>1e-6);
 assert(a==sebastian_tgamma_abs(3,.2,.3,4,.1,.5,10.2));
 assert(!std::isfinite(analysis_bin_t(-.2,0,.2,.3,4,.1,.5,10.2)));
 analysis_photon_t=false;
 analysis_merge_low_t=argc>3;
 return initialize_pass2_csv(argv[1],argv[2]) ? 0:1;
}
''')
        exe = d / 'check'
        subprocess.run(['c++', '-std=c++17', '-I'+str(SOURCE), str(harness),
                        str(SOURCE/'initialize_pass2_csv.cpp'), '-o', str(exe)], check=True)
        original = SOURCE/'imports/all_bin_v3.csv'
        for mode in ('nominal','merged'):
            subprocess.run([str(exe),str(original),str(d/f'{mode}.csv')]+
                           (['merge'] if mode=='merged' else []),check=True)
        read=lambda p:list(csv.DictReader(p.open()))
        valid=[r for r in read(original) if int(r['valid bin'])==1]
        nominal=read(d/'nominal.csv'); merged=read(d/'merged.csv')
        assert len(nominal)==len(valid)
        for a,b in zip(valid,nominal):
            for k in ('xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max','phimin','phimax'):
                assert a[k]==b[k]
        key=lambda r:tuple(float(r[k]) for k in ('xBmin','xBmax','Q2min','Q2max','phimin','phimax'))
        low={key(r):r for r in valid if float(r['t_abs_min'])==.11}
        high={key(r):r for r in valid if float(r['t_abs_min'])==.15}
        cells=[r for r in merged if float(r['t_abs_min'])==.11]
        assert len(cells)==len(low.keys()|high.keys())
        assert all(float(r['t_abs_max'])==.25 for r in cells)
        assert not any(float(r['t_abs_min'])==.15 for r in merged)
        assert len(merged)==len(nominal)-len(low.keys()&high.keys())
        for r in cells:
            donor=high.get(key(r),low.get(key(r)))
            # Cross-section extraction imports Frad/Fbin by this donor index.
            assert r['bin index']==donor['']
        print('PASS: nominal bins, merged union, donor factors, photon definition and beam energy')

if __name__=='__main__':
    main()
