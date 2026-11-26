#!/usr/bin/env python3
"""Repeat this study using the selected samples already on this Mac."""
from pathlib import Path
import os,sys,subprocess,importlib.util
HERE=Path(__file__).resolve().parent
if importlib.util.find_spec('uproot') is None:
 deps=Path('/private/tmp/eppi_diagnostics_deps')
 if deps.exists():os.environ['PYTHONPATH']=str(deps)+os.pathsep+os.environ.get('PYTHONPATH','')
repo=Path('/Users/tbhayward/Documents/Research.nosync/clas12_analysis_software/analysis_scripts/RGC_eppi-')
desktop=Path('/Users/tbhayward/Desktop')
command=[sys.executable,str(HERE/'run_study.py'),'--analysis-root',str(repo),'--selected-inputs',str(desktop/'channel_selection/output/selected_inputs.json'),'--coarse-dilution',str(desktop/'asymmetry_extraction/output_two_tprime_bins/coarse_dilution.json'),'--theory',str(desktop/'pion_theory_comparisons/data/theory'),'--baseline-fine',str(desktop/'asymmetry_extraction/output_fixed_UU'),'--baseline-coarse',str(desktop/'asymmetry_extraction/output_two_tprime_bins'),'--workers','3',*sys.argv[1:]]
if __name__=='__main__':raise SystemExit(subprocess.run(command).returncode)
