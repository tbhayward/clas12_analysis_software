#!/usr/bin/env python3
"""One-command plotting routine, optionally refreshing inputs from new fits."""
import argparse,json,subprocess,sys,os,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'pion_theory_mpl'))

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--output',type=Path,default=HERE/'output')
 p.add_argument('--measurements',type=Path,default=HERE/'measurements/measurements.csv')
 p.add_argument('--pip-results',type=Path)
 p.add_argument('--pim-fine',type=Path)
 p.add_argument('--pim-coarse',type=Path)
 p.add_argument('--theory',type=Path,default=HERE/'inputs')
 p.add_argument('--watermark-alpha',type=float,default=.10)
 p.add_argument('--dpi',type=int,default=300)
 p.add_argument('--skip-legacy',action='store_true')
 p.add_argument('--channels',nargs='+',choices=['pip','pim'],default=['pip','pim'])
 a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 if any([a.pip_results,a.pim_fine,a.pim_coarse]):
  defaults=json.loads((HERE/'data_sources.json').read_text())
  from prepare_measurements import prepare
  snapshot=a.output/'data/measurements'
  prepare(a.pip_results or HERE/'measurements/pip_source.csv',a.pim_fine or defaults['pim_fine'],a.pim_coarse or defaults['pim_coarse'],snapshot)
 else:snapshot=a.measurements.parent
 table=snapshot/'measurements.csv' if any([a.pip_results,a.pim_fine,a.pim_coarse]) else a.measurements
 command=[sys.executable,str(HERE/'plot_comparisons.py'),'--measurements',str(table),'--theory',str(a.theory),'--output',str(a.output),'--watermark-alpha',str(a.watermark_alpha),'--dpi',str(a.dpi),'--channels',*a.channels]
 subprocess.run(command,check=True)
 if not a.skip_legacy and 'pip' in a.channels:
  subprocess.run([sys.executable,str(HERE/'legacy_update/plot_kroll_comparisons.py'),'--input-csv',str(snapshot/'pip_source.csv'),'--cff-directory',str(a.theory/'pip'),'--output-dir',str(a.output/'.legacy_work'),'--watermark-alpha',str(a.watermark_alpha),'--dpi',str(a.dpi)],check=True)
 from chiral_sectors import run as run_chiral_sectors
 run_chiral_sectors(a.theory,table,a.output,a.dpi,a.watermark_alpha,a.channels)
 from analyze_cff_sensitivity import run as analyze_cffs
 analyze_cffs(a.theory,table,a.output)
 from summarize_inputs import summarize,extend_gallery
 summarize(a.theory,table,a.output)
 extend_gallery(a.output)
 print('Open '+str(a.output/'index.html'))
if __name__=='__main__':main()
