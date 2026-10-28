#!/usr/bin/env python3
"""Single-file pilot, not a production launcher. Run from CLAS12 checkout root.
Example: python3 processing_scripts/rebuild_pilot/pilot_paired.py --mode dvcs_mc --hipo /path/file.hipo --period fa18_inb --energy 10.604 --events 10000
"""
import argparse, os, subprocess, tempfile, shutil, sys, shlex
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['dvcs_mc','pi0_mc'],required=True);p.add_argument('--hipo',required=True);p.add_argument('--period',required=True,choices=['fa18_inb','fa18_out','sp18_inb','sp18_out','sp19_inb']);p.add_argument('--energy',type=float,required=True);p.add_argument('--events',type=int,default=10000);p.add_argument('--outdir',default='rebuild_pilot_output');a=p.parse_args()
root=Path.cwd(); scripts=root/'processing_scripts/rebuild_pilot'; source=Path(a.hipo).resolve();out=Path(a.outdir).resolve();out.mkdir(parents=True,exist_ok=True)
if not source.is_file():sys.exit('Input HIPO does not exist: '+str(source))
jar=root/'processing_classes/dist/processing_classes.jar'; assert jar.is_file(),jar
if not shutil.which('root-config'):sys.exit('root-config not on PATH')
name='processing_mc_rec_dvcs.groovy' if a.mode=='dvcs_mc' else 'processing_mc_rec_exclusive_pi0.groovy'
converter=scripts/'convert_paired_stream'; cpp=scripts/'convert_paired_stream.cpp'
flags=subprocess.check_output(['root-config','--cflags'],text=True).split();libs=subprocess.check_output(['root-config','--libs'],text=True).split()
subprocess.run(['g++','-O3','-std=c++17',*flags,str(cpp),'-o',str(converter),*libs],check=True)
prefix=out/(a.mode+'_'+a.period+'_'+source.stem)
fifo=Path(str(prefix)+'.txt.events.tsv');result=Path(str(prefix)+'.mc_events.root');selected=Path(str(prefix)+'.txt')
if any(f.exists() for f in [fifo,result,selected]):sys.exit('Refusing to overwrite existing pilot outputs: '+str(prefix))
os.mkfifo(fifo)
log=Path(str(prefix)+'.processing.log'); convlog=Path(str(prefix)+'.converter.log')
try:
 with log.open('w') as gl, convlog.open('w') as cl:
  converter_proc=subprocess.Popen([str(converter),str(fifo),str(result)],stdout=cl,stderr=subprocess.STDOUT)
  args=[str(root/'coatjava/bin/run-groovy'),'-cp',str(jar),str(scripts/name),str(source),str(selected),'1',str(a.energy),'11','1',str(a.events),a.period]
  with gl:
   # Load QADB in the worker's own tcsh environment, just as processing.csh does.
   # No manual `module load` is required in the user's interactive shell.
   quoted=' '.join(shlex.quote(x) for x in args)
   setup='source source_file.txt; if ($status != 0) exit 90; module load qadb/3.4.1; if ($status != 0) exit 91; exec '+quoted
   processor_proc=subprocess.Popen(['/bin/tcsh','-c',setup],stdout=gl,stderr=subprocess.STDOUT,cwd=root)
   code=processor_proc.wait()
  # If Groovy fails before opening FIFO, converter would otherwise wait forever.
  if code:converter_proc.terminate()
  try:convert_code=converter_proc.wait(timeout=120)
  except subprocess.TimeoutExpired:converter_proc.kill();convert_code=converter_proc.wait()
  if code or convert_code:sys.exit(f'Pilot FAILED: Groovy={code}, converter={convert_code}. Logs: {log}, {convlog}')
  print('Pilot ROOT:',result)
  print('Groovy log:',log)
  print('Converter log:',convlog)
  print('IMPORTANT: Validate tree entry count, generated/reconstructed values and parity before production.')
finally:
 if fifo.exists():fifo.unlink()
