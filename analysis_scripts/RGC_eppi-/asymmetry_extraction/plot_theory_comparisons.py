#!/usr/bin/env python3
"""Rebuild the saved pi- BM/CFF comparisons with the shared pion routine."""
from pathlib import Path
import subprocess,sys
routine=Path(__file__).resolve().parents[2]/'pion_theory_comparisons/run_comparisons.py'
if __name__=='__main__':
 raise SystemExit(subprocess.run([sys.executable,str(routine),'--channels','pim',*sys.argv[1:]]).returncode)
