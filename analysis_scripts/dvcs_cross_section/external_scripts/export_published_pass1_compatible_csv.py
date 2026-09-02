#!/usr/bin/env python3
import argparse
from pathlib import Path
from pass1_published_loader import load_published_pass1_dataframe

def main():
    p=argparse.ArgumentParser()
    p.add_argument('published',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--legacy-boundaries',type=Path,default=None)
    a=p.parse_args()
    df=load_published_pass1_dataframe(a.published,a.legacy_boundaries)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(a.output,index=False)
    print(f'[pass1-published] Wrote compatibility CSV: {a.output}')
if __name__=='__main__': main()
