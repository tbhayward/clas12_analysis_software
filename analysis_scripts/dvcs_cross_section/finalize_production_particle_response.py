#!/usr/bin/env python3
"""Associate the completed production CSV with its frozen response and final cuts."""
import argparse,hashlib,json,os
from pathlib import Path

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def finalize(origin,csv_path,cuts):
    origin=Path(origin).resolve();csv_path=Path(csv_path).resolve();cuts=Path(cuts).resolve();manifest=json.loads((origin/'input_manifest.json').read_text())
    result=dict(model_id=manifest['model_id'],csv=str(csv_path),csv_sha256=digest(csv_path),final_cuts=str(cuts),final_cuts_sha256=digest(cuts),response_manifest=str(origin/'input_manifest.json'),response_manifest_sha256=digest(origin/'input_manifest.json'),response_model=manifest['model'],data_and_generated_kinematics_modified=False,reconstructed_mc_modified=True)
    for target in (origin/'extraction_manifest.json',csv_path.with_suffix('.response.json')):
        temp=target.with_name(target.name+'.tmp');temp.write_text(json.dumps(result,indent=2)+'\n');os.replace(temp,target)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--response-output',required=True,type=Path);p.add_argument('--csv',required=True,type=Path);p.add_argument('--cuts',required=True,type=Path);a=p.parse_args();finalize(a.response_output,a.csv,a.cuts)
