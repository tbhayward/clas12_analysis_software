"""Replace only pi-minus measurements from the validated conference release table."""
import csv,hashlib,json,shutil
from pathlib import Path
import numpy as np


def prepare_release(pip_source,release,output):
 release=Path(release);output=Path(output);output.mkdir(parents=True,exist_ok=True)
 manifest=json.loads((release/'release_manifest.json').read_text())
 if manifest['LU_bins']!=24 or manifest['UL_bins']!=12 or manifest['polarized_constraints']!='none; LL and LLcosphi float as nuisance amplitudes':raise ValueError('Unexpected release definition')
 rows=[]
 with Path(pip_source).open() as f:pip=list(csv.DictReader(f))
 # The canonical pi+ source table is unchanged by a pi-minus refresh.
 from prepare_measurements import OBS,TP
 for r in pip:
  b=int(r['bin_number']);xi=(b-1)//6;ti=(b-1)%6
  if str(r['nominal_fit_valid']).lower() not in ['true','1']:raise ValueError('Invalid pi+ fit')
  for obs in OBS:rows.append(dict(channel='pip',observable=obs,x_index=xi,bin=b,x=float(r['mean_minus_tprime_gev2']),x_low=TP[ti],x_high=TP[ti+1],value=float(r[obs]),stat=float(r[obs+'_stat']),syst=float(r[obs+'_point_to_point_systematic']),status='measured'))
 with (release/'measurements.csv').open() as f:pim=list(csv.DictReader(f))
 for obs,n in [('lu1',24),('ul1',12),('ul2',12)]:
  subset=[r for r in pim if r['observable']==obs]
  if sorted(int(r['bin']) for r in subset)!=list(range(1,n+1)):raise ValueError('Incomplete pi-minus '+obs)
 if len(pim)!=48 or any(r['channel']!='pim' or r['status']!='measured' for r in pim):raise ValueError('Unexpected pi-minus rows')
 rows.extend(pim)
 fields=list(rows[0])
 with (output/'measurements.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)
 if Path(pip_source).resolve()!=(output/'pip_source.csv').resolve():shutil.copy2(pip_source,output/'pip_source.csv')
 provenance=dict(pip='24 bins, original pi+ final source including point-to-point systematic errors',pim='Preliminary release: 24 LU bins and 12 bins each for UL sinphi and UL sin2phi. LL amplitudes float but are not reported. UU fixed zero. Statistical uncertainties including dilution only.',cutset_id_pim=manifest['identity']['cutset_id'],release_signature=manifest['signature'],sources={str(Path(pip_source).resolve()):hashlib.sha256(Path(pip_source).read_bytes()).hexdigest(),str((release/'measurements.csv').resolve()):hashlib.sha256((release/'measurements.csv').read_bytes()).hexdigest()},rows=len(rows))
 (output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 from theory_io import load_measurements
 load_measurements(output/'measurements.csv')
 return rows
