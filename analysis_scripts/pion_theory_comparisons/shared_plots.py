"""Selected BM sector plots, optional GK overlays, and nine-page PDFs."""
from pathlib import Path
import shutil,zlib
import numpy as np
from PIL import Image

def make_pdf(files,path):
 objects=[b'',b''];pages=[]
 def add(content):
  objects.append(content);return len(objects)
 for file in files:
  with Image.open(file) as im:
   im=im.convert('RGB');w,h=im.size;raw=zlib.compress(im.tobytes())
  image=add(f'<< /Type /XObject /Subtype /Image /Width {w} /Height {h} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length {len(raw)} >>\nstream\n'.encode()+raw+b'\nendstream')
  pw,ph=w*72/300,h*72/300;draw=f'q {pw} 0 0 {ph} 0 0 cm /Figure Do Q'.encode()
  content=add(f'<< /Length {len(draw)} >>\nstream\n'.encode()+draw+b'\nendstream')
  pages.append(add(f'<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {pw} {ph}] /Resources << /XObject << /Figure {image} 0 R >> >> /Contents {content} 0 R >>'.encode()))
 objects[0]=b'<< /Type /Catalog /Pages 2 0 R >>'
 objects[1]=f'<< /Type /Pages /Count {len(pages)} /Kids [ '.encode()+b' '.join(f'{p} 0 R'.encode() for p in pages)+b' ] >>'
 with Path(path).open('wb') as f:
  f.write(b'%PDF-1.4\n%\xe2\xe3\xcf\xd3\n');offsets=[0]
  for n,obj in enumerate(objects,1):offsets.append(f.tell());f.write(f'{n} 0 obj\n'.encode()+obj+b'\nendobj\n')
  xref=f.tell();f.write(f'xref\n0 {len(objects)+1}\n0000000000 65535 f \n'.encode())
  for offset in offsets[1:]:f.write(f'{offset:010d} 00000 n \n'.encode())
  f.write(f'trailer\n<< /Size {len(objects)+1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n'.encode())

def build(theory,data,theory_path,output,dpi,alpha,channels):
 from chiral_sectors import draw
 from figure_layout import OBS_NAMES
 root=Path(output)/'figures/shared_plots';root.mkdir(parents=True,exist_ok=True)
 gk_path=Path(theory_path)/'gk_pip.csv';gk=np.loadtxt(gk_path,delimiter=',',skiprows=1) if gk_path.exists() else None
 selected=[('pip',o) for o in ['u1','u2','lu1','ul1','ul2','ll0','ll1']]+[('pim','lu1')]
 entries=[]
 for variant in ['bm','bm+gk']:
  folder=root/variant;folder.mkdir(exist_ok=True)
  for species,obs in selected:
   if species not in channels:continue
   name=f"{'pi_plus' if species=='pip' else 'pi_minus'}_{OBS_NAMES[obs]}_chiral_sector_sums_2x2.png"
   if variant=='bm':
    source=Path(output)/'figures'/('pi_plus/featured_stories' if species=='pip' else 'pi_minus/cff_breakdowns/chiral_sectors')/name
    shutil.copy2(source,folder/name)
    entries.append(dict(path=str((folder/name).relative_to(output)),title=name,category='shared plots BM'))
   else:entries.append(draw(theory,data,species,obs,Path(output),dpi,alpha,gk=gk,shared_variant=variant))
  ordered=[folder/f"{'pi_plus' if species=='pip' else 'pi_minus'}_{OBS_NAMES[obs]}_chiral_sector_sums_2x2.png" for species,obs in selected]
  comparison=Path(output)/'figures/charge_comparison/pi_plus_pi_minus_LU_sin_phi_with_pi0_theory_2x2.png'
  if comparison.exists():
   destination=folder/comparison.name;shutil.copy2(comparison,destination);ordered.append(destination)
   entries.append(dict(path=str(destination.relative_to(output)),title='π+ / π− LU sinφ with π0 theory',category='shared plots '+variant))
  if len(ordered)==9 and all(p.exists() for p in ordered):
   make_pdf(ordered,folder/'pion_chiral_sector_comparisons_9_pages.pdf')
   (folder/'pion_chiral_sector_comparisons_8_pages.pdf').unlink(missing_ok=True)
 # Remove the superseded flat collection after replacements have been written.
 for p in root.iterdir():
  if p.is_file() and (p.name.endswith('_chiral_sector_sums_2x2.png') or p.name=='pion_chiral_sector_comparisons_8_pages.pdf'):p.unlink()
 return entries
