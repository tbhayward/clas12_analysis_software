import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
spec=importlib.util.spec_from_file_location('audit_campaign',BASE/'audit_campaign.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
class AuditTest(unittest.TestCase):
 def test_complete_and_broken_campaign(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)
   for name in ['chunks','root','reports','hosts','claims']:(root/name).mkdir()
   def write(path,data):path.write_text(json.dumps(data))
   ds=dict(id='RGC_summer2022_C',target='C');inputfile=dict(path='/cache/test.hipo',size=100)
   job=dict(id='batch',dataset=ds['id'],target='C',run=16296,files=[inputfile])
   write(root/'manifest.json',dict(campaign_kind="loose_eppi_minus_v1",mode='test',datasets=[ds],jobs=[job]))
   write(root/'inventory.json',[dict(selected_runs=[16296],files={'16296':[inputfile]})])
   chunk=root/'chunks/batch.root';chunk.write_bytes(b'example batch')
   write(root/'chunks/batch.json',dict(entries=12,root_sha256=audit.digest(chunk),diagnostics=dict(counts=dict(rows=12))))
   final=root/'root'/('eppi_minus_loose_'+ds['id']+'.root');final.write_bytes(b'example final')
   write(root/'reports'/(ds['id']+'.json'),dict(entries=12,batches=1,file=str(final),sha256=audit.digest(final)))
   with contextlib.redirect_stdout(io.StringIO()):result=audit.inspect(root,True)
   self.assertEqual(result['errors'],[]);self.assertEqual(result['total_entries'],12)
   self.assertEqual(len(result['root_checks']),2)
   final.write_bytes(b'changed final')
   with contextlib.redirect_stdout(io.StringIO()):result=audit.inspect(root,True)
   self.assertTrue(any('checksum mismatch' in e for e in result['errors']))
   final.unlink();(root/'claims/failed_batch').mkdir()
   with contextlib.redirect_stdout(io.StringIO()):result=audit.inspect(root)
   self.assertTrue(any('Missing final ROOT' in e for e in result['errors']))
   self.assertTrue(any('Unresolved claims' in e for e in result['errors']))
if __name__=='__main__':unittest.main()
