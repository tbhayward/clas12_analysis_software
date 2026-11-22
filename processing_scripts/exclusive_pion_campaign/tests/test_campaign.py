import concurrent.futures
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('campaign',Path(__file__).resolve().parents[1]/'campaign.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
class CampaignTest(unittest.TestCase):
 def test_inventory(self):
  ds=c.datasets();self.assertEqual(len(ds),19);self.assertEqual(len({d['id'] for d in ds}),19)
  self.assertEqual(sum(d['target']=='CD2' for d in ds),1)
  for f,n in [('sidisdvcs_016137_00000.hipo',16137),('skim_006420.hipo',6420),('skim_011323_00001.hipo',11323)]:self.assertEqual(c.run_number(Path(f)),n)
  with self.assertRaises(ValueError):c.run_number(Path('mystery.hipo'))
 def test_batches(self):
  files=[dict(path=str(n),size=3) for n in range(10)];b=c.make_batches(files,10)
  self.assertEqual([len(x) for x in b],[3,3,3,1]);self.assertEqual(sum(map(len,b)),10)
 def test_atomic_claim_and_failure(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)
   for name in ['claims','chunks','work','logs']:(root/name).mkdir()
   inp=root/'input.hipo';inp.write_text('input');s=inp.stat()
   job=dict(id='one',dataset='RGB_spring2019_D2',target='D2',run=6420,bytes=5,files=[dict(path=str(inp),size=s.st_size,mtime_ns=s.st_mtime_ns)])
   calls=[]
   def fake_run(cmd,**kwargs):
    calls.append(cmd[0])
    if str(cmd[0])=='groovy':Path(cmd[-2]).write_text('events')
    else:
     Path(cmd[2]).write_text('root')
     return type('Result',(),dict(stdout='7\n'))()
   with patch.object(c,'run',side_effect=fake_run):
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as pool:list(pool.map(lambda _:c.process_job(root,root/'build',job),range(20)))
   self.assertEqual(calls.count('groovy'),1)
   self.assertEqual(c.read(root/'chunks/one.json')['entries'],7)
   self.assertEqual(list((root/'claims').iterdir()),[])
   (root/'chunks/one.json').unlink()
   with patch.object(c,'run',side_effect=RuntimeError('simulated failure')):c.process_job(root,root/'build',job)
   self.assertTrue((root/'claims/one/failed.json').exists());self.assertFalse((root/'chunks/one.json').exists())
   with patch.object(c,'run',side_effect=AssertionError('must not retry')):c.process_job(root,root/'build',job)
 def test_worker_finishes_and_merges_locally(self):
  from types import SimpleNamespace
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)
   for name in ['hosts','claims','chunks']:(root/name).mkdir()
   c.atomic_json(root/'manifest.json',dict(jobs=[]))
   args=SimpleNamespace(campaign=root,allow_other_host=False,workers=10)
   with patch.object(c.socket,'gethostname',return_value='ifarm2401'),patch.object(c,'verified',return_value=root/'build'),patch.object(c,'status',return_value=True),patch.object(c,'merge') as merger,patch.object(c.subprocess,'Popen',side_effect=AssertionError('No remote launch permitted')):
    c.worker(args)
    merger.assert_called_once_with(args)
   self.assertEqual(list((root/'hosts').iterdir()),[])
   self.assertEqual(list((root/'claims').iterdir()),[])
if __name__=='__main__':unittest.main()
