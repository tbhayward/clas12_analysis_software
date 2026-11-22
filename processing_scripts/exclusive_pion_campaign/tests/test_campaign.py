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
   (root/'build/source').mkdir(parents=True)
   (root/'build/source/processing_exclusive_piminus.groovy').write_text('if(args.length!=3) throw new Exception()')
   calls=[]
   def fake_run(cmd,**kwargs):
    calls.append(cmd[0])
    if str(cmd[0])=='groovy':
     self.assertEqual(len(cmd),7)
     self.assertIn('-Xmx2048m',kwargs['env']['JAVA_OPTS'])
     self.assertNotIn('ActiveProcessorCount',kwargs['env']['JAVA_OPTS'])
     Path(cmd[5]).write_text('events')
    else:
     Path(cmd[2]).write_text('root')
     return type('Result',(),dict(stdout='7\n',stderr='',returncode=0))()
   with patch.object(c,'run',side_effect=fake_run),patch.object(c.subprocess,'run',side_effect=fake_run):
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
 def test_recover_without_reprocessing(self):
  from types import SimpleNamespace
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)
   for name in ['hosts','claims/job','work/job','logs','chunks','build']:(root/name).mkdir(parents=True)
   c.atomic_json(root/'manifest.json',dict(jobs=[dict(id='job',target='D2',run=6265,bytes=5)]))
   c.atomic_json(root/'claims/job/failed.json',dict(error='old conversion failed'))
   (root/'logs/job.log').write_text('SUMMARY seen=4 rows=2\nold converter error\n')
   (root/'work/job/events.txt').write_text('saved rows')
   def command(cmd,**kwargs):
    if str(cmd[0])=='g++':Path(cmd[cmd.index('-o')+1]).write_text('binary')
   def convert(cmd,**kwargs):
    self.assertNotEqual(str(cmd[0]),'groovy')
    Path(cmd[2]).write_text('recovered root')
    return SimpleNamespace(returncode=0,stdout='2\n',stderr='')
   with patch.object(c,'verified',return_value=root/'build'),patch.object(c.subprocess,'check_output',return_value=''),patch.object(c,'run',side_effect=command),patch.object(c.subprocess,'run',side_effect=convert),patch.object(c,'merge') as merger:
    c.recover_conversion(SimpleNamespace(campaign=root))
    merger.assert_called_once()
   self.assertEqual(c.read(root/'chunks/job.json')['entries'],2)
   self.assertTrue((root/'work/job/events.txt').exists())
   self.assertFalse((root/'claims/job').exists())
 def test_progress_accounts_for_both_hosts(self):
  import io
  import contextlib
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)
   for name in ['chunks','claims/b','claims/c','work/b']:(root/name).mkdir(parents=True)
   jobs=[dict(id=n,bytes=1024) for n in 'abcd']
   c.atomic_json(root/'chunks/a.json',dict(entries=42))
   c.atomic_json(root/'claims/b/owner.json',dict(host='ifarm2402',started=0))
   c.atomic_json(root/'claims/c/failed.json',dict(error='failure'))
   c.atomic_json(root/'work/b/progress.json',dict(stage='processing',updated=c.time.time(),counts=dict(seen=200,rows=10,qa_reject=100,zero_helicity=20)))
   with contextlib.redirect_stdout(io.StringIO()) as output:state=c.progress_snapshot(root,dict(jobs=jobs))
   self.assertEqual(state,dict(done=1,active=1,queued=1,failed=1,entries=42))
   self.assertIn('ifarm2402',output.getvalue());self.assertIn('QA rejected 100',output.getvalue())
 def test_retest_selection_and_gate(self):
  import contextlib,io
  groups={n:[dict(size=n)] for n in range(1,11)}
  selected=c.retest_runs(groups,list(groups),[1])
  self.assertIn(1,selected);self.assertEqual(len(selected),3);self.assertNotIn(2,selected)
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);(root/'chunks').mkdir()
   c.atomic_json(root/'manifest.json',dict(retest_original_runs={'example':[1]},jobs=[dict(id='old',dataset='example',run=1),dict(id='new',dataset='example',run=5)]))
   c.atomic_json(root/'chunks/old.json',dict(entries=0,diagnostics=dict(counts=dict(qa_reject=100))))
   c.atomic_json(root/'chunks/new.json',dict(entries=0))
   with contextlib.redirect_stdout(io.StringIO()),self.assertRaises(RuntimeError):c.retest_report(root)
   c.atomic_json(root/'chunks/new.json',dict(entries=10))
   with contextlib.redirect_stdout(io.StringIO()) as output:c.retest_report(root)
   self.assertIn('PASS',output.getvalue())
if __name__=='__main__':unittest.main()
