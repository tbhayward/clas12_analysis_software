import os,re,subprocess,tempfile,unittest
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
BIN=Path(os.environ.get('PION_CONVERTER',BASE/'build/converter'))
@unittest.skipUnless(BIN.exists(),'Compile the converter first')
class ConverterTest(unittest.TestCase):
 def test_strict_conversion(self):
  source=(BASE/'src/convert_exclusive_piminus_txt_to_root.cpp').read_text();schema=re.search(r'const string schema="([^"]+)"',source)[1];names=schema.split()
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);csv=p/'run.csv';csv.write_text('# Mixed RGC polarization and legacy charge records\n5032,0,0,0,1250\n16137,10,5,5,-0.2,0.01\n6265,10,5,5,12345\n');text=p/'events.txt';root=p/'out.root'
   values={n:0 for n in names};values.update(runnum=16137,evnum=123,helicity=-1,detector1=2,detector2=1)
   row=' '.join(str(values[n]) for n in names)
   def invoke(target='ND3'):return subprocess.run([str(BIN),str(text),str(root),str(csv),target,'16137'],capture_output=True,text=True)
   text.write_text('# '+schema+'\n'+row+'\n');r=invoke();self.assertEqual(r.returncode,0,r.stderr)
   subprocess.run([str(BIN),'--validate',str(root),'1'],check=True)
   for bad in [row+' 9',' '.join(row.split()[:-1]),row.replace('123','nan'),row.replace('123','123.5')]:
    text.write_text('# '+schema+'\n'+bad+'\n');self.assertNotEqual(invoke().returncode,0)
   text.write_text('# '+schema+'\n'+row+'\n');csv.write_text('')
   self.assertNotEqual(invoke().returncode,0);self.assertEqual(invoke('ET').returncode,0)
   text.write_text('# '+schema+'\n');self.assertEqual(invoke('D2').returncode,0)
   subprocess.run([str(BIN),'--validate',str(root),'0'],check=True)
if __name__=='__main__':unittest.main()
