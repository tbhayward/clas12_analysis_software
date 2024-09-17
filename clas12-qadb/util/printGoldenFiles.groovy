// print list of golden files stored in QADB
// - "golden" means no defects

import clasqa.QADB
QADB qa = new QADB()

// open output file writer
def outFileObj = new File(System.getenv('QADB')+"/text/listOfGoldenFiles.txt")
def outFile = outFileObj.newWriter(false)

// run loop (sorted)
qa.getQaTree().sort{ a,b -> a.key.toInteger() <=> b.key.toInteger() }.each{
  runnum,runQA ->

  // file loop (sorted)
  runQA.sort{ c,d -> c.key.toInteger() <=> d.key.toInteger() }.each{
    filenum, fileQA ->
    if(fileQA.defect==0) outFile << "$runnum $filenum\n"
  }
} // end run loop

outFile.close()
