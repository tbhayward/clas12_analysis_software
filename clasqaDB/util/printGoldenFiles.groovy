// print list of golden files stored in QADB
// - "golden" means no defects

import clasqa.QADB

QADB qa = new QADB()

// run loop (sorted)
qa.getQaTree().sort{ a,b -> a.key.toInteger() <=> b.key.toInteger() }.each{
  runnum,runQA ->

  // file loop (sorted)
  runQA.sort{ c,d -> c.key.toInteger() <=> d.key.toInteger() }.each{
    filenum, fileQA ->
    if(fileQA.defect==0) println "$runnum $filenum"
  }
} // end run loop
