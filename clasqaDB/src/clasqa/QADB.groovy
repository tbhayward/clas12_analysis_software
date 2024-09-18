package clasqa

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import clasqa.Tools

class QADB {

  //..............
  // constructor
  //``````````````
  // arguments:
  // - runnumMin and runnumMax: if both are negative (default), then the
  //   entire QADB will be read; you can restrict to a specific range of
  //   runs to limit QADB, which may be more effecient
  // - verbose: if true, print (a lot) more information
  public QADB(int runnumMin=-1, runnumMax=-1, boolean verbose_=false) {

    // setup
    verbose = verbose_
    util = new Tools()
    nbits = util.bitDefinitions.size()
    dbDirN = System.getenv('QADB') + '/qadb'
    if(dbDirN==null) {
      System.err << "ERROR: env var QADB not set; source env.sh\n\n\n"
      return
    }
    if(verbose) println("QADB dir = ${dbDirN}/qadb")

    // concatenate trees
    qaTree = [:]
    chargeTree = [:]
    slurper = new JsonSlurper()
    def dbDir = new File(dbDirN)
    def dbFilter = ~/.*Tree.json$/
    def slurpAction = { tree,branch ->
      def rnum = branch.key.toInteger()
      if( ( runnumMin<0 && runnumMax<0 ) || 
          ( rnum>=runnumMin && rnum<=runnumMax)) {
        if(verbose) println "run $rnum"
        tree << branch
      }
    }
    dbDir.traverse( 
      type:groovy.io.FileType.FILES,
      maxDepth:1,
      nameFilter:dbFilter)
    { dbFile ->
      if(verbose) println "read " + dbFile.getAbsoluteFile()
      if(dbFile.name.contains("qaTree.json"))
        slurper.parse(dbFile).each{ obj -> slurpAction(qaTree,obj) }
      else if(dbFile.name.contains("chargeTree.json"))
        slurper.parse(dbFile).each{ obj -> slurpAction(chargeTree,obj) }
    }

    // defect mask used for asymmetry analysis
    asymMask = 0
    asymMask += 0x1 << Bit('TotalOutlier')
    asymMask += 0x1 << Bit('TerminalOutlier')
    asymMask += 0x1 << Bit('MarginalOutlier')
    asymMask += 0x1 << Bit('SectorLoss')

    // initialize local vars
    runnum = -1
    filenum = -1
    evnumMin = -1
    evnumMax = -1
    found = false
    mask = 0
    charge = 0
    chargeTotal = 0
    chargeCounted = false
    chargeCountedFiles = []
  }


  //...............................
  // golden QA cut
  //```````````````````````````````
  // returns false if the event is in a file with *any* defect
  public boolean golden(int runnum_, int evnum_) {
    def foundHere = query(runnum_,evnum_)
    return foundHere && defect==0
  }


  //.....................................
  // QA for spin asymmetry analysis
  //`````````````````````````````````````
  // if true, this event is good for a spin asymmetry analysis
  public boolean OkForAsymmetry(int runnum_, int evnum_) {

    // perform lookup
    def foundHere = query(runnum_,evnum_)
    if(!foundHere) return false;

    // check for bits which will always cause the file to be rejected 
    // (asymMask is defined in the constructor)
    if( defect & asymMask ) return false

    // special cases for `Misc` bit
    if(hasDefect('Misc')) {

      // check if this is a run on the list of runs with a large fraction of
      // events with undefined helicity; if so, accept this run, since none of
      // these files are marked with `Misc` for any other reasons
      if( runnum_ in [ 5128, 5129, 5130, 5158,
                       5159, 5160, 5163, 5165,
                       5166, 5167, 5168, 5169,
                       5180, 5181, 5182, 5183, 5567 ]) return true

      // check if this run had an FADC failure; there is no indication spin
      // asymmetries are impacted by this issue
      else if( runnum_ in 6736..6757 ) return true

      // otherwise, this file fails the QA
      else return false

    }

    // otherwise, this file passes the QA
    return true
  }

  
  //...............................
  // user-defined custom QA cuts
  //```````````````````````````````
  // first set which defect bits you want to filter out; by default
  // none are set; the variable `mask` will be applied as a mask
  // on the defect bits
  public void setMaskBit(String bitStr, boolean state=true) { 
    def defectBit = Bit(bitStr)
    if(defectBit<0 || defectBit>=nbits)
      System.err << "ERROR: QADB::setMaskBit called for unknown defectBit\n"
    else  {
      mask &= ~(0x1 << defectBit)
      if(state) mask |= (0x1 << defectBit)
    }
  }
  // access the custom mask, if you want to double-check it
  public int getMask() { return mask }
  // then call this method to check your custom QA cut for a given
  // run number and event number
  public boolean pass(int runnum_, int evnum_) {
    def foundHere = query(runnum_,evnum_)
    return foundHere && !(defect & mask)
  }


  //..............
  // accessors
  //``````````````
  // --- access this file's info
  public int getRunnum() { return found ? runnum.toInteger() : -1 }
  public int getFilenum() { return found ? filenum.toInteger() : -1 }
  public String getComment() { return found ? comment : "" }
  public int getEvnumMin() { return found ? evnumMin : -1 }
  public int getEvnumMax() { return found ? evnumMax : -1 }
  public double getCharge() { return found ? charge : -1 }
  // --- access QA info
  // check if the file has a particular defect
  // - if sector==0, checks the OR of all the sectors
  // - if an error is thrown, return true so file will be flagged
  public boolean hasDefect(String name_, int sector=0) {
    return hasDefectBit(Bit(name_),sector)
  }
  // - alternatively, check for defect by bit number
  public boolean hasDefectBit(int defect_, int sector=0) {
    return (getDefect(sector) >> defect_) & 0x1
  }
  // get defect bitmask; if sector==0, gets OR of all sectors' bitmasks
  public int getDefect(int sector=0) {
    if(!found) return -1
    if(sector==0) return defect
    else if(sector>=1 && sector<=6) return sectorDefect["$sector"]
    else {
      System.err << "ERROR: bad sector number in QADB::getDefect\n"
      return -1;
    }
  }
  // translate defect name to defect bit
  public int Bit(String name_) { return util.bit(name_) }
  // --- access the full tree
  public def getQaTree() { return qaTree }
  public def getChargeTree() { return chargeTree }


  //....................................................................
  // query qaTree to get QA info for this run number and event number
  // - a lookup is only performed if necessary: if the run number changes
  //   or if the event number goes outside of the range of the file which
  //   most recently queried
  // - this method is called automatically when evaluating QA cuts
  //````````````````````````````````````````````````````````````````````
  public boolean query(int runnum_, int evnum_) {

    // perform lookup, only if needed
    if( runnum_ != runnum ||
        ( runnum_ == runnum && (evnum_ < evnumMin || evnum_ > evnumMax ))
    ) {
      // reset vars
      runnum = runnum_
      filenum = -1
      evnumMin = -1
      evnumMax = -1
      charge = -1
      found = false

      // search for file which contains this event
      if(verbose) println "query qaTree..."
      qaFile = qaTree["$runnum"].find{
        if(verbose) println " search file "+it.key+" for event $evnum_"
        evnum_ >= it.value['evnumMin'] && evnum_ <= it.value['evnumMax']
      }

      // if file found, set variables
      if(qaFile!=null) {
        queryByFilenum(runnum_,qaFile.key.toInteger())
      }

      // print a warning if a file was not found for this event
      // - this warning is suppressed for 'tag1' events
      if(!found && runnum_!=0) {
        System.err << "WARNING: QADB::query could not find " <<
        "runnum=$runnum_ evnum=$evnum_\n"
      }
    }

    // result of query
    return found
  }

  //........................................
  // if you know the DST file number, you can call QueryByFilenum to perform
  // lookups via the file number, rather than via the event number
  // - you can subsequently call any QA cut method, such as `Golden()`;
  //   although QA cut methods require an event number, no additional lookup
  //   query will be performed since it already has been done in QueryByFilenum
  //````````````````````````````````````````
  public boolean queryByFilenum(int runnum_, int filenum_) {

    // if the run number or file number changed, perform new lookup
    if( runnum_ != runnum || filenum_ != filenum) {
      
      // reset vars
      runnum = runnum_
      filenum = filenum_
      evnumMin = -1
      evnumMax = -1
      charge = -1
      found = false

      if(qaTree["$runnum"]!=null) {
        if(qaTree["$runnum"]["$filenum"]!=null) {
          qaFileTree = qaTree["$runnum"]["$filenum"]
          evnumMin = qaFileTree['evnumMin']
          evnumMax = qaFileTree['evnumMax']
          comment = qaFileTree['comment']
          defect = qaFileTree['defect']
          sectorDefect = [:]
          qaFileTree['sectorDefects'].each{ sec,defs ->
            sectorDefect[sec] = 0
            defs.each{ 
              sectorDefect[sec] += 0x1 << it
            }
          }
          chargeMin = chargeTree["$runnum"]["$filenum"]["fcChargeMin"]
          chargeMax = chargeTree["$runnum"]["$filenum"]["fcChargeMax"]
          charge = chargeMax - chargeMin
          chargeCounted = false
          found = true
        }
      }

      // print a warning if a file was not found for this event
      // - this warning is suppressed for 'tag1' events
      if(!found && runnum_!=0) {
        System.err << "WARNING: QADB::queryByFilenum could not find " <<
          "runnum=$runnum_ filenum=$filenum_\n"
      }
    }

    // result of query
    return found
  }


  // get maximum file number for a given run (useful for QADB validation)
  public int getMaxFilenum(int runnum_) {
    int maxFilenum=0
    qaTree["$runnum_"].each{ 
      maxFilenum = it.key.toInteger() > maxFilenum ? 
                   it.key.toInteger() : maxFilenum
    }
    return maxFilenum
  }




  //.................................
  // Faraday Cup charge accumulator
  //`````````````````````````````````
  // -- accumulator
  // call this method after evaluating QA cuts (or at least after calling query())
  // to add the current file's charge to the total charge;
  // - charge is accumulated per DST file, since the QA filters per DST file
  // - a DST file's charge is only accounted for if we have not counted it before
  public void accumulateCharge() {
    if(!chargeCounted) {
      if(!( [runnum,filenum] in chargeCountedFiles )) {
        chargeTotal += charge
        chargeCountedFiles << [runnum,filenum]
      }
      chargeCounted = true
    }
  }
  // -- accessor
  // call this method at the end of your event loop
  public double getAccumulatedCharge() { return chargeTotal }
  // reset accumulated charge, if you ever need to
  public void resetAccumulatedCharge() { chargeTotal = 0 }



  //===============================================================

  public Tools util

  private def qaTree
  private def chargeTree
  private def qaFile
  private def qaFileTree

  private boolean verbose
  private def dbDirN
  private def slurper

  private def runnum,filenum,evnumMin,evnumMax,evnumMinTmp,evnumMaxTmp
  private double charge,chargeMin,chargeMax,chargeTotal
  private boolean chargeCounted
  private def chargeCountedFiles
  private int defect
  private def sectorDefect
  private String comment
  private boolean found
  private int nbits
  private def mask
  private def asymMask
}
