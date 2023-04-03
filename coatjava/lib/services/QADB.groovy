package clasqa

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import clasqa.Tools

class QADB {

  //..............
  // constructor
  //``````````````
  public QADB(boolean verbose_=false) {

    // setup
    verbose = verbose_
    util = new Tools()
    nbits = util.bitDefinitions.size()
    dbDirN = System.getenv('QADB')
    if(dbDirN==null) {
      System.err << "ERROR: env var QADB not set; source env.sh\n\n\n"
      return
    }
    if(verbose) println("QADB dir = $dbDirN")

    // concatenate qaTree.json files into one tree
    qaTree = [:]
    slurper = new JsonSlurper()
    def dbDir = new File(dbDirN)
    def dbFilter = ~/^qaTree.json$/
    dbDir.traverse( 
      type:groovy.io.FileType.FILES,
      maxDepth:1,
      nameFilter:dbFilter) { dbFile ->
      if(verbose) println "read " + dbFile.getAbsoluteFile()
      slurper.parse(dbFile).each{ obj -> qaTree << obj }
    }
    //if(verbose) println util.pPrint(qaTree) // pretty print qa tree

    // defect mask used for asymmetry analysis
    asymMask = 0
    asymMask += 0x1 << util.bit('TotalOutlier')
    asymMask += 0x1 << util.bit('TerminalOutlier')
    asymMask += 0x1 << util.bit('MarginalOutlier')
    asymMask += 0x1 << util.bit('SectorLoss')

    // initialize local vars
    runnum = -1
    filenum = -1
    evnumMin = -1
    evnumMax = -1
    found = false
    mask = 0

  }

  
  //...............................
  // methods for simple QA check
  //```````````````````````````````
  // first set which defect bits you want to filter out; by default
  // none are called
  public void setMaskBit(String bitStr, boolean state=true) { 
    def bit = util.bit(bitStr)
    if(bit<0 || bit>nbits)
      System.err << "ERROR: QADB::maskBit called for unknown bit\n"
    else  {
      mask &= ~(0x1 << bit)
      if(state) mask |= (0x1 << bit)
    }
  }
  public int getMask() { return mask }

  // then call this method for a run number and event number;
  // you do not need to call query;
  // if query fails, this method returns false
  public boolean pass(int runnum_, int evnum_) {
    def foundHere = query(runnum_,evnum_)
    return foundHere && !(defect & mask)
  }

  // alternatively, call this method to simply check if the file has a
  // defect
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
    if(hasDefectName('Misc')) {

      // if only the `Misc` defect bit is set, check if this is a run on the
      // list of runs with a large fraction of events with undefined helicity;
      // if so, accept this run, since none of these files are marked with
      // `Misc` for any other reasons
      if( defect == (0x1 << util.bit('Misc')) ) {
        if( runnum_ in [ 5128, 5129, 5130, 5158,
                         5159, 5160, 5163, 5165,
                         5166, 5167, 5168, 5169,
                         5180, 5181, 5182, 5183 ]) return true
        else return false
      }

      // otherwise, reject this file
      return false
    }

    return true
  }

    
      



  //..............
  // accessors
  //``````````````

  // access this file's info
  public int getFilenum() { return found ? filenum.toInteger() : -1 }
  public String getComment() { return found ? comment : "" }
  public int getEvnumMin() { return found ? evnumMin : -1 }
  public int getEvnumMax() { return found ? evnumMax : -1 }

  // check if the file has a particular defect; if sector==0, checks
  // the OR of all the sectors
  // - if an error is thrown, return true so file will be flagged
  public boolean hasDefect(int defect_, int sector=0) {
    if( ! (0..nbits).contains(defect_) ) {
      System.err << "ERROR: bad defect number for QADB::hasDefect\n"
      return true
    }
    if(!found) return true
    if(sector>0) return ( sectorDefect["$sector"] >> defect_ ) & 0x1
    else return ( defect >> defect_ ) & 0x1
  }
  // - alternatively, specify name of defect
  public boolean hasDefectName(String name_, int sector=0) {
    return hasDefect(util.bit(name_),sector)
  }


  // get defect bitmask; if sector==0, gets OR of all sectors' bitmasks
  public int getDefect(int sector=0) {
    if(!found) return -1
    if(sector>0) return sectorDefect["$sector"]
    else return defect
  }

  // access the full tree
  public def db() { return qaTree }




  //..............
  // query qaTree to get QA info for this run number and event number
  // - a lookup is only performed if necessary: if the run number changes
  //   or if the event number goes outside of the range of the file which
  //   most recently queried
  //``````````````
  public boolean query(int runnum_, int evnum_) {

    // perform lookup, only if needed
    if( runnum_ != runnum ||
        ( runnum_ == runnum && (evnum_ < evnumMin || evnum_ > evnumMax ))
    ) {
      runnum = runnum_
      filenum = -1
      evnumMin = -1
      evnumMax = -1
      found = false

      // search for file which contains this event
      qaTree["$runnum"].each{ f, qaFile ->
        evnumMinTmp = qaFile['evnumMin']
        evnumMaxTmp = qaFile['evnumMax']

        // found matching file, set all variables and stop the search
        if(evnum_>=evnumMinTmp && evnum_<=evnumMaxTmp) {
          filenum = f
          evnumMin = evnumMinTmp
          evnumMax = evnumMaxTmp
          comment = qaFile['comment']
          defect = qaFile['defect']
          sectorDefect = [:]

          qaFile['sectorDefects'].each{ sec,defs ->
            sectorDefect[sec] = 0
            defs.each{ 
              sectorDefect[sec] += 0x1 << it
            }
          }

          found = true
          return true
        }
      }

      // print a warning if a file was not found for this event
      // - this warning is suppressed for 'tag1' events
      if(!found) {
        if(runnum_!=0) {
          System.err << "WARNING: QADB::query could not find " <<
          "runnum=$runnum_ evnum=$evnum_\n"
        }
      }
    }

    // lookup was not needed, return current value of found
    return found
  }

  public Tools util

  private def qaTree

  private boolean verbose
  private def dbDirN
  private def slurper

  private def runnum,filenum,evnumMin,evnumMax,evnumMinTmp,evnumMaxTmp
  private def defect
  private def sectorDefect
  private def comment
  private boolean found
  private int nbits
  private def mask
  private def asymMask
}
