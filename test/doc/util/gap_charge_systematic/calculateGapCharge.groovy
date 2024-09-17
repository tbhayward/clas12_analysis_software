// determine the amount of charge *between* subsequent QA bins
// - for QA using DST 5-files, this may be nonzero
// - for QA using time bins, this should *always* be zero
// - since this is for estimating systematic for cross sections based on RG-A
//   Pass 1 data (Valerii, Sangbaek), we include a list of runs
//
// Run this script first, followed by `analyzeGapCharge.C` to make plots

import org.jlab.io.hipo.HipoDataSource
import clasqa.QADB

// arguments
if(args.length<1) {
  System.err.println "USAGE: run-groovy ${this.class.getSimpleName()}.groovy [RUN NUMBER]"
  System.err.println "- set [RUN NUMBER] to 0 to run this for hard-coded list of runs"
  System.exit(101)
}
def arg_runnum = args[0].toInteger()
QADB qa = new QADB()

// define output file
def datfileName   = "charge_gaps.dat"
def datfile       = new File(datfileName)
def datfileWriter = datfile.newWriter(false)
datfileWriter << '# Runs together with FC charge and gaps between bins\n'
datfileWriter << [
  'runnum/I',
  'binnum/I',
  'golden/I',
  'lastBin/I',
  'fcChargeMin/D',
  'fcChargeMax/D',
  'fcChargeGapToNextBin/D',
  'runset/I',
  'comment/C',
].join(':') << '\n'

// set run list(s)
def runListHash = []
if(arg_runnum == 0) {
  // Valerii's run lists
  runListHash = [
   [
     comment: 'Runs 45 nA',
     runs: [
       5032, 5036, 5038, 5039, 5040, 5041, 5043, 5045, 5052, 5053, 5116, 5117, 5119, 5120,
       5124, 5125, 5126, 5127, 5139, 5153, 5158, 5162, 5163, 5164, 5181, 5191, 5193, 5195,
       5196, 5197, 5198, 5199, 5200, 5201, 5202, 5203, 5204, 5205, 5206, 5208, 5211, 5212,
       5215, 5216, 5219, 5220, 5221, 5222, 5223, 5230, 5231, 5232, 5233, 5234, 5235, 5237,
       5238, 5248, 5345, 5346, 5347, 5349, 5351, 5354, 5355, 5367,
     ]
   ],
   [
     comment: 'Runs 50 - 55 nA',
     runs: [
       5342, 5343, 5344, 5356, 5357, 5358, 5359, 5360, 5361, 5362, 5366, 5368, 5369, 5372,
       5373, 5374, 5375, 5376, 5377, 5378, 5379, 5380, 5381, 5383, 5386, 5390, 5391, 5392,
       5393, 5398, 5401, 5403, 5404, 5406, 5407,
     ]
   ],
   [
     comment: 'Runs with wrong CCDB beam blocker constant 45 nA',
     runs: [
       5249, 5252, 5253, 5257, 5258, 5259, 5261, 5262, 5303, 5304, 5305, 5306, 5307, 5310,
       5311, 5315, 5317, 5318, 5319, 5320, 5323, 5324, 5335, 5339, 5340,
     ],
   ],
  ]
}
else {
  runListHash = [[comment: 'Single run', runs: [arg_runnum]]]
}

// loop over runs
runListHash.eachWithIndex{ hash, runset ->
  println "================================================================"
  println hash.comment
  println "================================================================"
  hash.runs.each{ runnum ->
    System.out.println runnum
    def qaTree     = qa.getQaTree()["$runnum"]
    def chargeTree = qa.getChargeTree()["$runnum"]
    qaTree.each{ binnum, bintree ->

      // the last bin has no subsequent bin, so mark it
      def nextBinnum = "${binnum.toInteger()+5}"
      def isLastBin  = qaTree[nextBinnum] == null

      // check GOLDEN bins
      def golden = bintree.defect == 0

      // calculate the gap charge
      // def gapCharge = nextBintree.fcChargeMin - bintree.fcChargeMax
      def thisBinChargeTree = chargeTree["$binnum"]
      def nextBinChargeTree = chargeTree["$nextBinnum"]
      def gapCharge = isLastBin ? -9999999.0 : nextBinChargeTree.fcChargeMin - thisBinChargeTree.fcChargeMax
      datfileWriter << [
        runnum,
        binnum,
        golden ? 1 : 0,
        isLastBin ? 1 : 0,
        thisBinChargeTree.fcChargeMin,
        thisBinChargeTree.fcChargeMax,
        gapCharge,
        runset,
        "\"${hash.comment.replace(' ','_')}\"",
      ].join(' ') << '\n'
    }
  }
}

// close
datfileWriter.flush()
datfileWriter.close()
System.out.println "Wrote $datfileName"
