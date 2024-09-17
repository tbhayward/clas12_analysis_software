// test a lookup to the QADB
// the arguments define the path to descend into the tree
// - to print the whole tree: run-groovy access.groovy
// - access a run or a file within a run
//   run-groovy access.groovy [runnum]
//   run-groovy access.groovy [runnum] [filenum]

import clasqa.QADB
import clasqa.Tools

Tools T = new Tools()
QADB qa = new QADB()

/*
println "=========== raw"
println T.jAccess(qa.getQaTree(),args)
println T.jAccess(qa.getChargeTree(),args)
*/

println "=========== pretty print"
println T.pPrint(T.jAccess(qa.getQaTree(),args))
println T.pPrint(T.jAccess(qa.getChargeTree(),args))
