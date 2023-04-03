// calculate total analyzed charge for an example event loop
// with QA cuts enabled
// - you must specify a hipo file as an argument

#include <iostream>

// clas12root headers
#include "reader.h"
#include "clas12reader.h"

// QADB header and namespace
#include "QADB.h"
using namespace QA;

using namespace clas12; // for clas12root
using namespace std;

int main(int argc, char** argv) {

  // instantiate clas12reader object for specified hipo file
  string infileN;
  if(argc<=1) {
    cerr << "USAGE: " << argv[0] << " [hipo file]" << endl;
    exit(0);
  };
  clas12reader * c12 = new clas12reader(string(argv[1]));


  // instantiate QADB
  QADB * qa = new QADB();
  // alternatively, specify run range to restrict QADB (may be more efficient)
  //QADB * qa = new QADB(5000,5500);


  // define variables
  int runnum,evnum;
  int evCount = 0;


  // event loop
  cout << "begin event loop..." << endl;
  while(c12->next()==true) {
    if(evCount%10000==0) cout << evCount << " events analyzed" << endl;
    
    // truncate event loop (for quick testing)
    if(evCount>1e5) { cout << "event loop truncated!" << endl; break; };

    // get run number and event number
    runnum = c12->runconfig()->getRun();
    evnum = c12->runconfig()->getEvent();

    // QA cuts
    if(qa->OkForAsymmetry(runnum,evnum)) {

      // accumulate charge; note that although the call to
      // QADB::accumulateCharge() charge happens for each
      // event within a DST file that passed the QA cuts, that
      // file's charge will only be accumulated once, so
      // overcounting is not possible 
      qa->AccumulateCharge();

      /* continue your analysis here */

    };

    evCount++;
  };

  // print charge
  cout << "\ntotal accumulated charge analyzed: " << endl;
  cout << "run=" << runnum << "  charge=" <<
    qa->GetAccumulatedCharge() << " nC" << endl;
};
