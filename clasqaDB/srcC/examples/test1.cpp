// dump QADB data for a specified file; this is for QADB validation and testing,
// but may still be useful as an example; see dumpQADB.cpp for a more
// user-friendly example
// - this program does requires a HIPO file and clas12root

#include <iostream>
#include <bitset>

// QADB header and namespace
#include "QADB.h"
using namespace QA;

// clas12root headers
#include "reader.h"
#include "clas12reader.h"

const int maxbit = 5; // max defect bit number

using namespace std;
using namespace clas12; // for clas12root

// print separator
void sep(string s,int n) {for(int k=0;k<n;k++) cout<<s; cout<<endl;};

// print error
void err(string s) { cerr << "ERROR: " << s << endl; };

// MAIN
int main(int argc, char ** argv) {

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


  
  // event loop
  int evnum,runnum,filenum;
  int filenumTmp = -1000;
  string defname;
  int chargeInt;
  while(c12->next()==true) {

    // get run number and event number
    runnum = c12->runconfig()->getRun();
    evnum = c12->runconfig()->getEvent();

    // query
    qa->Query(runnum,evnum);

    // skip tag1 events
    if(runnum==0) continue;

    // apply QA cut
    //if(!(qa->OkForAsymmetry(runnum,evnum))) continue;

    // get filenum; only print info if we queried a different file
    filenum = qa->GetFilenum();
    if(filenum==filenumTmp) continue;
    filenumTmp = filenum;
  
    sep("=",50);
    //err("test error print");

    // check run and file number accessors: make sure that they
    // are equal to what we asked for
    cout << "- run,evnum,evnumRange,filenum" << endl;
    cout << qa->GetRunnum() << " " << runnum << endl;
    if(qa->GetRunnum() != runnum) err("QADB::GetRunnum != runnum");

    // check event number: report an error if evnum min>=max
    cout << evnum << endl;
    cout << qa->GetEvnumMin() << " " << qa->GetEvnumMax() << endl;
    if(qa->GetEvnumMin() >= qa->GetEvnumMax())
      err("GetEvnumMin() >= GetEvnumMax()");
    cout << filenum << endl;

    // print charge (convert to pC and truncate, for easier comparison)
    chargeInt = (int) (1000*qa->GetCharge());
    cout << "- charge,comment" << endl;
    cout << chargeInt << endl;
    qa->AccumulateCharge();

    // print comment
    cout << "\"" << qa->GetComment() << "\"" << endl;


    // print overall defect info
    cout << "- defect" << endl;
    cout << qa->GetDefect() << endl;
    for(int sec=1; sec<=6; sec++) 
      cout << " " << qa->GetDefect(sec); cout << endl;

    // print defect bit info
    for(int bit=0; bit<=maxbit; bit++) {
      // translate bit number to name; check if QADB::Bit returns correct bit
      switch(bit) {
        case 0: defname="TotalOutlier"; break;
        case 1: defname="TerminalOutlier"; break;
        case 2: defname="MarginalOutlier"; break;
        case 3: defname="SectorLoss"; break;
        case 4: defname="LowLiveTime"; break;
        case 5: defname="Misc"; break;
      };
      if(qa->Bit(defname.c_str()) != bit) err("QADB::Bit problem");
      cout << qa->Bit(defname.c_str()) << " " << bit << " " << defname << endl;
      // print defect info
      cout << qa->HasDefect(defname.c_str()) << endl;
      for(int sec=1; sec<=6; sec++) 
        cout << " " << qa->HasDefect(defname.c_str(),sec); cout << endl;
      // print bit masking
      qa->SetMaskBit(defname.c_str());
      cout << qa->GetMask() << " " << qa->Pass(runnum,evnum) << endl;
      qa->SetMaskBit(defname.c_str(),false);
    };

    // print QA cuts (see above for custom cut check with mask)
    cout << "- cuts" << endl;
    cout << qa->Golden(runnum,evnum) << endl;
    cout << qa->OkForAsymmetry(runnum,evnum) << endl;
  };

  sep("=",50);

  // print accumulated charge
  cout << "- charge" << endl;
  cout << ((int)(1000*qa->GetAccumulatedCharge())) << endl;

  sep("=",50);
  
  return 0;
}

