// dump QADB data for a specified run; this is for QADB validation and testing,
// but may still be useful as an example; see dumpQADB.cpp for a more
// user-friendly example
// - this program does not require a HIPO file or clas12root; it only loops
//   through the QADB itself, performing lookup by file number
// - you can specify a run number as an argument

#include <iostream>
#include <bitset>

// QADB header and namespace
#include "QADB.h"
using namespace QA;

const int maxbit = 5; // max defect bit number

using namespace std;

// print separator
void sep(string s,int n) {for(int k=0;k<n;k++) cout<<s; cout<<endl;};

// print error
void err(string s) { cerr << "ERROR: " << s << endl; };

// MAIN
int main(int argc, char ** argv) {

  // specify run number
  int runnum = 5160; // default
  if(argc>1) runnum = (int) strtof(argv[1],NULL);
  cout << "test QADB for RUN NUMBER " << runnum << endl;


  // instantiate QADB
  QADB * qa = new QADB();
  // alternatively, specify run range to restrict QADB (may be more efficient)
  //QADB * qa = new QADB(5000,5500);


  
  // loop through files
  int evnum;
  string defname;
  int chargeInt;
  for(int filenum=0; filenum<=qa->GetMaxFilenum(runnum); filenum+=5) {
    sep("=",50);

    //err("test error print");

    // query by file number
    qa->QueryByFilenum(runnum,filenum);
    evnum = qa->GetEvnumMin(); // evnum needed for QA cut methods

    // check run and file number accessors: make sure that they
    // are equal to what we asked for
    cout << "- run,file,evnum" << endl;
    cout << qa->GetRunnum() << " " << runnum << endl;
    if(qa->GetRunnum() != runnum) err("QADB::GetRunnum != runnum");
    cout << qa->GetFilenum() << " " << filenum << endl;
    if(qa->GetFilenum() != filenum) err("QADB::GetFilenum != filenum");

    // check event number: report an error if evnum min>=max
    cout << qa->GetEvnumMin() << " " << qa->GetEvnumMax() << endl;
    if(qa->GetEvnumMin() >= qa->GetEvnumMax())
      err("GetEvnumMin() >= GetEvnumMax()");

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
  cout << "- charge, max filenum" << endl;
  cout << ((int)(1000*qa->GetAccumulatedCharge())) << endl;

  // print max file number
  cout << qa->GetMaxFilenum(runnum) << endl;

  sep("=",50);
  

  return 0;
}

