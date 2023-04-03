// dump QADB data for a specified run; this demonstrates access to the various
// data stored in QADB
// - this program does not require a HIPO file or clas12root; it only loops
//   through the QADB itself, performing lookup by file number
// - you can specify a run number as an argument

#include <iostream>
#include <bitset>

// QADB header and namespace
#include "QADB.h"
using namespace QA;

using namespace std;

// print separator
void sep(string s,int n) {for(int k=0;k<n;k++) cout<<s; cout<<endl;};

// MAIN
int main(int argc, char ** argv) {

  // specify run number
  int runnum = 5160; // default
  if(argc>1) runnum = (int) strtof(argv[1],NULL);
  cout << "dump QADB for RUN NUMBER " << runnum << endl;


  // instantiate QADB
  QADB * qa = new QADB();
  // alternatively, specify run range to restrict QADB (may be more efficient)
  //QADB * qa = new QADB(5000,5500);

  
  // loop through files
  int evnum;
  for(int filenum=0; filenum<=qa->GetMaxFilenum(runnum); filenum+=5) {
    sep("=",50);
    cout << "FILE NUMBER " << filenum << endl;

    // perform the lookup, by filenum
    if(qa->QueryByFilenum(runnum,filenum)) {

      // we need an event number within this file, to pass to QA criteria
      // checking methods, such as Golden; no additional Query will be called
      evnum = qa->GetEvnumMin();

      // print whether this file passes some QA cuts
      if(qa->Golden(runnum,evnum)) {
        cout << "- GOLDEN FILE!" << endl;
      } else {
        cout << "- not golden: file has defects" << endl;
        cout << (qa->OkForAsymmetry(runnum,evnum) ? "- OK" : "- NOT OK");
        cout << " for asymmetry analysis" << endl;
      };

      // print event number range
      sep("-",40);
      cout << "- event number range = " 
        << qa->GetEvnumMin() << " to "
        << qa->GetEvnumMax() << endl;

      // print charge (max accumulated charge minus min accumulated charge)
      printf("- charge (max-min) = %.3f nC\n",qa->GetCharge());

      // print defect bits (OR over all sectors)
      sep("-",40);
      cout << "- defect = " << qa->GetDefect() << " = 0b" <<
        bitset<16>(qa->GetDefect()) << endl;
      if(qa->HasDefect("TotalOutlier"))    cout << "   - TotalOutlier defect" << endl;
      if(qa->HasDefect("TerminalOutlier")) cout << "   - TerminalOutlier defect" << endl;
      if(qa->HasDefect("MarginalOutlier")) cout << "   - MarginalOutlier defect" << endl;
      if(qa->HasDefect("SectorLoss"))      cout << "   - SectorLoss defect" << endl;
      if(qa->HasDefect("LowLiveTime"))     cout << "   - LowLiveTime defect" << endl;
      if(qa->HasDefect("Misc"))            cout << "   - Misc defect" << endl;
      if(qa->HasDefect("TotalOutlierFT"))    cout << "   - TotalOutlierFT defect" << endl;
      if(qa->HasDefect("TerminalOutlierFT")) cout << "   - TerminalOutlierFT defect" << endl;
      if(qa->HasDefect("MarginalOutlierFT")) cout << "   - MarginalOutlierFT defect" << endl;
      if(qa->HasDefect("LossFT"))            cout << "   - LossFT defect" << endl;

      // print defect bits for each sector
      sep(".",40);
      for(int s=1; s<=6; s++) {
        cout << "  sector " << s << " defect = " <<
          qa->GetDefect(s) << " = 0b" <<
          bitset<16>(qa->GetDefect(s)) << endl;
      };

      // if there is a sector loss, print which sectors
      for(int s=1; s<=6; s++) {
        if(qa->HasDefect("SectorLoss",s)) {
          cout << "    -> sector " << s << " loss" << endl;
        };
      };

      // print comment
      sep(".",40);
      cout << "comment = \"" << qa->GetComment() << "\"" << endl;

    };

  };
  sep("=",50);
  

  return 0;
}

