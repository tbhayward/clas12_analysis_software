// print list of files, and whether they are Golden, OkForAsymmetry, etc.

#include "QADB.h"
#include <fstream>
using namespace std;

int main(int argc, char ** argv) {

  QA::QADB * qa = new QA::QADB();
  
  /// get sorted list of runs
  auto qaTree = qa->GetQaTree();
  vector<int> runList;
  for(auto it=qaTree->MemberBegin(); it!=qaTree->MemberEnd(); ++it)
    runList.push_back( atoi((it->name).GetString()) );
  sort(runList.begin(),runList.end());

  /// start output file
  string outFileN = getenv("QADB") ? getenv("QADB") : "";
  if(outFileN.compare("")==0) {
    cerr << "ERROR: QADB environment variable not set" << endl;
    return 1;
  };
  outFileN += "/text/summary.txt";
  ofstream outFile(outFileN);
  if(!outFile) {
    cerr << "ERROR: cannot create output file " << outFileN << endl;
    return 1;
  }
  outFile << "#Run File Golden OkForAsymmetry" << endl;

  /// loop over runs ----------------
  vector<int> fileList;
  for(auto runnum : runList) {

    /// get sorted list of files
    fileList.clear();
    auto runTree = (*qaTree)[to_string(runnum).c_str()].GetObject();
    for(auto it=runTree.MemberBegin(); it!=runTree.MemberEnd(); ++it)
      fileList.push_back( atoi((it->name).GetString()) );
    sort(fileList.begin(),fileList.end());

    /// loop over files -----------------
    for(auto filenum : fileList) {

      /// get some event number in the middle of the file, such as the average
      auto fileTree = runTree[to_string(filenum).c_str()].GetObject();
      int evnumMin = fileTree["evnumMin"].GetInt();
      int evnumMax = fileTree["evnumMax"].GetInt();
      int evnumAve = (int) ( (evnumMin+evnumMax)/2 );

      /// query by run number and event number
      outFile << runnum
        << " " << filenum
        << " " << qa->Golden(runnum,evnumAve)
        << " " << qa->OkForAsymmetry(runnum,evnumAve)
        << endl;

      /// use this output to check consistency with listOfGoldenFiles.txt
      // if(qa->Golden(runnum,evnumAve)) cout << runnum << " " << filenum << endl;

    } /// end file loop
  } /// end run loop

  outFile.close();
  return 0;
}

