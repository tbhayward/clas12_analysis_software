#ifndef QADB_H_
#define QADB_H_

#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <map>
#include <vector>
#include <algorithm>

namespace QA {
  class QADB {
    public:

      //.................
      // constructor
      //`````````````````
      //arguments:
      // - runnumMin and runnumMax: if both are negative (default), then the
      //   entire QADB will be read; you can restrict to a specific range of
      //   runs to limit QADB, which may be more effecient
      // - verbose: if true, print (a lot) more information
      inline QADB(int runnumMin_=-1, int runnumMax=-1, bool verbose_=false);


      //...............................
      // golden QA cut
      //```````````````````````````````
      // returns false if the event is in a file with *any* defect
      inline bool Golden(int runnum_, int evnum_) { 
        bool foundHere = this->Query(runnum_,evnum_);
        return foundHere && defect==0;
      };


      //.....................................
      // QA for spin asymmetry analysis
      //`````````````````````````````````````
      // if true, this event is good for a spin asymmetry analysis
      inline bool OkForAsymmetry(int runnum_, int evnum_);


      //...............................
      // user-defined custom QA cuts
      //```````````````````````````````
      // first set which defect bits you want to filter out; by default
      // none are set; the variable `mask` will be applied as a mask
      // on the defect bits
      inline void SetMaskBit(const char * defectName, bool state=true);
      // access the custom mask, if you want to double-check it
      inline int GetMask() { return mask; };
      // then call this method to check your custom QA cut for a given
      // run number and event number
      inline bool Pass(int runnum_, int evnum_);


      //.................
      // accessors
      //`````````````````
      // --- access this file's info
      inline int GetRunnum() { return found ? runnum : -1; };
      inline int GetFilenum() { return found ? filenum : -1; };
      inline std::string GetComment() { return found ? comment : ""; };
      inline int GetEvnumMin() { return found ? evnumMin : -1; };
      inline int GetEvnumMax() { return found ? evnumMax : -1; };
      inline double GetCharge() { return found ? charge : -1; };
      // --- access QA info
      // check if the file has a particular defect
      // - if sector==0, checks the OR of all the sectors
      // - if an error is thrown, return true so file will be flagged
      inline bool HasDefect(const char * defectName, int sector=0) {
        return this->HasDefectBit(Bit(defectName),sector);
      };
      // - alternatively, check for defect by bit number
      inline bool HasDefectBit(int defect_, int sector=0) {
        return (this->GetDefect(sector) >> defect_) & 0x1;
      };
      // get this file's defect bitmask;
      // - if sector==0, gets OR of all sectors' bitmasks
      inline int GetDefect(int sector=0);
      // translate defect name to defect bit
      inline int Bit(const char * defectName);


      //.................................................................
      // query qaTree to get QA info for this run number and event number
      // - a lookup is only performed if necessary: if the run number changes
      //   or if the event number goes outside of the range of the file which
      //   most recently queried
      // - this method is called automatically when evaluating QA cuts
      //`````````````````````````````````````````````````````````````````
      inline bool Query(int runnum_, int evnum_);
      // if you know the DST file number, you can call QueryByFilenum to perform
      // lookups via the file number, rather than via the event number
      // - you can subsequently call any QA cut method, such as `Golden()`;
      //   although QA cut methods require an event number, no additional lookup
      //   query will be performed since it already has been done in QueryByFilenum
      inline bool QueryByFilenum(int runnum_, int filenum_);
      // get maximum file number for a given run (useful for QADB validation)
      inline int GetMaxFilenum(int runnum_);


      //.................................
      // Faraday Cup charge
      //`````````````````````````````````
      // -- accumulator
      // call this method after evaluating QA cuts (or at least after calling Query())
      // to add the current file's charge to the total charge;
      // - charge is accumulated per DST file, since the QA filters per DST file
      // - a DST file's charge is only accounted for if we have not counted it before
      inline void AccumulateCharge();
      // -- accessor
      // returns total accumlated charge that passed your QA cuts; call this
      // method after your event loop
      inline double GetAccumulatedCharge() { return chargeTotal; };
      // reset accumulated charge, if you ever need to
      inline void ResetAccumulatedCharge() { chargeTotal = 0; };



    private:

      int runnumMin, runnumMax;
      bool verbose = true;
      std::vector<std::string> qaJsonList;
      std::vector<std::string> chargeJsonList;

      rapidjson::Document qaTree;
      rapidjson::Document chargeTree;
      char readBuffer[65536];
      inline void chainTrees(std::vector<std::string> jsonList, rapidjson::Document & outTree);

      char runnumStr[32];
      char filenumStr[32];
      int runnum,filenum,evnumMin,evnumMax,evnumMinTmp,evnumMaxTmp;
      int defect;
      int sectorDefect[6];
      char sectorStr[8];
      std::string comment;
      double charge, chargeMin, chargeMax, chargeTotal;
      bool chargeCounted;
      std::vector<std::pair<int,int>> chargeCountedFiles;

      std::map<std::string,int> defectNameMap;
      int nbits;

      bool found;
      int asymMask;
      int mask;
  };



  //...............
  // constructor
  //```````````````
  QADB::QADB(int runnumMin_, int runnumMax_, bool verbose_) {

    runnumMin = runnumMin_;
    runnumMax = runnumMax_;
    verbose = verbose_;

    // get QADB directory
    if(verbose) std::cout << "\n[+] find json files" << std::endl;
    std::string dbDirN = getenv("QADB") ? getenv("QADB") : "";
    if(dbDirN.compare("")==0) {
      std::cerr << "ERROR: QADB environment variable not set" << std::endl;
      return;
    };
    dbDirN += "/qadb";
    if(verbose) std::cout << "QADB at " << dbDirN << std::endl;

    // get list of json files
    DIR * dbDir = opendir(dbDirN.c_str());
    struct dirent * dbDirent;
    while((dbDirent=readdir(dbDir))) {
      std::string qaDirN = std::string(dbDirent->d_name);
      if(qaDirN.find("qa.")!=std::string::npos) {
        qaJsonList.push_back(dbDirN+"/"+qaDirN+"/qaTree.json");
        chargeJsonList.push_back(dbDirN+"/"+qaDirN+"/chargeTree.json");
      };
    };
    closedir(dbDir);
    if(verbose) {
      std::cout << "qaTree files:" << std::endl;
      for(std::string str : qaJsonList) std::cout << " - " << str << std::endl;
      std::cout << "chargeTree files:" << std::endl;
      for(std::string str : chargeJsonList) std::cout << " - " << str << std::endl;
    };

    // read json files and concatenate, including only runs within specified
    // range [runnumMin,runnumMax]
    if(verbose) std::cout << "\n[+] read specified runs from json files" << std::endl;
    qaTree.SetObject();
    chargeTree.SetObject();
    this->chainTrees(qaJsonList,qaTree);
    this->chainTrees(chargeJsonList,chargeTree);
    if(verbose) {
      std::cout << "full list of runs read from QADB:" << std::endl;
      for(auto it=qaTree.MemberBegin(); it!=qaTree.MemberEnd(); ++it)
        std::cout << (it->name).GetString() << std::endl;
      std::cout << "-----\n";
    };


    // define bits (order and names must match those in ../../src/clasqa/Tools.groovy)
    nbits=0;
    defectNameMap.insert(std::pair<std::string,int>("TotalOutlier",nbits++));
    defectNameMap.insert(std::pair<std::string,int>("TerminalOutlier",nbits++));
    defectNameMap.insert(std::pair<std::string,int>("MarginalOutlier",nbits++));
    defectNameMap.insert(std::pair<std::string,int>("SectorLoss",nbits++));
    defectNameMap.insert(std::pair<std::string,int>("LowLiveTime",nbits++));
    defectNameMap.insert(std::pair<std::string,int>("Misc",nbits++));
    defectNameMap.insert(std::pair<std::string,int>("TotalOutlierFT",nbits++));
    defectNameMap.insert(std::pair<std::string,int>("TerminalOutlierFT",nbits++));
    defectNameMap.insert(std::pair<std::string,int>("MarginalOutlierFT",nbits++));
    defectNameMap.insert(std::pair<std::string,int>("LossFT",nbits++));


    // defect mask used for asymmetry analysis
    asymMask = 0;
    asymMask += 0x1 << Bit("TotalOutlier");
    asymMask += 0x1 << Bit("TerminalOutlier");
    asymMask += 0x1 << Bit("MarginalOutlier");
    asymMask += 0x1 << Bit("SectorLoss");

    // initialize local vars
    runnum = -1;
    filenum = -1;
    evnumMin = -1;
    evnumMax = -1;
    found = false;
    mask = 0;
    charge = 0;
    chargeTotal = 0;
    chargeCounted = false;
  };


  //....................................................................
  // concatenate trees from JSON files in jsonList to tree outTree
  //````````````````````````````````````````````````````````````````````
  // - includes only runs within specified range [runnumMin,runnumMax]
  void QADB::chainTrees(std::vector<std::string> jsonList, rapidjson::Document & outTree) {

    // loop through list of json files
    for(std::string jsonFileN : jsonList) {

      // open json file stream
      if(verbose) std::cout << "read json stream " << jsonFileN << std::endl;
      FILE * jsonFile = fopen(jsonFileN.c_str(),"r");
      rapidjson::FileReadStream jsonStream(jsonFile,readBuffer,sizeof(readBuffer));
      
      // parse stream to tmpTree
      rapidjson::Document tmpTree(rapidjson::kObjectType);
      if(tmpTree.ParseStream(jsonStream).HasParseError()) {
        std::cerr << "ERROR: QADB could not parse " << jsonFileN << std::endl;
        return;
      };

      // append to tmpTree to outTree
      for(auto it=tmpTree.MemberBegin(); it!=tmpTree.MemberEnd(); ++it) {
        runnum = atoi((it->name).GetString());
        if( ( runnumMin<0 && runnumMax<0) ||
            ( runnum>=runnumMin && runnum<=runnumMax)
        ) {
          if(verbose) std::cout << "- add run " << runnum << std::endl;
          rapidjson::Value objKey, objVal;
          objKey.CopyFrom(it->name,outTree.GetAllocator());
          objVal.CopyFrom(it->value,outTree.GetAllocator());
          outTree.AddMember(objKey,objVal,outTree.GetAllocator());
        };
      };

      // close json file
      fclose(jsonFile);
    };

    // reset runnum
    runnum = -1;
  };


  //...................................
  // QA for spin asymmetry analysis
  //```````````````````````````````````
  bool QADB::OkForAsymmetry(int runnum_, int evnum_) {

    // perform lookup
    bool foundHere = this->Query(runnum_,evnum_);
    if(!foundHere) return false;

    // check for bits which will always cause the file to be rejected 
    // (asymMask is defined in the constructor)
    if( defect & asymMask ) return false;

    // special cases for `Misc` bit
    if(this->HasDefect("Misc")) {

      // check if this is a run on the list of runs with a large fraction of
      // events with undefined helicity; if so, accept this run, since none of
      // these files are marked with `Misc` for any other reasons
      if(runnum_==5128 ||
         runnum_==5129 ||
         runnum_==5130 ||
         runnum_==5158 ||
         runnum_==5159 ||
         runnum_==5160 ||
         runnum_==5163 ||
         runnum_==5165 ||
         runnum_==5166 ||
         runnum_==5167 ||
         runnum_==5168 ||
         runnum_==5169 ||
         runnum_==5180 ||
         runnum_==5181 ||
         runnum_==5182 ||
         runnum_==5183 ||
         runnum_==5567) return true;

      // check if this run had an FADC failure; there is no indication spin
      // asymmetries are impacted by this issue
      else if(runnum_>=6736 && runnum_<=6757) return true;

      // otherwise, this file fails the QA
      else return false;
    };

    // otherwise, this file passes the QA
    return true;
  };


  //...............................
  // user-defined custom QA cuts
  //```````````````````````````````
  void QADB::SetMaskBit(const char * defectName, bool state) {
    int defectBit = this->Bit(defectName);
    if(defectBit<0 || defectBit>=nbits)
      std::cerr << "ERROR: QADB::SetMaskBit called for unknown bit" << std::endl;
    else {
      mask &= ~(0x1 << defectBit);
      if(state) mask |= (0x1 << defectBit);
    };
  };
  bool QADB::Pass(int runnum_, int evnum_) {
    bool foundHere = this->Query(runnum_,evnum_);
    return foundHere && !(defect & mask);
  }




  //.................
  // accessors
  //`````````````````
  int QADB::GetDefect(int sector) {
    if(!found) return -1;
    if(sector==0) return defect;
    else if(sector>=1 && sector<=6) return sectorDefect[sector-1];
    else {
      std::cerr << "ERROR: bad sector number for QADB::GetDefect" << std::endl;
      return -1;
    };
  };
  int QADB::Bit(const char * defectName) {
    int defectBit;
    try { defectBit = defectNameMap.at(std::string(defectName)); }
    catch(const std::out_of_range & e) {
      std::cerr << "ERROR: QADB::Bit() unknown defectName" << std::endl;
      return -1;
    };
    return defectBit;
  };
    



  //.....................................................................
  // query qaTree to get QA info for this run number and event number
  // - a lookup is only performed if necessary: if the run number changes
  //   or if the event number goes outside of the range of the file which
  //   most recently queried
  //`````````````````````````````````````````````````````````````````````
  bool QADB::Query(int runnum_, int evnum_) {

    // if the run number changed, or if the event number is outside the range
    // of the previously queried file, perform a new lookup
    if( runnum_ != runnum ||
        ( runnum_ == runnum && (evnum_ < evnumMin || evnum_ > evnumMax ))
    ) {
      // reset vars
      runnum = runnum_;
      filenum = -1;
      evnumMin = -1;
      evnumMax = -1;
      charge = -1;
      found = false;

      // search for file which contains this event
      sprintf(runnumStr,"%d",runnum);
      if(qaTree.HasMember(runnumStr)) {
        auto runTree = qaTree[runnumStr].GetObject();

        for(auto it=runTree.MemberBegin(); it!=runTree.MemberEnd(); ++it) {
          auto fileTree = (it->value).GetObject();
          evnumMinTmp = fileTree["evnumMin"].GetInt();
          evnumMaxTmp = fileTree["evnumMax"].GetInt();
          if( evnum_ >= evnumMinTmp && evnum_ <= evnumMaxTmp ) {
            this->QueryByFilenum(runnum_,atoi((it->name).GetString()));
            break;
          };
        };
      };

      // print a warning if a file was not found for this event
      // - this warning is suppressed for 'tag1' events
      if(!found && runnum_!=0) {
        std::cerr << "WARNING: QADB::Query could not find runnum=" <<
          runnum_ << " evnum=" << evnum_ << std::endl;
      };
    };

    // result of query
    return found;
  };

  // ------ query by file number
  bool QADB::QueryByFilenum(int runnum_, int filenum_) {

    // if the run number or file number changed, perform new lookup
    if( runnum_ != runnum || filenum_ != filenum) {
      
      // reset vars
      runnum = runnum_;
      filenum = filenum_;
      evnumMin = -1;
      evnumMax = -1;
      charge = -1;
      found = false;

      sprintf(runnumStr,"%d",runnum);
      if(qaTree.HasMember(runnumStr)) {
        auto runTree = qaTree[runnumStr].GetObject();
        sprintf(filenumStr,"%d",filenum);
        if(runTree.HasMember(filenumStr)) {
          auto fileTree = runTree[filenumStr].GetObject();
          evnumMin = fileTree["evnumMin"].GetInt();
          evnumMax = fileTree["evnumMax"].GetInt();
          comment = fileTree["comment"].GetString();
          defect = fileTree["defect"].GetInt();
          auto sectorTree = fileTree["sectorDefects"].GetObject();
          for(int s=0; s<6; s++) {
            sprintf(sectorStr,"%d",s+1);
            const rapidjson::Value& defList = sectorTree[sectorStr];
            sectorDefect[s] = 0;
            for(rapidjson::SizeType i=0; i<defList.Size(); i++) {
              sectorDefect[s] += 0x1 << defList[i].GetInt();
            };
          };
          chargeMin = chargeTree[runnumStr][filenumStr]["fcChargeMin"].GetDouble();
          chargeMax = chargeTree[runnumStr][filenumStr]["fcChargeMax"].GetDouble();
          charge = chargeMax - chargeMin;
          chargeCounted = false;
          found = true;
        };
      };

      // print a warning if a file was not found for this event
      // - this warning is suppressed for 'tag1' events
      if(!found && runnum!=0) {
        std::cerr << "WARNING: QADB::QueryByFilenum could not find runnum=" <<
          runnum_ << " filenum=" << filenum_ << std::endl;
      };
    };

    // result of query
    return found;
  };

  // ----- return maximum filenum for a given runnum
  int QADB::GetMaxFilenum(int runnum_) {
    int maxFilenum=0;
    sprintf(runnumStr,"%d",runnum_);
    auto runTree = qaTree[runnumStr].GetObject();
    for(auto it=runTree.MemberBegin(); it!=runTree.MemberEnd(); ++it) {
      maxFilenum = atoi((it->name).GetString()) > maxFilenum ?
                   atoi((it->name).GetString()) : maxFilenum;
    };
    return maxFilenum;
  };




  //.................................
  // Faraday Cup charge accumulator
  //`````````````````````````````````
  void QADB::AccumulateCharge() {
    if(!chargeCounted) {
      if(
        find(
          chargeCountedFiles.begin(),
          chargeCountedFiles.end(),
          std::pair<int,int>(runnum,filenum)
        ) == chargeCountedFiles.end()
      ) {
        chargeTotal += charge;
        chargeCountedFiles.push_back(std::pair<int,int>(runnum,filenum));
      };
      chargeCounted = true;
    };
  };

};


#endif
