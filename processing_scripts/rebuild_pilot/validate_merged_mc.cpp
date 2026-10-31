#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include "validate_mc_companion.h"
int main(int argc,char**argv){
 if(argc!=4){std::cerr<<"usage: validate_merged_mc merged.root expected_entries expected_rec\n";return 2;}
 const long long expected=std::stoll(argv[2]), expectedRec=std::stoll(argv[3]);
 TFile f(argv[1],"READ");if(f.IsZombie()||f.TestBit(TFile::kRecovered)){std::cerr<<"BAD_ROOT\n";return 3;}
 auto*t=dynamic_cast<TTree*>(f.Get("MCEvents"));if(!t){std::cerr<<"MISSING_TREE\n";return 4;}
 for(const char* name:{"event_sequence","gen_kinematics_valid","rec_kinematics_valid","gen_t1","rec_t1"})if(!t->GetBranch(name)){std::cerr<<"MISSING_BRANCH "<<name<<"\n";return 5;}
 const long long n=t->GetEntries();if(n!=expected){std::cerr<<"BAD_ENTRIES actual="<<n<<" expected="<<expected<<"\n";return 6;}
 double seq=0,gen=0,rec=0;t->SetBranchAddress("event_sequence",&seq);t->SetBranchAddress("gen_kinematics_valid",&gen);t->SetBranchAddress("rec_kinematics_valid",&rec);
 long long nr=0,ng=0,previous=0;for(long long i=0;i<n;i++){
  if(t->GetEntry(i)<=0){std::cerr<<"READ_FAILURE entry="<<i<<"\n";return 7;}
  if(!std::isfinite(seq)||std::llround(seq)<1){std::cerr<<"INVALID_SEQUENCE entry="<<i<<"\n";return 8;}
  long long s=std::llround(seq);if(s==1){/* Sequence restarts at each merged input file. */}
  else if(s!=previous+1){std::cerr<<"BROKEN_SEQUENCE entry="<<i<<" seq="<<s<<" previous="<<previous<<"\n";return 9;}
  previous=s;
  if(gen>0.5)ng++;if(rec>0.5)nr++;
 }
 if(ng!=n||nr!=expectedRec){std::cerr<<"BAD_COUNTS gen="<<ng<<" rec="<<nr<<" expected_rec="<<expectedRec<<"\n";return 10;}
 if(!validate_mc_companion(f,t))return 11;
 std::cout<<"MERGED_VALID entries="<<n<<" gen="<<ng<<" rec="<<nr<<"\n";return 0;
}
