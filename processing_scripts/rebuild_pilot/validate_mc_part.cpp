#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <cmath>
int main(int argc,char**argv){
 if(argc!=3){std::cerr<<"usage: validate_mc_part file.root expected_entries (-1 unknown)\n";return 2;}
 TFile f(argv[1],"READ");if(f.IsZombie()||f.TestBit(TFile::kRecovered)){std::cerr<<"BAD_ROOT\n";return 3;}
 auto*t=dynamic_cast<TTree*>(f.Get("MCEvents"));if(!t){std::cerr<<"MISSING_MCEvents\n";return 4;}
 long long n=t->GetEntries(),expect=std::stoll(argv[2]);
 if(n<=0||(expect>=0&&n!=expect)){std::cerr<<"BAD_ENTRIES actual="<<n<<" expected="<<expect<<"\n";return 5;}
 for(auto name:{"event_sequence","gen_kinematics_valid","rec_kinematics_valid","gen_t1","rec_t1"})if(!t->GetBranch(name)){std::cerr<<"MISSING_BRANCH "<<name<<"\n";return 6;}
 double sequence=0,gen=0,rec=0; t->SetBranchAddress("event_sequence",&sequence);t->SetBranchAddress("gen_kinematics_valid",&gen);t->SetBranchAddress("rec_kinematics_valid",&rec);
 long long ng=0,nr=0;for(long long i=0;i<n;i++){if(t->GetEntry(i)<=0){std::cerr<<"READ_ERROR entry="<<i<<"\n";return 7;}if(!std::isfinite(sequence)||std::llround(sequence)!=i+1){std::cerr<<"BAD_SEQUENCE entry="<<i<<" value="<<sequence<<"\n";return 8;}if(gen>0.5)ng++;if(rec>0.5)nr++;}
 if(ng!=n){std::cerr<<"GEN_MISMATCH gen="<<ng<<" total="<<n<<"\n";return 9;}
 std::cout<<"VALID entries="<<n<<" gen="<<ng<<" rec="<<nr<<"\n";return 0;
}
