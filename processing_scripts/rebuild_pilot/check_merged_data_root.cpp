#include <TFile.h>
#include <TTree.h>
#include <cstdlib>
#include <iostream>
#include <memory>
int main(int argc,char**argv){
 if(argc!=4){std::cerr<<"usage: check_merged_data_root ROOT expected_rows expected_branches\n";return 2;}
 std::unique_ptr<TFile> f(TFile::Open(argv[1],"READ"));
 if(!f||f->IsZombie()){std::cerr<<"cannot open ROOT\n";return 2;}
 TTree* t=nullptr;f->GetObject("DataEvents",t);
 if(!t){std::cerr<<"missing DataEvents tree\n";return 2;}
 auto rows=t->GetEntries();auto branches=t->GetListOfBranches()->GetEntries();
 std::cout<<"[MERGED_CHECK] rows="<<rows<<" expected="<<argv[2]<<" branches="<<branches<<" expected_branches="<<argv[3]<<"\n";
 if(rows!=std::atoll(argv[2])||branches!=std::atoi(argv[3]))return 1;
 for(int i=1;i<=branches;i++){char name[20];std::snprintf(name,sizeof(name),"col%02d",i);if(!t->GetBranch(name)){std::cerr<<"missing "<<name<<"\n";return 1;}}
 return 0;
}
