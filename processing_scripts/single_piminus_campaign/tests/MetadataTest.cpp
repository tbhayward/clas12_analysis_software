#include "../src/RunMetadata.h"
#include <iostream>
int main(int argc,char** argv){try{
 double a=targetPolarization(argv[1],"ND3",16360),b=targetPolarization(argv[1],"ND3",16962),c=targetPolarization(argv[1],"ND3",17704);
 if(std::abs(a-.239)>1e-10||std::abs(b+.299)>1e-10||std::abs(c-.324)>1e-10)throw std::runtime_error("Wrong polarization values");
 if(targetPolarization(argv[1],"D2",6265)!=0)throw std::runtime_error("Unpolarized target");
 std::cout<<"Actual mixed-format CSV parsed correctly for all three selected ND3 runs and RGB\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
