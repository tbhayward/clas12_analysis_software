#!/usr/bin/env python3
"""Check response owners close at function exit, including failure returns."""
import os,subprocess,tempfile
from pathlib import Path
SOURCE=Path(__file__).resolve().parent
MOCK=r'''
#include <cassert>
#include <cstdlib>
#include <stdexcept>
static int live=0;
class TTree { public: virtual ~TTree()=default; };
#include "production_particle_response.h"
struct Chain: TTree {Chain(){++live;}~Chain(){--live;}};
static int run(bool fail){
 ResponseTreeOwners owners;
 owners.push_back(std::make_unique<Chain>());
 owners.push_back(std::make_unique<Chain>());
 assert(live==2);
 if(fail)return 1;
 return 0;
}
static void throws(){ResponseTreeOwners owners;owners.push_back(std::make_unique<Chain>());throw std::runtime_error("failure");}
int main(){
 std::atexit([]{assert(live==0);});
 assert(run(false)==0);assert(live==0);
 assert(run(true)==1);assert(live==0);
 try{throws();}catch(const std::runtime_error&){}assert(live==0);
}
'''
def main():
    text=(SOURCE/'production_particle_response.cpp').read_text();main=(SOURCE/'main.cpp').read_text()
    assert 'owned().push_back' not in text and 'static std::vector<std::unique_ptr<TChain>>' not in text
    assert 'ResponseTreeOwners particle_response_trees;' in main
    assert 'response_maps, particle_response_trees)' in main and 'owners.push_back(std::move(item.second))' in text
    with tempfile.TemporaryDirectory() as td:
        source=Path(td)/'lifetime.cpp';exe=Path(td)/'lifetime';source.write_text(MOCK)
        subprocess.run([os.environ.get('CXX','c++'),'-std=c++17','-I'+str(SOURCE),str(source),'-o',str(exe)],check=True)
        subprocess.run([str(exe)],check=True)
    print('PASS: scoped response-chain ownership, normal/error returns and exception cleanup before static shutdown')
if __name__=='__main__':main()
