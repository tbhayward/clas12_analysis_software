#!/usr/bin/env python3
"""Compile actual current-key parsing and fitting helpers without ROOT."""
from pathlib import Path
import re,subprocess,tempfile
ROOT=Path(__file__).resolve().parent

def definition(text,name,kind='function'):
 pattern=(r'static [^\n]*\b'+name+r'\(') if kind=='function' else (r'struct '+name+r'\s*\{')
 m=re.search(pattern,text);assert m,name
 brace=text.index('{',m.start());depth=1;i=brace+1
 while depth:depth+=(text[i]=='{')-(text[i]=='}');i+=1
 return text[m.start():i]+(';' if kind=='struct' else '')

def main():
 s=(ROOT/'current_dependence.cpp').read_text()
 code='''#include <map>
#include <set>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <limits>
#include <cmath>
#include <cassert>
#include <stdexcept>
struct PeriodTags {std::string display;};
static PeriodTags parse_period_from_key(const std::string& key){return {key.find("Sp18_out")!=std::string::npos?"Sp18 Out":"Fa18 Inb"};}
static void fatal(const std::string& s){throw std::runtime_error(s);}
'''
 for name in ['lower_ascii','has_substr','reference_current_nA','has_exact_current_token','parse_current_from_key']:code+=definition(s,name)+'\n'
 code+=definition(s,'FitResult','struct')+'\n'
 code+=definition(s,'CurrentPoint','struct')+'\n'
 code+=definition(s,'has_zero_and_nonzero_current_points')+'\n'
 for name in ['weighted_linear_fit','rel_at_current']:code+=definition(s,name)+'\n'
 loader=(ROOT/'load_trees.cpp').read_text()
 start=loader.index('const std::map<std::string, std::vector<std::string>> current_configs =')
 stop=loader.index('};',start)+2
 configs=loader[start:stop]
 assert 'current == "no_bkg" ? "0nA" : current' in loader
 code+='int main(){\n'+configs+'\n'
 code+='''assert(!has_zero_and_nonzero_current_points({}));
 assert(!has_zero_and_nonzero_current_points({{50,.8,.01}}));
 assert(!has_zero_and_nonzero_current_points({{45,.82,.01},{50,.8,.01}}));
 assert(has_zero_and_nonzero_current_points({{0,1.,.01},{50,.8,.01}}));
'''

 code+='''for(const auto& period:current_configs){
  std::set<int> settings;
  for(const auto& current:period.second){
   const std::string label=current=="no_bkg"?"0nA":current;
   settings.insert(parse_current_from_key("DVCS_"+period.first+"_"+label+"_rec_current"));
  }
  assert(settings.count(0)==1);
  assert(settings.size()==period.second.size());
  assert(settings.size()>=2);
 }

 for(int n:{0,5,35,40,45,50,55,75}){
  auto key=std::string("DVCS_Fa18_inb_")+std::to_string(n)+"nA_rec_current";
  assert(parse_current_from_key(key)==n);
 }
 assert(parse_current_from_key("DVCS_Fa18_inb_nobkg_rec")==0);
 assert(parse_current_from_key("DVCS_Fa18_inb_rec")==50);
 assert(parse_current_from_key("DVCS_Sp18_out_rec")==45);
 assert(!has_exact_current_token("45na",5));
 assert(!has_exact_current_token("50na",0));
 std::vector<double>x,y,sy;
 for(int n:{45,50,55}){
  x.push_back(parse_current_from_key(std::string("DVCS_Fa18_inb_")+std::to_string(n)+"nA_rec_current"));
  y.push_back(1.-.004*n);sy.push_back(.01);
 }
 auto f=weighted_linear_fit(x,y,sy);
 assert(std::abs(f.m+.004)<1e-10);
 assert(std::abs(rel_at_current(50,f)-.8)<1e-10);
}
'''
 with tempfile.TemporaryDirectory() as t:
  p=Path(t);(p/'test.cpp').write_text(code)
  subprocess.run(['c++','-std=c++17',str(p/'test.cpp'),'-o',str(p/'test')],check=True)
  subprocess.run([str(p/'test')],check=True)
 print('PASS: all five scans include distinct 0/nonzero currents; exact labels, legacy defaults, recovered nonzero slope')
if __name__=='__main__':main()
