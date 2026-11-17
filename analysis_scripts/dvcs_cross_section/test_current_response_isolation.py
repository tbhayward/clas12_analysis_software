#!/usr/bin/env python3
"""ROOT-free native regression checks for nominal current-response isolation."""
from pathlib import Path
import re, subprocess, tempfile, shutil
ROOT=Path(__file__).resolve().parent

def definition(text,name,kind='function'):
    pattern=(r'static [^\n]*\b'+name+r'\(') if kind=='function' else (r'struct '+name+r'\s*\{')
    m=re.search(pattern,text);assert m,name
    start=text.index('{',m.start());depth=1;i=start+1
    while depth:
        depth+=(text[i]=='{')-(text[i]=='}');i+=1
    return text[m.start():i]+(';' if kind=='struct' else '')+'\n'

def main():
    source=(ROOT/'current_systematics.cpp').read_text()
    header=(ROOT/'total_counts.h').read_text();variation=(ROOT/'cut_variation_runner.cpp').read_text()
    assert 'bool write_current_nuisance_responses = false;' in header
    assert 'count_opts.write_current_nuisance_responses = false;' in variation
    writer=definition((ROOT/'total_counts.cpp').read_text(),'write_current_nuisance_response_csv')
    assert writer.count('row < source_csv.rows.size()')==3 # gen, reconstructed and metadata
    assert 'const WeightedHelCounts zero{}' in writer
    code='''#include <map>
#include <unordered_map>
#include <vector>
#include <array>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include "current_response_provenance.h"
#include "current_systematics.h"
namespace fs=std::filesystem;
'''
    code+=definition(source,'Csv','struct')
    for name in ['trim','normalize_csv_field','split_csv','escape_csv','read_csv','num','tuple_first']:
        code+=definition(source,name)
    code+=definition(source,'YieldState','struct')
    code+='using SampleMap=std::map<std::string,YieldState>; using PeriodMap=std::map<std::string,SampleMap>; using ResponseByRow=std::map<int,PeriodMap>;\n'
    for name in ['read_responses','sample_nom','validate_response_provenance','validate_nominal_response_yields','sample_delta','nuisances_for_period','unfolded_for_period']:
        code+=definition(source,name)
    code+=definition(source,'BinResult','struct')
    for name in ['period_result','combined_10p6_result']:code+=definition(source,name)
    code+=r'''
int main(){
    auto fails=[](auto f){try{f();return false;}catch(const std::runtime_error&){return true;}};
    PeriodMap p;
    const std::array<std::string,5> periods={"Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out","Sp19 Inb"};
    for(auto per:periods){
        p[per]["epg_data"].nominal=100; p[per]["eppi0_data"].nominal=10;
        p[per]["eppi0_rec"].nominal=20; p[per]["misid_rec"].nominal=5;
        p[per]["epg_rec"].nominal=50; p[per]["epg_gen"].nominal=100;
        p[per]["epg_data"].delta["data:shared"]=2;
    }
    double u=0.; assert(unfolded_for_period(p,"Fa18 Inb",nullptr,0.,u)); assert(std::fabs(u-195.)<1e-12);
    auto one=period_result(p,"Fa18 Inb",0.);auto all=combined_10p6_result(p);
    assert(one.valid && all.valid);assert(std::fabs(one.scalar-all.scalar)<1e-12);
    auto sp=period_result(p,"Sp19 Inb",.0343);
    assert(sp.valid && std::fabs(sp.scalar-std::hypot(one.scalar,.0343))<1e-12);
    auto zero=p;zero["Fa18 Inb"]["eppi0_data"].nominal=0.;zero["Fa18 Inb"]["eppi0_rec"].nominal=0.;
    assert(unfolded_for_period(zero,"Fa18 Inb",nullptr,0.,u) && u==200.);
    zero["Fa18 Inb"].erase("eppi0_data");assert(!unfolded_for_period(zero,"Fa18 Inb",nullptr,0.,u));
    auto large=p;large["Fa18 Inb"]["epg_data"].delta["data:shared"]=300.;
    const std::string shared="data:shared"; assert(unfolded_for_period(large,"Fa18 Inb",&shared,-1.,u) && u<0.);
    assert(period_result(large,"Fa18 Inb",0.).scalar>1.);
    auto bad=p;bad["Fa18 Inb"]["epg_rec"].delta["mc:bad"]=60.;
    assert(fails([&]{period_result(bad,"Fa18 Inb",0.);}));
    Csv csv;auto column=[&](const std::string& name,double value){csv.index[name]=csv.header.size();csv.header.push_back(name);if(csv.rows.empty())csv.rows.emplace_back();csv.rows[0].push_back(std::to_string(value));};
    for(auto name:current_response_provenance::bin_columns())column(name,1.);
    for(auto per:periods){
        for(auto ch:{"ep->epg","ep->eppi0"})for(auto topo:{"(FD, FD)","(CD, FD)","(CD, FT)"})
            column(std::string("normalized raw yield, ")+ch+", "+topo+", exp, "+per+", unpol",std::string(topo)=="(FD, FD)"?(std::string(ch)=="ep->epg"?100.:10.):0.);
        column("reconstructed current corrected yield, ep->epg, mc, "+per,50.);
        column("reconstructed current corrected yield, ep->eppi0, mc, "+per,20.);
        column("reconstructed current corrected yield, ep->eppi0->epg, mc, "+per,5.);
        column("generated yield, ep->epg, mc, "+per,100.);
    }
    ResponseByRow response;response[0]=p;validate_nominal_response_yields(csv,response);
    response[0]["Fa18 Inb"]["epg_data"].nominal=101.;
    assert(fails([&]{validate_nominal_response_yields(csv,response);}));response[0]=p;
    CurrentSystematicsOptions opt;opt.nuisance_response_csv="responses.csv";opt.combined_cuts_json="cuts.json";opt.current_response_model_json="model.json";
    {std::ofstream f("responses.csv");f<<"response";}{std::ofstream f("cuts.json");f<<"cuts";}{std::ofstream f("model.json");f<<"model";}
    assert(fails([&]{validate_response_provenance(csv,opt);}));
    {std::ofstream f("responses.csv.provenance.csv");f<<"row";for(auto name:current_response_provenance::bin_columns())f<<','<<name;
      f<<",cuts_hash,calibration_hash,response_hash\n0";for(int i=0;i<8;++i)f<<",1";
      f<<','<<current_response_provenance::file_hash("cuts.json")<<','<<current_response_provenance::file_hash("model.json")<<','<<current_response_provenance::file_hash("responses.csv")<<'\n';}
    validate_response_provenance(csv,opt);
    csv.rows[0][csv.index.at("Q2min")]="2";assert(fails([&]{validate_response_provenance(csv,opt);}));csv.rows[0][csv.index.at("Q2min")]="1";
    {std::ofstream f("responses.csv");f<<"overwritten variation";}assert(fails([&]{validate_response_provenance(csv,opt);}));
    std::cout<<"PASS: export opt-in and variation isolation, bin/file provenance, nominal-yield agreement, explicit zero vs missing control, shared-period and Sp19 transfer responses, signed zero crossing and invalid excursion rejection\n";
}
'''
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp);(tmp/'test.cpp').write_text(code)
        subprocess.run([shutil.which('c++') or 'g++','-std=c++17','-O0','-I',str(ROOT),str(tmp/'test.cpp'),'-o',str(tmp/'test')],check=True)
        subprocess.run([str(tmp/'test')],cwd=tmp,check=True)
if __name__=='__main__':main()
