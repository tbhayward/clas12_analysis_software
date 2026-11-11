#!/usr/bin/env python3
"""Exercise the production C++ CSV header writer without loading ROOT."""
import csv,io,subprocess,tempfile
from pathlib import Path
SOURCE=Path(__file__).resolve().parent

def function(text,name):
    # Find the actual declaration rather than a preceding unrelated helper.
    import re
    match=re.search(r'static [^\n]+\b'+name+r'\(',text)
    assert match,name
    start=match.start();brace=text.index('{',match.end());depth=1;end=brace+1
    while depth:
        depth+=(text[end]=='{')-(text[end]=='}');end+=1
    return text[start:end]

def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--csv',type=Path);args=ap.parse_args()
    text=(SOURCE/'cross_sections.cpp').read_text()
    names=['split_csv_line','trim','unquote','quote_if_needed','join_csv_line','serialize_csv_header']
    implementation='\n'.join(function(text,n) for n in names)
    main_cpp='''
int main(){std::string line;while(std::getline(std::cin,line)){
 auto fields=split_csv_line(line);
 if(unquote(trim(fields.back()))!="cross section correction convention")fields.push_back("cross section correction convention");
 const auto output=serialize_csv_header(fields);
 if(serialize_csv_header(split_csv_line(output))!=output)return 2;
 std::cout<<output<<"\\n";
}}
'''
    headers=[['bin index','norm, Fa18 Inb','acceptance corrected yield, ep->epg, exp, 2018 (10.6 GeV), unpol','name with "quotes"','plain']]
    fixture=SOURCE/'imports/all_bin_v3.csv'
    if fixture.exists():headers.append(next(csv.reader(fixture.open())))
    # Production schema from the initializer, including every normalization name.
    init=SOURCE/'initialize_pass2_csv.cpp'
    if init.exists():
        import re
        header=re.findall(r'H\.push_back\("([^"\n]+)"\)',init.read_text())
        if header:headers.append(header)
    if args.csv:
        headers.append(next(csv.reader(args.csv.open())))
    encoded=[]
    for h in headers:
        stream=io.StringIO();csv.writer(stream,lineterminator='').writerow(h);encoded.append(stream.getvalue())
    with tempfile.TemporaryDirectory() as tmp:
        cpp=Path(tmp)/'headers.cpp';exe=Path(tmp)/'headers'
        cpp.write_text('#include <string>\n#include <vector>\n#include <sstream>\n#include <iostream>\n#include <cctype>\n'+implementation+main_cpp)
        subprocess.run(['c++','-std=c++17',str(cpp),'-o',str(exe)],check=True)
        outputs=subprocess.check_output([str(exe)],input='\n'.join(encoded)+'\n',text=True).splitlines()
    for old,new in zip(headers,outputs):
        actual=next(csv.reader([new]));assert actual==old+([] if old[-1]=='cross section correction convention' else ['cross section correction convention']),(old,actual)
    assert 'out_lines.push_back(serialize_csv_header(header));' in text
    print('PASS: production C++ header serialization preserves all column names, embedded quotes, normalization columns, and repeat-write identity')

if __name__=='__main__':main()
