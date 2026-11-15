#include "production_particle_response.h"
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TChainElement.h>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <memory>
#include <algorithm>
namespace {
using json=nlohmann::json;namespace fs=std::filesystem;
std::vector<std::unique_ptr<TChain>>& owned(){static std::vector<std::unique_ptr<TChain>> x;return x;}
std::string quote(const std::string& s){std::string x="'";for(char c:s)x+=c=='\''?"'\\''":std::string(1,c);return x+"'";}
json read(const fs::path& p){std::ifstream f(p);if(!f)throw std::runtime_error("Cannot read "+p.string());json x;f>>x;return x;}
std::vector<std::string> files(TTree* tree){std::vector<std::string> x;if(auto* c=dynamic_cast<TChain*>(tree)){for(auto* o:*c->GetListOfFiles())x.push_back(fs::absolute(o->GetTitle()).lexically_normal().string());}else if(tree->GetCurrentFile())x.push_back(fs::absolute(tree->GetCurrentFile()->GetName()).lexically_normal().string());else throw std::runtime_error("Response needs file-backed reconstructed trees");return x;}
bool command(const std::string& s){std::cout<<"[particle_response] "<<s<<std::endl;return std::system(s.c_str())==0;}
}
bool prepare_production_particle_response(const ProductionParticleResponseOptions& o,const std::vector<std::pair<std::string,ResponseTreeMap*>>& maps){try{
 if(!o.enabled())return true;if(o.workers<1||o.workers>10)throw std::runtime_error("Response workers must be in [1,10]");fs::create_directories(o.output_dir);json inv=json::array();
 for(const auto& category:maps)for(const auto& item:*category.second)inv.push_back({{"sample",category.first},{"tag",item.first},{"tree",item.second->GetName()},{"files",files(item.second)}});
 const fs::path inventory=fs::path(o.output_dir)/"input_inventory.json";{std::ofstream f(inventory);f<<inv.dump(2);}
 const std::string cmd="python3 prepare_production_particle_response.py --inventory "+quote(inventory.string())+(o.fit_response ? " --fit-response --initial-cuts "+quote((fs::path(o.output_dir)/"initial_cuts/combined_cuts.json").string()) : " --response-run "+quote(fs::absolute(o.fitted_run).string()))+" --global-cuts output/jsons/global_cuts_config.json --scratch-root "+quote(o.scratch_root)+" --output "+quote(o.output_dir)+" --workers "+std::to_string(o.workers);
 if(!command(cmd))return false;
 const fs::path manifest=fs::absolute(fs::path(o.output_dir)/"input_manifest.json");const auto m=read(manifest);std::vector<std::pair<TTree**,std::unique_ptr<TChain>>> prepared;
 for(const auto& category:maps){if(category.first=="data"||category.first=="pi0_data"||category.first.find("gen")!=std::string::npos)continue;
  for(auto& item:*category.second){auto chain=std::make_unique<TChain>(item.second->GetName());for(const auto& file:files(item.second)){if(!m.at("routes").contains(file))throw std::runtime_error("Unrouted reconstructed input: "+file);auto target=m.at("routes").at(file).at("target").get<std::string>();if(!fs::is_regular_file(target))throw std::runtime_error("Missing response file: "+target);if(chain->Add(target.c_str())==0)throw std::runtime_error("Cannot chain response file: "+target);}
   if(chain->GetEntries()!=item.second->GetEntries())throw std::runtime_error("Response changed entry count for "+item.first);prepared.emplace_back(&item.second,std::move(chain));
  }
 }
 for(auto& item:prepared){*item.first=item.second.get();owned().push_back(std::move(item.second));}
 if(setenv("DVCS_PARTICLE_RESPONSE_MANIFEST",manifest.c_str(),1)!=0)throw std::runtime_error("Cannot install response input manifest");
 std::cout<<"[particle_response] All reconstructed maps now use frozen full-statistics response; DATA and generated maps preserved.\n";return true;
}catch(const std::exception& e){std::cerr<<"[particle_response] FATAL: "<<e.what()<<'\n';return false;}}
bool write_production_particle_response_note(const ProductionParticleResponseOptions& o){if(!o.enabled())return true;return command("python3 plot_production_particle_response.py --response-output "+quote(o.output_dir)+" --cuts "+quote((fs::path(o.output_dir)/"final_cuts/combined_cuts.json").string())+" --output "+quote((fs::path(o.output_dir)/"analysis_note").string())+" --workers "+std::to_string(o.workers));}

bool finalize_production_particle_response(const ProductionParticleResponseOptions& o){if(!o.enabled())return true;return command("python3 finalize_production_particle_response.py --response-output "+quote(o.output_dir)+" --csv output/csvs/dvcs_pass2_analysis.csv --cuts "+quote((fs::path(o.output_dir)/"final_cuts/combined_cuts.json").string()));}
