#ifndef PRODUCTION_PARTICLE_RESPONSE_H
#define PRODUCTION_PARTICLE_RESPONSE_H
#include <map>
#include <string>
#include <vector>
class TTree;
using ResponseTreeMap=std::map<std::string,TTree*>;
struct ProductionParticleResponseOptions {
 bool fit_response=false;
 std::string fitted_run;
 bool enabled() const {return fit_response || !fitted_run.empty();}
 std::string scratch_root="/scratch/thayward/dvcs_particle_response/production";
 std::string output_dir="output/particle_response";
 int workers=10;
};
// Preparation preserves original maps until every cache file is validated.
// One canonical manifest is then installed for Python and all reconstructed maps.
bool prepare_production_particle_response(const ProductionParticleResponseOptions& options,
 const std::vector<std::pair<std::string,ResponseTreeMap*>>& maps);
bool write_production_particle_response_note(const ProductionParticleResponseOptions& options);
bool finalize_production_particle_response(const ProductionParticleResponseOptions& options);
#endif
