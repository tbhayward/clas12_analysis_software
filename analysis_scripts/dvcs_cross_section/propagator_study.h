#ifndef PROPAGATOR_STUDY_H
#define PROPAGATOR_STUDY_H

#include <map>
#include <string>
#include <vector>

class TTree;

// Build dvcsgen ycol-mirror efficiency canvases (phi-dependent) for standard DVCS bins.
//
// Efficiency is defined as:
//
//   eff(phi) = N(pass global cuts WITH ycol cut) / N(pass global cuts WITHOUT ycol cut)
//
// where the ycol cut is implemented inside global_cuts as:
//   P1_pos > cfg.dvcsgen_ycol_cut
// with P1_pos computed from beam energy and electron/photon kinematics.
//
// Inputs:
//   csv_main            : path to dvcs_pass2_analysis.csv (bin definitions)
//   dataTreesByLabel    : map from label -> vector<TTree*> (DATA trees)
//   combined_cuts_json  : path to output/jsons/combined_cuts.json (3-sigma cuts)
//   global_cuts_json    : path to output/jsons/global_cuts_config.json (provenance config)
//   out_root_dir        : e.g. "output/propagator_study"
//
// Returns true on success; throws std::runtime_error on fatal configuration/schema/IO issues.
bool run_propagator_study(const std::string &csv_main,
                          const std::map<std::string, std::vector<TTree *>> &dataTreesByLabel,
                          const std::string &combined_cuts_json,
                          const std::string &global_cuts_json,
                          const std::string &out_root_dir);

#endif