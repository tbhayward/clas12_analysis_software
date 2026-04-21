#include "branch_plot_config.h"

#include <TMath.h>
#include <map>
#include <string>

const std::map<std::string, BranchHistConfig>& getBranchPlotConfigs() {
    static const std::map<std::string, BranchHistConfig> cfg = {
        {"runnum",             {10000, 4000, 6000, true}},
        {"Delta_eta",          {100, -6.0, 1.0, true}},
        {"mc_Delta_eta",       {100, -6.0, 1.0, true}},
        {"Delta_phi",          {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_Delta_phi",       {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"Delta_phi12",        {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_Delta_phi12",     {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"Delta_phi13",        {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_Delta_phi13",     {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"Delta_phi23",        {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_Delta_phi23",     {100, 0.0, 2.0 * TMath::Pi(), true}},

        {"DepA",               {100, 0.0, 1.0, true}},
        {"mc_DepA",            {100, 0.0, 1.0, true}},
        {"DepB",               {100, 0.0, 1.0, true}},
        {"mc_DepB",            {100, 0.0, 1.0, true}},
        {"DepC",               {100, 0.0, 1.0, true}},
        {"mc_DepC",            {100, 0.0, 1.0, true}},
        {"DepV",               {100, 0.0, 2.0, true}},
        {"mc_DepV",            {100, 0.0, 2.0, true}},
        {"DepW",               {100, 0.0, 1.0, true}},
        {"mc_DepW",            {100, 0.0, 1.0, true}},

        {"e_p",                {100, 1.0, 9.0, true}},
        {"mc_e_p",             {100, 1.0, 9.0, true}},
        {"e_phi",              {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_e_phi",           {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"e_theta",            {100, 0.0, 90.0, true}},
        {"mc_e_theta",         {100, 0.0, 90.0, true}},

        {"eta",                {100, -3.0, 3.0, true}},
        {"mc_eta",             {100, -3.0, 3.0, true}},
        {"eta1",               {100, -3.0, 3.0, true}},
        {"mc_eta1",            {100, -3.0, 3.0, true}},
        {"eta2",               {100, -3.0, 3.0, true}},
        {"mc_eta2",            {100, -3.0, 3.0, true}},
        {"eta3",               {100, -3.0, 3.0, true}},
        {"mc_eta3",            {100, -3.0, 3.0, true}},
        {"eta12",              {100, -3.0, 3.0, true}},
        {"mc_eta12",           {100, -3.0, 3.0, true}},
        {"eta13",              {100, -3.0, 3.0, true}},
        {"mc_eta13",           {100, -3.0, 3.0, true}},
        {"eta23",              {100, -3.0, 3.0, true}},
        {"mc_eta23",           {100, -3.0, 3.0, true}},
        {"eta1_gN",            {100, -3.0, 3.0, true}},
        {"mc_eta1_gN",         {100, -3.0, 3.0, true}},
        {"eta2_gN",            {100, -3.0, 3.0, true}},
        {"mc_eta2_gN",         {100, -3.0, 3.0, true}},

        {"evnum",              {100, 0.0, 0.0, false}},
        {"mc_evnum",           {100, 0.0, 0.0, false}},
        {"helicity",           {2, -2.0, 2.0, true}},

        {"Mh",                 {100, 0.0, 4.0, true}},
        {"mc_Mh",              {100, 0.0, 4.0, true}},
        {"Mh12",               {100, 0.0, 2.5, true}},
        {"mc_Mh12",            {100, 0.0, 2.5, true}},
        {"Mh13",               {100, 0.0, 2.5, true}},
        {"mc_Mh13",            {100, 0.0, 2.5, true}},
        {"Mh23",               {100, 0.0, 2.5, true}},
        {"mc_Mh23",            {100, 0.0, 2.5, true}},

        {"Mx",                 {100, -2.0, 2.0, true}},
        {"mc_Mx",              {100, -2.0, 2.0, true}},
        {"Mx1",                {100, -2.0, 2.0, true}},
        {"mc_Mx1",             {100, -2.0, 2.0, true}},
        {"Mx2",                {100, 0.0, 2.0, true}},
        {"mc_Mx2",             {100, 0.0, 2.0, true}},
        {"Mx3",                {100, 0.0, 4.0, true}},
        {"mc_Mx3",             {100, 0.0, 4.0, true}},
        {"Mx12",               {100, 0.0, 4.0, true}},
        {"mc_Mx12",            {100, 0.0, 4.0, true}},
        {"Mx13",               {100, 0.0, 4.0, true}},
        {"mc_Mx13",            {100, 0.0, 4.0, true}},
        {"Mx23",               {100, 0.0, 4.0, true}},
        {"mc_Mx23",            {100, 0.0, 4.0, true}},

        {"phi",                {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_phi",             {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"phi1",               {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_phi1",            {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"phi2",               {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_phi2",            {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"phi3",               {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_phi3",            {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"phih",               {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_phih",            {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"phiR",               {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_phiR",            {100, 0.0, 2.0 * TMath::Pi(), true}},

        {"p_p",                {100, 0.0, 6.0, true}},
        {"mc_p_p",             {100, 0.0, 6.0, true}},
        {"p1_p",               {100, 0.0, 3.0, true}},
        {"mc_p1_p",            {100, 0.0, 3.0, true}},
        {"p2_p",               {100, 1.0, 9.0, true}},
        {"mc_p2_p",            {100, 1.0, 9.0, true}},
        {"p3_p",               {100, 0.0, 6.0, true}},
        {"mc_p3_p",            {100, 0.0, 6.0, true}},

        {"p_phi",              {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_p_phi",           {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"p1_phi",             {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_p1_phi",          {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"p2_phi",             {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_p2_phi",          {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"p3_phi",             {100, 0.0, 2.0 * TMath::Pi(), true}},
        {"mc_p3_phi",          {100, 0.0, 2.0 * TMath::Pi(), true}},

        {"pT",                 {100, 0.0, 0.3, true}},
        {"mc_pT",              {100, 0.0, 0.3, true}},
        {"pT1",                {100, 0.0, 1.2, true}},
        {"mc_pT1",             {100, 0.0, 1.2, true}},
        {"pT2",                {100, 0.0, 1.2, true}},
        {"mc_pT2",             {100, 0.0, 1.2, true}},
        {"pTpT",               {100, 0.0, 1.0, true}},
        {"mc_pTpT",            {100, 0.0, 1.0, true}},

        {"p_theta",            {100, 0.0, 90.0, true}},
        {"mc_p_theta",         {100, 0.0, 90.0, true}},
        {"p1_theta",           {100, 0.0, 90.0, true}},
        {"mc_p1_theta",        {100, 0.0, 90.0, true}},
        {"p2_theta",           {100, 6.0, 45.0, true}},
        {"mc_p2_theta",        {100, 0.0, 90.0, true}},
        {"p3_theta",           {100, 0.0, 90.0, true}},
        {"mc_p3_theta",        {100, 0.0, 90.0, true}},

        {"Q2",                 {100, 0.0, 9.0, true}},
        {"mc_Q2",              {100, 0.0, 9.0, true}},

        {"theta",              {100, 0.0, TMath::Pi(), true}},
        {"mc_theta",           {100, 0.0, TMath::Pi(), true}},

        {"t",                  {100, -15.0, 1.0, true}},
        {"mc_t",               {100, -15.0, 1.0, true}},
        {"tprime",                  {100, -15.0, 1.0, true}},
        {"mc_tprime",               {100, -15.0, 1.0, true}},
        {"t1",                 {100, -1.0, 1.0, true}},
        {"mc_t1",              {100, -1.0, 1.0, true}},
        {"t2",                 {100, -15.0, 1.0, true}},
        {"mc_t2",              {100, -15.0, 1.0, true}},
        {"tmin",               {100, -0.5, 0.0, true}},
        {"mc_tmin",            {100, -0.5, 0.0, true}},

        {"vz_e",               {100, -15.0, 15.0, true}},
        {"mc_vz_e",            {100, -15.0, 15.0, true}},
        {"vz_p",               {100, -15.0, 15.0, true}},
        {"mc_vz_p",            {100, -15.0, 15.0, true}},
        {"vz_p1",              {100, -15.0, 15.0, true}},
        {"mc_vz_p1",           {100, -15.0, 15.0, true}},
        {"vz_p2",              {100, -15.0, 15.0, true}},
        {"mc_vz_p2",           {100, -15.0, 15.0, true}},
        {"vz_p3",              {100, -15.0, 15.0, true}},
        {"mc_vz_p3",           {100, -15.0, 15.0, true}},

        {"W",                  {100, 2.0, 4.0, true}},
        {"mc_W",               {100, 2.0, 4.0, true}},
        {"x",                  {100, 0.0, 0.8, true}},
        {"mc_x",               {100, 0.0, 0.8, true}},

        {"xF",                 {100, -1.0, 1.0, true}},
        {"mc_xF",              {100, -1.0, 1.0, true}},
        {"xF1",                {100, -1.0, 1.0, true}},
        {"mc_xF1",             {100, -1.0, 1.0, true}},
        {"xF2",                {100, -1.0, 1.0, true}},
        {"mc_xF2",             {100, -1.0, 1.0, true}},
        {"xF3",                {100, -1.0, 1.0, true}},
        {"mc_xF3",             {100, -1.0, 1.0, true}},
        {"xF12",               {100, -1.0, 1.0, true}},
        {"mc_xF12",            {100, -1.0, 1.0, true}},
        {"xF13",               {100, -1.0, 1.0, true}},
        {"mc_xF13",            {100, -1.0, 1.0, true}},
        {"xF23",               {100, -1.0, 1.0, true}},
        {"mc_xF23",            {100, -1.0, 1.0, true}},

        {"y",                  {100, 0.0, 1.0, true}},
        {"mc_y",               {100, 0.0, 1.0, true}},
        {"z",                  {100, 0.0, 2.0, true}},
        {"mc_z",               {100, 0.0, 2.0, true}},
        {"z1",                 {100, 0.0, 2.0, true}},
        {"mc_z1",              {100, 0.0, 2.0, true}},
        {"z2",                 {100, 0.0, 2.0, true}},
        {"mc_z2",              {100, 0.0, 2.0, true}},

        {"zeta",               {100, 0.0, 2.0, true}},
        {"mc_zeta",            {100, 0.0, 2.0, true}},
        {"zeta1",              {100, 0.0, 2.0, true}},
        {"mc_zeta1",           {100, 0.0, 2.0, true}},
        {"zeta2",              {100, 0.0, 2.0, true}},
        {"mc_zeta2",           {100, 0.0, 2.0, true}},

        {"Emiss2",             {100, 0.0, 2.0, true}},
        {"theta_gamma_gamma",  {100, 0.0, 2.0, true}},
        {"pTmiss",             {100, 0.0, 0.5, true}},
        {"Mxgammasquared",     {100, 0.0, 2.0, true}}
    };
    return cfg;
}

bool hasBranchPlotConfig(const std::string& branch_name) {
    const auto& cfg = getBranchPlotConfigs();
    return cfg.find(branch_name) != cfg.end();
}

const BranchHistConfig* findBranchPlotConfig(const std::string& branch_name) {
    const auto& cfg = getBranchPlotConfigs();
    auto it = cfg.find(branch_name);
    if (it == cfg.end()) return nullptr;
    return &(it->second);
}