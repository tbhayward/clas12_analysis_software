// One streaming pass per reconstructed DATA/MC file; no event arrays in Python.
#include "global_cuts.h"
#include "exclusivity_resolution_core.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
using json=nlohmann::json;
using namespace resolution_study;
static json read_json(const std::string& path) {
    std::ifstream f(path); if(!f) throw std::runtime_error("Cannot read "+path);
    json j; f>>j; return j;
}
static GlobalCutConfig config(const json& j) {
    GlobalCutConfig c;
#define LOAD(name) c.name = j.at(#name).get<decltype(c.name)>()
    LOAD(t1_abs_max); LOAD(open_angle_min_deg);
    LOAD(enable_reconstructed_electron_p_cut); LOAD(reconstructed_electron_p_min_GeV);
    LOAD(enable_analysis_phase_space_cuts); LOAD(enable_analysis_y_min_cut);
    LOAD(analysis_y_min); LOAD(analysis_W_min_GeV); LOAD(enable_pTmiss_cut); LOAD(pTmiss_max);
    LOAD(enable_dvcsgen_ycol_cut); LOAD(dvcsgen_ycol_cut);
    LOAD(enable_auxiliary_fiducial_cuts); LOAD(auxiliary_require_distinct_fd_sectors);
    LOAD(auxiliary_e_theta_min_deg); LOAD(auxiliary_e_theta_max_deg);
    LOAD(auxiliary_fd_proton_theta_min_deg); LOAD(auxiliary_fd_proton_theta_max_deg);
    LOAD(auxiliary_fd_proton_p_min_inb_GeV); LOAD(auxiliary_fd_proton_p_min_out_GeV);
    LOAD(auxiliary_fd_photon_theta_min_deg); LOAD(auxiliary_fd_photon_theta_max_deg);
    LOAD(auxiliary_cd_proton_p_min_GeV); LOAD(auxiliary_cd_proton_theta_min_deg);
    LOAD(auxiliary_cd_proton_theta_max_deg); LOAD(auxiliary_ft_photon_p_min_GeV);
    LOAD(enable_topology_filter); LOAD(required_detector1); LOAD(required_detector2);
    LOAD(enable_electron_fd_sector_filter); LOAD(electron_fd_sector);
    LOAD(enable_proton_fd_sector_filter); LOAD(proton_fd_sector);
    LOAD(enable_proton_cd_sector_filter); LOAD(proton_cd_sector);
    LOAD(enable_photon_fd_sector_filter); LOAD(photon_fd_sector); LOAD(excluded_runs);
#undef LOAD
    return c;
}
struct Reader {
    struct Bound {TBranch* branch;TLeaf* conversion;double* value;};
    TTree* tree;
    std::map<std::string,double> v;
    std::vector<Bound> generated, reconstructed;
    Reader(TTree* t):tree(t) { tree->SetBranchStatus("*",0); }
    bool bind(const std::string& name, std::vector<Bound>& group, bool optional=false) {
        std::string physical=name;
        auto* b=tree->GetBranch(physical.c_str());
        if (!b) { if (optional) return false; throw std::runtime_error("Missing branch " + name); }
        auto* leaves=b->GetListOfLeaves();
        auto* leaf=leaves && leaves->GetEntriesFast()==1?dynamic_cast<TLeaf*>(leaves->At(0)):nullptr;
        const std::set<std::string> numeric={"Double_t","Float_t","Int_t","UInt_t","Long_t","ULong_t","Long64_t","ULong64_t","Short_t","UShort_t","Char_t","UChar_t","Bool_t"};
        if(!leaf || leaf->GetLeafCount() || leaf->GetLenStatic()!=1 || !numeric.count(leaf->GetTypeName()))
            throw std::runtime_error("Expected scalar numeric branch " + physical);
        auto& value=v[name]; tree->SetBranchStatus(physical.c_str(),1);
        bool direct=std::string(leaf->GetTypeName())=="Double_t";
        if (direct && tree->SetBranchAddress(physical.c_str(),&value)<0)
            throw std::runtime_error("Incompatible branch type for " + name);
        tree->AddBranchToCache(physical.c_str(),true); group.push_back({b,direct?nullptr:leaf,&value}); return true;
    }
    double get(const std::string& name) const { return v.at(name); }
    void read(const std::vector<Bound>& group, Long64_t i) {
        for (const auto& b:group) {
            if (b.branch->GetEntry(i)<=0) throw std::runtime_error("Failed branch read");
            if(b.conversion) *b.value=b.conversion->GetValue();
        }
    }
};
static Window window(const json& j) {
    Window w;w.active=true;w.upper=j.value("mode",std::string())=="upper_quantile";
    w.high=j.value("cut_high",std::numeric_limits<double>::quiet_NaN());
    w.low=j.value("cut_low",std::numeric_limits<double>::quiet_NaN());
    if(!w.upper && (!std::isfinite(w.low)||!std::isfinite(w.high)||w.high<=w.low)) {
        double mean=j.at("mean"),sd=j.at("std");
        if(!std::isfinite(mean)||!std::isfinite(sd)||sd<=0) throw std::runtime_error("Invalid fallback cut");
        w.low=mean-3*sd;w.high=mean+3*sd;
    }
    if(!std::isfinite(w.high)) throw std::runtime_error("Invalid saved cut");
    return w;
}
static void run(const json& spec,const json& job,const std::string& output) {
    auto cfg=config(spec.at("global_cuts"));
    auto no_pt=cfg;no_pt.enable_pTmiss_cut=false;
    std::string period=job.at("period"),sample=job.at("sample");
    std::string code=spec.at("period_codes").at(period);
    std::array<std::array<Window,NV>,3> cuts;
    std::array<std::string,3> names={"FD_FD","CD_FD","CD_FT"};
    for(int t=0;t<3;++t) {
        const auto& block=spec.at("exclusivity").at("DVCS_"+code+"_"+names[t]).at(sample);
        if(!block.contains("Mx2")) throw std::runtime_error("Missing mandatory Mx2 cut");
        for(auto it=block.begin();it!=block.end();++it) {
            auto found=std::find(VARIABLES.begin(),VARIABLES.end(),it.key());
            if(found==VARIABLES.end()) throw std::runtime_error("Unsupported saved cut "+it.key());
            cuts[t][found-VARIABLES.begin()]=window(it.value());
        }
    }
    std::vector<Stratum> strata;
    for(const auto& s:spec.at("strata")) strata.push_back({s.at("axis"),s.value("low",0.),s.value("high",0.)});
    const int ng=strata.size();
    auto axis_index=spec.at("strata_axis_index").get<std::vector<int>>();
    std::vector<int> groups;groups.reserve(7);
    std::vector<Histogram> hist;
    for(int t=0;t<3;++t) for(int g=0;g<ng;++g) for(const auto& v:spec.at("variables"))
        hist.emplace_back(v.at("bins"),v.at("low"),v.at("high"));
    std::array<uint64_t,3> global_count{},selected{};
    uint64_t entries=0;
    json files=json::array();
    for(const auto& input:job.at("files")) {
        std::string path=input;TFile f(path.c_str(),"READ");
        if(f.IsZombie()) throw std::runtime_error("Cannot open "+path);
        auto* tree=dynamic_cast<TTree*>(f.Get("PhysicsEvents"));
        if(!tree) throw std::runtime_error("Missing PhysicsEvents in "+path);
        tree->SetCacheSize(32*1024*1024);Reader b(tree);
        for(auto name:{"x","Q2","t1","detector1","detector2","e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi","open_angle_ep2","pTmiss"}) b.bind(name,b.reconstructed);
        bool runnum=b.bind("runnum",b.reconstructed,true);
        if(!runnum && !cfg.excluded_runs.empty()) throw std::runtime_error("Cannot enforce run blacklist without runnum");
        bool direct_delta=b.bind("Delta_phi",b.reconstructed,true);
        for(const auto& name:VARIABLES) if(name!="Delta_phi" && name!="pTmiss") b.bind(name,b.reconstructed);
        tree->StopCacheLearningPhase();
        // Store stable pointers outside the event loop; no per-event map lookups.
        std::vector<const double*> p;
        for(auto n:{"x","Q2","t1","detector1","detector2","e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi","open_angle_ep2","pTmiss"}) p.push_back(&b.v.at(n));
        std::array<const double*,NV> vp;
        for(size_t v=0;v<NV;++v) vp[v]=(v==0&&!direct_delta)?nullptr:&b.v.at(VARIABLES[v]);
        const double* rp=runnum?&b.v.at("runnum"):nullptr;
        auto global=[&](const GlobalCutConfig& c) {
            return passes_analysis_phase_space_cuts(period,*p[0],*p[1],*p[2],*p[5],c) &&
                passes_proton_momentum_cut(int(*p[3]),period,*p[8],c) &&
                passes_global_cuts(*p[2],*p[14],*p[15],int(*p[3]),int(*p[4]),period,
                    *p[5],*p[6],*p[7],*p[9],*p[10],*p[11],*p[12],*p[13],c);
        };
        for(Long64_t i=0;i<tree->GetEntries();++i) {
            ++entries;b.read(b.reconstructed,i);
            if(!std::isfinite(*p[3])||!std::isfinite(*p[4])) continue;
            int d1=int(*p[3]),d2=int(*p[4]);
            int t=d1==1&&d2==1?0:d1==2&&d2==1?1:d1==2&&d2==0?2:-1;
            if(t<0) continue;
            if(rp && is_excluded_run(int(*rp),cfg)) continue;
            bool gp=global(cfg);
            bool relaxed=gp || (cfg.enable_pTmiss_cut && global(no_pt));
            if(!relaxed) continue;
            if(gp) ++global_count[t];
            std::array<double,NV> values;
            for(size_t v=0;v<NV;++v) values[v]=vp[v]?*vp[v]:std::abs(std::remainder(*p[10]-*p[13],2*3.14159265358979323846));
            unsigned mask=failures(values,cuts[t]);
            if(gp && !mask) ++selected[t];
            std::array<double,6> axes={*p[5],*p[6]*DEG,*p[8],*p[9]*DEG,*p[11],*p[12]*DEG};
            groups.clear();groups.push_back(0);
            // Strata are ordered by six axes in the specification.
            for(int g=1;g<ng;++g) {
                int a=axis_index[g];
                if(strata[g].contains(axes[a])) groups.push_back(g);
            }
            for(size_t v=0;v<NV;++v) {
                if(!leave_one_out(mask,v) || (!gp && v!=3)) continue;
                bool own_pass=cuts[t][v].pass(values[v]);
                if(v==3 && cfg.enable_pTmiss_cut) own_pass=own_pass && values[v]<=cfg.pTmiss_max;
                for(int g:groups) hist[(t*ng+g)*NV+v].add(values[v],own_pass);
            }
        }
        files.push_back({{"path",path},{"entries",tree->GetEntries()},{"delta_phi_fallback",!direct_delta},{"has_runnum",runnum}});
    }
    json out={{"entries",entries},{"global_pass",global_count},{"selected",selected},{"files",files},{"histograms",json::array()}};
    for(const auto& h:hist) out["histograms"].push_back({{"counts",h.counts},{"nonfinite",h.nonfinite},{"pass_own_cut",h.pass_own_cut}});
    std::ofstream f(output);if(!f) throw std::runtime_error("Cannot write "+output);f<<out;
}
int main(int argc,char** argv) {
    try {if(argc!=4) throw std::runtime_error("usage: worker SPEC JOB OUTPUT");run(read_json(argv[1]),read_json(argv[2]),argv[3]);}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
