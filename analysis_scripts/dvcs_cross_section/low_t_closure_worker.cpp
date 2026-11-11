// Standalone paired-MC closure accumulator. Global selections call production code.
#include "global_cuts.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using json = nlohmann::json;
constexpr double DEG = 180.0 / 3.14159265358979323846;
static json read_json(const std::string& path) {
    std::ifstream f(path); if (!f) throw std::runtime_error("Cannot read " + path);
    json j; f >> j; return j;
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
static uint64_t mix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
static uint64_t path_hash(const std::string& s) {
    uint64_t h=1469598103934665603ULL;
    for (unsigned char c:s) { h^=c; h*=1099511628211ULL; }
    return h;
}
struct Reader {
    struct Bound {TBranch* branch;TLeaf* conversion;double* value;};
    TTree* tree;
    std::map<std::string,double> v;
    std::vector<Bound> generated, reconstructed;
    Reader(TTree* t):tree(t) { tree->SetBranchStatus("*",0); }
    bool bind(const std::string& name, std::vector<Bound>& group, bool optional=false) {
        std::string physical=name;
        if(!tree->GetBranch(physical.c_str()) &&
           (name=="rec_runnum" || name=="rec_detector1" || name=="rec_detector2")) physical=name.substr(4);
        auto* b=tree->GetBranch(physical.c_str());
        if (!b) { if (optional) return false; throw std::runtime_error("Missing paired branch " + name); }
        auto* leaves=b->GetListOfLeaves();
        auto* leaf=leaves && leaves->GetEntriesFast()==1?dynamic_cast<TLeaf*>(leaves->At(0)):nullptr;
        const std::set<std::string> numeric={"Double_t","Float_t","Int_t","UInt_t","Long_t","ULong_t","Long64_t","ULong64_t","Short_t","UShort_t","Char_t","UChar_t","Bool_t"};
        if(!leaf || leaf->GetLeafCount() || leaf->GetLenStatic()!=1 || !numeric.count(leaf->GetTypeName()))
            throw std::runtime_error("Expected scalar numeric paired branch " + physical);
        auto& value=v[name]; tree->SetBranchStatus(physical.c_str(),1);
        bool direct=std::string(leaf->GetTypeName())=="Double_t";
        if (direct && tree->SetBranchAddress(physical.c_str(),&value)<0)
            throw std::runtime_error("Incompatible branch type for " + name);
        tree->AddBranchToCache(physical.c_str(),true); group.push_back({b,direct?nullptr:leaf,&value}); return true;
    }
    double get(const std::string& name) const { return v.at(name); }
    void read(const std::vector<Bound>& group, Long64_t i) {
        for (const auto& b:group) {
            if (b.branch->GetEntry(i)<=0) throw std::runtime_error("Failed paired branch read");
            if(b.conversion) *b.value=b.conversion->GetValue();
        }
    }
};
struct Grid {
    std::vector<double> x,q,t,p;
    std::vector<int> lut;
    int n;
    explicit Grid(const json& j):x(j.at("edges")[0].get<std::vector<double>>()),
        q(j.at("edges")[1].get<std::vector<double>>()),t(j.at("edges")[2].get<std::vector<double>>()),
        p(j.at("edges")[3].get<std::vector<double>>()),lut(j.at("lookup").get<std::vector<int>>()),
        n(j.at("bins").size()) {}
    static int index(const std::vector<double>& e,double v) {
        if (!std::isfinite(v) || v<e.front() || v>=e.back()) return -1;
        return int(std::upper_bound(e.begin(),e.end(),v)-e.begin())-1;
    }
    int find(double xv,double qv,double tv,double pv) const {
        pv=std::fmod(pv*DEG,360.0); if (pv<0) pv+=360.0;
        int a=index(x,xv),b=index(q,qv),c=index(t,tv),d=index(p,pv);
        if (a<0 || b<0 || c<0 || d<0) return -1;
        return lut[(((a*int(q.size()-1)+b)*int(t.size()-1)+c)*int(p.size()-1)+d)];
    }
};
struct Scheme {
    Grid grid;
    int ng,ns,n;
    std::vector<double> train_gen,train_reco,train_reco2,train_diag,test_gen,test_gen2,test_reco,test_reco2,test_diag;
    Scheme(const json& j,int groups,int scenarios):grid(j),ng(groups),ns(scenarios),n(grid.n) {
        for(auto* v:{&train_gen,&train_reco,&train_reco2,&train_diag}) v->resize(ng*n);
        for(auto* v:{&test_gen,&test_gen2,&test_reco,&test_reco2,&test_diag}) v->resize(ns*ng*n);
    }
    void add(int g,int r,bool train,int angle,double response,const std::vector<double>& w) {
        for (int group:{0,angle}) {
            if (group<0) continue;
            int off=group*n;
            if (train) {
                if (g>=0) train_gen[off+g]+=1.;
                if (r>=0) { train_reco[off+r]+=response; train_reco2[off+r]+=response*response; }
                if (g>=0 && g==r) train_diag[off+g]+=response;
            } else for (int s=0;s<ns;++s) {
                int k=(s*ng+group)*n;
                if (g>=0) { test_gen[k+g]+=w[s]; test_gen2[k+g]+=w[s]*w[s]; }
                if (r>=0) { test_reco[k+r]+=response*w[s]; test_reco2[k+r]+=response*response*w[s]*w[s]; }
                if (g>=0 && g==r) test_diag[k+g]+=response*w[s]*w[s];
            }
        }
    }
};
struct Window {
    bool upper=false,skip=false;
    double low=0,high=0;
    explicit Window(const json& j) {
        upper=j.value("mode",std::string("symmetric_3sigma"))=="upper_quantile";
        high=j.value("cut_high",std::numeric_limits<double>::quiet_NaN());
        if(upper) { skip=!std::isfinite(high); return; }
        low=j.value("cut_low",std::numeric_limits<double>::quiet_NaN());
        if(!std::isfinite(low) || !std::isfinite(high) || high<=low) {
            double mean=j.value("mean",std::numeric_limits<double>::quiet_NaN());
            double sd=j.value("std",std::numeric_limits<double>::quiet_NaN());
            skip=!std::isfinite(mean) || !std::isfinite(sd) || sd<=0;
            low=mean-3*sd;high=mean+3*sd;
        }
    }
    bool pass(double value) const {
        return std::isfinite(value) && (skip || (upper?value<=high:value>=low && value<=high));
    }
};
struct PreparedCut { const double* value; bool delta; Window window; };
static int category(double t) {
    static const std::vector<double> e={0,.05,.08,.11,.15,.25,.4,.6,.8,1.};
    return int(std::upper_bound(e.begin(),e.end(),t)-e.begin());
}
static void weights(double t,const std::vector<std::array<double,2>>& scenarios,std::vector<double>& w) {
    const double bump=std::exp(-.5*std::pow((t-.13)/.025,2));
    for(size_t i=0;i<scenarios.size();++i) {
        const auto& s=scenarios[i];
        double value=s[0]==0?1.:std::exp(std::clamp(s[0]*(t-.2),-std::log(5.),std::log(5.)));
        value*=1+s[1]*bump;w[i]=value;
    }
}
static void run(const json& spec,const json& job,const std::string& destination) {
    auto cfg=config(spec.at("global_cuts"));
    std::string period=job.at("period"),tag=spec.at("period_tags").at(period);
    std::string code=spec.at("period_codes").at(period);
    const auto& cuts=spec.at("exclusivity");
    const auto& model=spec.at("current_model").at("mc").at("ep->epg").at(tag);
    std::array<double,7> response_weights;
    for(int i=0;i<7;++i) {
        std::string region=i==0?"FT":"S"+std::to_string(i);
        double factor=model.at(region).at("reference_factor");
        if(!std::isfinite(factor) || factor<=0) throw std::runtime_error("Invalid MC response for "+region);
        response_weights[i]=1./factor;
    }
    std::vector<std::array<double,2>> stress;
    for(const auto& scenario:spec.at("scenarios")) stress.push_back({scenario.at("slope"),scenario.at("bump")});
    std::vector<double> w(stress.size());
    std::array<std::string,3> names={"FD_FD","CD_FD","CD_FT"};
    std::set<std::string> variables;
    const std::set<std::string> allowed={"Delta_phi","theta","theta_gamma_gamma","pTmiss","Emiss2","Mx2","Mx2_2"};
    for (const auto& name:names) for (auto it=cuts.at("DVCS_"+code+"_"+name).at("mc").begin();it!=cuts.at("DVCS_"+code+"_"+name).at("mc").end();++it) {
        if (!allowed.count(it.key())) throw std::runtime_error("Unknown production exclusivity variable "+it.key());
        variables.insert(it.key());
    }
    auto angle_edges=spec.at("angle_edges").get<std::vector<double>>();
    int ng=angle_edges.size(),ns=spec.at("scenarios").size();
    std::vector<Scheme> schemes;
    for (const auto& grid:spec.at("schemes")) schemes.emplace_back(grid,ng,ns);
    std::vector<double> migration(12*11),topology(3*11*11),angles(8*11*11);
    struct Flow {uint64_t entries=0,valid_truth=0,generated_channel=0,valid_reco=0,
        global_pass=0,selected=0,reconstructed_without_valid_truth=0,nonfinite_truth_t=0;} flow;
    for (const auto& input:job.at("files")) {
        std::string path=input.get<std::string>(); TFile f(path.c_str(),"READ");
        if (f.IsZombie()) throw std::runtime_error("Cannot open "+path);
        auto* tree=dynamic_cast<TTree*>(f.Get("MCEvents"));
        if (!tree) throw std::runtime_error("Missing MCEvents in "+path);
        tree->SetCacheSize(32*1024*1024);
        Reader b(tree); std::vector<Reader::Bound> valid;
        b.bind("gen_kinematics_valid",valid);
        b.bind("rec_kinematics_valid",valid);
        for (auto name:{"x","Q2","t1","phi2","p1_theta"}) b.bind("gen_"+std::string(name),b.generated);
        for (auto name:{"x","Q2","t1","phi2","detector1","detector2","e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi","open_angle_ep2"}) b.bind("rec_"+std::string(name),b.reconstructed);
        bool has_run=b.bind("rec_runnum",b.reconstructed,true);
        bool delta_direct=true;
        if (cfg.enable_pTmiss_cut || variables.count("pTmiss")) variables.insert("pTmiss");
        for (const auto& v:variables) if (!b.bind("rec_"+v,b.reconstructed,v=="Delta_phi")) delta_direct=false;
        tree->StopCacheLearningPhase();
        std::array<std::vector<PreparedCut>,3> prepared;
        for(int i=0;i<3;++i) {
            const auto& cut=cuts.at("DVCS_"+code+"_"+names[i]).at("mc");
            for(auto it=cut.begin();it!=cut.end();++it) {
                bool delta=it.key()=="Delta_phi" && !delta_direct;
                prepared[i].push_back({delta?nullptr:&b.v.at("rec_"+it.key()),delta,Window(it.value())});
            }
        }
        const double zero=0.;
        const double& v_gen_Q2=b.v.at("gen_Q2");
        const double& v_gen_kinematics_valid=b.v.at("gen_kinematics_valid");
        const double& v_gen_p1_theta=b.v.at("gen_p1_theta");
        const double& v_gen_phi2=b.v.at("gen_phi2");
        const double& v_gen_t1=b.v.at("gen_t1");
        const double& v_gen_x=b.v.at("gen_x");
        const double& v_rec_Q2=b.v.at("rec_Q2");
        const double& v_rec_detector1=b.v.at("rec_detector1");
        const double& v_rec_detector2=b.v.at("rec_detector2");
        const double& v_rec_e_p=b.v.at("rec_e_p");
        const double& v_rec_e_phi=b.v.at("rec_e_phi");
        const double& v_rec_e_theta=b.v.at("rec_e_theta");
        const double& v_rec_kinematics_valid=b.v.at("rec_kinematics_valid");
        const double& v_rec_open_angle_ep2=b.v.at("rec_open_angle_ep2");
        const double& v_rec_p1_p=b.v.at("rec_p1_p");
        const double& v_rec_p1_phi=b.v.at("rec_p1_phi");
        const double& v_rec_p1_theta=b.v.at("rec_p1_theta");
        const double& v_rec_p2_p=b.v.at("rec_p2_p");
        const double& v_rec_p2_phi=b.v.at("rec_p2_phi");
        const double& v_rec_p2_theta=b.v.at("rec_p2_theta");
        const double& v_rec_pTmiss=b.v.count("rec_pTmiss")?b.v.at("rec_pTmiss"):zero;
        const double& v_rec_phi2=b.v.at("rec_phi2");
        const double& v_rec_runnum=has_run?b.v.at("rec_runnum"):zero;
        const double& v_rec_t1=b.v.at("rec_t1");
        const double& v_rec_x=b.v.at("rec_x");

        uint64_t seed=path_hash(path)^spec.at("seed").get<uint64_t>();
        for (Long64_t i=0;i<tree->GetEntries();++i) {
            ++flow.entries; b.read(valid,i);
            bool gv=v_gen_kinematics_valid>=.5,rv=v_rec_kinematics_valid>=.5;
            if (!gv) { if(rv) ++flow.reconstructed_without_valid_truth; continue; }
            b.read(b.generated,i);
            double gt=-v_gen_t1; if (!std::isfinite(gt)) { ++flow.nonfinite_truth_t; continue; }
            ++flow.valid_truth;
            bool generated=passes_generated_channel_cuts(v_gen_x,v_gen_Q2,cfg);
            if(generated) ++flow.generated_channel;
            bool selected=false; int top=-1; double rt=0,response=1.;
            if (rv) {
                ++flow.valid_reco; b.read(b.reconstructed,i);
                int d1=int(v_rec_detector1),d2=int(v_rec_detector2);
                top=d1==1&&d2==1?0:d1==2&&d2==1?1:d1==2&&d2==0?2:-1;
                double ptmiss=v_rec_pTmiss;
                bool blacklist=has_run && std::find(cfg.excluded_runs.begin(),cfg.excluded_runs.end(),int(v_rec_runnum))!=cfg.excluded_runs.end();
                bool global=top>=0 && !blacklist &&
                    passes_analysis_phase_space_cuts(period,v_rec_x,v_rec_Q2,v_rec_t1,v_rec_e_p,cfg) &&
                    passes_proton_momentum_cut(d1,period,v_rec_p1_p,cfg) &&
                    passes_global_cuts(v_rec_t1,v_rec_open_angle_ep2,ptmiss,d1,d2,period,
                        v_rec_e_p,v_rec_e_theta,v_rec_e_phi,v_rec_p1_theta,v_rec_p1_phi,
                        v_rec_p2_p,v_rec_p2_theta,v_rec_p2_phi,cfg);
                if(global) {
                    ++flow.global_pass; selected=true;
                    for(const auto& cut:prepared[top]) {
                        double v=cut.delta?std::abs(std::remainder(v_rec_p1_phi-v_rec_p2_phi,2*std::acos(-1.))):*cut.value;
                        if(!cut.window.pass(v)) {selected=false;break;}
                    }
                    if(selected) {
                        ++flow.selected;rt=-v_rec_t1;
                        int region=d2==0?0:fd_sector_from_phi_rad(v_rec_p2_phi);
                        if(region<0 || region>6) throw std::runtime_error("Invalid photon region");
                        response=response_weights[region];
                    }
                }
            }
            if(generated) {
                int gc=category(gt),rc=selected?category(rt):11; migration[rc*11+gc]+=1;
                if(selected) { topology[(top*11+rc)*11+gc]+=1; }
            }
            double theta=v_gen_p1_theta*DEG;
            int ai=Grid::index(angle_edges,theta); int angle=ai>=0?ai+1:-1;
            if(generated && selected && angle>=0) angles[(angle*11+category(rt))*11+category(gt)]+=1;
            bool train=(mix(uint64_t(i)^seed)&1)==0;
            bool weights_ready=false;
            for(auto& scheme:schemes) {
                int g=generated?scheme.grid.find(v_gen_x,v_gen_Q2,gt,v_gen_phi2):-1;
                int r=selected?scheme.grid.find(v_rec_x,v_rec_Q2,rt,v_rec_phi2):-1;
                if(g<0 && r<0) continue;
                if(!train && !weights_ready) {weights(gt,stress,w);weights_ready=true;}
                scheme.add(g,r,train,angle,response,w);
            }
        }
    }
    json output;output["flow"]["entries"] = flow.entries;output["flow"]["generated_channel"] = flow.generated_channel;output["flow"]["global_pass"] = flow.global_pass;output["flow"]["nonfinite_truth_t"] = flow.nonfinite_truth_t;output["flow"]["reconstructed_without_valid_truth"] = flow.reconstructed_without_valid_truth;output["flow"]["selected"] = flow.selected;output["flow"]["valid_reco"] = flow.valid_reco;output["flow"]["valid_truth"] = flow.valid_truth;output["migration"]=migration;output["topology"]=topology;output["angles"]=angles;
    for(const auto& scheme:schemes) {
        json a;
        a["train_gen"] = scheme.train_gen;
        a["train_reco"] = scheme.train_reco;
        a["train_reco2"] = scheme.train_reco2;
        a["train_diag"] = scheme.train_diag;
        a["test_gen"] = scheme.test_gen;
        a["test_gen2"] = scheme.test_gen2;
        a["test_reco"] = scheme.test_reco;
        a["test_reco2"] = scheme.test_reco2;
        a["test_diag"] = scheme.test_diag;
        output["schemes"].push_back(std::move(a));
    }
    std::ofstream out(destination);if(!out) throw std::runtime_error("Cannot write "+destination);out<<output;
}
int main(int argc,char** argv) {
    try {
        if(argc!=4) throw std::runtime_error("Usage: low_t_closure_worker spec.json job.json output.json");
        run(read_json(argv[1]),read_json(argv[2]),argv[3]); return 0;
    } catch(const std::exception& e) {std::cerr<<"[low-t-closure] "<<e.what()<<'\n';return 1;}
}
