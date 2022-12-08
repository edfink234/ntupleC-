#ifndef RDFOBJECTSH
#define RDFOBJECTSH

#include <vector>
#include <string>

#include "Math/Vector4D.h"
#include "ROOT/RDataFrame.hxx"

using namespace ROOT::Math;
using namespace ROOT::VecOps;

struct TruthParticle
{
    static const std::string PREFIX;
    
    RVec<int>   mc_pdg_id;
    RVec<int>   mc_barcode;
    RVec<int>   mc_parent_barcode;
    RVec<int>   mc_status;
    RVec<float> mc_pt;
    RVec<float> mc_charge;
    RVec<float> mc_eta;
    RVec<float> mc_phi;
    RVec<float> mc_e;
    RVec<float> mc_mass;
};

struct Electron final : public TruthParticle
{
    static const std::string PREFIX;
    static const int PDG_ID;
    
    RVec<float> electron_charge;
    RVec<float> electron_pt;
    RVec<float> electron_e;
    RVec<float> electron_eta;
    RVec<float> electron_phi;
    RVec<int>   electron_id;
    RVec<float> electron_isolation;
    RVec<float> electron_d0;
    RVec<float> electron_z0;
    RVec<int>   electron_id_medium;
};

struct Muon final : public TruthParticle
{
    const std::string PREFIX;
    static const int PDG_ID;
    
    RVec<float> muon_charge;
    RVec<float> muon_pt;
    RVec<float> muon_e;
    RVec<float> muon_eta;
    RVec<float> muon_phi;
};

struct Photon final : public TruthParticle
{
    static const int PDG_ID;

    RVec<float> photon_pt;
    RVec<float> photon_e;
    RVec<float> photon_eta;
    RVec<float> photon_phi;
    RVec<float> photon_etcone40;
    RVec<int>   photon_id;
    RVec<int>   photon_id_loose;
    RVec<int>   photon_id_tight;
    RVec<int>   photon_cluster_eta_be_2;
    RVec<int>   photon_id_nn;
};

struct Cluster final : public PhysicsObject
{
    static const std::string PREFIX;
    
    RVec<float> cluster_pt;
    RVec<float> cluster_eta;
    RVec<float> cluster_phi;
    RVec<float> cluster_e;
};

struct Track final : public PhysicsObject
{
    static const std::string PREFIX;
    
    RVec<float> track_pt;
    RVec<float> track_charge;
    RVec<float> track_eta;
    RVec<float> track_phi;
    RVec<float> track_e;
    RVec<int>   track_num_pixel_hits;
    RVec<int>   track_num_sct_hits;
};

#endif

