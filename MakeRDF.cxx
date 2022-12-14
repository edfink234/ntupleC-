#include "TChain.h"
#include "TH1F.h"
#include "TBranch.h"

#include "MakeRDF.h"

#include "RDFObjects.h"

TChain RDFTree::__chain{"physics"};
TChain RDFTree::__event_info_chain{"full_event_info"};

using namespace ROOT::VecOps;

SchottDataFrame MakeRDF(const std::vector<std::string>& files)
{
    for (const auto& f: files)
    {
        RDFTree::__chain.Add(f.c_str());
        RDFTree::__event_info_chain.Add(f.c_str());
    }
    RDFTree::__chain.AddFriend(&RDFTree::__event_info_chain);

    ROOT::RDataFrame df(RDFTree::__chain);
    
    auto NewDf = df.Define("truth_particles",[&](RVec<int>& mc_pdg_id, RVec<int>& mc_barcode, RVec<float>& mc_parent_barcode, RVec<float>& mc_status, RVec<float>& mc_pt, RVec<float>& mc_charge, RVec<float>& mc_eta, RVec<float>& mc_phi, RVec<float>& mc_e, RVec<float>& mc_mass)
    {
        RVec<TruthParticle> x;
        x.reserve(mc_pt.size());
        TruthParticle temp;
        for (size_t i = 0; i < mc_pt.size(); i++)
        {
            temp.mc_pdg_id =  mc_pdg_id[i];
            temp.mc_barcode =  mc_barcode[i];
            temp.mc_parent_barcode =  mc_parent_barcode[i];
            temp.mc_status =  mc_status[i];
            temp.mc_pt =  mc_pt[i];
            temp.mc_charge =  mc_charge[i];
            temp.mc_eta =  mc_eta[i];
            temp.mc_phi =  mc_phi[i];
            temp.mc_e =  mc_e[i];
            temp.mc_mass =  mc_mass[i];
            x.push_back(temp);
        }
        return x;
    }, {"mc_pdg_id", "mc_barcode", "mc_parent_barcode", "mc_status", "mc_pt", "mc_charge", "mc_eta", "mc_phi", "mc_e", "mc_mass"})
    .Define("electrons",[&](RVec<float>& electron_charge, RVec<float>& electron_pt, RVec<float>& electron_e, RVec<float>& electron_eta, RVec<float>& electron_phi, /*RVec<float>& electron_id,*/ RVec<float>& electron_isolation, RVec<float>& electron_d0, RVec<float>& electron_z0, /*RVec<float>& electron_id_medium,*/ RVec<TruthParticle>& truth_particles)
    {
        RVec<Electron> x;
        x.reserve(electron_pt.size());
        Electron temp;
        for (size_t i = 0; i < electron_pt.size(); i++)
        {
            temp.mc_pdg_id =  truth_particles[i].mc_pdg_id;
            temp.mc_barcode =  truth_particles[i].mc_barcode;
            temp.mc_parent_barcode =  truth_particles[i].mc_parent_barcode;
            temp.mc_status =  truth_particles[i].mc_status;
            temp.mc_pt =  truth_particles[i].mc_pt;
            temp.mc_charge =  truth_particles[i].mc_charge;
            temp.mc_eta =  truth_particles[i].mc_eta;
            temp.mc_phi =  truth_particles[i].mc_phi;
            temp.mc_e =  truth_particles[i].mc_e;
            temp.mc_mass =  truth_particles[i].mc_mass;
            temp.electron_charge =  electron_charge[i];
            temp.electron_pt =  electron_pt[i];
            temp.electron_e =  electron_e[i];
            temp.electron_eta =  electron_eta[i];
            temp.electron_phi =  electron_phi[i];
//            temp.electron_id =  electron_id[i];
            temp.electron_isolation =  electron_isolation[i];
            temp.electron_d0 =  electron_d0[i];
            temp.electron_z0 =  electron_z0[i];
//            temp.electron_id_medium =  electron_id_medium[i];
            x.push_back(temp);
        }
        return x;
    }, {"electron_charge", "electron_pt", "electron_e", "electron_eta", "electron_phi", /*"electron_id",*/ "electron_isolation", "electron_d0", "electron_z0", /*"electron_id_medium",*/ "truth_particles"})
    .Define("muons",[&](RVec<float>& muon_charge, RVec<float>& muon_pt, RVec<float>& muon_e, RVec<float>& muon_eta, RVec<float>& muon_phi, RVec<TruthParticle>& truth_particles)
    {
        RVec<Muon> x;
        x.reserve(muon_pt.size());
        Muon temp;
        for (size_t i = 0; i < muon_pt.size(); i++)
        {
            temp.mc_pdg_id =  truth_particles[i].mc_pdg_id;
            temp.mc_barcode =  truth_particles[i].mc_barcode;
            temp.mc_parent_barcode =  truth_particles[i].mc_parent_barcode;
            temp.mc_status =  truth_particles[i].mc_status;
            temp.mc_pt =  truth_particles[i].mc_pt;
            temp.mc_charge =  truth_particles[i].mc_charge;
            temp.mc_eta =  truth_particles[i].mc_eta;
            temp.mc_phi =  truth_particles[i].mc_phi;
            temp.mc_e =  truth_particles[i].mc_e;
            temp.mc_mass =  truth_particles[i].mc_mass;
            temp.muon_charge =  muon_charge[i];
            temp.muon_pt =  muon_pt[i];
            temp.muon_e =  muon_e[i];
            temp.muon_eta =  muon_eta[i];
            temp.muon_phi =  muon_phi[i];
            x.push_back(temp);
        }
        return x;
    }, {"muon_charge", "muon_pt", "muon_e", "muon_eta", "muon_phi", "truth_particles"})
    .Define("photons",[&](RVec<float>& photon_pt, RVec<float>& photon_e, RVec<float>& photon_eta, RVec<float>& photon_phi,  RVec<float>& photon_etcone40, RVec<int>& photon_id, RVec<int>& photon_id_loose, RVec<int>& photon_id_tight, RVec<int>& photon_cluster_eta_be_2, /*RVec<int>& photon_id_nn,*/ RVec<TruthParticle>& truth_particles)
    {
        RVec<Photon> x;
        x.reserve(photon_pt.size());
        Photon temp;
        for (size_t i = 0; i < photon_pt.size(); i++)
        {
            temp.mc_pdg_id =  truth_particles[i].mc_pdg_id;
            temp.mc_barcode =  truth_particles[i].mc_barcode;
            temp.mc_parent_barcode =  truth_particles[i].mc_parent_barcode;
            temp.mc_status =  truth_particles[i].mc_status;
            temp.mc_pt =  truth_particles[i].mc_pt;
            temp.mc_charge =  truth_particles[i].mc_charge;
            temp.mc_eta =  truth_particles[i].mc_eta;
            temp.mc_phi =  truth_particles[i].mc_phi;
            temp.mc_e =  truth_particles[i].mc_e;
            temp.mc_mass =  truth_particles[i].mc_mass;
            temp.photon_pt =  photon_pt[i];
            temp.photon_e =  photon_e[i];
            temp.photon_eta =  photon_eta[i];
            temp.photon_phi =  photon_phi[i];
            temp.photon_etcone40 =  photon_etcone40[i];
            temp.photon_id =  photon_id[i];
            temp.photon_id_loose =  photon_id_loose[i];
            temp.photon_id_tight =  photon_id_tight[i];
            temp.photon_cluster_eta_be_2 =  photon_cluster_eta_be_2[i];
//            temp.photon_id_nn = photon_id_nn[i];
            x.push_back(temp);
        }
        return x;
    }, {"photon_pt", "photon_e", "photon_eta", "photon_phi", "photon_etcone40", "photon_id", "photon_id_loose", "photon_id_tight", "photon_cluster_eta_be_2", /*"photon_id_nn",*/ "truth_particles"})
    .Define("clusters",[&](RVec<float>& cluster_pt, RVec<float>& cluster_eta, RVec<float>& cluster_phi, RVec<float>& cluster_e)
    {
        RVec<Cluster> x;
        x.reserve(cluster_pt.size());
        Cluster temp;
        for (size_t i = 0; i < cluster_pt.size(); i++)
        {
            temp.cluster_pt =  cluster_pt[i];
            temp.cluster_e =  cluster_e[i];
            temp.cluster_eta =  cluster_eta[i];
            temp.cluster_phi =  cluster_phi[i];
            x.push_back(temp);
        }
        return x;
    }, {"cluster_pt", "cluster_phi", "cluster_e", "cluster_eta"})
    .Define("tracks",[&](RVec<float>& track_pt, RVec<float>& track_eta, RVec<float>& track_phi, /*RVec<float>& track_e,*/ RVec<float>& track_charge, RVec<int>& track_num_pixel_hits, RVec<int>& track_num_sct_hits)
    {
        RVec<Track> x;
        x.reserve(track_pt.size());
        Track temp;
        for (size_t i = 0; i < track_pt.size(); i++)
        {
            temp.track_pt =  track_pt[i];
            temp.track_eta =  track_eta[i];
            temp.track_phi =  track_phi[i];
//            temp.track_e =  track_e[i];
            temp.track_charge =  track_charge[i];
            temp.track_num_pixel_hits =  track_num_pixel_hits[i];
            temp.track_num_sct_hits =  track_num_sct_hits[i];
            x.push_back(temp);
        }
        return x;
    }, {"track_pt", "track_charge", "track_eta", "track_phi", /*"track_e",*/ "track_num_pixel_hits", "track_num_sct_hits"});
    
    return NewDf;
}

