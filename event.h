#ifndef EVENTH
#define EVENTH

#include "objects.h"



struct Event
{
    UInt_t m_RunNumber;
    UInt_t m_RandomRunNumber;
    UInt_t ei_event_number;
    std::vector<double> *photon_pt= nullptr;
    std::vector<double> *photon_e = nullptr;
    std::vector<std::vector<string>> *photon_syst_name = nullptr;
    std::vector<std::vector<double>> *photon_syst_pt = nullptr;
    std::vector<std::vector<double>> *photon_syst_e = nullptr;
    std::vector<double> *electron_pt= nullptr;
    std::vector<double> *electron_e = nullptr;
    std::vector<std::vector<string>> *electron_syst_name = nullptr;
    std::vector<std::vector<double>> *electron_syst_pt = nullptr;
    std::vector<std::vector<double>> *electron_syst_e = nullptr;
    std::vector<double> *cluster_pt= nullptr;
    std::vector<double> *track_pt= nullptr;
    std::vector<string> *trigger_passed_triggers = nullptr;
    std::vector<double> *mc_pt= nullptr;
    std::vector<int> *mc_pdg_id = nullptr;
    std::vector<int> *mc_barcode = nullptr;
    std::vector<int> *mc_parent_barcode = nullptr;
    std::vector<int> *mc_status = nullptr;
    std::vector<int>* track_type = nullptr;

    int entry_number;
    int run_number;
    int random_run_number;
    int event_number;
    std::vector<Photon> photons;
    std::vector<Electron> electrons;
    std::vector<Cluster> clusters;
    std::vector<Track> tracks;
    std::vector<Track> pixel_tracks;
    std::vector<TruthParticle> truth_particles;
    std::vector<string> triggers;
    int muon_spectrometer_num_track_particles;

    static bool cache_truth;
    static bool load_reco;
    static bool load_photons;
    static bool load_electrons;
    static bool load_clusters;
    static bool load_tracks;
    static bool load_triggers;
    static string systematic;
};




#endif
