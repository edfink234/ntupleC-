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

//bool Event::cache_truth = true;
//bool Event::load_reco= true;
//bool Event::load_photons= true;
//bool Event::load_electrons= true;
//bool Event::load_clusters= true;
//bool Event::load_tracks= true;
//bool Event::load_triggers= true;
//string Event::systematic= "";


//class Event
//{
//private:
    
    
//    void __load_photons(TChain* entry, int entry_number);
//    void __load_electrons(TChain* entry, int entry_number);
//    void __load_clusters(TChain* entry, int entry_number);
//    void __load_tracks(TChain* entry, int entry_number);
//    void __load_truth_particles(TChain* entry, int entry_number);
//    void __load_triggers(TChain* entry, int entry_number);
    
//public:
//    Event(TChain* entry, TChain* event_info_entry = nullptr, int entry_number = 0);
//    Event();
//    ~Event();
//    std::vector<TruthParticle> find_truth_particles
//    (const std::vector<int>&& barcode = {},
//     const std::vector<int>&& parent_barcode = {},
//     const std::vector<int>&& pdg_id = {},
//     int* status_code = nullptr);
//    int run_number;
//    Event();
//    ~Event();
//    Event(const Event&);
//    Event& operator=(const Event&);
    
    
//    static bool cache_truth;
//    static bool load_reco;
//    static bool load_photons;
//    static bool load_electrons;
//    static bool load_clusters;
//    static bool load_tracks;
//    static bool load_triggers;
//    static string systematic;
    
//    TChain entry;
//    int entry_number;
//    int run_number;
//    int random_run_number;
//    int event_number;
//    std::vector<Photon> photons;
//    std::vector<Electron> electrons;
//    std::vector<Cluster> clusters;
//    std::vector<Track> tracks;
//    std::vector<Track> pixel_tracks;
//    std::vector<TruthParticle> truth_particles;
//    std::vector<string> triggers; 
//    int muon_spectrometer_num_track_particles;
//};


#endif
