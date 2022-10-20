#ifndef EVENTH
#define EVENTH

#include "objects.h"
#include "TInterpreter.h"




class Event
{
private:
    
//    TChain entry;
//    int run_number;
//    int random_run_number;
//    int event_number;
//    std::vector<Photon> photons;
//    std::vector<Electron> electrons;
//    std::vector<Cluster> clusters;
//    std::vector<Track> tracks;
//    std::vector<Track> pixel_tracks;
//    std::vector<TruthParticle> truth_particles;
//    std::vector<string> triggers; // std::vector<std::vector<string>> triggers; //maybe?
//    int muon_spectrometer_num_track_particles;
    
    void __load_photons(TChain* entry);
    void __load_electrons(TChain* entry);
    void __load_clusters(TChain* entry);
    void __load_tracks(TChain* entry);
    void __load_truth_particles(TChain* entry);
    void __load_triggers(TChain* entry);
    
public:
    Event(TChain* entry, TChain* event_info_entry = nullptr);
    Event();
    ~Event();
    std::vector<TruthParticle> find_truth_particles
    (const std::vector<int>&& barcode = {},
     const std::vector<int>&& parent_barcode = {},
     const std::vector<int>&& pdg_id = {},
     int* status_code = nullptr);
//    int run_number;
//    Event();
//    ~Event();
    Event(const Event&);
    Event& operator=(const Event&);
    
    
    static bool cache_truth;
    static bool load_reco;
    static bool load_photons;
    static bool load_electrons;
    static bool load_clusters;
    static bool load_tracks;
    static bool load_triggers;
    static string systematic;
    
    TChain entry;
    int run_number;
    int random_run_number;
    int event_number;
    std::vector<Photon> photons;
    std::vector<Electron> electrons;
    std::vector<Cluster> clusters;
    std::vector<Track> tracks;
    std::vector<Track> pixel_tracks;
    std::vector<TruthParticle> truth_particles;
    std::vector<string> triggers; // std::vector<std::vector<string>> triggers; //maybe?
    int muon_spectrometer_num_track_particles;
};


#endif
