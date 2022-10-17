#ifndef EVENTH
#define EVENTH

#include "objects.h"


class Event
{
private:
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
    vector<Photon> photons;
    vector<Electron> electrons;
    vector<Cluster> clusters;
    vector<Track> tracks;
    vector<Track> pixel_tracks;
    vector<TruthParticle> truth_particles;
    vector<string> triggers; // vector<vector<string>> triggers; //maybe?
    int muon_spectrometer_num_track_particles;
    
    void __load_photons(TChain* entry);
    void __load_electrons(TChain* entry);
    void __load_clusters(TChain* entry);
    void __load_tracks(TChain* entry);
    void __load_truth_particles(TChain* entry);
    void __load_triggers(TChain* entry);
    
public:
    Event(TChain* entry, TChain* event_info_entry = nullptr);
    vector<TruthParticle>& find_truth_particles(const vector<int>&&,
                                                const vector<int>&&, const vector<int>&&, int);
    
//    Event();
//    ~Event();
//    Event(const Event&);
//    Event& operator=(const Event&);
    
};


#endif
