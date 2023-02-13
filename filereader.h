#ifndef FILEREADERH
#define FILEREADERH

#include <vector>
#include <string>

#include "event.h"
#include "objects.h"

using FT = bool(TChain&); //alias for filter function type

/*
 FileReader: helper class for FileReaderRange
 */
class FileReader
{
private:
    friend class FileReaderRange;
    std::vector<std::string> __files;
    int __skip_first_events;
    int __current_index;
    std::vector<FT*> __event_filters;
    TChain __chain;
    TChain __event_info_chain;
    Long64_t __num_events;
    bool __has_event_info_chain;
    
    void __load_photon_addresses();
    void __load_electron_addresses();
    void __load_cluster_addresses();
    void __load_track_addresses();
    void __load_truth_particle_addresses(bool cache_truth);
    void __load_trigger_addresses();

    void __load_photons();
    void __load_electrons();
    void __load_muons();
    void __load_clusters();
    void __load_tracks();
    void __load_truth_particles();
    void __load_triggers();
    
public:
    FileReader();
    FileReader(const std::vector<std::string>& files, const char* tree_name = "physics", Long64_t num_events = -1, int skip_first_events = 0);
    ~FileReader();
    FileReader(const FileReader&);
    template <typename T>
    void add_event_filter(T&);
    bool __passes_event_filters();
    int current_index();
    Long64_t num_events();
    Event __current_event;
    
    std::vector<TruthParticle> find_truth_particles
    (const std::vector<int>&& barcode = {},
     const std::vector<int>&& parent_barcode = {},
     const std::vector<int>&& pdg_id = {},
     int* status_code = nullptr,
     bool inv = false);
};

/*
 FileReaderRange: class to easily iterate over a group of files with the same TTrees
 */
class FileReaderRange {
    FileReader f;
public:
     FileReaderRange(const std::vector<std::string>& files,
           const char* tree_name =
//           "../user.kschmied.28655874._000025.LGNTuple.root?#physics",
           "physics",
           Long64_t num_events = -1, int skip_first_events = 0);

     struct Iterator
     {
         Iterator(int i);
         Iterator (int i, FileReader& f);
         Iterator& operator++();

         FileReader& operator*();
         friend bool operator!=(const Iterator& a, const Iterator& b);
         private:
         
         int data;
         FileReader F;
     };

    Iterator begin();
    Iterator end();
};

#endif



