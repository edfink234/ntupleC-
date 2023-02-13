#include "filereader.h"

//                ************************
//                *      FileReader      *
//                ************************

/*
 FileReader constructor: Takes in a vector of files, and, optionally, the tree name,
 the number of events to process, and the number of initial events to skip. An additional
 TTree called full_event_info is always assumed to be inside all of the ROOT files.
 
 The constructor adds the files to the internal chain attributes, initializes
 the number of events to skip and to process.
 */

FileReader::FileReader(const std::vector<std::string>& files, const char* tree_name, Long64_t num_events, int skip_first_events) :

__files{move(files)},
__skip_first_events{skip_first_events},
__current_index{skip_first_events},
__chain{tree_name},
__event_info_chain{"full_event_info"},
__has_event_info_chain{true}

{
    for (const auto& f: __files)
    {
        __chain.Add(f.c_str());
        __event_info_chain.Add(f.c_str());
    }

    if ((num_events<0) || (num_events>__chain.GetEntries()))
    {
        __num_events = __chain.GetEntries();
    }
    else
    {
        __num_events = num_events;
    }
    __current_event.entry_number = skip_first_events;
}

FileReader::FileReader() = default;

FileReader::~FileReader() 
{
    __chain.Reset();
    __event_info_chain.Reset();
}

FileReader::FileReader(const FileReader& other)
{
    __files = other.__files;
    __skip_first_events = other.__skip_first_events;
    __current_index = other.__current_index;
    __current_event = other.__current_event;
    __event_filters = other.__event_filters;
    __has_event_info_chain = other.__has_event_info_chain;
    
    __chain.Reset();
    __chain.Add(const_cast<TChain*>(&other.__chain)); //😬
    __event_info_chain.Reset();
    __event_info_chain.Add(const_cast<TChain*>(&other.__event_info_chain)); //😬
    __num_events = other.__num_events;
}

/*
 Utility function to add bool(TChain&) functions
 to filter events
 */
template <typename T>
void FileReader::add_event_filter(T& filter)
{
    __event_filters.push_back(filter);
}

/*
 Checks if the event passes all event filters
 added by FileReader::add_event_filter
 */
bool FileReader::__passes_event_filters()
{
    if (__event_filters.empty())
    {
        return true;
    }
    for (auto& filter_func: __event_filters)
    {
        if (!((*filter_func)(__chain)))
        {
            return false;
        }
    }
    return true;
}

int FileReader::current_index()
{
    return __current_index;
}

Long64_t FileReader::num_events()
{
    return __num_events;
}

/*
 Sets the branch status of branches in the TChain `__chain', as
 well as set the Photon static attributes of Event
 to the addresses of those branches.
 */

void FileReader::__load_photon_addresses()
{
    __chain.SetBranchStatus("photon_syst_name",1);
    __chain.SetBranchStatus("photon_syst_pt",1);
    __chain.SetBranchStatus("photon_syst_e",1);
    
    __chain.SetBranchAddress("photon_syst_name",&(__current_event.photon_syst_name));
    __chain.SetBranchAddress("photon_syst_pt",&(__current_event.photon_syst_pt));
    __chain.SetBranchAddress("photon_syst_e",&(__current_event.photon_syst_e));
}

/*
 Sets the branch status of branches in the TChain `__chain', as
 well as set the Electron static attributes of Event
 to the addresses of those branches.
 */

void FileReader::__load_electron_addresses()
{
    __chain.SetBranchStatus("electron_syst_name",1);
    __chain.SetBranchStatus("electron_syst_pt",1);
    __chain.SetBranchStatus("electron_syst_e",1);
    
    __chain.SetBranchAddress("electron_syst_name",&(__current_event.electron_syst_name));
    __chain.SetBranchAddress("electron_syst_pt",&(__current_event.electron_syst_pt));
    __chain.SetBranchAddress("electron_syst_e",&(__current_event.electron_syst_e));
}

/*
 Sets the branch status of branches in the TChain `__chain', as
 well as set the Cluster static attributes of Event
 to the addresses of those branches. Currently only applies to branch
 cluster_pt 😅
 */

void FileReader::__load_cluster_addresses()
{
    __chain.SetBranchStatus("cluster_pt",1);
    __chain.SetBranchAddress("cluster_pt",&(__current_event.cluster_pt));
}

/*
 Sets the branch status of branches in the TChain `__chain', as
 well as set the Track static attributes of Event
 to the addresses of those branches. Currently only applies to branch
 track_type, but it's not as common so we check if it exists before
 trying to set it to avoid a possible run-time error.
 */

void FileReader::__load_track_addresses()
{
//    __chain.SetBranchStatus("track_pt",1);
//    __chain.SetBranchAddress("track_pt",&(__current_event.track_pt));
    if (__chain.GetListOfBranches()->FindObject("track_type"))
    {
        __chain.SetBranchStatus("track_type",1);
        __chain.SetBranchAddress("track_type",&(__current_event.track_type));
    }
}

/*
 Sets the branch status of branches in the TChain `__chain', as
 well as set the TruthParticle static attributes of Event
 to the addresses of those branches. Only applies if the
 cache_truth static variable is set to false
 */

void FileReader::__load_truth_particle_addresses(bool cache_truth)
{
    if (!cache_truth)
    {
        __chain.SetBranchStatus("mc_pdg_id",1);
        __chain.SetBranchStatus("mc_barcode",1);
        __chain.SetBranchStatus("mc_parent_barcode",1);
        __chain.SetBranchStatus("mc_status",1);
        
        __chain.SetBranchAddress("mc_pdg_id",&(__current_event.mc_pdg_id));
        __chain.SetBranchAddress("mc_barcode",&(__current_event.mc_barcode));
        __chain.SetBranchAddress("mc_parent_barcode",&(__current_event.mc_parent_barcode));
        __chain.SetBranchAddress("mc_status",&(__current_event.mc_status));
    }
}

/*
 if trigger_passed_triggers is in the list of branches, then turns the
 trigger_passed_triggers branch status to "on", and sets the FileReader's
 event's trigger_passed_triggers attribute's address to the
 trigger_passed_triggers branch.
 */

void FileReader::__load_trigger_addresses()
{
    
    if (__chain.GetListOfBranches()->FindObject("trigger_passed_triggers"))
    {
        __chain.SetBranchStatus("trigger_passed_triggers",1);
        __chain.SetBranchAddress("trigger_passed_triggers",&(__current_event.trigger_passed_triggers));
    }
}

/*
 Physically loads Photon objects into the FileReader's event
 attribute (in practice, for the current event). Will set the pt and
 energy of the Photon objects to the systematic values
 if Event::systematic is set to something besides "nominal" or "".
 */

void FileReader::__load_photons()
{
    auto length = (*(Photon::photon_pt)).size();
    double pt, energy;
    int index;
    
    for (size_t i = 0; i < length; i++)
    {
        if ((Event::systematic.empty()) || (Event::systematic == "nominal"))
        {
            pt = (*(Photon::photon_pt))[i];
            energy = (*(Photon::photon_e))[i];
        }
        else
        {
            index = -1;
            auto it = find(((*(__current_event.photon_syst_name))[i]).begin(), ((*(__current_event.photon_syst_name))[i]).end(), Event::systematic);
            
            if (it == ((*(__current_event.photon_syst_name))[i]).end())
            {
                continue;
            }
            index = it - ((*(__current_event.photon_syst_name))[i]).begin();
            pt = (*(__current_event.photon_syst_pt))[i][index];
            energy = (*(__current_event.photon_syst_e))[i][index];
        }
        if (pt < 0)
        {
            continue;
        }

        __current_event.photons.emplace_back(Photon(i, pt, energy));
    }
}

/*
 Physically loads Electron objects into the FileReader's event
 attribute (in practice, for the current event). Will set the pt and
 energy of the Electron objects to the systematic values
 if Event::systematic is set to something besides "nominal" or "".
 */

void FileReader::__load_electrons()
{
    auto length = (*(Electron::electron_pt)).size();
    double pt, energy;
    int index;
//    printf("length = %ld\n",length);
    for (size_t i = 0; i < length; i++)
    {
        if ((Event::systematic.empty()) || (Event::systematic == "nominal"))
        {
            pt = (*(Electron::electron_pt))[i];
            energy = (*(Electron::electron_e))[i];
        }
       
        else
        {
            index = -1;
            auto it = find(((*(__current_event.electron_syst_name))[i]).begin(), ((*(__current_event.electron_syst_name))[i]).end(), Event::systematic);
            
            if (it == ((*(__current_event.electron_syst_name))[i]).end())
            {
                continue;
            }
            index = it - ((*(__current_event.electron_syst_name))[i]).begin();
            pt = (*(__current_event.electron_syst_pt))[i][index];
            energy = (*(__current_event.electron_syst_e))[i][index];
        }
        if (pt < 0)
        {
            continue;
        }
        __current_event.electrons.emplace_back(Electron(i, pt, energy));
    }
}

/*
 Physically loads Muon objects into the FileReader's event
 attribute (in practice, for the current event). Will set the pt and
 energy of the Muon objects to the systematic values
 if Event::systematic is set to something besides "nominal" or "".
 However, currently there are no systematics implemented for Muons
 for our analysis, so that code is commented out.
 */

void FileReader::__load_muons()
{
    auto length = (*(Muon::muon_pt)).size();
    double pt, energy;
    int index;
    
    for (size_t i = 0; i < length; i++)
    {
        if ((Event::systematic.empty()) || (Event::systematic == "nominal"))
        {
            pt = (*(Muon::muon_pt))[i];
            energy = (*(Muon::muon_e))[i];
        }
//        else
//        {
//            index = -1;
//            auto it = find(((*(__current_event.electron_syst_name))[i]).begin(), ((*(__current_event.electron_syst_name))[i]).end(), Event::systematic);
//
//            if (it == ((*(__current_event.electron_syst_name))[i]).end())
//            {
//                continue;
//            }
//            index = it - ((*(__current_event.electron_syst_name))[i]).begin();
//            pt = (*(__current_event.electron_syst_pt))[i][index];
//            energy = (*(__current_event.electron_syst_e))[i][index];
//        }
//        if (pt < 0)
//        {
//            continue;
//        }
        __current_event.muons.emplace_back(Muon(i, pt, energy));
    }
}

/*
 Physically loads Cluster objects into the FileReader's event
 attribute (in practice, for the current event). Specifically, loops over
 the static cluster_pt attribute of cluster and pushes back Cluster
 objects from i = 0 to i = cluster_pt.size()-1 into the FileReader's
 Event's clusters vector attribute
 */

void FileReader::__load_clusters()
{
    auto length = (*(Cluster::cluster_pt)).size();
    
    for (size_t i = 0; i < length; i++)
    {
        __current_event.clusters.emplace_back(Cluster(i));
    }
}

/*
 Physically loads Track objects into the FileReader's event
 attribute (in practice, for the current event). Specifically,
 loops over the static track_pt attribute of cluster and pushes back Cluster
 objects from i = 0 to i = cluster_pt.size()-1 into the FileReader's
 Event's tracks or pixel_tracks vector attributes, depending on if
 Event.track_type is 0 or 1 respectively
 */

void FileReader::__load_tracks()
{
    auto length = (*(Track::track_pt)).size();
    int track_type;
    for (size_t i = 0; i < length; i++)
    {
        if (!__current_event.track_type)
        {
            track_type = 0;
        }
        else
        {
            track_type = (*__current_event.track_type)[i];
        }
        if (track_type==0)
        {
            __current_event.tracks.emplace_back(Track(i));
        }
        else if (track_type==1)
        {
            __current_event.pixel_tracks.emplace_back(Track(i));
        }
    }
}

/*
 Physically loads TruthParticle objects into the FileReader's event
 attribute (in practice, for the current event). Specifically, loops over
 the static mc_pt attribute of TruthParticle and pushes back TruthParticle
 objects from i = 0 to i = mc_pt.size()-1 into the FileReader's
 Event's truth_particles vector attribute
 */

void FileReader::__load_truth_particles()
{
    auto length = (*(TruthParticle::mc_pt)).size();
    for (size_t i = 0; i < length; i++)
    {
        __current_event.truth_particles.emplace_back(TruthParticle(i));
    }
}

/*
 Physically loads std::strings of triggers into the FileReader's event
 attribute (in practice, for the current event). Specifically, if the
 FileReader's Event's trigger_passed_triggers attribute was set by
 __load_trigger_addresses (because the branch trigger_passed_triggers
 exists in the __chain TTree), then __load_triggers loops
 over the static trigger_passed_triggers attribute
 of Event and pushes back std::strings of triggers from i = 0 to
 i = trigger_passed_triggers.size()-1 into the FileReader's
 Event's triggers vector attribute
 */

void FileReader::__load_triggers()
{
    if (__current_event.trigger_passed_triggers)
    {
        auto length = (*(__current_event.trigger_passed_triggers)).size();
        for (size_t i = 0; i < length; i++)
        {
            __current_event.triggers.push_back((*(__current_event.trigger_passed_triggers))[i]);
        }
    }
}

/*
 if Event::cache_truth (meaning if the user has specified loading truth_particles)
 ---------------------------------------------------------------------------------
 In this case, this function will return a list of truth particles whose barcodes,
 parent barcodes, and absolute values of the pdg ids are contained within the respective
 arguments passed to this function if they're not empty.
 Also in this case, the truth particles returned
 will all have status code equal to the status code passed to this function,
 unless none was passed. In the case inv = True, then truth particles returned
 will all have status code not equal to the one passed to this funciton
 
 if not Event::cache_truth (meaning if the user has specified to not load truth_particles)
 -----------------------------------------------------------------------------------------
 In this case, this function will loop over the mc_pdg_id branch stored in
 the FileReader's Event attribute (which is loaded by __load_truth_particle_addresses
 if Event::cache_truth = false), and will return a list of truth particles whose barcodes,
 parent barcodes, and absolute values of the pdg ids are contained by the corresponding mc
 branches and within the respective arguments passed to this function
 if they're not empty (or null in the case of status_code).
 */

std::vector<TruthParticle> FileReader::find_truth_particles
     (const std::vector<int>&& barcode,
      const std::vector<int>&& parent_barcode,
      const std::vector<int>&& pdg_id,
      int* status_code,
      bool inv)
{
    std::vector<TruthParticle> results;
    
    if (Event::cache_truth)
    {
        for (auto& tp : __current_event.truth_particles)
        {
            if ((!barcode.empty()) && (find(barcode.begin(),barcode.end(),tp.barcode()) == barcode.end()))
            {
                continue;
            }
            if ((!parent_barcode.empty()) && (find(parent_barcode.begin(),parent_barcode.end(),tp.parent_barcode()) == parent_barcode.end()))
            {
                continue;
            }
            if ((!pdg_id.empty()) && (find(pdg_id.begin(),pdg_id.end(),abs(tp.pdg_id)) == pdg_id.end()))
            {
                continue;
            }
            if (status_code && (tp.status_code() != *status_code) && (!inv))
            {
                continue;
            }
            if (status_code && (tp.status_code() == *status_code) && inv)
            {
                continue;
            }
            results.push_back(tp);
        }
    }
    else
    {
        auto length = (*__current_event.mc_pdg_id).size();

        for (size_t i = 0; i < length; i++)
        {
            if (!barcode.empty())
            {
                if (find(barcode.begin(),barcode.end(),(*(__current_event.mc_barcode))[i]) == barcode.end())
                {
                    continue;
                }
            }
            if (!parent_barcode.empty())
            {
                if (find(parent_barcode.begin(),parent_barcode.end(),(*(__current_event.mc_parent_barcode))[i]) == parent_barcode.end())
                {
                    continue;
                }
            }
            if (!pdg_id.empty())
            {
                if (find(pdg_id.begin(),pdg_id.end(),(*(__current_event.mc_pdg_id))[i]) == pdg_id.end())
                {
                    continue;
                }
            }
            if (status_code)
            {
                if ((*__current_event.mc_status)[i] != *status_code)
                {
                    continue;
                }
            }
            results.push_back(TruthParticle(i));
        }
    }
    return results;
}

//                *****************************
//                *      FileReaderRange      *
//                *****************************

/*
 FileReaderRange constructor: Takes in a vector of files, and, optionally, the tree name,
 the number of events to process, and the number of initial events to skip. An additional
 TTree called full_event_info is always assumed to be inside all of the ROOT files.
 
 The constructor calls the copy constructor of FileReader to initialize its
 internal FileReader attribute
 
 FileReaderRange allows the user to loop over TChains of files in an efficient and
 hassle-free manner.
 
 e.g.
 
 FileReaderRange reader(input_filenames); //create FileReaderRange object
 
 for (auto&& f: reader) //now iterate, easy as that!
 {
    //do something
 }
 */

FileReaderRange::FileReaderRange(const std::vector<std::string>& files,
       const char* tree_name,
       Long64_t num_events, int skip_first_events) :
f(files,tree_name, num_events, skip_first_events)
{
    
}

/*
 Iterator constructor called by end() method, takes in number of events to process
 and initalizes internal data attribute to this value.
 */
FileReaderRange::Iterator::Iterator(int i) : data{i} {}

/*
 Iterator constructor called by begin() method, takes in number of events to skip
 (event number to start at) and the internal FileReader object of FileReaderRange
 used to initialize the Iterator's data and FileReader attributes, respectively.
 
 Starts by deactivating all of the branches, then sets only the ones specified by
 the static attributes of Event. Then calls GetEntry once and subsequently loads
 the relevant objects into the FileReader's Event's vector attributes.
 */
FileReaderRange::Iterator::Iterator(int i, FileReader& f) : data{i}, F{f}
{
    F.__current_index = F.__current_event.entry_number = data;
    
    F.__current_event.truth_particles.reserve(6);
    F.__current_event.photons.reserve(6);
    F.__current_event.triggers.reserve(6);
    F.__current_event.tracks.reserve(6);
    F.__current_event.pixel_tracks.reserve(6);
    F.__current_event.clusters.reserve(6);
    F.__current_event.electrons.reserve(6);
    F.__current_event.muons.reserve(6);
    
    //Decativate all branches
    F.__chain.SetBranchStatus("*",0);
    F.__event_info_chain.SetBranchStatus("*",0);
    
    //SetAddresses
    if (F.__has_event_info_chain)
    {
        F.__event_info_chain.SetBranchStatus("m_RunNumber",1);
        F.__event_info_chain.SetBranchStatus("m_RandomRunNumber",1);
        F.__chain.SetBranchStatus("ei_event_number",1);

        F.__event_info_chain.SetBranchAddress("m_RunNumber",&(F.__current_event.m_RunNumber));
        F.__event_info_chain.SetBranchAddress("m_RandomRunNumber",&(F.__current_event.m_RandomRunNumber));
        F.__chain.SetBranchAddress("ei_event_number",&(F.__current_event.ei_event_number));
    }

    F.__load_truth_particle_addresses(Event::cache_truth);
        
    if (Event::cache_truth)
    {
        TruthParticle::SetTruthParticle(&F.__chain);
    }
        
    if (Event::load_reco)
    {
        if (Event::load_photons)
        {
            F.__load_photon_addresses();
            Photon::SetPhoton(&F.__chain);
        }
        
        if (Event::load_electrons)
        {
            F.__load_electron_addresses();
            Electron::SetElectron(&F.__chain);
        }
        
        if (Event::load_muons)
        {
            Muon::SetMuon(&F.__chain);
        }
        
        if (Event::load_clusters)
        {
//            F.__load_cluster_addresses();
            Cluster::SetCluster(&F.__chain);
        }
        
        if (Event::load_tracks)
        {
            F.__load_track_addresses();
            Track::SetTrack(&F.__chain);
        }
        
        if (Event::load_triggers)
        {
            F.__load_trigger_addresses();
        }
    }

    //GetAddresses
    if (F.__has_event_info_chain)
    {
        F.__event_info_chain.GetEntry(F.__current_index);
    }
    F.__chain.GetEntry(F.__current_index);
    
    if (Event::cache_truth)
    {
        F.__load_truth_particles();
    }
    
    if (Event::load_reco)
    {
        if (Event::load_photons)
        {
            F.__load_photons();
        }
        if (Event::load_electrons)
        {
            F.__load_electrons();
        }
        if (Event::load_muons)
        {
            F.__load_muons();
        }
        if (Event::load_clusters)
        {
            F.__load_clusters();
        }
        if (Event::load_tracks)
        {
            F.__load_tracks();
        }
        if (Event::load_triggers)
        {
            F.__load_triggers();
        }
    }
}

/*
 ++ operator of the custom iterator, describes how to procede to the
 next iteration, similar to Python's __next__ dunder method. First starts
 by clearing the vector attributes of the current event, increases the data
 attribute, calls GetEntry to load the addresses for the next event, and
 finally loads the relevant objects into the FileReader's Event's vector attributes.
 */
FileReaderRange::Iterator& FileReaderRange::Iterator::operator++()
{
    F.__current_event.truth_particles.clear();
    F.__current_event.photons.clear();
    F.__current_event.triggers.clear();
    F.__current_event.tracks.clear();
    F.__current_event.pixel_tracks.clear();
    F.__current_event.clusters.clear();
    F.__current_event.electrons.clear();
    F.__current_event.muons.clear();

    F.__current_index = F.__current_event.entry_number = ++data;
    
    if (F.__has_event_info_chain)
    {
        F.__event_info_chain.GetEntry(F.__current_index);
    }
    F.__chain.GetEntry(F.__current_index);
    
    if (Event::cache_truth)
    {
        F.__load_truth_particles();
    }
    
    if (Event::load_reco)
    {
        if (Event::load_photons)
        {
            F.__load_photons();
        }
        if (Event::load_electrons)
        {
            F.__load_electrons();
        }
        if (Event::load_muons)
        {
            F.__load_muons();
        }
        if (Event::load_clusters)
        {
            F.__load_clusters();
        }
        if (Event::load_tracks)
        {
            F.__load_tracks();
        }
        if (Event::load_triggers)
        {
            F.__load_triggers();
        }
    }
    
    return *this;
}
/*
 When you iterate over a FileReaderRange object, i.e.
 
 for (auto&& i: FileReaderRange(filenames))
 {
    //...
 }
 the variable 'i' is described by this derefencing operator,
 i.e., it's a FileReader reference.
 */
FileReader& FileReaderRange::Iterator::operator*()
{
    return F;
}

/*
 Determines when the iteration stops, i.e., when the data of
 the Iterator object constructed by begin() is equal to
 the Iterator object constructed by end()
 */
bool operator!=(const FileReaderRange::Iterator& a, const FileReaderRange::Iterator& b)
{
    return a.data != b.data;
}

/*
 Function that generates an Iterator object to describe how the iteration
 over a FileReaderRange object begins.
 */
FileReaderRange::Iterator FileReaderRange::begin()
{
    return Iterator(f.__skip_first_events, f);
}

/*
 Function that generates an Iterator object to describe how the iteration
 over a FileReaderRange object ends.
 */
FileReaderRange::Iterator FileReaderRange::end()
{
    return Iterator(f.__num_events);
//    return Iterator(100);
}


