#include "filereader.h"


FileReader::FileReader(std::vector<string>& files, const char* tree_name, Long64_t num_events, int skip_first_events) :

__files{move(files)},
__skip_first_events{skip_first_events},
__current_index{skip_first_events},
__chain{tree_name},
__event_info_chain{"full_event_info"},
__has_event_info_chain{true}


{
    for (auto& f: __files)
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
    
    __chain.Reset();
    __chain.Add(const_cast<TChain*>(&other.__chain)); //😬
    __event_info_chain.Reset();
    __event_info_chain.Add(const_cast<TChain*>(&other.__event_info_chain)); //😬
    __num_events = other.__num_events;
}

template <typename T>
void FileReader::add_event_filter(T& filter)
{
    __event_filters.push_back(filter);
}

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


//Event FileReader::event()
//{
//    return __current_event;
//}


int FileReader::current_index()
{
    return __current_index;
}
Long64_t FileReader::num_events()
{
    return __num_events;
}

void FileReader::__load_photon_addresses()
{

    __chain.SetBranchAddress("photon_pt",&(__current_event.photon_pt));
    __chain.SetBranchAddress("photon_e",&(__current_event.photon_e));
    __chain.SetBranchAddress("photon_syst_name",&(__current_event.photon_syst_name));
    __chain.SetBranchAddress("photon_syst_pt",&(__current_event.photon_syst_pt));
    __chain.SetBranchAddress("photon_syst_e",&(__current_event.photon_syst_e));
}



void FileReader::__load_electron_addresses()
{
    __chain.SetBranchAddress("electron_pt",&(__current_event.electron_pt));
    __chain.SetBranchAddress("electron_e",&(__current_event.electron_e));
    __chain.SetBranchAddress("electron_syst_name",&(__current_event.electron_syst_name));
    __chain.SetBranchAddress("electron_syst_pt",&(__current_event.electron_syst_pt));
    __chain.SetBranchAddress("electron_syst_e",&(__current_event.electron_syst_e));
}


void FileReader::__load_cluster_addresses()
{
    __chain.SetBranchAddress("cluster_pt",&(__current_event.cluster_pt));
}


void FileReader::__load_track_addresses()
{
    __chain.SetBranchAddress("track_pt",&(__current_event.track_pt));
    if (__chain.GetListOfBranches()->FindObject("track_type"))
    {
        __chain.SetBranchAddress("track_type",&(__current_event.track_type));
    }
}


void FileReader::__load_truth_particle_addresses()
{
    __chain.SetBranchAddress("mc_pt",&(__current_event.mc_pt));
}


void FileReader::__load_trigger_addresses()
{
    if (__chain.GetListOfBranches()->FindObject("trigger_passed_triggers"))
    {
        __chain.SetBranchAddress("trigger_passed_triggers",&(__current_event.trigger_passed_triggers));
    }
}

void FileReader::__load_photons()
{
    auto length = (*(__current_event.photon_pt)).size();
    double pt, energy;
    int index;
    
    for (size_t i = 0; i < length; ++i)
    {
        
        if ((Event::systematic.empty()) || (Event::systematic=="nominal"))
        {
            pt = (*(__current_event.photon_pt))[i];
            energy = (*(__current_event.photon_e))[i];
        }
        else
        {
            index = -1;
            auto it = find(((*(__current_event.photon_syst_name))[i]).begin(), ((*(__current_event.photon_syst_name))[i]).end(), Event::systematic);
            
            if (it==((*(__current_event.photon_syst_name))[i]).end())
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
        __current_event.photons.emplace_back(Photon(&__chain, i, pt, energy, __current_event.entry_number));
    }
}
void FileReader::__load_electrons()
{
    auto length = (*(__current_event.electron_pt)).size();
    double pt, energy;
    int index;
    
    for (size_t i = 0; i < length; ++i)
    {
        
        if ((Event::systematic.empty()) || (Event::systematic=="nominal"))
        {
            pt = (*(__current_event.electron_pt))[i];
            energy = (*(__current_event.electron_e))[i];
        }
        else
        {
            index = -1;
            auto it = find(((*(__current_event.electron_syst_name))[i]).begin(), ((*(__current_event.electron_syst_name))[i]).end(), Event::systematic);
            
            if (it==((*(__current_event.electron_syst_name))[i]).end())
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
        __current_event.electrons.emplace_back(Electron(&__chain, i, pt, energy, __current_event.entry_number));
    }
}
void FileReader::__load_clusters()
{
    auto length = (*(__current_event.cluster_pt)).size();
    
    for (size_t i = 0; i < length; ++i)
    {
        __current_event.clusters.emplace_back(Cluster(&__chain, i, __current_event.entry_number));
    }
}
void FileReader::__load_tracks()
{
    auto length = (*(__current_event.track_pt)).size();
    
    if (!__current_event.track_type)
    {
        for (size_t i = 0; i < length; ++i)
        {
            __current_event.tracks.emplace_back(Track(&__chain,i,__current_event.entry_number));
        }
    }
    else
    {
        for (size_t i = 0; i < length; ++i)
        {
            if ((*__current_event.track_type)[i]==1)
            {
                __current_event.pixel_tracks.emplace_back(Track(&__chain,i,__current_event.entry_number));
            }
        }
    }

}
void FileReader::__load_truth_particles()
{
    auto length = (*(__current_event.mc_pt)).size();
    for (size_t i = 0; i < length; ++i)
    {
        __current_event.truth_particles.emplace_back(TruthParticle(&__chain, i, __current_event.entry_number));
    }
}
void FileReader::__load_triggers()
{
    if (__current_event.trigger_passed_triggers)
    {
        auto length = (*(__current_event.trigger_passed_triggers)).size();
        for (size_t i = 0; i < length; ++i)
        {
            __current_event.triggers.push_back((*(__current_event.trigger_passed_triggers))[i]);
        }
    }
}

std::vector<TruthParticle> FileReader::find_truth_particles
     (const std::vector<int>&& barcode,
      const std::vector<int>&& parent_barcode,
      const std::vector<int>&& pdg_id,
      int* status_code)
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
            if (status_code && (tp.status_code() != *status_code))
            {
                continue;
            }
            
            results.push_back(tp);
        }
    }
    else
    {
        std::vector<int> *mc_pdg_id = nullptr;
        __chain.SetBranchAddress("mc_pdg_id",&mc_pdg_id);
        
        __chain.GetBranch("mc_pdg_id")->GetEntry(__current_index);
        
        auto length = (*mc_pdg_id).size();
        
        for (decltype(length) i = 0; i < length; ++i)
        {
            if (!barcode.empty())
            {
                std::vector<int> *mc_barcode = nullptr;
                __chain.SetBranchAddress("mc_barcode",&mc_barcode);
                
                __chain.GetBranch("mc_barcode")->GetEntry(__current_index);
                if (find(barcode.begin(),barcode.end(),(*mc_barcode)[i]) == barcode.end())
                {
                    continue;
                }
            }
            if (!parent_barcode.empty())
            {
                std::vector<int> *mc_parent_barcode = nullptr;
                __chain.SetBranchAddress("mc_parent_barcode",&mc_parent_barcode);
                
                __chain.GetBranch("mc_parent_barcode")->GetEntry(__current_index);
                if (find(parent_barcode.begin(),parent_barcode.end(),(*mc_parent_barcode)[i]) == parent_barcode.end())
                {
                    continue;
                }
            }
            if (!pdg_id.empty())
            {
                std::vector<int> *mc_pdg_id = nullptr;
                __chain.SetBranchAddress("mc_pdg_id",&mc_pdg_id);
                
                __chain.GetBranch("mc_pdg_id")->GetEntry(__current_index);
                if (find(pdg_id.begin(),pdg_id.end(),(*mc_pdg_id)[i]) == pdg_id.end())
                {
                    continue;
                }
            }
            if (status_code)
            {
                std::vector<int> *mc_status = nullptr;
                __chain.SetBranchAddress("mc_status",&mc_status);
                
                __chain.GetBranch("mc_status")->GetEntry(__current_index);
                if ((*mc_status)[i] != *status_code)
                {
                    continue;
                }
            }
           
            results.push_back(TruthParticle(&__chain, i, __current_index));
            
        }
    }
    return results;
}





FileReaderRange::FileReaderRange(std::vector<string>&& files,
       const char* tree_name,
       Long64_t num_events, int skip_first_events) :
f(files,tree_name, num_events, skip_first_events)
{}

FileReaderRange::Iterator::Iterator(int i) : data{i} {}

FileReaderRange::Iterator::Iterator(int i, FileReader& f) : data{i}, F{f}
{
    
    //SetAddresses
    
    if (F.__has_event_info_chain)
    {
        F.__event_info_chain.SetBranchAddress("m_RunNumber",&(F.__current_event.m_RunNumber));
        F.__event_info_chain.SetBranchAddress("m_RandomRunNumber",&(F.__current_event.m_RandomRunNumber));
        F.__event_info_chain.SetBranchAddress("ei_event_number",&(F.__current_event.ei_event_number));
    }
    if (Event::cache_truth)
    {
        F.__load_truth_particle_addresses();
    }
    if (Event::load_reco)
    {
        if (Event::load_photons)
        {
            
            F.__load_photon_addresses();
        }
        if (Event::load_electrons)
        {
            F.__load_electron_addresses();
        }
        if (Event::load_clusters)
        {
            F.__load_cluster_addresses();
        }
        if (Event::load_tracks)
        {
            F.__load_track_addresses();
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


FileReaderRange::Iterator& FileReaderRange::Iterator::operator++()
{
    F.__current_event = {}; //reset event
    
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
    

//    if (F.__passes_event_filters())
//    {
////        puts("hello?");
//        //create a new event on each loop -> instead:
//        //really just need to load all particles for a given entry
//        Event Temp(&(F.__chain), &(F.__event_info_chain), data);
//
//
//        F.__current_event = Temp;
////        F.__current_event.entry_number = F.__current_index;
//    }
    return *this;
}

FileReader& FileReaderRange::Iterator::operator*() {return F;}


bool operator!=(const FileReaderRange::Iterator& a, const FileReaderRange::Iterator& b)
{
    return a.data != b.data;
}

FileReaderRange::Iterator FileReaderRange::begin() {return Iterator(f.__skip_first_events, f);}
FileReaderRange::Iterator FileReaderRange::end() {return Iterator(30);}


