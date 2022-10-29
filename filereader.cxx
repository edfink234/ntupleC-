#include "filereader.h"

FileReader::FileReader(std::vector<std::string>& files, const char* tree_name, Long64_t num_events, int skip_first_events) :

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
    __chain.SetBranchStatus("photon_pt",1);
    __chain.SetBranchStatus("photon_e",1);
    __chain.SetBranchStatus("photon_syst_name",1);
    __chain.SetBranchStatus("photon_syst_pt",1);
    __chain.SetBranchStatus("photon_syst_e",1);
    
    __chain.SetBranchAddress("photon_pt",&(__current_event.photon_pt));
    __chain.SetBranchAddress("photon_e",&(__current_event.photon_e));
    __chain.SetBranchAddress("photon_syst_name",&(__current_event.photon_syst_name));
    __chain.SetBranchAddress("photon_syst_pt",&(__current_event.photon_syst_pt));
    __chain.SetBranchAddress("photon_syst_e",&(__current_event.photon_syst_e));
}

void FileReader::__load_electron_addresses()
{
    __chain.SetBranchStatus("electron_pt",1);
    __chain.SetBranchStatus("electron_e",1);
    __chain.SetBranchStatus("electron_syst_name",1);
    __chain.SetBranchStatus("electron_syst_pt",1);
    __chain.SetBranchStatus("electron_syst_e",1);
    
    __chain.SetBranchAddress("electron_pt",&(__current_event.electron_pt));
    __chain.SetBranchAddress("electron_e",&(__current_event.electron_e));
    __chain.SetBranchAddress("electron_syst_name",&(__current_event.electron_syst_name));
    __chain.SetBranchAddress("electron_syst_pt",&(__current_event.electron_syst_pt));
    __chain.SetBranchAddress("electron_syst_e",&(__current_event.electron_syst_e));
}

void FileReader::__load_cluster_addresses()
{
    __chain.SetBranchStatus("cluster_pt",1);
    __chain.SetBranchAddress("cluster_pt",&(__current_event.cluster_pt));
}

void FileReader::__load_track_addresses()
{
    __chain.SetBranchStatus("track_pt",1);
    __chain.SetBranchAddress("track_pt",&(__current_event.track_pt));
    if (__chain.GetListOfBranches()->FindObject("track_type"))
    {
        __chain.SetBranchStatus("track_type",1);
        __chain.SetBranchAddress("track_type",&(__current_event.track_type));
    }
}

void FileReader::__load_truth_particle_addresses(bool cache_truth)
{
    if (cache_truth)
    {
        __chain.SetBranchStatus("mc_pt",1);
        __chain.SetBranchAddress("mc_pt",&(__current_event.mc_pt));
    }
    else
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

void FileReader::__load_trigger_addresses()
{
    if (__chain.GetListOfBranches()->FindObject("trigger_passed_triggers"))
    {
        __chain.SetBranchStatus("trigger_passed_triggers",1);
        __chain.SetBranchAddress("trigger_passed_triggers",&(__current_event.trigger_passed_triggers));
    }
}

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

        __current_event.photons.emplace_back(Photon(i, pt, energy, __current_event.entry_number));
    }
}

void FileReader::__load_electrons()
{
    auto length = (*(Electron::electron_pt)).size();
    double pt, energy;
    int index;
    
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
        __current_event.electrons.emplace_back(Electron(i, pt, energy, __current_event.entry_number));
    }
}

void FileReader::__load_clusters()
{
    auto length = (*(Cluster::cluster_pt)).size();
    
    for (size_t i = 0; i < length; i++)
    {
        __current_event.clusters.emplace_back(Cluster(i, __current_event.entry_number));
    }
}

void FileReader::__load_tracks()
{
    auto length = (*(Track::track_pt)).size();
    
    if (!__current_event.track_type)
    {
        for (size_t i = 0; i < length; i++)
        {
            __current_event.tracks.emplace_back(Track(i,__current_event.entry_number));
        }
    }
    else
    {
        for (size_t i = 0; i < length; i++)
        {
            if ((*__current_event.track_type)[i] == 1)
            {
                __current_event.pixel_tracks.emplace_back(Track(i,__current_event.entry_number));
            }
        }
    }
}

void FileReader::__load_truth_particles()
{
    auto length = (*(TruthParticle::mc_pt)).size();
    for (size_t i = 0; i < length; i++)
    {

        __current_event.truth_particles.emplace_back(TruthParticle(i, __current_event.entry_number));
    }
}

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
        auto length = (*__current_event.mc_pdg_id).size();

        for (decltype(length) i = 0; i < length; i++)
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
            results.push_back(TruthParticle(i, __current_index));
        }
    }
    return results;
}

FileReaderRange::FileReaderRange(std::vector<std::string>&& files,
       const char* tree_name,
       Long64_t num_events, int skip_first_events) :
f(files,tree_name, num_events, skip_first_events)
{
    
}

FileReaderRange::Iterator::Iterator(int i) : data{i} {}

FileReaderRange::Iterator::Iterator(int i, FileReader& f) : data{i}, F{f}
{
    F.__current_index = F.__current_event.entry_number = data;
    F.__chain.SetBranchStatus("*",0);
    //SetAddresses

    if (F.__has_event_info_chain)
    {
        F.__chain.SetBranchStatus("m_RunNumber",1);
        F.__chain.SetBranchStatus("m_RandomRunNumber",1);
        F.__chain.SetBranchStatus("ei_event_number",1);
//        
        
        F.__event_info_chain.SetBranchAddress("m_RunNumber",&(F.__current_event.m_RunNumber));
        
        F.__event_info_chain.SetBranchAddress("m_RandomRunNumber",&(F.__current_event.m_RandomRunNumber));
        
        F.__event_info_chain.SetBranchAddress("ei_event_number",&(F.__current_event.ei_event_number));
        
    }

    F.__load_truth_particle_addresses(Event::cache_truth);

    if (Event::load_reco)
    {
        if (Event::load_photons)
        {
            
            F.__load_photon_addresses();
            
        }
        if (Event::load_electrons)
        {
            F.__load_electron_addresses();
            Photon::SetPhoton(&F.__chain);
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
        if (Event::load_clusters)
        {
            F.__load_cluster_addresses();
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
    F.__current_event.truth_particles.clear();
    F.__current_event.photons.clear();
    F.__current_event.triggers.clear();
    
//    F.__current_event = {}; //reset event
//    printf("data = %d\n",data);
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
    
    return *this;
}

FileReader& FileReaderRange::Iterator::operator*()
{
    return F;
}

bool operator!=(const FileReaderRange::Iterator& a, const FileReaderRange::Iterator& b)
{
    return a.data != b.data;
}

FileReaderRange::Iterator FileReaderRange::begin()
{
    return Iterator(f.__skip_first_events, f);
}

FileReaderRange::Iterator FileReaderRange::end()
{
    return Iterator(f.__num_events);
//    return Iterator(100);
}


