#include "event.h"
#include <algorithm>
using std::find;



bool Event::cache_truth = true;
bool Event::load_reco= true;
bool Event::load_photons= true;
bool Event::load_electrons= true;
bool Event::load_clusters= true;
bool Event::load_tracks= true;
bool Event::load_triggers= true;
string Event::systematic= "";

Event::Event() = default;
Event::~Event(){entry.Reset();};

Event::Event(TChain* entry, TChain* event_info_entry)
: run_number{-1}, random_run_number{-1}, event_number{-1}
{
    
//    static int count = 0;
    this->entry.Add(entry);
    if (event_info_entry)
    {
        UInt_t m_RunNumber;
        event_info_entry->SetBranchAddress("m_RunNumber",&m_RunNumber);
//        event_info_entry->GetBranch("m_RunNumber")->GetEntry();
        run_number = m_RunNumber;
//        if (count>4)
//        printf("run_number = %d\n",run_number);
        UInt_t m_RandomRunNumber;
        event_info_entry->SetBranchAddress("m_RandomRunNumber",&m_RandomRunNumber);
//        event_info_entry->GetBranch("m_RandomRunNumber")->GetEntry();
        random_run_number = m_RandomRunNumber;
        
        UInt_t ei_event_number;
        entry->SetBranchAddress("ei_event_number",&ei_event_number);
       //entry->GetBranch("ei_event_number")->GetEntry();
        event_number = ei_event_number;
        
        
    }
    if (cache_truth)
    {
        __load_truth_particles(entry);
    }
    if (load_reco)
    {
        if (load_photons)
        {
            __load_photons(entry);
        }
        if (load_electrons)
        {
            __load_electrons(entry);
        }
        if (load_clusters)
        {
            __load_clusters(entry);
        }
        if (load_tracks)
        {
            __load_tracks(entry);
        }
        if (load_triggers)
        {
            __load_triggers(entry);
        }
    }
//    ++count;
    
}
//TChain entry;
//int run_number;
//int random_run_number;
//int event_number;
//std::vector<Photon> photons;
//std::vector<Electron> electrons;
//std::vector<Cluster> clusters;
//std::vector<Track> tracks;
//std::vector<Track> pixel_tracks;
//std::vector<TruthParticle> truth_particles;
//std::vector<string> triggers; // std::vector<std::vector<string>> triggers; //maybe?
//int muon_spectrometer_num_track_particles;



Event::Event(const Event& other)
{
    entry.Reset();
    entry.Add(const_cast<TChain*>(&(other.entry))); //😬
    run_number = other.run_number;
    random_run_number = other.random_run_number;
    event_number = other.event_number;
    photons = other.photons;
    electrons = other.electrons;
    clusters = other.clusters;
    tracks = other.tracks;
    pixel_tracks = other.pixel_tracks;
    truth_particles = other.truth_particles;
    triggers = other.triggers;
    muon_spectrometer_num_track_particles = other.muon_spectrometer_num_track_particles;
}

Event& Event::operator=(const Event& other)
{
    entry.Reset();
    entry.Add(const_cast<TChain*>(&(other.entry))); //😬
    run_number = other.run_number;
    random_run_number = other.random_run_number;
    event_number = other.event_number;
    photons = other.photons;
    electrons = other.electrons;
    clusters = other.clusters;
    tracks = other.tracks;
    pixel_tracks = other.pixel_tracks;
    truth_particles = other.truth_particles;
    triggers = other.triggers;
    muon_spectrometer_num_track_particles = other.muon_spectrometer_num_track_particles;
    return *this;
}


void Event::__load_photons(TChain* entry)
{
    std::vector<double> *photon_pt= nullptr;
//    entry->SetBranchAddress("photon_pt",&photon_pt);
//    entry->GetBranch("photon_pt")->GetEntry();
//    
//    std::vector<double> *photon_e = nullptr;
//    entry->SetBranchAddress("photon_e",&photon_e);
//    //entry->GetBranch("photon_e")->GetEntry();
//    
//    std::vector<std::vector<string>> *photon_syst_name = nullptr;
//    entry->SetBranchAddress("photon_syst_name",&photon_syst_name);
//    //entry->GetBranch("photon_syst_name")->GetEntry();
//    
//    std::vector<std::vector<double>> *photon_syst_pt = nullptr;
//    entry->SetBranchAddress("photon_syst_pt",&photon_syst_pt);
//    //entry->GetBranch("photon_syst_pt")->GetEntry();
//    
//    std::vector<std::vector<double>> *photon_syst_e = nullptr;
//    entry->SetBranchAddress("photon_syst_e",&photon_syst_e);
//    //entry->GetBranch("photon_syst_e")->GetEntry();
//    
//    auto length = (*photon_pt).size();
//    
//    double pt, energy;
//    int index;
//    
//    for (decltype(length) i = 0; i < length; ++i)
//    {
//        if ((systematic.empty()) || (systematic=="nominal"))
//        {
//            pt = (*photon_pt)[i];
//            energy = (*photon_e)[i];
//        }
//        else
//        {
//            index = -1;
//            auto it = find(((*photon_syst_name)[i]).begin(), ((*photon_syst_name)[i]).end(), systematic);
//            
//            if (it==((*photon_syst_name)[i]).end())
//            {
//                continue;
//            }
//            index = it - ((*photon_syst_name)[i]).begin();
//            pt = (*photon_syst_pt)[i][index];
//            energy = (*photon_syst_e)[i][index];
//        }
//        if (pt < 0)
//        {
//            continue;
//        }
//        photons.emplace_back(Photon(entry, i, pt, energy));
//    }
    
}

void Event::__load_electrons(TChain* entry)
{
//    std::vector<double> *electron_pt= nullptr;
//    entry->SetBranchAddress("electron_pt",&electron_pt);
//   //entry->GetBranch("electron_pt")->GetEntry();
//    
//    std::vector<double> *electron_e = nullptr;
//    entry->SetBranchAddress("electron_e",&electron_e);
//   //entry->GetBranch("electron_e")->GetEntry();
//
//    std::vector<std::vector<string>> *electron_syst_name = nullptr;
//    entry->SetBranchAddress("electron_syst_name",&electron_syst_name);
//   //entry->GetBranch("electron_syst_name")->GetEntry();
//
//    std::vector<std::vector<double>> *electron_syst_pt = nullptr;
//    entry->SetBranchAddress("electron_syst_pt",&electron_syst_pt);
//   //entry->GetBranch("electron_syst_pt")->GetEntry();
//
//    std::vector<std::vector<double>> *electron_syst_e = nullptr;
//    entry->SetBranchAddress("electron_syst_e",&electron_syst_e);
//   //entry->GetBranch("electron_syst_e")->GetEntry();

//    auto length = (*electron_pt).size();
//
//    double pt, energy;
//    int index;
//
//    for (decltype(length) i = 0; i < length; ++i)
//    {
//        if ((systematic.empty()) || (systematic=="nominal"))
//        {
//            pt = (*electron_pt)[i];
//            energy = (*electron_e)[i];
//        }
//        else
//        {
//            index = -1;
//            auto it = find(((*electron_syst_name)[i]).begin(), ((*electron_syst_name)[i]).end(), systematic);
//
//            if (it==((*electron_syst_name)[i]).end())
//            {
//                continue;
//            }
//            index = it - ((*electron_syst_name)[i]).begin();
//            pt = (*electron_syst_pt)[i][index];
//            energy = (*electron_syst_e)[i][index];
//        }
//        if (pt < 0)
//        {
//            continue;
//        }
//        electrons.emplace_back(Electron(entry, i, pt, energy));
//    }

}

void Event::__load_clusters(TChain* entry)
{
    std::vector<double> *cluster_pt= nullptr;
    entry->SetBranchAddress("cluster_pt",&cluster_pt);
    
   //entry->GetBranch("cluster_pt")->GetEntry();
    
    auto length = (*cluster_pt).size();
    
    for (decltype(length) i = 0; i < length; ++i)
    {
        clusters.emplace_back(Cluster(entry, i));
    }
}

void Event::__load_tracks(TChain* entry)
{
    std::vector<double> *track_pt= nullptr;
    entry->SetBranchAddress("track_pt",&track_pt);
    
   //entry->GetBranch("track_pt")->GetEntry();
    
    auto length = (*track_pt).size();
    
    for (decltype(length) i = 0; i < length; ++i)
    {
        
        if (!static_cast<TBranch*>(entry->GetListOfBranches()->FindObject("track_type")))
        {
            tracks.emplace_back(Track(entry,i));
        }
        else
        {
            pixel_tracks.emplace_back(Track(entry,i));
        }
        
    }
}

void Event::__load_triggers(TChain* entry)
{
    puts("__load_triggers");
    if (entry->GetListOfBranches()->FindObject("trigger_passed_triggers"))
    {
        std::vector<string> *trigger_passed_triggers = nullptr;
        
        entry->SetBranchAddress("trigger_passed_triggers",&trigger_passed_triggers);
        
       entry->GetBranch("trigger_passed_triggers")->GetEntry();
        
        auto length = (*trigger_passed_triggers).size();
        printf("Length = %lu\n",length);
        for (decltype(length) i = 0; i < length; ++i)
        {
            triggers.push_back((*trigger_passed_triggers)[i]);
        }
    }
    
}

void Event::__load_truth_particles(TChain* entry)
{
    std::vector<double> *mc_pt= nullptr;
    entry->SetBranchAddress("mc_pt",&mc_pt);
    
   //entry->GetBranch("mc_pt")->GetEntry();
    
    auto length = (*mc_pt).size();
    
    for (decltype(length) i = 0; i < length; ++i)
    {
        truth_particles.emplace_back(TruthParticle(entry, i));
    }
}

std::vector<TruthParticle> Event::find_truth_particles
     (const std::vector<int>&& barcode,
      const std::vector<int>&& parent_barcode,
      const std::vector<int>&& pdg_id,
      int* status_code)
{
    std::vector<TruthParticle> results;
    if (cache_truth)
    {
        for (auto& tp : truth_particles)
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
        entry.SetBranchAddress("mc_pdg_id",&mc_pdg_id);
        
        entry.GetBranch("mc_pdg_id")->GetEntry();
        
        auto length = (*mc_pdg_id).size();
        
        for (decltype(length) i = 0; i < length; ++i)
        {
            if (!barcode.empty())
            {
                std::vector<int> *mc_barcode = nullptr;
                entry.SetBranchAddress("mc_barcode",&mc_barcode);
                
                entry.GetBranch("mc_barcode")->GetEntry();
                if (find(barcode.begin(),barcode.end(),(*mc_barcode)[i]) == barcode.end())
                {
                    continue;
                }
            }
            if (!parent_barcode.empty())
            {
                std::vector<int> *mc_parent_barcode = nullptr;
                entry.SetBranchAddress("mc_parent_barcode",&mc_parent_barcode);
                
                entry.GetBranch("mc_parent_barcode")->GetEntry();
                if (find(parent_barcode.begin(),parent_barcode.end(),(*mc_parent_barcode)[i]) == parent_barcode.end())
                {
                    continue;
                }
            }
            if (!pdg_id.empty())
            {
                std::vector<int> *mc_pdg_id = nullptr;
                entry.SetBranchAddress("mc_pdg_id",&mc_pdg_id);
                
                entry.GetBranch("mc_pdg_id")->GetEntry();
                if (find(pdg_id.begin(),pdg_id.end(),(*mc_pdg_id)[i]) == pdg_id.end())
                {
                    continue;
                }
            }
            if (status_code)
            {
                std::vector<int> *mc_status = nullptr;
                entry.SetBranchAddress("mc_status",&mc_status);
                
                entry.GetBranch("mc_status")->GetEntry();
                if ((*mc_status)[i] != *status_code)
                {
                    continue;
                }
            }
            
            results.push_back(TruthParticle(&entry,i));
            
        }
    }
    return results;
}

//Event::Event() = default;
//Event::~Event() = default;
//Event::Event(const Event& other)
//{
//    
//}
//Event& Event::operator=(const Event&)
//{
//    return *this;
//}
