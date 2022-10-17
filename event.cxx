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

Event::Event(TChain* entry, TChain* event_info_entry)
: run_number{-1}, random_run_number{-1}, event_number{-1}
{
    this->entry.Add(entry);
    if (event_info_entry)
    {
        UInt_t m_RunNumber;
        event_info_entry->SetBranchAddress("m_RunNumber",&m_RunNumber);
        event_info_entry->GetBranch("m_RunNumber")->GetEntry();
        run_number = m_RunNumber;
        
        UInt_t m_RandomRunNumber;
        event_info_entry->SetBranchAddress("m_RandomRunNumber",&m_RandomRunNumber);
        event_info_entry->GetBranch("m_RandomRunNumber")->GetEntry();
        random_run_number = m_RandomRunNumber;
        
        UInt_t ei_event_number;
        entry->SetBranchAddress("ei_event_number",&ei_event_number);
        entry->GetBranch("ei_event_number")->GetEntry();
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
    
}

void Event::__load_truth_particles(TChain* entry)
{
    vector<double> *mc_pt= nullptr;
    entry->SetBranchAddress("mc_pt",&mc_pt);
    
    entry->GetBranch("mc_pt")->GetEntry();
    
    auto length = (*mc_pt).size();
    
    for (decltype(length) i = 0; i < length; ++i)
    {
        truth_particles.emplace_back(TruthParticle(entry, i));
    }
}

void Event::__load_photons(TChain* entry)
{
    vector<double> *photon_pt= nullptr;
    entry->SetBranchAddress("photon_pt",&photon_pt);
    entry->GetBranch("photon_pt")->GetEntry();
    
    vector<double> *photon_e = nullptr;
    entry->SetBranchAddress("photon_e",&photon_e);
    entry->GetBranch("photon_e")->GetEntry();
    
    vector<vector<string>> *photon_syst_name = nullptr;
    entry->SetBranchAddress("photon_syst_name",&photon_syst_name);
    entry->GetBranch("photon_syst_name")->GetEntry();
    
    vector<vector<double>> *photon_syst_pt = nullptr;
    entry->SetBranchAddress("photon_syst_pt",&photon_syst_pt);
    entry->GetBranch("photon_syst_pt")->GetEntry();
    
    vector<vector<double>> *photon_syst_e = nullptr;
    entry->SetBranchAddress("photon_syst_e",&photon_syst_e);
    entry->GetBranch("photon_syst_e")->GetEntry();
    
    auto length = (*photon_pt).size();
    
    double pt, energy;
    int index;
    
    for (decltype(length) i = 0; i < length; ++i)
    {
        if ((systematic.empty()) || (systematic=="nominal"))
        {
            pt = (*photon_pt)[i];
            energy = (*photon_e)[i];
        }
        else
        {
            index = -1;
            auto it = find(((*photon_syst_name)[i]).begin(), ((*photon_syst_name)[i]).end(), systematic);
            
            if (it==((*photon_syst_name)[i]).end())
            {
                continue;
            }
            index = it - ((*photon_syst_name)[i]).begin();
            pt = (*photon_syst_pt)[i][index];
            energy = (*photon_syst_e)[i][index];
        }
        if (pt < 0)
        {
            continue;
        }
        photons.emplace_back(Photon(entry, i, pt, energy));
    }
    
//    entry->GetBranch("mc_pt")->GetEntry();
//    R__ASSERT((*mc_pt).size() > static_cast<size_t>(_index));
//    return (*mc_pt)[_index];
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
