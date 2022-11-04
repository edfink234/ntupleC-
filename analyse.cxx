/*
 Python version's time: 16.87959684530894 minutes, 13.25841596921285 minutes
 # events for EG_RESOLUTION_ALL__1down = 32227
 # events for nominal = 32215
 */


#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <utility>

#include "TH1F.h"
#include "TFile.h"

//https://root-forum.cern.ch/t/problem-in-accessing-vector-vector-float/27983/2
//#include "LinkDef.h"
#include "plotting.h"
#include "objects.h"
#include "event.h"
#include "filereader.h"
#include "candidate.h"

using Clock = std::chrono::high_resolution_clock;

//debugging: https://drake.mit.edu/profiling.html
constexpr std::array<const char*,4> CUTS =
{"truth","00_no_cuts","01_no_tracks","02_mass_cut"};
std::unordered_map<std::string,Plot> plots;
std::unordered_map<std::string,PlotGroup> plot_groups;

bool photon_selection(Photon& photon)
{
    if (!photon.id_())
    {
        return false;
    }
    
    if (photon.pt() < 2500)
    {
        return false;
    }
    
    if (abs(photon.eta()) > 2.37)
    {
        return false;
    }
    
    if ((1.37 < abs(photon.eta())) && (abs(photon.eta()) < 1.52))
    {
        return false;
    }
    return true;
}

bool truth_photon_selection(TruthParticle& photon)
{
    if (!photon.id_())
    {
        return false;
    }
    
    if (photon.pt() < 2500)
    {
        return false;
    }
    
    if (abs(photon.eta()) > 2.37)
    {
        return false;
    }
    
    if ((1.37 < abs(photon.eta())) && (abs(photon.eta()) < 1.52))
    {
        return false;
    }
    return true;
}


bool track_selection(Track& track)
{
    if (track.pt() < 100)
    {
        return false;
    }
    
    if (abs(track.eta()) > 2.5)
    {
        return false;
    }
    
    return true;
}

void fill_signal_hists(CandidateSet& candidate, std::string cutname, double weight = 1)
{
    // Single photon pT
    switch (candidate.CandidateSetType()) {
        case 'T':
            plot_groups.at(cutname+"/reco/photon_pt").fill(candidate.particle_a.pt()/1e3,weight);
            plot_groups.at(cutname+"/reco/photon_pt").fill(candidate.particle_b.pt()/1e3,weight);
            break;
        case 'P':
            plot_groups.at(cutname+"/reco/photon_pt").fill(candidate.particle_a_photon.pt()/1e3,weight);
            plot_groups.at(cutname+"/reco/photon_pt").fill(candidate.particle_b_photon.pt()/1e3,weight);
            break;
    }
    // Invariant diphoton mass
    plot_groups.at(cutname+"/candidate/mass").fill(candidate.four_momentum.M()/1e3,weight);
}

void fill_tracking_hists(Event& event, std::string cutname, double weight = 1)
{
    // Fill the tracking histograms.
    plots.at(cutname+"/tracking/num_tracks").fill(event.tracks.size(),weight);
    
    for (auto& track: event.tracks)
    {
        plots.at(cutname+"/tracking/track_pt").fill(track.pt()/1e3,weight);
        plots.at(cutname+"/tracking/track_eta").fill(track.eta(),weight);
    }
}

void run_analysis(const std::vector<std::string>& input_filenames, std::string systematic = "nominal", bool mc = false, TFile* output_file = nullptr)
{
    plots.clear();
    plot_groups.clear();
    
    std::cout <<"systematic = " << systematic << '\n';
    plots.emplace("cutflow",Plot(systematic+"/cutflow", "cutflow", 10, 0, 10));

    for (const std::string& cutname: CUTS)
    {
        plot_groups.emplace(
        cutname+"/reco/photon_pt",
        PlotGroup({
            Plot(systematic+"/"+cutname
            +"/reco/photon/pt",
            "pT / GeV", 20, 0, 25),
            Plot(systematic+"/"+cutname
            +"/reco/photon/pt_tight",
            "pT / GeV", 100, 0, 10)
        }));
        
        plot_groups.emplace(
        cutname+"/candidate/mass",
        PlotGroup({
            Plot(systematic+"/"+cutname
            +"/candidate/mass",
            "candidate mass", 250, 0, 50),
            Plot(systematic+"/"+cutname
            +"/candidate/mass_large",
            "candidate mass", 200, 0, 200),
            Plot(systematic+"/"+cutname
            +"/candidate/mass_large_fine",
            "candidate mass", 600, 0, 200)
        }));

        plots.emplace(cutname+"/tracking/num_tracks",
            Plot(systematic+"/"+cutname
            +"/tracking/num_tracks",
            "tracking num_tracks", 20, 0, 20));
        
        plots.emplace(cutname+"/tracking/track_pt",
            Plot(systematic+"/"+cutname
            +"/tracking/track_pt",
            "tracking track_pt", 140, 0, 7));
        
        plots.emplace(cutname+"/tracking/track_eta",
            Plot(systematic+"/"+cutname
            +"/tracking/track_eta",
            "tracking track_eta", 50, -2.5, 2.5));
    }
    
    FileReaderRange reader(input_filenames);
    Event::systematic = systematic;
    
    int num_passed_events = 0;
    int weight = 1;
    
    for (auto &&f: reader)
    {
        std::cout << "entry_number " << f.__current_event.entry_number  << '\n';
        
        if ((!mc) && (find(f.__current_event.triggers.begin(), f.__current_event.triggers.end(), "HLT_hi_upc_FgapAC3_hi_gg_upc_L1TAU1_TE4_VTE200") == f.__current_event.triggers.end()))
        {
            continue;
        }
     
        std::vector<Photon> photons;
        std::copy_if (f.__current_event.photons.begin(), f.__current_event.photons.end(), std::back_inserter(photons), photon_selection );
        f.__current_event.photons = move(photons);
        
        std::vector<Track> tracks;
        std::copy_if (f.__current_event.tracks.begin(), f.__current_event.tracks.end(), std::back_inserter(tracks), track_selection );
        f.__current_event.tracks = move(tracks);

        std::vector<TruthParticle>&& pre_truth_photons = f.find_truth_particles({},{},{Photon::PDG_ID});
        std::vector<TruthParticle> truth_photons;
        std::copy_if (pre_truth_photons.begin(), pre_truth_photons.end(), std::back_inserter(truth_photons), truth_photon_selection );

        if (truth_photons.size() == 2)
        {
            CandidateSet truth_candidate(std::make_pair(truth_photons[0],truth_photons[1]));
            if (truth_candidate.four_momentum.M() > 5000)
            {
                fill_signal_hists(truth_candidate,"truth",1);
            }
        }

        if (f.__current_event.photons.size()==2)
        {
            CandidateSet candidate(std::make_pair(f.__current_event.photons[0],f.__current_event.photons[1]));
            
            fill_signal_hists(candidate,"00_no_cuts",weight);
            fill_tracking_hists(f.__current_event,"00_no_cuts",weight);

            if (f.__current_event.tracks.empty() && candidate.four_momentum.M() > 5000)
            {
                fill_signal_hists(candidate,"02_mass_cut",weight);
                fill_tracking_hists(f.__current_event,"02_mass_cut",weight);
                num_passed_events += 1;
            }

            if (f.__current_event.tracks.empty())
            {
                fill_signal_hists(candidate,"01_no_tracks",weight);
                fill_tracking_hists(f.__current_event,"01_no_tracks",weight);
            }
        }
    }
    
    for (auto& plot_group: plot_groups)
    {
        plot_group.second.save(output_file);
    }
    
    for (auto& plot: plots)
    {
        plot.second.save(output_file);
    }

    std::cout << "# events for " << systematic << " = " << num_passed_events
    << '\n';
}

void analyse_haa()
{
    auto start_time = Clock::now();
    std::cout << "Run over MC\n";

    std::vector<std::string> input_filenames = {"../user.kschmied.28655874._000025.LGNTuple.root"};
    
//    std::vector<std::string> input_filenames = {"/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000025.LGNTuple.root"};
    
//    std::vector<std::string> input_filenames =
//    {
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000002.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000001.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000005.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000006.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000007.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000004.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000003.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000010.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000009.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000008.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000011.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000012.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000015.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000013.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000017.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000014.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000016.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000019.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000018.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000021.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000023.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000020.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000024.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000022.LGNTuple.root",
//        "/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000025.LGNTuple.root"
//    };
//    const char* output_filename = "example_mc_haa_out_test1cpp.root";
    const char* output_filename = "example_mc_haa_out_cppreg.root";
    TFile* output_file = TFile::Open(output_filename, "RECREATE");
    if (!output_file) {
        std::cout << "Error opening file\n";
       exit(-1);
    }
    
    std::vector<std::string> systematics = {"nominal", "EG_RESOLUTION_ALL__1down"};
    
    for (std::string& systematic: systematics)
    {
        run_analysis(input_filenames, systematic.c_str(), true, output_file);
    }
    

    
    output_file->Close();
    auto end_time = Clock::now();
    std::cout << "Time difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
}

int main()
{
    analyse_haa();
}

