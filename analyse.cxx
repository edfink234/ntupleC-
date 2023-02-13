/*
 Example analysis file for benchmarking different frameworks and for explaining how
 the framework works
 
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

#include "plotting.h"
#include "objects.h"
#include "event.h"
#include "filereader.h"
#include "candidate.h"

using Clock = std::chrono::high_resolution_clock;

constexpr std::array<const char*,4> CUTS =
{"truth","00_no_cuts","01_no_tracks","02_mass_cut"};
std::unordered_map<std::string,Plot> plots;
std::unordered_map<std::string,PlotGroup> plot_groups;

//reco-photon selection
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

//truth-photon selection
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

//track selection
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
/*
filters a vector based on a boolean function, if the function is false for a
particular element in the vector the original element is erased, not a copy
 */
template <typename T>
void filter(std::vector<T>& vec, bool (*func) (T&))
{
    for (auto i = vec.begin(); i != vec.end(); ++i)
    {
        if (!func(*i))
        {
            vec.erase(i);
            i--;
        }
    }
}

//fill signal histograms with an event
template <typename T>
void fill_signal_hists(CandidateSet<T>& candidate, std::string cutname, double weight = 1)
{
    plot_groups.at(cutname+"/reco/photon_pt").fill(candidate.particle_a.pt()/1e3,weight);
    plot_groups.at(cutname+"/reco/photon_pt").fill(candidate.particle_b.pt()/1e3,weight);

    // Invariant diphoton mass
    plot_groups.at(cutname+"/candidate/mass").fill(candidate.four_momentum.M()/1e3,weight);
}

//fill tracking histograms with an event
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

//driver function for the analysis
void run_analysis(const std::vector<std::string>& input_filenames, std::string systematic = "nominal", bool mc = false, TFile* output_file = nullptr)
{
    plots.clear();
    plot_groups.clear();
    
    std::cout <<"systematic = " << systematic << '\n';
    
    /*
    plots we want to create for this analysis. The constructors
    easily creates a directory structure to store the files in a
    single file!
     */
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
    
    FileReaderRange reader(input_filenames); //instantiate filereader
    Event::systematic = systematic; //run analysis for a given systematic
    
    int num_passed_events = 0;  //number of events that passed the selection below
    int weight = 1; //how to weight the histogram entries, here we set it to 1 (no weighting)
    
    for (auto &&f: reader) //syntactic sugar to easily read over files
    {
        std::cout << "entry_number " << f.__current_event.entry_number  << '\n';
        
        //if the event passes this trigger criteria, keep, else, reject
        if ((!mc) && (find(f.__current_event.triggers.begin(), f.__current_event.triggers.end(), "HLT_hi_upc_FgapAC3_hi_gg_upc_L1TAU1_TE4_VTE200") == f.__current_event.triggers.end()))
        {
            continue;
        }
     
        //filter reco-photons in event based on photon_selection, see function definition
        filter<Photon>(f.__current_event.photons,&photon_selection);
        
        //filter reco-tracks in event based on track_selection, see function definition
        filter<Track>(f.__current_event.tracks,&track_selection);
        
        //get truth particles with pdg_id = Photon::PDG_ID (22), i.e., photons, in current event
        std::vector<TruthParticle>&& truth_photons = f.find_truth_particles({},{},{Photon::PDG_ID});
        
        //filter truth photons in event based on truth_photon_selection, see function definition
        filter<TruthParticle>(truth_photons,&truth_photon_selection);
        
        //check if 2 truth photons remain after truth_photon_selection
        if (truth_photons.size() == 2)
        {
            //if so, make a CandidateSet out of them...
            CandidateSet<TruthParticle>
            truth_candidate(std::make_pair(truth_photons[0],truth_photons[1]));
            //to easily calculate the invariant mass and check if its greater than 5 GeV
            if (truth_candidate.four_momentum.M() > 5000)
            {
                //if so fill corresponding histos with these properties.
                fill_signal_hists(truth_candidate,"truth",1);
            }
        }

        //check if 2 reco photons remain after photon_selection
        if (f.__current_event.photons.size()==2)
        {
            //if so, make a CandidateSet out of them...
            CandidateSet<Photon>
            candidate(std::make_pair(f.__current_event.photons[0],f.__current_event.photons[1]));
            
            //fill corresponding histos with these properties.
            fill_signal_hists(candidate,"00_no_cuts",weight);
            fill_tracking_hists(f.__current_event,"00_no_cuts",weight);

            //check if invariant mass is greater than 5 GeV AND no tracks remain in event
            if (f.__current_event.tracks.empty() && candidate.four_momentum.M() > 5000)
            {
                //fill corresponding histos with these properties.
                fill_signal_hists(candidate,"02_mass_cut",weight);
                fill_tracking_hists(f.__current_event,"02_mass_cut",weight);
                num_passed_events += 1; //event passed :)
            }
            
            //check if no tracks remain in event
            if (f.__current_event.tracks.empty())
            {
                //fill corresponding histos with these properties.
                fill_signal_hists(candidate,"01_no_tracks",weight);
                fill_tracking_hists(f.__current_event,"01_no_tracks",weight);
            }
        }
    }
    
    //Save plots in plot_group to output_file
    for (auto& plot_group: plot_groups)
    {
        plot_group.second.save(output_file);
    }
    
    //Save plots to output_file
    for (auto& plot: plots)
    {
        plot.second.save(output_file);
    }

    std::cout << "# events for " << systematic << " = " << num_passed_events
    << '\n';
}

void analyse()
{
    auto start_time = Clock::now();
    std::cout << "Run over MC\n";

    std::vector<std::string> input_filenames = {"../user.kschmied.28655874._000025.LGNTuple.root"};
 
    const char* output_filename = "example_mc_haa_out_test1cppdummy.root";

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
    //close output_file, good practice
    output_file->Close();
    auto end_time = Clock::now();
    std::cout << "Time difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
}

int main()
{
    analyse();
}

