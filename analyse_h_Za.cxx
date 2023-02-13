/*
 Preselection for H->Za analysis.
 */

#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <iostream>

#include "TH1F.h"
#include "TFile.h"

//https://root-forum.cern.ch/t/problem-in-accessing-vector-vector-float/27983/2
#include "plotting.h"
#include "objects.h"
#include "event.h"
#include "filereader.h"

using Clock = std::chrono::high_resolution_clock;

//debugging: https://drake.mit.edu/profiling.html
constexpr std::array<const char*,3> CUTS = {"truth","reco","reco_2y"};
std::unordered_map<std::string,Plot> plots;
std::unordered_map<std::string,PlotGroup> plot_groups;

//ΔR between two leptons
double dR(TruthParticle& lepton1, TruthParticle& lepton2)
{
    return lepton1.delta_r(lepton2);
}

//preselection cuts
bool lepton_selection(std::vector<TruthParticle>& leptons)
{
    return
    (
     (leptons.size() == 2)
    &&
     (leptons[0].delta_r(leptons[1]) > 0.01)
    &&
    (leptons[0].charge() == -1*leptons[1].charge())
    &&
    ((((leptons[0].e()+leptons[1].e())*(leptons[0].e()+leptons[1].e())) -
    ((leptons[0].pt()+leptons[1].pt())*(leptons[0].pt()+leptons[1].pt()))) >= 6561e6)
    &&
    ((((leptons[0].e()+leptons[1].e())*(leptons[0].e()+leptons[1].e())) -
    ((leptons[0].pt()+leptons[1].pt())*(leptons[0].pt()+leptons[1].pt()))) <= 10201e6)
    &&
    ((leptons[0].pt() > 20e3 && leptons[1].pt() > 27e3)
     ||
    (leptons[1].pt() > 20e3 && leptons[0].pt() > 27e3))
     );
}

//photon selection cuts, not really part of preselection in the paper afaik
bool photon_selection(Photon& photon)
{
    return ((photon.id_())
            && (photon.pt() > 10000)
    && (abs(photon.eta()) < 2.37)
            && (!((1.37 < abs(photon.eta())) && (1.52 > abs(photon.eta()))))
            );
}

//track selection cuts, not really part of preselection in the paper afaik
bool track_selection(Track& track)
{
    return ((track.pt() > 1000) && (abs(track.eta()) < 2.5));
}

//function to fill histograms with lepton and photon pt and eta
void fill_signal_hists(std::vector<TruthParticle>& particles, std::string cutname, double weight = 1, std::string particleType = "photons")
{
    if (particleType == "leptons")
    {
        for (auto& lepton: particles)
        {
            plot_groups[cutname+"/leptons/pt"].fill(lepton.pt()/1e3,weight);
            plots[cutname+"/leptons/eta"].fill(lepton.eta(),weight);
        }
        return;
    }

    for (auto& photon: particles)
    {
        plot_groups.at(cutname+"/photons/pt").fill(photon.pt()/1e3,weight);
        plots.at(cutname+"/photons/eta").fill(photon.eta(),weight);
    }
}
/*
function overload to fill histograms with lepton and photon pt and eta for
a vector of Photon objects instead of Truthparticles
 */
void fill_signal_hists(std::vector<Photon>& particles, std::string cutname, double weight = 1, std::string particleType = "photons")
{
    if (particleType == "leptons")
    {
        for (auto& lepton: particles)
        {
            plot_groups[cutname+"/leptons/pt"].fill(lepton.pt()/1e3,weight);
            plots[cutname+"/leptons/eta"].fill(lepton.eta(),weight);
        }
        return;
    }

    for (auto& photon: particles)
    {
        plot_groups.at(cutname+"/photons/pt").fill(photon.pt()/1e3,weight);
        plots.at(cutname+"/photons/eta").fill(photon.eta(),weight);
    }
}

//driver function for the analysis
void run_analysis(std::vector<std::string>& input_filenames, std::string systematic = "nominal", bool mc = false)
{
//TODO: Make one container of a union or variant of Plot/PlotGroup instead of 2 containers, e.g.
//variant: https://en.cppreference.com/w/cpp/utility/variant
//union: https://www.sololearn.com/compiler-playground/cop9eIyns3c3
    
    /*
    plots we want to create for this analysis. The constructors
    easily creates a directory structure to store the files in a
    single file!
     */
    plots.emplace("cutflow",Plot(systematic+"/cutflow", "cutflow", 10, 0, 10));
    for (const std::string& cutname: CUTS)
    {
        plots.emplace(cutname+"/photons/eta",
            Plot(systematic+"/"+cutname
            +"/photons/eta",
            "eta", 40, -6, 6));

        plots.emplace(cutname+"/leptons/eta",
            Plot(systematic+"/"+cutname
            +"/leptons/eta",
            "eta", 40, -6, 6));

        plot_groups.emplace(
        cutname+"/photons/pt",
        PlotGroup({
            Plot(systematic+"/"+cutname
            +"/photons/pt",
            "pT / GeV", 20, 0, 200),
            Plot(systematic+"/"+cutname
            +"/photons/pt_tight",
            "pT / GeV", 100, 0, 100)
        }));

        plot_groups.emplace(
        cutname+"/leptons/pt",
        PlotGroup({
            Plot(systematic+"/"+cutname
            +"/leptons/pt",
            "pT / GeV", 20, 0, 200),
            Plot(systematic+"/"+cutname
            +"/leptons/pt_tight",
            "pT / GeV", 100, 0, 100)
        }));
    }
    
    FileReaderRange reader(input_filenames); //instantiate filereader
    Event::systematic = systematic; //run analysis for a given systematic
    Event::load_tracks = true; //it's always true anyway...
    int num_passed_events = 0; //number of events that passed the photon selection
    
    for (auto &&f: reader) // read over files
    {
        std::cout << "entry_number " << f.__current_event.entry_number  << '\n';
        
        int weight = 1;
        
        if (mc)
        {
            std::vector<TruthParticle>&& truth_higgs = f.find_truth_particles({},{},{35});
            if (!(truth_higgs.empty()))
            {
                std::vector<TruthParticle>&& truth_axions = f.find_truth_particles({},{truth_higgs[0].barcode()}, {36});
            }
            int temp = 1;
            //get truth particles with pdg_id = 22, i.e., photons, in current event
            std::vector<TruthParticle>&& truth_photons = f.find_truth_particles({},{},{22},&temp);
            //fill histograms with properties of these truth photons
            fill_signal_hists(truth_photons,"truth");
            
            //and the same for truth leptons and anti-leptons
            std::vector<TruthParticle>&& truth_leptons = f.find_truth_particles({},{},{11, 12, 13, 14, 15, 16, 17, 18},&temp);
            
            //check if lepton selection is passed...
            if (lepton_selection(truth_leptons))
            {
                //then fill corresponding histos with these properties.
                fill_signal_hists(truth_leptons, "truth",1,"leptons");
            }
        }
        //container to hold reco-photons that pass photon selection
        std::vector<Photon> photons;
        std::copy_if (f.__current_event.photons.begin(), f.__current_event.photons.end(), std::back_inserter(photons), photon_selection );
        //then fill the corresponding histograms with these selected photons
        fill_signal_hists(photons, "reco", weight);
        //if there are exactly two that passed the cuts
        if (photons.size() == 2)
        {
            //then fill the corresponding histograms with these selected photons
            fill_signal_hists(photons,"reco_2y",weight);
            num_passed_events++; //number of events that passed the photon selection
        }
    }

    std::cout << "# events for " << systematic << " = " << num_passed_events
    << '\n';
}

void analyse_haa()
{
    auto start_time = Clock::now();
    std::cout << "Run over MC\n";

    std::vector<std::string> input_filenames =
    {"../user.kschmied.28655874._000025.LGNTuple.root"};
    
    const char* output_filename = "example_mc_haa_out_testcpp.root";
    
    TFile* output_file = TFile::Open(output_filename, "RECREATE");
    
    if (!output_file) {
        std::cout << "Error opening file\n";
       exit(-1);
    }
    std::vector<std::string> systematics = {"nominal"};
    
    for (std::string& systematic: systematics)
    {
        run_analysis(input_filenames, systematic.c_str(), true);
    }
    //Save plots to output_file
    for (auto& plot: plots)
    {
        plot.second.save(output_file);
    }
    //Save plots in plot_group to output_file
    for (auto& plot_group: plot_groups)
    {
        plot_group.second.save(output_file);
    }
    //close output_file, good practice
    output_file->Close();
    
    auto end_time = Clock::now();
    
    std::cout << "Time difference: "
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()/1e9 << " seconds" << std::endl;
}

int main()
{
    analyse_haa();
}
