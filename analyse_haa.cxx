#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <chrono>
//https://root-forum.cern.ch/t/problem-in-accessing-vector-vector-float/27983/2
#include "LinkDef.h"
#include "TFile.h"
#include "plotting.h"
#include "objects.h"
#include "event.h"
#include "filereader.h"
#include "TH1F.h"
#include "TInterpreter.h"
#include "RtypesCore.h"
using Clock = std::chrono::high_resolution_clock;


//debugging: https://drake.mit.edu/profiling.html
const std::vector<string> CUTS = {"truth","reco","reco_2y"};
std::unordered_map<string,Plot> plots;
std::unordered_map<string,PlotGroup> plot_groups;


double dR(TruthParticle& lepton1, TruthParticle& lepton2)
{
    return
    pow(
        ((lepton1.eta() - lepton2.eta())*(lepton1.eta() - lepton2.eta())
    
         + (lepton1.phi() - lepton2.phi())*(lepton1.phi() - lepton2.phi())),0.5);
}

bool lepton_selection(vector<TruthParticle>& leptons)
{
    return
    (
     (leptons.size() == 2)
    &&
    (dR(leptons[0],leptons[1]) > 0.01)
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

bool photon_selection(Photon& photon)
{
//    cout << "photon_selection here\n";
    return ((photon.id_()) && (photon.pt() > 10000)
    && (abs(photon.eta()) < 2.37)
            && (!((1.37 < abs(photon.eta())) && (1.52 > abs(photon.eta())))));
}

bool track_selection(Track& track)
{
    return ((track.pt() > 1000) && (abs(track.eta()) < 2.5));
}

/*
 I realize I could just use templates here for fill_signal_hists,
 but it makes it ~10^-1 seconds slower...
 */

void fill_signal_hists(std::vector<TruthParticle>& particles, string cutname, double weight = 1, string particleType = "photons")
{
    
    if (particleType=="leptons")
    {
        for (auto lepton: particles)
        {
            plot_groups[cutname+string("/leptons/pt")].fill(lepton.pt()/1e3,weight);
            plots[cutname+string("/leptons/eta")].fill(lepton.eta(),weight);
        }
        return;
    }

    
    for (auto photon: particles)
    {

        plot_groups.at(cutname+string("/photons/pt")).fill(photon.pt()/1e3,weight);
        plots.at(cutname+string("/photons/eta")).fill(photon.eta(),weight);

    }
}

void fill_signal_hists(std::vector<Photon>& particles, string cutname, double weight = 1, string particleType = "photons")
{

    if (particleType=="leptons")
    {
        for (auto lepton: particles)
        {
            plot_groups[cutname+string("/leptons/pt")].fill(lepton.pt()/1e3,weight);
            plots[cutname+string("/leptons/eta")].fill(lepton.eta(),weight);
        }
        return;
    }

    for (auto photon: particles)
    {
        plot_groups.at(cutname+string("/photons/pt")).fill(photon.pt()/1e3,weight);
        plots.at(cutname+string("/photons/eta")).fill(photon.eta(),weight);
    }
}


void run_analysis(string& input_filename, string systematic = "nominal", bool mc = false)
{
//https://en.cppreference.com/w/cpp/utility/variant
    //or union
//https://www.sololearn.com/compiler-playground/cop9eIyns3c3
    

    
    plots.emplace("cutflow",Plot(systematic+string("/cutflow"), string("cutflow"), 10, 0, 10));
    

    for (const string& cutname: CUTS)
    {
        plots.emplace(cutname+string("/photons/eta"),
            Plot(systematic+string(1,'/')+cutname
            +string("/photons/eta"),
            string("eta"), 40, -6, 6));

        plots.emplace(cutname+string("/leptons/eta"),
            Plot(systematic+string(1,'/')+cutname
            +string("/leptons/eta"),
            string("eta"), 40, -6, 6));

        plot_groups.emplace(
        cutname+string("/photons/pt"),
        PlotGroup({
            Plot(systematic+string(1,'/')+cutname
            +string("/photons/pt"),
            string("pT / GeV"), 20, 0, 200),
            Plot(systematic+string(1,'/')+cutname
            +string("/photons/pt_tight"),
            string("pT / GeV"), 100, 0, 100)
        }));

        plot_groups.emplace(
        cutname+string("/leptons/pt"),
        PlotGroup({
            Plot(systematic+string(1,'/')+cutname
            +string("/leptons/pt"),
            string("pT / GeV"), 20, 0, 200),
            Plot(systematic+string(1,'/')+cutname
            +string("/leptons/pt_tight"),
            string("pT / GeV"), 100, 0, 100)
        }));

    }
    FileReaderRange reader({input_filename});
    Event::systematic = systematic;
    Event::load_tracks = true;
    int num_passed_events = 0;
    

    for (auto &&f: reader)
    {

        cout << "entry_number " << f.__current_event.entry_number  << '\n';
//
        int weight = 1;
        if (mc)
        {
            
            std::vector<TruthParticle>&& truth_higgs = f.find_truth_particles({},{},{35});
            if (!(truth_higgs.empty()))
            {
                cout << "found!";
                std::vector<TruthParticle>&& truth_axions = f.find_truth_particles({},{truth_higgs[0].barcode()}, {36});
            }
            int temp = 1;
            std::vector<TruthParticle>&& truth_photons = f.find_truth_particles({},{},{22},&temp);
//
//            for (auto i: truth_photons)
//            {
//                cout << string(i) << '\n';
//            }

//
            fill_signal_hists(truth_photons,"truth");
//
//
            
            std::vector<TruthParticle>&& truth_leptons = f.find_truth_particles({},{},{11, 12, 13, 14, 15, 16, 17, 18},&temp);
            
            

            if (lepton_selection(truth_leptons))
            {
                fill_signal_hists(truth_leptons, "truth",1,"leptons");
            }
           
//
        }
//
//
        
//        vector<TruthParticle> photons;
        
        vector<Photon> photons;
//        std::cout << photons.size() << ' ' << f.__current_event.photons.size()
//        << '\n';
        std::copy_if (f.__current_event.photons.begin(), f.__current_event.photons.end(), std::back_inserter(photons), photon_selection );
        fill_signal_hists(photons, "reco", weight);
//
        if (photons.size() == 2)
        {
//            for (auto i: photons)
//            {
//                cout << string(i) << '\n';
//            }
            fill_signal_hists(photons,"reco_2y",weight);
            num_passed_events++;
        }
       
    }
    
    cout << "# events for " << systematic << " = " << num_passed_events
    << '\n';
    
}


void analyse_haa()
{
    auto start_time = Clock::now();
    cout << "Run over MC\n";

    string input_filename("../user.kschmied.28655874._000025.LGNTuple.root");
    
    
//    const char* output_filename = "example_mc_haa_out_test1cpp.root";
    const char* output_filename = "example_mc_haa_out_testcpp.root";
    TFile* output_file = TFile::Open(output_filename, "RECREATE");
    if (!output_file) {
        std::cout << "Error opening file\n";
       exit(-1);
    }
    

    std::vector<string> systematics = {"nominal"};
    
    for (string& systematic: systematics)
    {
        run_analysis(input_filename, systematic.c_str(), true);
    }
    
    for (auto& plot: plots)
    {
        plot.second.save(output_file);
    }

    for (auto& plot_group: plot_groups)
    {
        plot_group.second.save(output_file);
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
