#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <string>
#include <vector>
#include <array>
#include <iostream>

#include "TH1F.h"
#include "TFile.h"

//https://root-forum.cern.ch/t/problem-in-accessing-vector-vector-float/27983/2
//#include "LinkDef.h"
#include "plotting.h"
#include "objects.h"
#include "event.h"
#include "filereader.h"

using Clock = std::chrono::high_resolution_clock;

//debugging: https://drake.mit.edu/profiling.html
constexpr std::array<const char*,3> CUTS = {"truth","reco","reco_2y"};
std::unordered_map<std::string,Plot> plots;
std::unordered_map<std::string,PlotGroup> plot_groups;

double dR(TruthParticle& lepton1, TruthParticle& lepton2)
{
    return
//    pow(((lepton1.eta() - lepton2.eta())*(lepton1.eta() - lepton2.eta()) + (lepton1.phi() - lepton2.phi())*(lepton1.phi() - lepton2.phi())),0.5);
//    sqrt((lepton1.eta() - lepton2.eta())*(lepton1.eta() - lepton2.eta()) + (lepton1.phi() - lepton2.phi())*(lepton1.phi() - lepton2.phi()));
    lepton1.delta_r(lepton2);
}


bool lepton_selection(std::vector<TruthParticle>& leptons)
{
//    std::cout << dR(leptons[0],leptons[1]) << '\n';
    return
    (
     (leptons.size() == 2)
    &&
//    (dR(leptons[0],leptons[1]) > 0.01)
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

bool photon_selection(Photon& photon)
{
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

void run_analysis(std::vector<std::string>& input_filenames, std::string systematic = "nominal", bool mc = false)
{
//https://en.cppreference.com/w/cpp/utility/variant
    //or union
//https://www.sololearn.com/compiler-playground/cop9eIyns3c3
    
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
    
    FileReaderRange reader(input_filenames);
    Event::systematic = systematic;
    Event::load_tracks = true;
    int num_passed_events = 0;
    
    for (auto &&f: reader)
    {
        std::cout << "entry_number " << f.__current_event.entry_number  << '\n';
        
        int weight = 1;
        
//        for (auto i: f.__current_event.triggers)
//        {
//            std::cout << std::string(i) << '\n';
//        }
        
        if (mc)
        {
            std::vector<TruthParticle>&& truth_higgs = f.find_truth_particles({},{},{35});
            if (!(truth_higgs.empty()))
            {
                std::cout << "found!";
                std::vector<TruthParticle>&& truth_axions = f.find_truth_particles({},{truth_higgs[0].barcode()}, {36});
            }
            int temp = 1;
            std::vector<TruthParticle>&& truth_photons = f.find_truth_particles({},{},{22},&temp);
//
//            for (auto i: truth_photons)
//            {
//                std::cout << std::string(i) << '\n';
//            }
            

            fill_signal_hists(truth_photons,"truth");
//
            std::vector<TruthParticle>&& truth_leptons = f.find_truth_particles({},{},{11, 12, 13, 14, 15, 16, 17, 18},&temp);
            
            if (lepton_selection(truth_leptons))
            {
                fill_signal_hists(truth_leptons, "truth",1,"leptons");
            }
        }
//        std::vector<TruthParticle> photons;
        std::vector<Photon> photons;
//        std::cout << photons.size() << ' ' << f.__current_event.photons.size()
//        << '\n';
        std::copy_if (f.__current_event.photons.begin(), f.__current_event.photons.end(), std::back_inserter(photons), photon_selection );
        
        
        fill_signal_hists(photons, "reco", weight);
        if (photons.size() == 2)
        {
//            for (auto i: photons)
//            {
//                std::cout << std::string(i) << '\n';
//            }
            fill_signal_hists(photons,"reco_2y",weight);
            num_passed_events++;
        }
    }
    
    std::cout << "# events for " << systematic << " = " << num_passed_events
    << '\n';
}

void analyse_haa()
{
    auto start_time = Clock::now();
    std::cout << "Run over MC\n";

//    std::vector<std::string> input_filenames = {"/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000025.LGNTuple.root"};
    
    std::vector<std::string> input_filenames = {"../user.kschmied.28655874._000025.LGNTuple.root"};//,"../user.kschmied.28655874._000024.LGNTuple.root"};
    
//    std::vector<std::string> input_filenames = {"/home/common/Haa/ntuples/Za/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726.root"};
    
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
