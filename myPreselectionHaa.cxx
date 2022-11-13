#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <utility>
#include <unordered_set>
#include <fstream>
#include <iomanip>
#include <regex>

#include "TH1F.h"
#include "TFile.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

//https://root-forum.cern.ch/t/problem-in-accessing-vector-vector-float/27983/2
//#include "LinkDef.h"
#include "plotting.h"
#include "objects.h"
#include "event.h"
#include "filereader.h"
#include "candidate.h"

using Clock = std::chrono::high_resolution_clock;

//debugging: https://drake.mit.edu/profiling.html

std::unordered_map<std::string,Plot> plots;

bool photon_selection(Photon& photon)
{
    if (!photon.id_loose())
    {
        return false;
    }
    
    if (photon.pt() <= 10000)
    {
        return false;
    }
    
    if (abs(photon.eta()) >= 2.37)
    {
        return false;
    }
    
    if ((1.37 < abs(photon.eta())) && (abs(photon.eta()) < 1.52))
    {
        return false;
    }
    return true;
}

bool lepton_selection(CandidateSet<TruthParticle>& leptons)
{
//    std::cout << dR(leptons[0],leptons[1]) << '\n';
    return
     (
    (leptons.particle_a.delta_r(leptons.particle_b) > 0.01)
    &&
    (leptons.particle_a.charge() == -1*leptons.particle_b.charge())
    &&
    ((leptons.four_momentum.M()/1e3 >= 81) &&
        (leptons.four_momentum.M()/1e3 <= 101))
    &&
    (leptons.four_momentum.Pt()/1e3 > 10)
    &&
    ((leptons.particle_a.pt() > 20e3 && leptons.particle_b.pt() > 27e3)
        ||
    (leptons.particle_b.pt() > 20e3 && leptons.particle_a.pt() > 27e3))
      );
}

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

void run_analysis(const std::vector<std::string>& input_filenames, std::string systematic = "nominal", bool mc = false, TFile* output_file = nullptr, std::string prefix = "")
{
    plots.clear();
//    ROOT::RVec<double> points;
    std::cout <<"systematic = " << systematic << '\n';
    const int MaxBins = 100, minBins = 100, inc = 100;
    
    for (int i=minBins; i<=MaxBins; i+=inc)
    {
        plots.emplace(prefix+"Preselectrion di-lepton mass distribution, bins = " + std::to_string(i),
            Plot(prefix+"Preselectrion di-lepton mass distribution, bins = " + std::to_string(i),
                 prefix+"Preselectrion di-lepton mass distribution, bins = " + std::to_string(i), i, 60, 120));
        
        plots.emplace(prefix+"Di-lepton p_T distribution before pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"Di-lepton p_T distribution before pre-selection, bins = " + std::to_string(i),
                 prefix+"Di-lepton p_T distribution before pre-selection, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i),
                 prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i), i, 0, 80));
        
        plots.emplace(prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i),
                 prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i), i, 0, 80));
        
        plots.emplace(prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i),
                 prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i), i, 0, 80));
        
        plots.emplace(prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i),
                 prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i), i, 0, 50));
        
        plots.emplace(prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i),
                 prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i), i, 0, 50));
        
        plots.emplace(prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i),
                 prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i), i, 0, 50));
        
        plots.emplace(prefix+"transverse momentum of photon pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"transverse momentum of photon pre-selection, bins = " + std::to_string(i),
                 prefix+"transverse momentum of photon pre-selection, bins = " + std::to_string(i), i, 0, 100));
        
        plots.emplace(prefix+"transverse momentum of Z-boson pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"transverse momentum of Z-boson pre-selection, bins = " + std::to_string(i),
                 prefix+"transverse momentum of Z-boson pre-selection, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"transverse momentum of axion pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"transverse momentum of axion pre-selection, bins = " + std::to_string(i),
                 prefix+"transverse momentum of axion pre-selection, bins = " + std::to_string(i), i, 0, 60));
        
        plots.emplace(prefix+"transverse momentum of higgs pre-selection, bins = " + std::to_string(i),
            Plot(prefix+"transverse momentum of higgs pre-selection, bins = " + std::to_string(i),
                 prefix+"transverse momentum of higgs pre-selection, bins = " + std::to_string(i), i, 0, 217));
        
        plots.emplace(prefix+"delta R between the reconstructed Z and a systems, bins = " + std::to_string(i),
            Plot(prefix+"delta R between the reconstructed Z and a systems, bins = " + std::to_string(i),
                 prefix+"delta R between the reconstructed Z and a systems, bins = " + std::to_string(i), i, 0, 6.5));
        
        plots.emplace(prefix+"delta phi between the reconstructed Z and a systems, bins = " + std::to_string(i),
            Plot(prefix+"delta phi between the reconstructed Z and a systems, bins = " + std::to_string(i),
                 prefix+"delta phi between the reconstructed Z and a systems, bins = " + std::to_string(i), i, 0, 6.25));
        
        plots.emplace(prefix+"delta eta between the reconstructed Z and a systems, bins = " + std::to_string(i),
            Plot(prefix+"delta eta between the reconstructed Z and a systems, bins = " + std::to_string(i),
                 prefix+"delta eta between the reconstructed Z and a systems, bins = " + std::to_string(i), i, 0, 6.25));
        
    }
    
    FileReaderRange reader(input_filenames);
    Event::systematic = systematic;
    
    int weight = 1;
    
//    std::unordered_set<std::string> all_triggers;
    
    int beforePreselection = 0,
        two_leptons = 0,
        opp_charge = 0,
        lep1_lep2_pt = 0,
        lep_same_flavor = 0,
        dilep_mass = 0,
        dilep_pt = 0;
    
    constexpr std::array<const char*,35> triggers = {
        "HLT_e26_lhtight_nod0_ivarloose",
        "HLT_e60_lhmedium_nod0",
        "HLT_e140_lhloose_nod0",
        "HLT_mu26_ivarmedium",
        "HLT_mu50",
        
        "HLT_e26_lhtight_nod0_ivarloose",
        "HLT_e60_lhmedium_nod0",
        "HLT_e140_lhloose_nod0",
        "HLT_mu26_ivarmedium",
        "HLT_mu50",
        
        "HLT_e24_lhtight_nod0_ivarloose",
        "HLT_e26_lhtight_nod0_ivarloose",
        "HLT_e60_lhmedium_nod0",
        "HLT_e60_medium",
        "HLT_e140_lhloose_nod0",
        "HLT_mu26_ivarmedium",
        "HLT_mu50",
        
        "HLT_e24_lhmedium_L1EM20VH",
        "HLT_e60_lhmedium",
        "HLT_e120_lhloose",
        "HLT_mu20_iloose_L1MU15",
        "HLT_mu50",
        
        "HLT_2e17_lhvloose_nod0_L12EM15VHI",
        "HLT_2e17_lhvloose_nod0",
        "HLT_2e24_lhvloose_nod0",
        "HLT_mu22_mu8noL1",
        
        "HLT_2e17_lhvloose_nod0_L12EM15VHI",
        "HLT_2e17_lhvloose_nod0",
        "HLT_2e24_lhvloose_nod0",
        "HLT_mu22_mu8noL1",
        
        "HLT_2e15_lhvloose_nod0_L12EM13VHI",
        "HLT_2e17_lhvloose_nod0",
        "HLT_mu22_mu8noL1",
        
        "HLT_2e12_lhvloose_L12EM10VH",
        "HLT_mu18_mu8noL1",
    };
    
    for (auto &&f: reader)
    {
        std::cout << "entry_number " << f.__current_event.entry_number  << '\n';
        //        std::copy(f.__current_event.triggers.begin(),f.__current_event.triggers.end(),std::inserter(all_triggers,all_triggers.end()));
        //        std::cout<<'\n';
        
        
        auto passBeforePreselection = [&](){
            const auto trigger_found = std::find_first_of(f.__current_event.triggers.begin(), f.__current_event.triggers.end(), triggers.begin(), triggers.end());
            if (trigger_found == f.__current_event.triggers.end())
            {
                return false;
            }
//            std::vector<TruthParticle>&& truth_muon_electrons = f.find_truth_particles({},{},{11,-11,13,-13},&weight);
//            for (auto& mu_elec: truth_muon_electrons)
//            {
//                if (abs(mu_elec.pdg_id)==11)
//                {
//                    if (abs(mu_elec.eta()) >= 2.37)
//                    {
//                        return false;
//                    }
//                    if ((1.37 < abs(mu_elec.eta())) && (abs(mu_elec.eta()) < 1.52))
//                    {
//                        return false;
//                    }
//                    if (mu_elec.pt() <= 20e3)
//                    {
//                        return false;
//                    }
//                }
//                else if (abs(mu_elec.pdg_id)==13)
//                {
//                    if (abs(mu_elec.eta()) >= 2.5)
//                    {
//                        return false;
//                    }
//                    if (mu_elec.pt() <= 15e3)
//                    {
//                        return false;
//                    }
//                }
//            }
            return true;
        };
        
        if (passBeforePreselection())
        {
            beforePreselection++;
            
            std::vector<TruthParticle>&& truth_leptons = f.find_truth_particles({},{},{11, 12, 13, 14, 15, 16, 17, 18},&weight);
            
            if (truth_leptons.size() == 2)
            {
                two_leptons++;
                CandidateSet<TruthParticle> candidate_dilepton(std::make_pair(truth_leptons[0],truth_leptons[1]));
                if (truth_leptons[0].charge() == -1*truth_leptons[1].charge())
                {
                    opp_charge++;
                    if ((truth_leptons[0].pt() > 20e3 && truth_leptons[1].pt() > 27e3)
                        ||
                        (truth_leptons[1].pt() > 20e3 && truth_leptons[0].pt() > 27e3))
                    {
                        lep1_lep2_pt++;
                        if (truth_leptons[0].same_flavour(truth_leptons[1]))
                        {
                            //https://particle.wiki/wiki/PDG_particle_numbering_scheme
                            lep_same_flavor++;
                            if ((candidate_dilepton.four_momentum.M()/1e3 >= 81) &&
                                (candidate_dilepton.four_momentum.M()/1e3 <= 101))
                            {
                                dilep_mass++;
                                if (candidate_dilepton.four_momentum.Pt()/1e3 > 10)
                                {
                                    dilep_pt++;
                                }
                            }
                        }
                    }
                    
                }
                
                for (int i=minBins; i<=MaxBins; i+=inc)
                {
                    plots.at(prefix+"Di-lepton p_T distribution before pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.four_momentum.Pt()/1e3,weight);
                }
                
                if (lepton_selection(candidate_dilepton))
                {
                    for (int i=minBins; i<=MaxBins; i+=inc)
                    {
                        plots.at(prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_a.pt()/1e3,weight);
                        plots.at(prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_b.pt()/1e3,weight);
                        if (candidate_dilepton.particle_b.pt() > candidate_dilepton.particle_a.pt())
                        {
                            plots.at(prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_b.pt()/1e3,weight);
                            plots.at(prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_a.pt()/1e3,weight);
                        }
                        else
                        {
                            plots.at(prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_a.pt()/1e3,weight);
                            plots.at(prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_b.pt()/1e3,weight);
                        }
                        
                        plots.at(prefix+"Preselectrion di-lepton mass distribution, bins = " + std::to_string(i)).fill(candidate_dilepton.four_momentum.M()/1e3,weight);
                        
                        plots.at(prefix+"transverse momentum of Z-boson pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.four_momentum.Pt()/1e3,weight);
                        
                        filter<Photon>(f.__current_event.photons,&photon_selection);
                        
                        if (f.__current_event.photons.size()==2)
                        {
                            plots.at(prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[0].pt()/1e3);
                            plots.at(prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[1].pt()/1e3);
                            if (f.__current_event.photons[0].pt() > f.__current_event.photons[1].pt())
                            {
                                plots.at(prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[0].pt()/1e3);
                                plots.at(prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[1].pt()/1e3);
                            }
                            else
                            {
                                plots.at(prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[1].pt()/1e3);
                                plots.at(prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[0].pt()/1e3);
                            }
                            
                            CandidateSet<Photon> candidate_diphoton(std::make_pair(f.__current_event.photons[0],f.__current_event.photons[1]));
                            if ((candidate_diphoton.four_momentum.Pt()/1e3) > 20)
                            {
                                plots.at(prefix+"transverse momentum of axion pre-selection, bins = " + std::to_string(i)).fill(candidate_diphoton.four_momentum.Pt()/1e3,weight);
                                
                                PtEtaPhiEVector tVecBoson = candidate_diphoton.four_momentum + candidate_dilepton.four_momentum;
                                
                                plots.at(prefix+"transverse momentum of higgs pre-selection, bins = " + std::to_string(i)).fill(tVecBoson.Pt()/1e3,weight);
                                
                                plots.at(prefix+"delta R between the reconstructed Z and a systems, bins = " + std::to_string(i)).fill(VectorUtil::DeltaR(candidate_diphoton.four_momentum,candidate_dilepton.four_momentum));
                                
                                plots.at(prefix+"delta phi between the reconstructed Z and a systems, bins = " + std::to_string(i)).fill(VectorUtil::DeltaPhi(candidate_diphoton.four_momentum,candidate_dilepton.four_momentum));
                                
                                plots.at(prefix+"delta eta between the reconstructed Z and a systems, bins = " + std::to_string(i)).fill(abs(candidate_diphoton.four_momentum.Eta() - candidate_dilepton.four_momentum.Eta()));
                                
                            }
                            
                        }
                        
                        else if (f.__current_event.photons.size()==1)
                        {
                            plots.at(prefix+"transverse momentum of photon pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[0].pt()/1e3);
                        }
                            
                        
                    }
                }
            }
        }
    }
            
    
    
//    std::ofstream out("triggers.txt", std::ios::app);
//
//    for (auto &&t: all_triggers)
//    {
//        out << t << '\n';
//    }
//    out.close();
//    std::cout << "Mean = " << Mean(points) << "\nStdDev = " << StdDev(points) << '\n';
    
    std::cout << "\n\n\n";
    std::cout << std::right << std::string(25,' ')
    << std::right << std::setw(25) << "beforePreselection"
    << std::right << std::setw(25) << "two_leptons"
    << std::right << std::setw(25) << "opp_charge"
    << std::right << std::setw(25) << "lep1_lep2_pt"
    << std::right << std::setw(25) << "lep_same_flavor"
    << std::right << std::setw(25) << "dilep_mass"
    << std::right << std::setw(25) << "dilep_pt\n";
    
    std::cout << std::right << std::setw(25) << "My Event numbers"
    << std::right << std::setw(25) << beforePreselection
    << std::right << std::setw(25) << two_leptons
    << std::right << std::setw(25) << opp_charge
    << std::right << std::setw(25) << lep1_lep2_pt
    << std::right << std::setw(25) << lep_same_flavor
    << std::right << std::setw(25) << dilep_mass
    << std::right << std::setw(25) << dilep_pt << '\n';
    
    std::cout << std::right << std::setw(25) << "My Event ratios"
    << std::right << std::setw(25) << beforePreselection/beforePreselection
    << std::right << std::setw(25) << static_cast<double>(two_leptons)/beforePreselection
    << std::right << std::setw(25) << static_cast<double>(opp_charge)/beforePreselection
    << std::right << std::setw(25) << static_cast<double>(lep1_lep2_pt)/beforePreselection
    << std::right << std::setw(25) << static_cast<double>(lep_same_flavor)/beforePreselection
    << std::right << std::setw(25) << static_cast<double>(dilep_mass)/beforePreselection
    << std::right << std::setw(25) << static_cast<double>(dilep_pt)/beforePreselection << '\n';
    
    //values from paper here
    std::cout << std::right << std::setw(25) << "Paper's Event ratios"
    << std::right << std::setw(25) << 23580.1/23580.1
    << std::right << std::setw(25) << 21606.75/23580.1
    << std::right << std::setw(25) << 21505.03/23580.1
    << std::right << std::setw(25) << 21375.48/23580.1
    << std::right << std::setw(25) << 21375.33/23580.1
    << std::right << std::setw(25) << 20543.06/23580.1
    << std::right << std::setw(25) << 19516.87/23580.1 << '\n';
    std::cout << "\n\n\n";
    
//    auto df1 = ROOT::RDataFrame(1).Define("beforePreselection", [&] { return beforePreselection; }).Define("two_leptons", [&] { return two_leptons; }).Define("opp_charge", [&] { return opp_charge; }).Define("lep1_lep2_pt", [&] { return lep1_lep2_pt; }).Define("lep_same_flavor", [&] { return lep_same_flavor; }).Define("dilep_mass", [&] { return dilep_mass; }).Define("dilep_pt", [&] { return dilep_pt; });
    
    
//    auto df1 = ROOT::RDataFrame(1).Define("beforePreselection", [&] { return 1; }).Define("two_leptons", [&] { return 2; }).Define("opp_charge", [&] { return 3; }).Define("lep1_lep2_pt", [&] { return 4; }).Define("lep_same_flavor", [&] { return 5; }).Define("dilep_mass", [&] { return 6; }).Define("dilep_pt", [&] { return 7; });
    
    
//    std::cout << df1.Display({"beforePreselection","two_leptons","opp_charge","lep1_lep2_pt","lep_same_flavor","dilep_mass","dilep_pt"},1)->AsString() << '\n';
    
//    std::cout << df1.Display({"beforePreselection","two_leptons","opp_charge","lep1_lep2_pt","lep_same_flavor","dilep_mass","dilep_pt"},1)->AsString() << '\n';
//
  
//    df1.Redefine("beforePreselection", [] {return 2;});
    
    for (auto& plot: plots)
    {
        plot.second.save(output_file);
    }

}

void myPreselectionHaa()
{
    auto start_time = Clock::now();
    std::cout << "Run over MC\n";
    
    std::vector<std::vector<std::string>> input_filenames = {{"/home/common/Haa/ntuples/Za/mc16_13TeV.600909.PhPy8EG_AZNLO_ggH125_mA5p0_Cyy0p01_Czh1p0.merge.AOD.e8324_e7400_s3126_r10724_r10726_v2.root"},
    {"/home/common/Haa/ntuples/Za/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v2.root"}};
    
    const char* output_filename = "mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA_p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v2_out.root";
    
    TFile* output_file = TFile::Open(output_filename, "RECREATE");
    if (!output_file) {
        std::cout << "Error opening file\n";
       exit(-1);
    }
    
    for (const auto& input_filename: input_filenames)
    {
        std::smatch m;
        regex_search(input_filename[0], m, std::regex("mA[0-9]+"));
        
        run_analysis(input_filename, "nominal", true, output_file, input_filename[0].substr(m.position(),m.length())+" ");
    }
        

    output_file->Close();
    auto end_time = Clock::now();
    std::cout << "Time difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
}

int main()
{
    myPreselectionHaa();
}
