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
#include <regex>
#include <cstdlib>

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

template <typename T>
bool lepton_selection(CandidateSet<T>& leptons)
{
    return
     (
    (leptons.particle_a.delta_r(leptons.particle_b) > 0.01)
    &&
//    ((leptons.particle_a.charge()*leptons.particle_b.charge()) < 0)
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
//        plots.emplace(prefix+"Preselectrion di-lepton mass distribution, bins = " + std::to_string(i),
//            Plot(prefix+"Preselectrion di-lepton mass distribution, bins = " + std::to_string(i),
//                 prefix+"Preselectrion di-lepton mass distribution, bins = " + std::to_string(i), i, 60, 120));
//
//        plots.emplace(prefix+"Di-lepton p_T distribution before pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"Di-lepton p_T distribution before pre-selection, bins = " + std::to_string(i),
//                 prefix+"Di-lepton p_T distribution before pre-selection, bins = " + std::to_string(i), i, 0, 200));
//
//        plots.emplace(prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i),
//                 prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i), i, 0, 80));
//
//        plots.emplace(prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i),
//                 prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i), i, 0, 80));
//
//        plots.emplace(prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i),
//                 prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i), i, 0, 80));
//
//        plots.emplace(prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i),
//                 prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i), i, 0, 50));
//
//        plots.emplace(prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i),
//                 prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i), i, 0, 50));
//
//        plots.emplace(prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i),
//                 prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i), i, 0, 50));
//
//        plots.emplace(prefix+"transverse momentum of photon pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"transverse momentum of photon pre-selection, bins = " + std::to_string(i),
//                 prefix+"transverse momentum of photon pre-selection, bins = " + std::to_string(i), i, 0, 100));
//
//        plots.emplace(prefix+"transverse momentum of Z-boson pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"transverse momentum of Z-boson pre-selection, bins = " + std::to_string(i),
//                 prefix+"transverse momentum of Z-boson pre-selection, bins = " + std::to_string(i), i, 0, 200));
//
//        plots.emplace(prefix+"transverse momentum of axion pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"transverse momentum of axion pre-selection, bins = " + std::to_string(i),
//                 prefix+"transverse momentum of axion pre-selection, bins = " + std::to_string(i), i, 0, 60));
//
//        plots.emplace(prefix+"transverse momentum of higgs pre-selection, bins = " + std::to_string(i),
//            Plot(prefix+"transverse momentum of higgs pre-selection, bins = " + std::to_string(i),
//                 prefix+"transverse momentum of higgs pre-selection, bins = " + std::to_string(i), i, 0, 217));
//
//        plots.emplace(prefix+"delta R between the reconstructed Z and a systems, bins = " + std::to_string(i),
//            Plot(prefix+"delta R between the reconstructed Z and a systems, bins = " + std::to_string(i),
//                 prefix+"delta R between the reconstructed Z and a systems, bins = " + std::to_string(i), i, 0, 6.5));
//
//        plots.emplace(prefix+"delta phi between the reconstructed Z and a systems, bins = " + std::to_string(i),
//            Plot(prefix+"delta phi between the reconstructed Z and a systems, bins = " + std::to_string(i),
//                 prefix+"delta phi between the reconstructed Z and a systems, bins = " + std::to_string(i), i, 0, 6.25));
//
//        plots.emplace(prefix+"delta eta between the reconstructed Z and a systems, bins = " + std::to_string(i),
//            Plot(prefix+"delta eta between the reconstructed Z and a systems, bins = " + std::to_string(i),
//                 prefix+"delta eta between the reconstructed Z and a systems, bins = " + std::to_string(i), i, 0, 6.25));
        
        plots.emplace(prefix+"pt distribution of all electron objects, bins = " + std::to_string(i),
            Plot(prefix+"pt distribution of all electron objects, bins = " + std::to_string(i),
                 prefix+"pt distribution of all electron objects, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"eta distribution of all electron objects, bins = " + std::to_string(i),
            Plot(prefix+"eta distribution of all electron objects, bins = " + std::to_string(i),
                 prefix+"eta distribution of all electron objects, bins = " + std::to_string(i), i, -0.0075, 0.0075));
        
        plots.emplace(prefix+"dilep pt before dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep pt before dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep pt before dilep cuts electrons objects, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"dilep mass before dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep mass before dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep mass before dilep cuts electrons objects, bins = " + std::to_string(i), i, 60, 120));
        
        plots.emplace(prefix+"dilep delta R before dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep delta R before dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep delta R before dilep cuts electrons objects, bins = " + std::to_string(i), i, 0, 6.5));
        
        plots.emplace(prefix+"dilep delta eta before dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep delta eta before dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep delta eta before dilep cuts electrons objects, bins = " + std::to_string(i), i, 0, 6.25));
        
        plots.emplace(prefix+"dilep pt after all dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep pt after all dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep pt after all dilep cuts electrons objects, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"dilep mass after all dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep mass after all dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep mass after all dilep cuts electrons objects, bins = " + std::to_string(i), i, 60, 120));
        
        plots.emplace(prefix+"dilep delta R after all dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep delta R after all dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep delta R after all dilep cuts electrons objects, bins = " + std::to_string(i), i, 0, 6.5));
        
        plots.emplace(prefix+"dilep delta eta after all dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep delta eta after all dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep delta eta after all dilep cuts electrons objects, bins = " + std::to_string(i), i, 0, 6.25));
    }
    
    FileReaderRange reader(input_filenames);
    
    Event::systematic = systematic;
    
    int weight = 1;
    
//    std::unordered_set<std::string> all_triggers;
    std::unordered_map <int,int> all_Z_products;
    
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
    
    
    int allEvents = 0, nZee = 0, nZeeFid = 0, nZmumu = 0, nZmumuFid = 0,
    elecObj = 0, exactlyTwo = 0, oppCharge = 0, leadingPt = 0, dR = 0,
    mLL = 0, pTLL = 0;
    
    auto findParentInChain = [](int targetBarcode, std::vector<TruthParticle>& startParticles, std::vector<TruthParticle>& truthChain)
    {
        std::vector<TruthParticle> truthSelected;
        bool foundParent;
        if (truthChain.size() >= 2)
        {
            TruthParticle tp;
            for (auto& tpe: startParticles)
            {
                tp = tpe;
                while (true)
                {
                    if (tp.parent_barcode() == targetBarcode)
                    {
                        truthSelected.push_back(tp);
                        break;
                    }
                    else
                    {
                        foundParent = false;
                        for (auto& tmp: truthChain)
                        {
                            if (tp.parent_barcode() == tmp.barcode())
                            {
                                tp = tmp;
                                foundParent = true;
                                break;
                            }
                        }
                        if (foundParent == false)
                        {
                            break;
                        }
                    }
                }
            }
        }
        return truthSelected;
    };
    
    auto sigElecCut = [](TruthParticle& tp)
    {
        return ((tp.pt()/1e3 > 20) && (abs(tp.eta()) < 2.37) &&
        !((1.37 < abs(tp.eta())) && (abs(tp.eta()) < 1.52)));
    };
    
    auto sigMuonCut = [](TruthParticle& tp)
    {
        return ((tp.pt()/1e3 > 20) && (abs(tp.eta()) < 2.4));
    };
    
    auto elecObjCuts = [](Electron& ep)
    {
        return ((ep.pt()/1e3 > 20) && (abs(ep.eta()) < 2.37) &&
        !((1.37 < abs(ep.eta())) && (abs(ep.eta()) < 1.52)) &&
                (ep.id_medium() == 1));
    };
    
    
    for (auto &&f: reader)
    {
        std::cout << "entry_number " << f.__current_event.entry_number  << '\n';
        //        std::copy(f.__current_event.triggers.begin(),f.__current_event.triggers.end(),std::inserter(all_triggers,all_triggers.end()));
        //        std::cout<<'\n';
        
        allEvents++;
        auto passBeforePreselection = [&](){
//            const auto trigger_found = std::find_first_of(f.__current_event.triggers.begin(), f.__current_event.triggers.end(), triggers.begin(), triggers.end());
//            if (trigger_found == f.__current_event.triggers.end())
//            {
//                return false;
//            }
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
        
        //cutflow
        if (passBeforePreselection())
        {
            std::vector<TruthParticle> truth_leptons;
            std::vector<TruthParticle> truth_electrons,
            truth_stableElectrons, truth_mu, truth_stableMu, truthSignalElectrons, truthSignalMu;
            std::vector<TruthParticle>&& truth_Zbosons = f.find_truth_particles({},{},{23},&weight, true);

            if (!(truth_Zbosons.empty()))
            {
                truth_leptons = f.find_truth_particles({},{truth_Zbosons[0].barcode()},{11, -11});
                
                truth_electrons = f.find_truth_particles({},{},{11, -11});
                truth_stableElectrons = f.find_truth_particles({},{},{11, -11},&weight);
                truth_mu = f.find_truth_particles({},{},{13, -13});
                truth_stableMu = f.find_truth_particles({},{},{13, -13},&weight);
            
                if (truth_electrons.size() >= 2)
                {
                    truthSignalElectrons = findParentInChain(truth_Zbosons[0].barcode(), truth_stableElectrons, truth_electrons);
                    if (truthSignalElectrons.size() == 2)
                    {
                        nZee++;
                    }
                    filter<TruthParticle>(truthSignalElectrons,sigElecCut);
                    if (truthSignalElectrons.size() == 2)
                    {
                        nZeeFid++;
                    }
                }
                if (truth_mu.size() >= 2)
                {
                    truthSignalMu = findParentInChain(truth_Zbosons[0].barcode(), truth_stableMu, truth_mu);
                    if (truthSignalMu.size() == 2)
                    {
                        nZmumu++;
                    }
                    filter<TruthParticle>(truthSignalMu,sigMuonCut);
                    if (truthSignalMu.size() == 2)
                    {
                        nZmumuFid++;
                    }
                }
            }
            
            std::vector<Electron>& reco_electrons = f.__current_event.electrons;

            for (int i=minBins; i<=MaxBins; i+=inc)
            {
                for (auto& j: reco_electrons)
                {
                    plots.at(prefix+"pt distribution of all electron objects, bins = " + std::to_string(i)).fill(j.pt()/1e3,weight);
                    plots.at(prefix+"eta distribution of all electron objects, bins = " + std::to_string(i)).fill(j.eta()/1e3,weight);
                }
            }
            
            filter<Electron>(reco_electrons, elecObjCuts);
            
            if (reco_electrons.size() >= 2)
            {
                elecObj++;
                if (reco_electrons.size() == 2)
                {
                    exactlyTwo++;
                    CandidateSet<Electron> candidateDiElectron(std::make_pair(reco_electrons[0],reco_electrons[1]));
                    for (int i=minBins; i<=MaxBins; i+=inc)
                    {
                        plots.at(prefix+"dilep pt before dilep cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.four_momentum.Pt()/1e3,weight);
                        plots.at(prefix+"dilep mass before dilep cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.four_momentum.M()/1e3,weight);
                        plots.at(prefix+"dilep delta R before dilep cuts electrons objects, bins = " + std::to_string(i)).fill(reco_electrons[0].delta_r(reco_electrons[1]),weight);
                        plots.at(prefix+"dilep delta eta before dilep cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.delta_eta(),weight);
                    }
                    
                    
                    if (reco_electrons[0].charge()*reco_electrons[1].charge() < 0)
                    {
                        oppCharge++;
                        if ((reco_electrons[0].pt() > 20e3 && reco_electrons[1].pt() > 27e3)
                            ||
                            (reco_electrons[1].pt() > 20e3 && reco_electrons[0].pt() > 27e3))
                        {
                            leadingPt++;
                            
                            if (reco_electrons[0].delta_r(reco_electrons[1]) > 0.01)
                            {
                                dR++;
                                
                                if ((candidateDiElectron.four_momentum.M()/1e3 >= 81) &&
                                    (candidateDiElectron.four_momentum.M()/1e3 <= 101))
                                {
                                    mLL++;
                                    if (candidateDiElectron.four_momentum.Pt()/1e3 > 10)
                                    {
                                        pTLL++;
                                        for (int i=minBins; i<=MaxBins; i+=inc)
                                        {
                                            plots.at(prefix+"dilep pt after all dilep cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.four_momentum.Pt()/1e3,weight);
                                            plots.at(prefix+"dilep mass after all dilep cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.four_momentum.M()/1e3,weight);
                                            plots.at(prefix+"dilep delta R after all dilep cuts electrons objects, bins = " + std::to_string(i)).fill(reco_electrons[0].delta_r(reco_electrons[1]),weight);
                                            plots.at(prefix+"dilep delta eta after all dilep cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.delta_eta(),weight);
                                        }
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }
            
            
//            if (truth_leptons.size() == 2)
//            {
//                two_leptons++;
//                CandidateSet<TruthParticle> candidate_dilepton(std::make_pair(truth_leptons[0],truth_leptons[1]));
////                CandidateSet<Electron> candidate_dilepton(std::make_pair(truth_leptons[0],truth_leptons[1]));
//                if (truth_leptons[0].charge() == -1*truth_leptons[1].charge())
//                {
//                    opp_charge++;
//                    if ((truth_leptons[0].pt() > 20e3 && truth_leptons[1].pt() > 27e3)
//                        ||
//                        (truth_leptons[1].pt() > 20e3 && truth_leptons[0].pt() > 27e3))
//                    {
//                        lep1_lep2_pt++;
//                        if (truth_leptons[0].same_flavour(truth_leptons[1]))
//                        {
//                            //https://particle.wiki/wiki/PDG_particle_numbering_scheme
//                            lep_same_flavor++;
//                            if ((candidate_dilepton.four_momentum.M()/1e3 >= 81) &&
//                                (candidate_dilepton.four_momentum.M()/1e3 <= 101))
//                            {
//                                dilep_mass++;
//                                if (candidate_dilepton.four_momentum.Pt()/1e3 > 10)
//                                {
//                                    dilep_pt++;
//                                }
//                            }
//                        }
//                    }
//
//                }
//
//                for (int i=minBins; i<=MaxBins; i+=inc)
//                {
//                    plots.at(prefix+"Di-lepton p_T distribution before pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.four_momentum.Pt()/1e3,weight);
//                }
//
////                if (lepton_selection<Electron>(candidate_dilepton))
//                if (lepton_selection<TruthParticle>(candidate_dilepton))
//                {
//                    for (int i=minBins; i<=MaxBins; i+=inc)
//                    {
//                        plots.at(prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_a.pt()/1e3,weight);
//                        plots.at(prefix+"Transverse momentum of leading and subleading leptons after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_b.pt()/1e3,weight);
//                        if (candidate_dilepton.particle_b.pt() > candidate_dilepton.particle_a.pt())
//                        {
//                            plots.at(prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_b.pt()/1e3,weight);
//                            plots.at(prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_a.pt()/1e3,weight);
//                        }
//                        else
//                        {
//                            plots.at(prefix+"Transverse momentum of leading lepton after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_a.pt()/1e3,weight);
//                            plots.at(prefix+"Transverse momentum of subleading lepton after pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.particle_b.pt()/1e3,weight);
//                        }
//
//                        plots.at(prefix+"Preselectrion di-lepton mass distribution, bins = " + std::to_string(i)).fill(candidate_dilepton.four_momentum.M()/1e3,weight);
//
//                        plots.at(prefix+"transverse momentum of Z-boson pre-selection, bins = " + std::to_string(i)).fill(candidate_dilepton.four_momentum.Pt()/1e3,weight);
//
//                        filter<Photon>(f.__current_event.photons,&photon_selection);
//
//                        if (f.__current_event.photons.size()==2)
//                        {
//                            plots.at(prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[0].pt()/1e3);
//                            plots.at(prefix+"Transverse momentum of leading and subleading photons after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[1].pt()/1e3);
//                            if (f.__current_event.photons[0].pt() > f.__current_event.photons[1].pt())
//                            {
//                                plots.at(prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[0].pt()/1e3);
//                                plots.at(prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[1].pt()/1e3);
//                            }
//                            else
//                            {
//                                plots.at(prefix+"Transverse momentum of leading photon after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[1].pt()/1e3);
//                                plots.at(prefix+"Transverse momentum of subleading photon after pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[0].pt()/1e3);
//                            }
//
//                            CandidateSet<Photon> candidate_diphoton(std::make_pair(f.__current_event.photons[0],f.__current_event.photons[1]));
//                            if ((candidate_diphoton.four_momentum.Pt()/1e3) > 20)
//                            {
//                                plots.at(prefix+"transverse momentum of axion pre-selection, bins = " + std::to_string(i)).fill(candidate_diphoton.four_momentum.Pt()/1e3,weight);
//
//                                PtEtaPhiEVector tVecBoson = candidate_diphoton.four_momentum + candidate_dilepton.four_momentum;
//
//                                plots.at(prefix+"transverse momentum of higgs pre-selection, bins = " + std::to_string(i)).fill(tVecBoson.Pt()/1e3,weight);
//
//                                plots.at(prefix+"delta R between the reconstructed Z and a systems, bins = " + std::to_string(i)).fill(VectorUtil::DeltaR(candidate_diphoton.four_momentum,candidate_dilepton.four_momentum));
//
//                                plots.at(prefix+"delta phi between the reconstructed Z and a systems, bins = " + std::to_string(i)).fill(VectorUtil::DeltaPhi(candidate_diphoton.four_momentum,candidate_dilepton.four_momentum));
//
//                                plots.at(prefix+"delta eta between the reconstructed Z and a systems, bins = " + std::to_string(i)).fill(abs(candidate_diphoton.four_momentum.Eta() - candidate_dilepton.four_momentum.Eta()));
//                            }
//                        }
//
//                        else if (f.__current_event.photons.size()==1)
//                        {
//                            plots.at(prefix+"transverse momentum of photon pre-selection, bins = " + std::to_string(i)).fill(f.__current_event.photons[0].pt()/1e3);
//                        }
//                    }
//                }
//            }
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
    
    std::ofstream out(prefix.substr(0,3)+".txt");
    
//    out <<
//    "TwoLeptons,OppCharge,Lep1Lep2Pt,LepSameFlavor,DilepMass,DilepPt\n"
    
    out << "allEvents,$\\geq 2 e^-$,$2e^-$,Opp Charge, $p_T^{\\mathrm{Leading}} > 27$ GeV,$\\Delta R > 0.01$,81 GeV $\\leq m_{ll} \\leq$ 101 GeV,$p_T^{ll}>$ 10 GeV\n"
    
    << allEvents << ',' << elecObj << ',' << exactlyTwo << ',' <<
    oppCharge << ',' << leadingPt << ',' << dR << ',' << mLL << ',' << pTLL << '\n';
    
//    << two_leptons << ',' << opp_charge << ','
//    << lep1_lep2_pt << ',' << lep_same_flavor << ',' << dilep_mass << ','
//    << dilep_pt << '\n'
    
//    << static_cast<double>(two_leptons)/static_cast<double>(two_leptons)
//    << ',' << static_cast<double>(opp_charge)/static_cast<double>(two_leptons) << ','
//    << static_cast<double>(lep1_lep2_pt)/static_cast<double>(two_leptons) << ','
//    << static_cast<double>(lep_same_flavor)/static_cast<double>(two_leptons) << ','
//    << static_cast<double>(dilep_mass)/static_cast<double>(two_leptons) << ','
//    << static_cast<double>(dilep_pt)/static_cast<double>(two_leptons) << '\n';
//
//    switch (prefix[2]) {
//        case '1':
//            out << 21606.75/21606.75 << ',' << 21505.03/21606.75
//            << ',' << 21375.48/21606.75 << ',' << 21375.33/21606.75 << ','
//            << 20543.06/21606.75 << ',' << 19516.87/21606.75 << '\n';
//            break;
//        case '5':
//            out << 21655.38/21655.38 << ',' << 21556.03/21655.38
//            << ',' << 21424.29/21655.38 << ',' << 21424.28/21655.38 << ','
//            << 20585.09/21655.38 << ',' << 19536.68/21655.38 << '\n';
//        default:
//            break;
//    }
    
    
    out.close();
    
    std::ofstream outfile(prefix.substr(0,3)+"_kristoff.txt");
    outfile <<
    "All Events,$Z\\rightarrow ee$,$Z\\rightarrow ee$ fiducial,$Z\\rightarrow \\mu\\mu$,$Z\\rightarrow \\mu\\mu$ fiducial\n"
    << allEvents << ',' << nZee << ',' << nZeeFid << ','
    << nZmumu << ',' << nZmumuFid << '\n';
        
    outfile.close();
    
    std::ofstream outfileZ("all_Z_products.txt");
    for (auto& i: all_Z_products)
    {
        outfileZ << i.first << '\t' << i.second << '\n';
    }
    
    outfileZ.close();
    
    for (auto& plot: plots)
    {
        plot.second.save(output_file);
    }

}

void myPreselectionHaa()
{
    auto start_time = Clock::now();
    std::cout << "Run over MC\n";
    
//    std::vector<std::vector<std::string>> input_filenames = {{"/home/common/Haa/ntuples/Za/mc16_13TeV.600909.PhPy8EG_AZNLO_ggH125_mA5p0_Cyy0p01_Czh1p0.merge.AOD.e8324_e7400_s3126_r10724_r10726_v2.root"},
//    {"/home/common/Haa/ntuples/Za/mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v2.root"}};
    
    
//    std::vector<std::vector<std::string>> input_filenames =
//    {{"mc16_13TeV.600909.PhPy8EG_AZNLO_ggH125_mA5p0_Cyy0p01_Czh1p0.merge.AOD.e8324_e7400_s3126_r10724_r10726_v2.root"},{"mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v2.root"}};
    
    std::vector<std::vector<std::string>> input_filenames = {{"mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0_allTruth_Test.root"}};
    
    const char* output_filename = "mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA_p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v2_out.root";
//    const char* output_filename = "mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0_allTruth_Test_out.root";
    
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
    std::ofstream out("someFile.txt",std::ios::app);
    out << "\\textbf{mA1} \\par \\hspace{-4cm} \\scalebox{0.9}{\n";
    out.close();
//    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA1.txt\").rename(index={0:r\"my events\",1:r\"my ratios\",2:r\"paper ratios\"}).style.format(precision=3).to_latex( hrules=True));' >> someFile.txt");
    
    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA1.txt\").rename(index={0:r\"events\"}).style.format(precision=3).to_latex( hrules=True));' >> someFile.txt");
    
    out.open("someFile.txt",std::ios::app);
    out << "}\n";
    out.close();
    
//    out.open("someFile.txt",std::ios::app);
//    out << "\\vspace{2cm} \\textbf{mA5} \\par \\hspace{-4cm}\n";
//    out.close();
//    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA5.txt\").rename(index={0:r\"my events\",1:r\"my ratios\",2:r\"paper ratios\"}).style.format(precision=3).to_latex( hrules=True));' >> someFile.txt");
    
    out.open("someFile.txt",std::ios::app);
    out << "\\vspace{2cm} \\textbf{mA1} \\par \\hspace{-4cm}\n";
    out.close();
    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA1_kristoff.txt\").rename(index={0:r\"events\"}).style.format(precision=3).to_latex(hrules=True));' >> someFile.txt");
    
//    out.open("someFile.txt",std::ios::app);
//    out << "\\vspace{2cm} \\textbf{mA5} \\par \\hspace{-4cm}\n";
//    out.close();
//    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA5_kristoff.txt\").rename(index={0:r\"events\"}).style.format(precision=3).to_latex(hrules=True));' >> someFile.txt");
    
    system("cat someFile.txt");
//    system("rm mA1.txt mA5.txt someFile.txt");
    system("rm mA1.txt mA1_kristoff.txt someFile.txt");

    output_file->Close();
    auto end_time = Clock::now();
    std::cout << "Time difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
}

int main()
{
    myPreselectionHaa();
}
