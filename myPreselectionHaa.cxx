#include <unordered_map>
#include <map>
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

template <typename T, typename U>
void filter(std::vector<T>& vec, bool (*func) (T&, U&), U& p)
{
    for (auto i = vec.begin(); i != vec.end(); ++i)
    {
        if (!func(*i, p))
        {
            vec.erase(i);
            i--;
        }
    }
}

template <typename T, typename U>
void filter(std::vector<T>& vec, bool (*func) (T&, U&), U p)
{
    for (auto i = vec.begin(); i != vec.end(); ++i)
    {
        if (!func(*i, p))
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
    const int MaxBins = 200, minBins = 200, inc = 100;
    
    static std::unordered_map<int,double> zeroPhotonsMatchedEff; //prefix[2]-'0'
    static std::unordered_map<int,double> twoPhotonsMatchedEff; //prefix[2]-'0'
    static std::unordered_map<int,double> onePhotonMatchedEff;  //prefix[2]-'0'
    static std::unordered_map<int,double> moreThanTwoPhotonMatchedEff;  //prefix[2]-'0'
    static std::unordered_map<int,double> zeroPhotonsMatchedEffLeptons; //prefix[2]-'0'
    static std::unordered_map<int,double> twoPhotonsMatchedEffLeptons; //prefix[2]-'0'
    static std::unordered_map<int,double> onePhotonMatchedEffLeptons;  //prefix[2]-'0'
    static std::unordered_map<int,double> moreThanTwoPhotonMatchedEffLeptons; //prefix[2]-'0'
    int truthPhotonsFiducialCount = 0;
    int truthPhotonsFiducialCountLepton = 0;
    
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
        
        plots.emplace(prefix+"dilep + photon mass before dilep and photon cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep + photon mass before dilep and photon cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep + photon mass before dilep and photon cuts electrons objects, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"dilep + photon mass after dilep cuts and before photon cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep + photon mass after dilep cuts and before photon cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep + photon mass after dilep cuts and before photon cuts electrons objects, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"dilep + photon mass after dilep cuts and photon cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep + photon mass after dilep cuts and photon cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep + photon mass after dilep cuts and photon cuts electrons objects, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"pt distribution of all electron objects, bins = " + std::to_string(i),
            Plot(prefix+"pt distribution of all electron objects, bins = " + std::to_string(i),
                 prefix+"pt distribution of all electron objects, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"eta distribution of all electron objects, bins = " + std::to_string(i),
            Plot(prefix+"eta distribution of all electron objects, bins = " + std::to_string(i),
                 prefix+"eta distribution of all electron objects, bins = " + std::to_string(i), i, -7.5, 7.5));
        
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
        
        plots.emplace(prefix+"dilep mass after all dilep and photon cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep mass after all dilep and photon cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep mass after all dilep and photon cuts electrons objects, bins = " + std::to_string(i), i, 60, 120));
    
        plots.emplace(prefix+"dilep delta R after all dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep delta R after all dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep delta R after all dilep cuts electrons objects, bins = " + std::to_string(i), i, 0, 6.5));
        
        plots.emplace(prefix+"dilep delta eta after all dilep cuts electrons objects, bins = " + std::to_string(i),
            Plot(prefix+"dilep delta eta after all dilep cuts electrons objects, bins = " + std::to_string(i),
                 prefix+"dilep delta eta after all dilep cuts electrons objects, bins = " + std::to_string(i), i, 0, 6.25));
        
        plots.emplace(prefix+"pt distribution of all electron objects after all dilep cuts, bins = " + std::to_string(i),
            Plot(prefix+"pt distribution of all electron objects after all dilep cuts, bins = " + std::to_string(i),
                 prefix+"pt distribution of all electron objects after all dilep cuts, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"eta distribution of all electron objects after all dilep cuts, bins = " + std::to_string(i),
            Plot(prefix+"eta distribution of all electron objects after all dilep cuts, bins = " + std::to_string(i),
                 prefix+"eta distribution of all electron objects after all dilep cuts, bins = " + std::to_string(i), i, -7.5, 7.5));
        
        plots.emplace(prefix+"pt distribution of all photon objects, bins = " + std::to_string(i),
            Plot(prefix+"pt distribution of all photon objects, bins = " + std::to_string(i),
                 prefix+"pt distribution of all photon objects, bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"eta distribution of all photon objects, bins = " + std::to_string(i),
            Plot(prefix+"eta distribution of all photon objects, bins = " + std::to_string(i),
                 prefix+"eta distribution of all photon objects, bins = " + std::to_string(i), i, -7.5, 7.5));
        
        plots.emplace(prefix+"delta R of reco photons with truth axion, bins = " + std::to_string(i),
            Plot(prefix+"delta R of reco photons with truth axion, bins = " + std::to_string(i),
                 prefix+"delta R of reco photons with truth axion, bins = " + std::to_string(i), i, 0, 6.5));
        
        plots.emplace(prefix+"delta R of reco-photon pairs with #Delta R_{#gamma_{1,2},a} < 0.2 and/or min#left(#Delta R_{#gamma_{1,2},a}#right), bins = " + std::to_string(i),
            Plot(prefix+"delta R of reco-photon pairs with #Delta R_{#gamma_{1,2},a} < 0.2 and/or min#left(#Delta R_{#gamma_{1,2},a}#right), bins = " + std::to_string(i),
                 prefix+"delta R of reco-photon pairs with #Delta R_{#gamma_{1,2},a} < 0.2 and/or min#left(#Delta R_{#gamma_{1,2},a}#right), bins = " + std::to_string(i), i, 0, 0.2));
        
        plots.emplace(prefix+"delta R of reco-photon pairs with #Delta R_{#gamma_{1,2},a} < 0.2 and/or min#left(#Delta R_{#gamma_{1,2},a}#right) after cuts, bins = " + std::to_string(i),
            Plot(prefix+"delta R of reco-photon pairs with #Delta R_{#gamma_{1,2},a} < 0.2 and/or min#left(#Delta R_{#gamma_{1,2},a}#right) after cuts, bins = " + std::to_string(i),
                 prefix+"delta R of reco-photon pairs with #Delta R_{#gamma_{1,2},a} < 0.2 and/or min#left(#Delta R_{#gamma_{1,2},a}#right) after cuts, bins = " + std::to_string(i), i, 0, 0.2));
        
        plots.emplace(prefix+"delta R of photon 1 and photon 2 after all cuts, bins = " + std::to_string(i),
            Plot(prefix+"delta R of photon 1 and photon 2 after all cuts, bins = " + std::to_string(i),
                 prefix+"delta R of photon 1 and photon 2 after all cuts, bins = " + std::to_string(i), i, 0, 0.2));
        
        plots.emplace(prefix+"Actual delta R of photon 1 and photon 2 decayed from axion, bins = " + std::to_string(i),
            Plot(prefix+"Actual delta R of photon 1 and photon 2 decayed from axion, bins = " + std::to_string(i),
                 prefix+"Actual delta R of photon 1 and photon 2 decayed from axion, bins = " + std::to_string(i), i, 0, 6.5));
        
        plots.emplace(prefix+"delta R of all reco photons with two truth photons that decayed from axion, bins = " + std::to_string(i),
            Plot(prefix+"delta R of all reco photons with two truth photons that decayed from axion, bins = " + std::to_string(i),
                 prefix+"delta R of all reco photons with two truth photons that decayed from axion, bins = " + std::to_string(i), i, 0, 1.75));
        
        plots.emplace(prefix+"P_{t}#left(#gamma_{1,a}^{t} + #gamma_{2,a}^{t}#right), bins = " + std::to_string(i),
            Plot(prefix+"P_{t}#left(#gamma_{1,a}^{t} + #gamma_{2,a}^{t}#right), bins = " + std::to_string(i),
                 prefix+"P_{t}#left(#gamma_{1,a}^{t} + #gamma_{2,a}^{t}#right), bins = " + std::to_string(i), i, 0, 200));
        
        plots.emplace(prefix+"P_{t}#left(#gamma_{1,a}^{t} + #gamma_{2,a}^{t} - #gamma_{reco}#right), bins = " + std::to_string(i),
            Plot(prefix+"P_{t}#left(#gamma_{1,a}^{t} + #gamma_{2,a}^{t} - #gamma_{reco}#right), bins = " + std::to_string(i),
                 prefix+"P_{t}#left(#gamma_{1,a}^{t} + #gamma_{2,a}^{t} - #gamma_{reco}#right), bins = " + std::to_string(i), i, 0, 200));
        
        
    }
    Event::cache_truth = false;
    
    FileReaderRange reader(input_filenames);
    
    Event::systematic = systematic;
    
    int weight = 1;
    
//    std::unordered_set<std::string> all_triggers;
//    std::map <int,int> all_pdg_ids;
//    std::unordered_map <int,int> all_Z_products;
    
    int beforePreselection = 0,
        two_leptons = 0,
        opp_charge = 0,
        lep1_lep2_pt = 0,
        lep_same_flavor = 0,
        dilep_mass = 0,
        dilep_pt = 0;
    
    int events_two_photons_in_direction_of_axion = 0, events_one_photons_in_direction_of_axion = 0;
    
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
        if (truthChain.size() >= 1)
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
    
    auto photonObjCuts = [](Photon& p)
    {
        return ((p.pt()/1e3 > 10) && (p.id_loose()==1) && (abs(p.eta()) < 2.37));
    };
//    int diphotonPairs=0, pdg_id35 = 0;
        
    auto truthPhotonsFiducial = [](std::vector<TruthParticle>& tps)
    {
        for (auto& i: tps)
        {
            if ((abs(i.eta()) >= 2.37) || (i.pt() <= 10))
            {
                return false;
            }
        }
        return true;
    };
    
    auto recoPhotonsFiducial = [](std::vector<Photon>& tps)
    {
        for (auto& i: tps)
        {
            if ((abs(i.eta()) >= 2.37) || (i.pt() <= 10) || (abs(i.eta()) > 1.37
                && abs(i.eta()) < 1.52))
            {
                return false;
            }
        }
        return true;
    };
    
    for (auto &&f: reader)
    {
        std::cout << "entry_number " << f.__current_event.entry_number  << '\n';
        //        std::copy(f.__current_event.triggers.begin(),f.__current_event.triggers.end(),std::inserter(all_triggers,all_triggers.end()));
        //        std::cout<<'\n';
        bool photons_from_axion_found;
        bool lepton_passed = false;
        
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
            //Electrons
            std::vector<Electron>& reco_electrons = f.__current_event.electrons;
            std::vector<Photon>& reco_photons = f.__current_event.photons;
            
            for (int i=minBins; i<=MaxBins; i+=inc)
            {
                for (auto& j: reco_electrons)
                {
                    plots.at(prefix+"pt distribution of all electron objects, bins = " + std::to_string(i)).fill(j.pt()/1e3,weight);
                    plots.at(prefix+"eta distribution of all electron objects, bins = " + std::to_string(i)).fill(j.eta(),weight);
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
                        if (reco_photons.size() == 2)
                        {
                            CandidateSet<Photon> candidateDiPhoton(std::make_pair(reco_photons[0],reco_photons[1]));
                            
                            plots.at(prefix+"dilep + photon mass before dilep and photon cuts electrons objects, bins = " + std::to_string(i)).fill((candidateDiElectron.four_momentum+candidateDiPhoton.four_momentum).M()/1e3,weight);
    //                        dilep + photon mass before dilep and photon cuts electrons objects, bins =
                        }
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
                                        lepton_passed=true;
                                        for (int i=minBins; i<=MaxBins; i+=inc)
                                        {
                                            plots.at(prefix+"dilep pt after all dilep cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.four_momentum.Pt()/1e3,weight);
                                            plots.at(prefix+"dilep mass after all dilep cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.four_momentum.M()/1e3,weight);
                                            plots.at(prefix+"dilep delta R after all dilep cuts electrons objects, bins = " + std::to_string(i)).fill(reco_electrons[0].delta_r(reco_electrons[1]),weight);
                                            plots.at(prefix+"dilep delta eta after all dilep cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.delta_eta(),weight);
                                            if (reco_photons.size() == 2)
                                            {
                                                // dilep + photon mass after dilep cuts and before photon cuts electrons objects, bins =
                                                CandidateSet<Photon> candidateDiPhoton(std::make_pair(reco_photons[0],reco_photons[1]));
                                                plots.at(prefix+"dilep + photon mass after dilep cuts and before photon cuts electrons objects, bins = " + std::to_string(i)).fill((candidateDiElectron.four_momentum+candidateDiPhoton.four_momentum).M()/1e3,weight);
                                            }
                                        }
                                        for (int i=minBins; i<=MaxBins; i+=inc)
                                        {
                                            for (auto& j: reco_electrons)
                                            {
                                                plots.at(prefix+"pt distribution of all electron objects after all dilep cuts, bins = " + std::to_string(i)).fill(j.pt()/1e3,weight);
                                                plots.at(prefix+"eta distribution of all electron objects after all dilep cuts, bins = " + std::to_string(i)).fill(j.eta(),weight);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            //Photons
//            std::vector<Photon> reco_photons = f.__current_event.photons;
            //            std::vector<TruthParticle> truth_axions = f.find_truth_particles({},{},{35});
            
            for (int i=minBins; i<=MaxBins; i+=inc)
            {
                for (auto& j: reco_photons)
                {
                    plots.at(prefix+"pt distribution of all photon objects, bins = " + std::to_string(i)).fill(j.pt()/1e3,weight);
                    plots.at(prefix+"eta distribution of all photon objects, bins = " + std::to_string(i)).fill(j.eta(),weight);
                }
                filter<Photon>(reco_photons,photonObjCuts);
                if (reco_photons.size()==2 && lepton_passed)
                {
                    plots.at(prefix+ "delta R of photon 1 and photon 2 after all cuts, bins = " + std::to_string(i)).fill(reco_photons[1].delta_r(reco_photons[0]),weight);
                    CandidateSet<Electron> candidateDiElectron(std::make_pair(reco_electrons[0],reco_electrons[1]));
                    
                    plots.at(prefix+ "dilep mass after all dilep and photon cuts electrons objects, bins = " + std::to_string(i)).fill(candidateDiElectron.four_momentum.M()/1e3,weight);
                    
                    //dilep + photon mass after dilep cuts and photon cuts electrons objects, bins =
                    CandidateSet<Photon> candidateDiPhoton(std::make_pair(reco_photons[0],reco_photons[1]));
                    plots.at(prefix+"dilep + photon mass after dilep cuts and photon cuts electrons objects, bins = " + std::to_string(i)).fill((candidateDiElectron.four_momentum+candidateDiPhoton.four_momentum).M()/1e3,weight);
                }
            }
        }
    }
    
    std::cout << "\n\n\n";
    
    std::ofstream out(prefix.substr(0,3)+".txt");
    
//    out <<
//    "TwoLeptons,OppCharge,Lep1Lep2Pt,LepSameFlavor,DilepMass,DilepPt\n"
    
    out << "allEvents,$\\geq 2 e^-$,$2e^-$,Opp Charge, $p_T^{\\mathrm{Leading}} > 27$ GeV,$\\Delta R > 0.01$,81 GeV $\\leq m_{ll} \\leq$ 101 GeV,$p_T^{ll}>$ 10 GeV\n"
    
    << allEvents << ',' << elecObj << ',' << exactlyTwo << ',' <<
    oppCharge << ',' << leadingPt << ',' << dR << ',' << mLL << ',' << pTLL << '\n';
    
    
    out.close();
    
    std::ofstream outfile(prefix.substr(0,3)+"_kristoff.txt");
    outfile <<
    "All Events,$Z\\rightarrow ee$,$Z\\rightarrow ee$ fiducial,$Z\\rightarrow \\mu\\mu$,$Z\\rightarrow \\mu\\mu$ fiducial\n"
    << allEvents << ',' << nZee << ',' << nZeeFid << ','
    << nZmumu << ',' << nZmumuFid << '\n';
        
    outfile.close();
    
    std::cout << "Events with atleast two reco-photons within ΔR < 0.2 of axion = " << events_two_photons_in_direction_of_axion
    << "\nEvents with exactly one reco-photon within ΔR < 0.2 of axion = " << events_one_photons_in_direction_of_axion << '\n';
    
    twoPhotonsMatchedEff[prefix[2]-'0'] /= (1.0*truthPhotonsFiducialCount); //prefix[2]-'0'
    onePhotonMatchedEff[prefix[2]-'0'] /= (1.0*truthPhotonsFiducialCount);
    moreThanTwoPhotonMatchedEff[prefix[2]-'0'] /= (1.0*truthPhotonsFiducialCount);
    zeroPhotonsMatchedEff[prefix[2]-'0'] /= (1.0*truthPhotonsFiducialCount);
    
    std::cout << "two photon reco-efficiency\n" << std::string(26,'=') << '\n';
    for (auto& i: twoPhotonsMatchedEff)
    {
        std::cout << i.first << " GeV \t" << i.second << '\n';
    }
    
    std::cout << '\n';
    
    std::cout << "one photon reco-efficiency\n" << std::string(26,'=') << '\n';
    for (auto& i: onePhotonMatchedEff)
    {
        std::cout << i.first << " GeV \t" << i.second << '\n';
    }
    
    std::cout << '\n';
    
    std::cout << "more than two photon reco-efficiency\n" << std::string(26,'=') << '\n';
    for (auto& i: moreThanTwoPhotonMatchedEff)
    {
        std::cout << i.first << " GeV \t" << i.second << '\n';
    }
    
    std::cout << '\n';
    
    std::cout << "fraction of signal-events within acceptance with no reco photons\n" << std::string(26,'=') << '\n';
    for (auto& i: zeroPhotonsMatchedEff)
    {
        std::cout << i.first << " GeV \t" << i.second << '\n';
    }
    
    std::cout << "\n\n";
    
    twoPhotonsMatchedEffLeptons[prefix[2]-'0'] /= (1.0*truthPhotonsFiducialCountLepton); //prefix[2]-'0'
    onePhotonMatchedEffLeptons[prefix[2]-'0'] /= (1.0*truthPhotonsFiducialCountLepton);
    moreThanTwoPhotonMatchedEffLeptons[prefix[2]-'0'] /= (1.0*truthPhotonsFiducialCountLepton);
    zeroPhotonsMatchedEffLeptons[prefix[2]-'0'] /= (1.0*truthPhotonsFiducialCountLepton);
    
    std::cout << "two photon reco-efficiency with leptons\n" << std::string(26,'=') << '\n';
    for (auto& i: twoPhotonsMatchedEffLeptons)
    {
        std::cout << i.first << " GeV \t" << i.second << '\n';
    }
    
    std::cout << '\n';
    
    std::cout << "one photon reco-efficiency with leptons\n" << std::string(26,'=') << '\n';
    for (auto& i: onePhotonMatchedEffLeptons)
    {
        std::cout << i.first << " GeV \t" << i.second << '\n';
    }
    
    std::cout << '\n';
    
    std::cout << "more than two photon reco-efficiency with leptons\n" << std::string(26,'=') << '\n';
    for (auto& i: moreThanTwoPhotonMatchedEffLeptons)
    {
        std::cout << i.first << " GeV \t" << i.second << '\n';
    }
    
    std::cout << '\n';
    
    std::cout << "fraction of signal-events within acceptance & w/ leptons with no reco photons\n" << std::string(26,'=') << '\n';
    for (auto& i: zeroPhotonsMatchedEffLeptons)
    {
        std::cout << i.first << " GeV \t" << i.second << '\n';
    }
    
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
    
    std::vector<std::vector<std::string>> input_filenames = {
//        {"mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA1p0_Cyy0p01_Czh1p0_allTruth_Test.root"},
        {"Ntuple_data_test.root"}
//        {"Ntuple_MC_Za_mA5p0_v4.root"}
        
    };
    
    const char* output_filename =
//    "mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA_p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v2_out.root";
    "Ntuple_data_test_out.root";
//    "Ntuple_MC_Za_mA5p0_v4_out.root";
    
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
//    std::ofstream out("someFile.txt",std::ios::app);
//    out << "\\textbf{mA1} \\par \\hspace{-4cm} \\scalebox{0.9}{\n";
//    out.close();
//    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA1.txt\").rename(index={0:r\"my events\",1:r\"my ratios\",2:r\"paper ratios\"}).style.format(precision=3).to_latex( hrules=True));' >> someFile.txt");
//
//    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA1.txt\").rename(index={0:r\"events\"}).style.format(precision=3).to_latex( hrules=True));' >> someFile.txt");
//
//    out.open("someFile.txt",std::ios::app);
//    out << "}\n\n";
//    out.close();
//
//    out.open("someFile.txt",std::ios::app);
//    out << "\\vspace{2cm} \\textbf{mA5} \\par \\hspace{-4cm}\n";
//    out.close();
//    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA5.txt\").rename(index={0:r\"my events\",1:r\"my ratios\",2:r\"paper ratios\"}).style.format(precision=3).to_latex( hrules=True));' >> someFile.txt");
//
//    out.open("someFile.txt",std::ios::app);
//    out << "\\vspace{2cm} \\textbf{mA1} \\par \\hspace{-4cm}\n";
//    out.close();
//    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA1_kristoff.txt\").rename(index={0:r\"events\"}).style.format(precision=3).to_latex(hrules=True));' >> someFile.txt");
//
//    out.open("someFile.txt",std::ios::app);
//    out << "\\vspace{2cm} \\textbf{mA5} \\par \\hspace{-4cm}\n";
//    out.close();
//    system("python3 -c 'import pandas as pd; print(pd.read_csv(r\"mA5_kristoff.txt\").rename(index={0:r\"events\"}).style.format(precision=3).to_latex(hrules=True));' >> someFile.txt");
//
//    system("cat someFile.txt");
//    system("rm mA1.txt mA5.txt someFile.txt");
//    system("rm mA1.txt mA1_kristoff.txt someFile.txt");

    output_file->Close();
    auto end_time = Clock::now();
    std::cout << "Time difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
}

int main()
{
    myPreselectionHaa();
}


//TODO: Match two reco-photons to two-truth photons
