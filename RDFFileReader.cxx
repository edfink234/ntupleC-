#include "ROOT/RDataFrame.hxx"
#include <ROOT/RLogger.hxx>
#include "TCanvas.h"
#include "TMath.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#include <vector>
#include <string>
#include <memory>
#include <tuple>

#include "RDFFileReader.h"
#include "RDFevent.h"

using namespace ROOT::VecOps;
using namespace ROOT::Math;
using namespace ROOT::Math::VectorUtil;

RDFFileReader::RDFFileReader(const std::vector<std::string>& files, const char* tree_name, unsigned int numThreads) :

__has_event_info_chain{true},
setThreads{numThreads},
__chain{tree_name},
__event_info_chain{"full_event_info"}
//df_chain(tree_name, files),
//df_event_info_chain("full_event_info", files)
{
    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);

    for (const auto& f: files)
    {
        __chain.Add(f.c_str());
        __event_info_chain.Add(f.c_str());
    }
    
    __chain.AddFriend(&__event_info_chain);
    
    df = std::make_unique<ROOT::RDataFrame>(__chain);
    
//    df->Describe().Print();
//    auto filtered_df = df.Filter([](RVec<int> pdg_ids, RVec<int> pdg_ids status_codes)
//    {
//        return std::
//    })
    
//    auto two = df->Define("two",[](){return 2;}, );
    
//    std::cout << two.Display({"two"},10)->AsString() << '\n';
    
//    df->ForeachSlot(
//    [&](unsigned int s, RVec<float>& mc_pt)
//    {
//        std::cout << mc_pt[0] << '\n';
//    } ,{"mc_pt"});
    
//    mc_barcode                                      ROOT::VecOps::RVec<int>                 Dataset
//    mc_charge                                       ROOT::VecOps::RVec<float>               Dataset
//    mc_decay_length                                 ROOT::VecOps::RVec<float>               Dataset
//    mc_decay_time                                   ROOT::VecOps::RVec<float>               Dataset
//    mc_e                                            ROOT::VecOps::RVec<float>               Dataset
//    mc_eta                                          ROOT::VecOps::RVec<float>               Dataset
//    mc_mass                                         ROOT::VecOps::RVec<float>               Dataset
//    mc_parent_barcode                               ROOT::VecOps::RVec<int>                 Dataset
//    mc_parent_pdg_id                                ROOT::VecOps::RVec<int>                 Dataset
//    mc_pdg_id                                       ROOT::VecOps::RVec<int>                 Dataset
//    mc_phi                                          ROOT::VecOps::RVec<float>               Dataset
//    mc_production_vertex_barcode                    ROOT::VecOps::RVec<int>                 Dataset
//    mc_production_vertex_status                     ROOT::VecOps::RVec<int>                 Dataset
//    mc_pt                                           ROOT::VecOps::RVec<float>               Dataset
//    mc_status                                       ROOT::VecOps::RVec<int>                 Dataset
//    mc_vertex_barcode                               ROOT::VecOps::RVec<int>                 Dataset
//    mc_vertex_incoming_barcodes                     ROOT::VecOps::RVec<vector<int>>         Dataset
//    mc_vertex_outgoing_barcodes                     ROOT::VecOps::RVec<vector<int>>         Dataset
//    mc_vertex_x                                     ROOT::VecOps::RVec<float>               Dataset
//    mc_vertex_y                                     ROOT::VecOps::RVec<float>               Dataset
//    mc_vertex_z                                     ROOT::VecOps::RVec<float>               Dataset
    
    //TODO: make tuple a struct so data members have names
    auto Tps = df->Define("truth_particles",[&](RVec<int>& mc_pdg_id_, RVec<int>& mc_status_, RVec<float>& mc_pt_, RVec<float>& mc_charge_, RVec<float>& mc_eta_, RVec<float>& mc_phi_, RVec<float>& mc_e_)
    {
        RVec<std::tuple<int, int, float, float, float, float, float>> x;
        x.reserve(mc_pt_.size());
        for (size_t i = 0; i < mc_pt_.size(); i++)
        {
            x.push_back(std::make_tuple(mc_pdg_id_[i],mc_status_[i],mc_pt_[i],mc_charge_[i],mc_eta_[i],mc_phi_[i],mc_e_[i]));
        }
        return x;

    }, {"mc_pdg_id", "mc_status", "mc_pt", "mc_charge", "mc_eta", "mc_phi", "mc_e"});
    
    //TODO: std::remove_if, "erase-remove idiom wiki"
    
    auto stable_truth_photons = Tps.Define("stable_truth_photons", [&](RVec<std::tuple<int, int, float, float, float, float, float>> truth_particles)
    {
        for (auto i = truth_particles.begin(); i != truth_particles.end(); ++i)
        {
            if (abs(std::get<0>(*i))!=22 || std::get<1>(*i)!=1)
            {
                //remove b/f erase???
                truth_particles.erase(i);
                i--;
            }
        }
        
        return truth_particles;

    }, {"truth_particles"});
    
    auto stable_truth_photons_pt = stable_truth_photons.Define("stable_truth_photons_pt", [&](RVec<std::tuple<int, int, float, float, float, float, float>>& stable_truth_photons)
    {
        RVec<float> pt;
        pt.reserve(stable_truth_photons.size());
        for (auto &i: stable_truth_photons)
        {
            pt.push_back(std::get<2>(i)/1e3);
        }
        return pt;

    }, {"stable_truth_photons"});
    
    auto stable_truth_photons_eta = stable_truth_photons.Define("stable_truth_photons_eta", [&](RVec<std::tuple<int, int, float, float, float, float, float>> stable_truth_photons)
    {
        RVec<float> eta;
        eta.reserve(stable_truth_photons.size());
        for (auto &i: stable_truth_photons)
        {
            eta.push_back(std::get<4>(i));
        }
        return eta;

    }, {"stable_truth_photons"});
    
    // Here we filter the truth_particles column in each event
    auto stable_truth_leptons = Tps.Define("stable_truth_leptons", [&](RVec<std::tuple<int, int, float, float, float, float, float>> truth_particles)
    {
        for (auto i = truth_particles.begin(); i != truth_particles.end(); ++i)
        {
            if ( ( (abs(std::get<0>(*i))!=11)
                && (abs(std::get<0>(*i))!=12)
                && (abs(std::get<0>(*i))!=13)
                && (abs(std::get<0>(*i))!=14)
                && (abs(std::get<0>(*i))!=15)
                && (abs(std::get<0>(*i))!=16)
                && (abs(std::get<0>(*i))!=17)
                && (abs(std::get<0>(*i))!=18))
                
                
                || std::get<1>(*i)!=1)
            {
                truth_particles.erase(i);
                i--;
            }
        }
        
        return truth_particles;

    }, {"truth_particles"});
    
//    Tps.Describe().Print();
//    std::cout << stable_truth_leptons.Display({"stable_truth_leptons"},100)->AsString();
//    stable_truth_leptons
    
//    {"mc_pdg_id", "mc_status", "mc_pt", "mc_charge", "mc_eta", "mc_phi", "mc_e"}
//          0            1          2          3          4          5        6
    
//    Now we actually "remove" certain events
    auto passed_truth_leptons = stable_truth_leptons.Filter(
    [&](RVec<std::tuple<int, int, float, float, float, float, float>> stable_truth_leptons)
    {
        return
                (
                    (stable_truth_leptons.size()==2)
                        &&
                    (DeltaR(
                            PtEtaPhiEVector(std::get<2>(stable_truth_leptons[0]),
                                            std::get<4>(stable_truth_leptons[0]),
                                            std::get<5>(stable_truth_leptons[0]),
                                            std::get<6>(stable_truth_leptons[0])),

                            PtEtaPhiEVector(std::get<2>(stable_truth_leptons[1]),
                                            std::get<4>(stable_truth_leptons[1]),
                                            std::get<5>(stable_truth_leptons[1]),
                                            std::get<6>(stable_truth_leptons[1]))
                           ) > 0.01)
                        &&
                    (std::get<3>(stable_truth_leptons[0]) == -1*std::get<3>(stable_truth_leptons[1]))

                        &&

                    ((((std::get<6>(stable_truth_leptons[0]) + std::get<6>(stable_truth_leptons[1]))
                    * (std::get<6>(stable_truth_leptons[0]) + std::get<6>(stable_truth_leptons[1])))
                    -
                    ((std::get<2>(stable_truth_leptons[0]) + std::get<2>(stable_truth_leptons[1]))
                    * (std::get<2>(stable_truth_leptons[0]) + std::get<2>(stable_truth_leptons[1]))))
                    >= 6561e6)

                        &&

                    ((((std::get<6>(stable_truth_leptons[0]) + std::get<6>(stable_truth_leptons[1]))
                    * (std::get<6>(stable_truth_leptons[0]) + std::get<6>(stable_truth_leptons[1])))
                    -
                    ((std::get<2>(stable_truth_leptons[0]) + std::get<2>(stable_truth_leptons[1]))
                    * (std::get<2>(stable_truth_leptons[0]) + std::get<2>(stable_truth_leptons[1]))))
                    <= 10201e6)

                        &&

                    (((std::get<2>(stable_truth_leptons[0]) > 20e3) && (std::get<2>(stable_truth_leptons[1]) > 27e3))
                        ||
                    ((std::get<2>(stable_truth_leptons[1]) > 20e3) && (std::get<2>(stable_truth_leptons[0]) > 27e3)))
////
                 );
        
    } , {"stable_truth_leptons"}, "passed_truth_leptons");
    
//    auto nEntriesAfterCuts = passed_truth_leptons.Count();
//    std::cout << *nEntriesAfterCuts << '\n';
         
//    std::cout << passed_truth_leptons.GetColumnNames().size();
    
//    for (auto &i: passed_truth_leptons.GetColumnNames())
//    {
//        std::cout << i << '\n';
//    }
    
    auto passed_truth_leptons_pt = passed_truth_leptons.Define("passed_truth_leptons_pt", [&](RVec<std::tuple<int, int, float, float, float, float, float>> passed_truth_leptons)
    {
        RVec<float> pt;
        pt.reserve(passed_truth_leptons.size());
        for (auto &i: passed_truth_leptons)
        {
            pt.push_back(std::get<2>(i)/1e3);
        }
        return pt;

    }, {"stable_truth_leptons"});

    auto passed_truth_leptons_eta = passed_truth_leptons.Define("passed_truth_leptons_eta", [&](RVec<std::tuple<int, int, float, float, float, float, float>> passed_truth_leptons)
    {
        RVec<float> eta;
        eta.reserve(passed_truth_leptons.size());
        for (auto &i: passed_truth_leptons)
        {
            eta.push_back(std::get<4>(i));
        }
        return eta;

    }, {"stable_truth_leptons"});

//    photon_cluster_BE1_R                            ROOT::VecOps::RVec<float>                                               Dataset
//    photon_cluster_BE1_Z                            ROOT::VecOps::RVec<float>                                               Dataset
//    photon_cluster_BE1_eta                          ROOT::VecOps::RVec<float>                                               Dataset
//    photon_cluster_e                                ROOT::VecOps::RVec<float>                                               Dataset
//    photon_cluster_eta                              ROOT::VecOps::RVec<float>                                               Dataset
//    photon_cluster_eta_be_2                         ROOT::VecOps::RVec<float>                                               Dataset
//    photon_cluster_phi                              ROOT::VecOps::RVec<float>                                               Dataset
//    photon_cluster_pt                               ROOT::VecOps::RVec<float>                                               Dataset
//    photon_e                                        ROOT::VecOps::RVec<float>                                               Dataset
//    photon_eta                                      ROOT::VecOps::RVec<float>                                               Dataset
//    photon_etcone40                                 ROOT::VecOps::RVec<float>                                               Dataset
//    photon_id                                       ROOT::VecOps::RVec<int>                                                 Dataset
//    photon_id_eff                                   ROOT::VecOps::RVec<float>                                               Dataset
//    photon_id_loose                                 ROOT::VecOps::RVec<int>                                                 Dataset
//    photon_id_tight                                 ROOT::VecOps::RVec<int>                                                 Dataset
//    photon_is_crack                                 ROOT::VecOps::RVec<int>                                                 Dataset
//    photon_iso_eff                                  ROOT::VecOps::RVec<float>                                               Dataset
//    photon_mc_barcode                               ROOT::VecOps::RVec<int>                                                 Dataset
//    photon_mc_origin                                ROOT::VecOps::RVec<int>                                                 Dataset
//    photon_phi                                      ROOT::VecOps::RVec<float>                                               Dataset
//    photon_pt                                       ROOT::VecOps::RVec<float>                                               Dataset
//    photon_pt_uncorrected                           ROOT::VecOps::RVec<float>                                               Dataset
//    photon_ptcone20                                 ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_delta_e                     ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_e_ratio                     ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_f_1                         ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_f_side                      ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_r_eta                       ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_r_phi                       ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_rhad                        ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_rhad1                       ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_w_eta_2                     ROOT::VecOps::RVec<float>                                               Dataset
//    photon_shower_shape_w_s_tot                     ROOT::VecOps::RVec<float>                                               Dataset
//    photon_syst_e                                   ROOT::VecOps::RVec<vector<float>>                                       Dataset
//    photon_syst_id_eff                              ROOT::VecOps::RVec<vector<float>>                                       Dataset
//    photon_syst_iso_eff                             ROOT::VecOps::RVec<vector<float>>                                       Dataset
//    photon_syst_name                                ROOT::VecOps::RVec<vector<string>>                                      Dataset
//    photon_syst_pt                                  ROOT::VecOps::RVec<vector<float>>                                       Dataset
//    photon_syst_trg_eff                             ROOT::VecOps::RVec<vector<float>>                                       Dataset
//    photon_topoetcone20                             ROOT::VecOps::RVec<float>                                               Dataset
//    photon_topoetcone40                             ROOT::VecOps::RVec<float>                                               Dataset
//    photon_trg_eff                                  ROOT::VecOps::RVec<float>                                               Dataset
    
    
    // photon_id, photon_pt, photon_eta
    auto photons = df->Define("photons",[&](RVec<int>& photon_id_, RVec<float>& photon_pt_, RVec<float>& photon_eta_)
    {
        RVec<std::tuple<int, float, float>> x;
        x.reserve(photon_pt_.size());
        for (size_t i = 0; i < photon_pt_.size(); i++)
        {
            if (photon_pt_[i] < 0)
            {
                continue;
            }
            x.push_back(std::make_tuple(photon_id_[i],photon_pt_[i],photon_eta_[i]));
        }
        return x;

    }, {"photon_id", "photon_pt", "photon_eta"});
    
    auto selected_photons = photons.Define("selected_photons", [&](RVec<std::tuple<int, float, float>> photons)
    {
        for (auto i = photons.begin(); i != photons.end(); ++i)
        {
            if ( (!std::get<0>(*i))|| (std::get<1>(*i)<=10000) || (std::get<2>(*i) >= 2.37)
                || (abs(std::get<2>(*i))>1.37 && abs(std::get<2>(*i))<1.52)
                )
            {
                photons.erase(i);
                i--;
            }
        }
        
        return photons;

    }, {"photons"});
    
    auto selected_photons_pt = selected_photons.Define("selected_photons_pt", [&](RVec<std::tuple<int, float, float>>& selected_photons)
    {
        RVec<float> pt;
        pt.reserve(selected_photons.size());
        for (auto &i: selected_photons)
        {
            pt.push_back(std::get<1>(i)/1e3);
        }
        return pt;

    }, {"selected_photons"});

    auto selected_photons_eta = selected_photons.Define("selected_photons_eta", [&](RVec<std::tuple<int, float, float>>& selected_photons)
    {
        RVec<float> eta;
        eta.reserve(selected_photons.size());
        for (auto &i: selected_photons)
        {
            eta.push_back(std::get<2>(i));
        }
        return eta;

    }, {"selected_photons"});

    auto passed_selected_photons = selected_photons.Filter(
    [&](RVec<std::tuple<int, float, float>>& selected_photons)
    {
        return (selected_photons.size()==2);
        
    } , {"selected_photons"}, "passed_selected_photons");
    
    auto passed_selected_photons_pt = passed_selected_photons.Define("passed_selected_photons_pt", [&](RVec<std::tuple<int, float, float>> passed_selected_photons)
    {
        RVec<float> pt;
        pt.reserve(passed_selected_photons.size());
        for (auto &i: passed_selected_photons)
        {
            pt.push_back(std::get<1>(i)/1e3);
        }
        return pt;

    }, {"selected_photons"});

    auto passed_selected_photons_eta = passed_selected_photons.Define("passed_selected_photons_eta", [&](RVec<std::tuple<int, float, float>> passed_selected_photons)
    {
        RVec<float> eta;
        eta.reserve(passed_selected_photons.size());
        for (auto &i: passed_selected_photons)
        {
            eta.push_back(std::get<2>(i));
        }
        return eta;

    }, {"selected_photons"});
    
    
    
    auto nEntriesAfterCuts = passed_selected_photons_eta.Count();
    std::cout << "# events for nominal = " << *nEntriesAfterCuts << '\n';
    
    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos;
    histos.reserve(10);
    
    histos = {stable_truth_photons_pt.Histo1D<RVec<float>>({"histName1", "histTitle", 20u, 0, 200}, "stable_truth_photons_pt"),
        stable_truth_photons_eta.Histo1D<RVec<float>>({"histName2", "histTitle", 40u, -6, 6}, "stable_truth_photons_eta"),
        stable_truth_photons_pt.Histo1D<RVec<float>>({"histName3", "histTitle", 100u, 0, 100}, "stable_truth_photons_pt"),
        passed_truth_leptons_pt.Histo1D<RVec<float>>({"histName4", "histTitle", 20u, 0, 200}, "passed_truth_leptons_pt"),
        passed_truth_leptons_eta.Histo1D<RVec<float>>({"histName5", "histTitle", 40u, -6, 6}, "passed_truth_leptons_eta"),
        passed_truth_leptons_pt.Histo1D<RVec<float>>({"histName6", "histTitle", 100u, 0, 100}, "passed_truth_leptons_pt"),
        selected_photons_pt.Histo1D<RVec<float>>({"histName7", "histTitle", 20u, 0, 200}, "selected_photons_pt"),
        selected_photons_eta.Histo1D<RVec<float>>({"histName8", "histTitle", 40u, -6, 6}, "selected_photons_eta"),
        passed_selected_photons_pt.Histo1D<RVec<float>>({"histName9", "histTitle", 100u, 0, 100}, "passed_selected_photons_pt"),
        passed_selected_photons_eta.Histo1D<RVec<float>>({"histName10", "histTitle", 40u, -6, 6}, "passed_selected_photons_eta"),
    };
    
    TCanvas *c1;
    for (auto& h: histos)
    {
        c1 = new TCanvas("","",800, 700);
        h->Draw("same");
        c1->SaveAs((h->GetName()+std::string(".png")).c_str());
    }
    
//    TCanvas *c1 = new TCanvas("","",800, 700);
//    auto myHist2 = stable_truth_photons_pt.Histo1D<RVec<float>>({"histName", "histTitle", 20u, 0, 200}, "stable_truth_photons_pt");
//    myHist2->Draw("same");
//    c1->SaveAs("stable_truth_photons_pt.png");
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = stable_truth_photons_eta.Histo1D<RVec<float>>({"histName", "histTitle", 40u, -6, 6}, "stable_truth_photons_eta");
//    myHist2->Draw("same");
//    c1->SaveAs("stable_truth_photons_eta.png");
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = stable_truth_photons_pt.Histo1D<RVec<float>>({"histName", "histTitle", 100u, 0, 100}, "stable_truth_photons_pt");
//    myHist2->Draw("same");
//    c1->SaveAs("stable_truth_photons_pt_tight.png");
//
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = passed_truth_leptons_pt.Histo1D<RVec<float>>({"histName", "histTitle", 20u, 0, 200}, "passed_truth_leptons_pt");
//    myHist2->Draw("same");
//    c1->SaveAs("passed_truth_leptons_pt.png");
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = passed_truth_leptons_eta.Histo1D<RVec<float>>({"histName", "histTitle", 40u, -6, 6}, "passed_truth_leptons_eta");
//    myHist2->Draw("same");
//    c1->SaveAs("passed_truth_leptons_eta.png");
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = passed_truth_leptons_pt.Histo1D<RVec<float>>({"histName", "histTitle", 100u, 0, 100}, "passed_truth_leptons_pt");
//    myHist2->Draw("same");
//    c1->SaveAs("passed_truth_leptons_pt_tight.png");
//
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = selected_photons_pt.Histo1D<RVec<float>>({"histName", "histTitle", 20u, 0, 200}, "selected_photons_pt");
//    myHist2->Draw("same");
//    c1->SaveAs("selected_photons_pt.png");
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = selected_photons_eta.Histo1D<RVec<float>>({"histName", "histTitle", 40u, -6, 6}, "selected_photons_eta");
//    myHist2->Draw("same");
//    c1->SaveAs("selected_photons_eta.png");
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = selected_photons_pt.Histo1D<RVec<float>>({"histName", "histTitle", 100u, 0, 100}, "selected_photons_pt");
//    myHist2->Draw("same");
//    c1->SaveAs("selected_photons_pt_tight.png");
//
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = passed_selected_photons_pt.Histo1D<RVec<float>>({"histName", "histTitle", 20u, 0, 200}, "passed_selected_photons_pt");
//    myHist2->Draw("same");
//    c1->SaveAs("passed_selected_photons_pt.png");
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = passed_selected_photons_eta.Histo1D<RVec<float>>({"histName", "histTitle", 40u, -6, 6}, "passed_selected_photons_eta");
//    myHist2->Draw("same");
//    c1->SaveAs("passed_selected_photons_eta.png");
//    c1 = new TCanvas("","",800, 700);
//    myHist2 = passed_selected_photons_pt.Histo1D<RVec<float>>({"histName", "histTitle", 100u, 0, 100}, "passed_selected_photons_pt");
//    myHist2->Draw("same");
//    c1->SaveAs("passed_selected_photons_pt_tight.png");
//
//
//    auto whatervver = makeRdat(filenames);
    
//    std::cout << df->Display({"calo_cluster_e","calo_cluster_em_probability","vertex_pos_z"},2)->AsString();
    
//    auto myHist2 = myDf.Histo1D<float>("myColumn");
    
//    df->ForeachSlot(
//    [&](unsigned int s, RVec<float>& mc_pt_)
//    {
//        for (auto& i: mc_pt_)
//        {
//            std::cout << i << ' ';
//        }
//        std::cout << '\n';
//    } ,{"mc_pt"});
    
    
//    stable_truth_photons.Display("stable_truth_photons",100)->Print();
    
//    Tps.Describe().Print();
    
    
                               
}

RDFFileReader::~RDFFileReader() = default;

RDFFileReader::__callEnableImplicitMT::__callEnableImplicitMT(unsigned int N)
{
    if (N==-1)
    {
        std::cout << "called\n";
        return;
    }
    std::cout << "not called\n";
    ROOT::EnableImplicitMT(N);
}
