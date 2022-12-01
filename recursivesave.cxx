#include <iostream>
#include <algorithm>

#include "TH1F.h"
#include "TKey.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"

std::string GetXTitle(std::string& title)
{
    if (title.find("pt distribution")!=std::string::npos)
    {
        return "p_{T} (GeV)";
    }
    else if (title.find("eta distribution")!=std::string::npos)
    {
        return "#eta";
    }
    else if (title.find("dilep pt")!=std::string::npos)
    {
        return "p_{T}^{ll} (GeV)";
    }
    else if (title.find("dilep mass")!=std::string::npos)
    {
        return "m_{ll} (GeV)";
    }
    else if (title.find("dilep delta R")!=std::string::npos)
    {
        return "#Delta R_{ll}";
    }
    else if (title.find("dilep delta eta")!=std::string::npos)
    {
        return "#Delta #eta_{ll}";
    }
    else if (title.find("delta R")!=std::string::npos)
    {
        return "#Delta R";
    }
    else if (title.find("P_t")!=std::string::npos)
    {
        return "p_{T} (GeV)";
    }
    return "";
}

void recursiveTH1Fsave(TList* f)
{
    for (auto i: *f)
    {
        if (static_cast<TKey*>(i)->GetClassName()==std::string("TH1F"))
        {
            TCanvas c1;
            std::string str = (static_cast<TKey*>(i)->GetName()+std::string(".pdf"));
            std::cout << str << '\n';
            str.erase(std::remove(str.begin(), str.end(), '/'), str.end());
            TH1F* h = static_cast<TH1F*>(static_cast<TKey*>(i)->ReadObj());
            
            h->SetXTitle(GetXTitle(str).c_str());
            h->Draw();
            c1.SaveAs((""+str).c_str());
        }
        else
        {
            recursiveTH1Fsave(static_cast<TDirectoryFile*>(((TKey*)(i))->ReadObj())->GetListOfKeys());
        }
    }
}

void recursivesave()
{
    TFile f("mc16_13TeV.600750.PhPy8EG_AZNLO_ggH125_mA_p0_Cyy0p01_Czh1p0.NTUPLE.e8324_e7400_s3126_r10724_r10726_v2_out.root");
//    TFile f("Ntuple_data_test_out.root");
//    TFile f("Ntuple_MC_Za_mA5p0_v4_out.root");
//    TFile f("example_mc_haa_out_testcpp.root");
    recursiveTH1Fsave(f.GetListOfKeys());
    system("convert *pdf -quality 100 file.pdf"); //Only if imagemagick is installed.
    system(R"--(ls *pdf | grep -xv "file.pdf" | parallel rm)--"); //Only if GNU parallel is installed
}

