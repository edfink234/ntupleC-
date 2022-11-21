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
    return "";
}

void recursiveTH1Fsave(TList* f)
{
    for (auto i: *f)
    {
        if (static_cast<TKey*>(i)->GetClassName()==std::string("TH1F"))
        {
            TCanvas c1;
            std::string str = (static_cast<TKey*>(i)->GetName()+std::string(".png"));
            str.erase(std::remove(str.begin(), str.end(), '/'), str.end());
            TH1F* h = static_cast<TH1F*>(static_cast<TKey*>(i)->ReadObj());
            
            h->SetXTitle(GetXTitle(str).c_str());
            h->Draw();
            c1.SaveAs(str.c_str());
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
    recursiveTH1Fsave(f.GetListOfKeys());
}

