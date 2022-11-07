#include <iostream>
#include <algorithm>

#include "TH1F.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"

void recursiveTH1Fsave(TList* f)
{
    for (auto i: *f)
    {
        if (static_cast<TKey*>(i)->GetClassName()==string("TH1F"))
        {
            TCanvas c1;
            ((TKey*)(i))->ReadObj()->Draw();
            string str = (static_cast<TKey*>(i)->GetName()+string(".pdf"));
            str.erase(std::remove(str.begin(), str.end(), '/'), str.end());
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
//    TFile f("ntupleC++_v2/example_mc_haa_out_testcpp_final.root");
    TFile f("ntupleC++_v2/example_mc_haa_out_cppreg_final.root");
    recursiveTH1Fsave(f.GetListOfKeys());
}

