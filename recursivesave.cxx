/*
 A utility file to help save arbitrarily nested TH1F objects into separate
 image files, specifically, anything that is allowed by TPad::SaveAs
 
 Simply run this in the ROOT prompt via
 
 root [0] .x recursivesave.cxx
 
 and it will save the
 */

#include <iostream>
#include <algorithm>

#include "TH1F.h"
#include "TKey.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"

//Function to infer x-axis title from histogram title
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
//function to recursively save TH1F objects from the given root file keys
void recursiveTH1Fsave(TList* f)
{
    for (auto i: *f)
    {
        if (static_cast<TKey*>(i)->GetClassName()==std::string("TH1F"))
        {
            TCanvas c1;
            std::string str = (static_cast<TKey*>(i)->GetName()+std::string(".png"));
            std::cout << str << '\n';
            str.erase(std::remove(str.begin(), str.end(), '/'), str.end());
            TH1F* h = static_cast<TH1F*>(static_cast<TKey*>(i)->ReadObj());
            
            h->SetXTitle(GetXTitle(str).c_str());
            h->Draw();
            c1.SaveAs((str).c_str());
        }
        else
        {
            recursiveTH1Fsave(static_cast<TDirectoryFile*>(((TKey*)(i))->ReadObj())->GetListOfKeys());
        }
    }
}

void recursivesave()
{
    TFile f("example_mc_haa_out_test1cppdummy.root"); //the file we want the TH1F's from
    recursiveTH1Fsave(f.GetListOfKeys()); //calling the function
//    system("convert *pdf -quality 100 file.pdf"); //Only if imagemagick is installed, converts all pngs into a single pdf
//    system(R"--(ls *pdf | grep -xv "file.pdf" | parallel rm)--"); //Only if GNU parallel is installed, removes all pdfs except file.pdf
}

