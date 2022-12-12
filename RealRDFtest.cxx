#include <iostream>
#include <vector>
#include <chrono>
#include <memory>

#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "MakeRDF.h"

using Clock = std::chrono::high_resolution_clock;

//https://root-forum.cern.ch/t/rdataframe-count-and-report-re-looping-over-whole-dataframe/46592
void RealRDFtest()
{
    auto start_time = Clock::now();

    std::vector<std::string> input_filenames = {"../user.kschmied.28655874._000025.LGNTuple.root"};
    
    SchottDataFrame df(MakeRDF(input_filenames));
    
    df.Describe().Print();

    
    
    auto end_time = Clock::now();
    std::cout << "Time difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
    
}


int main()
{
    RealRDFtest();
}
