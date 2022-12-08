#include <iostream>
#include <vector>
#include <chrono>

#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "RDFFileReader.h"

using Clock = std::chrono::high_resolution_clock;

//https://root-forum.cern.ch/t/rdataframe-count-and-report-re-looping-over-whole-dataframe/46592
void RdataFrametest()
{
    auto start_time = Clock::now();
//    std::vector<const char*> input_filenames = {"../user.kschmied.28655874._000025.LGNTuple.root","../user.kschmied.28655874._000024.LGNTuple.root"};
    std::vector<std::string> input_filenames = {"../user.kschmied.28655874._000025.LGNTuple.root"};
//    {"Ntuple_data_test.root"};
    
    RDFFileReader x(input_filenames);
//    RDFFileReader x(input_filenames);
    auto end_time = Clock::now();
    std::cout << "Time difference:"
       << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " nanoseconds" << std::endl;
    
}


int main()
{
    RdataFrametest();
}
