#ifndef RDFFILEREADERH
#define RDFFILEREADERH

#include "ROOT/RDataFrame.hxx"
#include "TChain.h"

#include <vector>
#include <string>
#include <memory>

#include "RDFevent.h"

class RDFFileReader
{
private:
    struct __callEnableImplicitMT
    {
        __callEnableImplicitMT(unsigned int);
    };
    __callEnableImplicitMT setThreads;
    TChain __chain;
    TChain __event_info_chain;
    std::unique_ptr<ROOT::RDataFrame> df;
    bool __has_event_info_chain;
    std::unique_ptr<Event> __current_event;

public:
    RDFFileReader(const std::vector<std::string>& files, const char* tree_name = "physics", unsigned int numThreads = -1);
    ~RDFFileReader();
    RDFFileReader(const RDFFileReader&) = delete;
};


#endif
