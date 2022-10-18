#include "filereader.h"

FileReader::FileReader(vector<string>& files, const char* tree_name, Long64_t num_events, int skip_first_events) :

__files{move(files)},
__skip_first_events{skip_first_events},
__current_index{skip_first_events},
__current_event{nullptr}

{
    __chain.Add(tree_name);
    __event_info_chain.Add("full_event_info");
    for (auto& f: files)
    {
        __chain.Add(f.c_str());
        __event_info_chain.Add(f.c_str());
    }
    if ((num_events<0) || (num_events>__chain.GetEntries()))
    {
        __num_events = __chain.GetEntries();
    }
    else
    {
        __num_events = num_events;
    }
}

FileReader::~FileReader()
{
    __chain.Reset();
    __event_info_chain.Reset();
}

FileReader::FileReader(const FileReader& other)
{
    __files = other.__files;
    __skip_first_events = other.__skip_first_events;
    __current_index = other.__current_index;
    __current_event = other.__current_event;
    __event_filters = other.__event_filters;
    
    __chain.Add(const_cast<TChain*>(&other.__chain)); //😬
    __event_info_chain.Add(const_cast<TChain*>(&other.__event_info_chain)); //😬
    __num_events = other.__num_events;
}

template <typename T>
void FileReader::add_event_filter(T& filter)
{
    __event_filters.push_back(filter);
}

bool FileReader::__passes_event_filters()
{
    if (__event_filters.empty())
    {
        return true;
    }
    for (auto& filter_func: __event_filters)
    {
        if (!((*filter_func)(__chain)))
        {
            return false;
        }
    }
    return true;
}


Event* FileReader::event()
{
    return __current_event;
}
int FileReader::current_index()
{
    return __current_index;
}
Long64_t FileReader::num_events()
{
    return __num_events;
}
