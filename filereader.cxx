#include "filereader.h"


FileReader::FileReader(std::vector<string>& files, const char* tree_name, Long64_t num_events, int skip_first_events) :

__files{move(files)},
__skip_first_events{skip_first_events},
__current_index{skip_first_events}
//__chain{tree_name},
//__event_info_chain{"full_event_info"}
//, __current_event{nullptr}

{
    
    __chain.Add(tree_name);
//    __event_info_chain.Add("/Users/edwardfinkelstein/mnt/droplet/user.kschmied.28655874._000025.LGNTuple.root?#full_event_info");
    __event_info_chain.Add("../user.kschmied.28655874._000025.LGNTuple.root?#full_event_info");
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

FileReader::FileReader() = default;

FileReader::~FileReader()
{
    __chain.Reset();
    __event_info_chain.Reset();
}

FileReader::FileReader(const FileReader& other)
{
    puts("called");
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


Event FileReader::event()
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

Range::Range(std::vector<string>&& files,
       const char* tree_name,
       Long64_t num_events, int skip_first_events) :
f(files,tree_name, num_events, skip_first_events)
{}

Range::Iterator::Iterator(int i) : data{i} {}

Range::Iterator::Iterator(int i, FileReader& f) : data{i}, F{f} {}

Range::Iterator& Range::Iterator::operator++()
{
    data++;
 Long64_t entryNumberWithinCurrentTree = F.__chain.LoadTree(F.__current_index);
 if (entryNumberWithinCurrentTree < 0)
 {
     // something went wrong.
     puts("something went wrong.");
     F.__current_index++;
     return *this;
 }
    F.__chain.GetEntry(F.__current_index);
    F.__event_info_chain.GetEntry(F.__current_index);
    F.__current_index++;
    if (F.__passes_event_filters())
    {

        Event Temp(&(F.__chain), &(F.__event_info_chain));
        F.__current_event = Temp;
    }
    return *this;
}

Event Range::Iterator::operator*() {return F.__current_event;}


bool operator!=(const Range::Iterator& a, const Range::Iterator& b)
{
    return a.data != b.data;
}

Range::Iterator Range::begin() {return Iterator(f.__skip_first_events, f);}
Range::Iterator Range::end() {return Iterator(20);}


