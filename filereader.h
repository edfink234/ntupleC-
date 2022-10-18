#ifndef FILEREADERH
#define FILEREADERH

#include "event.h"

using FT = bool(TChain&);

class FileReader
{
private:
    friend class Range;
    vector<string> __files;
    int __skip_first_events;
    int __current_index;
    Event* __current_event;
    vector<FT*> __event_filters;
    TChain __chain;
    TChain __event_info_chain;
    Long64_t __num_events;
public:
    FileReader(vector<string>&, const char* tree_name = "physics", Long64_t num_events = -1,
               int skip_first_events = 0);
    ~FileReader();
    FileReader(const FileReader&);
    template <typename T>
    void add_event_filter(T&);
    bool __passes_event_filters();
    Event* event();
    int current_index();
    Long64_t num_events();

};

class Range {
     int x, y;
//     const char* name;
//    vector<string>& files;
//    const char* tree_name;
//    Long64_t num_events;
//    int skip_first_events;
    FileReader f;
public:
//     Range() {}
     Range(vector<string>&& files, const char* tree_name = "/Users/edwardfinkelstein/mnt/droplet/user.kschmied.28655874._000025.LGNTuple.root?#physics", Long64_t num_events = -1,
           int skip_first_events = 0)
    : //files{files}, tree_name{tree_name}, num_events{num_events},
    //skip_first_events{skip_first_events}
    
    f(files,tree_name, num_events, skip_first_events)
    {}
//    : f{files, tree_name, num_events, skip_first_events}
//     {
//
//     }
    
     struct Iterator
     {
         Iterator (int i) : data{i} {}//{t.__setName("bye");}
         Iterator& operator++(){data++; return *this;}

         int operator*() const {return data;}
         friend bool operator!=(const Iterator& a, const Iterator& b){return a.data != b.data;}
         private:
         int data;
     };


     Iterator begin() {return Iterator(f.__skip_first_events);}
    Iterator end() {return Iterator(f.__num_events);}
};



#endif
