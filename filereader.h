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
//    Event* __current_event;
    Event __current_event;
    vector<FT*> __event_filters;
    TChain __chain;
    TChain __event_info_chain;
    Long64_t __num_events;
public:
    FileReader();
    FileReader(vector<string>&, const char* tree_name = "physics", Long64_t num_events = -1,
               int skip_first_events = 0);
    ~FileReader();
    FileReader(const FileReader&);
    template <typename T>
    void add_event_filter(T&);
    bool __passes_event_filters();
//    Event* event();
    Event event();
    int current_index();
    Long64_t num_events();

};

class Range {
    FileReader f;
public:

     Range(vector<string>&& files, const char* tree_name = "/Users/edwardfinkelstein/mnt/droplet/user.kschmied.28655874._000025.LGNTuple.root?#physics", Long64_t num_events = -1,
           int skip_first_events = 0)
    :
    f(files,tree_name, num_events, skip_first_events)
    {}
    
     struct Iterator
     {
         
         Iterator(int i) : data{i} {}
         Iterator (int i, FileReader& f) : data{i}, F{f} {}
         Iterator& operator++()
         {
//             data++;
             printf("%d\n", data++);
          Long64_t entryNumberWithinCurrentTree = F.__chain.LoadTree(F.__current_index); // The result can be used with TBranch::GetEntry
          if (entryNumberWithinCurrentTree < 0)
          {
              // something went wrong.
              puts("something went wrong.");
              F.__current_index++;
              return *this;
          }
             F.__chain.GetEntry(F.__current_index);
//             F.__event_info_chain.GetEntry(F.__current_index);
             F.__current_index++;
//             if (F.__passes_event_filters())
//             {
//                 F.__current_event = Event(&(F.__chain), &(F.__event_info_chain));
//             }
             return *this;
         }

//         Event operator*() const //{return data;}
//         {
//             return (F.__current_event);
//         }
         int operator*() const {return data;}
         friend bool operator!=(const Iterator& a, const Iterator& b){return a.data != b.data;}
         private:
         int data;
         FileReader F;
     };


    Iterator begin() {return Iterator(f.__skip_first_events, f);}
    Iterator end() {return Iterator(f.__num_events);}
};



#endif
