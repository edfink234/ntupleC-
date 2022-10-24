#ifndef FILEREADERH
#define FILEREADERH
#include <memory>
#include "event.h"



using FT = bool(TChain&);

class FileReader
{
private:
    friend class FileReaderRange;
    std::vector<string> __files;
    int __skip_first_events;
    int __current_index;
//    Event* __current_event;
    Event __current_event;
//    shared_ptr<Event> __current_event;
    std::vector<FT*> __event_filters;
    TChain __chain;
    TChain __event_info_chain;
    Long64_t __num_events;
public:
    FileReader();
    FileReader(std::vector<string>& files, const char* tree_name = "physics", Long64_t num_events = -1, int skip_first_events = 0);
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

class FileReaderRange {
    FileReader f;
public:

     FileReaderRange(std::vector<string>&& files,
           const char* tree_name =
//           "../user.kschmied.28655874._000025.LGNTuple.root?#physics",
           "physics",
           Long64_t num_events = -1, int skip_first_events = 0);

    
     struct Iterator
     {
         
         Iterator(int i);
         Iterator (int i, FileReader& f);
         Iterator& operator++();

         Event operator*();
         
         private:
         friend bool operator!=(const Iterator& a, const Iterator& b);
         int data;
         FileReader F;
     };


    Iterator begin();
    Iterator end();
//    Iterator end() {return Iterator(f.__num_events);}
};








#endif
