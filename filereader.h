#ifndef FILEREADERH
#define FILEREADERH
#include <memory>
#include "event.h"



using FT = bool(TChain&);

class FileReader
{
private:
    friend class Range;
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

class Range {
    FileReader f;
public:

     Range(std::vector<string>&& files,
           const char* tree_name =
           "../user.kschmied.28655874._000025.LGNTuple.root?#physics",
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


//class Range {
//    FileReader f;
//public:
//
//     Range(std::vector<string>&& files,
//           const char* tree_name =
//           "../user.kschmied.28655874._000025.LGNTuple.root?#physics",
////           "/Users/edwardfinkelstein/mnt/droplet/user.kschmied.28655874._000025.LGNTuple.root?#physics",
//           Long64_t num_events = -1, int skip_first_events = 0)
////    Range(std::vector<string>&& files, const char* tree_name = "physics", Long64_t num_events = -1, int skip_first_events = 0)
//    :
//    f(files,tree_name, num_events, skip_first_events)
//    {}
//
//     struct Iterator
//     {
//
//         Iterator(int i) : data{i} {}
//         Iterator (int i, FileReader& f) : data{i}, F{f}
//         {}
//         Iterator& operator++()
//         {
//             data++;
////             printf("%d\n", data++);
//          Long64_t entryNumberWithinCurrentTree = F.__chain.LoadTree(F.__current_index); // The result can be used with TBranch::GetEntry
//          if (entryNumberWithinCurrentTree < 0)
//          {
//              // something went wrong.
//              puts("something went wrong.");
//              F.__current_index++;
//              return *this;
//          }
//             F.__chain.GetEntry(F.__current_index);
//             F.__event_info_chain.GetEntry(F.__current_index);
//             F.__current_index++;
//             if (F.__passes_event_filters())
//             {
//                 //Very slow!!!
//                 Event Temp(&(F.__chain), &(F.__event_info_chain));
//                 F.__current_event = Temp;
//             }
//             return *this;
//         }
//
//         Event operator*() {return F.__current_event;}
//         friend bool operator!=(const Iterator& a, const Iterator& b){
////             printf("%d %d\n",a.data,b.data);
//             return a.data != b.data;}
//         private:
//         int data;
//         FileReader F;
//     };
//
//
//    Iterator begin() {return Iterator(f.__skip_first_events, f);}
//    Iterator end() {return Iterator(20);}
////    Iterator end() {return Iterator(f.__num_events);}
//};





#endif
