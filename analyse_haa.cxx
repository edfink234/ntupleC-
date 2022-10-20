#include <unordered_map>
//#include <utility>
#include "TFile.h"
#include "plotting.h"
#include "objects.h"
#include "event.h"
#include "filereader.h"
#include "TInterpreter.h"



const vector<string> CUTS = {"truth","reco","reco_2y"};
unordered_map<string,Plot> plots;
unordered_map<string,PlotGroup> plot_groups;


void run_analysis(string& input_filename, string systematic = "nominal", bool mc = false)
{
//https://en.cppreference.com/w/cpp/utility/variant
    //union, for now
//https://www.sololearn.com/compiler-playground/cop9eIyns3c3
    
    
    plots.emplace(string("cutflow"), Plot(systematic+string("/cutflow"), string("cutflow"), 10, 0, 10));
    
    
//    cout << systematic << mc << input_filename;
    for (const string& cutname: CUTS)
    {
        plots.emplace(cutname+string("/photons/eta"),
            Plot(systematic+string(1,'/')+cutname
            +string("/photons/eta"),
            string("eta"), 40, -6, 6));
        
        plots.emplace(cutname+string("/leptons/eta"),
            Plot(systematic+string(1,'/')+cutname
            +string("/leptons/eta"),
            string("eta"), 40, -6, 6));
        
        plot_groups.emplace(
        cutname+string("/photons/pt"),
        PlotGroup({
            Plot(systematic+string(1,'/')+cutname
            +string("/photons/pt"),
            string("pT / GeV"), 20, 0, 200),
            Plot(systematic+string(1,'/')+cutname
            +string("/photons/pt_tight"),
            string("pT / GeV"), 100, 0, 100)
        }));

        plot_groups.emplace(
        cutname+string("/leptons/pt"),
        PlotGroup({
            Plot(systematic+string(1,'/')+cutname
            +string("/leptons/pt"),
            string("pT / GeV"), 20, 0, 200),
            Plot(systematic+string(1,'/')+cutname
            +string("/leptons/pt_tight"),
            string("pT / GeV"), 100, 0, 100)
        }));
        
//        PlotGroup({
//                    Plot(systematic+string(1,'/')+cutname
//                    +string("/leptons/pt"),
//                    string("pT / GeV"), 20, 0, 200),
//                    Plot(systematic+string(1,'/')+cutname
//                    +string("/leptons/pt_tight"),
//                    string("pT / GeV"), 100, 0, 100)
//        });
        

    }
    
    Range x({input_filename});
    Event::systematic = systematic;
    Event::load_tracks = true;
    int num_passed_events = 0;
    
//    cout << mc;
    using std::endl;
    
    for (auto &&i: x)
    {
//        cout << 0 << endl;
        
        cout << "event_number " << i.event_number  << '\n';
        
        int weight = 1;
//        cout << mc << endl;
        if (mc)
        {
            vector<TruthParticle>&& truth_higgs = i.find_truth_particles({},{},{35});
            if (!(truth_higgs.empty()))
            {
                cout << "found!";
                vector<TruthParticle>&& truth_axions = i.find_truth_particles({},{truth_higgs[0].barcode()}, {36});
            }
            int temp = 1;
//            vector<TruthParticle>&& truth_photons = i.find_truth_particles({},{},{22},&temp);
//            
//            cout << "Size = " << truth_photons.size() << endl;
//            for (auto i: truth_photons)
//            {
//                cout << string(i) << '\n';
//            }
        }
//        exit(1);
    }
    
    for (auto& plot: plots)
    {
        plot.second.save();
    }
    
    for (auto& plot_group: plot_groups)
    {
        plot_group.second.save();
    }
    
    
}

struct stat buffer;

void analyse_haa()
{

    cout << "Run over MC\n";
//    string input_filename("/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000025.LGNTuple.root");
//    string input_filename("/Users/edwardfinkelstein/mnt/droplet/user.kschmied.28655874._000025.LGNTuple.root");
    string input_filename("../user.kschmied.28655874._000025.LGNTuple.root");
    
    
    
    
    //https://stackoverflow.com/a/12774387/18255427
    if ((stat (input_filename.c_str(), &buffer) != 0))
    {
        cout << "The MC file does not exist, please correct the path in analyse_haa.py\n";
        exit(1);
    }
    
    const char* output_filename = "example_mc_haa_out_test1cpp.root";
    TFile* output_file = TFile::Open(output_filename, "RECREATE");
    if (!output_file) {
       std::cout << "Error opening file" << std::endl;
       exit(-1);
    }
    

//
//    TFile output_file(output_file, "RECREATE");
    
    vector<string> systematics = {"nominal"};
    
    for (string& systematic: systematics)
    {
        run_analysis(input_filename, systematic.c_str(), true);
    }
//    output_file->Close();
    
    
}



int main()
{
    analyse_haa();
}
