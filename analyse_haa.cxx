#include <unordered_map>
//#include <utility>
#include "TFile.h"
#include "plotting.h"



const vector<string> CUTS = {"truth","reco","reco_2y"};
unordered_map<string,Plot> plots;
unordered_map<string,PlotGroup> plot_groups;


void run_analysis(string& input_filename, string systematic = "nominal", bool mc = false)
{
//https://en.cppreference.com/w/cpp/utility/variant
    //union, for now
//https://www.sololearn.com/compiler-playground/cop9eIyns3c3
    
    plots.emplace(string("cutflow"), Plot(systematic+string("/cutflow"), string("cutflow"), 10, 0, 10));
    
    
    
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
    
    for (auto& plot: plots)
    {
        plot.second.save();
    }
    
    for (auto& plot_group: plot_groups)
    {
        plot_group.second.save();
    }
    
    
}

//struct stat buffer;

void analyse_haa()
{
    cout << "Run over MC\n";
//    string input_filename("/home/common/Haa/ntuples/MC/background_v14/user.kschmied.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_v14_LGNTuple.root/user.kschmied.28655874._000025.LGNTuple.root");
    string input_filename("/Users/edwardfinkelstein/mnt/droplet/user.kschmied.28655874._000025.LGNTuple.root");
    
    
    
    
    //https://stackoverflow.com/a/12774387/18255427
//    if ((stat (input_filename.c_str(), &buffer) != 0))
//    {
//        cout << "The MC file does not exist, please correct the path in analyse_haa.py\n";
//        exit(1);
//    }
    
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
    
    
}



int main()
{
    analyse_haa();
}
