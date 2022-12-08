#include "RDFevent.h"

#include <string>

bool Event::cache_truth = true;
bool Event::load_reco = true;
bool Event::load_photons = true;
bool Event::load_electrons = true;
bool Event::load_muons = false;
bool Event::load_clusters = true;
bool Event::load_tracks = true;
bool Event::load_triggers = true;
std::string Event::systematic = "";

Event::Event()
{
    if (Event::cache_truth)
    {
        
    }
}
