#include "event.h"

/*
 Default values for static attributes of Event, can be changed by the
 user.
 */

bool Event::cache_truth = true;
bool Event::load_reco = true;
bool Event::load_photons = true;
bool Event::load_electrons = true;
bool Event::load_muons = true;
bool Event::load_clusters = true;
bool Event::load_tracks = true;
bool Event::load_triggers = true;
std::string Event::systematic = "";
