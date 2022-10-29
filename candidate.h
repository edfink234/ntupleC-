#ifndef CANDIDATEH
#define CANDIDATEH

#include "objects.h"

#include <utility>

class CandidateSet {
    TLorentzVector four_momentum;
private:
    TruthParticle particle_a;
    TruthParticle particle_b;

public:
    //Only 2 particles can be grouped together for now.
    CandidateSet(std::pair<TruthParticle, TruthParticle>);
    double z_mass();
    double delta_phi();
    double delta_eta();
    double acoplanarity();
};

#endif
