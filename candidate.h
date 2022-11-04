#ifndef CANDIDATEH
#define CANDIDATEH

#include "objects.h"

#include <utility>

class CandidateSet {
private:
    enum class ParticleType{TruthParticle, Photon};
    ParticleType particle_type;
public:
    //Only 2 particles can be grouped together for now.
    CandidateSet(std::pair<TruthParticle, TruthParticle>);
    CandidateSet(std::pair<Photon, Photon>);
    TLorentzVector four_momentum;
    TruthParticle particle_a;
    TruthParticle particle_b;
    Photon particle_a_photon;
    Photon particle_b_photon;
    double z_mass();
    double delta_phi();
    double delta_eta();
    double acoplanarity();
    char CandidateSetType();
};

#endif
