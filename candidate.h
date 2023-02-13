#ifndef CANDIDATEH
#define CANDIDATEH

#include "objects.h"

#include <utility>

#include "Math/Vector4D.h"

using namespace ROOT::Math;

/*
 CandidateSet: Templated class to easily store info about
 2 particle systems, either two TruthParticle,
 Electron, or Photon objects atm.
 */

template <typename T>
class CandidateSet {
public:
    CandidateSet(std::pair<T, T>);
//    TLorentzVector four_momentum;
    PtEtaPhiEVector four_momentum;
    T particle_a;
    T particle_b;
    double z_mass();
    double delta_phi();
    double delta_eta();
    double acoplanarity();
};

#endif
