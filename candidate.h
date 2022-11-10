#ifndef CANDIDATEH
#define CANDIDATEH

#include "objects.h"

#include <utility>

#include "Math/Vector4D.h"

using namespace ROOT::Math;

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
