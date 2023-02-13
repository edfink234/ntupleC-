#include "candidate.h"
#include "objects.h"

#include <utility>

#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

/*
 Constructor for a CandidateSet, i.e. a pair of exactly
 two PhysicsObjects of the same derived type. Only Electron,
 Photon, and TruthParticle are supported.
 
 e.g.
 
 CandidateSet<TruthParticle> truth_candidate(std::make_pair(TruthParticle(1),TruthParticle(2)));
 
 CandidateSet<Electron> electron_candidate(std::make_pair(Electron(1),Electron(2)));
 
 CandidateSet<Photon> photon_candidate(std::make_pair(Photon(1),Photon(2)));

 See more practical examples in analyse.cxx, analyse_haa.cxx, or myPreselectionHaa.cxx
 */

template <typename T>
CandidateSet<T>::CandidateSet(std::pair<T, T> particles)
: particle_a{particles.first},
  particle_b{particles.second}
{
    four_momentum = particle_a.Vector() + particle_b.Vector();
}

/*
 the invariant mass of the two-particle system
 */
template <typename T>
double CandidateSet<T>::z_mass()
{
    return four_momentum.M();
}

/*
 Δɸ of the two-particle system
 */
template <typename T>
double CandidateSet<T>::delta_phi()
{
//    return particle_a.Vector().DeltaPhi(particle_b.Vector());
    return VectorUtil::DeltaPhi(particle_a.Vector(),particle_b.Vector());
}

/*
 Δη of the two-particle system
 */
template <typename T>
double CandidateSet<T>::delta_eta()
{
    return particle_a.Vector().Eta() - particle_b.Vector().Eta();
}

/*
 acoplanarity of the two-particle system, defined as:
 
      acos(cos(ɸ_1 - ɸ_2))
 1 -  --------------------
               π
 */
template <typename T>
double CandidateSet<T>::acoplanarity()
{
//    return 1 - acos(cos(particle_a.phi() - particle_b.phi())) / TMath::Pi();
    return 1 - acos(cos(particle_a.phi() - particle_b.phi())) / Pi();
}

//forward declarations of allowed CandidateSet instances
template class CandidateSet<TruthParticle>;
template class CandidateSet<Photon>;
template class CandidateSet<Electron>;
