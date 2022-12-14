#include "candidate.h"
#include "objects.h"

#include <utility>

#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

template <typename T>
CandidateSet<T>::CandidateSet(std::pair<T, T> particles)
: particle_a{particles.first},
  particle_b{particles.second}
{
    four_momentum = particle_a.Vector() + particle_b.Vector();
}

template <typename T>
double CandidateSet<T>::z_mass()
{
    return four_momentum.M();
}

template <typename T>
double CandidateSet<T>::delta_phi()
{
//    return particle_a.Vector().DeltaPhi(particle_b.Vector());
    return VectorUtil::DeltaPhi(particle_a.Vector(),particle_b.Vector());
}

template <typename T>
double CandidateSet<T>::delta_eta()
{
    return particle_a.Vector().Eta() - particle_b.Vector().Eta();
}

template <typename T>
double CandidateSet<T>::acoplanarity()
{
//    return 1 - acos(cos(particle_a.phi() - particle_b.phi())) / TMath::Pi();
    return 1 - acos(cos(particle_a.phi() - particle_b.phi())) / Pi();
}

template class CandidateSet<TruthParticle>;
template class CandidateSet<Photon>;
template class CandidateSet<Electron>;
