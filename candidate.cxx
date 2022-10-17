#include "candidate.h"
#include <stdexcept>


CandidateSet::CandidateSet(std::pair<TruthParticle, TruthParticle> particles)
: particle_a{particles.first},
  particle_b{particles.second}
{
    four_momentum = particle_a.Vector() + particle_b.Vector();
}

double CandidateSet::z_mass()
{
    return four_momentum.M();
}
double CandidateSet::delta_phi()
{
    return particle_a.Vector().DeltaPhi(particle_b.Vector());
}
double CandidateSet::delta_eta()
{
    return particle_a.Vector().Eta() - particle_b.Vector().Eta();
}
double CandidateSet::acoplanarity()
{
    return 1 - acos(cos(particle_a.phi() - particle_b.phi())) / TMath::Pi();
}
