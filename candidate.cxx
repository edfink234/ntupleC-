#include "candidate.h"

#include <utility>

CandidateSet::CandidateSet(std::pair<TruthParticle, TruthParticle> particles)
: particle_a{particles.first},
  particle_b{particles.second},
  particle_type{ParticleType::TruthParticle}
{
    four_momentum = particle_a.Vector() + particle_b.Vector();
}

CandidateSet::CandidateSet(std::pair<Photon, Photon> particles)
: particle_a_photon{particles.first},
  particle_b_photon{particles.second},
  particle_type{ParticleType::Photon}
{
    four_momentum = particle_a_photon.Vector() + particle_b_photon.Vector();
}

double CandidateSet::z_mass()
{
    return four_momentum.M();
}

double CandidateSet::delta_phi()
{
    switch (particle_type)
    {
        case ParticleType::TruthParticle:
            return particle_a.Vector().DeltaPhi(particle_b.Vector());
            break;
        case ParticleType::Photon:
            return particle_a_photon.Vector().DeltaPhi(particle_b_photon.Vector());
            break;
    }
}

double CandidateSet::delta_eta()
{
    switch (particle_type)
    {
        case ParticleType::TruthParticle:
            return particle_a.Vector().Eta() - particle_b.Vector().Eta();
            break;
        case ParticleType::Photon:
            return particle_a_photon.Vector().Eta() - particle_b_photon.Vector().Eta();
            break;
    }
}

double CandidateSet::acoplanarity()
{
    switch (particle_type)
    {
        case ParticleType::TruthParticle:
            return 1 - acos(cos(particle_a.phi() - particle_b.phi())) / TMath::Pi();
            break;
        case ParticleType::Photon:
            return 1 - acos(cos(particle_a_photon.phi() - particle_b_photon.phi())) / TMath::Pi();
            break;
    }
}

char CandidateSet::CandidateSetType()
{
    switch (particle_type)
    {
        case ParticleType::TruthParticle:
            return 'T';
            break;
        case ParticleType::Photon:
            return 'P';
            break;
    }
}
