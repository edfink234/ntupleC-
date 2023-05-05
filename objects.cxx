#include "objects.h"
#include "event.h"

//#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#include <climits>

using namespace ROOT::Math;
              
//                *****************
//                * PhysicsObject *
//                *****************

/*
 Create a PhysicsObject that corresponds to an event
 in a TTree. The index argument refers to the element
 the static vector attributes to index/return, if
 one creates an object that derives from PhysicsObject,
 i.e., TruthParticle, Electron, Photon, Muon, Cluster, Track
 
 e.g.
 
 PhysicsObject p{1};
 */

PhysicsObject::PhysicsObject(int index)
:  _index{index}
{

}

PhysicsObject::PhysicsObject() = default;

PhysicsObject::~PhysicsObject() = default;

//TLorentzVector PhysicsObject::Vector()
//{
//    TLorentzVector vec;
//    vec.SetPtEtaPhiE(pt(),eta(),phi(),e());
//    return vec;
//}

/*
 As recommend by the docs, ROOT::Math::LorentzVector is the ``superior alternative''
 to TLorentzVector (commented out above)
 */
PtEtaPhiEVector PhysicsObject::Vector()
{
    PtEtaPhiEVector vec;
    vec.SetCoordinates(pt(),eta(),phi(),e());
    return vec;
}

double PhysicsObject::acoplanarity(TruthParticle& particle_b)
{
//    return 1 - acos(cos(this->phi() - particle_b.phi())) / TMath::Pi();
    return 1 - acos(cos(this->phi() - particle_b.phi())) / Pi();
}

double PhysicsObject::delta_r(TruthParticle& particle_b)
{
//    return this->Vector().DeltaR(particle_b.Vector());
    return VectorUtil::DeltaR(this->Vector(),particle_b.Vector());
}

PhysicsObject::operator std::string() //const
{
    return
    (
    "Particle(pT = "+std::to_string(pt())
    +", eta = " + std::to_string(eta())
    +", phi = " + std::to_string(phi())
    +") ");
}

const std::string PhysicsObject::PREFIX = "";

//                *****************
//                * TruthParticle *
//                *****************

/*
 Create a TruthParticle that corresponds to an event
 in a TTree. The index argument refers to the element
 the static vector attributes index/return.
 If the static variable Event::cache_truth
 is set to true, then its pdg_id is set to the _index'th
 element of the event's mc_pdg_id vector, else, it's
 set to INT_MIN, which no particle corresponds to.
 
 e.g.
 
 TruthParticle tp{1};
 */

TruthParticle::TruthParticle(int index)
:   PhysicsObject::PhysicsObject(index),
    Cluster_eta{0.0}
{
    //https://root-forum.cern.ch/t/looping-over-a-ttree/19899/5
    //https://root.cern/root/html534/tutorials/tree/hvector.C.html
    if (Event::cache_truth)
    {
        pdg_id = (*TruthParticle::mc_pdg_id)[_index];
    }
    else
    {
        pdg_id = INT_MIN; //so undefined :)
    }
    
}

TruthParticle::TruthParticle() = default;

TruthParticle::TruthParticle(const TruthParticle& particle_temp)
{
    pdg_id = particle_temp.pdg_id;
    _index = particle_temp._index;
    Cluster_eta = particle_temp.Cluster_eta;
}

TruthParticle& TruthParticle::operator=(const TruthParticle& particle_temp)
{
    pdg_id = particle_temp.pdg_id;
    _index = particle_temp._index;
    Cluster_eta = particle_temp.Cluster_eta;
    return *this;
}

TruthParticle::~TruthParticle() = default;

//TLorentzVector TruthParticle::Vector()
//{
//    TLorentzVector vec;
//    vec.SetPtEtaPhiE(pt(),eta(),phi(),e());
//    return vec;
//}
PtEtaPhiEVector TruthParticle::Vector()
{
    PtEtaPhiEVector vec;
    vec.SetCoordinates(pt(),eta(),phi(),e());
    return vec;
}

int TruthParticle::barcode()
{
    return (*TruthParticle::mc_barcode)[_index];
}

int TruthParticle::parent_barcode()
{
    return (*TruthParticle::mc_parent_barcode)[_index];
}

int TruthParticle::status_code()
{
    return (*TruthParticle::mc_status)[_index];
}

double TruthParticle::pt()
{
    return (*TruthParticle::mc_pt)[_index];
}

double TruthParticle::charge()
{
    if (Event::cache_truth)
    {
        return (*TruthParticle::mc_charge)[_index];
    }
    return 0; //photon
}

double TruthParticle::eta()
{
    return (*TruthParticle::mc_eta)[_index];
}

double TruthParticle::phi()
{
    return (*TruthParticle::mc_phi)[_index];
}

bool TruthParticle::same_flavour(const TruthParticle& other)
{
//    printf("%d %d\n",this->pdg_id,other.pdg_id);
    return (abs(this->pdg_id)==abs(other.pdg_id));
}

double TruthParticle::e()
{
    return (*TruthParticle::mc_e)[_index];
}

double TruthParticle::m()
{
    return (*TruthParticle::mc_mass)[_index];
}

TruthParticle::operator std::string()
{
    return
    (
    "TruthParticle(pdg_id = "+std::to_string(pdg_id)
    +", pT = " + std::to_string(pt())
    +", charge = " + std::to_string(charge())
    +", eta = " + std::to_string(eta())
    +", phi = " + std::to_string(phi())
    +") ");
}

int TruthParticle::id_(){return true;}
int TruthParticle::id_loose(){return true;}
int TruthParticle::id_tight(){return true;}

const std::string TruthParticle::PREFIX = "mc";
std::vector<int>* TruthParticle::mc_pdg_id = nullptr;
std::vector<int>* TruthParticle::mc_barcode = nullptr;
std::vector<int>* TruthParticle::mc_parent_barcode = nullptr;
std::vector<int>* TruthParticle::mc_status = nullptr;
std::vector<double>* TruthParticle::mc_pt = nullptr;
std::vector<double>* TruthParticle::mc_charge = nullptr;
std::vector<double>* TruthParticle::mc_eta = nullptr;
std::vector<double>* TruthParticle::mc_phi = nullptr;
std::vector<double>* TruthParticle::mc_e = nullptr;
std::vector<double>* TruthParticle::mc_mass = nullptr;

/*
 Sets the branch status of branches in the TChain `chain', as
 well as set the static attributes of Truthparticle
 to the addresses of those branches.
 */
void TruthParticle::SetTruthParticle(TChain* chain)
{
    if (chain->GetBranch("mc_pdg_id")) chain->SetBranchStatus("mc_pdg_id",1);
    if (chain->GetBranch("mc_barcode")) chain->SetBranchStatus("mc_barcode",1);
    if (chain->GetBranch("mc_parent_barcode")) chain->SetBranchStatus("mc_parent_barcode",1);
    if (chain->GetBranch("mc_status")) chain->SetBranchStatus("mc_status",1);
    if (chain->GetBranch("mc_pt")) chain->SetBranchStatus("mc_pt",1);
    if (chain->GetBranch("mc_charge")) chain->SetBranchStatus("mc_charge",1);
    if (chain->GetBranch("mc_eta")) chain->SetBranchStatus("mc_eta",1);
    if (chain->GetBranch("mc_phi")) chain->SetBranchStatus("mc_phi",1);
    if (chain->GetBranch("mc_e")) chain->SetBranchStatus("mc_e",1);
    if (chain->GetBranch("mc_mass")) chain->SetBranchStatus("mc_mass",1);
    
    if (chain->GetBranch("mc_pdg_id")) chain->SetBranchAddress("mc_pdg_id",&TruthParticle::mc_pdg_id);
    if (chain->GetBranch("mc_barcode")) chain->SetBranchAddress("mc_barcode",&TruthParticle::mc_barcode);
    if (chain->GetBranch("mc_parent_barcode")) chain->SetBranchAddress("mc_parent_barcode",&TruthParticle::mc_parent_barcode);
    if (chain->GetBranch("mc_status")) chain->SetBranchAddress("mc_status",&TruthParticle::mc_status);
    if (chain->GetBranch("mc_pt")) chain->SetBranchAddress("mc_pt",&TruthParticle::mc_pt);
    if (chain->GetBranch("mc_charge")) chain->SetBranchAddress("mc_charge",&TruthParticle::mc_charge);
    if (chain->GetBranch("mc_eta"))  chain->SetBranchAddress("mc_eta",&TruthParticle::mc_eta);
    if (chain->GetBranch("mc_phi")) chain->SetBranchAddress("mc_phi",&TruthParticle::mc_phi);
    if (chain->GetBranch("mc_e")) chain->SetBranchAddress("mc_e",&TruthParticle::mc_e);
    if (chain->GetBranch("mc_mass")) chain->SetBranchAddress("mc_mass",&TruthParticle::mc_mass);
}

//                *****************
//                *    Electron   *
//                *****************

/*
 Create an Electron that corresponds to an event
 in a TTree. The index argument refers to the element
 the static vector attributes index/return.
 The pt and the energy are in MeV.
 
 e.g.
 
 Electron e{1, 1.1, 1.1};
 */

Electron::Electron(int index, double pt, double energy)
: TruthParticle::TruthParticle(index),
  _systematic_pt{pt},
  _systematic_energy{energy}
{
    
}

Electron::Electron() = default;

Electron::Electron(const Electron& other)
{
    pdg_id = other.pdg_id;
    _index = other._index;
    Cluster_eta = other.Cluster_eta;
    _systematic_pt = other._systematic_pt;
    _systematic_energy = other._systematic_energy;
}

Electron& Electron::operator=(const Electron& other)
{
    pdg_id = other.pdg_id;
    _index = other._index;
    Cluster_eta = other.Cluster_eta;
    _systematic_pt = other._systematic_pt;
    _systematic_energy = other._systematic_energy;
    return *this;
}

int Electron::Pdg_id()
{
    return charge()*PDG_ID;
}

double Electron::pt()
{
    if (_systematic_pt) //if it's not 0
    {
        return _systematic_pt;
    }
    return (*Electron::electron_pt)[_index];
}

double Electron::charge()
{
    return (*Electron::electron_charge)[_index];
}

double Electron::e()
{
    if (_systematic_energy)
    {
        return _systematic_energy;
    }
    return (*Electron::electron_e)[_index];
}

double Electron::eta()
{
    return (*Electron::electron_eta)[_index];
}

double Electron::phi()
{
    return (*Electron::electron_phi)[_index];
}

int Electron::id_()
{
    return (*Electron::electron_id)[_index];
}

int Electron::id_medium()
{
    return (*Electron::electron_id_medium)[_index];
}

double Electron::isolation()
{
    return (*Electron::electron_isolation)[_index];
}

double Electron::d0()
{
    return (*Electron::electron_d0)[_index];
}

double Electron::z0()
{
    return (*Electron::electron_z0)[_index];
}

Electron::operator std::string()
{
    return
    (
    "Electron(pdg_id = "+std::to_string(Pdg_id())
    +", pT = " + std::to_string(pt())
    +", charge = " + std::to_string(charge())
    +", eta = " + std::to_string(eta())
    +", phi = " + std::to_string(phi())
    +") ");
}

const int Electron::PDG_ID = 11;
const std::string Electron::PREFIX = "electron";
std::vector<double>* Electron::electron_charge = nullptr;
std::vector<double>* Electron::electron_pt = nullptr;
std::vector<double>* Electron::electron_e = nullptr;
std::vector<double>* Electron::electron_eta = nullptr;
std::vector<double>* Electron::electron_phi = nullptr;
std::vector<int>* Electron::electron_id = nullptr;
std::vector<int>* Electron::electron_id_medium = nullptr;
std::vector<double>* Electron::electron_isolation = nullptr;
std::vector<double>* Electron::electron_d0 = nullptr;
std::vector<double>* Electron::electron_z0 = nullptr;

/*
 Sets the branch status of branches in the TChain `chain', as
 well as set the static attributes of Electron
 to the addresses of those branches.
 */

void Electron::SetElectron(TChain* chain)
{
    if (chain->GetBranch("electron_charge"))  chain->SetBranchStatus("electron_charge",1);
    if (chain->GetBranch("electron_pt"))  chain->SetBranchStatus("electron_pt",1);
    if (chain->GetBranch("electron_e"))  chain->SetBranchStatus("electron_e",1);
    if (chain->GetBranch("electron_eta"))  chain->SetBranchStatus("electron_eta",1);
    if (chain->GetBranch("electron_phi"))  chain->SetBranchStatus("electron_phi",1);
//    if (chain->GetBranch("electron_id"))  chain->SetBranchStatus("electron_id",1);
    if (chain->GetBranch("electron_isolation"))  chain->SetBranchStatus("electron_isolation",1);
    if (chain->GetBranch("electron_d0"))  chain->SetBranchStatus("electron_d0",1);
    if (chain->GetBranch("electron_z0"))  chain->SetBranchStatus("electron_z0",1);
//    if (chain->GetBranch("electron_id_medium"))  chain->SetBranchStatus("electron_id_medium",1);

    if (chain->GetBranch("electron_charge"))  chain->SetBranchAddress("electron_charge",&Electron::electron_charge);
    if (chain->GetBranch("electron_pt"))  chain->SetBranchAddress("electron_pt",&Electron::electron_pt);
    if (chain->GetBranch("electron_e"))  chain->SetBranchAddress("electron_e",&Electron::electron_e);
    if (chain->GetBranch("electron_eta"))  chain->SetBranchAddress("electron_eta",&Electron::electron_eta);
    if (chain->GetBranch("electron_phi"))  chain->SetBranchAddress("electron_phi",&Electron::electron_phi);
//    if (chain->GetBranch("electron_id"))  chain->SetBranchAddress("electron_id",&Electron::electron_id);
    if (chain->GetBranch("electron_isolation"))  chain->SetBranchAddress("electron_isolation",&Electron::electron_isolation);
    if (chain->GetBranch("electron_d0"))  chain->SetBranchAddress("electron_d0",&Electron::electron_d0);
    if (chain->GetBranch("electron_z0"))  chain->SetBranchAddress("electron_z0",&Electron::electron_z0);
//    if (chain->GetBranch("electron_id_medium"))  chain->SetBranchAddress("electron_id_medium",&Electron::electron_id_medium);
}

//                *****************
//                *      Muon     *
//                *****************

/*
 Create an Muon that corresponds to an event
 in a TTree. The index argument refers to the element
 the static vector attributes index/return.
 The pt and the energy are in MeV.
 
 e.g.
 
 Muon m{1, 1.1, 1.1};
 */

Muon::Muon(int index, double pt, double energy)
: TruthParticle::TruthParticle(index),
  _systematic_pt{pt},
  _systematic_energy{energy}
{
    
}

Muon::Muon() = default;

Muon::Muon(const Muon& other)
{
    pdg_id = other.pdg_id;
    _index = other._index;
    Cluster_eta = other.Cluster_eta;
    _systematic_pt = other._systematic_pt;
    _systematic_energy = other._systematic_energy;
}

Muon& Muon::operator=(const Muon& other)
{
    pdg_id = other.pdg_id;
    _index = other._index;
    Cluster_eta = other.Cluster_eta;
    _systematic_pt = other._systematic_pt;
    _systematic_energy = other._systematic_energy;
    return *this;
}

int Muon::Pdg_id()
{
    return charge()*PDG_ID;
}

double Muon::pt()
{
    if (_systematic_pt) //if it's not 0
    {
        return _systematic_pt;
    }
    return (*Muon::muon_pt)[_index];
}

double Muon::charge()
{
    return (*Muon::muon_charge)[_index];
}

double Muon::e()
{
    if (_systematic_energy)
    {
        return _systematic_energy;
    }
    return (*Muon::muon_e)[_index];
}

double Muon::eta()
{
    return (*Muon::muon_eta)[_index];
}

double Muon::phi()
{
    return (*Muon::muon_phi)[_index];
}

Muon::operator std::string()
{
    return
    (
    "Muon(pdg_id = "+std::to_string(Pdg_id())
    +", pT = " + std::to_string(pt())
    +", charge = " + std::to_string(charge())
    +", eta = " + std::to_string(eta())
    +", phi = " + std::to_string(phi())
    +") ");
}

const int Muon::PDG_ID = 13;
const std::string Muon::PREFIX = "muon";
std::vector<double>* Muon::muon_charge = nullptr;
std::vector<double>* Muon::muon_pt = nullptr;
std::vector<double>* Muon::muon_e = nullptr;
std::vector<double>* Muon::muon_eta = nullptr;
std::vector<double>* Muon::muon_phi = nullptr;

/*
 Sets the branch status of branches in the TChain `chain', as
 well as set the static attributes of Muon
 to the addresses of those branches.
 */

void Muon::SetMuon(TChain* chain)
{
    if (chain->GetBranch("muon_charge"))  chain->SetBranchStatus("muon_charge",1);
    if (chain->GetBranch("muon_pt"))  chain->SetBranchStatus("muon_pt",1);
    if (chain->GetBranch("muon_e"))  chain->SetBranchStatus("muon_e",1);
    if (chain->GetBranch("muon_eta"))  chain->SetBranchStatus("muon_eta",1);
    if (chain->GetBranch("muon_phi"))  chain->SetBranchStatus("muon_phi",1);

    if (chain->GetBranch("muon_charge"))  chain->SetBranchAddress("muon_charge",&Muon::muon_charge);
    if (chain->GetBranch("muon_pt"))  chain->SetBranchAddress("muon_pt",&Muon::muon_pt);
    if (chain->GetBranch("muon_e"))  chain->SetBranchAddress("muon_e",&Muon::muon_e);
    if (chain->GetBranch("muon_eta"))  chain->SetBranchAddress("muon_eta",&Muon::muon_eta);
    if (chain->GetBranch("muon_phi"))  chain->SetBranchAddress("muon_phi",&Muon::muon_phi);
}

//                *****************
//                *    Photon     *
//                *****************

/*
 Create a Photon that corresponds to an event
 in a TTree. The index argument refers to the element
 the static vector attributes index/return.
 The pt and the energy are in MeV.
 
 e.g.
 
 Photon p{1, 1.1, 1.1};
 */

Photon::Photon(int index, double pt, double energy)
: TruthParticle::TruthParticle(index),
  _systematic_pt{pt},
  _systematic_energy{energy}
{
    
}

Photon::Photon() = default;

Photon::Photon(const Photon& other)
{
    pdg_id = other.pdg_id;
    _index = other._index;
    Cluster_eta = other.Cluster_eta;
    _systematic_pt = other._systematic_pt;
    _systematic_energy = other._systematic_energy;
}

Photon& Photon::operator=(const Photon& other)
{
    pdg_id = other.pdg_id;
    _index = other._index;
    Cluster_eta = other.Cluster_eta;
    _systematic_pt = other._systematic_pt;
    _systematic_energy = other._systematic_energy;
    return *this;
}

double Photon::pt()
{
    if (_systematic_pt)
    {
        return _systematic_pt;
    }
    return (*Photon::photon_pt)[_index];
}

double Photon::e()
{
    if (_systematic_energy)
    {
        return _systematic_energy;
    }
    return (*Photon::photon_e)[_index];
}

double Photon::eta()
{
    return (*Photon::photon_eta)[_index];
}

double Photon::phi()
{
    return (*Photon::photon_phi)[_index];
}

double Photon::m()
{
    return 0;
}

double Photon::isolation()
{
    return (*Photon::photon_etcone40)[_index];
}

int Photon::id_()
{
    return (*Photon::photon_id)[_index];
}

int Photon::id_loose()
{
    return (*Photon::photon_id_loose)[_index];
}

int Photon::id_tight()
{
    return (*Photon::photon_id_tight)[_index];
}

double Photon::cluster_eta()
{
    return (*Photon::photon_cluster_eta_be_2)[_index];
}

int Photon::id_nn()
{
    return (*Photon::photon_id_nn)[_index];
}

Photon::operator std::string()
{
    return
    (
    "Photon(pdg_id = "+std::to_string(PDG_ID)
    +", pT = " + std::to_string(pt())
    +", charge = " + std::to_string(charge())
    +", eta = " + std::to_string(eta())
    +", phi = " + std::to_string(phi())
    +") ");
}

const int Photon::PDG_ID = 22;
const std::string Photon::PREFIX = "photon";

std::vector<double>* Photon::photon_pt = nullptr;
std::vector<double>* Photon::photon_e = nullptr;
std::vector<double>* Photon::photon_eta = nullptr;
std::vector<double>* Photon::photon_phi = nullptr;
std::vector<double>* Photon::photon_etcone40 = nullptr;
std::vector<int>* Photon::photon_id = nullptr;
std::vector<int>* Photon::photon_id_loose = nullptr;
std::vector<int>* Photon::photon_id_tight = nullptr;
std::vector<int>* Photon::photon_cluster_eta_be_2 = nullptr;
std::vector<int>* Photon::photon_id_nn = nullptr;

/*
 Sets the branch status of branches in the TChain `chain', as
 well as set the static attributes of Photon
 to the addresses of those branches.
 */

void Photon::SetPhoton(TChain* chain)
{
    if (chain->GetBranch("photon_pt")) chain->SetBranchStatus("photon_pt",1);
    if (chain->GetBranch("photon_e")) chain->SetBranchStatus("photon_e",1);
    if (chain->GetBranch("photon_eta")) chain->SetBranchStatus("photon_eta",1);
    if (chain->GetBranch("photon_phi")) chain->SetBranchStatus("photon_phi",1);
    if (chain->GetBranch("photon_etcone40")) chain->SetBranchStatus("photon_etcone40",1);
    if (chain->GetBranch("photon_id")) chain->SetBranchStatus("photon_id",1);
    if (chain->GetBranch("photon_id_loose")) chain->SetBranchStatus("photon_id_loose",1);
    if (chain->GetBranch("photon_id_tight")) chain->SetBranchStatus("photon_id_tight",1);
    if (chain->GetBranch("photon_cluster_eta_be_2")) chain->SetBranchStatus("photon_cluster_eta_be_2",1);
    
    if (chain->GetBranch("photon_pt")) chain->SetBranchAddress("photon_pt",&Photon::photon_pt);
    if (chain->GetBranch("photon_e")) chain->SetBranchAddress("photon_e",&Photon::photon_e);
    if (chain->GetBranch("photon_eta")) chain->SetBranchAddress("photon_eta",&Photon::photon_eta);
    if (chain->GetBranch("photon_phi")) chain->SetBranchAddress("photon_phi",&Photon::photon_phi);
    if (chain->GetBranch("photon_etcone40")) chain->SetBranchAddress("photon_etcone40",&Photon::photon_etcone40);
    if (chain->GetBranch("photon_id")) chain->SetBranchAddress("photon_id",&Photon::photon_id);
    if (chain->GetBranch("photon_id_loose")) chain->SetBranchAddress("photon_id_loose",&Photon::photon_id_loose);
    if (chain->GetBranch("photon_id_tight")) chain->SetBranchAddress("photon_id_tight",&Photon::photon_id_tight);
    if (chain->GetBranch("photon_cluster_eta_be_2")) chain->SetBranchAddress("photon_cluster_eta_be_2",&Photon::photon_cluster_eta_be_2);
//    chain->SetBranchAddress("photon_id_nn",&Photon::photon_id_nn);
}

//                *****************
//                *    Cluster    *
//                *****************

/*
 Create a Cluster that corresponds to an event
 in a TTree. The index argument refers to the element
 the static vector attributes index/return.
 
 e.g.
 
 Cluster c{1};
 */

Cluster::Cluster(int index)
:   PhysicsObject::PhysicsObject(index)
{
    
}

Cluster::Cluster() = default;

Cluster::Cluster(const Cluster& particle_temp)
{
    _index = particle_temp._index;
}

Cluster& Cluster::operator=(const Cluster& other)
{
    _index = other._index;
    return *this;
}

double Cluster::pt()
{
    return (*Cluster::cluster_pt)[_index];
}

double Cluster::eta()
{
    return (*Cluster::cluster_eta)[_index];
}

double Cluster::phi()
{
    return (*Cluster::cluster_phi)[_index];
}

double Cluster::e()
{
    return (*Cluster::cluster_e)[_index];
}

const std::string Cluster::PREFIX = "cluster";

std::vector<double>* Cluster::cluster_pt = nullptr;
std::vector<double>* Cluster::cluster_eta = nullptr;
std::vector<double>* Cluster::cluster_phi = nullptr;
std::vector<double>* Cluster::cluster_e = nullptr;

/*
 Sets the branch status of branches in the TChain `chain', as
 well as set the static attributes of Cluster
 to the addresses of those branches.
 */

void Cluster::SetCluster(TChain* chain)
{
    if (chain->GetBranch("cluster_pt"))  chain->SetBranchStatus("cluster_pt",1);
    if (chain->GetBranch("cluster_eta"))  chain->SetBranchStatus("cluster_eta",1);
    if (chain->GetBranch("cluster_phi"))  chain->SetBranchStatus("cluster_phi",1);
    if (chain->GetBranch("cluster_e"))  chain->SetBranchStatus("cluster_e",1);
    
    if (chain->GetBranch("cluster_pt"))  chain->SetBranchAddress("cluster_pt",&Cluster::cluster_pt);
    if (chain->GetBranch("cluster_eta"))  chain->SetBranchAddress("cluster_eta",&Cluster::cluster_eta);
    if (chain->GetBranch("cluster_phi"))  chain->SetBranchAddress("cluster_phi",&Cluster::cluster_phi);
    if (chain->GetBranch("cluster_e"))  chain->SetBranchAddress("cluster_e",&Cluster::cluster_e);
}

//                *****************
//                *    Track      *
//                *****************

/*
 Create a Track that corresponds to an event
 in a TTree. The index argument refers to the element
 the static vector attributes index/return.
 
 e.g.
 
 Track t{1};
 */

Track::Track(int index)
:   PhysicsObject::PhysicsObject(index)
{
    
}

Track::Track() = default;

Track::Track(const Track& particle_temp)
{
    _index = particle_temp._index;
}

Track& Track::operator=(const Track& other)
{
    _index = other._index;
    return *this;
}

double Track::pt()
{
    return (*Track::track_pt)[_index];
}

double Track::charge()
{
    return (*Track::track_charge)[_index];
}

double Track::eta()
{
    return (*Track::track_eta)[_index];
}

double Track::phi()
{
    return (*Track::track_phi)[_index];
}

double Track::e()
{
    return (*Track::track_e)[_index];
}

int Track::num_pixel_hits()
{
    return (*Track::track_num_pixel_hits)[_index];
}

int Track::num_sct_hits()
{
    return (*track_num_sct_hits)[_index];
}

Track::operator std::string()
{
    return
    (
     "Track(pT = "+std::to_string(pt())
    +", charge = " + std::to_string(charge())
    +", eta = " + std::to_string(eta())
    +", phi = " + std::to_string(phi())
    +") ");
}

const std::string Track::PREFIX = "track";

std::vector<double>* Track::track_pt = nullptr;
std::vector<double>* Track::track_charge = nullptr;
std::vector<double>* Track::track_eta = nullptr;
std::vector<double>* Track::track_phi = nullptr;
std::vector<double>* Track::track_e = nullptr;
std::vector<int>* Track::track_num_pixel_hits = nullptr;
std::vector<int>* Track::track_num_sct_hits = nullptr;

/*
 Sets the branch status of branches in the TChain `chain', as
 well as set the static attributes of Track
 to the addresses of those branches.
 */

void Track::SetTrack(TChain* chain)
{
    if (chain->GetBranch("track_pt"))  chain->SetBranchStatus("track_pt",1);
    if (chain->GetBranch("track_charge"))  chain->SetBranchStatus("track_charge",1);
    if (chain->GetBranch("track_eta"))  chain->SetBranchStatus("track_eta",1);
    if (chain->GetBranch("track_phi"))  chain->SetBranchStatus("track_phi",1);
//    if (chain->GetBranch("track_e"))  chain->SetBranchStatus("track_e",1);
    if (chain->GetBranch("track_num_pixel_hits"))  chain->SetBranchStatus("track_num_pixel_hits",1);
    if (chain->GetBranch("track_num_sct_hits"))  chain->SetBranchStatus("track_num_sct_hits",1);
    
    if (chain->GetBranch("track_pt"))  chain->SetBranchAddress("track_pt",&Track::track_pt);
    if (chain->GetBranch("track_charge"))  chain->SetBranchAddress("track_charge",&Track::track_charge);
    if (chain->GetBranch("track_eta"))  chain->SetBranchAddress("track_eta",&Track::track_eta);
    if (chain->GetBranch("track_phi"))  chain->SetBranchAddress("track_phi",&Track::track_phi);
//    if (chain->GetBranch("track_e"))  chain->SetBranchAddress("track_e",&Track::track_e);
    if (chain->GetBranch("track_num_pixel_hits"))  chain->SetBranchAddress("track_num_pixel_hits",&Track::track_num_pixel_hits);
    if (chain->GetBranch("track_num_sct_hits"))  chain->SetBranchAddress("track_num_sct_hits",&Track::track_num_sct_hits);
}
