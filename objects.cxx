#include "objects.h"
#include "TMath.h"
#include <cassert>

using std::to_string;

//TODO: Maybe add some namespaces?

                
//                *****************
//                * PhysicsObject *
//                *****************
                 
PhysicsObject::PhysicsObject() = default;

PhysicsObject::PhysicsObject(TChain* entry, int index, const char* name, const char* title)
:  _index{index}
{
    _entry.Add(entry);
    _entry.SetName(name);
    _entry.SetTitle(title);
}

PhysicsObject::~PhysicsObject() = default;

TLorentzVector PhysicsObject::Vector()
{
    TLorentzVector vec;
    vec.SetPtEtaPhiE(pt(),eta(),phi(),e());
    return vec;
}

double PhysicsObject::acoplanarity(TruthParticle& particle_b)
{
    return 1 - acos(cos(this->phi() - particle_b.phi())) / TMath::Pi();
}

double PhysicsObject::delta_r(TruthParticle& particle_b)
{
    return this->Vector().DeltaR(particle_b.Vector());
}

PhysicsObject::operator string() //const
{
    using s = string; //alias
    
    return
    (
    s("Particle(pT=")+to_string(pt())
    +s(", eta=") + to_string(eta())
    +s(", phi=") + to_string(phi())
    );
}

const string PhysicsObject::PREFIX = "";

//                *****************
//                * TruthParticle *
//                *****************

TruthParticle::TruthParticle() = default;

TruthParticle::TruthParticle(TChain* entry, int index, const char* name, const char* title)
:   PhysicsObject::PhysicsObject(entry, index, name, title),
    Cluster_eta{0.0}
{
    //https://root-forum.cern.ch/t/looping-over-a-ttree/19899/5
    //https://root.cern/root/html534/tutorials/tree/hvector.C.html
//    _entry.SetBranchAddress("mc_pdg_id",&pdg_id);
//    _entry.GetBranch("mc_pdg_id")->GetEntry(index);
    
    using std::vector;
    vector<int> *mc_pdg_id= nullptr;
    _entry.SetBranchAddress("mc_pdg_id",&mc_pdg_id);
    _entry.GetBranch("mc_pdg_id")->GetEntry(index);
    R__ASSERT((*mc_pdg_id).size() > static_cast<size_t>(index));
    pdg_id = (*mc_pdg_id)[index];
    
}

TruthParticle::TruthParticle(const TruthParticle& particle_temp)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&particle_temp._entry)); //😬
    pdg_id = particle_temp.pdg_id;
    _index = particle_temp._index;
    Cluster_eta = particle_temp.Cluster_eta;
}

TruthParticle& TruthParticle::operator=(const TruthParticle& particle_temp)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&particle_temp._entry)); //😬
    pdg_id = particle_temp.pdg_id;
    _index = particle_temp._index;
    Cluster_eta = particle_temp.Cluster_eta;
    return *this;
}

TruthParticle::~TruthParticle() = default;

TLorentzVector TruthParticle::Vector() 
{
    TLorentzVector vec;
    vec.SetPtEtaPhiE(pt(),eta(),phi(),e());
    return vec;
}

int TruthParticle::barcode()
{
    vector<int> *mc_barcode= nullptr;
    _entry.SetBranchAddress("mc_barcode",&mc_barcode);
    _entry.GetBranch("mc_barcode")->GetEntry(_index);
    R__ASSERT((*mc_barcode).size() > static_cast<size_t>(_index));
    return (*mc_barcode)[_index];
}

int TruthParticle::parent_barcode()
{
    vector<int> *mc_parent_barcode= nullptr;
    _entry.SetBranchAddress("mc_parent_barcode",&mc_parent_barcode);
    _entry.GetBranch("mc_parent_barcode")->GetEntry(_index);
    R__ASSERT((*mc_parent_barcode).size() > static_cast<size_t>(_index));
    return (*mc_parent_barcode)[_index];
}

int TruthParticle::status_code()
{
    vector<int> *mc_status= nullptr;
    _entry.SetBranchAddress("mc_status",&mc_status);
    _entry.GetBranch("mc_status")->GetEntry(_index);
    R__ASSERT((*mc_status).size() > static_cast<size_t>(_index));
    return (*mc_status)[_index];
}

double TruthParticle::pt()
{
    vector<double> *mc_pt= nullptr;
    _entry.SetBranchAddress("mc_pt",&mc_pt);
    _entry.GetBranch("mc_pt")->GetEntry(_index);
    R__ASSERT((*mc_pt).size() > static_cast<size_t>(_index));
    return (*mc_pt)[_index];
}

double TruthParticle::charge()
{
    vector<double> *mc_charge= nullptr;
    _entry.SetBranchAddress("mc_charge",&mc_charge);
    _entry.GetBranch("mc_charge")->GetEntry(_index);
    R__ASSERT((*mc_charge).size() > static_cast<size_t>(_index));
    return (*mc_charge)[_index];
}

double TruthParticle::eta()
{
    vector<double> *mc_eta= nullptr;
    _entry.SetBranchAddress("mc_eta",&mc_eta);
    _entry.GetBranch("mc_eta")->GetEntry(_index);
    R__ASSERT((*mc_eta).size() > static_cast<size_t>(_index));
    return (*mc_eta)[_index];
}

double TruthParticle::phi()
{
    vector<double> *mc_phi= nullptr;
    _entry.SetBranchAddress("mc_phi",&mc_phi);
    _entry.GetBranch("mc_phi")->GetEntry(_index);
    R__ASSERT((*mc_phi).size() > static_cast<size_t>(_index));
    return (*mc_phi)[_index];
}

double TruthParticle::e()
{
    vector<double> *mc_e= nullptr;
    _entry.SetBranchAddress("mc_e",&mc_e);
    _entry.GetBranch("mc_e")->GetEntry(_index);
    R__ASSERT((*mc_e).size() > static_cast<size_t>(_index));
    return (*mc_e)[_index];
}

double TruthParticle::m()
{
    vector<double> *mc_mass= nullptr;
    _entry.SetBranchAddress("mc_mass",&mc_mass);
    _entry.GetBranch("mc_mass")->GetEntry(_index);
    R__ASSERT((*mc_mass).size() > static_cast<size_t>(_index));
    return (*mc_mass)[_index];
}

TruthParticle::operator string() 
{
    using s = string; //alias
    
    return
    (
    s("TruthParticle(pdg_id=")+to_string(pdg_id)
    +s(", pT=") + to_string(pt())
    +s(", charge=") + to_string(charge())
    +s(", eta=") + to_string(eta())
    +s(", phi=") + to_string(phi())
    );
}

int TruthParticle::id_(){return true;}
int TruthParticle::id_loose(){return true;}
int TruthParticle::id_tight(){return true;}

const string TruthParticle::PREFIX = "mc";

//                *****************
//                *    Electron   *
//                *****************

Electron::Electron() = default;

Electron::Electron(TChain* entry, int index, double pt, double energy, const char* name, const char* title)
: TruthParticle::TruthParticle(entry, index, name, title),
  _systematic_pt{pt},
  _systematic_energy{energy}
{
    
}


Electron::Electron(const Electron& other)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&other._entry)); //😬
    pdg_id = other.pdg_id;
    _index = other._index;
    Cluster_eta = other.Cluster_eta;
    _systematic_pt = other._systematic_pt;
    _systematic_energy = other._systematic_energy;
}

Electron& Electron::operator=(const Electron& other)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&other._entry)); //😬
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
    if (_systematic_pt)
    {
        return _systematic_pt;
    }
    vector<double> *electron_pt= nullptr;
    _entry.SetBranchAddress("electron_pt",&electron_pt);
    _entry.GetBranch("electron_pt")->GetEntry(_index);
    R__ASSERT((*electron_pt).size() > static_cast<size_t>(_index));
    return (*electron_pt)[_index];
}

double Electron::e()
{
    if (_systematic_energy)
    {
        return _systematic_energy;
    }
    vector<double> *electron_e= nullptr;
    _entry.SetBranchAddress("electron_e",&electron_e);
    _entry.GetBranch("electron_e")->GetEntry(_index);
    R__ASSERT((*electron_e).size() > static_cast<size_t>(_index));
    return (*electron_e)[_index];
}

double Electron::eta()
{
    vector<double> *electron_eta= nullptr;
    _entry.SetBranchAddress("electron_eta",&electron_eta);
    _entry.GetBranch("electron_eta")->GetEntry(_index);
    R__ASSERT((*electron_eta).size() > static_cast<size_t>(_index));
    return (*electron_eta)[_index];
}

double Electron::phi()
{
    vector<double> *electron_phi= nullptr;
    _entry.SetBranchAddress("electron_phi",&electron_phi);
    _entry.GetBranch("electron_phi")->GetEntry(_index);
    R__ASSERT((*electron_phi).size() > static_cast<size_t>(_index));
    return (*electron_phi)[_index];
}

int Electron::id_()
{
    vector<int> *electron_id= nullptr;
    _entry.SetBranchAddress("electron_id",&electron_id);
    _entry.GetBranch("electron_id")->GetEntry(_index);
    R__ASSERT((*electron_id).size() > static_cast<size_t>(_index));
    return (*electron_id)[_index];
}

double Electron::isolation()
{
    vector<double> *electron_isolation= nullptr;
    _entry.SetBranchAddress("electron_isolation",&electron_isolation);
    _entry.GetBranch("electron_isolation")->GetEntry(_index);
    R__ASSERT((*electron_isolation).size() > static_cast<size_t>(_index));
    return (*electron_isolation)[_index];
}

double Electron::d0()
{
    vector<double> *electron_d0= nullptr;
    _entry.SetBranchAddress("electron_d0",&electron_d0);
    _entry.GetBranch("electron_d0")->GetEntry(_index);
    R__ASSERT((*electron_d0).size() > static_cast<size_t>(_index));
    return (*electron_d0)[_index];
}

double Electron::z0()
{
    vector<double> *electron_z0= nullptr;
    _entry.SetBranchAddress("electron_z0",&electron_z0);
    _entry.GetBranch("electron_z0")->GetEntry(_index);
    R__ASSERT((*electron_z0).size() > static_cast<size_t>(_index));
    return (*electron_z0)[_index];
}

Electron::operator string()
{
    using s = string; //alias
    
    return
    (
    s("TruthParticle(pdg_id=")+to_string(PDG_ID)
    +s(", pT=") + to_string(pt())
    +s(", charge=") + to_string(charge())
    +s(", eta=") + to_string(eta())
    +s(", phi=") + to_string(phi())
    );
}

const int Electron::PDG_ID = 11;
const string Electron::PREFIX = "electron";

//                *****************
//                *    Photon     *
//                *****************

Photon::Photon() = default;

Photon::Photon(TChain* entry, int index, double pt, double energy, const char* name, const char* title)
: TruthParticle::TruthParticle(entry, index, name, title),
  _systematic_pt{pt},
  _systematic_energy{energy}
{
    
}

Photon::Photon(const Photon& other)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&other._entry)); //😬
    pdg_id = other.pdg_id;
    _index = other._index;
    Cluster_eta = other.Cluster_eta;
    _systematic_pt = other._systematic_pt;
    _systematic_energy = other._systematic_energy;
}

Photon& Photon::operator=(const Photon& other)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&other._entry)); //😬
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
    vector<double> *photon_pt= nullptr;
    _entry.SetBranchAddress("photon_pt",&photon_pt);
    _entry.GetBranch("photon_pt")->GetEntry(_index);
    R__ASSERT((*photon_pt).size() > static_cast<size_t>(_index));
    return (*photon_pt)[_index];
}

double Photon::e()
{
    if (_systematic_energy)
    {
        return _systematic_energy;
    }
    vector<double> *photon_e= nullptr;
    _entry.SetBranchAddress("photon_e",&photon_e);
    _entry.GetBranch("photon_e")->GetEntry(_index);
    R__ASSERT((*photon_e).size() > static_cast<size_t>(_index));
    return (*photon_e)[_index];
}

double Photon::eta()
{
    vector<double> *photon_eta= nullptr;
    _entry.SetBranchAddress("photon_eta",&photon_eta);
    _entry.GetBranch("photon_eta")->GetEntry(_index);
    R__ASSERT((*photon_eta).size() > static_cast<size_t>(_index));
    return (*photon_eta)[_index];
}

double Photon::phi()
{
    vector<double> *photon_phi= nullptr;
    _entry.SetBranchAddress("photon_phi",&photon_phi);
    _entry.GetBranch("photon_phi")->GetEntry(_index);
    R__ASSERT((*photon_phi).size() > static_cast<size_t>(_index));
    return (*photon_phi)[_index];
}

double Photon::m()
{
    return 0;
}

double Photon::isolation()
{
    vector<double> *photon_etcone40= nullptr;
    _entry.SetBranchAddress("photon_etcone40",&photon_etcone40);
    _entry.GetBranch("photon_etcone40")->GetEntry(_index);
    R__ASSERT((*photon_etcone40).size() > static_cast<size_t>(_index));
    return (*photon_etcone40)[_index];
}

int Photon::id_()
{
    vector<int> *photon_id= nullptr;
    _entry.SetBranchAddress("photon_id",&photon_id);
    _entry.GetBranch("photon_id")->GetEntry(_index);
    R__ASSERT((*photon_id).size() > static_cast<size_t>(_index));
    return (*photon_id)[_index];
}

int Photon::id_loose()
{
    vector<int> *photon_id_loose= nullptr;
    _entry.SetBranchAddress("photon_id_loose",&photon_id_loose);
    _entry.GetBranch("photon_id_loose")->GetEntry(_index);
    R__ASSERT((*photon_id_loose).size() > static_cast<size_t>(_index));
    return (*photon_id_loose)[_index];
}

int Photon::id_tight()
{
    vector<int> *photon_id_tight= nullptr;
    _entry.SetBranchAddress("photon_id_tight",&photon_id_tight);
    _entry.GetBranch("photon_id_tight")->GetEntry(_index);
    R__ASSERT((*photon_id_tight).size() > static_cast<size_t>(_index));
    return (*photon_id_tight)[_index];
}

double Photon::cluster_eta()
{
    vector<int> *photon_cluster_eta_be_2= nullptr;
    _entry.SetBranchAddress("photon_cluster_eta_be_2",&photon_cluster_eta_be_2);
    _entry.GetBranch("photon_cluster_eta_be_2")->GetEntry(_index);
    R__ASSERT((*photon_cluster_eta_be_2).size() > static_cast<size_t>(_index));
    return (*photon_cluster_eta_be_2)[_index];
}

int Photon::id_nn()
{
    vector<int> *photon_id_nn= nullptr;
    _entry.SetBranchAddress("photon_id_nn",&photon_id_nn);
    _entry.GetBranch("photon_id_nn")->GetEntry(_index);
    R__ASSERT((*photon_id_nn).size() > static_cast<size_t>(_index));
    return (*photon_id_nn)[_index];
}

Photon::operator string()
{
    using s = string; //alias
    
    return
    (
    s("TruthParticle(pdg_id=")+to_string(PDG_ID)
    +s(", pT=") + to_string(pt())
    +s(", charge=") + to_string(charge())
    +s(", eta=") + to_string(eta())
    +s(", phi=") + to_string(phi())
    );
}

const int Photon::PDG_ID = 22;
const string Photon::PREFIX = "photon";

//                *****************
//                *    Cluster    *
//                *****************

Cluster::Cluster() = default;

Cluster::Cluster(TChain* entry, int index, const char* name, const char* title)
:   PhysicsObject::PhysicsObject(entry, index, name, title)
{
    
}

Cluster::Cluster(const Cluster& particle_temp)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&particle_temp._entry)); //😬
    _index = particle_temp._index;
}

Cluster& Cluster::operator=(const Cluster& other)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&other._entry)); //😬
    _index = other._index;
    return *this;
}

double Cluster::pt()
{
    vector<double> *cluster_pt= nullptr;
    _entry.SetBranchAddress("cluster_pt",&cluster_pt);
    _entry.GetBranch("cluster_pt")->GetEntry(_index);
    R__ASSERT((*cluster_pt).size() > static_cast<size_t>(_index));
    return (*cluster_pt)[_index];
}

double Cluster::eta()
{
    vector<double> *cluster_eta= nullptr;
    _entry.SetBranchAddress("cluster_eta",&cluster_eta);
    _entry.GetBranch("cluster_eta")->GetEntry(_index);
    R__ASSERT((*cluster_eta).size() > static_cast<size_t>(_index));
    return (*cluster_eta)[_index];
}

double Cluster::phi()
{
    vector<double> *cluster_phi= nullptr;
    _entry.SetBranchAddress("cluster_phi",&cluster_phi);
    _entry.GetBranch("cluster_phi")->GetEntry(_index);
    R__ASSERT((*cluster_phi).size() > static_cast<size_t>(_index));
    return (*cluster_phi)[_index];
}

double Cluster::e()
{
    vector<double> *cluster_e= nullptr;
    _entry.SetBranchAddress("cluster_e",&cluster_e);
    _entry.GetBranch("cluster_e")->GetEntry(_index);
    R__ASSERT((*cluster_e).size() > static_cast<size_t>(_index));
    return (*cluster_e)[_index];
}

const string Cluster::PREFIX = "cluster";

//                *****************
//                *    Track      *
//                *****************

Track::Track() = default;


Track::Track(TChain* entry, int index, const char* name, const char* title)
:   PhysicsObject::PhysicsObject(entry, index, name, title)
{
    
}

Track::Track(const Track& particle_temp)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&particle_temp._entry)); //😬
    _index = particle_temp._index;
}

Track& Track::operator=(const Track& other)
{
    _entry.Reset();
    _entry.Add(const_cast<TChain*>(&other._entry)); //😬
    _index = other._index;
    return *this;
}

double Track::pt()
{
    vector<double> *track_pt= nullptr;
    _entry.SetBranchAddress("track_pt",&track_pt);
    _entry.GetBranch("track_pt")->GetEntry(_index);
    R__ASSERT((*track_pt).size() > static_cast<size_t>(_index));
    return (*track_pt)[_index];
}

double Track::charge()
{
    vector<double> *track_charge= nullptr;
    _entry.SetBranchAddress("track_charge",&track_charge);
    _entry.GetBranch("track_charge")->GetEntry(_index);
    R__ASSERT((*track_charge).size() > static_cast<size_t>(_index));
    return (*track_charge)[_index];
}

double Track::eta()
{
    vector<double> *track_eta= nullptr;
    _entry.SetBranchAddress("track_eta",&track_eta);
    _entry.GetBranch("track_eta")->GetEntry(_index);
    R__ASSERT((*track_eta).size() > static_cast<size_t>(_index));
    return (*track_eta)[_index];
}

double Track::phi()
{
    vector<double> *track_phi= nullptr;
    _entry.SetBranchAddress("track_phi",&track_phi);
    _entry.GetBranch("track_phi")->GetEntry(_index);
    R__ASSERT((*track_phi).size() > static_cast<size_t>(_index));
    return (*track_phi)[_index];
}

double Track::e()
{
    vector<double> *track_e= nullptr;
    _entry.SetBranchAddress("track_e",&track_e);
    _entry.GetBranch("track_e")->GetEntry(_index);
    R__ASSERT((*track_e).size() > static_cast<size_t>(_index));
    return (*track_e)[_index];
}

int Track::num_pixel_hits()
{
    vector<int> *track_num_pixel_hits= nullptr;
    _entry.SetBranchAddress("track_num_pixel_hits",&track_num_pixel_hits);
    _entry.GetBranch("track_num_pixel_hits")->GetEntry(_index);
    R__ASSERT((*track_num_pixel_hits).size() > static_cast<size_t>(_index));
    return (*track_num_pixel_hits)[_index];
}

int Track::num_sct_hits()
{
    vector<int> *track_num_sct_hits= nullptr;
    _entry.SetBranchAddress("track_num_sct_hits",&track_num_sct_hits);
    _entry.GetBranch("track_num_sct_hits")->GetEntry(_index);
    R__ASSERT((*track_num_sct_hits).size() > static_cast<size_t>(_index));
    return (*track_num_sct_hits)[_index];
}

Track::operator string()
{
    using s = string; //alias
    
    return
    (
     s("Track(pT=")+to_string(pt())
    +s(", charge=") + to_string(charge())
    +s(", eta=") + to_string(eta())
    +s(", phi=") + to_string(phi())
    );
}

const string Track::PREFIX = "track";
