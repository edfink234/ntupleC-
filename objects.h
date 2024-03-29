#ifndef OBJECTSH
#define OBJECTSH

#include <vector>
#include <string>

//#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "TChain.h"

using namespace ROOT::Math;

class TruthParticle;

/*
 PhysicsObject: base class for all objects of actual interest
 */
class PhysicsObject
{
protected:
    static const std::string PREFIX;
    int _index;
public:
    PhysicsObject();
    PhysicsObject(int);
    virtual ~PhysicsObject() = 0;
//    virtual TLorentzVector Vector();
    virtual PtEtaPhiEVector Vector();
    double acoplanarity(TruthParticle&);
    double delta_r( TruthParticle&);
    operator std::string();

    virtual double pt() = 0;
    virtual double eta() = 0;
    virtual double phi() = 0;
    virtual double e() = 0;
};

/*
 TruthParticle: class that stores info for the particles that were actually involved
 in a particle physics event, only availabe in simulation
 */
class TruthParticle : public PhysicsObject
{
protected:
    double Cluster_eta;
    static const std::string PREFIX;
    friend class FileReader;
    static std::vector<int>* mc_pdg_id;
    static std::vector<int>* mc_barcode;
    static std::vector<int>* mc_parent_barcode;
    static std::vector<int>* mc_status;
    static std::vector<double>* mc_pt;
    static std::vector<double>* mc_charge;
    static std::vector<double>* mc_eta;
    static std::vector<double>* mc_phi;
    static std::vector<double>* mc_e;
    static std::vector<double>* mc_mass;
public:
    int pdg_id;
    static void SetTruthParticle(TChain*);
    TruthParticle();
    TruthParticle(int);
    TruthParticle(const TruthParticle&);
    TruthParticle& operator=(const TruthParticle&);
    virtual ~TruthParticle();
//    TLorentzVector Vector() override;
    PtEtaPhiEVector Vector() override;
    int barcode();
    int parent_barcode();
    int status_code();
    double pt() override;
    virtual double charge();
    double eta() override;
    double phi() override;
    double e() override;
    virtual bool same_flavour(const TruthParticle &) final;
    virtual double m();
    operator std::string();
    
    virtual int id_();
    virtual int id_loose();
    virtual int id_tight();
};

/*
 Electron: class that stores info for reconstructed electrons
 */
class Electron final : public TruthParticle
{
private:
    static const std::string PREFIX;
    double _systematic_pt;
    double _systematic_energy;
    friend class FileReader;
    static std::vector<double>* electron_charge;
    static std::vector<double>* electron_pt;
    static std::vector<double>* electron_e;
    static std::vector<double>* electron_eta;
    static std::vector<double>* electron_phi;
    static std::vector<int>* electron_id;
    static std::vector<double>* electron_isolation;
    static std::vector<double>* electron_d0;
    static std::vector<double>* electron_z0;
    static std::vector<int>* electron_id_medium;
public:
    static const int PDG_ID;
    static void SetElectron(TChain*);
    Electron();
    Electron(int, double pt=0, double energy = 0);
    Electron(const Electron&);
    Electron& operator=(const Electron&);
    int Pdg_id();
    double pt() override;
    double charge() override;
    double e() override;
    double eta() override;
    double phi() override;
    int id_() override;
    double isolation();
    int id_medium();
    double d0();
    double z0();
    operator std::string();
};

/*
 Muon: class that stores info for reconstructed muons
 */
class Muon final : public TruthParticle
{
private:
    static const std::string PREFIX;
    double _systematic_pt;
    double _systematic_energy;
    friend class FileReader;
    static std::vector<double>* muon_charge;
    static std::vector<double>* muon_pt;
    static std::vector<double>* muon_e;
    static std::vector<double>* muon_eta;
    static std::vector<double>* muon_phi;

public:
    static const int PDG_ID;
    static void SetMuon(TChain*);
    Muon();
    Muon(int, double pt=0, double energy = 0);
    Muon(const Muon&);
    Muon& operator=(const Muon&);
    int Pdg_id();
    double pt() override;
    double charge() override;
    double e() override;
    double eta() override;
    double phi() override;
    operator std::string();
};

/*
 Photon: class that stores info for reconstructed photons
 */
class Photon final : public TruthParticle
{
private:
    static const std::string PREFIX;
    double _systematic_pt;
    double _systematic_energy;
    friend class FileReader;
    static std::vector<double>* photon_pt;
    static std::vector<double>* photon_e;
    static std::vector<double>* photon_eta;
    static std::vector<double>* photon_phi;
    static std::vector<double>* photon_etcone40;
    static std::vector<int>* photon_id;
    static std::vector<int>* photon_id_loose;
    static std::vector<int>* photon_id_tight;
    static std::vector<int>* photon_cluster_eta_be_2;
    static std::vector<int>* photon_id_nn;
public:
    static const int PDG_ID;
    static void SetPhoton(TChain*);
    Photon();
    Photon(int, double pt=0, double energy = 0);
    Photon(const Photon&);
    Photon& operator=(const Photon&);
    double pt() override;
    double e() override;
    double eta() override;
    double phi() override;
    double m() override;
    double isolation();
    int id_() override;
    int id_loose() override;
    int id_tight() override;
    double cluster_eta();
    int id_nn();
    operator std::string();
};

/*
 Cluster: class that stores info for reconstructed clusters
 */
class Cluster final : public PhysicsObject
{
private:
    static const std::string PREFIX;
    friend class FileReader;
    static std::vector<double>* cluster_pt;
    static std::vector<double>* cluster_eta;
    static std::vector<double>* cluster_phi;
    static std::vector<double>* cluster_e;
public:
    static void SetCluster(TChain*);
    Cluster();
    Cluster(int);
    Cluster(const Cluster&);
    Cluster& operator=(const Cluster&);
    double pt() override;
    double eta() override;
    double phi() override;
    double e() override;
};

/*
 Track: class that stores info for reconstructed tracks
 */
class Track final : public PhysicsObject
{
private:
    static const std::string PREFIX;
    friend class FileReader;
    static std::vector<double>* track_pt;
    static std::vector<double>* track_charge;
    static std::vector<double>* track_eta;
    static std::vector<double>* track_phi;
    static std::vector<double>* track_e;
    static std::vector<int>* track_num_pixel_hits;
    static std::vector<int>* track_num_sct_hits;
public:
    static void SetTrack(TChain*);
    Track();
    Track(const Track&);
    Track(int);
    Track& operator=(const Track&);
    double pt() override;
    double charge();
    double eta() override;
    double phi() override;
    double e() override;
    int num_pixel_hits();
    int num_sct_hits();
    operator std::string();
};

#endif
