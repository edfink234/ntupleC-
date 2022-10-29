#ifndef OBJECTSH
#define OBJECTSH

#include <vector>
#include <string>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"

class TruthParticle;

class PhysicsObject
{
protected:
    static const std::string PREFIX;
    int _index;
    int _entry_number;
public:
    PhysicsObject();
    PhysicsObject(int, int entry_number = 0, const char* name = "", const char* title = "");
    virtual ~PhysicsObject() = 0;
    virtual TLorentzVector Vector(); //renamed jic it clashes w/ std::vector
    double acoplanarity(TruthParticle&);
    double delta_r( TruthParticle&);
    operator std::string();
    
    //TODO: getattr equivalent in C++;
    //TODO: Add property class?: https://en.wikipedia.org/wiki/Property_%28programming%29#C++
    //override in TruthParticle!
    //https://stackoverflow.com/questions/46446652/is-there-any-point-in-using-override-when-overriding-a-pure-virtual-function

    virtual double pt() = 0;
    virtual double eta() = 0;
    virtual double phi() = 0;
    virtual double e() = 0;
};

class TruthParticle : public PhysicsObject
{
protected:
    double Cluster_eta;
    int pdg_id;
    static const std::string PREFIX;
    friend class FileReader;
public:
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
    
    static void SetTruthParticle(TChain*);
    TruthParticle();
    TruthParticle(int, int entry_number = 0, const char* name = "", const char* title = "");
    TruthParticle(const TruthParticle&);
    TruthParticle& operator=(const TruthParticle&);
    virtual ~TruthParticle();
    TLorentzVector Vector() override;
    int barcode();
    int parent_barcode();
    int status_code();
    double pt() override;
    double charge();
    double eta() override;
    double phi() override;
    double e() override;
    virtual double m();
    operator std::string();
    
    virtual int id_();
    virtual int id_loose();
    virtual int id_tight();
};

class Electron final : public TruthParticle
{
private:
    static const std::string PREFIX;
    static const int PDG_ID;
    double _systematic_pt;
    double _systematic_energy;
    
public:
    static std::vector<double>* electron_pt;
    static std::vector<double>* electron_e;
    static std::vector<double>* electron_eta;
    static std::vector<double>* electron_phi;
    static std::vector<int>* electron_id;
    static std::vector<double>* electron_isolation;
    static std::vector<double>* electron_d0;
    static std::vector<double>* electron_z0;
    
    static void SetElectron(TChain*);
    Electron();
    Electron(int, double pt=0, double energy = 0, int entry_number = 0, const char* name = "", const char* title = "");
    Electron(const Electron&);
    Electron& operator=(const Electron&);
    int Pdg_id();
    double pt() override;
    double e() override;
    double eta() override;
    double phi() override;
    int id_() override;
    double isolation();
    double d0();
    double z0();
    operator std::string();
};

class Photon final : public TruthParticle
{
private:
    static const std::string PREFIX;
    static const int PDG_ID;
    double _systematic_pt;
    double _systematic_energy;
    
public:
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
    
    static void SetPhoton(TChain*);
    Photon();
    Photon(int, double pt=0, double energy = 0, int entry_number = 0, const char* name = "", const char* title = "");
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

class Cluster final : public PhysicsObject
{
private:
    static const std::string PREFIX;
    
public:
    static std::vector<double>* cluster_pt;
    static std::vector<double>* cluster_eta;
    static std::vector<double>* cluster_phi;
    static std::vector<double>* cluster_e;
    
    static void SetCluster(TChain*);
    Cluster();
    Cluster(int, int entry_number = 0, const char* name = "", const char* title = "");
    Cluster(const Cluster&);
    Cluster& operator=(const Cluster&);
    double pt() override;
    double eta() override;
    double phi() override;
    double e() override;
};

class Track final : public PhysicsObject
{
private:
    static const std::string PREFIX;
public:
    static std::vector<double>* track_pt;
    static std::vector<double>* track_charge;
    static std::vector<double>* track_eta;
    static std::vector<double>* track_phi;
    static std::vector<double>* track_e;
    static std::vector<int>* track_num_pixel_hits;
    static std::vector<int>* track_num_sct_hits;
    
    static void SetTrack(TChain*);
    Track();
    Track(const Track&);
    Track(int, int entry_number = 0, const char* name = "", const char* title = "");
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
