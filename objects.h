#ifndef OBJECTSH
#define OBJECTSH

#include <cmath>
#include <vector>
#include <string>
#include "TLorentzVector.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
using std::string;

class TruthParticle;

class PhysicsObject
{
protected:
    TChain _entry;
    static const string PREFIX;
    int _index;
public:
    PhysicsObject();
    PhysicsObject(TChain*, int,  const char* name = "", const char* title = "");
    virtual ~PhysicsObject() = 0;
    virtual TLorentzVector Vector(); //renamed jic it clashes w/ std::vector
    
    double acoplanarity(TruthParticle&);
    double delta_r( TruthParticle&);
    operator string();
    
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
    static const string PREFIX;
    friend class Event;
public:
    TruthParticle();
    TruthParticle(TChain*, int, const char* name = "", const char* title = "");
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
    operator string();
    
    virtual int id_();
    virtual int id_loose();
    virtual int id_tight();
};

class Electron final : public TruthParticle
{
private:
    static const string PREFIX;
    static const int PDG_ID;
    double _systematic_pt;
    double _systematic_energy;
public:
    Electron();
    Electron(TChain*, int, double pt=0, double energy = 0, const char* name = "", const char* title = "");
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
    operator string();
};

class Photon final : public TruthParticle
{
private:
    static const string PREFIX;
    static const int PDG_ID;
    double _systematic_pt;
    double _systematic_energy;
public:
    Photon();
    Photon(TChain*, int, double pt=0, double energy = 0, const char* name = "", const char* title = "");
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
    operator string();
};

class Cluster final : public PhysicsObject
{
private:
    static const string PREFIX;
public:
    Cluster();
    Cluster(TChain*, int, const char* name = "", const char* title = "");
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
    static const string PREFIX;
public:
    Track();
    Track(const Track&);
    Track(TChain*, int, const char* name = "", const char* title = "");
    Track& operator=(const Track&);
    double pt() override;
    double charge();
    double eta() override;
    double phi() override;
    double e() override;
    int num_pixel_hits();
    int num_sct_hits();
    operator string();
};



#endif
