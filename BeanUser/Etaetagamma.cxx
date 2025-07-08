//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  Etaetagamma                                                         //
//                                                                      //
// Select gamma pi0 pi0 and gamma eta eta events                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <set>
#include <cstdarg>

#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TMath.h>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"


#include "CLHEP/Units/PhysicalConstants.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/TwoVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
using CLHEP::pi;

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "AbsCor/AbsCor.h"
#include "VertexFit/VertexDbSvc.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include "MagneticField/MagneticFieldSvc.h"

#include "DstEvtRecTracks.h"
#include "TofHitStatus.h"
#include "ReadDst.h"
#include "EventTag/EventTagSvc.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif


#define TAG_SRC_SIZE 15

//const double twopi = 6.2831853;
//const double pi = 3.1415927;
const double mpi = 0.13957;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
//const double velc = 29.9792458;  tof_path unit in cm.
const double velc = 299.792458;   // tof path unit in mm


static double mpi0 = 0.1349766;
static double meta = 0.547853;
// static double mphi = 1.01941;
static double momega = 0.78257;

typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;

typedef std::vector<DstEvtRecTracks*> Vtracks;
/////////////////////////////////////////////////////////////////////////////

static std::vector<TNtuple* > m_tuple;

//ReadBeamParFromDb m_reader;
// Declare r0, z0 cut for charged tracks
static double m_vr0cut;
static double m_vz0cut;

//Declare energy, dphi, dthe cuts for fake gamma's
static double m_energyThreshold;
static double m_gammaPhiCut;
static double m_gammaThetaCut;
static double m_gammaAngleCut;

static int m_test4C;
static int m_test5C;

static int m_checkDedx;
static int m_checkTof;

static AbsCor* m_abscor = 0;
/////////////////////////////////////////////////////////////////////////////
static TH1D * hTotalMass;
static TH1I * hNGamma;

static TH1I * hNPi0Candidates;
static TH1I * hNEtaCandidates;

static TH1D * hEMCTime;

static TH1D * hMGammaGamma;
static TH1D * hMGammaGammaCut;


static TH1D * hPi0Chisq;
static TH1D * hEtaChisq;

static TH1D * hGood2EtaChisq;
static TH1D * hBad2EtaChisq;

static TH1D * hGood2EtaChisq4C;
static TH1D * hBad2EtaChisq4C;

static TH1D * hGood2EtaChisq2C;
static TH1D * hBad2EtaChisq2C;

static TH1D * hDeltaEta;
static TH1D * hDeltaPi;

static TH1D * hMRecoilEtaPair;
static TH1D * hMRecoilPi0Pair;
static TH1D * hMRecoilPi0Gamma;


static TH1D * hMPi0Pi0;
static TH1D * hMPi0Gamma;
static TH1D * hMEtaGamma;
static TH1D * hMEtaEta;


static TH1D * hDecayAngleCosPi0;
static TH1D * hDecayAngleCosEta;
static TH1D * hDecayAnglePi0;

static TH1D * hDecayAngleEta;
static TH1D * hGood2EtaDecayAngleCos;
static TH1D * hBad2EtaDecayAngleCos;


static TH1D * hDecayAnglePair;
static TH1D * hDecayAngleCosPair;

static TH2D * hDecayAngleCosEta_Tag;
static TH2D * hEtaChisq_Tag;
// static TH1D * hDecayAngleCosEtaPi0;

//DEBUG:
static TH1I * hNCut;
static TH1I * hNGoodTagCut;

static TH1D * hPi0RadPhotonEnergy;
static TH1D * hEtaRadPhotonEnergy;

static TH1D * hGood2EtaRadPhotonEnergy;
static TH1D * hBad2EtaRadPhotonEnergy;

static TH1D * hPhotonEnergy;

static TH1I *   hPi0Tags;
static TH1I *   hEtaTags;

static TH1I *   hEvrPi0Tags;
static TH1I *   hEvrEtaTags;

static TH1I *   hGood2EtaTags;
static TH1I *   hBad2EtaTags;


static TH1I *   hInitialTags;


static TH2D * hEGamma_Theta;

static TH1D * hGood2EtaDelta;
static TH1D * hBad2EtaDelta;

static TH1D * hBad2EtaEErr;
static TH1D * hGood2EtaEErr;

static TH1I * hGood2EtaNGamma;
static TH1I * hBad2EtaNGamma;

static TH1I * hGood2EtaNCombinations;
static TH1I * hBad2EtaNCombinations;


static TH1I * hGood2EtaNPi0Candidates;
static TH1I * hBad2EtaNPi0Candidates;

static TH1I * hGood2EtaNEtaCandidates;
static TH1I * hBad2EtaNEtaCandidates;


static TH1D * hGood2EtaMEtaEta;
static TH1D * hBad2EtaMEtaEta;

static TH1D * hGood2EtaMGammaPair;
static TH1D * hBad2EtaMGammaPair;

static TH1D * hGood2EtaMGammaGammaAll;
static TH1D * hBad2EtaMGammaGammaAll;

static TH1D * hGood2EtaDeltaMggMpi0;
static TH1D * hBad2EtaDeltaMggMpi0;

static TH2D * hBad2EtaChisq4C_NGamma;

//pi0
static TH1I *   hGood2Pi0Tags;
static TH1I *   hBad2Pi0Tags;

static EventTagSvc * m_EventTagSvc;

// static bool cut0 = 0;

//~ map<unsigned int, string> tagsInvolved;

//~ static TH2D * hMGammaGamma_chisq;

/////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
void EtaetagammaStartJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
    if( selector->Verbose() ) cout << " EtaEtagammaStartJob() " << endl;

    // init Absorption Correction
    m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

    // Initial values
    m_vr0cut=1.0;
    m_vz0cut=5.0;

    m_energyThreshold=0.04;
    m_gammaPhiCut=20.0;
    m_gammaThetaCut=20.0;
    m_gammaAngleCut=20.0;

    m_test4C = 1;
    m_test5C = 1;

    m_checkDedx = 1;
    m_checkTof = 1;

    if( selector->Verbose() ) {
        cout << endl;
        cout << " Rhopi internal parameters: " << endl;
        cout << " Rhopi.Vr0cut =            " << m_vr0cut << endl;
        cout << " Rhopi.Vz0cut =            " << m_vz0cut << endl;
        cout << " Rhopi.EnergyThreshold =   " << m_energyThreshold << endl;
        cout << " Rhopi.GammaPhiCut =       " << m_gammaPhiCut << endl;
        cout << " Rhopi.GammaThetaCut =     " << m_gammaThetaCut << endl;
        cout << " Rhopi.GammaAngleCut =     " << m_gammaAngleCut << endl;
        cout << " Rhopi.Test4C =            " << m_test4C << endl;
        cout << " Rhopi.Test5C =            " << m_test5C << endl;
        cout << " Rhopi.CheckDedx =         " << m_checkDedx << endl;
        cout << " Rhopi.CheckTof =          " << m_checkTof << endl;
        cout << endl;
    }

    // We have to initialize Magnetic field
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    MagneticFieldSvc* mf = MagneticFieldSvc::instance();
    if( !(mf->GetPath()).empty() ) {
    cerr << " RhopiStartJob::WARNING:"
         << " MagneticFieldSvc has already been initialized" << endl
         << "                         path = " << mf->GetPath() << endl;
    }
    // set path to directory with magnetic fields tables
    mf->SetPath( selector->AbsPath("Analysis/MagneticField") );
    mf->UseDBFlag(false); // like in the boss program
    mf->RunMode(3); // like in the boss program


    // Book n-tuples
    m_tuple.resize(15,0);

    // register in selector to save in given directory
    VecObj ntuples(m_tuple.begin(),m_tuple.end());
    selector->RegInDir(ntuples,"EtaEtagamma");


    VecObj his1;
    his1.resize(100);
    VecObj his2;
    his2.resize(100);
    int i = 0;

    hNCut = new TH1I("hNCut","N of cut(0-input, 1: neutral, 2: good photons, 3: 5or6 photons",11,0,10);

    hNGoodTagCut = new TH1I("hNGoodTagCut","N of good tags versus cut number",11,0,10); hNGoodTagCut->Fill("blah",0.0);


    hTotalMass =  new TH1D("total_mass", "Total energy of photons selected", 200, 0, 4);
    hNGamma =  new TH1I("hNGamma", "N gamma", 16, 0, 15);

    hGood2EtaNGamma =  new TH1I("hGood2EtaNGamma", "N gamma", 16, 0, 15);
    hBad2EtaNGamma =  new TH1I("hBad2EtaNGamma", "N gamma", 16, 0, 15);

    hBad2EtaChisq4C_NGamma  = new TH2D("hBad2EtaChisq4C_NGamma", "#chi_{4C}^2 vs N_#gamma",5,4,8,50,0,200);

    hGood2EtaNCombinations =  new TH1I("hGood2EtaNCombinations", "N combinations", 16, 0, 15);
    hBad2EtaNCombinations =  new TH1I("hBad2EtaNCombinations", "N combinations", 16, 0, 15);

    hGood2EtaNPi0Candidates =  new TH1I("hGood2EtaNPi0Candidates", "hGood2EtaNPi0Candidates", 6, 0, 5);
    hBad2EtaNPi0Candidates  =  new TH1I("hBad2EtaNPi0Candidates", "hBad2EtaNPi0Candidates", 6, 0, 5);
    hGood2EtaNEtaCandidates =  new TH1I("hGood2EtaNEtaCandidates", "hGood2EtaNEtaCandidates", 6, 0, 5);
    hBad2EtaNEtaCandidates =  new TH1I("hBad2EtaNEtaCandidates", "Bad2EtaNEtaCandidates", 6, 0, 5);

    hNPi0Candidates =  new TH1I("hNPi0Candidates", "N pi0 pairs", 16, 0, 15);
    hNEtaCandidates =  new TH1I("hNEtaCandidates", "N eta pairs", 16, 0, 15);



    hMGammaGamma =  new TH1D("hMGammaGamma", "M_{#gamma#gamma} (GeV)", 500, 0,1);
    hMGammaGammaCut =  new TH1D("hMGammaGammaCut", "cut M_{#gamma#gamma} (GeV)", 500, 0,1);

    hGood2EtaMGammaPair =  new TH1D("hGood2EtaMGammaPair", "M_{#gamma#gamma} (GeV)", 500, 0,1);
    hBad2EtaMGammaPair =  new TH1D("hBad2EtaMGammaPair", "M_{#gamma#gamma} (GeV)", 500, 0,1);

    hGood2EtaMGammaGammaAll =  new TH1D("hGood2EtaMGammaGammaAll", "M_{#gamma#gamma} (GeV) (all cobinations)", 500, 0,1);
    hBad2EtaMGammaGammaAll =  new TH1D("hBad2EtaMGammaGammaAll", "M_{#gamma#gamma} (GeV) (all combinations)", 500, 0,1);

    hGood2EtaDeltaMggMpi0 =  new TH1D("hGood2EtaDeltaMggMpi0", "min of M_{#gamma#gamma} - M_{#pi_0} (GeV) ", 100, -0.2,0.2);
    hBad2EtaDeltaMggMpi0 =  new TH1D("hBad2EtaDeltaMggMpi0", "min of M_{#gamma#gamma} - M_{#pi_0} (GeV) ", 100, -0.2,0.2);


    hPi0Chisq =  new TH1D("hPi0Chisq", "#chi^2", 200, 0,200);
    hEtaChisq =  new TH1D("hEtaChisq", "#chi^2", 200, 0,200);

    hGood2EtaChisq =  new TH1D("hGood2EtaChisq", "#chi^2 of right selected eta-eta-gamma", 200, 0,200);
    hBad2EtaChisq =  new TH1D("hBad2EtaChisq", "#chi^2 of misselected eta-eta-gamma ", 200, 0,200);

    hGood2EtaChisq4C =  new TH1D("hGood2EtaChisq4C", "4C #chi^2 of right selected eta-eta-gamma", 200, 0,200);
    hBad2EtaChisq4C =  new TH1D("hBad2EtaChisq4C", "%4C CHI^2 of misselected eta-eta-gamma ", 200, 0,200);

    hGood2EtaChisq2C =  new TH1D("hGood2EtaChisq2C", "2C #chi^2 of right selected eta-eta-gamma", 200, 0,200);
    hBad2EtaChisq2C =  new TH1D("hBad2EtaChisq2C", "%2C #chi^2 of misselected eta-eta-gamma ", 200, 0,200);



    hMRecoilEtaPair =  new TH1D("hMRecoilEtaPair", "Recoil mass of eta-eta pair(GeV)", 400, 0,4);
    hMRecoilPi0Pair =  new TH1D("hMRecoilPi0Pair", "Recoil mass of pi0-pi0 pair(GeV)", 400, 0,4);
    hMRecoilPi0Gamma = new TH1D("hMRecoilPi0Gamma", "Recoil mass of pi0-gammma (GeV)", 400, 0,4);

    //~ hPi0ChisqMin =  new TH1D("hPi0ChisqMin", "Min #chi^2 of 5 photons", 200, 0,200);
    hEMCTime = new TH1D("hEMCTime", "Time of EMC shower",101,0,100);
    hDeltaEta = new TH1D("hDeltaEta","%DELTA Eta",200,0,0.5);
    hDeltaPi = new TH1D("hDeltaPi","%DELTA Pi",200,0,0.5);

    hPhotonEnergy = new TH1D("hPhotonEnergy","Energy of photons",300,0,3);

    hGood2EtaEErr = new TH1D("hGood2EtaEErr","Error Energy of photons  of properly selected eta-eta-gamma",100,0,0.2);
    hBad2EtaEErr = new TH1D("hBad2EtaEErr","Erro Energy of wrongly selected eta-eta-gamma",100,0,0.2);

    hPi0RadPhotonEnergy = new TH1D("hPi0RadPhotonEnergy","Energy of radiative photon",300,0,3);
    hEtaRadPhotonEnergy = new TH1D("hEtaRadPhotonEnergy","Energy of radiative photon",300,0,3);

    hGood2EtaRadPhotonEnergy = new TH1D("hGood2EtaRadPhotonEnergy","Energy of radiative photon good",300,0,3);
    hBad2EtaRadPhotonEnergy = new TH1D("hBad2EtaRadPhotonEnergy","Energy of radiative photon bad",300,0,3);






    hEGamma_Theta = new TH2D("hEGamma_Theta", "EGamma vs Theta",200,0,3.14,100,0,2);
    //~ hMGammaGamma_chisq = new TH2D("hMGammaGamma_chisq", "M_{#gamma #gamma} vs #chi^2",200,0,500,200,0,1);

    hDecayAngleCosPi0 = new TH1D("hDecayAngleCosPi0", "cos DecayAngle of Pi0 candidates", 200, 0, 1);
    hDecayAngleCosEta = new TH1D("hDecayAngleCosEta", "cos DecayAngle of eta candidates", 200,0, 1);

    hGood2EtaDecayAngleCos = new TH1D("hGood2EtaDecayAngleCos", "cos DecayAngle of properly selected 2eta+gamma", 200,0, 1);
    hBad2EtaDecayAngleCos = new TH1D("hBad2EtaDecayAngleCos", "cos DecayAngle of wrongly selected 2eta+gamma", 200,0, 1);





    hDecayAngleCosEta_Tag = new TH2D("hDecayAngleCosEta_Tag", "cos DecayAngle of eta pair found versus tag", 50,0, 1, 10,1,0);
    hDecayAnglePi0 = new TH1D("hDecayAnglePi0", "DecayAngle of Pi0 candidates", 200, 0, 3.14);
    hDecayAngleEta = new TH1D("hDecayAngleEta", "DecayAngle of eta candidates", 200, 0, 3.14);

    hDecayAngleCosPair = new TH1D("hDecayAngleCosPair", "DecayAngle of pairs", 200, 0, 1);

    hDecayAnglePair = new TH1D("hDecayAnglePair", "DecayAngle of pairs", 200, -0, 3.14);

    hMPi0Pi0 = new TH1D("hMPi0Pi0", "Invariant mass of pi0-pi0 pair", 200, 0,4);
    hMEtaEta = new TH1D("hMEtaEta", "Invariant mass of eta-eta pair", 200, 0,4);

    hGood2EtaMEtaEta = new TH1D("hGood2EtaMEtaEta", "Invariant mass of eta-eta pair good ", 200, 0,4);
    hBad2EtaMEtaEta = new TH1D("hBad2EtaMEtaEta", "Invariant mass of eta-eta pair bad", 200, 0,4);

    hMPi0Gamma = new TH1D("hMPi0Gamma", "Invariant mass of pi0-gamma pair (6C-fit passed)", 200, 0,4);
    hMEtaGamma = new TH1D("hMEtaGamma", "Invariant mass of eta-gamma pair (6C-fit passed)", 200, 0,4);

    hPi0Tags = new TH1I("hPi0Tags", "Tags involved (6C-fit passed)", 11, 1,0); // lower bound greater than upper bound for automatically bound calculating (as stated in manual)
    hPi0Tags->Fill("blah",0.0);


    hGood2Pi0Tags = new TH1I("hGood2Pi0Tags", "Tags involved (6C-fit passed)", 11, 1,0);
    hGood2Pi0Tags->Fill("blah",0.0); //one MUST keep this for correct merging!

    hBad2Pi0Tags = new TH1I("hBad2Pi0Tags", "Tags involved (6C-fit passed)", 11, 1,0);
    hBad2Pi0Tags->Fill("blah",0.0); //one MUST keep this for correct merging!

    hEtaTags = new TH1I("hEtaTags", "Tags involved (6C-fit passed)", 11, 1,0); // lower bound greater than upper bound for automatically bound calculating (as stated in manual)
    hEtaTags->Fill("blah",0.0); //one MUST keep this for correct merging!

    hInitialTags = new TH1I("hInitialTags", "Tags involved (before any cuts done)", 11, 0,10); // lower bound greater than upper bound for automatically bound calculating (as stated in manual)
    hInitialTags->Fill("blah",0.0); //one MUST keep this for correct merging!

    hEvrPi0Tags = new TH1I("hEvrPi0Tags", "Possible tags for pi0, algorithm from mc", 11, 1,0); // lower bound greater than upper bound for automatically bound calculating (as stated in manual)
    hEvrPi0Tags->Fill("blah",0.0); //one MUST keep this for correct merging!

    hEvrEtaTags = new TH1I("hEvrEtaTags", "Possible tags for eta, algorithm from mc", 11, 1,0); // lower bound greater than upper bound for automatically bound calculating (as stated in manual)
    hEvrEtaTags->Fill("blah",0.0); //one MUST keep this for correct merging!

    hGood2EtaTags = new TH1I("hGood2EtaTags", "Possible tags for eta good", 11, 1,0); //
    hGood2EtaTags->Fill("blah",0.0); //one MUST keep this for correct merging!

    hBad2EtaTags = new TH1I("hBad2EtaTags", "Possible tags for eta bad", 11, 1,0); //
    hBad2EtaTags->Fill("blah",0.0); //one MUST keep this for correct merging!

    hEtaChisq_Tag = new TH2D("hEtaChisq_Tag", "eta chisq versus tag", 30,0, 200, 1,1,0);

    hGood2EtaDelta = new TH1D("hGood2EtaDelta","delta_eta for of properly selected 2eta+gamma ",200,0,0.5);
    hBad2EtaDelta = new TH1D("hBad2EtaDelta","delta_eta for of wrongly selected 2eta+gamma ",200,0,0.5);


    his1[i++] = hNCut;
    his1[i++] = hNGoodTagCut;
    his1[i++] = hTotalMass;
    his1[i++] = hNGamma;
    his1[i++] = hMGammaGamma;
    his1[i++] = hMGammaGammaCut;
    his1[i++] = hEGamma_Theta;
    his1[i++] = hEMCTime;
    his1[i++] = hPhotonEnergy;
    his1[i++] = hDecayAnglePair;
    his1[i++] = hDecayAngleCosPair;
    his1[i++] = hInitialTags;


    his1[i++] = hDeltaEta;
    his1[i++] = hDecayAngleCosEta;
    his1[i++] = hMEtaEta;
    his1[i++] = hEtaTags;

    his1[i++] = hEtaRadPhotonEnergy;
    his1[i++] = hBad2EtaRadPhotonEnergy;
    his1[i++] = hGood2EtaRadPhotonEnergy;

    his1[i++] = hDecayAngleEta;
    his1[i++] = hDecayAngleCosEta_Tag;
    his1[i++] = hEtaChisq;
    his1[i++] = hMRecoilEtaPair;
    his1[i++] = hEtaChisq_Tag;
    his1[i++] = hEvrEtaTags;
    his1[i++] = hGood2EtaChisq;
    his1[i++] = hBad2EtaChisq;

    his1[i++] = hGood2EtaChisq4C;
    his1[i++] = hBad2EtaChisq4C;

    his1[i++] = hGood2EtaChisq2C;
    his1[i++] = hBad2EtaChisq2C;

    his1[i++] = hGood2EtaDecayAngleCos;
    his1[i++] = hBad2EtaDecayAngleCos;
    his1[i++] = hGood2EtaDelta;
    his1[i++] = hBad2EtaDelta;

    his1[i++] = hGood2EtaEErr;
    his1[i++] = hBad2EtaEErr;

    his1[i++] = hGood2EtaNGamma;
    his1[i++] = hBad2EtaNGamma;
    his1[i++] = hGood2EtaNCombinations;
    his1[i++] = hBad2EtaNCombinations;
    his1[i++] = hGood2EtaTags;
    his1[i++] = hBad2EtaTags;

    his1[i++] = hGood2EtaMGammaPair;
    his1[i++] = hBad2EtaMGammaPair;
    his1[i++] = hGood2EtaMEtaEta;
    his1[i++] = hBad2EtaMEtaEta;
    his1[i++] = hMEtaGamma;

    his1[i++] = hGood2EtaMGammaGammaAll;
    his1[i++] = hBad2EtaMGammaGammaAll;
    his1[i++] = hBad2EtaChisq4C_NGamma;


    his1[i++] = hGood2EtaNPi0Candidates;
    his1[i++] = hBad2EtaNPi0Candidates ;
    his1[i++] = hGood2EtaNEtaCandidates;
    his1[i++] = hBad2EtaNEtaCandidates ;

    his1[i++] = hGood2EtaDeltaMggMpi0;
    his1[i++] = hBad2EtaDeltaMggMpi0;



    int j = 0;
    his2[j++] = hNPi0Candidates;
    his2[j++] = hDecayAngleCosPi0;
    his2[j++] = hMPi0Pi0;
    his2[j++] = hPi0RadPhotonEnergy;
    his2[j++] = hPi0Chisq;
    his2[j++] = hMPi0Gamma;
    his2[j++] = hPi0Tags;

    his2[j++] = hDeltaPi;
    his2[j++] = hDecayAnglePi0;
    his2[j++] = hMRecoilPi0Pair;
    his2[j++] = hMRecoilPi0Gamma;
    his2[j++] = hEvrPi0Tags;

    his2[j++] = hGood2Pi0Tags;
    his2[j++] = hBad2Pi0Tags;

    //~ his1[i++] = hMGammaGamma_chisq;
    selector->RegInDir(his1,"EtaEtagamma");


    selector->RegInDir(his2,"Pi0Pi0gamma");

    m_EventTagSvc = EventTagSvc::instance();
    m_EventTagSvc->setPdtFile(
        selector->AbsPath("Analysis/EventTag/share/pdt.table") );
    m_EventTagSvc->setDecayTabsFile(
        selector->AbsPath("Analysis/EventTag/share/decay.codes") );
    m_EventTagSvc->setUserDecayTabsFile(
        selector->AbsPath("BeanUser/mydecay.codes") );

    m_EventTagSvc->addUserChainTrig("f_4(2050)");
    m_EventTagSvc->addUserChainTrig("f_2");
    m_EventTagSvc->addUserChainTrig("a_20");
    m_EventTagSvc->addUserChainTrig("omega");
    m_EventTagSvc->addUserChainTrig("f_2(1950)");
    m_EventTagSvc->addUserChainTrig("eta'");
    m_EventTagSvc->addUserChainTrig("f'_2");
    m_EventTagSvc->addUserChainTrig("b_10");
    m_EventTagSvc->addUserChainTrig("eta");
    m_EventTagSvc->addUserChainTrig("phi");

    m_EventTagSvc->initialize();

}

int getFinalGammaCount(const TObjArray* arr) {
    int finals = 0;
    TObjArrayIter *it = (TObjArrayIter *) arr->MakeIterator();

    while (TMcParticle* particle = (TMcParticle *) it->Next()) {

        bool has_daughters = particle->getDaughters().size();

        if (!has_daughters) { //final state
            if ( particle->getParticleID()  == 22 ) {
                finals +=1;
            }
        }
    }

    delete it;
    return finals;
}

void dumpMcParticle(const TObjArray* particles,  int root, int shift) {
    TObjArrayIter *it = (TObjArrayIter *) particles->MakeIterator();
    while (TMcParticle* particle = (TMcParticle *) it->Next()) {
        if (particle->getMother() == root) { //daughter

            for (int i=0; i<shift; i++ ) {
                cout << "  ";
            }


            cout << m_EventTagSvc->pdg2name(particle->getParticleID())
                 //~ << "[" <<  particle->getParticleID()  << "]"

                << " " << particle->getInitialMomentumE()*1000 << " MeV"
                 << " p:(" << particle->getInitialMomentumX() << ","
                 << particle->getInitialMomentumY() << ","
                 << particle->getInitialMomentumZ() << ")"
            <<endl;

            dumpMcParticle(particles, particle->getTrackIndex(), shift+1);

        }

    }

    delete it;
}


void dumpMcParticles(const TObjArray* particles) {

    //~ TObjArrayIter *it = (TObjArrayIter *) particles->MakeIterator();
    //~ while (TMcParticle* particle = (TMcParticle *) it->Next()) {
        //~
        //~ cout << (TDatabasePDG::Instance())->GetParticle(particle->getParticleID())->GetTitle() <<","
             //~ << "index: " << particle->getTrackIndex() << ","
             //~ << "mother: " << particle->getMother() << endl;
        //~
    //~ }
    //~
    dumpMcParticle(particles, -99, 0);
}

static void formatTag(char * tag_str, unsigned int eventTag_) {
    unsigned int eventTag = eventTag_;

    unsigned char byte1 =   (eventTag & 0xFF); //last byte
    eventTag >>= 8;

    unsigned char byte2 =   (eventTag & 0xFF); //last byte
    eventTag >>= 8;

    unsigned char byte3 =   (eventTag & 0xFF); //last byte
    eventTag >>= 8;

    snprintf(tag_str, TAG_SRC_SIZE, "%d-%d-%d", byte1, byte2, byte3);

}

void     _handleMC(TMcEvent* m_TMcEvent,const TEvtRecEvent*  evtRecEvent,unsigned int eventTag) {
    TObjArrayIter *it = (TObjArrayIter *) m_TMcEvent->getMcParticleCol()->MakeIterator();

    //~ int motherIndex = -1;

    set<int> secondLevel;

    int finalGammaCount = 0;
    int finalOtherCount = 0;



    while (TMcParticle* particle = (TMcParticle *) it->Next()) {
        if (particle->getInitialMomentumE()*1000  > 1 ) {
            bool has_daughters = particle->getDaughters().size();
            if (!has_daughters) { //final state

                if ( particle->getParticleID()  == 22 ) {
                    finalGammaCount +=1;
                } else if (  abs(particle->getParticleID()) != 4 ) { // not c or c-bar quark
                    finalOtherCount += 1;
                }
            }


            if (particle->getMother() > -1) {
                secondLevel.insert(particle->getMother());
            }

        }


    }

    //~ if (finalOtherCount != 0) return;
    //~ if ( finalGammaCount != 5) return;

    if ( (finalGammaCount == 5) && (finalOtherCount == 0)) {

        int etaCount = 0;
        int pi0Count = 0;
        int secondLevelOtherCount = 0;


        for ( set<int>::iterator it=secondLevel.begin() ; it != secondLevel.end(); it++ ) {

            int pid = m_TMcEvent->getMcParticle(*it)->getParticleID()     ;

            if (pid == 111) { //pi0
                pi0Count+=1;
            } else if ( pid == 221) { //eta
                etaCount+=1;
            } else {
                secondLevelOtherCount+=1;
            }

        }


        if (((pi0Count==2) && (etaCount==0)) || ((etaCount==2) && (pi0Count==0))) {
            char tag_str[TAG_SRC_SIZE];
            unsigned int eventTagDecay =  (eventTag & 0xFFFFFF00) >> 8;

            formatTag(tag_str, eventTagDecay);

            cout << "(MC) !!! found 5 final gammas and 2 pi0/eta: " << tag_str << "(" << eventTagDecay<<")"  <<endl;
            dumpMcParticles(m_TMcEvent->getMcParticleCol());


            if ((etaCount==2) && (pi0Count==0)) {
                hEvrEtaTags->Fill(tag_str,1.0);
            }

            if ((pi0Count==2) && (etaCount==0))  {
                hEvrPi0Tags->Fill(tag_str,1.0);
            }

        }





          //~ cout <<   "pi0Count:" << pi0Count << endl;
    //~ cout <<   "etaCount:" << etaCount << endl;
    //~ cout <<   "secondLevelOtherCount:" << secondLevelOtherCount << endl;

    }

    //~ cout <<   "finalOtherCount:" << finalOtherCount << endl;
    //~ cout <<   "finalGammaCount:" << finalGammaCount << endl;


    delete it;


}



//constructs set, checks if there were intersections
bool intersection( set<int> & myset, int num, ... ) {
    va_list argptr;
    va_start( argptr, num );
    for( ; num > 0; num-- ) {
        int arg =  va_arg( argptr, int );
        if ( myset.insert(arg).second == false) return true;
    }

    va_end( argptr );

    return false;
}


typedef struct  {
    double chisq;
    bool kmfitOk;
    int radPhoton;
    pair<  pair<int,int>, pair<int,int>   > mesonPair;
    TLorentzVector pairMomentum;
    double totalE;
    int combinations;
    double chisq4C;
    double delta;
    double chisq2C;

    TLorentzVector radPhotonMomentum;
    pair< pair<TLorentzVector, TLorentzVector>, pair<TLorentzVector,TLorentzVector> > mesonPairMomenta;



} Combination;

void HepLorentzToRootLorentz(HepLorentzVector from, TLorentzVector& to) {
   to.SetXYZT( from.x(), from.y(), from.z(), from.t());
}
Combination findCombinationLIU(const Vtracks& photons,const vector<TLorentzVector>& photonsP, double mesonMass, double mesonWindow, double vetoMass, double vetoWindow) { // using C.Y. LIU cuts from PWA of J/\psi -> \gamma \eta \eta
    // bool pairFound;
    // int fits;
    Combination result;
    result.chisq   = -1;
    result.kmfitOk = false;
    result.combinations = 0;
    result.chisq4C = 9999.0;
    result.chisq2C = -100;
    result.delta = 9999.0;

    KinematicFit * kmfit = KinematicFit::instance();
    HepLorentzVector ecms(0.034,0,0,3.097); // j/psi
    int nGamma = photons.size();

    // if only 5 photons presented do fake 'for' loop from bad=-1 to bad=0 (1 value). All 5 tracks will be added to kmfit

    int forFrom, forTo;
    if (nGamma == 5) {
        forFrom = photons.size();
        forTo = photons.size()+1;
    } else {
        forFrom = 0;
        forTo = photons.size();
    }
    bool kmfit4COk = false;
    int badIndex = -1;


    for (int bad = forFrom; bad < forTo; bad ++) { //iterating over bad photon, fiting the rest 5 photons
        kmfit->init();

        // the rest 5 photons
        for (int i=0; i < bad; i++) {
            kmfit->AddTrack(i, 0.0, photons[i]->emcShower());
        }

        for (int i = bad + 1; i < int(photons.size()); i++) {
            kmfit->AddTrack(i-1, 0.0, photons[i]->emcShower());
        }

        kmfit->AddFourMomentum(0, ecms);

        bool oksq = kmfit->Fit();
        if (oksq) {
            double chi2 = kmfit->chisq();
            //~ result.combinations +=1;

            if(chi2 < result.chisq4C) {
                result.chisq4C = chi2;


                if (chi2 < 50)
                    kmfit4COk = true;

                badIndex = bad;

            }
        }
    }

    if (kmfit4COk) {
        Vtracks new_photons = photons;

        vector<TLorentzVector> new_photonsP = photonsP;

        if ( badIndex < int(new_photons.size()) ) { //there were 6 photons
            new_photons.erase(new_photons.begin() + badIndex);  //deletes bad photon from vector
            new_photonsP.erase(new_photonsP.begin() + badIndex);
        }



        for (int first = 0; first < int(new_photons.size()); first ++) {
            for (int second = first + 1; second < int(new_photons.size()); second ++) {
                double mass= (new_photonsP[first] + new_photonsP[second]).M();

                if (TMath::Abs(mass - vetoMass) < vetoWindow) {
                    return result; //veto failed
                }
            }
        }





        for (int pair1first=0; pair1first < int(photons.size()) - 1; pair1first++) {
            for (int pair1second=pair1first+1; pair1second < int(photons.size()); pair1second++) {  //select first two photons
                for (int rest = 0; rest < int(photons.size()) - 1; rest++) {
                    if ( (rest != pair1first) && (rest != pair1second)) { // rest(non-eta) photon selected
                        int pair2first = -1;
                        int pair2second = -1;


                        for (int i = 0; i < int(photons.size()); i ++)  {
                            if ( (i != pair1first) && (i != pair1second) && (i != rest) ) {
                                if (pair2first == -1) {
                                    pair2first = i;
                                } else {
                                    pair2second = i;
                                }
                            }
                        }

                        double MPair1 =   (photonsP[pair1first]+photonsP[pair1second]).M(); // mass
                        double MPair2 =   (photonsP[pair2first]+photonsP[pair2second]).M();

                        double delta = TMath::Sqrt(  ( MPair1 - mesonMass) *  ( MPair1 - mesonMass) +
                                                     ( MPair2 - mesonMass) *  ( MPair2 - mesonMass) );



                        //final cuts vvv
                        if ( TMath::Abs(MPair1 - mesonMass) > mesonWindow ) continue;
                        if ( TMath::Abs(MPair2 - mesonMass) > mesonWindow ) continue;



                        if (delta < result.delta) {
                            result.mesonPair.first.first = pair1first;
                            result.mesonPair.first.second = pair1second;

                            result.mesonPair.second.first = pair2first;
                            result.mesonPair.second.second = pair2second;

                            result.radPhoton = rest;
                            result.delta = delta;
                            result.kmfitOk = true;

                            if (delta < 0.05) {
                                result.combinations += 1;
                            }

                        }
                    }
                }
            }
        }





    }

    if (result.kmfitOk) {
        //~ if (goodEventTag) hNGoodTagCut->Fill(6);


        result.totalE = photons[result.mesonPair.first.first]->emcShower()->energy() +
                        photons[result.mesonPair.first.second]->emcShower()->energy() +
                        photons[result.mesonPair.second.first]->emcShower()->energy() +
                        photons[result.mesonPair.second.second]->emcShower()->energy() +
                        photons[result.radPhoton]->emcShower()->energy() ;


        result.pairMomentum =
            photonsP[result.mesonPair.first.first] +
            photonsP[result.mesonPair.first.second] +
            photonsP[result.mesonPair.second.first] +
            photonsP[result.mesonPair.second.second];

    }

    if (result.combinations > 1) {
        result.kmfitOk = 0;
    }

    return result;


}


Combination findCombination(const Vtracks& photons,const vector<TLorentzVector>& photonsP, const vector<int>& radPhotonCandidates, const vector<pair<int,int> >& mesonCandidates, double mesonMass,const vector<pair<int,int> >& vetoPairs) {
    // bool pairFound;
    int fits;
    Combination result;
    result.chisq   = 9999.;
    result.kmfitOk = false;
    result.combinations = 0;
    result.chisq4C = -100;
    result.chisq2C = -100;


    KinematicFit * kmfit = KinematicFit::instance();
    HepLorentzVector ecms(0.034,0,0,3.097); // j/psi
    //~ HepLorentzVector ecms(0.0405,0,0,3.686); // psi'


    for (int radIndex=0; radIndex < int(radPhotonCandidates.size()); radIndex++) {
        for (int first=0; first < int(mesonCandidates.size()); first ++ ){
            for (int second = first + 1; second < int(mesonCandidates.size()); second ++ ){
                // no intersections
                set<int> photons_set;
                if (intersection(photons_set, 5,
                                     radPhotonCandidates[radIndex],
                                     mesonCandidates[first].first,
                                     mesonCandidates[first].second,
                                     mesonCandidates[second].first,
                                     mesonCandidates[second].second ))   break ;

                //~ cout << "veto:";
                bool veto = false;
                for (int veto_index = 0; veto_index < int(vetoPairs.size()); veto_index ++){
                    if (photons_set.count(vetoPairs[veto_index].first) && photons_set.count(vetoPairs[veto_index].second)) {
                        veto = true;
                        break;
                    }


                    //~ cout << "count[" << veto_index << "]:"
                         //~ << photons_set.count(vetoPairs[veto_index].first)
                         //~ << ","
                         //~ << photons_set.count(vetoPairs[veto_index].second) << endl;

                    //~ cout << vetoPairs[veto_index].first << "," << vetoPairs[veto_index].second << ",";
                 }
                //~ cout << endl;

                // and there is not veto'ed pair among this 5 photons
                if (veto) {
                    break;
                }

                //~ if (!pairFound) {
                    //~ pairFound = true;
                    //~ if (goodEventTag) hNGoodTagCut->Fill(4);
                //~ }

                // do 6C fit
                fits+=1;
                kmfit->init();

                kmfit->AddTrack(0, 0.0, photons[radPhotonCandidates[radIndex]]->emcShower() );
                kmfit->AddTrack(1, 0.0, photons[mesonCandidates[first].first ]->emcShower());
                kmfit->AddTrack(2, 0.0, photons[mesonCandidates[first].second]->emcShower() );
                kmfit->AddTrack(3, 0.0, photons[mesonCandidates[second].first]->emcShower() );
                kmfit->AddTrack(4, 0.0, photons[mesonCandidates[second].second ]->emcShower());

                kmfit->AddFourMomentum(0, ecms);
                kmfit->AddResonance(1, mesonMass, 1, 2);
                kmfit->AddResonance(2, mesonMass, 3, 4);

                bool oksq = kmfit->Fit();

                if(oksq) {
                    double chi2 = kmfit->chisq();
                    result.combinations +=1;


                    if(chi2 < result.chisq) {
                        result.chisq = chi2;

                        if (chi2 < 70)  {
                            result.kmfitOk = true;
                        }


                        result.mesonPair.first = mesonCandidates[first];
                        result.mesonPair.second = mesonCandidates[second];
                        result.radPhoton = radPhotonCandidates[radIndex];

                        kmfit->init();

                        kmfit->AddTrack(0, 0.0, photons[radPhotonCandidates[radIndex]]->emcShower() );
                        kmfit->AddTrack(1, 0.0, photons[mesonCandidates[first].first ]->emcShower());
                        kmfit->AddTrack(2, 0.0, photons[mesonCandidates[first].second]->emcShower() );
                        kmfit->AddTrack(3, 0.0, photons[mesonCandidates[second].first]->emcShower() );
                        kmfit->AddTrack(4, 0.0, photons[mesonCandidates[second].second ]->emcShower());

                        kmfit->AddFourMomentum(0, ecms);
                        if (kmfit->Fit()) {
                            result.chisq4C = kmfit->chisq();
                        }

                        kmfit->init();
                        kmfit->AddTrack(0, 0.0, photons[mesonCandidates[first].first ]->emcShower());
                        kmfit->AddTrack(1, 0.0, photons[mesonCandidates[first].second]->emcShower() );
                        kmfit->AddTrack(2, 0.0, photons[mesonCandidates[second].first]->emcShower() );
                        kmfit->AddTrack(3, 0.0, photons[mesonCandidates[second].second ]->emcShower());
                        kmfit->AddResonance(0, mesonMass, 0, 1);
                        kmfit->AddResonance(1, mesonMass, 2, 3);

                        if (kmfit->Fit()) {
                            result.chisq2C = kmfit->chisq();
                        }
                    }
                }
            }
        }
    }


    //~ cout << "fits: " << fits << endl;
    if (result.kmfitOk) {
        //~ if (goodEventTag) hNGoodTagCut->Fill(6);
        HepLorentzToRootLorentz(kmfit->pfit(0), result.radPhotonMomentum);
        HepLorentzToRootLorentz(kmfit->pfit(1), result.mesonPairMomenta.first.first);
        HepLorentzToRootLorentz(kmfit->pfit(2), result.mesonPairMomenta.first.second);
        HepLorentzToRootLorentz(kmfit->pfit(3), result.mesonPairMomenta.second.first);
        HepLorentzToRootLorentz(kmfit->pfit(4), result.mesonPairMomenta.second.second);

        result.totalE = photons[result.mesonPair.first.first]->emcShower()->energy() +
                        photons[result.mesonPair.first.second]->emcShower()->energy() +
                        photons[result.mesonPair.second.first]->emcShower()->energy() +
                        photons[result.mesonPair.second.second]->emcShower()->energy() +
                        photons[result.radPhoton]->emcShower()->energy() ;


        result.pairMomentum =
            photonsP[result.mesonPair.first.first] +
            photonsP[result.mesonPair.first.second] +
            photonsP[result.mesonPair.second.first] +
            photonsP[result.mesonPair.second.second];


    }
    return result;
}


    //~ cout  << "TAG:"  << hex << eventTag << dec << endl;

//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
bool EtaetagammaEvent(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//-----------------------------------------------------------------------------
{
    if( selector->Verbose() ) cout << " RhopiEvent() " << endl;
    m_EventTagSvc->setMcEvent(m_TMcEvent);

    bool mc = 1;


    m_abscor->AbsorptionCorrection(selector);

    int event = m_TEvtHeader->getEventId();
    int runNo = m_TEvtHeader->getRunId();

    unsigned int eventTag = 0;
    if (mc) eventTag = m_EventTagSvc->getEventTag() ;

    if( selector->Verbose() )
        cout << "run, evtnum = " << runNo << " , " << event <<endl;





    // for MC
    unsigned int eventTagDecay =  (eventTag & 0xFFFFFF00) >> 8;


    char tag_str[TAG_SRC_SIZE];
    formatTag(tag_str, eventTagDecay);
    hInitialTags->Fill(tag_str,1.0);

    // GOOD eta decays:
    bool mc_good_2eta = 0;

    if ( (eventTagDecay == 459076) || //68-1-7
         (eventTagDecay == 67355) || //27-7-1   f_2(1950)   x
         (eventTagDecay == 67086) || //14-6-1   f_4(2050)    x
         (eventTagDecay == 66586) || // 26-4-1   f'_2         (1.525)          x
         (eventTagDecay == 65837) || //45-1-1  x
         (eventTagDecay == 67856) || //16-9-1  f_2(1275) x
         (eventTagDecay == 262467) || //67-1-4  eta omega; omega -> eta gamma
         (eventTagDecay == 65540) ) {// 4-0-1 gamma eta_c; eta_c -> 2eta
            mc_good_2eta = 1;
        }

    bool mc_good_2pi0 = 0;

    if ( (eventTagDecay == 526) ||  //14-2-0  f_4(2050)
       (eventTagDecay == 528) ||  //16-2-0 f_2
       (eventTagDecay == 42) ||  // 42-0-0 phsp
       (eventTagDecay == 581) ||  // 69-2-0  pi0-omega
       (eventTagDecay == 647) ||  //135-2-0 b_10-pi0
       (eventTagDecay == 2091) ||  //43-8-0 pi0+a_20
       (eventTagDecay == 1562) ){  //26-6-0 f'_2
            mc_good_2pi0 = 1;
        }


    // end  mc

    const TEvtRecEvent*  evtRecEvent = m_TEvtRecObject->getEvtRecEvent();



    if (evtRecEvent->totalCharged() > 0) {
      return false;
    }

    if (mc)  _handleMC(m_TMcEvent,evtRecEvent,eventTag);



    if( selector->Verbose() ) {
        cout <<"ncharg, nneu, tottks = "
             << evtRecEvent->totalCharged() << " , "
             << evtRecEvent->totalNeutral() << " , "
             << evtRecEvent->totalTracks() << endl;
    }

    // evtRecTrkCol
    const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

    double eTotal = 0.0;
    int nGamma = 0;

    Vtracks photons;
    photons.clear();
    vector<int> photonsIndicies;

    for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
        DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
        if(!itTrk->isEmcShowerValid())
            continue;

        RecEmcShower *emcTrk = itTrk->emcShower();

        //~ cout << "emctrk time=" << emcTrk->time() << endl;
        hEMCTime->Fill(emcTrk->time());
        if ( (emcTrk->time() < 0) || (emcTrk->time() >= 14))
            continue;

        double eraw = emcTrk->energy();

        if (selector->Verbose()) {
            cout << "gamma e=" << eraw << ", "
                 << "theta=" << emcTrk->theta() << ","
                 << "phi=" << emcTrk->phi() << ","
                 << "cell=" << emcTrk->cellId()
            << endl;
        }


        // bool good = 0;
        double absCosTheta = TMath::Abs(  TMath::Cos(emcTrk->theta()) );





        if  (!(   ( ( absCosTheta < 0.8 ) &&  ( eraw > 25E-3) ) ||
                  ( ( absCosTheta > 0.86 ) && ( absCosTheta < 0.92 ) && ( eraw > 50E-3)  )   )){
            continue;
        }


        hEGamma_Theta->Fill( emcTrk->theta(), emcTrk->energy());
        eTotal += eraw;
        nGamma+=1;
        photons.push_back(itTrk);
        photonsIndicies.push_back(i);
    }



    hNGamma->Fill(nGamma);

    if (nGamma < 5) return false;


    //~ return true;


    if (nGamma > 8) return false;
    //~ if (nGamma > 6) return false;
    if (nGamma > 5) return false;

    hTotalMass->Fill(eTotal);

    if (eTotal < 2.5)  return false;



    // Assign four-momentum to each photon selected
    vector<TLorentzVector> photonsP;
    photonsP.resize(photons.size());
    for (int i=0; i < int(photons.size()); i++) {
        RecEmcShower* shower = photons[i]->emcShower();
        double eraw  = shower->energy();
        double phi   = shower->phi();
        double the   = shower->theta();

        photonsP[i].SetPxPyPzE( eraw*sin(the)*cos(phi),
                                eraw*sin(the)*sin(phi),
                                eraw*cos(the),
                                eraw                 );



        // fill photon energy histogram
        hPhotonEnergy->Fill(eraw);
    }


    // Select pi0 and eta candidates


    vector<pair<int,int> > pi0Candidates, pi0VetoPairs;
    vector<pair<int,int> > etaCandidates;

    vector<int> radPhotonCandidates;

    for (int first=0; first < int(photons.size()); first++) {
        for (int second=first+1; second < int(photons.size()); second++) {
            TLorentzVector pair4Momentum = photonsP[first]+photonsP[second];
            double MGammaGamma =   pair4Momentum.M(); //invariant mass
            hMGammaGamma->Fill( MGammaGamma);


            // pi0 mass window from B.J LIU et al,"J/ψ and ψ(2S) Decays into γπ0π0, γηη, γηη0, and γη0η0 Final States"
            // LIU_collab.pdf

            pair<int, int> thisPair(first,second);

            TVector3 boost = -pair4Momentum.BoostVector();
            TLorentzVector gamma1Momentum = photonsP[first];
            gamma1Momentum.Boost(boost);

            double decayAngle = gamma1Momentum.Angle( pair4Momentum.Vect() );
            double decayAngleCos = TMath::Abs(TMath::Cos(decayAngle));


            hDecayAngleCosPair->Fill(decayAngleCos);
            hDecayAnglePair->Fill(decayAngle);


            if ( ( MGammaGamma > mpi0 - 0.06) &&  ( MGammaGamma < mpi0 + 0.04)) { //pi0 mass window
                if (decayAngleCos < 0.95)
                {
                    pi0Candidates.push_back(thisPair);
                    hMGammaGammaCut->Fill(MGammaGamma);
                }

                hDecayAngleCosPi0->Fill(decayAngleCos);
                hDecayAnglePi0->Fill(decayAngle);
            }


            //~ if ( ( MGammaGamma > mpi0 - 0.025) &&  ( MGammaGamma < mpi0 + 0.02)) { //pi0 veto mass window
            if ( ( MGammaGamma > mpi0 - 0.015) &&  ( MGammaGamma < mpi0 + 0.015)) { //pi0 veto mass window
                pi0VetoPairs.push_back(thisPair);
            }


            if ( (MGammaGamma > meta - 0.09) && ( MGammaGamma < meta + 0.06)) {
            //~ if ( (MGammaGamma > meta - 0.04) && ( MGammaGamma < meta + 0.04)) {
                if (decayAngleCos < 0.95)
                {
                    etaCandidates.push_back(thisPair);
                    hMGammaGammaCut->Fill(MGammaGamma);
                }
                hDecayAngleCosEta->Fill(decayAngleCos);
                hDecayAngleEta->Fill(decayAngle);
            }


            // decay angle


        }
    }

    for (int i=0; i < int(photons.size()); i++ ){
        //~ if (  photonsP[i].E()  > 0.3 ) {
            radPhotonCandidates.push_back(i);
         //~ }
    }


    //loop over all radiative Photons, and all pi0/eta pairs




    pair<  pair<int,int>, pair<int,int>   > pi0pair;




    // int fits = 0;


    // bool onePairFound = false;
    hNPi0Candidates->Fill(pi0Candidates.size());

    // fill histogram with tags
    formatTag(tag_str, eventTagDecay);

    TLorentzVector Pcms(0.034,0,0,3.097); //


    bool pi0_or_eta_found = 0;

    vector<pair<int,int> > nullvector;
    //~ Combination pi0result = findCombination(photons, photonsP, radPhotonCandidates, pi0Candidates, mpi0) ;
    Combination pi0result = findCombination(photons, photonsP, radPhotonCandidates, pi0Candidates, mpi0, nullvector) ; //MINE

    //~ Combination pi0result = findCombinationLIU(photons, photonsP, mpi0, 10E-3, meta, 1E-3) ;

    if (pi0result.kmfitOk) {
        cout << "Ngamma:" << nGamma << endl;
        hPi0Chisq->Fill(pi0result.chisq);
        //~ if (pi0result.chisq > 40) return false;

        // gamma-pi0 invariant mass
        double mPi0Gamma1 = (photonsP[pi0result.mesonPair.first.first] + photonsP[pi0result.mesonPair.first.second] +photonsP[pi0result.radPhoton]).M();
        double mPi0Gamma2 = (photonsP[pi0result.mesonPair.second.first] + photonsP[pi0result.mesonPair.second.second] +photonsP[pi0result.radPhoton]).M();

        hMPi0Gamma->Fill(mPi0Gamma1);
        hMPi0Gamma->Fill(mPi0Gamma2);

        if ( TMath::Abs(mPi0Gamma1 - momega)  < 10E-3) return false;
        if ( TMath::Abs(mPi0Gamma2 - momega)  < 10E-3) return false;


        hMPi0Pi0->Fill(pi0result.pairMomentum.M());
        hPi0RadPhotonEnergy->Fill(photonsP[pi0result.radPhoton].E());




        hMRecoilPi0Pair->Fill( (Pcms - pi0result.pairMomentum).M());

        hMRecoilPi0Gamma->Fill(( Pcms - (photonsP[pi0result.mesonPair.first.first] + photonsP[pi0result.mesonPair.first.second] +photonsP[pi0result.radPhoton])).M());
        hMRecoilPi0Gamma->Fill(( Pcms - (photonsP[pi0result.mesonPair.second.first] + photonsP[pi0result.mesonPair.second.second] +photonsP[pi0result.radPhoton])).M());

        hPi0Tags->Fill(tag_str,1.0);
        cout << "tag: " << tag_str << endl;
        if (mc) {
            cout << "pi0 selection passed. dump:" << endl;
            dumpMcParticles(m_TMcEvent->getMcParticleCol());


            if (mc_good_2pi0) {
                hGood2Pi0Tags->Fill(tag_str,1.0);
            } else {
                hBad2Pi0Tags->Fill(tag_str,1.0);
            }


        }

        pi0_or_eta_found = 1;
    }

    Combination etaresult = findCombination(photons, photonsP, radPhotonCandidates, etaCandidates, meta, pi0VetoPairs) ; //MINE
    //~ Combination etaresult = findCombinationLIU(photons, photonsP, meta, 40E-3, mpi0, 20E-3) ; // LIU
    //~ Combination etaresult = findCombinationLIU(photons, photonsP, meta, 40E-3, mpi0, 0E-3) ; // LIU


    if (etaresult.kmfitOk) {
    //~ if (etaresult.chisq > 40) return false;

        // DEBUG vvvv
        TLorentzVector pair4Momentum = photonsP[etaresult.mesonPair.first.first]   +   photonsP[etaresult.mesonPair.first.second]   ;
        TVector3 boost = -pair4Momentum.BoostVector();
        TLorentzVector gamma1Momentum = photonsP[etaresult.mesonPair.first.first];
        gamma1Momentum.Boost(boost);
        double decayAngle1 = gamma1Momentum.Angle( pair4Momentum.Vect() );
        double decayAngleCos1 = TMath::Abs(TMath::Cos(decayAngle1));
        hDecayAngleCosEta_Tag->Fill(decayAngle1, tag_str,1.0);


        pair4Momentum = photonsP[etaresult.mesonPair.second.first]   +   photonsP[etaresult.mesonPair.second.second]   ;
        boost = -pair4Momentum.BoostVector();
        gamma1Momentum = photonsP[etaresult.mesonPair.second.first];
        gamma1Momentum.Boost(boost);

        double decayAngle2 = gamma1Momentum.Angle( pair4Momentum.Vect() );
        double decayAngleCos2 = TMath::Abs(TMath::Cos(decayAngle2));
        hDecayAngleCosEta_Tag->Fill(decayAngleCos2, tag_str,1.0);





        double mEtaGamma1 = (photonsP[etaresult.mesonPair.first.first] + photonsP[etaresult.mesonPair.first.second] +photonsP[etaresult.radPhoton]).M();
        double mEtaGamma2 = (photonsP[etaresult.mesonPair.second.first] + photonsP[etaresult.mesonPair.second.second] +photonsP[etaresult.radPhoton]).M();

        hMEtaGamma->Fill(mEtaGamma1);
        hMEtaGamma->Fill(mEtaGamma2);
        //~
        //~ if ( TMath::Abs(mEtaGamma1 - mphi)  < 30E-3) return false;
        //~ if ( TMath::Abs(mEtaGamma2 - mphi)  < 30E-3) return false;



        hEtaChisq->Fill(etaresult.chisq);

        cout << "energy: " << etaresult.totalE << endl;
        cout << "chisq: " << etaresult.chisq << endl;


        hMEtaEta->Fill(etaresult.pairMomentum.M());
        hEtaRadPhotonEnergy->Fill(photonsP[etaresult.radPhoton].E());


        if (mc) {
            hEtaTags->Fill(tag_str,1.0);
            cout << "tag: " << tag_str << endl;

            hEtaChisq_Tag->Fill(etaresult.chisq, tag_str, 1.0);

            double MPair1 = (photonsP[etaresult.mesonPair.first.first]   +   photonsP[etaresult.mesonPair.first.second]).M()   ;
            double MPair2 = (photonsP[etaresult.mesonPair.second.first]   +   photonsP[etaresult.mesonPair.second.second]).M()   ;

            double delta_eta = TMath::Sqrt(  ( MPair1 - meta) *  ( MPair1 - meta) +
                                                    ( MPair2 - meta) *  ( MPair2 - meta) );




            vector<int> selectedPhotons;
            selectedPhotons.push_back(etaresult.mesonPair.first.first);
            selectedPhotons.push_back(etaresult.mesonPair.second.first);
            selectedPhotons.push_back(etaresult.mesonPair.first.second);
            selectedPhotons.push_back(etaresult.mesonPair.second.second);
            selectedPhotons.push_back(etaresult.radPhoton);




            if (mc_good_2eta) {
                hGood2EtaRadPhotonEnergy->Fill(photonsP[etaresult.radPhoton].E());
                hGood2EtaChisq->Fill(etaresult.chisq);
                hGood2EtaChisq4C->Fill(etaresult.chisq4C);
                hGood2EtaChisq2C->Fill(etaresult.chisq2C);

                hGood2EtaDecayAngleCos->Fill(decayAngleCos1);
                hGood2EtaDecayAngleCos->Fill(decayAngleCos2);
                hGood2EtaDelta->Fill(delta_eta);

                hGood2EtaEErr->Fill(photons[etaresult.mesonPair.first.first]->emcShower()->dE());
                hGood2EtaEErr->Fill(photons[etaresult.mesonPair.second.first]->emcShower()->dE());
                hGood2EtaEErr->Fill(photons[etaresult.mesonPair.first.second]->emcShower()->dE());
                hGood2EtaEErr->Fill(photons[etaresult.mesonPair.second.second]->emcShower()->dE());

                hGood2EtaNGamma->Fill(nGamma);

                hGood2EtaNCombinations->Fill(etaresult.combinations);

                hGood2EtaTags->Fill(tag_str,1.0);

                hGood2EtaMGammaPair->Fill( MPair1);
                hGood2EtaMGammaPair->Fill( MPair2);

                hGood2EtaMEtaEta->Fill(etaresult.pairMomentum.M());


                double delta_pi0_min_abs = 9999.0;
                double delta_pi0_min = 9999.0;

                for (int first = 0; first < int(selectedPhotons.size()); first ++)
                    for (int second = first + 1; second < int(selectedPhotons.size()); second ++) {
                        double mass =  (photonsP[selectedPhotons[first]] + photonsP[selectedPhotons[second]]).M();
                        hGood2EtaMGammaGammaAll->Fill( mass);

                        if ( TMath::Abs(mass - mpi0) < delta_pi0_min_abs) {
                            delta_pi0_min_abs = TMath::Abs(mass - mpi0);
                            delta_pi0_min = mass - mpi0;
                        }

                    }
                hGood2EtaDeltaMggMpi0->Fill(delta_pi0_min);



                hGood2EtaNPi0Candidates->Fill(pi0Candidates.size() );
                hGood2EtaNEtaCandidates->Fill(etaCandidates.size() );

            } else {
                hBad2EtaRadPhotonEnergy->Fill(photonsP[etaresult.radPhoton].E());
                hBad2EtaChisq4C->Fill(etaresult.chisq4C);
                hBad2EtaChisq2C->Fill(etaresult.chisq2C);
                hBad2EtaChisq->Fill(etaresult.chisq);
                hBad2EtaDecayAngleCos->Fill(decayAngleCos1);
                hBad2EtaDecayAngleCos->Fill(decayAngleCos2);
                hBad2EtaDelta->Fill(delta_eta);

                hBad2EtaEErr->Fill(photons[etaresult.mesonPair.first.first]->emcShower()->dE());
                hBad2EtaEErr->Fill(photons[etaresult.mesonPair.second.first]->emcShower()->dE());
                hBad2EtaEErr->Fill(photons[etaresult.mesonPair.first.second]->emcShower()->dE());
                hBad2EtaEErr->Fill(photons[etaresult.mesonPair.second.second]->emcShower()->dE());

                hBad2EtaNGamma->Fill(nGamma);
                hBad2EtaNCombinations->Fill(etaresult.combinations);
                hBad2EtaTags->Fill(tag_str,1.0);

                hBad2EtaMGammaPair->Fill( MPair1);
                hBad2EtaMGammaPair->Fill( MPair2);

                hBad2EtaMEtaEta->Fill(etaresult.pairMomentum.M());

                double delta_pi0_min_abs = 9999.0;
                double delta_pi0_min = 9999.0;

                for (int first = 0; first < int(selectedPhotons.size()); first ++)
                    for (int second = first + 1; second < int(selectedPhotons.size()); second ++) {
                        double mass =  (photonsP[selectedPhotons[first]] + photonsP[selectedPhotons[second]]).M();
                        hBad2EtaMGammaGammaAll->Fill( mass);

                        if ( TMath::Abs(mass - mpi0) < delta_pi0_min_abs) {
                            delta_pi0_min_abs = TMath::Abs(mass - mpi0);
                            delta_pi0_min = mass - mpi0;
                        }

                    }
                hBad2EtaDeltaMggMpi0->Fill(delta_pi0_min);

                hBad2EtaNPi0Candidates->Fill(pi0Candidates.size() );
                hBad2EtaNEtaCandidates->Fill(etaCandidates.size() );

                //~ hBad2EtaChisq4C_NGamma->Fill(nGamma*1.0,etaresult.chisq4C );
                //~ cout << "hBad2EtaChisq4C_NGamma->Fill(" << etaresult.chisq4C << "," << nGamma  <<")" << endl;
            }

            cout << "eta selection passed. dump:" << endl;
            dumpMcParticles(m_TMcEvent->getMcParticleCol());
        }


        hMRecoilEtaPair->Fill( (Pcms - etaresult.pairMomentum).M());





        pi0_or_eta_found = 1;
    }


    if (pi0_or_eta_found) {
        return true;
    }

    return false;

}

//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
void EtaetagammaEndJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
  if( selector->Verbose() ) cout << " RhopiEndJob() " << endl;

  //~ cout << "in finalize()" << endl;
  //~ cout << "in finalize()" << endl;
  hPi0Tags->GetXaxis()->LabelsOption(">");
  hEtaTags->GetXaxis()->LabelsOption(">");
}

#ifdef __cplusplus
}
#endif
