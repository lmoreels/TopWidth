#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <array>
#include "TRandom3.h"
#include "TNtuple.h"
#include <TFile.h>
#include <TLeaf.h>

// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
//#include "../macros/Style.C"

// user defined
#include "Tools/interface/ResolutionFunctions.h"
#include "Tools/interface/KinFitter.h"


using namespace std;
using namespace TopTree;


bool test = false;
bool testHistos = false;
bool testTTbarOnly = false;
bool makePlots = true;
bool calculateResolutionFunctions = false;
bool calculateAverageMass = false;
bool calculateLikelihood = true;
bool doKinFit = true;
bool useToys = false;
bool applyLeptonSF = true;
bool applyPU = true;
bool applyJER = true;
bool applyJEC = true;
bool applyBTagSF = true;
bool applyNloSF = false;

bool applyWidthSF = false;
double scaleWidth = 0.66;

bool useOldNtuples = true;
string systStr = "nominal";
string whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos )
  {
    if (useOldNtuples) return "160812";
    else return "170216";
  }
  else if ( syst.find("JERup") != std::string::npos ) return "160916";
  else if ( syst.find("JERdown") != std::string::npos ) return "160930";
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    return "160812";
  }
}
string ntupleDate = whichDate(systStr);
int verbose = 2;

string pathNtuples = "";
bool isData = false;
bool isTTbar = false;

int nofHardSelected = 0;
int nofMatchedEvents = 0;
int nofHadrMatchedEvents = 0;
int nofCorrectlyMatched = 0;
int nofNotCorrectlyMatched = 0;
int nofCP = 0, nofWP = 0, nofUP = 0;
int nofCP_TT = 0, nofWP_TT = 0, nofUP_TT = 0;
double Luminosity = 9999.;



///  Working points for b tagging  // Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
double CSVv2Loose =  0.460;
double CSVv2Medium = 0.800;
double CSVv2Tight = 0.935;

/// Top width
double genTopWidth = 1.366; // from fit
double genTopMass = 172.3; // from fit

// Temporarily, until calculated from TTbar sample
double chi2WMass = 80.385;
double sigmaChi2WMass = 10;
double chi2TopMass = 172.5; //180.0; //from mtop mass plot: 167.0
double sigmaChi2TopMass = 40;

/// Likelihood function
const int nCP = 356970;  //357134;
const int nWP = 127011;  //155459;
const int nUP = 259755;  //297689;

// Voigt/CB fit parameters
const double mu_CP = 1.01, sigma_CP = 0.0664, r_CP = 2.0298, norm_CP = 0.002557;
const double alpha_WP = -0.38, n_WP = 11.668, sigma_WP = 0.1634, mu_WP = 0.889, norm_WP = 0.003655;
const double alpha_UP = -0.967, n_UP = 2.45, sigma_UP = 0.1897, mu_UP = 0.97, norm_UP = 0.004305;
const double alpha_WPUP = -0.8207, n_WPUP = 3.05, sigma_WPUP = 0.2, mu_WPUP = 0.97, norm_WPUP = 0.003967;
const double norm_comb = 1.; //0.956039;
const double gammaConvConst = 0.0312447, gammaConvRico = 0.00792538;
//const double gammaConvConst = 0., gammaConvRico = 1.;

/// Average top mass
// TT gen match, TT reco match, TT reco wrongMatch WP/UP, TT reco noMatch, TT reco wrongPerm, TT reco wrongPerm W Ok, TT reco wrongPerm W Not Ok, TT reco, ST_t_top reco, ST_t_antitop reco, ST_tW_top reco, ST_tW_antitop reco, DYJets reco, WJets reco, data reco, all MC reco, all samples reco (data+MC) 
// also background in CP/WP/UP cats (unlike name suggests)
const int nofAveMasses = 17;
std::array<double, nofAveMasses> aveTopMass;
std::array<double, nofAveMasses> aveTopMass_noWidth    = {168.719, 167.105, 203.378, 204.724, 197.450, 201.182, 185.360, 193.207, 270.895, 267.167, 230.144, 229.649, 250.010, 242.091, 200.455, 193.762, 193.825};
std::array<double, nofAveMasses> aveTopMass_KFchi2cut5 = {168.719, 167.253, 192.093, 189.672, 196.716, 199.756, 165.839, 180.817, 249.629, 249.039, 227.992, 224.213, 221.995, 213.278, 184.884, 181.158, 181.191};
std::array<double, nofAveMasses> aveTopMass_KFchi2cut6 = {168.719, 167.206, 192.030, 189.733, 196.637, 199.866, 165.879, 180.978, 249.326, 251.815, 229.616, 224.083, 219.695, 208.480, 184.976, 181.327, 181.359};
std::array<double, nofAveMasses> aveTopMass_KFchi2cut7 = {168.719, 167.156, 191.966, 189.826, 196.450, 199.847, 166.001, 181.128, 248.543, 251.538, 228.486, 224.230, 221.079, 211.266, 185.483, 181.482, 181.517};
std::array<double, nofAveMasses> aveTopMass_KFchi2cut8 = {168.719, 167.111, 191.928, 189.896, 196.361, 199.930, 166.184, 181.282, 247.694, 250.639, 227.952, 224.751, 221.542, 210.461, 185.774, 181.641, 181.677};
std::array<double, nofAveMasses> aveTopMass_KFchi2cut9 = {168.719, 167.069, 191.902, 189.970, 196.277, 199.995, 166.503, 181.428, 248.440, 251.997, 227.283, 224.807, 223.771, 202.844, 185.814, 181.795, 181.830};
std::array<double, nofAveMasses> aveTopMass_KFchi2cut10= {168.719, 167.019, 191.879, 190.060, 196.152, 200.024, 166.695, 181.571, 248.206, 252.113, 226.179, 224.148, 223.302, 200.685, 186.035, 181.939, 181.975};
std::array<double, nofAveMasses> aveTopMass_KFchi2cut15= {168.719, 166.848, 192.016, 190.654, 195.740, 200.236, 167.958, 182.345, 250.376, 251.513, 223.667, 222.226, 222.102, 206.881, 187.017, 182.734, 182.772};
// No cut on KF chi2
std::array<double, nofAveMasses> aveTopMass_widthx0p5  = {166.680, 166.147, 202.979, 204.490, 196.533, 200.905, 179.901, 192.179, 259.186, 257.377, 232.275, 230.366, 241.172, 235.975, 200.818, 192.694, 192.770};
std::array<double, nofAveMasses> aveTopMass_widthx0p66 = {167.605, 167.067, 203.807, 205.278, 197.529, 201.867, 181.030, 193.041, 259.186, 257.377, 232.275, 230.366, 241.172, 235.975, 200.818, 193.547, 193.615};
std::array<double, nofAveMasses> aveTopMass_widthx0p75 = {167.986, 167.434, 204.120, 205.564, 197.958, 202.287, 181.491, 193.373, 259.186, 257.377, 232.275, 230.366, 241.172, 235.975, 200.818, 193.876, 193.940};
std::array<double, nofAveMasses> aveTopMass_widthx1    = {168.723, 168.111, 204.631, 205.986, 198.843, 203.183, 182.341, 193.935, 259.186, 257.377, 232.275, 230.366, 241.172, 235.975, 200.818, 194.432, 194.492};
std::array<double, nofAveMasses> aveTopMass_widthx2    = {169.970, 168.992, 204.813, 205.741, 200.850, 205.472, 183.275, 194.316, 259.186, 257.377, 232.275, 230.366, 241.172, 235.975, 200.818, 194.810, 194.866};
std::array<double, nofAveMasses> aveTopMass_widthx3    = {170.764, 169.371, 204.597, 205.079, 202.538, 207.556, 183.453, 194.266, 259.186, 257.377, 232.275, 230.366, 241.172, 235.975, 200.818, 194.760, 194.817};
std::array<double, nofAveMasses> aveTopMass_widthx4    = {169.941, 168.900, 202.590, 202.886, 201.317, 206.443, 181.884, 192.657, 259.186, 257.377, 232.275, 230.366, 241.172, 235.975, 200.818, 193.168, 193.240};
void getAveMasses(double width)
{
  for (int i = 0; i < nofAveMasses; i++)
  {
//     if ( width == 0.5 )       aveTopMass[i] = aveTopMass_widthx0p5[i];
//     else if ( width == 0.66 ) aveTopMass[i] = aveTopMass_widthx0p66[i];
//     else if ( width == 0.75 ) aveTopMass[i] = aveTopMass_widthx0p75[i];
//     else if ( width == 1. )   aveTopMass[i] = aveTopMass_widthx1[i];
//     else if ( width == 2. )   aveTopMass[i] = aveTopMass_widthx2[i];
//     else if ( width == 3. )   aveTopMass[i] = aveTopMass_widthx3[i];
//     else if ( width == 4. )   aveTopMass[i] = aveTopMass_widthx4[i];
//     else
//     {
//       if ( i == 0)
//         cout << "Average top mass for width " << width << " not found. Using average mass for width = 1..." << endl;
//       aveTopMass[i] = aveTopMass_widthx1[i];
//     }
    ///*if (! applyWidthSF)*/ aveTopMass[i] = aveTopMass_noWidth[i];
    aveTopMass[i] = aveTopMass_KFchi2cut5[i];
  }
}


// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;

vector < Dataset* > datasets;

ofstream txtMassGenMatched, txtMassRecoCP, txtMassRecoWPUP, txtMassRecoUP, txtMassRecoWP, txtMassRecoWPWOk, txtMassRecoWPWNotOk, txtMassReco;

/// Function prototypes
struct HighestPt
{
  bool operator()( TLorentzVector j1, TLorentzVector j2 ) const
  {
    return j1.Pt() > j2.Pt() ;
  }
};

struct greater
{
  template<class T>
  bool operator()(T const &a, T const &b) const { return a > b; }
};


string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();
void GetMetaData(TTree* tree, bool isData);
void InitTree(TTree* tree, bool isData);
void InitMSPlots();
void InitHisto1D();
void InitHisto2D();
void InitHisto1DMatch();
void InitHisto2DMatch();
void TruthMatching(vector<TLorentzVector> partons, vector<TLorentzVector> selectedJets, pair<unsigned int, unsigned int> *MCPermutation);
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearVars();
void ClearObjects();
void FillGeneralPlots(int d);
void FillMatchingPlots();
void FillKinFitPlots(int d);
long GetNEvents(TTree* fChain, string var, bool isData);
void GetHLTFraction(double* fractions);
double BreitWigner(double topPT, double scale);
double eventWeightCalculator(double topPT, double scale);
Double_t voigt(Double_t *x, Double_t *par);
Double_t crysBall_WP(Double_t *x);
Double_t crysBall_UP(Double_t *x);
Double_t crysBall_WPUP(Double_t *x);
Double_t logLikelihood(Double_t *x, Double_t *par);
Double_t fakeLikelihood(Double_t *x, Double_t *par);
Double_t widthToGammaTranslation(Double_t *x);



// Declaration of leaf types
Int_t           run_num;
Long64_t        evt_num;
Int_t           lumi_num;
Int_t           nvtx;
Int_t           npu;
Double_t        rho;
Bool_t          isTrigged;
Bool_t          hasExactly4Jets;
Bool_t          hasJetLeptonCleaning;
Int_t           nLeptons;
Int_t           nMuons;
Int_t           muon_charge[1];   //[nMuons]
Double_t        muon_pt[1];   //[nMuons]
Double_t        muon_phi[1];   //[nMuons]
Double_t        muon_eta[1];   //[nMuons]
Double_t        muon_E[1];   //[nMuons]
Double_t        muon_M[1];   //[nMuons]
Double_t        muon_d0[1];   //[nMuons]
Double_t        muon_chargedHadronIso[1];   //[nMuons]
Double_t        muon_neutralHadronIso[1];   //[nMuons]
Double_t        muon_photonIso[1];   //[nMuons]
Double_t        muon_puChargedHadronIso[1];   //[nMuons]
Double_t        muon_relIso[1];   //[nMuons]
Double_t        muon_pfIso[1];   //[nMuons]
Int_t           nJets;
Int_t           jet_charge[20];   //[nJets]
Double_t        jet_pt[20];   //[nJets]
Double_t        jet_phi[20];   //[nJets]
Double_t        jet_eta[20];   //[nJets]
Double_t        jet_E[20];   //[nJets]
Double_t        jet_M[20];   //[nJets]
Double_t        jet_bdiscr[20];   //[nJets]
Double_t        met_pt;
Double_t        met_phi;
Double_t        met_eta;
Double_t        met_Et;
Double_t        met_E;
Int_t           nMCParticles;
Int_t           mc_status[200];   //[nMCParticles]
Int_t           mc_pdgId[200];   //[nMCParticles]
Int_t           mc_mother[200];   //[nMCParticles]
Int_t           mc_granny[200];   //[nMCParticles]
Double_t        mc_pt[200];   //[nMCParticles]
Double_t        mc_phi[200];   //[nMCParticles]
Double_t        mc_eta[200];   //[nMCParticles]
Double_t        mc_E[200];   //[nMCParticles]
Double_t        mc_M[200];   //[nMCParticles]
Bool_t          mc_isLastCopy[200];   //[nMCParticles]
Bool_t          mc_isPromptFinalState[200];   //[nMCParticles]
Bool_t          mc_isHardProcess[200];   //[nMCParticles]
Bool_t          mc_fromHardProcessFinalState[200];   //[nMCParticles]
Double_t        nloWeight;
Double_t        puSF;
Double_t        btagSF;
Double_t        muonIdSF[1];   //[nMuons]
Double_t        muonIsoSF[1];   //[nMuons]
Double_t        muonTrigSFv2[1];   //[nMuons]
Double_t        muonTrigSFv3[1];   //[nMuons]

Long64_t        nEvents;
Long64_t        nEventsSel;
Int_t           cutFlow[10];
Int_t           appliedJER;
Int_t           appliedJES;
Int_t           appliedPU;
Long64_t        nofEventsHLTv2;
Long64_t        nofEventsHLTv3;
Long64_t        nofSelEventsHLTv2;
Long64_t        nofSelEventsHLTv3;

// List of branches
TBranch        *b_run_num;   //!
TBranch        *b_evt_num;   //!
TBranch        *b_lumi_num;   //!
TBranch        *b_nvtx;   //!
TBranch        *b_npu;   //!
TBranch        *b_rho;   //!
TBranch        *b_isTrigged;   //!
TBranch        *b_hasExactly4Jets;   //!
TBranch        *b_hasJetLeptonCleaning;   //!
TBranch        *b_nLeptons;   //!
TBranch        *b_nMuons;   //!
TBranch        *b_muon_charge;   //!
TBranch        *b_muon_pt;   //!
TBranch        *b_muon_phi;   //!
TBranch        *b_muon_eta;   //!
TBranch        *b_muon_E;   //!
TBranch        *b_muon_M;   //!
TBranch        *b_muon_d0;   //!
TBranch        *b_muon_chargedHadronIso;   //!
TBranch        *b_muon_neutralHadronIso;   //!
TBranch        *b_muon_photonIso;   //!
TBranch        *b_muon_puChargedHadronIso;   //!
TBranch        *b_muon_relIso;   //!
TBranch        *b_muon_pfIso;   //!
TBranch        *b_nJets;   //!
TBranch        *b_jet_charge;   //!
TBranch        *b_jet_pt;   //!
TBranch        *b_jet_phi;   //!
TBranch        *b_jet_eta;   //!
TBranch        *b_jet_E;   //!
TBranch        *b_jet_M;   //!
TBranch        *b_jet_bdiscr;   //!
TBranch        *b_met_pt;   //!
TBranch        *b_met_phi;   //!
TBranch        *b_met_eta;   //!
TBranch        *b_met_Et;   //!
TBranch        *b_met_E;   //!
TBranch        *b_nMCParticles;   //!
TBranch        *b_mc_status;   //!
TBranch        *b_mc_pdgId;   //!
TBranch        *b_mc_mother;   //!
TBranch        *b_mc_granny;   //!
TBranch        *b_mc_pt;   //!
TBranch        *b_mc_phi;   //!
TBranch        *b_mc_eta;   //!
TBranch        *b_mc_E;   //!
TBranch        *b_mc_M;   //!
TBranch        *b_mc_isLastCopy;   //!
TBranch        *b_mc_isPromptFinalState;   //!
TBranch        *b_mc_isHardProcess;   //!
TBranch        *b_mc_fromHardProcessFinalState;   //!
TBranch        *b_nloWeight;   //!
TBranch        *b_puSF;   //!
TBranch        *b_btagSF;   //!
TBranch        *b_muonIdSF;   //!
TBranch        *b_muonIsoSF;   //!
TBranch        *b_muonTrigSFv2;   //!
TBranch        *b_muonTrigSFv3;   //!

TBranch        *b_nEvents;   //!
TBranch        *b_nEventsSel;   //!
TBranch        *b_cutFlow;   //!
TBranch        *b_appliedJER;   //!
TBranch        *b_appliedJES;   //!
TBranch        *b_appliedPU;   //!
TBranch        *b_nofEventsHLTv2;   //!
TBranch        *b_nofEventsHLTv3;   //!
TBranch        *b_nofSelEventsHLTv2;   //!
TBranch        *b_nofSelEventsHLTv3;   //

double scaleFactor, widthSF;
vector<unsigned int> bJetId;
vector<int> selectedJetsCharge;
vector<double> selectedJetsBDiscr;
double bdiscrTop, bdiscrTop2, tempbdiscr;
int labelB1, labelB2;
int labelsReco[4];
double reco_minMlb, reco_minMl_nonb, reco_ttbarMass, reco_dRLepB_lep, reco_dRLepB_had;
double min_Mlb, dRLepB;
int labelMlb, labelMl_nonb;
double massForWidth;


/// Define TLVs
TLorentzVector muon, jet, mcpart;
TLorentzVector WCandidate;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
vector<TLorentzVector> selectedJetsKFcorrected;
vector<TLorentzVector> mcParticles;
vector<TLorentzVector> partons;
vector<TLorentzVector> partonsMatched;
vector<TLorentzVector> jetsMatched;

/// Matching
int pdgID_top = 6; //top quark

bool doMatching = true;
bool all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
bool all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
bool hadronicTopJetsMatched = false;
bool hadronicTopJetsMatched_MCdef_ = false;
pair<unsigned int, unsigned int> MCPermutation[4] = {pair<unsigned int,unsigned int>(9999,9999)};
int topQuark = -9999, antiTopQuark = -9999;
int genmuon = -9999;
bool muonmatched = false;
bool muPlusFromTop = false, muMinusFromTop = false;
vector<unsigned int> partonId;

double matchedWMass_reco, matchedTopMass_reco, matchedTopMass_gen;
double matchedMlb_corr, matchedMlb_wrong, matchedTtbarMass_corr, matchedTtbarMass_wrong;
double matchedDRLepB_corr, matchedDRLepB_wrong;

/// KinFitter
TKinFitter* kFitter;
bool addWMassKF = true;
bool addEqMassKF = false;
int kFitVerbosity = 0;
double kFitChi2 = 99.;
int nofAcceptedKFit = 0;
double Wmass_reco_orig, Wmass_reco_kf, topmass_reco_orig, topmass_reco_kf;
double toppt_reco_orig, toppt_reco_kf;

/// Likelihood
int nTot = 0;
Double_t f_CP = 1./3., f_WP = 1./3., f_UP = 1./3.;
Double_t widthArray[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2., 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.5, 2.55, 2.60, 2.65, 2.70, 2.75, 3., 3.25, 3.5, 3.75, 4.};
string widthArrayStr[] = {"0p5", "0p55", "0p6", "0p65", "0p7", "0p75", "0p8", "0p85", "0p9", "0p95", "1", "1p05", "1p1", "1p15", "1p2", "1p25", "1p3", "1p35", "1p4", "1p45", "1p5", "1p55", "1p6", "1p65", "1p7", "1p75", "1p8", "1p85", "1p9", "1p95", "2", "2p05", "2p1", "2p15", "2p2", "2p25", "2p3", "2p35", "2p4", "2p45", "2p5", "2p55", "2p6", "2p65", "2p7", "2p75", "3", "3p25", "3p5", "3p75", "4"};
const int nWidthsLL = sizeof(widthArray)/sizeof(widthArray[0]);
Double_t gammaArray[nWidthsLL];

Double_t loglike[nWidthsLL] = {0};
Double_t loglike_pd[10][nWidthsLL] = {{0}};
Double_t loglike_per_evt[nWidthsLL] = {0};
Double_t loglike_onlyGoodEvts[nWidthsLL] = {0};
Double_t fakelike_CP[nWidthsLL] = {0};
Double_t fakelike_CP_per_evt[nWidthsLL] = {0};
Double_t fakelike_onlyGoodEvts[nWidthsLL] = {0};
Double_t fakelike_CP_Res[nWidthsLL] = {0};

bool isGoodLL = false;
int nofGoodEvtsLL[10] = {0};
int nofBadEvtsLL[10] = {0};
Double_t aveTopMassLL = -1.;
Double_t maxMtDivAveMt = 0., minMtDivAveMt = 9999.;

ofstream txtLogLike, txtLogLikeTest, txtOutputLogLike;

/// Toys
TRandom3 random3;
double toy; // = random3.Uniform(0,1);
double toyMax;
int nToys[10] = {0}, nDataEvts;

/// Meta
string strSyst = "";
double eqLumi;
vector<int> vJER, vJES, vPU;

bool CharSearch( char str[], char substr[] )
{
  char *output = NULL;
  output = strstr(str, substr);
  if (output) return true;
  else return false;
}

int main(int argc, char* argv[])
{
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  string channel;
  if ( argc == 1) channel = "mu";
  //else if ( CharSearch(argv[1], "mu") || CharSearch(argv[1], "Mu") || CharSearch(argv[1], "MU") || CharSearch(argv[1], "m") ) channel = "mu";
  //else if ( CharSearch(argv[1], "el") || CharSearch(argv[1], "El") || CharSearch(argv[1], "EL") || CharSearch(argv[1], "e") ) channel = "el";
  //else if ( (argv[1]).find("all") != std::string::npos || (argv[1]).find("All") != std::string::npos || (argv[1]).find("ALL") != std::string::npos ) channel = "all";
  
  if (calculateAverageMass)
  {
    calculateLikelihood = false;
    useToys = false;
    makePlots = false;
  }
  if (test)
  {
    makePlots = false;
  }
  
  //string pathOutput = "test/";
  string pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (makePlots)
  {
    // Add channel to output path
    pathOutput += channel+"/";
    mkdir(pathOutput.c_str(),0777);
    if (useToys)
    {
      pathOutput+= "toys/";
      mkdir(pathOutput.c_str(),0777);
    }
    if (testHistos)
    {
      pathOutput += "test/";
      mkdir(pathOutput.c_str(),0777);
    }
    // Give timestamp to output path
    pathOutput += dateString+"/";
    mkdir(pathOutput.c_str(),0777);
  }
  
  pathNtuples = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate+"/";
  cout << "Using Ntuples from " << ntupleDate << ". This corresponds to systematics: " << systStr << endl;
  if ( applyWidthSF && scaleWidth != 1 ) cout << "TTbar sample width will be scaled by a factor " << scaleWidth << endl;
  if (calculateAverageMass) cout << "Calculating average mass values..." << endl;
  if (calculateLikelihood) cout << "Calculating -loglikelihood values..." << endl;
  if (testHistos) cout << "Testing histogram consistency..." << endl;
  
  
  /// xml file
  string xmlFileName ="config/topWidth.xml";
  
  const char *xmlfile = xmlFileName.c_str();
  
  cout << " - Using config file " << xmlfile << endl;
  
  
  
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  datasets.clear();
  TTreeLoader treeLoader; 
  treeLoader.LoadDatasets(datasets, xmlfile);
  
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    string dataSetName = datasets[d]->Name();
    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA") != std::string::npos )
      Luminosity = datasets[d]->EquivalentLumi();
    
    if ( dataSetName.find("QCD") != std::string::npos )
      datasets[d]->SetColor(kYellow);
    if ( dataSetName.find("TT") != std::string::npos )
    {
      datasets[d]->SetTitle("t#bar{t}");
      datasets[d]->SetColor(kRed+1);
    }
    //if ( dataSetName.find("TTbarJets_Other") != std::string::npos ) datasets[d]->SetColor(kRed-7);
    if ( dataSetName.find("WJets") != std::string::npos )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      datasets[d]->SetColor(kGreen-3);
    }
    if ( dataSetName.find("ZJets") != std::string::npos || dataSetName.find("DY") != std::string::npos )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      datasets[d]->SetColor(kAzure-2);
      //datasets[d]->SetColor(kMagenta);
    }
    if ( dataSetName.find("ST") != std::string::npos || dataSetName.find("SingleTop") != std::string::npos )
    {
      datasets[d]->SetTitle("ST");
      datasets[d]->SetColor(kBlue-2);
      //if ( dataSetName.find("tW") != std::string::npos )
      //{
      //  datasets[d]->SetTitle("ST tW");
      //  datasets[d]->SetColor(kBlue-4);
      //}
      //else
      //  datasets[d]->SetTitle("ST t");
    }
  }
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  ResolutionFunctions* rf = new ResolutionFunctions(calculateResolutionFunctions, true);
  KinFitter *kf = new KinFitter("PlotsForResolutionFunctions_testFit.root", addWMassKF, addEqMassKF);
    
  if (makePlots)
  {
    InitMSPlots();
    InitHisto1D();
    InitHisto2D();
  }
  
  vJER.clear(); vJES.clear(); vPU.clear();
  
  /// Load ave top mass
  getAveMasses(scaleWidth);
  
  if (calculateLikelihood)
  {
    txtLogLike.open(("likelihood_per_event_"+dateString+".txt").c_str());
    txtLogLike << "## -Log(likelihood) values per event" << endl;
    txtLogLike << "#  Widths : ";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      txtLogLike << widthArray[iWidth] << "  ";
    }
    txtLogLike << endl;
    
    txtLogLikeTest.open(("likelihood_per_event_test_"+dateString+".txt").c_str());
    txtLogLikeTest << "## ievt      recoTopMass     M_lb    dR(b,lep)    dR(b_h,lep)" << endl;
    txtLogLikeTest << "#  Gammas : ";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      gammaArray[iWidth] = widthToGammaTranslation(&widthArray[iWidth]);
      txtLogLikeTest << gammaArray[iWidth] << "  ";
    }
    
    /// Fraction of events that is CP, WP and UP
    nTot = nCP + nWP + nUP;
    f_CP = (double)nCP/(double)nTot;
    f_WP = (double)nWP/(double)nTot;
    f_UP = (double)nUP/(double)nTot;
    
    /// Average top mass for likelihood
    aveTopMassLL = aveTopMass[1];
    //aveTopMassLL = f_CP*aveTopMass[1] + f_WP*aveTopMass[4] + f_UP*aveTopMass[3];
    //aveTopMassLL = f_CP*aveTopMass[1] + (f_WP+f_UP)*aveTopMass[2];
    cout << "LogLike::Average top quark mass = " << aveTopMassLL << endl;
  }
  
  if (calculateAverageMass)
  {
    txtMassGenMatched.open(("averageMass/mass_gen_matched_TT_"+dateString+".txt").c_str());
    txtMassRecoCP.open(("averageMass/mass_reco_matched_TT_"+dateString+".txt").c_str());
    txtMassRecoWPUP.open(("averageMass/mass_reco_notCorrectMatch_TT_"+dateString+".txt").c_str());
    txtMassRecoUP.open(("averageMass/mass_reco_notMatched_TT_"+dateString+".txt").c_str());
    txtMassRecoWP.open(("averageMass/mass_reco_wrongPerm_TT_"+dateString+".txt").c_str());
    txtMassRecoWPWOk.open(("averageMass/mass_reco_wrongPerm_WOk_TT_"+dateString+".txt").c_str());
    txtMassRecoWPWNotOk.open(("averageMass/mass_reco_wrongPerm_WNotOk_TT_"+dateString+".txt").c_str());
  }
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName, slumi;
  double timePerDataSet[datasets.size()];
  
  int nEntries;
  double fracHLT[2] = {-1};
  GetHLTFraction(fracHLT);
  if ( fracHLT[0] == -1 || fracHLT[1] == -1 )
  {
    cout << "Something went wrong with the fraction calculation for trigger SFs!" << endl;
    exit(1);
  }
  cout << "The muon trigger scale factors will be scaled by " << fracHLT[0] << " for HLTv2 and " << fracHLT[1] << " for HLTv3." << endl;
  
  
  
  /// Loop over datasets
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    clock_t startDataSet = clock();
    
    ClearMetaData();
    
    dataSetName = datasets[d]->Name();
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
    }
    
    isData = false; isTTbar = false;
    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA") != std::string::npos )
    {
      isData = true;
    }
    else if ( dataSetName.find("TT") != std::string::npos )
    {
      isTTbar = true;
    }
    
    if (testTTbarOnly && ! isTTbar ) continue;
    
    
    if (calculateAverageMass)
    {
      txtMassReco.open(("averageMass/mass_reco_"+dataSetName+"_"+dateString+".txt").c_str());
    }
    
    if (calculateLikelihood)
    {
      txtLogLike << endl << "#  Dataset: " << dataSetName << endl;
      txtLogLikeTest << endl << "#  Dataset: " << dataSetName << endl;
    }
    
    string ntupleFileName = "Ntuples_"+dataSetName+".root";
    tFileMap[dataSetName.c_str()] = new TFile((pathNtuples+ntupleFileName).c_str(),"READ"); //create TFile for each dataset
    
    string tTreeName = "tree";
    string tStatsTreeName = "stats";
    
    /// Get meta data
    tStatsTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tStatsTreeName.c_str());
    GetMetaData(tStatsTree[dataSetName.c_str()], isData);
    
    tStatsTree[(dataSetName).c_str()]->GetEntry(0);
    vJER.push_back(appliedJER);
    vJES.push_back(appliedJES);
    vPU.push_back(appliedPU);
    
    if (useToys)
    {
      eqLumi = (double)GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData)/datasets[d]->Xsection();
      if (isData) nDataEvts = GetNEvents(tStatsTree[(dataSetName).c_str()], "nEventsSel", 1);
      else toyMax = Luminosity/eqLumi;
      
      cout << "TOYS::Number of selected data events: " << nDataEvts << endl;
      if (test)
        cout << "      Lumi : " << Luminosity << "; eqLumi: " << (double)datasets[d]->EquivalentLumi() << "; " << eqLumi << endl;
      cout << "TOYS::Lumi/eqLumi = " << toyMax << endl;
    }
    
    
    /// Get data
    tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
    nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
    cout << "                nEntries: " << nEntries << endl;
    
    
    // Set branch addresses and branch pointers
    InitTree(tTree[dataSetName.c_str()], isData);
    
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    int endEvent = nEntries;
    if (test || testHistos) endEvent = 2001;
    for (int ievt = 0; ievt < endEvent; ievt++)
    {
      ClearObjects();
      
      if (ievt%10000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
      
      //if (! isTTbar) continue;
      //if (ievt > 588000) break;
      //if (ievt != 2062 && ievt != 75831 && ievt != 113603 && ievt != 115687 && ievt != 155732 && ievt != 163161 && ievt != 186900 && ievt != 215759 && ievt != 233634 && ievt != 238021 && ievt != 243052 && ievt != 243674 && ievt != 266399 && ievt != 317190 && ievt != 317752 && ievt != 325854 && ievt != 330813 && ievt != 333620 && ievt != 347247 && ievt != 439571 && ievt != 450329 && ievt != 491328 && ievt != 510024 && ievt != 514196 && ievt != 538345 && ievt != 570225 && ievt != 576194 && ievt != 577278 && ievt != 587570) continue;
      
      
      /// Load event
      tTree[(dataSetName).c_str()]->GetEntry(ievt);
      
      
      /// Toys
      if (useToys)
      {
        toy = random3.Rndm();
        if ( toy > toyMax ) continue;
        (nToys[d])++;
      }
      
      
      /// Scale factors
      if (! isData)
      {
        if (applyLeptonSF) { scaleFactor *= muonIdSF[0] * muonIsoSF[0] * (fracHLT[0]*muonTrigSFv2[0] + fracHLT[1]*muonTrigSFv3[0]);}
        if (applyPU) { scaleFactor *= puSF;}
        if (applyBTagSF) { scaleFactor *= btagSF;}
        //if (applyNloSF) { scaleFactor *= nloWeight;}  // additional SF due to number of events with neg weight!!
        //cout << "Scalefactor: " << setw(6) << scaleFactor << "  btag SF: " << setw(6) << btagSF << "  pu SF: " << setw(6) << puSF << "  muonId: " << setw(6) << muonIdSF[0] << "  muonIso: " << setw(6) << muonIsoSF[0] << "  muonTrig: " << setw(6) << fracHLT[0]*muonTrigSFv2[0] + fracHLT[1]*muonTrigSFv3[0] << endl;
      }
      
      if (useToys) scaleFactor *= eqLumi;  // undo automatic scaling by eqLumi in MSPlots
      
      
      
      //////////////////////
      ///  Fill objects  ///
      //////////////////////
      
      muon.SetPtEtaPhiE(muon_pt[0], muon_eta[0], muon_phi[0], muon_E[0]);
      selectedLepton.push_back(muon);
      
      for (int iJet = 0; iJet < nJets; iJet++)
      {
        jet.Clear();
        jet.SetPtEtaPhiE(jet_pt[iJet], jet_eta[iJet], jet_phi[iJet], jet_E[iJet]);
        if ( jet.Pt() > 30. )
        {
          selectedJets.push_back(jet);
          selectedJetsCharge.push_back(jet_charge[iJet]);
          selectedJetsBDiscr.push_back(jet_bdiscr[iJet]);
        }
      }
      
      if ( selectedJets.size() > 4 ) continue;
      
      for (int iJet = 0; iJet < selectedJets.size(); iJet++)
      {
        if ( jet_bdiscr[iJet] > CSVv2Medium )
        {
          selectedBJets.push_back(selectedJets[iJet]);
          bJetId.push_back(iJet);  /// selectedBJets[j] = selectedJets[bJetId[j]]
                                   /// b discr of selectedBJets[j] = selectedJetsBDiscr[bJetId[j]]
        }
      }
      //std::sort(selectedBJets.begin(),selectedBJets.end(),HighestPt());  // already the case
      
      if ( selectedBJets.size() < 2 ) continue;
      
      /// Count nof events with exactly 4 jets with pT > 30 GeV of which 2 are b tagged
      nofHardSelected++;
      
      /// label jets with highest b discr
      for (int iJet = 0; iJet < selectedBJets.size(); iJet++)
      {
        tempbdiscr = selectedJetsBDiscr[bJetId[iJet]];
        if ( tempbdiscr > bdiscrTop )
        {
          bdiscrTop2 = bdiscrTop;
          bdiscrTop = tempbdiscr;

          labelB2 = labelB1;
          labelB1 = iJet;
        }
        else if ( tempbdiscr > bdiscrTop2 )
        {
          bdiscrTop2 = tempbdiscr;
          labelB2 = iJet;
        }
      }
        
      
      /// Make plots
      if (makePlots)
      {
        FillGeneralPlots(d);
      }
      
      
      
      /////////////////////////////
      ///  JET PARTON MATCHING  ///
      /////////////////////////////
      
      //if ( isTTbar || dataSetName.find("ST") != std::string::npos )  // no matches for ST
      if (isTTbar)
      {
        for (int iMC = 0; iMC < nMCParticles; iMC++)
        {
          mcpart.Clear();
          mcpart.SetPtEtaPhiE(mc_pt[iMC], mc_eta[iMC], mc_phi[iMC], mc_E[iMC]);
          mcParticles.push_back(mcpart);
        }
        
        
        for (unsigned int i = 0; i < mcParticles.size(); i++)
        {
          if ( test && verbose > 4 )
            cout << setw(3) << right << i << "  Status: " << setw(2) << mc_status[i] << "  pdgId: " << setw(3) << mc_pdgId[i] << "  Mother: " << setw(4) << mc_mother[i] << "  Granny: " << setw(4) << mc_granny[i] << "  Pt: " << setw(7) << left << mc_pt[i] << "  Eta: " << mc_eta[i] << endl;
          
          /// Find tops
          if ( mc_pdgId[i] == pdgID_top )  // isLastCopy() == status 62
          {
            if ( mc_status[i] == 22 ) topQuark = i;
            if ( topQuark == -9999 && mc_status[i] == 62 ) topQuark = i;
          }
          else if ( mc_pdgId[i] == -pdgID_top )
          {
            if ( mc_status[i] == 22 ) antiTopQuark = i;
            if ( antiTopQuark == -9999 && mc_status[i] == 62 ) antiTopQuark = i;
          }
          
          
          /// Status restriction: Final state particle or particle from hardest process
          if ( (mc_status[i] > 1 && mc_status[i] <= 20) || mc_status[i] >= 30 ) continue;
          
          /// Muons
          if ( mc_pdgId[i] == 13 && mc_mother[i] == -24 && mc_granny[i] == -pdgID_top )		// mu-, W-, tbar
          {
            muMinusFromTop = true;
            if ( mc_status[i] == 23 ) genmuon = i;
            else if ( mc_status[i] != 23 && genmuon == -9999 ) genmuon = i;
          }
          if ( mc_pdgId[i] == -13 && mc_mother[i] == 24 && mc_granny[i] == pdgID_top )		// mu+, W+, t
          {
            muPlusFromTop = true;
            if ( mc_status[i] == 23 ) genmuon = i;
            else if ( mc_status[i] != 23 && genmuon == -9999 ) genmuon = i;
          }
          
          /// Partons/gluons
          if ( abs(mc_pdgId[i]) < 6 || abs(mc_pdgId[i]) == 21 )  //light/b quarks, 6 should stay hardcoded, OR gluon
          {
            partons.push_back(mcParticles[i]);
            partonId.push_back(i);  /// partons[j] = mcParticles[partonId[j]]
          }
          
        }  // end loop mcParticles

        if (verbose > 3)
        {
          cout << "Size mcParticles:   " << mcParticles.size() << endl;
          cout << "Size partons:       " << partons.size() << endl;
          cout << "Size selectedJets:  " << selectedJets.size() << endl;
        }
        
        if ( muMinusFromTop && muPlusFromTop )
        {
          if (test) cout << "Both tops decay leptonically... Event " << ievt << " will not be matched." << endl;
          doMatching = false;
        }
        
        
        
        /////////////////////////////////////////
        ///  Scale factor ttbar sample width  ///
        /////////////////////////////////////////
        
        if (applyWidthSF)
        {
          if ( muon_charge[0] > 0 ) massForWidth = (mcParticles[antiTopQuark]).M();
          else if ( muon_charge[0] < 0 ) massForWidth = (mcParticles[topQuark]).M();
          
          widthSF = eventWeightCalculator(massForWidth, scaleWidth);
          
          if ( widthSF != widthSF )  // widthSF = NaN
          {
            cout << "Event " << ievt << ";  Event nb. " << evt_num << "; Run nb. " << run_num << endl;
            cout << "Top mass: " << (mcParticles[topQuark]).M() << "; antiTop mass: " << (mcParticles[antiTopQuark]).M() << "; Lepton charge: " << muon_charge[0] << "; width SF: " << widthSF << endl;
            
            continue;
          }
        }  // end applyWidthSF
        
        
        
        //////////////////
        ///  Matching  ///
        //////////////////
        
        if (doMatching)
        {
          TruthMatching(partons, selectedJets, MCPermutation);
          
          if (all4PartonsMatched && test && verbose > 3)
          {
            for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
            {
              //cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Status: " << setw(2) << mc_status[partonId[MCPermutation[iMatch].second]] << "  pdgId: " << setw(3) << mc_pdgId[partonId[MCPermutation[iMatch].second]] << "  Mother: " << setw(4) << mc_mother[partonId[MCPermutation[iMatch].second]] << "  Granny: " << setw(4) << mc_granny[partonId[MCPermutation[iMatch].second]] << endl;
              cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Pt: " << setw(7) << left << mc_pt[partonId[MCPermutation[iMatch].second]] << "  Eta: " << mc_eta[partonId[MCPermutation[iMatch].second]] << "  Phi: " << mc_phi[partonId[MCPermutation[iMatch].second]] << endl;
              cout << "Event  " << right << setw(4) << ievt << ";  Matched jet    " << iMatch << "  Pt: " << setw(7) << left << jet_pt[MCPermutation[iMatch].first] << "  Eta: " << jet_eta[MCPermutation[iMatch].first] << "  Phi: " << jet_phi[MCPermutation[iMatch].first] << endl;
            }
          }
          
          
          ///////////////////
          ///  Resolution functions
          ///////////////////
          
          if (hadronicTopJetsMatched)
          {
            
            for (unsigned int iMatch = 0; iMatch < 3; iMatch++)
            {
              /// MCPermutation[i].first  = jet number
              /// MCPermutation[i].second = parton number
              /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
              
              partonsMatched.push_back(partons[MCPermutation[iMatch].second]);
              jetsMatched.push_back(selectedJets[MCPermutation[iMatch].first]);
            }
            if (all4PartonsMatched)
            {
              partonsMatched.push_back(partons[MCPermutation[3].second]);
              jetsMatched.push_back(selectedJets[MCPermutation[3].first]);
            }
            
            if (all4PartonsMatched && calculateResolutionFunctions)
            {
              rf->fillJets(partonsMatched, jetsMatched);
              
              if (muonmatched) rf->fillMuon(mcParticles[genmuon], selectedLepton[0]);
              //if (electronmatched) rf->fillElectron(...)

            }  // end rf
          }
          
          if (hadronicTopJetsMatched)
          {  
            matchedWMass_reco = (jetsMatched[0] + jetsMatched[1]).M();
            matchedTopMass_reco = (jetsMatched[0] + jetsMatched[1] + jetsMatched[2]).M();
            matchedTopMass_gen = (partonsMatched[0] + partonsMatched[1] + partonsMatched[2]).M();
            
            if (calculateAverageMass) txtMassGenMatched << ievt << "  " << matchedWMass_reco << "  " << matchedTopMass_reco << "  " << widthSF << endl;
            
            if (makePlots)
            {
              FillMatchingPlots();
            }

          }  // end hadronicTopJetsMatched
        
        }  // end doMatching
        
        
      }  // end if TT
      
      
      
      /////////////////////////////////////
      ///  Reconstruction of Top quark  ///
      /////////////////////////////////////
      
      //int labelsReco[4] = {-9999, -9999, -9999, -9999};  // 0,1: light jets; 2: hadronic b; 3: leptonic b
      double deltaR;
      double minDeltaR = 9999.;
      
      for (int kjet = 0; kjet < selectedBJets.size(); kjet++)
      {
        if ( selectedBJets.size() > 2 && kjet != labelB1 && kjet != labelB2 ) continue;
        
        deltaR = ROOT::Math::VectorUtil::DeltaR(selectedBJets[kjet], selectedLepton[0]);
        
        if (deltaR < minDeltaR)
        {
          minDeltaR = deltaR;
          labelsReco[3] = bJetId[kjet];
        }
      }
      
      if ( labelsReco[3] == -9999 ) continue;
      
      if ( labelsReco[3] == bJetId[labelB1] ) labelsReco[2] = bJetId[labelB2];
      else if ( labelsReco[3] == bJetId[labelB2] ) labelsReco[2] = bJetId[labelB1];
      else cerr << endl << "Seems like something went wrong with the b jets..." << endl;
      
      for (int ijet = 0; ijet < selectedJets.size(); ijet++)
      {
        if ( ijet == labelsReco[2] || ijet == labelsReco[3] ) continue;
        
        if ( labelsReco[0] == -9999 ) labelsReco[0] = ijet;
        else if ( labelsReco[1] == -9999 ) labelsReco[1] = ijet;
        else cerr << endl << "Seems like there are too many jets..." << endl;
      }
      
      if ( labelsReco[0] == -9999 || labelsReco[1] == -9999 || labelsReco[2] == -9999 ) continue;
      
      
      
      ////////////////////////////
      ///   Kinematic Fit      ///
      ////////////////////////////
      
      if (doKinFit)
      {
        /*if (addEqMassKF)
          kFitter = kf->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]], selectedJets[labelsReco[2]], selectedJets[labelsReco[3]], selectedLepton[0], TLV NEUTRINO, kFitVerbosity);
        else if (addWMassKF)*/
          kFitter = kf->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]], kFitVerbosity);
        
        if ( kFitter->getStatus() != 0 )  // did not converge
        {
          if (test) cout << "Fit did not converge..." << endl;
          continue;
        }
        
        kFitChi2 = kFitter->getS();
        if (test) cout << "Fit converged: Chi2 = " << kFitChi2 << endl;
        
        if ( kFitChi2 > 5. ) continue;
        nofAcceptedKFit++;
        
        selectedJetsKFcorrected.clear();
        selectedJetsKFcorrected = kf->getCorrectedJets();
        
      }
      
      
      /// Reconstruct event
      //  Hadronic variables  // OBS: only W mass constraint ! Jet3 = selectedJets[labelsReco[2]] !
      if (! doKinFit)
      {
        selectedJetsKFcorrected.clear();
        selectedJetsKFcorrected.push_back(selectedJets[labelsReco[0]]);
        selectedJetsKFcorrected.push_back(selectedJets[labelsReco[1]]);
        selectedJetsKFcorrected.push_back(selectedJets[labelsReco[2]]);
      }
      if ( selectedJetsKFcorrected.size() == 2 ) selectedJetsKFcorrected.push_back(selectedJets[labelsReco[2]]);
      
      Wmass_reco_orig = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).M();
      Wmass_reco_kf = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).M();
      topmass_reco_orig = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).M();
      topmass_reco_kf = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).M();
      toppt_reco_orig = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).Pt();
      toppt_reco_kf = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).Pt();
      
      if ( topmass_reco_kf < 0. )
      {
        cout << endl <<"Event " << setw(5) << right << ievt << "   ";
        cout << "Top mass after kinFit is negative! I.e. " << topmass_reco_kf << " Before kinFit: " << topmass_reco_orig << "  ievt " << ievt << endl;
        cout << "Mass jet 1 & jet 2: " << (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).M() << "; Mass jet 1 & jet 3: " << (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[2]).M() << "; Mass jet 2 & jet 3: " << (selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).M() << endl;
        //if (test)
        //{
          cout << "Original:   Jet 1: pT " << selectedJets[labelsReco[0]].Pt() << "; Jet 2: pT " << selectedJets[labelsReco[1]].Pt() << "; Jet 3: pT " << selectedJets[labelsReco[2]].Pt() << endl;
          cout << "Corrected:  Jet 1: pT " << selectedJetsKFcorrected[0].Pt() << "; Jet 2: pT " << selectedJetsKFcorrected[1].Pt() << "; Jet 3: pT " << selectedJetsKFcorrected[2].Pt() << endl;
          cout << "Original:   Jet 1: px " << selectedJets[labelsReco[0]].Px() << "; py " << selectedJets[labelsReco[0]].Py() << "; pz " << selectedJets[labelsReco[0]].Pz() << "; E " << selectedJets[labelsReco[0]].E() << endl;
          cout << "Corrected:  Jet 1: px " << selectedJetsKFcorrected[0].Px() << "; py " << selectedJetsKFcorrected[0].Py() << "; pz " << selectedJetsKFcorrected[0].Pz() << "; E " << selectedJetsKFcorrected[0].E() << endl;
          cout << "Original:   Jet 2: px " << selectedJets[labelsReco[1]].Px() << "; py " << selectedJets[labelsReco[1]].Py() << "; pz " << selectedJets[labelsReco[1]].Pz() << "; E " << selectedJets[labelsReco[1]].E() << endl;
          cout << "Corrected:  Jet 2: px " << selectedJetsKFcorrected[1].Px() << "; py " << selectedJetsKFcorrected[1].Py() << "; pz " << selectedJetsKFcorrected[1].Pz() << "; E " << selectedJetsKFcorrected[1].E() << endl;
          cout << "Original:   Jet 1: pt " << selectedJets[labelsReco[0]].Pt() << "; eta " << selectedJets[labelsReco[0]].Eta() << "; phi " << selectedJets[labelsReco[0]].Phi() << "; M " << selectedJets[labelsReco[0]].M() << endl;
          cout << "Corrected:  Jet 1: pt " << selectedJetsKFcorrected[0].Pt() << "; eta " << selectedJetsKFcorrected[0].Eta() << "; phi " << selectedJetsKFcorrected[0].Phi() << "; M " << selectedJetsKFcorrected[0].M() << endl;
          cout << "Original:   Jet 2: pt " << selectedJets[labelsReco[1]].Pt() << "; eta " << selectedJets[labelsReco[1]].Eta() << "; phi " << selectedJets[labelsReco[1]].Phi() << "; M " << selectedJets[labelsReco[1]].M() << endl;
          cout << "Corrected:  Jet 2: pt " << selectedJetsKFcorrected[1].Pt() << "; eta " << selectedJetsKFcorrected[1].Eta() << "; phi " << selectedJetsKFcorrected[1].Phi() << "; M " << selectedJetsKFcorrected[1].M() << endl;
          cout << "Original:   Jet 1: x " << selectedJets[labelsReco[0]].X() << "; y " << selectedJets[labelsReco[0]].Y() << "; z " << selectedJets[labelsReco[0]].Z() << "; t " << selectedJets[labelsReco[0]].T() << endl;
          cout << "Corrected:  Jet 1: x " << selectedJetsKFcorrected[0].X() << "; y " << selectedJetsKFcorrected[0].Y() << "; z " << selectedJetsKFcorrected[0].Z() << "; t " << selectedJetsKFcorrected[0].T() << endl;
          cout << "Original:   Jet 2: x " << selectedJets[labelsReco[1]].X() << "; y " << selectedJets[labelsReco[1]].Y() << "; z " << selectedJets[labelsReco[1]].Z() << "; t " << selectedJets[labelsReco[1]].T() << endl;
          cout << "Corrected:  Jet 2: x " << selectedJetsKFcorrected[1].X() << "; y " << selectedJetsKFcorrected[1].Y() << "; z " << selectedJetsKFcorrected[1].Z() << "; t " << selectedJetsKFcorrected[1].T() << endl;
        //}
      }
      
      if ( doKinFit && makePlots )
      {
        FillKinFitPlots(d);
      }
      
      
      //  Leptonic variables
      reco_minMlb = (selectedLepton[0] + selectedJets[labelsReco[3]]).M();
      reco_dRLepB_lep = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[3]], selectedLepton[0] );  // deltaR between lepton and leptonic b jet
      reco_dRLepB_had = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
      
      //  Combination
      reco_ttbarMass = reco_minMlb + topmass_reco_kf;
      
      //Average mass
      if (calculateAverageMass) txtMassReco << ievt << "  " << Wmass_reco_kf << "  " << topmass_reco_kf << "  " << widthSF << endl;
      
      //Fill histos
      if (makePlots)
      {
        MSPlot["Reco_W_mass"]->Fill(Wmass_reco_kf, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Reco_hadTop_mass"]->Fill(topmass_reco_kf, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Reco_hadTop_pT"]->Fill(toppt_reco_kf, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Reco_mlb"]->Fill(reco_minMlb, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Reco_ttbar_mass"]->Fill(reco_ttbarMass, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Reco_dR_lep_b"]->Fill(reco_dRLepB_lep, datasets[d], true, Luminosity*scaleFactor);
        if (isTTbar)
        {
          histo1D["Reco_W_mass_reco"]->Fill(Wmass_reco_kf, widthSF);
          histo1D["Reco_top_mass_reco"]->Fill(topmass_reco_kf, widthSF);
          histo1D["dR_lep_b_reco"]->Fill(reco_dRLepB_lep, widthSF);
        }
        
        MSPlot["Reco_mTop_div_aveMTopMatch"]->Fill(topmass_reco_kf/aveTopMass[0], datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Reco_mTop_div_aveMTopTTReco"]->Fill(topmass_reco_kf/aveTopMass[7], datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Reco_mTop_div_aveMTopAllReco"]->Fill(topmass_reco_kf/aveTopMass[15], datasets[d], true, Luminosity*scaleFactor);
        
        if (isData)
        {
          histo1D[("mTop_div_aveMTop_"+dataSetName).c_str()]->Fill(topmass_reco_kf/aveTopMass[14]);
          MSPlot["Reco_mTop_div_aveMTop_stack"]->Fill(topmass_reco_kf/aveTopMass[14], datasets[d], true, Luminosity*scaleFactor*widthSF);
        }
        else
        {
          histo1D[("mTop_div_aveMTop_"+dataSetName).c_str()]->Fill(topmass_reco_kf/aveTopMass[d+6]);
          MSPlot["Reco_mTop_div_aveMTop_stack"]->Fill(topmass_reco_kf/aveTopMass[d+6], datasets[d], true, Luminosity*scaleFactor*widthSF);
          histo1D["mTop_div_aveMTop_allMCReco"]->Fill(topmass_reco_kf/aveTopMass[15]);
        }
        
      }  // end fill plots
      
      
      //Likelihood
      if (calculateLikelihood)
      {
        double tempAveMass = topmass_reco_kf/aveTopMassLL;
        if ( tempAveMass > maxMtDivAveMt ) maxMtDivAveMt = tempAveMass;
        if ( tempAveMass < minMtDivAveMt ) minMtDivAveMt = tempAveMass;
        for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
        {
          loglike_per_evt[iWidth] = logLikelihood(&tempAveMass, &gammaArray[iWidth]);
          loglike[iWidth] += loglike_per_evt[iWidth];
          loglike_pd[d][iWidth] += loglike_per_evt[iWidth];
        }

        /// make loglikelihood only with events that have minimum
        isGoodLL = false;
        for (int iWidth = 1; iWidth < nWidthsLL-1; iWidth++)
        {
          if ( loglike_per_evt[0] > loglike_per_evt[iWidth] && loglike_per_evt[iWidth] < loglike_per_evt[nWidthsLL-1] )
          {
            isGoodLL = true;
            break;
          }
        }
        if (isGoodLL)
        {
          nofGoodEvtsLL[d]++;
          for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
          {
            loglike_onlyGoodEvts[iWidth] += loglike_per_evt[iWidth];
          }
        }
        else
        {
          nofBadEvtsLL[d]++;
        }

        /// write debug file
        if ( ( isData && ievt%1000 == 0 )
             || ( isTTbar && ievt%100000 == 0 )
             || ( ! isData && ! isTTbar && ievt%100 == 0 ) )
        {
          txtLogLike << ievt << "  ";
          for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
          {
            txtLogLike << loglike_per_evt[iWidth] << ",  ";

          }
          txtLogLike << topmass_reco_kf << endl;
        }
      }
      
      
      
      ///////////////////////////////////
      ///  CHECK MATCHED COMBINATION  ///
      ///////////////////////////////////
      
      ///
      // 3 possibilities:
      // - correct top match: 3 jets selected with reco method correspond to the 3 matched jets (n.b. this is also true when the jets originating from the W boson and the b jet do not exactly correspond to the matched jets, because we are only interested in the reconstructed top quark.)
      // - wrong permutation: the correct jet combination exists in the selected jets, but is not chosen by the reco method.
      // - wrong (no) match:  the correct jet combination does not exist in the selected jets (e.g. when one jet is not selected.)
      
      double tempAveMass = topmass_reco_kf/aveTopMass[1];
      
      //if (isTTbar)
      if (! isData)
      {
        if (! applyWidthSF ) widthSF = 1.;
        
        if (hadronicTopJetsMatched)
        {
          /// Correct match
          if ( ( labelsReco[0] == MCPermutation[0].first || labelsReco[0] == MCPermutation[1].first || labelsReco[0] == MCPermutation[2].first ) && ( labelsReco[1] == MCPermutation[0].first || labelsReco[1] == MCPermutation[1].first || labelsReco[1] == MCPermutation[2].first ) && ( labelsReco[2] == MCPermutation[0].first || labelsReco[2] == MCPermutation[1].first || labelsReco[2] == MCPermutation[2].first ) )  // correct jets for top quark
          {
            nofCorrectlyMatched++;
            if (calculateAverageMass) txtMassRecoCP << ievt << "  " << Wmass_reco_kf << "  " << topmass_reco_kf << "  " << widthSF << endl;
            else if (makePlots)
            {
              histo1D["mTop_div_aveMTop_TT_reco_CP"]->Fill(topmass_reco_kf/aveTopMass[1], widthSF);
              histo1D["minMlb_reco_CP"]->Fill(reco_minMlb, widthSF);
              histo1D["dR_lep_b_lep_reco_CP"]->Fill(reco_dRLepB_lep, widthSF);
              histo1D["dR_lep_b_had_reco_CP"]->Fill(reco_dRLepB_had, widthSF);
              histo1D["ttbar_mass_reco_CP"]->Fill(reco_ttbarMass, widthSF);
              histo2D["dR_lep_b_lep_vs_had_CP"]->Fill(reco_dRLepB_lep, reco_dRLepB_had, widthSF);
              histo2D["ttbar_mass_vs_minMlb_CP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_lep < 3. )
                histo2D["ttbar_mass_vs_minMlb_dRlepCut_CP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_lep < 2. )
                histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_CP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_had > 1.2 )
                histo2D["ttbar_mass_vs_minMlb_dRhadCut_CP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_had > 2. )
                histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_CP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_lep < 3. && reco_dRLepB_had > 1.2 )
                histo2D["ttbar_mass_vs_minMlb_dRBothCuts_CP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_lep < 2. && reco_dRLepB_had > 2. )
                histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_CP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if (doKinFit) MSPlot["KF_top_mass_corr_CP"]->Fill(topmass_reco_kf, datasets[d], true, Luminosity*scaleFactor*widthSF);
              if (doKinFit)
              {
                histo1D["KF_top_mass_corr_CP"]->Fill(topmass_reco_kf, widthSF);
                histo1D["KF_Chi2_CP"]->Fill(kFitChi2);
              }
            }
            
            if ( tempAveMass > 0.5 && tempAveMass < 1.5)
            {
              nofCP++;
              if (isTTbar) nofCP_TT++;
            }
            
            if (calculateLikelihood)
            {
              
              for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
              {
                fakelike_CP_per_evt[iWidth] = fakeLikelihood(&tempAveMass, &gammaArray[iWidth]);
                fakelike_CP[iWidth] += fakelike_CP_per_evt[iWidth];
                if ( fabs(selectedJetsKFcorrected[2].Et()-mcParticles[partonId[MCPermutation[2].second]].Et()) < 0.005 * selectedJetsKFcorrected[2].Et() )
                  fakelike_CP_Res[iWidth] += fakelike_CP_per_evt[iWidth];
              }

              /// make loglikelihood only with events that have minimum
              isGoodLL = false;
              for (int iWidth = 1; iWidth < nWidthsLL-1; iWidth++)
              {
                if ( fakelike_CP_per_evt[0] > fakelike_CP_per_evt[iWidth] && fakelike_CP_per_evt[iWidth] < fakelike_CP_per_evt[nWidthsLL-1] )
                {
                  isGoodLL = true;
                  break;
                }
              }
              if (isGoodLL)
              {
                for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
                {
                  fakelike_onlyGoodEvts[iWidth] += fakelike_CP_per_evt[iWidth];
                }
              }
            }
            
            /// write debug file
            if ( calculateLikelihood && ! isGoodLL /*&& ( ( isData && nofGoodEvtsLL[d]%500 == 0 )
                 || ( isTTbar && nofGoodEvtsLL[d]%50000 == 0 )
                 || ( ! isData && ! isTTbar && nofGoodEvtsLL[d]%50 == 0 ) )*/ )
            {
              if (test)
              {
                txtLogLikeTest << setw(8) << right << ievt << "  CP  ";
                txtLogLikeTest << topmass_reco_kf << "  " << reco_minMlb << "  " << reco_dRLepB_lep << "  " << reco_dRLepB_had << endl;
              }
              if (makePlots)
              {
                histo1D["debugLL_hadr_top_mass_reco"]->Fill(topmass_reco_kf);
                histo1D["debugLL_minMlb_reco"]->Fill(reco_minMlb);
                histo1D["debugLL_ttbar_mass_reco"]->Fill(reco_ttbarMass);
                histo1D["debugLL_dR_lep_b_lep_reco"]->Fill(reco_dRLepB_lep);
                histo1D["debugLL_dR_lep_b_had_reco"]->Fill(reco_dRLepB_had);
              }
            }
            
          }  // end corr match
          else  // wrong permutation
          {
            nofNotCorrectlyMatched++;
            if (calculateAverageMass)
            {
              txtMassRecoWP << ievt << "  " << Wmass_reco_kf << "  " << topmass_reco_kf << "  " << widthSF << endl;
              txtMassRecoWPUP << ievt << "  " << Wmass_reco_kf << "  " << topmass_reco_kf << "  " << widthSF << endl;
            }
            else if (makePlots)
            {
              histo1D["mTop_div_aveMTop_TT_reco_WP"]->Fill(topmass_reco_kf/aveTopMass[1], widthSF);
              histo1D["mTop_div_aveMTop_TT_reco_WPUP"]->Fill(topmass_reco_kf/aveTopMass[1], widthSF);
              histo1D["minMlb_reco_WP"]->Fill(reco_minMlb, widthSF);
              histo1D["dR_lep_b_lep_reco_WP"]->Fill(reco_dRLepB_lep, widthSF);
              histo1D["dR_lep_b_had_reco_WP"]->Fill(reco_dRLepB_had, widthSF);
              histo1D["ttbar_mass_reco_WP"]->Fill(reco_ttbarMass, widthSF);
              histo2D["dR_lep_b_lep_vs_had_WP"]->Fill(reco_dRLepB_lep, reco_dRLepB_had, widthSF);
              histo2D["ttbar_mass_vs_minMlb_WP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_lep < 3. )
                histo2D["ttbar_mass_vs_minMlb_dRlepCut_WP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_lep < 2. )
                histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_WP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_had > 1.2 )
                histo2D["ttbar_mass_vs_minMlb_dRhadCut_WP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_had > 2. )
                histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_WP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_lep < 3. && reco_dRLepB_had > 1.2 )
                histo2D["ttbar_mass_vs_minMlb_dRBothCuts_WP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if ( reco_dRLepB_lep < 2. && reco_dRLepB_had > 2. )
                histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_WP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
              if (doKinFit) MSPlot["KF_top_mass_corr_WP"]->Fill(topmass_reco_kf, datasets[d], true, Luminosity*scaleFactor*widthSF);
              if (doKinFit)
              {
                histo1D["KF_top_mass_corr_WP"]->Fill(topmass_reco_kf, widthSF);
                histo1D["KF_Chi2_WP"]->Fill(kFitChi2);
              }
            }
            
            if ( tempAveMass > 0.5 && tempAveMass < 1.5)
            {
              nofWP++;
              if (isTTbar) nofWP_TT++;
            }
            
            /// write debug file
            if ( test && calculateLikelihood && ! isGoodLL /*&& ( ( isData && nofGoodEvtsLL[d]%500 == 0 )
                 || ( isTTbar && nofGoodEvtsLL[d]%50000 == 0 )
                 || ( ! isData && ! isTTbar && nofGoodEvtsLL[d]%50 == 0 ) )*/ )
            {
              txtLogLikeTest << setw(8) << right << ievt << "  WP  ";
              txtLogLikeTest << topmass_reco_kf << "  " << reco_minMlb << "  " << reco_dRLepB_lep << "  " << reco_dRLepB_had << endl;
            }

            if ( ( labelsReco[0] == MCPermutation[0].first || labelsReco[0] == MCPermutation[1].first || labelsReco[0] == MCPermutation[2].first ) && ( labelsReco[1] == MCPermutation[0].first || labelsReco[1] == MCPermutation[1].first || labelsReco[1] == MCPermutation[2].first ) )
            {
              if (calculateAverageMass) txtMassRecoWPWOk << ievt << "  " << Wmass_reco_kf << "  " << topmass_reco_kf << "  " << widthSF << endl;
              else if (makePlots)
                histo1D["mTop_div_aveMTop_TT_reco_WP_WOk"]->Fill(topmass_reco_kf/aveTopMass[5], widthSF);
            }
            else
            {
              if (calculateAverageMass) txtMassRecoWPWNotOk << ievt << "  " << Wmass_reco_kf << "  " << topmass_reco_kf << "  " << widthSF << endl;
              else if (makePlots)
                histo1D["mTop_div_aveMTop_TT_reco_WP_WNotOk"]->Fill(topmass_reco_kf/aveTopMass[6], widthSF);
            }
          }  // end wrong perm
        }  // end hadrTopMatch
        else  // no match
        {
          if (calculateAverageMass)
          {
            txtMassRecoUP << ievt << "  " << Wmass_reco_kf << "  " << topmass_reco_kf << "  " << widthSF << endl;
            txtMassRecoWPUP << ievt << "  " << Wmass_reco_kf << "  " << topmass_reco_kf << "  " << widthSF << endl;
          }
          else if (makePlots)
          {
            histo1D["mTop_div_aveMTop_TT_reco_UP"]->Fill(topmass_reco_kf/aveTopMass[1], widthSF);
            histo1D["mTop_div_aveMTop_TT_reco_WPUP"]->Fill(topmass_reco_kf/aveTopMass[1], widthSF);
            histo1D["minMlb_reco_UP"]->Fill(reco_minMlb, widthSF);
            histo1D["dR_lep_b_lep_reco_UP"]->Fill(reco_dRLepB_lep, widthSF);
            histo1D["dR_lep_b_had_reco_UP"]->Fill(reco_dRLepB_had, widthSF);
            histo1D["ttbar_mass_reco_UP"]->Fill(reco_ttbarMass, widthSF);
            histo2D["dR_lep_b_lep_vs_had_UP"]->Fill(reco_dRLepB_lep, reco_dRLepB_had, widthSF);
            histo2D["ttbar_mass_vs_minMlb_UP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
            if ( reco_dRLepB_lep < 3. )
              histo2D["ttbar_mass_vs_minMlb_dRlepCut_UP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
            if ( reco_dRLepB_lep < 2. )
              histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_UP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
            if ( reco_dRLepB_had > 1.2 )
              histo2D["ttbar_mass_vs_minMlb_dRhadCut_UP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
            if ( reco_dRLepB_had > 2. )
              histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_UP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
            if ( reco_dRLepB_lep < 3. && reco_dRLepB_had > 1.2 )
              histo2D["ttbar_mass_vs_minMlb_dRBothCuts_UP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
            if ( reco_dRLepB_lep < 2. && reco_dRLepB_had > 2. )
              histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_UP"]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
            if (doKinFit) MSPlot["KF_top_mass_corr_UP"]->Fill(topmass_reco_kf, datasets[d], true, Luminosity*scaleFactor*widthSF);
            if (doKinFit)
            {
              histo1D["KF_top_mass_corr_UP"]->Fill(topmass_reco_kf, widthSF);
              histo1D["KF_Chi2_UP"]->Fill(kFitChi2);
            }
          }
          
          if ( tempAveMass > 0.5 && tempAveMass < 1.5)
          {
            nofUP++;
            if (isTTbar) nofUP_TT++;
          }
          
          /// write debug file
          if ( test && calculateLikelihood && ! isGoodLL /*&& ( ( isData && nofGoodEvtsLL[d]%500 == 0 )
               || ( isTTbar && nofGoodEvtsLL[d]%50000 == 0 )
               || ( ! isData && ! isTTbar && nofGoodEvtsLL[d]%50 == 0 ) )*/ )
          {
            txtLogLikeTest << setw(8) << right << ievt << "  UP  ";
            txtLogLikeTest << topmass_reco_kf << "  " << reco_minMlb << "  " << reco_dRLepB_lep << "  " << reco_dRLepB_had << endl;
          }
            
        }  // end no match
        //if (test) cout << "checked match" << endl;
      }  // end TT / ! isData
      
      if (makePlots)
      {
        histo1D["mTop_div_aveMTop_bkgd"]->Fill(topmass_reco_kf/aveTopMass[15]);
      }
      
      
    }  // end loop events
    
    
    cout << endl;  /// Stronger selection in this analyser compared to Ntuples ==> endEvent --> nofHardSelected
    cout << "Number of events with exactly 4 jets with pT > 30 GeV: " << nofHardSelected << " (" << 100*((float)nofHardSelected/(float)endEvent) << "%)" << endl;
    cout << "Number of events accepted by kinFitter: " << nofAcceptedKFit << " (" << 100*((float)nofAcceptedKFit/(float)nofHardSelected) << "%)" << endl;
    
    //if ( isTTbar || dataSetName.find("ST") != std::string::npos )
    if (isTTbar)
    {
      cout << "Number of matched events: " << setw(8) << right << nofMatchedEvents << endl;
      cout << "Number of events with hadronic top matched: " << setw(8) << right << nofHadrMatchedEvents << endl;
      cout << "Correctly matched reconstructed events:     " << setw(8) << right << nofCorrectlyMatched << endl;
      cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatched << endl;
      if ( nofCorrectlyMatched != 0 || nofNotCorrectlyMatched != 0 )
        cout << "   ===> This means that " << 100*(float)nofCorrectlyMatched / (float)(nofCorrectlyMatched + nofNotCorrectlyMatched) << "% is correctly matched." << endl;
      
      /// Resolution functions
      if (calculateResolutionFunctions)
      {
        string rfFileName = "PlotsForResolutionFunctions.root";
        string rfFitFileName = "PlotsForResolutionFunctions_Fitted.root";
        TFile *foutRF = new TFile(rfFileName.c_str(), "RECREATE");
        foutRF->cd();

        rf->writeHistograms();

        foutRF->Close();
        
        rf->makeFit(rfFileName, rfFitFileName);
        rf->writeTable(rfFitFileName);

        delete foutRF;
      }
      
    }  // end TT
    
    if (calculateLikelihood)
    {
      cout << "Number of events with min in likelihood    " << setw(8) << right << nofGoodEvtsLL[d] << endl;
      cout << "Number of events without min in likelihood " << setw(8) << right << nofBadEvtsLL[d] << endl;
      if ( nofGoodEvtsLL[d] != 0 || nofBadEvtsLL[d] != 0 )
        cout << "   ===> " << 100*(float)nofGoodEvtsLL[d] / (float)(nofGoodEvtsLL[d] + nofBadEvtsLL[d]) << "% are 'good'." << endl;
    }
    
    if (useToys) cout << "TOYS::" << datasets[d]->Name() << ": " << nToys[d] << endl;
    
    if (calculateAverageMass) txtMassReco.close();
    
    tFileMap[dataSetName.c_str()]->Close();
    
    timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
    
  }  // end loop datasets
  
  if (calculateAverageMass)
  {
    txtMassGenMatched.close();
    txtMassRecoCP.close();
    txtMassRecoWPUP.close();
    txtMassRecoUP.close();
    txtMassRecoWP.close();
    txtMassRecoWPWOk.close();
    txtMassRecoWPWNotOk.close();
  }
  
  cout << "Number of events with 0.5 < mt/<mt> < 1.5 : CP: " << nofCP << " (" << (double)nofCP/((double)(nofCP+nofWP+nofUP)) << "%)   WP: " << nofWP << " (" << (double)nofWP/((double)(nofCP+nofWP+nofUP)) << "%)   UP: " << nofUP << " (" << (double)nofUP/((double)(nofCP+nofWP+nofUP)) << "%)" << endl;
  cout << "                               (TTbar only) CP: " << nofCP_TT << "  WP: " << nofWP_TT << "  UP: " << nofUP_TT << endl;
  
  if (calculateLikelihood)
  {
    txtLogLike.close();
    txtLogLikeTest.close();
    
    cout << "Minimum top mass divided by average top mass: " << minMtDivAveMt << endl;
    cout << "Maximum top mass divided by average top mass: " << maxMtDivAveMt << endl;
    
    cout << endl << "likelihood values (all samples) : {";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      cout << loglike[iWidth]/(1e+6);
      if ( iWidth != nWidthsLL-1 ) cout << ", ";
    }
    cout << "} *10^6;" << endl << "likelihood values (data-only) : ";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      cout << loglike_pd[0][iWidth];
      if ( iWidth != nWidthsLL-1 ) cout << ", ";
    }
    cout << endl << "likelihood values (only good events) : {";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      cout << loglike_onlyGoodEvts[iWidth]/(1e+3);
      if ( iWidth != nWidthsLL-1 ) cout << ", ";
    }
    cout << "} *10^3;" << endl << "fake likelihood values (CP) : {";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      cout << fakelike_CP[iWidth]/(1e+6);
      if ( iWidth != nWidthsLL-1 ) cout << ", ";
    }
    cout << "} *10^6;" << endl << "fake likelihood values (CP - only good events) : {";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      cout << fakelike_onlyGoodEvts[iWidth]/(1e+4);
      if ( iWidth != nWidthsLL-1 ) cout << ", ";
    }
    cout << "} *10^4;" << endl << "fake likelihood values (CP) with E_t,jet - E_t,q < 0.005*E_t,jet : {";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      cout << fakelike_CP_Res[iWidth];
      if ( iWidth != nWidthsLL-1 ) cout << ", ";
    }
    cout << "};" << endl << "widths: ";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      cout << widthArray[iWidth];
      if ( iWidth != nWidthsLL-1 ) cout << ", ";
    }
    cout << endl << "gammas: ";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      cout << gammaArray[iWidth];
      if ( iWidth != nWidthsLL-1 ) cout << ", ";
    }
    cout << endl << endl;
    
    /// Print output to file
    txtOutputLogLike.open("output_loglikelihood_gamma.txt");
    txtOutputLogLike << "Widths:      ";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      txtOutputLogLike << setw(5) << right << widthArray[iWidth] << "  ";
    }
    txtOutputLogLike << endl << "Gammas:      ";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      txtOutputLogLike << setw(5) << right << gammaArray[iWidth] << "  ";
    }
    txtOutputLogLike << endl << "LLikelihood: ";
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      txtOutputLogLike << setw(12) << right << loglike_onlyGoodEvts[iWidth] << "  ";
    }
    txtOutputLogLike << endl;
    txtOutputLogLike.close();
  }
  
  cout << "Processing time per dataset: " << endl;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
  }
  
  if (useToys)
  {
    cout << "TOYS::Number of selected data events: " << nDataEvts << endl;
    for (unsigned int d = 0; d < datasets.size(); d++)
    {
      cout << "TOYS::" << datasets[d]->Name() << ": " << nToys[d] << endl;
    }
  }
  
  //////////////////////////////////////////
  ///  Check Shape Changing Systematics  ///
  //////////////////////////////////////////
  
  int sumJER = 0, sumJES = 0, sumPU = 0;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    sumJER += vJER[d];
    sumJES += vJES[d];
    sumPU  += vPU[d];
  }
  
  if ( sumJER == 0 && sumJES == 0 && sumPU == 0 )
  {
    strSyst = "nominal";
  }
  else if ( abs(sumJER) == datasets.size() && sumJES == 0 && sumPU == 0 )
  {
    if ( sumJER > 0 ) { strSyst = "JERup";}
    else if ( sumJER < 0 ) { strSyst = "JERdown";}
  }
  else if ( sumJER == 0 && abs(sumJES) == datasets.size() && sumPU == 0 )
  {
    if ( sumJES > 0 ) { strSyst = "JESup";}
    else if ( sumJES < 0 ) { strSyst = "JESdown";}
  }
  else if ( sumJER == 0 && sumJES == 0 && abs(sumPU) == datasets.size() )
  {
    if ( sumPU > 0 ) { strSyst = "PUup";}
    else if ( sumPU < 0 ) { strSyst = "PUdown";}
  }
  else if (! testTTbarOnly)
  {
    cerr << "Shape changing systematics not consistent accross datasets or multiple applied at once" << endl;
    cerr << "Exiting...." << endl;
    exit(1);
  }
  if (! testTTbarOnly) cout << " - Systematics confirmed to be " << strSyst << endl;
  
  
  if (test)
  {
    cout << "Exiting because of test..." << endl;
    exit(1);
  }
  else if (calculateAverageMass)
  {
    cout << "Average mass calculated. Exiting..." << endl;
    exit(1);
  }
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  string rootFileName = "NtuplePlots_"+systStr+".root";
  mkdir((pathOutput+"MSPlot/").c_str(),0777);
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile ((pathOutput+rootFileName).c_str(), "RECREATE");
  cout << "   Output file is " << pathOutput+rootFileName << endl;
  
  ///Write histograms
  fout->cd();
  
  for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    //cout << "MSPlot: " << it->first << endl;
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(name, 1, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
    temp->Write(fout, name, true, pathOutput+"MSPlot", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
  }
  
  // 1D
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  th1dir->cd();
  gStyle->SetOptStat(1111);
  for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathOutput+it->first+".png").c_str() );
  }
  
  // 2D
  TDirectory* th2dir = fout->mkdir("2D_histograms");
  th2dir->cd();
  gStyle->SetPalette(55);
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
    TH2F *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first, "colz");
    tempCanvas->SaveAs( (pathOutput+it->first+".png").c_str() );
  }
  
  fout->Close();
  
  delete fout;
  
  
  
  double time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  if ( time >= 60 )
  {
    int mins = time/60;
    double secs = time - mins*60;
    
    if (mins >= 60 )
    {
      int hours = mins/60;
      mins = mins - hours*60;
      cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
    }
    else
      cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
  }
    
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  cout << " - Goodbye" << endl;
  
  return 0;
  
}  // end main



/// Functions

string ConvertIntToString(int Number, int pad)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  if ( pad > 1 ) { convert << std::setw(pad) << std::setfill('0');}
  convert << Number;
  return convert.str();
}

string MakeTimeStamp()
{
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  int year = now->tm_year - 100;  /// + 1900 to get current year
  int month = now->tm_mon + 1;
  int day = now->tm_mday;
  int hour = now->tm_hour;
  int min = now->tm_min;
  //int sec = now->tm_sec;
  
  string year_str = ConvertIntToString(year, 2);
  string month_str = ConvertIntToString(month, 2);
  string day_str = ConvertIntToString(day, 2);
  string hour_str = ConvertIntToString(hour, 2);
  string min_str = ConvertIntToString(min, 2);
  //string sec_str = ConvertIntToString(sec, 2);
  
  string date_str = year_str + month_str + day_str + "_" + hour_str + min_str;
  return date_str;
}


void GetMetaData(TTree* tree, bool isData)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  
  tree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
  tree->SetBranchAddress("nEventsSel", &nEventsSel, &b_nEventsSel);
  tree->SetBranchAddress("cutFlow", cutFlow, &b_cutFlow);
  tree->SetBranchAddress("appliedJER", &appliedJER, &b_appliedJER);
  tree->SetBranchAddress("appliedJES", &appliedJES, &b_appliedJES);
  tree->SetBranchAddress("appliedPU", &appliedPU, &b_appliedPU);
  if (isData)
  {
   tree->SetBranchAddress("nofEventsHLTv2", &nofEventsHLTv2, &b_nofEventsHLTv2);
   tree->SetBranchAddress("nofEventsHLTv3", &nofEventsHLTv3, &b_nofEventsHLTv3);
   tree->SetBranchAddress("nofSelEventsHLTv2", &nofSelEventsHLTv2, &b_nofSelEventsHLTv2);
   tree->SetBranchAddress("nofSelEventsHLTv3", &nofSelEventsHLTv3, &b_nofSelEventsHLTv3);
  }
}

void InitTree(TTree* tree, bool isData)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  
  tree->SetBranchAddress("run_num", &run_num, &b_run_num);
  tree->SetBranchAddress("evt_num", &evt_num, &b_evt_num);
  tree->SetBranchAddress("lumi_num", &lumi_num, &b_lumi_num);
  tree->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  tree->SetBranchAddress("npu", &npu, &b_npu);
  tree->SetBranchAddress("rho", &rho, &b_rho);
  tree->SetBranchAddress("isTrigged", &isTrigged, &b_isTrigged);
  tree->SetBranchAddress("hasExactly4Jets", &hasExactly4Jets, &b_hasExactly4Jets);
  tree->SetBranchAddress("hasJetLeptonCleaning", &hasJetLeptonCleaning, &b_hasJetLeptonCleaning);
  tree->SetBranchAddress("nLeptons", &nLeptons, &b_nLeptons);
  tree->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
  tree->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
  tree->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
  tree->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
  tree->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
  tree->SetBranchAddress("muon_E", muon_E, &b_muon_E);
  tree->SetBranchAddress("muon_M", muon_M, &b_muon_M);
  tree->SetBranchAddress("muon_d0", muon_d0, &b_muon_d0);
  tree->SetBranchAddress("muon_chargedHadronIso", muon_chargedHadronIso, &b_muon_chargedHadronIso);
  tree->SetBranchAddress("muon_neutralHadronIso", muon_neutralHadronIso, &b_muon_neutralHadronIso);
  tree->SetBranchAddress("muon_photonIso", muon_photonIso, &b_muon_photonIso);
  tree->SetBranchAddress("muon_puChargedHadronIso", muon_puChargedHadronIso, &b_muon_puChargedHadronIso);
  tree->SetBranchAddress("muon_relIso", muon_relIso, &b_muon_relIso);
  tree->SetBranchAddress("muon_pfIso", muon_pfIso, &b_muon_pfIso);
  tree->SetBranchAddress("nJets", &nJets, &b_nJets);
  tree->SetBranchAddress("jet_charge", jet_charge, &b_jet_charge);
  tree->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
  tree->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
  tree->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
  tree->SetBranchAddress("jet_E", jet_E, &b_jet_E);
  tree->SetBranchAddress("jet_M", jet_M, &b_jet_M);
  tree->SetBranchAddress("jet_bdiscr", jet_bdiscr, &b_jet_bdiscr);
  tree->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
  tree->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
  tree->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
  tree->SetBranchAddress("met_Et", &met_Et, &b_met_Et);
  tree->SetBranchAddress("met_E", &met_E, &b_met_E);
  if (! isData)
  {
    tree->SetBranchAddress("nMCParticles", &nMCParticles, &b_nMCParticles);
    tree->SetBranchAddress("mc_status", mc_status, &b_mc_status);
    tree->SetBranchAddress("mc_pdgId", mc_pdgId, &b_mc_pdgId);
    tree->SetBranchAddress("mc_mother", mc_mother, &b_mc_mother);
    tree->SetBranchAddress("mc_granny", mc_granny, &b_mc_granny);
    tree->SetBranchAddress("mc_pt", mc_pt, &b_mc_pt);
    tree->SetBranchAddress("mc_phi", mc_phi, &b_mc_phi);
    tree->SetBranchAddress("mc_eta", mc_eta, &b_mc_eta);
    tree->SetBranchAddress("mc_E", mc_E, &b_mc_E);
    tree->SetBranchAddress("mc_M", mc_M, &b_mc_M);
    if (! useOldNtuples)
    {
      tree->SetBranchAddress("mc_isLastCopy", mc_isLastCopy, &b_mc_isLastCopy);
      tree->SetBranchAddress("mc_isPromptFinalState", mc_isPromptFinalState, &b_mc_isPromptFinalState);
      tree->SetBranchAddress("mc_isHardProcess", mc_isHardProcess, &b_mc_isHardProcess);
      tree->SetBranchAddress("mc_fromHardProcessFinalState", mc_fromHardProcessFinalState, &b_mc_fromHardProcessFinalState);
    }
    tree->SetBranchAddress("nloWeight", &nloWeight, &b_nloWeight);
    tree->SetBranchAddress("puSF", &puSF, &b_puSF);
    tree->SetBranchAddress("btagSF", &btagSF, &b_btagSF);
    tree->SetBranchAddress("muonIdSF", muonIdSF, &b_muonIdSF);
    tree->SetBranchAddress("muonIsoSF", muonIsoSF, &b_muonIsoSF);
    tree->SetBranchAddress("muonTrigSFv2", muonTrigSFv2, &b_muonTrigSFv2);
    tree->SetBranchAddress("muonTrigSFv3", muonTrigSFv3, &b_muonTrigSFv3);
  }
  
}

void InitMSPlots()
{
  MSPlot["muon_pT"] = new MultiSamplePlot(datasets, "muon_pT", 22, 0, 440, "p_{T} [GeV]");
  MSPlot["muon_eta"] = new MultiSamplePlot(datasets, "muon_eta", 30, -3, 3, "Eta");
  MSPlot["muon_phi"] = new MultiSamplePlot(datasets, "muon_phi", 32, -3.2, 3.2, "Phi");
  MSPlot["muon_relIso"] = new MultiSamplePlot(datasets, "muon_relIso", 20, 0, 0.2, "relIso");
  MSPlot["muon_d0"] = new MultiSamplePlot(datasets, "muon_d0", 50, 0, 0.03, "d_{0}");
  MSPlot["leadingJet_pT"] = new MultiSamplePlot(datasets, "leadingJet_pT", 60, 0, 600, "p_{T} [GeV]");
  MSPlot["jet2_pT"] = new MultiSamplePlot(datasets, "jet2_pT", 40, 0, 400, "p_{T} [GeV]");
  MSPlot["jet3_pT"] = new MultiSamplePlot(datasets, "jet3_pT", 40, 0, 400, "p_{T} [GeV]");
  MSPlot["jet4_pT"] = new MultiSamplePlot(datasets, "jet4_pT", 40, 0, 400, "p_{T} [GeV]");
  MSPlot["Ht_4leadingJets"] = new MultiSamplePlot(datasets,"Ht_4leadingJets", 60, 0, 1200, "H_{T} [GeV]");
  MSPlot["met_pT"] = new MultiSamplePlot(datasets, "met_pT", 40, 0, 400, "p_{T} [GeV]");
  MSPlot["met_eta"] = new MultiSamplePlot(datasets, "met_eta", 30, -3, 3, "Eta");
  MSPlot["met_phi"] = new MultiSamplePlot(datasets, "met_phi", 32, -3.2, 3.2, "Phi");
  
  MSPlot["M3"] = new MultiSamplePlot(datasets, "M3", 40, 60, 460, "M [GeV]");

  MSPlot["nJets"] = new MultiSamplePlot(datasets, "nJets", 13, -0.5, 12.5, "# jets");
  MSPlot["nBJets"] = new MultiSamplePlot(datasets, "nBJets", 9, -0.5, 8.5, "# b jets");
  MSPlot["CSVv2Discr_allJets"] = new MultiSamplePlot(datasets, "CSVv2Discr_allJets", 48, 0.0, 1.2, "CSVv2 discriminant value");
  MSPlot["CSVv2Discr_leadingJet"] = new MultiSamplePlot(datasets, "CSVv2Discr_leadingJet", 48, 0.0, 1.2, "CSVv2 discriminant value of leading jet");
  MSPlot["CSVv2Discr_jet2"] = new MultiSamplePlot(datasets, "CSVv2Discr_jet2", 48, 0.0, 1.2, "CSVv2 discriminant value of jet2");
  MSPlot["CSVv2Discr_jet3"] = new MultiSamplePlot(datasets, "CSVv2Discr_jet3", 48, 0.0, 1.2, "CSVv2 discriminant value of jet3");
  MSPlot["CSVv2Discr_jet4"] = new MultiSamplePlot(datasets, "CSVv2Discr_jet4", 48, 0.0, 1.2, "CSVv2 discriminant value of jet4");
  MSPlot["CSVv2Discr_highest"] = new MultiSamplePlot(datasets, "CSVv2Discr_highest", 48, 0.0, 1.2, "Highest CSVv2 discriminant value");
  MSPlot["CSVv2Discr_jetNb"] = new MultiSamplePlot(datasets, "CSVv2Discr_jetNb", 8, -0.5, 7.5, "Jet number (in order of decreasing p_{T}) with highest CSVv2 discriminant value");
  
  MSPlot["min_Mlb"] = new MultiSamplePlot(datasets, "min_Mlb", 40, 0, 400, "M_{lb} [GeV]");
  MSPlot["dR_Lep_B"] = new MultiSamplePlot(datasets, "dR_Lep_B", 25, 0, 5, "#Delta R(l,b)");
  
  
  /// Reco
  MSPlot["Reco_W_mass"] = new MultiSamplePlot(datasets, "Reco_W_mass", 100, 0, 400, "M_{W} [GeV]");
  MSPlot["Reco_hadTop_mass"] = new MultiSamplePlot(datasets, "Reco_hadTop_mass", 80, 0, 800, "M_{t} [GeV]");
  MSPlot["Reco_hadTop_pT"] = new MultiSamplePlot(datasets, "Reco_hadTop_pT", 80, 0, 800, "p_{T} [GeV]");
  
  MSPlot["Reco_mlb"] = new MultiSamplePlot(datasets, "Reco_mlb", 80, 0, 800, "M_{lb} [GeV]");
  MSPlot["Reco_ttbar_mass"] = new MultiSamplePlot(datasets, "Reco_ttbar_mass", 50, 0, 1000, "M_{t#bar{t}} [GeV]");

  MSPlot["Reco_dR_lep_b"] = new MultiSamplePlot(datasets, "Reco_dR_lep_b", 25, 0, 5, "#Delta R(l,b)");
  
  
  MSPlot["Reco_mTop_div_aveMTopMatch"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass (from matched events)", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
  MSPlot["Reco_mTop_div_aveMTopTTReco"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass (of reconstructed TT events)", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
  MSPlot["Reco_mTop_div_aveMTopAllReco"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass (of reconstructed events from all datasets)", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
  
  MSPlot["Reco_mTop_div_aveMTop_stack"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
  
  /// KinFitter
  MSPlot["KF_W_mass_orig"] = new MultiSamplePlot(datasets, "W mass before kinFitter", 250, 0, 500, "m_{W} [GeV]");
  MSPlot["KF_top_mass_orig"] = new MultiSamplePlot(datasets, "Top mass before kinFitter", 400, 0, 800, "m_{t} [GeV]");
  MSPlot["KF_top_pt_orig"] = new MultiSamplePlot(datasets, "Top pt before kinFitter", 400, 0, 800, "p_{T} [GeV]");
  MSPlot["KF_W_mass_corr"] = new MultiSamplePlot(datasets, "W mass after kinFitter", 250, 0, 500, "m_{W,kf} [GeV]");
  MSPlot["KF_top_mass_corr"] = new MultiSamplePlot(datasets, "Top mass after kinFitter", 400, 0, 800, "m_{t,kf} [GeV]");
  MSPlot["KF_top_pt_corr"] = new MultiSamplePlot(datasets, "Top pt after kinFitter", 400, 0, 800, "p_{T,kf} [GeV]");
  
  MSPlot["KF_top_mass_corr_CP"] = new MultiSamplePlot(datasets, "Top mass after kinFitter for correct match (CP)", 400, 0, 800, "m_{t,kf} [GeV]");
  MSPlot["KF_top_mass_corr_WP"] = new MultiSamplePlot(datasets, "Top mass after kinFitter for wrong permutations (WP)", 400, 0, 800, "m_{t,kf} [GeV]");
  MSPlot["KF_top_mass_corr_UP"] = new MultiSamplePlot(datasets, "Top mass after kinFitter for no match (UP)", 400, 0, 800, "m_{t,kf} [GeV]");
}

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  InitHisto1DMatch();
  
  /// Reco
  histo1D["Reco_W_mass_reco"] = new TH1F("Reco_W_mass_reco","Reconstructed hadronic W mass; M_{W} [GeV]", 80, 0, 800);
  histo1D["Reco_top_mass_reco"] = new TH1F("Reco_top_mass_reco","Reconstructed top mass; M_{t} [GeV]", 80, 0, 800);
    
  /// m_t/<m_t>
  histo1D["mTop_div_aveMTop_TT_matched_jets"] = new TH1F("mTop_div_aveMTop_TT_matched_jets","Top mass divided by average top mass for matched TT sample (using jets from matched partons); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_TT_reco_CP"] = new TH1F("mTop_div_aveMTop_TT_reco_CP","Top mass divided by average top mass for matched TT sample (reco, correct top match); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_TT_reco_WPUP"] = new TH1F("mTop_div_aveMTop_TT_reco_WPUP","Top mass divided by average top mass for unmatched TT sample (reco, no top match & wrong permutations); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_TT_reco_UP"] = new TH1F("mTop_div_aveMTop_TT_reco_UP","Top mass divided by average top mass for unmatched TT sample (reco, no top match); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_TT_reco_WP"] = new TH1F("mTop_div_aveMTop_TT_reco_WP","Top mass divided by average top mass for matched TT sample (reco, wrong top match: wrong permutation); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_TT_reco_WP_WOk"] = new TH1F("mTop_div_aveMTop_TT_reco_WP_WOk","Top mass divided by average top mass for matched TT sample (reco, wrong top match: wrong permutation, W jets correct); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_TT_reco_WP_WNotOk"] = new TH1F("mTop_div_aveMTop_TT_reco_WP_WNotOk","Top mass divided by average top mass for matched TT sample (reco, wrong top match: wrong permutation, W jets not correct); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  
  histo1D["mTop_div_aveMTop_bkgd"] = new TH1F("mTop_div_aveMTop_bkgd","Top mass divided by average top mass for background samples; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_TT"] = new TH1F("mTop_div_aveMTop_TT","Top mass divided by average top mass for TT sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_ST_tW_top"] = new TH1F("mTop_div_aveMTop_ST_tW_top","Top mass divided by average top mass for ST tW top sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_ST_tW_antitop"] = new TH1F("mTop_div_aveMTop_ST_tW_antitop","Top mass divided by average top mass for ST tW antitop sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_ST_t_top"] = new TH1F("mTop_div_aveMTop_ST_t_top","Top mass divided by average top mass for ST t top sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_ST_t_antitop"] = new TH1F("mTop_div_aveMTop_ST_t_antitop","Top mass divided by average top mass for ST t antitop sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_DYJets"] = new TH1F("mTop_div_aveMTop_DYJets","Top mass divided by average top mass for DY+Jets sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_WJets"] = new TH1F("mTop_div_aveMTop_WJets","Top mass divided by average top mass for W+Jets sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_allMCReco"] = new TH1F("mTop_div_aveMTop_allMCReco","Top mass divided by average top mass for all MC samples; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["mTop_div_aveMTop_data"] = new TH1F("mTop_div_aveMTop_data","Top mass divided by average top mass for data sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  
  /// mlb
  histo1D["minMlb_reco_CP"]  = new TH1F("minMlb_reco_CP","Minimal reconstructed M_{lb} mass using events that have correct hadronic top match (CP); min(M_{lb}) [GeV]", 400, 0, 800);
  histo1D["minMlb_reco_WP"]  = new TH1F("minMlb_reco_WP","Minimal reconstructed M_{lb} mass using events that have wrong permutation hadronic top match (WP); min(M_{lb}) [GeV]", 400, 0, 800);
  histo1D["minMlb_reco_UP"]  = new TH1F("minMlb_reco_UP","Minimal reconstructed M_{lb} mass using events that have no hadronic top match (UP); min(M_{lb}) [GeV]", 400, 0, 800);
  
  /// dR
  //  lepton, b(lep)
  histo1D["dR_lep_b_reco"] = new TH1F("dR_lep_b_reco","Minimal delta R between the lepton and a b jet (for reconstructed events); #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_lep_b_lep_reco_CP"]  = new TH1F("dR_lep_b_lep_reco_CP","Minimal delta R between the lepton and the leptonic b jet (reco, correct match); #Delta R(l,b_{l})", 25, 0, 5);
  histo1D["dR_lep_b_lep_reco_WP"]  = new TH1F("dR_lep_b_lep_reco_WP","Minimal delta R between the lepton and the leptonic b jet (reco, wrong permutation); #Delta R(l,b_{l})", 25, 0, 5);
  histo1D["dR_lep_b_lep_reco_UP"]  = new TH1F("dR_lep_b_lep_reco_UP","Minimal delta R between the lepton and the leptonic b jet (reco, no match); #Delta R(l,b_{l})", 25, 0, 5);
  //  lepton, b(hadr)
  histo1D["dR_lep_b_had_reco_CP"]  = new TH1F("dR_lep_b_had_reco_CP","Minimal delta R between the lepton and the hadronic b jet (reco, correct match); #Delta R(l,b_{h})", 25, 0, 5);
  histo1D["dR_lep_b_had_reco_WP"]  = new TH1F("dR_lep_b_had_reco_WP","Minimal delta R between the lepton and the hadronic b jet (reco, wrong permutation); #Delta R(l,b_{h})", 25, 0, 5);
  histo1D["dR_lep_b_had_reco_UP"]  = new TH1F("dR_lep_b_had_reco_UP","Minimal delta R between the lepton and the hadronic b jet (reco, no match); #Delta R(l,b_{h})", 25, 0, 5);
  
  /// ttbar mass
  histo1D["ttbar_mass_reco_CP"] = new TH1F("ttbar_mass_reco_CP","Reconstructed mass of the top quark pair (reco, correct match); M_{t#bar{t}} [GeV]", 500, 0, 1000);
  histo1D["ttbar_mass_reco_WP"] = new TH1F("ttbar_mass_reco_WP","Reconstructed mass of the top quark pair (reco, wrong permutation); M_{t#bar{t}} [GeV]", 500, 0, 1000);
  histo1D["ttbar_mass_reco_UP"] = new TH1F("ttbar_mass_reco_UP","Reconstructed mass of the top quark pair (reco, no match); M_{t#bar{t}} [GeV]", 500, 0, 1000);
  
  // debug
  histo1D["debugLL_hadr_top_mass_reco"] = new TH1F("debugLL_hadr_top_mass_reco","Reconstructed mass of the hadronic top quark; M_{t} [GeV]", 400, 0, 800);
  histo1D["debugLL_minMlb_reco"]  = new TH1F("debugLL_minMlb_reco","Minimal reconstructed M_{lb} mass; min(M_{lb}) [GeV]", 400, 0, 800);
  histo1D["debugLL_ttbar_mass_reco"] = new TH1F("debugLL_ttbar_mass_reco","Reconstructed mass of the top quark pair; M_{t#bar{t}} [GeV]", 500, 0, 1000);
  histo1D["debugLL_dR_lep_b_lep_reco"]  = new TH1F("debugLL_dR_lep_b_lep_reco","Minimal delta R between the lepton and the leptonic b jet; #Delta R(l,b_{l})", 25, 0, 5);
  histo1D["debugLL_dR_lep_b_had_reco"]  = new TH1F("debugLL_dR_lep_b_had_reco","Minimal delta R between the lepton and the hadronic b jet; #Delta R(l,b_{h})", 25, 0, 5);
  
  /// KinFitter
  histo1D["KF_Chi2_TT"] = new TH1F("KF_Chi2_TT", "Chi2 value of kinFitter (TT); #chi^{2}", 200, 0, 20);
  histo1D["KF_W_mass_orig_TT"] = new TH1F("KF_W_mass_orig_TT", "W mass before kinFitter (TT); m_{W} [GeV]", 250, 0, 500);
  histo1D["KF_top_mass_orig_TT"] = new TH1F("KF_top_mass_orig_TT", "Top mass before kinFitter (TT); m_{t} [GeV]", 400, 0, 800);
  histo1D["KF_top_pt_orig_TT"] = new TH1F("KF_top_pt_orig_TT", "Top pt before kinFitter (TT); p_{T} [GeV]", 400, 0, 800);
  histo1D["KF_W_mass_corr_TT"] = new TH1F("KF_W_mass_corr_TT", "W mass after kinFitter (TT); m_{W,kf} [GeV]", 250, 0, 500);
  histo1D["KF_top_mass_corr_TT"] = new TH1F("KF_top_mass_corr_TT", "Top mass after kinFitter (TT); m_{t,kf} [GeV]", 400, 0, 800);
  histo1D["KF_top_pt_corr_TT"] = new TH1F("KF_top_pt_corr_TT", "Top pt after kinFitter (TT); p_{T,kf} [GeV]", 400, 0, 800);
  
  histo1D["KF_top_mass_orig_ex4j_chi2cut5_TT"] = new TH1F("KF_top_mass_orig_ex4j_chi2cut5_TT", "Top mass before kinFitter (TT) (exactly 4 jets - KF chi2 < 5); m_{t} [GeV]", 400, 0, 800);
  histo1D["KF_top_mass_corr_ex4j_chi2cut5_TT"] = new TH1F("KF_top_mass_corr_ex4j_chi2cut5_TT", "Top mass after kinFitter (TT) (exactly 4 jets - KF chi2 < 5); m_{t,kf} [GeV]", 400, 0, 800);
  
  histo1D["KF_top_mass_orig_ex4j_chi2cut2_TT"] = new TH1F("KF_top_mass_orig_ex4j_chi2cut2_TT", "Top mass before kinFitter (TT) (exactly 4 jets - KF chi2 < 2); m_{t} [GeV]", 400, 0, 800);
  histo1D["KF_top_mass_corr_ex4j_chi2cut2_TT"] = new TH1F("KF_top_mass_corr_ex4j_chi2cut2_TT", "Top mass after kinFitter (TT) (exactly 4 jets - KF chi2 < 2); m_{t,kf} [GeV]", 400, 0, 800);
  
  histo1D["KF_top_mass_orig_ex4j_chi2cut1p5_TT"] = new TH1F("KF_top_mass_orig_ex4j_chi2cut1p5_TT", "Top mass before kinFitter (TT) (exactly 4 jets - KF chi2 < 1.5); m_{t} [GeV]", 400, 0, 800);
  histo1D["KF_top_mass_corr_ex4j_chi2cut1p5_TT"] = new TH1F("KF_top_mass_corr_ex4j_chi2cut1p5_TT", "Top mass after kinFitter (TT) (exactly 4 jets - KF chi2 < 1.5); m_{t,kf} [GeV]", 400, 0, 800);
  
  histo1D["KF_jet0_Et_diff_TT"] = new TH1F("KF_jet0_Et_diff_TT", "Et difference after kinFitter for jet0 (TT); E_{T,kf} - E_{T,orig} [GeV]", 400, -50, 50);
  histo1D["KF_jet1_Et_diff_TT"] = new TH1F("KF_jet1_Et_diff_TT", "Et difference after kinFitter for jet1 (TT); E_{T,kf} - E_{T,orig} [GeV]", 400, -50, 50);
  
  histo1D["KF_top_mass_corr_CP"] = new TH1F("KF_top_mass_corr_CP", "Top mass after kinFitter for correct match (CP); m_{t,kf} [GeV]", 400, 0, 800);
  histo1D["KF_top_mass_corr_WP"] = new TH1F("KF_top_mass_corr_WP", "Top mass after kinFitter for wrong permutations (WP); m_{t,kf} [GeV]", 400, 0, 800);
  histo1D["KF_top_mass_corr_UP"] = new TH1F("KF_top_mass_corr_UP", "Top mass after kinFitter for no match (UP); m_{t,kf} [GeV]", 400, 0, 800);
  
  histo1D["KF_Chi2_CP"] = new TH1F("KF_Chi2_CP", "Chi2 value of kinFitter (CP); #chi^{2}", 200, 0, 20);
  histo1D["KF_Chi2_WP"] = new TH1F("KF_Chi2_WP", "Chi2 value of kinFitter (WP); #chi^{2}", 200, 0, 20);
  histo1D["KF_Chi2_UP"] = new TH1F("KF_Chi2_UP", "Chi2 value of kinFitter (UP); #chi^{2}", 200, 0, 20);
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
  InitHisto2DMatch();
  
  /// Reco
  histo2D["dR_lep_b_lep_vs_had_CP"] = new TH2F("dR_lep_b_lep_vs_had_CP","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, correct match); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_lep_b_lep_vs_had_WP"] = new TH2F("dR_lep_b_lep_vs_had_WP","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, wrong permutations); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_lep_b_lep_vs_had_UP"] = new TH2F("dR_lep_b_lep_vs_had_UP","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, no match); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  
  histo2D["ttbar_mass_vs_minMlb_CP"] = new TH2F("ttbar_mass_vs_minMlb_CP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_WP"] = new TH2F("ttbar_mass_vs_minMlb_WP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_UP"] = new TH2F("ttbar_mass_vs_minMlb_UP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  
  histo2D["ttbar_mass_vs_minMlb_dRlepCut_CP"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCut_CP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRlepCut_WP"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCut_WP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRlepCut_UP"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCut_UP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRhadCut_CP"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCut_CP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRhadCut_WP"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCut_WP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRhadCut_UP"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCut_UP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  
  histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_CP"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCutHard_CP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_WP"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCutHard_WP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_UP"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCutHard_UP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_CP"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCutHard_CP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_WP"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCutHard_WP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_UP"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCutHard_UP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  
  histo2D["ttbar_mass_vs_minMlb_dRBothCuts_CP"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCuts_CP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRBothCuts_WP"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCuts_WP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRBothCuts_UP"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCuts_UP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_CP"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCutsHard_CP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_WP"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCutsHard_WP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_UP"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCutsHard_UP","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  
  /// KinFitter
  histo2D["KF_W_mass_orig_vs_corr_TT"] = new TH2F("KF_W_mass_orig_vs_corr_TT", "W mass made with KF corrected jets vs. original jets (TT); m_{W,orig} [GeV]; m_{W,kf} [GeV]", 250, 0, 500, 250, 0, 500);
  histo2D["KF_top_mass_orig_vs_corr_TT"] = new TH2F("KF_top_mass_orig_vs_corr_TT", "Top mass made with KF corrected jets vs. original jets (TT); m_{t,orig} [GeV]; m_{t,kf} [GeV]", 400, 0, 800, 400, 0, 800);
  histo2D["KF_jets_Et_diff_TT"] = new TH2F("KF_jets_Et_diff_TT", "Et difference after kinFitter for jet1 in function of jet0 (TT); E_{T,0,kf} - E_{T,0,orig} [GeV]; E_{T,1,kf} - E_{T,1,orig} [GeV]", 400, -50, 50, 400, -50, 50);
  histo2D["KF_W_px_orig_vs_corr_TT"] = new TH2F("KF_W_px_orig_vs_corr_TT", "W p_{x} made with KF corrected jets vs. original jets (TT); p_{x,orig} [GeV]; p_{x,kf} [GeV]", 400, 0, 800, 400, 0, 800);
  histo2D["KF_W_py_orig_vs_corr_TT"] = new TH2F("KF_W_py_orig_vs_corr_TT", "W p_{y} made with KF corrected jets vs. original jets (TT); p_{y,orig} [GeV]; p_{y,kf} [GeV]", 400, 0, 800, 400, 0, 800);
  histo2D["KF_W_pz_orig_vs_corr_TT"] = new TH2F("KF_W_pz_orig_vs_corr_TT", "W p_{z} made with KF corrected jets vs. original jets (TT); p_{z,orig} [GeV]; p_{z,kf} [GeV]", 400, 0, 800, 400, 0, 800);
}

void InitHisto1DMatch()
{
  TH1::SetDefaultSumw2();
  
  histo1D["W_mass_reco_matched"] = new TH1F("W_mass_reco_matched","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 125, 0, 250);
  histo1D["top_mass_reco_matched"] = new TH1F("top_mass_reco_matched","Reconstructed top mass of matched events; M_{t} [GeV]", 175, 50, 400);
  histo1D["top_mass_gen_matched"] = new TH1F("top_mass_gen_matched","Generated top mass of matched events; M_{t} [GeV]", 1200, 100, 250);
  histo1D["mlb_matched_corr"]  = new TH1F("mlb_matched_corr","Reconstructed leptonic top mass using correctly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["mlb_matched_wrong"] = new TH1F("mlb_matched_wrong","Reconstructed leptonic top mass using wrongly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["ttbar_mass_matched_corr"] = new TH1F("ttbar_mass_matched_corr","Reconstructed mass of the top quark pair using correctly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong"] = new TH1F("ttbar_mass_matched_wrong","Reconstructed mass of the top quark pair using wrongly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["dR_lep_b_matched_corr"] = new TH1F("dR_lep_b_matched_corr","Delta R between the lepton and the leptonic b jet for matched events; #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_lep_b_matched_wrong"] = new TH1F("dR_lep_b_matched_wrong","Delta R between the lepton and the hadronic b jet for matched events; #Delta R(l,b_{had})", 25, 0, 5);
}

void InitHisto2DMatch()
{
  TH2::SetDefaultSumw2();
  
  /// Matched events
  histo2D["mlb_corr_mlb_wrong_matched"] = new TH2F("mlb_corr_mlb_wrong_matched","Wrongly constructed M_{lb} vs. correctly constructed M_{lb}; M_{lb_{lep}}; M_{lb_{had}}", 80, 0, 800, 80, 0, 800);
  histo2D["dR_lep_b_corr_dR_lep_b_wrong_matched"] = new TH2F("dR_lep_b_corr_dR_lep_b_wrong_matched","Wrongly constructed dR(l,b) vs. correctly constructed dR(l,b); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["mlb_dR_lep_b_corr_matched"] = new TH2F("mlb_dR_lep_b_corr_matched","dR(l,b) vs. M_{lb}; M_{lb_{lep}}; #Delta R(l,b_{lep})", 80, 0, 800, 25, 0, 5);
  histo2D["mlb_dR_lep_b_wrong_matched"] = new TH2F("mlb_dR_lep_b_wrong_matched","dR(l,b) vs. M_{lb}, both wrongly matched; M_{lb_{had}}; #Delta R(l,b_{had})", 80, 0, 800, 25, 0, 5);
}

void TruthMatching(vector<TLorentzVector> partons, vector<TLorentzVector> selectedJets, pair<unsigned int, unsigned int> *MCPermutation)  /// MCPermutation: 0,1 hadronic W jet; 2 hadronic b jet; 3 leptonic b jet
{
  JetPartonMatching matching = JetPartonMatching(partons, selectedJets, 2, true, true, 0.3);  // partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
  
  if (matching.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matching.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matching.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
  }

  if (verbose > 2)
    cout << "Matching done" << endl;
  
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int j = JetPartonPair[i].second;
    
    if ( fabs(mc_pdgId[partonId[j]]) < 6 )  //light/b quarks, 6 should stay hardcoded
    {
      if ( ( muPlusFromTop && mc_mother[partonId[j]] == -24 && mc_granny[partonId[j]] == -pdgID_top )
          || ( muMinusFromTop && mc_mother[partonId[j]] == 24 && mc_granny[partonId[j]] == pdgID_top ) )  // if mu+, check if mother of particle is W- and granny tbar --> then it is a quark from W- decay
      {
        if (verbose > 3)
          cout << "Light jet: " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;
        
        if (MCPermutation[0].first == 9999)
        {
          MCPermutation[0] = JetPartonPair[i];
        }
        else if (MCPermutation[1].first == 9999)
        {
          MCPermutation[1] = JetPartonPair[i];
        }
        else
        {
          cerr << "Found a third jet coming from a W boson which comes from a top quark..." << endl;
          cerr << " -- muMinusFromMtop: " << muMinusFromTop << " muPlusFromMtop: " << muPlusFromTop << endl;
          cerr << " -- pdgId: " << mc_pdgId[partonId[j]] << " mother: " << mc_mother[partonId[j]] << " granny: " << mc_granny[partonId[j]] << " Pt: " << mc_pt[partonId[j]] << endl;
//          cerr << " -- ievt: " << ievt << endl;
          exit(1);
        }
      }
    }
    if ( fabs(mc_pdgId[partonId[j]]) == 5 )
    {
      if ( ( muPlusFromTop && mc_mother[partonId[j]] == -pdgID_top )
          || ( muMinusFromTop && mc_mother[partonId[j]] == pdgID_top ) )  // if mu+ (top decay leptonic) and mother is antitop ---> hadronic b
      {
        if (verbose > 3)
          cout << "b jet:     " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;
        
        MCPermutation[2] = JetPartonPair[i];
      }
      else if ( ( muPlusFromTop && mc_mother[partonId[j]] == pdgID_top )
        || ( muMinusFromTop && mc_mother[partonId[j]] == -pdgID_top ) )  // if mu+ (top decay leptonic) and mother is top ---> leptonic b
      {
        if (verbose > 3)
          cout << "b jet:     " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;
        
        MCPermutation[3] = JetPartonPair[i];
      }
    }
  }  /// End loop over Jet Parton Pairs
  
  
  if ( MCPermutation[0].first != 9999 && MCPermutation[1].first != 9999 && MCPermutation[2].first != 9999 && MCPermutation[3].first != 9999 )
  {
    all4PartonsMatched = true;
    nofMatchedEvents++;
    if (MCPermutation[0].first < 4 && MCPermutation[1].first < 4 && MCPermutation[2].first < 4 && MCPermutation[3].first < 4)
      all4JetsMatched_MCdef_ = true;
  }
  else if (verbose > 3) cout << "Size JetPartonPair: " << JetPartonPair.size() << ". Not all partons matched!" << endl;
  
  if ( MCPermutation[0].first != 9999 && MCPermutation[1].first != 9999 && MCPermutation[2].first != 9999 )
  {
    hadronicTopJetsMatched = true;
    nofHadrMatchedEvents++;
    if ( MCPermutation[0].first < 4 && MCPermutation[1].first < 4 && MCPermutation[2].first < 4 )
      hadronicTopJetsMatched_MCdef_ = true;
  }
  
  if ( genmuon != -9999 && ROOT::Math::VectorUtil::DeltaR(mcParticles[genmuon], selectedLepton[0]) < 0.1 )
    muonmatched = true;
}

void ClearMetaData()
{
  nEvents = 0;
  nEventsSel = 0;
  for (int i = 0; i < 10; i++)
  {
    cutFlow[i] = 0;
  }
  appliedJER = 999;
  appliedJES = 999;
  appliedPU = 999;
  if (isData)
  {
    nofEventsHLTv2 = 0;
    nofEventsHLTv3 = 0;
    nofSelEventsHLTv2 = 0;
    nofSelEventsHLTv3 = 0;
  }
  
  strSyst = "";
  eqLumi = 1.;
  
  nofHardSelected = 0;
  nofMatchedEvents = 0;
  nofHadrMatchedEvents = 0;
  nofCorrectlyMatched = 0;
  nofNotCorrectlyMatched = 0;
  nofAcceptedKFit = 0;
  
  toyMax = 1.;
}

void ClearLeaves()
{
  run_num = -1;
  evt_num = -1;
  lumi_num = -1;
  nvtx = -1;
  npu = -1;
  rho = -999.;
  isTrigged = false;
  for (Int_t i = 0; i < 10; i++)
    cutFlow[i] = 0;
  hasExactly4Jets = false;
  hasJetLeptonCleaning = false;
  nLeptons = -1;
  nMuons = -1;
  muon_charge[0] = 0;
  muon_pt[0] = 0.;
  muon_phi[0] = 0.;
  muon_eta[0] = 0.;
  muon_E[0] = 0.;
  muon_M[0] = 0.;
  muon_d0[0] = 999.;
  muon_chargedHadronIso[0] = 999.;
  muon_neutralHadronIso[0] = 999.;
  muon_photonIso[0] = 999.;
  muon_puChargedHadronIso[0] = 999.;
  muon_relIso[0] = 999.;
  muon_pfIso[0] = 999.;
  nJets = -1;
  for (Int_t i = 0; i < 20; i++)
  {
    jet_charge[i] = 0;
    jet_pt[i] = 0.;
    jet_phi[i] = 0.;
    jet_eta[i] = 0.;
    jet_E[i] = 0.;
    jet_M[i] = 0.;
    jet_bdiscr[i] = -1.;
  }
  met_pt = 0.;
  met_phi = 0.;
  met_eta = 0.;
  met_Et = 0.;
  met_E = 0.;
  if (! isData)
  {
    nMCParticles = -1;
    for (Int_t i = 0; i < 200; i++)
    {
      mc_status[i] = -1;
      mc_pdgId[i] = 0;
      mc_mother[i] = 0;
      mc_granny[i] = 0;
      mc_pt[i] = 0.;
      mc_phi[i] = 0.;
      mc_eta[i] = 0.;
      mc_E[i] = 0.;
      mc_M[i] = 0.;
    }
  }
  nloWeight = 1.;
  puSF = 1.;
  btagSF = 1.;
  muonIdSF[0] = 1.;
  muonIsoSF[0] = 1.;
  muonTrigSFv2[0] = 1.;
  muonTrigSFv3[0] = 1.;
}

void ClearTLVs()
{
  muon.Clear();
  jet.Clear();
  mcpart.Clear();
  WCandidate.Clear();
  selectedLepton.clear();
  selectedJets.clear();
  selectedBJets.clear();
  selectedJetsKFcorrected.clear();
  mcParticles.clear();
  partons.clear();
  partonsMatched.clear();
  jetsMatched.clear();
}

void ClearMatching()
{
  partons.clear();
  partonsMatched.clear();
  jetsMatched.clear();
  
  doMatching = true;
  all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
  all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
  hadronicTopJetsMatched = false;
  hadronicTopJetsMatched_MCdef_ = false;
  for (int i = 0; i < 4; i++)
  {
    MCPermutation[i] = pair<unsigned int,unsigned int>(9999,9999);
  }
  topQuark = -9999;
  antiTopQuark = -9999;
  genmuon = -9999;
  muonmatched = false;
  muPlusFromTop = false;
  muMinusFromTop = false;
  partonId.clear();
  
  matchedWMass_reco = -1.;
  matchedTopMass_reco = -1.;
  matchedTopMass_gen = -1.;
  matchedMlb_corr = -1.;
  matchedMlb_wrong = -1.;
  matchedTtbarMass_corr = -1.;
  matchedTtbarMass_wrong = -1.;
  matchedDRLepB_corr = 999.;
  matchedDRLepB_wrong = 999.;
}

void ClearVars()
{
  ClearMatching();
  
  scaleFactor = 1.;
  widthSF = 1.;
  bJetId.clear();
  selectedJetsCharge.clear();
  selectedJetsBDiscr.clear();
  bdiscrTop = -99.;
  bdiscrTop2 = -99.;
  tempbdiscr = -99.;
  labelB1 = -9999;
  labelB2 = -9999;
  for (int i = 0; i < 4; i++)
  {
    labelsReco[i] = -9999;
  }
  reco_minMlb = 9999.;
  reco_minMl_nonb = 9999.;
  reco_ttbarMass = -1.;
  reco_dRLepB_lep = -1.;
  reco_dRLepB_had = -1.;
  min_Mlb = 9999.;
  dRLepB = -1.;
  labelMlb = -9999;
  labelMl_nonb = -9999;
  massForWidth = 0.01;
  kFitVerbosity = false;
  kFitChi2 = 99.;
  Wmass_reco_orig = -1.;
  Wmass_reco_kf = -1.;
  topmass_reco_orig = -1.;
  topmass_reco_kf = -1.;
  toy = -1.;
}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearVars();
}

void FillGeneralPlots(int d)
{
  double M3 = (selectedJets[0] + selectedJets[1] + selectedJets[2]).M();
  double Ht = selectedJets[0].Pt() + selectedJets[1].Pt() + selectedJets[2].Pt() + selectedJets[3].Pt();
  
  double Mlb_temp = -1;
  for (unsigned int i = 0; i < selectedBJets.size(); i++)
  {
    Mlb_temp = (selectedLepton[0] + selectedBJets[i]).M();
    if ( Mlb_temp < min_Mlb )
    {
      min_Mlb = Mlb_temp;
      dRLepB = ROOT::Math::VectorUtil::DeltaR( selectedBJets[i], selectedLepton[0]);
    }
  }
  
  
  MSPlot["muon_pT"]->Fill(selectedLepton[0].Pt(), datasets[d], true, Luminosity*scaleFactor);
  MSPlot["muon_eta"]->Fill(selectedLepton[0].Eta(), datasets[d], true, Luminosity*scaleFactor);
  MSPlot["muon_phi"]->Fill(selectedLepton[0].Phi(), datasets[d], true, Luminosity*scaleFactor);
  MSPlot["muon_relIso"]->Fill(muon_relIso[0], datasets[d], true, Luminosity*scaleFactor);
  MSPlot["muon_d0"]->Fill(muon_d0[0], datasets[d], true, Luminosity*scaleFactor);
  MSPlot["leadingJet_pT"]->Fill(selectedJets[0].Pt(), datasets[d], true, Luminosity*scaleFactor);
  MSPlot["jet2_pT"]->Fill(selectedJets[1].Pt(), datasets[d], true, Luminosity*scaleFactor);
  MSPlot["jet3_pT"]->Fill(selectedJets[2].Pt(), datasets[d], true, Luminosity*scaleFactor);
  MSPlot["jet4_pT"]->Fill(selectedJets[3].Pt(), datasets[d], true, Luminosity*scaleFactor);
  MSPlot["Ht_4leadingJets"]->Fill(Ht, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["met_pT"]->Fill(met_pt, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["met_eta"]->Fill(met_eta, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["met_phi"]->Fill(met_phi, datasets[d], true, Luminosity*scaleFactor);
  
  MSPlot["M3"]->Fill(M3, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["min_Mlb"]->Fill(min_Mlb, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["dR_Lep_B"]->Fill(dRLepB, datasets[d], true, Luminosity*scaleFactor);
  
  MSPlot["nJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
  MSPlot["nBJets"]->Fill(selectedBJets.size(), datasets[d], true, Luminosity*scaleFactor);
  
  MSPlot["CSVv2Discr_leadingJet"]->Fill(jet_bdiscr[0], datasets[d], true, Luminosity*scaleFactor);
  MSPlot["CSVv2Discr_jet2"]->Fill(jet_bdiscr[1], datasets[d], true, Luminosity*scaleFactor);
  MSPlot["CSVv2Discr_jet3"]->Fill(jet_bdiscr[2], datasets[d], true, Luminosity*scaleFactor);
  MSPlot["CSVv2Discr_jet4"]->Fill(jet_bdiscr[3], datasets[d], true, Luminosity*scaleFactor);
  
  int labelB = -1;
  double highestBDiscr = -999.;
  for (int iJet = 0; iJet < selectedJets.size(); iJet++)
  {
    MSPlot["CSVv2Discr_allJets"]->Fill(jet_bdiscr[iJet], datasets[d], true, Luminosity*scaleFactor);
    if ( jet_bdiscr[iJet] > highestBDiscr )
    {
      highestBDiscr = jet_bdiscr[iJet];
      labelB = iJet;
    }
  }
  MSPlot["CSVv2Discr_highest"]->Fill(highestBDiscr, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["CSVv2Discr_jetNb"]->Fill(labelB, datasets[d], true, Luminosity*scaleFactor);
}

void FillMatchingPlots()
{
  if (hadronicTopJetsMatched)
  {  
    matchedWMass_reco = (jetsMatched[0] + jetsMatched[1]).M();
    matchedTopMass_reco = (jetsMatched[0] + jetsMatched[1] + jetsMatched[2]).M();
    matchedTopMass_gen = (partonsMatched[0] + partonsMatched[1] + partonsMatched[2]).M();
    
    histo1D["W_mass_reco_matched"]->Fill(matchedWMass_reco, widthSF);
    histo1D["top_mass_reco_matched"]->Fill(matchedTopMass_reco, widthSF);
    histo1D["top_mass_gen_matched"]->Fill(matchedTopMass_gen, widthSF);
    
    histo1D["mTop_div_aveMTop_TT_matched_jets"]->Fill(matchedTopMass_reco/aveTopMass[0]);
    
    if (all4PartonsMatched && muonmatched)
    {
      matchedMlb_corr = (selectedLepton[0] + jetsMatched[3]).M();  // lept b
      matchedMlb_wrong = (selectedLepton[0] + jetsMatched[2]).M();  // hadr b
      matchedTtbarMass_corr = matchedTopMass_reco + matchedMlb_corr;
      matchedTtbarMass_wrong = matchedTopMass_reco + matchedMlb_wrong;
      matchedDRLepB_corr = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[3]);  // lept b
      matchedDRLepB_wrong = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[2]);  // hadr b
      
      histo1D["mlb_matched_corr"]->Fill(matchedMlb_corr);
      histo1D["mlb_matched_wrong"]->Fill(matchedMlb_wrong);
      histo1D["ttbar_mass_matched_corr"]->Fill(matchedTtbarMass_corr);
      histo1D["ttbar_mass_matched_wrong"]->Fill(matchedTtbarMass_wrong);
      histo1D["dR_lep_b_matched_corr"]->Fill(matchedDRLepB_corr);
      histo1D["dR_lep_b_matched_wrong"]->Fill(matchedDRLepB_wrong);
      
      histo2D["mlb_corr_mlb_wrong_matched"]->Fill(matchedMlb_corr, matchedMlb_wrong);
      histo2D["dR_lep_b_corr_dR_lep_b_wrong_matched"]->Fill(matchedDRLepB_corr, matchedDRLepB_wrong);
      histo2D["mlb_dR_lep_b_corr_matched"]->Fill(matchedMlb_corr, matchedDRLepB_corr);
      histo2D["mlb_dR_lep_b_wrong_matched"]->Fill(matchedMlb_wrong, matchedDRLepB_wrong);
    }  // end all4PartonsMatched && muonMatched
  }  // end hadronicTopJetsMatched
}

void FillKinFitPlots(int d)
{
  MSPlot["KF_W_mass_orig"]->Fill(Wmass_reco_orig, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["KF_W_mass_corr"]->Fill(Wmass_reco_kf, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["KF_top_mass_orig"]->Fill(topmass_reco_orig, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["KF_top_mass_corr"]->Fill(topmass_reco_kf, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["KF_top_pt_orig"]->Fill(toppt_reco_orig, datasets[d], true, Luminosity*scaleFactor);
  MSPlot["KF_top_pt_corr"]->Fill(toppt_reco_kf, datasets[d], true, Luminosity*scaleFactor);
  
  if (isTTbar)
  { 
    histo1D["KF_W_mass_orig_TT"]->Fill(Wmass_reco_orig);
    histo1D["KF_W_mass_corr_TT"]->Fill(Wmass_reco_kf);
    histo2D["KF_W_mass_orig_vs_corr_TT"]->Fill(Wmass_reco_orig, Wmass_reco_kf);
    histo2D["KF_W_px_orig_vs_corr_TT"]->Fill( (selectedJets[0] + selectedJets[1]).Px(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Px() );
    histo2D["KF_W_py_orig_vs_corr_TT"]->Fill( (selectedJets[0] + selectedJets[1]).Py(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Py() );
    histo2D["KF_W_pz_orig_vs_corr_TT"]->Fill( (selectedJets[0] + selectedJets[1]).Pz(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Pz() );
    
    histo1D["KF_top_mass_orig_TT"]->Fill(topmass_reco_orig);
    histo1D["KF_top_mass_corr_TT"]->Fill(topmass_reco_kf);
    histo2D["KF_top_mass_orig_vs_corr_TT"]->Fill(topmass_reco_orig, topmass_reco_kf);
    histo1D["KF_top_pt_orig_TT"]->Fill(toppt_reco_orig);
    histo1D["KF_top_pt_corr_TT"]->Fill(toppt_reco_kf);
    
    histo1D["KF_jet0_Et_diff_TT"]->Fill((selectedJetsKFcorrected[0] - selectedJets[labelsReco[0]]).Et());
    histo1D["KF_jet1_Et_diff_TT"]->Fill((selectedJetsKFcorrected[1] - selectedJets[labelsReco[1]]).Et());
    histo2D["KF_jets_Et_diff_TT"]->Fill((selectedJetsKFcorrected[0] - selectedJets[labelsReco[0]]).Et(), (selectedJetsKFcorrected[1] - selectedJets[labelsReco[1]]).Et());
    
    histo1D["KF_Chi2_TT"]->Fill(kFitChi2);
    
    if ( kFitChi2 < 5. )
    {
      histo1D["KF_top_mass_orig_ex4j_chi2cut5_TT"]->Fill(topmass_reco_orig);
      histo1D["KF_top_mass_corr_ex4j_chi2cut5_TT"]->Fill(topmass_reco_kf);
      
      if ( kFitChi2 < 2. )
      {
        histo1D["KF_top_mass_orig_ex4j_chi2cut2_TT"]->Fill(topmass_reco_orig);
        histo1D["KF_top_mass_corr_ex4j_chi2cut2_TT"]->Fill(topmass_reco_kf);
        
        if ( kFitChi2 < 1.5 )
        {
          histo1D["KF_top_mass_orig_ex4j_chi2cut1p5_TT"]->Fill(topmass_reco_orig);
          histo1D["KF_top_mass_corr_ex4j_chi2cut1p5_TT"]->Fill(topmass_reco_kf);
        }  // 1.5
      }  // 2
    }  // 5
  }
}

long GetNEvents(TTree* fChain, string var, bool isData)
{
  GetMetaData(fChain, isData);
  long varNew = 0;
  for (unsigned int iEntry = 0; iEntry < fChain->GetEntries(); iEntry++)
  {
    fChain->GetEntry(iEntry);
    varNew += (fChain->FindLeaf(var.c_str()))->GetValueLong64();
  }
  
  return varNew;
}

void GetHLTFraction(double* fractions)
{
  TFile* fileHLT = new TFile((pathNtuples+"Ntuples_data.root").c_str(),"READ");
  TTree* fChain = (TTree*) fileHLT->Get("stats");
  
  //fChain->SetBranchAddress("nofEventsHLTv2", &nofEventsHLTv2, &b_nofEventsHLTv2);
  //fChain->SetBranchAddress("nofEventsHLTv3", &nofEventsHLTv3, &b_nofEventsHLTv3);
  
  long nv2 = GetNEvents(fChain, "nofEventsHLTv2", 1);
  long nv3 = GetNEvents(fChain, "nofEventsHLTv3", 1);
  
  fractions[0] = ((double)nv2) / ((double)nv2 + (double)nv3);
  fractions[1] = ((double)nv3) / ((double)nv2 + (double)nv3);
  
  fileHLT->Close();
  ClearMetaData();
}

double BreitWigner(double topMass, double scale)
{
  double BWmass = genTopMass;
  double BWgamma = scale*genTopWidth;
  double gammaterm = sqrt( pow(BWmass, 4) + pow(BWmass*BWgamma, 2) );
  double numerator = 2*sqrt(2)*BWmass*BWgamma*gammaterm/( TMath::Pi()*sqrt( pow(BWmass, 2) + gammaterm ) );
  double denominator = pow(pow(topMass, 2) - pow(BWmass, 2), 2) + pow(BWmass*BWgamma, 2);
  
  return numerator/denominator;
}

double eventWeightCalculator(double topMass, double scale)
{
  return BreitWigner(topMass, scale)/BreitWigner(topMass, 1.);
}

Double_t voigt(Double_t *x, Double_t *par) {
  // Computation of Voigt function (normalised).
  // Voigt is a convolution of
  // gauss(xx) = 1/(sqrt(2*pi)*sigma) * exp(xx*xx/(2*sigma*sigma)
  // and
  // lorentz(xx) = (1/pi) * (lg/2) / (xx*xx + g*g/4)
  // functions.
  //
  // To calculate the Faddeeva function with relative error less than 10^(-r).
  // r can be set by the the user subject to the constraints 2 <= r <= 5.
  
  /// Voigt(x-mu, sigma, lg, r)
  return TMath::Voigt(x[0]-mu_CP, sigma_CP, par[0], r_CP)*norm_CP;
}

Double_t crysBall_WP(Double_t *x) {
  // params: alpha, n, sigma, mu
  if ( sigma_WP <= 0. ) return 0.;
  Double_t alpha = fabs(alpha_WP);
  Double_t A = pow( n_WP/alpha , n_WP) * exp(-alpha*alpha/2.);
  Double_t B = n_WP/alpha - alpha;
  
  Double_t ref = (x[0] - mu_WP)/sigma_WP;  // (x-mean)/sigma
  if ( alpha_WP < 0. ) ref = -ref;
  Double_t fitfunc = 1.;
  if ( ref > -alpha ) fitfunc = fitfunc * exp(-ref*ref/2.);
  else if (ref <= -alpha ) fitfunc = fitfunc * A * pow ( B - ref , -n_WP);
  return fitfunc*norm_WP;
}

Double_t crysBall_UP(Double_t *x) {
  // params: alpha, n, sigma, mu
  if ( sigma_UP <= 0. ) return 0.;
  Double_t alpha = fabs(alpha_UP);
  Double_t A = pow( n_UP/alpha , n_UP) * exp(-alpha*alpha/2.);
  Double_t B = n_UP/alpha - alpha;
  
  Double_t ref = (x[0] - mu_UP)/sigma_UP;  // (x-mean)/sigma
  if ( alpha_UP < 0. ) ref = -ref;
  Double_t fitfunc = 1.;
  if ( ref > -alpha ) fitfunc = fitfunc * exp(-ref*ref/2.);
  else if (ref <= -alpha ) fitfunc = fitfunc * A * pow ( B - ref , -n_UP);
  return fitfunc*norm_UP;
}

Double_t crysBall_WPUP(Double_t *x) {
  // params: alpha, n, sigma, mu
  if ( sigma_WPUP <= 0. ) return 0.;
  Double_t alpha = fabs(alpha_WPUP);
  Double_t A = pow( n_WPUP/alpha , n_WPUP) * exp(-alpha*alpha/2.);
  Double_t B = n_WPUP/alpha - alpha;
  
  Double_t ref = (x[0] - mu_WPUP)/sigma_WPUP;  // (x-mean)/sigma
  if ( alpha_WPUP < 0. ) ref = -ref;
  Double_t fitfunc = 1.;
  if ( ref > -alpha ) fitfunc = fitfunc * exp(-ref*ref/2.);
  else if (ref <= -alpha ) fitfunc = fitfunc * A * pow ( B - ref , -n_WPUP);
  return fitfunc*norm_WPUP;
}

Double_t logLikelihood(Double_t *x, Double_t *par) {
  return -TMath::Log( norm_comb*( f_CP*voigt(x, par) + f_WP*crysBall_WP(x) + f_UP*crysBall_UP(x) ) );
}

Double_t logLikelihood2(Double_t *x, Double_t *par) {
  return -TMath::Log( norm_comb*( f_CP*voigt(x, par) + (f_WP+f_UP)*crysBall_WPUP(x) ) );
}

Double_t fakeLikelihood(Double_t *x, Double_t *par) {
  return -TMath::Log( voigt(x, par) );
}

Double_t widthToGammaTranslation(Double_t *x) {
  return gammaConvConst + gammaConvRico * x[0];
}
