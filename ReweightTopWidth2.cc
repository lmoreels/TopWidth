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
#include <TGraph.h>
#include <TGraph2D.h>

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
#include "Tools/interface/EventReweighting.h"
#include "Tools/interface/KinFitter.h"
#include "Tools/interface/Likelihood.h"


using namespace std;
using namespace TopTree;


bool test = false;
bool testHistos = false;
bool testTTbarOnly = false;
bool doGenOnly = true;
bool makePlots = true;
bool makeReweightedPlots = true; 
bool makeMatchingPlots = false;
bool makeRecoPlots = false;
bool calculateResolutionFunctions = false;
bool calculateAverageMass = false;
bool doKinFit = true;
bool applyKinFitCut = true;
double kinFitCutValue = 5.;

bool doMETCleaning = true;
bool applyLeptonSF = true;
bool applyLeptonSFup = false;
bool applyLeptonSFdown = false;
bool applyPU = true;
bool applyPUup = false;
bool applyPUdown = false;
bool applyBTagSF = true;
bool applyBTagSFup = false;
bool applyBTagSFdown = false;
bool applyTopPtReweighting = false;
//bool applyJER = true;
//bool applyJEC = true;

bool doReweighting = false;
bool applyWidthSF = false;
double scaleWidth = 1.;

string systStr = "nominal";
pair<string,string> whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos )
  {
    return pair<string,string>("170607","170522");
  }
  else if ( syst.find("JECup") != std::string::npos ) return pair<string,string>("170602","170522");
  else if ( syst.find("JECdown") != std::string::npos ) return pair<string,string>("170606","170522");
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    return pair<string,string>("170607","170522");
  }
}
pair<string,string> ntupleDate = whichDate(systStr);
int verbose = 2;

string pathNtuples = "";
string pathNtuplesMC = "";
string pathNtuplesData = "";
string outputDirLL = "LikelihoodTemplates/";
string inputDirLL = "";
string inputDateLL = "170613_1657/";  // TT nominal
bool isData = false;
bool isTTbar = false;

int nofHardSelected = 0;
int nofMETCleaned = 0;
int nofMatchedEvents = 0;
int nofHadrMatchedEvents = 0;
int nofHadrMatchedEventsAKF = 0;
int nofCorrectlyMatched = 0;
int nofNotCorrectlyMatched = 0;
int nofCorrectlyMatchedAKF = 0;
int nofNotCorrectlyMatchedAKF = 0;
int nofCorrectlyMatchedAKFNoCut = 0;
int nofNotCorrectlyMatchedAKFNoCut = 0;
int nofNoMatchAKFNoCut = 0;
int nofCM = 0, nofWM = 0, nofNM = 0;
int nofCM_TT = 0, nofWM_TT = 0, nofNM_TT = 0;
double nofCMl = 0, nofWMl = 0, nofNMl = 0;
double nofCM_weighted = 0, nofWM_weighted = 0, nofNM_weighted = 0;

/// Lumi per data era
double lumi_runBCDEF = 19.67550334113;  // 1/fb
double lumi_runGH = 16.146177597883;  // 1/fb
double Luminosity = (lumi_runBCDEF + lumi_runGH)*1000;  // 1/pb

///  Working points for b tagging  // Updated 13/04/17, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
double CSVv2Loose  = 0.5426;
double CSVv2Medium = 0.8484;
double CSVv2Tight  = 0.9535;

/// Average top mass
// TT_genp_match, TT_genj_match, TT_reco_match, TT_reco_wrongMatch_WM/NM, TT_reco_noMatch, TT_reco_wrongPerm, TT_reco, ST_t_top, ST_t_antitop, ST_tW_top, ST_tW_antitop, DYJets, WJets, data, Reco, All, MC, Reco, All, Samples
// also background in CM/WM/NM cats (unlike name suggests)
const int nofAveMasses = 16;
//  KF chi2 < 5
//std::array<double, 14> aveTopMass = {171.826, 169.746, 167.556, 197.087, 196.662, 198.143, 182.150, 249.229, 246.893, 226.933, 225.681, 185.024, 184.880, 184.902};  // no DYJets, no WJets // Res 170608
std::array<double, 14> aveTopMass = {171.826, 169.746, 167.511, 197.053, 196.688, 197.911, 181.895, 249.468, 247.437, 227.529, 226.099, 184.794, 184.594, 184.62};  // Res 170608 Single Gaus
//std::array<double, 14> aveTopMass = {171.826, 169.746, 167.572, 196.603, 196.072, 197.919, 181.953, 247.003, 243.879, 226.505, 224.951, 184.717, 184.598, 184.616};  // Res 170515
//  no KF chi2 cut
//std::array<double, nofAveMasses> aveTopMass = {171.810, 168.728, 167.110, 203.721, 204.952, 198.233, 193.403, 270.895, 267.167, 230.144, 229.649, 250.010, 242.091, 200.455, 193.963, 194.025};

/// # events after kin fitter
//  KF chi2 < 5
int nEventsAKF[] = {100949, 544305, 8927, 8847, 7051, 4283, 3, 16, 55, 276, 1, 18, 123, 1415};

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TGraph*> graph;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;

vector < Dataset* > datasets;

ofstream txtMassGenPMatched, txtMassGenJMatched, txtMassRecoCM, txtMassRecoWMNM, txtMassRecoNM, txtMassRecoWM, txtMassReco;

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


string ConvertDoubleToString(double Number);
string DotReplace(double var);
string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();
bool fexists(const char *filename);
void GetMetaData(TTree* tree, bool isData);
void InitTree(TTree* tree, bool isData);
void InitHisto1D();
void InitHisto1DGen();
void InitHisto1DReweighted();
void InitHisto1DMatch();
void InitHisto1DReco();
void InitHisto2D();
void InitHisto2DGen();
void InitHisto2DReweighted();
void InitHisto2DMatch();
void TruthMatching(vector<TLorentzVector> partons, vector<TLorentzVector> selectedJets, pair<unsigned int, unsigned int> *MCPermutation);
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearVars();
void ClearObjects();
void FillControlPlots(vector<Dataset *> datasets, int d);
void FillControlPlots(vector<Dataset *> datasets, int d, string suffix);
void FillMatchingPlots();
void FillKinFitPlots(bool doneKinFit);
void FillCatsPlots(string catSuffix);
long GetNEvents(TTree* fChain, string var, bool isData);
void GetEraFraction(double* fractions);
void PrintKFDebug(int ievt);



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
Bool_t          hasErasedBadOrCloneMuon;
Bool_t          filterHBHENoise;
Bool_t          filterHBHEIso;
Bool_t          filterCSCTightHalo;
Bool_t          filterEcalDeadCell;
Bool_t          filterEEBadSc;
Bool_t          filterBadChCand;
Bool_t          filterBadMuon;
Bool_t          passedMETFilter;
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
Int_t           jet_nConstituents[20];   //[nJets]
Int_t           jet_nChConstituents[20];   //[nJets]
Int_t           jet_charge[20];   //[nJets]
Double_t        jet_pt[20];   //[nJets]
Double_t        jet_phi[20];   //[nJets]
Double_t        jet_eta[20];   //[nJets]
Double_t        jet_E[20];   //[nJets]
Double_t        jet_M[20];   //[nJets]
Double_t        jet_bdiscr[20];   //[nJets]
Double_t        met_px;
Double_t        met_py;
Double_t        met_pt;
Double_t        met_phi;
Double_t        met_eta;
Double_t        met_Et;
Double_t        met_E;
Double_t        met_corr_px;
Double_t        met_corr_py;
Double_t        met_corr_pt;
Double_t        met_corr_phi;
Double_t        met_corr_eta;
Double_t        met_corr_Et;
Double_t        met_corr_E;
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
Bool_t          hasGenTop;
Bool_t          hasGenTopWithStatus22;
Bool_t          hasGenTopWithStatus62;
Bool_t          hasGenAntiTop;
Bool_t          hasGenAntiTopWithStatus22;
Bool_t          hasGenAntiTopWithStatus62;

Long64_t        nEvents;
Long64_t        nEventsSel;
Long64_t        nofEventsWithGenTop;
Long64_t        nofEventsWithGenTopWithStatus22or62;
Long64_t        nofEventsWithGenAntiTop;
Long64_t        nofEventsWithGenAntiTopWithStatus22or62;
Long64_t        nofTTEventsWithoutBothGenTops;
Long64_t        nofTTEventsWithoutGenTop;
Long64_t        nofTTEventsWithoutGenAntiTop;
Long64_t        nofTTEventsWithoutBothGenTopsWithStatus22;
Long64_t        nofTTEventsWithoutGenTopWithStatus22;
Long64_t        nofTTEventsWithoutGenAntiTopWithStatus22;
Long64_t        nofTTEventsWithoutBothGenTopsWithStatus62;
Long64_t        nofTTEventsWithoutGenTopWithStatus62;
Long64_t        nofTTEventsWithoutGenAntiTopWithStatus62;

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
TBranch        *b_hasErasedBadOrCloneMuon;   //!
TBranch        *b_filterHBHENoise;   //!
TBranch        *b_filterHBHEIso;   //!
TBranch        *b_filterCSCTightHalo;   //!
TBranch        *b_filterEcalDeadCell;   //!
TBranch        *b_filterEEBadSc;   //!
TBranch        *b_filterBadChCand;   //!
TBranch        *b_filterBadMuon;   //!
TBranch        *b_passedMETFilter;   //!
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
TBranch        *b_jet_nConstituents;   //!
TBranch        *b_jet_nChConstituents;   //!
TBranch        *b_jet_charge;   //!
TBranch        *b_jet_pt;   //!
TBranch        *b_jet_phi;   //!
TBranch        *b_jet_eta;   //!
TBranch        *b_jet_E;   //!
TBranch        *b_jet_M;   //!
TBranch        *b_jet_bdiscr;   //!
TBranch        *b_met_px;   //!
TBranch        *b_met_py;   //!
TBranch        *b_met_pt;   //!
TBranch        *b_met_phi;   //!
TBranch        *b_met_eta;   //!
TBranch        *b_met_Et;   //!
TBranch        *b_met_E;   //!
TBranch        *b_met_corr_px;   //!
TBranch        *b_met_corr_py;   //!
TBranch        *b_met_corr_pt;   //!
TBranch        *b_met_corr_phi;   //!
TBranch        *b_met_corr_eta;   //!
TBranch        *b_met_corr_Et;   //!
TBranch        *b_met_corr_E;   //!
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
TBranch        *b_hasGenTop;   //!
TBranch        *b_hasGenTopWithStatus22;   //!
TBranch        *b_hasGenTopWithStatus62;   //!
TBranch        *b_hasGenAntiTop;   //!
TBranch        *b_hasGenAntiTopWithStatus22;   //!
TBranch        *b_hasGenAntiTopWithStatus62;   //!

TBranch        *b_nEvents;   //!
TBranch        *b_nEventsSel;   //!
TBranch        *b_nofEventsWithGenTop;   //!
TBranch        *b_nofEventsWithGenTopWithStatus22or62;   //!
TBranch        *b_nofEventsWithGenAntiTop;   //!
TBranch        *b_nofEventsWithGenAntiTopWithStatus22or62;   //!
TBranch        *b_nofTTEventsWithoutBothGenTops;   //!
TBranch        *b_nofTTEventsWithoutGenTop;   //!
TBranch        *b_nofTTEventsWithoutGenAntiTop;   //!
TBranch        *b_nofTTEventsWithoutBothGenTopsWithStatus22;   //!
TBranch        *b_nofTTEventsWithoutGenTopWithStatus22;   //!
TBranch        *b_nofTTEventsWithoutGenAntiTopWithStatus22;   //!
TBranch        *b_nofTTEventsWithoutBothGenTopsWithStatus62;   //!
TBranch        *b_nofTTEventsWithoutGenTopWithStatus62;   //!
TBranch        *b_nofTTEventsWithoutGenAntiTopWithStatus62;   //!


double lumiWeight, numWeight, scaleFactor, widthSF;
double topPtRewSF, topPtSF, antiTopPtSF;
bool foundTop62, foundAntiTop62;
vector<unsigned int> bJetId;
double bdiscrTop, bdiscrTop2, tempbdiscr;
int labelB1, labelB2;
int labelsReco[4];
double massForWidth;

string catSuffix = "";
string catSuffixList[] = {"_CM", "_WM", "_NM"};
bool isCM, isWM, isNM;


/// Define TLVs
TLorentzVector muon, jet, mcpart;
TLorentzVector WCandidate;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
vector<TLorentzVector> selectedJetsAKF;
vector<TLorentzVector> selectedJetsKFcorrected;
vector<TLorentzVector> selectedJetsKFMatched;
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
int genneutrino = -9999;
bool muonmatched = false;
bool muPlusFromTop = false, muMinusFromTop = false;
bool foundAllTTbarComponents = false;
int bjjDecay[] = {-9999, -9999, -9999}, blvDecay[] = {-9999, -9999, -9999};
vector<unsigned int> partonId;


/// KinFitter
bool doneKinFit = false;
TKinFitter* kFitter;
TKinFitter* kFitterMatched;
bool addWMassKF = true;
bool addEqMassKF = false;
int kFitVerbosity = 0;
double kFitChi2 = 99., kFitChi2Matched = 99.;
int nofAcceptedKFit = 0, nofAcceptedKFitMatched = 0;
bool passKFChi2MatchedCut = false;

ofstream txtDebugTopMass;


/// Variables
double M3, Ht, min_Mlb, dRLepB;
double M3_aKF, Ht_aKF;
double reco_W_mass_bKF, reco_top_mass_bKF, reco_top_pt_bKF, reco_mlb_bKF, reco_dRLepB_lep_bKF, reco_dRLepB_had_bKF, reco_ttbar_mass_bKF, redTopMass_bKF;
double reco_W_mass_aKF, reco_top_mass_aKF, reco_top_pt_aKF, reco_mlb_aKF, reco_dRLepB_lep_aKF, reco_dRLepB_had_aKF, reco_ttbar_mass_aKF, redTopMass;
double aveTopMassLL = aveTopMass[2];

double matched_W_mass_q, matched_top_mass_q;
double matched_W_mass_j, matched_top_mass_j, matched_top_mass_j_akF;
double matched_mlb_corr, matched_ttbarMass_corr, matched_dR_lep_b_corr;
double matched_mlb_wrong, matched_ttbarMass_wrong, matched_dR_lep_b_wrong;


/// Input samples & reweighting
int thisGenWidthId;
double genWidthArray[] = {0.2, 0.5, 1., 4., 8.};
string genWidthString[] = {"g0p2", "g0p5", "g1", "g4", "g8"};
const int nGenWidths = sizeof(genWidthArray)/sizeof(genWidthArray[0]);

int nEventsTT[] = {19285604,19579453,153843293,18940139,19524579};

double reweightArray[] = {0.2, 0.5, 0.75, 1., 1.5, 2., 3., 4., 8.};
string reweightString[] = {"s0p2","s0p5", "s0p75", "s1", "s1p5", "s2", "s3", "s4", "s8"};
const int nReweightings = sizeof(reweightArray)/sizeof(reweightArray[0]);

double m_top = -1.;
double m_antitop = -1.;
double m_bjj = -1.;
double m_blv = -1.;
double m_hadr = -1;
double m_lept = -1;
double evWeight_hadr = 1.;
double evWeight_lept = 1.;
double evWeight_prod = 1.;
double evWeight_prodsqrt = 1.;



/// Meta
string strSyst = "";
double eqLumi;

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
  cout << "* The following corrections are applied:    *" << endl;
  cout << "*   - Jet Energy Corrections                *" << endl;
  if ( systStr.find("JECdown") != std::string::npos )
    cout << "*   - Jet Energy Scale: scale down          *" << endl;
  else if ( systStr.find("JECup") != std::string::npos )
    cout << "*   - Jet Energy Scale: scale up            *" << endl;
  cout << "*   - Jet Energy Resolution: ";
  if ( systStr.find("JERdown") != std::string::npos )    cout << "scale down     *" << endl;
  else if ( systStr.find("JERup") != std::string::npos ) cout << "scale up       *" << endl;
  else                                                 cout << "nominal        *" << endl;
  cout << "*   - Jet/lepton Cleaning                   *" << endl;
  if (doMETCleaning) cout << "*   - MET Cleaning                          *" << endl;
  if (applyTopPtReweighting) cout << "*   - Top pT reweighting                    *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  string channel;
  if ( argc == 1 ) channel = "mu";
  //else if ( CharSearch(argv[1], "mu") || CharSearch(argv[1], "Mu") || CharSearch(argv[1], "MU") || CharSearch(argv[1], "m") ) channel = "mu";
  //else if ( CharSearch(argv[1], "el") || CharSearch(argv[1], "El") || CharSearch(argv[1], "EL") || CharSearch(argv[1], "e") ) channel = "el";
  //else if ( (argv[1]).find("all") != std::string::npos || (argv[1]).find("All") != std::string::npos || (argv[1]).find("ALL") != std::string::npos ) channel = "all";
  
  if (calculateAverageMass)
  {
    makePlots = false;
    doGenOnly = false;
  }
  if (calculateResolutionFunctions)
  {
    testTTbarOnly = true;
    doGenOnly = true;
    makePlots = false;
    calculateAverageMass = false;
    doKinFit = false;
  }
  if (test) makePlots = false;
  if (testHistos)
  {
    makePlots = true;
    doGenOnly = false;
  }
  if (! makePlots)
  {
    makeReweightedPlots = false; 
    makeMatchingPlots = false;
    makeRecoPlots = false;
  }
  
  string pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (makePlots)
  {
    pathOutput += "reweighting/";
    mkdir(pathOutput.c_str(),0777);
    // Add channel to output path
    pathOutput += channel+"/";
    mkdir(pathOutput.c_str(),0777);
    if (testHistos)
    {
      pathOutput += "test/";
      mkdir(pathOutput.c_str(),0777);
    }
    // Give timestamp to output path
    pathOutput += dateString+"/";
    mkdir(pathOutput.c_str(),0777);
  }
  
  pathNtuples = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.first+"/";
  cout << "Using Ntuples from " << ntupleDate.first << ". This corresponds to systematics: " << systStr << endl;
  if (calculateAverageMass) cout << "Calculating average mass values..." << endl;
  if (testHistos) cout << "Testing histogram consistency..." << endl;
  if (doGenOnly) cout << "Running only matching..." << endl;
  if (applyWidthSF) cout << "TTbar sample width will be scaled by a factor " << scaleWidth << endl;
  else scaleWidth = 1.;
  
  /// xml file
  string xmlFileName ="config/topWidth_ttbarOnly.xml";
  
  const char *xmlfile = xmlFileName.c_str();
  
  cout << " - Using config file " << xmlfile << endl;
  
  
  
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  datasets.clear();
  TTreeLoader treeLoader; 
  treeLoader.LoadDatasets(datasets, xmlfile);
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  EventReweighting *rew = new EventReweighting(false);  // no correction for number of events
  ResolutionFunctions* rf = new ResolutionFunctions(calculateResolutionFunctions, true);
  KinFitter *kf;
  KinFitter *kfMatched;
  Likelihood *like;
  
  if (! calculateResolutionFunctions)
  {
    kf = new KinFitter("PlotsForResolutionFunctions_testFit_170608_S.root", addWMassKF, addEqMassKF);
    kfMatched = new KinFitter("PlotsForResolutionFunctions_testFit_170608_S.root", addWMassKF, addEqMassKF);
  }
  
  if (makePlots)
  {
    InitHisto1D();
    InitHisto2D();
  }  
  
  if (calculateAverageMass)
  {
    mkdir("averageMass/",0777);
    txtMassGenPMatched.open(("averageMass/mass_genp_matched_TT_"+dateString+".txt").c_str());
    txtMassGenJMatched.open(("averageMass/mass_genj_matched_TT_"+dateString+".txt").c_str());
    txtMassRecoCM.open(("averageMass/mass_reco_matched_TT_"+dateString+".txt").c_str());
    txtMassRecoWMNM.open(("averageMass/mass_reco_notCorrectMatch_TT_"+dateString+".txt").c_str());
    txtMassRecoNM.open(("averageMass/mass_reco_notMatched_TT_"+dateString+".txt").c_str());
    txtMassRecoWM.open(("averageMass/mass_reco_wrongPerm_TT_"+dateString+".txt").c_str());
  }
  
  txtDebugTopMass.open("debug_missing_topQ.txt");
  txtDebugTopMass << "## Events where top quark(s) not found in genParticle collection" << endl;
  txtDebugTopMass << "#  If lepton charge > 0 : leptonically decaying top, hadronically decaying antitop" << endl;
  txtDebugTopMass << "#  If lepton charge < 0 : hadronically decaying top, leptonically decaying antitop" << endl;
  
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName, slumi;
  double timePerDataSet[datasets.size()];
  
  int nEntries;
  // double fracDataEras[2] = {-1.};
//   //GetEraFraction(fracDataEras);
//   fracDataEras[0] = lumi_runBCDEF/(lumi_runBCDEF + lumi_runGH);
//   fracDataEras[1] = lumi_runGH/(lumi_runBCDEF + lumi_runGH);
//   if ( fracDataEras[0] == -1. || fracDataEras[1] == -1. )
//   {
//     cerr << "Something went wrong with the fraction calculation for muon SFs!" << endl;
//     exit(1);
//   }
//   cout << "The muon scale factors will be scaled by " << fracDataEras[0] << " for eras B-F and " << fracDataEras[1] << " for eras G-H." << endl;
  
  
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
    
    isTTbar = false;
    if ( dataSetName.find("TT") != std::string::npos )
    {
      isTTbar = true;
      doReweighting = true;
      thisGenWidthId = 2;
      if ( dataSetName.find("width") != std::string::npos || dataSetName.find("Width") != std::string::npos )
      {
        if ( dataSetName.find("x0p2") != std::string::npos )
        {
          thisGenWidthId = 0;
          scaleWidth = 0.2;
        }
        else if ( dataSetName.find("x0p5") != std::string::npos )
        {
          thisGenWidthId = 1;
          scaleWidth = 0.5;
        }
        else if ( dataSetName.find("x4") != std::string::npos )
        {
          thisGenWidthId = 3;
          scaleWidth = 4.;
        }
        else if ( dataSetName.find("x8") != std::string::npos )
        {
          thisGenWidthId = 4;
          scaleWidth = 8.;
        }
        applyWidthSF = false;
        doReweighting = false;
      }
    }
    
    
    if (calculateAverageMass)
    {
      txtMassReco.open(("averageMass/mass_reco_"+dataSetName+"_"+dateString+".txt").c_str());
    }
    
    
    string ntupleFileName = "Ntuples_"+dataSetName+".root";
    tFileMap[dataSetName.c_str()] = new TFile((pathNtuples+ntupleFileName).c_str(),"READ"); //create TFile for each dataset
    
    string tTreeName = "tree";
    string tStatsTreeName = "stats";
    
    /// Get meta data
    tStatsTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tStatsTreeName.c_str());
    GetMetaData(tStatsTree[dataSetName.c_str()], false);
    
    tStatsTree[(dataSetName).c_str()]->GetEntry(0);
    
    /// eqLumi calculation
    if (isData) lumiWeight = 1.;
    else
    {
      eqLumi = (double)GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", false)/datasets[d]->Xsection();  // 1/pb
      lumiWeight = Luminosity/eqLumi;
    }
    
    // Scale number of events
    //numWeight = TMath::MinElement(nEventsTT)/nEventsTT[d];  // Don't scale up events
    numWeight = nEventsTT[2]/nEventsTT[thisGenWidthId];  // Don't put extra SFs on reweighted samples
    
    /// Get data
    tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
    nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
    cout << "                nEntries  : " << nEntries << endl;
    cout << "                eqLumi    : " << eqLumi << "/pb = " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", false) << " / " << datasets[d]->Xsection() << " pb" << endl;
    cout << "                lumiWeight: " << lumiWeight << endl;
    cout << "                numWeight : " << numWeight << /*" = " << TMath::MinElement(nEventsTT) << " / " << nEventsTT[d] <<*/ endl;
    
    // Set branch addresses and branch pointers
    InitTree(tTree[dataSetName.c_str()], false);
    
    
    
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
      
      //if (ievt > 588000) break;
      //if (ievt != 2062 && ievt != 75831 && ievt != 113603 && ievt != 115687 && ievt != 155732 && ievt != 163161 && ievt != 186900 && ievt != 215759 && ievt != 233634 && ievt != 238021 && ievt != 243052 && ievt != 243674 && ievt != 266399 && ievt != 317190 && ievt != 317752 && ievt != 325854 && ievt != 330813 && ievt != 333620 && ievt != 347247 && ievt != 439571 && ievt != 450329 && ievt != 491328 && ievt != 510024 && ievt != 514196 && ievt != 538345 && ievt != 570225 && ievt != 576194 && ievt != 577278 && ievt != 587570) continue;
      
      
      /// Load event
      tTree[(dataSetName).c_str()]->GetEntry(ievt);
      
      
      
      //////////////////////
      ///  Fill objects  ///
      //////////////////////
      
      muon.SetPtEtaPhiE(muon_pt[0], muon_eta[0], muon_phi[0], muon_E[0]);
      selectedLepton.push_back(muon);
      
      for (int iJet = 0; iJet < nJets; iJet++)
      {
        jet.Clear();
        jet.SetPtEtaPhiE(jet_pt[iJet], jet_eta[iJet], jet_phi[iJet], jet_E[iJet]);
        selectedJets.push_back(jet);
      }
      
      if ( ! calculateResolutionFunctions && selectedJets.size() > 4 ) continue;
      nofHardSelected++;
      
      if (! passedMETFilter) continue;
      nofMETCleaned++;
      
      for (int iJet = 0; iJet < selectedJets.size(); iJet++)
      {
        if ( jet_bdiscr[iJet] > CSVv2Medium )
        {
          selectedBJets.push_back(selectedJets[iJet]);
          bJetId.push_back(iJet);  /// selectedBJets[j] = selectedJets[bJetId[j]]
        }
      }
      //std::sort(selectedBJets.begin(),selectedBJets.end(),HighestPt());  // already the case
      
      /// label jets with highest b discr
      for (int iJet = 0; iJet < selectedBJets.size(); iJet++)
      {
        tempbdiscr = jet_bdiscr[bJetId[iJet]];
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
      
      
      if (! applyWidthSF ) widthSF = 1.;
      else if ( applyWidthSF && ! isTTbar ) widthSF = 1.;  // also for data
      
      
      
      /////////////////////////////
      ///  JET PARTON MATCHING  ///
      /////////////////////////////
      
      for (int iMC = 0; iMC < nMCParticles; iMC++)
      {
        mcpart.Clear();
        mcpart.SetPtEtaPhiE(mc_pt[iMC], mc_eta[iMC], mc_phi[iMC], mc_E[iMC]);
        mcParticles.push_back(mcpart);
      }
      
      foundTop62 = false; foundAntiTop62 = false;
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
        
        if ( mc_pdgId[i] == pdgID_top )  // top
        {
          if ( mc_status[i] == 62 )
          {
            foundTop62 = true;
            if ( mcParticles[i].Pt() < 400 ) topPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
            else topPtSF = TMath::Exp(0.0615-0.0005*400.);
          }
          else if (! foundTop62 )
          {
            if ( mcParticles[i].Pt() < 400 ) topPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
            else topPtSF = TMath::Exp(0.0615-0.0005*400.);
          }
        }
        else if ( mc_pdgId[i] == -pdgID_top )  // antitop
        {
          if ( mc_status[i] == 62 )
          {
            foundAntiTop62 = true;
            if ( mcParticles[i].Pt() < 400 ) antiTopPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
            else antiTopPtSF = TMath::Exp(0.0615-0.0005*400.);
          }
          else if (! foundAntiTop62)
          {
            if ( mcParticles[i].Pt() < 400 ) antiTopPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
            else antiTopPtSF = TMath::Exp(0.0615-0.0005*400.);
          }
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
        
        /// Neutrinos
        if ( mc_pdgId[i] == -14 && mc_mother[i] == -24 && mc_granny[i] == -pdgID_top )    // vbar, W-, tbar
        {
          if ( mc_status[i] == 23 ) genneutrino = i;
          else if ( mc_status[i] != 23 && genneutrino == -9999 ) genneutrino = i;
        }

        if ( mc_pdgId[i] == 14 && mc_mother[i] == 24 && mc_granny[i] == pdgID_top )    // v, W+, t
        {
          if ( mc_status[i] == 23 ) genneutrino = i;
          else if ( mc_status[i] != 23 && genneutrino == -9999 ) genneutrino = i;
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
      
      if ( isTTbar && (topQuark == -9999 || antiTopQuark == -9999) )
      {
        txtDebugTopMass << "Event " << ievt << ";  Event nb. " << evt_num << "; Run nb. " << run_num << endl;
        txtDebugTopMass << "Top mass id: " << topQuark << "; antiTop mass id: " << antiTopQuark << "; Lepton charge: " << muon_charge[0] << endl;
        continue;
      }
      
      if (isTTbar) topPtRewSF = TMath::Sqrt(topPtSF*antiTopPtSF);
      else topPtRewSF = 1.;
      if (applyTopPtReweighting) scaleFactor *= topPtRewSF;
      
      
      
      /////////////////////////////////////////
      ///  Scale factor ttbar sample width  ///
      /////////////////////////////////////////
      
      if ( muon_charge[0] > 0 ) massForWidth = (mcParticles[antiTopQuark]).M();
      else if ( muon_charge[0] < 0 ) massForWidth = (mcParticles[topQuark]).M();
      
      if ( applyWidthSF && isTTbar )
      {
        widthSF = rew->EventWeightCalculator(massForWidth, scaleWidth);
        
        if ( widthSF != widthSF )  // widthSF = NaN
        {
          cout << "Event " << ievt << ";  Event nb. " << evt_num << "; Run nb. " << run_num << endl;
          cout << "Top mass: " << (mcParticles[topQuark]).M() << "; antiTop mass: " << (mcParticles[antiTopQuark]).M() << "; Lepton charge: " << muon_charge[0] << "; width SF: " << widthSF << endl;
          
          continue;
        }
        
        if (makePlots) histo1D["width_SF"]->Fill(widthSF);
      }  // end applyWidthSF
      
      
      
      //////////////////
      ///  Matching  ///
      //////////////////
      
      if (doMatching)
      {
        TruthMatching(partons, selectedJets, MCPermutation);
        
        if (test && verbose > 3 && all4PartonsMatched)
        {
          for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
          {
            cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Pt: " << setw(7) << left << mc_pt[partonId[MCPermutation[iMatch].second]] << "  Eta: " << mc_eta[partonId[MCPermutation[iMatch].second]] << "  Phi: " << mc_phi[partonId[MCPermutation[iMatch].second]] << endl;
            cout << "Event  " << right << setw(4) << ievt << ";  Matched jet    " << iMatch << "  Pt: " << setw(7) << left << jet_pt[MCPermutation[iMatch].first] << "  Eta: " << jet_eta[MCPermutation[iMatch].first] << "  Phi: " << jet_phi[MCPermutation[iMatch].first] << endl;
          }
        }
        
        ///////////////////
        ///  Reconstruct gen event
        ///////////////////
        
        for (unsigned int i = 0; i < partons.size(); i++)  /// partons[j] = mcParticles[partonId[j]]
        {
          if ( muMinusFromTop )
          {
            // hadronic decay  --> top
            if ( mc_pdgId[partonId[i]] == 5 && mc_mother[partonId[i]] == pdgID_top ) bjjDecay[0] = partonId[i];
            else if ( abs(mc_pdgId[partonId[i]]) < 5 && mc_mother[partonId[i]] == 24 && mc_granny[partonId[i]] == pdgID_top )
            {
              //cout << "In hadronic loop: " << setw(3) << right << partonId[i] << "  pdgId: " << setw(3) << mc_pdgId[partonId[i]] << "  Mother: " << setw(4) << mc_mother[partonId[i]] << "  Granny: " << setw(4) << mc_granny[partonId[i]] << endl;

              if ( bjjDecay[1] == -9999 ) bjjDecay[1] = partonId[i];
              else if ( bjjDecay[2] == -9999 ) bjjDecay[2] = partonId[i];
              else { cerr << "ERROR: Too many partons for hadronic top decay!  Event: " << ievt << endl; continue; }
            }
            
            if ( mc_pdgId[partonId[i]] == -5 && mc_mother[partonId[i]] == -pdgID_top ) blvDecay[0] = partonId[i];  // leptonic b
          }
          else if ( muPlusFromTop )
          {
            // hadronic decay  --> antitop
            if ( mc_pdgId[partonId[i]] == -5 && mc_mother[partonId[i]] == -pdgID_top ) bjjDecay[0] = partonId[i];
            else if ( abs(mc_pdgId[partonId[i]]) < 5 && mc_mother[partonId[i]] == -24 && mc_granny[partonId[i]] == -pdgID_top )
            {
              //cout << "In hadronic loop: " << setw(3) << right << partonId[i] << "  pdgId: " << setw(3) << mc_pdgId[partonId[i]] << "  Mother: " << setw(4) << mc_mother[partonId[i]] << "  Granny: " << setw(4) << mc_granny[partonId[i]] << endl;
              
              if ( bjjDecay[1] == -9999 ) bjjDecay[1] = partonId[i];
              else if ( bjjDecay[2] == -9999 ) bjjDecay[2] = partonId[i];
              else { cerr << "ERROR: Too many partons for hadronic antitop decay!  Event: " << ievt << endl; continue; }
            }
            
            if ( mc_pdgId[partonId[i]] == 5 && mc_mother[partonId[i]] == pdgID_top ) blvDecay[0] = partonId[i];  // leptonic b
          }
        }  // end loop partons
        blvDecay[1] = genmuon;
        blvDecay[2] = genneutrino;
        
        if ( bjjDecay[0] != -9999 && bjjDecay[1] != -9999 && bjjDecay[2] != -9999 && blvDecay[0] != -9999 && blvDecay[1] != -9999 && blvDecay[2] != -9999)
        {
          foundAllTTbarComponents = true;
          if ( test && verbose > 2 ) cout << "All ttbar components found... Skipping event " << ievt << "..." << endl;
        }
        
        
        ///////////////////
        ///  
        ///////////////////
        
        if (foundAllTTbarComponents)
        {
          m_top = (mcParticles[topQuark]).M();
          m_antitop = (mcParticles[antiTopQuark]).M();
          m_bjj = (mcParticles[bjjDecay[0]] + mcParticles[bjjDecay[1]] + mcParticles[bjjDecay[2]]).M();
          m_blv = (mcParticles[blvDecay[0]] + mcParticles[blvDecay[1]] + mcParticles[blvDecay[2]]).M();
          
          m_hadr = -1.;
          m_lept = -1.;
          if ( muMinusFromTop )
          {
            m_hadr = m_top;
            m_lept = m_antitop;
          }
          else if ( muPlusFromTop )
          {
            m_hadr = m_antitop;
            m_lept = m_top;
          }
          
          if (applyWidthSF && doReweighting && makeReweightedPlots)
          {
            for (int s = 0; s < nReweightings; s++)
            {
              evWeight_hadr = rew->EventWeightCalculator(m_hadr, reweightArray[s]);
              evWeight_lept = rew->EventWeightCalculator(m_lept, reweightArray[s]);
              evWeight_prod = evWeight_hadr*evWeight_lept;
              evWeight_prodsqrt = TMath::Sqrt(evWeight_prod);
              
              histo1D[("top_mass_hadr_gen_"+reweightString[s]).c_str()]->Fill(m_hadr, evWeight_hadr);
              histo1D[("top_mass_lept_gen_"+reweightString[s]).c_str()]->Fill(m_lept, evWeight_hadr);
              histo1D[("bjj_mass_gen_"+reweightString[s]).c_str()]->Fill(m_bjj, evWeight_hadr);
              histo1D[("blv_mass_gen_"+reweightString[s]).c_str()]->Fill(m_blv, evWeight_hadr);
              
              histo1D[("top_mass_hadr_gen_prodSF_"+reweightString[s]).c_str()]->Fill(m_hadr, evWeight_prod);
              histo1D[("top_mass_lept_gen_prodSF_"+reweightString[s]).c_str()]->Fill(m_lept, evWeight_prod);
              histo1D[("bjj_mass_gen_prodSF_"+reweightString[s]).c_str()]->Fill(m_bjj, evWeight_prod);
              histo1D[("blv_mass_gen_prodSF_"+reweightString[s]).c_str()]->Fill(m_blv, evWeight_prod);
              
              histo1D[("top_mass_hadr_gen_prodsqrtSF_"+reweightString[s]).c_str()]->Fill(m_hadr, evWeight_prodsqrt);
              histo1D[("top_mass_lept_gen_prodsqrtSF_"+reweightString[s]).c_str()]->Fill(m_lept, evWeight_prodsqrt);
              histo1D[("bjj_mass_gen_prodsqrtSF_"+reweightString[s]).c_str()]->Fill(m_bjj, evWeight_prodsqrt);
              histo1D[("blv_mass_gen_prodsqrtSF_"+reweightString[s]).c_str()]->Fill(m_blv, evWeight_prodsqrt);
              
              histo1D[("Width_SF_hadr_"+reweightString[s]).c_str()]->Fill(evWeight_hadr);
              histo1D[("Width_SF_lept_"+reweightString[s]).c_str()]->Fill(evWeight_lept);
              histo1D[("Width_SF_prod_"+reweightString[s]).c_str()]->Fill(evWeight_prod);
              histo1D[("Width_SF_prodsqrt_"+reweightString[s]).c_str()]->Fill(evWeight_prodsqrt);
              
              histo2D[("Width_SF_hadr_lept_"+reweightString[s]).c_str()]->Fill(evWeight_hadr, evWeight_lept);
              histo2D[("Width_SF_hadr_prod_"+reweightString[s]).c_str()]->Fill(evWeight_hadr, evWeight_prod);
              histo2D[("Width_SF_hadr_prodsqrt_"+reweightString[s]).c_str()]->Fill(evWeight_hadr, evWeight_prodsqrt);
              histo2D[("top_mass_hadr_gen_vs_weight_"+reweightString[s]).c_str()]->Fill(m_hadr, evWeight_hadr);
              //histo2D[("top_mass_lept_gen_vs_weight_"+reweightString[s]).c_str()]->Fill(m_lept, evWeight_lept);
            }
          }
          else if (makePlots)
          {
            histo1D[("top_mass_hadr_gen_"+genWidthString[thisGenWidthId]).c_str()]->Fill(m_hadr, numWeight);
            histo1D[("top_mass_hadr_gen_"+genWidthString[thisGenWidthId]+"_fewerBins").c_str()]->Fill(m_hadr, numWeight);
            histo1D[("top_mass_lept_gen_"+genWidthString[thisGenWidthId]).c_str()]->Fill(m_lept, numWeight);
            histo1D[("bjj_mass_gen_"+genWidthString[thisGenWidthId]).c_str()]->Fill(m_bjj, numWeight);
            histo1D[("blv_mass_gen_"+genWidthString[thisGenWidthId]).c_str()]->Fill(m_blv, numWeight);
            
            histo2D["top_mass_hadr_lept_gen_"+genWidthString[thisGenWidthId]]->Fill(m_hadr, m_lept, numWeight);
            histo2D["top_mass_hadr_bjj_gen_"+genWidthString[thisGenWidthId]]->Fill(m_hadr, m_bjj, numWeight);
            histo2D["top_mass_lept_blv_gen_"+genWidthString[thisGenWidthId]]->Fill(m_lept, m_blv, numWeight);
            histo2D["top_mass_bjj_blv_gen_"+genWidthString[thisGenWidthId]]->Fill(m_bjj, m_blv, numWeight);
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
          
          if (isTTbar && all4PartonsMatched && calculateResolutionFunctions)
          {
            rf->fillJets(partonsMatched, jetsMatched);
            
            if (muonmatched) rf->fillMuon(mcParticles[genmuon], selectedLepton[0]);
            //if (electronmatched) rf->fillElectron(...)
            
          }  // end rf
          
          
          matched_W_mass_j = (jetsMatched[0] + jetsMatched[1]).M();
          matched_W_mass_q = (partonsMatched[0] + partonsMatched[1]).M();
          matched_top_mass_j = (jetsMatched[0] + jetsMatched[1] + jetsMatched[2]).M();
          matched_top_mass_q = (partonsMatched[0] + partonsMatched[1] + partonsMatched[2]).M();
          
          if (calculateAverageMass)
          {
            txtMassGenPMatched << ievt << "  " << matched_top_mass_q << endl;
            txtMassGenJMatched << ievt << "  " << matched_top_mass_j << endl;
          }
          
          
          /// KF for matched jets
          if (doKinFit)
          {
            kFitterMatched = kfMatched->doFit(jetsMatched[0], jetsMatched[1], kFitVerbosity);
            
            if ( kFitterMatched->getStatus() != 0 )  // did not converge
            {
              if (test && verbose > 2) cout << "Event " << ievt << ": Fit for matched events did not converge..." << endl;
              continue;
            }             
            kFitChi2Matched = kFitterMatched->getS();
            if (test && verbose > 4) cout << "Fit converged: Chi2 = " << kFitChi2Matched << endl;
            
            if ( applyKinFitCut && kFitChi2Matched < kinFitCutValue ) passKFChi2MatchedCut = true;
            if (passKFChi2MatchedCut)
            {
              nofAcceptedKFitMatched++;
              
              selectedJetsKFMatched.clear();
              selectedJetsKFMatched = kfMatched->getCorrectedJets();
              
              if ( selectedJetsKFMatched.size() == 2 ) selectedJetsKFMatched.push_back(jetsMatched[2]);
              
              matched_top_mass_j_akF = (selectedJetsKFMatched[0] + selectedJetsKFMatched[1] + selectedJetsKFMatched[2]).M();
              
            }  // passKFChi2MatchedCut
          }  // end KF
          
          
          if (isTTbar && makeMatchingPlots)
          {
            FillMatchingPlots();
          }
          
        }  // end hadronicTopJetsMatched
        
      }  // end doMatching
      
      
      
      if (doGenOnly) continue;
      
      
      
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
      
      
      
      ///////////////////////////////////
      ///  CHECK MATCHED COMBINATION  ///
      ///////////////////////////////////
      
      ///
      // 3 possibilities:
      // - correct top match: 3 jets selected with reco method correspond to the 3 matched jets (n.b. this is also true when the jets originating from the W boson and the b jet do not exactly correspond to the matched jets, because we are only interested in the reconstructed top quark.)
      // - wrong permutation: the correct jet combination exists in the selected jets, but is not chosen by the reco method.
      // - wrong (no) match:  the correct jet combination does not exist in the selected jets (e.g. when one jet is not selected.)
      
      
      if (! isData)
      {
        if (hadronicTopJetsMatched)
        {
          /// Correct match
          if ( ( labelsReco[0] == MCPermutation[0].first || labelsReco[0] == MCPermutation[1].first || labelsReco[0] == MCPermutation[2].first ) && ( labelsReco[1] == MCPermutation[0].first || labelsReco[1] == MCPermutation[1].first || labelsReco[1] == MCPermutation[2].first ) && ( labelsReco[2] == MCPermutation[0].first || labelsReco[2] == MCPermutation[1].first || labelsReco[2] == MCPermutation[2].first ) )  // correct jets for top quark
          {
            isCM = true;
            nofCorrectlyMatched++;
          }
          else  // wrong permutation
          {
            isWM = true;
            nofNotCorrectlyMatched++;
          }
        }  // end hadrTopMatch
        else  // no match
        {
          isNM = true;
        }
        
        
        if ( (! isCM && ! isWM && ! isNM) || (isCM && isWM) || (isCM && isNM) || (isWM && isNM) )
          cerr << "Something wrong with trigger logic CM/WM/NM !! " << endl;
        
      }  // not Data
      
      
      if (isCM) catSuffix = catSuffixList[0];
      else if (isWM) catSuffix = catSuffixList[1];
      else if (isNM) catSuffix = catSuffixList[2];
      
      
      /// Fill variables before performing kinFit
      reco_W_mass_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).M();
      reco_top_mass_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).M();
      reco_top_pt_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).Pt();
      reco_mlb_bKF = (selectedLepton[0] + selectedJets[labelsReco[3]]).M();
      reco_dRLepB_lep_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[3]], selectedLepton[0] );  // deltaR between lepton and leptonic b jet
      reco_dRLepB_had_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
      reco_ttbar_mass_bKF = reco_mlb_bKF + reco_top_mass_bKF;
      redTopMass_bKF = reco_top_mass_bKF/aveTopMassLL;
      
      if (makeRecoPlots)
      {
        if (isTTbar && doKinFit) FillKinFitPlots(doneKinFit);
      }
      
      
      
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
          if (test && verbose > 2) cout << "Event " << ievt << ": Fit did not converge..." << endl;
          continue;
        }
        
        kFitChi2 = kFitter->getS();
        if (test && verbose > 4) cout << "Fit converged: Chi2 = " << kFitChi2 << endl;
        
        doneKinFit = true;
        if (makeRecoPlots)
        {
          if (! isData) histo1D["KF_Chi2"+catSuffix+"_wide"]->Fill(kFitChi2);
          if (isTTbar) histo1D["KF_Chi2_TT"]->Fill(kFitChi2);
        }
        
        if (isCM) nofCorrectlyMatchedAKFNoCut++;
        else if (isWM) nofNotCorrectlyMatchedAKFNoCut++;
        else if (isNM) nofNoMatchAKFNoCut++;
        
        if ( applyKinFitCut && kFitChi2 > kinFitCutValue ) continue;
        nofAcceptedKFit++;
        if (hadronicTopJetsMatched) nofHadrMatchedEventsAKF++;
        if (isCM) nofCorrectlyMatchedAKF++;
        else if (isWM) nofNotCorrectlyMatchedAKF++;
        
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
      else if ( selectedJetsKFcorrected.size() == 2 ) selectedJetsKFcorrected.push_back(selectedJets[labelsReco[2]]);
      
      /// Make pT ordered jet collection after KF
      selectedJetsAKF = selectedJetsKFcorrected;
      selectedJetsAKF.push_back(selectedJets[labelsReco[3]]);
      std::sort(selectedJetsAKF.begin(),selectedJetsAKF.end(),HighestPt());
      
      /// Define variables
      reco_W_mass_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).M();
      reco_top_mass_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).M();
      reco_top_pt_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).Pt();
      reco_mlb_aKF = (selectedLepton[0] + selectedJets[labelsReco[3]]).M();
      reco_dRLepB_lep_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[3]], selectedLepton[0] );  // deltaR between lepton and leptonic b jet
      reco_dRLepB_had_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
      reco_ttbar_mass_aKF = reco_mlb_aKF + reco_top_mass_aKF;
      redTopMass = reco_top_mass_aKF/aveTopMassLL;
      
      if ( reco_top_mass_aKF < 0. )
        PrintKFDebug(ievt);
      
      if (calculateAverageMass) txtMassReco << ievt << "  " << reco_top_mass_aKF << endl;
      
      if ( doKinFit && makeRecoPlots )
      {
        if (isTTbar) FillKinFitPlots(doneKinFit);
        if (! isData)
        {
          histo1D["allSim_top_mass"]->Fill(reco_top_mass_aKF, lumiWeight*scaleFactor*widthSF);
          histo1D["allSim_red_top_mass"]->Fill(redTopMass, lumiWeight*scaleFactor*widthSF);
        }
      }
      
      
      
      if (calculateAverageMass && ! isData)
      {
        if (isCM) txtMassRecoCM << ievt << "  " << reco_top_mass_aKF << endl;
        else
        {
          txtMassRecoWMNM << ievt << "  " << reco_top_mass_aKF << endl;
          if (isWM)
            txtMassRecoWM << ievt << "  " << reco_top_mass_aKF << endl;
          else if (isNM)
            txtMassRecoNM << ievt << "  " << reco_top_mass_aKF << endl;
        }
      }  // end aveMassCalc
      
            
      //Fill histos
      if (makeRecoPlots)
      {
        histo1D[("red_top_mass_"+dataSetName).c_str()]->Fill(redTopMass);
        
        FillCatsPlots(catSuffix);
        
      }  // end makePlots
      
      
      
      
    }  // end loop events
    
    
    cout << endl;  /// Stronger selection in this analyser compared to Ntuples ==> endEvent --> nofHardSelected
    cout << "Number of events with exactly 4 jets with pT > 30 GeV: " << nofHardSelected << " (" << 100*((float)nofHardSelected/(float)endEvent) << "%)" << endl;
    cout << "Number of events with clean MET: " << nofMETCleaned << " (" << 100*((float)nofMETCleaned/(float)nofHardSelected) << "%)" << endl;
    if (doKinFit) cout << "Number of clean events accepted by kinFitter: " << nofAcceptedKFit << " (" << 100*((float)nofAcceptedKFit/(float)nofMETCleaned) << "%)" << endl;
    
    //if ( isTTbar || dataSetName.find("ST") != std::string::npos )
    if (! isData && nofHadrMatchedEvents > 0 )
    {
      cout << "Number of matched events: " << setw(8) << right << nofMatchedEvents << endl;
      cout << "Number of events with hadronic top matched (before KF): " << setw(8) << right << nofHadrMatchedEvents << " (" << 100*((float)nofHadrMatchedEvents/(float)nofMETCleaned) << "%)" << endl;
      if (doKinFit) cout << "Number of events with hadronic top matched (after KF):  " << setw(8) << right << nofHadrMatchedEventsAKF << " (" << 100*((float)nofHadrMatchedEventsAKF/(float)nofAcceptedKFit) << "%)" << endl;
      if (! doGenOnly)
      {
        cout << "Correctly matched reconstructed events:     " << setw(8) << right << nofCorrectlyMatched << endl;
        cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatched << endl;
        if ( nofCorrectlyMatched != 0 || nofNotCorrectlyMatched != 0 )
          cout << "   ===> This means that " << 100*(float)nofCorrectlyMatched / (float)(nofCorrectlyMatched + nofNotCorrectlyMatched) << "% of matched events is correctly matched." << endl;
        
        if (doKinFit)
        {
          cout << "                        " << 100*(float)nofCorrectlyMatched / (float)nofMETCleaned << "% of all events is correctly matched before kinfitter." << endl;
          cout << " --- Kinematic fit" << endl;
          cout << " --- Before chi2 cut --- " << endl;
          cout << "Correctly matched reconstructed events    : " << setw(8) << right << nofCorrectlyMatchedAKFNoCut << endl;
          cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatchedAKFNoCut << endl;
          if ( nofCorrectlyMatchedAKFNoCut != 0 || nofNotCorrectlyMatchedAKFNoCut != 0 )
            cout << "   ===> This means that " << 100*(float)nofCorrectlyMatchedAKFNoCut / (float)(nofCorrectlyMatchedAKFNoCut + nofNotCorrectlyMatchedAKFNoCut) << "% of matched events is correctly matched after KF." << endl;
          
          cout << "                        " << 100*(float)nofCorrectlyMatchedAKFNoCut / (float)(nofNotCorrectlyMatchedAKFNoCut+nofNoMatchAKFNoCut) << "% of all events accepted by kinfitter is correctly matched." << endl;
          
          cout << " --- After chi2 cut --- " << endl;
          cout << "Correctly matched reconstructed events (after KF): " << setw(8) << right << nofCorrectlyMatchedAKF << endl;
          cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatchedAKF << endl;
          if ( nofCorrectlyMatchedAKF != 0 || nofNotCorrectlyMatchedAKF != 0 )
            cout << "   ===> This means that " << 100*(float)nofCorrectlyMatchedAKF / (float)(nofCorrectlyMatchedAKF + nofNotCorrectlyMatchedAKF) << "% of matched events is correctly matched after KF." << endl;
          
          cout << "                        " << 100*(float)nofCorrectlyMatchedAKF / (float)nofAcceptedKFit << "% of all events accepted by kinfitter is correctly matched." << endl;
        }
        else cout << "                        " << 100*(float)nofCorrectlyMatched / (float)nofMETCleaned << "% of all events is correctly matched." << endl;
      }
      
      if (doKinFit)
      {
        cout << " --- Kinematic fit for gen events" << endl;
        cout << "Number of generated matched events accepted by kinFitter: " << nofAcceptedKFitMatched << " (" << 100*((float)nofAcceptedKFitMatched/(float)nofHardSelected) << "%)" << endl;
      }
      
      /// Resolution functions
      if (isTTbar && calculateResolutionFunctions)
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
      
    }  // end ! isData
    
//     if (calculateLikelihood)
//     {
//       cout << "Number of events with min in likelihood    " << setw(8) << right << nofGoodEvtsLL[d] << endl;
//       cout << "Number of events without min in likelihood " << setw(8) << right << nofBadEvtsLL[d] << endl;
//       if ( nofGoodEvtsLL[d] != 0 || nofBadEvtsLL[d] != 0 )
//         cout << "   ===> " << 100*(float)nofGoodEvtsLL[d] / (float)(nofGoodEvtsLL[d] + nofBadEvtsLL[d]) << "% are 'good'." << endl;
//       if (! isData) cout << "Number of events with min in parton likelihood    " << setw(8) << right << nofGoodEvtsLL_gen[d] << endl;
//     }
    
    
    if (calculateAverageMass) txtMassReco.close();
    
    tFileMap[dataSetName.c_str()]->Close();
    
    timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
    
  }  // end loop datasets
  
  
  if (calculateAverageMass)
  {
    txtMassGenPMatched.close();
    txtMassGenJMatched.close();
    txtMassRecoCM.close();
    txtMassRecoWMNM.close();
    txtMassRecoNM.close();
    txtMassRecoWM.close();
  }
  
  if (applyWidthSF) txtDebugTopMass.close();
  
  
  cout << "Processing time per dataset: " << endl;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
  }
  
  
  
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
  
  if (makePlots)
  {
    string rootFileName = "ReweightingPlots_"+systStr+".root";
    
    cout << " - Recreate output file ..." << endl;
    TFile *fout = new TFile ((pathOutput+rootFileName).c_str(), "RECREATE");
    cout << "   Output file is " << pathOutput+rootFileName << endl;
    
    ///Write histograms
    fout->cd();
    
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
  }
  
  
  
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

string ConvertDoubleToString(double Number)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

string DotReplace(double var)
{
  string str = ConvertDoubleToString(var);
  replace(str.begin(), str.end(), '.', 'p');
  return str;
}

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

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile.good();
}

void GetMetaData(TTree* tree, bool isData)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  
  tree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
  tree->SetBranchAddress("nEventsSel", &nEventsSel, &b_nEventsSel);
  tree->SetBranchAddress("nofEventsWithGenTop", &nofEventsWithGenTop, &b_nofEventsWithGenTop);
  tree->SetBranchAddress("nofEventsWithGenTopWithStatus22or62", &nofEventsWithGenTopWithStatus22or62, &b_nofEventsWithGenTopWithStatus22or62);
  tree->SetBranchAddress("nofEventsWithGenAntiTop", &nofEventsWithGenAntiTop, &b_nofEventsWithGenAntiTop);
  tree->SetBranchAddress("nofEventsWithGenAntiTopWithStatus22or62", &nofEventsWithGenAntiTopWithStatus22or62, &b_nofEventsWithGenAntiTopWithStatus22or62);
  tree->SetBranchAddress("nofTTEventsWithoutBothGenTops", &nofTTEventsWithoutBothGenTops, &b_nofTTEventsWithoutBothGenTops);
  tree->SetBranchAddress("nofTTEventsWithoutGenTop", &nofTTEventsWithoutGenTop, &b_nofTTEventsWithoutGenTop);
  tree->SetBranchAddress("nofTTEventsWithoutGenAntiTop", &nofTTEventsWithoutGenAntiTop, &b_nofTTEventsWithoutGenAntiTop);
  tree->SetBranchAddress("nofTTEventsWithoutBothGenTopsWithStatus22", &nofTTEventsWithoutBothGenTopsWithStatus22, &b_nofTTEventsWithoutBothGenTopsWithStatus22);
  tree->SetBranchAddress("nofTTEventsWithoutGenTopWithStatus22", &nofTTEventsWithoutGenTopWithStatus22, &b_nofTTEventsWithoutGenTopWithStatus22);
  tree->SetBranchAddress("nofTTEventsWithoutGenAntiTopWithStatus22", &nofTTEventsWithoutGenAntiTopWithStatus22, &b_nofTTEventsWithoutGenAntiTopWithStatus22);
  tree->SetBranchAddress("nofTTEventsWithoutBothGenTopsWithStatus62", &nofTTEventsWithoutBothGenTopsWithStatus62, &b_nofTTEventsWithoutBothGenTopsWithStatus62);
  tree->SetBranchAddress("nofTTEventsWithoutGenTopWithStatus62", &nofTTEventsWithoutGenTopWithStatus62, &b_nofTTEventsWithoutGenTopWithStatus62);
  tree->SetBranchAddress("nofTTEventsWithoutGenAntiTopWithStatus62", &nofTTEventsWithoutGenAntiTopWithStatus62, &b_nofTTEventsWithoutGenAntiTopWithStatus62);
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
  tree->SetBranchAddress("hasErasedBadOrCloneMuon", &hasErasedBadOrCloneMuon, &b_hasErasedBadOrCloneMuon);
  tree->SetBranchAddress("filterHBHENoise", &filterHBHENoise, &b_filterHBHENoise);
  tree->SetBranchAddress("filterHBHEIso", &filterHBHEIso, &b_filterHBHEIso);
  tree->SetBranchAddress("filterCSCTightHalo", &filterCSCTightHalo, &b_filterCSCTightHalo);
  tree->SetBranchAddress("filterEcalDeadCell", &filterEcalDeadCell, &b_filterEcalDeadCell);
  tree->SetBranchAddress("filterEEBadSc", &filterEEBadSc, &b_filterEEBadSc);
  tree->SetBranchAddress("filterBadChCand", &filterBadChCand, &b_filterBadChCand);
  tree->SetBranchAddress("filterBadMuon", &filterBadMuon, &b_filterBadMuon);
  tree->SetBranchAddress("passedMETFilter", &passedMETFilter, &b_passedMETFilter);
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
  tree->SetBranchAddress("jet_nConstituents", jet_nConstituents, &b_jet_nConstituents);
  tree->SetBranchAddress("jet_nChConstituents", jet_nChConstituents, &b_jet_nChConstituents);
  tree->SetBranchAddress("jet_charge", jet_charge, &b_jet_charge);
  tree->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
  tree->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
  tree->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
  tree->SetBranchAddress("jet_E", jet_E, &b_jet_E);
  tree->SetBranchAddress("jet_M", jet_M, &b_jet_M);
  tree->SetBranchAddress("jet_bdiscr", jet_bdiscr, &b_jet_bdiscr);
  tree->SetBranchAddress("met_px", &met_px, &b_met_px);
  tree->SetBranchAddress("met_py", &met_py, &b_met_py);
  tree->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
  tree->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
  tree->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
  tree->SetBranchAddress("met_Et", &met_Et, &b_met_Et);
  tree->SetBranchAddress("met_E", &met_E, &b_met_E);
  tree->SetBranchAddress("met_corr_px", &met_corr_px, &b_met_corr_px);
  tree->SetBranchAddress("met_corr_py", &met_corr_py, &b_met_corr_py);
  tree->SetBranchAddress("met_corr_pt", &met_corr_pt, &b_met_corr_pt);
  tree->SetBranchAddress("met_corr_phi", &met_corr_phi, &b_met_corr_phi);
  tree->SetBranchAddress("met_corr_eta", &met_corr_eta, &b_met_corr_eta);
  tree->SetBranchAddress("met_corr_Et", &met_corr_Et, &b_met_corr_Et);
  tree->SetBranchAddress("met_corr_E", &met_corr_E, &b_met_corr_E);
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
    tree->SetBranchAddress("mc_isLastCopy", mc_isLastCopy, &b_mc_isLastCopy);
    tree->SetBranchAddress("mc_isPromptFinalState", mc_isPromptFinalState, &b_mc_isPromptFinalState);
    tree->SetBranchAddress("mc_isHardProcess", mc_isHardProcess, &b_mc_isHardProcess);
    tree->SetBranchAddress("mc_fromHardProcessFinalState", mc_fromHardProcessFinalState, &b_mc_fromHardProcessFinalState);
    tree->SetBranchAddress("hasGenTop", &hasGenTop, &b_hasGenTop);
    tree->SetBranchAddress("hasGenTopWithStatus22", &hasGenTopWithStatus22, &b_hasGenTopWithStatus22);
    tree->SetBranchAddress("hasGenTopWithStatus62", &hasGenTopWithStatus62, &b_hasGenTopWithStatus62);
    tree->SetBranchAddress("hasGenAntiTop", &hasGenAntiTop, &b_hasGenAntiTop);
    tree->SetBranchAddress("hasGenAntiTopWithStatus22", &hasGenAntiTopWithStatus22, &b_hasGenAntiTopWithStatus22);
    tree->SetBranchAddress("hasGenAntiTopWithStatus62", &hasGenAntiTopWithStatus62, &b_hasGenAntiTopWithStatus62);
  }
}

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  InitHisto1DGen();
  if (makeReweightedPlots) InitHisto1DReweighted();
  if (makeMatchingPlots) InitHisto1DMatch();
  if (makeRecoPlots) InitHisto1DReco();
  
  
  /// SFs
  histo1D["width_SF"] = new TH1F("width_SF", "Scale factor to change the ttbar distribution width; width SF", 500, 0, 5);
  
  
}

void InitHisto1DGen()
{
  TH1::SetDefaultSumw2();
  
  for (int g = 0; g < nGenWidths; g++)
  {
    histo1D[("top_mass_hadr_gen_"+genWidthString[g]).c_str()] = new TH1F(("top_mass_hadr_gen_"+genWidthString[g]).c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 4000, 120, 220);
    histo1D[("top_mass_hadr_gen_"+genWidthString[g]+"_fewerBins").c_str()] = new TH1F(("top_mass_hadr_gen_"+genWidthString[g]+"_fewerBins").c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 800, 120, 220);
    histo1D[("top_mass_lept_gen_"+genWidthString[g]).c_str()] = new TH1F(("top_mass_lept_gen_"+genWidthString[g]).c_str(), "Mass of generated top quark with leptonic decay; M_{t_{lept}} [GeV]", 4000, 120, 220);
    histo1D[("bjj_mass_gen_"+genWidthString[g]).c_str()] = new TH1F(("bjj_mass_gen_"+genWidthString[g]).c_str(), "Mass of generated bqq quarks; M_{bqq} [GeV]", 2000, 50, 300);
    histo1D[("blv_mass_gen_"+genWidthString[g]).c_str()] = new TH1F(("blv_mass_gen_"+genWidthString[g]).c_str(), "Mass of generated b, lepton and neutrino; M_{blv} [GeV]", 2000, 50, 300);
    
    /// Matching
//     histo1D[("top_mass_gen_matched_"+genWidthString[g]).c_str()] = new TH1F(("top_mass_gen_matched_"+genWidthString[g]).c_str(), "Generated top mass of matched events; M_{t} [GeV]", 800, 50, 300);
//     histo1D[("top_mass_reco_matched_"+genWidthString[g]).c_str()] = new TH1F(("top_mass_reco_matched_"+genWidthString[g]).c_str(), "Reconstructed top mass of matched events; M_{t} [GeV]", 400, 0, 400);
//     histo1D[("reduced_top_mass_gen_matched_"+genWidthString[g]).c_str()] = new TH1F(("reduced_top_mass_gen_matched_"+genWidthString[g]).c_str(), "Reduced reconstructed top mass (m_bqq) of matched events; M_{t}/<M_{t}>", 400, 0.8, 1.2);
//     histo1D[("reduced_top_mass_reco_matched_"+genWidthString[g]).c_str()] = new TH1F(("reduced_top_mass_reco_matched_"+genWidthString[g]).c_str(), "Reduced reconstructed top mass (m_bjj) of matched events; M_{t}/<M_{t}>", 400, 0, 2.4);
    
  }
}

void InitHisto1DReweighted()
{
  TH1::SetDefaultSumw2();
  
  for (int s = 0; s < nReweightings; s++)
  {
    histo1D[("top_mass_hadr_gen_"+reweightString[s]).c_str()] = new TH1F(("top_mass_hadr_gen_"+reweightString[s]).c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 4000, 120, 220);
    histo1D[("top_mass_lept_gen_"+reweightString[s]).c_str()] = new TH1F(("top_mass_lept_gen_"+reweightString[s]).c_str(), "Mass of generated top quark with leptonic decay; M_{t_{lept}} [GeV]", 4000, 120, 220);
    histo1D[("bjj_mass_gen_"+reweightString[s]).c_str()] = new TH1F(("bjj_mass_gen_"+reweightString[s]).c_str(), "Mass of generated bqq quarks; M_{bqq} [GeV]", 2000, 50, 300);
    histo1D[("blv_mass_gen_"+reweightString[s]).c_str()] = new TH1F(("blv_mass_gen_"+reweightString[s]).c_str(), "Mass of generated b, lepton and neutrino; M_{blv} [GeV]", 2000, 50, 300);
    
    histo1D[("top_mass_hadr_gen_prodSF_"+reweightString[s]).c_str()] = new TH1F(("top_mass_hadr_gen_prodSF_"+reweightString[s]).c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 4000, 120, 220);
    histo1D[("top_mass_lept_gen_prodSF_"+reweightString[s]).c_str()] = new TH1F(("top_mass_lept_gen_prodSF_"+reweightString[s]).c_str(), "Mass of generated top quark with leptonic decay; M_{t_{lept}} [GeV]", 4000, 120, 220);
    histo1D[("bjj_mass_gen_prodSF_"+reweightString[s]).c_str()] = new TH1F(("bjj_mass_gen_prodSF_"+reweightString[s]).c_str(), "Mass of generated bqq quarks; M_{bqq} [GeV]", 2000, 50, 300);
    histo1D[("blv_mass_gen_prodSF_"+reweightString[s]).c_str()] = new TH1F(("blv_mass_gen_prodSF_"+reweightString[s]).c_str(), "Mass of generated b, lepton and neutrino; M_{blv} [GeV]", 2000, 50, 300);
    
    histo1D[("top_mass_hadr_gen_prodsqrtSF_"+reweightString[s]).c_str()] = new TH1F(("top_mass_hadr_gen_prodsqrtSF_"+reweightString[s]).c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 4000, 120, 220);
    histo1D[("top_mass_lept_gen_prodsqrtSF_"+reweightString[s]).c_str()] = new TH1F(("top_mass_lept_gen_prodsqrtSF_"+reweightString[s]).c_str(), "Mass of generated top quark with leptonic decay; M_{t_{lept}} [GeV]", 4000, 120, 220);
    histo1D[("bjj_mass_gen_prodsqrtSF_"+reweightString[s]).c_str()] = new TH1F(("bjj_mass_gen_prodsqrtSF_"+reweightString[s]).c_str(), "Mass of generated bqq quarks; M_{bqq} [GeV]", 2000, 50, 300);
    histo1D[("blv_mass_gen_prodsqrtSF_"+reweightString[s]).c_str()] = new TH1F(("blv_mass_gen_prodsqrtSF_"+reweightString[s]).c_str(), "Mass of generated b, lepton and neutrino; M_{blv} [GeV]", 2000, 50, 300);
    
    /// Matching
//     histo1D[("top_mass_gen_matched_"+reweightString[s]).c_str()] = new TH1F(("top_mass_gen_matched_"+reweightString[s]).c_str(), "Generated top mass of matched events; M_{t} [GeV]", 800, 50, 300);
//     histo1D[("top_mass_reco_matched_"+reweightString[s]).c_str()] = new TH1F(("top_mass_reco_matched_"+reweightString[s]).c_str(), "Reconstructed top mass of matched events; M_{t} [GeV]", 400, 0, 400);
//     histo1D[("reduced_top_mass_gen_matched_"+reweightString[s]).c_str()] = new TH1F(("reduced_top_mass_gen_matched_"+reweightString[s]).c_str(), "Reduced reconstructed top mass (m_bqq) of matched events; M_{t}/<M_{t}>", 400, 0.8, 1.2);
//     histo1D[("reduced_top_mass_reco_matched_"+reweightString[s]).c_str()] = new TH1F(("reduced_top_mass_reco_matched_"+reweightString[s]).c_str(), "Reduced reconstructed top mass (m_bjj) of matched events; M_{t}/<M_{t}>", 400, 0, 2.4);
    
    /// SFs
    histo1D[("Width_SF_hadr_"+reweightString[s]).c_str()] = new TH1F(("Width_SF_hadr_"+reweightString[s]).c_str(), "Hadronic scale factor to change the ttbar distribution width; width SF", 5001, -0.0005, 5.0005);
    histo1D[("Width_SF_lept_"+reweightString[s]).c_str()] = new TH1F(("Width_SF_lept_"+reweightString[s]).c_str(), "Leptonic scale factor to change the ttbar distribution width; width SF", 5001, -0.0005, 5.0005);
    histo1D[("Width_SF_prod_"+reweightString[s]).c_str()] = new TH1F(("Width_SF_prod_"+reweightString[s]).c_str(), "Product of the hadronic and leptonic scale factor to change the ttbar distribution width; width SF", 5001, -0.0005, 5.0005);
    histo1D[("Width_SF_prodsqrt_"+reweightString[s]).c_str()] = new TH1F(("Width_SF_prodsqrt_"+reweightString[s]).c_str(), "Sqrt of the product of the hadronic and leptonic scale factor to change the ttbar distribution width; width SF", 5001, -0.0005, 5.0005);
  }
}

void InitHisto1DMatch()
{
  TH1::SetDefaultSumw2();
  
  histo1D["matched_W_mass_reco"] = new TH1F("matched_W_mass_reco","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 125, 0, 250);
  histo1D["matched_top_mass_reco"] = new TH1F("matched_top_mass_reco","Reconstructed top mass of matched events; M_{t} [GeV]", 175, 50, 400);
  histo1D["matched_top_mass_gen"] = new TH1F("matched_top_mass_gen","Generated top mass of matched events; M_{t} [GeV]", 1200, 150, 190);
  
  histo1D["matched_red_top_mass_TT_partons"] = new TH1F("matched_red_top_mass_TT_partons","Reduced top mass for matched TT sample (using matched partons); M_{t}/<M_{t}>", 800, 0.8, 1.2);
  histo1D["matched_red_top_mass_TT_jets"] = new TH1F("matched_red_top_mass_TT_jets","Reduced top mass for matched TT sample (using jets from matched partons); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  
  histo1D["matched_mlb_corr"]  = new TH1F("matched_mlb_corr","Reconstructed leptonic top mass using correctly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["matched_mlb_wrong"] = new TH1F("matched_mlb_wrong","Reconstructed leptonic top mass using wrongly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["matched_ttbar_mass_corr"] = new TH1F("matched_ttbar_mass_corr","Reconstructed mass of the top quark pair using correctly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["matched_ttbar_mass_wrong"] = new TH1F("matched_ttbar_mass_wrong","Reconstructed mass of the top quark pair using wrongly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["matched_dR_lep_b_corr"] = new TH1F("matched_dR_lep_b_corr","Delta R between the lepton and the leptonic b jet for matched events; #Delta R(l,b)", 25, 0, 5);
  histo1D["matched_dR_lep_b_wrong"] = new TH1F("matched_dR_lep_b_wrong","Delta R between the lepton and the hadronic b jet for matched events; #Delta R(l,b_{had})", 25, 0, 5);
  
  if (doKinFit) histo1D["matched_top_mass_reco_aKF"] = new TH1F("matched_top_mass_reco_aKF", "Top mass after kinFitter for matched events; m_{t,kf} [GeV]", 80, 0, 400);
}

void InitHisto1DReco()
{
  TH1::SetDefaultSumw2();
  
  /// m_t/<m_t>
  histo1D["red_top_mass_TT_CM"] = new TH1F("red_top_mass_TT_CM","Reduced top mass for matched TT sample (reco, correct top match); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_TT_WMNM"] = new TH1F("red_top_mass_TT_WMNM","Reduced top mass for unmatched TT sample (reco, no & wrong top match); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_TT_NM"] = new TH1F("red_top_mass_TT_NM","Reduced top mass for unmatched TT sample (reco, no top match); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_TT_WM"] = new TH1F("red_top_mass_TT_WM","Reduced top mass for matched TT sample (reco, wrong top match: wrong permutation); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  
  histo1D["red_top_mass_bkgd"] = new TH1F("red_top_mass_bkgd","Reduced top mass for background samples; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_TT"] = new TH1F("red_top_mass_TT","Reduced top mass for TT sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_ST_tW_top"] = new TH1F("red_top_mass_ST_tW_top","Reduced top mass for ST tW top sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_ST_tW_antitop"] = new TH1F("red_top_mass_ST_tW_antitop","Reduced top mass for ST tW antitop sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_ST_t_top"] = new TH1F("red_top_mass_ST_t_top","Reduced top mass for ST t top sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_ST_t_antitop"] = new TH1F("red_top_mass_ST_t_antitop","Reduced top mass for ST t antitop sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_DYJets"] = new TH1F("red_top_mass_DYJets","Reduced top mass for DY+Jets sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_WJets"] = new TH1F("red_top_mass_WJets","Reduced top mass for W+Jets sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["red_top_mass_data"] = new TH1F("red_top_mass_data","Reduced top mass for data sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  
  /// mlb
  histo1D["mlb_CM"]  = new TH1F("minMlb_CM","m_{lb} mass for events that have a correct hadronic top match (CM); m_{lb} [GeV]", 400, 0, 800);
  histo1D["mlb_WM"]  = new TH1F("minMlb_WM","m_{lb} for events that have a wrong hadronic top match (WM); m_{lb} [GeV]", 400, 0, 800);
  histo1D["mlb_NM"]  = new TH1F("minMlb_NM","m_{lb} for events that have no hadronic top match (NM); m_{lb} [GeV]", 400, 0, 800);
  
  /// dR
  //  lepton, b(lep)
  histo1D["dR_lep_b_lep_CM"]  = new TH1F("dR_lep_b_lep_CM","Minimal delta R between the lepton and the leptonic b jet (reco, correct match); #Delta R(l,b_{l})", 25, 0, 5);
  histo1D["dR_lep_b_lep_WM"]  = new TH1F("dR_lep_b_lep_WM","Minimal delta R between the lepton and the leptonic b jet (reco, wrong permutation); #Delta R(l,b_{l})", 25, 0, 5);
  histo1D["dR_lep_b_lep_NM"]  = new TH1F("dR_lep_b_lep_NM","Minimal delta R between the lepton and the leptonic b jet (reco, no match); #Delta R(l,b_{l})", 25, 0, 5);
  //  lepton, b(hadr)
  histo1D["dR_lep_b_had_CM"]  = new TH1F("dR_lep_b_had_CM","Minimal delta R between the lepton and the hadronic b jet (reco, correct match); #Delta R(l,b_{h})", 25, 0, 5);
  histo1D["dR_lep_b_had_WM"]  = new TH1F("dR_lep_b_had_WM","Minimal delta R between the lepton and the hadronic b jet (reco, wrong permutation); #Delta R(l,b_{h})", 25, 0, 5);
  histo1D["dR_lep_b_had_NM"]  = new TH1F("dR_lep_b_had_NM","Minimal delta R between the lepton and the hadronic b jet (reco, no match); #Delta R(l,b_{h})", 25, 0, 5);
  
  /// ttbar mass
  histo1D["ttbar_mass_CM"] = new TH1F("ttbar_mass_CM","Reconstructed mass of the top quark pair (reco, correct match); m_{t#bar{t}} [GeV]", 500, 0, 1000);
  histo1D["ttbar_mass_WM"] = new TH1F("ttbar_mass_WM","Reconstructed mass of the top quark pair (reco, wrong permutation); m_{t#bar{t}} [GeV]", 500, 0, 1000);
  histo1D["ttbar_mass_NM"] = new TH1F("ttbar_mass_NM","Reconstructed mass of the top quark pair (reco, no match); m_{t#bar{t}} [GeV]", 500, 0, 1000);
  
  
  /// KinFitter
  if (doKinFit)
  {
    histo1D["KF_Chi2_TT"] = new TH1F("KF_Chi2_TT", "Chi2 value of kinFitter (TT); #chi^{2}", 200, 0, 20);
    histo1D["KF_W_mass_orig_TT"] = new TH1F("KF_W_mass_orig_TT", "W mass before kinFitter (TT); m_{W} [GeV]", 250, 0, 500);
    histo1D["KF_top_mass_orig_TT"] = new TH1F("KF_top_mass_orig_TT", "Top mass before kinFitter (TT); m_{t} [GeV]", 400, 0, 800);
    histo1D["KF_top_pt_orig_TT"] = new TH1F("KF_top_pt_orig_TT", "Top pt before kinFitter (TT); p_{T} [GeV]", 400, 0, 800);
    histo1D["KF_dR_lep_b_orig_TT"] = new TH1F("KF_dR_lep_b_orig_TT","Minimal delta R between the lepton and a b jet before kinFitter (for reconstructed TT events); #Delta R(l,b)", 25, 0, 5);
    histo1D["KF_W_mass_corr_TT"] = new TH1F("KF_W_mass_corr_TT", "W mass after kinFitter (TT); m_{W,kf} [GeV]", 250, 0, 500);
    histo1D["KF_top_mass_corr_TT"] = new TH1F("KF_top_mass_corr_TT", "Top mass after kinFitter (TT); m_{t,kf} [GeV]", 400, 0, 800);
    histo1D["KF_top_pt_corr_TT"] = new TH1F("KF_top_pt_corr_TT", "Top pt after kinFitter (TT); p_{T,kf} [GeV]", 400, 0, 800);
    histo1D["KF_dR_lep_b_corr_TT"] = new TH1F("KF_dR_lep_b_corr_TT","Minimal delta R between the lepton and a b jet after kinFitter (for reconstructed TT events); #Delta R(l,b)", 25, 0, 5);
    
    histo1D["KF_top_mass_orig_ex4j_chi2cut5_TT"] = new TH1F("KF_top_mass_orig_ex4j_chi2cut5_TT", "Top mass before kinFitter (TT) (exactly 4 jets - KF chi2 < 5); m_{t} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_ex4j_chi2cut5_TT"] = new TH1F("KF_top_mass_corr_ex4j_chi2cut5_TT", "Top mass after kinFitter (TT) (exactly 4 jets - KF chi2 < 5); m_{t,kf} [GeV]", 400, 0, 800);
    
    histo1D["KF_top_mass_orig_ex4j_chi2cut2_TT"] = new TH1F("KF_top_mass_orig_ex4j_chi2cut2_TT", "Top mass before kinFitter (TT) (exactly 4 jets - KF chi2 < 2); m_{t} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_ex4j_chi2cut2_TT"] = new TH1F("KF_top_mass_corr_ex4j_chi2cut2_TT", "Top mass after kinFitter (TT) (exactly 4 jets - KF chi2 < 2); m_{t,kf} [GeV]", 400, 0, 800);
    
    histo1D["KF_top_mass_orig_ex4j_chi2cut1p5_TT"] = new TH1F("KF_top_mass_orig_ex4j_chi2cut1p5_TT", "Top mass before kinFitter (TT) (exactly 4 jets - KF chi2 < 1.5); m_{t} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_ex4j_chi2cut1p5_TT"] = new TH1F("KF_top_mass_corr_ex4j_chi2cut1p5_TT", "Top mass after kinFitter (TT) (exactly 4 jets - KF chi2 < 1.5); m_{t,kf} [GeV]", 400, 0, 800);
    
    histo1D["KF_jet0_Et_diff_TT"] = new TH1F("KF_jet0_Et_diff_TT", "Et difference after kinFitter for jet0 (TT); E_{T,kf} - E_{T,orig} [GeV]", 400, -50, 50);
    histo1D["KF_jet1_Et_diff_TT"] = new TH1F("KF_jet1_Et_diff_TT", "Et difference after kinFitter for jet1 (TT); E_{T,kf} - E_{T,orig} [GeV]", 400, -50, 50);
    
    histo1D["KF_top_mass_corr_CM"] = new TH1F("KF_top_mass_corr_CM", "Top mass after kinFitter for correct match (CM); m_{t,kf} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_WM"] = new TH1F("KF_top_mass_corr_WM", "Top mass after kinFitter for wrong permutations (WM); m_{t,kf} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_NM"] = new TH1F("KF_top_mass_corr_NM", "Top mass after kinFitter for no match (NM); m_{t,kf} [GeV]", 400, 0, 800);
    
    histo1D["KF_Chi2_CM"] = new TH1F("KF_Chi2_CM", "Chi2 value of kinFitter (CM); #chi^{2}", 200, 0, 5);
    histo1D["KF_Chi2_WM"] = new TH1F("KF_Chi2_WM", "Chi2 value of kinFitter (WM); #chi^{2}", 200, 0, 5);
    histo1D["KF_Chi2_NM"] = new TH1F("KF_Chi2_NM", "Chi2 value of kinFitter (NM); #chi^{2}", 200, 0, 5);
    histo1D["KF_Chi2_CM_wide"] = new TH1F("KF_Chi2_CM_wide", "Chi2 value of kinFitter (CM); #chi^{2}", 200, 0, 20);
    histo1D["KF_Chi2_WM_wide"] = new TH1F("KF_Chi2_WM_wide", "Chi2 value of kinFitter (WM); #chi^{2}", 200, 0, 20);
    histo1D["KF_Chi2_NM_wide"] = new TH1F("KF_Chi2_NM_wide", "Chi2 value of kinFitter (NM); #chi^{2}", 200, 0, 20);
  }
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
  InitHisto2DGen();
  if (makeReweightedPlots) InitHisto2DReweighted();
  if (makeMatchingPlots) InitHisto2DMatch();
  if (makeRecoPlots)
  {
    histo2D["dR_lep_b_lep_vs_had_CM"] = new TH2F("dR_lep_b_lep_vs_had_CM","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, correct match); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
    histo2D["dR_lep_b_lep_vs_had_WM"] = new TH2F("dR_lep_b_lep_vs_had_WM","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, wrong permutations); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
    histo2D["dR_lep_b_lep_vs_had_NM"] = new TH2F("dR_lep_b_lep_vs_had_NM","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, no match); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
    
    histo2D["ttbar_mass_vs_minMlb_CM"] = new TH2F("ttbar_mass_vs_minMlb_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
    histo2D["ttbar_mass_vs_minMlb_WM"] = new TH2F("ttbar_mass_vs_minMlb_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
    histo2D["ttbar_mass_vs_minMlb_NM"] = new TH2F("ttbar_mass_vs_minMlb_NM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
    
    /// KinFitter
    if (doKinFit)
    {
      histo2D["KF_W_mass_orig_vs_corr_TT"] = new TH2F("KF_W_mass_orig_vs_corr_TT", "W mass made with KF corrected jets vs. original jets (TT); m_{W,orig} [GeV]; m_{W,kf} [GeV]", 250, 0, 500, 250, 0, 500);
      histo2D["KF_top_mass_orig_vs_corr_TT"] = new TH2F("KF_top_mass_orig_vs_corr_TT", "Top mass made with KF corrected jets vs. original jets (TT); m_{t,orig} [GeV]; m_{t,kf} [GeV]", 400, 0, 800, 400, 0, 800);
      histo2D["KF_jets_Et_diff_TT"] = new TH2F("KF_jets_Et_diff_TT", "Et difference after kinFitter for jet1 in function of jet0 (TT); E_{T,0,kf} - E_{T,0,orig} [GeV]; E_{T,1,kf} - E_{T,1,orig} [GeV]", 400, -50, 50, 400, -50, 50);
      histo2D["KF_W_px_orig_vs_corr_TT"] = new TH2F("KF_W_px_orig_vs_corr_TT", "W p_{x} made with KF corrected jets vs. original jets (TT); p_{x,orig} [GeV]; p_{x,kf} [GeV]", 400, 0, 800, 400, 0, 800);
      histo2D["KF_W_py_orig_vs_corr_TT"] = new TH2F("KF_W_py_orig_vs_corr_TT", "W p_{y} made with KF corrected jets vs. original jets (TT); p_{y,orig} [GeV]; p_{y,kf} [GeV]", 400, 0, 800, 400, 0, 800);
      histo2D["KF_W_pz_orig_vs_corr_TT"] = new TH2F("KF_W_pz_orig_vs_corr_TT", "W p_{z} made with KF corrected jets vs. original jets (TT); p_{z,orig} [GeV]; p_{z,kf} [GeV]", 400, 0, 800, 400, 0, 800);
    }
  }
}

void InitHisto2DGen()
{
  TH2::SetDefaultSumw2();
  
  for (int g = 0; g < nGenWidths; g++)
  {
    histo2D["top_mass_hadr_lept_gen_"+genWidthString[g]] = new TH2F(("top_mass_hadr_lept_gen_"+genWidthString[g]).c_str(),("Mass of generated top quark with leptonic decay vs. hadronic decay ("+genWidthString[g]+"); M_{t_{hadr}} [GeV]; M_{t_{lept}} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
    histo2D["top_mass_hadr_bjj_gen_"+genWidthString[g]] = new TH2F(("top_mass_hadr_bjj_gen_"+genWidthString[g]).c_str(),("Generated mass of bjj quarks vs. hadronically decaying top quark mass ("+genWidthString[g]+"); M_{t_{hadr}} [GeV]; M_{bjj} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
    histo2D["top_mass_lept_blv_gen_"+genWidthString[g]] = new TH2F(("top_mass_lept_blv_gen_"+genWidthString[g]).c_str(),("Mass of generated b, lepton and neutrino vs. leptonically decaying top quark mass ("+genWidthString[g]+"); M_{t_{lept}} [GeV]; M_{blv} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
    histo2D["top_mass_bjj_blv_gen_"+genWidthString[g]] = new TH2F(("top_mass_bjj_blv_gen_"+genWidthString[g]).c_str(),("Mass of generated b, lepton and neutrino vs. mass of bjj quarks ("+genWidthString[g]+"); M_{bjj} [GeV]; M_{blv} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
  }
}

void InitHisto2DReweighted()
{
  TH2::SetDefaultSumw2();
  
  for (int s = 0; s < nReweightings; s++)
  {
    histo2D[("Width_SF_hadr_lept_"+reweightString[s]).c_str()] = new TH2F(("Width_SF_hadr_lept_"+reweightString[s]).c_str(), "Leptonic vs. hadronic scale factor to change the ttbar distribution width; SF_{h}; SF_{l}", 5001, -0.0005, 5.0005, 5001, -0.0005, 5.0005);
    histo2D[("Width_SF_hadr_prod_"+reweightString[s]).c_str()] = new TH2F(("Width_SF_hadr_prod_"+reweightString[s]).c_str(), "Product of hadronic and leptonic SF vs. hadronic scale factor to change the ttbar distribution width; SF_{h}; SF_{h} #times SF_{l}", 5001, -0.0005, 5.0005, 5001, -0.0005, 5.0005);
    histo2D[("Width_SF_hadr_prodsqrt_"+reweightString[s]).c_str()] = new TH2F(("Width_SF_hadr_prodsqrt99_"+reweightString[s]).c_str(), "Sqrt of the product of hadronic and leptonic SF vs. hadronic scale factor to change the ttbar distribution width; SF_{h}; #sqrt{SF_{h} #times SF_{l}}", 5001, -0.0005, 5.0005, 5001, -0.0005, 5.0005);
    
    histo2D[("top_mass_hadr_gen_vs_weight_"+reweightString[s]).c_str()] = new TH2F(("top_mass_hadr_gen_vs_weight_"+reweightString[s]).c_str(), "Weights vs. mass of generated top quark (hadronic decay); M_{t_{hadr}} [GeV]; SF_{h}", 1000, 120, 220, 5001, -0.0005, 5.0005);
    //histo2D[("top_mass_lept_gen_vs_weight_"+reweightString[s]).c_str()] = new TH2F(("top_mass_lept_gen_vs_weight_"+reweightString[s]).c_str(), "Weights vs. mass of generated top quark (leptonic decay); M_{t_{lept}} [GeV]; SF_{l}", 1000, 120, 220, 5001, -0.0005, 5.0005);
  }
  }

void InitHisto2DMatch()
{
  TH2::SetDefaultSumw2();
  
  /// Matched events
  histo2D["matched_mlb_corr_mlb_wrong"] = new TH2F("matched_mlb_corr_mlb_wrong","Wrongly constructed M_{lb} vs. correctly constructed M_{lb}; M_{lb_{lep}}; M_{lb_{had}}", 80, 0, 800, 80, 0, 800);
  histo2D["matched_dR_lep_b_corr_dR_lep_b_wrong"] = new TH2F("matched_dR_lep_b_corr_dR_lep_b_wrong","Wrongly constructed dR(l,b) vs. correctly constructed dR(l,b); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["matched_mlb_dR_lep_b_corr"] = new TH2F("matched_mlb_dR_lep_b_corr","dR(l,b) vs. M_{lb}; M_{lb_{lep}}; #Delta R(l,b_{lep})", 80, 0, 800, 25, 0, 5);
  histo2D["matched_mlb_dR_lep_b_wrong"] = new TH2F("matched_mlb_dR_lep_b_wrong","dR(l,b) vs. M_{lb}, both wrongly matched; M_{lb_{had}}; #Delta R(l,b_{had})", 80, 0, 800, 25, 0, 5);
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
  nofEventsWithGenTop = 0;
  nofEventsWithGenTopWithStatus22or62 = 0;
  nofEventsWithGenAntiTop = 0;
  nofEventsWithGenAntiTopWithStatus22or62 = 0;
  nofTTEventsWithoutBothGenTops = 0;
  nofTTEventsWithoutGenTop = 0;
  nofTTEventsWithoutGenAntiTop = 0;
  nofTTEventsWithoutBothGenTopsWithStatus22 = 0;
  nofTTEventsWithoutGenTopWithStatus22 = 0;
  nofTTEventsWithoutGenAntiTopWithStatus22 = 0;
  nofTTEventsWithoutBothGenTopsWithStatus62 = 0;
  nofTTEventsWithoutGenTopWithStatus62 = 0;
  nofTTEventsWithoutGenAntiTopWithStatus62 = 0;

  
  strSyst = "";
  eqLumi = 1.;
  lumiWeight = 1.;
  
  nofHardSelected = 0;
  nofMETCleaned = 0;
  nofMatchedEvents = 0;
  nofHadrMatchedEvents = 0;
  nofHadrMatchedEventsAKF = 0;
  nofCorrectlyMatched = 0;
  nofNotCorrectlyMatched = 0;
  nofCorrectlyMatchedAKF = 0;
  nofNotCorrectlyMatchedAKF = 0;
  nofCorrectlyMatchedAKFNoCut = 0;
  nofNotCorrectlyMatchedAKFNoCut = 0;
  nofNoMatchAKFNoCut = 0;
  nofAcceptedKFit = 0;
  nofAcceptedKFitMatched = 0;
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
  hasExactly4Jets = false;
  hasJetLeptonCleaning = false;
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
    jet_nConstituents[i] = -1;
    jet_nChConstituents[i] = -1;
    jet_pt[i] = 0.;
    jet_phi[i] = 0.;
    jet_eta[i] = 0.;
    jet_E[i] = 0.;
    jet_M[i] = 0.;
    jet_bdiscr[i] = -1.;
  }
  met_px = 0.;
  met_py = 0.;
  met_pt = 0.;
  met_phi = 0.;
  met_eta = 0.;
  met_Et = 0.;
  met_E = 0.;
  met_corr_px = 0.;
  met_corr_py = 0.;
  met_corr_pt = 0.;
  met_corr_phi = 0.;
  met_corr_eta = 0.;
  met_corr_Et = 0.;
  met_corr_E = 0.;
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
    hasGenTop = false;
    hasGenTopWithStatus22 = false;
    hasGenTopWithStatus62 = false;
    hasGenAntiTop = false;
    hasGenAntiTopWithStatus22 = false;
    hasGenAntiTopWithStatus62 = false;
  }
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
  selectedJetsAKF.clear();
  selectedJetsKFcorrected.clear();
  selectedJetsKFMatched.clear();
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
  foundAllTTbarComponents = false;
  for (int i = 0; i < 3; i++)
  {
    bjjDecay[i] = -9999;
    blvDecay[i] = -9999;
  }
  partonId.clear();
  
  matched_W_mass_q = -1.;
  matched_top_mass_q = -1.;
  matched_W_mass_j = -1.;
  matched_top_mass_j = -1.;
  matched_mlb_corr = -1.;
  matched_ttbarMass_corr = -1.;
  matched_dR_lep_b_corr = 999.;
  matched_mlb_wrong = -1.;
  matched_ttbarMass_wrong = -1.;
  matched_dR_lep_b_wrong = 999.;
}

void ClearVars()
{
  ClearMatching();
  
  scaleFactor = 1.;
  widthSF = 1.;
  topPtRewSF = 1.;
  topPtSF = 1.;
  antiTopPtSF = 1.;
  m_top = -1.;
  m_antitop = -1.;
  m_bjj = -1.;
  m_blv = -1.;
  m_hadr = -1;
  m_lept = -1;
  evWeight_hadr = 1.;
  evWeight_lept = 1.;
  evWeight_prod = 1.;
  evWeight_prodsqrt = 1.;
  bJetId.clear();
  bdiscrTop = -99.;
  bdiscrTop2 = -99.;
  tempbdiscr = -99.;
  labelB1 = -9999;
  labelB2 = -9999;
  for (int i = 0; i < 4; i++)
  {
    labelsReco[i] = -9999;
  }
  massForWidth = 0.01;
  catSuffix = "";
  isCM = false;
  isWM = false;
  isNM = false;
  doneKinFit = false;
  kFitVerbosity = false;
  kFitChi2 = 99.;
  kFitChi2Matched = 99.;
  passKFChi2MatchedCut = false;
  
  M3 = -1.;
  M3_aKF = -1.;
  Ht = -1.;
  Ht_aKF = -1;
  min_Mlb = 9999.;
  dRLepB = -1.;
  reco_W_mass_bKF = -1.;
  reco_top_mass_bKF = -1.;
  reco_top_pt_bKF = -1.;
  reco_mlb_bKF = -1.;
  reco_dRLepB_lep_bKF = -1.;
  reco_dRLepB_had_bKF = -1.;
  reco_ttbar_mass_bKF = -1.;
  redTopMass_bKF = -1.;
  reco_W_mass_aKF = -1.;
  reco_top_mass_aKF = -1.;
  reco_top_pt_aKF = -1.;
  reco_mlb_aKF = -1.;
  reco_dRLepB_lep_aKF = -1.;
  reco_dRLepB_had_aKF = -1.;
  reco_ttbar_mass_aKF = -1.;
  redTopMass = -1.;

}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearVars();
}

void FillMatchingPlots()
{
  if (hadronicTopJetsMatched)
  {  
    histo1D["matched_W_mass_reco"]->Fill(matched_W_mass_j, widthSF);
    histo1D["matched_top_mass_reco"]->Fill(matched_top_mass_j, widthSF);
    histo1D["matched_top_mass_gen"]->Fill(matched_top_mass_q, widthSF);
    
    histo1D["matched_red_top_mass_TT_partons"]->Fill(matched_top_mass_q/aveTopMass[0], widthSF);
    histo1D["matched_red_top_mass_TT_jets"]->Fill(matched_top_mass_j/aveTopMass[1], widthSF);
    
    if (doKinFit)
    {
      histo1D["matched_top_mass_reco_aKF"]->Fill(matched_top_mass_j_akF, widthSF);
    }
    
    if (all4PartonsMatched && muonmatched)
    {
      matched_mlb_corr = (selectedLepton[0] + jetsMatched[3]).M();  // lept b
      matched_mlb_wrong = (selectedLepton[0] + jetsMatched[2]).M();  // hadr b
      matched_ttbarMass_corr = matched_top_mass_j + matched_mlb_corr;
      matched_ttbarMass_wrong = matched_top_mass_j + matched_mlb_wrong;
      matched_dR_lep_b_corr = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[3]);  // lept b
      matched_dR_lep_b_wrong = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[2]);  // hadr b
      
      histo1D["matched_mlb_corr"]->Fill(matched_mlb_corr);
      histo1D["matched_mlb_wrong"]->Fill(matched_mlb_wrong);
      histo1D["matched_ttbar_mass_corr"]->Fill(matched_ttbarMass_corr);
      histo1D["matched_ttbar_mass_wrong"]->Fill(matched_ttbarMass_wrong);
      histo1D["matched_dR_lep_b_corr"]->Fill(matched_dR_lep_b_corr);
      histo1D["matched_dR_lep_b_wrong"]->Fill(matched_dR_lep_b_wrong);
      
      histo2D["matched_mlb_corr_mlb_wrong"]->Fill(matched_mlb_corr, matched_mlb_wrong);
      histo2D["matched_dR_lep_b_corr_dR_lep_b_wrong"]->Fill(matched_dR_lep_b_corr, matched_dR_lep_b_wrong);
      histo2D["matched_mlb_dR_lep_b_corr"]->Fill(matched_mlb_corr, matched_dR_lep_b_corr);
      histo2D["matched_mlb_dR_lep_b_wrong"]->Fill(matched_mlb_wrong, matched_dR_lep_b_wrong);
    }  // end all4PartonsMatched && muonMatched
  }  // end hadronicTopJetsMatched
}

void FillKinFitPlots(bool doneKinFit)
{
  if (isTTbar)
  {
    if (! doneKinFit)
    {
      histo1D["KF_W_mass_orig_TT"]->Fill(reco_W_mass_bKF);
      histo1D["KF_top_mass_orig_TT"]->Fill(reco_top_mass_bKF);
      histo1D["KF_top_pt_orig_TT"]->Fill(reco_top_pt_bKF);
      histo1D["KF_dR_lep_b_orig_TT"]->Fill(reco_dRLepB_lep_bKF);
    }
    else
    {
      histo1D["KF_W_mass_corr_TT"]->Fill(reco_W_mass_bKF);
      histo2D["KF_W_mass_orig_vs_corr_TT"]->Fill(reco_W_mass_bKF, reco_W_mass_aKF);
      histo2D["KF_W_px_orig_vs_corr_TT"]->Fill( (selectedJets[0] + selectedJets[1]).Px(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Px() );
      histo2D["KF_W_py_orig_vs_corr_TT"]->Fill( (selectedJets[0] + selectedJets[1]).Py(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Py() );
      histo2D["KF_W_pz_orig_vs_corr_TT"]->Fill( (selectedJets[0] + selectedJets[1]).Pz(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Pz() );
      
      histo1D["KF_top_mass_corr_TT"]->Fill(reco_top_mass_aKF);
      histo2D["KF_top_mass_orig_vs_corr_TT"]->Fill(reco_top_mass_bKF, reco_top_mass_aKF);
      histo1D["KF_top_pt_corr_TT"]->Fill(reco_top_pt_aKF);
      histo1D["KF_dR_lep_b_corr_TT"]->Fill(reco_dRLepB_lep_aKF);
      
      histo1D["KF_jet0_Et_diff_TT"]->Fill((selectedJetsKFcorrected[0] - selectedJets[labelsReco[0]]).Et());
      histo1D["KF_jet1_Et_diff_TT"]->Fill((selectedJetsKFcorrected[1] - selectedJets[labelsReco[1]]).Et());
      histo2D["KF_jets_Et_diff_TT"]->Fill((selectedJetsKFcorrected[0] - selectedJets[labelsReco[0]]).Et(), (selectedJetsKFcorrected[1] - selectedJets[labelsReco[1]]).Et());
    }
    
    if ( kFitChi2 < 5. )
    {
      if (! doneKinFit) histo1D["KF_top_mass_orig_ex4j_chi2cut5_TT"]->Fill(reco_top_mass_bKF);
      else histo1D["KF_top_mass_corr_ex4j_chi2cut5_TT"]->Fill(reco_top_mass_aKF);
      
      if ( kFitChi2 < 2. )
      {
        if (! doneKinFit) histo1D["KF_top_mass_orig_ex4j_chi2cut2_TT"]->Fill(reco_top_mass_bKF);
        else histo1D["KF_top_mass_corr_ex4j_chi2cut2_TT"]->Fill(reco_top_mass_aKF);
        
        if ( kFitChi2 < 1.5 )
        {
          if (! doneKinFit) histo1D["KF_top_mass_orig_ex4j_chi2cut1p5_TT"]->Fill(reco_top_mass_bKF);
          else histo1D["KF_top_mass_corr_ex4j_chi2cut1p5_TT"]->Fill(reco_top_mass_aKF);
        }  // 1.5
      }  // 2
    }  // 5
  }
}

void FillCatsPlots(string catSuffix)
{
  if (! isData)
  {
    histo1D["red_top_mass_bkgd"]->Fill(redTopMass);
    if ( isWM || isNM ) histo1D["red_top_mass_TT_WMNM"]->Fill(redTopMass, widthSF);
    histo1D[("red_top_mass_TT"+catSuffix).c_str()]->Fill(redTopMass, widthSF);
    histo1D[("mlb"+catSuffix).c_str()]->Fill(reco_mlb_aKF, widthSF);
    histo1D[("dR_lep_b_lep"+catSuffix).c_str()]->Fill(reco_dRLepB_lep_aKF, widthSF);
    histo1D[("dR_lep_b_had"+catSuffix).c_str()]->Fill(reco_dRLepB_had_aKF, widthSF);
    histo1D[("ttbar_mass"+catSuffix).c_str()]->Fill(reco_ttbar_mass_aKF, widthSF);
    histo2D[("dR_lep_b_lep_vs_had"+catSuffix).c_str()]->Fill(reco_dRLepB_lep_aKF, reco_dRLepB_had_aKF, widthSF);
    histo2D[("ttbar_mass_vs_minMlb"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//    if ( reco_dRLepB_lep < 3. )
//    {
//      histo2D[("ttbar_mass_vs_minMlb_dRlepCut"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//      if ( reco_dRLepB_had > 1.2 )
//        histo2D[("ttbar_mass_vs_minMlb_dRBothCuts"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//      
//      if ( reco_dRLepB_lep < 2. )
//      {
//        histo2D[("ttbar_mass_vs_minMlb_dRlepCutHard"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//        if ( reco_dRLepB_had > 2. )
//          histo2D[("ttbar_mass_vs_minMlb_dRBothCutsHard"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//      }
//    }
//    if ( reco_dRLepB_had > 1.2 )
//    {
//      histo2D[("ttbar_mass_vs_minMlb_dRhadCut"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//      if ( reco_dRLepB_had > 2. )
//        histo2D[("ttbar_mass_vs_minMlb_dRhadCutHard"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//    }
    if (doKinFit)
    {
      histo1D["KF_Chi2"+catSuffix]->Fill(kFitChi2);
      histo1D["KF_top_mass_corr"+catSuffix]->Fill(reco_top_mass_aKF, widthSF);
    }
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

void GetEraFraction(double* fractions)
{
  TFile* file = new TFile((pathNtuplesData+"Ntuples_data.root").c_str(),"READ");
  TTree* fChain = (TTree*) file->Get("stats");
  
  //fChain->SetBranchAddress("nofEventsHLTv2", &nofEventsHLTv2, &b_nofEventsHLTv2);
  //fChain->SetBranchAddress("nofEventsHLTv3", &nofEventsHLTv3, &b_nofEventsHLTv3);
  
  long nBCDEF = GetNEvents(fChain, "nofEventsRunB", 1) + GetNEvents(fChain, "nofEventsRunCD", 1) + GetNEvents(fChain, "nofEventsRunEF", 1);
  long nGH = GetNEvents(fChain, "nofEventsRunG", 1) + GetNEvents(fChain, "nofEventsRunH", 1);
  
  fractions[0] = ((double)nBCDEF) / ((double)nBCDEF + (double)nGH);
  fractions[1] = ((double)nGH) / ((double)nBCDEF + (double)nGH);
  
  file->Close();
  ClearMetaData();
}

void PrintKFDebug(int ievt)
{
  cout << endl <<"Event " << setw(5) << right << ievt << "   ";
  cout << "Top mass after kinFit is negative! I.e. " << reco_top_mass_aKF << " Before kinFit: " << reco_top_mass_bKF << "  ievt " << ievt << endl;
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

