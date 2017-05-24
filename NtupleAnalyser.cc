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
#include "Tools/interface/KinFitter.h"
#include "Tools/interface/EventReweighting.h"
#include "Tools/interface/Likelihood.h"


using namespace std;
using namespace TopTree;


bool test = false;
bool testHistos = false;
bool testTTbarOnly = false;
bool doGenOnly = false;
bool makePlots = true;
bool calculateResolutionFunctions = false;
bool calculateAverageMass = false;
bool makeTGraphs = false;
bool calculateLikelihood = true;
bool doPseudoExps = true;
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


bool applyWidthSF = false;
double scaleWidth = 1.;

string systStr = "nominal";
pair<string,string> whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos )
  {
    return pair<string,string>("170517","170522");
  }
  else if ( syst.find("JECup") != std::string::npos ) return pair<string,string>("","170522");
  else if ( syst.find("JECdown") != std::string::npos ) return pair<string,string>("","170522");
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    return pair<string,string>("170517","170522");
  }
}
pair<string,string> ntupleDate = whichDate(systStr);
int verbose = 2;

string pathNtuples = "";
string pathNtuplesMC = "";
string pathNtuplesData = "";
string outputDirLL = "LikelihoodTemplates/";
string inputDirLL = "";
string inputDateLL = "170519_1140/";
bool isData = false;
bool isTTbar = false;

int nofHardSelected = 0;
int nofMETCleaned = 0;
int nofMatchedEvents = 0;
int nofHadrMatchedEvents = 0;
int nofCorrectlyMatched = 0;
int nofNotCorrectlyMatched = 0;
int nofCP = 0, nofWP = 0, nofUP = 0;
int nofCP_TT = 0, nofWP_TT = 0, nofUP_TT = 0;
double nofCP_weighted = 0, nofWP_weighted = 0, nofUP_weighted = 0;

/// Lumi per data era
double lumi_runBCDEF = 19.67550334113;  // 1/fb
double lumi_runGH = 16.146177597883;  // 1/fb
double Luminosity = (lumi_runBCDEF + lumi_runGH)*1000;  // 1/pb

///  Working points for b tagging  // Updated 13/04/17, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
double CSVv2Loose  = 0.5426;
double CSVv2Medium = 0.8484;
double CSVv2Tight  = 0.9535;

/// Average top mass
// TT_genp_match, TT_genj_match, TT_reco_match, TT_reco_wrongMatch_WP/UP, TT_reco_noMatch, TT_reco_wrongPerm, TT_reco, ST_t_top, ST_t_antitop, ST_tW_top, ST_tW_antitop, DYJets, WJets, data, Reco, All, MC, Reco, All, Samples
// also background in CP/WP/UP cats (unlike name suggests)
const int nofAveMasses = 16;
//  KF chi2 < 5
std::array<double, nofAveMasses> aveTopMass = {171.825, 169.751, 167.558, 193.812, 191.823, 198.213, 182.175, 252.030, 252.837, 233.056, 230.495, 214.340, 212.065, 185.199, 182.561, 182.876};
//  no KF chi2 cut
//std::array<double, nofAveMasses> aveTopMass = {171.810, 168.728, 167.110, 203.721, 204.952, 198.233, 193.403, 270.895, 267.167, 230.144, 229.649, 250.010, 242.091, 200.455, 193.963, 194.025};

/// # events after kin fitter
//  KF chi2 < 5
int nEventsAKF[] = {106005, 573511, 1401, 1366, 653, 412, 158, 49};

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;
map<string,TGraph*> graph;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;

vector < Dataset* > datasets, datasetsMSP, datasetsTemp;

ofstream txtMassGenPMatched, txtMassGenJMatched, txtMassRecoCP, txtMassRecoWPUP, txtMassRecoUP, txtMassRecoWP, txtMassReco;

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
void FillGeneralPlots(vector<Dataset *> datasets, int d);
void FillMatchingPlots();
void FillKinFitPlots();
void FillCatsPlots(string catSuffix);
void FillMSPlots(int d);
long GetNEvents(TTree* fChain, string var, bool isData);
void GetEraFraction(double* fractions);
void CheckSystematics(vector<int> vJER, vector<int> vJES, vector<int> vPU);
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
Bool_t          isDataRunB;
Bool_t          isDataRunC;
Bool_t          isDataRunD;
Bool_t          isDataRunE;
Bool_t          isDataRunF;
Bool_t          isDataRunG;
Bool_t          isDataRunH;
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
Double_t        met_pt;
Double_t        met_phi;
Double_t        met_eta;
Double_t        met_Et;
Double_t        met_E;
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
Double_t        btagSF;
Double_t        btagSF_up;
Double_t        btagSF_down;
Double_t        puSF;
Double_t        puSF_up;
Double_t        puSF_down;
Double_t        muonIdSF_BCDEF[1];   //[nMuons]
Double_t        muonIdSF_GH[1];   //[nMuons]
Double_t        muonIdSF_up_BCDEF[1];   //[nMuons]
Double_t        muonIdSF_up_GH[1];   //[nMuons]
Double_t        muonIdSF_down_BCDEF[1];   //[nMuons]
Double_t        muonIdSF_down_GH[1];   //[nMuons]
Double_t        muonIsoSF_BCDEF[1];   //[nMuons]
Double_t        muonIsoSF_GH[1];   //[nMuons]
Double_t        muonIsoSF_up_BCDEF[1];   //[nMuons]
Double_t        muonIsoSF_up_GH[1];   //[nMuons]
Double_t        muonIsoSF_down_BCDEF[1];   //[nMuons]
Double_t        muonIsoSF_down_GH[1];   //[nMuons]
Double_t        muonTrigSF_BCDEF[1];   //[nMuons]
Double_t        muonTrigSF_GH[1];   //[nMuons]
Double_t        muonTrigSF_up_BCDEF[1];   //[nMuons]
Double_t        muonTrigSF_up_GH[1];   //[nMuons]
Double_t        muonTrigSF_down_BCDEF[1];   //[nMuons]
Double_t        muonTrigSF_down_GH[1];   //[nMuons]
Double_t        muonTrackSF_eta[1];   //[nMuons]
Double_t        muonTrackSF_aeta[1];   //[nMuons]
Double_t        muonTrackSF_nPV[1];   //[nMuons]

Long64_t        nEvents;
Long64_t        nEventsSel;
Int_t           cutFlow[10];
Int_t           cutFlow2[10];
Int_t           appliedJER;
Int_t           appliedJES;
Int_t           appliedPU;
Long64_t        nofEventsRunB;
Long64_t        nofEventsRunCD;
Long64_t        nofEventsRunEF;
Long64_t        nofEventsRunG;
Long64_t        nofEventsRunH;
Long64_t        nofSelEventsRunB;
Long64_t        nofSelEventsRunCD;
Long64_t        nofSelEventsRunEF;
Long64_t        nofSelEventsRunG;
Long64_t        nofSelEventsRunH;
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
TBranch        *b_isDataRunB;   //!
TBranch        *b_isDataRunC;   //!
TBranch        *b_isDataRunD;   //!
TBranch        *b_isDataRunE;   //!
TBranch        *b_isDataRunF;   //!
TBranch        *b_isDataRunG;   //!
TBranch        *b_isDataRunH;   //!
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
TBranch        *b_met_pt;   //!
TBranch        *b_met_phi;   //!
TBranch        *b_met_eta;   //!
TBranch        *b_met_Et;   //!
TBranch        *b_met_E;   //!
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
TBranch        *b_btagSF;   //!
TBranch        *b_btagSF_up;   //!
TBranch        *b_btagSF_down;   //!
TBranch        *b_puSF;   //!
TBranch        *b_puSF_up;   //!
TBranch        *b_puSF_down;   //!
TBranch        *b_muonIdSF_BCDEF;   //!
TBranch        *b_muonIdSF_GH;   //!
TBranch        *b_muonIdSF_up_BCDEF;   //!
TBranch        *b_muonIdSF_up_GH;   //!
TBranch        *b_muonIdSF_down_BCDEF;   //!
TBranch        *b_muonIdSF_down_GH;   //!
TBranch        *b_muonIsoSF_BCDEF;   //!
TBranch        *b_muonIsoSF_GH;   //!
TBranch        *b_muonIsoSF_up_BCDEF;   //!
TBranch        *b_muonIsoSF_up_GH;   //!
TBranch        *b_muonIsoSF_down_BCDEF;   //!
TBranch        *b_muonIsoSF_down_GH;   //!
TBranch        *b_muonTrigSF_BCDEF;   //!
TBranch        *b_muonTrigSF_GH;   //!
TBranch        *b_muonTrigSF_up_BCDEF;   //!
TBranch        *b_muonTrigSF_up_GH;   //!
TBranch        *b_muonTrigSF_down_BCDEF;   //!
TBranch        *b_muonTrigSF_down_GH;   //!
TBranch        *b_muonTrackSF_eta;   //!
TBranch        *b_muonTrackSF_aeta;   //!
TBranch        *b_muonTrackSF_nPV;   //!

TBranch        *b_nEvents;   //!
TBranch        *b_nEventsSel;   //!
TBranch        *b_cutFlow;   //!
TBranch        *b_cutFlow2;   //!
TBranch        *b_appliedJER;   //!
TBranch        *b_appliedJES;   //!
TBranch        *b_appliedPU;   //!
TBranch        *b_nofEventsRunB;   //!
TBranch        *b_nofEventsRunCD;   //!
TBranch        *b_nofEventsRunEF;   //!
TBranch        *b_nofEventsRunG;   //!
TBranch        *b_nofEventsRunH;   //!
TBranch        *b_nofSelEventsRunB;   //!
TBranch        *b_nofSelEventsRunCD;   //!
TBranch        *b_nofSelEventsRunEF;   //!
TBranch        *b_nofSelEventsRunG;   //!
TBranch        *b_nofSelEventsRunH;   //!
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


double scaleFactor, widthSF;
double topPtRewSF, topPtSF, antiTopPtSF;
vector<unsigned int> bJetId;
double bdiscrTop, bdiscrTop2, tempbdiscr;
int labelB1, labelB2;
int labelsReco[4];
double reco_minMlb, reco_minMl_nonb, reco_ttbarMass, reco_dRLepB_lep, reco_dRLepB_had;
double min_Mlb, dRLepB;
int labelMlb, labelMl_nonb;
double massForWidth;

string catSuffix = "";
string catSuffixList[] = {"_CP", "_WP", "_UP"};
bool isCP, isWP, isUP;


/// Define TLVs
TLorentzVector muon, jet, mcpart;
TLorentzVector WCandidate;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
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
bool muonmatched = false;
bool muPlusFromTop = false, muMinusFromTop = false;
vector<unsigned int> partonId;

double matchedWMass_reco, matchedWMass_gen, matchedTopMass_reco, matchedTopMass_gen;
double matchedMlb_corr, matchedMlb_wrong, matchedTtbarMass_corr, matchedTtbarMass_wrong;
double matchedDRLepB_corr, matchedDRLepB_wrong;

/// KinFitter
TKinFitter* kFitter;
TKinFitter* kFitterMatched;
bool addWMassKF = true;
bool addEqMassKF = false;
int kFitVerbosity = 0;
double kFitChi2 = 99., kFitChi2Matched = 99.;
int nofAcceptedKFit = 0, nofAcceptedKFitMatched = 0;
double Wmass_reco_orig, Wmass_reco_kf, topmass_reco_orig, topmass_reco_kf, topmass_reco_kf_matched;
double toppt_reco_orig, toppt_reco_kf;

/// Likelihood
//Double_t widthArray[] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2., 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.5, 2.55, 2.60, 2.65, 2.70, 2.75, 3., 3.25, 3.5, 3.75, 4., 4.25, 4.5, 4.75, 5., 5.25, 5.5, 5.75, 6., 6.25, 6.5, 6.75, 7., 7.25, 7.5, 7.75, 8.};
Double_t widthArray[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.};
const int nWidthsLL = sizeof(widthArray)/sizeof(widthArray[0]);
string stringWidthArray[nWidthsLL];


//bool isGoodLL = false, isGoodLL_gen = false;
//int nofGoodEvtsLL[20] = {0};
//int nofGoodEvtsLL_gen[20] = {0};
//int nofBadEvtsLL[20] = {0};
Double_t aveTopMassLL = aveTopMass[2];
Double_t redTopMass = -1.;
Double_t maxRedTopMass = 0., minRedTopMass = 9999.;
Double_t minCutRedTopMass = 0.6, maxCutRedTopMass = 1.4;

ofstream txtDebugTopMass, txtDebugPUSF;

/// Pseudo experiments
const int nPseudoExps = 1000;
int nPsExps = nPseudoExps;
TRandom3 random3;
double toyValues[nPseudoExps]; // = random3.Uniform(0,1);
double toyMax;
int nEvtsInPseudoExp[nPseudoExps][15] = {0}, nDataEvts;

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
  cout << "*********************************************" << endl;
  cout << "* The following scale factors are applied:  *" << endl;
  cout << "*   - Lepton scale factors: ";
  if (applyLeptonSF)          cout << "nominal         *" << endl;
  else if (applyLeptonSFup)   cout << "scale down      *" << endl;
  else if (applyLeptonSFdown) cout << "scale up        *" << endl;
  cout << "*   - Pile up: ";
  if (applyPU)          cout << "nominal                      *" << endl;
  else if (applyPUup)   cout << "scale up                     *" << endl;
  else if (applyPUdown) cout << "scale down                   *" << endl;
  cout << "*   - B tag scale factors: ";
  if (applyBTagSF)     cout << "nominal          *" << endl;
  if (applyBTagSFup)   cout << "scale down       *" << endl;
  if (applyBTagSFdown) cout << "scale up         *" << endl;
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
    calculateLikelihood = false;
    doPseudoExps = false;
    makePlots = false;
    doGenOnly = false;
  }
  if (calculateResolutionFunctions)
  {
    testTTbarOnly = true;
    doGenOnly = true;
    makePlots = false;
    calculateAverageMass = false;
    calculateLikelihood = false;
    doPseudoExps = false;
    doKinFit = false;
  }
  if (test) makePlots = false;
  if (testHistos) makePlots = true;
  if (doGenOnly)
  {
    doPseudoExps = false;
  }
  if (makeTGraphs) calculateLikelihood = false;
  if (calculateLikelihood) makeTGraphs = false;
  else doPseudoExps = false;
  if (doPseudoExps) makePlots = false;
  //string pathOutput = "test/";
  string pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (makePlots)
  {
    // Add channel to output path
    pathOutput += channel+"/";
    mkdir(pathOutput.c_str(),0777);
    if (doPseudoExps)
    {
      pathOutput+= "PseudoExp/";
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
  
  if (makeTGraphs)
  {
    mkdir(outputDirLL.c_str(),0777);
    outputDirLL += dateString+"/";
    mkdir(outputDirLL.c_str(),0777);
  }
  else if (calculateLikelihood)
  {
    inputDirLL = outputDirLL+inputDateLL;
  }
  
  pathNtuplesMC = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.first+"/";
  pathNtuplesData = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.second+"/";
  cout << "Using Ntuples from " << ntupleDate.first << " for MC and " << ntupleDate.second << " for data. This corresponds to systematics: " << systStr << endl;
  if (calculateAverageMass) cout << "Calculating average mass values..." << endl;
  if (calculateLikelihood) cout << "Calculating -loglikelihood values..." << endl;
  if (doPseudoExps)       cout << "              for pseudo experiments" << endl;
  if (testHistos) cout << "Testing histogram consistency..." << endl;
  if (doGenOnly) cout << "Running only matching..." << endl;
  if (applyWidthSF) cout << "TTbar sample width will be scaled by a factor " << scaleWidth << endl;
  else scaleWidth = 1.;
  
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
  treeLoader.LoadDatasets(datasetsMSP, xmlfile);
  
  int dTT = -1;
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    string dataSetName = datasets[d]->Name();
//    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA") != std::string::npos )
//      Luminosity = datasets[d]->EquivalentLumi();
    
    if ( dataSetName.find("QCD") != std::string::npos )
    {
      datasets[d]->SetColor(kYellow);
      datasetsMSP[d]->SetColor(kYellow);
    }
    if ( dataSetName.find("TT") != std::string::npos )
    {
      datasets[d]->SetTitle("t#bar{t}");
      datasets[d]->SetColor(kRed+1);
      datasetsMSP[dTT+2]->SetName("TT_UP");
      datasetsMSP[dTT+2]->SetTitle("t#bar{t} UP");
      datasetsMSP[dTT+2]->SetColor(kRed-10);
      
      dTT = d;
    }
    //if ( dataSetName.find("TTbarJets_Other") != std::string::npos ) datasets[d]->SetColor(kRed-7);
    if ( dataSetName.find("WJets") != std::string::npos )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      datasets[d]->SetColor(kGreen-3);
      datasetsMSP[d]->SetTitle("W#rightarrowl#nu");
      datasetsMSP[d]->SetColor(kGreen-3);
    }
    if ( dataSetName.find("ZJets") != std::string::npos || dataSetName.find("DY") != std::string::npos )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      //datasets[d]->SetColor(kAzure-2);
      datasets[d]->SetColor(kMagenta);
      datasetsMSP[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      //datasetsMSP[d]->SetColor(kAzure-2);
      datasetsMSP[d]->SetColor(kMagenta);
    }
    if ( dataSetName.find("ST") != std::string::npos || dataSetName.find("SingleTop") != std::string::npos )
    {
      datasets[d]->SetTitle("ST");
      datasets[d]->SetColor(kBlue-2);
      datasetsMSP[d]->SetTitle("ST");
      datasetsMSP[d]->SetColor(kBlue-2);
      //if ( dataSetName.find("tW") != std::string::npos )
      //{
      //  datasets[d]->SetTitle("ST tW");
      //  datasets[d]->SetColor(kBlue-4);
      //}
      //else
      //  datasets[d]->SetTitle("ST t");
    }
  }
  
  treeLoader.LoadDatasets(datasetsTemp, xmlfile);
  datasetsTemp[dTT]->SetName("TT_WP");
  datasetsTemp[dTT]->SetTitle("t#bar{t} WP");
  datasetsTemp[dTT]->SetColor(kRed-9);
  datasetsMSP.insert(datasetsMSP.begin()+dTT, 1, datasetsTemp[dTT]);
  datasetsTemp.clear();
  treeLoader.LoadDatasets(datasetsTemp, xmlfile);
  datasetsTemp[dTT]->SetName("TT_CP");
  datasetsTemp[dTT]->SetTitle("t#bar{t} CP");
  datasetsTemp[dTT]->SetColor(kRed-7);
  datasetsMSP.insert(datasetsMSP.begin()+dTT, 1, datasetsTemp[dTT]);
  
//   cout << "Content datasets  " << endl;
//   for (int d = 0; d < datasets.size(); d++)
//     cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
//   cout << "Content datasetsMSP  " << endl;
//   for (int d = 0; d < datasetsMSP.size(); d++)
//     cout << "   Dataset " << d << ": " << datasetsMSP[d]->Name() << " / title : " << datasetsMSP[d]->Title() << endl;
  
  
  
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
    kf = new KinFitter("PlotsForResolutionFunctions_testFit.root", addWMassKF, addEqMassKF);
    kfMatched = new KinFitter("PlotsForResolutionFunctions_testFit.root", addWMassKF, addEqMassKF);
  }
  
  if (makeTGraphs)
  {
    like = new Likelihood(minCutRedTopMass, maxCutRedTopMass, outputDirLL, dateString, makeTGraphs, false, false);  // calculateGoodEvtLL, verbose
  }
  if (calculateLikelihood)
  {
    like = new Likelihood(minCutRedTopMass, maxCutRedTopMass, inputDirLL, dateString, makeTGraphs, true, false);  // calculateGoodEvtLL, verbose
    calculateLikelihood = like->ConstructTGraphsFromFile();
    calculateLikelihood = like->ConstructTGraphsFromFile("CorrectMatchLikelihood_");
  }
  if (doPseudoExps)
  {
    nPsExps = like->InitPull(nPseudoExps);
  }
  
  if (makePlots)
  {
    if (! doGenOnly) InitMSPlots();
    InitHisto1D();
    InitHisto2D();
  }
  
  vJER.clear(); vJES.clear(); vPU.clear();
  
  
  if (calculateAverageMass)
  {
    mkdir("averageMass/",0777);
    txtMassGenPMatched.open(("averageMass/mass_genp_matched_TT_"+dateString+".txt").c_str());
    txtMassGenJMatched.open(("averageMass/mass_genj_matched_TT_"+dateString+".txt").c_str());
    txtMassRecoCP.open(("averageMass/mass_reco_matched_TT_"+dateString+".txt").c_str());
    txtMassRecoWPUP.open(("averageMass/mass_reco_notCorrectMatch_TT_"+dateString+".txt").c_str());
    txtMassRecoUP.open(("averageMass/mass_reco_notMatched_TT_"+dateString+".txt").c_str());
    txtMassRecoWP.open(("averageMass/mass_reco_wrongPerm_TT_"+dateString+".txt").c_str());
  }
  
  txtDebugTopMass.open("debug_missing_topQ.txt");
  txtDebugTopMass << "## Events where top quark(s) not found in genParticle collection" << endl;
  txtDebugTopMass << "#  If lepton charge > 0 : leptonically decaying top, hadronically decaying antitop" << endl;
  txtDebugTopMass << "#  If lepton charge < 0 : hadronically decaying top, leptonically decaying antitop" << endl;
  
  txtDebugPUSF.open("debug_pu_sf.txt");
  txtDebugPUSF << "## PU SF = 0" << endl;
  txtDebugPUSF << "#  nVtx    nTruePU" << endl;
  
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName, slumi;
  double timePerDataSet[datasets.size()];
  
  int nEntries;
  double fracDataEras[2] = {-1.};
  //GetEraFraction(fracDataEras);
  fracDataEras[0] = lumi_runBCDEF/(lumi_runBCDEF + lumi_runGH);
  fracDataEras[1] = lumi_runGH/(lumi_runBCDEF + lumi_runGH);
  if ( fracDataEras[0] == -1. || fracDataEras[1] == -1. )
  {
    cerr << "Something went wrong with the fraction calculation for muon SFs!" << endl;
    exit(1);
  }
  cout << "The muon scale factors will be scaled by " << fracDataEras[0] << " for eras B-F and " << fracDataEras[1] << " for eras G-H." << endl;
  
  
  bool hasFoundTTbar = false;
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
      pathNtuples = pathNtuplesData;
    }
    else if ( dataSetName.find("TT") != std::string::npos )
    {
      isTTbar = true;
      hasFoundTTbar = true;
    }
    
    if (! isData)
    {
      pathNtuples = pathNtuplesMC;
    }
    
    if (testTTbarOnly && ! isTTbar)
    {
      cout << "Skipping dataset..." << endl;
      continue;
    }
    
    if (calculateAverageMass)
    {
      txtMassReco.open(("averageMass/mass_reco_"+dataSetName+"_"+dateString+".txt").c_str());
    }
    
    
    string ntupleFileName = "Ntuples_"+dataSetName+".root";
    
    
    /// Change name of ttbar dataset to TT (whether it is nominal or widthxX)
    if (isTTbar) dataSetName = "TT";
    
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
    
    /// eqLumi calculation
    if (isData) eqLumi = Luminosity;
    else eqLumi = (double)GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData)/datasets[d]->Xsection();  // 1/pb
    
    if (doPseudoExps)
    {
      if (! isData) toyMax = Luminosity/eqLumi;  // CHECK: better nDataEvts/nEvtsPassKinFit ??
      //else nDataEvts = GetNEvents(tStatsTree[(dataSetName).c_str()], "nEventsSel", 1);
      else nDataEvts = nEventsAKF[d];
      cout << "PseudoExperiments::Number of selected data events: " << nDataEvts << endl;
      
      if (test)
        cout << "      Lumi : " << Luminosity << "/pb; eqLumi: " << eqLumi << "/pb." << endl;
      cout << "PseudoExperiments::Lumi/eqLumi = " << toyMax;
      if (! isData)
      {
        toyMax *= 0.93;  // small overshoot in MC --> scale down
        cout << " x 0.92 = " << toyMax;
      }
      cout << endl;
    }
    
    
    /// Get data
    tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
    nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
    cout << "                nEntries: " << nEntries << endl;
    if (isData) cout << "                Lumi    : " << Luminosity << "/pb" << endl;
    else cout << "                eqLumi  : " << eqLumi << "/pb = " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData) << " / " << datasets[d]->Xsection() << " pb" << endl;
    
    if (isData)
    {
      cout << "NEvents data:   " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData) << endl;
      cout << "NEvents Run B:  " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunB", isData) << endl;
      cout << "NEvents Run CD: " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunCD", isData) << endl;
      cout << "NEvents Run EF: " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunEF", isData) << endl;
      cout << "NEvents Run G:  " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunG", isData) << endl;
      cout << "NEvents Run H:  " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunH", isData) << endl;
    }
    
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
      
      
      /// Scale factors
      if (! isData)
      {
        if (applyLeptonSF) { scaleFactor *= muonTrackSF_eta[0] * (fracDataEras[0]*muonIdSF_BCDEF[0] + fracDataEras[1]*muonIdSF_GH[0]) * (fracDataEras[0]*muonIsoSF_BCDEF[0] + fracDataEras[1]*muonIsoSF_GH[0]) * (fracDataEras[0]*muonTrigSF_BCDEF[0] + fracDataEras[1]*muonTrigSF_GH[0]);}
        else if (applyLeptonSFup) { scaleFactor *= muonTrackSF_eta[0]*1.01 * (fracDataEras[0]*muonIdSF_up_BCDEF[0] + fracDataEras[1]*muonIdSF_up_GH[0]) * (fracDataEras[0]*muonIsoSF_up_BCDEF[0] + fracDataEras[1]*muonIsoSF_up_GH[0]) * (fracDataEras[0]*muonTrigSF_up_BCDEF[0] + fracDataEras[1]*muonTrigSF_up_GH[0]);}  // CHECK syst muon track SF
        else if (applyLeptonSFdown) { scaleFactor *= muonTrackSF_eta[0]*0.99 * (fracDataEras[0]*muonIdSF_down_BCDEF[0] + fracDataEras[1]*muonIdSF_down_GH[0]) * (fracDataEras[0]*muonIsoSF_down_BCDEF[0] + fracDataEras[1]*muonIsoSF_down_GH[0]) * (fracDataEras[0]*muonTrigSF_down_BCDEF[0] + fracDataEras[1]*muonTrigSF_down_GH[0]);}  // CHECK syst muon track SF
        if (applyPU) { scaleFactor *= puSF;}
        else if (applyPUup) { scaleFactor *= puSF_up;}
        else if (applyPUdown) { scaleFactor *= puSF_down;}
        if (applyBTagSF) { scaleFactor *= btagSF;}
        else if (applyBTagSFup) { scaleFactor *= btagSF_up;}
        else if (applyBTagSFdown) { scaleFactor *= btagSF_down;}
//        cout << "Scalefactor: " << setw(6) << scaleFactor << "  btag SF: " << setw(6) << btagSF << "  pu SF: " << setw(6) << puSF << "  muonId: " << setw(6) << fracDataEras[0]*muonIdSF_BCDEF[0] + fracDataEras[1]*muonIdSF_GH[0] << "  muonIso: " << setw(6) << fracDataEras[0]*muonIsoSF_BCDEF[0] + fracDataEras[1]*muonIsoSF_GH[0] << "  muonTrig: " << setw(6) << fracDataEras[0]*muonTrigSF_BCDEF[0] + fracDataEras[1]*muonTrigSF_GH[0] << endl;
        if (applyPU && puSF == 0) txtDebugPUSF << nvtx << "    " << npu << endl;
      }
      
      //if (useToys) scaleFactor *= eqLumi;  // undo eqLumi scaling in MSPlots
      
      
      
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
      
      if ( selectedJets.size() > 4 ) continue;
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
      
      //if ( isTTbar || dataSetName.find("ST") != std::string::npos )
      if (! isData)
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
          
          if ( mc_pdgId[i] == pdgID_top && mc_status[i] == 62 )  // top
          {
            if ( mcParticles[i].Pt() < 400 ) topPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
            else topPtSF = TMath::Exp(0.0615-0.0005*400.);
          }
          else if ( mc_pdgId[i] == -pdgID_top && mc_status[i] == 62 )  // antitop
          {
            if ( mcParticles[i].Pt() < 400 ) antiTopPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
            else antiTopPtSF = TMath::Exp(0.0615-0.0005*400.);
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
          
          if (makePlots) histo1D["Width_SF"]->Fill(widthSF);
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
            
            
            matchedWMass_reco = (jetsMatched[0] + jetsMatched[1]).M();
            matchedWMass_gen = (partonsMatched[0] + partonsMatched[1]).M();
            matchedTopMass_reco = (jetsMatched[0] + jetsMatched[1] + jetsMatched[2]).M();
            matchedTopMass_gen = (partonsMatched[0] + partonsMatched[1] + partonsMatched[2]).M();
            
            if (calculateAverageMass)
            {
              txtMassGenPMatched << ievt << "  " << matchedTopMass_gen << endl;
              txtMassGenJMatched << ievt << "  " << matchedTopMass_reco << endl;
            }
            
            
            /// KF for matched jets
            if (doKinFit && ! doPseudoExps)
            {
              kFitterMatched = kfMatched->doFit(jetsMatched[0], jetsMatched[1], kFitVerbosity);
              
              if ( kFitterMatched->getStatus() != 0 )  // did not converge
              {
                if (test && verbose > 2) cout << "Event " << ievt << ": Fit for matched events did not converge..." << endl;
                continue;
              }
              
              kFitChi2Matched = kFitterMatched->getS();
              if (test && verbose > 4) cout << "Fit converged: Chi2 = " << kFitChi2Matched << endl;
              
              if ( applyKinFitCut && kFitChi2Matched > kinFitCutValue ) continue;
              nofAcceptedKFitMatched++;
              
              selectedJetsKFMatched.clear();
              selectedJetsKFMatched = kfMatched->getCorrectedJets();
              
              if ( selectedJetsKFMatched.size() == 2 ) selectedJetsKFMatched.push_back(jetsMatched[2]);
              
              topmass_reco_kf_matched = (selectedJetsKFMatched[0] + selectedJetsKFMatched[1] + selectedJetsKFMatched[2]).M();
              
              if (calculateLikelihood)
              {
                double temp = topmass_reco_kf_matched/aveTopMassLL;
                like->CalculateGenLikelihood(temp, massForWidth, scaleWidth, isTTbar, isData);
              }
            }  // end KF
            
            
            if (isTTbar && makePlots)
            {
              FillMatchingPlots();
            }

          }  // end hadronicTopJetsMatched
        
        }  // end doMatching
        
        
      }  // end if not data
      
      
      
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
            isCP = true;
            nofCorrectlyMatched++;
          }
          else  // wrong permutation
          {
            isWP = true;
            nofNotCorrectlyMatched++;
          }
        }  // end hadrTopMatch
        else  // no match
        {
          isUP = true;
        }
        
        
        if ( (! isCP && ! isWP && ! isUP) || (isCP && isWP) || (isCP && isUP) || (isWP && isUP) )
          cerr << "Something wrong with trigger logic CP/WP/UP !! " << endl;
        
      }  // not Data
      
      
      if (isCP) catSuffix = catSuffixList[0];
      else if (isWP) catSuffix = catSuffixList[1];
      else if (isUP) catSuffix = catSuffixList[2];
      
      
      
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
        
        if ( applyKinFitCut && kFitChi2 > kinFitCutValue ) continue;
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
      else if ( selectedJetsKFcorrected.size() == 2 ) selectedJetsKFcorrected.push_back(selectedJets[labelsReco[2]]);
      
      Wmass_reco_orig = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).M();
      Wmass_reco_kf = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).M();
      topmass_reco_orig = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).M();
      topmass_reco_kf = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).M();
      toppt_reco_orig = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).Pt();
      toppt_reco_kf = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).Pt();
      
      if ( topmass_reco_kf < 0. )
        PrintKFDebug(ievt);
      
      if (calculateAverageMass) txtMassReco << ievt << "  " << topmass_reco_kf << endl;
      
      if ( doKinFit && makePlots )
        FillKinFitPlots();
      
      
      //  Leptonic variables
      reco_minMlb = (selectedLepton[0] + selectedJets[labelsReco[3]]).M();
      reco_dRLepB_lep = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[3]], selectedLepton[0] );  // deltaR between lepton and leptonic b jet
      reco_dRLepB_had = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
      
      //  Combination
      reco_ttbarMass = reco_minMlb + topmass_reco_kf;
      
      
      /// Reduced top mass
      redTopMass = topmass_reco_kf/aveTopMassLL;
      if (makeTGraphs) like->FillHistograms(redTopMass, massForWidth, isTTbar, isData, catSuffix);
      if (calculateLikelihood)
      {
        like->CalculateLikelihood(redTopMass, massForWidth, scaleWidth, isTTbar, isData);
        if (isCP && ! doPseudoExps)  // isCP ensures ! isData
          like->CalculateCPLikelihood(redTopMass, massForWidth, scaleWidth, isTTbar, isData);
      }
      
      if ( redTopMass > maxRedTopMass ) maxRedTopMass = redTopMass;
      if ( redTopMass < minRedTopMass ) minRedTopMass = redTopMass;
      if ( redTopMass > minCutRedTopMass && redTopMass < maxCutRedTopMass )
      {
        if (isCP)
        {
          nofCP++;
          nofCP_weighted += widthSF;
          if (isTTbar) nofCP_TT++;
        }
        else if (isWP)
        {
          nofWP++;
          nofWP_weighted += widthSF;
          if (isTTbar) nofWP_TT++;
        }
        else if (isUP)
        {
          nofUP++;
          nofUP_weighted += widthSF;
          if (isTTbar) nofUP_TT++;
        }
      }
      
      if (calculateAverageMass && ! isData)
      {
        if (isCP) txtMassRecoCP << ievt << "  " << topmass_reco_kf << endl;
        else
        {
          txtMassRecoWPUP << ievt << "  " << topmass_reco_kf << endl;
          if (isWP)
            txtMassRecoWP << ievt << "  " << topmass_reco_kf << endl;
          else if (isUP)
            txtMassRecoUP << ievt << "  " << topmass_reco_kf << endl;
        }
      }  // end aveMassCalc
      
      
      /// Pseudo experiments
      if (doPseudoExps && ! isData)
      {
        random3.RndmArray(nPsExps, toyValues);
        for (int iPsExp = 0; iPsExp < nPsExps; iPsExp++)
        {
          //toyValue = random3.Rndm();
          //if ( toyValue > toyMax ) continue;
          if ( toyValues[iPsExp] > toyMax ) continue;
          (nEvtsInPseudoExp[iPsExp][d])++;
          like->AddPsExp(iPsExp, isData);
        }
      }
      
      //Fill histos
      if (makePlots)
      {
        histo1D[("Reduced_top_mass_"+dataSetName).c_str()]->Fill(redTopMass);
        if (isTTbar)
        {
          histo1D["Reco_W_mass_reco"]->Fill(Wmass_reco_kf, widthSF);
          histo1D["Reco_top_mass_reco"]->Fill(topmass_reco_kf, widthSF);
          histo1D["dR_lep_b_reco"]->Fill(reco_dRLepB_lep, widthSF);
        }
        
        FillCatsPlots(catSuffix);
        
        int dMSP = d;
        if (hasFoundTTbar && ! isTTbar) dMSP = d+2;
        else if (isTTbar && isWP) dMSP = d+1;
        else if (isTTbar && isUP) dMSP = d+2;
        
        FillMSPlots(dMSP);
        
      }  // end makePlots
      
      
      
      
    }  // end loop events
    
    
    cout << endl;  /// Stronger selection in this analyser compared to Ntuples ==> endEvent --> nofHardSelected
    cout << "Number of events with exactly 4 jets with pT > 30 GeV: " << nofHardSelected << " (" << 100*((float)nofHardSelected/(float)endEvent) << "%)" << endl;
    cout << "Number of events with clean MET: " << nofMETCleaned << " (" << 100*((float)nofMETCleaned/(float)nofHardSelected) << "%)" << endl;
    if (doKinFit) cout << "Number of clean events accepted by kinFitter: " << nofAcceptedKFit << " (" << 100*((float)nofAcceptedKFit/(float)nofMETCleaned) << "%)" << endl;
    
    //if ( isTTbar || dataSetName.find("ST") != std::string::npos )
    if (! isData && nofMatchedEvents > 0 )
    {
      cout << "Number of matched events: " << setw(8) << right << nofMatchedEvents << endl;
      cout << "Number of events with hadronic top matched: " << setw(8) << right << nofHadrMatchedEvents << endl;
      if (! doGenOnly)
      {
        cout << "Correctly matched reconstructed events:     " << setw(8) << right << nofCorrectlyMatched << endl;
        cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatched << endl;
        if ( nofCorrectlyMatched != 0 || nofNotCorrectlyMatched != 0 )
          cout << "   ===> This means that " << 100*(float)nofCorrectlyMatched / (float)(nofCorrectlyMatched + nofNotCorrectlyMatched) << "% is correctly matched." << endl;
      }
      if (doKinFit) cout << "Number of events accepted by kinFitter: " << nofAcceptedKFitMatched << " (" << 100*((float)nofAcceptedKFitMatched/(float)nofHardSelected) << "%)" << endl;
      
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
    txtMassRecoCP.close();
    txtMassRecoWPUP.close();
    txtMassRecoUP.close();
    txtMassRecoWP.close();
  }
  
  if (applyWidthSF) txtDebugTopMass.close();
  
  if (! doGenOnly && ! testTTbarOnly)
  {
    cout << "Number of events with " << minCutRedTopMass << " < mt/<mt> < " << maxCutRedTopMass << " : CP: " << nofCP << " (" << 100*(double)nofCP/((double)(nofCP+nofWP+nofUP)) << "%)   WP: " << nofWP << " (" << 100*(double)nofWP/((double)(nofCP+nofWP+nofUP)) << "%)   UP: " << nofUP << " (" << 100*(double)nofUP/((double)(nofCP+nofWP+nofUP)) << "%)   Total: " << nofCP+nofWP+nofUP << endl;
    cout << "                                  weighted: CP: " << nofCP_weighted << " (" << 100*nofCP_weighted/(nofCP_weighted+nofWP_weighted+nofUP_weighted) << "%)   WP: " << nofWP_weighted << " (" << 100*nofWP_weighted/(nofCP_weighted+nofWP_weighted+nofUP_weighted) << "%)   UP: " << nofUP_weighted << " (" << 100*nofUP_weighted/(nofCP_weighted+nofWP_weighted+nofUP_weighted) << "%)   Total: " << (int)(nofCP_weighted+nofWP_weighted+nofUP_weighted) << endl;
    cout << "                               (TTbar only) CP: " << nofCP_TT << "              WP: " << nofWP_TT << "              UP: " << nofUP_TT << endl;
  }
  
  
  if (makeTGraphs)
  {
    like->WriteHistograms("ReducedTopMassPlots.root");
    like->ConstructTGraphsFromHisto("TGraphFunctions.root");
  }
  
  if (calculateLikelihood)
  {
    cout << "Minimum reduced top mass: " << minRedTopMass << endl;
    cout << "Maximum reduced top mass: " << maxRedTopMass << endl;
    
    /// Print output to file
    string llFileName = "output_loglikelihood_widthx1";
    if (applyWidthSF) llFileName = "output_loglikelihood_widthx"+DotReplace(scaleWidth);
    if (doGenOnly) llFileName = "output_loglikelihood_parton_widthx"+DotReplace(scaleWidth);
    //if (useToys) llFileName = "output_loglikelihood_toys";
    like->PrintLikelihoodOutput(llFileName+".txt");
    //if (unblind) like->PrintLikelihoodOutputData(llFileName+"_data.txt");
    like->PrintMtmLikelihoodOutput(llFileName+"_Mtm.txt");
    
    /// Calculate output width
    cout << "Standard output width: " << endl;
    like->GetOutputWidth(scaleWidth);
    if (! doPseudoExps)
    {
      cout << "Output width for correctly matched events (using likelihood with only CP template): " << endl;
      like->GetOutputWidth(scaleWidth, "CP");
      cout << "Output width for generated evetns (using likelihood with only CP template): " << endl;
      like->GetOutputWidth(scaleWidth, "gen");
    }
    cout << "Output width from file (standard calculation): " << endl;
    like->GetOutputWidth(llFileName+".txt", scaleWidth);
  }
  
  if (doPseudoExps)
  {
    like->CalculatePull(scaleWidth);
  }
  
  
  cout << "Processing time per dataset: " << endl;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
  }
  
  if (doPseudoExps)
  {
    cout << "PseudoExperiments::Number of selected data events: " << nDataEvts << endl;
    double totMCPsExp;
    for (int i = 0; i < 3; i++)
    {
      totMCPsExp = 0;
      cout << "PseudoExperiment " << std::setw(3) << std::right << i << endl;
      for (unsigned int d = 1; d < datasets.size(); d++)
      {
        totMCPsExp += nEvtsInPseudoExp[i][d];
        cout << "                      " << datasets[d]->Name() << ": " << nEvtsInPseudoExp[i][d] << endl;
      }
      cout << "                                            " << "Total MC: " << totMCPsExp << endl;
    }
  }
  
  
  ///  Check Shape Changing Systematics
  if (! testTTbarOnly) CheckSystematics(vJER, vJES, vPU);
  
  
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
  if (! doGenOnly) mkdir((pathOutput+"MSPlot/").c_str(),0777);
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile ((pathOutput+rootFileName).c_str(), "RECREATE");
  cout << "   Output file is " << pathOutput+rootFileName << endl;
  
  ///Write histograms
  fout->cd();
  
  if (! doGenOnly)
  {
    for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
      //cout << "MSPlot: " << it->first << endl;
      MultiSamplePlot *temp = it->second;
      string name = it->first;
      temp->Draw(name, 1, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
      temp->Write(fout, name, true, pathOutput+"MSPlot", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
    }
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
  tree->SetBranchAddress("cutFlow", cutFlow, &b_cutFlow);
  tree->SetBranchAddress("cutFlow2", cutFlow2, &b_cutFlow2);
  tree->SetBranchAddress("appliedJER", &appliedJER, &b_appliedJER);
  tree->SetBranchAddress("appliedJES", &appliedJES, &b_appliedJES);
  tree->SetBranchAddress("appliedPU", &appliedPU, &b_appliedPU);
  if (isData)
  {
   tree->SetBranchAddress("nofEventsRunB", &nofEventsRunB, &b_nofEventsRunB);
   tree->SetBranchAddress("nofEventsRunCD", &nofEventsRunCD, &b_nofEventsRunCD);
   tree->SetBranchAddress("nofEventsRunEF", &nofEventsRunEF, &b_nofEventsRunEF);
   tree->SetBranchAddress("nofEventsRunG", &nofEventsRunG, &b_nofEventsRunG);
   tree->SetBranchAddress("nofEventsRunH", &nofEventsRunH, &b_nofEventsRunH);
   tree->SetBranchAddress("nofSelEventsRunB", &nofSelEventsRunB, &b_nofEventsRunB);
   tree->SetBranchAddress("nofSelEventsRunCD", &nofSelEventsRunCD, &b_nofSelEventsRunCD);
   tree->SetBranchAddress("nofSelEventsRunEF", &nofSelEventsRunEF, &b_nofSelEventsRunEF);
   tree->SetBranchAddress("nofSelEventsRunG", &nofSelEventsRunG, &b_nofSelEventsRunG);
   tree->SetBranchAddress("nofSelEventsRunH", &nofSelEventsRunH, &b_nofSelEventsRunH);
  }
  else if (isTTbar)
  {
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
  tree->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
  tree->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
  tree->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
  tree->SetBranchAddress("met_Et", &met_Et, &b_met_Et);
  tree->SetBranchAddress("met_E", &met_E, &b_met_E);
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
    tree->SetBranchAddress("btagSF", &btagSF, &b_btagSF);
    tree->SetBranchAddress("btagSF_up", &btagSF_up, &b_btagSF_up);
    tree->SetBranchAddress("btagSF_down", &btagSF_down, &b_btagSF_down);
    tree->SetBranchAddress("puSF", &puSF, &b_puSF);
    tree->SetBranchAddress("puSF_up", &puSF_up, &b_puSF_up);
    tree->SetBranchAddress("puSF_down", &puSF_down, &b_puSF_down);
    tree->SetBranchAddress("muonIdSF_BCDEF", muonIdSF_BCDEF, &b_muonIdSF_BCDEF);
    tree->SetBranchAddress("muonIdSF_GH", muonIdSF_GH, &b_muonIdSF_GH);
    tree->SetBranchAddress("muonIdSF_up_BCDEF", muonIdSF_up_BCDEF, &b_muonIdSF_up_BCDEF);
    tree->SetBranchAddress("muonIdSF_up_GH", muonIdSF_up_GH, &b_muonIdSF_up_GH);
    tree->SetBranchAddress("muonIdSF_down_BCDEF", muonIdSF_down_BCDEF, &b_muonIdSF_down_BCDEF);
    tree->SetBranchAddress("muonIdSF_down_GH", muonIdSF_down_GH, &b_muonIdSF_down_GH);
    tree->SetBranchAddress("muonIsoSF_BCDEF", muonIsoSF_BCDEF, &b_muonIsoSF_BCDEF);
    tree->SetBranchAddress("muonIsoSF_GH", muonIsoSF_GH, &b_muonIsoSF_GH);
    tree->SetBranchAddress("muonIsoSF_up_BCDEF", muonIsoSF_up_BCDEF, &b_muonIsoSF_up_BCDEF);
    tree->SetBranchAddress("muonIsoSF_up_GH", muonIsoSF_up_GH, &b_muonIsoSF_up_GH);
    tree->SetBranchAddress("muonIsoSF_down_BCDEF", muonIsoSF_down_BCDEF, &b_muonIsoSF_down_BCDEF);
    tree->SetBranchAddress("muonIsoSF_down_GH", muonIsoSF_down_GH, &b_muonIsoSF_down_GH);
    tree->SetBranchAddress("muonTrigSF_BCDEF", muonTrigSF_BCDEF, &b_muonTrigSF_BCDEF);
    tree->SetBranchAddress("muonTrigSF_GH", muonTrigSF_GH, &b_muonTrigSF_GH);
    tree->SetBranchAddress("muonTrigSF_up_BCDEF", muonTrigSF_up_BCDEF, &b_muonTrigSF_up_BCDEF);
    tree->SetBranchAddress("muonTrigSF_up_GH", muonTrigSF_up_GH, &b_muonTrigSF_up_GH);
    tree->SetBranchAddress("muonTrigSF_down_BCDEF", muonTrigSF_down_BCDEF, &b_muonTrigSF_down_BCDEF);
    tree->SetBranchAddress("muonTrigSF_down_GH", muonTrigSF_down_GH, &b_muonTrigSF_down_GH);
    tree->SetBranchAddress("muonTrackSF_eta", muonTrackSF_eta, &b_muonTrackSF_eta);
    tree->SetBranchAddress("muonTrackSF_aeta", muonTrackSF_aeta, &b_muonTrackSF_aeta);
    tree->SetBranchAddress("muonTrackSF_nPV", muonTrackSF_nPV, &b_muonTrackSF_nPV);
  }
  else
  {
    tree->SetBranchAddress("isDataRunB", &isDataRunB, &b_isDataRunB);
    tree->SetBranchAddress("isDataRunC", &isDataRunC, &b_isDataRunC);
    tree->SetBranchAddress("isDataRunD", &isDataRunD, &b_isDataRunD);
    tree->SetBranchAddress("isDataRunE", &isDataRunE, &b_isDataRunE);
    tree->SetBranchAddress("isDataRunF", &isDataRunF, &b_isDataRunF);
    tree->SetBranchAddress("isDataRunG", &isDataRunG, &b_isDataRunG);
    tree->SetBranchAddress("isDataRunH", &isDataRunH, &b_isDataRunH);
  }
}

void InitMSPlots()
{
  MSPlot["muon_pT"] = new MultiSamplePlot(datasetsMSP, "muon_pT", 22, 0, 440, "p_{T} [GeV]");
  MSPlot["muon_eta"] = new MultiSamplePlot(datasetsMSP, "muon_eta", 30, -3, 3, "Eta");
  MSPlot["muon_phi"] = new MultiSamplePlot(datasetsMSP, "muon_phi", 32, -3.2, 3.2, "Phi");
  MSPlot["muon_relIso"] = new MultiSamplePlot(datasetsMSP, "muon_relIso", 20, 0, 0.2, "relIso");
  MSPlot["muon_d0"] = new MultiSamplePlot(datasetsMSP, "muon_d0", 50, 0, 0.03, "d_{0}");
  MSPlot["leadingJet_pT"] = new MultiSamplePlot(datasetsMSP, "leadingJet_pT", 60, 0, 600, "p_{T} [GeV]");
  MSPlot["jet2_pT"] = new MultiSamplePlot(datasetsMSP, "jet2_pT", 40, 0, 400, "p_{T} [GeV]");
  MSPlot["jet3_pT"] = new MultiSamplePlot(datasetsMSP, "jet3_pT", 40, 0, 400, "p_{T} [GeV]");
  MSPlot["jet4_pT"] = new MultiSamplePlot(datasetsMSP, "jet4_pT", 40, 0, 400, "p_{T} [GeV]");
  MSPlot["Ht_4leadingJets"] = new MultiSamplePlot(datasetsMSP,"Ht_4leadingJets", 60, 0, 1200, "H_{T} [GeV]");
  MSPlot["met_pT"] = new MultiSamplePlot(datasetsMSP, "met_pT", 40, 0, 400, "E_{T}^{miss} p_{T} [GeV]");
  MSPlot["met_eta"] = new MultiSamplePlot(datasetsMSP, "met_eta", 30, -3, 3, "E_{T}^{miss} Eta");
  MSPlot["met_phi"] = new MultiSamplePlot(datasetsMSP, "met_phi", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
  MSPlot["met_corr_pT"] = new MultiSamplePlot(datasetsMSP, "met_corr_pT", 40, 0, 400, "E_{T}^{miss} p_{T} [GeV]");
  MSPlot["met_corr_eta"] = new MultiSamplePlot(datasetsMSP, "met_corr_eta", 30, -3, 3, "E_{T}^{miss} Eta");
  MSPlot["met_corr_phi"] = new MultiSamplePlot(datasetsMSP, "met_corr_phi", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
  
  MSPlot["M3"] = new MultiSamplePlot(datasetsMSP, "M3", 40, 60, 460, "M_{3} [GeV]");

  MSPlot["nJets"] = new MultiSamplePlot(datasetsMSP, "nJets", 13, -0.5, 12.5, "# jets");
  MSPlot["nBJets"] = new MultiSamplePlot(datasetsMSP, "nBJets", 9, -0.5, 8.5, "# b jets");
  MSPlot["CSVv2Discr_allJets"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_allJets", 48, 0.0, 1.2, "CSVv2 discriminant value");
  MSPlot["CSVv2Discr_leadingJet"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_leadingJet", 48, 0.0, 1.2, "CSVv2 discriminant value of leading jet");
  MSPlot["CSVv2Discr_jet2"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet2", 48, 0.0, 1.2, "CSVv2 discriminant value of jet2");
  MSPlot["CSVv2Discr_jet3"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet3", 48, 0.0, 1.2, "CSVv2 discriminant value of jet3");
  MSPlot["CSVv2Discr_jet4"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet4", 48, 0.0, 1.2, "CSVv2 discriminant value of jet4");
  MSPlot["CSVv2Discr_highest"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_highest", 48, 0.0, 1.2, "Highest CSVv2 discriminant value");
  MSPlot["CSVv2Discr_jetNb"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jetNb", 8, -0.5, 7.5, "Jet number (in order of decreasing p_{T}) with highest CSVv2 discriminant value");
  
  MSPlot["min_Mlb"] = new MultiSamplePlot(datasetsMSP, "min_Mlb", 40, 0, 400, "M_{lb} [GeV]");
  MSPlot["dR_Lep_B"] = new MultiSamplePlot(datasetsMSP, "dR_Lep_B", 25, 0, 5, "#Delta R(l,b)");
  
  
  /// Reco
  MSPlot["Reco_W_mass"] = new MultiSamplePlot(datasetsMSP, "Reco_W_mass", 100, 0, 400, "M_{W} [GeV]");
  MSPlot["Reco_hadTop_mass"] = new MultiSamplePlot(datasetsMSP, "Reco_hadTop_mass", 80, 0, 800, "M_{t} [GeV]");
  MSPlot["Reco_hadTop_pT"] = new MultiSamplePlot(datasetsMSP, "Reco_hadTop_pT", 80, 0, 800, "p_{T} [GeV]");
  
  MSPlot["Reco_mlb"] = new MultiSamplePlot(datasetsMSP, "Reco_mlb", 80, 0, 800, "M_{lb} [GeV]");
  MSPlot["Reco_ttbar_mass"] = new MultiSamplePlot(datasetsMSP, "Reco_ttbar_mass", 50, 0, 1000, "M_{t#bar{t}} [GeV]");

  MSPlot["Reco_dR_lep_b_min"] = new MultiSamplePlot(datasetsMSP, "Reco_dR_lep_b_min", 25, 0, 5, "#Delta R(l,b_{l})");
  MSPlot["Reco_dR_lep_b_max"] = new MultiSamplePlot(datasetsMSP, "Reco_dR_lep_b_max", 25, 0, 5, "#Delta R(l,b_{h})");
  
  MSPlot["Reco_Red_top_mass"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
  MSPlot["Reco_Red_top_mass_fewerBins"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass ", 44, 0.2, 2.4, "M_{t}/<M_{t}>");
  
  /// KinFitter
  if (doKinFit)
  {
    MSPlot["KF_Chi2"] = new MultiSamplePlot(datasetsMSP, "Chi2 value of kinFitter", 200, 0, 5, "#chi^{2}");
    MSPlot["KF_Chi2_narrow"] = new MultiSamplePlot(datasetsMSP, "Chi2 value of kinFitter  ", 200, 0, 2, "#chi^{2}");
    MSPlot["KF_Chi2_wide"] = new MultiSamplePlot(datasetsMSP, "Chi2 value of kinFitter ", 200, 0, 20, "#chi^{2}");
    MSPlot["KF_W_mass_orig"] = new MultiSamplePlot(datasetsMSP, "W mass before kinFitter", 250, 0, 500, "m_{W} [GeV]");
    MSPlot["KF_top_mass_orig"] = new MultiSamplePlot(datasetsMSP, "Top mass before kinFitter", 400, 0, 800, "m_{t} [GeV]");
    MSPlot["KF_top_pt_orig"] = new MultiSamplePlot(datasetsMSP, "Top pt before kinFitter", 400, 0, 800, "p_{T} [GeV]");
    MSPlot["KF_W_mass_corr"] = new MultiSamplePlot(datasetsMSP, "W mass after kinFitter", 250, 0, 500, "m_{W,kf} [GeV]");
    MSPlot["KF_top_mass_corr"] = new MultiSamplePlot(datasetsMSP, "Top mass after kinFitter", 400, 0, 800, "m_{t,kf} [GeV]");
    MSPlot["KF_top_pt_corr"] = new MultiSamplePlot(datasetsMSP, "Top pt after kinFitter", 400, 0, 800, "p_{T,kf} [GeV]");
    
    MSPlot["KF_top_mass_orig_short"] = new MultiSamplePlot(datasetsMSP, "KF_top_mass_orig_short", 80, 0, 400, "m_{t} [GeV]");
    MSPlot["KF_top_mass_corr_short"] = new MultiSamplePlot(datasetsMSP, "KF_top_mass_corr_short", 80, 0, 400, "m_{t,kf} [GeV]");
    MSPlot["KF_top_mass_corr_fewerBins"] = new MultiSamplePlot(datasetsMSP, "Top mass after kinFitter ", 100, 90, 250, "m_{t,kf} [GeV]");
  }
}

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  InitHisto1DMatch();
//  if (makePlotsTGraph) InitHisto1DTGraph();
  
  /// SFs
  histo1D["Width_SF"] = new TH1F("Width_SF", "Scale factor to change the ttbar distribution width; width SF", 500, 0, 5);
  
  /// Reco
  histo1D["Reco_W_mass_reco"] = new TH1F("Reco_W_mass_reco","Reconstructed hadronic W mass; M_{W} [GeV]", 80, 0, 800);
  histo1D["Reco_top_mass_reco"] = new TH1F("Reco_top_mass_reco","Reconstructed top mass; M_{t} [GeV]", 80, 0, 800);
    
  /// m_t/<m_t>
  histo1D["Reduced_top_mass_TT_reco_CP"] = new TH1F("Reduced_top_mass_TT_reco_CP","Reduced top mass for matched TT sample (reco, correct top match); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_TT_reco_WPUP"] = new TH1F("Reduced_top_mass_TT_reco_WPUP","Reduced top mass for unmatched TT sample (reco, no top match & wrong permutations); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_TT_reco_UP"] = new TH1F("Reduced_top_mass_TT_reco_UP","Reduced top mass for unmatched TT sample (reco, no top match); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_TT_reco_WP"] = new TH1F("Reduced_top_mass_TT_reco_WP","Reduced top mass for matched TT sample (reco, wrong top match: wrong permutation); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  
  histo1D["Reduced_top_mass_bkgd"] = new TH1F("Reduced_top_mass_bkgd","Reduced top mass for background samples; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_TT"] = new TH1F("Reduced_top_mass_TT","Reduced top mass for TT sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_ST_tW_top"] = new TH1F("Reduced_top_mass_ST_tW_top","Reduced top mass for ST tW top sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_ST_tW_antitop"] = new TH1F("Reduced_top_mass_ST_tW_antitop","Reduced top mass for ST tW antitop sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_ST_t_top"] = new TH1F("Reduced_top_mass_ST_t_top","Reduced top mass for ST t top sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_ST_t_antitop"] = new TH1F("Reduced_top_mass_ST_t_antitop","Reduced top mass for ST t antitop sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_DYJets"] = new TH1F("Reduced_top_mass_DYJets","Reduced top mass for DY+Jets sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_WJets"] = new TH1F("Reduced_top_mass_WJets","Reduced top mass for W+Jets sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  histo1D["Reduced_top_mass_data"] = new TH1F("Reduced_top_mass_data","Reduced top mass for data sample; M_{t}/<M_{t}>", 880, 0.2, 2.4);
  
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
//  histo1D["debugLL_hadr_top_mass_reco"] = new TH1F("debugLL_hadr_top_mass_reco","Reconstructed mass of the hadronic top quark; M_{t} [GeV]", 400, 0, 800);
//  histo1D["debugLL_minMlb_reco"]  = new TH1F("debugLL_minMlb_reco","Minimal reconstructed M_{lb} mass; min(M_{lb}) [GeV]", 400, 0, 800);
//  histo1D["debugLL_ttbar_mass_reco"] = new TH1F("debugLL_ttbar_mass_reco","Reconstructed mass of the top quark pair; M_{t#bar{t}} [GeV]", 500, 0, 1000);
//  histo1D["debugLL_dR_lep_b_lep_reco"]  = new TH1F("debugLL_dR_lep_b_lep_reco","Minimal delta R between the lepton and the leptonic b jet; #Delta R(l,b_{l})", 25, 0, 5);
//  histo1D["debugLL_dR_lep_b_had_reco"]  = new TH1F("debugLL_dR_lep_b_had_reco","Minimal delta R between the lepton and the hadronic b jet; #Delta R(l,b_{h})", 25, 0, 5);
  
  /// KinFitter
  if (doKinFit)
  {
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
    
    histo1D["KF_Chi2_CP"] = new TH1F("KF_Chi2_CP", "Chi2 value of kinFitter (CP); #chi^{2}", 200, 0, 5);
    histo1D["KF_Chi2_WP"] = new TH1F("KF_Chi2_WP", "Chi2 value of kinFitter (WP); #chi^{2}", 200, 0, 5);
    histo1D["KF_Chi2_UP"] = new TH1F("KF_Chi2_UP", "Chi2 value of kinFitter (UP); #chi^{2}", 200, 0, 5);
    histo1D["KF_Chi2_CP_wide"] = new TH1F("KF_Chi2_CP_wide", "Chi2 value of kinFitter (CP); #chi^{2}", 200, 0, 20);
    histo1D["KF_Chi2_WP_wide"] = new TH1F("KF_Chi2_WP_wide", "Chi2 value of kinFitter (WP); #chi^{2}", 200, 0, 20);
    histo1D["KF_Chi2_UP_wide"] = new TH1F("KF_Chi2_UP_wide", "Chi2 value of kinFitter (UP); #chi^{2}", 200, 0, 20);
    
    histo1D["KF_top_mass_orig_match"] = new TH1F("KF_top_mass_orig_match", "Top mass before kinFitter for matched events; m_{t} [GeV]", 80, 0, 400);
    histo1D["KF_top_mass_corr_match"] = new TH1F("KF_top_mass_corr_match", "Top mass after kinFitter for matched events; m_{t,kf} [GeV]", 80, 0, 400);
  }
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

void InitHisto1DMatch()
{
  TH1::SetDefaultSumw2();
  
  histo1D["W_mass_reco_matched"] = new TH1F("W_mass_reco_matched","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 125, 0, 250);
  histo1D["top_mass_reco_matched"] = new TH1F("top_mass_reco_matched","Reconstructed top mass of matched events; M_{t} [GeV]", 175, 50, 400);
  histo1D["top_mass_gen_matched"] = new TH1F("top_mass_gen_matched","Generated top mass of matched events; M_{t} [GeV]", 1200, 150, 190);
  
  histo1D["Reduced_top_mass_TT_matched_partons"] = new TH1F("Reduced_top_mass_TT_matched_partons","Reduced top mass for matched TT sample (using matched partons); M_{t}/<M_{t}>", 800, 0.8, 1.2);
  histo1D["Reduced_top_mass_TT_matched_jets"] = new TH1F("Reduced_top_mass_TT_matched_jets","Reduced top mass for matched TT sample (using jets from matched partons); M_{t}/<M_{t}>", 880, 0.2, 2.4);
  
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
    cutFlow2[i] = 0;
  }
  appliedJER = 999;
  appliedJES = 999;
  appliedPU = 999;
  if (isData)
  {
    nofEventsRunB = 0;
    nofEventsRunCD = 0;
    nofEventsRunEF = 0;
    nofEventsRunG = 0;
    nofEventsRunH = 0;
    nofSelEventsRunB = 0;
    nofSelEventsRunCD = 0;
    nofSelEventsRunEF = 0;
    nofSelEventsRunG = 0;
    nofSelEventsRunH = 0;
  }
  else
  {
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
  }
  
  strSyst = "";
  eqLumi = 1.;
  
  nofHardSelected = 0;
  nofMETCleaned = 0;
  nofMatchedEvents = 0;
  nofHadrMatchedEvents = 0;
  nofCorrectlyMatched = 0;
  nofNotCorrectlyMatched = 0;
  nofAcceptedKFit = 0;
  nofAcceptedKFitMatched = 0;
  
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
  met_pt = 0.;
  met_phi = 0.;
  met_eta = 0.;
  met_Et = 0.;
  met_E = 0.;
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
    
    btagSF = 1.;
    btagSF_up = 1.;
    btagSF_down = 1.;
    puSF = 1.;
    puSF_up = 1.;
    puSF_down = 1.;
    muonIdSF_BCDEF[0] = 1.;
    muonIdSF_GH[0] = 1.;
    muonIdSF_up_BCDEF[0] = 1.;
    muonIdSF_up_GH[0] = 1.;
    muonIdSF_down_BCDEF[0] = 1.;
    muonIdSF_down_GH[0] = 1.;
    muonIsoSF_BCDEF[0] = 1.;
    muonIsoSF_GH[0] = 1.;
    muonIsoSF_up_BCDEF[0] = 1.;
    muonIsoSF_up_GH[0] = 1.;
    muonIsoSF_down_BCDEF[0] = 1.;
    muonIsoSF_down_GH[0] = 1.;
    muonTrigSF_BCDEF[0] = 1.;
    muonTrigSF_GH[0] = 1.;
    muonTrigSF_up_BCDEF[0] = 1.;
    muonTrigSF_up_GH[0] = 1.;
    muonTrigSF_down_BCDEF[0] = 1.;
    muonTrigSF_down_GH[0] = 1.;
    muonTrackSF_eta[0] = 1.;
    muonTrackSF_aeta[0] = 1.;
    muonTrackSF_nPV[0] = 1.;
  }
  else
  {
    isDataRunB = false;
    isDataRunC = false;
    isDataRunD = false;
    isDataRunE = false;
    isDataRunF = false;
    isDataRunG = false;
    isDataRunH = false;
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
  topPtRewSF = 1.;
  topPtSF = 1.;
  antiTopPtSF = 1.;
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
  catSuffix = "";
  isCP = false;
  isWP = false;
  isUP = false;
  kFitVerbosity = false;
  kFitChi2 = 99.;
  kFitChi2Matched = 99.;
  Wmass_reco_orig = -1.;
  Wmass_reco_kf = -1.;
  topmass_reco_orig = -1.;
  topmass_reco_kf = -1.;
  topmass_reco_kf_matched = -1.;
  redTopMass = -1.;
  //toyValue = -1.;
  for (int i = 0; i < nPsExps; i++)
  {
    toyValues[i] = 1.;
  }
}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearVars();
}

void FillGeneralPlots(vector<Dataset *> datasets, int d)
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
  
  
  MSPlot["muon_pT"]->Fill(selectedLepton[0].Pt(), datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["muon_eta"]->Fill(selectedLepton[0].Eta(), datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["muon_phi"]->Fill(selectedLepton[0].Phi(), datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["muon_relIso"]->Fill(muon_relIso[0], datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["muon_d0"]->Fill(muon_d0[0], datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["leadingJet_pT"]->Fill(selectedJets[0].Pt(), datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["jet2_pT"]->Fill(selectedJets[1].Pt(), datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["jet3_pT"]->Fill(selectedJets[2].Pt(), datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["jet4_pT"]->Fill(selectedJets[3].Pt(), datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["Ht_4leadingJets"]->Fill(Ht, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["met_pT"]->Fill(met_pt, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["met_eta"]->Fill(met_eta, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["met_phi"]->Fill(met_phi, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["met_corr_pT"]->Fill(met_corr_pt, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["met_corr_eta"]->Fill(met_corr_eta, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["met_corr_phi"]->Fill(met_corr_phi, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  
  MSPlot["M3"]->Fill(M3, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["min_Mlb"]->Fill(min_Mlb, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["dR_Lep_B"]->Fill(dRLepB, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  
  MSPlot["nJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["nBJets"]->Fill(selectedBJets.size(), datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  
  MSPlot["CSVv2Discr_leadingJet"]->Fill(jet_bdiscr[0], datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["CSVv2Discr_jet2"]->Fill(jet_bdiscr[1], datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["CSVv2Discr_jet3"]->Fill(jet_bdiscr[2], datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["CSVv2Discr_jet4"]->Fill(jet_bdiscr[3], datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  
  int labelB = -1;
  double highestBDiscr = -999.;
  for (int iJet = 0; iJet < selectedJets.size(); iJet++)
  {
    MSPlot["CSVv2Discr_allJets"]->Fill(jet_bdiscr[iJet], datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    if ( jet_bdiscr[iJet] > highestBDiscr )
    {
      highestBDiscr = jet_bdiscr[iJet];
      labelB = iJet;
    }
  }
  MSPlot["CSVv2Discr_highest"]->Fill(highestBDiscr, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["CSVv2Discr_jetNb"]->Fill(labelB, datasets[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
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
    
    histo1D["Reduced_top_mass_TT_matched_partons"]->Fill(matchedTopMass_gen/aveTopMass[0], widthSF);
    histo1D["Reduced_top_mass_TT_matched_jets"]->Fill(matchedTopMass_reco/aveTopMass[1], widthSF);
    
    if (doKinFit)
    {
      histo1D["KF_top_mass_orig_match"]->Fill(matchedTopMass_reco, widthSF);
      histo1D["KF_top_mass_corr_match"]->Fill(topmass_reco_kf_matched, widthSF);
    }
    
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

void FillKinFitPlots()
{
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

void FillCatsPlots(string catSuffix)
{
  if (! isData)
  {
    histo1D["Reduced_top_mass_bkgd"]->Fill(redTopMass);
    if ( isWP || isUP ) histo1D["Reduced_top_mass_TT_reco_WPUP"]->Fill(redTopMass, widthSF);
    histo1D[("Reduced_top_mass_TT_reco"+catSuffix).c_str()]->Fill(redTopMass, widthSF);
    histo1D[("minMlb_reco"+catSuffix).c_str()]->Fill(reco_minMlb, widthSF);
    histo1D[("dR_lep_b_lep_reco"+catSuffix).c_str()]->Fill(reco_dRLepB_lep, widthSF);
    histo1D[("dR_lep_b_had_reco"+catSuffix).c_str()]->Fill(reco_dRLepB_had, widthSF);
    histo1D[("ttbar_mass_reco"+catSuffix).c_str()]->Fill(reco_ttbarMass, widthSF);
    histo2D[("dR_lep_b_lep_vs_had"+catSuffix).c_str()]->Fill(reco_dRLepB_lep, reco_dRLepB_had, widthSF);
    histo2D[("ttbar_mass_vs_minMlb"+catSuffix).c_str()]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
    if ( reco_dRLepB_lep < 3. )
    {
      histo2D[("ttbar_mass_vs_minMlb_dRlepCut"+catSuffix).c_str()]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
      if ( reco_dRLepB_had > 1.2 )
        histo2D[("ttbar_mass_vs_minMlb_dRBothCuts"+catSuffix).c_str()]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
      
      if ( reco_dRLepB_lep < 2. )
      {
        histo2D[("ttbar_mass_vs_minMlb_dRlepCutHard"+catSuffix).c_str()]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
        if ( reco_dRLepB_had > 2. )
          histo2D[("ttbar_mass_vs_minMlb_dRBothCutsHard"+catSuffix).c_str()]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
      }
    }
    if ( reco_dRLepB_had > 1.2 )
    {
      histo2D[("ttbar_mass_vs_minMlb_dRhadCut"+catSuffix).c_str()]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
      if ( reco_dRLepB_had > 2. )
        histo2D[("ttbar_mass_vs_minMlb_dRhadCutHard"+catSuffix).c_str()]->Fill(reco_minMlb, reco_ttbarMass, widthSF);
    }
    if (doKinFit)
    {
      histo1D[("KF_Chi2"+catSuffix).c_str()]->Fill(kFitChi2);
      histo1D[("KF_Chi2"+catSuffix+"_wide").c_str()]->Fill(kFitChi2);
      histo1D[("KF_top_mass_corr"+catSuffix).c_str()]->Fill(topmass_reco_kf, widthSF);
    }
  }
}

void FillMSPlots(int d)
{
  FillGeneralPlots(datasetsMSP, d);
  
  MSPlot["Reco_W_mass"]->Fill(Wmass_reco_kf, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["Reco_hadTop_mass"]->Fill(topmass_reco_kf, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["Reco_hadTop_pT"]->Fill(toppt_reco_kf, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["Reco_mlb"]->Fill(reco_minMlb, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["Reco_ttbar_mass"]->Fill(reco_ttbarMass, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["Reco_dR_lep_b_min"]->Fill(reco_dRLepB_lep, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  MSPlot["Reco_dR_lep_b_max"]->Fill(reco_dRLepB_had, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  if ( redTopMass < 2.4 )
  {
    MSPlot["Reco_Red_top_mass"]->Fill(redTopMass, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    MSPlot["Reco_Red_top_mass_fewerBins"]->Fill(redTopMass, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
  }
  
  if (doKinFit)
  {
    if ( kFitChi2 < 5 )
    {
      MSPlot["KF_Chi2"]->Fill(kFitChi2, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
      if ( kFitChi2 < 2 )
        MSPlot["KF_Chi2_narrow"]->Fill(kFitChi2, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    }
    MSPlot["KF_Chi2_wide"]->Fill(kFitChi2, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    MSPlot["KF_W_mass_orig"]->Fill(Wmass_reco_orig, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    MSPlot["KF_W_mass_corr"]->Fill(Wmass_reco_kf, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    MSPlot["KF_top_mass_orig"]->Fill(topmass_reco_orig, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    MSPlot["KF_top_mass_corr"]->Fill(topmass_reco_kf, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    if ( topmass_reco_kf < 250 && topmass_reco_kf > 90 )
      MSPlot["KF_top_mass_corr_fewerBins"]->Fill(topmass_reco_kf, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    MSPlot["KF_top_pt_orig"]->Fill(toppt_reco_orig, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    MSPlot["KF_top_pt_corr"]->Fill(toppt_reco_kf, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    
    MSPlot["KF_top_mass_orig_short"]->Fill(topmass_reco_orig, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
    MSPlot["KF_top_mass_corr_short"]->Fill(topmass_reco_kf, datasetsMSP[d], true, Luminosity/eqLumi*scaleFactor*widthSF);
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

void CheckSystematics(vector<int> vJER, vector<int> vJES, vector<int> vPU)
{
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
  else
  {
    cerr << "Shape changing systematics not consistent accross datasets or multiple applied at once" << endl;
    cerr << "Exiting...." << endl;
    exit(1);
  }
  cout << " - Systematics confirmed to be " << strSyst << endl;
}

void PrintKFDebug(int ievt)
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

