#include <TStyle.h>
#include <TPaveText.h>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <array>
#include <TRandom3.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TGraph.h>
#include <TGraph2D.h>

// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
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
bool makePlots = false;
bool fourJetsOnly = false;

bool doKinFit = true;
bool applyKinFitCut = true;
double kinFitCutValue = 5.;

bool doMETCleaning = true;

string systStr = "nominal";
pair<string,string> whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos )
  {
    return pair<string,string>("171015","170921");
  }
  else if ( syst.find("JESup") != std::string::npos ) return pair<string,string>("170922","170921");
  else if ( syst.find("JESdown") != std::string::npos ) return pair<string,string>("170923","170921");
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    return pair<string,string>("171015","171021");
  }
}
pair<string,string> ntupleDate = whichDate(systStr);
int verbose = 2;

string pathNtuples = "";
string pathNtuplesMC = "";
string pathNtuplesData = "";
string pathNtuplesSyst = "";
string pathOutput = "";

int nofHardSelected = 0;
int nofMETCleaned = 0;
int nofAfterDRmincut = 0;
int nofMatchedEvents = 0;
int nofHadrMatchedEvents = 0;
int nofHadrMatchedEventsAKF = 0;
int nofCorrectlyMatched = 0;
int nofNotCorrectlyMatched = 0;
int nofUnmatched = 0;
int nofCorrectlyMatchedAKF = 0;
int nofNotCorrectlyMatchedAKF = 0;
int nofUnmatchedAKF = 0;
int nofCorrectlyMatchedAKFNoCut = 0;
int nofNotCorrectlyMatchedAKFNoCut = 0;
int nofUnmatchedAKFNoCut = 0;

int corrMatchHadrB = 0;

int nofCMnJets[4] = {0};
int nofWMnJets[4] = {0};
int nofUMnJets[4] = {0};
int nofCMnJetsAKF[4] = {0};
int nofWMnJetsAKF[4] = {0};
int nofUMnJetsAKF[4] = {0};

/// Lumi per data era
double lumi_runBCDEF = 19.67550334113;  // 1/fb
double lumi_runGH = 16.146177597883;  // 1/fb
double Luminosity = (lumi_runBCDEF + lumi_runGH)*1000;  // 1/pb

///  Working points for b tagging  // Updated 13/04/17, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
double CSVv2Loose  = 0.5426;
double CSVv2Medium = 0.8484;
double CSVv2Tight  = 0.9535;


// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TGraph*> graph;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;

vector < Dataset* > datasets;


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
void GetMetaData(TTree* tree);
void InitTree(TTree* tree);
void InitHisto1D();
void InitHisto2D();
void TruthMatching(vector<TLorentzVector> partons, vector<TLorentzVector> selectedJets, pair<unsigned int, unsigned int> *MCPermutation);
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearVars();
void ClearObjects();
long GetNEvents(TTree* fChain, string var);
long GetNEvents(TTree* fChain, string var, unsigned int index);
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
// Bool_t          filterHBHENoise;
// Bool_t          filterHBHEIso;
// Bool_t          filterCSCTightHalo;
// Bool_t          filterEcalDeadCell;
// Bool_t          filterEEBadSc;
// Bool_t          filterBadChCand;
// Bool_t          filterBadMuon;
Bool_t          passedMETFilter = true;
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
// Bool_t          hasGenTopWithStatus22;
// Bool_t          hasGenTopWithStatus62;
Bool_t          hasGenAntiTop;
// Bool_t          hasGenAntiTopWithStatus22;
// Bool_t          hasGenAntiTopWithStatus62;
Double_t        weight1001;
Double_t        weight1002;
Double_t        weight1003;
Double_t        weight1004;
Double_t        weight1005;
Double_t        weight1007;
Double_t        weight1009;
// Double_t        upFragWeight;
// Double_t        centralFragWeight;
// Double_t        downFragWeight;
// Double_t        petersonFragWeight;
// Double_t        semilepbrUp;
// Double_t        semilepbrDown;
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
Double_t        cutFlowWeighted[10];
Double_t        cutFlow2Weighted[10];
Int_t           appliedJER;
Int_t           appliedJES;
Int_t           appliedPU;
// Long64_t        nofEventsWithGenTop;
// Long64_t        nofEventsWithGenTopWithStatus22or62;
// Long64_t        nofEventsWithGenAntiTop;
// Long64_t        nofEventsWithGenAntiTopWithStatus22or62;
// Long64_t        nofTTEventsWithoutBothGenTops;
Long64_t        nofTTEventsWithoutAGenTop;
// Long64_t        nofTTEventsWithoutGenTop;
// Long64_t        nofTTEventsWithoutGenAntiTop;
// Long64_t        nofTTEventsWithoutBothGenTopsWithStatus22;
// Long64_t        nofTTEventsWithoutGenTopWithStatus22;
// Long64_t        nofTTEventsWithoutGenAntiTopWithStatus22;
// Long64_t        nofTTEventsWithoutBothGenTopsWithStatus62;
// Long64_t        nofTTEventsWithoutGenTopWithStatus62;
// Long64_t        nofTTEventsWithoutGenAntiTopWithStatus62;
Double_t        sumWeight1001;
Double_t        sumWeight1002;
Double_t        sumWeight1003;
Double_t        sumWeight1004;
Double_t        sumWeight1005;
Double_t        sumWeight1007;
Double_t        sumWeight1009;

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
// TBranch        *b_filterHBHENoise;   //!
// TBranch        *b_filterHBHEIso;   //!
// TBranch        *b_filterCSCTightHalo;   //!
// TBranch        *b_filterEcalDeadCell;   //!
// TBranch        *b_filterEEBadSc;   //!
// TBranch        *b_filterBadChCand;   //!
// TBranch        *b_filterBadMuon;   //!
// TBranch        *b_passedMETFilter;   //!
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
// TBranch        *b_hasGenTopWithStatus22;   //!
// TBranch        *b_hasGenTopWithStatus62;   //!
TBranch        *b_hasGenAntiTop;   //!
// TBranch        *b_hasGenAntiTopWithStatus22;   //!
// TBranch        *b_hasGenAntiTopWithStatus62;   //!
TBranch        *b_weight1001;   //!
TBranch        *b_weight1002;   //!
TBranch        *b_weight1003;   //!
TBranch        *b_weight1004;   //!
TBranch        *b_weight1005;   //!
TBranch        *b_weight1007;   //!
TBranch        *b_weight1009;   //!
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
TBranch        *b_cutFlowWeighted;   //!
TBranch        *b_cutFlow2Weighted;   //!
TBranch        *b_appliedJER;   //!
TBranch        *b_appliedJES;   //!
TBranch        *b_appliedPU;   //!
// TBranch        *b_nofEventsWithGenTop;   //!
// TBranch        *b_nofEventsWithGenTopWithStatus22or62;   //!
// TBranch        *b_nofEventsWithGenAntiTop;   //!
// TBranch        *b_nofEventsWithGenAntiTopWithStatus22or62;   //!
// TBranch        *b_nofTTEventsWithoutBothGenTops;   //!
TBranch        *b_nofTTEventsWithoutAGenTop;   //!
// TBranch        *b_nofTTEventsWithoutGenTop;   //!
// TBranch        *b_nofTTEventsWithoutGenAntiTop;   //!
// TBranch        *b_nofTTEventsWithoutBothGenTopsWithStatus22;   //!
// TBranch        *b_nofTTEventsWithoutGenTopWithStatus22;   //!
// TBranch        *b_nofTTEventsWithoutGenAntiTopWithStatus22;   //!
// TBranch        *b_nofTTEventsWithoutBothGenTopsWithStatus62;   //!
// TBranch        *b_nofTTEventsWithoutGenTopWithStatus62;   //!
// TBranch        *b_nofTTEventsWithoutGenAntiTopWithStatus62;   //!
TBranch        *b_sumWeight1001;   //!
TBranch        *b_sumWeight1002;   //!
TBranch        *b_sumWeight1003;   //!
TBranch        *b_sumWeight1004;   //!
TBranch        *b_sumWeight1005;   //!
TBranch        *b_sumWeight1007;   //!
TBranch        *b_sumWeight1009;   //!


long nEventsDataSet;
double xSection;
double lumiWeight, scaleFactor;
bool foundTop22, foundAntiTop22;
bool foundTop62, foundAntiTop62;
bool foundLastCopyTop, foundLastCopyAntitop;
vector<unsigned int> bJetId;
double bdiscrTop, bdiscrTop2, tempbdiscr;
int labelB1, labelB2;
int labelsReco[4], labelsRecoKF[4];
double massHadTopQ, massLepTopQ;

string catSuffix = "";
string catSuffixList[] = {"_CM", "_WM", "_UM"};
bool isCM, isWM, isUM;


/// Define TLVs
TLorentzVector muon, jet, mcpart;
TLorentzVector WCandidate;
TLorentzVector neutrino;
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
bool muonmatched = false;
bool foundMuPlus = false, foundMuMinus = false;
bool muPlusFromTop = false, muMinusFromTop = false;
vector<unsigned int> partonId;


/// KinFitter
bool doneKinFit = false;
TKinFitter* kFitter01;
TKinFitter* kFitter02;
TKinFitter* kFitter12;
bool addWMassKF = true;
bool addEqMassKF = false;
int kFitVerbosity = 0;
double kFitChi2_01 = 99., kFitChi2_02 = 99., kFitChi2_12 = 99., kFitChi2_min = 99.;
int labelsKF[3] = {-9999, -9999, -9999};
bool min01, min02, min12;
int nofAcceptedKFit = 0;
double nofAcceptedKFitWeighted = 0.;

Double_t minTopMass = 100., maxTopMass = 245.;


/// Variables
double M3, Ht, min_Mlb, dRLepB;
double M3_aKF, Ht_aKF;
double reco_W_mass_bKF, reco_top_mass_bKF, reco_top_pt_bKF, reco_mlb_bKF, reco_dRLepB_lep_bKF, reco_dRLepB_had_bKF, reco_ttbar_mass_bKF;
double reco_W_mass_aKF, reco_top_mass_aKF, reco_top_mass_alt_aKF, reco_top_pt_aKF, reco_mlb_aKF, reco_dRLepB_lep_aKF, reco_dRLepB_had_aKF, reco_ttbar_mass_aKF;
double tempDR;

double matched_W_mass_q, matched_top_mass_q;
double matched_W_mass_j, matched_top_mass_j, matched_top_mass_j_akF;
double matched_mlb_corr, matched_ttbarMass_corr, matched_dR_lep_b_corr;
double matched_mlb_wrong, matched_ttbarMass_wrong, matched_dR_lep_b_wrong;


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
  else                                                   cout << "nominal        *" << endl;
  cout << "*   - Jet/lepton Cleaning                   *" << endl;
  if (doMETCleaning) cout << "*   - MET Cleaning                          *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  string channel;
  if ( argc == 1 ) channel = "mu";
  //else if ( CharSearch(argv[1], "mu") || CharSearch(argv[1], "Mu") || CharSearch(argv[1], "MU") || CharSearch(argv[1], "m") ) channel = "mu";
  //else if ( CharSearch(argv[1], "el") || CharSearch(argv[1], "El") || CharSearch(argv[1], "EL") || CharSearch(argv[1], "e") ) channel = "el";
  //else if ( (argv[1]).find("all") != std::string::npos || (argv[1]).find("All") != std::string::npos || (argv[1]).find("ALL") != std::string::npos ) channel = "all";
  
  if (test) makePlots = false;
  if (testHistos) makePlots = true;
  
  pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (makePlots)
  {
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
  cout << "Using Ntuples from " << ntupleDate.first << " for MC. This corresponds to systematics: " << systStr << endl;
  if (testHistos) cout << "Testing histogram consistency..." << endl;
  
  
  
  /// xml file
  string xmlFileName ="config/topWidth.xml";
  
  const char *xmlfile = xmlFileName.c_str();
  
  cout << " - Using config file " << xmlfile << endl;
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  EventReweighting *rew = new EventReweighting(false);  // no correction for number of events
  ResolutionFunctions* rf = new ResolutionFunctions(false, true);
  KinFitter *kf01 = new KinFitter("input/PlotsForResolutionFunctions_testFit_171023.root", addWMassKF, addEqMassKF);
  KinFitter *kf02 = new KinFitter("input/PlotsForResolutionFunctions_testFit_171023.root", addWMassKF, addEqMassKF);
//  KinFitter *kf12 = new KinFitter("input/PlotsForResolutionFunctions_testFit_170915.root", addWMassKF, addEqMassKF);
  
  
  if (makePlots)
  {
    InitHisto1D();
    InitHisto2D();
  }
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName;
  int nEntries;
  bool skipEvent = false;
  
  /// Load original Ntuple
  cout << "  - Loading original Ntuples..." << endl;
  TFile* origNtuple = new TFile((pathNtuples+"Ntuples_TT_nominal.root").c_str(), "READ");
  
  string tTreeName = "tree";
  string tStatsTreeName = "stats";
  TTree* origTree = (TTree*) origNtuple->Get(tTreeName.c_str());
  TTree* origStatsTree = (TTree*) origNtuple->Get(tStatsTreeName.c_str());
  
  // Set branch addresses and branch pointers
  InitTree(origTree);
  GetMetaData(origStatsTree);
  
  nEntries = (int)origTree->GetEntries();
  cout << "                nEntries: " << nEntries << endl;
  
  
  
  ////////////////////////////////////
  ///  Loop on events
  ////////////////////////////////////
  
  int endEvent = nEntries;
  if (test || testHistos) endEvent = 100;
  for (int ievt = 0; ievt < endEvent; ievt++)
  {
    ClearObjects();
    
    if (ievt%10000 == 0)
      std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
    
    
    /// Load event
    origTree->GetEntry(ievt);
    
    
    
    //////////////////////
    ///  Fill objects  ///
    //////////////////////
    
    muon.SetPtEtaPhiE(muon_pt[0], muon_eta[0], muon_phi[0], muon_E[0]);
    selectedLepton.push_back(muon);
    
    for (int iJet = 0; iJet < nJets; iJet++)
    {
      jet.Clear();
      jet.SetPtEtaPhiE(jet_pt[iJet], jet_eta[iJet], jet_phi[iJet], jet_E[iJet]);
      if ( jet_pt[iJet] < 250. ) selectedJets.push_back(jet);
    }
    
    if ( selectedJets.size() < 4 ) continue;
    
    for (int iJet = 0; iJet < selectedJets.size(); iJet++)
    {
      if ( jet_bdiscr[iJet] > CSVv2Medium )
      {
        selectedBJets.push_back(selectedJets[iJet]);
        bJetId.push_back(iJet);  /// selectedBJets[j] = selectedJets[bJetId[j]]
      }
    }
    
    if ( selectedBJets.size() < 2 ) continue;
    
    
    if ( fourJetsOnly && selectedJets.size() > 4 ) continue;
    //if ( selectedJets.size() > 5 ) continue;
    nofHardSelected++;
    
    if (! passedMETFilter) continue;
    nofMETCleaned++;
    
    skipEvent = false;
    for (int iJet = 0; iJet < selectedJets.size(); iJet++)
    {
      for (int jJet = iJet+1; jJet < selectedJets.size(); jJet++)
      {
        tempDR = ROOT::Math::VectorUtil::DeltaR(selectedJets[iJet], selectedJets[jJet]);
        if ( tempDR < 0.6 )
        {
          skipEvent = true;
          break;
        }
      }
      if (skipEvent) break;
    }
    if (skipEvent) continue;
    nofAfterDRmincut++;
    
    
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
    
    
    
    
    /////////////////////////////
    ///  JET PARTON MATCHING  ///
    /////////////////////////////
    
    for (int iMC = 0; iMC < nMCParticles; iMC++)
    {
      mcpart.Clear();
      mcpart.SetPtEtaPhiE(mc_pt[iMC], mc_eta[iMC], mc_phi[iMC], mc_E[iMC]);
      mcParticles.push_back(mcpart);
    }
    
    foundTop22 = false; foundAntiTop22 = false; foundTop62 = false; foundAntiTop62 = false;
    foundLastCopyTop = false; foundLastCopyAntitop = false;
    for (unsigned int i = 0; i < mcParticles.size(); i++)
    {
      if ( test && verbose > 4 )
        cout << setw(3) << right << i << "  Status: " << setw(2) << mc_status[i] << "  pdgId: " << setw(3) << mc_pdgId[i] << "  Mother: " << setw(4) << mc_mother[i] << "  Granny: " << setw(4) << mc_granny[i] << "  Pt: " << setw(7) << left << mc_pt[i] << "  Eta: " << mc_eta[i] << endl;
      
      /// Find tops
      if ( mc_pdgId[i] == pdgID_top )  // isLastCopy() == status 62
      {
        if ( mc_status[i] == 22 )
        {
          topQuark = i;
          foundTop22 = true;
        }
        if ( topQuark == -9999 && mc_status[i] == 62 ) topQuark = i;
      }
      else if ( mc_pdgId[i] == -pdgID_top )
      {
        if ( mc_status[i] == 22 )
        {
          antiTopQuark = i;
          foundAntiTop22 = true;
        }
        if ( antiTopQuark == -9999 && mc_status[i] == 62 ) antiTopQuark = i;
      }
      
    
      
      
      /// Status restriction: Final state particle or particle from hardest process
      if ( (mc_status[i] > 1 && mc_status[i] <= 20) || mc_status[i] >= 30 ) continue;
      
      /// Muons   // from top or from W (in ST tW channel)
      if ( mc_pdgId[i] == 13 && mc_mother[i] == -24 )		// mu-, W-
      {
        foundMuMinus = true;
        if ( mc_granny[i] == -pdgID_top ) muMinusFromTop = true;  // t~
        
        if ( mc_status[i] == 23 ) genmuon = i;
        else if ( genmuon == -9999 ) genmuon = i;
      }
      if ( mc_pdgId[i] == -13 && mc_mother[i] == 24 )		// mu+, W+
      {
        foundMuPlus = true;
        if ( mc_granny[i] == pdgID_top ) muPlusFromTop = true;  // t
        
        if ( mc_status[i] == 23 ) genmuon = i;
        else if ( genmuon == -9999 ) genmuon = i;
      }
      
      /// Partons/gluons
      if ( abs(mc_pdgId[i]) < 6 /*|| abs(mc_pdgId[i]) == 21*/ )  //light/b quarks, 6 should stay hardcoded, OR gluon
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
    else if ( foundMuMinus && foundMuPlus )
    {
      if (test) cout << "Found fully leptonic decay of ST tW... Event " << ievt << " will not be matched." << endl;
      doMatching = false;
    }
    
    
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
        
      }  // end hadronicTopJetsMatched
      
    }  // end doMatching
    
    
    
    
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
    
//    if ( labelsReco[3] == -9999 ) continue;
    
    if ( labelsReco[3] == bJetId[labelB1] ) labelsReco[2] = bJetId[labelB2];
    else if ( labelsReco[3] == bJetId[labelB2] ) labelsReco[2] = bJetId[labelB1];
    else cerr << endl << "Seems like something went wrong with the b jets..." << endl;
    
    for (int ijet = 0; ijet < selectedJets.size(); ijet++)
    {
      if ( ijet == labelsReco[2] || ijet == labelsReco[3] ) continue;
      
      if ( labelsReco[0] == -9999 ) labelsReco[0] = ijet;
      else if ( labelsReco[1] == -9999 ) labelsReco[1] = ijet;
//      else cerr << endl << "Seems like there are too many jets..." << endl;
    }
    
//    if ( labelsReco[0] == -9999 || labelsReco[1] == -9999 || labelsReco[2] == -9999 ) continue;
    
//     if ( selectedJets.size() > 4)
//     {
//       for (int ijet = 0; ijet < selectedJets.size(); ijet++)
//       {
//         for (int jjet = ijet; jjet < selectedJets.size(); jjet++)
//         {
//           for (int kjet = 0; kjet < selectedJets.size(); kjet++)
//           {
//             if ( kjet == ijet || kjet == jjet ) continue;
//             //if ( jet_bdiscr[kjet] < CSVv2Medium ) continue;
// 
//             for (int ljet = 0; ljet < selectedJets.size(); ljet++)
//             {
//               if ( ljet == ijet || ljet == jjet || ljet == kjet ) continue;
//               //if ( jet_bdiscr[ljet] < CSVv2Medium ) continue;
// 
//               deltaR = ROOT::Math::VectorUtil::DeltaR(selectedJets[ijet], selectedJets[jjet]) + ROOT::Math::VectorUtil::DeltaR(selectedJets[ijet], selectedJets[kjet]) + ROOT::Math::VectorUtil::DeltaR(selectedJets[jjet], selectedJets[kjet]) + ROOT::Math::VectorUtil::DeltaR(selectedJets[ljet], selectedLepton[0]);
// 
//               if (deltaR < minDeltaR)
//               {
//                 minDeltaR = deltaR;
//                 labelsReco[0] = ijet;
//                 labelsReco[1] = jjet;
//                 labelsReco[2] = kjet;
//                 labelsReco[3] = ljet;
//               }
//             }
//           }
//         }
//       }
//     }
    
    if ( labelsReco[0] == -9999 || labelsReco[1] == -9999 || labelsReco[2] == -9999 || labelsReco[3] == -9999 ) continue;
    
    /// Fill variables before performing kinFit
    reco_W_mass_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).M();
    reco_top_mass_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).M();
    reco_top_pt_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).Pt();
    reco_mlb_bKF = (selectedLepton[0] + selectedJets[labelsReco[3]]).M();
    reco_dRLepB_lep_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[3]], selectedLepton[0] );  // deltaR between lepton and leptonic b jet
    reco_dRLepB_had_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
    reco_ttbar_mass_bKF = reco_mlb_bKF + reco_top_mass_bKF;
    
    //if ( fabs(reco_top_mass_bKF - reco_mlb_bKF) < 20 ) continue;
    
    
    
    ///////////////////////////////////
    ///  CHECK MATCHED COMBINATION  ///
    ///////////////////////////////////
    
    ///
    // 3 possibilities:
    // - correct top match: 3 jets selected with reco method correspond to the 3 matched jets (n.b. this is also true when the jets originating from the W boson and the b jet do not exactly correspond to the matched jets, because we are only interested in the reconstructed top quark.)
    // - wrong permutation: the correct jet combination exists in the selected jets, but is not chosen by the reco method.
    // - wrong (no) match:  the correct jet combination does not exist in the selected jets (e.g. when one jet is not selected.)
    
    
    if (hadronicTopJetsMatched)
    {
      /// Correct match
      if ( ( labelsReco[0] == MCPermutation[0].first || labelsReco[0] == MCPermutation[1].first || labelsReco[0] == MCPermutation[2].first ) && ( labelsReco[1] == MCPermutation[0].first || labelsReco[1] == MCPermutation[1].first || labelsReco[1] == MCPermutation[2].first ) && ( labelsReco[2] == MCPermutation[0].first || labelsReco[2] == MCPermutation[1].first || labelsReco[2] == MCPermutation[2].first ) )  // correct jets for top quark
      {
        isCM = true;
        nofCorrectlyMatched++;
        if ( labelsReco[2] == MCPermutation[2].first ) corrMatchHadrB++;
        if ( selectedJets.size() == 4 ) nofCMnJets[0]++;
        else if ( selectedJets.size() == 5 ) nofCMnJets[1]++;
        else if ( selectedJets.size() == 6 ) nofCMnJets[2]++;
        else nofCMnJets[3]++;
      }
      else  // wrong permutation
      {
        isWM = true;
        nofNotCorrectlyMatched++;
        if ( selectedJets.size() == 4 ) nofWMnJets[0]++;
        else if ( selectedJets.size() == 5 ) nofWMnJets[1]++;
        else if ( selectedJets.size() == 6 ) nofWMnJets[2]++;
        else nofWMnJets[3]++;
        if (test) cout << "Event " << ievt << ": (" << MCPermutation[0].first << "," << MCPermutation[1].first << "," << MCPermutation[2].first << "," << MCPermutation[3].first << ")" << endl;
      }
    }  // end hadrTopMatch
    else  // no match
    {
      isUM = true;
      nofUnmatched++;
      if ( selectedJets.size() == 4 ) nofUMnJets[0]++;
      else if ( selectedJets.size() == 5 ) nofUMnJets[1]++;
      else if ( selectedJets.size() == 6 ) nofUMnJets[2]++;
      else nofUMnJets[3]++;
    }
    
    
    if ( (! isCM && ! isWM && ! isUM) || (isCM && isWM) || (isCM && isUM) || (isWM && isUM) )
      cerr << "Something wrong with trigger logic CM/WM/UM !! " << endl;
    
    
    if (isCM) catSuffix = catSuffixList[0];
    else if (isWM) catSuffix = catSuffixList[1];
    else if (isUM) catSuffix = catSuffixList[2];
    
    
    
    ////////////////////////////
    ///   Kinematic Fit      ///
    ////////////////////////////
    
    for (int i = 0; i < 4; i++)
    {
      labelsRecoKF[i] = labelsReco[i];
    }
    
    if (doKinFit)
    {
      if (addEqMassKF)
      {
        cout << "Event " << ievt << ": original combination = (" << labelsReco[0] << "," << labelsReco[1] << "," << labelsReco[2] << "," << labelsReco[3] << ") ====> ";
        if (isCM) cout << "CM";
        else if (isWM) cout << "WM";
        else cout << "UM";
        cout << endl;
        neutrino.SetPtEtaPhiE(met_corr_pt, 1e-10, met_corr_phi, met_corr_E);
        //cout << "Neutrino pt = " << met_corr_pt << " and in TLV: " << neutrino.Pt() << "    Et = " << met_corr_Et << " and in TLV: " << neutrino.Et() << "   E = " << met_corr_E << " and in TLV: " << neutrino.E() << endl;
        double tempProb, maxProb = 999.;
        for ( int ijet = 0; ijet < selectedJets.size(); ijet++)
        {
          for (int jjet = ijet+1; jjet < selectedJets.size(); jjet++)
          {
            //if ( jjet == ijet ) continue;
            for (int kjet = 0; kjet < selectedJets.size(); kjet++)
            {
              if ( kjet == ijet || kjet == jjet ) continue;
              for (int ljet = 0; ljet < selectedJets.size(); ljet++)
              {
                if ( ljet == ijet || ljet == jjet || ljet == kjet ) continue;
                
                kFitter02 = kf02->doFit(selectedJets[ijet], selectedJets[jjet], selectedJets[kjet], selectedJets[ljet], selectedLepton[0], neutrino, kFitVerbosity);
                
                if ( kFitter02->getStatus() != 0 ) continue;
                //tempProb = TMath::Log(TMath::Prob(kFitter02->getS(),kFitter02->getNDF()));
                //if ( tempProb > maxProb)
                tempProb = fabs(kFitter02->getS()/kFitter02->getNDF() - 1.);
                if ( tempProb < maxProb )
                {
                  maxProb = tempProb;
                  labelsRecoKF[0] = ijet;
                  labelsRecoKF[1] = jjet;
                  labelsRecoKF[2] = kjet;
                  labelsRecoKF[3] = ljet;
                  if (test)
                    cout << "Event " << ievt << ": for combination (" << labelsRecoKF[0] << "," << labelsRecoKF[1] << "," << labelsRecoKF[2] << "," << labelsRecoKF[3] << "), the chi2 value is " << kFitter02->getS() << " and the ndf is " << kFitter02->getNDF() << " ===> " << TMath::Prob(kFitter02->getS(),kFitter02->getNDF()) << endl;
                }
              }
            }
          }
        }
        
        if ( maxProb == 999. ) continue;
        kFitter01 = kf01->doFit(selectedJets[labelsRecoKF[0]], selectedJets[labelsRecoKF[1]], selectedJets[labelsRecoKF[2]], selectedJets[labelsRecoKF[3]], selectedLepton[0], neutrino, kFitVerbosity);
        
        //kFitter01 = kf01->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]], selectedJets[labelsReco[2]], selectedJets[labelsReco[3]], selectedLepton[0], neutrino, kFitVerbosity);
        //kFitter02 = kf02->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]], selectedJets[labelsReco[3]], selectedJets[labelsReco[2]], selectedLepton[0], neutrino, kFitVerbosity);
      }
      else if (addWMassKF)
      {
        kFitter01 = kf01->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]], kFitVerbosity);
//       kFitter02 = kf02->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[2]], kFitVerbosity);
//       kFitter12 = kf12->doFit(selectedJets[labelsReco[1]], selectedJets[labelsReco[2]], kFitVerbosity);
      
       if ( kFitter01->getStatus() != 0 /*&& kFitter02->getStatus() != 0 && kFitter12->getStatus() != 0 */)  // did not converge
       {
         if (test && verbose > 2) cout << "Event " << ievt << ": Fit did not converge..." << endl;
         continue;
       }
     }
     
//       if (kFitter01->getStatus() == 0 )
//       {
//         kFitChi2_01 = kFitter01->getS();
//         if (test)
//         {
//           cout << "Event " << ievt << " (";
//           if (isCM) cout << "CM";
//           else if (isWM) cout << "WM";
//           else cout << "UM";
//           cout << "): for combination 1, the chi2 value is " << kFitChi2_01 << " and the ndf is " << kFitter01->getNDF() << " ===> " << TMath::Prob(kFitChi2_01,kFitter01->getNDF()) << endl;
//         }
//       }
//       if (kFitter02->getStatus() == 0 )
//       {
//         kFitChi2_02 = kFitter02->getS();
//         if (test) cout << "Event " << ievt << "     : for combination 2, the chi2 value is " << kFitChi2_02 << " and the ndf is " << kFitter02->getNDF() << " ===> " << TMath::Prob(kFitChi2_02,kFitter02->getNDF()) << endl;
//       }
// //       kFitChi2_02 = kFitter02->getS();
// //       kFitChi2_12 = kFitter12->getS();
//       min01 = false; min02 = false; min12 = false;
// //       if ( kFitter01->getStatus() == 0 && kFitChi2_01 < kFitChi2_min )
// //       {
// //         min01 = true;
//       if (kFitter01->getStatus() == 0 && kFitter02->getStatus() == 0)
//       {
//         kFitChi2_min = kFitChi2_01;
//         if ( kFitChi2_02 < kFitChi2_01 )
//         {
//           min02 = true;
//           kFitChi2_min = kFitChi2_02;
//         }
//       }
//       else if (kFitter02->getStatus() == 0)
//       {
//         min02 = true;
//         kFitChi2_min = kFitChi2_02;
//       }
//       else kFitChi2_min = kFitChi2_01;

//         for (int i = 0; i < 3; i++) labelsKF[i] = labelsReco[i];
//       }
//       if ( kFitter02->getStatus() == 0 && kFitChi2_02 < kFitChi2_min )
//       {
//         min01 = false; min02 = true;
//         kFitChi2_min = kFitChi2_02;
//         labelsKF[0] = labelsReco[0];
//         labelsKF[1] = labelsReco[2];
//         labelsKF[2] = labelsReco[1];
//       }
//       if ( kFitter12->getStatus() == 0 && kFitChi2_12 < kFitChi2_min )
//       {
//         min01 = false; min02 = false; min12 = true;
//         kFitChi2_min = kFitChi2_12;
//         labelsKF[0] = labelsReco[1];
//         labelsKF[1] = labelsReco[2];
//         labelsKF[2] = labelsReco[0];
//       }
      
      /// Put correct order back into labelsReco
      //for (int i = 0; i < 3; i++) labelsReco[i] = labelsKF[i];
      if (test && verbose > 4) cout << "Fit converged: Chi2 = " << kFitChi2_min << endl;
      
      doneKinFit = true;
      
      if (isCM) nofCorrectlyMatchedAKFNoCut++;
      else if (isWM) nofNotCorrectlyMatchedAKFNoCut++;
      else if (isUM) nofUnmatchedAKFNoCut++;
      
      kFitChi2_min = kFitter01->getS();
      if ( applyKinFitCut && kFitChi2_min > kinFitCutValue ) continue;
      nofAcceptedKFit++;
      nofAcceptedKFitWeighted += scaleFactor;
      if (hadronicTopJetsMatched) nofHadrMatchedEventsAKF++;
      if (isCM) nofCorrectlyMatchedAKF++;
      else if (isWM) nofNotCorrectlyMatchedAKF++;
      else if (isUM) nofUnmatchedAKF++;
      
      selectedJetsKFcorrected.clear();
//       if (min12) selectedJetsKFcorrected = kf12->getCorrectedJets();
//       else if (min02) selectedJetsKFcorrected = kf02->getCorrectedJets();
//       else 
      selectedJetsKFcorrected = kf01->getCorrectedJets();
      
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
    if (addEqMassKF)
    {
      for (int i = 0; i < 4; i++)
        selectedJetsAKF.push_back(selectedJetsKFcorrected.at(i));
    }
    else
    {
      selectedJetsAKF = selectedJetsKFcorrected;
      selectedJetsAKF.push_back(selectedJets[labelsReco[3]]);
      std::sort(selectedJetsAKF.begin(),selectedJetsAKF.end(),HighestPt());
    }
    
    /// Define variables
    reco_W_mass_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).M();
    reco_top_mass_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).M();
    reco_top_mass_alt_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJets[labelsReco[3]]).M();
    reco_top_pt_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).Pt();
    reco_mlb_aKF = (selectedLepton[0] + selectedJets[labelsReco[3]]).M();
    reco_dRLepB_lep_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[3]], selectedLepton[0] );  // deltaR between lepton and leptonic b jet
    reco_dRLepB_had_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
    reco_ttbar_mass_aKF = reco_mlb_aKF + reco_top_mass_aKF;
    
    //if ( fabs(reco_top_mass_aKF - reco_mlb_aKF) < 25 ) continue;
    
    if ( reco_top_mass_aKF < 0. )
      PrintKFDebug(ievt);
    
    
    if ( reco_top_mass_aKF > maxTopMass || reco_top_mass_aKF < minTopMass ) continue;
    if ( reco_top_mass_alt_aKF < maxTopMass-5 ) continue;
    
    
    if (hadronicTopJetsMatched)
    {
      /// Correct match
      if ( ( labelsRecoKF[0] == MCPermutation[0].first || labelsRecoKF[0] == MCPermutation[1].first || labelsRecoKF[0] == MCPermutation[2].first ) && ( labelsRecoKF[1] == MCPermutation[0].first || labelsRecoKF[1] == MCPermutation[1].first || labelsRecoKF[1] == MCPermutation[2].first ) && ( labelsRecoKF[2] == MCPermutation[0].first || labelsRecoKF[2] == MCPermutation[1].first || labelsRecoKF[2] == MCPermutation[2].first ) )  // correct jets for top quark
      {
        if ( selectedJets.size() == 4 ) nofCMnJetsAKF[0]++;
        else if ( selectedJets.size() == 5 ) nofCMnJetsAKF[1]++;
        else if ( selectedJets.size() == 6 ) nofCMnJetsAKF[2]++;
        else nofCMnJetsAKF[3]++;
      }
      else  // wrong permutation
      {
        if ( selectedJets.size() == 4 ) nofWMnJetsAKF[0]++;
        else if ( selectedJets.size() == 5 ) nofWMnJetsAKF[1]++;
        else if ( selectedJets.size() == 6 ) nofWMnJetsAKF[2]++;
        else nofWMnJetsAKF[3]++;
      }
    }
    else  // unmatched
    {
      if ( selectedJets.size() == 4 ) nofUMnJetsAKF[0]++;
      else if ( selectedJets.size() == 5 ) nofUMnJetsAKF[1]++;
      else if ( selectedJets.size() == 6 ) nofUMnJetsAKF[2]++;
      else nofUMnJetsAKF[3]++;
    }
    
  }  // end loop events
  
  
  cout << endl;  /// Stronger selection in this analyser compared to Ntuples ==> endEvent --> nofHardSelected
  cout << "Number of events with exactly 4 jets with pT > 30 GeV: " << nofHardSelected << " (" << 100*((float)nofHardSelected/(float)endEvent) << "%)" << endl;
  cout << "Number of events with clean MET: " << nofMETCleaned << " (" << 100*((float)nofMETCleaned/(float)nofHardSelected) << "%)" << endl;
  cout << "Number of events after min dR cut: " << nofAfterDRmincut << " (" << 100*((float)nofAfterDRmincut/(float)nofMETCleaned) << "%)" << endl;
  if (doKinFit) cout << "Number of clean events accepted by kinFitter: " << nofAcceptedKFit << " (" << 100*((float)nofAcceptedKFit/(float)nofAfterDRmincut) << "%)" << endl;
  
  if (nofHadrMatchedEvents > 0 )
  {
    cout << "Number of matched events: " << setw(8) << right << nofMatchedEvents << endl;
    cout << "Number of events with hadronic top matched (before KF): " << setw(8) << right << nofHadrMatchedEvents << " (" << 100*((float)nofHadrMatchedEvents/(float)nofMETCleaned) << "%)" << endl;
    if (doKinFit) cout << "Number of events with hadronic top matched (after KF):  " << setw(8) << right << nofHadrMatchedEventsAKF << " (" << 100*((float)nofHadrMatchedEventsAKF/(float)nofAcceptedKFit) << "%)" << endl;
    
    
    cout << "Correctly matched reconstructed events:     " << setw(8) << right << nofCorrectlyMatched << endl;
    cout << "                                     hereof " << setw(8) << right  << corrMatchHadrB << " (" << 100*((float)corrMatchHadrB/(float)nofCorrectlyMatched) << "%) have the b jet correctly matched" << endl;
    cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatched << endl;
    if ( nofCorrectlyMatched != 0 || nofNotCorrectlyMatched != 0 )
      cout << "   ===> This means that " << 100*(float)nofCorrectlyMatched / (float)(nofCorrectlyMatched + nofNotCorrectlyMatched) << "% of matched events is correctly matched." << endl;
    
    if (doKinFit)
    {
      cout << "                        " << 100*(float)nofCorrectlyMatched / (float)nofMETCleaned << "% of all events is correctly matched before kinfitter." << endl;
      cout << " --- Kinematic fit --- Before chi2 cut --- " << endl;
      cout << "Correctly matched reconstructed events    : " << setw(8) << right << nofCorrectlyMatchedAKFNoCut << endl;
      cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatchedAKFNoCut << endl;
      if ( nofCorrectlyMatchedAKFNoCut != 0 || nofNotCorrectlyMatchedAKFNoCut != 0 )
        cout << "   ===> This means that " << 100*(float)nofCorrectlyMatchedAKFNoCut / (float)(nofCorrectlyMatchedAKFNoCut + nofNotCorrectlyMatchedAKFNoCut) << "% of matched events is correctly matched after KF." << endl;
      
      cout << "                        " << 100*(float)nofCorrectlyMatchedAKFNoCut / (float)(nofCorrectlyMatchedAKFNoCut+nofNotCorrectlyMatchedAKFNoCut+nofUnmatchedAKFNoCut) << "% of all events accepted by kinfitter is correctly matched." << endl;
      
      cout << " --- Kinematic fit --- After chi2 cut --- " << endl;
      cout << "Correctly matched reconstructed events (after KF): " << setw(8) << right << nofCorrectlyMatchedAKF << endl;
      cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatchedAKF << endl;
      if ( nofCorrectlyMatchedAKF != 0 || nofNotCorrectlyMatchedAKF != 0 )
        cout << "   ===> This means that " << 100*(float)nofCorrectlyMatchedAKF / (float)(nofCorrectlyMatchedAKF + nofNotCorrectlyMatchedAKF) << "% of matched events is correctly matched after KF." << endl;
      
      cout << "                        " << 100*(float)nofCorrectlyMatchedAKF / (float)nofAcceptedKFit << "% of all events accepted by kinfitter is correctly matched." << endl;
    }
    else cout << "                        " << 100*(float)nofCorrectlyMatched / (float)nofMETCleaned << "% of all events is correctly matched." << endl;
  }
  
  cout << endl;
  cout << " 4 jets:   CM " << setw(6) << right << nofCMnJets[0] << " (" << 100*((float)nofCMnJets[0]/(float)(nofCMnJets[0]+nofWMnJets[0]+nofUMnJets[0])) << "%)     WM " << setw(6) << right << nofWMnJets[0] << " (" << 100*((float)nofWMnJets[0]/(float)(nofCMnJets[0]+nofWMnJets[0]+nofUMnJets[0])) << "%)     UM " << setw(6) << right << nofUMnJets[0] << " (" << 100*((float)nofUMnJets[0]/(float)(nofCMnJets[0]+nofWMnJets[0]+nofUMnJets[0])) << "%)" << endl;
  if (! fourJetsOnly)
  {
    cout << " 5 jets:   CM " << setw(6) << right << nofCMnJets[1] << " (" << 100*((float)nofCMnJets[1]/(float)(nofCMnJets[1]+nofWMnJets[1]+nofUMnJets[1])) << "%)     WM " << setw(6) << right << nofWMnJets[1] << " (" << 100*((float)nofWMnJets[1]/(float)(nofCMnJets[1]+nofWMnJets[1]+nofUMnJets[1])) << "%)     UM " << setw(6) << right << nofUMnJets[1] << " (" << 100*((float)nofUMnJets[1]/(float)(nofCMnJets[1]+nofWMnJets[1]+nofUMnJets[1])) << "%)" << endl;
    cout << " 6 jets:   CM " << setw(6) << right << nofCMnJets[2] << " (" << 100*((float)nofCMnJets[2]/(float)(nofCMnJets[2]+nofWMnJets[2]+nofUMnJets[2])) << "%)     WM " << setw(6) << right << nofWMnJets[2] << " (" << 100*((float)nofWMnJets[2]/(float)(nofCMnJets[2]+nofWMnJets[2]+nofUMnJets[2])) << "%)     UM " << setw(6) << right << nofUMnJets[2] << " (" << 100*((float)nofUMnJets[2]/(float)(nofCMnJets[2]+nofWMnJets[2]+nofUMnJets[2])) << "%)" << endl;
    cout << ">6 jets:   CM " << setw(6) << right << nofCMnJets[3] << " (" << 100*((float)nofCMnJets[3]/(float)(nofCMnJets[3]+nofWMnJets[3]+nofUMnJets[3])) << "%)     WM " << setw(6) << right << nofWMnJets[3] << " (" << 100*((float)nofWMnJets[3]/(float)(nofCMnJets[3]+nofWMnJets[3]+nofUMnJets[3])) << "%)     UM " << setw(6) << right << nofUMnJets[3] << " (" << 100*((float)nofUMnJets[3]/(float)(nofCMnJets[3]+nofWMnJets[3]+nofUMnJets[3])) << "%)" << endl;
    
    int total = nofCMnJets[0]+nofWMnJets[0]+nofUMnJets[0] + nofCMnJets[1]+nofWMnJets[1]+nofUMnJets[1] + nofCMnJets[2]+nofWMnJets[2]+nofUMnJets[2] + nofCMnJets[3]+nofWMnJets[3]+nofUMnJets[3];
    cout << "  total:   CM " << setw(6) << right << nofCMnJets[0]+nofCMnJets[1]+nofCMnJets[2]+nofCMnJets[3] << " (" << 100*((float)(nofCMnJets[0]+nofCMnJets[1]+nofCMnJets[2]+nofCMnJets[3])/(float)(total)) << "%)     WM " << setw(6) << right << nofWMnJets[0]+nofWMnJets[1]+nofWMnJets[2]+nofWMnJets[3] << " (" << 100*((float)(nofWMnJets[0]+nofWMnJets[1]+nofWMnJets[2]+nofWMnJets[3])/(float)(total)) << "%)     UM " << setw(6) << right << nofUMnJets[0]+nofUMnJets[1]+nofUMnJets[2]+nofUMnJets[3] << " (" << 100*((float)(nofUMnJets[0]+nofUMnJets[1]+nofUMnJets[2]+nofUMnJets[3])/(float)(total)) << "%)" << endl;
  }
  
  
  cout << endl;
  cout << "After KF" << endl;
  cout << " 4 jets:   CM " << setw(6) << right << nofCMnJetsAKF[0] << " (" << 100*((float)nofCMnJetsAKF[0]/(float)(nofCMnJetsAKF[0]+nofWMnJetsAKF[0]+nofUMnJetsAKF[0])) << "%)     WM " << setw(6) << right << nofWMnJetsAKF[0] << " (" << 100*((float)nofWMnJetsAKF[0]/(float)(nofCMnJetsAKF[0]+nofWMnJetsAKF[0]+nofUMnJetsAKF[0])) << "%)     UM " << setw(6) << right << nofUMnJetsAKF[0] << " (" << 100*((float)nofUMnJetsAKF[0]/(float)(nofCMnJetsAKF[0]+nofWMnJetsAKF[0]+nofUMnJetsAKF[0])) << "%)" << endl;
  if (! fourJetsOnly)
  {
    cout << " 5 jets:   CM " << setw(6) << right << nofCMnJetsAKF[1] << " (" << 100*((float)nofCMnJetsAKF[1]/(float)(nofCMnJetsAKF[1]+nofWMnJetsAKF[1]+nofUMnJetsAKF[1])) << "%)     WM " << setw(6) << right << nofWMnJetsAKF[1] << " (" << 100*((float)nofWMnJetsAKF[1]/(float)(nofCMnJetsAKF[1]+nofWMnJetsAKF[1]+nofUMnJetsAKF[1])) << "%)     UM " << setw(6) << right << nofUMnJetsAKF[1] << " (" << 100*((float)nofUMnJetsAKF[1]/(float)(nofCMnJetsAKF[1]+nofWMnJetsAKF[1]+nofUMnJetsAKF[1])) << "%)" << endl;
    cout << " 6 jets:   CM " << setw(6) << right << nofCMnJetsAKF[2] << " (" << 100*((float)nofCMnJetsAKF[2]/(float)(nofCMnJetsAKF[2]+nofWMnJetsAKF[2]+nofUMnJetsAKF[2])) << "%)     WM " << setw(6) << right << nofWMnJetsAKF[2] << " (" << 100*((float)nofWMnJetsAKF[2]/(float)(nofCMnJetsAKF[2]+nofWMnJetsAKF[2]+nofUMnJetsAKF[2])) << "%)     UM " << setw(6) << right << nofUMnJetsAKF[2] << " (" << 100*((float)nofUMnJetsAKF[2]/(float)(nofCMnJetsAKF[2]+nofWMnJetsAKF[2]+nofUMnJetsAKF[2])) << "%)" << endl;
    cout << ">6 jets:   CM " << setw(6) << right << nofCMnJetsAKF[3] << " (" << 100*((float)nofCMnJetsAKF[3]/(float)(nofCMnJetsAKF[3]+nofWMnJetsAKF[3]+nofUMnJetsAKF[3])) << "%)     WM " << setw(6) << right << nofWMnJetsAKF[3] << " (" << 100*((float)nofWMnJetsAKF[3]/(float)(nofCMnJetsAKF[3]+nofWMnJetsAKF[3]+nofUMnJetsAKF[3])) << "%)     UM " << setw(6) << right << nofUMnJetsAKF[3] << " (" << 100*((float)nofUMnJetsAKF[3]/(float)(nofCMnJetsAKF[3]+nofWMnJetsAKF[3]+nofUMnJetsAKF[3])) << "%)" << endl;
    
    int total = nofCMnJetsAKF[0]+nofWMnJetsAKF[0]+nofUMnJetsAKF[0] + nofCMnJetsAKF[1]+nofWMnJetsAKF[1]+nofUMnJetsAKF[1] + nofCMnJetsAKF[2]+nofWMnJetsAKF[2]+nofUMnJetsAKF[2] + nofCMnJetsAKF[3]+nofWMnJetsAKF[3]+nofUMnJetsAKF[3];
    cout << "  total:   CM " << setw(6) << right << nofCMnJetsAKF[0]+nofCMnJetsAKF[1]+nofCMnJetsAKF[2]+nofCMnJetsAKF[3] << " (" << 100*((float)(nofCMnJetsAKF[0]+nofCMnJetsAKF[1]+nofCMnJetsAKF[2]+nofCMnJetsAKF[3])/(float)(total)) << "%)     WM " << setw(6) << right << nofWMnJetsAKF[0]+nofWMnJetsAKF[1]+nofWMnJetsAKF[2]+nofWMnJetsAKF[3] << " (" << 100*((float)(nofWMnJetsAKF[0]+nofWMnJetsAKF[1]+nofWMnJetsAKF[2]+nofWMnJetsAKF[3])/(float)(total)) << "%)     UM " << setw(6) << right << nofUMnJetsAKF[0]+nofUMnJetsAKF[1]+nofUMnJetsAKF[2]+nofUMnJetsAKF[3] << " (" << 100*((float)(nofUMnJetsAKF[0]+nofUMnJetsAKF[1]+nofUMnJetsAKF[2]+nofUMnJetsAKF[3])/(float)(total)) << "%)" << endl;
  }
  
  
  origNtuple->Close();
  
  
  if (test)
  {
    cout << "Exiting because of test..." << endl;
    exit(1);
  }
  
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  if (makePlots)
  {
    string rootFileName = "NtuplePlots_"+systStr+".root";
    
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

void GetMetaData(TTree* tree)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  
  tree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
  tree->SetBranchAddress("nEventsSel", &nEventsSel, &b_nEventsSel);
  tree->SetBranchAddress("cutFlow", cutFlow, &b_cutFlow);
  tree->SetBranchAddress("cutFlow2", cutFlow2, &b_cutFlow2);
  tree->SetBranchAddress("cutFlowWeighted", cutFlowWeighted, &b_cutFlowWeighted);     // temporarily!
  tree->SetBranchAddress("cutFlow2Weighted", cutFlow2Weighted, &b_cutFlow2Weighted);  // temporarily!
  tree->SetBranchAddress("appliedJER", &appliedJER, &b_appliedJER);
  tree->SetBranchAddress("appliedJES", &appliedJES, &b_appliedJES);
  tree->SetBranchAddress("appliedPU", &appliedPU, &b_appliedPU);
  
//   tree->SetBranchAddress("nofEventsWithGenTop", &nofEventsWithGenTop, &b_nofEventsWithGenTop);
//   tree->SetBranchAddress("nofEventsWithGenTopWithStatus22or62", &nofEventsWithGenTopWithStatus22or62, &b_nofEventsWithGenTopWithStatus22or62);
//   tree->SetBranchAddress("nofEventsWithGenAntiTop", &nofEventsWithGenAntiTop, &b_nofEventsWithGenAntiTop);
//   tree->SetBranchAddress("nofEventsWithGenAntiTopWithStatus22or62", &nofEventsWithGenAntiTopWithStatus22or62, &b_nofEventsWithGenAntiTopWithStatus22or62);
//   tree->SetBranchAddress("nofTTEventsWithoutBothGenTops", &nofTTEventsWithoutBothGenTops, &b_nofTTEventsWithoutBothGenTops);
  tree->SetBranchAddress("nofTTEventsWithoutAGenTop", &nofTTEventsWithoutAGenTop, &b_nofTTEventsWithoutAGenTop);
//   tree->SetBranchAddress("nofTTEventsWithoutGenTop", &nofTTEventsWithoutGenTop, &b_nofTTEventsWithoutGenTop);
//   tree->SetBranchAddress("nofTTEventsWithoutGenAntiTop", &nofTTEventsWithoutGenAntiTop, &b_nofTTEventsWithoutGenAntiTop);
//   tree->SetBranchAddress("nofTTEventsWithoutBothGenTopsWithStatus22", &nofTTEventsWithoutBothGenTopsWithStatus22, &b_nofTTEventsWithoutBothGenTopsWithStatus22);
//   tree->SetBranchAddress("nofTTEventsWithoutGenTopWithStatus22", &nofTTEventsWithoutGenTopWithStatus22, &b_nofTTEventsWithoutGenTopWithStatus22);
//   tree->SetBranchAddress("nofTTEventsWithoutGenAntiTopWithStatus22", &nofTTEventsWithoutGenAntiTopWithStatus22, &b_nofTTEventsWithoutGenAntiTopWithStatus22);
//   tree->SetBranchAddress("nofTTEventsWithoutBothGenTopsWithStatus62", &nofTTEventsWithoutBothGenTopsWithStatus62, &b_nofTTEventsWithoutBothGenTopsWithStatus62);
//   tree->SetBranchAddress("nofTTEventsWithoutGenTopWithStatus62", &nofTTEventsWithoutGenTopWithStatus62, &b_nofTTEventsWithoutGenTopWithStatus62);
//   tree->SetBranchAddress("nofTTEventsWithoutGenAntiTopWithStatus62", &nofTTEventsWithoutGenAntiTopWithStatus62, &b_nofTTEventsWithoutGenAntiTopWithStatus62);
  tree->SetBranchAddress("sumWeight1001", &sumWeight1001, &b_sumWeight1001);
  tree->SetBranchAddress("sumWeight1002", &sumWeight1002, &b_sumWeight1002);
  tree->SetBranchAddress("sumWeight1003", &sumWeight1003, &b_sumWeight1003);
  tree->SetBranchAddress("sumWeight1004", &sumWeight1004, &b_sumWeight1004);
  tree->SetBranchAddress("sumWeight1005", &sumWeight1005, &b_sumWeight1005);
  tree->SetBranchAddress("sumWeight1007", &sumWeight1007, &b_sumWeight1007);
  tree->SetBranchAddress("sumWeight1009", &sumWeight1009, &b_sumWeight1009);
}

void InitTree(TTree* tree)
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
//   tree->SetBranchAddress("filterHBHENoise", &filterHBHENoise, &b_filterHBHENoise);
//   tree->SetBranchAddress("filterHBHEIso", &filterHBHEIso, &b_filterHBHEIso);
//   tree->SetBranchAddress("filterCSCTightHalo", &filterCSCTightHalo, &b_filterCSCTightHalo);
//   tree->SetBranchAddress("filterEcalDeadCell", &filterEcalDeadCell, &b_filterEcalDeadCell);
//   tree->SetBranchAddress("filterEEBadSc", &filterEEBadSc, &b_filterEEBadSc);
//   tree->SetBranchAddress("filterBadChCand", &filterBadChCand, &b_filterBadChCand);
//   tree->SetBranchAddress("filterBadMuon", &filterBadMuon, &b_filterBadMuon);
//   tree->SetBranchAddress("passedMETFilter", &passedMETFilter, &b_passedMETFilter);
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
//   tree->SetBranchAddress("hasGenTopWithStatus22", &hasGenTopWithStatus22, &b_hasGenTopWithStatus22);
//   tree->SetBranchAddress("hasGenTopWithStatus62", &hasGenTopWithStatus62, &b_hasGenTopWithStatus62);
  tree->SetBranchAddress("hasGenAntiTop", &hasGenAntiTop, &b_hasGenAntiTop);
//   tree->SetBranchAddress("hasGenAntiTopWithStatus22", &hasGenAntiTopWithStatus22, &b_hasGenAntiTopWithStatus22);
//   tree->SetBranchAddress("hasGenAntiTopWithStatus62", &hasGenAntiTopWithStatus62, &b_hasGenAntiTopWithStatus62);
  tree->SetBranchAddress("weight1001", &weight1001, &b_weight1001);
  tree->SetBranchAddress("weight1002", &weight1002, &b_weight1002);
  tree->SetBranchAddress("weight1003", &weight1003, &b_weight1003);
  tree->SetBranchAddress("weight1004", &weight1004, &b_weight1004);
  tree->SetBranchAddress("weight1005", &weight1005, &b_weight1005);
  tree->SetBranchAddress("weight1007", &weight1007, &b_weight1007);
  tree->SetBranchAddress("weight1009", &weight1009, &b_weight1009);
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

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
  
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
      //if (isTTbar)
      //{
        if ( ( muPlusFromTop && mc_mother[partonId[j]] == -24 /*&& mc_granny[partonId[j]] == -pdgID_top*/ )
            || ( muMinusFromTop && mc_mother[partonId[j]] == 24 /*&& mc_granny[partonId[j]] == pdgID_top*/ ) )  // if mu+, check if mother of particle is W- and granny tbar --> then it is a quark from W- decay
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
      //}
      //       else   /// normally only ST tW should give results
      //       {
      //         if ( ( foundMuPlus && mc_mother[partonId[j]] == -24 && mc_granny[partonId[j]] == -pdgID_top )
      //           || ( foundMuMinus && mc_mother[partonId[j]] == 24 && mc_granny[partonId[j]] == pdgID_top ) )  // if mu+, check if mother of particle is W- and granny tbar --> then it is a quark from W- decay
      //         {
      //           if (MCPermutation[0].first == 9999)
      //           {
      //             MCPermutation[0] = JetPartonPair[i];
      //           }
      //           else if (MCPermutation[1].first == 9999)
      //           {
      //             MCPermutation[1] = JetPartonPair[i];
      //           }
      //           else
      //           {
      //             cerr << "Found a third jet coming from a W boson which comes from a top quark..." << endl;
      //             cerr << " -- muMinusFromMtop: " << muMinusFromTop << " muPlusFromMtop: " << muPlusFromTop << endl;
      //             cerr << " -- pdgId: " << mc_pdgId[partonId[j]] << " mother: " << mc_mother[partonId[j]] << " granny: " << mc_granny[partonId[j]] << " Pt: " << mc_pt[partonId[j]] << endl;
      //   //          cerr << " -- ievt: " << ievt << endl;
      //             exit(1);
      //           }
      //         }
      //       }
    }
    if ( fabs(mc_pdgId[partonId[j]]) == 5 )
    {
      if ( ( (muPlusFromTop || foundMuPlus) && mc_mother[partonId[j]] == -pdgID_top )
          || ( (muMinusFromTop || foundMuMinus) && mc_mother[partonId[j]] == pdgID_top ) )  // if mu+ (top decay leptonic) and mother is antitop ---> hadronic b
      {
        if (verbose > 3)
          cout << "b jet:     " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;
        
        MCPermutation[2] = JetPartonPair[i];
      }
      else if ( ( (muPlusFromTop || foundMuPlus) && mc_mother[partonId[j]] == pdgID_top )
               || ( (muMinusFromTop || foundMuMinus) && mc_mother[partonId[j]] == -pdgID_top ) )  // if mu+ (top decay leptonic) and mother is top ---> leptonic b
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
    cutFlowWeighted[i] = 0;
    cutFlow2Weighted[i] = 0;
  }
  appliedJER = 999;
  appliedJES = 999;
  appliedPU = 999;
//   nofEventsWithGenTop = 0;
//   nofEventsWithGenTopWithStatus22or62 = 0;
//   nofEventsWithGenAntiTop = 0;
//   nofEventsWithGenAntiTopWithStatus22or62 = 0;
//   nofTTEventsWithoutBothGenTops = 0;
  nofTTEventsWithoutAGenTop = 0;
//   nofTTEventsWithoutGenTop = 0;
//   nofTTEventsWithoutGenAntiTop = 0;
//   nofTTEventsWithoutBothGenTopsWithStatus22 = 0;
//   nofTTEventsWithoutGenTopWithStatus22 = 0;
//   nofTTEventsWithoutGenAntiTopWithStatus22 = 0;
//   nofTTEventsWithoutBothGenTopsWithStatus62 = 0;
//   nofTTEventsWithoutGenTopWithStatus62 = 0;
//   nofTTEventsWithoutGenAntiTopWithStatus62 = 0;
  sumWeight1001 = 0;
  sumWeight1002 = 0;
  sumWeight1003 = 0;
  sumWeight1004 = 0;
  sumWeight1005 = 0;
  sumWeight1007 = 0;
  sumWeight1009 = 0;
  
  strSyst = "";
  nEventsDataSet = 0;
  xSection = 1.;
  eqLumi = 1.;
  lumiWeight = 1.;
  
  nofHardSelected = 0;
  nofMETCleaned = 0;
  nofAfterDRmincut = 0;
  nofMatchedEvents = 0;
  nofHadrMatchedEvents = 0;
  nofHadrMatchedEventsAKF = 0;
  nofCorrectlyMatched = 0;
  nofNotCorrectlyMatched = 0;
  nofUnmatched = 0;
  nofCorrectlyMatchedAKF = 0;
  nofNotCorrectlyMatchedAKF = 0;
  nofUnmatchedAKF = 0;
  nofCorrectlyMatchedAKFNoCut = 0;
  nofNotCorrectlyMatchedAKFNoCut = 0;
  nofUnmatchedAKFNoCut = 0;
  nofAcceptedKFit = 0;
  nofAcceptedKFitWeighted = 0.;
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
//   hasGenTopWithStatus22 = false;
//   hasGenTopWithStatus62 = false;
  hasGenAntiTop = false;
//   hasGenAntiTopWithStatus22 = false;
//   hasGenAntiTopWithStatus62 = false;
  
  weight1001 = 1.;
  weight1002 = 1.;
  weight1003 = 1.;
  weight1004 = 1.;
  weight1005 = 1.;
  weight1007 = 1.;
  weight1009 = 1.;
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

void ClearTLVs()
{
  muon.Clear();
  jet.Clear();
  mcpart.Clear();
  WCandidate.Clear();
  neutrino.Clear();
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
  foundMuPlus = false;
  foundMuMinus = false;
  muPlusFromTop = false;
  muMinusFromTop = false;
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
  bJetId.clear();
  bdiscrTop = -99.;
  bdiscrTop2 = -99.;
  tempbdiscr = -99.;
  labelB1 = -9999;
  labelB2 = -9999;
  for (int i = 0; i < 4; i++)
  {
    labelsReco[i] = -9999;
    labelsRecoKF[i] = -9999;
  }
  massHadTopQ = 0.01;
  massLepTopQ = 0.01;
  catSuffix = "";
  isCM = false;
  isWM = false;
  isUM = false;
  doneKinFit = false;
  kFitVerbosity = false;
  kFitChi2_01 = 99.;
  kFitChi2_02 = 99.;
  kFitChi2_12 = 99.;
  kFitChi2_min = 99.;
  min01 = false;
  min02 = false;
  min12 = false;
  
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
  reco_W_mass_aKF = -1.;
  reco_top_mass_aKF = -1.;
  reco_top_mass_alt_aKF = -1.;
  reco_top_pt_aKF = -1.;
  reco_mlb_aKF = -1.;
  reco_dRLepB_lep_aKF = -1.;
  reco_dRLepB_had_aKF = -1.;
  reco_ttbar_mass_aKF = -1.;
}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearVars();
}

long GetNEvents(TTree* fChain, string var)
{
  GetMetaData(fChain);
  long varNew = 0;
  for (unsigned int iEntry = 0; iEntry < fChain->GetEntries(); iEntry++)
  {
    fChain->GetEntry(iEntry);
    varNew += (fChain->FindLeaf(var.c_str()))->GetValueLong64();
  }
  
  return varNew;
}

long GetNEvents(TTree* fChain, string var, unsigned int index)
{
  GetMetaData(fChain);
  long varNew = 0;
  for (unsigned int iEntry = 0; iEntry < fChain->GetEntries(); iEntry++)
  {
    fChain->GetEntry(iEntry);
    varNew += (fChain->FindLeaf(var.c_str()))->GetValueLong64(index);
  }
  
  return varNew;
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

