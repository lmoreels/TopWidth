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
bool unblind = false;
bool doGenOnly = false;
bool makePlots = false;
bool makeControlPlots = false;
bool makeLikelihoodPlots = false;
bool calculateResolutionFunctions = false;
bool calculateAverageMass = false;
bool calculateFractions = false;
bool makeTGraphs = false;
bool useTTTemplates = false;
bool calculateLikelihood = true;
bool doPseudoExps = false;
bool doKinFit = true;
bool applyKinFitCut = true;
double kinFitCutValue = 5.;

bool doMETCleaning = true;
bool applyLeptonSF = true;
bool applyPU = true;
bool applyBTagSF = true;

bool runTTbar = true;
bool runSTtW = true;
bool runSTt = true;
bool runOther = true;


bool rewHadTopOnly = false;
bool applyWidthSF = false;
double scaleWidth = 0.6;


bool runListWidths = true;
double listWidths[] = {0.2, 0.4, 0.5, 0.6, 0.8, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.};
int nWidths = sizeof(listWidths)/sizeof(listWidths[0]);
double thisWidth;

bool runSystematics = false;
bool runRateSystematics = false;
bool runSampleSystematics = false;
int nSystematics;
string thisSystematic;

string listRateSyst[] = {"nominal", "leptonIdSFup", "leptonIdSFdown", "leptonIsoSFup", "leptonIsoSFdown", "leptonTrigSFup", "leptonTrigSFdown", "leptonTrkSFup", "leptonTrkSFdown", "puSFup", "puSFdown", "btagSFup", "btagSFdown", "topPtReweighting", "lumiup", "lumidown", "renFac1002", "renFac1003", "renFac1004", "renFac1005", "renFac1007", "renFac1009"};
int nRateSystematics = sizeof(listRateSyst)/sizeof(listRateSyst[0]);
string listSampleSyst[] = {"tuneup", "tunedown", "isrup", "isrdown", "fsrup", "fsrdown", "hdampup", "hdampdown", "mpiERD", "qcdERD", "gluonMove", "gluonMoveERD", "mass169p5", "mass175p5", "herwig"};
string dataSetNameSyst[] = {"TT_tune_up", "TT_tune_down", "TT_isr_up", "TT_isr_down", "TT_fsr_up", "TT_fsr_down", "TT_hdamp_up", "TT_hdamp_down", "TT_erdOn", "TT_QCD_erdOn", "TT_gluon_move", "TT_gluon_move_erdOn", "TT_mass169p5", "TT_mass175p5", "TT_herwigpp"};
int nSampleSystematics = sizeof(listSampleSyst)/sizeof(listSampleSyst[0]);
const char *xmlSyst = "config/topWidth_extra.xml";


TFile *fileWidths;


string systStr = "nominal";
pair<string,string> whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos )
  {
    return pair<string,string>("170712","170522");
  }
  else if ( syst.find("JECup") != std::string::npos ) return pair<string,string>("170602","170522");
  else if ( syst.find("JECdown") != std::string::npos ) return pair<string,string>("170606","170522");
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    return pair<string,string>("170712","170522");
  }
}
pair<string,string> ntupleDate = whichDate(systStr);
string ntupleSystDate = "170803";
int verbose = 2;

string pathNtuples = "";
string pathNtuplesMC = "";
string pathNtuplesData = "";
string pathNtuplesSyst = "";
string pathOutput = "";
string outputDirLL = "LikelihoodTemplates/";
string inputDirLL = "";
string inputDateLL = "170829_1010/";  // use relativeSF instead of lumiWeight
//string inputDateLL = "170814_1625/";  // all SFs, all samples; widthSF had&lep; redtopmass [0.6, 1.4]; other widths
//string inputDateLL = "170802_1317/";  // all SFs, ttbar-only; widthSF had&lep; redtopmass [0.6, 1.4]
//string inputDateLL = "170802_1222/";  // all SFs, all samples; widthSF had&lep; redtopmass [0.6, 1.4]
//string inputDateLL = "170801_1451/";  // all SFs, ttbar-only; mlb cut 160; widthSF had-only; redtopmass [0.7, 1.3]
//string inputDateLL = "170801_1252/";  // all SFs, ttbar-only; mlb cut 160; widthSF had & lep; redtopmass [0.7, 1.3]
//string inputDateLL = "170731_1824/";  // all SFs, ttbar-only; no smoothing; mlb cut 160; widthSF had-only
//string inputDateLL = "170731_1823/";  // all SFs, ttbar-only; no smoothing; mlb cut 160; widthSF had & lep
//string inputDateLL = "170728_1418/";  // w/o W+1/2jets, all SFs, ttbar-only; no smoothing; widthSF had & lep
//string inputDateLL = "170728_1324/";  // w/o W+1/2jets, all SFs, ttbar-only; no smoothing; widthSF had-only
string whichTemplates()
{
  if      (rewHadTopOnly) return "170728_1324/";
  else return "170802_1317/";
}

bool isData = false;
bool isTTbar = false;
bool isST = false;
bool isSTtW = false;
bool isSTt = false;
bool isOther = false;
bool isHerwig = false;

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
double nofCMl = 0., nofWMl = 0., nofNMl = 0.;
double nofCM_weighted = 0., nofWM_weighted = 0., nofNM_weighted = 0.;
double nofCMout_weighted = 0., nofWMout_weighted = 0., nofNMout_weighted = 0.;

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
std::array<double, 14> aveTopMass = {171.833, 169.809, 167.636, 197.975, 197.718, 198.582, 182.317, 252.174, 249.964, 229.383, 227.814, 184.794, 185.096, 185.046};  // with SFs
//std::array<double, 14> aveTopMass = {171.826, 169.746, 167.511, 197.053, 196.687, 197.911, 181.895, 249.468, 247.437, 227.530, 226.099, 184.794, 184.594, 184.624};  // Res 170608 Single Gaus
//std::array<double, 14> aveTopMass = {171.826, 169.746, 167.556, 197.087, 196.662, 198.143, 182.150, 249.229, 246.893, 226.933, 225.681, 185.024, 184.880, 184.902};  // no DYJets, no WJets // Res 170608

/// # events after kin fitter
//  KF chi2 < 5
int nEventsAKF[] = {100949, 547108, 8927, 8847, 7051, 4283, 3, 16, 55, 276, 1, 18, 123, 1415};

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlotCP;
map<string,TH1D*> histo1DLike;
map<string,TGraph*> graph;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;

vector < Dataset* > datasets, datasetsMSP, datasetsTemp, datasetsSyst;
vector<string> dataSetNames;
vector<int> includeDataSets;

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
void InitMSPlots();
void InitHisto1D();
void InitHisto2D();
void InitHisto1DMatch();
void InitHisto2DMatch();
void InitLikelihoodPlots();
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
void FillMSPlots(int d, bool doneKinFit);
void FillLikelihoodPlots();
void WriteLikelihoodPlots();
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
Double_t        weight1001;
Double_t        weight1002;
Double_t        weight1003;
Double_t        weight1004;
Double_t        weight1005;
Double_t        weight1007;
Double_t        weight1009;
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
Long64_t        nofTTEventsWithoutAGenTop;
Long64_t        nofTTEventsWithoutGenTop;
Long64_t        nofTTEventsWithoutGenAntiTop;
Long64_t        nofTTEventsWithoutBothGenTopsWithStatus22;
Long64_t        nofTTEventsWithoutGenTopWithStatus22;
Long64_t        nofTTEventsWithoutGenAntiTopWithStatus22;
Long64_t        nofTTEventsWithoutBothGenTopsWithStatus62;
Long64_t        nofTTEventsWithoutGenTopWithStatus62;
Long64_t        nofTTEventsWithoutGenAntiTopWithStatus62;
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
TBranch        *b_nofTTEventsWithoutAGenTop;   //!
TBranch        *b_nofTTEventsWithoutGenTop;   //!
TBranch        *b_nofTTEventsWithoutGenAntiTop;   //!
TBranch        *b_nofTTEventsWithoutBothGenTopsWithStatus22;   //!
TBranch        *b_nofTTEventsWithoutGenTopWithStatus22;   //!
TBranch        *b_nofTTEventsWithoutGenAntiTopWithStatus22;   //!
TBranch        *b_nofTTEventsWithoutBothGenTopsWithStatus62;   //!
TBranch        *b_nofTTEventsWithoutGenTopWithStatus62;   //!
TBranch        *b_nofTTEventsWithoutGenAntiTopWithStatus62;   //!
TBranch        *b_sumWeight1001;   //!
TBranch        *b_sumWeight1002;   //!
TBranch        *b_sumWeight1003;   //!
TBranch        *b_sumWeight1004;   //!
TBranch        *b_sumWeight1005;   //!
TBranch        *b_sumWeight1007;   //!
TBranch        *b_sumWeight1009;   //!


long nEventsDataSet;
double xSection;
double lumiWeight, scaleFactor, widthSF, relativeSF, eqLumi_TT;
double thisLeptonSF, thisLeptonIdSF, thisLeptonIsoSF, thisLeptonTrigSF;
double renFacSumNom, renFacSum1002, renFacSum1003, renFacSum1004, renFacSum1005, renFacSum1007, renFacSum1009;
double topPtRewSF, topPtSF, antiTopPtSF;
bool foundTop22, foundAntiTop22;
bool foundTop62, foundAntiTop62;
bool foundLastCopyTop, foundLastCopyAntitop;
vector<unsigned int> bJetId;
double bdiscrTop, bdiscrTop2, tempbdiscr;
int labelB1, labelB2;
int labelsReco[4];
double massHadTopQ, massLepTopQ;

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
bool muonmatched = false;
bool foundMuPlus = false, foundMuMinus = false;
bool muPlusFromTop = false, muMinusFromTop = false;
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


/// Likelihood
Double_t aveTopMassLL = aveTopMass[2];
Double_t maxRedTopMass = 0., minRedTopMass = 9999.;
Double_t minCutRedTopMass = 0.6, maxCutRedTopMass = 1.4;
int nWidthsLike = 0;
vector<double> widthsLike;
vector<double> loglike_per_evt;
double redTopMassArray[] = {0.65, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1., 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.4};
int nLikeMasses = sizeof(redTopMassArray)/sizeof(redTopMassArray[0]) - 1;
Likelihood *like;


ofstream txtDebugTopMass, txtDebugPUSF;

/// Pseudo experiments
const int nPseudoExps = 400;
int nPsExps = nPseudoExps;
TRandom3 random3;
double toyValues[nPseudoExps]; // = random3.Uniform(0,1);
double toyMax;
int nEvtsInPseudoExp[nPseudoExps][15] = {0}, nDataEvts;
double nEvtsInPseudoExpW[nPseudoExps][15] = {0.};


/// Variables
double M3, Ht, min_Mlb, dRLepB;
double M3_aKF, Ht_aKF;
double reco_W_mass_bKF, reco_top_mass_bKF, reco_top_pt_bKF, reco_mlb_bKF, reco_dRLepB_lep_bKF, reco_dRLepB_had_bKF, reco_ttbar_mass_bKF, redTopMass_bKF;
double reco_W_mass_aKF, reco_top_mass_aKF, reco_top_pt_aKF, reco_mlb_aKF, reco_dRLepB_lep_aKF, reco_dRLepB_had_aKF, reco_ttbar_mass_aKF, redTopMass;

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
  cout << "* The following scale factors are applied:  *" << endl;
  if (applyLeptonSF) cout << "*   - Lepton scale factors: nominal         *" << endl;
  if (applyPU)       cout << "*   - Pile up: nominal                      *" << endl;
  if (applyBTagSF)   cout << "*   - B tag scale factors: nominal          *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  string channel;
  if ( argc == 1 ) channel = "mu";
  //else if ( CharSearch(argv[1], "mu") || CharSearch(argv[1], "Mu") || CharSearch(argv[1], "MU") || CharSearch(argv[1], "m") ) channel = "mu";
  //else if ( CharSearch(argv[1], "el") || CharSearch(argv[1], "El") || CharSearch(argv[1], "EL") || CharSearch(argv[1], "e") ) channel = "el";
  //else if ( (argv[1]).find("all") != std::string::npos || (argv[1]).find("All") != std::string::npos || (argv[1]).find("ALL") != std::string::npos ) channel = "all";
  
  if (calculateAverageMass)
  {
    makeTGraphs = false;
    calculateLikelihood = false;
    doPseudoExps = false;
    makePlots = false;
    doGenOnly = false;
    runListWidths = false;
    runSystematics = false;
  }
  if (calculateResolutionFunctions)
  {
    testTTbarOnly = true;
    doGenOnly = true;
    makePlots = false;
    calculateAverageMass = false;
    calculateLikelihood = false;
    makeTGraphs = false;
    doPseudoExps = false;
    doKinFit = false;
    runListWidths = false;
    runSystematics = false;
  }
  if (test) makePlots = false;
  if (testHistos)
  {
    makePlots = true;
    makeControlPlots = true;
    doGenOnly = false;
    makeTGraphs = false;
    runListWidths = false;
    runSystematics = false;
  }
  if (doGenOnly)
  {
    doPseudoExps = false;
    runSystematics = false;
  }
  if (calculateFractions)
  {
    calculateLikelihood = false;
    doPseudoExps = false;
    runListWidths = false;
    runSystematics = false;
  }
  if (makeTGraphs)
  {
    runListWidths = false;
    calculateLikelihood = false;
  }
  if (calculateLikelihood) makeTGraphs = false;
  else
  {
    makeLikelihoodPlots = false;
    doPseudoExps = false;
  }
  if (doPseudoExps)
  {
    makePlots = false;
    runListWidths = false;
  }
  if (runListWidths)
  {
    makePlots = false;
    unblind = false;
  }
  if (runSystematics)
  {
    runListWidths = false;
    unblind = false;
    makeTGraphs = false;
    calculateLikelihood = true;
    doPseudoExps = false;
  }
  else
  {
    runRateSystematics = false;
    runSampleSystematics = false;
  }
  if (! runRateSystematics && ! runSampleSystematics) runSystematics = false;
  
  if (! makePlots) makeControlPlots = false;
  
  pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (makePlots || makeLikelihoodPlots)
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
  pathNtuplesSyst = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleSystDate+"/";
  if (calculateAverageMass) cout << "Calculating average mass values..." << endl;
  if (calculateLikelihood) cout << "Calculating -loglikelihood values using templates from " << inputDateLL << endl;
  if (doPseudoExps)       cout << "              for pseudo experiments" << endl;
  if (testHistos) cout << "Testing histogram consistency..." << endl;
  if (doGenOnly) cout << "Running only matching..." << endl;
  if (applyWidthSF) cout << "TTbar sample width will be scaled by a factor " << scaleWidth << endl;
  else scaleWidth = 1.;
  
  bool runAll = false;
  if (runTTbar && runSTtW && runSTt && runOther) runAll = true;
  
  if (testTTbarOnly || (runTTbar && ! runSTtW && ! runSTt && ! runOther) )
    cout << "Only running ttbar events..." << endl;
  else if (runAll) cout << "Running all datasets" << endl;
  else if (runTTbar) cout << "Running ttbar events..." << endl;
  else if (runSTtW || runSTt)
  {
    if (runSTtW && runSTt && ! runTTbar && ! runOther) cout << "Only running ST events..." << endl;
    else if (runSTtW && ! runSTt && ! runTTbar && ! runOther) cout << "Only running ST tW events..." << endl;
    else if (! runSTtW && runSTt && ! runTTbar && ! runOther) cout << "Only running ST t events..." << endl;
    else if (runSTtW && runSTt) cout << "Running ST events..." << endl;
    else if (runSTtW) cout << "Running ST tW events..." << endl;
    else if (runSTt) cout << "Running ST t events..." << endl;
  }
  else if (runOther)
  {
    if (! runTTbar && ! runSTtW && ! runSTt) cout << "Only running DY+jets and W+jets events... " << endl;
    else cout << "Running DY+jets and W+jets events... " << endl;
  }
  
  
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
    
    if ( dataSetName.find("TT") != std::string::npos )
      dataSetNames.push_back("TT_nominal");   // Temporarily ! Change to "TT" when making new (final?) templates
    else
      dataSetNames.push_back(dataSetName);
    
    
    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA") != std::string::npos )
      includeDataSets.push_back(0);
//      Luminosity = datasets[d]->EquivalentLumi();
    else if ( dataSetName.find("TT") != std::string::npos )
      includeDataSets.push_back(1);
    else if ( dataSetName.find("ST") != std::string::npos )
      includeDataSets.push_back(1);
    else
      includeDataSets.push_back(1);
    
    if ( dataSetName.find("QCD") != std::string::npos )
    {
      datasets[d]->SetColor(kYellow);
      datasetsMSP[d]->SetColor(kYellow);
    }
    if ( dataSetName.find("TT") != std::string::npos )
    {
      datasets[d]->SetTitle("t#bar{t}");
      datasets[d]->SetColor(kRed+1);
      datasetsMSP[dTT+2]->SetName("TT_NM");
      datasetsMSP[dTT+2]->SetTitle("t#bar{t} unmatched");
      datasetsMSP[dTT+2]->SetColor(kRed-10);
      
      dTT = d;
    }
    //if ( dataSetName.find("TTbarJets_Other") != std::string::npos ) datasets[d]->SetColor(kRed-7);
    if ( dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      //datasets[d]->SetColor(kGreen-3);
      datasets[d]->SetColor(kBlue-2);
      datasetsMSP[d]->SetTitle("W#rightarrowl#nu");
      //datasetsMSP[d]->SetColor(kGreen-3);
      datasetsMSP[d]->SetColor(kBlue-2);
    }
    if ( dataSetName.find("ZJets") != std::string::npos || dataSetName.find("DY") != std::string::npos )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{#font[122]{\55}}");
      //datasets[d]->SetColor(kAzure-2);
      //datasets[d]->SetColor(kMagenta);
      datasets[d]->SetColor(kAzure+6);
      datasetsMSP[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{#font[122]{\55}}");
      //datasetsMSP[d]->SetColor(kAzure-2);
      //datasetsMSP[d]->SetColor(kMagenta);
      datasetsMSP[d]->SetColor(kAzure+6);
    }
    if ( dataSetName.find("ST") != std::string::npos || dataSetName.find("SingleTop") != std::string::npos )
    {
      datasets[d]->SetTitle("ST");
      //datasets[d]->SetColor(kBlue-2);
      datasets[d]->SetColor(kOrange-4);  // 595, 615, 800
      datasetsMSP[d]->SetTitle("ST");
      //datasetsMSP[d]->SetColor(kBlue-2);
      datasetsMSP[d]->SetColor(kOrange-4);
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
  datasetsTemp[dTT]->SetName("TT_WM");
  datasetsTemp[dTT]->SetTitle("t#bar{t} wrong");
  datasetsTemp[dTT]->SetColor(kRed-9);
  datasetsMSP.insert(datasetsMSP.begin()+dTT, 1, datasetsTemp[dTT]);
  datasetsTemp.clear();
  treeLoader.LoadDatasets(datasetsTemp, xmlfile);
  datasetsTemp[dTT]->SetName("TT_CM");
  datasetsTemp[dTT]->SetTitle("t#bar{t} correct");
  datasetsTemp[dTT]->SetColor(kRed-7);
  datasetsMSP.insert(datasetsMSP.begin()+dTT, 1, datasetsTemp[dTT]);
  
//   cout << "Content datasets  " << endl;
//   for (int d = 0; d < datasets.size(); d++)
//     cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
//   cout << "Content datasetsMSP  " << endl;
//   for (int d = 0; d < datasetsMSP.size(); d++)
//     cout << "   Dataset " << d << ": " << datasetsMSP[d]->Name() << " / title : " << datasetsMSP[d]->Title() << endl;
  
  
  /// Load systematic samples
  if (runSampleSystematics)
  {
    datasetsSyst.clear();
    treeLoader.LoadDatasets(datasetsSyst, xmlSyst);
    cout << "Found " << datasetsSyst.size() << " extra systematics samples..." << endl;
  }
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  EventReweighting *rew = new EventReweighting(false);  // no correction for number of events
  ResolutionFunctions* rf = new ResolutionFunctions(calculateResolutionFunctions, true);
  KinFitter *kf;
  KinFitter *kfMatched;
//  Likelihood *like;
  
  if (! calculateResolutionFunctions)
  {
    kf = new KinFitter("PlotsForResolutionFunctions_testFit_170608_S.root", addWMassKF, addEqMassKF);
    kfMatched = new KinFitter("PlotsForResolutionFunctions_testFit_170608_S.root", addWMassKF, addEqMassKF);
  }
  
  if (makeTGraphs || calculateFractions)
  {
    like = new Likelihood(minCutRedTopMass, maxCutRedTopMass, outputDirLL, dateString, rewHadTopOnly, makeTGraphs, false, true);  // calculateGoodEvtLL, verbose
  }
  if (calculateLikelihood)
  {
    like = new Likelihood(minCutRedTopMass, maxCutRedTopMass, inputDirLL, dateString, rewHadTopOnly, makeTGraphs, false, true);  // calculateGoodEvtLL, verbose
    if (useTTTemplates)
    {
      calculateLikelihood = like->ConstructTGraphsFromFile(dataSetNames, includeDataSets);
    }
    else
    {
      calculateLikelihood = like->ConstructTGraphsFromFile();
      if (! runSystematics) calculateLikelihood = like->ConstructTGraphsFromFile("CorrectMatchLikelihood_");
      if (! runSystematics) calculateLikelihood = like->ConstructTGraphsFromFile("MatchLikelihood_");
    }
    widthsLike.clear();
    widthsLike = like->GetWidths();
    nWidthsLike = widthsLike.size();
    
    if (runSystematics) fileWidths = new TFile(("OutputLikelihood/"+dateString+"/OutputWidths_syst.root").c_str(), "RECREATE");
    else if (runListWidths) fileWidths = new TFile(("OutputLikelihood/"+dateString+"/OutputWidths.root").c_str(), "RECREATE");
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
    MSPlot["nPVs_beforePU_"] = new MultiSamplePlot(datasets, "nPVs_beforePU_", 46, -0.5, 45.5, "# PVs");
    MSPlot["nPVs_afterPU_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_", 46, -0.5, 45.5, "# PVs");
    MSPlot["nPVs_afterPU_up_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_up_", 46, -0.5, 45.5, "# PVs");
    MSPlot["nPVs_afterPU_down_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_down_", 46, -0.5, 45.5, "# PVs");
    MSPlot["rho_"] = new MultiSamplePlot(datasets, "#rho", 41, -0.5, 40.5, "#rho");
    MSPlot["nJets_"] = new MultiSamplePlot(datasets, "nJets_", 13, -0.5, 12.5, "# jets");
    MSPlot["leadingJet_pT_"] = new MultiSamplePlot(datasets, "leadingJet_pT_", 40, 0, 400, "p_{T}", "GeV");
    MSPlot["jet_pT_allJets_"] = new MultiSamplePlot(datasets, "jet_pT_allJets_", 40, 0, 400, "p_{T}", "GeV");
    MSPlot["leadingJet_pT_aKF_"] = new MultiSamplePlot(datasets, "leadingJet_pT_aKF_", 40, 0, 400, "p_{T}", "GeV");
    MSPlot["jet_pT_allJets_aKF_"] = new MultiSamplePlot(datasets, "jet_pT_allJets_aKF_", 40, 0, 400, "p_{T}", "GeV");
    MSPlot["btag_SF_"] = new MultiSamplePlot(datasets, "btag_SF_", 80, 0., 2., "btag SF");
  }
  if (makeLikelihoodPlots)
  {
    InitLikelihoodPlots();
  }
  
  vJER.clear(); vJES.clear(); vPU.clear();
  
  
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
  
  int dMSP;
  bool hasFoundTTbar = false;
  bool doReweighting = false;
  
  /// Loop over systematics or widths
  if (runSystematics)
  {
    if (runRateSystematics) nSystematics = nRateSystematics;
    else if (runSampleSystematics ) nSystematics = nSampleSystematics;
  }
  int endSys = nSystematics;
  if (runListWidths) endSys = nWidths;
  if (! runSystematics && ! runListWidths) endSys = 1;
  for (int iSys = 0; iSys < endSys; iSys++)
  {
    if (runListWidths)
    {
      thisWidth = listWidths[iSys];
      cout << endl << "Running over widths... Now at width " << thisWidth << " x SM width... " << endl;
    }
    else
    {
      thisWidth = scaleWidth;
      cout << endl << "Width for ttbar sample is " << thisWidth << " x SM width..." << endl;
    }
    
    if (runRateSystematics)
    {
      thisSystematic = listRateSyst[iSys];
      cout << endl << "Running over systematics...  " << thisSystematic << endl;
      if (! applyLeptonSF && thisSystematic.find("lepton") != std::string::npos ) continue;
      else if (! applyBTagSF && thisSystematic.find("btag") != std::string::npos ) continue;
      else if (! applyPU && thisSystematic.find("pu") != std::string::npos ) continue;
      else if ( thisSystematic.find("lumiup") != std::string::npos ) Luminosity *= 1.025;
      else if ( thisSystematic.find("lumidown") != std::string::npos ) Luminosity *= 0.975;
    }
    else if (runSampleSystematics)
    {
      thisSystematic = listSampleSyst[iSys];
      cout << endl << "Running over systematics...  " << thisSystematic << endl;
      
      // Temporarily!
      if ( thisSystematic.find("fsr") != std::string::npos || thisSystematic.find("FSR") != std::string::npos || thisSystematic.find("Fsr") != std::string::npos ) continue;
      else if ( thisSystematic.find("herwig") != std::string::npos || thisSystematic.find("HERWIG") != std::string::npos || thisSystematic.find("Herwig") != std::string::npos ) continue;
      else if ( thisSystematic.find("gluonMoveERD") != std::string::npos ) continue;
      else if ( thisSystematic.find("JES") != std::string::npos ) continue;
      else if ( thisSystematic.find("JER") != std::string::npos ) continue;
    }
    
    
    /// Clear counters and likelihood
    nofCM = 0; nofWM = 0; nofNM = 0;
    nofCM_TT = 0; nofWM_TT = 0; nofNM_TT = 0;
    nofCMl = 0; nofWMl = 0; nofNMl = 0;
    nofCM_weighted = 0; nofWM_weighted = 0; nofNM_weighted = 0;
    nofCMout_weighted = 0., nofWMout_weighted = 0., nofNMout_weighted = 0.;
    
    if (calculateLikelihood) like->ClearLikelihoods();
    
    hasFoundTTbar = false;
    doReweighting = false;
    
    eqLumi_TT = -1.;
    
    
    
    ////////////////////////////////////
    /// Loop over datasets
    ////////////////////////////////////
    
    //for (int d = 0; d < 1; d++)
    for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
    {
      clock_t startDataSet = clock();
      
      ClearMetaData();
      
      dataSetName = datasets[d]->Name();
      isData = false; isTTbar = false; isST = false; isSTtW = false; isSTt = false; isOther = false; isHerwig = false;
      if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA") != std::string::npos )
      {
        isData = true;
        pathNtuples = pathNtuplesData;
      }
      else if ( dataSetName.find("TT") != std::string::npos )
      {
        isTTbar = true;
        hasFoundTTbar = true;
        if ( dataSetName.find("width") != std::string::npos || dataSetName.find("Width") != std::string::npos )
        {
          if ( dataSetName.find("x0p2") != std::string::npos ) thisWidth = 0.2;
          else if ( dataSetName.find("x0p5") != std::string::npos ) thisWidth = 0.5;
          else if ( dataSetName.find("x4") != std::string::npos ) thisWidth = 4.;
          else if ( dataSetName.find("x8") != std::string::npos ) thisWidth = 8.;
          applyWidthSF = false;
          if (runListWidths)
          {
            runListWidths = false;
            endSys = 1;
          }
          else if (runSystematics)
          {
            cerr << "ERROR: Trying to run systematics on ttbar sample with non-nominal width... Exiting..." << endl;
            exit(1);
          }
        }
        if (runSampleSystematics) dataSetName = dataSetNameSyst[iSys];
        if ( dataSetName.find("herwig") != std::string::npos || dataSetName.find("HERWIG") != std::string::npos || dataSetName.find("Herwig") != std::string::npos ) isHerwig = true;
        if ( dataSetName.find("isr") != std::string::npos || dataSetName.find("fsr") != std::string::npos || dataSetName.find("hdamp") != std::string::npos || dataSetName.find("erdOn") != std::string::npos || dataSetName.find("ERD") != std::string::npos || dataSetName.find("gluon") != std::string::npos || dataSetName.find("tune") != std::string::npos || dataSetName.find("mass") != std::string::npos || dataSetName.find("mtop") != std::string::npos || isHerwig )
        {
          applyWidthSF = false;
          runListWidths = false;
        }
      }
      else if ( dataSetName.find("ST") != std::string::npos )
      {
        isST = true;
        if ( dataSetName.find("tW") != std::string::npos ) isSTtW = true;
        else isSTt = true;
      }
      else
      {
        isOther = true;
      }
      
      if (verbose > 1)
      {
        cout << "   Dataset " << d << ": " << dataSetName << " / title : " << datasets[d]->Title() << endl;
      }
      
      if (isData && (runListWidths || runSystematics) )
      {
        cout << "Skipping data";
        if (runSystematics) cout << " when running systematics";
        cout << "..." << endl;
        continue;
      }
      
      if (! isData && useTTTemplates && includeDataSets[d] == 0 )
      {
        cout << "Skipping dataset, because not included in likelihood calculation" << endl;
        continue;
      }
      
      doReweighting = false;
      if ( isTTbar && (applyWidthSF || runListWidths) ) doReweighting = true;
      
      if (! isData)
      {
        pathNtuples = pathNtuplesMC;
      }
      
      if (testTTbarOnly && ! isTTbar)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      
      if (! runTTbar && isTTbar)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      if (! runSTtW && isSTtW)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      if (! runSTt && isSTt)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      if (! runOther && isOther)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      
//       if (doPseudoExps && isOther)
//       {
//         cout << "Skipping DY+jets & W+jets datasets when doing pseudo experiments." << endl;
//         continue;
//       }
      
      
      string ntupleFileName = pathNtuples+"Ntuples_"+dataSetName+".root";
      if (runSampleSystematics && isTTbar) ntupleFileName = pathNtuplesSyst+"Ntuples_"+dataSetName+".root";
      
      /// Change name of ttbar dataset to TT (whether it is nominal or widthxX)
      if (isTTbar) dataSetName = "TT";
      
      tFileMap[dataSetName.c_str()] = new TFile(ntupleFileName.c_str(),"READ"); //create TFile for each dataset
      
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
      if (isData) lumiWeight = 1.;
      else
      {
        nEventsDataSet = GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData);
        if (isTTbar)
        {
          nEventsDataSet -= GetNEvents(tStatsTree[(dataSetName).c_str()], "nofTTEventsWithoutAGenTop", isData);
        }
        xSection = datasets[d]->Xsection();  // pb
        eqLumi = (double)nEventsDataSet/xSection;  // 1/pb
        if (isTTbar) eqLumi_TT = eqLumi;
        
        lumiWeight = Luminosity/eqLumi;
        if ( eqLumi_TT != -1. )
          relativeSF = eqLumi_TT/eqLumi;
        else if ( makeTGraphs || calculateLikelihood )
        {
          cerr << "ERROR: Cannot find eqLumi of ttbar sample for relative scale factor... " << endl;
          exit(1);
        }
      }
      
      if (doPseudoExps)
      {
        lumiWeight = 1.;  // select same amount as in data, so no reweighting necessary
        
        if (! isData) toyMax = Luminosity/eqLumi;  // CHECK: better nDataEvts/nEvtsPassKinFit ??
                                                   //else nDataEvts = GetNEvents(tStatsTree[(dataSetName).c_str()], "nEventsSel", 1);
        else nDataEvts = nEventsAKF[d];
        cout << "PseudoExperiments::Number of selected data events: " << nDataEvts << endl;
        
        if (test)
          cout << "      Lumi : " << Luminosity << "/pb; eqLumi: " << eqLumi << "/pb." << endl;
        cout << "PseudoExperiments::Lumi/eqLumi = " << toyMax;
//         if (! isData)
//         {
//           toyMax *= 0.908;  // small overshoot in MC --> scale down --> fixed now
//           cout << " x 0.908 = " << toyMax;
//         }
        cout << endl;
      }
      
      if (calculateAverageMass)
      {
        txtMassReco.open(("averageMass/mass_reco_"+dataSetName+"_"+dateString+".txt").c_str());
      }
      
      
      /// Get data
      tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
      nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
      cout << "                nEntries  : " << nEntries << endl;
      if (isData) cout << "                Lumi    : " << Luminosity << "/pb" << endl;
      else
      {
        cout << "                eqLumi    : " << eqLumi << "/pb = " << nEventsDataSet << " / " << xSection << " pb" << endl;
        if (isTTbar)
          cout << "                              ( " << nEventsDataSet << " = " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData) << " - " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofTTEventsWithoutAGenTop", isData) << " )" << endl;
        cout << "                lumiWeight: " << lumiWeight << endl;
      }
      
      if (test && isData)
      {
        cout << "NEvents data:   " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData) << endl;
        cout << "NEvents Run B:  " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunB", isData) << endl;
        cout << "NEvents Run CD: " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunCD", isData) << endl;
        cout << "NEvents Run EF: " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunEF", isData) << endl;
        cout << "NEvents Run G:  " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunG", isData) << endl;
        cout << "NEvents Run H:  " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunH", isData) << endl;
      }
      
      if (runRateSystematics && isTTbar && thisSystematic.find("renFac") != std::string::npos )
      {
        renFacSumNom = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1001", isData);
        renFacSum1002 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1002", isData);
        renFacSum1003 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1003", isData);
        renFacSum1004 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1004", isData);
        renFacSum1005 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1005", isData);
        renFacSum1007 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1007", isData);
        renFacSum1009 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1009", isData);
        cout << endl << "            Sum weight1001: " << (long)renFacSumNom << endl;
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
          thisLeptonIdSF = fracDataEras[0]*muonIdSF_BCDEF[0] + fracDataEras[1]*muonIdSF_GH[0];
          thisLeptonIsoSF = fracDataEras[0]*muonIsoSF_BCDEF[0] + fracDataEras[1]*muonIsoSF_GH[0];
          thisLeptonTrigSF = fracDataEras[0]*muonTrigSF_BCDEF[0] + fracDataEras[1]*muonTrigSF_GH[0];
          thisLeptonSF = muonTrackSF_eta[0] * thisLeptonIdSF * thisLeptonIsoSF * thisLeptonTrigSF;
          if (runRateSystematics)
          {
            if ( thisSystematic.find("lepton") != std::string::npos )  // apply some lepton SF systematic
            {
              if (applyBTagSF) { scaleFactor *= btagSF;}
              if (applyPU) { scaleFactor *= puSF;}
              if ( thisSystematic.find("IdSFup") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * (fracDataEras[0]*muonIdSF_up_BCDEF[0] + fracDataEras[1]*muonIdSF_up_GH[0]) * thisLeptonIsoSF * thisLeptonTrigSF;}
              else if ( thisSystematic.find("IdSFdown") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * (fracDataEras[0]*muonIdSF_down_BCDEF[0] + fracDataEras[1]*muonIdSF_down_GH[0]) * thisLeptonIsoSF * thisLeptonTrigSF;}
              else if ( thisSystematic.find("IsoSFup") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * thisLeptonIdSF * (fracDataEras[0]*muonIsoSF_up_BCDEF[0] + fracDataEras[1]*muonIsoSF_up_GH[0]) * thisLeptonTrigSF;}
              else if ( thisSystematic.find("IsoSFdown") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * thisLeptonIdSF * (fracDataEras[0]*muonIsoSF_down_BCDEF[0] + fracDataEras[1]*muonIsoSF_down_GH[0]) * thisLeptonTrigSF;}
              else if ( thisSystematic.find("TrigSFup") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * thisLeptonIdSF * thisLeptonIsoSF * (fracDataEras[0]*muonTrigSF_up_BCDEF[0] + fracDataEras[1]*muonTrigSF_up_GH[0]);}
              else if ( thisSystematic.find("TrigSFdown") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * thisLeptonIdSF * thisLeptonIsoSF * (fracDataEras[0]*muonTrigSF_down_BCDEF[0] + fracDataEras[1]*muonTrigSF_down_GH[0]);}
              else if ( thisSystematic.find("TrkSFup") != std::string::npos ) { scaleFactor *= 1.01*thisLeptonSF;}
              else if ( thisSystematic.find("TrkSFdown") != std::string::npos ) { scaleFactor *= 0.99*thisLeptonSF;}
            }
            else if ( thisSystematic.find("btag") != std::string::npos )
            {
              if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
              if (applyPU) { scaleFactor *= puSF;}
              if ( thisSystematic.find("up") != std::string::npos ) { scaleFactor *= btagSF_up;}
              else if ( thisSystematic.find("down") != std::string::npos ) { scaleFactor *= btagSF_down;}
            }
            else if ( thisSystematic.find("pu") != std::string::npos )
            {
              if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
              if (applyBTagSF) { scaleFactor *= btagSF;}
              if ( thisSystematic.find("up") != std::string::npos ) { scaleFactor *= puSF_up;}
              else if ( thisSystematic.find("down") != std::string::npos ) { scaleFactor *= puSF_down;}
            }
            else if ( thisSystematic.find("renFac") != std::string::npos && isTTbar )
            {
              if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
              if (applyBTagSF) { scaleFactor *= btagSF;}
              if (applyPU) { scaleFactor *= puSF;}
              
              if ( thisSystematic.find("1002") != std::string::npos ) scaleFactor *= weight1002 * renFacSumNom/renFacSum1002;
              else if ( thisSystematic.find("1003") != std::string::npos ) scaleFactor *= weight1003 * renFacSumNom/renFacSum1003;
              else if ( thisSystematic.find("1004") != std::string::npos ) scaleFactor *= weight1004 * renFacSumNom/renFacSum1004;
              else if ( thisSystematic.find("1005") != std::string::npos ) scaleFactor *= weight1005 * renFacSumNom/renFacSum1005;
              else if ( thisSystematic.find("1007") != std::string::npos ) scaleFactor *= weight1007 * renFacSumNom/renFacSum1007;
              else if ( thisSystematic.find("1009") != std::string::npos ) scaleFactor *= weight1009 * renFacSumNom/renFacSum1009;
            }
            else
            {
              if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
              if (applyBTagSF) { scaleFactor *= btagSF;}
              if (applyPU) { scaleFactor *= puSF;}
            }
          }
          else if (runSampleSystematics)
          {
            if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
            if (applyBTagSF) { scaleFactor *= btagSF;}
            if (applyPU) { scaleFactor *= puSF;}
          }
          else
          {
            if (makePlots && passedMETFilter)
            {
              MSPlot["nPVs_beforePU_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF);
              MSPlot["nPVs_afterPU_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF*puSF);
              MSPlot["nPVs_afterPU_up_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF*puSF_up);
              MSPlot["nPVs_afterPU_down_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF*puSF_down);
            }
            
            if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
            if (applyBTagSF) { scaleFactor *= btagSF;}
            if (applyPU) { scaleFactor *= puSF;}
            
            if (makePlots)
            {
              MSPlot["btag_SF_"]->Fill(btagSF, datasets[d], true, lumiWeight*scaleFactor);
            }
            if ( thisWidth == 1 && applyPU && puSF == 0 ) txtDebugPUSF << nvtx << "    " << npu << endl;
          }
        }
        else if (makePlots && passedMETFilter)
        {
          MSPlot["nPVs_beforePU_"]->Fill(nvtx, datasets[d], false, lumiWeight);
          MSPlot["nPVs_afterPU_"]->Fill(nvtx, datasets[d], false, lumiWeight);
          MSPlot["nPVs_afterPU_up_"]->Fill(nvtx, datasets[d], false, lumiWeight);
          MSPlot["nPVs_afterPU_down_"]->Fill(nvtx, datasets[d], false, lumiWeight);
        }
        
        
        
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
        
        if (makePlots && passedMETFilter)
        {
          MSPlot["nJets_"]->Fill(selectedJets.size(), datasets[d], true, lumiWeight*scaleFactor);
          MSPlot["rho_"]->Fill(rho, datasets[d], true, lumiWeight*scaleFactor);
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
        
        
        if (! doReweighting ) widthSF = 1.;
        //else if ( applyWidthSF && ! isTTbar ) widthSF = 1.;  // also for data
        
        
        
        /////////////////////////////
        ///  JET PARTON MATCHING  ///
        /////////////////////////////
        
        //if ( isTTbar || dataSetName.find("ST") != std::string::npos )
        if (! isData && ! isHerwig)
        {
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
              
              if (makePlots && isTTbar)
              {
                histo1D["genTop_status"]->Fill(mc_status[i]);
                if ( mc_status[i] == 22 )
                  histo1D["genTop_status22_pT"]->Fill(mc_pt[i]);
                else if ( mc_status[i] == 62 )
                  histo1D["genTop_status62_pT"]->Fill(mc_pt[i]);
              }
              if (isTTbar && mc_isLastCopy[i])
              {
                foundLastCopyTop = true;
                if (makePlots) histo1D["genTop_isLastCopy_status"]->Fill(mc_status[i]);
              }
            }
            else if ( mc_pdgId[i] == -pdgID_top )
            {
              if ( mc_status[i] == 22 )
              {
                antiTopQuark = i;
                foundAntiTop22 = true;
              }
              if ( antiTopQuark == -9999 && mc_status[i] == 62 ) antiTopQuark = i;
              
              if (makePlots && isTTbar)
              {
                histo1D["genAntitop_status"]->Fill(mc_status[i]);
                if ( mc_status[i] == 22 )
                  histo1D["genAntitop_status22_pT"]->Fill(mc_pt[i]);
                else if ( mc_status[i] == 62 )
                  histo1D["genAntitop_status62_pT"]->Fill(mc_pt[i]);
              }
              if (isTTbar && mc_isLastCopy[i])
              {
                foundLastCopyAntitop = true;
                if (makePlots) histo1D["genAntitop_isLastCopy_status"]->Fill(mc_status[i]);
              }
            }
            
            if ( mc_pdgId[i] == pdgID_top && mc_status[i] == 62 )  // top
            {
              foundTop62 = true;
              if ( mcParticles[i].Pt() < 800. ) topPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
              else topPtSF = TMath::Exp(0.0615-0.0005*800.);
            }
            else if ( mc_pdgId[i] == -pdgID_top && mc_status[i] == 62 )  // antitop
            {
              foundAntiTop62 = true;
              if ( mcParticles[i].Pt() < 800. ) antiTopPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
              else antiTopPtSF = TMath::Exp(0.0615-0.0005*800.);
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
            if ( abs(mc_pdgId[i]) < 6 || abs(mc_pdgId[i]) == 21 )  //light/b quarks, 6 should stay hardcoded, OR gluon
            {
              partons.push_back(mcParticles[i]);
              partonId.push_back(i);  /// partons[j] = mcParticles[partonId[j]]
            }
            
          }  // end loop mcParticles
          
          if (makePlots && isTTbar)
          {
            if (foundLastCopyTop) histo1D["genTop_hasLastCopy"]->Fill(1.);
            else histo1D["genTop_hasLastCopy"]->Fill(0.);
            if (foundLastCopyAntitop) histo1D["genAntitop_hasLastCopy"]->Fill(1.);
            else histo1D["genAntitop_hasLastCopy"]->Fill(0.);
            if (foundTop22) histo1D["genTop_hasStatus22"]->Fill(1.);
            else histo1D["genTop_hasStatus22"]->Fill(0.);
            if (foundAntiTop22) histo1D["genAntitop_hasStatus22"]->Fill(1.);
            else histo1D["genAntitop_hasStatus22"]->Fill(0.);
          }
          
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
          
          if ( isTTbar && (topQuark == -9999 || antiTopQuark == -9999) )
          {
            if ( thisWidth == 1 )
            {
              txtDebugTopMass << "Event " << ievt << ";  Event nb. " << evt_num << "; Run nb. " << run_num << endl;
              txtDebugTopMass << "Top mass id: " << topQuark << "; antiTop mass id: " << antiTopQuark << "; Lepton charge: " << muon_charge[0] << endl;
            }
            continue;
          }
//           else if ( dataSetName.find("ST_tW") != std::string::npos && ( (foundMuMinus && topQuark == -9999) || (foundMuPlus && antiTopQuark == -9999) ) )
//           {
//             if ( thisWidth == 1 )
//             {
//               txtDebugTopMass << "Event " << ievt << ";  Event nb. " << evt_num << "; Run nb. " << run_num << endl;
//               txtDebugTopMass << "Top mass id: " << topQuark << "; antiTop mass id: " << antiTopQuark << "; Lepton charge: " << muon_charge[0] << endl;
//             }
//             continue;
//           }
          
          if (runSystematics)
          {
            doMatching = false;
            if ( thisSystematic.find("topPtRew") != std::string::npos )
            {
              if (isTTbar) topPtRewSF = TMath::Sqrt(topPtSF*antiTopPtSF);
              else topPtRewSF = 1.;
              scaleFactor *= topPtRewSF;
            }
          }
          
          
          
          /////////////////////////////////////////
          ///  Scale factor ttbar sample width  ///
          /////////////////////////////////////////
          
          if ( muon_charge[0] > 0 )
          {
            massHadTopQ = (mcParticles[antiTopQuark]).M();
            massLepTopQ = (mcParticles[topQuark]).M();
          }
          else if ( muon_charge[0] < 0 )
          {
            massHadTopQ = (mcParticles[topQuark]).M();
            massLepTopQ =  (mcParticles[antiTopQuark]).M();
          }
          
          if ( doReweighting )
          {
            if (rewHadTopOnly) widthSF = rew->EventWeightCalculatorNonRel(massHadTopQ, thisWidth);
            else widthSF = rew->EventWeightCalculatorNonRel(massHadTopQ, thisWidth) * rew->EventWeightCalculatorNonRel(massLepTopQ, thisWidth);
            
            if ( widthSF != widthSF )  // widthSF = NaN
            {
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
                txtMassGenPMatched << ievt << "  " << matched_top_mass_q << "  " << scaleFactor << "  " << lumiWeight << endl;
                txtMassGenJMatched << ievt << "  " << matched_top_mass_j << "  " << scaleFactor << "  " << lumiWeight << endl;
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
                
                if ( applyKinFitCut && kFitChi2Matched < kinFitCutValue ) passKFChi2MatchedCut = true;
                if (passKFChi2MatchedCut)
                {
                  nofAcceptedKFitMatched++;
                  
                  selectedJetsKFMatched.clear();
                  selectedJetsKFMatched = kfMatched->getCorrectedJets();
                  
                  if ( selectedJetsKFMatched.size() == 2 ) selectedJetsKFMatched.push_back(jetsMatched[2]);
                  
                  matched_top_mass_j_akF = (selectedJetsKFMatched[0] + selectedJetsKFMatched[1] + selectedJetsKFMatched[2]).M();
                  
                  if (calculateLikelihood)
                  {
                    double temp = matched_top_mass_j_akF/aveTopMassLL;
                    if (! useTTTemplates && ! runSystematics)
                      like->CalculateGenLikelihood(temp, massHadTopQ, massLepTopQ, thisWidth, doReweighting, isData);
                  }
                }  // passKFChi2MatchedCut
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
        
        dMSP = d;
        if (hasFoundTTbar && ! isTTbar) dMSP = d+2;
        else if (isTTbar && isWM) dMSP = d+1;
        else if (isTTbar && isNM) dMSP = d+2;
        
        
        /// Fill variables before performing kinFit
        reco_W_mass_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).M();
        reco_top_mass_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).M();
        reco_top_pt_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).Pt();
        reco_mlb_bKF = (selectedLepton[0] + selectedJets[labelsReco[3]]).M();
        reco_dRLepB_lep_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[3]], selectedLepton[0] );  // deltaR between lepton and leptonic b jet
        reco_dRLepB_had_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
        reco_ttbar_mass_bKF = reco_mlb_bKF + reco_top_mass_bKF;
        redTopMass_bKF = reco_top_mass_bKF/aveTopMassLL;
        
        
        
        if (makePlots)
        {
          if (isTTbar && doKinFit) FillKinFitPlots(doneKinFit);
          
          FillMSPlots(dMSP, doneKinFit);
          //FillControlPlots(datasetsMSP, dMSP);
          
          MSPlot["leadingJet_pT_"]->Fill(selectedJets[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
          for (int iJet = 0; iJet < selectedJets.size(); iJet++)
          {
            MSPlot["jet_pT_allJets_"]->Fill(selectedJets[iJet].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
          }
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
          if (makePlots)
          {
            MSPlot["KF_Chi2_wide"]->Fill(kFitChi2, datasetsMSP[dMSP], true, lumiWeight*scaleFactor*widthSF);
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
        
        if (calculateAverageMass) txtMassReco << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
        
        if ( doKinFit && makePlots )
        {
          if (isTTbar) FillKinFitPlots(doneKinFit);
          if (! isData)
          {
            histo1D["allSim_top_mass"]->Fill(reco_top_mass_aKF, lumiWeight*scaleFactor*widthSF);
            histo1D["allSim_red_top_mass"]->Fill(redTopMass, lumiWeight*scaleFactor*widthSF);
          }
        }
        
        
        
        ////////////////////
        ///  Likelihood  ///
        ////////////////////
        
        if (makeTGraphs) like->FillHistograms(redTopMass, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isData, catSuffix);
        if (calculateLikelihood)
        {
          loglike_per_evt = like->CalculateLikelihood(redTopMass, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, thisWidth, doReweighting, isData);
          if ( ! doPseudoExps && ! useTTTemplates && ! runSystematics && isCM )  // isCM ensures ! isData
            like->CalculateCMLikelihood(redTopMass, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, doReweighting, isData);
          if ( ! doPseudoExps && ! useTTTemplates && ! runSystematics && ( isCM || isWM ) )
            like->CalculateTempLikelihood(redTopMass, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, doReweighting, isData);
        }
        
        /// Calculate fraction of events in category outside interval
        if (! isData)
        {
          if (isCM) nofCMout_weighted += lumiWeight*scaleFactor*widthSF;
          else if (isWM) nofWMout_weighted += lumiWeight*scaleFactor*widthSF;
          else if (isNM) nofNMout_weighted += lumiWeight*scaleFactor*widthSF;
        }
        
        if ( redTopMass > maxRedTopMass ) maxRedTopMass = redTopMass;
        if ( redTopMass < minRedTopMass ) minRedTopMass = redTopMass;
        if ( ! isData && redTopMass > minCutRedTopMass && redTopMass < maxCutRedTopMass )
        {
          if (isCM)
          {
            nofCM++;
            nofCMl += lumiWeight*scaleFactor;
            nofCM_weighted += lumiWeight*scaleFactor*widthSF;
            if (isTTbar) nofCM_TT++;
          }
          else if (isWM)
          {
            nofWM++;
            nofWMl += lumiWeight*scaleFactor;
            nofWM_weighted += lumiWeight*scaleFactor*widthSF;
            if (isTTbar) nofWM_TT++;
          }
          else if (isNM)
          {
            nofNM++;
            nofNMl += lumiWeight*scaleFactor;
            nofNM_weighted += lumiWeight*scaleFactor*widthSF;
            if (isTTbar) nofNM_TT++;
          }
          
          if (calculateFractions)
          {
            like->AddToFraction(d, lumiWeight*scaleFactor, massHadTopQ, massLepTopQ, doReweighting, isCM, isWM, isNM);
          }
          
//           if (isTTbar && makeLikelihoodPlots)
//           {
//             FillLikelihoodPlots();
//           }
        }
        
        if (calculateAverageMass && ! isData)
        {
          if (isCM) txtMassRecoCM << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
          else
          {
            txtMassRecoWMNM << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
            if (isWM)
              txtMassRecoWM << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
            else if (isNM)
              txtMassRecoNM << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
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
            nEvtsInPseudoExpW[iPsExp][d] += scaleFactor;
            if ( redTopMass > minCutRedTopMass && redTopMass < maxCutRedTopMass )
              like->AddPsExp(iPsExp, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, doReweighting, isData);
            
            /// Fill plots only for first pseudo experiment
            if ( makePlots && iPsExp == 0 )
            {
              /// Combine DY & W+jets
              if ( dataSetName.find("DY") != std::string::npos ) histo1D["red_top_mass_DYJets"]->Fill(redTopMass);
              else if ( dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos )
                histo1D["red_top_mass_WJets"]->Fill(redTopMass);
              else histo1D[("red_top_mass_"+dataSetName).c_str()]->Fill(redTopMass);
              
              FillCatsPlots(catSuffix);
              
              int dMSP = d;
              if (hasFoundTTbar && ! isTTbar) dMSP = d+2;
              else if (isTTbar && isWM) dMSP = d+1;
              else if (isTTbar && isNM) dMSP = d+2;
              
              FillMSPlots(dMSP, doneKinFit);
            }
          }
        }
        
        //Fill histos
        if ( makePlots && (! doPseudoExps || isData) )
        {
          /// Combine DY & W+jets
          if ( dataSetName.find("DY") != std::string::npos ) histo1D["red_top_mass_DYJets"]->Fill(redTopMass);
          else if ( dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos )
            histo1D["red_top_mass_WJets"]->Fill(redTopMass);
          else histo1D[("red_top_mass_"+dataSetName).c_str()]->Fill(redTopMass);
          
          FillCatsPlots(catSuffix);
          
          FillMSPlots(dMSP, doneKinFit);
          
          MSPlot["leadingJet_pT_aKF_"]->Fill(selectedJetsAKF[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
          for (int iJet = 0; iJet < selectedJetsAKF.size(); iJet++)
          {
            MSPlot["jet_pT_allJets_aKF_"]->Fill(selectedJetsAKF[iJet].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
          }
          
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
            cout << " --- Kinematic fit --- Before chi2 cut --- " << endl;
            cout << "Correctly matched reconstructed events    : " << setw(8) << right << nofCorrectlyMatchedAKFNoCut << endl;
            cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatchedAKFNoCut << endl;
            if ( nofCorrectlyMatchedAKFNoCut != 0 || nofNotCorrectlyMatchedAKFNoCut != 0 )
              cout << "   ===> This means that " << 100*(float)nofCorrectlyMatchedAKFNoCut / (float)(nofCorrectlyMatchedAKFNoCut + nofNotCorrectlyMatchedAKFNoCut) << "% of matched events is correctly matched after KF." << endl;
            
            cout << "                        " << 100*(float)nofCorrectlyMatchedAKFNoCut / (float)(nofNotCorrectlyMatchedAKFNoCut+nofNoMatchAKFNoCut) << "% of all events accepted by kinfitter is correctly matched." << endl;
            
            cout << " --- Kinematic fit --- After chi2 cut --- " << endl;
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
          cout << " --- Kinematic fit --- Gen events" << endl;
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
      
      
      
      if (calculateAverageMass) txtMassReco.close();
      
      tFileMap[dataSetName.c_str()]->Close();
      
      timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
      
    }  // end loop datasets
    
    if (! doGenOnly && ! testTTbarOnly)
    {
      cout << "Number of events with " << minCutRedTopMass << " < mt/<mt> < " << maxCutRedTopMass << " : CM: " << nofCM << " (" << 100*(double)nofCM/((double)(nofCM+nofWM+nofNM)) << "%)   WM: " << nofWM << " (" << 100*(double)nofWM/((double)(nofCM+nofWM+nofNM)) << "%)   NM: " << nofNM << " (" << 100*(double)nofNM/((double)(nofCM+nofWM+nofNM)) << "%)   Total: " << nofCM+nofWM+nofNM << endl;
      cout << "Number of events with " << minCutRedTopMass << " < mt/<mt> < " << maxCutRedTopMass << " : CM: " << nofCMl << " (" << 100*nofCMl/(nofCMl+nofWMl+nofNMl) << "%)   WM: " << nofWMl << " (" << 100*nofWMl/(nofCMl+nofWMl+nofNMl) << "%)   NM: " << nofNMl << " (" << 100*nofNMl/(nofCMl+nofWMl+nofNMl) << "%)   Total: " << nofCMl+nofWMl+nofNMl << endl;
      cout << "                                  weighted: CM: " << nofCM_weighted << " (" << 100*nofCM_weighted/(nofCM_weighted+nofWM_weighted+nofNM_weighted) << "%)   WM: " << nofWM_weighted << " (" << 100*nofWM_weighted/(nofCM_weighted+nofWM_weighted+nofNM_weighted) << "%)   NM: " << nofNM_weighted << " (" << 100*nofNM_weighted/(nofCM_weighted+nofWM_weighted+nofNM_weighted) << "%)   Total: " << (int)(nofCM_weighted+nofWM_weighted+nofNM_weighted) << endl;
      cout << "                               (TTbar only) CM: " << nofCM_TT << "               WM: " << nofWM_TT << "               NM: " << nofNM_TT << endl;
      cout << endl << "Number of events outside interval:          CM: " << nofCMout_weighted << " (" << 100*nofCMout_weighted/(nofCMout_weighted+nofWMout_weighted+nofNMout_weighted) << "%)   WM: " << nofWMout_weighted << " (" << 100*nofWMout_weighted/(nofCMout_weighted+nofWMout_weighted+nofNMout_weighted) << "%)   NM: " << nofNMout_weighted << " (" << 100*nofNMout_weighted/(nofCMout_weighted+nofWMout_weighted+nofNMout_weighted) << "%)   Total: " << (int)(nofCMout_weighted+nofWMout_weighted+nofNMout_weighted) << endl;
    }
    
    
    if (calculateFractions)
    {
      like->CalculateFractions(dataSetNames);
    }
    if (makeTGraphs)
    {
      like->WriteHistograms("ReducedTopMassPlots.root");
      like->ConstructTGraphsFromHisto("TGraphFunctions.root", dataSetNames, includeDataSets);
    }
    
    if (calculateLikelihood)
    {
      cout << "Minimum reduced top mass: " << minRedTopMass << endl;
      cout << "Maximum reduced top mass: " << maxRedTopMass << endl;
      
      /// Print output to file
      string llFileName = "output_loglikelihood_widthx"+DotReplace(thisWidth);
      if (doGenOnly) llFileName = "output_loglikelihood_parton_widthx"+DotReplace(thisWidth);
      //if (useToys) llFileName = "output_loglikelihood_toys";
      if (runSystematics) llFileName = "output_loglikelihood_"+thisSystematic;
      like->PrintLikelihoodOutput(llFileName+".txt");
      if (unblind) like->PrintLikelihoodOutputData(llFileName+"_data.txt");
      like->PrintMtmLikelihoodOutput(llFileName+"_Mtm.txt");
      
      /// Calculate output width
      if ( runSystematics || runListWidths)
      {
        if (runSystematics) cout << "Output width for " << thisSystematic << ": " << endl;
        else cout << "Standard output width: " << endl;
        fileWidths->cd();
        if (runSystematics) like->GetOutputWidth(thisWidth, thisSystematic, true, false);
        else like->GetOutputWidth(thisWidth, true, false);
        if (! useTTTemplates && ! runSystematics)
        {
          cout << "Output width for correctly matched events (using likelihood with only CM template): " << endl;
          like->GetOutputWidth(thisWidth, "CM", true, false);
          cout << "Output width for correctly & wrongly matched events (using likelihood with only CM & WM templates): " << endl;
          like->GetOutputWidth(thisWidth, "matched", true, false);
        }
      }
      else
      {
        cout << "Standard output width: " << endl;
        like->GetOutputWidth(thisWidth, true, true);
        if (unblind) like->GetOutputWidth(thisWidth, "data", true, true);
        if (! doPseudoExps && ! useTTTemplates)
        {
          cout << "Output width for correctly matched events (using likelihood with only CM template): " << endl;
          like->GetOutputWidth(thisWidth, "CM", true, true);
          cout << "Output width for correctly & wrongly matched events (using likelihood with only CM & WM templates): " << endl;
          like->GetOutputWidth(thisWidth, "matched", true, true);
          //cout << "Output width for generated events (using likelihood with only CM template): " << endl;
          //like->GetOutputWidth(scaleWidth, "gen", true);
        }
        //cout << "Output width from file (standard calculation): " << endl;
        //like->GetOutputWidth(llFileName+".txt", scaleWidth, true);
      }
    }
    if (calculateFractions)
    {
      like->CalculateFractions(dataSetNames);
    }
    
    if (doPseudoExps)
    {
      like->CalculatePull(thisWidth);
    }
    
    
    cout << endl << "Processing time per dataset: " << endl;
    for (unsigned int d = 0; d < datasets.size(); d++)
    {
      cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
    }
    
    if (doPseudoExps)
    {
      cout << endl << "PseudoExperiments::Number of selected data events: " << nDataEvts << endl;
      double totMCPsExp, totMCPsExpW;
      for (int i = 0; i < 3; i++)
      {
        totMCPsExp = 0;
        totMCPsExpW = 0.;
        cout << "PseudoExperiment " << std::setw(3) << std::right << i << endl;
        for (unsigned int d = 1; d < datasets.size(); d++)
        {
          totMCPsExp += nEvtsInPseudoExp[i][d];
          totMCPsExpW += nEvtsInPseudoExpW[i][d];
          cout << "                      " << datasets[d]->Name() << ": " << nEvtsInPseudoExp[i][d] << endl;
        }
        cout << "                                            Total MC: " << totMCPsExp << endl;
        cout << "                                           (weighted: " << totMCPsExpW << ")" << endl;
      }
    }
    
    
    ///  Check Shape Changing Systematics
    if (! testTTbarOnly && runAll && ! runListWidths && ! runSystematics && ! useTTTemplates && ! doPseudoExps)
      CheckSystematics(vJER, vJES, vPU);
    
    
  }  // end loop systematics/widths
  
  if (runSystematics || runListWidths)
  {
    fileWidths->Close();
    delete fileWidths;
  }
  
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
    string rootFileName = "NtuplePlots_"+systStr+".root";
    if (! doGenOnly)
    {
      mkdir((pathOutput+"MSPlot/").c_str(),0777);
      if (makeControlPlots) mkdir((pathOutput+"MSPlot/ControlPlots/").c_str(),0777);
    }
    
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
      
      for (map<string,MultiSamplePlot*>::const_iterator it = MSPlotCP.begin(); it != MSPlotCP.end(); it++)
      {
        //cout << "MSPlot: " << it->first << endl;
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        temp->Draw(name, 1, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
        temp->Write(fout, name, true, pathOutput+"MSPlot/ControlPlots", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
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
  }
  if (makeLikelihoodPlots) WriteLikelihoodPlots();
  
  
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
  if (! isData) tree->SetBranchAddress("cutFlowWeighted", cutFlowWeighted, &b_cutFlowWeighted);     // temporarily!
  if (! isData) tree->SetBranchAddress("cutFlow2Weighted", cutFlow2Weighted, &b_cutFlow2Weighted);  // temporarily!
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
    tree->SetBranchAddress("nofTTEventsWithoutAGenTop", &nofTTEventsWithoutAGenTop, &b_nofTTEventsWithoutAGenTop);
    tree->SetBranchAddress("nofTTEventsWithoutGenTop", &nofTTEventsWithoutGenTop, &b_nofTTEventsWithoutGenTop);
    tree->SetBranchAddress("nofTTEventsWithoutGenAntiTop", &nofTTEventsWithoutGenAntiTop, &b_nofTTEventsWithoutGenAntiTop);
    tree->SetBranchAddress("nofTTEventsWithoutBothGenTopsWithStatus22", &nofTTEventsWithoutBothGenTopsWithStatus22, &b_nofTTEventsWithoutBothGenTopsWithStatus22);
    tree->SetBranchAddress("nofTTEventsWithoutGenTopWithStatus22", &nofTTEventsWithoutGenTopWithStatus22, &b_nofTTEventsWithoutGenTopWithStatus22);
    tree->SetBranchAddress("nofTTEventsWithoutGenAntiTopWithStatus22", &nofTTEventsWithoutGenAntiTopWithStatus22, &b_nofTTEventsWithoutGenAntiTopWithStatus22);
    tree->SetBranchAddress("nofTTEventsWithoutBothGenTopsWithStatus62", &nofTTEventsWithoutBothGenTopsWithStatus62, &b_nofTTEventsWithoutBothGenTopsWithStatus62);
    tree->SetBranchAddress("nofTTEventsWithoutGenTopWithStatus62", &nofTTEventsWithoutGenTopWithStatus62, &b_nofTTEventsWithoutGenTopWithStatus62);
    tree->SetBranchAddress("nofTTEventsWithoutGenAntiTopWithStatus62", &nofTTEventsWithoutGenAntiTopWithStatus62, &b_nofTTEventsWithoutGenAntiTopWithStatus62);
    tree->SetBranchAddress("sumWeight1001", &sumWeight1001, &b_sumWeight1001);
    tree->SetBranchAddress("sumWeight1002", &sumWeight1002, &b_sumWeight1002);
    tree->SetBranchAddress("sumWeight1003", &sumWeight1003, &b_sumWeight1003);
    tree->SetBranchAddress("sumWeight1004", &sumWeight1004, &b_sumWeight1004);
    tree->SetBranchAddress("sumWeight1005", &sumWeight1005, &b_sumWeight1005);
    tree->SetBranchAddress("sumWeight1007", &sumWeight1007, &b_sumWeight1007);
    tree->SetBranchAddress("sumWeight1009", &sumWeight1009, &b_sumWeight1009);
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
    if (isTTbar)
    {
      tree->SetBranchAddress("weight1001", &weight1001, &b_weight1001);
      tree->SetBranchAddress("weight1002", &weight1002, &b_weight1002);
      tree->SetBranchAddress("weight1003", &weight1003, &b_weight1003);
      tree->SetBranchAddress("weight1004", &weight1004, &b_weight1004);
      tree->SetBranchAddress("weight1005", &weight1005, &b_weight1005);
      tree->SetBranchAddress("weight1007", &weight1007, &b_weight1007);
      tree->SetBranchAddress("weight1009", &weight1009, &b_weight1009);
    }
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
  /// Control plots
  if (makeControlPlots)
  {
    MSPlotCP["muon_pT"] = new MultiSamplePlot(datasetsMSP, "muon_pT", 22, 0, 440, "p_{T}", "GeV");
    MSPlotCP["muon_eta"] = new MultiSamplePlot(datasetsMSP, "muon_eta", 30, -3, 3, "Eta");
    MSPlotCP["muon_phi"] = new MultiSamplePlot(datasetsMSP, "muon_phi", 32, -3.2, 3.2, "Phi");
    MSPlotCP["muon_relIso"] = new MultiSamplePlot(datasetsMSP, "muon_relIso", 20, 0, 0.2, "relIso");
    MSPlotCP["muon_d0"] = new MultiSamplePlot(datasetsMSP, "muon_d0", 60, 0, 0.003, "d_{0}");
    MSPlotCP["leadingJet_pT"] = new MultiSamplePlot(datasetsMSP, "leadingJet_pT", 60, 0, 600, "p_{T}", "GeV");
    MSPlotCP["jet2_pT"] = new MultiSamplePlot(datasetsMSP, "jet2_pT", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet3_pT"] = new MultiSamplePlot(datasetsMSP, "jet3_pT", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet4_pT"] = new MultiSamplePlot(datasetsMSP, "jet4_pT", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet_pT_allJets"] = new MultiSamplePlot(datasetsMSP, "jet_pT_allJets", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["Ht_4leadingJets"] = new MultiSamplePlot(datasetsMSP,"Ht_4leadingJets", 60, 0, 1200, "H_{T}", "GeV");
    MSPlotCP["met_pT"] = new MultiSamplePlot(datasetsMSP, "met_pT", 40, 0, 400, "E_{T}^{miss} p_{T}", "GeV");
    MSPlotCP["met_eta"] = new MultiSamplePlot(datasetsMSP, "met_eta", 30, -3, 3, "E_{T}^{miss} Eta");
    MSPlotCP["met_phi"] = new MultiSamplePlot(datasetsMSP, "met_phi", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
    MSPlotCP["met_corr_pT"] = new MultiSamplePlot(datasetsMSP, "met_corr_pT", 40, 0, 400, "E_{T}^{miss} p_{T}", "GeV");
    MSPlotCP["met_corr_eta"] = new MultiSamplePlot(datasetsMSP, "met_corr_eta", 30, -3, 3, "E_{T}^{miss} Eta");
    MSPlotCP["met_corr_phi"] = new MultiSamplePlot(datasetsMSP, "met_corr_phi", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
    
    MSPlotCP["muon_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_pT_aKF", 22, 0, 440, "p_{T}", "GeV");
    MSPlotCP["muon_eta_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_eta_aKF", 30, -3, 3, "Eta");
    MSPlotCP["muon_phi_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_phi_aKF", 32, -3.2, 3.2, "Phi");
    MSPlotCP["muon_relIso_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_relIso_aKF", 20, 0, 0.2, "relIso");
    MSPlotCP["muon_d0_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_d0_aKF", 60, 0, 0.003, "d_{0}");
    MSPlotCP["leadingJet_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "leadingJet_pT_aKF", 60, 0, 600, "p_{T}", "GeV");
    MSPlotCP["jet2_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "jet2_pT_aKF", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet3_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "jet3_pT_aKF", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet4_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "jet4_pT_aKF", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet_pT_allJets_aKF"] = new MultiSamplePlot(datasetsMSP, "jet_pT_allJets_aKF", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["Ht_4leadingJets_aKF"] = new MultiSamplePlot(datasetsMSP,"Ht_4leadingJets_aKF", 60, 0, 1200, "H_{T}", "GeV");
    MSPlotCP["met_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "met_pT_aKF", 40, 0, 400, "E_{T}^{miss} p_{T}", "GeV");
    MSPlotCP["met_eta_aKF"] = new MultiSamplePlot(datasetsMSP, "met_eta_aKF", 30, -3, 3, "E_{T}^{miss} Eta");
    MSPlotCP["met_phi_aKF"] = new MultiSamplePlot(datasetsMSP, "met_phi_aKF", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
    MSPlotCP["met_corr_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "met_corr_pT_aKF", 40, 0, 400, "E_{T}^{miss} p_{T}", "GeV");
    MSPlotCP["met_corr_eta_aKF"] = new MultiSamplePlot(datasetsMSP, "met_corr_eta_aKF", 30, -3, 3, "E_{T}^{miss} Eta");
    MSPlotCP["met_corr_phi_aKF"] = new MultiSamplePlot(datasetsMSP, "met_corr_phi_aKF", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
    
    
    MSPlotCP["nJets"] = new MultiSamplePlot(datasetsMSP, "nJets", 13, -0.5, 12.5, "# jets");
    MSPlotCP["nBJets"] = new MultiSamplePlot(datasetsMSP, "nBJets", 9, -0.5, 8.5, "# b jets");
    MSPlotCP["CSVv2Discr_allJets"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_allJets", 48, 0.0, 1.2, "CSVv2 discriminant value");
    MSPlotCP["CSVv2Discr_leadingJet"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_leadingJet", 48, 0.0, 1.2, "CSVv2 discriminant value of leading jet");
    MSPlotCP["CSVv2Discr_jet2"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet2", 48, 0.0, 1.2, "CSVv2 discriminant value of jet2");
    MSPlotCP["CSVv2Discr_jet3"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet3", 48, 0.0, 1.2, "CSVv2 discriminant value of jet3");
    MSPlotCP["CSVv2Discr_jet4"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet4", 48, 0.0, 1.2, "CSVv2 discriminant value of jet4");
    MSPlotCP["CSVv2Discr_highest"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_highest", 48, 0.0, 1.2, "Highest CSVv2 discriminant value");
    MSPlotCP["CSVv2Discr_jetNb"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jetNb", 8, -0.5, 7.5, "Jet number (in order of decreasing p_{T}) with highest CSVv2 discriminant value");
    
    MSPlotCP["nJets_aKF"] = new MultiSamplePlot(datasetsMSP, "nJets_aKF", 13, -0.5, 12.5, "# jets");
    MSPlotCP["nBJets_aKF"] = new MultiSamplePlot(datasetsMSP, "nBJets_aKF", 9, -0.5, 8.5, "# b jets");
    MSPlotCP["CSVv2Discr_allJets_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_allJets_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value");
    MSPlotCP["CSVv2Discr_leadingJet_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_leadingJet_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value of leading jet");
    MSPlotCP["CSVv2Discr_jet2_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet2_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value of jet2");
    MSPlotCP["CSVv2Discr_jet3_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet3_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value of jet3");
    MSPlotCP["CSVv2Discr_jet4_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet4_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value of jet4");
    MSPlotCP["CSVv2Discr_highest_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_highest_aKF", 48, 0.0, 1.2, "Highest CSVv2 discriminant value");
    MSPlotCP["CSVv2Discr_jetNb_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jetNb_aKF", 8, -0.5, 7.5, "Jet number (in order of decreasing p_{T}) with highest CSVv2 discriminant value");
    
    MSPlotCP["M3"] = new MultiSamplePlot(datasetsMSP, "M3", 40, 60, 460, "M_{3}", "GeV");
    MSPlotCP["min_Mlb"] = new MultiSamplePlot(datasetsMSP, "min_Mlb", 40, 0, 400, "M_{lb}", "GeV");
    MSPlotCP["dR_Lep_B"] = new MultiSamplePlot(datasetsMSP, "dR_Lep_B", 25, 0, 5, "#Delta R(l,b)");
    
    MSPlotCP["M3_aKF"] = new MultiSamplePlot(datasetsMSP, "M3_aKF", 40, 60, 460, "M_{3}", "GeV");
    MSPlotCP["min_Mlb_aKF"] = new MultiSamplePlot(datasetsMSP, "min_Mlb_aKF", 40, 0, 400, "M_{lb}", "GeV");
    MSPlotCP["dR_Lep_B_aKF"] = new MultiSamplePlot(datasetsMSP, "dR_Lep_B_aKF", 25, 0, 5, "#Delta R(l,b)");
  }
  
  /// SFs
  MSPlot["scaleFactor"] = new MultiSamplePlot(datasetsMSP, "scaleFactor", 80, 0., 2., "SF");
  MSPlot["btag_SF"] = new MultiSamplePlot(datasetsMSP, "btag_SF", 80, 0., 2., "btag SF");
  MSPlot["pu_SF"] = new MultiSamplePlot(datasetsMSP, "pu_SF", 80, 0., 2., "pu SF");
  
  /// Reco
  MSPlot["W_mass"] = new MultiSamplePlot(datasetsMSP, "W mass before kinFitter", 50, 0, 200, "m_{W}", "GeV");
  MSPlot["top_mass"] = new MultiSamplePlot(datasetsMSP, "Top mass before kinFitter", 40, 0, 400, "m_{t}", "GeV");
  MSPlot["top_mass_zoom"] = new MultiSamplePlot(datasetsMSP, "Top mass before kinFitter (zoomed)", 40, 130, 210, "m_{t}", "GeV");
  MSPlot["red_top_mass_manyBins"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass (many bins)", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
  MSPlot["red_top_mass"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass", 44, 0.2, 2.4, "M_{t}/<M_{t}>");
  MSPlot["top_pT"] = new MultiSamplePlot(datasetsMSP, "Top pt before kinFitter", 80, 0, 400, "p_{T}", "GeV");
  
  MSPlot["mlb"] = new MultiSamplePlot(datasetsMSP, "mlb before kinFitter", 80, 0, 800, "m_{lb}", "GeV");
  MSPlot["ttbar_mass"] = new MultiSamplePlot(datasetsMSP, "ttbar mass before kinFitter", 50, 0, 1000, "m_{t#bar{t}}", "GeV");
  
  MSPlot["dR_lep_b_min"] = new MultiSamplePlot(datasetsMSP, "Minimum dR(lep,b) before kinFitter", 25, 0, 5, "#Delta R(l,b_{l})");
  MSPlot["dR_lep_b_max"] = new MultiSamplePlot(datasetsMSP, "Maximum dR(lep,b) before kinFitter", 25, 0, 5, "#Delta R(l,b_{h})");
  
  if (doKinFit)
  {
    MSPlot["scaleFactor_aKF"] = new MultiSamplePlot(datasetsMSP, "scaleFactor_aKF", 80, 0., 2., "SF");
    MSPlot["btag_SF_aKF"] = new MultiSamplePlot(datasetsMSP, "btag_SF_aKF", 80, 0., 2., "btag SF");
    MSPlot["pu_SF_aKF"] = new MultiSamplePlot(datasetsMSP, "pu_SF_aKF", 80, 0., 2., "pu SF");
    MSPlot["W_mass_aKF"] = new MultiSamplePlot(datasetsMSP, "W mass after kinFitter", 50, 0, 200, "m_{W,kf}", "GeV");
    MSPlot["top_mass_aKF"] = new MultiSamplePlot(datasetsMSP, "Top mass after kinFitter", 40, 0, 400, "m_{t,kf}", "GeV");
    MSPlot["top_mass_aKF_zoom"] = new MultiSamplePlot(datasetsMSP, "Top mass after kinFitter (zoomed)", 40, 130, 210, "m_{t,kf}", "GeV");
    MSPlot["red_top_mass_aKF_manyBins"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass after KF (many bins)", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
    MSPlot["red_top_mass_aKF"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass after KF", 44, 0.2, 2.4, "M_{t}/<M_{t}>");
    MSPlot["top_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "Top pt after kinFitter", 80, 0, 400, "p_{T}", "GeV");
    
    MSPlot["mlb_aKF"] = new MultiSamplePlot(datasetsMSP, "mlb after kinFitter", 80, 0, 800, "m_{lb}", "GeV");
    MSPlot["ttbar_mass_aKF"] = new MultiSamplePlot(datasetsMSP, "ttbar mass after kinFitter", 50, 0, 1000, "m_{t#bar{t}}", "GeV");
    
    MSPlot["dR_lep_b_min_aKF"] = new MultiSamplePlot(datasetsMSP, "Minimum dR(lep,b) after kinFitter", 25, 0, 5, "#Delta R(l,b_{l})");
    MSPlot["dR_lep_b_max_aKF"] = new MultiSamplePlot(datasetsMSP, "Maximum dR(lep,b) after kinFitter", 25, 0, 5, "#Delta R(l,b_{h})");
    
    MSPlot["KF_Chi2"] = new MultiSamplePlot(datasetsMSP, "Chi2 value of kinFitter", 50, 0, 5, "#chi^{2}");
    MSPlot["KF_Chi2_narrow"] = new MultiSamplePlot(datasetsMSP, "Chi2 value of kinFitter (zoomed)", 50, 0, 2, "#chi^{2}");
    MSPlot["KF_Chi2_wide"] = new MultiSamplePlot(datasetsMSP, "Chi2 value of kinFitter (large)", 50, 0, 20, "#chi^{2}");
  }
}

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  InitHisto1DMatch();
  
  histo1D["genTop_status"] = new TH1F("genTop_status", "Status code of the generated top quark; status", 100, -0.5, 99.5);
  histo1D["genTop_isLastCopy_status"] = new TH1F("genTop_isLastCopy_status", "Status code of the last copy of the generated top quark; status", 100, -0.5, 99.5);
  histo1D["genTop_hasLastCopy"] = new TH1F("genTop_hasLastCopy", "Generated top quark has isLastCopy tag; has isLastCopy", 2, -0.5, 1.5);
  histo1D["genTop_hasStatus22"] = new TH1F("genTop_hasStatus22", "Generated top quark has status 22; has status 22", 2, -0.5, 1.5);
  histo1D["genAntitop_status"] = new TH1F("genAntitop_status", "Status code of the generated antitop quark; status", 100, -0.5, 99.5);
  histo1D["genAntitop_isLastCopy_status"] = new TH1F("genAntitop_isLastCopy_status", "Status code of the last copy of the generated antitop quark; status", 100, -0.5, 99.5);
  histo1D["genAntitop_hasLastCopy"] = new TH1F("genAntitop_hasLastCopy", "Generated antitop quark has isLastCopy tag; has isLastCopy", 2, -0.5, 1.5);
  histo1D["genAntitop_hasStatus22"] = new TH1F("genAntitop_hasStatus22", "Generated antitop quark has status 22; has status 22", 2, -0.5, 1.5);
  
  histo1D["genTop_status22_pT"] = new TH1F("genTop_status22_pT", "pT of the generated top quark with status 22; p_{T} (GeV)", 250, 0., 250.);
  histo1D["genTop_status62_pT"] = new TH1F("genTop_status62_pT", "pT of the generated top quark with status 62; p_{T} (GeV)", 250, 0., 250.);
  histo1D["genAntitop_status22_pT"] = new TH1F("genAntitop_status22_pT", "pT of the generated antitop quark with status 22; p_{T} (GeV)", 250, 0., 250.);
  histo1D["genAntitop_status62_pT"] = new TH1F("genAntitop_status62_pT", "pT of the generated antitop quark with status 62; p_{T} (GeV)", 250, 0., 250.);
  
  /// SFs
  histo1D["width_SF"] = new TH1F("width_SF", "Scale factor to change the ttbar distribution width; width SF", 500, 0, 5);
  //histo1D["btag_SF"] = new TH1F("btag_SF", "b tag scale factor; btag SF", 80, 0, 2);
  
  if (makeControlPlots)
  {
    histo1D["dR_jets_min"]  = new TH1F("dR_jets_min","Minimal delta R between two jets; #Delta R(j_{1},j_{2})", 35, 0, 3.5);
//    histo1D["dR_jets_sum"]  = new TH1F("dR_jets_sum","Sum of delta R between all jets; #Delta R(j_{1},j_{2})", 60, 0, 12);
  }
  
  
  /// Systematic comparison
  histo1D["allSim_top_mass"] = new TH1F("allSim_top_mass","Reconstructed top mass for all simulated samples; m_{t}", 32, 130, 210);
  histo1D["allSim_red_top_mass"] = new TH1F("allSim_red_top_mass","Reduced top mass for all simulated samples; m_{t}/<m_{t}>", 110, 0.2, 2.4);
  
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
  
  InitHisto2DMatch();
  
  histo2D["puSF_vs_nTruePU"] = new TH2F("puSF_vs_nTruePU","nTruePU vs. pu SF (TT); pu SF; nTruePU", 80, 0, 2, 60, 0, 60);
  histo2D["puSF_vs_nVtx"] = new TH2F("puSF_vs_nVtx","nVtx vs. pu SF (TT); pu SF; nVtx", 80, 0, 2, 60, 0, 60);
  histo2D["nVtx_vs_nTruePU"] = new TH2F("nVtx_vs_nTruePU","nTruePU vs. nVtx (TT); nVtx; nTruePU", 60, 0, 60, 60, 0, 60);
  
  /// Reco
  histo2D["dR_lep_b_lep_vs_had_CM"] = new TH2F("dR_lep_b_lep_vs_had_CM","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, correct match); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_lep_b_lep_vs_had_WM"] = new TH2F("dR_lep_b_lep_vs_had_WM","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, wrong permutations); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_lep_b_lep_vs_had_NM"] = new TH2F("dR_lep_b_lep_vs_had_NM","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, no match); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  
  histo2D["ttbar_mass_vs_minMlb_CM"] = new TH2F("ttbar_mass_vs_minMlb_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_WM"] = new TH2F("ttbar_mass_vs_minMlb_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_NM"] = new TH2F("ttbar_mass_vs_minMlb_NM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCut_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCut_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCut_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCut_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCut_NM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCut_NM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCut_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCut_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCut_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCut_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCut_NM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCut_NM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCutHard_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCutHard_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_NM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCutHard_NM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCutHard_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCutHard_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_NM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCutHard_NM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCuts_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCuts_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCuts_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCuts_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCuts_NM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCuts_NM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCutsHard_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCutsHard_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_NM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCutsHard_NM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  
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

void InitHisto2DMatch()
{
  TH2::SetDefaultSumw2();
  
  /// Matched events
  histo2D["matched_mlb_corr_mlb_wrong"] = new TH2F("matched_mlb_corr_mlb_wrong","Wrongly constructed M_{lb} vs. correctly constructed M_{lb}; M_{lb_{lep}}; M_{lb_{had}}", 80, 0, 800, 80, 0, 800);
  histo2D["matched_dR_lep_b_corr_dR_lep_b_wrong"] = new TH2F("matched_dR_lep_b_corr_dR_lep_b_wrong","Wrongly constructed dR(l,b) vs. correctly constructed dR(l,b); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["matched_mlb_dR_lep_b_corr"] = new TH2F("matched_mlb_dR_lep_b_corr","dR(l,b) vs. M_{lb}; M_{lb_{lep}}; #Delta R(l,b_{lep})", 80, 0, 800, 25, 0, 5);
  histo2D["matched_mlb_dR_lep_b_wrong"] = new TH2F("matched_mlb_dR_lep_b_wrong","dR(l,b) vs. M_{lb}, both wrongly matched; M_{lb_{had}}; #Delta R(l,b_{had})", 80, 0, 800, 25, 0, 5);
}

void InitLikelihoodPlots()
{
  for (int iMass = 0; iMass < nLikeMasses; iMass++)
  {
    histo1DLike["loglike_vs_width_m"+DotReplace(redTopMassArray[iMass])] = new TH1D(("loglike_vs_width_m"+DotReplace(redTopMassArray[iMass])).c_str(), ("loglike_vs_width_m"+DotReplace(redTopMassArray[iMass])+"; width (a.u.); -log(likelihood)").c_str(), nWidthsLike, -0.5, nWidthsLike-0.5);
//     histo1DLike["loglike_vs_width_CM_m"+DotReplace(redTopMassArray[iMass])] = new TH1D(("loglike_vs_width_CM_m"+DotReplace(redTopMassArray[iMass])).c_str(), ("loglike_vs_width_CM_m"+DotReplace(redTopMassArray[iMass])+"; width (a.u.); -log(likelihood)").c_str(), nWidthsLike, -0.5, nWidthsLike-0.5);
//     histo1DLike["loglike_vs_width_WM_m"+DotReplace(redTopMassArray[iMass])] = new TH1D(("loglike_vs_width_WM_m"+DotReplace(redTopMassArray[iMass])).c_str(), ("loglike_vs_width_WM_m"+DotReplace(redTopMassArray[iMass])+"; width (a.u.); -log(likelihood)").c_str(), nWidthsLike, -0.5, nWidthsLike-0.5);
//     histo1DLike["loglike_vs_width_NM_m"+DotReplace(redTopMassArray[iMass])] = new TH1D(("loglike_vs_width_NM_m"+DotReplace(redTopMassArray[iMass])).c_str(), ("loglike_vs_width_NM_m"+DotReplace(redTopMassArray[iMass])+"; width (a.u.); -log(likelihood)").c_str(), nWidthsLike, -0.5, nWidthsLike-0.5);
  }
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
      if (isTTbar)
      {
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
      }
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
    nofTTEventsWithoutAGenTop = 0;
    nofTTEventsWithoutGenTop = 0;
    nofTTEventsWithoutGenAntiTop = 0;
    nofTTEventsWithoutBothGenTopsWithStatus22 = 0;
    nofTTEventsWithoutGenTopWithStatus22 = 0;
    nofTTEventsWithoutGenAntiTopWithStatus22 = 0;
    nofTTEventsWithoutBothGenTopsWithStatus62 = 0;
    nofTTEventsWithoutGenTopWithStatus62 = 0;
    nofTTEventsWithoutGenAntiTopWithStatus62 = 0;
    sumWeight1001 = 0;
    sumWeight1002 = 0;
    sumWeight1003 = 0;
    sumWeight1004 = 0;
    sumWeight1005 = 0;
    sumWeight1007 = 0;
    sumWeight1009 = 0;
  }
  
  strSyst = "";
  nEventsDataSet = 0;
  xSection = 1.;
  eqLumi = 1.;
  lumiWeight = 1.;
  relativeSF = 1.;
  
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
  widthSF = 1.;
  thisLeptonSF = 1.;
  thisLeptonIdSF = 1.;
  thisLeptonIsoSF = 1.;
  thisLeptonTrigSF = 1.;
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
  massHadTopQ = 0.01;
  massLepTopQ = 0.01;
  catSuffix = "";
  isCM = false;
  isWM = false;
  isNM = false;
  doneKinFit = false;
  kFitVerbosity = false;
  kFitChi2 = 99.;
  kFitChi2Matched = 99.;
  passKFChi2MatchedCut = false;
  //toyValue = -1.;
  for (int i = 0; i < nPsExps; i++)
  {
    toyValues[i] = 1.;
  }
  
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
  loglike_per_evt.clear();
}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearVars();
}

void FillControlPlots(vector<Dataset *> datasets, int d)
{
  FillControlPlots(datasets, d, "");
}

void FillControlPlots(vector<Dataset *> datasets, int d, string suffix)
{
  M3 = (selectedJets[0] + selectedJets[1] + selectedJets[2]).M();
  Ht = selectedJets[0].Pt() + selectedJets[1].Pt() + selectedJets[2].Pt() + selectedJets[3].Pt();
  
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
  
  MSPlotCP["muon_pT"+suffix]->Fill(selectedLepton[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["muon_eta"+suffix]->Fill(selectedLepton[0].Eta(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["muon_phi"+suffix]->Fill(selectedLepton[0].Phi(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["muon_relIso"+suffix]->Fill(muon_relIso[0], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["muon_d0"+suffix]->Fill(muon_d0[0], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  if (suffix.find("aKF") != std::string::npos )
  {
    MSPlotCP["leadingJet_pT"+suffix]->Fill(selectedJetsAKF[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet2_pT"+suffix]->Fill(selectedJetsAKF[1].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet3_pT"+suffix]->Fill(selectedJetsAKF[2].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet4_pT"+suffix]->Fill(selectedJetsAKF[3].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    Ht_aKF = selectedJetsAKF[0].Pt() + selectedJetsAKF[1].Pt() + selectedJetsAKF[2].Pt() + selectedJetsAKF[3].Pt();
    MSPlotCP["Ht_4leadingJets"+suffix]->Fill(Ht_aKF, datasets[d], true, lumiWeight*scaleFactor*widthSF);
    M3_aKF = (selectedJetsAKF[0] + selectedJetsAKF[1] + selectedJetsAKF[2]).M();
    MSPlotCP["M3"+suffix]->Fill(M3_aKF, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  }
  else
  {
    MSPlotCP["leadingJet_pT"+suffix]->Fill(selectedJets[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet2_pT"+suffix]->Fill(selectedJets[1].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet3_pT"+suffix]->Fill(selectedJets[2].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet4_pT"+suffix]->Fill(selectedJets[3].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["Ht_4leadingJets"+suffix]->Fill(Ht, datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["M3"+suffix]->Fill(M3, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  }
  MSPlotCP["met_pT"+suffix]->Fill(met_pt, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_eta"+suffix]->Fill(met_eta, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_phi"+suffix]->Fill(met_phi, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_corr_pT"+suffix]->Fill(met_corr_pt, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_corr_eta"+suffix]->Fill(met_corr_eta, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_corr_phi"+suffix]->Fill(met_corr_phi, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  MSPlotCP["min_Mlb"+suffix]->Fill(min_Mlb, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["dR_Lep_B"+suffix]->Fill(dRLepB, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  MSPlotCP["nJets"+suffix]->Fill(selectedJets.size(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["nBJets"+suffix]->Fill(selectedBJets.size(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  MSPlotCP["CSVv2Discr_leadingJet"+suffix]->Fill(jet_bdiscr[0], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["CSVv2Discr_jet2"+suffix]->Fill(jet_bdiscr[1], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["CSVv2Discr_jet3"+suffix]->Fill(jet_bdiscr[2], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["CSVv2Discr_jet4"+suffix]->Fill(jet_bdiscr[3], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  int labelB = -1;
  double highestBDiscr = -999.;
  for (int iJet = 0; iJet < selectedJets.size(); iJet++)
  {
    if ( suffix.find("aKF") != std::string::npos )
      MSPlotCP["jet_pT_allJets"+suffix]->Fill(selectedJetsAKF[iJet].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    else
      MSPlotCP["jet_pT_allJets"+suffix]->Fill(selectedJets[iJet].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    
    MSPlotCP["CSVv2Discr_allJets"+suffix]->Fill(jet_bdiscr[iJet], datasets[d], true, lumiWeight*scaleFactor*widthSF);
    if ( jet_bdiscr[iJet] > highestBDiscr )
    {
      highestBDiscr = jet_bdiscr[iJet];
      labelB = iJet;
    }
  }
  MSPlotCP["CSVv2Discr_highest"+suffix]->Fill(highestBDiscr, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["CSVv2Discr_jetNb"+suffix]->Fill(labelB, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  if (isTTbar)
  {
    pair<unsigned int,unsigned int> hasMinDR = pair<unsigned int,unsigned int>(9999,9999);
    double tempMin = 999.;
    double tempDR;
    for (unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
    {
      for (unsigned int jJet = iJet+1; jJet < selectedJets.size(); jJet++)
      {
        tempDR = ROOT::Math::VectorUtil::DeltaR( selectedJets[iJet], selectedJets[jJet]);
        if ( tempDR < tempMin)
        {
          tempMin = tempDR;
          hasMinDR = pair<unsigned int,unsigned int>(iJet,jJet);
        }
      }
    }
    tempDR = ROOT::Math::VectorUtil::DeltaR( selectedJets[hasMinDR.first], selectedJets[hasMinDR.second]);
    histo1D["dR_jets_min"]->Fill(tempDR);
  }
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
  if (isTTbar)
  {
    histo2D["puSF_vs_nTruePU"]->Fill(puSF,npu);
    histo2D["puSF_vs_nVtx"]->Fill(puSF,nvtx);
    histo2D["nVtx_vs_nTruePU"]->Fill(nvtx,npu);
  }
  
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

void FillMSPlots(int d, bool doneKinFit)
{
  string suffix = "";
  if (doneKinFit) suffix = "_aKF";
  
  if (makeControlPlots) FillControlPlots(datasetsMSP, d, suffix);
  
  if (! isData)
  {
    MSPlot["scaleFactor"+suffix]->Fill(scaleFactor, datasetsMSP[d], false, lumiWeight*scaleFactor*widthSF);
    MSPlot["btag_SF"+suffix]->Fill(btagSF, datasetsMSP[d], false, lumiWeight*scaleFactor*widthSF);
    MSPlot["pu_SF"+suffix]->Fill(puSF, datasetsMSP[d], false, lumiWeight*scaleFactor*widthSF);
  }
  
  if (! doneKinFit)
  {
    MSPlot["W_mass"]->Fill(reco_W_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_mass"]->Fill(reco_top_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    if ( reco_top_mass_bKF < 210 && reco_top_mass_bKF > 130 )
    {
      MSPlot["top_mass_zoom"]->Fill(reco_top_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    if ( redTopMass_bKF < 2.4 )
    {
      MSPlot["red_top_mass_manyBins"]->Fill(redTopMass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      MSPlot["red_top_mass"]->Fill(redTopMass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    MSPlot["top_pT"]->Fill(reco_top_pt_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["mlb"]->Fill(reco_mlb_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["ttbar_mass"]->Fill(reco_ttbar_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_lep_b_min"]->Fill(reco_dRLepB_lep_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_lep_b_max"]->Fill(reco_dRLepB_had_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
  }
  else if (doKinFit)
  {
    if ( kFitChi2 < 5 )
    {
      MSPlot["KF_Chi2"]->Fill(kFitChi2, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      if ( kFitChi2 < 2 )
        MSPlot["KF_Chi2_narrow"]->Fill(kFitChi2, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    
    MSPlot["W_mass"+suffix]->Fill(reco_W_mass_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_mass"+suffix]->Fill(reco_top_mass_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    if ( reco_top_mass_aKF < 210 && reco_top_mass_aKF > 130 )
    {
      MSPlot["top_mass"+suffix+"_zoom"]->Fill(reco_top_mass_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    if ( redTopMass < 2.4 )
    {
      MSPlot["red_top_mass"+suffix+"_manyBins"]->Fill(redTopMass, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      MSPlot["red_top_mass"+suffix]->Fill(redTopMass, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    MSPlot["top_pT"+suffix]->Fill(reco_top_pt_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["mlb"+suffix]->Fill(reco_mlb_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["ttbar_mass"+suffix]->Fill(reco_ttbar_mass_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_lep_b_min"+suffix]->Fill(reco_dRLepB_lep_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_lep_b_max"+suffix]->Fill(reco_dRLepB_had_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
  }
}

void FillLikelihoodPlots()
{
  for (int iMass = 0; iMass < nLikeMasses; iMass++)
  {
    loglike_per_evt.clear();
    if ( redTopMassArray[iMass] > minCutRedTopMass && maxCutRedTopMass > redTopMassArray[iMass])
    {
      loglike_per_evt = like->CalculateLikelihood(redTopMassArray[iMass], 1., 172.5, 172.5, 1., false, false);
      for (int iWidth = 0; iWidth < nWidthsLike; iWidth++)  // Deliberately filling iWidth and not actual width
      {
//        if (redTopMass > redTopMassArray[iMass] && redTopMass < redTopMassArray[iMass+1])
//        {
          histo1DLike["loglike_vs_width_m"+DotReplace(redTopMassArray[iMass])]->SetBinContent(iWidth, loglike_per_evt[iWidth]);
//          if (isCM) histo1DLike["loglike_vs_width_CM_m"+DotReplace(redTopMassArray[iMass])]->SetBinContent(iWidth, loglike_per_evt[iWidth]);
//          if (isWM) histo1DLike["loglike_vs_width_WM_m"+DotReplace(redTopMassArray[iMass])]->SetBinContent(iWidth, loglike_per_evt[iWidth]);
//          if (isNM) histo1DLike["loglike_vs_width_NM_m"+DotReplace(redTopMassArray[iMass])]->SetBinContent(iWidth, loglike_per_evt[iWidth]);
//        }
      }  // end width

//      break;  // only fill once per event (because only one redTopMass)
    }
  }  // end masses
}

void WriteLikelihoodPlots()
{
  FillLikelihoodPlots();
  
  mkdir((pathOutput+"Likelihood/").c_str(),0777);
  string fileName = pathOutput+"Likelihood/LikelihoodPlots.root";
  cout << " - Recreate output file ..." << endl;
  TFile *foutLike = new TFile (fileName.c_str(), "RECREATE");
  cout << "   Output file is " << fileName << endl;
  foutLike->cd();
  
  //gStyle->SetOptStat(1111);
  for (std::map<std::string,TH1D*>::const_iterator it = histo1DLike.begin(); it != histo1DLike.end(); it++)
  {
    TH1D *temp = it->second;
    for (int iBin = 0; iBin < temp->GetNbinsX(); iBin++)
    {
      temp->GetXaxis()->SetBinLabel(iBin+1, (DotReplace(widthsLike[iBin])).c_str());
    }
    temp->SetStats(0);
    temp->Write();
//     TCanvas* tempCanvas = new TCanvas((it->first).c_str(), (it->first).c_str());
//     tempCanvas->cd();
//     temp->Draw();
//     tempCanvas->SaveAs( (pathOutput+"Likelihood/"+it->first+".png").c_str() );
  }
  
  foutLike->Close();
  delete foutLike;
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

