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
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
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
using namespace reweight;
using namespace TopTree;


bool test = false;
bool testHistos = false;
bool testTTbarOnly = false;
bool makePlots = true;
bool doKinFit = true;
bool applyKinFitCut = true;
double kinFitCutValue = 5.;

bool doMETCleaning = true;
bool applyLeptonSF = true;
bool applyPU = true;
bool applyBTagSF = true;




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
int verbose = 2;

string pathNtuples = "";
string pathNtuplesMC = "";
string pathNtuplesData = "";
bool isData = false;
bool isTTbar = false;

int nofHardSelected = 0;
int nofMETCleaned = 0;

// Pile-up
string pathCalPileup = "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/";
LumiReWeighting LumiWeights = LumiReWeighting(pathCalPileup+"MCPileup_Summer16.root", pathCalPileup+"pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet.root", "pileup", "pileup");
LumiReWeighting LumiWeights_up = LumiReWeighting(pathCalPileup+"MCPileup_Summer16.root", pathCalPileup+"pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysPlus.root", "pileup", "pileup");
LumiReWeighting LumiWeights_down = LumiReWeighting(pathCalPileup+"MCPileup_Summer16.root", pathCalPileup+"pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysMinus.root", "pileup", "pileup");


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
std::array<double, 14> aveTopMass = {171.826, 169.746, 167.511, 197.053, 196.687, 197.911, 181.895, 249.468, 247.437, 227.530, 226.099, 184.794, 184.594, 184.624};  // Res 170608 Single Gaus
//std::array<double, 14> aveTopMass = {171.826, 169.746, 167.572, 196.603, 196.072, 197.919, 181.953, 247.003, 243.879, 226.505, 224.951, 184.717, 184.598, 184.616};  // Res 170515
//  no KF chi2 cut
//std::array<double, nofAveMasses> aveTopMass = {171.810, 168.728, 167.110, 203.721, 204.952, 198.233, 193.403, 270.895, 267.167, 230.144, 229.649, 250.010, 242.091, 200.455, 193.963, 194.025};

/// # events after kin fitter
//  KF chi2 < 5
int nEventsAKF[] = {100949, 547108, 8927, 8847, 7051, 4283, 3, 16, 55, 276, 1, 18, 123, 1415};

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
void GetMetaData(TTree* tree, bool isData);
void InitTree(TTree* tree, bool isData);
void InitHisto1D();
void InitHisto2D();
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearVars();
void ClearObjects();
void FillCatsPlots();
long GetNEvents(TTree* fChain, string var, bool isData);
void GetEraFraction(double* fractions);



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


long nEventsDataSet;
double xSection;
double lumiWeight, scaleFactor, widthSF;
double thisLeptonSF, thisLeptonIdSF, thisLeptonIsoSF, thisLeptonTrigSF;
vector<unsigned int> bJetId;
double bdiscrTop, bdiscrTop2, tempbdiscr;
int labelB1, labelB2;
int labelsReco[4];


/// Define TLVs
TLorentzVector muon, jet, mcpart;
TLorentzVector WCandidate;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
vector<TLorentzVector> selectedJetsAKF;
vector<TLorentzVector> selectedJetsKFcorrected;
vector<TLorentzVector> mcParticles;
vector<TLorentzVector> partons;


/// KinFitter
bool doneKinFit = false;
TKinFitter* kFitter;
bool addWMassKF = true;
bool addEqMassKF = false;
int kFitVerbosity = 0;
double kFitChi2 = 99.;
int nofAcceptedKFit = 0;



/// Variables
double min_Mlb, dRLepB;
double M3_aKF, Ht_aKF;
double reco_W_mass_aKF, reco_top_mass_aKF, reco_top_pt_aKF, reco_mlb_aKF, reco_dRLepB_lep_aKF, reco_dRLepB_had_aKF, reco_ttbar_mass_aKF, redTopMass;

double puSF_nTruePU, puSF_nTruePU_up, puSF_nTruePU_down, puSF_nVtx, puSF_nVtx_up, puSF_nVtx_down;


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
  
  if (test) makePlots = false;
  if (testHistos)
  {
    makePlots = true;
  }
  
  string pathOutput = "OutputPlots/";
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
  
  
  pathNtuplesMC = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.first+"/";
  pathNtuplesData = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.second+"/";
  cout << "Using Ntuples from " << ntupleDate.first << " for MC and " << ntupleDate.second << " for data. This corresponds to systematics: " << systStr << endl;
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
    
    if ( dataSetName.find("TT") != std::string::npos )
    {
      datasets[d]->SetTitle("t#bar{t}");
      datasets[d]->SetColor(kRed+1);
    }
    //if ( dataSetName.find("TTbarJets_Other") != std::string::npos ) datasets[d]->SetColor(kRed-7);
    if ( dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      //datasets[d]->SetColor(kGreen-3);
      datasets[d]->SetColor(kBlue-2);
    }
    if ( dataSetName.find("ZJets") != std::string::npos || dataSetName.find("DY") != std::string::npos )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{#font[122]{\55}}");
      //datasets[d]->SetColor(kAzure-2);
      //datasets[d]->SetColor(kMagenta);
      datasets[d]->SetColor(kAzure+6);
    }
    if ( dataSetName.find("ST") != std::string::npos || dataSetName.find("SingleTop") != std::string::npos )
    {
      datasets[d]->SetTitle("ST");
      //datasets[d]->SetColor(kBlue-2);
      datasets[d]->SetColor(kOrange-4);  // 595, 615, 800
    }
  }
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  KinFitter *kf = new KinFitter("PlotsForResolutionFunctions_testFit_170608_S.root", addWMassKF, addEqMassKF);
  
  if (makePlots)
  {
    InitHisto1D();
    InitHisto2D();
  }
  
  
  
  
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
  
  /// No Loop over systematics or widths
  for (int iSys = 0; iSys < 1; iSys++)
  {
   
    hasFoundTTbar = false;
    
    
    ////////////////////////////////////
    /// Loop over datasets
    ////////////////////////////////////
    
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
        
        lumiWeight = Luminosity/eqLumi;
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
          
          if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
          if (applyBTagSF) { scaleFactor *= btagSF;}
          if (applyPU) { scaleFactor *= puSF;}
          
          
          puSF_nTruePU = LumiWeights.ITweight( npu );
          puSF_nTruePU_up = LumiWeights_up.ITweight( npu );
          puSF_nTruePU_down = LumiWeights_down.ITweight( npu );
          puSF_nVtx = LumiWeights.ITweight( nvtx );
          puSF_nVtx_up = LumiWeights_up.ITweight( nvtx );
          puSF_nVtx_down = LumiWeights_down.ITweight( nvtx );
          
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
          kFitter = kf->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]], kFitVerbosity);
          
          if ( kFitter->getStatus() != 0 )  // did not converge
          {
            if (test && verbose > 2) cout << "Event " << ievt << ": Fit did not converge..." << endl;
            continue;
          }
          
          kFitChi2 = kFitter->getS();
          if (test && verbose > 4) cout << "Fit converged: Chi2 = " << kFitChi2 << endl;
          
          doneKinFit = true;
          
          
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
        redTopMass = reco_top_mass_aKF/aveTopMass[2];
        
        
        
        //Fill histos
        if ( makePlots )
        {
          
          FillCatsPlots();
          
        }  // end makePlots
        
        
        
        
      }  // end loop events
      
      
      cout << endl;  /// Stronger selection in this analyser compared to Ntuples ==> endEvent --> nofHardSelected
      cout << "Number of events with exactly 4 jets with pT > 30 GeV: " << nofHardSelected << " (" << 100*((float)nofHardSelected/(float)endEvent) << "%)" << endl;
      cout << "Number of events with clean MET: " << nofMETCleaned << " (" << 100*((float)nofMETCleaned/(float)nofHardSelected) << "%)" << endl;
      if (doKinFit) cout << "Number of clean events accepted by kinFitter: " << nofAcceptedKFit << " (" << 100*((float)nofAcceptedKFit/(float)nofMETCleaned) << "%)" << endl;
      
      
      tFileMap[dataSetName.c_str()]->Close();
      
      timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
      
    }  // end loop datasets
    
    
    
    cout << endl << "Processing time per dataset: " << endl;
    for (unsigned int d = 0; d < datasets.size(); d++)
    {
      cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
    }
    
    
    
  }  // end loop systematics/widths
  
  
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
  if (isTTbar)
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
  if (! isData)
  {
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
}

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  /// SFs
  histo1D["width_SF"] = new TH1F("width_SF", "Scale factor to change the ttbar distribution width; width SF", 500, 0, 5);
  histo1D["btag_SF"] = new TH1F("btag_SF", "b tag scale factor; btag SF", 80, 0, 2);
  

}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
  histo2D["puSF_vs_nTruePU"] = new TH2F("puSF_vs_nTruePU","nTruePU vs. pu SF (TT); pu SF; nTruePU", 80, 0, 2, 60, 0, 60);
  histo2D["puSF_vs_nVtx"] = new TH2F("puSF_vs_nVtx","nVtx vs. pu SF (TT); pu SF; nVtx", 80, 0, 2, 60, 0, 60);
  histo2D["nVtx_vs_nTruePU"] = new TH2F("nVtx_vs_nTruePU","nTruePU vs. nVtx (TT); nVtx; nTruePU", 60, 0, 60, 60, 0, 60);
  
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
  if (! isData)
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
  }
  
  strSyst = "";
  nEventsDataSet = 0;
  xSection = 1.;
  eqLumi = 1.;
  lumiWeight = 1.;
  
  nofHardSelected = 0;
  nofMETCleaned = 0;
  nofAcceptedKFit = 0;

  

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
  if (! isData)
  {
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
}

void ClearTLVs()
{
  muon.Clear();
  jet.Clear();
  WCandidate.Clear();
  selectedLepton.clear();
  selectedJets.clear();
  selectedBJets.clear();
  selectedJetsAKF.clear();
  selectedJetsKFcorrected.clear();
}


void ClearVars()
{
  scaleFactor = 1.;
  widthSF = 1.;
  thisLeptonSF = 1.;
  thisLeptonIdSF = 1.;
  thisLeptonIsoSF = 1.;
  thisLeptonTrigSF = 1.;
  puSF_nTruePU = 1.;
  puSF_nTruePU_up = 1.;
  puSF_nTruePU_down = 1.;
  puSF_nVtx = 1.;
  puSF_nVtx_up = 1.;
  puSF_nVtx_down = 1.;
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
  doneKinFit = false;
  kFitVerbosity = false;
  kFitChi2 = 99.;
  
  M3_aKF = -1.;
  Ht_aKF = -1;
  min_Mlb = 9999.;
  dRLepB = -1.;
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

void FillCatsPlots()
{
  if (isTTbar)
  {
    histo2D["puSF_vs_nTruePU"]->Fill(puSF_nTruePU,npu);
    histo2D["puSF_vs_nVtx"]->Fill(puSF_nVtx,nvtx);
    histo2D["nVtx_vs_nTruePU"]->Fill(nvtx,npu);
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
