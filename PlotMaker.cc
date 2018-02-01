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
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
//#include "../macros/Style.C"





using namespace std;
using namespace TopTree;

bool runLocally = false;
bool runFromPnfs = false;

bool test = false;
bool testHistos = false;
bool testTTbarOnly = false;
bool skipData = false;
bool makePlots = true;


bool applyLeptonSF = true;
bool applyPU = true;
bool applyBTagSF = true;

bool runTTbar = true;
bool runSTtW = true;
bool runSTt = true;
bool runOther = true;


bool runGenWidth = false;
double listGenWidths[] = {0.2, 0.5, /*0.8, 2.,*/ 4./*, 8.*/};
int nGenWidths = sizeof(listGenWidths)/sizeof(listGenWidths[0]);


string suffixes[] = {"_mb", "", "_2b", "_2b_4j"};
int nSuffixes = sizeof(suffixes)/sizeof(suffixes[0]);


bool newTrees = true;
string systStr = "nominal";
pair<string,string> whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos )
  {
    if (newTrees)
      return pair<string,string>("171125","171125");
    else
      return pair<string,string>("171015","171021");
  }
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    return pair<string,string>("171125","171125");
  }
}
pair<string,string> ntupleDate = whichDate(systStr);
int verbose = 2;

string pathNtuples = "";
string pathNtuplesMC = "";
string pathNtuplesData = "";
string pathOutput = "";


bool isData = false;
bool isTTbar = false;
bool isTTsemilep = false;
bool isTTother = false;
bool isST = false;
bool isSTtW = false;
bool isSTt = false;
bool isOther = false;

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
map<string,TH2F*> histo2DT;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlotCP;
map<string,MultiSamplePlot*> MSPlotT;
map<string,TH1D*> histo1DLike;
map<string,TGraph*> graph;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;

vector < Dataset* > datasets;

double sumTopMass, sumMlb, sumEvents, sumCMTopMass, sumCMMlb, sumCMhadEvents, sumCMlepEvents;
Double_t aveMlbMass = 97.1331, aveMlbMassCM = 97.3438;  //180110_1159


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
string ConvertIntToString(int nb, int pad = 0);
string MakeTimeStamp();
bool fexists(const char *filename);
void GetMetaData(TTree* tree, bool isData);
void InitTree(TTree* tree, bool isData);
void InitMSPlots();
void InitHisto1D();
void InitHisto2D();
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearVars();
void ClearObjects();
void FillPlots(int d, string suffix);
long GetNEvents(TTree* fChain, string var, bool isData);
long GetNEvents(TTree* fChain, string var, unsigned int index, bool isData);
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
Bool_t          isDataRunB;
Bool_t          isDataRunC;
Bool_t          isDataRunD;
Bool_t          isDataRunE;
Bool_t          isDataRunF;
Bool_t          isDataRunG;
Bool_t          isDataRunH;
Int_t           nMuons;
Double_t        muon_charge[1];   //[nMuons]
Double_t        muon_pt[1];   //[nMuons]
Double_t        muon_phi[1];   //[nMuons]
Double_t        muon_eta[1];   //[nMuons]
Double_t        muon_E[1];   //[nMuons]
Double_t        muon_M[1];   //[nMuons]
Double_t        muon_d0[1];   //[nMuons]
Double_t        muon_dz[1];   //[nMuons]
Double_t        muon_chargedHadronIso[1];   //[nMuons]
Double_t        muon_neutralHadronIso[1];   //[nMuons]
Double_t        muon_photonIso[1];   //[nMuons]
Double_t        muon_puChargedHadronIso[1];   //[nMuons]
Double_t        muon_relIso[1];   //[nMuons]
Double_t        muon_pfIso[1];   //[nMuons]
Double_t        muon_chi2[1];   //[nMuons]
Int_t           muon_nofValidHits[1];   //[nMuons]
Int_t           muon_nofValidMuHits[1];   //[nMuons]
Int_t           muon_nofValidPixelHits[1];   //[nMuons]
Int_t           muon_nofMatchedStations[1];   //[nMuons]
Int_t           muon_nofTrackerLayersWithMeasurement[1];   //[nMuons]
Int_t           nJets;
Int_t           jet_nConstituents[20];   //[nJets]
Int_t           jet_nChConstituents[20];   //[nJets]
Double_t        jet_charge[20];   //[nJets]
Double_t        jet_pt[20];   //[nJets]
Double_t        jet_phi[20];   //[nJets]
Double_t        jet_eta[20];   //[nJets]
Double_t        jet_E[20];   //[nJets]
Double_t        jet_M[20];   //[nJets]
Double_t        jet_bdiscr[20];   //[nJets]
Double_t        jet_hadronFlavour[20];   //[nJets]
Int_t           nGenJets;
Int_t           genJet_pdgId[40];   //[nGenJets]
Double_t        genJet_charge[40];   //[nGenJets]
Double_t        genJet_pt[40];   //[nGenJets]
Double_t        genJet_phi[40];   //[nGenJets]
Double_t        genJet_eta[40];   //[nGenJets]
Double_t        genJet_E[40];   //[nGenJets]
Double_t        genJet_M[40];   //[nGenJets]
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
Bool_t          hasGenAntiTop;
Double_t        btagSF;
Double_t        puSF;
Double_t        muonIdSF_BCDEF[1];   //[nMuons]
Double_t        muonIdSF_GH[1];   //[nMuons]
Double_t        muonIsoSF_BCDEF[1];   //[nMuons]
Double_t        muonIsoSF_GH[1];   //[nMuons]
Double_t        muonTrigSF_BCDEF[1];   //[nMuons]
Double_t        muonTrigSF_GH[1];   //[nMuons]
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
Long64_t        nofEventsWithGenAntiTop;
Long64_t        nofTTEventsWithoutAGenTop;

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
TBranch        *b_muon_dz;   //!
TBranch        *b_muon_chargedHadronIso;   //!
TBranch        *b_muon_neutralHadronIso;   //!
TBranch        *b_muon_photonIso;   //!
TBranch        *b_muon_puChargedHadronIso;   //!
TBranch        *b_muon_relIso;   //!
TBranch        *b_muon_pfIso;   //!
TBranch        *b_muon_chi2;   //!
TBranch        *b_muon_nofValidHits;   //!
TBranch        *b_muon_nofValidMuHits;   //!
TBranch        *b_muon_nofValidPixelHits;   //!
TBranch        *b_muon_nofMatchedStations;   //!
TBranch        *b_muon_nofTrackerLayersWithMeasurement;   //!
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
TBranch        *b_jet_hadronFlavour;   //!
TBranch        *b_nGenJets;   //!
TBranch        *b_genJet_pdgId;   //!
TBranch        *b_genJet_charge;   //!
TBranch        *b_genJet_pt;   //!
TBranch        *b_genJet_phi;   //!
TBranch        *b_genJet_eta;   //!
TBranch        *b_genJet_E;   //!
TBranch        *b_genJet_M;   //!
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
TBranch        *b_hasGenAntiTop;   //!
TBranch        *b_btagSF;   //!
TBranch        *b_puSF;   //!
TBranch        *b_muonIdSF_BCDEF;   //!
TBranch        *b_muonIdSF_GH;   //!
TBranch        *b_muonIsoSF_BCDEF;   //!
TBranch        *b_muonIsoSF_GH;   //!
TBranch        *b_muonTrigSF_BCDEF;   //!
TBranch        *b_muonTrigSF_GH;   //!
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
TBranch        *b_nofEventsWithGenAntiTop;   //!
TBranch        *b_nofTTEventsWithoutAGenTop;   //!


long nEventsDataSet;
double xSection;
double lumiWeight, scaleFactor;
double thisLeptonSF, thisLeptonIdSF, thisLeptonIsoSF, thisLeptonTrigSF;



/// Define TLVs
TLorentzVector muon, jet;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
vector<unsigned int> bJetId;

double tempDR;



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

string dateString = "";

int main(int argc, char* argv[])
{
  dateString = MakeTimeStamp();
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
  cout << "*   - MET Cleaning                          *" << endl;
  cout << "*********************************************" << endl;
  cout << "* The following scale factors are applied:  *" << endl;
  if (applyLeptonSF) cout << "*   - Lepton scale factors: nominal         *" << endl;
  if (applyPU)       cout << "*   - Pile up: nominal                      *" << endl;
  if (applyBTagSF)   cout << "*   - B tag scale factors: nominal          *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  string channel = "mu";
  
  
  
  if (testHistos) makePlots = true;
    
  pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (makePlots)
  {
    pathOutput += "general/";
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
  
  
  ntupleDate = whichDate(systStr);
  if (runLocally)
  {
    pathNtuplesMC   = "/Volumes/LaCie/Ntuples/"+ntupleDate.first+"/";
    pathNtuplesData = "/Volumes/LaCie/Ntuples/"+ntupleDate.second+"/";
  }
  else if (runFromPnfs)
  {
    pathNtuplesMC   = "/pnfs/iihe/cms/store/user/lmoreels/Ntuples/CMSSW_80X/"+ntupleDate.first+"/";
    pathNtuplesData = "/pnfs/iihe/cms/store/user/lmoreels/Ntuples/CMSSW_80X/"+ntupleDate.second+"/";
  }
  else
  {
    pathNtuplesMC   = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.first+"/";
    pathNtuplesData = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.second+"/";
  }
  cout << "Using Ntuples from " << ntupleDate.first << " for MC and " << ntupleDate.second << " for data. This corresponds to systematics: " << systStr << endl;
  
  if (testHistos) cout << "Testing histogram consistency..." << endl;
  
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
  
  int dTT = -1;
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    string dataSetName = datasets[d]->Name();
    
    
    if ( dataSetName.find("QCD") != std::string::npos )
    {
      datasets[d]->SetColor(kYellow);
    }
    if ( dataSetName.find("TT") != std::string::npos )
    {
      datasets[d]->SetTitle("t#bar{t}");
      datasets[d]->SetColor(kRed+1);
    }
    //if ( dataSetName.find("TTbarJets_Other") != std::string::npos ) datasets[d]->SetColor(kRed-7);
    if ( dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      //datasets[d]->SetTitle("Other");
      datasets[d]->SetColor(kGreen+3);  // kGreen+2
      //datasets[d]->SetColor(kBlue-2);
    }
    if ( dataSetName.find("ZJets") != std::string::npos || dataSetName.find("DY") != std::string::npos )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{#font[122]{\55}}");
      //datasets[d]->SetTitle("Other");
      datasets[d]->SetColor(kGreen-3);  // kGreen+2
      //datasets[d]->SetColor(kAzure-2);
      //datasets[d]->SetColor(kMagenta);
      //datasets[d]->SetColor(kAzure+6);
    }
    if ( dataSetName.find("ST") != std::string::npos || dataSetName.find("SingleTop") != std::string::npos )
    {
      datasets[d]->SetTitle("ST");
      //datasets[d]->SetColor(kBlue-2);
      //datasets[d]->SetColor(kOrange-4);  // 595, 615, 800
      datasets[d]->SetColor(kGreen-10);  //kGreen-7,9
    }
  }
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
    
  if (makePlots)
  {
    InitMSPlots();
    InitHisto1D();
    InitHisto2D();
  }
  
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName;
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
  
  
  
  ////////////////////////////////////
  /// Loop over datasets
  ////////////////////////////////////
  
  //for (int d = 0; d < 1; d++)
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    clock_t startDataSet = clock();
    
    ClearMetaData();
    
    dataSetName = datasets[d]->Name();
    isData = false; isTTbar = false; isST = false; isSTtW = false; isSTt = false; isOther = false;
    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA") != std::string::npos )
    {
      isData = true;
      pathNtuples = pathNtuplesData;
    }
    else if ( dataSetName.find("TT") != std::string::npos )
    {
      isTTbar = true;
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
    
    if (isData && skipData)
    {
      cout << "Skipping data..." << endl;
      continue;
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
    
    
    
    string ntupleFileName = pathNtuples+"Ntuples_"+dataSetName+".root";
    
    /// Change name of ttbar dataset to TT (whether it is nominal or widthxX)
    if (isTTbar) dataSetName = "TT";
    
    tFileMap[dataSetName.c_str()] = new TFile(ntupleFileName.c_str(),"READ"); //create TFile for each dataset
    
    string tTreeName = "tree";
    string tStatsTreeName = "stats";
    
    /// Get meta data
    tStatsTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tStatsTreeName.c_str());
    GetMetaData(tStatsTree[dataSetName.c_str()], isData);
    
    tStatsTree[(dataSetName).c_str()]->GetEntry(0);
    
    /// eqLumi calculation
    if (isData)
    {
      lumiWeight = 1.;
    }
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
    
    if (test && isData)
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
    if (test || testHistos) endEvent = 1000;
    for (int ievt = 0; ievt < endEvent; ievt++)
    {
      ClearObjects();
      
      if (ievt%10000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
      
      
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
        if (makePlots)
        {
          MSPlot["nPVs_beforePU_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF);
          MSPlot["nPVs_afterPU_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF*puSF);
        }
      }
      else if (makePlots)
      {
        MSPlot["nPVs_beforePU_"]->Fill(nvtx, datasets[d], false, lumiWeight);
        MSPlot["nPVs_afterPU_"]->Fill(nvtx, datasets[d], false, lumiWeight);
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
        if ( jet_pt[iJet] > 30. && jet_pt[iJet] < 250. ) selectedJets.push_back(jet);
      }
      
      for (int iJet = 0; iJet < selectedJets.size(); iJet++)
      {
        if ( jet_bdiscr[iJet] > CSVv2Medium )
        {
          selectedBJets.push_back(selectedJets[iJet]);
        }
      }
      
      FillPlots(d, "_mb");
      
      if ( selectedJets.size() < 4 ) continue;
      
      FillPlots(d, "");
      
      selectedBJets.clear();
      for (int iJet = 0; iJet < selectedJets.size(); iJet++)
      {
        if ( jet_bdiscr[iJet] > CSVv2Medium )
        {
          selectedBJets.push_back(selectedJets[iJet]);
          bJetId.push_back(iJet);  /// selectedBJets[j] = selectedJets[bJetId[j]]
        }
      }
      
      
      if ( selectedBJets.size() < 2 ) continue;
      
      FillPlots(d, "_2b");
      
      
      if ( selectedJets.size() > 4 ) continue;
      
      FillPlots(d, "_2b_4j");
      
      
      
      
      
      
    }  // end loop events
    
    
    
    tFileMap[dataSetName.c_str()]->Close();
    
    timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
    
  }  // end loop datasets
  
  
  cout << endl << "Processing time per dataset: " << endl;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
  }
  
  
  
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
  tree->SetBranchAddress("cutFlowWeighted", cutFlowWeighted, &b_cutFlowWeighted);
  tree->SetBranchAddress("cutFlow2Weighted", cutFlow2Weighted, &b_cutFlow2Weighted);
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
    if (newTrees)
    {
      tree->SetBranchAddress("nofEventsWithGenTop", &nofEventsWithGenTop, &b_nofEventsWithGenTop);
      tree->SetBranchAddress("nofEventsWithGenAntiTop", &nofEventsWithGenAntiTop, &b_nofEventsWithGenAntiTop);
    }
    tree->SetBranchAddress("nofTTEventsWithoutAGenTop", &nofTTEventsWithoutAGenTop, &b_nofTTEventsWithoutAGenTop);
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
  tree->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
  tree->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
  tree->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
  tree->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
  tree->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
  tree->SetBranchAddress("muon_E", muon_E, &b_muon_E);
  tree->SetBranchAddress("muon_M", muon_M, &b_muon_M);
  tree->SetBranchAddress("muon_d0", muon_d0, &b_muon_d0);
  if (newTrees) tree->SetBranchAddress("muon_dz", muon_dz, &b_muon_dz);
  tree->SetBranchAddress("muon_chargedHadronIso", muon_chargedHadronIso, &b_muon_chargedHadronIso);
  tree->SetBranchAddress("muon_neutralHadronIso", muon_neutralHadronIso, &b_muon_neutralHadronIso);
  tree->SetBranchAddress("muon_photonIso", muon_photonIso, &b_muon_photonIso);
  tree->SetBranchAddress("muon_puChargedHadronIso", muon_puChargedHadronIso, &b_muon_puChargedHadronIso);
  tree->SetBranchAddress("muon_relIso", muon_relIso, &b_muon_relIso);
  tree->SetBranchAddress("muon_pfIso", muon_pfIso, &b_muon_pfIso);
  if (newTrees)
  {
    tree->SetBranchAddress("muon_chi2", muon_chi2, &b_muon_chi2);
    tree->SetBranchAddress("muon_nofValidHits", muon_nofValidHits, &b_muon_nofValidHits);
    tree->SetBranchAddress("muon_nofValidMuHits", muon_nofValidMuHits, &b_muon_nofValidMuHits);
    tree->SetBranchAddress("muon_nofValidPixelHits", muon_nofValidPixelHits, &b_muon_nofValidPixelHits);
    tree->SetBranchAddress("muon_nofMatchedStations", muon_nofMatchedStations, &b_muon_nofMatchedStations);
    tree->SetBranchAddress("muon_nofTrackerLayersWithMeasurement", muon_nofTrackerLayersWithMeasurement, &b_muon_nofTrackerLayersWithMeasurement);
  }
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
  if (newTrees && ! isData) tree->SetBranchAddress("jet_hadronFlavour", jet_hadronFlavour, &b_jet_hadronFlavour);   // TEMPORARILY !
  if (! isData)
  {
    if (newTrees)
    {
      tree->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets);
      tree->SetBranchAddress("genJet_pdgId", genJet_pdgId, &b_genJet_pdgId);
      tree->SetBranchAddress("genJet_charge", genJet_charge, &b_genJet_charge);
      tree->SetBranchAddress("genJet_pt", genJet_pt, &b_genJet_pt);
      tree->SetBranchAddress("genJet_phi", genJet_phi, &b_genJet_phi);
      tree->SetBranchAddress("genJet_eta", genJet_eta, &b_genJet_eta);
      tree->SetBranchAddress("genJet_E", genJet_E, &b_genJet_E);
      tree->SetBranchAddress("genJet_M", genJet_M, &b_genJet_M);
    }
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
    tree->SetBranchAddress("hasGenAntiTop", &hasGenAntiTop, &b_hasGenAntiTop);
    tree->SetBranchAddress("btagSF", &btagSF, &b_btagSF);
    tree->SetBranchAddress("puSF", &puSF, &b_puSF);
    tree->SetBranchAddress("muonIdSF_BCDEF", muonIdSF_BCDEF, &b_muonIdSF_BCDEF);
    tree->SetBranchAddress("muonIdSF_GH", muonIdSF_GH, &b_muonIdSF_GH);
    tree->SetBranchAddress("muonIsoSF_BCDEF", muonIsoSF_BCDEF, &b_muonIsoSF_BCDEF);
    tree->SetBranchAddress("muonIsoSF_GH", muonIsoSF_GH, &b_muonIsoSF_GH);
    tree->SetBranchAddress("muonTrigSF_BCDEF", muonTrigSF_BCDEF, &b_muonTrigSF_BCDEF);
    tree->SetBranchAddress("muonTrigSF_GH", muonTrigSF_GH, &b_muonTrigSF_GH);
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
  /// Do not split ttbar dataset
  MSPlot["nPVs_beforePU_"] = new MultiSamplePlot(datasets, "nPVs_beforePU_", 46, -0.5, 45.5, "# PVs");
  MSPlot["nPVs_afterPU_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_", 46, -0.5, 45.5, "# PVs");
  
  
  for (int i = 0; i < nSuffixes; i++)
  {
    MSPlot["rho_bPU"+suffixes[i]] = new MultiSamplePlot(datasets, ("rho_bPU"+suffixes[i]).c_str(), 41, -0.5, 40.5, "#rho");
    MSPlot["rho"+suffixes[i]] = new MultiSamplePlot(datasets, ("rho"+suffixes[i]).c_str(), 41, -0.5, 40.5, "#rho");
    MSPlot["nJets"+suffixes[i]] = new MultiSamplePlot(datasets, ("nJets"+suffixes[i]).c_str(), 13, -0.5, 12.5, "# jets");
    MSPlot["nBJets"+suffixes[i]] = new MultiSamplePlot(datasets, ("nBJets"+suffixes[i]).c_str(), 9, -0.5, 8.5, "# b jets");
    MSPlot["allDR"+suffixes[i]] = new MultiSamplePlot(datasets, ("allDR"+suffixes[i]).c_str(), 35, 0., 3.5, "#Delta R (j,j)");
    MSPlot["minDR"+suffixes[i]] = new MultiSamplePlot(datasets, ("minDR"+suffixes[i]).c_str(), 35, 0., 3.5, "min #Delta R (j,j)");
    MSPlot["maxDR"+suffixes[i]] = new MultiSamplePlot(datasets, ("maxDR"+suffixes[i]).c_str(), 45, 1., 5.5, "max #Delta R (j,j)");
    
    MSPlot["muon_relIso"+suffixes[i]] = new MultiSamplePlot(datasets, ("muon_relIso"+suffixes[i]).c_str(), 16, 0., 0.16, "relIso");
    MSPlot["muon_d0"+suffixes[i]] = new MultiSamplePlot(datasets, ("muon_d0"+suffixes[i]).c_str(), 60, 0., 0.003, "d_{0}");
    MSPlot["muon_dz"+suffixes[i]] = new MultiSamplePlot(datasets, ("muon_dz"+suffixes[i]).c_str(), 50, 0., 0.001, "d_{z}");
    MSPlot["muon_chi2"+suffixes[i]] = new MultiSamplePlot(datasets, ("muon_chi2"+suffixes[i]).c_str(), 40, 0., 10, "#chi^{2}");
    MSPlot["muon_nTrackLayers"+suffixes[i]] = new MultiSamplePlot(datasets, ("muon_nTrackLayers"+suffixes[i]).c_str(), 21, 4.5, 25.5, "# meas. in tracker layers");
    MSPlot["muon_nValidHits"+suffixes[i]] = new MultiSamplePlot(datasets, ("muon_nValidHits"+suffixes[i]).c_str(), 41, -0.5, 40.5, "# valid hits");
    MSPlot["muon_nValidPixelHits"+suffixes[i]] = new MultiSamplePlot(datasets, ("muon_nValidPixelHits"+suffixes[i]).c_str(), 11, -0.5, 10.5, "# valid pixel hits");
    MSPlot["muon_nValidMuHits"+suffixes[i]] = new MultiSamplePlot(datasets, ("muon_nValidMuHits"+suffixes[i]).c_str(), 56, -0.5, 55.5, "# valid muon hits");
    MSPlot["muon_nMatchedStations"+suffixes[i]] = new MultiSamplePlot(datasets, ("muon_nMatchedStations"+suffixes[i]).c_str(), 9, -0.5, 8.5, "# matched stations");
  }
}

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  for (int i = 0; i < nSuffixes; i++)
  {
    histo1D["allDR"+suffixes[i]] = new TH1F(("allDR"+suffixes[i]).c_str(),"Delta R between any two jets; #Delta R(j,j)", 35, 0, 3.5);
    histo1D["minDR"+suffixes[i]] = new TH1F(("minDR"+suffixes[i]).c_str(),"Minimal delta R between two jets; #Delta R(j,j)", 35, 0, 3.5);
    histo1D["maxDR"+suffixes[i]] = new TH1F(("maxDR"+suffixes[i]).c_str(),"Maximal delta R between two jets; #Delta R(j,j)", 35, 0, 3.5);
    
    histo1D["CSVv2Discr_b"+suffixes[i]] = new TH1F(("CSVv2Discr_b"+suffixes[i]).c_str(),"CSVv2 discriminator value for b partons; CSVv2 discriminant value", 48, 0., 1.2);
    histo1D["CSVv2Discr_c"+suffixes[i]] = new TH1F(("CSVv2Discr_c"+suffixes[i]).c_str(),"CSVv2 discriminator value for c partons; CSVv2 discriminant value", 48, 0., 1.2);
    histo1D["CSVv2Discr_udsg"+suffixes[i]] = new TH1F(("CSVv2Discr_udsg"+suffixes[i]).c_str(),"CSVv2 discriminator value for udsg partons; CSVv2 discriminant value", 48, 0., 1.2);
  }
  
  
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
  
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
    nofTTEventsWithoutAGenTop = 0;
  }
  
  strSyst = "";
  nEventsDataSet = 0;
  xSection = 1.;
  eqLumi = 1.;
  lumiWeight = 1.;
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
  muon_charge[0] = 0.;
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
    jet_charge[i] = 0.;
    jet_nConstituents[i] = -1;
    jet_nChConstituents[i] = -1;
    jet_pt[i] = 0.;
    jet_phi[i] = 0.;
    jet_eta[i] = 0.;
    jet_E[i] = 0.;
    jet_M[i] = 0.;
    jet_bdiscr[i] = -1.;
    jet_hadronFlavour[i] = -1.;
  }
  
  if (! isData)
  {
    nGenJets = -1;
    for (Int_t i = 0; i < 40; i++)
    {
      genJet_pdgId[i] = 0;
      genJet_charge[i] = 0.;
      genJet_pt[i] = 0.;
      genJet_phi[i] = 0.;
      genJet_eta[i] = 0.;
      genJet_E[i] = 0.;
      genJet_M[i] = 0.;
    }
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
    hasGenAntiTop = false;
    
    btagSF = 1.;
    puSF = 1.;
    muonIdSF_BCDEF[0] = 1.;
    muonIdSF_GH[0] = 1.;
    muonIsoSF_BCDEF[0] = 1.;
    muonIsoSF_GH[0] = 1.;
    muonTrigSF_BCDEF[0] = 1.;
    muonTrigSF_GH[0] = 1.;
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
  selectedLepton.clear();
  selectedJets.clear();
  selectedBJets.clear();
}

void ClearVars()
{
  scaleFactor = 1.;
  thisLeptonSF = 1.;
  thisLeptonIdSF = 1.;
  thisLeptonIsoSF = 1.;
  thisLeptonTrigSF = 1.;
  
  bJetId.clear();
}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearVars();
}

void FillPlots(int d, string suffix)
{
  if (isData) MSPlot["rho_bPU"+suffix]->Fill(rho, datasets[d], false, lumiWeight*thisLeptonSF*btagSF);
  else MSPlot["rho_bPU"+suffix]->Fill(rho, datasets[d], true, lumiWeight*thisLeptonSF*btagSF);
  MSPlot["rho"+suffix]->Fill(rho, datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["nJets"+suffix]->Fill(selectedJets.size(), datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["nBJets"+suffix]->Fill(selectedBJets.size(), datasets[d], true, lumiWeight*scaleFactor);
  
  double minDR = 999., maxDR = -999.;
  for (int iJet = 0; iJet < selectedJets.size(); iJet++)
  {
    if (isTTbar)
    {
      if ( jet_hadronFlavour[iJet] == 5 ) histo1D["CSVv2Discr_b"+suffix]->Fill(jet_bdiscr[iJet]);
      else if ( jet_hadronFlavour[iJet] == 4 ) histo1D["CSVv2Discr_c"+suffix]->Fill(jet_bdiscr[iJet]);
      else histo1D["CSVv2Discr_udsg"+suffix]->Fill(jet_bdiscr[iJet]);
    }
    
    for (int jJet = iJet+1; jJet < selectedJets.size(); jJet++)
    {
      tempDR = ROOT::Math::VectorUtil::DeltaR(selectedJets[iJet], selectedJets[jJet]);
      MSPlot["allDR"+suffix]->Fill(tempDR, datasets[d], true, lumiWeight*scaleFactor);
      if (isTTbar) histo1D["allDR"+suffix]->Fill(tempDR, lumiWeight*scaleFactor);
      if ( tempDR < minDR ) minDR = tempDR;
      if ( tempDR > maxDR ) maxDR = tempDR;
    }
  }
  
  MSPlot["minDR"+suffix]->Fill(minDR, datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["maxDR"+suffix]->Fill(maxDR, datasets[d], true, lumiWeight*scaleFactor);
  if (isTTbar)
  {
    histo1D["minDR"+suffix]->Fill(minDR, lumiWeight*scaleFactor);
    histo1D["maxDR"+suffix]->Fill(maxDR, lumiWeight*scaleFactor);
  }
  
  MSPlot["muon_relIso"+suffix]->Fill(muon_relIso[0], datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["muon_d0"+suffix]->Fill(muon_d0[0], datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["muon_dz"+suffix]->Fill(muon_dz[0], datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["muon_chi2"+suffix]->Fill(muon_chi2[0], datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["muon_nTrackLayers"+suffix]->Fill(muon_nofTrackerLayersWithMeasurement[0], datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["muon_nValidHits"+suffix]->Fill(muon_nofValidHits[0], datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["muon_nValidPixelHits"+suffix]->Fill(muon_nofValidPixelHits[0], datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["muon_nValidMuHits"+suffix]->Fill(muon_nofValidMuHits[0], datasets[d], true, lumiWeight*scaleFactor);
  MSPlot["muon_nMatchedStations"+suffix]->Fill(muon_nofMatchedStations[0], datasets[d], true, lumiWeight*scaleFactor);
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

long GetNEvents(TTree* fChain, string var, unsigned int index, bool isData)
{
  GetMetaData(fChain, isData);
  long varNew = 0;
  for (unsigned int iEntry = 0; iEntry < fChain->GetEntries(); iEntry++)
  {
    fChain->GetEntry(iEntry);
    varNew += (fChain->FindLeaf(var.c_str()))->GetValueLong64(index);
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
