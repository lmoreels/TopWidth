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
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
//#include "../macros/Style.C"

// user defined
#include "Tools/interface/EventReweighting.h"


using namespace std;
using namespace TopTree;


bool test = false;
bool testHistos = false;
bool testTTbarOnly = false;
bool doGenOnly = false;
bool makePlots = true;
bool makeReweightedPlots = true; 
bool calculateAverageMass = false;

bool doMETCleaning = true;
bool applyLeptonSF = true;
bool applyPU = true;
bool applyBTagSF = true;

bool doReweighting = false;
bool applyMassSF = false;
double scaleMass = 172.5.;

string systStr = "nominal";
string ntupleSystDate = "170803";
pair<string,string> whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos )
  {
    return pair<string,string>("170712","170920");
  }
  else if ( syst.find("JESup") != std::string::npos ) return pair<string,string>("170904","170904");
  else if ( syst.find("JESdown") != std::string::npos ) return pair<string,string>("170905","170905");
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    return pair<string,string>("170712","170920");
  }
}
pair<string,string> ntupleDate = whichDate(systStr);
int verbose = 2;

string pathNtuples = "";
string pathNtuplesMC = "";
string pathNtuplesData = "";
string pathNtuplesSyst = "";
string dateString = "";
string outputDirLL = "LikelihoodTemplates/";
string inputDirLL = "";
string inputDateLL = "170918_1300/";  // TT nominal
bool isData = false;
bool isTTbar = false;

int nofHardSelected = 0;
int nofMETCleaned = 0;
int nofMatchedEvents = 0;
int nofHadrMatchedEvents = 0;
int nofCorrectlyMatched = 0;
int nofNotCorrectlyMatched = 0;

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
std::array<double, 5> aveTopMass = {171.791, 169.728, 186.218, 182.657, 167.919};  // Res 170915
//std::array<double, 5> aveTopMass = {171.791, 169.728, 185.155, 181.914, 167.723};  // Res 170608 Single Gaus

/// # events
vector<long> nEventsTot, nEventsSelc, nEventsHardSel, nEventsAKF;


// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TGraph*> graph;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;

vector < Dataset* > datasets;

ofstream txtMassGenPMatched, txtMassGenJMatched, txtMassRecoBKF, txtMassRecoAKF, txtMassRecoAKFMatched;

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
void InitHisto2D();
void InitHisto2DGen();
void InitHisto2DReweighted();
void TruthMatching(vector<TLorentzVector> partons, vector<TLorentzVector> selectedJets, pair<unsigned int, unsigned int> *MCPermutation);
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearVars();
void ClearObjects();
long GetNEvents(TTree* fChain, string var, bool isData);
void GetEraFraction(double* fractions);
void CalculateAverageTopMass();



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
TBranch        *b_nofTTEventsWithoutAGenTop;   //!
TBranch        *b_nofTTEventsWithoutGenTop;   //!
TBranch        *b_nofTTEventsWithoutGenAntiTop;   //!
TBranch        *b_nofTTEventsWithoutBothGenTopsWithStatus22;   //!
TBranch        *b_nofTTEventsWithoutGenTopWithStatus22;   //!
TBranch        *b_nofTTEventsWithoutGenAntiTopWithStatus22;   //!
TBranch        *b_nofTTEventsWithoutBothGenTopsWithStatus62;   //!
TBranch        *b_nofTTEventsWithoutGenTopWithStatus62;   //!
TBranch        *b_nofTTEventsWithoutGenAntiTopWithStatus62;   //!


double lumiWeight, numWeight, scaleFactor, massSF;
bool foundTop62, foundAntiTop62;
vector<unsigned int> bJetId;
double bdiscrTop, bdiscrTop2, tempbdiscr;
int labelB1, labelB2;
double massHadTopQ, massLepTopQ;


/// Define TLVs
TLorentzVector muon, jet, mcpart;
TLorentzVector WCandidate;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
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
int bqqDecay[] = {-9999, -9999, -9999}, blvDecay[] = {-9999, -9999, -9999};
vector<unsigned int> partonId;


ofstream txtDebugTopMass;


/// Variables
double matched_W_mass_q, matched_top_mass_q;
double matched_W_mass_j, matched_top_mass_j, matched_top_mass_j_akF;
double matched_mlb_corr, matched_ttbarMass_corr, matched_dR_lep_b_corr;
double matched_mlb_wrong, matched_ttbarMass_wrong, matched_dR_lep_b_wrong;


/// Input samples & reweighting
int thisGenMassId;
double genMassArray[] = {169.5, 172.5, 175.5};
string genMassString[] = {"g169p5", "g172p5", "g175p5"};
const int nGenMasses = sizeof(genMassArray)/sizeof(genMassArray[0]);

int nEventsTT[] = {56273062, 153843293, 39649356};

double reweightArray[] = {169.5, 171.5, 172.5, 173.5, 175.5};
string reweightString[] = {"s169p5", "s171p5", "s172p5", "s173p5", "s175p5"};
const int nReweightings = sizeof(reweightArray)/sizeof(reweightArray[0]);

double m_top = -1.;
double m_antitop = -1.;
double m_bqq = -1.;
double m_blv = -1.;
double m_hadr = -1;
double m_lept = -1;
double evWeight_hadr = 1.;
double evWeight_lept = 1.;
double evWeight_prod = 1.;



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
  else                                                 cout << "nominal        *" << endl;
  cout << "*   - Jet/lepton Cleaning                   *" << endl;
  if (doMETCleaning) cout << "*   - MET Cleaning                          *" << endl;
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
  if (test) makePlots = false;
  if (testHistos)
  {
    makePlots = true;
    doGenOnly = false;
  }
  if (! makePlots)
  {
    makeReweightedPlots = false; 
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
  
  pathNtuplesMC = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.first+"/";
  pathNtuplesSyst = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleSystDate+"/";
  cout << "Using Ntuples from " << ntupleDate.first << ". This corresponds to systematics: " << systStr << endl;
  if (calculateAverageMass) cout << "Calculating average mass values..." << endl;
  if (testHistos) cout << "Testing histogram consistency..." << endl;
  if (doGenOnly) cout << "Running only matching..." << endl;
  if (applyMassSF) cout << "TTbar sample mass will be scaled by a factor " << scaleMass << endl;
  else scaleMass = 172.5;
  
  /// xml file
  string xmlFileName ="config/topWidth_ttbarOnlyMass.xml";
  
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
  
  if (makePlots)
  {
    InitHisto1D();
    InitHisto2D();
  }  
  
  if (calculateAverageMass)
  {
    mkdir("averageMass/",0777);
    txtMassGenPMatched.open(("averageMass/mass_genp_matched_"+dateString+".txt").c_str());
    txtMassGenJMatched.open(("averageMass/mass_genj_matched_"+dateString+".txt").c_str());
    txtMassRecoBKF.open(("averageMass/mass_recoBKF_"+dateString+".txt").c_str());
    txtMassRecoAKF.open(("averageMass/mass_recoAKF_"+dateString+".txt").c_str());
    txtMassRecoAKFMatched.open(("averageMass/mass_recoAKF_matched_"+dateString+".txt").c_str());
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
      thisGenMassId = 1;
      pathNtuples = pathNtuplesMC;
      if ( dataSetName.find("mass") != std::string::npos || dataSetName.find("Mass") != std::string::npos )
      {
        pathNtuples = pathNtuplesSyst;
        if ( dataSetName.find("169") != std::string::npos )
        {
          thisGenMassId = 0;
          scaleMass = 169.5;
        }
        else if ( dataSetName.find("175") != std::string::npos )
        {
          thisGenMassId = 2;
          scaleMass = 175.5;
        }
        applyMassSF = false;
        doReweighting = false;
      }
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
    numWeight = (double)nEventsTT[1]/(double)nEventsTT[thisGenMassId];  // Don't put extra SFs on reweighted samples
    
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
      
      
      if (! applyMassSF ) massSF = 1.;
      else if ( applyMassSF && ! isTTbar ) massSF = 1.;  // also for data
      
      
      
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
      
      
      
      /////////////////////////////////////////
      ///  Scale factor ttbar sample mass  ///
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
      
      if ( applyMassSF && isTTbar )
      {
        massSF = rew->MassEventWeightCalculatorNonRel(massHadTopQ, scaleMass);
        
        if ( massSF != massSF )  // massSF = NaN
        {
          cout << "Event " << ievt << ";  Event nb. " << evt_num << "; Run nb. " << run_num << endl;
          cout << "Top mass: " << (mcParticles[topQuark]).M() << "; antiTop mass: " << (mcParticles[antiTopQuark]).M() << "; Lepton charge: " << muon_charge[0] << "; mass SF: " << massSF << endl;
          
          continue;
        }
        
        if (makePlots) histo1D["mass_SF"]->Fill(massSF);
      }  // end applyMassSF
      
      
      
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
            if ( mc_pdgId[partonId[i]] == 5 && mc_mother[partonId[i]] == pdgID_top ) bqqDecay[0] = partonId[i];
            else if ( abs(mc_pdgId[partonId[i]]) < 5 && mc_mother[partonId[i]] == 24 && mc_granny[partonId[i]] == pdgID_top )
            {
              //cout << "In hadronic loop: " << setw(3) << right << partonId[i] << "  pdgId: " << setw(3) << mc_pdgId[partonId[i]] << "  Mother: " << setw(4) << mc_mother[partonId[i]] << "  Granny: " << setw(4) << mc_granny[partonId[i]] << endl;

              if ( bqqDecay[1] == -9999 ) bqqDecay[1] = partonId[i];
              else if ( bqqDecay[2] == -9999 ) bqqDecay[2] = partonId[i];
              else { cerr << "ERROR: Too many partons for hadronic top decay!  Event: " << ievt << endl; continue; }
            }
            
            if ( mc_pdgId[partonId[i]] == -5 && mc_mother[partonId[i]] == -pdgID_top ) blvDecay[0] = partonId[i];  // leptonic b
          }
          else if ( muPlusFromTop )
          {
            // hadronic decay  --> antitop
            if ( mc_pdgId[partonId[i]] == -5 && mc_mother[partonId[i]] == -pdgID_top ) bqqDecay[0] = partonId[i];
            else if ( abs(mc_pdgId[partonId[i]]) < 5 && mc_mother[partonId[i]] == -24 && mc_granny[partonId[i]] == -pdgID_top )
            {
              //cout << "In hadronic loop: " << setw(3) << right << partonId[i] << "  pdgId: " << setw(3) << mc_pdgId[partonId[i]] << "  Mother: " << setw(4) << mc_mother[partonId[i]] << "  Granny: " << setw(4) << mc_granny[partonId[i]] << endl;
              
              if ( bqqDecay[1] == -9999 ) bqqDecay[1] = partonId[i];
              else if ( bqqDecay[2] == -9999 ) bqqDecay[2] = partonId[i];
              else { cerr << "ERROR: Too many partons for hadronic antitop decay!  Event: " << ievt << endl; continue; }
            }
            
            if ( mc_pdgId[partonId[i]] == 5 && mc_mother[partonId[i]] == pdgID_top ) blvDecay[0] = partonId[i];  // leptonic b
          }
        }  // end loop partons
        blvDecay[1] = genmuon;
        blvDecay[2] = genneutrino;
        
        if ( bqqDecay[0] != -9999 && bqqDecay[1] != -9999 && bqqDecay[2] != -9999 && blvDecay[0] != -9999 && blvDecay[1] != -9999 && blvDecay[2] != -9999)
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
          m_bqq = (mcParticles[bqqDecay[0]] + mcParticles[bqqDecay[1]] + mcParticles[bqqDecay[2]]).M();
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
          
          if (doReweighting && makeReweightedPlots)
          {
            for (int s = 0; s < nReweightings; s++)
            {
              evWeight_hadr = rew->MassEventWeightCalculatorNonRel(m_hadr, reweightArray[s]);
              evWeight_lept = rew->MassEventWeightCalculatorNonRel(m_lept, reweightArray[s]);
              evWeight_prod = evWeight_hadr*evWeight_lept;
              
              histo1D[("top_mass_hadr_gen_"+reweightString[s]).c_str()]->Fill(m_hadr, evWeight_hadr);
              histo1D[("top_mass_hadr_gen_"+reweightString[s]+"_fewerBins").c_str()]->Fill(m_hadr, evWeight_hadr);
              histo1D[("top_mass_lept_gen_"+reweightString[s]).c_str()]->Fill(m_lept, evWeight_hadr);
              histo1D[("bqq_mass_gen_"+reweightString[s]).c_str()]->Fill(m_bqq, evWeight_hadr);
              histo1D[("blv_mass_gen_"+reweightString[s]).c_str()]->Fill(m_blv, evWeight_hadr);
              
              histo1D[("top_mass_hadr_gen_prodSF_"+reweightString[s]).c_str()]->Fill(m_hadr, evWeight_prod);
              histo1D[("top_mass_hadr_gen_prodSF_"+reweightString[s]+"_fewerBins").c_str()]->Fill(m_hadr, evWeight_prod);
//               for (int i = 0; i < nGenMasses; i++)
//               {
//                 if ( reweightArray[s] == genMassArray[i] )
//                   histo1D[("top_mass_hadr_gen_prodSF_"+reweightString[s]+"_th").c_str()]->Fill(m_hadr, evWeight_prod);
//               }
              histo1D[("top_mass_lept_gen_prodSF_"+reweightString[s]).c_str()]->Fill(m_lept, evWeight_prod);
              histo1D[("bqq_mass_gen_prodSF_"+reweightString[s]).c_str()]->Fill(m_bqq, evWeight_prod);
              histo1D[("blv_mass_gen_prodSF_"+reweightString[s]).c_str()]->Fill(m_blv, evWeight_prod);
              
              histo1D[("Mass_SF_hadr_"+reweightString[s]).c_str()]->Fill(evWeight_hadr);
              histo1D[("Mass_SF_lept_"+reweightString[s]).c_str()]->Fill(evWeight_lept);
              histo1D[("Mass_SF_prod_"+reweightString[s]).c_str()]->Fill(evWeight_prod);
              
              histo2D[("Mass_SF_hadr_lept_"+reweightString[s]).c_str()]->Fill(evWeight_hadr, evWeight_lept);
              histo2D[("Mass_SF_hadr_prod_"+reweightString[s]).c_str()]->Fill(evWeight_hadr, evWeight_prod);
              histo2D[("top_mass_hadr_gen_vs_weight_"+reweightString[s]).c_str()]->Fill(m_hadr, evWeight_hadr);
              //histo2D[("top_mass_lept_gen_vs_weight_"+reweightString[s]).c_str()]->Fill(m_lept, evWeight_lept);
              histo2D["top_mass_hadr_lept_gen_"+reweightString[s]]->Fill(m_hadr, m_lept, evWeight_hadr);
              histo2D["top_mass_hadr_lept_gen_prodSF_"+reweightString[s]]->Fill(m_hadr, m_lept, evWeight_prod);
            }
          }
          if (makePlots)
          {
            histo1D[("top_mass_hadr_gen_"+genMassString[thisGenMassId]).c_str()]->Fill(m_hadr, numWeight);
            histo1D[("top_mass_hadr_gen_"+genMassString[thisGenMassId]+"_fewerBins").c_str()]->Fill(m_hadr, numWeight);
//            histo1D[("top_mass_hadr_gen_"+genMassString[thisGenMassId]+"_th").c_str()]->Fill(m_hadr, numWeight);
            histo1D[("top_mass_lept_gen_"+genMassString[thisGenMassId]).c_str()]->Fill(m_lept, numWeight);
            histo1D[("bqq_mass_gen_"+genMassString[thisGenMassId]).c_str()]->Fill(m_bqq, numWeight);
            histo1D[("blv_mass_gen_"+genMassString[thisGenMassId]).c_str()]->Fill(m_blv, numWeight);
            
            histo2D["top_mass_hadr_lept_gen_"+genMassString[thisGenMassId]]->Fill(m_hadr, m_lept, numWeight);
            histo2D["top_mass_hadr_bqq_gen_"+genMassString[thisGenMassId]]->Fill(m_hadr, m_bqq, numWeight);
            histo2D["top_mass_lept_blv_gen_"+genMassString[thisGenMassId]]->Fill(m_lept, m_blv, numWeight);
            histo2D["top_mass_bqq_blv_gen_"+genMassString[thisGenMassId]]->Fill(m_bqq, m_blv, numWeight);
          }
        }
      }
      
      
    }  // end loop events
    
    
    
    cout << endl;  /// Stronger selection in this analyser compared to Ntuples ==> endEvent --> nofHardSelected
    cout << "Number of events with exactly 4 jets with pT > 30 GeV: " << nofHardSelected << " (" << 100*((float)nofHardSelected/(float)endEvent) << "%)" << endl;
    cout << "Number of events with clean MET: " << nofMETCleaned << " (" << 100*((float)nofMETCleaned/(float)nofHardSelected) << "%)" << endl;
    
    //if ( isTTbar || dataSetName.find("ST") != std::string::npos )
    if (! isData && nofHadrMatchedEvents > 0 )
    {
      cout << "Number of matched events: " << setw(8) << right << nofMatchedEvents << endl;
      cout << "Number of events with hadronic top matched (before KF): " << setw(8) << right << nofHadrMatchedEvents << " (" << 100*((float)nofHadrMatchedEvents/(float)nofMETCleaned) << "%)" << endl;
    }  // end ! isData
    
    
    
    tFileMap[dataSetName.c_str()]->Close();
    
    timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
    
  }  // end loop datasets
  
  
  if (calculateAverageMass)
  {
    txtMassGenPMatched.close();
    txtMassGenJMatched.close();
    txtMassRecoBKF.close();
    txtMassRecoAKF.close();
    txtMassRecoAKFMatched.close();
  }
  
  if (applyMassSF) txtDebugTopMass.close();
  
  cout << endl;
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
    CalculateAverageTopMass();
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
  if (makeReweightedPlots)
  {
    InitHisto1DReweighted();
  }
  
  
  /// SFs
  histo1D["mass_SF"] = new TH1F("mass_SF", "Scale factor to change the ttbar distribution mass; mass SF", 500, 0, 5);
  
}

void InitHisto1DGen()
{
  for (int g = 0; g < nGenMasses; g++)
  {
    histo1D[("top_mass_hadr_gen_"+genMassString[g]).c_str()] = new TH1F(("top_mass_hadr_gen_"+genMassString[g]).c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 4000, 120, 220);
    histo1D[("top_mass_hadr_gen_"+genMassString[g]+"_fewerBins").c_str()] = new TH1F(("top_mass_hadr_gen_"+genMassString[g]+"_fewerBins").c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 800, 120, 220);
    histo1D[("top_mass_lept_gen_"+genMassString[g]).c_str()] = new TH1F(("top_mass_lept_gen_"+genMassString[g]).c_str(), "Mass of generated top quark with leptonic decay; M_{t_{lept}} [GeV]", 4000, 120, 220);
    histo1D[("bqq_mass_gen_"+genMassString[g]).c_str()] = new TH1F(("bqq_mass_gen_"+genMassString[g]).c_str(), "Mass of generated bqq quarks; M_{bqq} [GeV]", 2000, 50, 300);
    histo1D[("blv_mass_gen_"+genMassString[g]).c_str()] = new TH1F(("blv_mass_gen_"+genMassString[g]).c_str(), "Mass of generated b, lepton and neutrino; M_{blv} [GeV]", 2000, 50, 300);
    
  }
  
//   histo1D["top_mass_hadr_gen_g0p2_th"] = new TH1F("top_mass_hadr_gen_g0p2_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 4000, 120, 220);
//   histo1D["top_mass_hadr_gen_g0p5_th"] = new TH1F("top_mass_hadr_gen_g0p5_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 2000, 120, 220);
//   histo1D["top_mass_hadr_gen_g1_th"] = new TH1F("top_mass_hadr_gen_g1_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 1000, 120, 220);
//   histo1D["top_mass_hadr_gen_g4_th"] = new TH1F("top_mass_hadr_gen_g4_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 320, 120, 220);
//   histo1D["top_mass_hadr_gen_g8_th"] = new TH1F("top_mass_hadr_gen_g8_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 200, 120, 220);
}

void InitHisto1DReweighted()
{
  for (int s = 0; s < nReweightings; s++)
  {
    histo1D[("top_mass_hadr_gen_"+reweightString[s]).c_str()] = new TH1F(("top_mass_hadr_gen_"+reweightString[s]).c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 4000, 120, 220);
    histo1D[("top_mass_hadr_gen_"+reweightString[s]+"_fewerBins").c_str()] = new TH1F(("top_mass_hadr_gen_"+reweightString[s]+"_fewerBins").c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 800, 120, 220);
    histo1D[("top_mass_lept_gen_"+reweightString[s]).c_str()] = new TH1F(("top_mass_lept_gen_"+reweightString[s]).c_str(), "Mass of generated top quark with leptonic decay; M_{t_{lept}} [GeV]", 4000, 120, 220);
    histo1D[("bqq_mass_gen_"+reweightString[s]).c_str()] = new TH1F(("bqq_mass_gen_"+reweightString[s]).c_str(), "Mass of generated bqq quarks; M_{bqq} [GeV]", 2000, 50, 300);
    histo1D[("blv_mass_gen_"+reweightString[s]).c_str()] = new TH1F(("blv_mass_gen_"+reweightString[s]).c_str(), "Mass of generated b, lepton and neutrino; M_{blv} [GeV]", 2000, 50, 300);
    
    histo1D[("top_mass_hadr_gen_prodSF_"+reweightString[s]).c_str()] = new TH1F(("top_mass_hadr_gen_prodSF_"+reweightString[s]).c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 4000, 120, 220);
    histo1D[("top_mass_hadr_gen_prodSF_"+reweightString[s]+"_fewerBins").c_str()] = new TH1F(("top_mass_hadr_gen_prodSF_"+reweightString[s]+"_fewerBins").c_str(), "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 800, 120, 220);
    histo1D[("top_mass_lept_gen_prodSF_"+reweightString[s]).c_str()] = new TH1F(("top_mass_lept_gen_prodSF_"+reweightString[s]).c_str(), "Mass of generated top quark with leptonic decay; M_{t_{lept}} [GeV]", 4000, 120, 220);
    histo1D[("bqq_mass_gen_prodSF_"+reweightString[s]).c_str()] = new TH1F(("bqq_mass_gen_prodSF_"+reweightString[s]).c_str(), "Mass of generated bqq quarks; M_{bqq} [GeV]", 2000, 50, 300);
    histo1D[("blv_mass_gen_prodSF_"+reweightString[s]).c_str()] = new TH1F(("blv_mass_gen_prodSF_"+reweightString[s]).c_str(), "Mass of generated b, lepton and neutrino; M_{blv} [GeV]", 2000, 50, 300);
    
    /// SFs
    histo1D[("Mass_SF_hadr_"+reweightString[s]).c_str()] = new TH1F(("Mass_SF_hadr_"+reweightString[s]).c_str(), "Hadronic scale factor to change the ttbar distribution mass; mass SF", 5001, -0.0005, 5.0005);
    histo1D[("Mass_SF_lept_"+reweightString[s]).c_str()] = new TH1F(("Mass_SF_lept_"+reweightString[s]).c_str(), "Leptonic scale factor to change the ttbar distribution mass; mass SF", 5001, -0.0005, 5.0005);
    histo1D[("Mass_SF_prod_"+reweightString[s]).c_str()] = new TH1F(("Mass_SF_prod_"+reweightString[s]).c_str(), "Product of the hadronic and leptonic scale factor to change the ttbar distribution mass; mass SF", 5001, -0.0005, 5.0005);
  }
  
//   histo1D["top_mass_hadr_gen_prodSF_s0p2_th"] = new TH1F("top_mass_hadr_gen_prodSF_s0p2_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 4000, 120, 220);
//   histo1D["top_mass_hadr_gen_prodSF_s0p5_th"] = new TH1F("top_mass_hadr_gen_prodSF_s0p5_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 2000, 120, 220);
//   histo1D["top_mass_hadr_gen_prodSF_s1_th"] = new TH1F("top_mass_hadr_gen_prodSF_s1_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 1000, 120, 220);
//   histo1D["top_mass_hadr_gen_prodSF_s4_th"] = new TH1F("top_mass_hadr_gen_prodSF_s4_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 320, 120, 220);
//   histo1D["top_mass_hadr_gen_prodSF_s8_th"] = new TH1F("top_mass_hadr_gen_prodSF_s8_th", "Mass of generated top quark with hadronic decay; M_{t_{hadr}} [GeV]", 200, 120, 220);
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
  InitHisto2DGen();
  if (makeReweightedPlots) InitHisto2DReweighted();
}

void InitHisto2DGen()
{
  for (int g = 0; g < nGenMasses; g++)
  {
    histo2D["top_mass_hadr_lept_gen_"+genMassString[g]] = new TH2F(("top_mass_hadr_lept_gen_"+genMassString[g]).c_str(),("Mass of generated top quark with leptonic decay vs. hadronic decay ("+genMassString[g]+"); m_{t_{hadr}} [GeV]; m_{t_{lept}} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
    histo2D["top_mass_hadr_bqq_gen_"+genMassString[g]] = new TH2F(("top_mass_hadr_bqq_gen_"+genMassString[g]).c_str(),("Generated mass of bqq quarks vs. hadronically decaying top quark mass ("+genMassString[g]+"); m_{t_{hadr}} [GeV]; m_{bqq} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
    histo2D["top_mass_lept_blv_gen_"+genMassString[g]] = new TH2F(("top_mass_lept_blv_gen_"+genMassString[g]).c_str(),("Mass of generated b, lepton and neutrino vs. leptonically decaying top quark mass ("+genMassString[g]+"); m_{t_{lept}} [GeV]; m_{blv} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
    histo2D["top_mass_bqq_blv_gen_"+genMassString[g]] = new TH2F(("top_mass_bqq_blv_gen_"+genMassString[g]).c_str(),("Mass of generated b, lepton and neutrino vs. mass of bqq quarks ("+genMassString[g]+"); m_{bqq} [GeV]; m_{blv} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
  }
}

void InitHisto2DReweighted()
{
  for (int s = 0; s < nReweightings; s++)
  {
    histo2D[("Mass_SF_hadr_lept_"+reweightString[s]).c_str()] = new TH2F(("Mass_SF_hadr_lept_"+reweightString[s]).c_str(), "Leptonic vs. hadronic scale factor to change the ttbar distribution mass; SF_{h}; SF_{l}", 5001, -0.0005, 5.0005, 5001, -0.0005, 5.0005);
    histo2D[("Mass_SF_hadr_prod_"+reweightString[s]).c_str()] = new TH2F(("Mass_SF_hadr_prod_"+reweightString[s]).c_str(), "Product of hadronic and leptonic SF vs. hadronic scale factor to change the ttbar distribution mass; SF_{h}; SF_{h} #times SF_{l}", 5001, -0.0005, 5.0005, 5001, -0.0005, 5.0005);
    
    histo2D[("top_mass_hadr_gen_vs_weight_"+reweightString[s]).c_str()] = new TH2F(("top_mass_hadr_gen_vs_weight_"+reweightString[s]).c_str(), "Weights vs. mass of generated top quark (hadronic decay); m_{t_{hadr}} [GeV]; SF_{h}", 1000, 120, 220, 5001, -0.0005, 5.0005);
    //histo2D[("top_mass_lept_gen_vs_weight_"+reweightString[s]).c_str()] = new TH2F(("top_mass_lept_gen_vs_weight_"+reweightString[s]).c_str(), "Weights vs. mass of generated top quark (leptonic decay); m_{t_{lept}} [GeV]; SF_{l}", 1000, 120, 220, 5001, -0.0005, 5.0005);
    
    histo2D["top_mass_hadr_lept_gen_"+reweightString[s]] = new TH2F(("top_mass_hadr_lept_gen_"+reweightString[s]).c_str(),("Mass of generated top quark with leptonic decay vs. hadronic decay ("+reweightString[s]+"); m_{t_{hadr}} [GeV]; m_{t_{lept}} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
    histo2D["top_mass_hadr_lept_gen_prodSF_"+reweightString[s]] = new TH2F(("top_mass_hadr_lept_gen_prodSF_"+reweightString[s]).c_str(),("Mass of generated top quark with leptonic decay vs. hadronic decay ("+reweightString[s]+"); m_{t_{hadr}} [GeV]; m_{t_{lept}} [GeV]").c_str(), 1000, 50, 300, 1000, 50, 300);
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
  nofCorrectlyMatched = 0;
  nofNotCorrectlyMatched = 0;
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
    bqqDecay[i] = -9999;
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
  massSF = 1.;
  m_top = -1.;
  m_antitop = -1.;
  m_bqq = -1.;
  m_blv = -1.;
  m_hadr = -1;
  m_lept = -1;
  evWeight_hadr = 1.;
  evWeight_lept = 1.;
  evWeight_prod = 1.;
  bJetId.clear();
  bdiscrTop = -99.;
  bdiscrTop2 = -99.;
  tempbdiscr = -99.;
  labelB1 = -9999;
  labelB2 = -9999;
  massHadTopQ = 0.01;
  massLepTopQ = 0.01;
}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearVars();
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

void CalculateAverageTopMass()
{
  ifstream fileIn;
  streampos currentPosition;
  ofstream fileOut;
  
  string pathInput = "averageMass/";
  string outputFileName = pathInput+"averageMass_"+dateString+".txt";
  fileOut.open(outputFileName.c_str());
  string inputFiles[] = {"genp_matched", "genj_matched", "recoBKF", "recoAKF", "recoAKF_matched"};
  int nInputs = sizeof(inputFiles)/sizeof(inputFiles[0]);
  
  char dataLine[1024];
  int nEntries, eventId;
  double massTop, sumTop, meanTop;
  
  for (int iFile = 0; iFile < nInputs; iFile++)
  {
    nEntries = 0; sumTop = 0.; meanTop = 0.;
    
    string inputFileName = pathInput+"mass_"+inputFiles[iFile]+"_"+dateString+".txt";
    fileIn.open(inputFileName.c_str());
    cout << "Opening " << inputFileName << "..." << endl;
    fileIn.getline(dataLine,sizeof(dataLine));
    currentPosition = fileIn.tellg();
    
    /// Loop over input file
    while ( fileIn.good() )
    {
      eventId = -1; massTop = 0.;
      
      fileIn.seekg(currentPosition);
      fileIn.getline(dataLine,sizeof(dataLine));
      istringstream iss(dataLine);
      iss >> eventId >> massTop;
      currentPosition = fileIn.tellg();
      
      nEntries++;
      sumTop += massTop;
    }
    
    /// Calculate mean
    meanTop = sumTop/((double)nEntries);
    
    fileOut << left << setw(18) << inputFiles[iFile];
    cout.setf(ios::fixed,ios::floatfield);
    fileOut << "   " << fixed << showpoint << setprecision(3) << meanTop << endl;
    
    /// Close input file
    fileIn.close();
  }
  
  fileOut.close();
  cout << "Average masses can be found in " << outputFileName << endl;
}
