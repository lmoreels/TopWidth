#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
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


bool test = true;
bool testHistos = false;
bool testTTbarOnly = true;
bool calculateResolutionFunctions = false;
bool calculateAverageMass = false;
bool calculateLikelihood = false;
bool useToys = false;
bool applyLeptonSF = true;
bool applyPU = true;
bool applyJER = true;
bool applyJEC = true;
bool applyBTagSF = true;
bool applyNloSF = false;

bool applyWidthSF = false;
float scaleWidth = 0.66;

string systStr = "nominal";
string whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos ) return "160812";
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

int nofMatchedEvents = 0;
int nofHadrMatchedEvents = 0;
int nofChi2First4 = 0;
int nofCorrectlyMatched_chi2 = 0;
int nofNotCorrectlyMatched_chi2 = 0;
double Luminosity = 9999.;



///  Working points for b tagging  // Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
float CSVv2Loose =  0.460;
float CSVv2Medium = 0.800;
float CSVv2Tight = 0.935;

/// Top width
double genTopWidth = 1.33; // from fit
double genTopMass = 172.5; // from fit

// Temporarily, until calculated from TTbar sample
float chi2WMass = 80.385;
float sigmaChi2WMass = 10;
float chi2TopMass = 172.5; //180.0; //from mtop mass plot: 167.0
float sigmaChi2TopMass = 40;

/// Likelihood function
const int nCP =  593345;
const int nWP = 1085957;
const int nUP = 1710501; // ttbar only: 1648250;

// normalisation factor from integral, not from fit
const double mu_CP = 0.9984, sigma_CP = 0.08913, r_CP = 2.47, norm_CP = 0.002558;
const double alpha_WP = -0.3614, n_WP = 20, sigma_WP = 0.1278, mu_WP = 0.7436, norm_WP = 0.004463;
const double alpha_UP = -0.3639, n_UP = 20, sigma_UP = 0.1439, mu_UP = 0.7167, norm_UP = 0.003993;
const double norm_comb = 1.; //0.956039;

/// Average top mass
// TT gen match, TT reco match, TT reco wrongMatch WP/UP, TT reco noMatch, TT reco wrongPerm, TT reco wrongPerm W Ok, TT reco wrongPerm W Not Ok, TT reco, ST_t_top reco, ST_t_antitop reco, ST_tW_top reco, ST_tW_antitop reco, DYJets reco, WJets reco, data reco, all MC reco, all samples reco (data+MC) 
/// also background in CP/WP/UP cats (unlike name suggests)
//  no cuts [[nominal]]
float aveTopMass[] = {166.933, 168.186, 202.938, 210.914, 190.375, 204.165, 182.950, 196.562, 240.086, 232.843, 219.273, 221.202, 213.004, 200.023, 199.332, 196.855, 196.877};
//  with cut on dR
//float aveTopMass[] = {166.933, 168.193, 201.842, 209.812, 189.347, 202.179, 182.512, 195.553, 238.844, 232.204, 218.916, 220.941, 212.830, 199.355, 198.380, 195.846, 195.870};
/// TT only for cats
/// no cut on chi2
//float aveTopMass[] = {166.933, 168.186, 202.651, 210.589, 190.375, 204.165, 182.950, 196.562, 240.086, 232.843, 219.273, 221.202, 213.004, 200.023, 199.332, 196.855};
/// cut: chi2 < 2
//float aveTopMass[] = {166.933, 167.531, 192.382, 183.662, 200.401, 174.695, 190.354, 217.753, 214.542, 210.731, 212.884, 200.177, 184.198, 192.049, 190.539};
/// cut: chi2 < 4
//float aveTopMass[] = {166.933, 167.876, 195.023, 184.251, 201.242, 175.062, 192.004, 221.724, 218.145, 212.657, 214.665, 201.905, 189.114, 193.777, 192.206};
/// cut: chi2 < 5
//float aveTopMass[] = {166.933, 167.963, 195.902, 184.367, 201.420, 175.136, 192.538, 223.597, 219.594, 213.119, 215.540, 202.621, 189.915, 194.406, 192.748};

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

string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();
void GetMetaData(TTree* tree, bool isData);
void InitTree(TTree* tree, bool isData);
void InitMSPlots();
void InitHisto1D();
void InitHisto2D();
void TruthMatching(vector<TLorentzVector> partons, vector<TLorentzVector> selectedJets, pair<unsigned int, unsigned int> *MCPermutation);
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearObjects();
void FillGeneralPlots(int d);
long GetNEvents(TTree* fChain, string var, bool isData);
void GetHLTFraction(double* fractions);
double BreitWigner(double topPT, double scale);
double eventWeightCalculator(double topPT, double scale);
Double_t voigt(Double_t *x, Double_t *par);
Double_t crysBall_WP(Double_t *x);
Double_t crysBall_UP(Double_t *x);
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
double reco_hadWMass, reco_hadTopMass, reco_hadTopPt;
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

/// KinFitter
TKinFitter* kFitter;
double kFitChi2 = 99.;
int nofAcceptedKFit = 0;

/// Likelihood
int nTot = 0;
double f_CP = 1./3., f_WP = 1./3., f_UP = 1./3.;
double widthArray[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4.};
string widthArrayStr[] = {"0p5", "0p55", "0p6", "0p65", "0p7", "0p75", "0p8", "0p85", "0p9", "0p95", "1", "1p05", "1p1", "1p15", "1p2", "1p25", "1p3", "1p35", "1p4", "1p45", "1p5", "1p55", "1p6", "1p65", "1p7", "1p75", "1p8", "1p85", "1p9", "1p95", "2", "2p25", "2p5", "2p75", "3", "3p25", "3p5", "3p75", "4"};
const int nWidths = sizeof(widthArray)/sizeof(widthArray[0]);
double gammaArray[nWidths];

double loglike[nWidths] = {0};
double loglike_pd[10][nWidths] = {{0}};
double loglike2[nWidths] = {0};
double loglike2_pd[10][nWidths] = {{0}};
double loglike_per_evt[nWidths] = {0};
double loglike_onlyGoodEvts[nWidths] = {0};
double fakelike_CP[nWidths] = {0};
double fakelike_CP_per_evt[nWidths] = {0};
double fakelike_onlyGoodEvts[nWidths] = {0};
double fakelike_CP_Res[nWidths] = {0};
double fakelike_CP_Res2[nWidths] = {0};

bool isGoodLL = false;
int nofGoodEvtsLL[10] = {0};
int nofBadEvtsLL[10] = {0};
double maxMtDivAveMt = 0., minMtDivAveMt = 9999.;

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
  }
  
  //string pathOutput = "test/";
  string pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (! test && ! calculateAverageMass)
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
  if (calculateLikelihood) cout << "Calculating -loglikelihood values..." << endl;
  
  
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
    if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
      Luminosity = datasets[d]->EquivalentLumi();
    
    if ( dataSetName.find("QCD") == 0 )
      datasets[d]->SetColor(kYellow);
    if ( dataSetName.find("TT") == 0 )
    {
      datasets[d]->SetTitle("t#bar{t}");
      datasets[d]->SetColor(kRed+1);
    }
    //if ( dataSetName.find("TTbarJets_Other") == 0 ) datasets[d]->SetColor(kRed-7);
    if ( dataSetName.find("WJets") == 0 )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      datasets[d]->SetColor(kGreen-3);
    }
    if ( dataSetName.find("ZJets") == 0 || dataSetName.find("DY") == 0 )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      datasets[d]->SetColor(kAzure-2);
      //datasets[d]->SetColor(kMagenta);
    }
    if ( dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") == 0 )
    {
      datasets[d]->SetTitle("ST");
      datasets[d]->SetColor(kBlue-2);
      //if ( dataSetName.find("tW") == 0 )
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
  KinFitter *kf = new KinFitter("PlotsForResolutionFunctions_testFit.root");
  
  if (! test && ! calculateAverageMass)
  {
    InitMSPlots();
    InitHisto1D();
    InitHisto2D();
  }
  
  vJER.clear(); vJES.clear(); vPU.clear();
  
  if (calculateLikelihood)
  {
    txtLogLike.open(("likelihood_per_event_"+dateString+".txt").c_str());
    txtLogLike << "## -Log(likelihood) values per event" << endl;
    txtLogLike << "#  Widths : ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      txtLogLike << widthArray[iWidth] << "  ";
    }
    txtLogLike << endl;
    
    txtLogLikeTest.open(("likelihood_per_event_test_"+dateString+".txt").c_str());
    txtLogLikeTest << "## ievt      recoTopMass     M_lb    dR(b,lep)    dR(b_h,lep)" << endl;
    txtLogLikeTest << "#  Gammas : ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      gammaArray[iWidth] = widthToGammaTranslation(&widthArray[iWidth]);
      txtLogLikeTest << gammaArray[iWidth] << "  ";
    }
    
    /// Fraction of events that is CP, WP and UP
    nTot = nCP + nWP + nUP;
    f_CP = (double)nCP/(double)nTot;
    f_WP = (double)nWP/(double)nTot;
    f_UP = (double)nUP/(double)nTot;
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
    if (testTTbarOnly && (datasets[d]->Name()).find("TT") == std::string::npos ) continue;
    
    clock_t startDataSet = clock();
    
    ClearMetaData();
    
    dataSetName = datasets[d]->Name();
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
    }
    
    isData = false;
    if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
    {
      isData = true;
    }
    
    if (calculateAverageMass)
    {
//       if (dataSetName.find("TT") == 0 )
//       {
//         txtMassGenMatched.open(("averageMass/mass_gen_matched_TT_"+dateString+".txt").c_str());
//         txtMassRecoCP.open(("averageMass/mass_reco_matched_TT_"+dateString+".txt").c_str());
//         txtMassRecoWPUP.open(("averageMass/mass_reco_notCorrectMatch_TT_"+dateString+".txt").c_str());
//         txtMassRecoUP.open(("averageMass/mass_reco_notMatched_TT_"+dateString+".txt").c_str());
//         txtMassRecoWP.open(("averageMass/mass_reco_wrongPerm_TT_"+dateString+".txt").c_str());
//         txtMassRecoWPWOk.open(("averageMass/mass_reco_wrongPerm_WOk_TT_"+dateString+".txt").c_str());
//         txtMassRecoWPWNotOk.open(("averageMass/mass_reco_wrongPerm_WNotOk_TT_"+dateString+".txt").c_str());
//       }
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
    if (test || testHistos) endEvent = 21;
    for (int ievt = 0; ievt < endEvent; ievt++)
    {
      ClearObjects();
      
      
      if (ievt%10000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
      
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
      
      /// Fill objects
      muon.SetPtEtaPhiE(muon_pt[0], muon_eta[0], muon_phi[0], muon_E[0]);
      selectedLepton.push_back(muon);
      
      for (int iJet = 0; iJet < nJets; iJet++)
      {
        jet.Clear();
        jet.SetPtEtaPhiE(jet_pt[iJet], jet_eta[iJet], jet_phi[iJet], jet_E[iJet]);
        selectedJets.push_back(jet);
      }
      
      for (int iJet = 0; iJet < selectedJets.size(); iJet++)
      {
        if ( jet_bdiscr[iJet] > CSVv2Medium )
        {
          selectedBJets.push_back(selectedJets[iJet]);
          bJetId.push_back(iJet);  /// selectedBJets[j] = selectedJets[bJetId[j]]
        }
      }
      //std::sort(selectedBJets.begin(),selectedBJets.end(),HighestPt());  // already the case
      
      
      
      
      /////////////////////////////
      ///  JET PARTON MATCHING  ///
      /////////////////////////////
      
      //if ( dataSetName.find("TT") == 0 || dataSetName.find("ST") == 0 )  // no matches for ST
      if ( dataSetName.find("TT") == 0 )
      {
        for (int iMC = 0; iMC < nMCParticles; iMC++)
        {
          mcpart.Clear();
          mcpart.SetPtEtaPhiE(mc_pt[iMC], mc_eta[iMC], mc_phi[iMC], mc_E[iMC]);
          mcParticles.push_back(mcpart);
        }
        
        
        for (unsigned int i = 0; i < mcParticles.size(); i++)
        {
          if (verbose > 4)
            cout << setw(3) << right << i << "  Status: " << setw(2) << mc_status[i] << "  pdgId: " << setw(3) << mc_pdgId[i] << "  Mother: " << setw(4) << mc_mother[i] << "  Granny: " << setw(4) << mc_granny[i] << "  Pt: " << setw(7) << left << mc_pt[i] << "  Eta: " << mc_eta[i] << endl;
          
          
          if ( (mc_status[i] > 1 && mc_status[i] <= 20) || mc_status[i] >= 30 ) continue;  /// Final state particle or particle from hardest process
          
          
          if ( mc_pdgId[i] == pdgID_top && mc_granny[i] == 2212 )
            topQuark = i;
          else if( mc_pdgId[i] == -pdgID_top && mc_granny[i] == 2212 )
            antiTopQuark = i;
          
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
          
          if (all4PartonsMatched)
          {
            
            for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
            {
              /// MCPermutation[i].first  = jet number
              /// MCPermutation[i].second = parton number
              /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
              
              partonsMatched.push_back(partons[MCPermutation[iMatch].second]);
              jetsMatched.push_back(selectedJets[MCPermutation[iMatch].first]);
            }
            
            if (calculateResolutionFunctions)
            {
              rf->fillJets(partonsMatched, jetsMatched);
              
              if (muonmatched) rf->fillMuon(mcParticles[genmuon], selectedLepton[0]);
              //if (electronmatched) rf->fillElectron(...)

            }  // end rf
          }
          
          /// Plot variables for matched events
          if (hadronicTopJetsMatched)
          {  
            float matchedWMass_reco = (jetsMatched[0] + jetsMatched[1]).M();
            float matchedTopMass_reco = (jetsMatched[0] + jetsMatched[1] + jetsMatched[2]).M();
            float matchedTopMass_gen = (partonsMatched[0] + partonsMatched[1] + partonsMatched[2]).M();
            
            if (calculateAverageMass) txtMassGenMatched << ievt << "  " << matchedWMass_reco << "  " << matchedTopMass_reco << endl;
            
            if (! test && ! calculateAverageMass)
            {
              histo1D["W_mass_reco_matched"]->Fill(matchedWMass_reco, widthSF);
              histo1D["top_mass_reco_matched"]->Fill(matchedTopMass_reco, widthSF);
              histo1D["top_mass_gen_matched"]->Fill(matchedTopMass_gen, widthSF);
              
              histo1D["mTop_div_aveMTop_TT_matched_jets"]->Fill(matchedTopMass_reco/aveTopMass[0]);
              
              if ( hadronicTopJetsMatched_MCdef_ )
              {
                histo1D["W_mass_reco_first4matched"]->Fill(matchedWMass_reco, widthSF);
                histo1D["top_mass_reco_first4matched"]->Fill(matchedTopMass_reco, widthSF);
                histo1D["top_mass_gen_first4matched"]->Fill(matchedTopMass_gen, widthSF);
              }
              if (hasExactly4Jets)
              {
                histo1D["W_mass_reco_matched_4jets"]->Fill(matchedWMass_reco, widthSF);
                histo1D["top_mass_reco_matched_4jets"]->Fill(matchedTopMass_reco, widthSF);
                histo1D["top_mass_gen_matched_4jets"]->Fill(matchedTopMass_gen, widthSF);
              }
              
              if (all4PartonsMatched && muonmatched)
              {
                float matchedMlb_corr = (selectedLepton[0] + jetsMatched[3]).M();  // lept b
                float matchedMlb_wrong = (selectedLepton[0] + jetsMatched[2]).M();  // hadr b
                float matchedTtbarMass_corr = matchedTopMass_reco + matchedMlb_corr;
                float matchedTtbarMass_wrong = matchedTopMass_reco + matchedMlb_wrong;
                double matchedDRLepB_corr = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[3]);  // lept b
                double matchedDRLepB_wrong = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[2]);  // hadr b
                
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
                
                if (hasExactly4Jets)
                {
                  histo1D["mlb_matched_corr_4jets"]->Fill(matchedMlb_corr);
                  histo1D["mlb_matched_wrong_4jets"]->Fill(matchedMlb_wrong);
                  histo1D["ttbar_mass_matched_corr_4jets"]->Fill(matchedTtbarMass_corr);
                  histo1D["ttbar_mass_matched_wrong_4jets"]->Fill(matchedTtbarMass_wrong);
                  histo1D["dR_lep_b_matched_corr_4jets"]->Fill(matchedDRLepB_corr);
                  histo1D["dR_lep_b_matched_wrong_4jets"]->Fill(matchedDRLepB_wrong);
                }
                
              }  // end all4PartonsMatched && muonMatched
              
            }  // end fill plots

          }  // end hadronicTopJetsMatched
        
        }  // end doMatching
        
        
      }  // end if TT
      
      
      
      /////////////////////////////////////
      ///  Reconstruction of Top quark  ///
      /////////////////////////////////////
      
      int labelsReco[3] = {-9999, -9999, -9999};		// 0,1 = light jets, 2 = hadronic b-jet.
      float recoWMass, chi2;
      float smallestChi2 = 9999.;
      double deltaR;
      double minDeltaR = 9999.;
      
      for (int ijet = 0; ijet < selectedJets.size(); ijet++)
      {
        for (int jjet = ijet+1; jjet < selectedJets.size(); jjet++)
        {
          recoWMass = (selectedJets[ijet] + selectedJets[jjet]).M();
          chi2 = pow( (recoWMass - chi2WMass)/sigmaChi2WMass, 2);
          
          if ( chi2 < smallestChi2)
          {
            smallestChi2 = chi2;
            labelsReco[0] = ijet;
            labelsReco[1] = jjet;
          }
        }
      }
      //if ( smallestChi2 > 2 ) continue;
      
      WCandidate = selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]];
      for (int kjet = 0; kjet < selectedJets.size(); kjet++)
      {
        if ( kjet == labelsReco[0] || kjet == labelsReco[1] ) continue;
        if ( selectedJets.size() > 4 && jet_bdiscr[kjet] < CSVv2Medium ) continue; // if exactly 4 jets, use all remaining jets; if more than 4 jets, only use remaining b jets
        
        deltaR = ROOT::Math::VectorUtil::DeltaR(selectedJets[kjet], WCandidate);
        
        if (deltaR < minDeltaR)
        {
          minDeltaR = deltaR;
          labelsReco[2] = kjet;
        }
      }
      if ( selectedJets.size() > 4 && labelsReco[2] == -9999 )
      {
        for (int kjet = 0; kjet < selectedJets.size(); kjet++)
        {
          if ( kjet == labelsReco[0] || kjet == labelsReco[1] ) continue;

          deltaR = ROOT::Math::VectorUtil::DeltaR(selectedJets[kjet], WCandidate);

          if (deltaR < minDeltaR)
          {
            minDeltaR = deltaR;
            labelsReco[2] = kjet;
          }
        }
      }
      
      ///-----Test KinFit------///
      cout << "Event " << setw(5) << right << ievt << "   ";
      kFitter = kf->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]], selectedJets[labelsReco[2]]);
      if ( kFitter->getStatus() == 0 )
      {
        kFitChi2 = kFitter->getS();
        cout << "Fit converged: Chi2 = " << kFitChi2 << endl;
      }
      //if ( kFitChi2 > 10. ) continue;
      nofAcceptedKFit++;
      
      
      
      if ( labelsReco[0] < 4 && labelsReco[1] < 4 && labelsReco[2] < 4 )
      {
        nofChi2First4++;
      }
      
      if ( labelsReco[0] != -9999 && labelsReco[1] != -9999 && labelsReco[2] != -9999 )
      {
        reco_hadWMass = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).M();
        reco_hadTopMass = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).M();
        reco_hadTopPt = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).Pt();
        
        /// Leptonic top mass
        double reco_Mlb_temp = 99999.;
        for (unsigned int i = 0; i < selectedJets.size(); i++)
        {
          if ( i == labelsReco[0] || i == labelsReco[1] || i == labelsReco[2] ) continue;
          
          reco_Mlb_temp = (selectedLepton[0] + selectedJets[i]).M();
          if ( jet_bdiscr[i] > CSVv2Medium && reco_Mlb_temp < reco_minMlb )
          {
            reco_minMlb = reco_Mlb_temp;
            reco_dRLepB_lep = ROOT::Math::VectorUtil::DeltaR( selectedJets[i], selectedLepton[0] );
          }
          if ( jet_bdiscr[i] < CSVv2Medium && reco_Mlb_temp < reco_minMl_nonb )
          {
            reco_minMl_nonb = reco_Mlb_temp;
            labelMl_nonb = i;
          }
        }
        
        if ( reco_minMlb == 9999. )  // if no Mlb combination with b-tagged jet, take jet that is not b tagged
        {
          reco_minMlb = reco_minMl_nonb;
          reco_dRLepB_lep = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelMl_nonb], selectedLepton[0] );
        }
        reco_dRLepB_had = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
        
        reco_ttbarMass = reco_minMlb + reco_hadTopMass;
        
        
        /// Test cut on dR
//        if ( reco_dRLepB_lep > 3. && reco_dRLepB_had < 1.2 ) continue;  // skip event
        
        
        if (calculateAverageMass) txtMassReco << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
        
        if (calculateLikelihood)
        {
          if (isData)
          { 
            double tempAveMass = reco_hadTopMass/aveTopMass[14];
            for (int iWidth = 0; iWidth < nWidths; iWidth++)
            {
              loglike[iWidth] += logLikelihood(&tempAveMass, &gammaArray[iWidth]);
              loglike_pd[d][iWidth] += logLikelihood(&tempAveMass, &gammaArray[iWidth]);
            }
          }
          else
          {
            double tempAveMass = reco_hadTopMass/aveTopMass[15];  // d+6 or [15]? (All MC together)
            for (int iWidth = 0; iWidth < nWidths; iWidth++)
            {
              loglike[iWidth] += logLikelihood(&tempAveMass, &gammaArray[iWidth]);
              loglike_pd[d][iWidth] += logLikelihood(&tempAveMass, &gammaArray[iWidth]);
            }
          }
          double tempAveMass = reco_hadTopMass/aveTopMass[16];
          if ( tempAveMass > maxMtDivAveMt ) maxMtDivAveMt = tempAveMass;
          if ( tempAveMass < minMtDivAveMt ) minMtDivAveMt = tempAveMass;
          for (int iWidth = 0; iWidth < nWidths; iWidth++)
          {
            loglike_per_evt[iWidth] = logLikelihood(&tempAveMass, &gammaArray[iWidth]);
            loglike2[iWidth] += loglike_per_evt[iWidth];
            loglike2_pd[d][iWidth] += loglike_per_evt[iWidth];
          }
          
          /// make loglikelihood only with events that have minimum
          isGoodLL = false;
          for (int iWidth = 1; iWidth < nWidths-1; iWidth++)
          {
            if ( loglike_per_evt[0] > loglike_per_evt[iWidth] && loglike_per_evt[iWidth] < loglike_per_evt[nWidths-1] )
            {
              isGoodLL = true;
              break;
            }
          }
          if (isGoodLL)
          {
            nofGoodEvtsLL[d]++;
            for (int iWidth = 0; iWidth < nWidths; iWidth++)
            {
              loglike_onlyGoodEvts[iWidth] += loglike_per_evt[iWidth];
            }
          }
          else
          {
            nofBadEvtsLL[d]++;
          }
          
          /// write debug file
          if ( ( dataSetName.find("data") == 0 && ievt%1000 == 0 ) 
               || ( dataSetName.find("TT") == 0 && ievt%100000 == 0 )
               || ( dataSetName.find("data") != 0 && dataSetName.find("TT") != 0 && ievt%100 == 0 ) )
          {
            txtLogLike << ievt << "  ";
            for (int iWidth = 0; iWidth < nWidths; iWidth++)
            {
              txtLogLike << loglike_per_evt[iWidth] << ",  ";
              
            }
            txtLogLike << reco_hadTopMass << endl;
          }
        }
        
        
        //Fill histos
        if (! test && ! calculateAverageMass)
        {
          /// Hadronic part only
          MSPlot["Chi2_value"]->Fill(smallestChi2, datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_W_mass"]->Fill(reco_hadWMass, datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_hadTop_mass"]->Fill(reco_hadTopMass, datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_hadTop_pT"]->Fill(reco_hadTopPt, datasets[d], true, Luminosity*scaleFactor);
          
          MSPlot["Chi2_mTop_div_aveMTopMatch"]->Fill(reco_hadTopMass/aveTopMass[0], datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_mTop_div_aveMTopTTChi2"]->Fill(reco_hadTopMass/aveTopMass[7], datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_mTop_div_aveMTopAllChi2"]->Fill(reco_hadTopMass/aveTopMass[15], datasets[d], true, Luminosity*scaleFactor);
          
          if (isData)
          {
            histo1D[("mTop_div_aveMTop_"+dataSetName).c_str()]->Fill(reco_hadTopMass/aveTopMass[14]);
            MSPlot["Chi2_mTop_div_aveMTop_stack"]->Fill(reco_hadTopMass/aveTopMass[14], datasets[d], true, Luminosity*scaleFactor*widthSF);
          }
          else
          {
            histo1D[("mTop_div_aveMTop_"+dataSetName).c_str()]->Fill(reco_hadTopMass/aveTopMass[d+6]);
            MSPlot["Chi2_mTop_div_aveMTop_stack"]->Fill(reco_hadTopMass/aveTopMass[d+6], datasets[d], true, Luminosity*scaleFactor*widthSF);
            histo1D["mTop_div_aveMTop_allMCReco"]->Fill(reco_hadTopMass/aveTopMass[15]);
          }
          
          if (hasExactly4Jets)
          {
            MSPlot["Chi2_hadTop_mass_4jets"]->Fill(reco_hadTopMass, datasets[d], true, Luminosity*scaleFactor);
          }
          
          if ( dataSetName.find("TT") == 0 )
          {
            histo1D["Chi2_W_mass_reco"]->Fill(reco_hadWMass, widthSF);
            histo1D["Chi2_top_mass_reco"]->Fill(reco_hadTopMass, widthSF);
            if (hasExactly4Jets)
            {
              histo1D["Chi2_top_mass_reco_4jets"]->Fill(reco_hadTopMass, widthSF);
            }
          }
          
          
          /// Plots with lepton
          MSPlot["Chi2_mlb"]->Fill(reco_minMlb, datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_ttbar_mass"]->Fill(reco_ttbarMass, datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_dR_lep_b"]->Fill(reco_dRLepB_lep, datasets[d], true, Luminosity*scaleFactor);
          if ( dataSetName.find("TT") == 0 )
          {
            histo1D["dR_lep_b_unmatched_chi2"]->Fill(reco_dRLepB_lep, widthSF);
          }
          
          if (hasExactly4Jets)
          {
            MSPlot["Chi2_mlb_4jets"]->Fill(reco_minMlb, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["Chi2_ttbar_mass_4jets"]->Fill(reco_ttbarMass, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["Chi2_dR_lep_b_4jets"]->Fill(reco_dRLepB_lep, datasets[d], true, Luminosity*scaleFactor);
          }
          
          if ( reco_minMlb < 200 )
          {
            MSPlot["Chi2_hadTop_mass_mlb_cut"]->Fill(reco_hadTopMass, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["Chi2_ttbar_mass_mlb_cut"]->Fill(reco_ttbarMass, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["Chi2_dR_lep_b_mlb_cut"]->Fill(reco_dRLepB_lep, datasets[d], true, Luminosity*scaleFactor);
            if (hasExactly4Jets)
            {
              MSPlot["Chi2_hadTop_mass_4jets_mlb_cut"]->Fill(reco_hadTopMass, datasets[d], true, Luminosity*scaleFactor);
              MSPlot["Chi2_ttbar_mass_4jets_mlb_cut"]->Fill(reco_ttbarMass, datasets[d], true, Luminosity*scaleFactor);
              MSPlot["Chi2_dR_lep_b_4jets_mlb_cut"]->Fill(reco_dRLepB_lep, datasets[d], true, Luminosity*scaleFactor);
            }
          }
          
        }  // end fill plots
        
      }  // end labels
      
      
      
      ///////////////////////////////////
      ///  CHECK MATCHED COMBINATION  ///
      ///////////////////////////////////
      
      ///
      // 3 possibilities:
      // - correct top match: 3 jets selected with reco method correspond to the 3 matched jets (n.b. this is also true when the jets originating from the W boson and the b jet do not exactly correspond to the matched jets, because we are only interested in the reconstructed top quark.)
      // - wrong permutation: the correct jet combination exists in the selected jets, but is not chosen by the reco method.
      // - wrong (no) match:  the correct jet combination does not exist in the selected jets (e.g. when one jet is not selected.)
      
      
      //if ( dataSetName.find("TT") == 0 )
      if (! isData)
      {
        if (! applyWidthSF ) widthSF = 1.;
        //if (test) cout << "Reco done, about to check match..." << endl;
        
        if (hadronicTopJetsMatched)
        { 
          if ( test && reco_hadTopMass == -1 )
          {
            cout << "Event: " << ievt << "; mass for width: " << massForWidth << ";  widthSF " << widthSF << endl;
            cout << "Size selectedJets:  " << selectedJets.size() << endl;
            cout << "Labels MCParticles: " << MCPermutation[0].first << "  " << MCPermutation[1].first << "  " << MCPermutation[2].first << endl;
            cout << "Labels jets:        " << labelsReco[0] << "  " << labelsReco[1] << "  " << labelsReco[2] << endl;
          }
          
          if ( ( labelsReco[0] == MCPermutation[0].first || labelsReco[0] == MCPermutation[1].first || labelsReco[0] == MCPermutation[2].first ) && ( labelsReco[1] == MCPermutation[0].first || labelsReco[1] == MCPermutation[1].first || labelsReco[1] == MCPermutation[2].first ) && ( labelsReco[2] == MCPermutation[0].first || labelsReco[2] == MCPermutation[1].first || labelsReco[2] == MCPermutation[2].first ) )  // correct jets for top quark
          {
            nofCorrectlyMatched_chi2++;
            if (calculateAverageMass) txtMassRecoCP << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
            else if (! test)
            {
              histo1D["mTop_div_aveMTop_TT_reco_CP"]->Fill(reco_hadTopMass/aveTopMass[1], widthSF);
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
            }
            
            if (calculateLikelihood)
            {
              float sumJetEnergies = selectedJets[MCPermutation[0].first].E() + selectedJets[MCPermutation[1].first].E() + selectedJets[MCPermutation[2].first].E();
              float sumJetMinusPartonEnergies = sumJetEnergies - mcParticles[partonId[MCPermutation[0].second]].E() - mcParticles[partonId[MCPermutation[1].second]].E() - mcParticles[partonId[MCPermutation[2].second]].E();
              
              double tempAveMass = reco_hadTopMass/aveTopMass[16];
              for (int iWidth = 0; iWidth < nWidths; iWidth++)
              {
                fakelike_CP_per_evt[iWidth] = fakeLikelihood(&tempAveMass, &gammaArray[iWidth]);
                fakelike_CP[iWidth] += fakelike_CP_per_evt[iWidth];
                if ( fabs(sumJetMinusPartonEnergies) < 0.005 * sumJetEnergies )
                  fakelike_CP_Res[iWidth] += fakelike_CP_per_evt[iWidth];
                if ( fabs(selectedJets[MCPermutation[0].first].E()-mcParticles[partonId[MCPermutation[0].second]].E()) < 0.005 * selectedJets[MCPermutation[0].first].E() 
                    && fabs(selectedJets[MCPermutation[1].first].E()-mcParticles[partonId[MCPermutation[1].second]].E()) < 0.005 * selectedJets[MCPermutation[1].first].E()
                    && fabs(selectedJets[MCPermutation[2].first].E()-mcParticles[partonId[MCPermutation[2].second]].E()) < 0.005 * selectedJets[MCPermutation[2].first].E() )
                  fakelike_CP_Res2[iWidth] += fakelike_CP_per_evt[iWidth];
              }

              /// make loglikelihood only with events that have minimum
              isGoodLL = false;
              for (int iWidth = 1; iWidth < nWidths-1; iWidth++)
              {
                if ( fakelike_CP_per_evt[0] > fakelike_CP_per_evt[iWidth] && fakelike_CP_per_evt[iWidth] < fakelike_CP_per_evt[nWidths-1] )
                {
                  isGoodLL = true;
                  break;
                }
              }
              if (isGoodLL)
              {
                for (int iWidth = 0; iWidth < nWidths; iWidth++)
                {
                  fakelike_onlyGoodEvts[iWidth] += fakelike_CP_per_evt[iWidth];
                }
              }
            }
            
            /// write debug file
            if ( calculateLikelihood && ! isGoodLL /*&& ( ( dataSetName.find("data") == 0 && nofGoodEvtsLL[d]%500 == 0 ) 
                 || ( dataSetName.find("TT") == 0 && nofGoodEvtsLL[d]%50000 == 0 )
                 || ( dataSetName.find("data") != 0 && dataSetName.find("TT") != 0 && nofGoodEvtsLL[d]%50 == 0 ) )*/ )
            {
              if (test)
              {
                txtLogLikeTest << setw(8) << right << ievt << "  CP  ";
                txtLogLikeTest << reco_hadTopMass << "  " << reco_minMlb << "  " << reco_dRLepB_lep << "  " << reco_dRLepB_had << endl;
              }
              histo1D["debugLL_hadr_top_mass_reco"]->Fill(reco_hadTopMass);
              histo1D["debugLL_minMlb_reco"]->Fill(reco_minMlb);
              histo1D["debugLL_ttbar_mass_reco"]->Fill(reco_ttbarMass);
              histo1D["debugLL_dR_lep_b_lep_reco"]->Fill(reco_dRLepB_lep);
              histo1D["debugLL_dR_lep_b_had_reco"]->Fill(reco_dRLepB_had);
            }
            
          }  // end corr match
          else  // wrong permutation
          {
            nofNotCorrectlyMatched_chi2++;
            if (calculateAverageMass)
            {
              txtMassRecoWP << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
              txtMassRecoWPUP << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
            }
            else if (! test)
            {
              histo1D["mTop_div_aveMTop_TT_reco_WP"]->Fill(reco_hadTopMass/aveTopMass[4], widthSF);
              histo1D["mTop_div_aveMTop_TT_reco_WPUP"]->Fill(reco_hadTopMass/aveTopMass[2], widthSF);
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
            }
            
            /// write debug file
            if ( test && calculateLikelihood && ! isGoodLL /*&& ( ( dataSetName.find("data") == 0 && nofGoodEvtsLL[d]%500 == 0 ) 
                 || ( dataSetName.find("TT") == 0 && nofGoodEvtsLL[d]%50000 == 0 )
                 || ( dataSetName.find("data") != 0 && dataSetName.find("TT") != 0 && nofGoodEvtsLL[d]%50 == 0 ) )*/ )
            {
              txtLogLikeTest << setw(8) << right << ievt << "  WP  ";
              txtLogLikeTest << reco_hadTopMass << "  " << reco_minMlb << "  " << reco_dRLepB_lep << "  " << reco_dRLepB_had << endl;
            }

            if ( ( labelsReco[0] == MCPermutation[0].first || labelsReco[0] == MCPermutation[1].first || labelsReco[0] == MCPermutation[2].first ) && ( labelsReco[1] == MCPermutation[0].first || labelsReco[1] == MCPermutation[1].first || labelsReco[1] == MCPermutation[2].first ) )
            {
              if (calculateAverageMass) txtMassRecoWPWOk << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
              else if (! test)
                histo1D["mTop_div_aveMTop_TT_reco_WP_WOk"]->Fill(reco_hadTopMass/aveTopMass[5], widthSF);
            }
            else
            {
              if (calculateAverageMass) txtMassRecoWPWNotOk << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
              else if (! test)
                histo1D["mTop_div_aveMTop_TT_reco_WP_WNotOk"]->Fill(reco_hadTopMass/aveTopMass[6], widthSF);
            }
          }  // end wrong perm
        }  // end hadrTopMatch
        else  // no match
        {
          if (calculateAverageMass)
          {
            txtMassRecoUP << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
            txtMassRecoWPUP << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
          }
          else if (! test)
          {
            histo1D["mTop_div_aveMTop_TT_reco_UP"]->Fill(reco_hadTopMass/aveTopMass[3], widthSF);
            histo1D["mTop_div_aveMTop_TT_reco_WPUP"]->Fill(reco_hadTopMass/aveTopMass[2], widthSF);
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
          }
          
          /// write debug file
          if ( test && calculateLikelihood && ! isGoodLL /*&& ( ( dataSetName.find("data") == 0 && nofGoodEvtsLL[d]%500 == 0 ) 
               || ( dataSetName.find("TT") == 0 && nofGoodEvtsLL[d]%50000 == 0 )
               || ( dataSetName.find("data") != 0 && dataSetName.find("TT") != 0 && nofGoodEvtsLL[d]%50 == 0 ) )*/ )
          {
            txtLogLikeTest << setw(8) << right << ievt << "  UP  ";
            txtLogLikeTest << reco_hadTopMass << "  " << reco_minMlb << "  " << reco_dRLepB_lep << "  " << reco_dRLepB_had << endl;
          }
            
        }  // end no match
        //if (test) cout << "checked match" << endl;
      }  // end TT / ! isData
      
      if (! calculateAverageMass && ! test)
      {
        histo1D["mTop_div_aveMTop_bkgd"]->Fill(reco_hadTopMass/aveTopMass[15]);
      }
      
      
      ////////////////////
      ///  Make plots  ///
      ////////////////////
      
      if (! test && ! calculateAverageMass)
      {
        FillGeneralPlots(d);
      }
      
    }  // end loop events
    
    
    cout << endl;
    cout << "Number of chi2 made with the 4 most energetic jets: " << nofChi2First4 << " (" << 100*((float)nofChi2First4/(float)endEvent) << "%)" << endl;
    cout << "Number of events accepted by kinFitter: " << nofAcceptedKFit << " (" << 100*((float)nofAcceptedKFit/(float)endEvent) << "%)" << endl;
    
    //if ( dataSetName.find("TT") == 0 || dataSetName.find("ST") == 0 )
    if ( dataSetName.find("TT") == 0 )
    {
      cout << "Number of matched events: " << setw(8) << right << nofMatchedEvents << endl;
      cout << "Number of events with hadronic top matched: " << setw(8) << right << nofHadrMatchedEvents << endl;
      cout << "Correctly matched reconstructed events:     " << setw(8) << right << nofCorrectlyMatched_chi2 << endl;
      cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatched_chi2 << endl;
      if ( nofCorrectlyMatched_chi2 != 0 || nofNotCorrectlyMatched_chi2 != 0 )
        cout << "   ===> This means that " << 100*(float)nofCorrectlyMatched_chi2 / (float)(nofCorrectlyMatched_chi2 + nofNotCorrectlyMatched_chi2) << "% is correctly matched." << endl;
      
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
      
//       if (calculateAverageMass)
//       {
//         txtMassGenMatched.close();
//         txtMassRecoCP.close();
//         txtMassRecoWPUP.close();
//         txtMassRecoUP.close();
//         txtMassRecoWP.close();
//         txtMassRecoWPWOk.close();
//         txtMassRecoWPWNotOk.close();
//       }
      
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
  
  if (calculateLikelihood)
  {
    txtLogLike.close();
    txtLogLikeTest.close();
    
    cout << "Minimum top mass divided by average top mass: " << minMtDivAveMt << endl;
    cout << "Maximum top mass divided by average top mass: " << maxMtDivAveMt << endl;
    
    cout << endl << "likelihood values (all samples) : {";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << loglike[iWidth]/(1e+6);
      else cout << loglike[iWidth]/(1e+6) << ", ";
    }
    cout << "} *10^6 " << endl << "likelihood values (data-only) : ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << loglike_pd[0][iWidth];
      else cout << loglike_pd[0][iWidth] << ", ";
    }
    cout << endl << "likelihood values 2 (all samples) : {";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << loglike2[iWidth]/(1e+6);
      else cout << loglike2[iWidth]/(1e+6) << ", ";
    }
    cout << "} *10^6 " << endl << "likelihood values 2 (data-only) : ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << loglike2_pd[0][iWidth];
      else cout << loglike2_pd[0][iWidth] << ", ";
    }
    cout << endl << "likelihood values (only good events) : {";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << loglike_onlyGoodEvts[iWidth]/(1e+6);
      else cout << loglike_onlyGoodEvts[iWidth]/(1e+6) << ", ";
    }
    cout << "} *10^6 " << endl << "fake likelihood values (CP) : ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << fakelike_CP[iWidth]/(1e+6);
      else cout << fakelike_CP[iWidth]/(1e+6) << ", ";
    }
    cout << "} *10^6" << endl << "fake likelihood values (CP - only good events) : {";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << fakelike_onlyGoodEvts[iWidth]/(1e+6);
      else cout << fakelike_onlyGoodEvts[iWidth]/(1e+6) << ", ";
    }
    cout << "} *10^6 " << endl << "fake likelihood values (CP) with E_jet - E_q < 0.01*E_jet (sum): ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << fakelike_CP_Res[iWidth]/(1e+6);
      else cout << fakelike_CP_Res[iWidth]/(1e+6) << ", ";
    }
    cout << "} *10^6 " << endl << "fake likelihood values (CP) with E_jet - E_q < 0.01*E_jet (seperate): ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << fakelike_CP_Res2[iWidth]/(1e+6);
      else cout << fakelike_CP_Res2[iWidth]/(1e+6) << ", ";
    }
    cout << "} *10^6" << endl << "widths: ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << widthArray[iWidth];
      else cout << widthArray[iWidth] << ", ";
    }
    cout << endl << "gammas: ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      if ( iWidth == nWidths-1 ) cout << gammaArray[iWidth];
      else cout << gammaArray[iWidth] << ", ";
    }
    cout << endl << endl;
    
    /// Print output to file
    txtOutputLogLike.open("output_loglikelihood.txt");
    txtOutputLogLike << "Widths:      ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      txtOutputLogLike << setw(5) << right << widthArray[iWidth] << "  ";
    }
    txtOutputLogLike << "Gammas:      ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      txtOutputLogLike << setw(5) << right << gammaArray[iWidth] << "  ";
    }
    txtOutputLogLike << endl << "LLikelihood: ";
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
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
    float secs = time - mins*60;
    
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
  
  
  /// Chi2
  MSPlot["Chi2_value"] = new MultiSamplePlot(datasets, "Chi2_value", 800, 0, 200, "#chi^{2} value");
  MSPlot["Chi2_W_mass"] = new MultiSamplePlot(datasets, "Chi2_W_mass", 100, 0, 400, "M_{W} [GeV]");
  MSPlot["Chi2_hadTop_mass"] = new MultiSamplePlot(datasets, "Chi2_hadTop_mass", 80, 0, 800, "M_{t} [GeV]");
  MSPlot["Chi2_hadTop_mass_4jets"] = new MultiSamplePlot(datasets, "Chi2_hadTop_mass_4jets", 80, 0, 800, "M_{t} [GeV]");
  MSPlot["Chi2_hadTop_pT"] = new MultiSamplePlot(datasets, "Chi2_hadTop_pT", 80, 0, 800, "p_{T} [GeV]");
  
  MSPlot["Chi2_mlb"] = new MultiSamplePlot(datasets, "Chi2_mlb", 80, 0, 800, "M_{lb} [GeV]");
  MSPlot["Chi2_mlb_4jets"] = new MultiSamplePlot(datasets, "Chi2_mlb_4jets", 80, 0, 800, "M_{lb} [GeV]");

  MSPlot["Chi2_hadTop_mass_mlb_cut"] = new MultiSamplePlot(datasets, "Chi2_hadTop_mass_mlb_cut", 40, 0, 800, "M_{t} [GeV]");
  MSPlot["Chi2_hadTop_mass_4jets_mlb_cut"] = new MultiSamplePlot(datasets, "Chi2_hadTop_mass_4jets_mlb_cut", 40, 0, 800, "M_{t} [GeV]");

  MSPlot["Chi2_ttbar_mass"] = new MultiSamplePlot(datasets, "Chi2_ttbar_mass", 50, 0, 1000, "M_{t#bar{t}} [GeV]");
  MSPlot["Chi2_ttbar_mass_4jets"] = new MultiSamplePlot(datasets, "Chi2_ttbar_mass_4jets", 50, 0, 1000, "M_{t#bar{t}} [GeV]");
  MSPlot["Chi2_ttbar_mass_mlb_cut"] = new MultiSamplePlot(datasets, "Chi2_ttbar_mass_mlb_cut", 50, 0, 1000, "M_{t#bar{t}} [GeV]");
  MSPlot["Chi2_ttbar_mass_4jets_mlb_cut"] = new MultiSamplePlot(datasets, "Chi2_ttbar_mass_4jets_mlb_cut", 50, 0, 1000, "M_{t#bar{t}} [GeV]");

  MSPlot["Chi2_dR_lep_b"] = new MultiSamplePlot(datasets, "Chi2_dR_lep_b", 25, 0, 5, "#Delta R(l,b)");
  MSPlot["Chi2_dR_lep_b_4jets"] = new MultiSamplePlot(datasets, "Chi2_dR_lep_b_4jets", 25, 0, 5, "#Delta R(l,b)");
  MSPlot["Chi2_dR_lep_b_mlb_cut"] = new MultiSamplePlot(datasets, "Chi2_dR_lep_b_mlb_cut", 25, 0, 5, "#Delta R(l,b)");
  MSPlot["Chi2_dR_lep_b_4jets_mlb_cut"] = new MultiSamplePlot(datasets, "Chi2_dR_lep_b_4jets_mlb_cut", 25, 0, 5, "#Delta R(l,b)");
  
  
  MSPlot["Chi2_mTop_div_aveMTopMatch"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass (from matched events)", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
  MSPlot["Chi2_mTop_div_aveMTopTTChi2"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass (from chi2 of reco TT events)", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
  MSPlot["Chi2_mTop_div_aveMTopAllChi2"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass (from chi2 of reco events from all datasets)", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
  
  MSPlot["Chi2_mTop_div_aveMTop_stack"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass", 880, 0.2, 2.4, "M_{t}/<M_{t}>");
}

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  /// Chi2
  histo1D["Chi2_W_mass_reco"] = new TH1F("Chi2_W_mass_reco","Reconstructed hadronic W mass using a Chi2; M_{W} [GeV]", 80, 0, 800);
  histo1D["Chi2_top_mass_reco"] = new TH1F("Chi2_top_mass_reco","Reconstructed top mass using a Chi2; M_{t} [GeV]", 80, 0, 800);
  histo1D["Chi2_top_mass_reco_4jets"] = new TH1F("Chi2_top_mass_reco_4jets","Reconstructed top mass using a Chi2 (exactly 4 jets); M_{t} [GeV]", 80, 0, 800);
  
  /// Matching
  histo1D["W_mass_reco_matched"] = new TH1F("W_mass_reco_matched","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 125, 0, 250);
  histo1D["top_mass_reco_matched"] = new TH1F("top_mass_reco_matched","Reconstructed top mass of matched events; M_{t} [GeV]", 175, 50, 400);
  histo1D["top_mass_gen_matched"] = new TH1F("top_mass_gen_matched","Generated top mass of matched events; M_{t} [GeV]", 175, 0, 350);
  histo1D["mlb_matched_corr"]  = new TH1F("mlb_matched_corr","Reconstructed leptonic top mass using correctly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["mlb_matched_wrong"] = new TH1F("mlb_matched_wrong","Reconstructed leptonic top mass using wrongly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["ttbar_mass_matched_corr"] = new TH1F("ttbar_mass_matched_corr","Reconstructed mass of the top quark pair using correctly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong"] = new TH1F("ttbar_mass_matched_wrong","Reconstructed mass of the top quark pair using wrongly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["dR_lep_b_matched_corr"] = new TH1F("dR_lep_b_matched_corr","Delta R between the lepton and the leptonic b jet for matched events; #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_lep_b_matched_wrong"] = new TH1F("dR_lep_b_matched_wrong","Delta R between the lepton and the hadronic b jet for matched events; #Delta R(l,b_{had})", 25, 0, 5);
  
  histo1D["W_mass_reco_matched_4jets"] = new TH1F("W_mass_reco_matched_4jets","Reconstructed hadronic W mass of matched events with exactly 4 jets; M_{W} [GeV]", 125, 0, 250);
  histo1D["top_mass_reco_matched_4jets"] = new TH1F("top_mass_reco_matched_4jets","Reconstructed top mass of matched events with exactly 4 jets; M_{t} [GeV]", 175, 50, 400);
  histo1D["top_mass_gen_matched_4jets"] = new TH1F("top_mass_gen_matched_4jets","Generated top mass of matched events with exactly 4 jets; M_{t} [GeV]", 175, 0, 350);
  histo1D["mlb_matched_corr_4jets"] = new TH1F("mlb_matched_corr_4jets","Reconstructed leptonic top mass using correctly matched events with exactly 4 jets; M_{lb} [GeV]", 80, 0, 800);
  histo1D["mlb_matched_wrong_4jets"] = new TH1F("mlb_matched_wrong_4jets","Reconstructed leptonic top mass using wrongly matched events with exactly 4 jets; M_{lb} [GeV]", 80, 0, 800);
  histo1D["ttbar_mass_matched_corr_4jets"] = new TH1F("ttbar_mass_matched_corr_4jets","Reconstructed mass of the top quark pair using correctly matched events with exactly 4 jets; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong_4jets"] = new TH1F("ttbar_mass_matched_wrong_4jets","Reconstructed mass of the top quark pair using wrongly matched events with exactly 4 jets; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["dR_lep_b_matched_corr_4jets"] = new TH1F("dR_lep_b_matched_corr_4jets","Delta R between the lepton and the leptonic b jet for matched events with exactly 4 jets; #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_lep_b_matched_wrong_4jets"] = new TH1F("dR_lep_b_matched_wrong_4jets","Delta R between the lepton and the hadronic b jet for matched events with exactly 4 jets; #Delta R(l,b_{had})", 25, 0, 5);
  
  histo1D["W_mass_reco_first4matched"] = new TH1F("W_mass_reco_first4matched","Reconstructed hadronic W mass of events where 4 hardest jets are matched; M_{W} [GeV]", 125, 0, 250);
  histo1D["top_mass_reco_first4matched"] = new TH1F("top_mass_reco_first4matched","Reconstructed top mass of events where 4 hardest jets are matched; M_{t} [GeV]", 175, 50, 400);
  histo1D["top_mass_gen_first4matched"]= new TH1F("top_mass_gen_first4matched","Generated top mass of events where partons are matched to 4 hardest jets; M_{t} [GeV]", 175, 0, 350);
  
  
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
  histo1D["dR_lep_b_unmatched_chi2"] = new TH1F("dR_lep_b_unmatched_chi2","Minimal delta R between the lepton and a b jet (looking at b jets not in chi2); #Delta R(l,b)", 25, 0, 5);
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
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
  /// Matched events
  histo2D["mlb_corr_mlb_wrong_matched"] = new TH2F("mlb_corr_mlb_wrong_matched","Wrongly constructed M_{lb} vs. correctly constructed M_{lb}; M_{lb_{lep}}; M_{lb_{had}}", 80, 0, 800, 80, 0, 800);
  histo2D["dR_lep_b_corr_dR_lep_b_wrong_matched"] = new TH2F("dR_lep_b_corr_dR_lep_b_wrong_matched","Wrongly constructed dR(l,b) vs. correctly constructed dR(l,b); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["mlb_dR_lep_b_corr_matched"] = new TH2F("mlb_dR_lep_b_corr_matched","dR(l,b) vs. M_{lb}; M_{lb_{lep}}; #Delta R(l,b_{lep})", 80, 0, 800, 25, 0, 5);
  histo2D["mlb_dR_lep_b_wrong_matched"] = new TH2F("mlb_dR_lep_b_wrong_matched","dR(l,b) vs. M_{lb}, both wrongly matched; M_{lb_{had}}; #Delta R(l,b_{had})", 80, 0, 800, 25, 0, 5);
  
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
  
  nofMatchedEvents = 0;
  nofHadrMatchedEvents = 0;
  nofChi2First4 = 0;
  nofCorrectlyMatched_chi2 = 0;
  nofNotCorrectlyMatched_chi2 = 0;
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
  
  scaleFactor = 1.;
  widthSF = 1.;
  bJetId.clear();
  reco_hadWMass = -1.;
  reco_hadTopMass = -1.;
  reco_hadTopPt = -1.;
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
  kFitChi2 = 99.;
  toy = -1.;
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
  partonId.clear();
}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearMatching();
}

void FillGeneralPlots(int d)
{
  float M3 = (selectedJets[0] + selectedJets[1] + selectedJets[2]).M();
  float Ht = selectedJets[0].Pt() + selectedJets[1].Pt() + selectedJets[2].Pt() + selectedJets[3].Pt();
  
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
  float highestBDiscr = -999.;
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
  return BreitWigner(topMass, scale)/BreitWigner(topMass, 1);
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
  if ( alpha_WP < 0 ) ref = -ref;
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
  if ( alpha_UP < 0 ) ref = -ref;
  Double_t fitfunc = 1.;
  if ( ref > -alpha ) fitfunc = fitfunc * exp(-ref*ref/2.);
  else if (ref <= -alpha ) fitfunc = fitfunc * A * pow ( B - ref , -n_UP);
  return fitfunc*norm_UP;
}

Double_t logLikelihood(Double_t *x, Double_t *par) {
  return -TMath::Log( norm_comb*( f_CP*voigt(x, par) + f_WP*crysBall_WP(x) + f_UP*crysBall_UP(x) ) );
}

Double_t fakeLikelihood(Double_t *x, Double_t *par) {
  return -TMath::Log( voigt(x, par) );
}

Double_t widthToGammaTranslation(Double_t *x) {
  return 0.0247249 + 0.00689775 * x[0];
}
