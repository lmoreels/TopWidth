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
#include "TopTreeAnalysisBase/MCInformation/interface/TransferFunctions.h"


using namespace std;
using namespace TopTree;


bool test = false;
bool testHistos = false;
bool calculateTransferFunctions = false;
bool calculateAverageMass = false;
bool applyLeptonSF = true;
bool applyPU = true;
bool applyJER = true;
bool applyJEC = true;
bool applyBTagSF = true;
bool applyNloSF = false;

bool applyWidthSF = false;
float scaleWidth = 2.;

string systStr = "nominal";
string whichDate(string syst)
{
  if ( syst.find("nominal") == 0 ) return "160812";
  else if ( syst.find("JERup") == 0 ) return "160916";
  else if ( syst.find("JERdown") == 0 ) return "160930";
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
int nofChi2First4 = 0;
int nofCorrectlyMatched_chi2 = 0;
int nofNotCorrectlyMatched_chi2 = 0;
float Luminosity = 9999.;



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

// Average top mass
// TT gen match, TT reco match, TT reco noMatch, TT reco wrongPerm, TT reco wrongJets, TT reco, ST_t_top reco, ST_t_antitop reco, ST_tW_top reco, ST_tW_antitop reco, DYJets reco, WJets reco, data reco, all MC reco
float aveTopMass[] = {166.931, 169.290, 183.759, 162.863, 185.982, 197.667, 242.341, 235.626, 220.731, 222.694, 214.982, 198.189, 200.452, 197.965};

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;

vector < Dataset* > datasets;

ofstream txtMassGenMatched, txtMassRecoMatched, txtMassRecoNotMatched, txtMassRecoWrongPerm, txtMassRecoWrongJets, txtMassReco;

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
long GetNEvents(TTree* fChain, string var, bool isData);
void GetHLTFraction(double* fractions);
double BreitWigner(double topPT, double scale);
double eventWeightCalculator(double topPT, double scale);
//void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string pathPNG);
//void MSPCreator (string pathPNG);
//void TH2FPlotter (int nBinsX,float lowX, float highX, string sVarofinterestX );



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
double reco_minMlb, reco_ttbarMass, reco_dRLepB;
double min_Mlb, dRLepB;
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
vector<TLorentzVector> bJetsAfterChi2;

/// Matching
int pdgID_top = 6; //top quark

bool doMatching = true;
bool all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
bool all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
bool hadronictopJetsMatched_MCdef_ = false;
pair<unsigned int, unsigned int> MCPermutation[4] = {pair<unsigned int,unsigned int>(9999,9999)};
int topQuark = -9999, antiTopQuark = -9999;
int genmuon = -9999;
bool muonmatched = false;
bool muPlusFromTop = false, muMinusFromTop = false;
vector<unsigned int> partonId;

  
/// Meta
string strSyst = "";
vector<int> vJER, vJES, vPU;


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
  else if ( argv[1] == "mu" || argv[1] == "Mu" || argv[1] == "MU" || argv[1] == "m" ) channel = "mu";
  else if ( argv[1] == "el" || argv[1] == "El" || argv[1] == "EL" || argv[1] == "e" ) channel = "el";
  //else if ( argv[1] == "all" || argv[1] == "All" || argv[1] == "ALL" ) channel = "all";
  
  //string pathOutput = "test/";
  string pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (! test && ! calculateAverageMass)
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
  
  pathNtuples = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate+"/";
  cout << "Using Ntuples from " << ntupleDate << ". This corresponds to systematics: " << systStr << endl;
  if ( applyWidthSF && scaleWidth != 1 ) cout << "TTbar sample width will be scaled by a factor " << scaleWidth << endl;
  
  
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
  
  TransferFunctions* tf = new TransferFunctions(calculateTransferFunctions);
  
  if (! test && ! calculateAverageMass)
  {
    InitMSPlots();
    InitHisto1D();
    InitHisto2D();
  }
  
  vJER.clear(); vJES.clear(); vPU.clear();
  
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName, slumi;
  double timePerDataSet[datasets.size()] = {0};
  
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
    
    isData = false;
    if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
    {
      isData = true;
    }
    
    if (calculateAverageMass)
    {
      if (dataSetName.find("TT") == 0 )
      {
        txtMassGenMatched.open(("averageMass/mass_gen_matched_TT_"+dateString+".txt").c_str());
        txtMassRecoMatched.open(("averageMass/mass_reco_matched_TT_"+dateString+".txt").c_str());
        txtMassRecoNotMatched.open(("averageMass/mass_reco_notMatched_TT_"+dateString+".txt").c_str());
        txtMassRecoWrongPerm.open(("averageMass/mass_reco_wrongPerm_TT_"+dateString+".txt").c_str());
        txtMassRecoWrongJets.open(("averageMass/mass_reco_wrongJets_TT_"+dateString+".txt").c_str());
      }
      txtMassReco.open(("averageMass/mass_reco_"+dataSetName+"_"+dateString+".txt").c_str());
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
    if (test || testHistos) endEvent = 1001;
    for (int ievt = 0; ievt < endEvent; ievt++)
    {
      ClearObjects();
      
      
      if (ievt%1000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
      
      /// Load event
      tTree[(dataSetName).c_str()]->GetEntry(ievt);
      
      /// Scale factors
      if (! isData)
      {
        if (applyLeptonSF) { scaleFactor *= muonIdSF[0] * muonIsoSF[0] * (fracHLT[0]*muonTrigSFv2[0] + fracHLT[1]*muonTrigSFv3[0]);}
        if (applyPU) { scaleFactor *= puSF;}
        if (applyBTagSF) { scaleFactor *= btagSF;}
        //if (applyNloSF) { scaleFactor *= nloWeight;}  // additional SF due to number of events with neg weight!!
        //cout << "Scalefactor: " << setw(6) << scaleFactor << "  btag SF: " << setw(6) << btagSF << "  pu SF: " << setw(6) << puSF << "  muonId: " << setw(6) << muonIdSF[0] << "  muonIso: " << setw(6) << muonIsoSF[0] << "  muonTrig: " << setw(6) << fracHLT[0]*muonTrigSFv2[0] + fracHLT[1]*muonTrigSFv3[0] << endl;
      }
      
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
        
        if (doMatching)
        {
          TruthMatching(partons, selectedJets, MCPermutation);
          
          if (all4PartonsMatched && test)
          {
            for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
            {
              //cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Status: " << setw(2) << mc_status[partonId[MCPermutation[iMatch].second]] << "  pdgId: " << setw(3) << mc_pdgId[partonId[MCPermutation[iMatch].second]] << "  Mother: " << setw(4) << mc_mother[partonId[MCPermutation[iMatch].second]] << "  Granny: " << setw(4) << mc_granny[partonId[MCPermutation[iMatch].second]] << endl;
              cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Pt: " << setw(7) << left << mc_pt[partonId[MCPermutation[iMatch].second]] << "  Eta: " << mc_eta[partonId[MCPermutation[iMatch].second]] << "  Phi: " << mc_phi[partonId[MCPermutation[iMatch].second]] << endl;
              cout << "Event  " << right << setw(4) << ievt << ";  Matched jet    " << iMatch << "  Pt: " << setw(7) << left << jet_pt[MCPermutation[iMatch].first] << "  Eta: " << jet_eta[MCPermutation[iMatch].first] << "  Phi: " << jet_phi[MCPermutation[iMatch].first] << endl;
            }
          }
          
          
          ///////////////////
          ///  Transfer functions
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
            
            if (calculateTransferFunctions)
            {
              tf->fillJets(partonsMatched, jetsMatched);
              
              if (muonmatched) tf->fillMuon(mcParticles[genmuon], selectedLepton[0]);
              //if (electronmatched) tf->fillElectron(...)

            }  // end tf
            
            
            
            /// Plot variables for matched events
            float matchedWMass_reco = (jetsMatched[0] + jetsMatched[1]).M();
            float matchedTopMass_reco = (jetsMatched[0] + jetsMatched[1] + jetsMatched[2]).M();
            float matchedTopMass_gen = (partonsMatched[0] + partonsMatched[1] + partonsMatched[2]).M();
            
            if (calculateAverageMass) txtMassGenMatched << ievt << "  " << matchedWMass_reco << "  " << matchedTopMass_reco << "  " << matchedTopMass_gen << endl;
            
            if (! test && ! calculateAverageMass)
            {
              histo1D["W_mass_reco_matched"]->Fill(matchedWMass_reco);
              histo1D["top_mass_reco_matched"]->Fill(matchedTopMass_reco);
              histo1D["top_mass_gen_matched"]->Fill(matchedTopMass_gen);
              
              histo1D["mTop_div_aveMTop_TT_matched_reco"]->Fill(matchedTopMass_reco/aveTopMass[0]);
              
              if ( all4JetsMatched_MCdef_ )
              {
                histo1D["W_mass_reco_first4matched"]->Fill(matchedWMass_reco);
                histo1D["top_mass_reco_first4matched"]->Fill(matchedTopMass_reco);
                histo1D["top_mass_gen_first4matched"]->Fill(matchedTopMass_gen);
              }
              if (hasExactly4Jets)
              {
                histo1D["W_mass_reco_matched_4jets"]->Fill(matchedWMass_reco);
                histo1D["top_mass_reco_matched_4jets"]->Fill(matchedTopMass_reco);
                histo1D["top_mass_gen_matched_4jets"]->Fill(matchedTopMass_gen);
              }
              
              if (muonmatched)
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
                
              }  // end muonMatched
              
            }  // end fill plots

          }  // end all4PartonsMatched
        
        }  // end doMatching
        
        
        /////////////////////////////////////////
        ///  Scale factor ttbar sample width  ///
        /////////////////////////////////////////
        
        if (applyWidthSF)
        {
          if ( muon_charge > 0 ) massForWidth = (mcParticles[antiTopQuark]).M();
          else if ( muon_charge < 0 ) massForWidth = (mcParticles[topQuark]).M();
          
          widthSF = eventWeightCalculator(massForWidth, scaleWidth);
          
        }  // end applyWidthSF
        
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
      
      
      if ( labelsReco[0] < 4 && labelsReco[1] < 4 && labelsReco[2] < 4 )
      {
        nofChi2First4++;
      }
      
      if ( labelsReco[0] != -9999 && labelsReco[1] != -9999 && labelsReco[2] != -9999 )
      {
        reco_hadWMass = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).M();
        reco_hadTopMass = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).M();
        reco_hadTopPt = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).Pt();
        
        if (calculateAverageMass) txtMassReco << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
        
        //Fill histos
        if (! test && ! calculateAverageMass)
        {
          MSPlot["Chi2_value"]->Fill(smallestChi2, datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_W_mass"]->Fill(reco_hadWMass, datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_hadTop_mass"]->Fill(reco_hadTopMass, datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_hadTop_pT"]->Fill(reco_hadTopPt, datasets[d], true, Luminosity*scaleFactor);
          
          MSPlot["Chi2_mTop_div_aveMTopMatch"]->Fill(reco_hadTopMass/aveTopMass[0], datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_mTop_div_aveMTopTTChi2"]->Fill(reco_hadTopMass/aveTopMass[5], datasets[d], true, Luminosity*scaleFactor);
          MSPlot["Chi2_mTop_div_aveMTopAllChi2"]->Fill(reco_hadTopMass/aveTopMass[13], datasets[d], true, Luminosity*scaleFactor);
          
          if (isData) histo1D[("mTop_div_aveMTop_"+dataSetName).c_str()]->Fill(reco_hadTopMass/aveTopMass[12]);
          else histo1D[("mTop_div_aveMTop_"+dataSetName).c_str()]->Fill(reco_hadTopMass/aveTopMass[d+4]);
          
          if (hasExactly4Jets)
          {
            MSPlot["Chi2_hadTop_mass_4jets"]->Fill(reco_hadTopMass, datasets[d], true, Luminosity*scaleFactor);
          }
          
          if ( dataSetName.find("TT") == 0 )
          {
            histo1D["Chi2_W_mass_reco"]->Fill(reco_hadWMass);
            histo1D["Chi2_top_mass_reco"]->Fill(reco_hadTopMass);
            if (hasExactly4Jets)
            {
              histo1D["Chi2_top_mass_reco_4jets"]->Fill(reco_hadTopMass);
            }
          }
          
          /// Leptonic top mass
          for (int i = 0; i < selectedJets.size(); i++)
          {
            if ( i != labelsReco[0] && i != labelsReco[1] && i != labelsReco[2] 
                && jet_bdiscr[i] > CSVv2Medium )
              bJetsAfterChi2.push_back(selectedJets[i]);
          }
          
          if ( bJetsAfterChi2.size() > 0 )
          {
            double reco_Mlb_temp = 99999.;
            for (unsigned int i = 0; i < bJetsAfterChi2.size(); i++)
            {
              reco_Mlb_temp = (selectedLepton[0] + bJetsAfterChi2[i]).M();
              if ( reco_Mlb_temp < reco_minMlb )
              {
                reco_minMlb = reco_Mlb_temp;
                reco_dRLepB = ROOT::Math::VectorUtil::DeltaR( bJetsAfterChi2[i], selectedLepton[0]);
              }
            }
            
            reco_ttbarMass = reco_minMlb + reco_hadTopMass;
            
            MSPlot["Chi2_mlb"]->Fill(reco_minMlb, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["Chi2_ttbar_mass"]->Fill(reco_ttbarMass, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["Chi2_dR_lep_b"]->Fill(reco_dRLepB, datasets[d], true, Luminosity*scaleFactor);
            if ( dataSetName.find("TT") == 0 )
            {
              histo1D["dR_lep_b_unmatched_chi2"]->Fill(reco_dRLepB);
            }
            
            if (hasExactly4Jets)
            {
              MSPlot["Chi2_mlb_4jets"]->Fill(reco_minMlb, datasets[d], true, Luminosity*scaleFactor);
              MSPlot["Chi2_ttbar_mass_4jets"]->Fill(reco_ttbarMass, datasets[d], true, Luminosity*scaleFactor);
              MSPlot["Chi2_dR_lep_b_4jets"]->Fill(reco_dRLepB, datasets[d], true, Luminosity*scaleFactor);
            }
            
            if ( reco_minMlb < 200 )
            {
              MSPlot["Chi2_hadTop_mass_mlb_cut"]->Fill(reco_hadTopMass, datasets[d], true, Luminosity*scaleFactor);
              MSPlot["Chi2_ttbar_mass_mlb_cut"]->Fill(reco_ttbarMass, datasets[d], true, Luminosity*scaleFactor);
              MSPlot["Chi2_dR_lep_b_mlb_cut"]->Fill(reco_dRLepB, datasets[d], true, Luminosity*scaleFactor);
              if (hasExactly4Jets)
              {
                MSPlot["Chi2_hadTop_mass_4jets_mlb_cut"]->Fill(reco_hadTopMass, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["Chi2_ttbar_mass_4jets_mlb_cut"]->Fill(reco_ttbarMass, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["Chi2_dR_lep_b_4jets_mlb_cut"]->Fill(reco_dRLepB, datasets[d], true, Luminosity*scaleFactor);
              }
            }
            
          }  // end has bjet not in chi2
          
        }  // end fill plots
        
      }  // end labels
      
      
      ///////////////////////////////////
      ///  CHECK MATCHED COMBINATION  ///
      ///////////////////////////////////
      
      if ( dataSetName.find("TT") == 0 && all4PartonsMatched )
      {
        if ( ( (labelsReco[0] == MCPermutation[0].first && labelsReco[1] == MCPermutation[1].first) 
            || (labelsReco[0] == MCPermutation[1].first && labelsReco[1] == MCPermutation[0].first) )
            && labelsReco[2] == MCPermutation[2].first ) // correct jets, correct permutation
        {
          nofCorrectlyMatched_chi2++;
          if (calculateAverageMass) txtMassRecoMatched << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
          else
          {
            histo1D["dR_lep_b_reco_and_corr_match_chi2"]->Fill(reco_dRLepB);
            histo1D["mTop_div_aveMTop_TT_corr_match_chi2_reco"]->Fill(reco_hadTopMass/aveTopMass[1]);
          }
        }
        else
        {
          nofNotCorrectlyMatched_chi2++;
          if (calculateAverageMass) txtMassRecoNotMatched << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
          else
          {
            histo1D["dR_lep_b_reco_and_wrong_match_chi2"]->Fill(reco_dRLepB);
            histo1D["mTop_div_aveMTop_TT_wrong_match_chi2_reco"]->Fill(reco_hadTopMass/aveTopMass[2]);
          }
          
          if ( ( labelsReco[0] == MCPermutation[0].first || labelsReco[0] == MCPermutation[1].first || labelsReco[0] == MCPermutation[2].first ) && ( labelsReco[1] == MCPermutation[0].first || labelsReco[1] == MCPermutation[1].first || labelsReco[1] == MCPermutation[2].first ) && ( labelsReco[2] == MCPermutation[0].first || labelsReco[2] == MCPermutation[1].first || labelsReco[2] == MCPermutation[2].first ) ) // correct jets, wrong permutation
          {
            if (calculateAverageMass) txtMassRecoWrongPerm << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
            else histo1D["mTop_div_aveMTop_TT_wrong_perm_match_chi2_reco"]->Fill(reco_hadTopMass/aveTopMass[3]);
          }
          else // wrong jets
          {
            if (calculateAverageMass) txtMassRecoWrongJets << ievt << "  " << reco_hadWMass << "  " << reco_hadTopMass << endl;
            else histo1D["mTop_div_aveMTop_TT_wrong_jets_match_chi2_reco"]->Fill(reco_hadTopMass/aveTopMass[4]);
          }
        }
      }
      else if (! isData && dataSetName.find("TT") != 0  && ! calculateAverageMass)
      {
        histo1D["mTop_div_aveMTop_bkgd"]->Fill(reco_hadTopMass/aveTopMass[13]);
      }
      
      
      ////////////////////
      ///  Make plots  ///
      ////////////////////
      
      if (! test && ! calculateAverageMass)
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
        if ( dataSetName.find("TT") == 0 )
        {
          histo1D["dR_lep_b_unmatched_all"]->Fill(dRLepB);
        }
        
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
      }  // end fill plots
      
    }  // end loop events
    
    
    cout << endl;
    cout << "Number of chi2 made with the 4 most energetic jets: " << nofChi2First4 << " (" << 100*((float)nofChi2First4/(float)endEvent) << "%)" << endl;
    
    if ( dataSetName.find("TT") == 0 )
    {
      cout << "Number of matched events: " << nofMatchedEvents << endl;
      cout << "Correctly matched for chi2Free:     " << setw(8) << right << nofCorrectlyMatched_chi2 << endl;
      cout << "Not correctly matched for chi2Free: " << setw(8) << right << nofNotCorrectlyMatched_chi2 << endl;
      if ( nofCorrectlyMatched_chi2 != 0 || nofNotCorrectlyMatched_chi2 != 0 )
        cout << "   ===> This means that " << 100*(float)nofCorrectlyMatched_chi2 / (float)(nofCorrectlyMatched_chi2 + nofNotCorrectlyMatched_chi2) << "% is correctly matched." << endl;
      
      /// Transfer functions
      if (calculateTransferFunctions)
      {
        string tfFileName = "PlotsForTransferFunctions.root";
        TFile *foutTF = new TFile(tfFileName.c_str(), "RECREATE");
        foutTF->cd();

        tf->writeHistograms();

        foutTF->Close();

        tf->writeTable(tfFileName);

        delete foutTF;
      }
      
      if (calculateAverageMass)
      {
        txtMassGenMatched.close();
        txtMassRecoMatched.close();
        txtMassRecoNotMatched.close();
        txtMassRecoWrongPerm.close();
        txtMassRecoWrongJets.close();
      }
      
    }  // end TT
    
    
    if (calculateAverageMass) txtMassReco.close();
    
    tFileMap[dataSetName.c_str()]->Close();
    
    timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
    
  }  // end loop datasets
  
  
  cout << "Processing time per dataset: " << endl;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
  }
  
  
  //////////////////////////////////////////
  ///  Check Shape Changing Systematics  ///
  //////////////////////////////////////////
  
  int sumJER = 0, sumJES = 0, sumPU = 0;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    sumJER += vJER[d];
    sumJES += vJES[d];
    sumPU += vPU[d];
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
  cout << "  - Systematics confirmed to be " << strSyst << endl;
  
  
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
  MSPlot["Chi2_value"] = new MultiSamplePlot(datasets, "Chi2_value", 200, 0, 200, "#chi^{2} value");
  MSPlot["Chi2_W_mass"] = new MultiSamplePlot(datasets, "Chi2_W_mass", 40, 0, 800, "M_{W} [GeV]");
  MSPlot["Chi2_hadTop_mass"] = new MultiSamplePlot(datasets, "Chi2_hadTop_mass", 40, 0, 800, "M_{t} [GeV]");
  MSPlot["Chi2_hadTop_mass_4jets"] = new MultiSamplePlot(datasets, "Chi2_hadTop_mass_4jets", 40, 0, 800, "M_{t} [GeV]");
  MSPlot["Chi2_hadTop_pT"] = new MultiSamplePlot(datasets, "Chi2_hadTop_pT", 40, 0, 800, "p_{T} [GeV]");
  
  MSPlot["Chi2_mlb"] = new MultiSamplePlot(datasets, "Chi2_mlb", 40, 0, 800, "M_{lb} [GeV]");
  MSPlot["Chi2_mlb_4jets"] = new MultiSamplePlot(datasets, "Chi2_mlb_4jets", 40, 0, 800, "M_{lb} [GeV]");

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
  
  
  MSPlot["Chi2_mTop_div_aveMTopMatch"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass (from matched events)", 50, 0, 2, "M_{t}/<M_{t}>");
  MSPlot["Chi2_mTop_div_aveMTopTTChi2"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass (from chi2 of reco TT events)", 50, 0, 2, "M_{t}/<M_{t}>");
  MSPlot["Chi2_mTop_div_aveMTopAllChi2"] = new MultiSamplePlot(datasets, "Top quark mass divided by average top mass (from chi2 of reco events from all datasets)", 50, 0, 2, "M_{t}/<M_{t}>");
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
  
  /// Test dR
  histo1D["dR_lep_b_unmatched_all"] = new TH1F("dR_lep_b_unmatched_all","Minimal delta R between the lepton and a b jet (looking at all b jets); #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_lep_b_unmatched_chi2"] = new TH1F("dR_lep_b_unmatched_chi2","Minimal delta R between the lepton and a b jet (looking at b jets not in chi2); #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_lep_b_reco_and_corr_match_chi2"]  = new TH1F("dR_lep_b_reco_and_corr_match_chi2","Minimal delta R between the lepton and a b jet (where jet combination chi2 equals matching); #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_lep_b_reco_and_wrong_match_chi2"]  = new TH1F("dR_lep_b_reco_and_wrong_match_chi2","Minimal delta R between the lepton and a b jet (where jet combination chi2 differs from matching); #Delta R(l,b)", 25, 0, 5);
  
  /// m_t/<m_t>
  histo1D["mTop_div_aveMTop_TT_matched_reco"] = new TH1F("mTop_div_aveMTop_TT_matched_reco","Top mass divided by average top mass for matched TT sample (reco); M_{t}/<M_{t}>", 200, 0.4, 2.4);
  histo1D["mTop_div_aveMTop_TT_corr_match_chi2_reco"] = new TH1F("mTop_div_aveMTop_TT_corr_match_chi2_reco","Top mass divided by average top mass for matched TT sample (reco via chi2 and correct match); M_{t}/<M_{t}>", 200, 0.4, 2.4);
  histo1D["mTop_div_aveMTop_TT_wrong_match_chi2_reco"] = new TH1F("mTop_div_aveMTop_TT_wrong_match_chi2_reco","Top mass divided by average top mass for matched TT sample (reco via chi2 and wrong match); M_{t}/<M_{t}>", 200, 0.4, 2.4);
  histo1D["mTop_div_aveMTop_TT_wrong_perm_match_chi2_reco"] = new TH1F("mTop_div_aveMTop_TT_wrong_perm_chi2_reco","Top mass divided by average top mass for matched TT sample (reco via chi2 and wrong match: correct jets, wrong permutation); M_{t}/<M_{t}>", 200, 0.4, 2.4);
  histo1D["mTop_div_aveMTop_TT_wrong_jets_match_chi2_reco"] = new TH1F("mTop_div_aveMTop_TT_wrong_jets_chi2_reco","Top mass divided by average top mass for matched TT sample (reco via chi2 and wrong match: wrong jets); M_{t}/<M_{t}>", 200, 0.4, 2.4);
  histo1D["mTop_div_aveMTop_bkgd"] = new TH1F("mTop_div_aveMTop_bkgd","Top mass divided by average top mass for background samples; M_{t}/<M_{t}>", 50, 0, 2);
  histo1D["mTop_div_aveMTop_TT"] = new TH1F("mTop_div_aveMTop_TT","Top mass divided by average top mass for TT sample; M_{t}/<M_{t}>", 50, 0, 2);
  histo1D["mTop_div_aveMTop_ST_tW_top"] = new TH1F("mTop_div_aveMTop_ST_tW_top","Top mass divided by average top mass for ST tW top sample; M_{t}/<M_{t}>", 50, 0, 2);
  histo1D["mTop_div_aveMTop_ST_tW_antitop"] = new TH1F("mTop_div_aveMTop_ST_tW_antitop","Top mass divided by average top mass for ST tW antitop sample; M_{t}/<M_{t}>", 50, 0, 2);
  histo1D["mTop_div_aveMTop_ST_t_top"] = new TH1F("mTop_div_aveMTop_ST_t_top","Top mass divided by average top mass for ST t top sample; M_{t}/<M_{t}>", 50, 0, 2);
  histo1D["mTop_div_aveMTop_ST_t_antitop"] = new TH1F("mTop_div_aveMTop_ST_t_antitop","Top mass divided by average top mass for ST t antitop sample; M_{t}/<M_{t}>", 50, 0, 2);
  histo1D["mTop_div_aveMTop_DYJets"] = new TH1F("mTop_div_aveMTop_DYJets","Top mass divided by average top mass for DY+Jets sample; M_{t}/<M_{t}>", 50, 0, 2);
  histo1D["mTop_div_aveMTop_WJets"] = new TH1F("mTop_div_aveMTop_WJets","Top mass divided by average top mass for W+Jets sample; M_{t}/<M_{t}>", 50, 0, 2);
  histo1D["mTop_div_aveMTop_data"] = new TH1F("mTop_div_aveMTop_data","Top mass divided by average top mass for data sample; M_{t}/<M_{t}>", 50, 0, 2);
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
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

  if ( MCPermutation[0].first < 4 && MCPermutation[1].first < 4 && MCPermutation[2].first < 4 )
    hadronictopJetsMatched_MCdef_ = true;
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
  
  nofMatchedEvents = 0;
  nofChi2First4 = 0;
  nofCorrectlyMatched_chi2 = 0;
  nofNotCorrectlyMatched_chi2 = 0;
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
  muon_charge[1] = 0;
  muon_pt[1] = 0.;
  muon_phi[1] = 0.;
  muon_eta[1] = 0.;
  muon_E[1] = 0.;
  muon_M[1] = 0.;
  muon_d0[1] = 999.;
  muon_chargedHadronIso[1] = 999.;
  muon_neutralHadronIso[1] = 999.;
  muon_photonIso[1] = 999.;
  muon_puChargedHadronIso[1] = 999.;
  muon_relIso[1] = 999.;
  muon_pfIso[1] = 999.;
  nJets = -1;
  for (Int_t i = 0; i < 20; i++)
  {
    jet_charge[20] = 0;
    jet_pt[20] = 0.;
    jet_phi[20] = 0.;
    jet_eta[20] = 0.;
    jet_E[20] = 0.;
    jet_M[20] = 0.;
    jet_bdiscr[20] = -1.;
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
      mc_status[200] = -1;
      mc_pdgId[200] = 0;
      mc_mother[200] = 0;
      mc_granny[200] = 0;
      mc_pt[200] = 0.;
      mc_phi[200] = 0.;
      mc_eta[200] = 0.;
      mc_E[200] = 0.;
      mc_M[200] = 0.;
    }
  }
  nloWeight = 1.;
  puSF = 1.;
  btagSF = 1.;
  muonIdSF[1] = 1.;
  muonIsoSF[1] = 1.;
  muonTrigSFv2[1] = 1.;
  muonTrigSFv3[1] = 1.;
  
  scaleFactor = 1.;
  widthSF = 1.;
  bJetId.clear();
  reco_hadWMass = -1.;
  reco_hadTopMass = -1.;
  reco_hadTopPt = -1.;
  reco_minMlb = 9999.;
  reco_ttbarMass = -1.;
  reco_dRLepB = -1.;
  min_Mlb = 9999.;
  dRLepB = -1.;
  massForWidth = 0.01;
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
  bJetsAfterChi2.clear();
}

void ClearMatching()
{
  partons.clear();
  partonsMatched.clear();
  jetsMatched.clear();
  
  doMatching = true;
  all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
  all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
  hadronictopJetsMatched_MCdef_ = false;
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

