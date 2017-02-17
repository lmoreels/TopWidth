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


bool test = false;

bool applyWidthSF = false;
float scaleWidth = 1.;

string systStr = "nominal";
string whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos ) return "170216";
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
bool makePlots = true;

int nofHardSelected = 0;
int nofTopNoStatus22 = 0;
int nofAntiTopNoStatus22 = 0;
int nofTopGrannyProton = 0;
int nofAntiTopGrannyProton = 0;
int nofTopMotherGluonGrannyParton = 0;
int nofAntiTopMotherGluonGrannyParton = 0;
int nofNoTops = 0;
int nofTopST = 0;
int nofAntiTopST = 0;
double Luminosity = 9999.;



///  Working points for b tagging  // Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
float CSVv2Loose =  0.460;
float CSVv2Medium = 0.800;
float CSVv2Tight = 0.935;

/// Top width
double genTopWidth = 1.33; // from fit
double genTopMass = 172.5; // from fit


map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

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


string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();
void GetMetaData(TTree* tree, bool isData);
void InitTree(TTree* tree, bool isData);
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearVars();
void ClearObjects();
long GetNEvents(TTree* fChain, string var, bool isData);
void GetHLTFraction(double* fractions);
double BreitWigner(double topPT, double scale);
double eventWeightCalculator(double topPT, double scale);



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
double massForWidth;


/// Define TLVs
TLorentzVector muon, jet, mcpart;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
vector<TLorentzVector> mcParticles;

/// Matching
int pdgID_top = 6; //top quark
int topQuark = -9999, antiTopQuark = -9999;
int genmuon = -9999;
bool muPlusFromTop = false, muMinusFromTop = false;

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
  
  clock_t start = clock();
  
  string channel;
  if ( argc == 1) channel = "mu";
  //else if ( CharSearch(argv[1], "mu") || CharSearch(argv[1], "Mu") || CharSearch(argv[1], "MU") || CharSearch(argv[1], "m") ) channel = "mu";
  //else if ( CharSearch(argv[1], "el") || CharSearch(argv[1], "El") || CharSearch(argv[1], "EL") || CharSearch(argv[1], "e") ) channel = "el";
  //else if ( (argv[1]).find("all") != std::string::npos || (argv[1]).find("All") != std::string::npos || (argv[1]).find("ALL") != std::string::npos ) channel = "all";
  
  
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
  for (int d = 1; d < 2; d++)   //Only ttbar
  {
    clock_t startDataSet = clock();
    
    ClearMetaData();
    
    dataSetName = datasets[d]->Name();
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
    }
    
    isData = false;
    if ( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
    {
      isData = true;
    }
    
    
    string ntupleFileName = "Ntuples_"+dataSetName+".root";
    tFileMap[dataSetName.c_str()] = new TFile((pathNtuples+ntupleFileName).c_str(),"READ"); //create TFile for each dataset
    
    string tTreeName = "tree";
    string tStatsTreeName = "stats";
    
    /// Get meta data
    tStatsTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tStatsTreeName.c_str());
    GetMetaData(tStatsTree[dataSetName.c_str()], isData);
    
    tStatsTree[(dataSetName).c_str()]->GetEntry(0);
    
    
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
    int startEvent = 0;
    if (test) 
    {
      startEvent = 1503;
      endEvent = startEvent + 41;
    }
    for (int ievt = startEvent; ievt < endEvent; ievt++)
    {
      ClearObjects();
      
      //if ( ievt != 35 && ievt != 283 && ievt != 378 && ievt != 479 && ievt != 501 && ievt != 572 && ievt != 876 && ievt != 985 && ievt != 1071 && ievt != 1191 && ievt != 1230 && ievt != 1311 && ievt != 1704 && ievt != 1721 && ievt != 1734 && ievt != 1754 && ievt != 1814 && ievt != 2008 && ievt != 2067 && ievt != 2317 && ievt != 2481 && ievt != 2497 && ievt != 2695 && ievt != 2851 && ievt != 2911 && ievt != 2963 && ievt != 3003 && ievt != 3043 && ievt != 3294 && ievt != 3308 && ievt != 3425 && ievt != 3444 && ievt != 3508 && ievt != 3562 && ievt != 3732 && ievt != 3740 && ievt != 3851 && ievt != 4017 && ievt != 4073 && ievt != 4201 && ievt != 4373 && ievt != 4433 && ievt != 4783 && ievt != 4878 && ievt != 4879 && ievt != 4960 && ievt != 4966 && ievt != 5077 && ievt != 5345 && ievt != 5422 ) continue;
      
      if (ievt%10000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
      
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
        
      
      
      ////////////////////////////
      ///  LOAD GEN PARTICLES  ///
      ////////////////////////////
      
      if ( dataSetName.find("TT") == 0 )
      {
        for (int iMC = 0; iMC < nMCParticles; iMC++)
        {
          mcpart.Clear();
          mcpart.SetPtEtaPhiE(mc_pt[iMC], mc_eta[iMC], mc_phi[iMC], mc_E[iMC]);
          mcParticles.push_back(mcpart);
        }
        
        if (test) cout << "Event:  " << ievt << endl;
        for (unsigned int i = 0; i < mcParticles.size(); i++)
        {
          //if (test && verbose > 4)
          if ( test && abs(mc_pdgId[i]) == pdgID_top )
            cout << setw(3) << right << i << "  Status: " << setw(2) << mc_status[i] << "  pdgId: " << setw(3) << mc_pdgId[i] << "  Mother: " << setw(4) << mc_mother[i] << "  Granny: " << setw(4) << mc_granny[i] << "  Mass: " << setw(7) << left << mc_M[i] << "  Pt: " << setw(7) << left << mc_pt[i] << "  Eta: " << setw(7) << left << mc_eta[i] << "  iLC: " << mc_isLastCopy[i] << "  isHP: " << mc_isHardProcess[i] << endl;
          
          
          if ( mc_pdgId[i] == pdgID_top )  // isLastCopy() == status 62
          {
            if ( mc_status[i] == 22 ) topQuark = i;
            if ( topQuark == -9999 && mc_isLastCopy[i] )
            {
              topQuark = i;
              nofTopNoStatus22++;
            }
          }
          else if ( mc_pdgId[i] == -pdgID_top )
          {
            if ( mc_status[i] == 22 ) antiTopQuark = i;
            if ( antiTopQuark == -9999 && mc_isLastCopy[i] )
            {
              antiTopQuark = i;
              nofAntiTopNoStatus22++;
            }
          }
          
          if ( mc_pdgId[i] == pdgID_top)
          {
            if ( mc_granny[i] == 2212 )
              nofTopGrannyProton++;
            else if ( mc_mother[i] == 21 && abs(mc_granny[i]) < 6 )
              nofTopMotherGluonGrannyParton++;
          }
          else if( mc_pdgId[i] == -pdgID_top )
          {
            if ( mc_granny[i] == 2212 )
              nofAntiTopGrannyProton++;
            else if ( mc_mother[i] == 21 && abs(mc_granny[i]) < 6 )
              nofAntiTopMotherGluonGrannyParton++;
          }
          
          
          if ( (mc_status[i] > 1 && mc_status[i] <= 20) || mc_status[i] >= 30 ) continue;  /// Final state particle or particle from hardest process
          
          
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
            
          }
          
        }  // end loop mcParticles
        
        
        if ( muMinusFromTop && muPlusFromTop )
        {
          if (test) cout << "Both tops decay leptonically... Event " << ievt << " will not be matched." << endl;
        }
        
        if ( topQuark == -9999 && antiTopQuark == -9999 )
          nofNoTops++;
        else if ( topQuark == -9999 )
          nofAntiTopST++;
        if ( antiTopQuark == -9999 )
          nofTopST++;
        
        
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
        
        
      }  // end if TT
      
      
      
      
    }  // end loop events
    
    
    cout << endl;  /// Stronger selection in this analyser compared to Ntuples ==> endEvent --> nofHardSelected
    cout << "Number of events with exactly 4 jets with pT > 30 GeV: " << nofHardSelected << " (" << 100*((float)nofHardSelected/(float)endEvent) << "%)" << endl << endl;
    
    cout << "Nb. of events without tops   " << nofNoTops << endl;
    cout << "Nb. of events with single top     " << nofTopST << endl;
    cout << "Nb. of events with single antitop " << nofAntiTopST << endl << endl;
    cout << "Nb. events without status 22 top      " << nofTopNoStatus22 << endl;
    cout << "Nb. events without status 22 antitop  " << nofAntiTopNoStatus22 << endl << endl;
    cout << "Nb. events with top granny proton        " << nofTopGrannyProton << endl;
    cout << "Nb. events with antitop granny proton    " << nofAntiTopGrannyProton << endl;
    cout << "Nb. events with top from glu from part     " << nofTopMotherGluonGrannyParton << endl;
    cout << "Nb. events with antitop from glu from part " << nofAntiTopMotherGluonGrannyParton << endl << endl;
    
    
    tFileMap[dataSetName.c_str()]->Close();
    
    timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
    
  }  // end loop datasets
  
  
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
    tree->SetBranchAddress("mc_isLastCopy", mc_isLastCopy, &b_mc_isLastCopy);
    tree->SetBranchAddress("mc_isPromptFinalState", mc_isPromptFinalState, &b_mc_isPromptFinalState);
    tree->SetBranchAddress("mc_isHardProcess", mc_isHardProcess, &b_mc_isHardProcess);
    tree->SetBranchAddress("mc_fromHardProcessFinalState", mc_fromHardProcessFinalState, &b_mc_fromHardProcessFinalState);
    tree->SetBranchAddress("nloWeight", &nloWeight, &b_nloWeight);
    tree->SetBranchAddress("puSF", &puSF, &b_puSF);
    tree->SetBranchAddress("btagSF", &btagSF, &b_btagSF);
    tree->SetBranchAddress("muonIdSF", muonIdSF, &b_muonIdSF);
    tree->SetBranchAddress("muonIsoSF", muonIsoSF, &b_muonIsoSF);
    tree->SetBranchAddress("muonTrigSFv2", muonTrigSFv2, &b_muonTrigSFv2);
    tree->SetBranchAddress("muonTrigSFv3", muonTrigSFv3, &b_muonTrigSFv3);
  }
  
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
  selectedLepton.clear();
  selectedJets.clear();
  selectedBJets.clear();
  mcParticles.clear();
}

void ClearMatching()
{
  topQuark = -9999;
  antiTopQuark = -9999;
  genmuon = -9999;
  muPlusFromTop = false;
  muMinusFromTop = false;
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
