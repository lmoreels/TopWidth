#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"
#include <TFile.h>

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


bool calculateTransferFunctions = false;
bool applyLeptonSF = true;
bool applyPU = true;
bool applyJER = true;
bool applyJEC = true;
bool applyBTagSF = true;
bool applyNloSF = false;

string ntupleDate = "160729";
int verbose = 2;

string pathNtuples = "";
bool isData = false;

int nofMatchedEvents = 0;
float Luminosity = 9999.;

///  Working points for b tagging  // Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
float CSVv2Loose =  0.460;
float CSVv2Medium = 0.800;
float CSVv2Tight = 0.935;

// Temporarily, until calculated from TTbar sample
float chi2WMass = 80.385;
float sigmaChi2WMass = 15;
float chi2TopMass = 172.5;
float sigmaChi2TopMass = 40;

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;

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

string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();
void InitTree(TTree* tree, bool isData);
void InitMSPlots();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearObjects();
void GetHLTFraction(double* fractions);
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
Bool_t          cutFlow[10];
Bool_t          hasExactly4Jets;
Bool_t          hasJetLeptonCleaning;
Int_t           appliedJER;
Int_t           appliedJES;
Int_t           appliedPU;
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

// List of branches
TBranch        *b_run_num;   //!
TBranch        *b_evt_num;   //!
TBranch        *b_lumi_num;   //!
TBranch        *b_nvtx;   //!
TBranch        *b_npu;   //!
TBranch        *b_rho;   //!
TBranch        *b_isTrigged;   //!
TBranch        *b_cutFlow;   //!
TBranch        *b_hasExactly4Jets;   //!
TBranch        *b_hasJetLeptonCleaning;   //!
TBranch        *b_appliedJER;   //!
TBranch        *b_appliedJES;   //!
TBranch        *b_appliedPU;   //!
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

float scaleFactor, normFactor;
vector<unsigned int> bJetId;
double recoWMass, recoTopMass, recoTopPt;

/// Define TLVs
TLorentzVector muon, jet, mcpart;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
vector<TLorentzVector> mcParticles;
vector<TLorentzVector> partons;
vector<TLorentzVector> partonsMatched;
vector<TLorentzVector> jetsMatched;

/// Matching
int pdgID_top = 6; //top quark

bool all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
bool all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
bool hadronictopJetsMatched_MCdef_ = false;
pair<unsigned int, unsigned int> leptonicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
pair<unsigned int, unsigned int> hadronicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
pair<unsigned int, unsigned int> hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
pair<unsigned int, unsigned int> hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
pair<unsigned int, unsigned int> MCPermutation[4] = {pair<unsigned int,unsigned int>(9999,9999)};
vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
int topQuark = -9999, antiTopQuark = -9999;
int genmuon = -9999;
bool muonmatched = false;
bool muPlusFromTop = false, muMinusFromTop = false;
vector<unsigned int> partonId;


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
  // Add channel to output path
  pathOutput += channel+"/";
  mkdir(pathOutput.c_str(),0777);
  // Give timestamp to output path
  pathOutput += dateString+"/";
  mkdir(pathOutput.c_str(),0777);
  
  
  pathNtuples = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate+"/";
  cout << "Using Ntuples from " << ntupleDate << endl;
  
  
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
      datasets[d]->SetColor(kCyan);
      //datasets[d]->SetColor(kMagenta);
    }
    if ( dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") == 0 )
    {
      datasets[d]->SetColor(kBlue-2);
      if ( dataSetName.find("tW") == 0 )
        datasets[d]->SetTitle("ST tW");
      else
        datasets[d]->SetTitle("ST t");
    }
  }
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  TransferFunctions* tf = new TransferFunctions(calculateTransferFunctions);
  
  
  InitMSPlots();
  
  
  
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
    
    
    string ntupleFileName = "Ntuples_"+dataSetName+".root";
    tFileMap[dataSetName.c_str()] = new TFile((pathNtuples+ntupleFileName).c_str(),"READ"); //create TFile for each dataset
    
    string tTreeName = "tree";
    //string tStatsTreeName = "stats";
    tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
    //tStatsTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tStatsTreeName.c_str());
    
    nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
    cout << "                     nEntries: " << nEntries << endl;
    
    
    // Set branch addresses and branch pointers
    InitTree(tTree[dataSetName.c_str()], isData);
    
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    
    for (int ievt = 0; ievt < nEntries; ievt++)
    //for (int ievt = 0; ievt < 1000; ievt++)
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
          bJetId.push_back(iJet);
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
          
          
          if ( mc_pdgId[i] == pdgID_top )
            topQuark = i;
          else if( mc_pdgId[i] == -pdgID_top )
            antiTopQuark = i;
          
          if ( mc_status[i] == 23 && mc_pdgId[i] == 13 && mc_mother[i] == -24 && mc_granny[i] == -pdgID_top )		// mu-, W-, tbar
          {
            muMinusFromTop = true;
            genmuon = i;
          }
          if ( mc_status[i] == 23 && mc_pdgId[i] == -13 && mc_mother[i] == 24 && mc_granny[i] == pdgID_top )		// mu+, W+, t
          {
            muPlusFromTop = true;
            genmuon = i;
          }
          
          if ( abs(mc_pdgId[i]) < 6 || abs(mc_pdgId[i]) == 21 )  //light/b quarks, 6 should stay hardcoded, OR gluon
          {
            partons.push_back(mcParticles[i]);
            partonId.push_back(i);
          }
          
        }  // end loop mcParticles

        if (verbose > 3)
        {
          cout << "Size mcParticles:   " << mcParticles.size() << endl;
          cout << "Size partons:       " << partons.size() << endl;
          cout << "Size selectedJets:  " << selectedJets.size() << endl;
        }
        
        JetPartonMatching matching = JetPartonMatching(partons, selectedJets, 2, true, true, 0.3);  // partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
        
        if (matching.getNumberOfAvailableCombinations() != 1)
          cerr << "matching.getNumberOfAvailableCombinations() = " << matching.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
        
        
        /// Fill match in JetPartonPair; // vector< pair<unsigned int, unsigned int> > 
                                         // First one is jet number, second one is mcParticle number
        
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
              
              if (hadronicWJet1_.first == 9999)
              {
                hadronicWJet1_ = JetPartonPair[i];
                MCPermutation[0] = JetPartonPair[i];
              }
              else if (hadronicWJet2_.first == 9999)
              {
                hadronicWJet2_ = JetPartonPair[i];
                MCPermutation[1] = JetPartonPair[i];
              }
              else
              {
                cerr << "Found a third jet coming from a W boson which comes from a top quark..." << endl;
                cerr << " -- muMinusFromMtop: " << muMinusFromTop << " muPlusFromMtop: " << muPlusFromTop << endl;
                cerr << " -- pdgId: " << mc_pdgId[partonId[j]] << " mother: " << mc_mother[partonId[j]] << " granny: " << mc_granny[partonId[j]] << " Pt: " << mc_pt[partonId[j]] << endl;
                cerr << " -- ievt: " << ievt << endl;
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
              
              hadronicBJet_ = JetPartonPair[i];
              MCPermutation[2] = JetPartonPair[i];
            }
            else if ( ( muPlusFromTop && mc_mother[partonId[j]] == pdgID_top )
              || ( muMinusFromTop && mc_mother[partonId[j]] == -pdgID_top ) )  // if mu+ (top decay leptonic) and mother is top ---> leptonic b
            {
              if (verbose > 3)
                cout << "b jet:     " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;
              
              leptonicBJet_ = JetPartonPair[i];
              MCPermutation[3] = JetPartonPair[i];
            }
          }
        }  /// End loop over Jet Parton Pairs
        
        
        if ( hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999 && hadronicBJet_.first != 9999 && leptonicBJet_.first != 9999 )
        {
          all4PartonsMatched = true;
          nofMatchedEvents++;
          if (hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4 && leptonicBJet_.first < 4)
            all4JetsMatched_MCdef_ = true;
        }
        else if (verbose > 3) cout << "Size JetPartonPair: " << JetPartonPair.size() << ". Not all partons matched!" << endl;
        
        if ( hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4 )
          hadronictopJetsMatched_MCdef_ = true;
        if ( genmuon != -9999 && ROOT::Math::VectorUtil::DeltaR(mcParticles[genmuon], selectedLepton[0]) < 0.1 )
          muonmatched = true;
        
        
        
        ///////////////////
        ///  Transfer functions
        ///////////////////
        
        if (all4PartonsMatched && calculateTransferFunctions)
        {
          
          for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
          {
            /// MCPermutation[i].first  = jet number
            /// MCPermutation[i].second = parton number
            /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
            
            partonsMatched.push_back(partons[MCPermutation[iMatch].second]);
            jetsMatched.push_back(selectedJets[MCPermutation[iMatch].first]);
          }
          
          tf->fillJets(partonsMatched, jetsMatched);
          
          if (muonmatched) tf->fillMuon(mcParticles[genmuon], selectedLepton[0]);
          //if (electronmatched) tf->fillElectron(...)
          
        }  // end tf
        
        
      }  // end if TT
      
      
      
      ///////////////
      ///  CHI 2  ///
      ///////////////
      
      int labelsReco[3] = {-9999, -9999, -9999};		// 0,1 = light jets, 2 = hadronic b-jet.
      float recoWMass, recoTopMass, WTerm, topTerm, chi2;
      float smallestChi2 = 9999.;
      
      for (int ijet = 0; ijet < 4; ijet++)
      {
        for (int jjet = ijet+1; jjet < 4; jjet++)
        {
          for (int kjet = 0; kjet < 4; kjet++)
          {
            if ( ijet != kjet && jjet != kjet )
            {
              recoWMass = (selectedJets[ijet] + selectedJets[jjet]).M();
              recoTopMass = (selectedJets[ijet] + selectedJets[jjet] + selectedJets[kjet]).M();
              
              WTerm = pow( (recoWMass - chi2WMass)/sigmaChi2WMass, 2);
              topTerm = pow ( (recoTopMass - chi2TopMass)/sigmaChi2TopMass, 2);
              
              chi2 = WTerm + topTerm;
              
              if ( chi2 < smallestChi2)
              {
                smallestChi2 = chi2;
                labelsReco[0] = ijet;
                labelsReco[1] = jjet;
                labelsReco[2] = kjet;
              }
            }
          }
        }
      }
      
      if (labelsReco[0] != -9999 && labelsReco[1] != -9999 && labelsReco[2] != -9999)
      {
        recoWMass = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).M();
        recoTopMass = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).M();
        recoTopPt = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).Pt();
      }
      
      
      /// Make plots
      float Ht = selectedJets[0].Pt() + selectedJets[1].Pt() + selectedJets[2].Pt() + selectedJets[3].Pt();
      
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
      
    }  // end loop events
    
    cout << endl;
    
    if ( dataSetName.find("TT") == 0 )
    {
      cout << "Number of matched events: " << nofMatchedEvents << endl;
      
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
    }  // end TT
    
    
    tFileMap[dataSetName.c_str()]->Close();
    
    timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
    
  }  // end loop datasets
  
  
  cout << "Processing time per dataset: " << endl;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
  }
  
  
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  string rootFileName = "NtuplePlots_nominal.root";
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
   tree->SetBranchAddress("cutFlow", cutFlow, &b_cutFlow);
   tree->SetBranchAddress("hasExactly4Jets", &hasExactly4Jets, &b_hasExactly4Jets);
   tree->SetBranchAddress("hasJetLeptonCleaning", &hasJetLeptonCleaning, &b_hasJetLeptonCleaning);
   tree->SetBranchAddress("appliedJER", &appliedJER, &b_appliedJER);
   tree->SetBranchAddress("appliedJES", &appliedJES, &b_appliedJES);
   tree->SetBranchAddress("appliedPU", &appliedPU, &b_appliedPU);
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

  MSPlot["nJets"] = new MultiSamplePlot(datasets, "nJets", 13, -0.5, 12.5, "# jets");
  MSPlot["nBJets"] = new MultiSamplePlot(datasets, "nBJets", 9, -0.5, 8.5, "# b jets");
  MSPlot["CSVv2Discr_allJets"] = new MultiSamplePlot(datasets, "CSVv2Discr_allJets", 48, 0.0, 1.2, "CSVv2 discriminant value");
  MSPlot["CSVv2Discr_leadingJet"] = new MultiSamplePlot(datasets, "CSVv2Discr_leadingJet", 48, 0.0, 1.2, "CSVv2 discriminant value of leading jet");
  MSPlot["CSVv2Discr_jet2"] = new MultiSamplePlot(datasets, "CSVv2Discr_jet2", 48, 0.0, 1.2, "CSVv2 discriminant value of jet2");
  MSPlot["CSVv2Discr_jet3"] = new MultiSamplePlot(datasets, "CSVv2Discr_jet3", 48, 0.0, 1.2, "CSVv2 discriminant value of jet3");
  MSPlot["CSVv2Discr_jet4"] = new MultiSamplePlot(datasets, "CSVv2Discr_jet4", 48, 0.0, 1.2, "CSVv2 discriminant value of jet4");
  MSPlot["CSVv2Discr_highest"] = new MultiSamplePlot(datasets, "CSVv2Discr_highest", 48, 0.0, 1.2, "Highest CSVv2 discriminant value");
  MSPlot["CSVv2Discr_jetNb"] = new MultiSamplePlot(datasets, "CSVv2Discr_jetNb", 48, 0.0, 1.2, "Jet number (in order of decreasing p_{T}) with highest CSVv2 discriminant value");
  
//   MSPlot["min_M_lb"] = new MultiSamplePlot(datasets, "min_M_lb", 40, 0, 400, "M_{lb} [GeV]");
//   MSPlot["dR_Lep_B"] = new MultiSamplePlot(datasets, "dR_Lep_B", 50, 0, 10, "#Delta R(l,b)");
  
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
  normFactor = 1.;
  bJetId.clear();
  recoWMass = -1.;
  recoTopMass = -1.;
  recoTopPt = -1.;
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
  partons.clear();
  partonsMatched.clear();
  jetsMatched.clear();
}

void ClearMatching()
{
  partons.clear();
  partonsMatched.clear();
  jetsMatched.clear();
  
  all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
  all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
  hadronictopJetsMatched_MCdef_ = false;
  leptonicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
  hadronicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
  hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
  hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
  for (int i = 0; i < 4; i++)
  {
    MCPermutation[i] = pair<unsigned int,unsigned int>(9999,9999);
  }
  JetPartonPair.clear();
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

void GetHLTFraction(double* fractions)
{
  TFile* fileHLT = new TFile((pathNtuples+"Ntuples_data.root").c_str(),"READ");
  TTree* fChain = (TTree*) fileHLT->Get("stats");
  
  Long64_t        nofEventsHLTv2;
  Long64_t        nofEventsHLTv3;
  TBranch        *b_nofEventsHLTv2;   //!
  TBranch        *b_nofEventsHLTv3;   //!
  
  fChain->SetBranchAddress("nofEventsHLTv2", &nofEventsHLTv2, &b_nofEventsHLTv2);
  fChain->SetBranchAddress("nofEventsHLTv3", &nofEventsHLTv3, &b_nofEventsHLTv3);
  
//  Long64_t        nofSelEventsHLTv2;
//  Long64_t        nofSelEventsHLTv3;
//  TBranch        *b_nofSelEventsHLTv2;   //!
//  TBranch        *b_nofSelEventsHLTv3;   //!
//  fChain->SetBranchAddress("nofSelEventsHLTv2", &nofSelEventsHLTv2, &b_nofSelEventsHLTv2);
//  fChain->SetBranchAddress("nofSelEventsHLTv3", &nofSelEventsHLTv3, &b_nofSelEventsHLTv3);
  
  fractions[0] = ((double)nofEventsHLTv2) / ((double)nofEventsHLTv2 + (double)nofEventsHLTv3);
  fractions[1] = ((double)nofEventsHLTv3) / ((double)nofEventsHLTv2 + (double)nofEventsHLTv3);
  
  fileHLT->Close();
}

