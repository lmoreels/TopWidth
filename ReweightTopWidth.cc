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
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"  // needed for ROOT::Math::VectorUtil::DeltaR()
//#include "../macros/Style.C"



using namespace std;
using namespace TopTree;


bool test = false;
bool makePlots = true;
string systStr = "nominal";

string ntupleDate = "160812"; // nominal
//string ntupleDate = "160916"; // JERup
int verbose = 2;

string pathNtuples = "";

int nofHardSelected = 0;
int nofSkippedEvents = 0;
int nofMatchedEvents = 0;
double nofEvents_hadr[3] = {0.};
double nofEvents_lept[3] = {0.};
double nofEvents_test[3] = {0.};


///  Working points for b tagging  // Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
double CSVv2Loose =  0.460;
double CSVv2Medium = 0.800;
double CSVv2Tight = 0.935;

// Temporarily, until calculated from TTbar sample
double chi2WMass = 80.385;
double sigmaChi2WMass = 15;
double chi2TopMass = 172.5; //180.0; //from mtop mass plot: 167.0
double sigmaChi2TopMass = 40;

// Average top mass
// TT gen match, TT reco match, TT reco wrongMatch WP/UP, TT reco noMatch, TT reco wrongPerm, TT reco wrongPerm W Ok, TT reco wrongPerm W Not Ok, TT reco, ST_t_top reco, ST_t_antitop reco, ST_tW_top reco, ST_tW_antitop reco, DYJets reco, WJets reco, data reco, all MC reco, all samples reco (data+MC) 
double aveTopMass[] = {168.719, 167.253, 192.093, 189.672, 196.716, 199.756, 165.839, 180.817, 249.629, 249.039, 227.992, 224.213, 221.995, 213.278, 184.884, 181.158, 181.191};


/// Top width
double genTopWidth = 1.31; // gen  //1.363; // from fit  //4.015; // from reco fit
double genTopMass = 172.5; // gen  //172.3; // from fit  //168.6; // from reco fit


// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

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
Double_t BreitWigner2Func(Double_t *x, Double_t *scale);
double BreitWigner(double topPT, double scale);
double BreitWigner2(double topMass, double scale);
double BreitWignerNonRel(double topMass, double scale);
double eventWeightCalculator(double topPT, double scale);
void DrawBW(TF1* function, string name, std::vector<double> scales, double min, double max, bool writeToFile);
void FitBW(TH1F* histoIn, string name, double scale);
void DrawBWonHisto(TH1F* histo, TF1* function, string name, double scale, double norm, double min, double max, bool writeToFile);



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


double scaleFactor;
vector<unsigned int> bJetId;
vector<int> selectedJetsCharge;
vector<double> selectedJetsBDiscr;
double reco_hadWMass, reco_hadTopMass, reco_hadTopPt;
double reco_minMlb, reco_ttbarMass, reco_dRLepB;
double min_Mlb, dRLepB;

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
pair<unsigned int, unsigned int> MCPermutation[4] = {pair<unsigned int,unsigned int>(9999,9999)};
int topQuark = -9999, antiTopQuark = -9999;
int genmuon = -9999;
int genneutrino = -9999;
bool muonmatched = false;
bool muPlusFromTop = false, muMinusFromTop = false;
int bjjDecay[] = {-9999, -9999, -9999}, blvDecay[] = {-9999, -9999, -9999};
vector<unsigned int> partonId;


/// Scale width
Double_t eventSF_gen = 1.;
Double_t eventSF_reco = 1.;
Double_t evWeight_hadr = 1.;
Double_t evWeight_lept = 1.;
Double_t evWeight_bjj = 1.;
Double_t evWeight_blv = 1.;
Double_t m_top = -1.;
Double_t m_antitop = -1.;
Double_t m_bjj = -1.;
Double_t m_blv = -1.;
Double_t m_hadr = -1;
Double_t m_lept = -1;
Double_t scaling[] = {1., 2., 0.5};
string scalingString[] = {"orig", "s2", "s0p5"};
int nScalings = sizeof(scaling)/sizeof(scaling[0]);





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
  if ( argc == 1 ) channel = "mu";
  //else if ( CharSearch(argv[1], "mu") || CharSearch(argv[1], "Mu") || CharSearch(argv[1], "MU") || CharSearch(argv[1], "m") ) channel = "mu";
  //else if ( CharSearch(argv[1], "el") || CharSearch(argv[1], "El") || CharSearch(argv[1], "EL") || CharSearch(argv[1], "e") ) channel = "el";
  //else if ( (argv[1]).find("all") != std::string::npos || (argv[1]).find("All") != std::string::npos || (argv[1]).find("ALL") != std::string::npos ) channel = "all";
  
  if (test) makePlots = false;
  
  string pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (makePlots)
  {
    pathOutput += "reweighting/";
    mkdir(pathOutput.c_str(),0777);
    // Add channel to output path
    pathOutput += channel+"/";
    mkdir(pathOutput.c_str(),0777);
    // Give timestamp to output path
    pathOutput += dateString+"/";
    mkdir(pathOutput.c_str(),0777);
  }
  
  
  pathNtuples = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate+"/";
  cout << "Using Ntuples from " << ntupleDate << ". This corresponds to systematics: " << systStr << endl;
  
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  if (makePlots)
  {
    InitHisto1D();
    InitHisto2D();
  }
  ClearMetaData();
  
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  int nEntries;
  
  /// Load original Ntuple
  cout << "  - Loading original Ntuples..." << endl;
  TFile* origNtuple = new TFile((pathNtuples+"Ntuples_TT.root").c_str(), "READ");

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
  if (test) endEvent = nEntries;
  for (int ievt = 0; ievt < endEvent; ievt++)
  {
    ClearObjects();
    
    
    if (ievt%10000 == 0)
      std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
    
    /// Load event
    origTree->GetEntry(ievt);
    
    
    /// Fill objects
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
    
    if ( selectedJets.size() != 4 ) continue;
    
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
    
    
    
    /////////////////////////////
    ///  JET PARTON MATCHING  ///
    /////////////////////////////

    for (int iMC = 0; iMC < nMCParticles; iMC++)
    {
      mcpart.Clear();
      mcpart.SetPtEtaPhiE(mc_pt[iMC], mc_eta[iMC], mc_phi[iMC], mc_E[iMC]);
      mcParticles.push_back(mcpart);
    }
    
    
    for (unsigned int i = 0; i < mcParticles.size(); i++)
    {
      if (verbose > 4)
        cout << setw(3) << right << i << "  Status: " << setw(2) << mc_status[i] << "  pdgId: " << setw(3) << mc_pdgId[i] << "  Mother: " << setw(4) << mc_mother[i] << "  Granny: " << setw(4) << mc_granny[i] << "  Pt: " << setw(7) << left << mc_pt[i] << "  Eta: " << mc_eta[i] << "  Mass: " << mc_M[i] << endl;
      
      
      /// Find tops
      if ( mc_pdgId[i] == pdgID_top )  // isLastCopy() == status 62
      {
        if ( mc_status[i] == 22 ) topQuark = i;
        else if ( topQuark == -9999 && mc_status[i] == 62 ) topQuark = i;
      }
      else if ( mc_pdgId[i] == -pdgID_top )
      {
        if ( mc_status[i] == 22 ) antiTopQuark = i;
        else if ( antiTopQuark == -9999 && mc_status[i] == 62 ) antiTopQuark = i;
      }
      
      
      /// Status restriction: Final state particle or particle from hardest process
      if ( (mc_status[i] > 1 && mc_status[i] <= 20) || mc_status[i] >= 30 ) continue;
      
      if ( mc_pdgId[i] == 13 && mc_mother[i] == -24 && mc_granny[i] == -pdgID_top )    // mu-, W-, tbar
      {
        muMinusFromTop = true;
        if ( mc_status[i] == 23 ) genmuon = i;
        else if ( mc_status[i] != 23 && genmuon == -9999 ) genmuon = i;
      }
      if ( mc_pdgId[i] == -13 && mc_mother[i] == 24 && mc_granny[i] == pdgID_top )    // mu+, W+, t
      {
        muPlusFromTop = true;
        if ( mc_status[i] == 23 ) genmuon = i;
        else if ( mc_status[i] != 23 && genmuon == -9999 ) genmuon = i;
      }
      
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
    
    if ( topQuark == -9999 || antiTopQuark == -9999 )
    {
      if ( test && verbose > 2 && topQuark == -9999 && antiTopQuark == -9999 ) cout << "WARNING: Both top and antitop quark not found";
      else if ( test && verbose > 2 && topQuark == -9999 ) cout << "WARNING: Top quark not found";
      else if ( test && verbose > 2 && antiTopQuark == -9999 ) cout << "WARNING: Antitop quark not found";
      if ( test && verbose > 2 ) cout << "... Skipping event " << ievt << "... " << endl;
      nofSkippedEvents++;
      continue;
    }
    
    if ( muMinusFromTop && muPlusFromTop )
    {
      if (test && verbose > 2 ) cout << "Both tops decay leptonically... Something went wrong... Event " << ievt << endl;
      nofSkippedEvents++;
      continue;
    }
    
    
    TruthMatching(partons, selectedJets, MCPermutation);
    
    if (! all4PartonsMatched) continue;
    
    if ( test && verbose > 3 )
    {
      for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
      {
        //cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Status: " << setw(2) << mc_status[partonId[MCPermutation[iMatch].second]] << "  pdgId: " << setw(3) << mc_pdgId[partonId[MCPermutation[iMatch].second]] << "  Mother: " << setw(4) << mc_mother[partonId[MCPermutation[iMatch].second]] << "  Granny: " << setw(4) << mc_granny[partonId[MCPermutation[iMatch].second]] << endl;
        cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Pt: " << setw(7) << left << mc_pt[partonId[MCPermutation[iMatch].second]] << "  Eta: " << mc_eta[partonId[MCPermutation[iMatch].second]] << "  Phi: " << mc_phi[partonId[MCPermutation[iMatch].second]] << endl;
        cout << "Event  " << right << setw(4) << ievt << ";  Matched jet    " << iMatch << "  Pt: " << setw(7) << left << jet_pt[MCPermutation[iMatch].first] << "  Eta: " << jet_eta[MCPermutation[iMatch].first] << "  Phi: " << jet_phi[MCPermutation[iMatch].first] << endl;
      }
    }
    
    
    for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
    {
      /// MCPermutation[i].first  = jet number
      /// MCPermutation[i].second = parton number
      /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
      
      partonsMatched.push_back(partons[MCPermutation[iMatch].second]);
      jetsMatched.push_back(selectedJets[MCPermutation[iMatch].first]);
    }
    
    
    
    ///////////////////
    ///  Calculate event weight
    ///////////////////
    
    
    /// Calculate masses before matching
    for (unsigned int i = 0; i < partons.size(); i++)  /// partons[j] = mcParticles[partonId[j]]
    {
      if ( muMinusFromTop )
      {
        // hadronic decay  --> top
        if ( mc_pdgId[partonId[i]] == 5 && mc_mother[partonId[i]] == pdgID_top ) bjjDecay[0] = partonId[i];
        else if ( abs(mc_pdgId[partonId[i]]) < 5 && mc_mother[partonId[i]] == 24 && mc_granny[partonId[i]] == pdgID_top )
        {
          //cout << "In hadronic loop: " << setw(3) << right << partonId[i] << "  pdgId: " << setw(3) << mc_pdgId[partonId[i]] << "  Mother: " << setw(4) << mc_mother[partonId[i]] << "  Granny: " << setw(4) << mc_granny[partonId[i]] << endl;
          
          if ( bjjDecay[1] == -9999 ) bjjDecay[1] = partonId[i];
          else if ( bjjDecay[2] == -9999 ) bjjDecay[2] = partonId[i];
          else { cerr << "ERROR: Too many partons for hadronic top decay!  Event: " << ievt << endl; continue; }
        }
        
        if ( mc_pdgId[partonId[i]] == -5 && mc_mother[partonId[i]] == -pdgID_top ) blvDecay[0] = partonId[i];  // leptonic b
      }
      else if ( muPlusFromTop )
      {
        // hadronic decay  --> antitop
        if ( mc_pdgId[partonId[i]] == -5 && mc_mother[partonId[i]] == -pdgID_top ) bjjDecay[0] = partonId[i];
        else if ( abs(mc_pdgId[partonId[i]]) < 5 && mc_mother[partonId[i]] == -24 && mc_granny[partonId[i]] == -pdgID_top )
        {
          //cout << "In hadronic loop: " << setw(3) << right << partonId[i] << "  pdgId: " << setw(3) << mc_pdgId[partonId[i]] << "  Mother: " << setw(4) << mc_mother[partonId[i]] << "  Granny: " << setw(4) << mc_granny[partonId[i]] << endl;
          
          if ( bjjDecay[1] == -9999 ) bjjDecay[1] = partonId[i];
          else if ( bjjDecay[2] == -9999 ) bjjDecay[2] = partonId[i];
          else { cerr << "ERROR: Too many partons for hadronic antitop decay!  Event: " << ievt << endl; continue; }
        }
        
        if ( mc_pdgId[partonId[i]] == 5 && mc_mother[partonId[i]] == pdgID_top ) blvDecay[0] = partonId[i];  // leptonic b
      }
    }  // end loop partons
    blvDecay[1] = genmuon;
    blvDecay[2] = genneutrino;
    
    if ( test && verbose > 3 )
    {
      cout << "bjj " << bjjDecay[0] << "  " << bjjDecay[1] << "  " << bjjDecay[2] << endl;
      cout << "blv " << blvDecay[0] << "  " << blvDecay[1] << "  " << blvDecay[2] << endl;
    }
    
    if ( bjjDecay[0] == -9999 || bjjDecay[1] == -9999 || bjjDecay[2] == -9999 || blvDecay[0] == -9999 || blvDecay[1] == -9999 || blvDecay[2] == -9999)
    {
      if ( test && verbose > 2 ) cout << "Not all ttbar components found... Skipping event " << ievt << "..." << endl;
      continue;
    }
    
    if ( test && verbose > 3 )
    {
      cout << "Hadronic b: " << setw(3) << right << bjjDecay[0] << "  pdgId: " << setw(3) << mc_pdgId[bjjDecay[0]] << "  Mother: " << setw(4) << mc_mother[bjjDecay[0]] << "  Granny: " << setw(4) << mc_granny[bjjDecay[0]] << endl;
      cout << "Light jet1: " << setw(3) << right << bjjDecay[1] << "  pdgId: " << setw(3) << mc_pdgId[bjjDecay[1]] << "  Mother: " << setw(4) << mc_mother[bjjDecay[1]] << "  Granny: " << setw(4) << mc_granny[bjjDecay[1]] << endl;
      cout << "Light jet2: " << setw(3) << right << bjjDecay[2] << "  pdgId: " << setw(3) << mc_pdgId[bjjDecay[2]] << "  Mother: " << setw(4) << mc_mother[bjjDecay[2]] << "  Granny: " << setw(4) << mc_granny[bjjDecay[2]] << endl;
      cout << "Leptonic b: " << setw(3) << right << blvDecay[0] << "  pdgId: " << setw(3) << mc_pdgId[blvDecay[0]] << "  Mother: " << setw(4) << mc_mother[blvDecay[0]] << "  Granny: " << setw(4) << mc_granny[blvDecay[0]] << endl;
      cout << "Muon:       " << setw(3) << right << blvDecay[1] << "  pdgId: " << setw(3) << mc_pdgId[blvDecay[1]] << "  Mother: " << setw(4) << mc_mother[blvDecay[1]] << "  Granny: " << setw(4) << mc_granny[blvDecay[1]] << endl;
      cout << "Neutrino:   " << setw(3) << right << blvDecay[2] << "  pdgId: " << setw(3) << mc_pdgId[blvDecay[2]] << "  Mother: " << setw(4) << mc_mother[blvDecay[2]] << "  Granny: " << setw(4) << mc_granny[blvDecay[2]] << endl;
    }
      
    m_top = (mcParticles[topQuark]).M();
    m_antitop = (mcParticles[antiTopQuark]).M();
    m_bjj = (mcParticles[bjjDecay[0]] + mcParticles[bjjDecay[1]] + mcParticles[bjjDecay[2]]).M();
    m_blv = (mcParticles[blvDecay[0]] + mcParticles[blvDecay[1]] + mcParticles[blvDecay[2]]).M();
    
    //cout << "m_top: " << m_top << ";  m_antitop: " << m_antitop << endl;
    
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
    
    
    /// Loop over different top width scalings
    for (int s = 0; s < nScalings; s++)
    {
      //cout << "Scale " << s << ": " << scaling[s] << " times top width" << endl;
      
      /// Calculate event weights when not nominal
      if ( s > 0 )
      {
        evWeight_hadr = eventWeightCalculator(m_hadr, scaling[s]);
        evWeight_lept = eventWeightCalculator(m_lept, scaling[s]);
        evWeight_bjj = eventWeightCalculator(m_bjj, scaling[s]);
        evWeight_blv = eventWeightCalculator(m_blv, scaling[s]);
      }
      
      nofEvents_hadr[s] += evWeight_hadr;
      nofEvents_lept[s] += evWeight_lept;
      nofEvents_test[s] += (evWeight_hadr+evWeight_lept)/2.;
      
      /// Fill plots before matching
      if (makePlots)
      {
        histo1D[("top_mass_hadr_gen_"+scalingString[s]).c_str()]->Fill(m_hadr, evWeight_hadr);
        histo1D[("top_mass_lept_gen_"+scalingString[s]).c_str()]->Fill(m_lept, evWeight_lept);
        histo1D[("bjj_mass_gen_"+scalingString[s]).c_str()]->Fill(m_bjj, evWeight_bjj);
        histo1D[("blv_mass_gen_"+scalingString[s]).c_str()]->Fill(m_blv, evWeight_blv);
        histo1D[("Width_SF_hadr_"+scalingString[s]).c_str()]->Fill(evWeight_hadr);
        histo1D[("Width_SF_lept_"+scalingString[s]).c_str()]->Fill(evWeight_lept);
        histo2D[("Width_SF_hadr_lept_"+scalingString[s]).c_str()]->Fill(evWeight_hadr, evWeight_lept);
        histo2D[("top_mass_hadr_gen_vs_weight_"+scalingString[s]).c_str()]->Fill(m_hadr, evWeight_hadr);
        histo2D[("top_mass_lept_gen_vs_weight_"+scalingString[s]).c_str()]->Fill(m_lept, evWeight_lept);
        if ( s == 0 )
        {
          histo2D["top_mass_hadr_lept_gen_orig"]->Fill(m_hadr, m_lept);
          histo2D["top_mass_hadr_bjj_gen_orig"]->Fill(m_hadr, m_bjj);
          histo2D["top_mass_lept_blv_gen_orig"]->Fill(m_lept, m_blv);
          histo2D["top_mass_bjj_blv_gen_orig"]->Fill(m_bjj, m_blv);
        }
      }
      
      
      /// Outliers???
      if ( s == 0 && ( m_hadr > 10 && m_hadr < 110 ) && ( m_bjj > 163 && m_bjj < 178 ) )
      {
        cout << endl;
        cout << "Event " << ievt << ": Hadronic top mass much smaller than reconstructed mass from quarks!" << endl;
        cout << "Top quark :    " << "  Status: " << setw(2) << mc_status[topQuark] << "  pdgId: " << setw(3) << mc_pdgId[topQuark] << "  Mother: " << setw(4) << mc_mother[topQuark] << "  Granny: " << setw(4) << mc_granny[topQuark] << "  Pt: " << setw(7) << left << mc_pt[topQuark] << "  Eta: " << mc_eta[topQuark] << "  Mass: " << mc_M[topQuark] << endl;
        cout << "Antitop quark :" << "  Status: " << setw(2) << mc_status[antiTopQuark] << "  pdgId: " << setw(3) << mc_pdgId[antiTopQuark] << "  Mother: " << setw(4) << mc_mother[antiTopQuark] << "  Granny: " << setw(4) << mc_granny[antiTopQuark] << "  Pt: " << setw(7) << left << mc_pt[antiTopQuark] << "  Eta: " << mc_eta[antiTopQuark] << "  Mass: " << mc_M[antiTopQuark] << endl;
        if (muMinusFromTop) cout << "Top decays hadronically, antitop decays leptonically" << endl;
        else if (muPlusFromTop) cout << "Top decays leptonically, antitop decays hadronically" << endl;
        cout << endl;
      }
      
      if ( s == 0 && ( m_lept > 10 && m_lept < 110 ) && ( m_blv > 163 && m_blv < 178 ) )
      {
        cout << endl;
        cout << "Event " << ievt << ": Leptonic top mass much smaller than reconstructed mass from b quark, lepton & neutrino!" << endl;
        cout << "Top quark :    " << "  Status: " << setw(2) << mc_status[topQuark] << "  pdgId: " << setw(3) << mc_pdgId[topQuark] << "  Mother: " << setw(4) << mc_mother[topQuark] << "  Granny: " << setw(4) << mc_granny[topQuark] << "  Pt: " << setw(7) << left << mc_pt[topQuark] << "  Eta: " << mc_eta[topQuark] << "  Mass: " << mc_M[topQuark] << endl;
        cout << "Antitop quark :" << "  Status: " << setw(2) << mc_status[antiTopQuark] << "  pdgId: " << setw(3) << mc_pdgId[antiTopQuark] << "  Mother: " << setw(4) << mc_mother[antiTopQuark] << "  Granny: " << setw(4) << mc_granny[antiTopQuark] << "  Pt: " << setw(7) << left << mc_pt[antiTopQuark] << "  Eta: " << mc_eta[antiTopQuark] << "  Mass: " << mc_M[antiTopQuark] << endl;
        if (muMinusFromTop) cout << "Top decays hadronically, antitop decays leptonically" << endl;
        else if (muPlusFromTop) cout << "Top decays leptonically, antitop decays hadronically" << endl;
        cout << endl;
      }
      
      
      
      if (all4PartonsMatched)
      {
        /// Plot variables for matched events
        double matchedTopMass_gen = (partonsMatched[0] + partonsMatched[1] + partonsMatched[2]).M();
        double matchedTopMass_reco = (jetsMatched[0] + jetsMatched[1] + jetsMatched[2]).M();
        
        if ( s > 0 )
        {
          eventSF_gen = eventWeightCalculator(matchedTopMass_gen, scaling[s]);
          eventSF_reco = eventWeightCalculator(matchedTopMass_gen, scaling[s]);
        }
        else
        {
          eventSF_gen = 1.;
          eventSF_reco = 1.;
        }
        
        if (makePlots)
        {
          histo1D[("top_mass_gen_matched_"+scalingString[s]).c_str()]->Fill(matchedTopMass_gen, eventSF_gen);
          histo1D[("top_mass_reco_matched_"+scalingString[s]).c_str()]->Fill(matchedTopMass_reco, eventSF_reco);
        }
        
        
      }  // end all4PartonsMatched
    
    
    }  // end loop scalings
    
    
    
  }  // end loop events
  
  
  cout << endl;
  cout << "Number of events passing hard selection: " << nofHardSelected << " (" << 100*((float)nofHardSelected/(float)endEvent) << "%)" << endl;
  cout << "Number of skipped events: " << nofSkippedEvents << " (" << 100*((float)nofSkippedEvents/(float)nofHardSelected) << "%)" << endl;
  cout << "Number of matched events: " << nofMatchedEvents << endl;
  cout << endl;
  cout << "Number of events without scaling hadronic width: " << nofEvents_hadr[0] << "; hadr width x 2: " << nofEvents_hadr[1] << "; hadr width x 0.5: " << nofEvents_hadr[2] << endl;
  cout << "Number of events without scaling leptonic width: " << nofEvents_lept[0] << "; lept width x 2: " << nofEvents_lept[1] << "; lept width x 0.5: " << nofEvents_lept[2] << endl;
  cout << "Number of events without scaling width:          " << nofEvents_test[0] << "; test width x 2: " << nofEvents_test[1] << "; test width x 0.5: " << nofEvents_test[2] << endl;
  
  
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
    //TDirectory* th1dir = fout->mkdir("1D_histograms");
    //th1dir->cd();
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
    
    fout->cd();
    TF1 *BW = new TF1("Breit_Wigner", BreitWigner2Func, 120., 220., 2);
    std::vector<double> scales = {0.5, 1., 2.};
    DrawBW(BW, pathOutput+"Breit_Wigner", scales, 120, 220, true);
    
    for (int s = 0; s < nScalings; s++)
    {
      string histoName = "top_mass_gen_matched_"+scalingString[s];
      cout << "Integral for SM width times " << scaling[s] << ":  "  << histo1D[histoName.c_str()]->Integral(0, histo1D[histoName.c_str()]->GetNbinsX()+1) << endl;
      FitBW(histo1D[histoName], pathOutput+histoName+"_withBW", scaling[s]);
      //DrawBWonHisto(histo1D[histoName], BW, pathOutput+histoName+"_withBW", scaling[s], 1., 120., 220., true);
    }
    
    
    fout->Close();
    
    delete fout;
  }
  
  
  
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


void GetMetaData(TTree* tree)
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
  
  tree->SetBranchStatus("*",1);
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
  
  tree->SetBranchStatus("*",1);
}


void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  /// SFs
  histo1D["Width_SF_hadr_orig"] = new TH1F("Width_SF_hadr_orig", "Hadronic scale factor to change the ttbar distribution width (no scaling); width SF", 501, 0.4985, 2.0015);
  histo1D["Width_SF_hadr_s2"] = new TH1F("Width_SF_hadr_s2", "Hadronic scale factor to change the ttbar distribution width (scaled by factor 2); width SF", 501, 0.4985, 2.0015);
  histo1D["Width_SF_hadr_s0p5"] = new TH1F("Width_SF_hadr_s0p5", "Hadronic scale factor to change the ttbar distribution width (scaled by factor 0.5); width SF", 501, 0.4985, 2.0015);
  histo1D["Width_SF_lept_orig"] = new TH1F("Width_SF_lept_orig", "Leptonic scale factor to change the ttbar distribution width (no scaling); width SF", 501, 0.4985, 2.0015);
  histo1D["Width_SF_lept_s2"] = new TH1F("Width_SF_lept_s2", "Leptonic scale factor to change the ttbar distribution width (scaled by factor 2); width SF", 501, 0.4985, 2.00152);
  histo1D["Width_SF_lept_s0p5"] = new TH1F("Width_SF_lept_s0p5", "Leptonic scale factor to change the ttbar distribution width (scaled by factor 0.5); width SF", 501, 0.4985, 2.0015);
  
  histo1D["top_mass_hadr_gen_orig"] = new TH1F("top_mass_hadr_gen_orig","Mass of generated top quark with hadronic decay (no scaling); M_{t_{hadr}} [GeV]", 4000, 120, 220);
  histo1D["top_mass_hadr_gen_s2"] = new TH1F("top_mass_hadr_gen_s2","Mass of generated top quark with hadronic decay (scaled by factor 2); M_{t_{hadr}} [GeV]", 4000, 120, 220);
  histo1D["top_mass_hadr_gen_s0p5"] = new TH1F("top_mass_hadr_gen_s0p5","Mass of generated top quark with hadronic decay (scaled by factor 0.5); M_{t_{hadr}} [GeV]", 4000, 120, 220);
  
  histo1D["top_mass_lept_gen_orig"] = new TH1F("top_mass_lept_gen_orig","Mass of generated top quark with leptonic decay (no scaling); M_{t_{lept}} [GeV]", 4000, 120, 220);
  histo1D["top_mass_lept_gen_s2"] = new TH1F("top_mass_lept_gen_s2","Mass of generated top quark with leptonic decay (scaled by factor 2); M_{t_{lept}} [GeV]", 4000, 120, 220);
  histo1D["top_mass_lept_gen_s0p5"] = new TH1F("top_mass_lept_gen_s0p5","Mass of generated top quark with leptonic decay (scaled by factor 0.5); M_{t_{lept}} [GeV]", 4000, 120, 220);
  
  histo1D["bjj_mass_gen_orig"] = new TH1F("bjj_mass_gen_orig","Mass of generated bjj quarks (no scaling); M_{bjj} [GeV]", 2000, 50, 300);
  histo1D["bjj_mass_gen_s2"] = new TH1F("bjj_mass_gen_s2","Mass of generated bjj quarks (scaled by factor 2); M_{bjj} [GeV]", 2000, 50, 300);
  histo1D["bjj_mass_gen_s0p5"] = new TH1F("bjj_mass_gen_s0p5","Mass of generated bjj quarks (scaled by factor 0.5); M_{bjj} [GeV]", 2000, 50, 300);
  
  histo1D["blv_mass_gen_orig"] = new TH1F("blv_mass_gen_orig","Mass of generated b, lepton and neutrino (no scaling); M_{blv} [GeV]", 2000, 50, 300);
  histo1D["blv_mass_gen_s2"] = new TH1F("blv_mass_gen_s2","Mass of generated b, lepton and neutrino (scaled by factor 2); M_{blv} [GeV]", 2000, 50, 300);
  histo1D["blv_mass_gen_s0p5"] = new TH1F("blv_mass_gen_s0p5","Mass of generated b, lepton and neutrino (scaled by factor 0.5); M_{blv} [GeV]", 2000, 50, 300);
  
  
  /// Matching
  histo1D["top_mass_reco_matched_orig"] = new TH1F("top_mass_reco_matched_orig","Reconstructed top mass of matched events (no scaling); M_{t} [GeV]", 400, 0, 400);
  histo1D["top_mass_reco_matched_s2"] = new TH1F("top_mass_reco_matched_s2","Reconstructed top mass of matched events (scaled by factor 2); M_{t} [GeV]", 400, 0, 400);
  histo1D["top_mass_reco_matched_s0p5"] = new TH1F("top_mass_reco_matched_s0p5","Reconstructed top mass of matched events (scaled by factor 0.5); M_{t} [GeV]", 400, 0, 400);
  
  histo1D["top_mass_gen_matched_orig"] = new TH1F("top_mass_gen_matched_orig","Generated top mass of matched events (no scaling); M_{t} [GeV]", 800, 50, 300);
  histo1D["top_mass_gen_matched_s2"] = new TH1F("top_mass_gen_matched_s2","Generated top mass of matched events (scaled by factor 2); M_{t} [GeV]", 800, 50, 300);
  histo1D["top_mass_gen_matched_s0p5"] = new TH1F("top_mass_gen_matched_s0p5","Generated top mass of matched events (scaled by factor 0.5); M_{t} [GeV]", 800, 50, 300);
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
  /// SFs
  histo2D["Width_SF_hadr_lept_orig"] = new TH2F("Width_SF_hadr_lept_orig", "Leptonic vs. hadronic scale factor to change the ttbar distribution width; hadr. width SF (no scaling); lept. width SF", 501, 0.4985, 2.0015, 501, 0.4985, 2.0015);
  histo2D["Width_SF_hadr_lept_s2"] = new TH2F("Width_SF_hadr_lept_s2", "Leptonic vs. hadronic scale factor to change the ttbar distribution width (scaled by factor 2); hadr. width SF; lept. width SF", 501, 0.4985, 2.0015, 501, 0.4985, 2.0015);
  histo2D["Width_SF_hadr_lept_s0p5"] = new TH2F("Width_SF_hadr_lept_s0p5", "Leptonic vs. hadronic scale factor to change the ttbar distribution width (scaled by factor 0.5); hadr. width SF; lept. width SF", 501, 0.4985, 2.0015, 501, 0.4985, 2.0015);
  
  histo2D["top_mass_hadr_lept_gen_orig"] = new TH2F("top_mass_hadr_lept_gen_orig","Mass of generated top quark with leptonic decay vs. hadronic decay (no scaling); M_{t_{hadr}} [GeV]; M_{t_{lept}} [GeV]", 1000, 50, 300, 1000, 50, 300);
  histo2D["top_mass_hadr_bjj_gen_orig"] = new TH2F("top_mass_hadr_bjj_gen_orig","Generated mass of bjj quarks vs. hadronically decaying top quark mass (no scaling); M_{t_{hadr}} [GeV]; M_{bjj} [GeV]", 1000, 50, 300, 1000, 50, 300);
  histo2D["top_mass_lept_blv_gen_orig"] = new TH2F("top_mass_lept_blv_gen_orig","Mass of generated b, lepton and neutrino vs. leptonically decaying top quark mass (no scaling); M_{t_{lept}} [GeV]; M_{blv} [GeV]", 1000, 50, 300, 1000, 50, 300);
  histo2D["top_mass_bjj_blv_gen_orig"] = new TH2F("top_mass_bjj_blv_gen_orig","Mass of generated b, lepton and neutrino vs. mass of bjj quarks (no scaling); M_{bjj} [GeV]; M_{blv} [GeV]", 1000, 50, 300, 1000, 50, 300);
  
  // Weights
  histo2D["top_mass_hadr_gen_vs_weight_orig"] = new TH2F("top_mass_hadr_gen_vs_weight_orig","Weights vs. mass of generated top quark (hadronic decay); M_{t_{hadr}} [GeV]; weight", 1000, 120, 220, 501, -0.0025, 2.5025);
  histo2D["top_mass_hadr_gen_vs_weight_s2"] = new TH2F("top_mass_hadr_gen_vs_weight_s2","Weights vs. mass of generated top quark (hadronic decay); M_{t_{hadr}} [GeV]; weight", 1000, 120, 220, 501, -0.0025, 2.5025);
  histo2D["top_mass_hadr_gen_vs_weight_s0p5"] = new TH2F("top_mass_hadr_gen_vs_weight_s0p5","Weights vs. mass of generated top quark (hadronic decay); M_{t_{hadr}} [GeV]; weight", 1000, 120, 220, 501, -0.0025, 2.5025);
  
  histo2D["top_mass_lept_gen_vs_weight_orig"] = new TH2F("top_mass_lept_gen_vs_weight_orig","Weights vs. mass of generated top quark (leptonic decay); M_{t_{lept}} [GeV]; weight", 1000, 120, 220, 501, -0.0025, 2.5025);
  histo2D["top_mass_lept_gen_vs_weight_s2"] = new TH2F("top_mass_lept_gen_vs_weight_s2","Weights vs. mass of generated top quark (leptonic decay); M_{t_{lept}} [GeV]; weight", 1000, 120, 220, 501, -0.0025, 2.5025);
  histo2D["top_mass_lept_gen_vs_weight_s0p5"] = new TH2F("top_mass_lept_gen_vs_weight_s0p5","Weights vs. mass of generated top quark (leptonic decay); M_{t_{lept}} [GeV]; weight", 1000, 120, 220, 501, -0.0025, 2.5025);
}

void TruthMatching(vector<TLorentzVector> partons, vector<TLorentzVector> selectedJets, pair<unsigned int, unsigned int> *MCPermutation)  /// MCPermutation: 0,1 hadronic W jet; 2 hadronic b jet; 3 leptonic b jet
{
  JetPartonMatching matching = JetPartonMatching(partons, selectedJets, 2, true, true, 0.3);  // partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
  
  if (matching.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matching.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle/parton number
  JetPartonPair.clear();
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matching.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
  }

  if (verbose > 2)
    cout << "Matching done" << endl;
  
  
  for (unsigned int iPair = 0; iPair < JetPartonPair.size(); iPair++)
  {
    unsigned int j = JetPartonPair[iPair].second;
    
    if ( fabs(mc_pdgId[partonId[j]]) < 6 )  //light/b quarks, 6 should stay hardcoded
    {
      if ( ( muPlusFromTop && mc_mother[partonId[j]] == -24 && mc_granny[partonId[j]] == -pdgID_top )
          || ( muMinusFromTop && mc_mother[partonId[j]] == 24 && mc_granny[partonId[j]] == pdgID_top ) )  // if mu+, check if mother of particle is W- and granny tbar --> then it is a quark from W- decay
      {
        if (verbose > 3)
          cout << "Light jet: " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;
        
        if (MCPermutation[0].first == 9999)
        {
          MCPermutation[0] = JetPartonPair[iPair];
        }
        else if (MCPermutation[1].first == 9999)
        {
          MCPermutation[1] = JetPartonPair[iPair];
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
        
        MCPermutation[2] = JetPartonPair[iPair];
      }
      else if ( ( muPlusFromTop && mc_mother[partonId[j]] == pdgID_top )
        || ( muMinusFromTop && mc_mother[partonId[j]] == -pdgID_top ) )  // if mu+ (top decay leptonic) and mother is top ---> leptonic b
      {
        if (verbose > 3)
          cout << "b jet:     " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;
        
        MCPermutation[3] = JetPartonPair[iPair];
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
  
  nofMatchedEvents = 0;
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
  for (int i = 0; i < 4; i++)
  {
    MCPermutation[i] = pair<unsigned int,unsigned int>(9999,9999);
  }
  topQuark = -9999;
  antiTopQuark = -9999;
  genmuon = -9999;
  genneutrino = -9999;
  muonmatched = false;
  muPlusFromTop = false;
  muMinusFromTop = false;
  for (int i = 0; i < 3; i++)
  {
    bjjDecay[i] = -9999;
    blvDecay[i] = -9999;
  }
  partonId.clear();
}

void ClearVars()
{
  ClearMatching();
  
  scaleFactor = 1.;
  eventSF_gen = 1.;
  eventSF_reco = 1.;
  evWeight_hadr = 1.;
  evWeight_lept = 1.;
  evWeight_bjj = 1.;
  evWeight_blv = 1.;
  m_top = -1.;
  m_antitop = -1.;
  m_bjj = -1.;
  m_blv = -1.;
  m_hadr = -1;
  m_lept = -1;
  bJetId.clear();
  selectedJetsCharge.clear();
  selectedJetsBDiscr.clear();
  reco_hadWMass = -1.;
  reco_hadTopMass = -1.;
  reco_hadTopPt = -1.;
  reco_minMlb = 9999.;
  reco_ttbarMass = -1.;
  reco_dRLepB = -1.;
  min_Mlb = 9999.;
  dRLepB = -1.;
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

Double_t BreitWigner2Func(Double_t *x, Double_t *scale)
{
  double BWmass = genTopMass;
  double BWgamma = scale[0]*genTopWidth;
  //double numerator = 2.*pow(BWmass*BWgamma, 2)/TMath::Pi();
  double gammaterm = sqrt( pow(BWmass, 4) + pow(BWmass*BWgamma, 2) );
  double numerator = 2.*sqrt(2.)*BWmass*BWgamma*gammaterm/( TMath::Pi()*sqrt( pow(BWmass, 2) + gammaterm ) );
  double denominator = pow(pow(x[0], 2) - pow(BWmass, 2), 2) + pow(x[0], 4)*pow(BWgamma/BWmass, 2);
  
  return numerator/denominator*scale[1];
}

double BreitWigner(double topMass, double scale)
{
  double BWmass = genTopMass; // or topMass (reco)?
  double BWgamma = scale*genTopWidth;
  double gammaterm = sqrt( pow(BWmass, 4) + pow(BWmass*BWgamma, 2) );
  double numerator = 2.*sqrt(2.)*BWmass*BWgamma*gammaterm/( TMath::Pi()*sqrt( pow(BWmass, 2) + gammaterm ) );
  double denominator = pow(pow(topMass, 2) - pow(BWmass, 2), 2) + pow(BWmass*BWgamma, 2);
  
  return numerator/denominator;
}

double BreitWigner2(double topMass, double scale)
{
  double BWmass = genTopMass;
  double BWgamma = scale*genTopWidth;
  //double numerator = 2.*pow(BWmass*BWgamma, 2)/TMath::Pi();
  double gammaterm = sqrt( pow(BWmass, 4) + pow(BWmass*BWgamma, 2) );
  double numerator = 2.*sqrt(2.)*BWmass*BWgamma*gammaterm/( TMath::Pi()*sqrt( pow(BWmass, 2) + gammaterm ) );
  double denominator = pow(pow(topMass, 2) - pow(BWmass, 2), 2) + pow(topMass, 4)*pow(BWgamma/BWmass, 2);
  
  return numerator/denominator;
}

double BreitWignerNonRel(double topMass, double scale)
{
  double BWmass = genTopMass;
  double BWgamma = scale*genTopWidth/2.;

  double bw = BWgamma/( pow(topMass - BWmass, 2) + pow(BWgamma, 2) );
  
  return bw/(TMath::Pi());
}

double eventWeightCalculator(double topMass, double scale)
{
  return BreitWigner2(topMass, scale)/BreitWigner2(topMass, 1.);
}

void DrawBW(TF1* function, string name, std::vector<double> scales, double min, double max, bool writeToFile)
{
  Color_t colours[] = {kRed, kBlue+2, kGreen-7};
  //Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+1, kCyan+1, kBlue+2, kMagenta, kViolet-5, kPink+10};
  TCanvas* c2 = new TCanvas(name.c_str(), name.c_str());
  c2->cd();
  function->FixParameter(0, scales[0]);
  function->FixParameter(1, 1.);
  function->SetLineColor(colours[0]);
  function->SetRange(min, max);
  function->DrawCopy();
  c2->Update();
  for (int i = 1; i < scales.size(); i++)
  {
    function->FixParameter(0,scales[i]);
    function->SetLineColor(colours[i]);
    function->DrawCopy("same");
    c2->Update();
  }
  if (writeToFile) c2->Write();
  c2->SaveAs((name+".png").c_str());
  c2->Close();
  
  delete c2;
}

void FitBW(TH1F* histoIn, string name, double scale)
{
  cout << "Fitting " << histoIn->GetName() << " for widthx" << scale << endl;
  
  // Normalise histogram
  TH1F *histo = (TH1F*) histoIn->Clone(histoIn->GetName());
  histo->Integral(0, histo->GetNbinsX()+1);
  histo->Scale(scale);
  
  TF1 *BWfit = new TF1("Breit_Wigner", BreitWigner2Func, 170., 174., 2);
  BWfit->SetParNames("widthx", "norm");
  BWfit->FixParameter(0, scale);
  
  std::string func_title = std::string(histo->GetName())+"_Fitted";
  BWfit->SetName(func_title.c_str());
  histo->Fit(BWfit, "R");
  gStyle->SetOptFit(0111);
  histo->SetStats(1);
  histo->Write();
  BWfit->Write();
  
  DrawBWonHisto(histoIn, BWfit, name, BWfit->GetParameter(0), BWfit->GetParameter(1), 120., 220., true);
}

void DrawBWonHisto(TH1F* histo, TF1* function, string name, double scale, double norm, double min, double max, bool writeToFile)
{
  Color_t colours[] = {kBlue+2, kRed};
  TCanvas* c3 = new TCanvas(name.c_str(), name.c_str());
  c3->cd();
  histo->GetXaxis()->SetRangeUser(min, max);
  histo->SetLineStyle(1);
  histo->SetLineColor(colours[0]);
  histo->Draw();
  
  function->FixParameter(0, scale);
  function->FixParameter(1, norm);
  function->SetLineColor(colours[1]);
  function->SetRange(min, max);
  function->DrawCopy("same");
  c3->Update();
  
  if (writeToFile) c3->Write();
  c3->SaveAs((name+".png").c_str());
  c3->Close();
  
  delete c3;
}
