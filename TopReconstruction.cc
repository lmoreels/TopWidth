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
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
//#include "../macros/Style.C"



using namespace std;
using namespace TopTree;


bool test = false;
string systStr = "nominal";

string ntupleDate = "160812"; // nominal
//string ntupleDate = "160916"; // JERup
int verbose = 2;

string pathNtuples = "";

int nofMatchedEvents = 0;

const int nofCats = 5;
int nofCorrectlyMatched_deltaRAll[nofCats] = {0, 0, 0, 0, 0};  // 0: All, 1: == 4jets, 2: == 5jets, 3: == 6jets, 4: > 6jets
int nofNotCorrectlyMatched_deltaRAll[nofCats] = {0, 0, 0, 0, 0};
int nofCorrectlyMatched_deltaR1B[nofCats] = {0, 0, 0, 0, 0};
int nofNotCorrectlyMatched_deltaR1B[nofCats] = {0, 0, 0, 0, 0};
int nofCorrectlyMatched_chi2All[nofCats] = {0, 0, 0, 0, 0};
int nofNotCorrectlyMatched_chi2All[nofCats] = {0, 0, 0, 0, 0};
int nofCorrectlyMatched_chi2W[nofCats] = {0, 0, 0, 0, 0};
int nofNotCorrectlyMatched_chi2W[nofCats] = {0, 0, 0, 0, 0};
int nofCorrectlyMatched_chi2W1B[nofCats] = {0, 0, 0, 0, 0};
int nofNotCorrectlyMatched_chi2W1B[nofCats] = {0, 0, 0, 0, 0};
int nofCorrectlyMatched_chi2WDeltaRW[nofCats] = {0, 0, 0, 0, 0};
int nofNotCorrectlyMatched_chi2WDeltaRW[nofCats] = {0, 0, 0, 0, 0};
int nofCorrectlyMatched_chi2WDeltaRW1B[nofCats] = {0, 0, 0, 0, 0};
int nofNotCorrectlyMatched_chi2WDeltaRW1B[nofCats] = {0, 0, 0, 0, 0};

int nofCorrectlyMatched_comb = 0;
int nofNotCorrectlyMatched_comb = 0;


///  Working points for b tagging  // Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
float CSVv2Loose =  0.460;
float CSVv2Medium = 0.800;
float CSVv2Tight = 0.935;

// Temporarily, until calculated from TTbar sample
float chi2WMass = 80.385;
float sigmaChi2WMass = 10;
float chi2TopMass = 172.5; //180.0; //from mtop mass plot: 167.0
float sigmaChi2TopMass = 40;

// Average top mass
// TT match, TT chi2, ST_t_top, ST_t_antitop, ST_tW_top, ST_tW_antitop, DYJets, WJets, data, allChi2
float aveTopMass[] = {166.922, 170.433, 178.150, 203.137, 201.759, 185.668, 186.755, 189.150, 182.008, 179.775, 178.315};

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;


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
void ClearObjects();



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


float scaleFactor;
vector<unsigned int> bJetId;
double reco_hadWMass, reco_hadTopMass, reco_hadTopPt;
double reco_minMlb, reco_ttbarMass, reco_dRLepB;
double min_Mlb, dRLepB;

/// Define TLVs
TLorentzVector muon, jet, mcpart;
TLorentzVector topCandidate, WCandidate;
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

bool all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
bool all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
bool hadronictopJetsMatched_MCdef_ = false;
pair<unsigned int, unsigned int> MCPermutation[4] = {pair<unsigned int,unsigned int>(9999,9999)};
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
  if ( argc == 1 ) channel = "mu";
  else if ( argv[1] == "mu" || argv[1] == "Mu" || argv[1] == "MU" || argv[1] == "m" ) channel = "mu";
  else if ( argv[1] == "el" || argv[1] == "El" || argv[1] == "EL" || argv[1] == "e" ) channel = "el";
  //else if ( argv[1] == "all" || argv[1] == "All" || argv[1] == "ALL" ) channel = "all";
  
//   //string pathOutput = "test/";
//   string pathOutput = "OutputPlots/";
//   mkdir(pathOutput.c_str(),0777);
//   pathOutput += "topreco/";
//   mkdir(pathOutput.c_str(),0777);
//   // Add channel to output path
//   pathOutput += channel+"/";
//   mkdir(pathOutput.c_str(),0777);
//   if (test)
//   {
//     pathOutput += "test/";
//     mkdir(pathOutput.c_str(),0777);
//   }
//   // Give timestamp to output path
//   pathOutput += dateString+"/";
//   mkdir(pathOutput.c_str(),0777);
  
  
  pathNtuples = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate+"/";
  cout << "Using Ntuples from " << ntupleDate << ". This corresponds to systematics: " << systStr << endl;
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  InitHisto1D();
  InitHisto2D();
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
  if (test) endEvent = 1001;
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
      if (test) cout << "Both tops decay leptonically... Something went wrong... Event " << ievt << endl;
      continue;
    }
    
    
    TruthMatching(partons, selectedJets, MCPermutation);
    
    if ( all4PartonsMatched && test && verbose > 3 )
    {
      for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
      {
        //cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Status: " << setw(2) << mc_status[partonId[MCPermutation[iMatch].second]] << "  pdgId: " << setw(3) << mc_pdgId[partonId[MCPermutation[iMatch].second]] << "  Mother: " << setw(4) << mc_mother[partonId[MCPermutation[iMatch].second]] << "  Granny: " << setw(4) << mc_granny[partonId[MCPermutation[iMatch].second]] << endl;
        cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Pt: " << setw(7) << left << mc_pt[partonId[MCPermutation[iMatch].second]] << "  Eta: " << mc_eta[partonId[MCPermutation[iMatch].second]] << "  Phi: " << mc_phi[partonId[MCPermutation[iMatch].second]] << endl;
        cout << "Event  " << right << setw(4) << ievt << ";  Matched jet    " << iMatch << "  Pt: " << setw(7) << left << jet_pt[MCPermutation[iMatch].first] << "  Eta: " << jet_eta[MCPermutation[iMatch].first] << "  Phi: " << jet_phi[MCPermutation[iMatch].first] << endl;
      }
    }
    
    
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


      /// Plot variables for matched events
      float matchedWMass_reco = (jetsMatched[0] + jetsMatched[1]).M();
      float matchedTopMass_reco = (jetsMatched[0] + jetsMatched[1] + jetsMatched[2]).M();
      float matchedTopMass_gen = (partonsMatched[0] + partonsMatched[1] + partonsMatched[2]).M();

//       histo1D["W_mass_reco_matched"]->Fill(matchedWMass_reco);
//       histo1D["top_mass_reco_matched"]->Fill(matchedTopMass_reco);
//       histo1D["top_mass_gen_matched"]->Fill(matchedTopMass_gen);

      if (muonmatched)
      {
        float matchedMlb_corr = (selectedLepton[0] + jetsMatched[3]).M();  // lept b
        float matchedMlb_wrong = (selectedLepton[0] + jetsMatched[2]).M();  // hadr b
        float matchedTtbarMass_corr = matchedTopMass_reco + matchedMlb_corr;
        float matchedTtbarMass_wrong = matchedTopMass_reco + matchedMlb_wrong;
        double matchedDRLepB_corr = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[3]);  // lept b
        double matchedDRLepB_wrong = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[2]);  // hadr b
        
//         histo1D["mlb_matched_corr"]->Fill(matchedMlb_corr);
//         histo1D["mlb_matched_wrong"]->Fill(matchedMlb_wrong);
//         histo1D["ttbar_mass_matched_corr"]->Fill(matchedTtbarMass_corr);
//         histo1D["ttbar_mass_matched_wrong"]->Fill(matchedTtbarMass_wrong);
//         histo1D["dR_lep_b_matched_corr"]->Fill(matchedDRLepB_corr);
//         histo1D["dR_lep_b_matched_wrong"]->Fill(matchedDRLepB_wrong);
        
      }  // end muonMatched

    }  // end all4PartonsMatched



    /////////////////////////////////////
    ///  Reconstruction of Top quark  ///
    /////////////////////////////////////
    
    // Labels: 0,1 = light jets, 2 = hadronic b-jet
    int labelsRecoDeltaRAll[3] = {-9999, -9999, -9999};
    int labelsRecoDeltaR1B[3] = {-9999, -9999, -9999};
    int labelsRecoChi2All[3] = {-9999, -9999, -9999};
    int labelsRecoChi2W[3] = {-9999, -9999, -9999};
    int labelsRecoChi2W1B[3] = {-9999, -9999, -9999};
    int labelsRecoChi2WDeltaRW[3] = {-9999, -9999, -9999};
    int labelsRecoChi2WDeltaRW1B[3] = {-9999, -9999, -9999};
    
    double deltaRAll, deltaR1B, deltaR;
    double minDeltaRAll = 9999., minDeltaR1B = 9999., minDeltaR = 9999.;
    
    
    ///------DELTA R ALL JETS------///
    
    for (int ijet = 0; ijet < selectedJets.size(); ijet++)
    {
      for (int jjet = ijet+1; jjet < selectedJets.size(); jjet++)
      {
        for (int kjet = jjet+1; kjet < selectedJets.size(); kjet++)
        {
          topCandidate = selectedJets[ijet] + selectedJets[jjet] + selectedJets[kjet];
          deltaRAll = sqrt( pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[ijet], topCandidate), 2) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[jjet], topCandidate), 2) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[kjet], topCandidate), 2) );
          
          if (deltaRAll < minDeltaRAll)
          {
            minDeltaRAll = deltaRAll;
            labelsRecoDeltaRAll[0] = ijet;
            labelsRecoDeltaRAll[1] = jjet;
            labelsRecoDeltaRAll[2] = kjet;
          }
        }
      }
    }
    
    
    
    ///------DELTA R 1 B JET------///
    
    for (int ijet = 0; ijet < selectedJets.size(); ijet++)
    {
      for (int jjet = ijet+1; jjet < selectedJets.size(); jjet++)
      {
        for (int kjet = jjet+1; kjet < selectedJets.size(); kjet++)
        {
          if ( jet_bdiscr[ijet] < CSVv2Medium && jet_bdiscr[jjet] < CSVv2Medium && jet_bdiscr[kjet] < CSVv2Medium ) continue;  // no b jet
          
          topCandidate = selectedJets[ijet] + selectedJets[jjet] + selectedJets[kjet];
          deltaR1B = sqrt( pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[ijet], topCandidate), 2) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[jjet], topCandidate), 2) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[kjet], topCandidate), 2) );
          
          if (deltaR1B < minDeltaR1B)
          {
            minDeltaR1B = deltaR1B;
            labelsRecoDeltaR1B[0] = ijet;
            labelsRecoDeltaR1B[1] = jjet;
            labelsRecoDeltaR1B[2] = kjet;
          }
        }
      }
    }
    
    
    
    ///------CHI2 ALL JETS---------///
    
    float recoWMass, recoTopMass, WTerm, topTerm, chi2;
    float smallestChi2 = 9999.;

    for (int ijet = 0; ijet < selectedJets.size(); ijet++)
    {
      for (int jjet = ijet+1; jjet < selectedJets.size(); jjet++)
      {
        for (int kjet = 0; kjet < selectedJets.size(); kjet++)
        {
          if ( ijet == kjet || jjet == kjet ) continue;
          
          recoWMass = (selectedJets[ijet] + selectedJets[jjet]).M();
          recoTopMass = (selectedJets[ijet] + selectedJets[jjet] + selectedJets[kjet]).M();
          
          WTerm = (recoWMass - chi2WMass)/sigmaChi2WMass;
          topTerm = (recoTopMass - chi2TopMass)/sigmaChi2TopMass;
          
          chi2 = pow(WTerm, 2) + pow(topTerm, 2);
          
          if ( chi2 < smallestChi2)
          {
            smallestChi2 = chi2;
            labelsRecoChi2All[0] = ijet;
            labelsRecoChi2All[1] = jjet;
            labelsRecoChi2All[2] = kjet;
          }
        }
      }
    }


//     if ( labelsRecoChi2All[0] != -9999 && labelsRecoChi2All[1] != -9999 && labelsRecoChi2All[2] != -9999 )
//     {
//       reco_hadWMass = (selectedJets[labelsRecoChi2All[0]] + selectedJets[labelsRecoChi2All[1]]).M();
//       reco_hadTopMass = (selectedJets[labelsRecoChi2All[0]] + selectedJets[labelsRecoChi2All[1]] + selectedJets[labelsRecoChi2All[2]]).M();
//       reco_hadTopPt = (selectedJets[labelsRecoChi2All[0]] + selectedJets[labelsRecoChi2All[1]] + selectedJets[labelsRecoChi2All[2]]).Pt();
// 
//       
//       histo1D["Chi2_W_mass_reco"]->Fill(reco_hadWMass);
//       histo1D["Chi2_top_mass_reco"]->Fill(reco_hadTopMass);
//       
// 
//       /// Leptonic top mass
//       for (int i = 0; i < selectedJets.size(); i++)
//       {
//         if ( i != labelsRecoChi2All[0] && i != labelsRecoChi2All[1] && i != labelsRecoChi2All[2] 
//             && jet_bdiscr[i] > CSVv2Medium )
//           bJetsAfterChi2.push_back(selectedJets[i]);
//       }
// 
//       if ( bJetsAfterChi2.size() > 0 )
//       {
//         double reco_Mlb_temp = 99999.;
//         for (unsigned int i = 0; i < bJetsAfterChi2.size(); i++)
//         {
//           reco_Mlb_temp = (selectedLepton[0] + bJetsAfterChi2[i]).M();
//           if ( reco_Mlb_temp < reco_minMlb )
//           {
//             reco_minMlb = reco_Mlb_temp;
//             reco_dRLepB = ROOT::Math::VectorUtil::DeltaR( bJetsAfterChi2[i], selectedLepton[0]);
//           }
//         }
//         reco_ttbarMass = reco_minMlb + reco_hadTopMass;
//         histo1D["dR_lep_b_unmatched_chi2"]->Fill(reco_dRLepB);
//         
// 
//       }  // end has bjet not in chi2
// 
//     }  // end labels
    
    
    
    ///------CHI2 W BOSON----------///
    
    smallestChi2 = 9999.;
    for (int ijet = 0; ijet < selectedJets.size(); ijet++)
    {
      for (int jjet = ijet+1; jjet < selectedJets.size(); jjet++)
      {
        recoWMass = (selectedJets[ijet] + selectedJets[jjet]).M();
        chi2 = pow( (recoWMass - chi2WMass)/sigmaChi2WMass, 2);
        
        if ( chi2 < smallestChi2)
        {
          smallestChi2 = chi2;
          labelsRecoChi2W[0] = ijet;
          labelsRecoChi2W[1] = jjet;
        }
      }
    }
    
    for (int kjet = 0; kjet < selectedJets.size(); kjet++)
    {
      if ( kjet == labelsRecoChi2W[0] || kjet == labelsRecoChi2W[1] ) continue;
      
      topCandidate = selectedJets[labelsRecoChi2W[0]] + selectedJets[labelsRecoChi2W[1]] + selectedJets[kjet];
      deltaR = sqrt( pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[labelsRecoChi2W[0]], topCandidate), 2) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[labelsRecoChi2W[1]], topCandidate), 2) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[kjet], topCandidate), 2) );

      if (deltaR < minDeltaR)
      {
        minDeltaR = deltaR;
        labelsRecoChi2W[2] = kjet;
      }
    }
    
    ///------CHI2 W BOSON 1B-------///
    labelsRecoChi2W1B[0] = labelsRecoChi2W[0];
    labelsRecoChi2W1B[1] = labelsRecoChi2W[1];
    
    minDeltaR = 9999.;
    for (int kjet = 0; kjet < selectedBJets.size(); kjet++)
    {
      if ( bJetId[kjet] == labelsRecoChi2W1B[0] || bJetId[kjet] == labelsRecoChi2W1B[1] ) continue;
      
      topCandidate = selectedJets[labelsRecoChi2W[0]] + selectedJets[labelsRecoChi2W[1]] + selectedJets[bJetId[kjet]];
      deltaR = sqrt( pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[labelsRecoChi2W[0]], topCandidate), 2) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[labelsRecoChi2W[1]], topCandidate), 2) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[bJetId[kjet]], topCandidate), 2) );

      if (deltaR < minDeltaR)
      {
        minDeltaR = deltaR;
        labelsRecoChi2W1B[2] = bJetId[kjet];
      }
    }
    
    ///------CHI2 W BOSON DELTA R W-------///
    labelsRecoChi2WDeltaRW[0] = labelsRecoChi2W[0];
    labelsRecoChi2WDeltaRW[1] = labelsRecoChi2W[1];
    
    minDeltaR = 9999.;
    WCandidate = selectedJets[labelsRecoChi2WDeltaRW[0]] + selectedJets[labelsRecoChi2WDeltaRW[1]];
    for (int kjet = 0; kjet < selectedJets.size(); kjet++)
    {
      if ( kjet == labelsRecoChi2WDeltaRW[0] || kjet == labelsRecoChi2WDeltaRW[1] ) continue;
      
      deltaR = ROOT::Math::VectorUtil::DeltaR(selectedJets[kjet], WCandidate);

      if (deltaR < minDeltaR)
      {
        minDeltaR = deltaR;
        labelsRecoChi2WDeltaRW[2] = kjet;
      }
    }
    
    ///------CHI2 W BOSON DELTA R W 1 B JET-------///
    labelsRecoChi2WDeltaRW1B[0] = labelsRecoChi2W[0];
    labelsRecoChi2WDeltaRW1B[1] = labelsRecoChi2W[1];
    
    minDeltaR = 9999.;
    WCandidate = selectedJets[labelsRecoChi2WDeltaRW1B[0]] + selectedJets[labelsRecoChi2WDeltaRW1B[1]];
    for (int kjet = 0; kjet < selectedBJets.size(); kjet++)
    {
      if ( bJetId[kjet] == labelsRecoChi2WDeltaRW1B[0] || bJetId[kjet] == labelsRecoChi2WDeltaRW1B[1] ) continue;
      
      deltaR = ROOT::Math::VectorUtil::DeltaR(selectedJets[bJetId[kjet]], WCandidate);

      if (deltaR < minDeltaR)
      {
        minDeltaR = deltaR;
        labelsRecoChi2WDeltaRW1B[2] = bJetId[kjet];
      }
    }
    
    
    
    ///////////////////////////////////
    ///  CHECK MATCHED COMBINATION  ///
    ///////////////////////////////////

    if (all4PartonsMatched)
    {
      if ( labelsRecoDeltaRAll[0] == MCPermutation[0].first 
          && ( (labelsRecoDeltaRAll[1] == MCPermutation[1].first && labelsRecoDeltaRAll[2] == MCPermutation[2].first) 
            || (labelsRecoDeltaRAll[1] == MCPermutation[2].first && labelsRecoDeltaRAll[2] == MCPermutation[1].first) ) )
      {
        nofCorrectlyMatched_deltaRAll[0]++;
        if ( selectedJets.size() == 4 ) nofCorrectlyMatched_deltaRAll[1]++;
        if ( selectedJets.size() == 5 ) nofCorrectlyMatched_deltaRAll[2]++;
        if ( selectedJets.size() == 6 ) nofCorrectlyMatched_deltaRAll[3]++;
        if ( selectedJets.size() > 6 )  nofCorrectlyMatched_deltaRAll[4]++;
      }
      else
      {
        nofNotCorrectlyMatched_deltaRAll[0]++;
        if ( selectedJets.size() == 4 ) nofNotCorrectlyMatched_deltaRAll[1]++;
        if ( selectedJets.size() == 5 ) nofNotCorrectlyMatched_deltaRAll[2]++;
        if ( selectedJets.size() == 6 ) nofNotCorrectlyMatched_deltaRAll[3]++;
        if ( selectedJets.size() > 6 )  nofNotCorrectlyMatched_deltaRAll[4]++;
      }
      
      if ( labelsRecoDeltaR1B[0] == MCPermutation[0].first 
          && ( (labelsRecoDeltaR1B[1] == MCPermutation[1].first && labelsRecoDeltaR1B[2] == MCPermutation[2].first) 
            || (labelsRecoDeltaR1B[1] == MCPermutation[2].first && labelsRecoDeltaR1B[2] == MCPermutation[1].first) ) )
      {
        nofCorrectlyMatched_deltaR1B[0]++;
        if ( selectedJets.size() == 4 ) nofCorrectlyMatched_deltaR1B[1]++;
        if ( selectedJets.size() == 5 ) nofCorrectlyMatched_deltaR1B[2]++;
        if ( selectedJets.size() == 6 ) nofCorrectlyMatched_deltaR1B[3]++;
        if ( selectedJets.size() > 6 )  nofCorrectlyMatched_deltaR1B[4]++;
      }
      else
      {
        nofNotCorrectlyMatched_deltaR1B[0]++;
        if ( selectedJets.size() == 4 ) nofNotCorrectlyMatched_deltaR1B[1]++;
        if ( selectedJets.size() == 5 ) nofNotCorrectlyMatched_deltaR1B[2]++;
        if ( selectedJets.size() == 6 ) nofNotCorrectlyMatched_deltaR1B[3]++;
        if ( selectedJets.size() > 6 )  nofNotCorrectlyMatched_deltaR1B[4]++;
      }
      
      if ( labelsRecoChi2All[0] == MCPermutation[0].first 
          && ( (labelsRecoChi2All[1] == MCPermutation[1].first && labelsRecoChi2All[2] == MCPermutation[2].first) 
            || (labelsRecoChi2All[1] == MCPermutation[2].first && labelsRecoChi2All[2] == MCPermutation[1].first) ) )
      {
        nofCorrectlyMatched_chi2All[0]++;
        if ( selectedJets.size() == 4 ) nofCorrectlyMatched_chi2All[1]++;
        if ( selectedJets.size() == 5 ) nofCorrectlyMatched_chi2All[2]++;
        if ( selectedJets.size() == 6 ) nofCorrectlyMatched_chi2All[3]++;
        if ( selectedJets.size() > 6 )  nofCorrectlyMatched_chi2All[4]++;
      }
      else
      {
        nofNotCorrectlyMatched_chi2All[0]++;
        if ( selectedJets.size() == 4 ) nofNotCorrectlyMatched_chi2All[1]++;
        if ( selectedJets.size() == 5 ) nofNotCorrectlyMatched_chi2All[2]++;
        if ( selectedJets.size() == 6 ) nofNotCorrectlyMatched_chi2All[3]++;
        if ( selectedJets.size() > 6 )  nofNotCorrectlyMatched_chi2All[4]++;
      }
      
      if ( labelsRecoChi2W[0] == MCPermutation[0].first 
          && ( (labelsRecoChi2W[1] == MCPermutation[1].first && labelsRecoChi2W[2] == MCPermutation[2].first) 
            || (labelsRecoChi2W[1] == MCPermutation[2].first && labelsRecoChi2W[2] == MCPermutation[1].first) ) )
      {
        nofCorrectlyMatched_chi2W[0]++;
        if ( selectedJets.size() == 4 )
        {
          nofCorrectlyMatched_chi2W[1]++;
          nofCorrectlyMatched_comb++;
        }
        if ( selectedJets.size() == 5 ) nofCorrectlyMatched_chi2W[2]++;
        if ( selectedJets.size() == 6 ) nofCorrectlyMatched_chi2W[3]++;
        if ( selectedJets.size() > 6 )  nofCorrectlyMatched_chi2W[4]++;
      }
      else
      {
        nofNotCorrectlyMatched_chi2W[0]++;
        if ( selectedJets.size() == 4 )
        {
          nofNotCorrectlyMatched_chi2W[1]++;
          nofNotCorrectlyMatched_comb++;
        }
        if ( selectedJets.size() == 5 ) nofNotCorrectlyMatched_chi2W[2]++;
        if ( selectedJets.size() == 6 ) nofNotCorrectlyMatched_chi2W[3]++;
        if ( selectedJets.size() > 6 )  nofNotCorrectlyMatched_chi2W[4]++;
      }
      
      if ( labelsRecoChi2W1B[0] == MCPermutation[0].first 
          && ( (labelsRecoChi2W1B[1] == MCPermutation[1].first && labelsRecoChi2W1B[2] == MCPermutation[2].first) 
            || (labelsRecoChi2W1B[1] == MCPermutation[2].first && labelsRecoChi2W1B[2] == MCPermutation[1].first) ) )
      {
        nofCorrectlyMatched_chi2W1B[0]++;
        if ( selectedJets.size() > 4 )  nofCorrectlyMatched_comb++;
        if ( selectedJets.size() == 4 ) nofCorrectlyMatched_chi2W1B[1]++;
        if ( selectedJets.size() == 5 ) nofCorrectlyMatched_chi2W1B[2]++;
        if ( selectedJets.size() == 6 ) nofCorrectlyMatched_chi2W1B[3]++;
        if ( selectedJets.size() > 6 )  nofCorrectlyMatched_chi2W1B[4]++;
      }
      else
      {
        nofNotCorrectlyMatched_chi2W1B[0]++;
        if ( selectedJets.size() > 4 )  nofNotCorrectlyMatched_comb++;
        if ( selectedJets.size() == 4 ) nofNotCorrectlyMatched_chi2W1B[1]++;
        if ( selectedJets.size() == 5 ) nofNotCorrectlyMatched_chi2W1B[2]++;
        if ( selectedJets.size() == 6 ) nofNotCorrectlyMatched_chi2W1B[3]++;
        if ( selectedJets.size() > 6 )  nofNotCorrectlyMatched_chi2W1B[4]++;
      }
      
      if ( labelsRecoChi2WDeltaRW[0] == MCPermutation[0].first 
          && ( (labelsRecoChi2WDeltaRW[1] == MCPermutation[1].first && labelsRecoChi2WDeltaRW[2] == MCPermutation[2].first) 
            || (labelsRecoChi2WDeltaRW[1] == MCPermutation[2].first && labelsRecoChi2WDeltaRW[2] == MCPermutation[1].first) ) )
      {
        nofCorrectlyMatched_chi2WDeltaRW[0]++;
        if ( selectedJets.size() == 4 )
        {
          nofCorrectlyMatched_chi2WDeltaRW[1]++;
          nofCorrectlyMatched_comb++;
        }
        if ( selectedJets.size() == 5 ) nofCorrectlyMatched_chi2WDeltaRW[2]++;
        if ( selectedJets.size() == 6 ) nofCorrectlyMatched_chi2WDeltaRW[3]++;
        if ( selectedJets.size() > 6 )  nofCorrectlyMatched_chi2WDeltaRW[4]++;
      }
      else
      {
        nofNotCorrectlyMatched_chi2WDeltaRW[0]++;
        if ( selectedJets.size() == 4 )
        {
          nofNotCorrectlyMatched_chi2WDeltaRW[1]++;
          nofNotCorrectlyMatched_comb++;
        }
        if ( selectedJets.size() == 5 ) nofNotCorrectlyMatched_chi2WDeltaRW[2]++;
        if ( selectedJets.size() == 6 ) nofNotCorrectlyMatched_chi2WDeltaRW[3]++;
        if ( selectedJets.size() > 6 )  nofNotCorrectlyMatched_chi2WDeltaRW[4]++;
      }
      
      if ( labelsRecoChi2WDeltaRW1B[0] == MCPermutation[0].first 
          && ( (labelsRecoChi2WDeltaRW1B[1] == MCPermutation[1].first && labelsRecoChi2WDeltaRW1B[2] == MCPermutation[2].first) 
            || (labelsRecoChi2WDeltaRW1B[1] == MCPermutation[2].first && labelsRecoChi2WDeltaRW1B[2] == MCPermutation[1].first) ) )
      {
        nofCorrectlyMatched_chi2WDeltaRW1B[0]++;
        if ( selectedJets.size() > 4 )  nofCorrectlyMatched_comb++;
        if ( selectedJets.size() == 4 ) nofCorrectlyMatched_chi2WDeltaRW1B[1]++;
        if ( selectedJets.size() == 5 ) nofCorrectlyMatched_chi2WDeltaRW1B[2]++;
        if ( selectedJets.size() == 6 ) nofCorrectlyMatched_chi2WDeltaRW1B[3]++;
        if ( selectedJets.size() > 6 )  nofCorrectlyMatched_chi2WDeltaRW1B[4]++;
      }
      else
      {
        nofNotCorrectlyMatched_chi2WDeltaRW1B[0]++;
        if ( selectedJets.size() > 4 )  nofNotCorrectlyMatched_comb++;
        if ( selectedJets.size() == 4 ) nofNotCorrectlyMatched_chi2WDeltaRW1B[1]++;
        if ( selectedJets.size() == 5 ) nofNotCorrectlyMatched_chi2WDeltaRW1B[2]++;
        if ( selectedJets.size() == 6 ) nofNotCorrectlyMatched_chi2WDeltaRW1B[3]++;
        if ( selectedJets.size() > 6 )  nofNotCorrectlyMatched_chi2WDeltaRW1B[4]++;
      }
    }  // end all4PartonsMatched


  }  // end loop events


  cout << endl;
  cout << "Number of matched events: " << nofMatchedEvents << endl << endl;
  
  
  float percentage_deltaRAll[nofCats] = {0.}, percentage_deltaR1B[nofCats] = {0.}, percentage_chi2All[nofCats] = {0.}, percentage_chi2W[nofCats] = {0.}, percentage_chi2W1B[nofCats] = {0.}, percentage_chi2WDeltaRW[nofCats] = {0.}, percentage_chi2WDeltaRW1B[nofCats] = {0.};
  for (int iCat = 0; iCat < nofCats; iCat++)
  {
    if ( nofCorrectlyMatched_deltaRAll[iCat] != 0 || nofNotCorrectlyMatched_deltaRAll[iCat] != 0 )
    {
      percentage_deltaRAll[iCat] = 100*(float)nofCorrectlyMatched_deltaRAll[iCat] / (float)(nofCorrectlyMatched_deltaRAll[iCat] + nofNotCorrectlyMatched_deltaRAll[iCat]);
    }
    
    if ( nofCorrectlyMatched_deltaR1B[iCat] != 0 || nofNotCorrectlyMatched_deltaR1B[iCat] != 0 )
    {
      percentage_deltaR1B[iCat] = 100*(float)nofCorrectlyMatched_deltaR1B[iCat] / (float)(nofCorrectlyMatched_deltaR1B[iCat] + nofNotCorrectlyMatched_deltaR1B[iCat]);
    }
    
    if ( nofCorrectlyMatched_chi2All[iCat] != 0 || nofNotCorrectlyMatched_chi2All[iCat] != 0 )
    {
      percentage_chi2All[iCat] = 100*(float)nofCorrectlyMatched_chi2All[iCat] / (float)(nofCorrectlyMatched_chi2All[iCat] + nofNotCorrectlyMatched_chi2All[iCat]);
    }
    
    if ( nofCorrectlyMatched_chi2W[iCat] != 0 || nofNotCorrectlyMatched_chi2W[iCat] != 0 )
    {
      percentage_chi2W[iCat] = 100*(float)nofCorrectlyMatched_chi2W[iCat] / (float)(nofCorrectlyMatched_chi2W[iCat] + nofNotCorrectlyMatched_chi2W[iCat]);
    }
    
    if ( nofCorrectlyMatched_chi2W1B[iCat] != 0 || nofNotCorrectlyMatched_chi2W1B[iCat] != 0 )
    {
      percentage_chi2W1B[iCat] = 100*(float)nofCorrectlyMatched_chi2W1B[iCat] / (float)(nofCorrectlyMatched_chi2W1B[iCat] + nofNotCorrectlyMatched_chi2W1B[iCat]);
    }
    
    if ( nofCorrectlyMatched_chi2WDeltaRW[iCat] != 0 || nofNotCorrectlyMatched_chi2WDeltaRW[iCat] != 0 )
    {
      percentage_chi2WDeltaRW[iCat] = 100*(float)nofCorrectlyMatched_chi2WDeltaRW[iCat] / (float)(nofCorrectlyMatched_chi2WDeltaRW[iCat] + nofNotCorrectlyMatched_chi2WDeltaRW[iCat]);
    }
    
    if ( nofCorrectlyMatched_chi2WDeltaRW1B[iCat] != 0 || nofNotCorrectlyMatched_chi2WDeltaRW1B[iCat] != 0 )
    {
      percentage_chi2WDeltaRW1B[iCat] = 100*(float)nofCorrectlyMatched_chi2WDeltaRW1B[iCat] / (float)(nofCorrectlyMatched_chi2WDeltaRW1B[iCat] + nofNotCorrectlyMatched_chi2WDeltaRW1B[iCat]);
    }
  }  // end for nofCats
  
  
  cout << "///------DELTA R ALL JETS------///" << endl;
  cout << "Correctly matched for deltaRAll:     " << setw(8) << right << nofCorrectlyMatched_deltaRAll[0] << endl;
  cout << "Not correctly matched for deltaRAll: " << setw(8) << right << nofNotCorrectlyMatched_deltaRAll[0] << endl;
  if ( percentage_deltaRAll[0] != 0. )
    cout << "   ===> This means that " << percentage_deltaRAll[0] << "% is correctly matched." << endl << endl;
  cout << "///------DELTA R 1 B JET------///" << endl;
  cout << "Correctly matched for deltaR1B:     " << setw(8) << right << nofCorrectlyMatched_deltaR1B[0] << endl;
  cout << "Not correctly matched for deltaR1B: " << setw(8) << right << nofNotCorrectlyMatched_deltaR1B[0] << endl;
  if ( percentage_deltaR1B[0] != 0. )
    cout << "   ===> This means that " << percentage_deltaR1B[0] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 ALL JETS---------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Used top mass: " << chi2TopMass << "; sigma: " << sigmaChi2TopMass << endl;
  cout << "Correctly matched for chi2All:     " << setw(8) << right << nofCorrectlyMatched_chi2All[0] << endl;
  cout << "Not correctly matched for chi2All: " << setw(8) << right << nofNotCorrectlyMatched_chi2All[0] << endl;
  if ( percentage_chi2All[0] != 0. )
    cout << "   ===> This means that " << percentage_chi2All[0] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON----------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2W:     " << setw(8) << right << nofCorrectlyMatched_chi2W[0] << endl;
  cout << "Not correctly matched for chi2W: " << setw(8) << right << nofNotCorrectlyMatched_chi2W[0] << endl;
  if ( percentage_chi2W[0] != 0. )
    cout << "   ===> This means that " << percentage_chi2W[0] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON 1B-------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2W1B:     " << setw(8) << right << nofCorrectlyMatched_chi2W1B[0] << endl;
  cout << "Not correctly matched for chi2W1B: " << setw(8) << right << nofNotCorrectlyMatched_chi2W1B[0] << endl;
  if ( percentage_chi2W1B[0] != 0. )
    cout << "   ===> This means that " << percentage_chi2W1B[0] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON DELTA R W----------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2WDeltaRW:     " << setw(8) << right << nofCorrectlyMatched_chi2WDeltaRW[0] << endl;
  cout << "Not correctly matched for chi2WDeltaRW: " << setw(8) << right << nofNotCorrectlyMatched_chi2WDeltaRW[0] << endl;
  if ( percentage_chi2WDeltaRW[0] != 0. )
    cout << "   ===> This means that " << percentage_chi2WDeltaRW[0] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON DELTA R W 1B-------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2WDeltaRW1B:     " << setw(8) << right << nofCorrectlyMatched_chi2WDeltaRW1B[0] << endl;
  cout << "Not correctly matched for chi2WDeltaRW1B: " << setw(8) << right << nofNotCorrectlyMatched_chi2WDeltaRW1B[0] << endl;
  if ( percentage_chi2WDeltaRW1B[0] != 0. )
    cout << "   ===> This means that " << percentage_chi2WDeltaRW1B[0] << "% is correctly matched." << endl << endl;
  
  cout << "---------EVENTS WITH EXACTLY 4 JETS---------" << endl;
  cout << "///------DELTA R ALL JETS------///" << endl;
  cout << "Correctly matched for deltaRAll:     " << setw(8) << right << nofCorrectlyMatched_deltaRAll[1] << endl;
  cout << "Not correctly matched for deltaRAll: " << setw(8) << right << nofNotCorrectlyMatched_deltaRAll[1] << endl;
  if ( percentage_deltaRAll[1] != 0. )
    cout << "   ===> This means that " << percentage_deltaRAll[1] << "% is correctly matched." << endl << endl;
  cout << "///------DELTA R 1 B JET------///" << endl;
  cout << "Correctly matched for deltaR1B:     " << setw(8) << right << nofCorrectlyMatched_deltaR1B[1] << endl;
  cout << "Not correctly matched for deltaR1B: " << setw(8) << right << nofNotCorrectlyMatched_deltaR1B[1] << endl;
  if ( percentage_deltaR1B[1] != 0. )
    cout << "   ===> This means that " << percentage_deltaR1B[1] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 ALL JETS---------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Used top mass: " << chi2TopMass << "; sigma: " << sigmaChi2TopMass << endl;
  cout << "Correctly matched for chi2All:     " << setw(8) << right << nofCorrectlyMatched_chi2All[1] << endl;
  cout << "Not correctly matched for chi2All: " << setw(8) << right << nofNotCorrectlyMatched_chi2All[1] << endl;
  if ( percentage_chi2All[1] != 0. )
    cout << "   ===> This means that " << percentage_chi2All[1] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON----------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2W:     " << setw(8) << right << nofCorrectlyMatched_chi2W[1] << endl;
  cout << "Not correctly matched for chi2W: " << setw(8) << right << nofNotCorrectlyMatched_chi2W[1] << endl;
  if ( percentage_chi2W[1] != 0. )
    cout << "   ===> This means that " << percentage_chi2W[1] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON 1B-------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2W1B:     " << setw(8) << right << nofCorrectlyMatched_chi2W1B[1] << endl;
  cout << "Not correctly matched for chi2W1B: " << setw(8) << right << nofNotCorrectlyMatched_chi2W1B[1] << endl;
  if ( percentage_chi2W1B[1] != 0. )
    cout << "   ===> This means that " << percentage_chi2W1B[1] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON DELTA R W----------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2WDeltaRW:     " << setw(8) << right << nofCorrectlyMatched_chi2WDeltaRW[1] << endl;
  cout << "Not correctly matched for chi2WDeltaRW: " << setw(8) << right << nofNotCorrectlyMatched_chi2WDeltaRW[1] << endl;
  if ( percentage_chi2WDeltaRW[1] != 0. )
    cout << "   ===> This means that " << percentage_chi2WDeltaRW[1] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON DELTA R W 1B-------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2WDeltaRW1B:     " << setw(8) << right << nofCorrectlyMatched_chi2WDeltaRW1B[1] << endl;
  cout << "Not correctly matched for chi2WDeltaRW1B: " << setw(8) << right << nofNotCorrectlyMatched_chi2WDeltaRW1B[1] << endl;
  if ( percentage_chi2WDeltaRW1B[1] != 0. )
    cout << "   ===> This means that " << percentage_chi2WDeltaRW1B[1] << "% is correctly matched." << endl << endl;
  
  cout << "---------EVENTS WITH EXACTLY 5 JETS---------" << endl;
  cout << "///------DELTA R ALL JETS------///" << endl;
  cout << "Correctly matched for deltaRAll:     " << setw(8) << right << nofCorrectlyMatched_deltaRAll[2] << endl;
  cout << "Not correctly matched for deltaRAll: " << setw(8) << right << nofNotCorrectlyMatched_deltaRAll[2] << endl;
  if ( percentage_deltaRAll[2] != 0. )
    cout << "   ===> This means that " << percentage_deltaRAll[2] << "% is correctly matched." << endl << endl;
  cout << "///------DELTA R 1 B JET------///" << endl;
  cout << "Correctly matched for deltaR1B:     " << setw(8) << right << nofCorrectlyMatched_deltaR1B[2] << endl;
  cout << "Not correctly matched for deltaR1B: " << setw(8) << right << nofNotCorrectlyMatched_deltaR1B[2] << endl;
  if ( percentage_deltaR1B[2] != 0. )
    cout << "   ===> This means that " << percentage_deltaR1B[2] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 ALL JETS---------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Used top mass: " << chi2TopMass << "; sigma: " << sigmaChi2TopMass << endl;
  cout << "Correctly matched for chi2All:     " << setw(8) << right << nofCorrectlyMatched_chi2All[2] << endl;
  cout << "Not correctly matched for chi2All: " << setw(8) << right << nofNotCorrectlyMatched_chi2All[2] << endl;
  if ( percentage_chi2All[2] != 0. )
    cout << "   ===> This means that " << percentage_chi2All[2] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON----------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2W:     " << setw(8) << right << nofCorrectlyMatched_chi2W[2] << endl;
  cout << "Not correctly matched for chi2W: " << setw(8) << right << nofNotCorrectlyMatched_chi2W[2] << endl;
  if ( percentage_chi2W[2] != 0. )
    cout << "   ===> This means that " << percentage_chi2W[2] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON 1B-------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2W1B:     " << setw(8) << right << nofCorrectlyMatched_chi2W1B[2] << endl;
  cout << "Not correctly matched for chi2W1B: " << setw(8) << right << nofNotCorrectlyMatched_chi2W1B[2] << endl;
  if ( percentage_chi2W1B[2] != 0. )
    cout << "   ===> This means that " << percentage_chi2W1B[2] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON DELTA R W----------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2WDeltaRW:     " << setw(8) << right << nofCorrectlyMatched_chi2WDeltaRW[2] << endl;
  cout << "Not correctly matched for chi2WDeltaRW: " << setw(8) << right << nofNotCorrectlyMatched_chi2WDeltaRW[2] << endl;
  if ( percentage_chi2WDeltaRW[2] != 0. )
    cout << "   ===> This means that " << percentage_chi2WDeltaRW[2] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON DELTA R W 1B-------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2WDeltaRW1B:     " << setw(8) << right << nofCorrectlyMatched_chi2WDeltaRW1B[2] << endl;
  cout << "Not correctly matched for chi2WDeltaRW1B: " << setw(8) << right << nofNotCorrectlyMatched_chi2WDeltaRW1B[2] << endl;
  if ( percentage_chi2WDeltaRW1B[2] != 0. )
    cout << "   ===> This means that " << percentage_chi2WDeltaRW1B[2] << "% is correctly matched." << endl << endl;
  
  cout << "---------EVENTS WITH EXACTLY 6 JETS---------" << endl;
  cout << "///------DELTA R ALL JETS------///" << endl;
  cout << "Correctly matched for deltaRAll:     " << setw(8) << right << nofCorrectlyMatched_deltaRAll[3] << endl;
  cout << "Not correctly matched for deltaRAll: " << setw(8) << right << nofNotCorrectlyMatched_deltaRAll[3] << endl;
  if ( percentage_deltaRAll[3] != 0. )
    cout << "   ===> This means that " << percentage_deltaRAll[3] << "% is correctly matched." << endl << endl;
  cout << "///------DELTA R 1 B JET------///" << endl;
  cout << "Correctly matched for deltaR1B:     " << setw(8) << right << nofCorrectlyMatched_deltaR1B[3] << endl;
  cout << "Not correctly matched for deltaR1B: " << setw(8) << right << nofNotCorrectlyMatched_deltaR1B[3] << endl;
  if ( percentage_deltaR1B[3] != 0. )
    cout << "   ===> This means that " << percentage_deltaR1B[3] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 ALL JETS---------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Used top mass: " << chi2TopMass << "; sigma: " << sigmaChi2TopMass << endl;
  cout << "Correctly matched for chi2All:     " << setw(8) << right << nofCorrectlyMatched_chi2All[3] << endl;
  cout << "Not correctly matched for chi2All: " << setw(8) << right << nofNotCorrectlyMatched_chi2All[3] << endl;
  if ( percentage_chi2All[3] != 0. )
    cout << "   ===> This means that " << percentage_chi2All[3] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON----------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2W:     " << setw(8) << right << nofCorrectlyMatched_chi2W[3] << endl;
  cout << "Not correctly matched for chi2W: " << setw(8) << right << nofNotCorrectlyMatched_chi2W[3] << endl;
  if ( percentage_chi2W[3] != 0. )
    cout << "   ===> This means that " << percentage_chi2W[3] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON 1B-------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2W1B:     " << setw(8) << right << nofCorrectlyMatched_chi2W1B[3] << endl;
  cout << "Not correctly matched for chi2W1B: " << setw(8) << right << nofNotCorrectlyMatched_chi2W1B[3] << endl;
  if ( percentage_chi2W1B[3] != 0. )
    cout << "   ===> This means that " << percentage_chi2W1B[3] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON DELTA R W----------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2WDeltaRW:     " << setw(8) << right << nofCorrectlyMatched_chi2WDeltaRW[3] << endl;
  cout << "Not correctly matched for chi2WDeltaRW: " << setw(8) << right << nofNotCorrectlyMatched_chi2WDeltaRW[3] << endl;
  if ( percentage_chi2WDeltaRW[3] != 0. )
    cout << "   ===> This means that " << percentage_chi2WDeltaRW[3] << "% is correctly matched." << endl << endl;
  cout << "///------CHI2 W BOSON DELTA R W 1B-------///" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched for chi2WDeltaRW1B:     " << setw(8) << right << nofCorrectlyMatched_chi2WDeltaRW1B[3] << endl;
  cout << "Not correctly matched for chi2WDeltaRW1B: " << setw(8) << right << nofNotCorrectlyMatched_chi2WDeltaRW1B[3] << endl;
  if ( percentage_chi2WDeltaRW1B[3] != 0. )
    cout << "   ===> This means that " << percentage_chi2WDeltaRW1B[3] << "% is correctly matched." << endl << endl;
  
  
  cout << "---------(Attempt at) Best Combination---------" << endl;
  cout << "Use Chi2W when exactly 4 jets && Chi2W1B when more than 4 jets" << endl;
  cout << "Used W mass:   " << chi2WMass << "; sigma: " << sigmaChi2WMass << endl;
  cout << "Correctly matched:     " << setw(8) << right << nofCorrectlyMatched_comb << endl;
  cout << "Not correctly matched: " << setw(8) << right << nofNotCorrectlyMatched_comb << endl;
  if ( nofCorrectlyMatched_comb != 0 || nofNotCorrectlyMatched_comb != 0 )
    cout << "   ===> This means that " << 100*(float)nofCorrectlyMatched_comb / (float)(nofCorrectlyMatched_comb + nofNotCorrectlyMatched_comb) << "% is correctly matched." << endl << endl;
  
  
  cout << endl;
  cout << "              All jets      4jets         5jets         6jets         6+ jets" << endl;
  cout << setw(6) << right << fixed << showpoint << setprecision(3);
//  cout << "test          33.333%       33.333%       33.333%       33.333%       33.333%" << endl;
  cout << "deltaRAll        " << percentage_deltaRAll[0] << "%       "  << percentage_deltaRAll[1] << "%       "  << percentage_deltaRAll[2] << "%       "  << percentage_deltaRAll[3] << "%       "  << percentage_deltaRAll[4] << "%       "  << endl;
  cout << "deltaR1B         " << percentage_deltaR1B[0] << "%       "  << percentage_deltaR1B[1] << "%       "  << percentage_deltaR1B[2] << "%       "  << percentage_deltaR1B[3] << "%       "  << percentage_deltaR1B[4] << "%       "  << endl;
  cout << "chi2All          " << percentage_chi2All[0] << "%       "  << percentage_chi2All[1] << "%       "  << percentage_chi2All[2] << "%       "  << percentage_chi2All[3] << "%       "  << percentage_chi2All[4] << "%       "  << endl;
  cout << "chi2W            " << percentage_chi2W[0] << "%       "  << percentage_chi2W[1] << "%       "  << percentage_chi2W[2] << "%       "  << percentage_chi2W[3] << "%       "  << percentage_chi2W[4] << "%       "  << endl;
  cout << "chi2W1B          " << percentage_chi2W1B[0] << "%       "  << percentage_chi2W1B[1] << "%       "  << percentage_chi2W1B[2] << "%       "  << percentage_chi2W1B[3] << "%       "  << percentage_chi2W1B[4] << "%       "  << endl;
  cout << "chi2WDeltaRW     " << percentage_chi2WDeltaRW[0] << "%       "  << percentage_chi2WDeltaRW[1] << "%       "  << percentage_chi2WDeltaRW[2] << "%       "  << percentage_chi2WDeltaRW[3] << "%       "  << percentage_chi2WDeltaRW[4] << "%       "  << endl;
  cout << "chi2WDeltaRW1B   " << percentage_chi2WDeltaRW1B[0] << "%       "  << percentage_chi2WDeltaRW1B[1] << "%       "  << percentage_chi2WDeltaRW1B[2] << "%       "  << percentage_chi2WDeltaRW1B[3] << "%       "  << percentage_chi2WDeltaRW1B[4] << "%       "  << endl;
  cout << endl;
  
  
  origNtuple->Close();

  
  if (test)
  {
    cout << "Exiting because of test..." << endl;
    exit(1);
  }
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
//   string rootFileName = "NtuplePlots_"+systStr+".root";
//   
//   cout << " - Recreate output file ..." << endl;
//   TFile *fout = new TFile ((pathOutput+rootFileName).c_str(), "RECREATE");
//   cout << "   Output file is " << pathOutput+rootFileName << endl;
//   
//   ///Write histograms
//   fout->cd();
//   
//   // 1D
//   TDirectory* th1dir = fout->mkdir("1D_histograms");
//   th1dir->cd();
//   gStyle->SetOptStat(1111);
//   for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
//   {
//     TH1F *temp = it->second;
//     int N = temp->GetNbinsX();
//     temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
//     temp->SetBinContent(N+1,0);
//     temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
//     temp->Write();
//     TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
//     tempCanvas->SaveAs( (pathOutput+it->first+".png").c_str() );
//   }
//   
//   // 2D
//   TDirectory* th2dir = fout->mkdir("2D_histograms");
//   th2dir->cd();
//   gStyle->SetPalette(55);
//   for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
//   {
//     TH2F *temp = it->second;
//     temp->Write();
//     TCanvas* tempCanvas = TCanvasCreator(temp, it->first, "colz");
//     tempCanvas->SaveAs( (pathOutput+it->first+".png").c_str() );
//   }
//   
//   fout->Close();
//   
//   delete fout;
  
  
  
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
}


void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
//   /// Chi2
//   histo1D["Chi2_W_mass_reco"] = new TH1F("Chi2_W_mass_reco","Reconstructed hadronic W mass using a Chi2; M_{W} [GeV]", 80, 0, 800);
//   histo1D["Chi2_top_mass_reco"] = new TH1F("Chi2_top_mass_reco","Reconstructed top mass using a Chi2; M_{t} [GeV]", 80, 0, 800);
//   histo1D["Chi2_top_mass_reco_4jets"] = new TH1F("Chi2_top_mass_reco_4jets","Reconstructed top mass using a Chi2 (exactly 4 jets); M_{t} [GeV]", 80, 0, 800);
//   
//   /// Matching
//   histo1D["W_mass_reco_matched"] = new TH1F("W_mass_reco_matched","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 125, 0, 250);
//   histo1D["top_mass_reco_matched"] = new TH1F("top_mass_reco_matched","Reconstructed top mass of matched events; M_{t} [GeV]", 175, 50, 400);
//   histo1D["top_mass_gen_matched"] = new TH1F("top_mass_gen_matched","Generated top mass of matched events; M_{t} [GeV]", 175, 0, 350);
//   histo1D["mlb_matched_corr"]  = new TH1F("mlb_matched_corr","Reconstructed leptonic top mass using correctly matched events; M_{lb} [GeV]", 80, 0, 800);
//   histo1D["mlb_matched_wrong"] = new TH1F("mlb_matched_wrong","Reconstructed leptonic top mass using wrongly matched events; M_{lb} [GeV]", 80, 0, 800);
//   histo1D["ttbar_mass_matched_corr"] = new TH1F("ttbar_mass_matched_corr","Reconstructed mass of the top quark pair using correctly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
//   histo1D["ttbar_mass_matched_wrong"] = new TH1F("ttbar_mass_matched_wrong","Reconstructed mass of the top quark pair using wrongly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
//   histo1D["dR_lep_b_matched_corr"] = new TH1F("dR_lep_b_matched_corr","Delta R between the lepton and the leptonic b jet for matched events; #Delta R(l,b)", 25, 0, 5);
//   histo1D["dR_lep_b_matched_wrong"] = new TH1F("dR_lep_b_matched_wrong","Delta R between the lepton and the hadronic b jet for matched events; #Delta R(l,b_{had})", 25, 0, 5);
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
  
  nofMatchedEvents = 0;
  for (int i = 0; i < nofCats; i++)
  {
    nofCorrectlyMatched_deltaRAll[i] = 0;
    nofNotCorrectlyMatched_deltaRAll[i] = 0;
    nofCorrectlyMatched_deltaR1B[i] = 0;
    nofNotCorrectlyMatched_deltaR1B[i] = 0;
    nofCorrectlyMatched_chi2All[i] = 0;
    nofNotCorrectlyMatched_chi2All[i] = 0;
    nofCorrectlyMatched_chi2W[i] = 0;
    nofNotCorrectlyMatched_chi2W[i] = 0;
    nofCorrectlyMatched_chi2W1B[i] = 0;
    nofNotCorrectlyMatched_chi2W1B[i] = 0;
    nofCorrectlyMatched_chi2WDeltaRW[i] = 0; 
    nofNotCorrectlyMatched_chi2WDeltaRW[i] = 0; 
    nofCorrectlyMatched_chi2WDeltaRW1B[i] = 0;
    nofNotCorrectlyMatched_chi2WDeltaRW1B[i] = 0;
  }
  nofCorrectlyMatched_comb = 0;
  nofNotCorrectlyMatched_comb = 0;
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
  nloWeight = 1.;
  puSF = 1.;
  btagSF = 1.;
  muonIdSF[1] = 1.;
  muonIsoSF[1] = 1.;
  muonTrigSFv2[1] = 1.;
  muonTrigSFv3[1] = 1.;
  
  scaleFactor = 1.;
  bJetId.clear();
  reco_hadWMass = -1.;
  reco_hadTopMass = -1.;
  reco_hadTopPt = -1.;
  reco_minMlb = 9999.;
  reco_ttbarMass = -1.;
  reco_dRLepB = -1.;
  min_Mlb = 9999.;
  dRLepB = -1.;
}

void ClearTLVs()
{
  muon.Clear();
  jet.Clear();
  mcpart.Clear();
  topCandidate.Clear();
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

