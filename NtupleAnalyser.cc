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
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
//#include "../macros/Style.C"

// user defined
#include "TopTreeAnalysisBase/MCInformation/interface/TransferFunctions.h"


using namespace std;
using namespace TopTree;


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

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;

map<string,TFile*> tFileMap;
map<string,TFile*> globalTFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

map<string,TNtuple*> ntuple;
map<string,TNtuple*> ntree;
map<string,TNtuple*> otree;


/// Function prototypes
string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();
void InitTree(TTree* tree, bool isData);
void ClearLeaves();
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
  
  string pathOutput = "test/";
  //string pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  // Add channel to output path
  pathOutput += channel+"/";
  mkdir(pathOutput.c_str(),0777);
  // Give timestamp to output path
  pathOutput += dateString+"/";
  mkdir(pathOutput.c_str(),0777);
  
  
  pathNtuples = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate+"/";
  
  
  /// xml file
  string xmlFileName ="config/topWidth.xml";
  
  const char *xmlfile = xmlFileName.c_str();
  
  cout << " - Using config file " << xmlfile << endl;
  
  
  
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets(datasets, xmlfile);
  
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName, filepath, slumi;
  
  int nEntries;
  double scaleFactor, normFactor;
  double fracHLT[2] = {-1};
  GetHLTFraction(fracHLT);
  if ( fracHLT[0] == -1 || fracHLT[1] == -1 )
  {
    cout << "Something went wrong with the fraction calculation for trigger SFs!" << endl;
    exit(1);
  }
  
  
  /// Loop over datasets
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    isData = false;
    
    dataSetName = datasets[d]->Name();
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
    }
    
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
    cout << "                 nEntries: " << nEntries << endl;
    
    
    // Set branch addresses and branch pointers
    InitTree(tTree[dataSetName.c_str()], isData);
    
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    /// Define TLVs
    vector<TLorentzVector*> selectedLepton;
    vector<TLorentzVector*> selectedJets;
    vector<TLorentzVector*> mcParticles;
    
    
    for (int ievt = 0; ievt < nEntries; ievt++)
    {
      selectedLepton.clear();
      selectedJets.clear();
      mcParticles.clear();
      ClearLeaves();
      scaleFactor = 1.;
      normFactor = 1.;
      
      
      if (ievt%1000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
      
      /// Load event
      tTree[(dataSetName).c_str()]->GetEntry(ievt);
      
      // Scale factors
      if (! isData)
      {
        if (applyLeptonSF) { scaleFactor *= muonIdSF[0] * muonIsoSF[0] * (fracHLT[0]*muonTrigSFv2[0] + fracHLT[1]*muonTrigSFv3[0]);}
        if (applyPU) { scaleFactor *= puSF;}
        if (applyBTagSF) { scaleFactor *= btagSF;}
        //if (applyNloSF) { scaleFactor *= nloWeight;}  // additional SF due to number of events with neg weight!!
      }
      
      
      
      
    }  // end loop events
    
    
    
    
    
    tFileMap[dataSetName.c_str()]->Close();
    
  }  // end loop datasets
  
  
  
  
  
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

