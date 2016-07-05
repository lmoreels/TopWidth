///////////////////
///   Ntupler   ///
///////////////////


#include "TStyle.h"
#include <ctime>
#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <map>
#include <array>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TTree.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include "TVectorD.h"
#include <errno.h>
#include "TRandom3.h"
#include "TRandom.h"
#include "TProfile.h"
#include <cstdlib>

//used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagCalibrationStandalone.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

// user defined
#include "TopTreeAnalysisBase/Tools/interface/Trigger.h"
#include "TopTreeAnalysisBase/MCInformation/interface/TransferFunctions.h"


using namespace std;
using namespace reweight;
using namespace TopTree;


map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;


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



int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  pathOutput = "Ntupler_"+dateString+"/";
  
  
  
  /////////////////////
  ///  Configuration
  /////////////////////
  
  bool testTTbarOnly = false;  
  bool calculateTransferFunctions = false;
  bool printTriggers = false;
  bool applyTriggers = true;
  bool applyLeptonSF = true;
  bool applyPU = true;
  bool applyPUup = false;
  bool applyPUdown = false;
  bool applyJER = true;
  bool applyJERup = false;
  bool applyJERdown = false;
  bool applyJEC = true;
  bool applyJESup = false;
  bool applyJESdown = false;
  bool calculateBTagSF = false;
  bool applyBTagSF = true;
  bool applyJetLeptonCleaning = true;
  
  if (  (applyPUup    && (             applyPUdown || applyJERup || applyJERdown || applyJESup || applyJESdown))
     || (applyPUdown  && (applyPUup ||                applyJERup || applyJERdown || applyJESup || applyJESdown))
     || (applyJERup   && (applyPUup || applyPUdown ||               applyJERdown || applyJESup || applyJESdown))
     || (applyJERdown && (applyPUup || applyPUdown || applyJERup ||                 applyJESup || applyJESdown))
     || (applyJESup   && (applyPUup || applyPUdown || applyJERup || applyJERdown ||               applyJESdown))
     || (applyJESdown && (applyPUup || applyPUdown || applyJERup || applyJERdown || applyJESup                )) )
  {
    cerr << "SCALE FACTORS: Cannot scale up/down multiple scale factors at once!" << endl;
    cerr << "  - Stopping the program... " << endl;
    exit(1);
  }
  if (calculateBTagSF && applyBTagSF)
  {
    cerr << "SCALE FACTORS: Cannot calculate & apply b-tag scale factors at the same time!" << endl;
    cerr << "  - Stopping the program... " << endl;
    exit(1);
  }
  
  cout << "* The following scale factors are applied:  *" << endl;
  if (applyTriggers) cout << "*   - Triggers                              *" << endl;
  if (applyLeptonSF) cout << "*   - Lepton scale factors                  *" << endl;
  if (applyPU)       cout << "*   - Pile up                               *" << endl;
  if (applyJER)
  {
    cout << "*   - Jet Energy Resolution: ";
    if (applyJERdown)    cout << "scale down     *" << endl;
    else if (applyJERup) cout << "scale up       *" << endl;
    else                 cout << "nominal        *" << endl;
  }
  if (applyJEC)         cout << "*   - Jet Energy Corrections                *" << endl;
  if (applyJESdown)     cout << "*   - Jet Energy Scale: scale down          *" << endl;
  else if (applyJESup)  cout << "*   - Jet Energy Scale: scale up            *" << endl;
  if (calculateBTagSF)  cout << "*   - Preparing histos for b tag SFs...     *" << endl;
  else if (applyBTagSF) cout << "*   - B tag scale factors                   *" << endl;
  
  if (applyJetLeptonCleaning) cout << "*   - Jet/lepton Cleaning                   *" << endl;
  cout << "*********************************************" << endl;
  
  
  
  bool hasNegWeight = false;
  bool eventSelected = false;
  bool hasExactly4Jets = false;
  bool has1bjet = false;
  bool has2bjets = false;
  int nofSelectedEvents = 0;
  int nofMatchedEvents = 0;
  int nofEventsWith1BJet = 0;
  int nofEventsWith2BJets = 0;
  int nofNegWeights = 0;
  int nofPosWeights = 0;
  int nofEventsHLTv2 = 0;
  int nofEventsHLTv3 = 0;
  int nofEventsJetLeptonCleaned = 0;
  
  /// xml file
  string xmlFileName ="config/topWidth_skimmed.xml";
  
  if (argc > 1)
  {
    xmlFileName = (string)argv[1];
  }
  
  const char *xmlfile = xmlFileName.c_str();
  
  cout << " - Using config file " << xmlfile << endl;
  
  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
  
  
  ////////////////////////////////////
  ///  AnalysisEnvironment
  ////////////////////////////////////
  
  AnalysisEnvironment anaEnv;
  cout << " - Loading environment ..." << endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  
  cout << anaEnv.JetCollection << " " <<  anaEnv.METCollection << " "
    << anaEnv.ElectronCollection << " " << anaEnv.MuonCollection << " "
    << anaEnv.PrimaryVertexCollection << " " << anaEnv.GenJetCollection << " "
    << endl;
  
//   anaEnv.PrimaryVertexCollection = "PrimaryVertex";
//   anaEnv.JetCollection = "PFJets_slimmedJets";
//   anaEnv.METCollection = "PFMET_slimmedMETs";
//   anaEnv.MuonCollection = "Muons_slimmedMuons";
//   anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
//   anaEnv.GenJetCollection = "GenJets_slimmedGenJets";
//   anaEnv.MCParticlesCollection = "MCParticles";
//   anaEnv.loadFatJetCollection = false;
//   anaEnv.loadGenJetCollection = true;
//   anaEnv.loadNPGenEventCollection = false;
//   anaEnv.loadMCParticles = true;
//   anaEnv.JetType = 2;
//   anaEnv.METType = 2;
  
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
//  verbose = 2;
  float oldLuminosity = anaEnv.Luminosity;  // in 1/pb
  cout << "Analysis environment luminosity for rescaling " << oldLuminosity << endl;
  
  
  
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;
  
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets(datasets, xmlfile); cout << "Number of datasets: " << datasets.size() << endl;
  for (unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
//   const string dName, dTitle;
//   const int color, ls, lw;
//   const float normf, EqLumi, xSect;
//   vector<string> vecfileNames;
//   Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
//   theDataset->SetEquivalentLuminosity(EqLumi);
//   datasets.push_back(theDataset);
  
  float Luminosity = oldLuminosity;
  
  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    string dataSetName = datasets[d]->Name();
    if ( (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) && Luminosity > datasets[d]->EquivalentLumi() ) { Luminosity = datasets[d]->EquivalentLumi();}
  }
  
  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
  
  
  //Global variable
  //TRootEvent* event = 0;
  TRootRun *runInfos = new TRootRun();
  
  //nof selected events
  //double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  Double_t nloweight = 0;
  
  
  
  ////////////////////////////
  ///  Initialise trigger  ///
  ////////////////////////////
  
  //Trigger* trigger = new Trigger(hasMuon, hasElectron, trigSingleLep, trigDoubleLep);
  Trigger* trigger = new Trigger(1, 0, 1, 0);
  
  
  
  //////////////////////////////////
  ///  Initialise scale factors  ///
  //////////////////////////////////
  
  string pathCalLept = "../TopTreeAnalysisBase/Calibrations/LeptonSF/";
  string pathCalBTag = "../TopTreeAnalysisBase/Calibrations/BTagging/";
  string pathCalPileup = "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/";
  string pathCalJEC = "../TopTreeAnalysisBase/Calibrations/JECFiles/";
  
  /// Leptons
  cout << " - Loading lepton scale factors ...";
  if (! applyLeptonSF) { cout << "    --- At the moment these are not used in the analysis";}
  cout << endl;
  
  double muonSFID, muonSFIso, muonSFTrig;
  //MuonSFWeight *muonSFWeight_ = new MuonSFWeight(pathCalLept+"Muon_SF_TopEA.root","SF_totErr", true, false, false); // (... , ... , extendRange, debug, print warning)
  MuonSFWeight *muonSFWeightID_T = new MuonSFWeight(pathCalLept+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);
  MuonSFWeight *muonSFWeightID_M = new MuonSFWeight(pathCalLept+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);
//  MuonSFWeight *muonSFWeightID_L = new MuonSFWeight(pathCalLept+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);
//  MuonSFWeight *muonSFWeightID_S = new MuonSFWeight(pathCalLept+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_SoftID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Soft ID
  
  MuonSFWeight *muonSFWeightIso_TT = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
//   MuonSFWeight *muonSFWeightIso_TM = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Tight RelIso, Medium ID
//   MuonSFWeight *muonSFWeightIso_LT = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Loose RelIso, Tight ID
//   MuonSFWeight *muonSFWeightIso_LM = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Loose RelIso, Medium ID
//   MuonSFWeight *muonSFWeightIso_LT = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseRelIso_DEN_LooseID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Loose RelIso, Loose ID
  
  double weightMuonHLTv2, weightMuonHLTv3;
  MuonSFWeight *muonSFWeightTrigHLTv4p2 = new MuonSFWeight(pathCalLept+"SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins/abseta_pt_ratio", true, false, false);
  MuonSFWeight *muonSFWeightTrigHLTv4p3 = new MuonSFWeight(pathCalLept+"SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio", true, false, false);
  
  //ElectronSFWeight *electronSFWeight_ = new ElectronSFWeight(pathCalLept+"Elec_SF_TopEA.root","GlobalSF", true, false, false); // (... , ... , extendRange, debug, print warning)
  
  
  /// B tag
  // documentation at http://mon.iihe.ac.be/%7Esmoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees_v4.pdf
  cout << " - Loading b tag scale factors ...";
  if (! applyBTagSF) { cout << "     --- At the moment these are not used in the analysis";}
  BTagCalibration *bTagCalib = new BTagCalibration("CSVv2", pathCalBTag+"CSVv2_76X_combToMujets.csv"); 
  BTagCalibrationReader *bTagReader_M = new BTagCalibrationReader(bTagCalib, BTagEntry::OP_MEDIUM, "mujets","central");
  BTagWeightTools *bTagHistoTool_M = new BTagWeightTools(bTagReader_M,"PlotsForBTagSFs.root",30., 999., 2.4);
  
  
  /// Pile-up
  cout << " - Loading pile-up scale factors ...";
  if (! applyPU) { cout << "   --- At the moment these are not used in the analysis";}
  cout << endl;
  
  LumiReWeighting LumiWeights(pathCalPileup+"pileup_MC_RunIIFall15DR76-Asympt25ns.root",pathCalPileup+"pileup_2015Data76X_25ns-Run246908-260627Cert.root","pileup","pileup");
  
  
  /// JEC
  vector<JetCorrectorParameters> vCorrParam;
  
  
  
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
  
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    nofSelectedEvents = 0;
    nofEventsWith1BJet = 0;
    nofEventsWith2BJets = 0;
    nofEventsJetLeptonCleaned = 0;
    nofNegWeights = 0;
    nofPosWeights = 0;
    float sumWeights = 0.;
    double nloSF = 1.;
    int iFile = -1;
    string previousFilename = "";
    bool nlo = false;
    bool isData = false;
    
    string dataSetName = datasets[d]->Name();
    
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << "/ title : " << datasets[d]->Title() << endl;
      cout << "      -> Equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
      cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    }
    
    //open files and load
    cout << "Load Dataset" << endl;
    treeLoader.LoadDataset(datasets[d], anaEnv);
    cout << "Load Dataset" << endl;
    
    if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
    {
      isData = true;
    }
    
    if ( (datasets[d]->Title()).find("amc@nlo") == 0 )
    {
      nlo = true;
      cout << "         This is an amc@nlo sample." << endl;
    }
    
    /// book triggers
    if (applyTriggers) { trigger->bookTriggers(isData);}
    
    
    
    ///////////////////////////////////////////
    ///  Initialise Jet Energy Corrections  ///
    ///////////////////////////////////////////
    
    vCorrParam.clear();
    
    if (isData)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
    }
    else
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt");
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); //true means redo also L1
    
    ////////////////////////////
    ///  Create output file  ///
    ////////////////////////////
    
    string rootFileName = "Ntuples_output_"+dataSetName+".root";
    rootFileName = "Ntuples_output_"+dataSetName+"_"+ConvertIntToStr(JobNum,0)+".root";
    
    cout << " - Recreate output file ..." << endl;
    TFile *fout = new TFile ((pathOutput+rootFileName).c_str(), "RECREATE");
    fout->cd();
    
    TTree* myTree = new TTree("tree","tree");
    TTree* basicEvent = new TTree("event","event");
    TTree* globalTree = new TTree("globaltree","globaltree");  // no selection applied
    
    
    
    ////////////////////////////////////
    ///  Define variables for trees  ///
    ////////////////////////////////////
    
    // event related variables
    Int_t run_num;
    Int_t evt_num;
    Int_t lumi_num;
    Int_t nvtx;
    Int_t npu;
    
    Bool_t isTrigged;
    Bool_t isSelected;
//    Bool_t passedMETFilter;
    
    Int_t appliedJER;
    Int_t appliedJES;
    
    Int_t nEvents;
    Int_t nEventsSel;
    Int_t nofPosWeights;
    Int_t nofNegWeights;
    Int_t sumW;
    
    Int_t nofEventsHLTv2;
    Int_t nofEventsHLTv3;
    
    Double_t puSF;
    Double_t btagSF;
    Double_t muonIDSF[10];
    Double_t muonIsoSF[10];
    Double_t muonTrigSFv2[10];
    Double_t muonTrigSFv3[10];
//    Double_t electronSF[10];
    Double_t nloWeight; // for amc@nlo samples
    
    
    Int_t nLeptons;
    
//     /// Variables for electrons
//     Int_t nElectrons;
//     Int_t electron_charge[10];
//     Double_t electron_pt[10];
//     Double_t electron_phi[10];
//     Double_t electron_eta[10];
//     Double_t electron_eta_superCluster[10];
//     Double_t electron_E[10];
//     Double_t electron_M[10];
//     Double_t electron_d0[10];
//     Double_t electron_d0BeamSpot[10];
//     Double_t electron_chargedHadronIso[10];
//     Double_t electron_neutralHadronIso[10];
//     Double_t electron_photonIso[10];
//     Double_t electron_pfIso[10];
//     
//     Double_t electron_sigmaIEtaIEta[10];
//     Double_t electron_deltaEtaIn[10];
//     Double_t electron_deltaPhiIn[10];
//     Double_t electron_hadronicOverEm[10];
//     Int_t electron_missingHits[10];
//     Bool_t electron_passConversion[10];
//     Bool_t electron_isId[10];
//     Bool_t electron_isIso[10];
//     Bool_t electron_isEBEEGap[10];
    
    /// Variables for muons
    Int_t nMuons;
    Int_t muon_charge[10];
    Double_t muon_pt[10];
    Double_t muon_phi[10];
    Double_t muon_eta[10];
    Double_t muon_E[10];
    Double_t muon_M[10];
    Double_t muon_d0[10];
    Double_t muon_d0BeamSpot[10];
    Double_t muon_chargedHadronIso[10];
    Double_t muon_neutralHadronIso[10];
    Double_t muon_photonIso[10];
    Double_t muon_relIso[10];
    Double_t muon_pfIso[10];
    Bool_t muon_isId[10];
    Bool_t muon_isIso[10];
    
    /// Variables for jets
    Int_t nJets;
    Int_t jet_charge[20];
    Double_t jet_pt[20];
    Double_t jet_phi[20];
    Double_t jet_eta[20];
    Double_t jet_E[20];
    Double_t jet_M[20];
    Double_t jet_bdiscr[20];
    
    /// met
    Double_t met_pt;
    Double_t met_phi;
    Double_t met_eta;
    Double_t met_E;
    
    /// mcparticles
    Int_t nMCParticles;
    Int_t mc_status[200];
    Int_t mc_pdgId[200];
    Int_t mc_mother[200];
    Int_t mc_granny[200];
    Double_t mc_pt[200];
    Double_t mc_phi[200];
    Double_t mc_eta[200];
    Double_t mc_E[200];
    Double_t mc_M[200];
    
    
    
    /////////////////////////
    ///  Define branches  ///
    /////////////////////////
    
    basicEvent->Branch("run_num",&run_num,"run_num/I");
    basicEvent->Branch("evt_num",&evt_num,"evt_num/I");
    basicEvent->Branch("lumi_num",&lumi_num,"lumi_num/I");
    basicEvent->Branch("nvtx",&nvtx,"nvtx/I");
    basicEvent->Branch("npu",&npu,"npu/I");
    basicEvent->Branch("nEv" , &nEv, "nEv/I");
    basicEvent->Branch("isTrigged",&isTrigged,"isTrigged/O");
    basicEvent->Branch("isSelected",&isSelected,"isSelected/O");
//    basicEvent->Branch("passedMETFilter", &passedMETFilter,"passedMETFilter/O");
    basicEvent->Branch("appliedJER",&appliedJER,"appliedJER/I");
    basicEvent->Branch("appliedJES", &appliedJES, "appliedJES/I");
    
    globalTree->Branch("run_num",&run_num,"run_num/I");
    globalTree->Branch("evt_num",&evt_num,"evt_num/I");
    globalTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
    globalTree->Branch("nvtx",&nvtx,"nvtx/I");
    globalTree->Branch("npu",&npu,"npu/I");
    globalTree->Branch("nEv" , &nEv, "nEv/I");
    globalTree->Branch("isTrigged",&isTrigged,"isTrigged/O");
    glocalTree->Branch("isSelected",&isSelected,"isSelected/O");
//    globalTree->Branch("passedMETFilter", &passedMETFilter,"passedMETFilter/O");
    globalTree->Branch("appliedJER",&appliedJER,"appliedJER/I");
    globalTree->Branch("appliedJES", &appliedJES, "appliedJES/I");
    if (isData)
    {
      globalTree->Branch("nofEventsHLTv2",&nofEventsHLTv2,"nofEventsHLTv2/I");
      globalTree->Branch("nofEventsHLTv3",&nofEventsHLTv3,"nofEventsHLTv3/I");
    }
    if (nlo)
    {
      globalTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
      globalTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
      globalTree->Branch("sumW", &sumW, "sumW/I");
    }
    
    myTree->Branch("run_num",&run_num,"run_num/I");
    myTree->Branch("evt_num",&evt_num,"evt_num/I");
    myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
    myTree->Branch("nvtx",&nvtx,"nvtx/I");
    myTree->Branch("npu",&npu,"npu/I");
    myTree->Branch("nEvSel" , &nEvSel, "nEvSel/I");
    myTree->Branch("isTrigged",&isTrigged,"isTrigged/O");
    myTree->Branch("isSelected",&isSelected,"isSelected/O");  // normally all events should be selected here
//    myTree->Branch("passedMETFilter", &passedMETFilter,"passedMETFilter/O");
    myTree->Branch("appliedJER",&appliedJER,"appliedJER/I");
    myTree->Branch("appliedJES", &appliedJES, "appliedJES/I");
    if (isData)
    {
      myTree->Branch("nofEventsHLTv2",&nofEventsHLTv2,"nofEventsHLTv2/I");
      myTree->Branch("nofEventsHLTv3",&nofEventsHLTv3,"nofEventsHLTv3/I");
    }
    if (nlo)
    {
      myTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
      myTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
      myTree->Branch("sumW", &sumW, "sumW/I");
    }
    
    
    /// SFs
    globalTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
    globalTree->Branch("puSF",&puSF,"puSF/D");
    globalTree->Branch("btagSF",&btagSF,"btagSF/D");
    globalTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
    globalTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
    globalTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
    globalTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
//    globalTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
    
    myTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
    myTree->Branch("puSF",&puSF,"puSF/D");
    myTree->Branch("btagSF",&btagSF,"btagSF/D");
    myTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
    myTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
    myTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
    myTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
//    myTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
    
    
    globalTree->Branch("nLeptons",&nLeptons, "nLeptons/I");
    myTree->Branch("nLeptons",&nLeptons, "nLeptons/I");

    // electrons
//     globalTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
//     globalTree->Branch("electon_charge",&electon_charge,"electon_charge[nElectrons]/I");
//     globalTree->Branch("electon_pt",&electon_pt,"electon_pt[nElectrons]/D");
//     globalTree->Branch("electon_phi",&electon_phi,"electon_phi[nElectrons]/D");
//     globalTree->Branch("electon_eta",&electon_eta,"electon_eta[nElectrons]/D");
//     globalTree->Branch("electon_eta_superCluster",&electon_eta_superCluster,"electon_eta_superCluster[nElectrons]/D");
//     globalTree->Branch("electon_E",&electon_E,"electon_E[nElectrons]/D");
//     globalTree->Branch("electon_M",&electon_M,"electon_M[nElectrons]/D");
//     globalTree->Branch("electon_d0",&electon_d0,"electon_d0[nElectrons]/D");
//     globalTree->Branch("electon_d0BeamSpot",&electon_d0BeamSpot,"electon_d0BeamSpot[nElectrons]/D");
//     globalTree->Branch("electon_chargedHadronIso",&electon_chargedHadronIso,"electon_chargedHadronIso[nElectrons]/D");
//     globalTree->Branch("electon_neutralHadronIso",&electon_neutralHadronIso,"electon_neutralHadronIso[nElectrons]/D");
//     globalTree->Branch("electon_photonIso",&electon_photonIso,"electon_photonIso[nElectrons]/D");
//     globalTree->Branch("electon_pfIso",&electon_pfIso,"electon_pfIso[nElectrons]/D");
//     globalTree->Branch("electon_sigmaIEtaIEta",&electon_sigmaIEtaIEta,"electon_sigmaIEtaIEta[nElectrons]/D");
//     globalTree->Branch("electon_deltaEtaIn",&electon_deltaEtaIn,"electon_deltaEtaIn[nElectrons]/D");
//     globalTree->Branch("electon_deltaPhiIn",&electon_deltaPhiIn,"electon_deltaPhiIn[nElectrons]/D");
//     globalTree->Branch("electon_hadronicOverEm",&electon_hadronicOverEm,"electon_hadronicOverEm[nElectrons]/D");
//     globalTree->Branch("electon_missingHits",&electon_missingHits,"electon_missingHits[nElectrons]/I");
//     globalTree->Branch("electon_passConversion",&electon_passConversion,"electon_passConversion[nElectrons]/O)");
//     globalTree->Branch("electon_isId",&electon_isId,"electon_isId[nElectrons]/O)");
//     globalTree->Branch("electon_isIso",&electon_isIso,"electon_isIso[nElectrons]/O)");
//     globalTree->Branch("electon_isEBEEGap",&electon_electon_isEBEEGap,"electon_isEBEEGap[nElectrons]/O)");
//     
//     myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
//     myTree->Branch("electon_charge",&electon_charge,"electon_charge[nElectrons]/I");
//     myTree->Branch("electon_pt",&electon_pt,"electon_pt[nElectrons]/D");
//     myTree->Branch("electon_phi",&electon_phi,"electon_phi[nElectrons]/D");
//     myTree->Branch("electon_eta",&electon_eta,"electon_eta[nElectrons]/D");
//     myTree->Branch("electon_eta_superCluster",&electon_eta_superCluster,"electon_eta_superCluster[nElectrons]/D");
//     myTree->Branch("electon_E",&electon_E,"electon_E[nElectrons]/D");
//     myTree->Branch("electon_M",&electon_M,"electon_M[nElectrons]/D");
//     myTree->Branch("electon_d0",&electon_d0,"electon_d0[nElectrons]/D");
//     myTree->Branch("electon_d0BeamSpot",&electon_d0BeamSpot,"electon_d0BeamSpot[nElectrons]/D");
//     myTree->Branch("electon_chargedHadronIso",&electon_chargedHadronIso,"electon_chargedHadronIso[nElectrons]/D");
//     myTree->Branch("electon_neutralHadronIso",&electon_neutralHadronIso,"electon_neutralHadronIso[nElectrons]/D");
//     myTree->Branch("electon_photonIso",&electon_photonIso,"electon_photonIso[nElectrons]/D");
//     myTree->Branch("electon_pfIso",&electon_pfIso,"electon_pfIso[nElectrons]/D");
//     myTree->Branch("electon_sigmaIEtaIEta",&electon_sigmaIEtaIEta,"electon_sigmaIEtaIEta[nElectrons]/D");
//     myTree->Branch("electon_deltaEtaIn",&electon_deltaEtaIn,"electon_deltaEtaIn[nElectrons]/D");
//     myTree->Branch("electon_deltaPhiIn",&electon_deltaPhiIn,"electon_deltaPhiIn[nElectrons]/D");
//     myTree->Branch("electon_hadronicOverEm",&electon_hadronicOverEm,"electon_hadronicOverEm[nElectrons]/D");
//     myTree->Branch("electon_missingHits",&electon_missingHits,"electon_missingHits[nElectrons]/I");
//     myTree->Branch("electon_passConversion",&electon_passConversion,"electon_passConversion[nElectrons]/O)");
//     myTree->Branch("electon_isId",&electon_isId,"electon_isId[nElectrons]/O)");
//     myTree->Branch("electon_isIso",&electon_isIso,"electon_isIso[nElectrons]/O)");
//     myTree->Branch("electon_isEBEEGap",&electon_electon_isEBEEGap,"electon_isEBEEGap[nElectrons]/O)");
    
    
    // muons
    globalTree->Branch("nMuons",&nMuons, "nMuons/I");
    globalTree->Branch("muon_charge",&muon_charge,"muon_charge[nMuons]/I");
    globalTree->Branch("muon_pt",&muon_pt,"muon_pt[nMuons]/D");
    globalTree->Branch("muon_phi",&muon_phi,"muon_phi[nMuons]/D");
    globalTree->Branch("muon_eta",&muon_eta,"muon_eta[nMuons]/D");
    globalTree->Branch("muon_E",&muon_E,"muon_E[nMuons]/D");
    globalTree->Branch("muon_M",&muon_M,"muon_M[nMuons]/D");
    globalTree->Branch("muon_d0",&muon_d0,"muon_d0[nMuons]/D");
    globalTree->Branch("muon_d0BeamSpot",&muon_d0BeamSpot,"muon_d0BeamSpot[nMuons]/D");
    globalTree->Branch("muon_chargedHadronIso",&muon_chargedHadronIso,"muon_chargedHadronIso[nMuons]/D");
    globalTree->Branch("muon_neutralHadronIso",&muon_neutralHadronIso,"muon_neutralHadronIso[nMuons]/D");
    globalTree->Branch("muon_photonIso",&muon_photonIso,"muon_photonIso[nMuons]/D");
    globalTree->Branch("muon_relIso",&muon_relIso,"muon_relIso[nMuons]/D");
    globalTree->Branch("muon_pfIso",&muon_pfIso,"muon_pfIso[nMuons]/D");
    globalTree->Branch("muon_isId",&muon_isId,"muon_isId[nMuons]/O");
    globalTree->Branch("muon_isIso",&muon_isIso,"muon_isIso[nMuons]/O");
    
    myTree->Branch("nMuons",&nMuons, "nMuons/I");
    myTree->Branch("muon_charge",&muon_charge,"muon_charge[nMuons]/I");
    myTree->Branch("muon_pt",&muon_pt,"muon_pt[nMuons]/D");
    myTree->Branch("muon_phi",&muon_phi,"muon_phi[nMuons]/D");
    myTree->Branch("muon_eta",&muon_eta,"muon_eta[nMuons]/D");
    myTree->Branch("muon_E",&muon_E,"muon_E[nMuons]/D");
    myTree->Branch("muon_M",&muon_M,"muon_M[nMuons]/D");
    myTree->Branch("muon_d0",&muon_d0,"muon_d0[nMuons]/D");
    myTree->Branch("muon_d0BeamSpot",&muon_d0BeamSpot,"muon_d0BeamSpot[nMuons]/D");
    myTree->Branch("muon_chargedHadronIso",&muon_chargedHadronIso,"muon_chargedHadronIso[nMuons]/D");
    myTree->Branch("muon_neutralHadronIso",&muon_neutralHadronIso,"muon_neutralHadronIso[nMuons]/D");
    myTree->Branch("muon_photonIso",&muon_photonIso,"muon_photonIso[nMuons]/D");
    myTree->Branch("muon_relIso",&muon_relIso,"muon_relIso[nMuons]/D");
    myTree->Branch("muon_pfIso",&muon_pfIso,"muon_pfIso[nMuons]/D");
    myTree->Branch("muon_isId",&muon_isId,"muon_isId[nMuons]/O");
    myTree->Branch("muon_isIso",&muon_isIso,"muon_isIso[nMuons]/O");
    
    
    // jets
    globalTree->Branch("nJets",&nJets,"nJets/I");
    globalTree->Branch("jet_charge",&jet_charge,"jet_charge[nJets]/I");
    globalTree->Branch("jet_pt",&jet_pt,"jet_pt[nJets]/D");
    globalTree->Branch("jet_phi",&jet_phi,"jet_phi[nJets]/D");
    globalTree->Branch("jet_eta",&jet_eta,"jet_eta[nJets]/D");
    globalTree->Branch("jet_E",&jet_E,"jet_E[nJets]/D");
    globalTree->Branch("jet_M",&jet_M,"jet_M[nJets]/D");
    globalTree->Branch("jet_bdisc",&jet_bdisc,"jet_bdisc[nJets]/D");
    
    myTree->Branch("nJets",&nJets,"nJets/I");
    myTree->Branch("jet_charge",&jet_charge,"jet_charge[nJets]/I");
    myTree->Branch("jet_pt",&jet_pt,"jet_pt[nJets]/D");
    myTree->Branch("jet_phi",&jet_phi,"jet_phi[nJets]/D");
    myTree->Branch("jet_eta",&jet_eta,"jet_eta[nJets]/D");
    myTree->Branch("jet_E",&jet_E,"jet_E[nJets]/D");
    myTree->Branch("jet_M",&jet_M,"jet_M[nJets]/D");
    myTree->Branch("jet_bdisc",&jet_bdisc,"jet_bdisc[nJets]/D");
    
    
    // met
    globalTree->Branch("met_pt", &met_pt, "met_pt/D");
    globalTree->Branch("met_phi", &met_phi, "met_phi/D");
    globalTree->Branch("met_eta", &met_eta,"met_eta/D");
    globalTree->Branch("met_E", &met_E,"met_E/D");
    
    myTree->Branch("met_pt", &met_pt, "met_pt/D");
    myTree->Branch("met_phi", &met_phi, "met_phi/D");
    myTree->Branch("met_eta", &met_eta,"met_eta/D");
    myTree->Branch("met_E", &met_E,"met_E/D");
    
    
    // mcparticles
    if (! isData)
    {
      globalTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
      globalTree->Branch("mc_status",&mc_status,"mc_status[nMCParticles]/I");
      globalTree->Branch("mc_pdgId",&mc_pdgId,"mc_pdgId[nMCParticles]/I");
      globalTree->Branch("mc_mother",&mc_mother,"mc_mother[nMCParticles]/I");
      globalTree->Branch("mc_granny",&mc_granny,"mc_granny[]/I");
      globalTree->Branch("mc_pt",&mc_pt,"mc_pt[nMCParticles]/D");
      globalTree->Branch("mc_phi",&mc_phi,"mc_phi[nMCParticles]/D");
      globalTree->Branch("mc_eta",&mc_eta,"mc_eta[nMCParticles]/D");
      globalTree->Branch("mc_E",&mc_E,"mc_E[nMCParticles]/D");
      globalTree->Branch("mc_M",&mc_M,"mc_M[nMCParticles]/D");
      
      myTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
      myTree->Branch("mc_status",&mc_status,"mc_status[nMCParticles]/I");
      myTree->Branch("mc_pdgId",&mc_pdgId,"mc_pdgId[nMCParticles]/I");
      myTree->Branch("mc_mother",&mc_mother,"mc_mother[nMCParticles]/I");
      myTree->Branch("mc_granny",&mc_granny,"mc_granny[nMCParticles]/I");
      myTree->Branch("mc_pt",&mc_pt,"mc_pt[nMCParticles]/D");
      myTree->Branch("mc_phi",&mc_phi,"mc_phi[nMCParticles]/D");
      myTree->Branch("mc_eta",&mc_eta,"mc_eta[nMCParticles]/D");
      myTree->Branch("mc_E",&mc_E,"mc_E[nMCParticles]/D");
      myTree->Branch("mc_M",&mc_M,"mc_M[nMCParticles]/D");
    }
    
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    nEvents = 0;
    nEventsSel = 0;
    
    /// Get run information
    datasets[d]->runTree()->SetBranchStatus("runInfos*",1);
    datasets[d]->runTree()->SetBranchAddress("runInfos",&runInfos);
    
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    //for (unsigned int ievt = 0; ievt < 201; ievt++)
    {
      nEvents++;
      
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      
      
      
      
      
      
      
      
      
      
      
      
    }  // end loop events
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
    fout->Close();
    delete fout;
    
  }  // end loop datasets
  
  /// To write plots b tagging:
  delete bTagHistoTool_M;
  
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;
  
  return 0;
}
