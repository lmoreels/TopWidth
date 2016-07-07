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
  
  
  
  ///////////////////
  ///  Selection  ///
  ///////////////////
  
  float muonPTSel = 26.; // GeV
  float muonEtaSel = 2.1;
  float muonRelIsoSel = 0.15;  // Tight muon
  string muonWP = "Tight";
  
  float muonPTVeto = 10.; // GeV
  float muonEtaVeto = 2.5;
  float muonRelIsoVeto = 0.25;  // Loose muon
  
  float electronPTSel = 24.; // GeV
  float electronEtaSel = 2.5;
  string electronWP = "Tight";
  
  float electronPTVeto = 15.; // GeV
  float electronEtaVeto = 2.5;
  
  float jetPT = 20.; // GeV, 30 GeV for selected jets
  float jetEta = 2.4;  // to allow b tagging
  
  
  
  //////////////////////////////////////
  ///  Working points for b tagging  ///
  //////////////////////////////////////
  
  /// Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
  
  float CSVv2Loose =  0.460;
  float CSVv2Medium = 0.800;
  float CSVv2Tight = 0.935;
  
  
  
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
    
    if ( (datasets[d]->Title()).find("amc") == 0 )
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
    TTree* statTree = new TTree("stats","stats");
//    TTree* globalTree = new TTree("globaltree","globaltree");  // no selection applied
    
    
    
    ////////////////////////////////////
    ///  Define variables for trees  ///
    ////////////////////////////////////
    
    // event related variables
    Int_t run_num;
    Int_t evt_num;
    Int_t lumi_num;
    Int_t nvtx;
    Int_t npu;
    Double_t rho;
    
    Bool_t isTrigged;
    Bool_t isSelected;
    Bool_t hasExactly4Jets;
//    Bool_t passedMETFilter;
    Bool_t cutFlow[10];
    
    Int_t appliedJER;
    Int_t appliedJES;
    
    Int_t nEvents;
    Int_t nEventsSel;
    Int_t nofPosWeights;
    Int_t nofNegWeights;
    Double_t sumW;
    
    Int_t nofEventsHLTv2;
    Int_t nofEventsHLTv3;
    
    Double_t puSF;
    Double_t btagSF;
    Double_t muonIdSF[10];
    Double_t muonIsoSF[10];
    Double_t muonTrigSFv2[10];
    Double_t muonTrigSFv3[10];
//    Double_t electronSF[10];
    Double_t nloWeight; // for amc@nlo samples
    
    
    Int_t nLeptons;
    
    /// Variables for electrons
    Int_t nElectrons;
//    Int_t electron_charge[10];
//    Double_t electron_pt[10];
//    Double_t electron_phi[10];
//    Double_t electron_eta[10];
//    Double_t electron_eta_superCluster[10];
//    Double_t electron_E[10];
//    Double_t electron_M[10];
//    Double_t electron_d0[10];
//    Double_t electron_chargedHadronIso[10];
//    Double_t electron_neutralHadronIso[10];
//    Double_t electron_photonIso[10];
//    Double_t electron_pfIso[10];
//
//    Double_t electron_sigmaIEtaIEta[10];
//    Double_t electron_deltaEtaIn[10];
//    Double_t electron_deltaPhiIn[10];
//    Double_t electron_hadronicOverEm[10];
//    Int_t electron_missingHits[10];
//    Bool_t electron_passConversion[10];
//    Bool_t electron_isEBEEGap[10];
    
    /// Variables for muons
    Int_t nMuons;
    Int_t muon_charge[10];
    Double_t muon_pt[10];
    Double_t muon_phi[10];
    Double_t muon_eta[10];
    Double_t muon_E[10];
    Double_t muon_M[10];
    Double_t muon_d0[10];
    Double_t muon_chargedHadronIso[10];
    Double_t muon_neutralHadronIso[10];
    Double_t muon_photonIso[10];
    Double_t muon_puChargedHadronIso[10];
    Double_t muon_relIso[10];
    Double_t muon_pfIso[10];
    
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
    Double_t met_Et;
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
    basicEvent->Branch("rho",&rho,"rho/D");
    basicEvent->Branch("isTrigged",&isTrigged,"isTrigged/O");
    basicEvent->Branch("cutFlow",&cutFlow,"cutFlow[10]/O");
    basicEvent->Branch("appliedJER",&appliedJER,"appliedJER/I");
    basicEvent->Branch("appliedJES", &appliedJES, "appliedJES/I");
    
    statTree->Branch("nEvents" , &nEvents, "nEvents/I");
    statTree->Branch("nEventsSel" , &nEventsSel, "nEventsSel/I");
    if (isData)
    {
      statTree->Branch("nofEventsHLTv2",&nofEventsHLTv2,"nofEventsHLTv2/I");
      statTree->Branch("nofEventsHLTv3",&nofEventsHLTv3,"nofEventsHLTv3/I");
    }
    if (nlo)
    {
      statTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
      statTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
      statTree->Branch("sumW", &sumW, "sumW/D");
    }
    
//    globalTree->Branch("run_num",&run_num,"run_num/I");
//    globalTree->Branch("evt_num",&evt_num,"evt_num/I");
//    globalTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
//    globalTree->Branch("nvtx",&nvtx,"nvtx/I");
//    globalTree->Branch("npu",&npu,"npu/I");
//    globalTree->Branch("rho",&rho,"rho/D");
//    globalTree->Branch("isTrigged",&isTrigged,"isTrigged/O");
//    globalTree->Branch("appliedJER",&appliedJER,"appliedJER/I");
//    globalTree->Branch("appliedJES", &appliedJES, "appliedJES/I");
    
    myTree->Branch("run_num",&run_num,"run_num/I");
    myTree->Branch("evt_num",&evt_num,"evt_num/I");
    myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
    myTree->Branch("nvtx",&nvtx,"nvtx/I");
    myTree->Branch("npu",&npu,"npu/I");
    myTree->Branch("rho",&rho,"rho/D");
    myTree->Branch("isTrigged",&isTrigged,"isTrigged/O");
    myTree->Branch("cutFlow",&cutFlow,"cutFlow[10]/O");
    myTree->Branch("hasExactly4Jets",&hasExactly4Jets,"hasExactly4Jets/O");
//    myTree->Branch("passedMETFilter", &passedMETFilter,"passedMETFilter/O");
    myTree->Branch("appliedJER",&appliedJER,"appliedJER/I");
    myTree->Branch("appliedJES", &appliedJES, "appliedJES/I");
    
    
    
    /// SFs
//    globalTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
//    globalTree->Branch("puSF",&puSF,"puSF/D");
//    globalTree->Branch("btagSF",&btagSF,"btagSF/D");
//    globalTree->Branch("muonIdSF",&muonIdSF,"muonIdSF[nMuons]/D");
//    globalTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
//    globalTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
//    globalTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
//    globalTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
    
    myTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
    myTree->Branch("puSF",&puSF,"puSF/D");
    myTree->Branch("btagSF",&btagSF,"btagSF/D");
    myTree->Branch("muonIdSF",&muonIdSF,"muonIdSF[nMuons]/D");
    myTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
    myTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
    myTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
//    myTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
    
    
//    globalTree->Branch("nLeptons",&nLeptons, "nLeptons/I");
    myTree->Branch("nLeptons",&nLeptons, "nLeptons/I");

    // electrons
//    globalTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
//    globalTree->Branch("electon_charge",&electon_charge,"electon_charge[nElectrons]/I");
//    globalTree->Branch("electon_pt",&electon_pt,"electon_pt[nElectrons]/D");
//    globalTree->Branch("electon_phi",&electon_phi,"electon_phi[nElectrons]/D");
//    globalTree->Branch("electon_eta",&electon_eta,"electon_eta[nElectrons]/D");
//    globalTree->Branch("electon_eta_superCluster",&electon_eta_superCluster,"electon_eta_superCluster[nElectrons]/D");
//    globalTree->Branch("electon_E",&electon_E,"electon_E[nElectrons]/D");
//    globalTree->Branch("electon_M",&electon_M,"electon_M[nElectrons]/D");
//    globalTree->Branch("electon_d0",&electon_d0,"electon_d0[nElectrons]/D");
//    globalTree->Branch("electon_chargedHadronIso",&electon_chargedHadronIso,"electon_chargedHadronIso[nElectrons]/D");
//    globalTree->Branch("electon_neutralHadronIso",&electon_neutralHadronIso,"electon_neutralHadronIso[nElectrons]/D");
//    globalTree->Branch("electon_photonIso",&electon_photonIso,"electon_photonIso[nElectrons]/D");
//    globalTree->Branch("electon_pfIso",&electon_pfIso,"electon_pfIso[nElectrons]/D");
//    globalTree->Branch("electon_sigmaIEtaIEta",&electon_sigmaIEtaIEta,"electon_sigmaIEtaIEta[nElectrons]/D");
//    globalTree->Branch("electon_deltaEtaIn",&electon_deltaEtaIn,"electon_deltaEtaIn[nElectrons]/D");
//    globalTree->Branch("electon_deltaPhiIn",&electon_deltaPhiIn,"electon_deltaPhiIn[nElectrons]/D");
//    globalTree->Branch("electon_hadronicOverEm",&electon_hadronicOverEm,"electon_hadronicOverEm[nElectrons]/D");
//    globalTree->Branch("electon_missingHits",&electon_missingHits,"electon_missingHits[nElectrons]/I");
//    globalTree->Branch("electon_passConversion",&electon_passConversion,"electon_passConversion[nElectrons]/O)");
//    globalTree->Branch("electon_isEBEEGap",&electon_electon_isEBEEGap,"electon_isEBEEGap[nElectrons]/O)");
    
//    myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
//    myTree->Branch("electon_charge",&electon_charge,"electon_charge[nElectrons]/I");
//    myTree->Branch("electon_pt",&electon_pt,"electon_pt[nElectrons]/D");
//    myTree->Branch("electon_phi",&electon_phi,"electon_phi[nElectrons]/D");
//    myTree->Branch("electon_eta",&electon_eta,"electon_eta[nElectrons]/D");
//    myTree->Branch("electon_eta_superCluster",&electon_eta_superCluster,"electon_eta_superCluster[nElectrons]/D");
//    myTree->Branch("electon_E",&electon_E,"electon_E[nElectrons]/D");
//    myTree->Branch("electon_M",&electon_M,"electon_M[nElectrons]/D");
//    myTree->Branch("electon_d0",&electon_d0,"electon_d0[nElectrons]/D");
//    myTree->Branch("electon_chargedHadronIso",&electon_chargedHadronIso,"electon_chargedHadronIso[nElectrons]/D");
//    myTree->Branch("electon_neutralHadronIso",&electon_neutralHadronIso,"electon_neutralHadronIso[nElectrons]/D");
//    myTree->Branch("electon_photonIso",&electon_photonIso,"electon_photonIso[nElectrons]/D");
//    myTree->Branch("electon_pfIso",&electon_pfIso,"electon_pfIso[nElectrons]/D");
//    myTree->Branch("electon_sigmaIEtaIEta",&electon_sigmaIEtaIEta,"electon_sigmaIEtaIEta[nElectrons]/D");
//    myTree->Branch("electon_deltaEtaIn",&electon_deltaEtaIn,"electon_deltaEtaIn[nElectrons]/D");
//    myTree->Branch("electon_deltaPhiIn",&electon_deltaPhiIn,"electon_deltaPhiIn[nElectrons]/D");
//    myTree->Branch("electon_hadronicOverEm",&electon_hadronicOverEm,"electon_hadronicOverEm[nElectrons]/D");
//    myTree->Branch("electon_missingHits",&electon_missingHits,"electon_missingHits[nElectrons]/I");
//    myTree->Branch("electon_passConversion",&electon_passConversion,"electon_passConversion[nElectrons]/O)");
//    myTree->Branch("electon_isEBEEGap",&electon_electon_isEBEEGap,"electon_isEBEEGap[nElectrons]/O)");
    
    
    // muons
//    globalTree->Branch("nMuons",&nMuons, "nMuons/I");
//    globalTree->Branch("muon_charge",&muon_charge,"muon_charge[nMuons]/I");
//    globalTree->Branch("muon_pt",&muon_pt,"muon_pt[nMuons]/D");
//    globalTree->Branch("muon_phi",&muon_phi,"muon_phi[nMuons]/D");
//    globalTree->Branch("muon_eta",&muon_eta,"muon_eta[nMuons]/D");
//    globalTree->Branch("muon_E",&muon_E,"muon_E[nMuons]/D");
//    globalTree->Branch("muon_M",&muon_M,"muon_M[nMuons]/D");
//    globalTree->Branch("muon_d0",&muon_d0,"muon_d0[nMuons]/D");
//    globalTree->Branch("muon_chargedHadronIso",&muon_chargedHadronIso,"muon_chargedHadronIso[nMuons]/D");
//    globalTree->Branch("muon_neutralHadronIso",&muon_neutralHadronIso,"muon_neutralHadronIso[nMuons]/D");
//    globalTree->Branch("muon_photonIso",&muon_photonIso,"muon_photonIso[nMuons]/D");
//    globalTree->Branch("muon_relIso",&muon_relIso,"muon_relIso[nMuons]/D");
//    globalTree->Branch("muon_pfIso",&muon_pfIso,"muon_pfIso[nMuons]/D");
    
    myTree->Branch("nMuons",&nMuons, "nMuons/I");
    myTree->Branch("muon_charge",&muon_charge,"muon_charge[nMuons]/I");
    myTree->Branch("muon_pt",&muon_pt,"muon_pt[nMuons]/D");
    myTree->Branch("muon_phi",&muon_phi,"muon_phi[nMuons]/D");
    myTree->Branch("muon_eta",&muon_eta,"muon_eta[nMuons]/D");
    myTree->Branch("muon_E",&muon_E,"muon_E[nMuons]/D");
    myTree->Branch("muon_M",&muon_M,"muon_M[nMuons]/D");
    myTree->Branch("muon_d0",&muon_d0,"muon_d0[nMuons]/D");
    myTree->Branch("muon_chargedHadronIso",&muon_chargedHadronIso,"muon_chargedHadronIso[nMuons]/D");
    myTree->Branch("muon_neutralHadronIso",&muon_neutralHadronIso,"muon_neutralHadronIso[nMuons]/D");
    myTree->Branch("muon_photonIso",&muon_photonIso,"muon_photonIso[nMuons]/D");
    myTree->Branch("muon_puChargedHadronIso",&muon_puChargedHadronIso,"muon_puChargedHadronIso[nMuons]/D");
    myTree->Branch("muon_relIso",&muon_relIso,"muon_relIso[nMuons]/D");
    myTree->Branch("muon_pfIso",&muon_pfIso,"muon_pfIso[nMuons]/D");    
    
    // jets
//    globalTree->Branch("nJets",&nJets,"nJets/I");
//    globalTree->Branch("jet_charge",&jet_charge,"jet_charge[nJets]/I");
//    globalTree->Branch("jet_pt",&jet_pt,"jet_pt[nJets]/D");
//    globalTree->Branch("jet_phi",&jet_phi,"jet_phi[nJets]/D");
//    globalTree->Branch("jet_eta",&jet_eta,"jet_eta[nJets]/D");
//    globalTree->Branch("jet_E",&jet_E,"jet_E[nJets]/D");
//    globalTree->Branch("jet_M",&jet_M,"jet_M[nJets]/D");
//    globalTree->Branch("jet_bdiscr",&jet_bdiscr,"jet_bdiscr[nJets]/D");
    
    myTree->Branch("nJets",&nJets,"nJets/I");
    myTree->Branch("jet_charge",&jet_charge,"jet_charge[nJets]/I");
    myTree->Branch("jet_pt",&jet_pt,"jet_pt[nJets]/D");
    myTree->Branch("jet_phi",&jet_phi,"jet_phi[nJets]/D");
    myTree->Branch("jet_eta",&jet_eta,"jet_eta[nJets]/D");
    myTree->Branch("jet_E",&jet_E,"jet_E[nJets]/D");
    myTree->Branch("jet_M",&jet_M,"jet_M[nJets]/D");
    myTree->Branch("jet_bdiscr",&jet_bdiscr,"jet_bdiscr[nJets]/D");
    
    
    // met
//    globalTree->Branch("met_pt", &met_pt, "met_pt/D");
//    globalTree->Branch("met_phi", &met_phi, "met_phi/D");
//    globalTree->Branch("met_eta", &met_eta,"met_eta/D");
//    globalTree->Branch("met_E", &met_E,"met_E/D");
    
    myTree->Branch("met_pt", &met_pt, "met_pt/D");
    myTree->Branch("met_phi", &met_phi, "met_phi/D");
    myTree->Branch("met_eta", &met_eta,"met_eta/D");
    myTree->Branch("met_Et", &met_Et,"met_Et/D");
    myTree->Branch("met_E", &met_E,"met_E/D");
    
    
    // mcparticles
    if (! isData)
    {
//      globalTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
//      globalTree->Branch("mc_status",&mc_status,"mc_status[nMCParticles]/I");
//      globalTree->Branch("mc_pdgId",&mc_pdgId,"mc_pdgId[nMCParticles]/I");
//      globalTree->Branch("mc_mother",&mc_mother,"mc_mother[nMCParticles]/I");
//      globalTree->Branch("mc_granny",&mc_granny,"mc_granny[]/I");
//      globalTree->Branch("mc_pt",&mc_pt,"mc_pt[nMCParticles]/D");
//      globalTree->Branch("mc_phi",&mc_phi,"mc_phi[nMCParticles]/D");
//      globalTree->Branch("mc_eta",&mc_eta,"mc_eta[nMCParticles]/D");
//      globalTree->Branch("mc_E",&mc_E,"mc_E[nMCParticles]/D");
//      globalTree->Branch("mc_M",&mc_M,"mc_M[nMCParticles]/D");
      
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
    nofPosWeights = 0;
    nofNegWeights = 0;
    sumW = 0.;
    nofEventsHLTv2 = 0;
    nofEventsHLTv3 = 0;
    
    /// Get run information
    datasets[d]->runTree()->SetBranchStatus("runInfos*",1);
    datasets[d]->runTree()->SetBranchAddress("runInfos",&runInfos);
    
    /// Define objects
    vector < TRootVertex* > vertex;
    vector < TRootMuon* > init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* > init_jets_corrected;
    vector < TRootJet* > init_jets;
    vector < TRootMET* > mets;
    vector < TRootGenJet* > genjets;
    vector < TRootMCParticle* > mcParticles;
    
    vector < TRootPFJet* > selectedJets;
    vector < TRootPFJet* > selectedBJets;
    vector < TRootMuon* > selectedMuons;
    vector < TRootMuon* > vetoMuons;
    vector < TRootElectron* > selectedElectrons;
    vector < TRootElectron* > vetoElectrons;
    
    
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    //for (unsigned int ievt = 0; ievt < 201; ievt++)
    {
      nEvents++;
      
      /// Clear objects
      vertex.clear();
      init_muons.clear();
      init_electrons.clear();
      init_jets_corrected.clear();
      init_jets.clear();
      mets.clear();
      genjets.clear();
      mcParticles.clear();
      
      selectedJets.clear();
      selectedBJets.clear();
      selectedMuons.clear();
      vetoMuons.clear();
      selectedElectrons.clear();
      vetoElectrons.clear();
      
      /// Reset other stuff
      isTrigged = false;
      isSelected = false;
      hasExactly4Jets = false;
      //passedMETFilter = false;
      cutFlow[] = {0};
      appliedJER = 0;
      appliedJES = 0;
      puSF = 1.;
      btagSF = 1.;
      muonIdSF[] = {1.};
      muonIsoSF[] = {1.};
      muonTrigSFv2[] = {1.};
      muonTrigSFv3[] = {1.};
      //electronSF[] = {1.};
      nloWeight = 1.; // for amc@nlo samples
      
      nLeptons = -1;
      nElectrons = -1;
      nMuons = -1;
      nJets = -1;
      
//      electron_charge[] = {0};
//      electron_pt[] = {0.};
//      electron_phi[] = {0.};
//      electron_eta[] = {0.};
//      electron_eta_superCluster[] = {0.};
//      electron_E[] = {0.};
//      electron_M[] = {0.};
//      electron_d0[] = {-1.};
//      electron_chargedHadronIso[] = {-1.};
//      electron_neutralHadronIso[] = {-1.};
//      electron_photonIso[] = {-1.};
//      electron_pfIso[] = {-1.};
//      electron_sigmaIEtaIEta[] = {-1.};
//      electron_deltaEtaIn[] = {-1.};
//      electron_deltaPhiIn[] = {-1.};
//      electron_hadronicOverEm[] = {-1.};
//      electron_missingHits[] = {-1};
//      electron_passConversion[] = {0};
//      electron_isEBEEGap[] = {0};

      muon_charge[] = {0};
      muon_pt[] = {0.};
      muon_phi[] = {0.};
      muon_eta[] = {0.};
      muon_E[] = {0.};
      muon_M[] = {0.};
      muon_d0[] = {-1.};
      muon_chargedHadronIso[] = {-1.};
      muon_neutralHadronIso[] = {-1.};
      muon_photonIso[] = {-1.};
      muon_puChargedHadronIso = {-1.};
      muon_relIso[] = {-1.};
      muon_pfIso[] = {-1.};
      
      jet_charge[] = {0};
      jet_pt[] = {0.};
      jet_phi[] = {0.};
      jet_eta[] = {0.};
      jet_E[] = {0.};
      jet_M[] = {0.};
      jet_bdiscr[] = {-1.};
      
      met_pt = 0;
      met_phi = 0;
      met_eta = 0;
      met_Et = 0;
      met_E = 0;
      
      /// mcparticles
      nMCParticles = -1;
      mc_status[] = {-1};
      mc_pdgId[] = {0};
      mc_mother[] = {0};
      mc_granny[] = {0};
      mc_pt[] = {0.};
      mc_phi[] = {0.};
      mc_eta[] = {0.};
      mc_E[] = {0.};
      mc_M[] = {0.};
      
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets, mets);
      init_jets_corrected = init_jets;
      
      datasets[d]->eventTree()->LoadTree(ievt);
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      run_num = event->runId();
      evt_num = event->eventId();
      lumi_num = event->lumiBlockId();
      nvtx = vertex.size();
      npu = (int)event->nTruePU();
      rho = event->fixedGridRhoFastjetAll();
      
      if (isData)
      {
        if ( run_num < 256630 )
        {
          cerr << "-- Dataset 2015C included..." << endl;
          exit(1);
        }
        else if ( run_num >= 256630 && run_num <= 257819 )
        {
          nofEventsHLTv2++;
        }
        else
        {
          nofEventsHLTv3++;
        }
      }

      if (! isData )
      {
        genjets = treeLoader.LoadGenJet(ievt,false);
        treeLoader.LoadMCEvent(ievt, 0, mcParticles, false);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      
      
      //////////////////////////////////////
      ///  SCALEFACTORS AND CORRECTIONS  ///
      //////////////////////////////////////
      
      if (nlo)
      {
        if ( event->getWeight(1001) != -9999. )
        {
          nloWeight = event->getWeight(1001)/abs(event->originalXWGTUP());
          //mc_scaleupweight = event->getWeight(1005)/abs(event->originalXWGTUP());
          //mc_scaledownweight = event->getWeight(1009)/abs(event->originalXWGTUP());
          if ( nloWeight >= 0. ) 
          {
            nofPosWeights++;
          }
          else
          {
            nofNegWeights++;
          }
        }
        if ( event->getWeight(1) != -9999. )
        {
          nloWeight = event->getWeight(1)/abs(event->originalXWGTUP());
          //mc_scaleupweight = event->getWeight(5)/abs(event->originalXWGTUP());
          //mc_scaledownweight = event->getWeight(9)/abs(event->originalXWGTUP());
          if ( nloWeight >= 0. )
          {
            nofPosWeights++;
          }
          else
          {
            nofNegWeights++;
          }
        }
        
        sumW += nloWeight;
      }
      
      jetTools->correctJets(init_jets_corrected, rho, isData);
      
      if (! isData)
      {
        puSF = LumiWeights.ITweight( (int)event->nTruePU() );
        // up syst -> puSF = LumiWeightsUp.ITweight( (int)event->nTruePU() );
        // down syst -> puSF = LumiWeightsDown.ITweight( (int)event->nTruePU() );
        
        if (applyJERdown)    jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus", false);
        else if (applyJERup) jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus", false);
        else                 jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal", false);
        
        /// Example how to apply JES systematics
        //jetTools->correctJetJESUnc(init_jets_corrected, "minus", 1);
        //jetTools->correctJetJESUnc(init_jets_corrected, "plus", 1);
      }
      
      
      
      /////////////////
      ///  Trigger  ///
      /////////////////
      
      trigger->checkAvail(run_num, datasets, d, &treeLoader, event, false);
      isTrigged = trigger->checkIfFired();
      
      
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      //Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets, rho);
      
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2.);
      
      selectedJets = selection.GetSelectedJets(jetPT, jetEta, true, "Tight");  // PtThr, EtaThr, applyJetID, TightLoose
      selectedMuons = selection.GetSelectedMuons(muonPTSel, muonEtaSel, muonRelIsoSel, muonWP, "Spring15");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      selectedElectrons = selection.GetSelectedElectrons(electronPTSel, electronEtaSel, electronWP, "Spring15_25ns", true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased
      vetoMuons = selection.GetSelectedMuons(muonPTVeto, muonEtaVeto, muonRelIsoVeto, "Loose", "Spring15");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      vetoElectronsSemiMu = selection.GetSelectedElectrons(electronPTVeto, electronEtaVeto, "Veto", "Spring15_25ns", true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased
      
      
      if (applyJetLeptonCleaning)
      {
        if(verbose > 3) cout << "  - Applying jet/lepton cleaning... " << endl; 
        
        vector<TRootPFJet*> selectedJetsBC;
        selectedJetsBC = selectedJets;
        selectedJets.clear();
        
        for (int iOrigJet = 0; iOrigJet < selectedJetsBC.size(); iOrigJet++)
        {
          bool toBeErased = false;
          for (int iMuon = 0; iMuon < selectedMuons.size(); iMuon++)
          {
            if ( selectedJetsBC[iOrigJet]->DeltaR(*selectedMuons[iMuon]) < 0.4 )
            {
              toBeErased = true;
              break;
            }
          }
          if (toBeErased) continue;
          for (int iElectron = 0; iElectron < selectedElectrons.size(); iElectron++)
          {
            if ( selectedJetsBC[iOrigJet]->DeltaR(*selectedElectrons[iElectron]) < 0.3 )
            {
              toBeErased = true;
              break;
            }
          }
          if (! toBeErased)
          {
            selectedJets.push_back(selectedJetsBC[iOrigJet]);
          }
        }
        if ( verbose > 3 )
        {
          if ( selectedJetsBC.size() != selectedJets.size() ) cout << "--> original = " << selectedJetsBC.size()  << " after cleaning = " << selectedJets.size() << endl;
        }
        nofEventsJetLeptonCleaned++;
        
      }  // end jet cleaning
      
      
      for (int i = 0; i < selectedJets.size(); i++)
      {
        if ( selectedJets[i]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium )
          selectedBJets.push_back(selectedJets[i]);
      }
      
      
      
      /// Fill variables for tree
      nJets = selectedJets.size();
      nMuons = selectedMuons.size();
      nElectrons = selectedElectrons.size();
      nLeptons = nMuons + nElectrons;
      
      for(Int_t iJet = 0; iJet < nJets; iJet++)
      {
        jet_charge[iJet] = selectedJets[iJet]->charge();
        jet_pt[iJet] = selectedJets[iJet]->Pt();
        jet_phi[iJet] = selectedJets[iJet]->Phi();
        jet_eta[iJet] = selectedJets[iJet]->Eta();
        jet_E[iJet] = selectedJets[iJet]->E();
        jet_M[iJet] = selectedJets[iJet]->M();
        jet_bdiscr[iJet] = selectedJets[iJet]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
      }
      
      for (Int_t iMuon = 0; iMuon < nMuons; iMuon++)
      {
        muon_charge[iMuon] = selectedMuons[iMuon]->charge();
        muon_pt[iMuon] = selectedMuons[iMuon]->Pt();
        muon_phi[iMuon] = selectedMuons[iMuon]->Phi();
        muon_eta[iMuon] = selectedMuons[iMuon]->Eta();
        muon_E[iMuon] = selectedMuons[iMuon]->E();
        muon_M[iMuon] = selectedMuons[iMuon]->M();
        muon_d0[iMuon] = selectedMuons[iMuon]->d0();
        muon_chargedHadronIso[iMuon] = selectedMuons[iMuon]->chargedHadronIso(4);
        muon_neutralHadronIso[iMuon] = selectedMuons[iMuon]->neutralHadronIso(4);
        muon_photonIso[iMuon] = selectedMuons[iMuon]->->photonIso(4);
        muon_puChargedHadronIso[iMuon] = selectedMuons[iMuon]->puChargedHadronIso(4);
        muon_relIso[iMuon] = ( muon_chargedHadronIso[iMuon] + max( 0.0, muon_neutralHadronIso[iMuon] + muon_photonIso[iMuon] - 0.5*muon_puChargedHadronIso[iMuon] ) ) / muon_pt[iMuon];  // dR = 0.4, dBeta corrected
        muon_pfIso[iMuon] = selectedMuons[iMuon]->relPfIso(4,0);
      }
      
//      for (Int_t iElectron = 0; iElectron < nElectrons; iElectron++)
//      {
//        electron_charge[iElectron] = selectedElectrons[iElectron]->charge();
//        electron_pt[iElectron] = selectedElectrons[iElectron]->Pt();
//        electron_phi[iElectron] = selectedElectrons[iElectron]->Phi();
//        electron_eta[iElectron] = selectedElectrons[iElectron]->Eta();
//        electron_eta_superCluster[iElectron] = selectedElectrons[iElectron]->superClusterEta();
//        electron_E[iElectron] = selectedElectrons[iElectron]->E();
//        electron_M[iElectron] = selectedElectrons[iElectron]->M();
//        electron_d0[iElectron] = selectedElectrons[iElectron]->d0();
//        electron_chargedHadronIso[iElectron] = selectedElectrons[iElectron]->;
//        electron_neutralHadronIso[iElectron] = selectedElectrons[iElectron]->;
//        electron_photonIso[iElectron] = selectedElectrons[iElectron]->;
//        electron_pfIso[iElectron] = selectedElectrons[iElectron]->->relPfIso(3,0);
//        electron_sigmaIEtaIEta[iElectron] = selectedElectrons[iElectron]->;
//        electron_deltaEtaIn[iElectron] = selectedElectrons[iElectron]->;
//        electron_deltaPhiIn[iElectron] = selectedElectrons[iElectron]->;
//        electron_hadronicOverEm[iElectron] = selectedElectrons[iElectron]->;
//        electron_missingHits[iElectron] = selectedElectrons[iElectron]->;
//        electron_passConversion[iElectron] = selectedElectrons[iElectron]->;
//        electron_isEBEEGap[iElectron] = selectedElectrons[iElectron]->;
//      }
      
      met_pt = mets[0]->Pt();
      met_phi = mets[0]->Phi();
      met_eta = mets[0]->Eta();
      met_Et = mets[0]->Et();
      met_E = mets[0]->E();
      
      if (! isData)
      {
        nMCParticles = mcParticles.size();
        for (Int_t iMC = 0; iMC < nMCParticles; iMC)
        {
          mc_status[iMC] = mcParticles[iMC]->status();
          mc_pdgId[iMC] = mcParticles[iMC]->type();
          mc_mother[iMC] = mcParticles[iMC]->motherType();
          mc_granny[iMC] = mcParticles[iMC]->grannyType();
          mc_pt[iMC] = mcParticles[iMC]->Pt();
          mc_phi[iMC] = mcParticles[iMC]->Phi();
          mc_eta[iMC] = mcParticles[iMC]->Eta();
          mc_E[iMC] = mcParticles[iMC]->E();
          mc_M[iMC] = mcParticles[iMC]->M();
        }
      }
      
      /// Fill scalefactors
      if (! isData)
      {
        if (applyBTagSF) btagSF = bTagHistoTool_M->getMCEventWeight(selectedJets);
        
        for (int iMuon = 0; iMuon < selectedMuons.size(); iMuon++)
        {
          muonIdSF[iMuon] = muonSFWeightID_T->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);  // eta, pt, shiftUpDown;
          muonIsoSF[iMuon] = muonSFWeightIso_TT->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
          muonTrigSFv2[iMuon] = muonSFWeightTrigHLTv4p2->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
          muonTrigSFv3[iMuon] = muonSFWeightTrigHLTv4p3->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
        }
        
        //for (int iElectron = 0; iElectron < selectedElectrons(); iElectron++)
        //{
        //  electronSF[iElectron] = ;
        //}
      }
      
//      globalTree->Fill();
      
      
      ////// Selection
      cutFlow[0] = 1;
      if (isTrigged)
      {
        cutFlow[1] = 1;
        if (isGoodPV)
        {
          cutFlow[2] = 1;
          if (selectedMuons.size() == 1)
          {
            cutFlow[3] = 1;
            if (vetoMuons.size() == 1)
            {
              cutFlow[4] = 1;
              if (vetoElectronsSemiMu.size() == 0)
              {
                cutFlow[5] = 1;
                
                /// First 4 jets need pT > 30 GeV
                if (selectedJets.size() >= 4)
                {
                  if (selectedJets[3]->Pt() < 30) selectedJets.clear();
                }
                
                if ( selectedJets.size() >= 4 )
                {
                  cutFlow[6] = 1;
                  if ( selectedJets.size() == 4 ) hasExactly4Jets = true;
                  
                  if ( selectedBJets.size() > 0 )
                  {
                    cutFlow[7] = 1;
                    if ( selectedBJets.size() > 1 )
                    {
                      cutFlow[8] = 1;
                      isSelected = true;
                    }  // at least 2 b-tagged jets
                  }  // at least 1 b-tagged jet
                  
                }  // at least 4 jets
              }  // no veto electrons
            }  // no additional loose muons (tight muon is also loose muon)
          }  // 1 good muon
        }  // good PV
      }  // trigged
      
      eventTree->Fill;
      
      if (! isSelected)
      {
        continue;
      }
      
      nEventsSel++;
      myTree->Fill();
      
      /// B-tagging
      if (calculateBTagSF && ! isData)
      {
        bTagHistoTool_M->FillMCEfficiencyHistos(selectedJets);
      }
      
      
      
      
      
      
      
      
      
    }  // end loop events
    
    statTree->Fill();
    
    /// Write to file
    fout->cd();
    eventTree->Write();
    myTree->Write();
    statTree->Write();
    fout->Close();
    delete fout;
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  // end loop datasets
  
  /// To write plots b tagging:
  delete bTagHistoTool_M;
  
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;
  
  double time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  cout << " - Goodbye" << endl;
  
  return 0;
}
