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
#include "TFile.h"
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
#include "Tools/interface/ResolutionFunctions.h"
#include "Tools/interface/Trigger.h"


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
  
  string pathOutput = "NtupleOutput/";
  mkdir(pathOutput.c_str(),0777);
  
  string xmlFileName ="config/topWidth.xml";
  int maxMCParticles = -1;
  
  
  ///////////////////////////
  ///  Process arguments  ///
  ///////////////////////////
  
  string dName, dTitle, channel;
  int color, ls, lw, jobNum = 0, startEvent = 0, endEvent = 200, JES, JER, fillBtagHisto;
  float normf, eqLumi, xSect, preselEff;
  string fileName;
  vector<string> vecfileNames;
  int ndatasets;
  
  bool localgridSubmission = false;
  if ( argc > 2 ) localgridSubmission = true;
  if ( argc == 2 && ((string)argv[1]).find(".xml") != std::string::npos )
  {
    cout << "One argument added for xml file, no localgrid submission" << endl;
    xmlFileName = (string)argv[1];
  }
  else if ( argc > 2 && argc < 17 )
  {
    cerr << "Too few input arguments from script. Check again." << endl;
    return 1;
  }
  
  if (localgridSubmission)
  {
    //xmlFileName = "topWidth_localgrid.xml";
    dName         = argv[1];
    dTitle        = argv[2];
    color         = strtol(argv[4], NULL, 10);
    ls            = strtol(argv[5], NULL, 10);
    lw            = strtol(argv[6], NULL, 10);
    normf         = strtod(argv[7], NULL);
    eqLumi        = strtod(argv[8], NULL);
    xSect         = strtod(argv[9], NULL);
    preselEff     = strtod(argv[10], NULL);
    fileName      = argv[11];
    // if there only two arguments after the fileName, the jobNum will be set to 0 by default as an integer is expected and it will get a string (lastfile of the list) 
    JES           = strtol(argv[argc-7], NULL,10);
    JER           = strtol(argv[argc-6], NULL,10);
    fillBtagHisto = strtol(argv[argc-5], NULL,10);
    channel       = argv[argc-4];
    jobNum        = strtol(argv[argc-3], NULL, 10);
    startEvent    = strtol(argv[argc-2], NULL, 10);
    endEvent      = strtol(argv[argc-1], NULL, 10);
    
    // all the files are stored from arg 11 to argc-7
    vecfileNames.clear();
    for(int args = 11; args < argc-7; args++) 
    {
      vecfileNames.push_back(argv[args]);
    }
    
    // Update output path according to channel
    mkdir((pathOutput+channel).c_str(),0777);
    pathOutput += channel+"/";
  }
  // Give timestamp to output path
  mkdir((pathOutput+dateString).c_str(),0777);
  pathOutput += dateString+"/";
  
  
  
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
  bool applyJESup = false;  // NOT IMPLEMENTED YET !!
  bool applyJESdown = false;
  bool calculateBTagSF = false;
  bool applyBTagSF = true;
  bool applyJetLeptonCleaning = true;
  
  if (localgridSubmission)
  {
    if ( JES == 0 ) { applyJEC = true; applyJESup = false; applyJESdown = false;}
    else if ( JES == 1 ) { applyJEC = false; applyJESup = true; applyJESdown = false;}
    else if ( JES == -1 ) { applyJEC = false; applyJESup = false; applyJESdown = true;}
    
    if ( JER == 0 ) { applyJER = true; applyJERup = false; applyJERdown = false;}
    else if ( JER == 1 ) { applyJER = false; applyJERup = true; applyJERdown = false;}
    else if ( JER == -1 ) { applyJER = false; applyJERup = false; applyJERdown = true;}
    
    //if ( fillBtagHisto == 0 ) { applyBTagSF = true; calculateBTagSF = false;}
    //else if ( fillBtagHisto == 1 ) { applyBTagSF = false; calculateBTagSF = true;}
  }
  
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
  
  if (localgridSubmission)
  {
    cout << "Using localgrid submission" << endl;
    cout << "---Dataset accepted from command line---" << endl;
    cout << "Dataset Name: " << dName << endl;
    cout << "Dataset Title: " << dTitle << endl;
    cout << "Dataset color: " << color << endl;
    cout << "Dataset ls: " << ls << endl;
    cout << "Dataset lw: " << lw << endl;
    cout << "Dataset normf: " << normf << endl;
    cout << "Dataset EqLumi: " << eqLumi << endl;
    cout << "Dataset xSect: " << xSect << endl;
    cout << "Dataset File Name: " << vecfileNames[0] << endl;
    if ( vecfileNames.size() > 1 )
    {
      for (unsigned int i = 1; i < vecfileNames.size(); i++)
      {
        cout << "                   " << vecfileNames[i] << endl;
      }
    }
    cout << "Channel is " << channel << endl;
    cout << "Beginning Event: " << startEvent << endl;
    cout << "Ending Event: " << endEvent << endl;
    cout << "JobNum: " << jobNum << endl;
    cout << "----------------------------------------" << endl << endl;
  }
  
  
  /// xml file
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
  Dataset* theDataset;
  
  cout << " - Loading datasets ..." << endl;
  if (localgridSubmission)
  {
    theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
    theDataset->SetEquivalentLuminosity(eqLumi);
    datasets.push_back(theDataset);
    //ndatasets = datasets.size() - 1;
    ndatasets = datasets.size();
  }
  else
  {
    treeLoader.LoadDatasets(datasets, xmlfile); cout << "Number of datasets: " << datasets.size() << endl;
    for (unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
    ndatasets = datasets.size();
  }
  
  float Luminosity = oldLuminosity;
  
  
  for (unsigned int d = 0; d < ndatasets; d++)
  {
    string dataSetName = datasets[d]->Name();
    if ( (dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA")!= std::string::npos) && Luminosity > datasets[d]->EquivalentLumi() ) { Luminosity = datasets[d]->EquivalentLumi();}
  }
  
  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
  
  
  //Global variable
  //TRootEvent* event = 0;
  TRootRun *runInfos = new TRootRun();
  
  //nof selected events
  //double NEvtsData = 0;
  //Double_t *nEvents = new Double_t[ndatasets];
  //Double_t nloweight = 0;
  
  
  
  ////////////////////////////
  ///  Initialise trigger  ///
  ////////////////////////////
  
  //Trigger* trigger = new Trigger(hasMuon, hasElectron, trigSingleLep, trigDoubleLep);
  Trigger* trigger = new Trigger(1, 0, 1, 0);
  
  
  
  //////////////////////////////////
  ///  Initialise scale factors  ///
  //////////////////////////////////
  
  string pathCalLept = "../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/20170413/";
  string pathCalBTag = "../TopTreeAnalysisBase/Calibrations/BTagging/";
  string pathCalPileup = "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/";
  string pathCalJEC = "../TopTreeAnalysisBase/Calibrations/JECFiles/";
  
  /// Leptons
  cout << " - Loading lepton scale factors ...";
  if (! applyLeptonSF) { cout << "    --- At the moment these are not used in the analysis";}
  cout << endl;
  
  double muonSFID, muonSFIso, muonSFTrigBCDEF, muonSFTrigGH, muonSFTrack;
  MuonSFWeight* muonSFWeightID_T_BCDEF = new MuonSFWeight(pathCalLept+"ID_EfficienciesAndSF_BCDEF.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false); // (... , ... , extendRange, debug, print warning)
  MuonSFWeight* muonSFWeightID_T_GH = new MuonSFWeight(pathCalLept+"ID_EfficienciesAndSF_GH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false);
  
  MuonSFWeight* muonSFWeightIso_TT_BCDEF = new MuonSFWeight(pathCalLept+"Iso_EfficienciesAndSF_BCDEF.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
  MuonSFWeight* muonSFWeightIso_TT_GH = new MuonSFWeight(pathCalLept+"Iso_EfficienciesAndSF_GH.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
  
  MuonSFWeight *muonSFWeightTrig_BCDEF = new MuonSFWeight(pathCalLept+"TrigEfficienciesAndSF_RunBtoF.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
  MuonSFWeight *muonSFWeightTrig_GH = new MuonSFWeight(pathCalLept+"TrigEfficienciesAndSF_GH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
  
  TFile *muontrackfile = new TFile((pathCalLept+"Tracking_EfficienciesAndSF_BCDEFGH.root").c_str(),"read");
  TGraph* h_muonSFWeightTrackEta = (TGraph*) muontrackfile->Get("ratio_eff_eta3_dr030e030_corr")->Clone();//Tracking efficiency as function of eta
  TGraph* h_muonSFWeightTrackPV = (TGraph*) muontrackfile->Get("ratio_eff_vtx_dr030e030_corr")->Clone();//Tracking efficiency as function of nPV
  
  
  /// B tag
  // documentation at http://mon.iihe.ac.be/%7Esmoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees_v4.pdf
  cout << " - Loading b tag scale factors ...";
  if (! applyBTagSF) { cout << "     --- At the moment these are not used in the analysis";}
  BTagCalibration *bTagCalib = new BTagCalibration("CSVv2", pathCalBTag+"CSVv2_76X_combToMujets.csv"); 
  BTagCalibrationReader *bTagReader_M = new BTagCalibrationReader(bTagCalib, BTagEntry::OP_MEDIUM, "mujets","central");
  BTagWeightTools *bTagHistoTool_M;
  
  
  /// Pile-up
  cout << " - Loading pile-up scale factors ...";
  if (! applyPU) { cout << "   --- At the moment these are not used in the analysis";}
  cout << endl;
  
  LumiReWeighting LumiWeights(pathCalPileup+"MCPileup_Summer16.root", pathCalPileup+"pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet.root", "pileup", "pileup");
  LumiReWeighting LumiWeights_up(pathCalPileup+"MCPileup_Summer16.root", pathCalPileup+"pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysPlus.root", "pileup", "pileup");
  LumiReWeighting LumiWeights_down(pathCalPileup+"MCPileup_Summer16.root", pathCalPileup+"pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysMinus.root", "pileup", "pileup");
  
  /// JEC
  vector<JetCorrectorParameters> vCorrParam;
  
  
  
  ///////////////////
  ///  Selection  ///
  ///////////////////
  
  float muonPTSel = 26.; // GeV
  float muonEtaSel = 2.4;
  float muonRelIsoSel = 0.15;  // Tight muon
  string muonWP = "Tight";
  
  float muonPTVeto = 10.; // GeV
  float muonEtaVeto = 2.5;
  float muonRelIsoVeto = 0.25;  // Loose muon
  
  float electronPTVeto = 15.; // GeV
  float electronEtaVeto = 2.5;
  
  float jetPT = 30.; // GeV
  float jetEta = 2.4;  // to allow b tagging
  
  
  
  //////////////////////////////////////
  ///  Working points for b tagging  ///
  //////////////////////////////////////
  
  /// Updated 13/04/17, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
  
  double CSVv2Loose =  0.5426;
  double CSVv2Medium = 0.8484;
  double CSVv2Tight = 0.9535;
  
  
  
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << ndatasets << " datasets !" << endl;
  
  for (unsigned int d = 0; d < ndatasets; d++)
  //for (unsigned int d = 1; d < 2; d++)
  {
    bool nlo = false;
    bool isData = false;
    
    string dataSetName = datasets[d]->Name();
    
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
      cout << "      -> Equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
      cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    }
    
    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA") != std::string::npos )
    {
      isData = true;
    }
    
    if ( (datasets[d]->Title()).find("amc") != std::string::npos || (datasets[d]->Title()).find("AMC") != std::string::npos || (datasets[d]->Title()).find("Amc") != std::string::npos || (datasets[d]->Title()).find("aMC") != std::string::npos )
    {
      nlo = true;
      cout << "         This is an amc@nlo sample." << endl;
    }
    
    if (calculateBTagSF && isData)
    {
      cout << "  Calculating btag scale factors.... Skipping data..." << endl;
      continue;
    }
    
    anaEnv.METCollection = "PFMET_slimmedMETs";
    if (isData) anaEnv.METCollection = "PFMET_slimmedMETsMuEGClean";
    
    //open files and load
    cout << "Load Dataset" << endl;
    treeLoader.LoadDataset(datasets[d], anaEnv);
    cout << "Load Dataset" << endl;
    
    /// book triggers
    trigger->bookTriggers(isData);
    
    
    
    /////////////////
    ///  BTag SF  ///
    /////////////////
    
    if (applyBTagSF)
    {
      /// Use seperate per data set?? (We have these...)
      //  string pathBTagHistos = BTagHistos/160729/merged/";
      //  bTagHistoTool_M = new BTagWeightTools(bTagReader_M, pathBTagHistos+"BTagSFs_"+dataSetName+"_mujets_central.root", false, 20., 600., 2.4);
      bTagHistoTool_M = new BTagWeightTools(bTagReader_M, "PlotsForBTagSFs.root", false, 20., 600., 2.4);
    }
    else if (calculateBTagSF && ! isData)
    {
      mkdir(("BTagHistos/"+dateString).c_str(),0777);
      bTagHistoTool_M = new BTagWeightTools(bTagReader_M,"BTagHistos/"+dateString+"/BTagSFs_"+dataSetName+"_"+ConvertIntToString(jobNum,0)+"_mujets_central.root", false, 20., 600., 2.4);
    }
    
    
    
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
    
    string rootFileName = "Ntuples_output_"+dataSetName+"_"+ConvertIntToString(jobNum,0)+".root";
    
    cout << " - Recreate output file ..." << rootFileName << endl;
    TFile *fout = new TFile ((pathOutput+rootFileName).c_str(), "RECREATE");
    fout->cd();
    
    TTree* myTree = new TTree("tree","tree");
    TTree* statTree = new TTree("stats","stats");
//    TTree* globalTree = new TTree("globaltree","globaltree");  // no selection applied
    
    
    
    ////////////////////////////////////
    ///  Define variables for trees  ///
    ////////////////////////////////////
    
    // stats of dataset
    Long64_t nEvents;
    Long64_t nEventsSel;
    Int_t nofPosWeights;
    Int_t nofNegWeights;
    Double_t sumW;
    
    Long64_t nofEventsHLTv2;
    Long64_t nofEventsHLTv3;
    Long64_t nofSelEventsHLTv2;
    Long64_t nofSelEventsHLTv3;
    
    // event related variables
    Int_t run_num;
    Long64_t evt_num;
    Int_t lumi_num;
    Int_t nvtx;
    Int_t npu;
    Double_t rho;
    
    Bool_t isTrigged;
    Bool_t isSelected;
    Bool_t hasExactly4Jets;
    Bool_t hasJetLeptonCleaning;
//    Bool_t passedMETFilter;
    Int_t cutFlow[10];
    
    Int_t appliedJER;
    Int_t appliedJES;
    Int_t appliedPU;
    
    Bool_t hasHLTv2;
    Bool_t hasHLTv3;
    Bool_t hasPosWeight;
    Bool_t hasNegWeight;
    
    Double_t puSF;
    Double_t puSF_up;
    Double_t puSF_down;
    Double_t btagSF;
    Double_t muonIdSF[10];
    Double_t muonIsoSF[10];
    Double_t muonTrigBCDEF[10];
    Double_t muonTrigGH[10];
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
    Bool_t mc_isLastCopy[200];
    Bool_t mc_isPromptFinalState[200];
    Bool_t mc_isHardProcess[200];
    Bool_t mc_fromHardProcessFinalState[200];
    
    
    
    /////////////////////////
    ///  Define branches  ///
    /////////////////////////
    
    
    statTree->Branch("nEvents" , &nEvents, "nEvents/L");
    statTree->Branch("nEventsSel" , &nEventsSel, "nEventsSel/L");
    statTree->Branch("cutFlow",&cutFlow,"cutFlow[10]/I");
    statTree->Branch("appliedJER",&appliedJER,"appliedJER/I");
    statTree->Branch("appliedJES", &appliedJES, "appliedJES/I");
    statTree->Branch("appliedPU", &appliedPU, "appliedPU/I");
    if (isData)
    {
      statTree->Branch("nofEventsHLTv2",&nofEventsHLTv2,"nofEventsHLTv2/L");
      statTree->Branch("nofEventsHLTv3",&nofEventsHLTv3,"nofEventsHLTv3/L");
      statTree->Branch("nofSelEventsHLTv2",&nofSelEventsHLTv2,"nofSelEventsHLTv2/L");
      statTree->Branch("nofSelEventsHLTv3",&nofSelEventsHLTv3,"nofSelEventsHLTv3/L");
    }
    if (nlo)
    {
      statTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
      statTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
      statTree->Branch("sumW", &sumW, "sumW/D");
      statTree->Branch("hasPosWeight",&hasPosWeight,"hasPosWeight/O");
      statTree->Branch("hasNegWeight",&hasNegWeight,"hasNegWeight/O");
    }
    
//    globalTree->Branch("run_num",&run_num,"run_num/I");
//    globalTree->Branch("evt_num",&evt_num,"evt_num/L");
//    globalTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
//    globalTree->Branch("nvtx",&nvtx,"nvtx/I");
//    globalTree->Branch("npu",&npu,"npu/I");
//    globalTree->Branch("rho",&rho,"rho/D");
//    globalTree->Branch("isTrigged",&isTrigged,"isTrigged/O");
//    globalTree->Branch("hasJetLeptonCleaning",&hasJetLeptonCleaning,"hasJetLeptonCleaning/O");
//    globalTree->Branch("appliedJER",&appliedJER,"appliedJER/I");
//    globalTree->Branch("appliedJES", &appliedJES, "appliedJES/I");
//    globalTree->Branch("appliedPU", &appliedPU, "appliedPU/I");
    
    myTree->Branch("run_num",&run_num,"run_num/I");
    myTree->Branch("evt_num",&evt_num,"evt_num/L");
    myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
    myTree->Branch("nvtx",&nvtx,"nvtx/I");
    myTree->Branch("npu",&npu,"npu/I");
    myTree->Branch("rho",&rho,"rho/D");
    myTree->Branch("isTrigged",&isTrigged,"isTrigged/O");
    myTree->Branch("hasExactly4Jets",&hasExactly4Jets,"hasExactly4Jets/O");
    myTree->Branch("hasJetLeptonCleaning",&hasJetLeptonCleaning,"hasJetLeptonCleaning/O");
//    myTree->Branch("passedMETFilter", &passedMETFilter,"passedMETFilter/O");
    if (isData)
    {
      myTree->Branch("hasHLTv2",&hasHLTv2,"hasHLTv2/O");
      myTree->Branch("hasHLTv3",&hasHLTv3,"hasHLTv3/O");
    }
    
    
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
//      globalTree->Branch("mc_isLastCopy", &mc_isLastCopy, "mc_isLastCopy[nMCParticles]/O");
//      globalTree->Branch("mc_isPromptFinalState", &mc_isPromptFinalState, "mc_isPromptFinalState[nMCParticles]/O");
//      globalTree->Branch("mc_isHardProcess", &mc_isHardProcess, "mc_isHardProcess[nMCParticles]/O");
//      globalTree->Branch("mc_fromHardProcessFinalState", &mc_fromHardProcessFinalState, "mc_fromHardProcessFinalState[nMCParticles]/O");
      
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
      myTree->Branch("mc_isLastCopy", &mc_isLastCopy, "mc_isLastCopy[nMCParticles]/O");
      myTree->Branch("mc_isPromptFinalState", &mc_isPromptFinalState, "mc_isPromptFinalState[nMCParticles]/O");
      myTree->Branch("mc_isHardProcess", &mc_isHardProcess, "mc_isHardProcess[nMCParticles]/O");
      myTree->Branch("mc_fromHardProcessFinalState", &mc_fromHardProcessFinalState, "mc_fromHardProcessFinalState[nMCParticles]/O");
    }
    
    
    /// SFs
//    globalTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
//    globalTree->Branch("puSF",&puSF,"puSF/D");
//    globalTree->Branch("btagSF",&btagSF,"btagSF/D");
//    globalTree->Branch("muonIdSF",&muonIdSF,"muonIdSF[nMuons]/D");
//    globalTree->Branch("muonIsoSF",&muonIsoSF, "muonIsoSF[nMuons]/D");
//    globalTree->Branch("muonTrigSFv2",&muonTrigSFv2,"muonTrigSFv2[nMuons]/D");
//    globalTree->Branch("muonTrigSFv3",&muonTrigSFv3,"muonTrigSFv3[nMuons]/D");
//    globalTree->Branch("electronSF",&electronSF,"electronSF[nElectrons]/D");
    
    if (! isData)
    {
      myTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
      myTree->Branch("puSF",&puSF,"puSF/D");
      myTree->Branch("btagSF",&btagSF,"btagSF/D");
      myTree->Branch("muonIdSF",&muonIdSF,"muonIdSF[nMuons]/D");
      myTree->Branch("muonIsoSF",&muonIsoSF, "muonIsoSF[nMuons]/D");
      //myTree->Branch("muonTrigSFv2",&muonTrigSFv2,"muonTrigSFv2[nMuons]/D");
      //myTree->Branch("muonTrigSFv3",&muonTrigSFv3,"muonTrigSFv3[nMuons]/D");
//      myTree->Branch("electronSF",&electronSF,"electronSF[nElectrons]/D");
    }
    
    
    
    ////////////////////////////
    ///  Determine range of events to run over
    ////////////////////////////
    
    
    unsigned int ending = datasets[d]->NofEvtsToRunOver();
    double end_d = ending;
    if ( localgridSubmission && endEvent < ending )
      end_d = endEvent;
    //end_d = 2000;  // for testing
    
    if ( end_d < startEvent )
    {
	    cout << "Starting event larger than number of events. Exiting..." << endl;
	    exit(1);
    }
    if ( verbose > 1 )
    {
      cout << "Number of events in total dataset = " << ending << endl;
      cout << "Will run over " << (end_d - startEvent) << " events..." << endl;
      cout << "Starting event = = = = " << startEvent  << endl;
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
    nofSelEventsHLTv2 = 0;
    nofSelEventsHLTv3 = 0;
    for (Int_t i = 0; i < 10; i++)
      cutFlow[i] = 0;
    
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
    vector < TRootPFJet* > selectedJetsBC;
    vector < TRootMuon* > selectedMuons;
    vector < TRootMuon* > vetoMuons;
    vector < TRootElectron* > selectedElectrons;
    vector < TRootElectron* > vetoElectrons;
    
    /// Systematics
    if (applyJEC) { appliedJES = 0;}
    else if (applyJESup) { appliedJES = 1;}
    else if (applyJESdown) { appliedJES = -1;}
    if (applyJER) { appliedJER = 0;}
    else if (applyJERup) { appliedJER = 1;}
    else if (applyJERdown) { appliedJER = -1;}
    if (applyPU) { appliedPU = 0;}
    else if (applyPUup) { appliedPU = 1;}
    else if (applyPUdown) { appliedPU = -1;}
    
    
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
    for (unsigned int ievt = startEvent; ievt < end_d; ievt++)
    //for (unsigned int ievt = 0; ievt < 10000000; ievt++)
    {
      nEvents++;
      
      if (ievt%100000 == 0)
        cout << "Processing event " << ievt << "..." << endl;
      
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
      selectedJetsBC.clear();
      selectedMuons.clear();
      vetoMuons.clear();
      selectedElectrons.clear();
      vetoElectrons.clear();
      
      /// Reset other stuff
      isTrigged = false;
      isSelected = false;
      hasExactly4Jets = false;
      hasJetLeptonCleaning = false;
      //passedMETFilter = false;
      puSF = 1.;
      btagSF = 1.;
      for (Int_t i = 0; i < 10; i++)
      {
        muonIdSF[i] = 1.;
        muonIsoSF[i] = 1.;
        //muonTrigSFv2[i] = 1.;
        //muonTrigSFv3[i] = 1.;
        //electronSF[i] = 1.;
      }
      nloWeight = 1.; // for amc@nlo samples
      
      hasHLTv2 = false;
      hasHLTv3 = false;
      hasPosWeight = false;
      hasNegWeight = false;
    
      
      nLeptons = -1;
      nElectrons = -1;
      nMuons = -1;
      nJets = -1;
      
      for (Int_t i = 0; i < 10; i++)
      {
//        electron_charge[i] = 0;
//        electron_pt[i] = 0.;
//        electron_phi[i] = 0.;
//        electron_eta[i] = 0.;
//        electron_eta_superCluster[i] = 0.;
//        electron_E[i] = 0.;
//        electron_M[i] = 0.;
//        electron_d0[i] = -1.;
//        electron_chargedHadronIso[i] = -1.;
//        electron_neutralHadronIso[i] = -1.;
//        electron_photonIso[i] = -1.;
//        electron_pfIso[i] = -1.;
//        electron_sigmaIEtaIEta[i] = -1.;
//        electron_deltaEtaIn[i] = -1.;
//        electron_deltaPhiIn[i] = -1.;
//        electron_hadronicOverEm[i] = -1.;
//        electron_missingHits[i] = -1;
//        electron_passConversion[i] = 0;
//        electron_isEBEEGap[i] = 0;
        
        muon_charge[i] = 0;
        muon_pt[i] = 0.;
        muon_phi[i] = 0.;
        muon_eta[i] = 0.;
        muon_E[i] = 0.;
        muon_M[i] = 0.;
        muon_d0[i] = 999.;
        muon_chargedHadronIso[i] = 999.;
        muon_neutralHadronIso[i] = 999.;
        muon_photonIso[i] = 999.;
        muon_puChargedHadronIso[i] = 999.;
        muon_relIso[i] = 999.;
        muon_pfIso[i] = 999.;
      }
      
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
      
      /// mcparticles
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
        mc_isLastCopy[i] = false;
        mc_isPromptFinalState[i] = false;
        mc_isHardProcess[i] = false;
        mc_fromHardProcessFinalState[i] = false;
      }
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets, mets);
      init_jets_corrected = init_jets;
      
      datasets[d]->eventTree()->LoadTree(ievt);
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
          cerr << "   Run number is " << run_num << endl;
          //cerr << "   File name is " << datasets[d]->eventTree()->GetFile()->GetName() << endl;
          exit(1);
        }
        else if ( run_num >= 256630 && run_num <= 257819 )
        {
          nofEventsHLTv2++;
          hasHLTv2 = true;
        }
        else
        {
          nofEventsHLTv3++;
          hasHLTv3 = true;
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
            hasPosWeight = true;
          }
          else
          {
            nofNegWeights++;
            hasNegWeight = true;
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
            hasPosWeight = true;
          }
          else
          {
            nofNegWeights++;
            hasNegWeight = true;
          }
        }
        
        sumW += nloWeight;
      }
      
      jetTools->correctJets(init_jets_corrected, rho, isData);
      
      if (! isData)
      {
        puSF = LumiWeights.ITweight( (int)event->nTruePU() );
        puSF_up = LumiWeights_up.ITweight( (int)event->nTruePU() );
        puSF_down = LumiWeights_down.ITweight( (int)event->nTruePU() );
        
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
      vetoMuons = selection.GetSelectedMuons(muonPTVeto, muonEtaVeto, muonRelIsoVeto, "Loose", "Spring15");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      //vetoElectrons = selection.GetSelectedElectrons(electronPTVeto, electronEtaVeto, "Veto", "Spring15_25ns", true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased
      
      if (applyJetLeptonCleaning)
      {
        if(verbose > 3) cout << "  - Applying jet/lepton cleaning... " << endl; 
        
        selectedJetsBC.clear();
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
              hasJetLeptonCleaning = true;
              break;
            }
          }
          if (toBeErased) continue;
          for (int iElectron = 0; iElectron < selectedElectrons.size(); iElectron++)
          {
            if ( selectedJetsBC[iOrigJet]->DeltaR(*selectedElectrons[iElectron]) < 0.3 )
            {
              toBeErased = true;
              hasJetLeptonCleaning = true;
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
        muon_photonIso[iMuon] = selectedMuons[iMuon]->photonIso(4);
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
        if (nMCParticles > maxMCParticles) maxMCParticles = nMCParticles;
        for (Int_t iMC = 0; iMC < nMCParticles; iMC++)
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
          mc_isLastCopy[iMC] = mcParticles[iMC]->isLastCopy();
          mc_isPromptFinalState[iMC] = mcParticles[iMC]->isPromptFinalState();
          mc_isHardProcess[iMC] = mcParticles[iMC]->isHardProcess();
          mc_fromHardProcessFinalState[iMC] = mcParticles[iMC]->fromHardProcessFinalState();
        }
      }
      
      
      /// Fill scalefactors
      if (! isData)
      {
        if (applyBTagSF) btagSF = bTagHistoTool_M->getMCEventWeight(selectedJets);
        
        for (int iMuon = 0; iMuon < selectedMuons.size(); iMuon++)
        {
          //muonIdSF[iMuon] = muonSFWeightID_T->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);  // eta, pt, shiftUpDown;
          //muonIsoSF[iMuon] = muonSFWeightIso_TT->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
          //muonTrigSFv2[iMuon] = muonSFWeightTrigHLTv4p2->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
          //muonTrigSFv3[iMuon] = muonSFWeightTrigHLTv4p3->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
        }
        
        //for (int iElectron = 0; iElectron < selectedElectrons(); iElectron++)
        //{
        //  electronSF[iElectron] = ;
        //}
      }
      
//      globalTree->Fill();
      
      
      ////// Selection
      cutFlow[0]++;
      if (isTrigged)
      {
        cutFlow[1]++;
        if (isGoodPV)
        {
          cutFlow[2]++;
          if (selectedMuons.size() == 1)
          {
            cutFlow[3]++;
            if (vetoMuons.size() == 1)
            {
              cutFlow[4]++;
              if (vetoElectrons.size() == 0)
              {
                cutFlow[5]++;
                
//                 /// First 4 jets need pT > 30 GeV
//                 if (selectedJets.size() >= 4)
//                 {
//                   if (selectedJets[3]->Pt() < 30) selectedJets.clear();
//                 }
                
                if ( selectedJets.size() >= 4 )
                {
                  cutFlow[6]++;
                  if ( selectedJets.size() == 4 ) hasExactly4Jets = true;
                  
                  if ( selectedBJets.size() > 0 )
                  {
                    cutFlow[7]++;
                    if ( selectedBJets.size() > 1 )
                    {
                      cutFlow[8]++;
                      isSelected = true;
                    }  // at least 2 b-tagged jets
                  }  // at least 1 b-tagged jet
                  
                }  // at least 4 jets
              }  // no veto electrons
            }  // no additional loose muons (tight muon is also loose muon)
          }  // 1 good muon
        }  // good PV
      }  // trigged
      
      
      if (! isSelected)
      {
        continue;
      }
      
      nEventsSel++;
      
      if (isData)
      {
        if ( run_num >= 256630 && run_num <= 257819 )
        {
          nofSelEventsHLTv2++;
        }
        else
        {
          nofSelEventsHLTv3++;
        }
      }
      
      if (!calculateBTagSF)
      {
        myTree->Fill();
      }
      
      
      /// B-tagging
      if (calculateBTagSF && ! isData)
      {
        bTagHistoTool_M->FillMCEfficiencyHistos(selectedJets);
      }
      
      
      
      
      
      
      
      
      
    }  // end loop events
    cout << "Max MCParticles: " << maxMCParticles << endl << endl;
    
    cout << "Fill trees..." << endl;
    if (! calculateBTagSF)
    {
      statTree->Fill();
      
      /// Write to file
      fout->cd();
      myTree->Write();
      statTree->Write();
    }
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
