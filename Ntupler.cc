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
#include "TFile.h"
#include "TTree.h"
#include <errno.h>
#include "TRandom3.h"
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
#include "Tools/interface/Trigger.h"


using namespace std;
using namespace reweight;
using namespace TopTree;


map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;


bool test = false;
bool fillLooseTree = false;


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
  
  string xmlFileName ="config/topWidth_data.xml";
//  string xmlFileName ="config/topWidth_MC.xml";
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
    test = false;
    
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
  bool applyJESup = false;  // Check implementation
  bool applyJESdown = false;
  bool calculateBTagSF = false;
  bool applyBTagSF = true;
  bool applyJetLeptonCleaning = true;
  
  if (localgridSubmission)
  {
    if ( JES == 0 ) { applyJEC = true; /*applyJESup = false; applyJESdown = false;*/}  // Check
    else if ( JES == 1 ) { /*applyJEC = false;*/ applyJESup = true; applyJESdown = false;}
    else if ( JES == -1 ) { /*applyJEC = false;*/ applyJESup = false; applyJESdown = true;}
    
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
  if (! fillLooseTree)  cout << "* Not filling tree with loose objects       *" << endl;
  else                  cout << "* Filling loose objects tree for TT nominal *" << endl;
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
  TRootEvent* event = 0;
  TRootRun *runInfos = new TRootRun();
  
  
  
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
  
  MuonSFWeight* muonSFWeightID_T_BCDEF = new MuonSFWeight(pathCalLept+"IDEfficienciesAndSF_BCDEF.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false); // (... , ... , extendRange, debug, print warning)
  MuonSFWeight* muonSFWeightID_T_GH = new MuonSFWeight(pathCalLept+"IDEfficienciesAndSF_GH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false);
  
  MuonSFWeight* muonSFWeightIso_TT_BCDEF = new MuonSFWeight(pathCalLept+"IsoEfficienciesAndSF_BCDEF.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
  MuonSFWeight* muonSFWeightIso_TT_GH = new MuonSFWeight(pathCalLept+"IsoEfficienciesAndSF_GH.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
  
  MuonSFWeight *muonSFWeightTrig_BCDEF = new MuonSFWeight(pathCalLept+"TrigEfficienciesAndSF_RunBtoF.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
  MuonSFWeight *muonSFWeightTrig_GH = new MuonSFWeight(pathCalLept+"TrigEfficienciesAndSF_GH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
  
  TFile *muontrackfile = new TFile((pathCalLept+"Tracking_EfficienciesAndSF_BCDEFGH.root").c_str(),"read");
  TGraph* h_muonSFWeightTrackEta = (TGraph*) muontrackfile->Get("ratio_eff_eta3_dr030e030_corr")->Clone();//Tracking efficiency as function of eta
  TGraph* h_muonSFWeightTrackAEta = (TGraph*) muontrackfile->Get("ratio_eff_aeta_dr030e030_corr")->Clone();//Tracking efficiency as function of abs(eta)
  TGraph* h_muonSFWeightTrackPV = (TGraph*) muontrackfile->Get("ratio_eff_vtx_dr030e030_corr")->Clone();//Tracking efficiency as function of nPV
  
  
  /// B tag
  // documentation at http://mon.iihe.ac.be/%7Esmoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees_v4.pdf
  cout << " - Loading b tag scale factors ...";
  if (! applyBTagSF) { cout << "     --- At the moment these are not used in the analysis";}
  BTagCalibration *bTagCalib = new BTagCalibration("CSVv2", pathCalBTag+"CSVv2Moriond17_2017_1_26_BtoH.csv");
  BTagCalibrationReader *bTagReader_M = new BTagCalibrationReader(bTagCalib, BTagEntry::OP_MEDIUM, "comb", "central");
  BTagCalibrationReader *bTagReader_M_up = new BTagCalibrationReader(bTagCalib, BTagEntry::OP_MEDIUM, "comb", "up");
  BTagCalibrationReader *bTagReader_M_down = new BTagCalibrationReader(bTagCalib, BTagEntry::OP_MEDIUM, "comb", "down");
  BTagWeightTools *bTagHistoTool_M;
  BTagWeightTools *bTagHistoTool_M_up;
  BTagWeightTools *bTagHistoTool_M_down;
  
  
  /// Pile-up
  cout << " - Loading pile-up scale factors ...";
  if (! applyPU) { cout << "   --- At the moment these are not used in the analysis";}
  cout << endl;
  
  LumiReWeighting LumiWeights(pathCalPileup+"MCPileup_Summer16.root", pathCalPileup+"pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet.root", "pileup", "pileup");
  LumiReWeighting LumiWeights_up(pathCalPileup+"MCPileup_Summer16.root", pathCalPileup+"pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysPlus.root", "pileup", "pileup");
  LumiReWeighting LumiWeights_down(pathCalPileup+"MCPileup_Summer16.root", pathCalPileup+"pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysMinus.root", "pileup", "pileup");
  
  /// JEC
  vector<JetCorrectorParameters> vCorrParam;
  JetCorrectionUncertainty *jecUnc;
  
  
  
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
  
  float electronPTSel = 34.; // GeV
  float electronEtaSel = 2.1;
  string electronWP = "Tight";
  
  float electronPTVeto = 15.; // GeV
  float electronEtaVeto = 2.5;
  
  float jetPT = 30.; // GeV
  float jetEta = 2.4;  // to allow b tagging
  
  
  
  //////////////////////////////////////
  ///  Working points for b tagging  ///
  //////////////////////////////////////
  
  /// Updated 13/04/17, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
  
  double CSVv2Loose  = 0.5426;
  double CSVv2Medium = 0.8484;
  double CSVv2Tight  = 0.9535;
  
  
  
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
    bool isTTbar = false;
    
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
    else if ( dataSetName.find("TT") != std::string::npos )
    {
      isTTbar = true;
    }
    
    if (! isTTbar || dataSetName.find("nominal") == std::string::npos )
    {
      fillLooseTree = false;
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
    //if (isData) anaEnv.METCollection = "PFMET_slimmedMETsMuEGClean";
    
    //open files and load
    cout << "Load Dataset" << endl;
    treeLoader.LoadDataset(datasets[d], anaEnv);
    cout << "Load Dataset" << endl;
    
    /// book triggers
    trigger->bookTriggers(false);  // only for MuEG triggers: isDataRunH or not
    
    
    
    /////////////////
    ///  BTag SF  ///
    /////////////////
    
    if (applyBTagSF)
    {
      /// Use seperate per data set?? (We have these...)
      //  string pathBTagHistos = BTagHistos/160729/Merged/";
      //  bTagHistoTool_M = new BTagWeightTools(bTagReader_M, pathBTagHistos+"BTagSFs_"+dataSetName+"_comb_central.root", false, 20., 600., 2.4);
      bTagHistoTool_M = new BTagWeightTools(bTagReader_M, "PlotsForBTagSFs_nominal.root", false, 30., 600., 2.4);
      bTagHistoTool_M_up = new BTagWeightTools(bTagReader_M_up, "PlotsForBTagSFs_up.root", false, 30., 600., 2.4);
      bTagHistoTool_M_down = new BTagWeightTools(bTagReader_M_down, "PlotsForBTagSFs_down.root", false, 30., 600., 2.4);
    }
    else if (calculateBTagSF && ! isData)
    {
      mkdir(("BTagHistos/"+dateString).c_str(),0777);
      bTagHistoTool_M = new BTagWeightTools(bTagReader_M, "BTagHistos/"+dateString+"/BTagSFs_"+dataSetName+"_"+ConvertIntToString(jobNum,0)+"_comb_central.root", false, 30., 600., 2.4);
      bTagHistoTool_M_up = new BTagWeightTools(bTagReader_M_up, "BTagHistos/"+dateString+"/BTagSFs_"+dataSetName+"_"+ConvertIntToString(jobNum,0)+"_comb_up.root", false, 30., 600., 2.4);
      bTagHistoTool_M_down = new BTagWeightTools(bTagReader_M_down, "BTagHistos/"+dateString+"/BTagSFs_"+dataSetName+"_"+ConvertIntToString(jobNum,0)+"_comb_down.root", false, 30., 600., 2.4);
    }
    
    
    
    ///////////////////////////////////////////
    ///  Initialise Jet Energy Corrections  ///
    ///////////////////////////////////////////
    
    vCorrParam.clear();
    
    if (! isData)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
    }
    else if ( dataSetName.find("Run2016B") != std::string::npos || dataSetName.find("Run2016C") != std::string::npos || dataSetName.find("Run2016D") != std::string::npos )
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if ( dataSetName.find("Run2016E") != std::string::npos || dataSetName.find("Run2016F") != std::string::npos )
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if ( dataSetName.find("Run2016G") != std::string::npos )
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if ( dataSetName.find("Run2016H") != std::string::npos )
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt");
    }
    
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
    TTree* looseTree = new TTree("looseObj","looseObj");
    
    
    
    ////////////////////////////////////
    ///  Define variables for trees  ///
    ////////////////////////////////////
    
    /// stats of dataset
    Long64_t nEvents;
    Long64_t nEventsSel;
    Int_t nofPosWeights;
    Int_t nofNegWeights;
    Double_t sumW;
    Int_t cutFlow[10];
    Int_t cutFlow2[10];
    Int_t appliedJER;
    Int_t appliedJES;
    Int_t appliedPU;
    
    Long64_t nofEventsRunB;
    Long64_t nofEventsRunCD;
    Long64_t nofEventsRunEF;
    Long64_t nofEventsRunG;
    Long64_t nofEventsRunH;
    Long64_t nofSelEventsRunB;
    Long64_t nofSelEventsRunCD;
    Long64_t nofSelEventsRunEF;
    Long64_t nofSelEventsRunG;
    Long64_t nofSelEventsRunH;
    
    Long64_t nofEventsWithGenTop;
    Long64_t nofEventsWithGenTopWithStatus22or62;
    Long64_t nofEventsWithGenAntiTop;
    Long64_t nofEventsWithGenAntiTopWithStatus22or62;
    Long64_t nofTTEventsWithoutBothGenTops;
    Long64_t nofTTEventsWithoutGenTop;
    Long64_t nofTTEventsWithoutGenAntiTop;
    Long64_t nofTTEventsWithoutBothGenTopsWithStatus22;
    Long64_t nofTTEventsWithoutGenTopWithStatus22;
    Long64_t nofTTEventsWithoutGenAntiTopWithStatus22;
    Long64_t nofTTEventsWithoutBothGenTopsWithStatus62;
    Long64_t nofTTEventsWithoutGenTopWithStatus62;
    Long64_t nofTTEventsWithoutGenAntiTopWithStatus62;
    
    /// event related variables
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
    Bool_t hasLooseJetLeptonCleaning;
    Bool_t hasErasedBadOrCloneMuon;
    Bool_t hasErasedBadOrCloneLooseMuon;
    
    // filters
    Bool_t filterHBHENoise;
    Bool_t filterHBHEIso;
    Bool_t filterCSCTightHalo;
    Bool_t filterEcalDeadCell;
    Bool_t filterEEBadSc;  // recommended for data-only
    Bool_t filterBadChCand;
    Bool_t filterBadMuon;
    Bool_t passedMETFilter;
    
    Bool_t isDataRunB;
    Bool_t isDataRunC;
    Bool_t isDataRunD;
    Bool_t isDataRunE;
    Bool_t isDataRunF;
    Bool_t isDataRunG;
    Bool_t isDataRunH;
    Bool_t hasPosWeight;
    Bool_t hasNegWeight;
    Bool_t hasGenTop;
    Bool_t hasGenTopWithStatus22;
    Bool_t hasGenTopWithStatus62;
    Bool_t hasGenAntiTop;
    Bool_t hasGenAntiTopWithStatus22;
    Bool_t hasGenAntiTopWithStatus62;
    
    Double_t nloWeight; // for amc@nlo samples
    Double_t btagSF;
    Double_t btagSF_up;
    Double_t btagSF_down;
    Double_t puSF;
    Double_t puSF_up;
    Double_t puSF_down;
    Double_t muonIdSF_BCDEF[10];
    Double_t muonIdSF_GH[10];
    Double_t muonIdSF_up_BCDEF[10];
    Double_t muonIdSF_up_GH[10];
    Double_t muonIdSF_down_BCDEF[10];
    Double_t muonIdSF_down_GH[10];
    Double_t muonIsoSF_BCDEF[10];
    Double_t muonIsoSF_GH[10];
    Double_t muonIsoSF_up_BCDEF[10];
    Double_t muonIsoSF_up_GH[10];
    Double_t muonIsoSF_down_BCDEF[10];
    Double_t muonIsoSF_down_GH[10];
    Double_t muonTrigSF_BCDEF[10];
    Double_t muonTrigSF_GH[10];
    Double_t muonTrigSF_up_BCDEF[10];
    Double_t muonTrigSF_up_GH[10];
    Double_t muonTrigSF_down_BCDEF[10];
    Double_t muonTrigSF_down_GH[10];
    Double_t muonTrackSF_eta[10];
    Double_t muonTrackSF_aeta[10];
    Double_t muonTrackSF_nPV[10];
    
    
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
    
    Int_t nLooseMuons;
    Int_t muon_loose_charge[20];
    Double_t muon_loose_pt[20];
    Double_t muon_loose_phi[20];
    Double_t muon_loose_eta[20];
    Double_t muon_loose_E[20];
    Double_t muon_loose_d0[20];
    Double_t muon_loose_relIso[20];
    Bool_t isGlobalLooseMuon[20];
    Bool_t isTrackerLooseMuon[20];
    
    /// Variables for jets
    Int_t nJets;
    Int_t jet_nConstituents[20];
    Int_t jet_nChConstituents[20];
    Int_t jet_charge[20];
    Double_t jet_pt[20];
    Double_t jet_phi[20];
    Double_t jet_eta[20];
    Double_t jet_E[20];
    Double_t jet_M[20];
    Double_t jet_bdiscr[20];
    
    Int_t nLooseJets;
    Int_t jet_loose_nConstituents[30];
    Int_t jet_loose_nChConstituents[30];
    Int_t jet_loose_charge[30];
    Double_t jet_loose_pt[30];
    Double_t jet_loose_phi[30];
    Double_t jet_loose_eta[30];
    Double_t jet_loose_E[30];
    Double_t jet_loose_M[30];
    Double_t jet_loose_bdiscr[30];
    
    /// met
    Double_t met_px;
    Double_t met_py;
    Double_t met_pt;
    Double_t met_phi;
    Double_t met_eta;
    Double_t met_Et;
    Double_t met_E;
    
    Double_t met_corr_px;
    Double_t met_corr_py;
    Double_t met_corr_pt;
    Double_t met_corr_phi;
    Double_t met_corr_eta;
    Double_t met_corr_Et;
    Double_t met_corr_E;
    
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
    statTree->Branch("cutFlow2",&cutFlow2,"cutFlow2[10]/I");
    statTree->Branch("appliedJER",&appliedJER,"appliedJER/I");
    statTree->Branch("appliedJES", &appliedJES, "appliedJES/I");
    statTree->Branch("appliedPU", &appliedPU, "appliedPU/I");
    if (isData)
    {
      statTree->Branch("nofEventsRunB",&nofEventsRunB,"nofEventsRunB/L");
      statTree->Branch("nofEventsRunCD",&nofEventsRunCD,"nofEventsRunCD/L");
      statTree->Branch("nofEventsRunEF",&nofEventsRunEF,"nofEventsRunEF/L");
      statTree->Branch("nofEventsRunG",&nofEventsRunG,"nofEventsRunG/L");
      statTree->Branch("nofEventsRunH",&nofEventsRunH,"nofEventsRunH/L");
      statTree->Branch("nofSelEventsRunB",&nofSelEventsRunB,"nofSelEventsRunB/L");
      statTree->Branch("nofSelEventsRunCD",&nofSelEventsRunCD,"nofSelEventsRunCD/L");
      statTree->Branch("nofSelEventsRunEF",&nofSelEventsRunEF,"nofSelEventsRunEF/L");
      statTree->Branch("nofSelEventsRunG",&nofSelEventsRunG,"nofSelEventsRunG/L");
      statTree->Branch("nofSelEventsRunH",&nofSelEventsRunH,"nofSelEventsRunH/L");
    }
    else
    {
      statTree->Branch("nofEventsWithGenTop", &nofEventsWithGenTop, "nofEventsWithGenTop/L");
      statTree->Branch("nofEventsWithGenTopWithStatus22or62", &nofEventsWithGenTopWithStatus22or62, "nofEventsWithGenTopWithStatus22or62/L");
      statTree->Branch("nofEventsWithGenAntiTop", &nofEventsWithGenAntiTop, "nofEventsWithGenAntiTop/L");
      statTree->Branch("nofEventsWithGenAntiTopWithStatus22or62", &nofEventsWithGenAntiTopWithStatus22or62, "nofEventsWithGenAntiTopWithStatus22or62/L");
      if (isTTbar)
      {
        statTree->Branch("nofTTEventsWithoutBothGenTops", &nofTTEventsWithoutBothGenTops, "nofTTEventsWithoutBothGenTops/L");
        statTree->Branch("nofTTEventsWithoutGenTop", &nofTTEventsWithoutGenTop, "nofTTEventsWithoutGenTop/L");
        statTree->Branch("nofTTEventsWithoutGenAntiTop", &nofTTEventsWithoutGenAntiTop, "nofTTEventsWithoutGenAntiTop/L");
        statTree->Branch("nofTTEventsWithoutBothGenTopsWithStatus22", &nofTTEventsWithoutBothGenTopsWithStatus22, "nofTTEventsWithoutBothGenTopsWithStatus22/L");
        statTree->Branch("nofTTEventsWithoutGenTopWithStatus22", &nofTTEventsWithoutGenTopWithStatus22, "nofTTEventsWithoutGenTopWithStatus22/L");
        statTree->Branch("nofTTEventsWithoutGenAntiTopWithStatus22", &nofTTEventsWithoutGenAntiTopWithStatus22, "nofTTEventsWithoutGenAntiTopWithStatus22/L");
        statTree->Branch("nofTTEventsWithoutBothGenTopsWithStatus62", &nofTTEventsWithoutBothGenTopsWithStatus62, "nofTTEventsWithoutBothGenTopsWithStatus62/L");
        statTree->Branch("nofTTEventsWithoutGenTopWithStatus62", &nofTTEventsWithoutGenTopWithStatus62, "nofTTEventsWithoutGenTopWithStatus62/L");
        statTree->Branch("nofTTEventsWithoutGenAntiTopWithStatus62", &nofTTEventsWithoutGenAntiTopWithStatus62, "nofTTEventsWithoutGenAntiTopWithStatus62/L");
      }
    }
    if (nlo)
    {
      statTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
      statTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
      statTree->Branch("sumW", &sumW, "sumW/D");
      myTree->Branch("hasPosWeight",&hasPosWeight,"hasPosWeight/O");
      myTree->Branch("hasNegWeight",&hasNegWeight,"hasNegWeight/O");
    }
    
    myTree->Branch("run_num",&run_num,"run_num/I");
    myTree->Branch("evt_num",&evt_num,"evt_num/L");
    myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
    myTree->Branch("nvtx",&nvtx,"nvtx/I");
    myTree->Branch("npu",&npu,"npu/I");
    myTree->Branch("rho",&rho,"rho/D");
    myTree->Branch("isTrigged",&isTrigged,"isTrigged/O");
    myTree->Branch("hasExactly4Jets",&hasExactly4Jets,"hasExactly4Jets/O");
    myTree->Branch("hasJetLeptonCleaning",&hasJetLeptonCleaning,"hasJetLeptonCleaning/O");
    myTree->Branch("hasErasedBadOrCloneMuon",&hasErasedBadOrCloneMuon,"hasErasedBadOrCloneMuon/O");
    if (fillLooseTree)
    {
      looseTree->Branch("hasJetLeptonCleaning",&hasJetLeptonCleaning,"hasJetLeptonCleaning/O");
      looseTree->Branch("hasLooseJetLeptonCleaning",&hasLooseJetLeptonCleaning,"hasLooseJetLeptonCleaning/O");
      looseTree->Branch("hasErasedBadOrCloneMuon",&hasErasedBadOrCloneMuon,"hasErasedBadOrCloneMuon/O");
      looseTree->Branch("hasErasedBadOrCloneLooseMuon",&hasErasedBadOrCloneLooseMuon,"hasErasedBadOrCloneLooseMuon/O");
    }
    
    myTree->Branch("filterHBHENoise",&filterHBHENoise,"filterHBHENoise/O");
    myTree->Branch("filterHBHEIso",&filterHBHEIso,"filterHBHEIso/O");
    myTree->Branch("filterCSCTightHalo",&filterCSCTightHalo,"filterCSCTightHalo/O");
    myTree->Branch("filterEcalDeadCell",&filterEcalDeadCell,"filterEcalDeadCell/O");
    myTree->Branch("filterEEBadSc",&filterEEBadSc,"filterEEBadSc/O");  // recommended for data-only
    myTree->Branch("filterBadChCand",&filterBadChCand,"filterBadChCand/O");
    myTree->Branch("filterBadMuon",&filterBadMuon,"filterBadMuon/O");
    myTree->Branch("passedMETFilter", &passedMETFilter,"passedMETFilter/O");
    
    if (isData)
    {
      myTree->Branch("isDataRunB",&isDataRunB,"isDataRunB/O");
      myTree->Branch("isDataRunC",&isDataRunC,"isDataRunC/O");
      myTree->Branch("isDataRunD",&isDataRunD,"isDataRunD/O");
      myTree->Branch("isDataRunE",&isDataRunE,"isDataRunE/O");
      myTree->Branch("isDataRunF",&isDataRunF,"isDataRunF/O");
      myTree->Branch("isDataRunG",&isDataRunG,"isDataRunG/O");
      myTree->Branch("isDataRunH",&isDataRunH,"isDataRunH/O");
    }
    
    // muons
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
    
    if (fillLooseTree)
    {
      looseTree->Branch("nMuons",&nMuons, "nMuons/I");
      looseTree->Branch("muon_charge",&muon_charge,"muon_charge[nMuons]/I");
      looseTree->Branch("muon_pt",&muon_pt,"muon_pt[nMuons]/D");
      looseTree->Branch("muon_phi",&muon_phi,"muon_phi[nMuons]/D");
      looseTree->Branch("muon_eta",&muon_eta,"muon_eta[nMuons]/D");
      looseTree->Branch("muon_E",&muon_E,"muon_E[nMuons]/D");
      looseTree->Branch("muon_d0",&muon_d0,"muon_d0[nMuons]/D");
      looseTree->Branch("muon_relIso",&muon_relIso,"muon_relIso[nMuons]/D");
      
      looseTree->Branch("nLooseMuons",&nLooseMuons, "nLooseMuons/I");
      looseTree->Branch("muon_loose_charge",&muon_loose_charge,"muon_loose_charge[nLooseMuons]/I");
      looseTree->Branch("muon_loose_pt",&muon_loose_pt,"muon_loose_pt[nLooseMuons]/D");
      looseTree->Branch("muon_loose_phi",&muon_loose_phi,"muon_loose_phi[nLooseMuons]/D");
      looseTree->Branch("muon_loose_eta",&muon_loose_eta,"muon_loose_eta[nLooseMuons]/D");
      looseTree->Branch("muon_loose_E",&muon_loose_E,"muon_loose_E[nLooseMuons]/D");
      looseTree->Branch("muon_loose_d0",&muon_loose_d0,"muon_loose_d0[nLooseMuons]/D");
      looseTree->Branch("muon_loose_relIso",&muon_loose_relIso,"muon_loose_relIso[nLooseMuons]/D");
      looseTree->Branch("isGlobalLooseMuon",&isGlobalLooseMuon,"isGlobalLooseMuon[nLooseMuons]/O");
      looseTree->Branch("isTrackerLooseMuon",&isTrackerLooseMuon,"isTrackerLooseMuon[nLooseMuons]/O");
    }
    
    // jets
    myTree->Branch("nJets",&nJets,"nJets/I");
    myTree->Branch("jet_nConstituents",&jet_nConstituents,"jet_nConstituents[nJets]/I");
    myTree->Branch("jet_nChConstituents",&jet_nChConstituents,"jet_nChConstituents[nJets]/I");
    myTree->Branch("jet_charge",&jet_charge,"jet_charge[nJets]/I");
    myTree->Branch("jet_pt",&jet_pt,"jet_pt[nJets]/D");
    myTree->Branch("jet_phi",&jet_phi,"jet_phi[nJets]/D");
    myTree->Branch("jet_eta",&jet_eta,"jet_eta[nJets]/D");
    myTree->Branch("jet_E",&jet_E,"jet_E[nJets]/D");
    myTree->Branch("jet_M",&jet_M,"jet_M[nJets]/D");
    myTree->Branch("jet_bdiscr",&jet_bdiscr,"jet_bdiscr[nJets]/D");
    
    if (fillLooseTree)
    {
      looseTree->Branch("nJets",&nJets,"nJets/I");
      looseTree->Branch("jet_nConstituents",&jet_nConstituents,"jet_nConstituents[nJets]/I");
      looseTree->Branch("jet_nChConstituents",&jet_nChConstituents,"jet_nChConstituents[nJets]/I");
      looseTree->Branch("jet_charge",&jet_charge,"jet_charge[nJets]/I");
      looseTree->Branch("jet_pt",&jet_pt,"jet_pt[nJets]/D");
      looseTree->Branch("jet_phi",&jet_phi,"jet_phi[nJets]/D");
      looseTree->Branch("jet_eta",&jet_eta,"jet_eta[nJets]/D");
      looseTree->Branch("jet_E",&jet_E,"jet_E[nJets]/D");
      looseTree->Branch("jet_M",&jet_M,"jet_M[nJets]/D");
      looseTree->Branch("jet_bdiscr",&jet_bdiscr,"jet_bdiscr[nJets]/D");
      
      looseTree->Branch("nLooseJets",&nLooseJets,"nLooseJets/I");
      looseTree->Branch("jet_loose_nConstituents",&jet_loose_nConstituents,"jet_loose_nConstituents[nLooseJets]/I");
      looseTree->Branch("jet_loose_nChConstituents",&jet_loose_nChConstituents,"jet_loose_nChConstituents[nLooseJets]/I");
      looseTree->Branch("jet_loose_charge",&jet_loose_charge,"jet_loose_charge[nLooseJets]/I");
      looseTree->Branch("jet_loose_pt",&jet_loose_pt,"jet_loose_pt[nLooseJets]/D");
      looseTree->Branch("jet_loose_phi",&jet_loose_phi,"jet_loose_phi[nLooseJets]/D");
      looseTree->Branch("jet_loose_eta",&jet_loose_eta,"jet_loose_eta[nLooseJets]/D");
      looseTree->Branch("jet_loose_E",&jet_loose_E,"jet_loose_E[nLooseJets]/D");
      looseTree->Branch("jet_loose_M",&jet_loose_M,"jet_loose_M[nLooseJets]/D");
      looseTree->Branch("jet_loose_bdiscr",&jet_loose_bdiscr,"jet_loose_bdiscr[nLooseJets]/D");
    }
    
    // met
    myTree->Branch("met_px", &met_px, "met_px/D");
    myTree->Branch("met_py", &met_py, "met_py/D");
    myTree->Branch("met_pt", &met_pt, "met_pt/D");
    myTree->Branch("met_phi", &met_phi, "met_phi/D");
    myTree->Branch("met_eta", &met_eta,"met_eta/D");
    myTree->Branch("met_Et", &met_Et,"met_Et/D");
    myTree->Branch("met_E", &met_E,"met_E/D");
    
    myTree->Branch("met_corr_px", &met_corr_px, "met_corr_px/D");
    myTree->Branch("met_corr_py", &met_corr_py, "met_corr_py/D");
    myTree->Branch("met_corr_pt", &met_corr_pt, "met_corr_pt/D");
    myTree->Branch("met_corr_phi", &met_corr_phi, "met_corr_phi/D");
    myTree->Branch("met_corr_eta", &met_corr_eta,"met_corr_eta/D");
    myTree->Branch("met_corr_Et", &met_corr_Et,"met_corr_Et/D");
    myTree->Branch("met_corr_E", &met_corr_E,"met_corr_E/D");
    
    
    // mcparticles
    if (! isData)
    {
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
      myTree->Branch("hasGenTop", &hasGenTop, "hasGenTop/O");
      myTree->Branch("hasGenTopWithStatus22", &hasGenTopWithStatus22, "hasGenTopWithStatus22/O");
      myTree->Branch("hasGenTopWithStatus62", &hasGenTopWithStatus62, "hasGenTopWithStatus62/O");
      myTree->Branch("hasGenAntiTop", &hasGenAntiTop, "hasGenAntiTop/O");
      myTree->Branch("hasGenAntiTopWithStatus22", &hasGenAntiTopWithStatus22, "hasGenAntiTopWithStatus22/O");
      myTree->Branch("hasGenAntiTopWithStatus62", &hasGenAntiTopWithStatus62, "hasGenAntiTopWithStatus62/O");
    }
    
    
    /// SFs
    if (! isData)
    {
      if (nlo) myTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
      myTree->Branch("btagSF",&btagSF,"btagSF/D");
      myTree->Branch("btagSF_up",&btagSF_up,"btagSF_up/D");
      myTree->Branch("btagSF_down",&btagSF_down,"btagSF_down/D");
      myTree->Branch("puSF",&puSF,"puSF/D");
      myTree->Branch("puSF_up",&puSF_up,"puSF_up/D");
      myTree->Branch("puSF_down",&puSF_down,"puSF_down/D");
      
      myTree->Branch("muonIdSF_BCDEF",&muonIdSF_BCDEF,"muonIdSF_BCDEF[nMuons]/D");
      myTree->Branch("muonIdSF_GH",&muonIdSF_GH,"muonIdSF_GH[nMuons]/D");
      myTree->Branch("muonIdSF_up_BCDEF",&muonIdSF_up_BCDEF,"muonIdSF_up_BCDEF[nMuons]/D");
      myTree->Branch("muonIdSF_up_GH",&muonIdSF_up_GH,"muonIdSF_up_GH[nMuons]/D");
      myTree->Branch("muonIdSF_down_BCDEF",&muonIdSF_down_BCDEF,"muonIdSF_down_BCDEF[nMuons]/D");
      myTree->Branch("muonIdSF_down_GH",&muonIdSF_down_GH,"muonIdSF_down_GH[nMuons]/D");
      myTree->Branch("muonIsoSF_BCDEF",&muonIsoSF_BCDEF, "muonIsoSF_BCDEF[nMuons]/D");
      myTree->Branch("muonIsoSF_GH",&muonIsoSF_GH, "muonIsoSF_GH[nMuons]/D");
      myTree->Branch("muonIsoSF_up_BCDEF",&muonIsoSF_up_BCDEF, "muonIsoSF_up_BCDEF[nMuons]/D");
      myTree->Branch("muonIsoSF_up_GH",&muonIsoSF_up_GH, "muonIsoSF_up_GH[nMuons]/D");
      myTree->Branch("muonIsoSF_down_BCDEF",&muonIsoSF_down_BCDEF, "muonIsoSF_down_BCDEF[nMuons]/D");
      myTree->Branch("muonIsoSF_down_GH",&muonIsoSF_down_GH, "muonIsoSF_down_GH[nMuons]/D");
      myTree->Branch("muonTrigSF_BCDEF",&muonTrigSF_BCDEF,"muonTrigSF_BCDEF[nMuons]/D");
      myTree->Branch("muonTrigSF_GH",&muonTrigSF_GH,"muonTrigSF_GH[nMuons]/D");
      myTree->Branch("muonTrigSF_up_BCDEF",&muonTrigSF_up_BCDEF,"muonTrigSF_up_BCDEF[nMuons]/D");
      myTree->Branch("muonTrigSF_up_GH",&muonTrigSF_up_GH,"muonTrigSF_up_GH[nMuons]/D");
      myTree->Branch("muonTrigSF_down_BCDEF",&muonTrigSF_down_BCDEF,"muonTrigSF_down_BCDEF[nMuons]/D");
      myTree->Branch("muonTrigSF_down_GH",&muonTrigSF_down_GH,"muonTrigSF_down_GH[nMuons]/D");
      myTree->Branch("muonTrackSF_eta",&muonTrackSF_eta,"muonTrackSF_eta[nMuons]/D");
      myTree->Branch("muonTrackSF_aeta",&muonTrackSF_aeta,"muonTrackSF_aeta[nMuons]/D");
      myTree->Branch("muonTrackSF_nPV",&muonTrackSF_nPV,"muonTrackSF_nPV[nMuons]/D");
    }
    
    
    
    ////////////////////////////
    ///  Determine range of events to run over
    ////////////////////////////
    
    
    unsigned int ending = datasets[d]->NofEvtsToRunOver();
    double end_d = ending;
    if ( localgridSubmission && endEvent < ending )
      end_d = endEvent;
    if (test) end_d = 2000;  // for testing
    
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
    
    
    
    /////////////////////////
    ///  Define internal variables
    /////////////////////////
    
    bool isGoodPV;
    bool isBadMuon, isCloneMuon;
    bool toBeErased;
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    nEvents = 0;
    nEventsSel = 0;
    nofPosWeights = 0;
    nofNegWeights = 0;
    sumW = 0.;
    nofEventsRunB = 0;
    nofEventsRunCD = 0;
    nofEventsRunEF = 0;
    nofEventsRunG = 0;
    nofEventsRunH = 0;
    nofSelEventsRunB = 0;
    nofSelEventsRunCD = 0;
    nofSelEventsRunEF = 0;
    nofSelEventsRunG = 0;
    nofSelEventsRunH = 0;
    nofEventsWithGenTop = 0;
    nofEventsWithGenTopWithStatus22or62 = 0;
    nofEventsWithGenAntiTop = 0;
    nofEventsWithGenAntiTopWithStatus22or62 = 0;
    nofTTEventsWithoutBothGenTops = 0;
    nofTTEventsWithoutGenTop = 0;
    nofTTEventsWithoutGenAntiTop = 0;
    nofTTEventsWithoutBothGenTopsWithStatus22 = 0;
    nofTTEventsWithoutGenTopWithStatus22 = 0;
    nofTTEventsWithoutGenAntiTopWithStatus22 = 0;
    nofTTEventsWithoutBothGenTopsWithStatus62 = 0;
    nofTTEventsWithoutGenTopWithStatus62 = 0;
    nofTTEventsWithoutGenAntiTopWithStatus62 = 0;
    
    for (Int_t i = 0; i < 10; i++)
    {
      cutFlow[i] = 0;
      cutFlow2[i] = 0;
    }
    
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
    vector < TRootMET* > mets_corrected;
    vector < TRootGenJet* > genjets;
    vector < TRootMCParticle* > mcParticles;
    
    vector < TRootPFJet* > selectedJets;
    vector < TRootPFJet* > selectedBJets;
    vector < TRootPFJet* > selectedJetsBC;
    vector < TRootPFJet* > selectedLooseJets;
    vector < TRootPFJet* > selectedLooseJetsBC;
    vector < TRootMuon* > selectedMuons;
    vector < TRootMuon* > selectedMuonsBC;
    vector < TRootMuon* > vetoMuons;
    vector < TRootMuon* > selectedLooseMuons;
    vector < TRootMuon* > selectedLooseMuonsBC;
    vector < TRootElectron* > selectedElectrons;
    vector < TRootElectron* > vetoElectrons;
    vector < TRootElectron* > selectedLooseElectrons;
    
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
      mets_corrected.clear();
      genjets.clear();
      mcParticles.clear();
      
      selectedJets.clear();
      selectedBJets.clear();
      selectedJetsBC.clear();
      selectedLooseJets.clear();
      selectedLooseJetsBC.clear();
      selectedMuons.clear();
      selectedMuonsBC.clear();
      vetoMuons.clear();
      selectedLooseMuons.clear();
      selectedLooseMuonsBC.clear();
      selectedElectrons.clear();
      vetoElectrons.clear();
      selectedLooseElectrons.clear();
      
      /// Reset other stuff
      isTrigged = false;
      isSelected = false;
      hasExactly4Jets = false;
      hasJetLeptonCleaning = false;
      hasLooseJetLeptonCleaning = false;
      hasErasedBadOrCloneMuon = false;
      hasErasedBadOrCloneLooseMuon = false;
      
      filterHBHENoise = false;
      filterHBHEIso = false;
      filterCSCTightHalo = false;
      filterEcalDeadCell = false;
      filterEEBadSc = false;  // recommended for data-only
      filterBadChCand = false;
      filterBadMuon = false;
      passedMETFilter = false;
      
      isDataRunB = false;
      isDataRunC = false;
      isDataRunD = false;
      isDataRunE = false;
      isDataRunF = false;
      isDataRunG = false;
      isDataRunH = false;
      hasPosWeight = false;
      hasNegWeight = false;
      hasGenTop = false;
      hasGenTopWithStatus22 = false;
      hasGenTopWithStatus62 = false;
      hasGenAntiTop = false;
      hasGenAntiTopWithStatus22 = false;
      hasGenAntiTopWithStatus62 = false;
      
      nloWeight = 1.; // for amc@nlo samples
      btagSF = 1.;
      btagSF_up = 1.;
      btagSF_down = 1.;
      puSF = 1.;
      puSF_up = 1.;
      puSF_down = 1.;
      for (Int_t i = 0; i < 10; i++)
      {
        muonIdSF_BCDEF[i] = 1.;
        muonIdSF_GH[i] = 1.;
        muonIdSF_up_BCDEF[i] = 1.;
        muonIdSF_up_GH[i] = 1.;
        muonIdSF_down_BCDEF[i] = 1.;
        muonIdSF_down_GH[i] = 1.;
        muonIsoSF_BCDEF[i] = 1.;
        muonIsoSF_GH[i] = 1.;
        muonIsoSF_up_BCDEF[i] = 1.;
        muonIsoSF_up_GH[i] = 1.;
        muonIsoSF_down_BCDEF[i] = 1.;
        muonIsoSF_down_GH[i] = 1.;
        muonTrigSF_BCDEF[i] = 1.;
        muonTrigSF_GH[i] = 1.;
        muonTrigSF_up_BCDEF[i] = 1.;
        muonTrigSF_up_GH[i] = 1.;
        muonTrigSF_down_BCDEF[i] = 1.;
        muonTrigSF_down_GH[i] = 1.;
        muonTrackSF_eta[i] = 1.;
        muonTrackSF_aeta[i] = 1.;
        muonTrackSF_nPV[i] = 1.;
      }
      
      
      nMuons = -1;
      nLooseMuons = -1;
      nJets = -1;
      nLooseJets = -1;
      
      for (Int_t i = 0; i < 10; i++)
      {
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
        muon_loose_charge[i] = 0;
        muon_loose_pt[i] = 0.;
        muon_loose_phi[i] = 0.;
        muon_loose_eta[i] = 0.;
        muon_loose_E[i] = 0.;
        muon_loose_d0[i] = 999.;
        muon_loose_relIso[i] = 999.;
        isGlobalLooseMuon[i] = false;
        isTrackerLooseMuon[i] = false;
      }
      
      for (Int_t i = 0; i < 20; i++)
      {
        jet_nConstituents[i] = 0.;
        jet_nChConstituents[i] = 0.;
        jet_charge[i] = 0;
        jet_pt[i] = 0.;
        jet_phi[i] = 0.;
        jet_eta[i] = 0.;
        jet_E[i] = 0.;
        jet_M[i] = 0.;
        jet_bdiscr[i] = -1.;
      }
      
      for (Int_t i = 0; i < 30; i++)
      {
        jet_loose_nConstituents[i] = 0;
        jet_loose_nChConstituents[i] = 0;
        jet_loose_charge[i] = 0;
        jet_loose_pt[i] = 0.;
        jet_loose_phi[i] = 0.;
        jet_loose_eta[i] = 0.;
        jet_loose_E[i] = 0.;
        jet_loose_M[i] = 0.;
        jet_loose_bdiscr[i] = 0.;
      }
      
      met_px = 0.;
      met_py = 0.;
      met_pt = 0.;
      met_phi = 0.;
      met_eta = 0.;
      met_Et = 0.;
      met_E = 0.;
      
      met_corr_px = 0.;
      met_corr_py = 0.;
      met_corr_pt = 0.;
      met_corr_phi = 0.;
      met_corr_eta = 0.;
      met_corr_Et = 0.;
      met_corr_E = 0.;
      
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
      
      /// Clear vars
      isGoodPV = false;
      isBadMuon = false;
      isCloneMuon = false;
      toBeErased = false;
      
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets, mets, false);
      init_jets_corrected = init_jets;
      mets_corrected = mets;
      
      datasets[d]->eventTree()->LoadTree(ievt);
      run_num = event->runId();
      evt_num = event->eventId();
      lumi_num = event->lumiBlockId();
      nvtx = vertex.size();
      npu = (int)event->nTruePU();
      rho = event->fixedGridRhoFastjetAll();
      
      if (isData)
      {
        if ( run_num < 272007 )
        {
          cerr << "-- Dataset 2016A included..." << endl;
          cerr << "   Run number is " << run_num << endl;
          //cerr << "   File name is " << datasets[d]->eventTree()->GetFile()->GetName() << endl;
          exit(1);
        }
        else if ( run_num >= 272007 && run_num <= 275376 )
        {
          nofEventsRunB++;
          isDataRunB = true;
        }
        else if ( run_num >= 275657 && run_num <= 276283 )
        {
          nofEventsRunCD++;
          isDataRunC = true;
        }
        else if ( run_num >= 276315 && run_num <= 276811 )
        {
          nofEventsRunCD++;
          isDataRunD = true;
        }
        else if ( run_num >= 276831 && run_num <= 277420 )
        {
          nofEventsRunEF++;
          isDataRunE = true;
        }
        else if ( run_num >= 277772 && run_num <= 278808 )
        {
          nofEventsRunEF++;
          isDataRunF = true;
        }
        else if ( run_num >= 278820 && run_num <= 280385 )
        {
          nofEventsRunG++;
          isDataRunG = true;
        }
        else if ( run_num >= 280919 && run_num <= 284044 )
        {
          nofEventsRunH++;
          isDataRunH = true;
        }
        else
        {
          cerr << "-- This run is not recognised..." << endl;
          cerr << "   Run number is " << run_num << endl;
          //cerr << "   File name is " << datasets[d]->eventTree()->GetFile()->GetName() << endl;
          exit(1);
        }
      }
      
      if (! isData )
      {
        genjets = treeLoader.LoadGenJet(ievt,false);
        treeLoader.LoadMCEvent(ievt, 0, mcParticles, false);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      
      
      ////////////////////////////
      ///  ENERGY CORRECTIONS  ///
      ////////////////////////////
      
      jetTools->correctJets(init_jets_corrected, rho, isData);
      jetTools->correctMETTypeOne(init_jets_corrected, mets_corrected[0], isData);
      
      if (applyJESdown)
        jetTools->correctJetJESUnc(init_jets_corrected, mets_corrected[0], "minus", 1); // with or without met?
      else if (applyJESup)
        jetTools->correctJetJESUnc(init_jets_corrected, mets_corrected[0], "plus", 1);
      
      if (! isData)
      {
        if (applyJERdown)
          jetTools->correctJetJER(init_jets_corrected, genjets, mets_corrected[0], "minus", false);
        else if (applyJERup)
          jetTools->correctJetJER(init_jets_corrected, genjets, mets_corrected[0], "plus", false);
        else
          jetTools->correctJetJER(init_jets_corrected, genjets, mets_corrected[0], "nominal", false);
      }
      
      
      
      /////////////////
      ///  Trigger  ///
      /////////////////
      
      trigger->checkAvail(run_num, datasets, d, &treeLoader, event, false);
      isTrigged = trigger->checkIfFired();
      
      
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      /// Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets_corrected, rho);
      
      isGoodPV = selection.isPVSelected(vertex, 4, 24., 2.);
      
      selectedJets = selection.GetSelectedJets(jetPT, jetEta, true, "Tight");  // PtThr, EtaThr, applyJetID, TightLoose
      selectedLooseJets = selection.GetSelectedJets(jetPT, jetEta, true, "Loose");
      selectedMuons = selection.GetSelectedMuons(muonPTSel, muonEtaSel, muonRelIsoSel, muonWP, "Summer16");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      selectedElectrons = selection.GetSelectedElectrons(electronPTSel, electronEtaSel, electronWP, "Spring16_80X", true, true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased, VID
      vetoMuons = selection.GetSelectedMuons(muonPTVeto, muonEtaVeto, muonRelIsoVeto, "Loose", "Summer16");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      vetoElectrons = selection.GetSelectedElectrons(electronPTVeto, electronEtaVeto, "Veto", "Spring16_80X", true, true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased
      selectedLooseMuons = selection.GetSelectedMuons(muonPTSel, muonEtaSel, muonRelIsoVeto, "Loose", "Summer16");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      selectedLooseElectrons = selection.GetSelectedElectrons(electronPTSel, electronEtaSel, "Loose", "Spring16_80X", true, true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased, VID
      
      /// Reject bad muons
      isBadMuon = false, isCloneMuon = false;
      selectedMuonsBC.clear();
      selectedMuonsBC = selectedMuons;
      selectedMuons.clear();
      for (int iOrigMuon = 0; iOrigMuon < selectedMuonsBC.size(); iOrigMuon++)
      {
        isBadMuon = selectedMuonsBC[iOrigMuon]->isBad80X();
        isCloneMuon = selectedMuonsBC[iOrigMuon]->isClone80X();
        if (! isBadMuon && ! isCloneMuon) selectedMuons.push_back(selectedMuonsBC[iOrigMuon]);
        else hasErasedBadOrCloneMuon = true;
      }
      
      isBadMuon = false; isCloneMuon = false;
      selectedLooseMuonsBC.clear();
      selectedLooseMuonsBC = selectedLooseMuons;
      selectedLooseMuons.clear();
      for (int iOrigMuon = 0; iOrigMuon < selectedLooseMuonsBC.size(); iOrigMuon++)
      {
        isBadMuon = selectedLooseMuonsBC[iOrigMuon]->isBad80X();
        isCloneMuon = selectedLooseMuonsBC[iOrigMuon]->isClone80X();
        if (! isBadMuon && ! isCloneMuon) selectedLooseMuons.push_back(selectedLooseMuonsBC[iOrigMuon]);
        else hasErasedBadOrCloneLooseMuon = true;
      }
      
      
      if (applyJetLeptonCleaning)
      {
        if(verbose > 3) cout << "  - Applying jet/lepton cleaning... " << endl; 
        
        selectedJetsBC.clear();
        selectedJetsBC = selectedJets;
        selectedJets.clear();
        
        for (int iOrigJet = 0; iOrigJet < selectedJetsBC.size(); iOrigJet++)
        {
          toBeErased = false;
          for (int iMuon = 0; iMuon < selectedLooseMuons.size(); iMuon++)
          {
            if ( selectedJetsBC[iOrigJet]->DeltaR(*selectedLooseMuons[iMuon]) < 0.4 )
            {
              toBeErased = true;
              hasJetLeptonCleaning = true;
              break;
            }
          }
          if (toBeErased) continue;
          for (int iElectron = 0; iElectron < selectedLooseElectrons.size(); iElectron++)
          {
            if ( selectedJetsBC[iOrigJet]->DeltaR(*selectedLooseElectrons[iElectron]) < 0.3 )
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
        
        
        /// Loose jets
        selectedLooseJetsBC.clear();
        selectedLooseJetsBC = selectedLooseJets;
        selectedLooseJets.clear();
        
        for (int iOrigJet = 0; iOrigJet < selectedLooseJetsBC.size(); iOrigJet++)
        {
          bool toBeErased = false;
          for (int iMuon = 0; iMuon < selectedLooseMuons.size(); iMuon++)
          {
            if ( selectedLooseJetsBC[iOrigJet]->DeltaR(*selectedLooseMuons[iMuon]) < 0.4 )
            {
              toBeErased = true;
              hasLooseJetLeptonCleaning = true;
              break;
            }
          }
          if (toBeErased) continue;
          for (int iElectron = 0; iElectron < selectedLooseElectrons.size(); iElectron++)
          {
            if ( selectedLooseJetsBC[iOrigJet]->DeltaR(*selectedLooseElectrons[iElectron]) < 0.3 )
            {
              toBeErased = true;
              hasLooseJetLeptonCleaning = true;
              break;
            }
          }
          if (! toBeErased)
          {
            selectedLooseJets.push_back(selectedLooseJetsBC[iOrigJet]);
          }
        }
        if ( verbose > 3 )
        {
          if ( selectedLooseJetsBC.size() != selectedLooseJets.size() ) cout << "--> original = " << selectedLooseJetsBC.size()  << " after cleaning = " << selectedLooseJets.size() << endl;
        }
                
      }  // end jet cleaning
      
      
      for (int i = 0; i < selectedJets.size(); i++)
      {
        if ( selectedJets[i]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium )
          selectedBJets.push_back(selectedJets[i]);
      }
      
      
      /// Fill tree with loose objects before event selection
      //  Define jet & muon variables
      if ( ! calculateBTagSF )
      {
        nJets = selectedJets.size();
        nLooseJets = selectedLooseJets.size();
        nMuons = selectedMuons.size();
        nLooseMuons = selectedLooseMuons.size();
        
        for(Int_t iJet = 0; iJet < nJets; iJet++)
        {
          jet_charge[iJet] = selectedJets[iJet]->charge();
          jet_nConstituents[iJet] = selectedJets[iJet]->nConstituents();
          jet_nChConstituents[iJet] = selectedJets[iJet]->chargedMultiplicity();
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
        
        // TO DO: Include nominal event scale factors in loose tree (muon, PU)
        
        if ( fillLooseTree && nLooseMuons > 0 )
        {
          for(Int_t iJet = 0; iJet < nLooseJets; iJet++)
          {
            jet_loose_charge[iJet] = selectedLooseJets[iJet]->charge();
            jet_loose_nConstituents[iJet] = selectedLooseJets[iJet]->nConstituents();
            jet_loose_nChConstituents[iJet] = selectedLooseJets[iJet]->chargedMultiplicity();
            jet_loose_pt[iJet] = selectedLooseJets[iJet]->Pt();
            jet_loose_phi[iJet] = selectedLooseJets[iJet]->Phi();
            jet_loose_eta[iJet] = selectedLooseJets[iJet]->Eta();
            jet_loose_E[iJet] = selectedLooseJets[iJet]->E();
            jet_loose_M[iJet] = selectedLooseJets[iJet]->M();
            jet_loose_bdiscr[iJet] = selectedLooseJets[iJet]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
          }
          
          for (Int_t iMuon = 0; iMuon < nLooseMuons; iMuon++)
          {
            muon_loose_charge[iMuon] = selectedLooseMuons[iMuon]->charge();
            muon_loose_pt[iMuon] = selectedLooseMuons[iMuon]->Pt();
            muon_loose_phi[iMuon] = selectedLooseMuons[iMuon]->Phi();
            muon_loose_eta[iMuon] = selectedLooseMuons[iMuon]->Eta();
            muon_loose_E[iMuon] = selectedLooseMuons[iMuon]->E();
            muon_loose_d0[iMuon] = selectedLooseMuons[iMuon]->d0();
            muon_loose_relIso[iMuon] = ( selectedLooseMuons[iMuon]->chargedHadronIso(4) + max( 0.0, selectedLooseMuons[iMuon]->neutralHadronIso(4) + selectedLooseMuons[iMuon]->photonIso(4) - 0.5*selectedLooseMuons[iMuon]->puChargedHadronIso(4) ) ) / muon_loose_pt[iMuon];  // dR = 0.4, dBeta corrected
            isGlobalLooseMuon[iMuon] = selectedLooseMuons[iMuon]->isGlobalMuon();
            isTrackerLooseMuon[iMuon] = selectedLooseMuons[iMuon]->isTrackerMuon();
          }
          
          looseTree->Fill();
        }
      }
      
      
      ////// Selection
      cutFlow[0]++;
      cutFlow2[0]++;
      if (isTrigged)
      {
        cutFlow[1]++;
        cutFlow2[1]++;
        if (isGoodPV)
        {
          cutFlow[2]++;
          cutFlow2[2]++;
          if (selectedMuons.size() == 1)
          {
            cutFlow[3]++;
            cutFlow2[3]++;
            if (vetoMuons.size() == 1)
            {
              cutFlow[4]++;
              cutFlow2[4]++;
              if (vetoElectrons.size() == 0)
              {
                cutFlow[5]++;
                cutFlow2[5]++;
                if ( selectedJets.size() >= 4 )
                {
                  cutFlow[6]++;
                  cutFlow2[6]++;
                  if ( selectedJets.size() == 4 )
                  {
                    cutFlow2[7]++;
                    hasExactly4Jets = true;
                  }
                  
                  if ( selectedBJets.size() > 0 )
                  {
                    cutFlow[7]++;
                    if (hasExactly4Jets) cutFlow2[8]++;
                    if ( selectedBJets.size() > 1 )
                    {
                      cutFlow[8]++;
                      isSelected = true;
                      
                      if (hasExactly4Jets)
                      {
                        cutFlow[9]++;
                        cutFlow2[9]++;
                      }
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
        if (isDataRunB) nofSelEventsRunB++;
        else if (isDataRunC || isDataRunD) nofSelEventsRunCD++;
        else if (isDataRunE || isDataRunF) nofSelEventsRunEF++;
        else if (isDataRunG) nofSelEventsRunG++;
        else if (isDataRunH) nofSelEventsRunH++;
      }
      
      
      /// B-tagging
      if (calculateBTagSF && ! isData)
      {
        bTagHistoTool_M->FillMCEfficiencyHistos(selectedJets);
        bTagHistoTool_M_up->FillMCEfficiencyHistos(selectedJets);
        bTagHistoTool_M_down->FillMCEfficiencyHistos(selectedJets);
        continue;
      }
      
      /// The next part is not run when b-tag SFs are calculated
      /// The tree is not written
      
      
      
      ///////////////////////////////
      ///  Fill Object Variables  ///
      ///////////////////////////////
      
      
      met_px = mets[0]->Px();
      met_py = mets[0]->Py();
      met_pt = sqrt(met_px*met_px + met_py*met_py);
      met_phi = mets[0]->Phi();
      met_eta = mets[0]->Eta();
      met_Et = mets[0]->Et();
      met_E = mets[0]->E();
      
      met_corr_px = mets_corrected[0]->Px();
      met_corr_py = mets_corrected[0]->Py();
      met_corr_pt = sqrt(met_corr_px*met_corr_px + met_corr_py*met_corr_py);
      met_corr_phi = mets_corrected[0]->Phi();
      met_corr_eta = mets_corrected[0]->Eta();
      met_corr_Et = mets_corrected[0]->Et();
      met_corr_E = mets_corrected[0]->E();
      
      filterHBHENoise = event->getHBHENoiseFilter();
      filterHBHEIso = event->getHBHENoiseIsoFilter();
      filterCSCTightHalo = event->getglobalTightHalo2016Filter();
      filterEcalDeadCell = event->getEcalDeadCellTriggerPrimitiveFilter();
      if (isData) filterEEBadSc = event->getEEBadScFilter();  // recommended for data-only
      else filterEEBadSc = true;
      filterBadChCand = event->getBadChCandFilter();
      filterBadMuon = event->getBadPFMuonFilter();
      
      if ( isGoodPV && filterHBHENoise && filterHBHEIso && filterCSCTightHalo && filterEcalDeadCell && filterEEBadSc && filterBadChCand && filterBadMuon ) passedMETFilter = true;
      
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
          
          if ( mc_pdgId[iMC] == 6 )
          {
            hasGenTop = true;
            nofEventsWithGenTop++;
            if ( mc_status[iMC] == 22 || mc_status[iMC] == 62 )
            {
              nofEventsWithGenTopWithStatus22or62++;
              if ( mc_status[iMC] == 22 ) hasGenTopWithStatus22 = true;
              else if ( mc_status[iMC] == 62 ) hasGenTopWithStatus62 = true;
            }
          }
          else if ( mc_pdgId[iMC] == -6 )
          {
            hasGenAntiTop = true;
            nofEventsWithGenAntiTop++;
            if ( mc_status[iMC] == 22 || mc_status[iMC] == 62 )
            {
              nofEventsWithGenAntiTopWithStatus22or62++;
              if ( mc_status[iMC] == 22 ) hasGenAntiTopWithStatus22 = true;
              else if ( mc_status[iMC] == 62 ) hasGenAntiTopWithStatus62 = true;
            }
          }
        }
        
        if (isTTbar)
        {
          if (! hasGenTop || ! hasGenAntiTop) nofTTEventsWithoutBothGenTops++;
          else if (! hasGenTop) nofTTEventsWithoutGenTop++;
          else if (! hasGenAntiTop) nofTTEventsWithoutGenAntiTop++;
          
          if (! hasGenTopWithStatus22 || ! hasGenAntiTopWithStatus22) nofTTEventsWithoutBothGenTopsWithStatus22++;
          else if (! hasGenTopWithStatus22) nofTTEventsWithoutGenTopWithStatus22++;
          else if (! hasGenAntiTopWithStatus22) nofTTEventsWithoutGenAntiTopWithStatus22++;
          
          if (! hasGenTopWithStatus62 || ! hasGenAntiTopWithStatus62) nofTTEventsWithoutBothGenTopsWithStatus62++;
          else if (! hasGenTopWithStatus62) nofTTEventsWithoutGenTopWithStatus62++;
          else if (! hasGenAntiTopWithStatus62) nofTTEventsWithoutGenAntiTopWithStatus62++;
        }
      }
      
      
      /// Fill scalefactors
      if (! isData)
      {
        if (applyBTagSF)
        {
          btagSF = bTagHistoTool_M->getMCEventWeight(selectedJets, false);
          btagSF_up = bTagHistoTool_M_up->getMCEventWeight(selectedJets, false);
          btagSF_down = bTagHistoTool_M_down->getMCEventWeight(selectedJets, false);
        }
        
        puSF = LumiWeights.ITweight( (int)event->nTruePU() );
        puSF_up = LumiWeights_up.ITweight( (int)event->nTruePU() );
        puSF_down = LumiWeights_down.ITweight( (int)event->nTruePU() );
        
        for (int iMuon = 0; iMuon < selectedMuons.size(); iMuon++)
        {
          muonIdSF_BCDEF[iMuon] = muonSFWeightID_T_BCDEF->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);  // eta, pt, shiftUpDown;
          muonIdSF_GH[iMuon] = muonSFWeightID_T_GH->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
          muonIdSF_up_BCDEF[iMuon] = muonSFWeightID_T_BCDEF->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 1);
          muonIdSF_up_GH[iMuon] = muonSFWeightID_T_GH->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 1);
          muonIdSF_down_BCDEF[iMuon] = muonSFWeightID_T_BCDEF->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), -1);
          muonIdSF_down_GH[iMuon] = muonSFWeightID_T_GH->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), -1);
          
          muonIsoSF_BCDEF[iMuon] = muonSFWeightIso_TT_BCDEF->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
          muonIsoSF_GH[iMuon] = muonSFWeightIso_TT_GH->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
          muonIsoSF_up_BCDEF[iMuon] = muonSFWeightIso_TT_BCDEF->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 1);
          muonIsoSF_up_GH[iMuon] = muonSFWeightIso_TT_GH->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 1);
          muonIsoSF_down_BCDEF[iMuon] = muonSFWeightIso_TT_BCDEF->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), -1);
          muonIsoSF_down_GH[iMuon] = muonSFWeightIso_TT_GH->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), -1);
          
          muonTrigSF_BCDEF[iMuon] = muonSFWeightTrig_BCDEF->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
          muonTrigSF_GH[iMuon] =    muonSFWeightTrig_GH->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 0);
          muonTrigSF_up_BCDEF[iMuon] = muonSFWeightTrig_BCDEF->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 1);
          muonTrigSF_up_GH[iMuon] = muonSFWeightTrig_GH->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), 1);
          muonTrigSF_down_BCDEF[iMuon] = muonSFWeightTrig_BCDEF->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), -1);
          muonTrigSF_down_GH[iMuon] = muonSFWeightTrig_GH->at(selectedMuons[iMuon]->Eta(), selectedMuons[iMuon]->Pt(), -1);
          
          muonTrackSF_eta[iMuon] = h_muonSFWeightTrackEta->Eval(selectedMuons[iMuon]->Eta());
          muonTrackSF_aeta[iMuon] = h_muonSFWeightTrackAEta->Eval(fabs(selectedMuons[iMuon]->Eta()));
          muonTrackSF_nPV[iMuon] = h_muonSFWeightTrackPV->Eval(npu);
        }
        
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
      }
      
      
      if (! calculateBTagSF)
      {
        myTree->Fill();
      }
      
      
    }  // end loop events
    if (! isData && ! calculateBTagSF) cout << "Max MCParticles: " << maxMCParticles << endl << endl;
    
    if (! calculateBTagSF)
    {
      cout << "Fill trees..." << endl;
      statTree->Fill();
      
      /// Write to file
      fout->cd();
      myTree->Write();
      statTree->Write();
      if (fillLooseTree) looseTree->Write();
    }
    fout->Close();
    delete fout;
    
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  // end loop datasets
  
  /// To write plots b tagging:
  delete bTagHistoTool_M;
  delete bTagHistoTool_M_up;
  delete bTagHistoTool_M_down;
  
  
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
