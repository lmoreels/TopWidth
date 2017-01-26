////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////



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
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

// used TopTreeAnalysis classes
#include "../TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

// user defined
#include "Tools/interface/Trigger.h"


using namespace std;
using namespace reweight;
using namespace TopTree;


string ConvertIntToString(int Number, bool pad)
{
  ostringstream convert;
  convert.clear();
  if ( pad && Number < 10 ) { convert << std::setw(2) << std::setfill('0');}
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
  
  string year_str = ConvertIntToString(year, true);
  string month_str = ConvertIntToString(month, true);
  string day_str = ConvertIntToString(day, true);
  string hour_str = ConvertIntToString(hour, true);
  string min_str = ConvertIntToString(min, true);
  //string sec_str = ConvertIntToString(sec, true);
  
  string date_str = year_str + month_str + day_str + "_" + hour_str + min_str;
  return date_str;
}


int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "***********************************" << endl;
  cout << "***   Beginning of program      ***" << endl;
  cout << "***********************************" << endl;
  cout << "Current time: " << dateString << endl;
  
  clock_t start = clock();
  
  bool useOneFourthOfDataSets = false;
  bool useTestSample = true;
  
  
  string rootFileName = "testWeights_output_FullDataSet_"+dateString+".root";
  string selectiontableMu = "SelectionTable_testWeightsFull_SemiMu_"+dateString+".tex";
  string pathPNG = "PlotsWeights_"+dateString+"/";
  int iReducedDataSets = 1;
  
  if (useOneFourthOfDataSets)
  {
    rootFileName = "testWeights_output_oneFourthOfDataSets_"+dateString+".root";
    selectiontableMu = "SelectionTable_testWeightsOneFourth_SemiMu_"+dateString+".tex";
    pathPNG = "PlotsWeightsOneFourth_"+dateString+"/";
    iReducedDataSets = 4;
  }
  if (useTestSample)
  {
    rootFileName = "testWeights_output_testSample_"+dateString+".root";
    selectiontableMu = "SelectionTable_testWeightsSample_SemiMu_"+dateString+".tex";
    pathPNG = "PlotsWeightsTestSample_"+dateString+"/";
    iReducedDataSets = 200;
  }
  
  
  
  /////////////////////
  ///  Configuration
  /////////////////////
  
  bool printTriggers = false;
  bool applyTriggers = true;
  bool applyJER = true;
  bool applyJEC = true;
  bool applyLeptonSF = false;
  bool applyPU = true;
  bool nlo = false;
  bool hasNegWeight = false;
  bool eventSelected = false;
  bool has1bjet = false;
  bool has2bjets = false;
  int nofSelectedEvents = 0;
  int nb_bTaggedJets = 0;
  int nofEventsWith1BJet = 0;
  int nofEventsWith2BJets = 0;
  int nofNegWeights = 0;
  
  
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
  
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  float oldLuminosity = anaEnv.Luminosity;  // in 1/pb
  cout << "Analysis environment luminosity for rescaling " << oldLuminosity << endl;
  
  
  
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;
  //vector < Dataset* > datasetsMu;
  //vector < Dataset* > datasetsEl;
  
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets(datasets, xmlfile); cout << "Number of datasets: " << datasets.size() << endl;
  for (unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = 566.487920485;
  
  
  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    string dataSetName = datasets[d]->Name();
    
    if ( (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) && Luminosity > datasets[d]->EquivalentLumi() ) { Luminosity = datasets[d]->EquivalentLumi();}
    
    if ( dataSetName.find("QCD") == 0 ) { datasets[d]->SetColor(kYellow);}
    if ( dataSetName.find("TT") == 0 ) { datasets[d]->SetColor(kRed+1);}
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
    if ( dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") ==0 )
    { 
      datasets[d]->SetColor(kBlue-2);
      if ( dataSetName.find("tW") == 0 ) { datasets[d]->SetTitle("ST tW");}
    }
    //if (dataSetName.find("NP") == 0 )
    //{
    //	datasets[d]->SetTitle("Signal");
    //	datasets[d]->SetColor(kGreen+4);
    //}
  }
  
  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
  if ( iReducedDataSets != 1 )
  {
    Luminosity = Luminosity/((double) iReducedDataSets);
    cout << "Running over 1/" << iReducedDataSets << " of the dataset, so luminosity changed to " << Luminosity << endl;
  }
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  
  //Global variable
  //TRootEvent* event = 0;
  TRootRun *runInfos = new TRootRun();
  
  //nof selected events
  //double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  Double_t mc_baseweight = 0, mc_baseweight1 = 0, mc_baseweight2 = 0;
  
  
  
  ////////////////////////////////////
  ///  Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
	
  map<string,TH2F*> histo2D;
    
  histo2D["muon_SF_ID"] = new TH2F("muon_SF_ID", "Muon ID scale factors in function of #eta (x) and p_{T} (y); #eta, p_{T} [GeV]", 21, 0, 2.1, 15, 0, 150);
  histo2D["muon_SF_Iso"] = new TH2F("muon_SF_Iso", "Muon relIso scale factors in function of #eta (x) and p_{T} (y); #eta, p_{T} [GeV]", 21, 0, 2.1, 15, 0, 150);
  
  
  
  ////////////////////////////////////
  ///  MultiSamplePlot
  ////////////////////////////////////
  
  map<string,MultiSamplePlot*> MSPlot;
  
  /// Plots before event selection
  MSPlot["init_nJets"] = new MultiSamplePlot(datasets, "init_nJets", 13, -0.5, 12.5, "# jets");
  MSPlot["init_nMuons"] = new MultiSamplePlot(datasets, "init_nMuons", 13, -0.5, 12.5, "# muons");
  MSPlot["init_nElectrons"] = new MultiSamplePlot(datasets, "init_nElectrons", 13, -0.5, 12.5, "# electrons");
  MSPlot["init_nPVs_before"] = new MultiSamplePlot(datasets, "init_nPVs_before", 41, -0.5, 40.5, "# PVs before reweighting");
  MSPlot["init_nPVs_after"] = new MultiSamplePlot(datasets, "init_nPVs_after", 41, -0.5, 40.5, "# PVs after reweighting");
  
  MSPlot["init_leadingJet_pT"] = new MultiSamplePlot(datasets, "init_leadingJet_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["init_leadingJet_eta"] = new MultiSamplePlot(datasets, "init_leadingJet_eta", 30, -3, 3, "Eta");
  MSPlot["init_leadingJet_CSVv2Discr"] = new MultiSamplePlot(datasets, "init_leadingJet_CSVv2Discr", 80, -0.5, 1.5, "CSVv2 discriminant value");
  MSPlot["init_leadingMuon_pT"] = new MultiSamplePlot(datasets, "init_leadingMuon_pT", 22, 0, 440, "p_{T} [GeV]");
  MSPlot["init_leadingMuon_eta"] = new MultiSamplePlot(datasets, "init_leadingMuon_eta", 30, -3, 3, "Eta");
  MSPlot["init_leadingMuon_phi"] = new MultiSamplePlot(datasets, "init_leadingMuon_phi", 32, -3.2, 3.2, "Phi");
  MSPlot["init_muon_relIso"] = new MultiSamplePlot(datasets, "init_muon_relIso", 30, 0, 0.3, "relIso");
  MSPlot["init_muon_d0"] = new MultiSamplePlot(datasets, "init_muon_d0", 50, 0, 5, "d_{0}");
  MSPlot["init_leadingElectron_pT"] = new MultiSamplePlot(datasets, "init_leadingElectron_pT", 22, 0, 440, "p_{T} [GeV]");
  MSPlot["init_leadingElectron_eta"] = new MultiSamplePlot(datasets, "init_leadingElectron_eta", 60, -3, 3, "Eta");
  MSPlot["init_met_pT"] = new MultiSamplePlot(datasets, "init_met_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["init_met_eta"] = new MultiSamplePlot(datasets, "init_met_eta", 30, -3, 3, "Eta");
  MSPlot["init_met_phi"] = new MultiSamplePlot(datasets, "init_met_phi", 32, -3.2, 3.2, "Phi");
  
  /// Event Selection
  MSPlot["Selection"] = new MultiSamplePlot(datasets, "Selection", 13, -0.5, 12.5, "Cutflow");
  
  /// Plots after event selection
  MSPlot["muon_pT"] = new MultiSamplePlot(datasets, "muon_pT", 22, 0, 440, "p_{T} [GeV]");
  MSPlot["muon_eta"] = new MultiSamplePlot(datasets, "muon_eta", 30, -3, 3, "Eta");
  MSPlot["muon_phi"] = new MultiSamplePlot(datasets, "muon_phi", 32, -3.2, 3.2, "Phi");
  MSPlot["muon_relIso"] = new MultiSamplePlot(datasets, "muon_relIso", 30, 0, 0.3, "relIso");
  MSPlot["muon_d0"] = new MultiSamplePlot(datasets, "muon_d0", 50, 0, 5, "d_{0}");
  MSPlot["leadingJet_pT"] = new MultiSamplePlot(datasets, "leadingJet_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["jet2_pT"] = new MultiSamplePlot(datasets, "jet2_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["jet3_pT"] = new MultiSamplePlot(datasets, "jet3_pT", 25, 0, 500, "p_{T} [GeV]");
  MSPlot["jet4_pT"] = new MultiSamplePlot(datasets, "jet4_pT", 25, 0, 500, "p_{T} [GeV]");
  MSPlot["Ht_4leadingJets"] = new MultiSamplePlot(datasets,"Ht_4leadingJets", 60, 0, 1200, "H_{T} [GeV]");
  MSPlot["met_pT"] = new MultiSamplePlot(datasets, "met_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["met_eta"] = new MultiSamplePlot(datasets, "met_eta", 30, -3, 3, "Eta");
  MSPlot["met_phi"] = new MultiSamplePlot(datasets, "met_phi", 32, -3.2, 3.2, "Phi");
  
  MSPlot["nJets"] = new MultiSamplePlot(datasets, "nJets", 13, -0.5, 12.5, "# jets");
  MSPlot["nBJets"] = new MultiSamplePlot(datasets, "nBJets", 13, -0.5, 12.5, "# b jets");
  MSPlot["bJet1_pT"] = new MultiSamplePlot(datasets, "bJet1_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["bJet2_pT"] = new MultiSamplePlot(datasets, "bJet2_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["bJet1_CSVv2Discr"] = new MultiSamplePlot(datasets, "bJet1_CSVv2Discr", 20, 0.8, 1.3, "CSVv2 discriminant value");
  MSPlot["bJet2_CSVv2Discr"] = new MultiSamplePlot(datasets, "bJet2_CSVv2Discr", 20, 0.8, 1.3, "CSVv2 discriminant value");
  
  MSPlot["1b_muon_pT"] = new MultiSamplePlot(datasets, "1b_muon_pT", 22, 0, 440, "p_{T} [GeV]");
  MSPlot["1b_muon_eta"] = new MultiSamplePlot(datasets, "1b_muon_eta", 30, -3, 3, "Eta");
  MSPlot["1b_muon_phi"] = new MultiSamplePlot(datasets, "1b_muon_phi", 32, -3.2, 3.2, "Phi");
  MSPlot["1b_muon_relIso"] = new MultiSamplePlot(datasets, "1b_muon_relIso", 30, 0, 0.3, "relIso");
  MSPlot["1b_muon_d0"] = new MultiSamplePlot(datasets, "1b_muon_d0", 50, 0, 5, "d_{0}");
  MSPlot["1b_leadingJet_pT"] = new MultiSamplePlot(datasets, "1b_leadingJet_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["1b_jet2_pT"] = new MultiSamplePlot(datasets, "1b_jet2_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["1b_jet3_pT"] = new MultiSamplePlot(datasets, "1b_jet3_pT", 25, 0, 500, "p_{T} [GeV]");
  MSPlot["1b_jet4_pT"] = new MultiSamplePlot(datasets, "1b_jet4_pT", 25, 0, 500, "p_{T} [GeV]");
  MSPlot["1b_Ht_4leadingJets"] = new MultiSamplePlot(datasets,"1b_Ht_4leadingJets", 60, 0, 1200, "H_{T} [GeV]");
  
  MSPlot["nPVs_before"] = new MultiSamplePlot(datasets, "nPVs_before", 41, -0.5, 40.5, "# PVs before reweighting");
  MSPlot["nPVs_after"] = new MultiSamplePlot(datasets, "nPVs_after", 41, -0.5, 40.5, "# PVs after reweighting");
  
  /// Scale factors
  MSPlot["pileup_SF"] = new MultiSamplePlot(datasets,"pileup_SF", 80, 0, 4, "lumiWeight");
  MSPlot["weightIndex"] = new MultiSamplePlot(datasets,"weightIndex", 5, -2.5, 2.5, "0: None; 1: Central scale variation 1; 2: scale_variation 1");
  MSPlot["weightIndex_diffs"] = new MultiSamplePlot(datasets,"weightIndex_diffs", 2, -0.5, 1.5, "Central scale variation 1 equals scale_variation 1");
  MSPlot["nloWeight_good"] = new MultiSamplePlot(datasets,"nloWeight_good", 400, -2.0, 2.0, "weights for amc@nlo samples");
  MSPlot["nloWeight_bad"] = new MultiSamplePlot(datasets,"nloWeight_bad", 420, -21.0, 21.0, "weights for amc@nlo samples");
  
  
  
  ////////////////////////////////////
  ///  Selection table
  ////////////////////////////////////
  
  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  
  //CutsSelecTableSemiMu.push_back("Missing $E_T$");
  //CutsSelecTableSemiMu.push_back("$H_T$ cut");
  CutsSelecTableSemiMu.push_back("$\\geq$ 1 b-jet (CSVMv2)");
  CutsSelecTableSemiMu.push_back("$\\geq$ 2 b-jets (CSVMv2)");
  //CutsSelecTableSemiMu.push_back("actually trigged");
  
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  //selecTableSemiMu.SetLuminosity(LuminosityMu);
  selecTableSemiMu.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
  
  
  ////////////////////////////
  ///  Initialise trigger  ///
  ////////////////////////////
  
  //Trigger* trigger = new Trigger(hasMuon, hasElectron, trigSingleLep, trigDoubleLep);
  Trigger* trigger = new Trigger(1, 0, 1, 0);
  
  
  
  //////////////////////////////////
  ///  Initialise scale factors  ///
  //////////////////////////////////
  
  string pathCalLept = "../TopTreeAnalysisBase/Calibrations/LeptonSF/";
  string pathCalPileup = "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/";
  
  /// Leptons
  cout << " - Loading lepton scale factors ...";
  if (! applyLeptonSF) { cout << "    --- At the moment these are not used in the analysis";}
  cout << endl;
  
  double muonSFID, muonSFIso, muonSFTrig;
  //MuonSFWeight *muonSFWeight_ = new MuonSFWeight(pathCalLept+"Muon_SF_TopEA.root","SF_totErr", true, false, false); // (... , ... , extendRange, debug, print warning)
  MuonSFWeight *muonSFWeightID_T = new MuonSFWeight(pathCalLept+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);
  
  MuonSFWeight *muonSFWeightIso_TT = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
  
  MuonSFWeight *muonSFWeightTrigHLTv4p2 = new MuonSFWeight(pathCalLept+"SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins/abseta_pt_ratio", true, false, false);
  MuonSFWeight *muonSFWeightTrigHLTv4p3 = new MuonSFWeight(pathCalLept+"SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio", true, false, false);
  
  
  /// Pile-up
  cout << " - Loading pile-up scale factors ...";
  if (! applyPU) { cout << "   --- At the moment these are not used in the analysis";}
  cout << endl;
  
  //LumiReWeighting LumiWeights(pathCalPileup+"pileup_MC_RunIISpring15DR74-Asympt25ns.root",pathCalPileup+"pileup_2015Data74X_25ns-Run254231-258750Cert/nominal.root","pileup60","pileup");  // old PU
  LumiReWeighting LumiWeights(pathCalPileup+"pileup_MC_RunIISpring15DR74-Asympt25ns.root",pathCalPileup+"pileup_2015Data74X_25ns-Run246908-260627Cert.root","pileup50","pileup");  // new PU
  
  
  
  ///////////////////////////////
  ///  Single Muon Selection  ///
  ///////////////////////////////
  
  /// Updated 19/01/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO
  
  float muonPTSel = 26.; // GeV
  float muonEtaSel = 2.1;
  float muonRelIsoSel = 0.15;  // Tight muon
  string muonWP = "Tight";
  
  float muonPTVeto = 10.; // GeV
  float muonEtaVeto = 2.5;
  float muonRelIsoVeto = 0.25;  // Loose muon
  
  
  
  ///////////////////////////////////
  ///  Single Electron Selection  ///
  ///////////////////////////////////
  
  // To do
  float electronPTSel = 24.; // GeV
  float electronEtaSel = 2.1;  // because of electron trigger
  string electronWP = "Tight";
  
  float electronPTVeto = 15.; // GeV
  float electronEtaVeto = 2.5;
  
  
  
  ///////////////////////
  ///  Jet Selection  ///
  ///////////////////////
  
  // Updated 19/01/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopJME
  float jetPT = 20.; // GeV
  float jetEta = 2.4;  // to allow b tagging
  
  
  
  //////////////////////////////////////
  ///  Working points for b tagging  ///
  //////////////////////////////////////
  
  /// Updated 14/01/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
  
  float CSVv2Loose =  0.605;
  float CSVv2Medium = 0.890;
  float CSVv2Tight = 0.970;
  
  
  
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
    nofNegWeights = 0;
    nlo = false;
    bool isData = false;
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
    
    if ( dataSetName.find("DY") == 0 || dataSetName.find("ZJets") == 0 || dataSetName.find("Zjets") == 0 || dataSetName.find("WJets") == 0 || dataSetName.find("Wjets") == 0 || dataSetName.find("ST_tch") == 0 )
    {
      nlo = true;
    }
    else
    {
      cout << "Skipping dataset " << dataSetName << endl;
      continue;
    }
    
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << "/ title : " << datasets[d]->Title() << endl;
      cout << "      -> Equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
      cout << "      -> Norm factor equals " << datasets[d]->NormFactor() << endl;
      cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    }
    
    //open files and load
    cout << "LoadEvent" << endl;
    treeLoader.LoadDataset(datasets[d], anaEnv);
    cout << "LoadEvent" << endl;
    
    if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
    {
      isData = true;
    }
    
    
    /// book triggers
    if (applyTriggers) { trigger->bookTriggers(isData);}
    
    
    
    ///////////////////////////////////////////
    ///  Initialise Jet Energy Corrections  ///
    ///////////////////////////////////////////
    
    vector<JetCorrectorParameters> vCorrParam;
    string pathCalJEC = "../TopTreeAnalysisBase/Calibrations/JECFiles/";

    if (isData)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
    }
    else
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer15_25nsV6_MC_Uncertainty_AK4PFchs.txt");
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); //true means redo also L1
    
    
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    
    /// Get run information
    datasets[d]->runTree()->SetBranchStatus("runInfos*",1);
    datasets[d]->runTree()->SetBranchAddress("runInfos",&runInfos);
    
    if (verbose > 1)
      //cout << "	Loop over events " << endl;
      cout << "	Loop over events  (" << ((int)((double)datasets[d]->NofEvtsToRunOver())/((double)iReducedDataSets)) << "/" << datasets[d]->NofEvtsToRunOver() << ")" << endl;
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    //for (unsigned int ievt = 0; ievt < 1000; ievt++)
    {
      
      if ( ievt%iReducedDataSets != 0 ) { continue;}
      
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
      
//      if (ievt%1000 == 0)
//        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
      if (((int)nEvents[d])%1000 == 0)
        std::cout << "Processing the " << ((int)nEvents[d]) << "th event (" << (nEvents[d]*((double)iReducedDataSets)/((double)datasets[d]->NofEvtsToRunOver()))*100  << "%)" << flush << "\r";
        //std::cout << "Processing the " << ((int)nEvents[d]) << "th event (" << (nEvents[d]*((double)iReducedDataSets)/((double)datasets[d]->NofEvtsToRunOver()))*100  << "%)" << std::endl;
      
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
      
      datasets[d]->eventTree()->LoadTree(ievt);
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      int currentRun = event->runId();
      
      if (! isData ) {
        genjets = treeLoader.LoadGenJet(ievt,false);
        //sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      
      // BE CAREFUL: TRootGenEvent is now obsolete!
      
      
      /// Fix negative event weights for amc@nlo
      if ( nlo )
      {
        mc_baseweight1 = event->getWeight(1001)/(abs(event->originalXWGTUP()));
        mc_baseweight2 = event->getWeight(1)/(abs(event->originalXWGTUP()));
        if ( ievt < 1000 )
        {
          cout << "scale_variation:         getWeight(1001) " << setw(5) << right << event->getWeight(1001) << "; originalXWGTUP " << setw(5) << event->originalXWGTUP() << endl;
          cout << "central scale variation: getWeight(1001) " << setw(5) << right << event->getWeight(1) << "; originalXWGTUP " << setw(5) << event->originalXWGTUP() << endl;
        }
        if ( mc_baseweight2 >= 0. )
        {
          MSPlot["weightIndex"]->Fill(1, datasets[d], false, Luminosity);
        }
        if ( mc_baseweight1 >= 0. )
        {
          MSPlot["weightIndex"]->Fill(2, datasets[d], false, Luminosity);
        }
        if ( mc_baseweight2 < 0. && mc_baseweight2 != -9999. )
        {
          MSPlot["weightIndex"]->Fill(-1, datasets[d], false, Luminosity);
        }
        if ( mc_baseweight1 < 0. && mc_baseweight1 != -9999. )
        {
          MSPlot["weightIndex"]->Fill(-2, datasets[d], false, Luminosity);
        }
        if ( mc_baseweight1 == -9999 && mc_baseweight2 == -9999. )
        {
          MSPlot["weightIndex"]->Fill(0, datasets[d], false, Luminosity);
        }
        
        if ( mc_baseweight1 != -9999 && mc_baseweight2 != -9999. )
        {
          if (mc_baseweight1 == mc_baseweight2)
          {
            MSPlot["weightIndex_diffs"]->Fill(1, datasets[d], false, Luminosity);
          }
          else MSPlot["weightIndex_diffs"]->Fill(0, datasets[d], false, Luminosity);
        }
        
        MSPlot["nloWeight_good"]->Fill(mc_baseweight1, datasets[d], false, Luminosity);
        MSPlot["nloWeight_bad"]->Fill(mc_baseweight2, datasets[d], false, Luminosity);
        
      }
      
      
      
      /////////////////////////////////////
      ///  DETERMINE EVENT SCALEFACTOR  ///
      /////////////////////////////////////
      
      // scale factor for the event
      double scaleFactor = 1.;
      
      hasNegWeight = false;
      if ( nlo && mc_baseweight < 0.0 )
      {
        //scaleFactor = -1.;
        hasNegWeight = true;
      }
      
      
      /// Plot number of primary vertices before PU reweighting
      MSPlot["init_nPVs_before"]->Fill(vertex.size(), datasets[d], true, Luminosity);
      
      
      /// PU reweighting
      double lumiWeight = 1.;
      
      if ( applyPU && ! isData )
      {
        lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
        //lumiWeight = LumiWeights.ITweight( vertex.size() );
        /// Outdated syst up/down !
        // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
        // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );
      }
      
      scaleFactor = scaleFactor*lumiWeight;
      
      
      /// Plot number of primary vertices after PU reweighting
      MSPlot["init_nPVs_after"]->Fill(vertex.size(), datasets[d], true, Luminosity*lumiWeight);
      
      
      
      ////////////////////////////
      ///  Include trigger set up here when using data
      ////////////////////////////
      
      bool trigged = false;
      bool fileChanged = false;
      bool runChanged = false;
      
      /// Fill selection table before trigger
      selecTableSemiMu.Fill(d,0,scaleFactor);
      MSPlot["Selection"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
      
      
      if ( ! applyTriggers && previousFilename != currentFilename )
      {
        fileChanged = true;
        previousFilename = currentFilename;
        iFile++;
        cout << "File changed!!! => iFile = " << iFile << endl;
      }
      
      if (applyTriggers)
      {
        trigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTriggers);
        trigged = trigger->checkIfFired();
        
        if (! trigged ) { continue;}
      }
      
      /// Fill selection table after trigger
      selecTableSemiMu.Fill(d,1,scaleFactor);
      MSPlot["Selection"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
      
      
      
      //////////////////////////////////////
      ///  Jet Energy Scale Corrections  ///
      //////////////////////////////////////
      
      if (applyJER && ! isData)
      {
        jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal", false);
        //jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus", false);
        //jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus", false);
        
        /// Example how to apply JES systematics
        //jetTools->correctJetJESUnc(init_jets_corrected, "minus", 1);
        //jetTools->correctJetJESUnc(init_jets_corrected, "plus", 1);
        //cout << "JER smeared!!! " << endl;
      }
      
      
      if (applyJEC)
      {
        jetTools->correctJets(init_jets_corrected, event->fixedGridRhoFastjetAll(), isData);
      }
      
      
      /// Fill control plots
      MSPlot["init_nJets"]->Fill(init_jets_corrected.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["init_nMuons"]->Fill(init_muons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["init_nElectrons"]->Fill(init_electrons.size(), datasets[d], true, Luminosity*scaleFactor);
      
      if ( init_jets_corrected.size() > 0 )
      {
        MSPlot["init_leadingJet_pT"]->Fill(init_jets_corrected[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["init_leadingJet_eta"]->Fill(init_jets_corrected[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["init_leadingJet_CSVv2Discr"]->Fill(init_jets_corrected[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
      }
      if ( init_muons.size() > 0 )
      {
        MSPlot["init_leadingMuon_pT"]->Fill(init_muons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["init_leadingMuon_eta"]->Fill(init_muons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["init_leadingMuon_phi"]->Fill(init_muons[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["init_muon_relIso"]->Fill( ( init_muons[0]->chargedHadronIso(4) + max( 0.0, init_muons[0]->neutralHadronIso(4) + init_muons[0]->photonIso(4) - 0.5*init_muons[0]->puChargedHadronIso(4) ) ) / init_muons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["init_muon_d0"]->Fill(init_muons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
      }
      if ( init_electrons.size() > 0 )
      {
        MSPlot["init_leadingElectron_pT"]->Fill(init_electrons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["init_leadingElectron_eta"]->Fill(init_electrons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
      }
      if ( mets.size() > 0 )
      {
         MSPlot["init_met_pT"]->Fill(mets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
         MSPlot["init_met_eta"]->Fill(mets[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
         MSPlot["init_met_phi"]->Fill(mets[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
      }
      
      
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      //Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
      
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2.);
      
      vector<TRootPFJet*> selectedJets = selection.GetSelectedJets(jetPT, jetEta, true, "Tight");  // PtThr, EtaThr, applyJetID, TightLoose
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(muonPTSel, muonEtaSel, muonRelIsoSel, muonWP, "Spring15");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(electronPTSel, electronEtaSel, electronWP, "Spring15_25ns", true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased
      vector<TRootMuon*> vetoMuons = selection.GetSelectedMuons(muonPTVeto, muonEtaVeto, muonRelIsoVeto, "Loose", "Spring15");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedElectrons(electronPTVeto, electronEtaVeto, "Veto", "Spring15_25ns", true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased
      
      
      if (selectedJets.size() >= 4)
      {
        if (selectedJets[3]->Pt() < 30) selectedJets.clear();
      }
      
      
      vector<TRootMCParticle*> mcParticles;
      
      if ( dataSetName.find("TT") == 0 )
      {
        treeLoader.LoadMCEvent(ievt, 0, mcParticles, false);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      
      eventSelected = false;
      has1bjet = false;
      has2bjets = false;
      nb_bTaggedJets = 0;
      
      /// Continue with selection table
      if (isGoodPV)
      {
        selecTableSemiMu.Fill(d,2,scaleFactor);
        MSPlot["Selection"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
        if (selectedMuons.size() == 1)
        {
          selecTableSemiMu.Fill(d,3,scaleFactor);
          MSPlot["Selection"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
          /// Apply muon scale factor
          if (applyLeptonSF && ! isData )
          {
            muonSFID = muonSFWeightID_T->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);  // eta, pt, shiftUpDown
            muonSFIso = muonSFWeightIso_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);  // eta, pt, shiftUpDown
            histo2D["muon_SF_ID"]->Fill(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), muonSFID);
            histo2D["muon_SF_Iso"]->Fill(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), muonSFIso);
            scaleFactor = scaleFactor*muonSFID*muonSFIso;
          }
          if (vetoMuons.size() == 1) {
            selecTableSemiMu.Fill(d,4,scaleFactor);
            MSPlot["Selection"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
            if (vetoElectronsSemiMu.size() == 0) {
              selecTableSemiMu.Fill(d,5,scaleFactor);
              MSPlot["Selection"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
              if ( selectedJets.size() >= 4 )
              {
                selecTableSemiMu.Fill(d,6,scaleFactor);
                MSPlot["Selection"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
                eventSelected = true;
                
                for (unsigned int i = 0; i < selectedJets.size(); i++)
                {
                  if (selectedJets[i]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium) nb_bTaggedJets++;
                }
                		
                if ( nb_bTaggedJets >= 1 )
                {
                  selecTableSemiMu.Fill(d,7,scaleFactor);
                  MSPlot["Selection"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
                  has1bjet = true;
                  
                  if ( nb_bTaggedJets >= 2 )
                  {
                    selecTableSemiMu.Fill(d,8,scaleFactor);
                    MSPlot["Selection"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
                    has2bjets = true;
                  }  // at least 2 b-tagged jets
                }  // at least 1 b-tagged jets

              }  // at least 4 jets
            }  // no loose electrons
          }  // no additional loose muons (tight muon is also loose muon)
        }  // 1 good muon
      }  // good PV
      
      
//       if ( applyTriggers && ! trigged ) { continue;}
//       selecTableSemiMu.Fill(d,9,scaleFactor);
//       MSPlot["Selection"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
      
      
      /// Do some stuff with selected events
      
      if (! eventSelected )
      {
        //cout << "Event no. " << ievt << " was not selected. " << endl;
        continue;
      }
      
      nofSelectedEvents++;
      if (hasNegWeight)
      {
        nofNegWeights++;
      }
      
      if (verbose > 3)
        cout << endl << "  Event " << ievt << " is selected" << endl;
      if (verbose > 4)
        cout << "Event Id: " << event->eventId() << "  Run Id: " << event->runId() << "  Lumi block Id: " << event->lumiBlockId() << endl;
      
      
      /// Pile-up
      MSPlot["nPVs_before"]->Fill(vertex.size(), datasets[d], true, Luminosity);
      MSPlot["nPVs_after"]->Fill(vertex.size(), datasets[d], true, Luminosity*lumiWeight);
      
      
      
      ////////////////////////////
      ///  Find b-tagged jets  ///
      ////////////////////////////
      
      int label_bJet1 = -9999;
      int label_bJet2 = -9999;
      float pT_bJet1 = -9999.;
      float pT_bJet2 = -9999.;
      for (unsigned int i = 0; i < selectedJets.size(); i++)
      {
        if ( ! has1bjet ) break;
        if (selectedJets[i]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium)
        {
          if ( ! has2bjets )
          {
            label_bJet1 = i;
            pT_bJet1 = selectedJets[label_bJet1]->Pt();
            break;
          }
          else
          {
            if (selectedJets[i]->Pt() > pT_bJet1)
            {
              // Save previous as second best
              if(label_bJet1 >= 0)
              {
                label_bJet2 = label_bJet1;
                pT_bJet2 = pT_bJet1;
              }
              // Keep new one
              label_bJet1 = i;
              pT_bJet1 = selectedJets[label_bJet1]->Pt();
            }
            else if (selectedJets[i]->Pt() > pT_bJet2)
            {
              label_bJet2 = i;
              pT_bJet2 = selectedJets[label_bJet2]->Pt();
            }
          }
        }
      }
      
      
      
      ////////////////////
      ///  FILL PLOTS  ///
      ////////////////////
      
      double HT = selectedJets[0]->Pt()+selectedJets[1]->Pt()+selectedJets[2]->Pt()+selectedJets[3]->Pt();
      double relIsoMu = ( selectedMuons[0]->chargedHadronIso(4) + max( 0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4) ) ) / selectedMuons[0]->Pt();  // dR = 0.4, dBeta corrected
      
      
      MSPlot["muon_pT"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["muon_eta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["muon_phi"]->Fill(selectedMuons[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["muon_relIso"]->Fill(relIsoMu, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["muon_d0"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["leadingJet_pT"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["jet2_pT"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["jet3_pT"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["jet4_pT"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Ht_4leadingJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
      if ( mets.size() > 0 )
      {
         MSPlot["met_pT"]->Fill(mets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
         MSPlot["met_eta"]->Fill(mets[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
         MSPlot["met_phi"]->Fill(mets[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
      }
      
      MSPlot["nJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nBJets"]->Fill(nb_bTaggedJets, datasets[d], true, Luminosity*scaleFactor);
      if (has1bjet)
      {
        MSPlot["bJet1_pT"]->Fill(selectedJets[label_bJet1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["bJet1_CSVv2Discr"]->Fill(selectedJets[label_bJet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
        
        if (has2bjets)
        {
          MSPlot["bJet2_pT"]->Fill(selectedJets[label_bJet2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
          MSPlot["bJet2_CSVv2Discr"]->Fill(selectedJets[label_bJet2]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
        }
        
        MSPlot["1b_muon_pT"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["1b_muon_eta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["1b_muon_phi"]->Fill(selectedMuons[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["1b_muon_relIso"]->Fill(relIsoMu, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["1b_muon_d0"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["1b_leadingJet_pT"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["1b_jet2_pT"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["1b_jet3_pT"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["1b_jet4_pT"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["1b_Ht_4leadingJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
        
      }  /// end 1b
      
      MSPlot["pileup_SF"]->Fill(lumiWeight, datasets[d], true, Luminosity*scaleFactor);
      
      
      
      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  /// Loop on events
    
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofEventsWith1BJet << " events with 1 b tagged jet." << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofEventsWith2BJets << " events with 2 b tagged jets." << endl;
    
    if ( nlo )
      cout << "Data set " << datasets[d]->Title() << " has " << nofNegWeights << " events with negative weights." << endl;
    
    
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  /// Loop on datasets
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  //string pathPNG = "PlotsOneFourth/";
  mkdir(pathPNG.c_str(),0777);
  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  
  ///Write histograms
  fout->cd();
  for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    //cout << "MSPlot: " << it->first << endl;
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(name, 0, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
    temp->Write(fout, name, true, pathPNG+"MSPlot/", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
  }
  

  // 2D
  TDirectory* th2dir = fout->mkdir("2D_histograms");
  th2dir->cd();
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
    TH2F *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  
  
  
  ///Selection tables
  selecTableSemiMu.TableCalculator(false, true, true, true, true);
  selecTableSemiMu.Write(selectiontableMu.c_str(), true, true, true, true, true, true, false);
  
  fout->Close();
  
  delete fout;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;
  
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
  
}
