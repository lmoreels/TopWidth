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

//used TopTreeAnalysis classes
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
#include "../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "../TopTreeAnalysisBase/Tools/interface/Trigger.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/TransferFunctions.h"


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
  bool useOneTenthOfDataSets = false;
  bool useOneFiftiethOfDataSets = false;
  bool useTestSample = true;
  
  
  string rootFileName = "testAnalyser_output_FullDataSet_"+dateString+".root";
  string selectiontableMu = "SelectionTable_testFull_SemiMu_"+dateString+".tex";
  string pathPNG = "Plots_"+dateString+"/";
  int iReducedDataSets = 1;
  
  if (useOneFourthOfDataSets)
  {
    rootFileName = "testAnalyser_output_oneFourthOfDataSets_"+dateString+".root";
    selectiontableMu = "SelectionTable_testOneFourth_SemiMu_"+dateString+".tex";
    pathPNG = "PlotsOneFourth_"+dateString+"/";
    iReducedDataSets = 4;
  }
  if (useOneTenthOfDataSets)
  {
    rootFileName = "testAnalyser_output_oneTenthOfDataSets_"+dateString+".root";
    selectiontableMu = "SelectionTable_testOneTenth_SemiMu_"+dateString+".tex";
    pathPNG = "PlotsOneTenth_"+dateString+"/";
    iReducedDataSets = 10;
  }
  if (useOneFiftiethOfDataSets)
  {
    rootFileName = "testAnalyser_output_oneFiftiethOfDataSets_"+dateString+".root";
    selectiontableMu = "SelectionTable_testOneFiftieth_SemiMu_"+dateString+".tex";
    pathPNG = "PlotsOneFiftieth_"+dateString+"/";
    iReducedDataSets = 50;
  }
  if (useTestSample)
  {
    rootFileName = "testAnalyser_output_testSample_"+dateString+".root";
    selectiontableMu = "SelectionTable_testSample_SemiMu_"+dateString+".tex";
    pathPNG = "PlotsTestSample_"+dateString+"/";
    iReducedDataSets = 200;
  }
  
  
  
  /////////////////////
  ///  Configuration
  /////////////////////
  
  bool testTTbarOnly = false;  
  bool calculateTransferFunctions = true;
  bool printTriggers = false;
  bool applyTriggers = true;
  bool applyJER = true;
  bool applyJEC = true;
  bool applyLeptonSF = true;
  bool applyPU = true;
  
  bool hasNegWeight = false;
  bool eventSelected = false;
  bool has1bjet = false;
  bool has2bjets = false;
  int nofSelectedEvents = 0;
  int nofMatchedEvents = 0;
  int nb_bTaggedJets = 0;
  int nofEventsWith1BJet = 0;
  int nofEventsWith2BJets = 0;
  int nofNegWeights = 0;
  int nofPosWeights = 0;
  int nofEventsHLTv2 = 0;
  int nofEventsHLTv3 = 0;
  
  
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
  
  float Luminosity = oldLuminosity;
  //float LuminosityMu = oldLuminosity;
  //float LuminosityEl = oldLuminosity;
  
  //bool foundMu = false;
  //bool foundEl = false;
  
  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    //if ( Luminosity > datasets[d]->EquivalentLumi() ) { Luminosity = datasets[d]->EquivalentLumi();}
    string dataSetName = datasets[d]->Name();
    
    if ( (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) && Luminosity > datasets[d]->EquivalentLumi() ) { Luminosity = datasets[d]->EquivalentLumi();}
    
    //if (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
    //  LuminosityMu = datasets[d]->EquivalentLumi();
    //  foundMu=true;
    //}
    //if (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
    //  LuminosityEl = datasets[d]->EquivalentLumi();
    //  foundEl=true;
    //}
    
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
//   if ( iReducedDataSets != 1 )
//   {
//     Luminosity = Luminosity/((double) iReducedDataSets);
//     cout << "Running over 1/" << iReducedDataSets << " of the dataset, so luminosity changed to " << Luminosity << endl;
//   }
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  
  //Global variable
  //TRootEvent* event = 0;
  TRootRun *runInfos = new TRootRun();
  
  //nof selected events
  //double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  Double_t mc_baseweight = 0;
  
  
  
  ////////////////////////////////////
  ///  Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
	
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;
  
  histo1D["muon_pT"] = new TH1F("muon_pT","Transverse momentum of the muon; p_{T} [GeV]", 200, 0, 200);
  histo1D["muon_eta"] = new TH1F("muon_eta","Pseudorapidity of the muon; #eta", 60, -3, 3);
  histo1D["leadingJet_pT"] = new TH1F("leadingJet_pT","Transverse momentum of the leading jet; p_{T} [GeV]", 800, 0, 800);
  histo1D["Ht_4leadingJets"] = new TH1F("Ht_4leadingJets","Scalar sum of transverse momenta of the 4 leading jets; H_{T} [GeV]", 120, 0, 1200);
  
  histo1D["WMass_reco_matched"] = new TH1F("WMass_reco_matched","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 500, 0, 500);
  histo1D["topMass_reco_matched"] = new TH1F("topMass_reco_matched","Reconstructed top mass of matched events; M_{t} [GeV]", 500, 0, 500);
  histo1D["topMass_gen_matched"] = new TH1F("topMass_gen_matched","Generated top mass of matched events; M_{t} [GeV]", 500, 0, 500);
  histo1D["WMass_reco_first4matched"] = new TH1F("WMass_reco_first4matched","Reconstructed hadronic W mass of events where 4 hardest jets are matched; M_{W} [GeV]", 500, 0, 500);
  histo1D["topMass_reco_first4matched"] = new TH1F("topMass_reco_first4matched","Reconstructed top mass of events where 4 hardest jets are matched; M_{t} [GeV]", 500, 0, 500);
  histo1D["topMass_gen_first4matched"] = new TH1F("topMass_gen_first4matched","Generated top mass of events where partons are matched to 4 hardest jets; M_{t} [GeV]", 500, 0, 500);
  histo1D["WMass_reco_2b_notMatched"] = new TH1F("WMass_reco_2b_notMatched","Reconstructed hadronic W mass of unmatched events with 2 b-tagged jets; M_{W} [GeV]", 500, 0, 500);
  histo1D["topMass_reco_2b_notMatched"] = new TH1F("topMass_reco_2b_notMatched","Reconstructed top mass of unmatched events with 2 b-tagged jets; M_{t} [GeV]", 500, 0, 500);
  histo1D["WMass_reco_1b_notMatched"] = new TH1F("WMass_reco_1b_notMatched","Reconstructed hadronic W mass of unmatched events with 1 b-tagged jet; M_{W} [GeV]", 500, 0, 500);
  histo1D["topMass_reco_1b_notMatched"] = new TH1F("topMass_reco_1b_notMatched","Reconstructed top mass of unmatched events with 1 b-tagged jet; M_{t} [GeV]", 500, 0, 500);
  histo1D["WMass_reco_0b_notMatched"] = new TH1F("WMass_reco_0b_notMatched","Reconstructed hadronic W mass of unmatched events without b tagging; M_{W} [GeV]", 500, 0, 500);
  histo1D["topMass_reco_0b_notMatched"] = new TH1F("topMass_reco_0b_notMatched","Reconstructed top mass of unmatched events without b tagging; M_{t} [GeV]", 500, 0, 500);
  
  histo2D["logLikeWidthMass_reco_matched"] = new TH2F("logLikeWidthMass_reco_matched", "-Log Likelihood of reconstructed matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 62, 144.75, 175.75, 595, 0.55, 60.05);
  histo2D["logLikeWidthMass_gen_matched"] = new TH2F("logLikeWidthMass_gen_matched", "-Log Likelihood of generated matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 62, 144.75, 175.75, 295, 0.55, 30.05);
  histo2D["logLikeWidthMass_reco_matched_zoom"] = new TH2F("logLikeWidthMass_reco_matched_zoom", "-Log Likelihood of reconstructed matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 62, 144.75, 175.75, 400, 12.05, 52.05);
  //histo2D["logLikeWidthMass_reco"] = new TH2F("logLikeWidthMass_reco", "-Log Likelihood of reconstructed matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 10, 167.25, 172.25, 35, 0.55, 4.05);  // sample with mt = 169.5
  //histo2D["logLikeWidthMass_gen"] = new TH2F("logLikeWidthMass_gen", "-Log Likelihood of generated matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 10, 167.25, 172.25, 35, 0.55, 4.05);  // sample with mt = 169.5
  
  //histo2D["muon_SF_ID"] = new TH2F("muon_SF_ID", "Muon ID scale factors in function of #eta (x) and p_{T} (y); #eta, p_{T} [GeV]", 21, 0, 2.1, 15, 0, 150);
  //histo2D["muon_SF_Iso"] = new TH2F("muon_SF_Iso", "Muon relIso scale factors in function of #eta (x) and p_{T} (y); #eta, p_{T} [GeV]", 21, 0, 2.1, 15, 0, 150);
  //histo2D["muon_SF_Trig"] = new TH2F("muon_SF_Trig", "Muon trigger scale factors in function of #eta (x) and p_{T} (y); #eta, p_{T} [GeV]", 21, 0, 2.1, 15, 0, 150);
  
  
  
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
  MSPlot["init_muon_d0"] = new MultiSamplePlot(datasets, "init_muon_d0", 50, 0, 0.03, "d_{0}");
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
  MSPlot["muon_d0"] = new MultiSamplePlot(datasets, "muon_d0", 50, 0, 0.03, "d_{0}");
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
  MSPlot["1b_muon_d0"] = new MultiSamplePlot(datasets, "1b_muon_d0", 50, 0, 0.03, "d_{0}");
  MSPlot["1b_leadingJet_pT"] = new MultiSamplePlot(datasets, "1b_leadingJet_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["1b_jet2_pT"] = new MultiSamplePlot(datasets, "1b_jet2_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["1b_jet3_pT"] = new MultiSamplePlot(datasets, "1b_jet3_pT", 25, 0, 500, "p_{T} [GeV]");
  MSPlot["1b_jet4_pT"] = new MultiSamplePlot(datasets, "1b_jet4_pT", 25, 0, 500, "p_{T} [GeV]");
  MSPlot["1b_Ht_4leadingJets"] = new MultiSamplePlot(datasets,"1b_Ht_4leadingJets", 60, 0, 1200, "H_{T} [GeV]");
  MSPlot["1b_nJets"] = new MultiSamplePlot(datasets, "1b_nJets", 13, -0.5, 12.5, "# jets");
  MSPlot["1b_met_pT"] = new MultiSamplePlot(datasets, "1b_met_pT", 40, 0, 800, "p_{T} [GeV]");
  
  MSPlot["muon_pT_noJetCut"] = new MultiSamplePlot(datasets, "muon_pT_noJetCut", 22, 0, 440, "p_{T} [GeV]");
  MSPlot["muon_relIso_noJetCut"] = new MultiSamplePlot(datasets, "muon_relIso_noJetCut", 30, 0, 0.3, "relIso");
  MSPlot["met_pT_noJetCut"] = new MultiSamplePlot(datasets, "met_pT_noJetCut", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["nJets_noJetCut"] = new MultiSamplePlot(datasets, "nJets_noJetCut", 13, -0.5, 12.5, "# jets");
      
  MSPlot["nPVs_before"] = new MultiSamplePlot(datasets, "nPVs_before", 41, -0.5, 40.5, "# PVs before reweighting");
  MSPlot["nPVs_after"] = new MultiSamplePlot(datasets, "nPVs_after", 41, -0.5, 40.5, "# PVs after reweighting");
  
  /// Scale factors
  MSPlot["pileup_SF"] = new MultiSamplePlot(datasets,"pileup_SF", 80, 0, 4, "lumiWeight");
  MSPlot["nloWeight"] = new MultiSamplePlot(datasets,"nloWeight", 200, -2.0, 2.0, "weights for amc@nlo samples");
  MSPlot["weightIndex"] = new MultiSamplePlot(datasets,"weightIndex", 5, -2.5, 2.5, "0: None; 1: scale_variation 1; 2: Central scale variation 1");
  
  
  
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
  
  
  
  ///////////////////////////////////////
  ///  Initialise Transfer Functions  ///
  ///////////////////////////////////////
  
  TransferFunctions* tf = new TransferFunctions(calculateTransferFunctions);
  
  
  
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
  
  /// Leptons
  cout << " - Loading lepton scale factors ...";
  if (! applyLeptonSF) { cout << "    --- At the moment these are not used in the analysis";}
  cout << endl;
  
  double muonSFID, muonSFIso, muonSFTrig;
  //MuonSFWeight *muonSFWeight_ = new MuonSFWeight(pathCalLept+"Muon_SF_TopEA.root","SF_totErr", true, false, false); // (... , ... , extendRange, debug, print warning)
  MuonSFWeight *muonSFWeightID_T = new MuonSFWeight(pathCalLept+"MuonID_Z_RunCD_Reco74X_Dec1.root", "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);
  MuonSFWeight *muonSFWeightID_M = new MuonSFWeight(pathCalLept+"MuonID_Z_RunCD_Reco74X_Dec1.root", "NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);
//  MuonSFWeight *muonSFWeightID_L = new MuonSFWeight(pathCalLept+"MuonID_Z_RunCD_Reco74X_Dec1.root", "NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);
  
  MuonSFWeight *muonSFWeightIso_TT = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco74X_Dec1.root", "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
//   MuonSFWeight *muonSFWeightIso_TM = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco74X_Dec1.root", "NUM_TightRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Tight RelIso, Medium ID
//   MuonSFWeight *muonSFWeightIso_LT = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco74X_Dec1.root", "NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Loose RelIso, Tight ID
//   MuonSFWeight *muonSFWeightIso_LM = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco74X_Dec1.root", "NUM_LooseRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Loose RelIso, Medium ID
//   MuonSFWeight *muonSFWeightIso_LT = new MuonSFWeight(pathCalLept+"MuonIso_Z_RunCD_Reco74X_Dec1.root", "NUM_LooseRelIso_DEN_LooseID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Loose RelIso, Loose ID
  
  double weightMuonHLTv2, weightMuonHLTv3;
  MuonSFWeight *muonSFWeightTrigHLTv4p2 = new MuonSFWeight(pathCalLept+"SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins/abseta_pt_ratio", true, false, false);
  MuonSFWeight *muonSFWeightTrigHLTv4p3 = new MuonSFWeight(pathCalLept+"SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio", true, false, false);
  
  //ElectronSFWeight *electronSFWeight_ = new ElectronSFWeight(pathCalLept+"Elec_SF_TopEA.root","GlobalSF", true, false, false); // (... , ... , extendRange, debug, print warning)
  
  
  /// B tag
  
  
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
  
  
  
  ////////////////////////////////
  ///  Define TLorentzVectors  ///
  ////////////////////////////////
  
  // Matching
  vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
  TLorentzVector topQuark, antiTopQuark;
  // Transfer functions
  vector<TLorentzVector*> partonTLV, jetTLV;
  //TLorentzVector genMuTLV, selMuTLV, genElTLV, selElTLV;
  
  
  
  /////////////////////////////////////////////
  ///  Define variables for top propagator  ///
  /////////////////////////////////////////////
  
  //float genTopMass = 172.5;  // Check!
  //float genTopWidth = 1.3;  // Check!
  float listTopMass[] = {145.0, 145.5, 146.0, 146.5, 147.0, 147.5, 148.0, 148.5, 149.0, 149.5, 150.0, 150.5, 151.0, 151.5, 152.0, 152.5, 153.0, 153.5, 154.0, 154.5, 155.0, 155.5, 156.0, 156.5, 157.0, 157.5, 158.0, 158.5, 159.0, 159.5, 160.0, 160.5, 161.0, 161.5, 162.0, 162.5, 163.0, 163.5, 164.0, 164.5, 165.0, 165.5, 166.0, 166.5, 167.0, 167.5, 168.0, 168.5, 169.0, 169.5, 170.0, 170.5, 171.0, 171.5, 172.0, 172.5, 173.0, 173.5, 174.0, 174.5, 175.0, 175.5};
  //float listTopMass[] = {161.0, 161.5, 162.0, 162.5, 163.0, 163.5, 164.0, 164.5, 165.0, 165.5, 166.0, 166.5, 167.0, 167.5, 168.0, 168.5, 169.0, 169.5, 170.0, 170.5, 171.0, 171.5, 172.0, 172.5, 173.0, 173.5, 174.0, 174.5, 175.0, 175.5};
  float listTopWidth[] = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8, 21.9, 22.0, 22.1, 22.2, 22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9, 23.0, 23.1, 23.2, 23.3, 23.4, 23.5, 23.6, 23.7, 23.8, 23.9, 24.0, 24.1, 24.2, 24.3, 24.4, 24.5, 24.6, 24.7, 24.8, 24.9, 25.0, 25.1, 25.2, 25.3, 25.4, 25.5, 25.6, 25.7, 25.8, 25.9, 26.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 26.9, 27.0, 27.1, 27.2, 27.3, 27.4, 27.5, 27.6, 27.7, 27.8, 27.9, 28.0, 28.1, 28.2, 28.3, 28.4, 28.5, 28.6, 28.7, 28.8, 28.9, 29.0, 29.1, 29.2, 29.3, 29.4, 29.5, 29.6, 29.7, 29.8, 29.9, 30.0, 30.1, 30.2, 30.3, 30.4, 30.5, 30.6, 30.7, 30.8, 30.9, 31.0, 31.1, 31.2, 31.3, 31.4, 31.5, 31.6, 31.7, 31.8, 31.9, 32.0, 32.1, 32.2, 32.3, 32.4, 32.5, 32.6, 32.7, 32.8, 32.9, 33.0, 33.1, 33.2, 33.3, 33.4, 33.5, 33.6, 33.7, 33.8, 33.9, 34.0, 34.1, 34.2, 34.3, 34.4, 34.5, 34.6, 34.7, 34.8, 34.9, 35.0, 35.1, 35.2, 35.3, 35.4, 35.5, 35.6, 35.7, 35.8, 35.9, 36.0, 36.1, 36.2, 36.3, 36.4, 36.5, 36.6, 36.7, 36.8, 36.9, 37.0, 37.1, 37.2, 37.3, 37.4, 37.5, 37.6, 37.7, 37.8, 37.9, 38.0, 38.1, 38.2, 38.3, 38.4, 38.5, 38.6, 38.7, 38.8, 38.9, 39.0, 39.1, 39.2, 39.3, 39.4, 39.5, 39.6, 39.7, 39.8, 39.9, 40.0, 40.1, 40.2, 40.3, 40.4, 40.5, 40.6, 40.7, 40.8, 40.9, 41.0, 41.1, 41.2, 41.3, 41.4, 41.5, 41.6, 41.7, 41.8, 41.9, 42.0, 42.1, 42.2, 42.3, 42.4, 42.5, 42.6, 42.7, 42.8, 42.9, 43.0, 43.1, 43.2, 43.3, 43.4, 43.5, 43.6, 43.7, 43.8, 43.9, 44.0, 44.1, 44.2, 44.3, 44.4, 44.5, 44.6, 44.7, 44.8, 44.9, 45.0, 45.1, 45.2, 45.3, 45.4, 45.5, 45.6, 45.7, 45.8, 45.9, 46.0, 46.1, 46.2, 46.3, 46.4, 46.5, 46.6, 46.7, 46.8, 46.9, 47.0, 47.1, 47.2, 47.3, 47.4, 47.5, 47.6, 47.7, 47.8, 47.9, 48.0, 48.1, 48.2, 48.3, 48.4, 48.5, 48.6, 48.7, 48.8, 48.9, 49.0, 49.1, 49.2, 49.3, 49.4, 49.5, 49.6, 49.7, 49.8, 49.9, 50.0, 50.1, 50.2, 50.3, 50.4, 50.5, 50.6, 50.7, 50.8, 50.9, 51.0, 51.1, 51.2, 51.3, 51.4, 51.5, 51.6, 51.7, 51.8, 51.9, 52.0, 52.1, 52.2, 52.3, 52.4, 52.5, 52.6, 52.7, 52.8, 52.9, 53.0, 53.1, 53.2, 53.3, 53.4, 53.5, 53.6, 53.7, 53.8, 53.9, 54.0, 54.1, 54.2, 54.3, 54.4, 54.5, 54.6, 54.7, 54.8, 54.9, 55.0, 55.1, 55.2, 55.3, 55.4, 55.5, 55.6, 55.7, 55.8, 55.9, 56.0, 56.1, 56.2, 56.3, 56.4, 56.5, 56.6, 56.7, 56.8, 56.9, 57.0, 57.1, 57.2, 57.3, 57.4, 57.5, 57.6, 57.7, 57.8, 57.9, 58.0, 58.1, 58.2, 58.3, 58.4, 58.5, 58.6, 58.7, 58.8, 58.9, 59.0, 59.1, 59.2, 59.3, 59.4, 59.5, 59.6, 59.7, 59.8, 59.9, 60.0};
  const int sizeListTopMass = sizeof(listTopMass)/sizeof(listTopMass[0]);
  const int sizeListTopWidth = sizeof(listTopWidth)/sizeof(listTopWidth[0]);
  
  double gammaProp_matched, numTopPropagator_matched;
  double denomTopPropagator_gen_matched, topPropagator_gen_matched;
  double denomTopPropagator_reco_matched, topPropagator_reco_matched;
  double likelihood_gen_matched[sizeListTopMass][sizeListTopWidth] = {{0}};
  double likelihood_reco_matched[sizeListTopMass][sizeListTopWidth] = {{0}};
  
  
  
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  double timePerDataSet[datasets.size()] = {0};
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
  
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    clock_t startDataSet = clock();
    
    nofSelectedEvents = 0;
    nofEventsWith1BJet = 0;
    nofEventsWith2BJets = 0;
    nofNegWeights = 0;
    nofPosWeights = 0;
    float sumWeights = 0.;
    double nloSF = 1.;
    int iFile = -1;
    string previousFilename = "";
    mcParticlesTLV.clear(); selectedJetsTLV.clear();  // vector<TLV>
    partonTLV.clear(); jetTLV.clear();  // vector<TLV*>
    topQuark.Clear(); antiTopQuark.Clear(); //genMuTLV.Clear(); selMuTLV.Clear(); genElTLV.Clear(); selElTLV.Clear();
    
    bool nlo = false;
    bool isData = false;
    
    string dataSetName = datasets[d]->Name();
    if ( testTTbarOnly && dataSetName.find("TT") != 0 )
    {
      cout << "Skipping data set " << dataSetName << " ..." << endl;
      continue;
    }
    
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << "/ title : " << datasets[d]->Title() << endl;
      cout << "      -> Equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
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
    
    if ( dataSetName.find("DY") == 0 || dataSetName.find("ZJets") == 0 || dataSetName.find("Zjets") == 0 || dataSetName.find("Z+jets") == 0 || dataSetName.find("WJets") == 0 || dataSetName.find("Wjets") == 0 || dataSetName.find("W+jets") == 0 || dataSetName.find("ST_tch") == 0 )
    {
      nlo = true;
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
      
      if (isData)
      {
        if ( currentRun < 256630 )
        {
          cerr << "-- Dataset 2015C included..." << endl;
          exit(1);
        }
        else if ( currentRun >= 256630 && currentRun <= 257819 )
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
        //sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      
      
      /////////////////////////////////////
      ///  DETERMINE EVENT SCALEFACTOR  ///
      /////////////////////////////////////
      
      // scale factor for the event
      double scaleFactor = 1.;
      
      
      /// Fix negative event weights for amc@nlo
      hasNegWeight = false;
      if (! isData && dataSetName.find("ST_tW") != 0 )  // not data & not ST tW channel
      {
        if ( event->getWeight(1001) != -9999. )
        {
          mc_baseweight = event->getWeight(1001)/abs(event->originalXWGTUP());
          //mc_scaleupweight = event->getWeight(1005)/abs(event->originalXWGTUP());
          //mc_scaledownweight = event->getWeight(1009)/abs(event->originalXWGTUP());
          if ( mc_baseweight >= 0. ) 
          {
            nofPosWeights++;
            MSPlot["weightIndex"]->Fill(1., datasets[d], false, Luminosity);
          }
          else
          {
            if (nlo) hasNegWeight = true;
            nofNegWeights++;
            MSPlot["weightIndex"]->Fill(-1., datasets[d], false, Luminosity);
          }
        }
        if ( event->getWeight(1) != -9999. )
        {
          mc_baseweight = event->getWeight(1)/abs(event->originalXWGTUP());
          //mc_scaleupweight = event->getWeight(5)/abs(event->originalXWGTUP());
          //mc_scaledownweight = event->getWeight(9)/abs(event->originalXWGTUP());
          if ( mc_baseweight >= 0. )
          {
            nofPosWeights++;
            MSPlot["weightIndex"]->Fill(2., datasets[d], false, Luminosity);
          }
          else
          {
            if (nlo) hasNegWeight = true;
            nofNegWeights++;
            MSPlot["weightIndex"]->Fill(-2., datasets[d], false, Luminosity);
          }
        }
        if ( event->getWeight(1001) == -9999. && event->getWeight(1) == -9999. )
        {
          cout << "WARNING: No weight found for event " << ievt << " in dataset " << dataSetName << endl;
          cout << "         Event Id: " << event->eventId() << "  Run Id: " << event->runId() << "  Lumi block Id: " << event->lumiBlockId() << endl;
          cout << "         Weight type is different from 'scale_variation' (1001) or 'Central scale variation' (1)." << endl;
        }
        if ( event->getWeight(1001) != -9999. && event->getWeight(1) != -9999. )
        {
          cout << "WARNING: Two weight types found for event " << ievt << " in dataset " << dataSetName << endl;
          cout << "         Event Id: " << event->eventId() << "  Run Id: " << event->runId() << "  Lumi block Id: " << event->lumiBlockId() << endl;
          cout << "         Check which weight type should be used when." << endl;
        }
        
        MSPlot["nloWeight"]->Fill(mc_baseweight, datasets[d], false, Luminosity);
        sumWeights += mc_baseweight;
        
        scaleFactor *= mc_baseweight;
      }
      
      
      
      /////////////////
      ///  Pile-up  ///
      /////////////////
      
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
      MSPlot["Selection"]->Fill(1., datasets[d], true, Luminosity*scaleFactor);
      
      
      
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
      
      
//       if (selectedJets.size() >= 4)
//       {
//         if (selectedJets[3]->Pt() < 30) selectedJets.clear();
//       }
      
      
      vector<TRootMCParticle*> mcParticles;
      
      if ( dataSetName.find("TT") == 0 )
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles, false);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      
      eventSelected = false;
      has1bjet = false;
      has2bjets = false;
      nb_bTaggedJets = 0;
      muonSFID = muonSFIso = muonSFTrig = 1.;
      
      
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
            muonSFTrig = weightMuonHLTv2 * muonSFWeightTrigHLTv4p2->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0) + weightMuonHLTv3 * muonSFWeightTrigHLTv4p3->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
            
            //histo2D["muon_SF_ID"]->Fill(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), muonSFID);
            //histo2D["muon_SF_Iso"]->Fill(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), muonSFIso);
            //histo2D["muon_SF_Trig"]->Fill(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), muonSFTrig);
            scaleFactor = scaleFactor*muonSFID*muonSFIso*muonSFTrig;
          }
          if (vetoMuons.size() == 1) {
            selecTableSemiMu.Fill(d,4,scaleFactor);
            MSPlot["Selection"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
            if (vetoElectronsSemiMu.size() == 0) {
              selecTableSemiMu.Fill(d,5,scaleFactor);
              MSPlot["Selection"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
              
              double testRelIso = ( selectedMuons[0]->chargedHadronIso(4) + max( 0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4) ) ) / selectedMuons[0]->Pt();  // dR = 0.4, dBeta corrected
              MSPlot["muon_pT_noJetCut"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["muon_relIso_noJetCut"]->Fill(testRelIso, datasets[d], true, Luminosity*scaleFactor);
              MSPlot["met_pT_noJetCut"]->Fill(mets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              MSPlot["nJets_noJetCut"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
              
              /// First 4 jets need pT > 30 GeV
              if (selectedJets.size() >= 4)
              {
                if (selectedJets[3]->Pt() < 30) selectedJets.clear();
              }
              
              
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
      
      if (verbose > 3)
        cout << endl << "  Event " << ievt << " is selected" << endl;
      if (verbose > 4)
        cout << "Event Id: " << event->eventId() << "  Run Id: " << event->runId() << "  Lumi block Id: " << event->lumiBlockId() << endl;
      
      
      /// Pile-up
      MSPlot["nPVs_before"]->Fill(vertex.size(), datasets[d], true, Luminosity);
      MSPlot["nPVs_after"]->Fill(vertex.size(), datasets[d], true, Luminosity*lumiWeight);
      
      
      
      /////////////////////////////
      ///  JET PARTON MATCHING  ///
      /////////////////////////////
      
      int MCPermutation[4];
			
      bool all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
      bool all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
      bool hadronictopJetsMatched_MCdef_ = false;
      
      pair<unsigned int, unsigned int> leptonicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
      pair<unsigned int, unsigned int> hadronicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
      pair<unsigned int, unsigned int> hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
      pair<unsigned int, unsigned int> hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
      
      int pdgID_top = 6; //top quark
      
      vector<TRootMCParticle*> mcParticlesMatching_;
      int genmuon = -9999;
      bool muonmatched = false;
      
      if ( dataSetName.find("TT") == 0 )
      {
        mcParticlesTLV.clear(); selectedJetsTLV.clear();  // vector<TLV>
        topQuark.Clear(); antiTopQuark.Clear();
        
        bool muPlusFromTop = false, muMinusFromTop = false;
        mcParticlesMatching_.clear();
        
        
        for (unsigned int i = 0; i < mcParticles.size(); i++)
        {
          if (verbose > 4)
            cout << setw(3) << right << i << "  Status: " << setw(2) << mcParticles[i]->status() << "  pdgId: " << setw(3) << mcParticles[i]->type() << "  Mother: " << setw(4) << mcParticles[i]->motherType() << "  Granny: " << setw(4) << mcParticles[i]->grannyType() << "  Pt: " << setw(7) << left << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
          
          
          if ( (mcParticles[i]->status() > 1 && mcParticles[i]->status() <= 20) || mcParticles[i]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
          
          
//           if (verbose > 4 && ( (mcParticles[i]->status() > 20 && mcParticles[i]->status() < 30) 
//               || ( mcParticles[i]->status() == 1 && (abs(mcParticles[i]->type()) == 13 || abs(mcParticles[i]->type()) == 14) ) ) )
//             cout << setw(3) << right << i << "  Status: " << setw(2) << mcParticles[i]->status() << "  pdgId: " << setw(3) << mcParticles[i]->type() << "  Mother: " << setw(4) << mcParticles[i]->motherType() << "  Granny: " << setw(4) << mcParticles[i]->grannyType() << "  Pt: " << setw(7) << left << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
          
          
          if ( mcParticles[i]->type() == pdgID_top )
            topQuark = *mcParticles[i];
          else if( mcParticles[i]->type() == -pdgID_top )
            antiTopQuark = *mcParticles[i];
					
          if ( mcParticles[i]->status() == 23 && mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )		// mu-, W-, tbar
          {
            muMinusFromTop = true;
            genmuon = i;
          }
          if ( mcParticles[i]->status() == 23 && mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )		// mu+, W+, t
          {
            muPlusFromTop = true;
            genmuon = i;
	    		}
          
          if ( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 )  //light/b quarks, 6 should stay hardcoded, OR gluon
          {
            mcParticlesTLV.push_back(*mcParticles[i]);
            mcParticlesMatching_.push_back(mcParticles[i]);
          }
          
        }
        
        // take all the selectedJets_ to study the radiation stuff, selectedJets_ are already ordened in decreasing Pt()
        for (unsigned int i = 0; i < selectedJets.size(); i++)
          selectedJetsTLV.push_back(*selectedJets[i]);
        
        if (verbose > 3)
        {
          cout << "Size mcParticles:          " << mcParticles.size() << endl;
          cout << "Size mcParticlesTLV:       " << mcParticlesTLV.size() << endl;
          cout << "Size mcParticlesMatching_: " << mcParticlesMatching_.size() << endl;
          cout << "Size selectedJetsTLV:      " << selectedJetsTLV.size() << endl;
        }
        
        
        JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);		// partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
        
        if (matching.getNumberOfAvailableCombinations() != 1)
          cerr << "matching.getNumberOfAvailableCombinations() = " << matching.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
        
        
        vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
        
        for (unsigned int i = 0; i < mcParticlesTLV.size(); i++)
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
          
          if ( fabs(mcParticlesMatching_[j]->type()) < 6 )  //light/b quarks, 6 should stay hardcoded
          {
            if ( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -pdgID_top )
                || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == pdgID_top ) )  // if mu+, check if mother of particle is W- and granny tbar --> then it is a quark from W- decay
            {
              if (verbose > 3)
                cout << "Light jet: " << j << "  Status: " << mcParticlesMatching_[j]->status() << "  pdgId: " << mcParticlesMatching_[j]->type() << "  Mother: " << mcParticlesMatching_[j]->motherType() << "  Granny: " << mcParticlesMatching_[j]->grannyType() << "  Pt: " << mcParticlesMatching_[j]->Pt() << "  Eta: " << mcParticlesMatching_[j]->Eta() << "  Phi: " << mcParticlesMatching_[j]->Phi() << "  Mass: " << mcParticlesMatching_[j]->M() << endl;
              if (hadronicWJet1_.first == 9999)
              {
                hadronicWJet1_ = JetPartonPair[i];
                MCPermutation[0] = JetPartonPair[i].first;
              }
              else if (hadronicWJet2_.first == 9999)
              {
                hadronicWJet2_ = JetPartonPair[i];
                MCPermutation[1] = JetPartonPair[i].first;
              }
              else
              {
                cerr << "Found a third jet coming from a W boson which comes from a top quark..." << endl;
                cerr << " -- muMinusFromMtop: " << muMinusFromTop << " muPlusFromMtop: " << muPlusFromTop << endl;
                cerr << " -- pdgId: " << mcParticlesMatching_[j]->type() << " mother: " << mcParticlesMatching_[j]->motherType() << " granny: " << mcParticlesMatching_[j]->grannyType() << " Pt: " << mcParticlesMatching_[j]->Pt() << endl;
                cerr << " -- ievt: " << ievt << endl;
                exit(1);
              }
            }
          }
          if ( fabs(mcParticlesMatching_[j]->type()) == 5 )
          {
            if ( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -pdgID_top )
                || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == pdgID_top ) )  // if mu+ (top decay leptonic) and mother is antitop ---> hadronic b
            {
              if (verbose > 3)
                cout << "b jet:     " << j << "  Status: " << mcParticlesMatching_[j]->status() << "  pdgId: " << mcParticlesMatching_[j]->type() << "  Mother: " << mcParticlesMatching_[j]->motherType() << "  Granny: " << mcParticlesMatching_[j]->grannyType() << "  Pt: " << mcParticlesMatching_[j]->Pt() << "  Eta: " << mcParticlesMatching_[j]->Eta() << "  Phi: " << mcParticlesMatching_[j]->Phi() << "  Mass: " << mcParticlesMatching_[j]->M() << endl;
              hadronicBJet_ = JetPartonPair[i];
              MCPermutation[2] = JetPartonPair[i].first;
            }
            else if ( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == pdgID_top )
              || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == -pdgID_top ) )
            {
              if (verbose > 3)
                cout << "b jet:     " << j << "  Status: " << mcParticlesMatching_[j]->status() << "  pdgId: " << mcParticlesMatching_[j]->type() << "  Mother: " << mcParticlesMatching_[j]->motherType() << "  Granny: " << mcParticlesMatching_[j]->grannyType() << "  Pt: " << mcParticlesMatching_[j]->Pt() << "  Eta: " << mcParticlesMatching_[j]->Eta() << "  Phi: " << mcParticlesMatching_[j]->Phi() << "  Mass: " << mcParticlesMatching_[j]->M() << endl;
              leptonicBJet_ = JetPartonPair[i];
              MCPermutation[3] = JetPartonPair[i].first;
            }
          }
        }  /// End loop over Jet Parton Pairs
        
        
        if (hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999 && hadronicBJet_.first != 9999 && leptonicBJet_.first != 9999)
        {
          
          all4PartonsMatched = true;
          nofMatchedEvents++;
          if (hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4 && leptonicBJet_.first < 4)
            all4JetsMatched_MCdef_ = true;
	  		}
        else if (verbose > 3) cout << "Size JetPartonPair: " << JetPartonPair.size() << ". Not all partons matched!" << endl;
        
        if (hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4)
          hadronictopJetsMatched_MCdef_ = true;
        if (genmuon != -9999 && ROOT::Math::VectorUtil::DeltaR( (TLorentzVector)*mcParticles[genmuon], (TLorentzVector)*selectedMuons[0]) < 0.1)
          muonmatched = true;
        
        
        
        ///////////////////
        ///  Transfer functions
        ///////////////////
        
        if (all4PartonsMatched && calculateTransferFunctions)
        {
          partonTLV.clear(); jetTLV.clear();  // vector<TLV*>
          //genMuTLV.Clear(); selMuTLV.Clear(); genElTLV.Clear(); selElTLV.Clear();
          
          for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
          {
            /// JetPartonPair[i].first  = jet number
            /// JetPartonPair[i].second = parton number
            /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
            
            partonTLV.push_back(mcParticles[JetPartonPair[iMatch].second]);
            jetTLV.push_back(selectedJets[JetPartonPair[iMatch].first]);
          }
          
          tf->fillJets(partonTLV, jetTLV);
          
          if (muonmatched) tf->fillMuon((TLorentzVector*) mcParticles[genmuon], (TLorentzVector*) selectedMuons[0]);
          //if (electronmatched) tf->fillElectron(...)
          
        }  // end tf
        
        
      }  /// End matching
      
      
      
      //////////////////////////////////
      ///  TOP PROPAGATOR (MATCHED)  ///
      //////////////////////////////////
      
      if ( dataSetName.find("TT") == 0 && all4PartonsMatched )
      {
        /// MCPermutation = JetPartonPair[i].first  = jet number
        ///                 JetPartonPair[i].second = parton number
        /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
        
        
        float WMassReco_matched = (*selectedJets[MCPermutation[0]] + *selectedJets[MCPermutation[1]]).M();
        float topMassReco_matched = (*selectedJets[MCPermutation[0]] + *selectedJets[MCPermutation[1]] + *selectedJets[MCPermutation[2]]).M();
        float topMassGen_matched = (*mcParticlesMatching_[hadronicWJet1_.second] + *mcParticlesMatching_[hadronicWJet2_.second] + *mcParticlesMatching_[hadronicBJet_.second]).M();
        
        for (unsigned int jMass = 0; jMass < sizeListTopMass; jMass++)
        {
          for (unsigned int jWidth = 0; jWidth < sizeListTopWidth; jWidth++)
          {
            gammaProp_matched = sqrt( pow( listTopMass[jMass], 4 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 ) );
            numTopPropagator_matched = ( 2 * sqrt(2) * listTopMass[jMass] * listTopWidth[jWidth] * gammaProp_matched ) / ( TMath::Pi() * sqrt( pow(listTopMass[jMass], 2) + gammaProp_matched ) );
            
            /// Generated mass
            denomTopPropagator_gen_matched = pow( pow(topMassGen_matched, 2) - pow(listTopMass[jMass], 2), 2 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 );
            
            topPropagator_gen_matched = numTopPropagator_matched/denomTopPropagator_gen_matched;
            
            likelihood_gen_matched[jMass][jWidth] += -TMath::Log10(topPropagator_gen_matched);
            
            /// Reconstructed mass
            denomTopPropagator_reco_matched = pow( pow(topMassReco_matched, 2) - pow(listTopMass[jMass], 2), 2 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 );
            
            topPropagator_reco_matched = numTopPropagator_matched/denomTopPropagator_reco_matched;
            
            likelihood_reco_matched[jMass][jWidth] += -TMath::Log10(topPropagator_reco_matched);
            
          }  /// End loop jWidth
        }  /// End loop jMass
        
        
        /// Fill plots
        histo1D["WMass_reco_matched"]->Fill(WMassReco_matched);
        histo1D["topMass_reco_matched"]->Fill(topMassReco_matched);
        histo1D["topMass_gen_matched"]->Fill(topMassGen_matched);
        if ( all4JetsMatched_MCdef_ )
        {
          histo1D["WMass_reco_first4matched"]->Fill(WMassReco_matched);
          histo1D["topMass_reco_first4matched"]->Fill(topMassReco_matched);
          histo1D["topMass_gen_first4matched"]->Fill(topMassGen_matched);
        }
        
      }  // end TT && matched
      
      
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
      
      
      
      //////////////////////
      // CHI2 FOR 2 BTAGS //				// !!! chi2WMass, sigmaChi2WMass, chi2TopMass, sigmaChi2TopMass
      //////////////////////
      
      int labelsReco[4] = {-9999,-9999,-9999,-9999};		// 0 = leptonic b-jet, 1 = hadronic b-jet, 2,3 = light jets.
      float recoWMass, recoTopMass, chi2_mass;
      if ( label_bJet1 != -9999 && label_bJet2 != -9999 )
      {
        nofEventsWith2BJets++;
        float recoTopMass_bJet1, recoTopMass_bJet2, WTerm, topTerm_bJet1, topTerm_bJet2, chi2_bJet1, chi2_bJet2;
        float smallestChi2 = 999999.;
        for (int ijet = 0; ijet < 4; ijet++)
        {
          for (int jjet = ijet+1; jjet < 4; jjet++)
          {
            if ( ijet != label_bJet1 && ijet != label_bJet2 && jjet != label_bJet1 && jjet != label_bJet2 )
            {
              recoWMass = ( *selectedJets[ijet] + *selectedJets[jjet]).M();
              recoTopMass_bJet1 = (*selectedJets[ijet] + *selectedJets[jjet] + *selectedJets[label_bJet1]).M();
              recoTopMass_bJet2 = (*selectedJets[ijet] + *selectedJets[jjet] + *selectedJets[label_bJet2]).M();
              
//               WTerm = pow( (recoWMass - chi2Wmass)/sigmaChi2Wmass, 2);
//               topTerm_bJet1 = pow( (recoTopMass_bJet1 - chi2Topmass)/sigmaChi2Topmass, 2);
//               topTerm_bJet2 = pow( (recoTopMass_bJet2 - chi2Topmass)/sigmaChi2Topmass, 2);
//               
//               chi2_bJet1 = WTerm + topTerm_bJet1;
//               chi2_bJet2 = WTerm + topTerm_bJet2;
//               
//               
//               if (chi2_bJet1 < smallestChi2) {
//                 smallestChi2 = chi2_bJet1;
//                 labelsReco[0] = label_bJet2;
//                 labelsReco[1] = label_bJet1;
//                 labelsReco[2] = ijet;
//                 labelsReco[3] = jjet;
//                 recoTopMass = recoTopMass_bJet1;
//               }
//               if (chi2_bJet2 < smallestChi2) {
//                 smallestChi2 = chi2_bJet2;
//                 labelsReco[0] = label_bJet1;
//                 labelsReco[1] = label_bJet2;
//                 labelsReco[2] = ijet;
//                 labelsReco[3] = jjet;
//                 recoTopMass = recoTopMass_bJet2;
//               }
//               chi2_mass = smallestChi2;
            }
          }
        }
        
        
 	  		//Fill histos
// 	  		if (labelsReco[0] != -9999 && labelsReco[1] != -9999 && labelsReco[2] != -9999 && labelsReco[3] != -9999)
//        {
 	    		//if (useMassesAndResolutions && eventselectedSemiMu) selecTableSemiMu.Fill(d,11,scaleFactor);
 	    		//if (useMassesAndResolutions && eventselectedSemiEl) selecTableSemiEl.Fill(d,12,scaleFactor);
          
          if ( dataSetName.find("TT") == 0 )
          {
            histo1D["WMass_reco_2b_notMatched"]->Fill(recoWMass);
            histo1D["topMass_reco_2b_notMatched"]->Fill(recoTopMass_bJet1, 0.5);
            histo1D["topMass_reco_2b_notMatched"]->Fill(recoTopMass_bJet2, 0.5);
          }
// 
// 	    		MSPlot["Chi2_2btags"]->Fill(smallestChi2, datasets[d], true, Luminosity*scaleFactor);
// 
// 					float Wmass_2btags = (*selectedJets[labelsReco[2]] + *selectedJets[labelsReco[3]]).M();
// 					MSPlot["W_Mass_2btags"+Flav]->Fill(Wmass_2btags, datasets[d], true, Luminosity*scaleFactor);
// 	    		float HtTop_2btags = selectedJets[labelsReco[1]]->Pt() + selectedJets[labelsReco[2]]->Pt() + selectedJets[labelsReco[3]]->Pt();
// 	    		MSPlot["hadTop_Ht_2btags"+Flav]->Fill(HtTop_2btags, datasets[d], true, Luminosity*scaleFactor);
// 	    		float hadtopmass_2btags = ( *selectedJets[labelsReco[1]] + *selectedJets[labelsReco[2]] +  *selectedJets[labelsReco[3]]).M();
// 	    		float hadtoppt_2btags = ( *selectedJets[labelsReco[1]] + *selectedJets[labelsReco[2]] +  *selectedJets[labelsReco[3]]).Pt();
// 	    		MSPlot["hadTop_Mass_2btags"+Flav]->Fill(hadtopmass_2btags, datasets[d], true, Luminosity*scaleFactor);
// 	    		MSPlot["hadTop_Pt_2btags"+Flav]->Fill(hadtoppt_2btags, datasets[d], true, Luminosity*scaleFactor);
// 	  		}
// 	  		else
//        {
// 	    		//When no Chi2 combination is found:
// 	    		cout << "Eventnr. " << ievt << ": no Chi2 found." << endl;
// 	  		}
      
      
      }  // end 2 b tags
      
      
      
      /////////////////////
      // CHI2 FOR 1 BTAG //
      /////////////////////
      
      else if ( label_bJet1 != -9999 && label_bJet2 == -9999 )
      {
        nofEventsWith1BJet++;
        float WTerm, topTerm, chi2;
        float smallestChi2 = 999999.;
        for (int ijet = 0; ijet < 4; ijet++)
        {
          for (int jjet = ijet+1; jjet < 4; jjet++)
          {
            if ( ijet != label_bJet1 && jjet != label_bJet1 )
            {
              recoWMass = ( *selectedJets[ijet] + *selectedJets[jjet]).M();
              recoTopMass = (*selectedJets[ijet] + *selectedJets[jjet] + *selectedJets[label_bJet1]).M();
              
//               WTerm = pow( (recoWMass - chi2Wmass)/sigmaChi2Wmass, 2);
//               topTerm = pow( (recoTopMass_bJet1 - chi2Topmass)/sigmaChi2Topmass, 2);
//               
//               chi2 = WTerm + topTerm;
//               
//               
//               if (chi2 < smallestChi2) {
//                 smallestChi2 = chi2;
//                 //labelsReco[0] = ;
//                 labelsReco[1] = label_bJet1;
//                 labelsReco[2] = ijet;
//                 labelsReco[3] = jjet;
//               }
//               chi2_mass = smallestChi2;
            }
          }
        }
        
        //Fill histos
// 	  		if (labelsReco[0] != -9999 && labelsReco[1] != -9999 && labelsReco[2] != -9999 && labelsReco[3] != -9999)
//        { 
          if ( dataSetName.find("TT") == 0 )
          {
            histo1D["WMass_reco_1b_notMatched"]->Fill(recoWMass);
            histo1D["topMass_reco_1b_notMatched"]->Fill(recoTopMass);
          }
// 	  		}
        
      }  // end 1 b tag
      
      
      
      ///////////////////////
      // CHI2 FOR NO BTAGS //
      ///////////////////////
      
      else
      {
        float WTerm, topTerm, chi2;
        float smallestChi2 = 999999.;
        for (int ijet = 0; ijet < 4; ijet++)
        {
          for (int jjet = ijet+1; jjet < 4; jjet++)
          {
            for (int kjet = 0; kjet < 4; kjet++)
            {
              if ( ijet != kjet && jjet != kjet )
              {
                recoWMass = ( *selectedJets[ijet] + *selectedJets[jjet]).M();
                recoTopMass = (*selectedJets[ijet] + *selectedJets[jjet] + *selectedJets[kjet]).M();
                
//                 WTerm = pow( (recoWMass - chi2Wmass)/sigmaChi2Wmass, 2);
//                 topTerm = pow( (recoTopMass - chi2Topmass)/sigmaChi2Topmass, 2);
// 
//                 chi2 = WTerm + topTerm;
// 
// 
//                 if (chi2 < smallestChi2) {
//                   smallestChi2 = chi2;
//                   //labelsReco[0] = ;
//                   labelsReco[1] = kjet;
//                   labelsReco[2] = ijet;
//                   labelsReco[3] = jjet;
//                 }
//                 chi2_mass = smallestChi2;
              }
            }
          }
        }
        
        //Fill histos
// 	  		if (labelsReco[0] != -9999 && labelsReco[1] != -9999 && labelsReco[2] != -9999 && labelsReco[3] != -9999)
//        { 
          if ( dataSetName.find("TT") == 0 )
          {
            histo1D["WMass_reco_0b_notMatched"]->Fill(recoWMass);
            histo1D["topMass_reco_0b_notMatched"]->Fill(recoTopMass);
          }
// 	  		}
//       histo1D["WMass_reco_0b_notMatched"]
//       histo1D["topMass_reco_0b_notMatched"]
         
      }  // end no b tags
      
      
      
      ////////////////////
      ///  FILL PLOTS  ///
      ////////////////////
      
      double HT = selectedJets[0]->Pt()+selectedJets[1]->Pt()+selectedJets[2]->Pt()+selectedJets[3]->Pt();
      double relIsoMu = ( selectedMuons[0]->chargedHadronIso(4) + max( 0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4) ) ) / selectedMuons[0]->Pt();  // dR = 0.4, dBeta corrected
      
      if ( dataSetName.find("TT") == 0 )
      {
        histo1D["muon_pT"]->Fill(selectedMuons[0]->Pt());
        histo1D["muon_eta"]->Fill(selectedMuons[0]->Eta());
        histo1D["leadingJet_pT"]->Fill(selectedJets[0]->Pt());
        histo1D["Ht_4leadingJets"]->Fill(HT);
      }
      
      
      /// Fill MSPlots
      
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
        MSPlot["1b_nJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
        if ( mets.size() > 0 ) MSPlot["1b_met_pT"]->Fill(mets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        
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
    if ( dataSetName.find("TT") == 0 )
      cout << "Number of matched events: " << nofMatchedEvents << endl;
    
    //if (nlo)
    if (! isData && dataSetName.find("ST") != 0 )  // not data & not ST tW channel
    {
      cout << "Data set " << datasets[d]->Title() << " has " << nofPosWeights << " events with positive weights and " << nofNegWeights << " events with negative weights." << endl;
      cout << "         Pos - neg is " << nofPosWeights - nofNegWeights << ", pos + neg is " << nofPosWeights + nofNegWeights << endl;
      cout << "The sum of the weights is " << ((int)sumWeights) << ", whereas the total number of events is " << ((int)nEvents[d]) << endl;
      
      /// Determine scale factor due to negative weights
      nloSF = ((double) (nofPosWeights - nofNegWeights))/((double) (nofPosWeights + nofNegWeights));
      cout << "This corresponds to an event scale factor of " << nloSF << " (NB: not necessarily representative for whole sample, this is skimmed set!)" << endl;
    }
    
    if (isData)
    {
      weightMuonHLTv2 = ((double) nofEventsHLTv2) / ((double) (nofEventsHLTv2 + nofEventsHLTv3));
      weightMuonHLTv3 = ((double) nofEventsHLTv3) / ((double) (nofEventsHLTv2 + nofEventsHLTv3));
      cout << "The muon trigger scale factors will be scaled by " << weightMuonHLTv2 << " for HLTv2 and " << weightMuonHLTv3 << " for HLTv3." << endl;
    }
    cout << endl;
    
    
    /// Fill histogram log likelihood && Transfer functions
    if ( dataSetName.find("TT") == 0 )
    {
      for (unsigned int jMass = 0; jMass < sizeListTopMass; jMass++)
      {
        for (unsigned int jWidth = 0; jWidth < sizeListTopWidth; jWidth++)
        {
          histo2D["logLikeWidthMass_reco_matched"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_reco_matched[jMass][jWidth]);
          histo2D["logLikeWidthMass_reco_matched_zoom"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_reco_matched[jMass][jWidth]);
          histo2D["logLikeWidthMass_gen_matched"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_gen_matched[jMass][jWidth]);
        }
      }
      
      
      /// Transfer functions
      TFile *foutTF = new TFile("PlotsForTransferFunctions.root", "RECREATE");
      foutTF->cd();
      
      tf->writeOutputFiles();
      
      foutTF->Close();
      delete foutTF;
      
    }  // end TT
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
    timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
    
  }  /// Loop on datasets
  
  
  cout << "Processing time per dataset: " << endl;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
  }
  
  
  
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
    temp->Draw(name, 1, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
    temp->Write(fout, name, true, pathPNG+"MSPlot/", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
  }
  
  // 1D
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  th1dir->cd();
  for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
    //TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    //tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
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
