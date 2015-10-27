////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////



#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <map>
#include <TH1.h>
#include <TH2.h>
#include <array>
#include <vector>
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
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"


using namespace std;
using namespace reweight;
using namespace TopTree;


int main (int argc, char *argv[])
{
  
  string rootFileName = "testAnalyser_output.root";
  
  clock_t start = clock();
  
  
  /////////////////////
  ///  Configuration
  /////////////////////
  
  bool eventSelected = false;
  bool has1bjet = false;
  bool has2bjets = false;
  int nofSelectedEvents = 0;
  int nofMatchedEvents = 0;
  int nb_bTaggedJets = 0;
  int nofEventsWith1BJet = 0;
  int nofEventsWith2BJets = 0;
  
  
  /// xml file
  string xmlFileName ="config/test.xml";
  
  if (argc > 1)
    xmlFileName = (string)argv[1];
  
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
    if (Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
    
    //if (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
    //  LuminosityMu = datasets[d]->EquivalentLumi();
    //  foundMu=true;
    //}
    //if (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
    //  LuminosityEl = datasets[d]->EquivalentLumi();
    //  foundEl=true;
    //}
    
    if ( dataSetName.find("QCD") == 0 ) datasets[d]->SetColor(kYellow);
    if ( dataSetName.find("TT") == 0 ) datasets[d]->SetColor(kRed+1);
    //if ( dataSetName.find("TTbarJets_Other") == 0 ) datasets[d]->SetColor(kRed-7);
    if ( dataSetName.find("WJets") == 0 )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      datasets[d]->SetColor(kGreen-3);
    }
    if ( dataSetName.find("ZJets") == 0 )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      datasets[d]->SetColor(kMagenta);
    }
    if ( dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") ==0 )
      datasets[d]->SetColor(kBlue-2);
    //if (dataSetName.find("NP") == 0 )
    //{
    //	datasets[d]->SetTitle("Signal");
    //	datasets[d]->SetColor(kGreen+4);
    //}
  }
  
  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
  
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  
  //Global variable
  //TRootEvent* event = 0;
  
  //nof selected events
  //double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  
  
  ////////////////////////////////////
  ///  Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
	
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;
  
  histo1D["muon_pT"] = new TH1F("muon_pT","Transverse momentum of the muon; p_{T} [GeV]", 200, 0, 200);
  histo1D["muon_eta"] = new TH1F("muon_eta","Pseudorapidity of the muon; #eta", 60, -3, 3);
  histo1D["leadingJet_pT"] = new TH1F("leadingJet_pT","Transverse momentum of the leading jet; p_{T} [GeV]", 800, 0, 800);
  histo1D["Ht_4leadingJets"] = new TH1F("Ht_4leadingJets","Scalar sum of transverse momenta of the 4 leading jets; H_{T} [GeV]", 120, 0, 1200);
  
  histo1D["W_Mass_Reco_matched"] = new TH1F("W_Mass_Reco_matched","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 500, 0, 500);
  histo1D["top_Mass_Reco_matched"] = new TH1F("top_Mass_Reco_matched","Reconstructed top mass of matched events; M_{t} [GeV]", 500, 0, 500);
  histo1D["top_Mass_Gen_matched"] = new TH1F("top_Mass_Gen_matched","Generated top mass of matched events; M_{t} [GeV]", 500, 0, 500);
  histo1D["W_Mass_Reco_first4matched"] = new TH1F("W_Mass_Reco_first4matched","Reconstructed hadronic W mass of events where 4 hardest jets are matched; M_{W} [GeV]", 500, 0, 500);
  histo1D["top_Mass_Reco_first4matched"] = new TH1F("top_Mass_Reco_first4matched","Reconstructed top mass of events where 4 hardest jets are matched; M_{t} [GeV]", 500, 0, 500);
  histo1D["top_Mass_Gen_first4matched"] = new TH1F("top_Mass_Gen_first4matched","Generated top mass of events where partons are matched to 4 hardest jets; M_{t} [GeV]", 500, 0, 500);
  histo1D["W_Mass_Reco_notMatched"] = new TH1F("W_Mass_Reco_notMatched","Reconstructed hadronic W mass of unmatched events; M_{W} [GeV]", 500, 0, 500);
  histo1D["top_Mass_Reco_notMatched"] = new TH1F("top_Mass_Reco_notMatched","Reconstructed top mass of unmatched events; M_{t} [GeV]", 500, 0, 500);
  
  histo2D["LogLikeWidthMass_Reco"] = new TH2F("LogLikeWidthMass_Reco", "-Log Likelihood of reconstructed matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 60, 144.75, 174.75, 595, 0.55, 60.05);
  histo2D["LogLikeWidthMass_Gen"] = new TH2F("LogLikeWidthMass_Gen", "-Log Likelihood of generated matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 60, 144.75, 174.75, 295, 0.55, 30.05);
  histo2D["LogLikeWidthMass_Reco_zoom"] = new TH2F("LogLikeWidthMass_Reco_zoom", "-Log Likelihood of reconstructed matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 60, 144.75, 174.75, 400, 12.05, 52.05);
  //histo2D["LogLikeWidthMass_Reco"] = new TH2F("LogLikeWidthMass_Reco", "-Log Likelihood of reconstructed matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 10, 167.25, 172.25, 35, 0.55, 4.05);  // sample with mt = 169.5
  //histo2D["LogLikeWidthMass_Gen"] = new TH2F("LogLikeWidthMass_Gen", "-Log Likelihood of generated matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 10, 167.25, 172.25, 35, 0.55, 4.05);  // sample with mt = 169.5
  
  
  
  ////////////////////////////////////
  ///  MultiSamplePlot
  ////////////////////////////////////
  
  map<string,MultiSamplePlot*> MSPlot;
  
  MSPlot["muon_pT"] = new MultiSamplePlot(datasets, "muon_pT", 22, 0, 440, "p_{T} [GeV]");
  MSPlot["muon_eta"] = new MultiSamplePlot(datasets, "muon_eta", 60, -3, 3, "Eta");
  MSPlot["leadingJet_pT"] = new MultiSamplePlot(datasets, "leadingJet_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["jet2_pT"] = new MultiSamplePlot(datasets, "jet2_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["jet3_pT"] = new MultiSamplePlot(datasets, "jet3_pT", 25, 0, 500, "p_{T} [GeV]");
  MSPlot["jet4_pT"] = new MultiSamplePlot(datasets, "jet4_pT", 25, 0, 500, "p_{T} [GeV]");
  MSPlot["Ht_4leadingJets"] = new MultiSamplePlot(datasets,"Ht_4leadingJets", 120, 0, 1200, "H_{T} [GeV]");
  
  MSPlot["nJets"] = new MultiSamplePlot(datasets, "nJets", 11, -0.5, 10.5, "# jets");
  MSPlot["nBJets"] = new MultiSamplePlot(datasets, "nBJets", 11, -0.5, 10.5, "# b jets");
  MSPlot["bJet1_pT"] = new MultiSamplePlot(datasets, "bJet1_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["bJet2_pT"] = new MultiSamplePlot(datasets, "bJet2_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["bJet1_CSVv2Discr"] = new MultiSamplePlot(datasets, "bJet1_CSVv2Discr", 80, -0.5, 1.5, "CSVv2 discriminant value");
  MSPlot["bJet2_CSVv2Discr"] = new MultiSamplePlot(datasets, "bJet2_CSVv2Discr", 80, -0.5, 1.5, "CSVv2 discriminant value");
  
  
  
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
  CutsSelecTableSemiMu.push_back("$\\geq$ 1 b-jet (CSVM)");
  CutsSelecTableSemiMu.push_back("$\\geq$ 2 b-jets (CSVM)");
  
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  //selecTableSemiMu.SetLuminosity(LuminosityMu);
  selecTableSemiMu.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
  
  
  ///////////////////////////////
  ///  Single Muon Selection  ///
  ///////////////////////////////
  
  /// Updated 27/10/15, https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO
  
  float muonPTSel = 26; // GeV
  float muonEtaSel = 2.1;
  float muonRelIsoSel = 0.15;  // Tight muon
  
  float muonPTVeto = 10; // GeV
  float muonEtaVeto = 2.5;
  float muonRelIsoVeto = 0.25;  // Loose muon
  
  
  
  ///////////////////////////////////
  ///  Single Electron Selection  ///
  ///////////////////////////////////
  
  // To do
  float electronEtaSel = 2.5;
  float electronEtaVeto = 2.5;
  
  
  
  ///////////////////////
  ///  Jet Selection  ///
  ///////////////////////
  
  // To do
  
  
  
  //////////////////////////////////////
  ///  Working points for b tagging  ///
  //////////////////////////////////////
  
  /// Updated 27/10/15, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
  
  float CSVMv2Loose =  0.605;
  float CSVMv2Medium = 0.890;
  float CSVMv2Tight = 0.970;
  
  
  
  /////////////////////////////////////////////
  ///  Define variables for top propagator  ///
  /////////////////////////////////////////////
  
  //float genTopMass = 172.5;  // Check!
  //float genTopWidth = 1.3;  // Check!
  float listTopMass[] = {145.0, 145.5, 146.0, 146.5, 147.0, 147.5, 148.0, 148.5, 149.0, 149.5, 150.0, 150.5, 151.0, 151.5, 152.0, 152.5, 153.0, 153.5, 154.0, 154.5, 155.0, 155.5, 156.0, 156.5, 157.0, 157.5, 158.0, 158.5, 159.0, 159.5, 160.0, 160.5, 161.0, 161.5, 162.0, 162.5, 163.0, 163.5, 164.0, 164.5, 165.0, 165.5, 166.0, 166.5, 167.0, 167.5, 168.0, 168.5, 169.0, 169.5, 170.0, 170.5, 171.0, 171.5, 172.0, 172.5, 173.0, 173.5, 174.0, 174.5};
  //float listTopMass[] = {167.5, 168.0, 168.5, 169.0, 169.5, 170.0, 170.5, 171.0, 171.5, 172.0};  // sample with mt = 169.5
  float listTopWidth[] = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8, 21.9, 22.0, 22.1, 22.2, 22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9, 23.0, 23.1, 23.2, 23.3, 23.4, 23.5, 23.6, 23.7, 23.8, 23.9, 24.0, 24.1, 24.2, 24.3, 24.4, 24.5, 24.6, 24.7, 24.8, 24.9, 25.0, 25.1, 25.2, 25.3, 25.4, 25.5, 25.6, 25.7, 25.8, 25.9, 26.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 26.9, 27.0, 27.1, 27.2, 27.3, 27.4, 27.5, 27.6, 27.7, 27.8, 27.9, 28.0, 28.1, 28.2, 28.3, 28.4, 28.5, 28.6, 28.7, 28.8, 28.9, 29.0, 29.1, 29.2, 29.3, 29.4, 29.5, 29.6, 29.7, 29.8, 29.9, 30.0, 30.1, 30.2, 30.3, 30.4, 30.5, 30.6, 30.7, 30.8, 30.9, 31.0, 31.1, 31.2, 31.3, 31.4, 31.5, 31.6, 31.7, 31.8, 31.9, 32.0, 32.1, 32.2, 32.3, 32.4, 32.5, 32.6, 32.7, 32.8, 32.9, 33.0, 33.1, 33.2, 33.3, 33.4, 33.5, 33.6, 33.7, 33.8, 33.9, 34.0, 34.1, 34.2, 34.3, 34.4, 34.5, 34.6, 34.7, 34.8, 34.9, 35.0, 35.1, 35.2, 35.3, 35.4, 35.5, 35.6, 35.7, 35.8, 35.9, 36.0, 36.1, 36.2, 36.3, 36.4, 36.5, 36.6, 36.7, 36.8, 36.9, 37.0, 37.1, 37.2, 37.3, 37.4, 37.5, 37.6, 37.7, 37.8, 37.9, 38.0, 38.1, 38.2, 38.3, 38.4, 38.5, 38.6, 38.7, 38.8, 38.9, 39.0, 39.1, 39.2, 39.3, 39.4, 39.5, 39.6, 39.7, 39.8, 39.9, 40.0, 40.1, 40.2, 40.3, 40.4, 40.5, 40.6, 40.7, 40.8, 40.9, 41.0, 41.1, 41.2, 41.3, 41.4, 41.5, 41.6, 41.7, 41.8, 41.9, 42.0, 42.1, 42.2, 42.3, 42.4, 42.5, 42.6, 42.7, 42.8, 42.9, 43.0, 43.1, 43.2, 43.3, 43.4, 43.5, 43.6, 43.7, 43.8, 43.9, 44.0, 44.1, 44.2, 44.3, 44.4, 44.5, 44.6, 44.7, 44.8, 44.9, 45.0, 45.1, 45.2, 45.3, 45.4, 45.5, 45.6, 45.7, 45.8, 45.9, 46.0, 46.1, 46.2, 46.3, 46.4, 46.5, 46.6, 46.7, 46.8, 46.9, 47.0, 47.1, 47.2, 47.3, 47.4, 47.5, 47.6, 47.7, 47.8, 47.9, 48.0, 48.1, 48.2, 48.3, 48.4, 48.5, 48.6, 48.7, 48.8, 48.9, 49.0, 49.1, 49.2, 49.3, 49.4, 49.5, 49.6, 49.7, 49.8, 49.9, 50.0, 50.1, 50.2, 50.3, 50.4, 50.5, 50.6, 50.7, 50.8, 50.9, 51.0, 51.1, 51.2, 51.3, 51.4, 51.5, 51.6, 51.7, 51.8, 51.9, 52.0, 52.1, 52.2, 52.3, 52.4, 52.5, 52.6, 52.7, 52.8, 52.9, 53.0, 53.1, 53.2, 53.3, 53.4, 53.5, 53.6, 53.7, 53.8, 53.9, 54.0, 54.1, 54.2, 54.3, 54.4, 54.5, 54.6, 54.7, 54.8, 54.9, 55.0, 55.1, 55.2, 55.3, 55.4, 55.5, 55.6, 55.7, 55.8, 55.9, 56.0, 56.1, 56.2, 56.3, 56.4, 56.5, 56.6, 56.7, 56.8, 56.9, 57.0, 57.1, 57.2, 57.3, 57.4, 57.5, 57.6, 57.7, 57.8, 57.9, 58.0, 58.1, 58.2, 58.3, 58.4, 58.5, 58.6, 58.7, 58.8, 58.9, 59.0, 59.1, 59.2, 59.3, 59.4, 59.5, 59.6, 59.7, 59.8, 59.9, 60.0};
  const int sizeListTopMass = sizeof(listTopMass)/sizeof(listTopMass[0]);
  const int sizeListTopWidth = sizeof(listTopWidth)/sizeof(listTopWidth[0]);
  
  double gammaProp, numTopPropagator;
  double denomTopPropagator_gen, topPropagator_gen;
  double denomTopPropagator_reco, topPropagator_reco;
  double likelihood_gen[sizeListTopMass][sizeListTopWidth] = {{0}};
  double likelihood_reco[sizeListTopMass][sizeListTopWidth] = {{0}};
  
  
  
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
  
  for (unsigned int d = 0; d < datasets.size(); d++)
  { 
    cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
    nofSelectedEvents = 0;
    nofEventsWith1BJet = 0;
    nofEventsWith2BJets = 0;
    
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << "/ title : " << datasets[d]->Title() << endl;
    if (verbose > 1)
      cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    //open files and load
    cout << "LoadEvent" << endl;
    treeLoader.LoadDataset(datasets[d], anaEnv);
    cout << "LoadEvent" << endl;
    
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    //for (unsigned int ievt = 0; ievt < 1000; ievt++)
    {
      
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      //vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
      
      if (ievt%1000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
      
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
      
      //if (! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
      //  genjets = treeLoader.LoadGenJet(ievt,false);
      //  //sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      //}
      
      
      // BE CAREFUL: TRootGenEvent is now obsolete!
      
      
      
      /////////////////////////////////////
      ///  DETERMINE EVENT SCALEFACTOR  ///
      /////////////////////////////////////
      
      // scale factor for the event
      float scaleFactor = 1.;
      
      
      // PU reweighting
      
      // old method
      //cout << "scalefactor " << scaleFactor << endl;
      double lumiWeight = 1; //LumiWeights.ITweight( (int)event->nTruePU() ); // currently no pile-up reweighting applied
      
      if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
        lumiWeight=1;
      
      // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
      // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );
      
      scaleFactor = scaleFactor*lumiWeight;
      
      
      
      ////////////////////////////
      ///  Include trigger set up here when using data
      ////////////////////////////
      
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if (previousFilename != currentFilename)
      {
        previousFilename = currentFilename;
        iFile++;
        cout << "File changed!!! => iFile = " << iFile << endl;
      }
      
      int currentRun = event->runId();
      
      if (previousRun != currentRun)
        previousRun = currentRun;
      
      
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      //Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
      
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2.);
      
      vector<TRootPFJet*> selectedJets = selection.GetSelectedJets(20, 2.5, true, "Tight");  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(muonPTSel, muonEtaSel, muonRelIsoSel, "Tight", "Spring15");  // GetSelectedMuons(float PtThr, float etaThr, float relIso, string WorkingPoint, string ProductionCampaign)
      vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(26, electronEtaSel, "Tight", "Spring15_25ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
      vector<TRootMuon*> vetoMuons = selection.GetSelectedMuons(muonPTVeto, muonEtaVeto, muonRelIsoVeto, "Loose", "Spring15");  // GetSelectedMuons(float PtThr, float etaThr, float relIso, string WorkingPoint, string ProductionCampaign)
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedElectrons(15, electronEtaVeto, "Loose", "Spring15_25ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
      
      
      /// Sort objects according to pT
      sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)
        
      
      //if (selectedJets.size() >= 4)
      //  if (selectedJets[3]->Pt() < 30) selectedJets.clear();
      
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
      
      selecTableSemiMu.Fill(d,0,scaleFactor);
      /// At the moment do not use trigger
      selecTableSemiMu.Fill(d,1,scaleFactor);
      if (isGoodPV)
      {
        selecTableSemiMu.Fill(d,2,scaleFactor);
        if (selectedMuons.size() == 1)
        {
          selecTableSemiMu.Fill(d,3,scaleFactor);
          if (vetoMuons.size() == 1) {
            selecTableSemiMu.Fill(d,4,scaleFactor);
            if (vetoElectronsSemiMu.size() == 0) {
              selecTableSemiMu.Fill(d,5,scaleFactor);
              if ( selectedJets.size() >= 4 )
              {
                selecTableSemiMu.Fill(d,6,scaleFactor);
                eventSelected = true;
                
                for (unsigned int i = 0; i < selectedJets.size(); i++)
                {
                  if (selectedJets[i]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVMv2Medium) nb_bTaggedJets++;
                }
                		
                if ( nb_bTaggedJets >= 1 )
                {
                  selecTableSemiMu.Fill(d,7,scaleFactor);
                  has1bjet = true;
                  
                  if ( nb_bTaggedJets >= 2 )
                  {
                    selecTableSemiMu.Fill(d,8,scaleFactor);
                    has2bjets = true;
                  }  // at least 2 b-tagged jets
                }  // at least 1 b-tagged jets

              }  // at least 4 jets
            }  // no loose electrons
          }  // no additional loose muons (tight muon is also loose muon)
        }  // 1 good muon
      }  // good PV
      
      
      
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
        vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
        TLorentzVector topQuark, antiTopQuark;
        
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
					
          if ( mcParticles[i]->status() == 1 && mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )		// mu-, W-, tbar
          {
            muMinusFromTop = true;
            genmuon = i;
          }
          if ( mcParticles[i]->status() == 1 && mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )		// mu+, W+, t
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
          cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<" .  This should be equal to 1 !!!"<<endl;
        
        
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
        if (genmuon != -9999 && ROOT::Math::VectorUtil::DeltaR( (TLorentzVector)*mcParticles[genmuon], (TLorentzVector)*selectedMuons[0]) < 0.3)
          muonmatched = true;
        
        
      }  /// End matching
      
      
      
      //////////////////////////////////
      ///  TOP PROPAGATOR (MATCHED)  ///
      //////////////////////////////////
      
      if ( dataSetName.find("TT") == 0 && all4PartonsMatched )
      {
        /// MCPermutation = JetPartonPair[i].first  = jet number
        ///                 JetPartonPair[i].second = parton number
        /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
        
        
        float WMassReco = (*selectedJets[MCPermutation[0]] + *selectedJets[MCPermutation[1]]).M();
        float topMassReco = (*selectedJets[MCPermutation[0]] + *selectedJets[MCPermutation[1]] + *selectedJets[MCPermutation[2]]).M();
        float topMassGen = (*mcParticlesMatching_[hadronicWJet1_.second] + *mcParticlesMatching_[hadronicWJet2_.second] + *mcParticlesMatching_[hadronicBJet_.second]).M();
        
        for (unsigned int jMass = 0; jMass < sizeListTopMass; jMass++)
        {
          for (unsigned int jWidth = 0; jWidth < sizeListTopWidth; jWidth++)
          {
            gammaProp = sqrt( pow( listTopMass[jMass], 4 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 ) );
            numTopPropagator = ( 2 * sqrt(2) * listTopMass[jMass] * listTopWidth[jWidth] * gammaProp ) / ( TMath::Pi() * sqrt( pow(listTopMass[jMass], 2) + gammaProp ) );
            
            /// Generated mass
            denomTopPropagator_gen = pow( pow(topMassGen, 2) - pow(listTopMass[jMass], 2), 2 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 );
            
            topPropagator_gen = numTopPropagator/denomTopPropagator_gen;
            
            likelihood_gen[jMass][jWidth] += -TMath::Log10(topPropagator_gen);
            
            /// Reconstructed mass
            denomTopPropagator_reco = pow( pow(topMassReco, 2) - pow(listTopMass[jMass], 2), 2 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 );
            
            topPropagator_reco = numTopPropagator/denomTopPropagator_reco;
            
            likelihood_reco[jMass][jWidth] += -TMath::Log10(topPropagator_reco);
            
          }  /// End loop jWidth
        }  /// End loop jMass
        
        
        /// Fill plots
        histo1D["W_Mass_Reco_matched"]->Fill(WMassReco);
        histo1D["top_Mass_Reco_matched"]->Fill(topMassReco);
        histo1D["top_Mass_Gen_matched"]->Fill(topMassGen);
        if ( all4JetsMatched_MCdef_ )
        {
          histo1D["W_Mass_Reco_first4matched"]->Fill(WMassReco);
          histo1D["top_Mass_Reco_first4matched"]->Fill(topMassReco);
          histo1D["top_Mass_Gen_first4matched"]->Fill(topMassGen);
        }
      }
      
      
      ////////////////////////////
      ///  Find b-tagged jets  ///
      ////////////////////////////
      
      int label_bJet1 = -9999;
      int label_bJet2 = -9999;
      float pT_bJet1 = -9999.;
      float pT_bJet2 = -9999.;
      for (unsigned int i = 0; i < selectedJets.size(); i++) {
        if ( ! has1bjet ) break;
        if (selectedJets[i]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVMv2Medium) {
          if ( ! has2bjets ) {
            label_bJet1 = i;
            pT_bJet1 = selectedJets[label_bJet1]->Pt();
            break;
          }
          else {
            if (selectedJets[i]->Pt() > pT_bJet1) {
              // Save previous as second best
              if(label_bJet1 >= 0){
                label_bJet2 = label_bJet1;
                pT_bJet2 = pT_bJet1;
              }
              // Keep new one
              label_bJet1 = i;
              pT_bJet1 = selectedJets[label_bJet1]->Pt();
            }
            else if (selectedJets[i]->Pt() > pT_bJet2) {
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
//      BTagCosThetaCalculation* bTagCosThetaCalculation = new BTagCosThetaCalculation();
//      float CosThetaCalculation;
      if (label_bJet1 != -9999 && label_bJet2 != -9999) {
        nofEventsWith2BJets++;
        float recoWMass, recoTopMass_bJet1, recoTopMass_bJet2, WTerm, topTerm_bJet1, topTerm_bJet2, chi2_bJet1, chi2_bJet2;
        float smallestChi2 = 999999.;
        for (int ijet=0; ijet<4; ijet++) {
          for (int jjet=ijet+1; jjet<4; jjet++) {
            if (ijet != label_bJet1 && ijet != label_bJet2 && jjet != label_bJet1 && jjet != label_bJet2) {
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
//               }
//               if (chi2_bJet2 < smallestChi2) {
//                 smallestChi2 = chi2_bJet2;
//                 labelsReco[0] = label_bJet1;
//                 labelsReco[1] = label_bJet2;
//                 labelsReco[2] = ijet;
//                 labelsReco[3] = jjet;
//               }
            }
          }
        }
        
        if ( dataSetName.find("TT") == 0 )
          {
            histo1D["W_Mass_Reco_notMatched"]->Fill(recoWMass);
            histo1D["top_Mass_Reco_notMatched"]->Fill(recoTopMass_bJet1, 0.5);
            histo1D["top_Mass_Reco_notMatched"]->Fill(recoTopMass_bJet2, 0.5);
          }
        
//         //Cos Theta*
//         TLorentzVector lepton = *selectedMuons[0];
//         vector<TLorentzVector> jets;
//         jets.clear();
//         for(int i=0; i< selectedJets.size(); i++) jets.push_back(*selectedJets[i]);
//         
//         TLorentzVector Neutrino(999,999,999,999);
//         
//         CosThetaCalculation = bTagCosThetaCalculation->CalcOrigKins(labelsReco[0],labelsReco[1],lepton,jets,80.4,172.5);
//         if(CosThetaCalculation != 999)
//         {
//           MSPlot["Cos_Theta*"]->Fill(CosThetaCalculation, datasets[d], true, Luminosity*scaleFactor);
//           if (CosThetaCalculation > 1)
//             cout << "Eventnr: " << ievt << "; Cos Theta: " << CosThetaCalculation << endl;
//           
//           Neutrino =  bTagCosThetaCalculation->GetNeutrino();
//         }
        
        
        
					
        
// 	  		//Fill histos
// 	  		if (labelsReco[0] != -9999 && labelsReco[1] != -9999 && labelsReco[2] != -9999 && labelsReco[3] != -9999) {
// 	    		labelsChanged++;
// 	    		if (useMassesAndResolutions && eventselectedSemiMu) selecTableSemiMu.Fill(d,11,scaleFactor);
// 	    		//if (useMassesAndResolutions && eventselectedSemiEl) selecTableSemiEl.Fill(d,12,scaleFactor);
// 
// 	    		MSPlot["Chi2_2btags"+Flav]->Fill(smallestChi2, datasets[d], true, Luminosity*scaleFactor);
// 
// 					float Wmass_2btags = (*selectedJets[labelsReco[2]] + *selectedJets[labelsReco[3]]).M();
// 					MSPlot["W_Mass_2btags"+Flav]->Fill(Wmass_2btags, datasets[d], true, Luminosity*scaleFactor);
// 	    		float HtTop_2btags = selectedJets[labelsReco[1]]->Pt() + selectedJets[labelsReco[2]]->Pt() + selectedJets[labelsReco[3]]->Pt();
// 	    		MSPlot["hadTop_Ht_2btags"+Flav]->Fill(HtTop_2btags, datasets[d], true, Luminosity*scaleFactor);
// 	    		float hadtopmass_2btags = ( *selectedJets[labelsReco[1]] + *selectedJets[labelsReco[2]] +  *selectedJets[labelsReco[3]]).M();
// 	    		float hadtoppt_2btags = ( *selectedJets[labelsReco[1]] + *selectedJets[labelsReco[2]] +  *selectedJets[labelsReco[3]]).Pt();
// 	    		MSPlot["hadTop_Mass_2btags"+Flav]->Fill(hadtopmass_2btags, datasets[d], true, Luminosity*scaleFactor);
// 	    		MSPlot["hadTop_Pt_2btags"+Flav]->Fill(hadtoppt_2btags, datasets[d], true, Luminosity*scaleFactor);
// 
// 
// 	    		if(CosThetaCalculation != 999) {
// 	      		float leptonicTopMass_2btags = (lepton + Neutrino + jets[labelsReco[0]]).M();
// 	      		float TTbarMass_2btags = leptonicTopMass_2btags + hadtopmass_2btags;
// 	      		MSPlot["TTbar_Mass_2btags"+Flav]->Fill(TTbarMass_2btags, datasets[d], true, Luminosity*scaleFactor);
// 
// 	      		float leptonicTopAngle_2btags = (lepton + Neutrino + jets[labelsReco[0]]).Phi();
// 	      		float hadronicTopAngle_2btags = (jets[labelsReco[1]] + jets[labelsReco[2]] + jets[labelsReco[3]]).Phi();
// 	      		float TTbarAngle_2btags = leptonicTopAngle_2btags - hadronicTopAngle_2btags;
// 	      		if (TTbarAngle_2btags < 0) {
// 							TTbarAngle_2btags = - TTbarAngle_2btags;
// 	      		}
// 	      		MSPlot["TTbar_Angle_2btags"+Flav]->Fill(TTbarAngle_2btags, datasets[d], true, Luminosity*scaleFactor);
// 
// 	    		}
// 	  		}
// 	  		else {
// 	    		//When no Chi2 combination is found:
// 	    		cout << "Eventnr. " << ievt << ": no Chi2 found." << endl;
// 	  		}
      
      
      }
      
      
      
      ////////////////////
      ///  FILL PLOTS  ///
      ////////////////////
      
      double HT = selectedJets[0]->Pt()+selectedJets[1]->Pt()+selectedJets[2]->Pt()+selectedJets[3]->Pt();
      
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
      MSPlot["leadingJet_pT"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["jet2_pT"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["jet3_pT"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["jet4_pT"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Ht_4leadingJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
      
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
      }
      
      
      
      
      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  /// Loop on events
    
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofEventsWith1BJet << " events with 1 b tagged jet." << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofEventsWith2BJets << " events with 2 b tagged jets." << endl;
    cout << "Number of matched events: " << nofMatchedEvents << endl;
    
    
    /// Fill histogram log likelihood
    if ( dataSetName.find("TT") == 0 )
    {
      for (unsigned int jMass = 0; jMass < sizeListTopMass; jMass++)
        {
          for (unsigned int jWidth = 0; jWidth < sizeListTopWidth; jWidth++)
          {
            histo2D["LogLikeWidthMass_Reco"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_reco[jMass][jWidth]);
            histo2D["LogLikeWidthMass_Reco_zoom"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_reco[jMass][jWidth]);
            histo2D["LogLikeWidthMass_Gen"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_gen[jMass][jWidth]);
          }
        }
    }
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  /// Loop on datasets
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  string pathPNG = "Plots/";
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
  string selectiontableMu = "SelectionTable_test_SemiMu.tex";
  //selecTableSemiMu.Write(selectiontableMu.c_str());
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
