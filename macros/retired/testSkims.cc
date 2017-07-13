////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////



#include "TStyle.h"
#include <ctime>
#include <cmath>
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
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/BTagCalibrationStandalone.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"

// user defined
#include "Tools/interface/ResolutionFunctions.h"
#include "Tools/interface/Trigger.h"


using namespace std;
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
  
  string rootFileName = "testSkims_output_1file_"+dateString+".root";
  string selectiontableMu = "SelectionTable_testSkims_SemiMu_"+dateString+".tex";
  string selectionTableSkimTest = "SelectionTable_testSkimStepByStep_1file_"+dateString+".tex";
  string pathPNG = "Plots_Skims_1file_"+dateString+"/";
  int iReducedDataSets = 1;
  
  
  /////////////////////
  ///  Configuration
  /////////////////////
  
  bool eventSelected = false;
  bool has1bjet = false;
  bool has2bjets = false;
  int nofSelectedEvents = 0;
  int nb_bTaggedJets = 0;
  int nofEventsWith1BJet = 0;
  int nofEventsWith2BJets = 0;
  
  
  /// xml file
  string xmlFileName ="config/testSkims.xml";
  
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
  
  float Luminosity = 1;
  //float LuminosityMu = oldLuminosity;
  //float LuminosityEl = oldLuminosity;
  
  //bool foundMu = false;
  //bool foundEl = false;
  
  
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
  
  histo1D["muon_pT_00"] = new TH1F("muon_pT_00","Transverse momentum of the muon; p_{T} [GeV]", 200, 0, 200);
  histo1D["muon_eta_00"] = new TH1F("muon_eta_00","Pseudorapidity of the muon; #eta", 60, -3, 3);
  histo1D["leadingJet_pT_00"] = new TH1F("leadingJet_pT_00","Transverse momentum of the leading jet; p_{T} [GeV]", 800, 0, 800);
  histo1D["Ht_4leadingJets_00"] = new TH1F("Ht_4leadingJets_00","Scalar sum of transverse momenta of the 4 leading jets; H_{T} [GeV]", 120, 0, 1200);
  histo1D["nJets_00"] = new TH1F("nJets_00","Number of jets", 13, -0.5, 12.5);
  histo1D["nBJets_00"] = new TH1F("nBJets_00","Number of b-tagged jets", 13, -0.5, 12.5);
  
  histo1D["muon_pT_01"] = new TH1F("muon_pT_01","Transverse momentum of the muon; p_{T} [GeV]", 200, 0, 200);
  histo1D["muon_eta_01"] = new TH1F("muon_eta_01","Pseudorapidity of the muon; #eta", 60, -3, 3);
  histo1D["leadingJet_pT_01"] = new TH1F("leadingJet_pT_01","Transverse momentum of the leading jet; p_{T} [GeV]", 800, 0, 800);
  histo1D["Ht_4leadingJets_01"] = new TH1F("Ht_4leadingJets_01","Scalar sum of transverse momenta of the 4 leading jets; H_{T} [GeV]", 120, 0, 1200);
  histo1D["nJets_01"] = new TH1F("nJets_01","Number of jets", 13, -0.5, 12.5);
  histo1D["nBJets_01"] = new TH1F("nBJets_01","Number of b-tagged jets", 13, -0.5, 12.5);
  
//   histo1D["muon_pT_02"] = new TH1F("muon_pT_02","Transverse momentum of the muon; p_{T} [GeV]", 200, 0, 200);
//   histo1D["muon_eta_02"] = new TH1F("muon_eta_02","Pseudorapidity of the muon; #eta", 60, -3, 3);
//   histo1D["leadingJet_pT_02"] = new TH1F("leadingJet_pT_02","Transverse momentum of the leading jet; p_{T} [GeV]", 800, 0, 800);
//   histo1D["Ht_4leadingJets_02"] = new TH1F("Ht_4leadingJets_02","Scalar sum of transverse momenta of the 4 leading jets; H_{T} [GeV]", 120, 0, 1200);
//   histo1D["nJets_02"] = new TH1F("nJets_02","Number of jets", 13, -0.5, 12.5);
//   histo1D["nBJets_02"] = new TH1F("nBJets_02","Number of b-tagged jets", 13, -0.5, 12.5);
//   
//   histo1D["muon_pT_03"] = new TH1F("muon_pT_03","Transverse momentum of the muon; p_{T} [GeV]", 200, 0, 200);
//   histo1D["muon_eta_03"] = new TH1F("muon_eta_03","Pseudorapidity of the muon; #eta", 60, -3, 3);
//   histo1D["leadingJet_pT_03"] = new TH1F("leadingJet_pT_03","Transverse momentum of the leading jet; p_{T} [GeV]", 800, 0, 800);
//   histo1D["Ht_4leadingJets_03"] = new TH1F("Ht_4leadingJets_03","Scalar sum of transverse momenta of the 4 leading jets; H_{T} [GeV]", 120, 0, 1200);
//   histo1D["nJets_03"] = new TH1F("nJets_03","Number of jets", 13, -0.5, 12.5);
//   histo1D["nBJets_03"] = new TH1F("nBJets_03","Number of b-tagged jets", 13, -0.5, 12.5);
  
  
  
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
  
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  //selecTableSemiMu.SetLuminosity(LuminosityMu);
  selecTableSemiMu.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
  vector<string> CutsSelecTableSkimTest;
  CutsSelecTableSkimTest.push_back(string("preselected"));
  CutsSelecTableSkimTest.push_back(string("1 muon"));
  CutsSelecTableSkimTest.push_back(string("4 jets"));
  CutsSelecTableSkimTest.push_back(string("muon pT > 18 GeV"));
  CutsSelecTableSkimTest.push_back(string("jet 1 pT > 18 GeV"));
  CutsSelecTableSkimTest.push_back(string("jet 2 pT > 18 GeV"));
  CutsSelecTableSkimTest.push_back(string("jet 3 pT > 18 GeV"));
  CutsSelecTableSkimTest.push_back(string("jet 4 pT > 18 GeV"));
  CutsSelecTableSkimTest.push_back(string("muon eta < 2.1"));
  CutsSelecTableSkimTest.push_back(string("jet eta < 2.4"));
  SelectionTable selecTableSkimTest(CutsSelecTableSkimTest, datasets);
  selecTableSkimTest.SetLuminosity(Luminosity);
  
  ///////////////////////////////
  ///  Single Muon Selection  ///
  ///////////////////////////////
  
  /// Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO
  
  float muonPTSel = 26; // GeV
  float muonEtaSel = 2.1;
  float muonRelIsoSel = 0.15;  // Tight muon
  string muonWP = "Tight";
  
  float muonPTVeto = 10; // GeV
  float muonEtaVeto = 2.5;
  float muonRelIsoVeto = 0.25;  // Loose muon
  
  
  
  ///////////////////////////////////
  ///  Single Electron Selection  ///
  ///////////////////////////////////
  
  // To do
  float electronPTSel = 24; // GeV
  float electronEtaSel = 2.5;
  string electronWP = "Tight";
  
  float electronPTVeto = 15; // GeV
  float electronEtaVeto = 2.5;
  
  
  
  ///////////////////////
  ///  Jet Selection  ///
  ///////////////////////
  
  /// Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopJME
  
  float jetPT = 20; // GeV
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
    bool isData = false;
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
    
    cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
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
      vector < TRootJet* > init_fatjets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
      
//      if (ievt%1000 == 0)
//        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
      if (((int)nEvents[d])%10000 == 0)
        std::cout << "Processing the " << ((int)nEvents[d]) << "th event (" << (nEvents[d]*((double)iReducedDataSets)/((double)datasets[d]->NofEvtsToRunOver()))*100  << "%)" << flush << "\r";
      
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
      
      if (! isData ) {
        genjets = treeLoader.LoadGenJet(ievt,false);
        //sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      
      // BE CAREFUL: TRootGenEvent is now obsolete!
      
      
      
      /////////////////////////////////////
      ///  DETERMINE EVENT SCALEFACTOR  ///
      /////////////////////////////////////
      
      // scale factor for the event
      float scaleFactor = 1.;
      
      
      
      ////////////////////////////
      ///  Include trigger set up here when using data
      ////////////////////////////
      
      bool trigged = false;
      
      /// Fill selection table before trigger
      selecTableSemiMu.Fill(d,0,scaleFactor);
      
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if (previousFilename != currentFilename)
      {
        previousFilename = currentFilename;
        iFile++;
        cout << "File changed!!! => iFile = " << iFile << endl;
      }
      
      
      /// Fill selection table after trigger (hypothetically..)
      selecTableSemiMu.Fill(d,1,scaleFactor);
      
      
      
      /////////////////
      ///  Skim tests
      /////////////////
      
      selecTableSkimTest.Fill(d,0,scaleFactor);
      if ( init_muons.size() > 0 )
      {
        selecTableSkimTest.Fill(d,1,scaleFactor);
        if ( init_jets_corrected.size() >= 4 )
        {
          selecTableSkimTest.Fill(d,2,scaleFactor);
          if ( init_muons[0]->Pt() >= 18. )
          {
            selecTableSkimTest.Fill(d,3,scaleFactor);
            if ( init_jets_corrected[0]->Pt() >= 18. )
            {
              selecTableSkimTest.Fill(d,4,scaleFactor);
              if ( init_jets_corrected[1]->Pt() >= 18. )
              {
                selecTableSkimTest.Fill(d,5,scaleFactor);
                if ( init_jets_corrected[2]->Pt() >= 18. )
                {
                  selecTableSkimTest.Fill(d,6,scaleFactor);
                  if ( init_jets_corrected[3]->Pt() >= 18. )
                  {
                    selecTableSkimTest.Fill(d,7,scaleFactor);
                    if ( fabs(init_muons[0]->Eta()) < 2.1 )
                    {
                      selecTableSkimTest.Fill(d,8,scaleFactor);
                      if ( fabs(init_jets_corrected[0]->Eta()) < 2.4 && fabs(init_jets_corrected[1]->Eta()) < 2.4 && fabs(init_jets_corrected[2]->Eta()) < 2.4 && fabs(init_jets_corrected[3]->Eta()) < 2.4 )
                      {
                        selecTableSkimTest.Fill(d,9,scaleFactor);
                      }
                    }
                  }  /// jet 4 pT > 18 GeV
                }  /// jet 3 pT > 18 GeV
              }  /// jet 2 pT > 18 GeV
            }  /// jet 1 pT > 18 GeV
          }  /// muon pT > 18 GeV
        }  /// 4 jets
      }  /// 1 muon
      
      
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      //Declare selection instance
      Run2Selection selection(init_jets_corrected, init_fatjets, init_muons, init_electrons, mets, event->fixedGridRhoFastjetAll());
      
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2.);
      
      vector<TRootPFJet*> selectedJets = selection.GetSelectedJets(jetPT, jetEta, true, "Tight");  // PtThr, EtaThr, applyJetID, TightLoose
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(muonPTSel, muonEtaSel, muonRelIsoSel, muonWP, "Spring15");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(electronPTSel, electronEtaSel, electronWP, "Spring15_25ns", true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased
      vector<TRootMuon*> vetoMuons = selection.GetSelectedMuons(muonPTVeto, muonEtaVeto, muonRelIsoVeto, "Loose", "Spring15");  // PtThr, etaThr, relIso, WorkingPoint, ProductionCampaign
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedElectrons(electronPTVeto, electronEtaVeto, "Veto", "Spring15_25ns", true);  // PtThr, etaThr, WorkingPoint, ProductionCampaign, CutsBased
      
      
//      /// Sort objects according to pT
//      sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)
        
      
      //if (selectedJets.size() >= 4)
      //  if (selectedJets[3]->Pt() < 30) selectedJets.clear();
      
      
      vector<TRootPFJet*> selectedBJets;
      for (int i = 0; i < selectedJets.size(); i++)
      {
        if ( selectedJets[i]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium )
          selectedBJets.push_back(selectedJets[i]);
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
        if (selectedMuons.size() == 1)
        {
          selecTableSemiMu.Fill(d,3,scaleFactor);
          if (vetoMuons.size() == 1) {
            selecTableSemiMu.Fill(d,4,scaleFactor);
            if (vetoElectronsSemiMu.size() == 0) {
              selecTableSemiMu.Fill(d,5,scaleFactor);
              
              /// First 4 jets need pT > 30 GeV
              if (selectedJets.size() >= 4)
              {
                if (selectedJets[3]->Pt() < 30) selectedJets.clear();
              }
              
              if ( selectedJets.size() >= 4 )
              {
                selecTableSemiMu.Fill(d,6,scaleFactor);
                eventSelected = true;
                
                for (unsigned int i = 0; i < selectedJets.size(); i++)
                {
                  if (selectedJets[i]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium) nb_bTaggedJets++;
                }
                		
                if ( nb_bTaggedJets >= 1 )
                {
                  selecTableSemiMu.Fill(d,7,scaleFactor);
                  has1bjet = true;
                  
                  if ( nb_bTaggedJets >= 2 )
                  {
                    selecTableSemiMu.Fill(d,8,scaleFactor);
                    has2bjets = true;
                    nofEventsWith2BJets++;
                  }  // at least 2 b-tagged jets
                  else { nofEventsWith1BJet++;}
                }  // at least 1 b-tagged jets

              }  // at least 4 jets
            }  // no veto electrons
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
      
      
      
      ////////////////////////////
      ///  Find b-tagged jets  ///
      ////////////////////////////
      
      int label_bJet1 = -9999;
      int label_bJet2 = -9999;
      float pT_bJet1 = -9999.;
      float pT_bJet2 = -9999.;
      for (unsigned int i = 0; i < selectedBJets.size(); i++)
      {
        if ( label_bJet1 != -9999 )
        {
          if ( selectedJets[i]->Pt() > pT_bJet1 )
          {
            // Save previous as second best
            label_bJet2 = label_bJet1;
            pT_bJet2 = pT_bJet1;
            
            // Keep new one
            label_bJet1 = i;
            pT_bJet1 = selectedJets[label_bJet1]->Pt();
          }
          else if ( selectedJets[i]->Pt() > pT_bJet2 )
          {
            label_bJet2 = i;
            pT_bJet2 = selectedJets[label_bJet2]->Pt();
          }
        }
        else
        {
          label_bJet1 = i;
          pT_bJet1 = selectedJets[label_bJet1]->Pt();
        }
      }
      
      
      
      ////////////////////
      ///  FILL PLOTS  ///
      ////////////////////
      
      double HT = selectedJets[0]->Pt()+selectedJets[1]->Pt()+selectedJets[2]->Pt()+selectedJets[3]->Pt();
      
      histo1D[("muon_pT_"+ConvertIntToString(d,true)).c_str()]->Fill(selectedMuons[0]->Pt());
      histo1D[("muon_eta_"+ConvertIntToString(d,true)).c_str()]->Fill(selectedMuons[0]->Eta());
      histo1D[("leadingJet_pT_"+ConvertIntToString(d,true)).c_str()]->Fill(selectedJets[0]->Pt());
      histo1D[("Ht_4leadingJets_"+ConvertIntToString(d,true)).c_str()]->Fill(HT);
      histo1D[("nJets_"+ConvertIntToString(d,true)).c_str()]->Fill(selectedJets.size());
      histo1D[("nBJets_"+ConvertIntToString(d,true)).c_str()]->Fill(selectedBJets.size());
      
      
      
      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  /// Loop on events
    
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofEventsWith1BJet << " events with 1 b tagged jet." << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofEventsWith2BJets << " events with 2 b tagged jets." << endl;
    
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  /// Loop on datasets
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  //string pathPNG = "PlotsOneFourth/";
  //mkdir(pathPNG.c_str(),0777);
  
  ///Write histograms
  fout->cd();  
  // 1D
  //TDirectory* th1dir = fout->mkdir("1D_histograms");
  //th1dir->cd();
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

//   // 2D
//   TDirectory* th2dir = fout->mkdir("2D_histograms");
//   th2dir->cd();
//   for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
//   {
//     TH2F *temp = it->second;
//     temp->Write();
//     TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
//     tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
//   }
  
  
  
  ///Selection tables
  selecTableSemiMu.TableCalculator(false, true, true, true, true);
  selecTableSemiMu.Write(selectiontableMu.c_str(), true, true, true, true, true, true, false);
  
  selecTableSkimTest.TableCalculator(false, true, true, true, true);
  selecTableSkimTest.Write(selectionTableSkimTest.c_str(), true, true, true, true, true, true, false);
  
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
