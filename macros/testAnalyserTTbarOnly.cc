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
#include "../TopTreeAnalysisBase/Tools/interface/Trigger.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFunctions.h"


using namespace std;
using namespace reweight;
using namespace TopTree;


string ConvertIntToString(int Number, bool pad)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
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
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  bool useOneHalfOfDataSets = true;
  bool useOneFourthOfDataSets = false;
  bool useOneTenthOfDataSets = false;
  bool useTestSample = false;
  
  
  string rootFileName = "testAnalyser_output_TTbarOnly_"+dateString+".root";
  string selectiontableMu = "SelectionTable_testTTbarOnly_"+dateString+".tex";
  string pathPNG = "Plots_TTbarOnly_"+dateString+"/";
  int iReducedDataSets = 1;
  
  if (useOneHalfOfDataSets)
  {
    rootFileName = "testAnalyser_output_oneHalfTTbarOnly_"+dateString+".root";
    selectiontableMu = "SelectionTable_testOneHalfTTbarOnly_"+dateString+".tex";
    pathPNG = "PlotsOneHalf_TTbarOnly_"+dateString+"/";
    iReducedDataSets = 2;
  }
  if (useOneFourthOfDataSets)
  {
    rootFileName = "testAnalyser_output_TTbarOnly_"+dateString+".root";
    selectiontableMu = "SelectionTable_testOneFourthTTbarOnly_"+dateString+".tex";
    pathPNG = "PlotsOneFourth_TTbarOnly_"+dateString+"/";
    iReducedDataSets = 4;
  }
  if (useOneTenthOfDataSets)
  {
    rootFileName = "testAnalyser_output_oneTenthTTbarOnly_"+dateString+".root";
    selectiontableMu = "SelectionTable_testOneTenthTTbarOnly_"+dateString+".tex";
    pathPNG = "PlotsOneTenth_TTbarOnly_"+dateString+"/";
    iReducedDataSets = 10;
  }
  if (useTestSample)
  {
    rootFileName = "testAnalyser_output_testSampleTTbarOnly_"+dateString+".root";
    selectiontableMu = "SelectionTable_testSampleTTbarOnly_"+dateString+".tex";
    pathPNG = "PlotsTestSample_TTbarOnly_"+dateString+"/";
    iReducedDataSets = 200;
  }
  
  
  
  /////////////////////
  ///  Configuration
  /////////////////////
  
  bool testTTbarOnly = true;  
  bool calculateResolutionFunctions = true;
  bool doChi2 = false;
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
//  int nofCorrectlyMatched = 0;
//  int nofNotCorrectlyMatched = 0;
  int nofEventsWith1BJet = 0;
  int nofEventsWith2BJets = 0;
  int nofNegWeights = 0;
  int nofPosWeights = 0;
  int nofEventsHLTv2 = 0;
  int nofEventsHLTv3 = 0;
  int nofEventsJetLeptonCleaned = 0;
  
  
  /// xml file
  string xmlFileName ="config/testTTbarOnly.xml";
  
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
    if ( dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") == 0 )
    { 
      datasets[d]->SetColor(kBlue-2);
      if ( dataSetName.find("tW") == 0 )
      { 
        datasets[d]->SetTitle("ST tW");
        datasets[d]->SetColor(kBlue-4);
      }
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
  
  
  histo1D["WMass_reco_matched"] = new TH1F("WMass_reco_matched","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 125, 0, 250);
  histo1D["topMass_reco_matched"] = new TH1F("topMass_reco_matched","Reconstructed top mass of matched events; M_{t} [GeV]", 300, 50, 400);
  histo1D["topMass_gen_matched"] = new TH1F("topMass_gen_matched","Generated top mass of matched events; M_{t} [GeV]", 300, 100, 250);
  histo1D["WMass_reco_first4matched"] = new TH1F("WMass_reco_first4matched","Reconstructed hadronic W mass of events where 4 hardest jets are matched; M_{W} [GeV]", 125, 0, 250);
  histo1D["topMass_reco_first4matched"] = new TH1F("topMass_reco_first4matched","Reconstructed top mass of events where 4 hardest jets are matched; M_{t} [GeV]", 350, 50, 400);
  histo1D["topMass_gen_first4matched"] = new TH1F("topMass_gen_first4matched","Generated top mass of events where partons are matched to 4 hardest jets; M_{t} [GeV]", 300, 50, 400);
  
  
  // Chi2 free
  // 2b
  histo1D["2b_Chi2Free_W_mass_reco_notMatched"] = new TH1F("2b_Chi2Free_W_mass_reco_notMatched","Reconstructed hadronic W mass of unmatched events with 2 b-tagged jets (free Chi2); M_{W} [GeV]", 80, 0, 800);
  histo1D["2b_Chi2Free_top_mass_reco_notMatched"] = new TH1F("2b_Chi2Free_top_mass_reco_notMatched","Reconstructed top mass of unmatched events with 2 b-tagged jets (free Chi2); M_{t} [GeV]", 80, 0, 800);
  histo1D["2b_Chi2Free_top_mass_reco_notMatched_4jets"] = new TH1F("2b_Chi2Free_top_mass_reco_notMatched_4jets","Reconstructed top mass of unmatched events with 2 b-tagged jets out of 4 in total (free Chi2); M_{t} [GeV]", 80, 0, 800);
  
  histo1D["2b_Chi2Free_lepTop_mass_notMatched"] = new TH1F("2b_lepTop_mass_notMatched","Reconstructed leptonic top mass of unmatched events with 2 b-tagged jets; M_{lb} [GeV]", 80, 0, 800);
  histo1D["2b_Chi2Free_lepTop_mass_notMatched_4jets"] = new TH1F("2b_lepTop_mass_notMatched_4jets","Reconstructed leptonic top mass of unmatched events with 2 b-tagged jets out of 4 in total; M_{lb} [GeV]", 80, 0, 800);

  histo1D["2b_Chi2Free_hadTop_mass_notMatched_mlb_cut"] = new TH1F("2b_Chi2Free_hadTop_mass_notMatched_mlb_cut","Reconstructed hadronic top mass of unmatched events with 2 b-tagged jets and a cut on M_{lb} (free Chi2); M_{t} [GeV]", 80, 0, 800);
  histo1D["2b_Chi2Free_hadTop_mass_notMatched_4jets_mlb_cut"] = new TH1F("2b_Chi2Free_hadTop_mass_notMatched_4jets_mlb_cut","Reconstructed hadronic top mass of unmatched events with 2 b-tagged jets out of 4 in total and a cut on M_{lb} (free Chi2); M_{t} [GeV]", 80, 0, 800);
  
  histo2D["2b_2D_unmatched_hadTopMass_lepTopMass"] = new TH2F("2b_2D_unmatched_hadTopMass_lepTopMass", "Mass of the leptonic top mass vs. the mass of the hadronic top mass (unmatched); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);

  histo1D["2b_Chi2Free_ttbar_mass_notMatched"] = new TH1F("2b_Chi2Free_ttbar_mass_notMatched","Reconstructed ttbar mass of unmatched events with 2 b-tagged jets (free Chi2); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["2b_Chi2Free_ttbar_mass_notMatched_4jets"] = new TH1F("2b_Chi2Free_ttbar_mass_notMatched_4jets","Reconstructed ttbar mass of unmatched events with 2 b-tagged jets out of 4 in total (free Chi2); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["2b_Chi2Free_ttbar_mass_notMatched_mlb_cut"] = new TH1F("2b_Chi2Free_ttbar_mass_notMatched_mlb_cut","Reconstructed ttbar mass of unmatched events with 2 b-tagged jets and a cut on M_{lb} (free Chi2); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["2b_Chi2Free_ttbar_mass_notMatched_4jets_mlb_cut"] = new TH1F("2b_Chi2Free_ttbar_mass_notMatched_4jets_mlb_cut","Reconstructed ttbar mass of unmatched events with 2 b-tagged jets out of 4 in total and a cut on M_{lb} (free Chi2); M_{t#bar{t}} [GeV]", 100, 0, 1000);

  histo1D["2b_Chi2Free_dR_lep_b_notMatched"] = new TH1F("2b_Chi2Free_dR_Lep_b_notMatched","Minimal delta R between the lepton and a b jet in events with 2 b-tagged jets (free Chi2); #Delta R(l,b)", 52, 0, 5);
  histo1D["2b_Chi2Free_dR_lep_b_notMatched_4jets"] = new TH1F("2b_Chi2Free_dR_Lep_b_notMatched_4jets","Minimal delta R between the lepton and a b jet in events with 2 b-tagged jets out of 4 in total (free Chi2); #Delta R(l,b)", 25, 0, 5);
  histo1D["2b_Chi2Free_dR_lep_b_notMatched_mlb_cut"] = new TH1F("2b_Chi2Free_dR_Lep_b_notMatched_mlb_cut","Minimal delta R between the lepton and a b jet in events with 2 b-tagged jets and a cut on M_{lb} (free Chi2); #Delta R(l,b)", 25, 0, 5);
  histo1D["2b_Chi2Free_dR_lep_b_notMatched_4jets_mlb_cut"] = new TH1F("2b_Chi2Free_dR_Lep_b_notMatched_4jets_mlb_cut","Minimal delta R between the lepton and a b jet in events with 2 b-tagged jets out of 4 in total and a cut on M_{lb} (free Chi2); #Delta R(l,b)", 25, 0, 5);
  
  // 1b
  histo1D["1b_Chi2Free_W_mass_reco_notMatched"] = new TH1F("1b_Chi2Free_W_mass_reco_notMatched","Reconstructed hadronic W mass of unmatched events with 1 b-tagged jet (free Chi2); M_{W} [GeV]", 80, 0, 800);
  histo1D["1b_Chi2Free_top_mass_reco_notMatched"] = new TH1F("1b_Chi2Free_top_mass_reco_notMatched","Reconstructed top mass of unmatched events with 1 b-tagged jet (free Chi2); M_{t} [GeV]", 80, 0, 800);
  histo1D["1b_Chi2Free_top_mass_reco_notMatched_4jets"] = new TH1F("1b_Chi2Free_top_mass_reco_notMatched_4jets","Reconstructed top mass of unmatched events with 1 b-tagged jet out of 4 in total (free Chi2); M_{t} [GeV]", 80, 0, 800);
  
  histo1D["1b_Chi2Free_lepTop_mass_notMatched"] = new TH1F("1b_lepTop_mass_notMatched","Reconstructed leptonic top mass of unmatched events with 1 b-tagged jet; M_{lb} [GeV]", 80, 0, 800);
  histo1D["1b_Chi2Free_lepTop_mass_notMatched_4jets"] = new TH1F("1b_lepTop_mass_notMatched_4jets","Reconstructed leptonic top mass of unmatched events with 1 b-tagged jet out of 4 in total; M_{lb} [GeV]", 80, 0, 800);

  histo1D["1b_Chi2Free_hadTop_mass_notMatched_mlb_cut"] = new TH1F("1b_Chi2Free_hadTop_mass_notMatched_mlb_cut","Reconstructed hadronic top mass of unmatched events with 1 b-tagged jet and a cut on M_{lb} (free Chi2); M_{t} [GeV]", 80, 0, 800);
  histo1D["1b_Chi2Free_hadTop_mass_notMatched_4jets_mlb_cut"] = new TH1F("1b_Chi2Free_hadTop_mass_notMatched_4jets_mlb_cut","Reconstructed hadronic top mass of unmatched events with 1 b-tagged jet out of 4 in total and a cut on M_{lb} (free Chi2); M_{t} [GeV]", 80, 0, 800);
  
  histo2D["1b_2D_unmatched_hadTopMass_lepTopMass"] = new TH2F("1b_2D_unmatched_hadTopMass_lepTopMass", "Mass of the leptonic top mass vs. the mass of the hadronic top mass (unmatched); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  
  histo1D["1b_Chi2Free_ttbar_mass_notMatched"] = new TH1F("1b_Chi2Free_ttbar_mass_notMatched","Reconstructed ttbar mass of unmatched events with 1 b-tagged jet (free Chi2); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["1b_Chi2Free_ttbar_mass_notMatched_4jets"] = new TH1F("1b_Chi2Free_ttbar_mass_notMatched_4jets","Reconstructed ttbar mass of unmatched events with 1 b-tagged jet out of 4 in total (free Chi2); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["1b_Chi2Free_ttbar_mass_notMatched_mlb_cut"] = new TH1F("1b_Chi2Free_ttbar_mass_notMatched_mlb_cut","Reconstructed ttbar mass of unmatched events with 1 b-tagged jet and a cut on M_{lb} (free Chi2); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["1b_Chi2Free_ttbar_mass_notMatched_4jets_mlb_cut"] = new TH1F("1b_Chi2Free_ttbar_mass_notMatched_4jets_mlb_cut","Reconstructed ttbar mass of unmatched events with 1 b-tagged jet out of 4 in total and a cut on M_{lb} (free Chi2); M_{t#bar{t}} [GeV]", 100, 0, 1000);

  histo1D["1b_Chi2Free_dR_lep_b_notMatched"] = new TH1F("1b_Chi2Free_dR_Lep_b_notMatched","Minimal delta R between the lepton and a b jet in events with 1 b-tagged jet (free Chi2); #Delta R(l,b)", 52, 0, 5);
  histo1D["1b_Chi2Free_dR_lep_b_notMatched_4jets"] = new TH1F("1b_Chi2Free_dR_Lep_b_notMatched_4jets","Minimal delta R between the lepton and a b jet in events with 1 b-tagged jet out of 4 in total (free Chi2); #Delta R(l,b)", 25, 0, 5);
  histo1D["1b_Chi2Free_dR_lep_b_notMatched_mlb_cut"] = new TH1F("1b_Chi2Free_dR_Lep_b_notMatched_mlb_cut","Minimal delta R between the lepton and a b jet in events with 1 b-tagged jet and a cut on M_{lb} (free Chi2); #Delta R(l,b)", 25, 0, 5);
  histo1D["1b_Chi2Free_dR_lep_b_notMatched_4jets_mlb_cut"] = new TH1F("1b_Chi2Free_dR_Lep_b_notMatched_4jets_mlb_cut","Minimal delta R between the lepton and a b jet in events with 1 b-tagged jet out of 4 in total and a cut on M_{lb} (free Chi2); #Delta R(l,b)", 25, 0, 5);
  
  
  /// Leptonic top quark
  histo1D["lepTop_mass_matched_corr"] = new TH1F("lepTop_mass_matched_corr","Reconstructed leptonic top mass using correctly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_wrong"] = new TH1F("lepTop_mass_matched_wrong","Reconstructed leptonic top mass using wrongly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_wrong_hadrB"] = new TH1F("lepTop_mass_matched_wrong_hadrB","Reconstructed leptonic top mass using wrongly matched events (using hadronic b jet); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_wrong_noB"] = new TH1F("lepTop_mass_matched_wrong_noB","Reconstructed leptonic top mass using wrongly matched events (b jet is mistagged); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_corr_4jets"] = new TH1F("lepTop_mass_matched_corr_4jets","Reconstructed leptonic top mass using correctly matched events with exactly 4 jets; M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_wrong_4jets"] = new TH1F("lepTop_mass_matched_wrong_4jets","Reconstructed leptonic top mass using wrongly matched events with exactly 4 jets; M_{lb} [GeV]", 80, 0, 800);
  
  
  /// ttbar mass
  histo1D["ttbar_mass_matched_corr"] = new TH1F("ttbar_mass_matched_corr","Reconstructed mass of the top quark pair using correctly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong"] = new TH1F("ttbar_mass_matched_wrong","Reconstructed mass of the top quark pair using wrongly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong_hadrB"] = new TH1F("ttbar_mass_matched_wrong_hadrB","Reconstructed mass of the top quark pair using wrongly matched events (using hadronic b jet); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong_noB"] = new TH1F("ttbar_mass_matched_wrong_noB","Reconstructed mass of the top quark pair using wrongly matched events (b jet is mistagged); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_corr_4jets"] = new TH1F("ttbar_mass_matched_corr_4jets","Reconstructed mass of the top quark pair using correctly matched events with exactly 4 jets; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong_4jets"] = new TH1F("ttbar_mass_matched_wrong_4jets","Reconstructed mass of the top quark pair using wrongly matched events with exactly 4 jets; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  
  /// delta R (lep, b)
  histo1D["dR_Lep_B"] = new TH1F("dR_Lep_B","Delta R between the lepton and the leptonic b jet; #Delta R(l,b)", 25, 0, 5);
  
  /// 2D plots: lepTopMass, ttbarMass, dR_Lep_B
  histo2D["2D_matched_lepTopMass_corr_ttbarMass_corr"] = new TH2F("2D_matched_lepTopMass_corr_ttbarMass_corr", "ttbarMass vs. leptonic top mass using correctly matched events; M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_lepTopMass_wrong_ttbarMass_wrong_hadrB"] = new TH2F("2D_matched_lepTopMass_wrong_ttbarMass_wrong_hadrB", "ttbarMass vs. leptonic top mass using wrongly matched events (using hadronic b jet); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  
  histo2D["2D_matched_lepTopMass_corr_dR_Lep_B"] = new TH2F("2D_matched_lepTopMass_corr_dR_Lep_B", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using correctly matched events; M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_lepTopMass_wrong_hadrB_dR_Lep_B"] = new TH2F("2D_matched_lepTopMass_wrong_hadrB_dR_Lep_B", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using wrongly matched events (using hadronic b jet); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  
  histo2D["2D_matched_ttbarMass_corr_dR_Lep_B"] = new TH2F("2D_matched_ttbarMass_corr_dR_Lep_B", "delta R between the lepton and the leptonic b jet vs. ttbarMass using correctly matched events; M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_ttbarMass_wrong_hadrB_dR_Lep_B"] = new TH2F("2D_matched_ttbarMass_wrong_hadrB_dR_Lep_B", "delta R between the lepton and the leptonic b jet vs. ttbarMass using wrongly matched events (using hadronic b jet); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  
  histo2D["2D_matched_hadTopMass_lepTopMass_corr"] = new TH2F("2D_matched_hadTopMass_lepTopMass_corr", "Mass of the leptonic top mass (correctly matched) vs. the mass of the hadronic top mass; M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_hadTopMass_lepTopMass_wrong_hadrB"] = new TH2F("2D_matched_hadTopMass_lepTopMass_wrong_hadrB", "Mass of the leptonic top mass (wrongly matched) vs. the mass of the hadronic top mass; M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
//  histo1D["2b_lepTop_mass_notMatched"] = new TH1F("2b_lepTop_mass_notMatched","Reconstructed leptonic top mass of unmatched events with 2 b-tagged jets; M_{lb} [GeV]", 80, 0, 800);
//  histo1D["2b_lepTop_mass_notMatched_4jets"] = new TH1F("2b_lepTop_mass_notMatched","Reconstructed leptonic top mass of unmatched events with 2 b-tagged jets out of 4 in total; M_{lb} [GeV]", 80, 0, 800);
  
  
  
  /// log likelihood
  histo2D["logLikeWidthMass_gen_matched"] = new TH2F("logLikeWidthMass_gen_matched", "-Log Likelihood of generated matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 62, 144.75, 175.75, 295, 0.55, 30.05);
  histo2D["logLikeWidthMass_reco_matched"] = new TH2F("logLikeWidthMass_reco_matched", "-Log Likelihood of reconstructed matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 62, 144.75, 175.75, 595, 0.55, 60.05);
  histo2D["logLikeWidthMass_reco_unmatched"] = new TH2F("logLikeWidthMass_reco_unmatched", "-Log Likelihood of reconstructed unmatched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 62, 144.75, 175.75, 595, 0.55, 60.05);
  histo2D["logLikeWidthMass_gen_matched_zoom"] = new TH2F("logLikeWidthMass_gen_matched_zoom", "-Log Likelihood of generated matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 22, 164.75, 175.75, 295, 0.55, 30.05);
  histo2D["logLikeWidthMass_reco_matched_zoom"] = new TH2F("logLikeWidthMass_reco_matched_zoom", "-Log Likelihood of reconstructed matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 22, 164.75, 175.75, 400, 12.05, 52.05);
  histo2D["logLikeWidthMass_reco_unmatched_zoom"] = new TH2F("logLikeWidthMass_reco_unmatched_zoom", "-Log Likelihood of reconstructed unmatched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 22, 164.75, 175.75, 400, 12.05, 52.05);
  //histo2D["logLikeWidthMass_reco"] = new TH2F("logLikeWidthMass_reco", "-Log Likelihood of reconstructed matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 10, 167.25, 172.25, 35, 0.55, 4.05);  // sample with mt = 169.5
  //histo2D["logLikeWidthMass_gen"] = new TH2F("logLikeWidthMass_gen", "-Log Likelihood of generated matched events VS top mass and top width; M_{t} [GeV]; #Gamma_{t} [GeV]", 10, 167.25, 172.25, 35, 0.55, 4.05);  // sample with mt = 169.5
  
  
  //histo2D["muon_SF_ID"] = new TH2F("muon_SF_ID", "Muon ID scale factors in function of #eta (x) and p_{T} (y); #eta, p_{T} [GeV]", 21, 0, 2.1, 15, 0, 150);
  //histo2D["muon_SF_Iso"] = new TH2F("muon_SF_Iso", "Muon relIso scale factors in function of #eta (x) and p_{T} (y); #eta, p_{T} [GeV]", 21, 0, 2.1, 15, 0, 150);
  //histo2D["muon_SF_Trig"] = new TH2F("muon_SF_Trig", "Muon trigger scale factors in function of #eta (x) and p_{T} (y); #eta, p_{T} [GeV]", 21, 0, 2.1, 15, 0, 150);
  
  
  histo1D["topMass_reco_matched_bin_Infto1p0"] = new TH1F("topMass_reco_matched_bin_Infto1p0","Reconstructed top mass of matched events (genMass < -1#sigma); M_{t} [GeV]", 350, 50, 400);
  histo1D["topMass_reco_matched_bin_1p0to0p5"] = new TH1F("topMass_reco_matched_bin_1p0to0p5","Reconstructed top mass of matched events (-1#sigma < genMass < -0.5#sigma); M_{t} [GeV]", 350, 50, 400);
  histo1D["topMass_reco_matched_bin_0p5"] = new TH1F("topMass_reco_matched_bin_0p5","Reconstructed top mass of matched events (-0.5#sigma < genMass < 0.5#sigma); M_{t} [GeV]", 350, 50, 400);
  histo1D["topMass_reco_matched_bin_0p5to1p0"] = new TH1F("topMass_reco_matched_bin_0p5to1p0","Reconstructed top mass of matched events (0.5#sigma < genMass < 1#sigma); M_{t} [GeV]", 350, 50, 400);
  histo1D["topMass_reco_matched_bin_1p0toInf"] = new TH1F("topMass_reco_matched_bin_1p0toInf","Reconstructed top mass of matched events (genMass > 1#sigma); M_{t} [GeV]", 350, 50, 400);
  
  histo1D["lepTop_mass_matched_corr_bin_Infto1p0"] = new TH1F("lepTop_mass_matched_corr_bin_Infto1p0", "Reconstructed leptonic top mass using correctly matched events (genMass < -1#sigma); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_corr_bin_1p0to0p5"] = new TH1F("lepTop_mass_matched_corr_bin_1p0to0p5", "Reconstructed leptonic top mass using correctly matched events (-1#sigma < genMass < -0.5#sigma); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_corr_bin_0p5"] = new TH1F("lepTop_mass_matched_corr_bin_0p5", "Reconstructed leptonic top mass using correctly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_corr_bin_0p5to1p0"] = new TH1F("lepTop_mass_matched_corr_bin_0p5to1p0", "Reconstructed leptonic top mass using correctly matched events (0.5#sigma < genMass < 1#sigma); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_corr_bin_1p0toInf"] = new TH1F("lepTop_mass_matched_corr_bin_1p0toInf", "Reconstructed leptonic top mass using correctly matched events (genMass > 1#sigma); M_{lb} [GeV]", 80, 0, 800);
  
  histo1D["lepTop_mass_matched_wrong_bin_Infto1p0"] = new TH1F("lepTop_mass_matched_wrong_bin_Infto1p0", "Reconstructed leptonic top mass using wrongly matched events (genMass < -1#sigma); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_wrong_bin_1p0to0p5"] = new TH1F("lepTop_mass_matched_wrong_bin_1p0to0p5", "Reconstructed leptonic top mass using wrongly matched events (-1#sigma < genMass < -0.5#sigma); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_wrong_bin_0p5"] = new TH1F("lepTop_mass_matched_wrong_bin_0p5", "Reconstructed leptonic top mass using wrongly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_wrong_bin_0p5to1p0"] = new TH1F("lepTop_mass_matched_wrong_bin_0p5to1p0", "Reconstructed leptonic top mass using wrongly matched events (0.5#sigma < genMass < 1#sigma); M_{lb} [GeV]", 80, 0, 800);
  histo1D["lepTop_mass_matched_wrong_bin_1p0toInf"] = new TH1F("lepTop_mass_matched_wrong_bin_1p0toInf", "Reconstructed leptonic top mass using wrongly matched events (genMass > 1#sigma); M_{lb} [GeV]", 80, 0, 800);
  
  histo1D["ttbar_mass_matched_corr_bin_Infto1p0"] = new TH1F("ttbar_mass_matched_corr_bin_Infto1p0","Reconstructed mass of the top quark pair using correctly matched events (genMass < -1#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_corr_bin_1p0to0p5"] = new TH1F("ttbar_mass_matched_corr_bin_1p0to0p5","Reconstructed mass of the top quark pair using correctly matched events (-1#sigma < genMass < -0.5#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_corr_bin_0p5"] = new TH1F("ttbar_mass_matched_corr_bin_0p5","Reconstructed mass of the top quark pair using correctly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_corr_bin_0p5to1p0"] = new TH1F("ttbar_mass_matched_corr_bin_0p5to1p0","Reconstructed mass of the top quark pair using correctly matched events (0.5#sigma < genMass < 1#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_corr_bin_1p0toInf"] = new TH1F("ttbar_mass_matched_corr_bin_1p0toInf","Reconstructed mass of the top quark pair using correctly matched events (genMass > 1#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  
  histo1D["ttbar_mass_matched_wrong_bin_Infto1p0"] = new TH1F("ttbar_mass_matched_wrong_bin_Infto1p0","Reconstructed mass of the top quark pair using wrongly matched events (genMass < -1#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong_bin_1p0to0p5"] = new TH1F("ttbar_mass_matched_wrong_bin_1p0to0p5","Reconstructed mass of the top quark pair using wrongly matched events (-1#sigma < genMass < -0.5#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong_bin_0p5"] = new TH1F("ttbar_mass_matched_wrong_bin_0p5","Reconstructed mass of the top quark pair using wrongly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong_bin_0p5to1p0"] = new TH1F("ttbar_mass_matched_wrong_bin_0p5to1p0","Reconstructed mass of the top quark pair using wrongly matched events (0.5#sigma < genMass < 1#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["ttbar_mass_matched_wrong_bin_1p0toInf"] = new TH1F("ttbar_mass_matched_wrong_bin_1p0toInf","Reconstructed mass of the top quark pair using wrongly matched events (genMass > 1#sigma); M_{t#bar{t}} [GeV]", 100, 0, 1000);
  
  histo1D["dR_Lep_B_bin_Infto1p0"] = new TH1F("dR_Lep_B_bin_Infto1p0","Delta R between the lepton and the leptonic b jet (genMass < -1#sigma); #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_Lep_B_bin_1p0to0p5"] = new TH1F("dR_Lep_B_bin_1p0to0p5","Delta R between the lepton and the leptonic b jet (-1#sigma < genMass < -0.5#sigma); #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_Lep_B_bin_0p5"] = new TH1F("dR_Lep_B_bin_0p5","Delta R between the lepton and the leptonic b jet (-0.5#sigma < genMass < 0.5#sigma); #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_Lep_B_bin_0p5to1p0"] = new TH1F("dR_Lep_B_bin_0p5to1p0","Delta R between the lepton and the leptonic b jet (0.5#sigma < genMass < 1#sigma); #Delta R(l,b)", 25, 0, 5);
  histo1D["dR_Lep_B_bin_1p0toInf"] = new TH1F("dR_Lep_B_bin_1p0toInf","Delta R between the lepton and the leptonic b jet (genMass > 1#sigma); #Delta R(l,b)", 25, 0, 5);
  
  
  histo2D["2D_matched_lepTopMass_corr_ttbarMass_corr_bin_Infto1p0"] = new TH2F("2D_matched_lepTopMass_corr_ttbarMass_corr_bin_Infto1p0", "ttbarMass vs. leptonic top mass using correctly matched events (genMass < -1#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_lepTopMass_corr_ttbarMass_corr_bin_1p0to0p5"] = new TH2F("2D_matched_lepTopMass_corr_ttbarMass_corr_bin_1p0to0p5", "ttbarMass vs. leptonic top mass using correctly matched events (-1#sigma < genMass < -0.5#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_lepTopMass_corr_ttbarMass_corr_bin_0p5"] = new TH2F("2D_matched_lepTopMass_corr_ttbarMass_corr_bin_0p5", "ttbarMass vs. leptonic top mass using correctly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_lepTopMass_corr_ttbarMass_corr_bin_0p5to1p0"] = new TH2F("2D_matched_lepTopMass_corr_ttbarMass_corr_bin_0p5to1p0", "ttbarMass vs. leptonic top mass using correctly matched events (0.5#sigma < genMass < 1#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_lepTopMass_corr_ttbarMass_corr_bin_1p0toInf"] = new TH2F("2D_matched_lepTopMass_corr_ttbarMass_corr_bin_1p0toInf", "ttbarMass vs. leptonic top mass using correctly matched events (genMass > 1#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  
  histo2D["2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_Infto1p0"] = new TH2F("2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_Infto1p0", "ttbarMass vs. leptonic top mass using wrongly matched events (genMass < -1#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_1p0to0p5"] = new TH2F("2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_1p0to0p5", "ttbarMass vs. leptonic top mass using wrongly matched events (-1#sigma < genMass < -0.5#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_0p5"] = new TH2F("2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_0p5", "ttbarMass vs. leptonic top mass using wrongly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_0p5to1p0"] = new TH2F("2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_0p5to1p0", "ttbarMass vs. leptonic top mass using wrongly matched events (0.5#sigma < genMass < 1#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_1p0toInf"] = new TH2F("2D_matched_lepTopMass_wrong_ttbarMass_wrong_bin_1p0toInf", "ttbarMass vs. leptonic top mass using wrongly matched events (genMass > 1#sigma); M_{lb} [GeV]; M_{t#bar{t}} [GeV]", 100, 0, 1000, 100, 0, 1000);
      
  histo2D["2D_matched_lepTopMass_corr_dR_Lep_B_bin_Infto1p0"] = new TH2F("2D_matched_lepTopMass_corr_dR_Lep_B_bin_Infto1p0", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using correctly matched events (genMass < -1#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_lepTopMass_corr_dR_Lep_B_bin_1p0to0p5"] = new TH2F("2D_matched_lepTopMass_corr_dR_Lep_B_bin_1p0to0p5", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using correctly matched events (-1#sigma < genMass < -0.5#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_lepTopMass_corr_dR_Lep_B_bin_0p5"] = new TH2F("2D_matched_lepTopMass_corr_dR_Lep_B_bin_0p5", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using correctly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_lepTopMass_corr_dR_Lep_B_bin_0p5to1p0"] = new TH2F("2D_matched_lepTopMass_corr_dR_Lep_B_bin_0p5to1p0", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using correctly matched events (0.5#sigma < genMass < 1#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_lepTopMass_corr_dR_Lep_B_bin_1p0toInf"] = new TH2F("2D_matched_lepTopMass_corr_dR_Lep_B_bin_1p0toInf", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using correctly matched events (genMass > 1#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  
  histo2D["2D_matched_lepTopMass_wrong_dR_Lep_B_bin_Infto1p0"] = new TH2F("2D_matched_lepTopMass_wrong_dR_Lep_B_bin_Infto1p0", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using wrongly matched events (genMass < -1#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_lepTopMass_wrong_dR_Lep_B_bin_1p0to0p5"] = new TH2F("2D_matched_lepTopMass_wrong_dR_Lep_B_bin_1p0to0p5", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using wrongly matched events (-1#sigma < genMass < -0.5#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_lepTopMass_wrong_dR_Lep_B_bin_0p5"] = new TH2F("2D_matched_lepTopMass_wrong_dR_Lep_B_bin_0p5", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using wrongly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_lepTopMass_wrong_dR_Lep_B_bin_0p5to1p0"] = new TH2F("2D_matched_lepTopMass_wrong_dR_Lep_B_bin_0p5to1p0", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using wrongly matched events (0.5#sigma < genMass < 1#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_lepTopMass_wrong_dR_Lep_B_bin_1p0toInf"] = new TH2F("2D_matched_lepTopMass_wrong_dR_Lep_B_bin_1p0toInf", "delta R between the lepton and the leptonic b jet vs. leptonic top mass using wrongly matched events (genMass > 1#sigma); M_{lb} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  
  histo2D["2D_matched_ttbarMass_corr_dR_Lep_B_bin_Infto1p0"] = new TH2F("2D_matched_ttbarMass_corr_dR_Lep_B_bin_Infto1p0", "delta R between the lepton and the leptonic b jet vs. ttbarMass using correctly matched events (genMass < -1#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_ttbarMass_corr_dR_Lep_B_bin_1p0to0p5"] = new TH2F("2D_matched_ttbarMass_corr_dR_Lep_B_bin_1p0to0p5", "delta R between the lepton and the leptonic b jet vs. ttbarMass using correctly matched events (-1#sigma < genMass < -0.5#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_ttbarMass_corr_dR_Lep_B_bin_0p5"] = new TH2F("2D_matched_ttbarMass_corr_dR_Lep_B_bin_0p5", "delta R between the lepton and the leptonic b jet vs. ttbarMass using correctly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_ttbarMass_corr_dR_Lep_B_bin_0p5to1p0"] = new TH2F("2D_matched_ttbarMass_corr_dR_Lep_B_bin_0p5to1p0", "delta R between the lepton and the leptonic b jet vs. ttbarMass using correctly matched events (0.5#sigma < genMass < 1#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_ttbarMass_corr_dR_Lep_B_bin_1p0toInf"] = new TH2F("2D_matched_ttbarMass_corr_dR_Lep_B_bin_1p0toInf", "delta R between the lepton and the leptonic b jet vs. ttbarMass using correctly matched events (genMass > 1#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  
  histo2D["2D_matched_ttbarMass_wrong_dR_Lep_B_bin_Infto1p0"] = new TH2F("2D_matched_ttbarMass_wrong_dR_Lep_B_bin_Infto1p0", "delta R between the lepton and the leptonic b jet vs. ttbarMass using wrongly matched events (genMass < -1#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_ttbarMass_wrong_dR_Lep_B_bin_1p0to0p5"] = new TH2F("2D_matched_ttbarMass_wrong_dR_Lep_B_bin_1p0to0p5", "delta R between the lepton and the leptonic b jet vs. ttbarMass using wrongly matched events (-1#sigma < genMass < -0.5#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_ttbarMass_wrong_dR_Lep_B_bin_0p5"] = new TH2F("2D_matched_ttbarMass_wrong_dR_Lep_B_bin_0p5", "delta R between the lepton and the leptonic b jet vs. ttbarMass using wrongly matched events (-0.5#sigma < genMass < 0.5#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_ttbarMass_wrong_dR_Lep_B_bin_0p5to1p0"] = new TH2F("2D_matched_ttbarMass_wrong_dR_Lep_B_bin_0p5to1p0", "delta R between the lepton and the leptonic b jet vs. ttbarMass using wrongly matched events (0.5#sigma < genMass < 1#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  histo2D["2D_matched_ttbarMass_wrong_dR_Lep_B_bin_1p0toInf"] = new TH2F("2D_matched_ttbarMass_wrong_dR_Lep_B_bin_1p0toInf", "delta R between the lepton and the leptonic b jet vs. ttbarMass using wrongly matched events (genMass > 1#sigma); M_{t#bar{t}} [GeV]; #Delta R(l,b)", 100, 0, 1000, 25, 0, 5);
  
  histo2D["2D_matched_hadTopMass_lepTopMass_corr_bin_Infto1p0"] = new TH2F("2D_matched_hadTopMass_lepTopMass_corr_bin_Infto1p0", "Mass of the leptonic top mass (correctly matched) vs. the mass of the hadronic top mass (genMass < -1#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_hadTopMass_lepTopMass_corr_bin_1p0to0p5"] = new TH2F("2D_matched_hadTopMass_lepTopMass_corr_bin_1p0to0p5", "Mass of the leptonic top mass (correctly matched) vs. the mass of the hadronic top mass (-1#sigma < genMass < -0.5#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_hadTopMass_lepTopMass_corr_bin_0p5"] = new TH2F("2D_matched_hadTopMass_lepTopMass_corr_bin_0p5", "Mass of the leptonic top mass (correctly matched) vs. the mass of the hadronic top mass (-0.5#sigma < genMass < 0.5#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_hadTopMass_lepTopMass_corr_bin_0p5to1p0"] = new TH2F("2D_matched_hadTopMass_lepTopMass_corr_bin_0p5to1p0", "Mass of the leptonic top mass (correctly matched) vs. the mass of the hadronic top mass (0.5#sigma < genMass < 1#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_hadTopMass_lepTopMass_corr_bin_1p0toInf"] = new TH2F("2D_matched_hadTopMass_lepTopMass_corr_bin_1p0toInf", "Mass of the leptonic top mass (correctly matched) vs. the mass of the hadronic top mass (genMass > 1#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  
  histo2D["2D_matched_hadTopMass_lepTopMass_wrong_bin_Infto1p0"] = new TH2F("2D_matched_hadTopMass_lepTopMass_wrong_bin_Infto1p0", "Mass of the leptonic top mass (wrongly matched) vs. the mass of the hadronic top mass (genMass < -1#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_hadTopMass_lepTopMass_wrong_bin_1p0to0p5"] = new TH2F("2D_matched_hadTopMass_lepTopMass_wrong_bin_1p0to0p5", "Mass of the leptonic top mass (wrongly matched) vs. the mass of the hadronic top mass (-1#sigma < genMass < -0.5#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_hadTopMass_lepTopMass_wrong_bin_0p5"] = new TH2F("2D_matched_hadTopMass_lepTopMass_wrong_bin_0p5", "Mass of the leptonic top mass (wrongly matched) vs. the mass of the hadronic top mass (-0.5#sigma < genMass < 0.5#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_hadTopMass_lepTopMass_wrong_bin_0p5to1p0"] = new TH2F("2D_matched_hadTopMass_lepTopMass_wrong_bin_0p5to1p0", "Mass of the leptonic top mass (wrongly matched) vs. the mass of the hadronic top mass (0.5#sigma < genMass < 1#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  histo2D["2D_matched_hadTopMass_lepTopMass_wrong_bin_1p0toInf"] = new TH2F("2D_matched_hadTopMass_lepTopMass_wrong_bin_1p0toInf", "Mass of the leptonic top mass (wrongly matched) vs. the mass of the hadronic top mass (genMass > 1#sigma); M_{t} [GeV]; M_{lb} [GeV]", 100, 0, 1000, 100, 0, 1000);
  
  
  
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
  //CutsSelecTableSemiMu.push_back(string("muon dR(jet) > 0.4"));
  
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
  
  ResolutionFunctions* rf = new ResolutionFunctions(calculateResolutionFunctions);
  
  
  
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
  // documentation at http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
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
  
  
  
  ///////////////////////////////
  ///  Single Muon Selection  ///
  ///////////////////////////////
  
  /// Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO
  
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
  float electronEtaSel = 2.5;
  string electronWP = "Tight";
  
  float electronPTVeto = 15.; // GeV
  float electronEtaVeto = 2.5;
  
  
  
  ///////////////////////
  ///  Jet Selection  ///
  ///////////////////////
  
  /// Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopJME
  
  float jetPT = 20.; // GeV, 30 GeV for selected jets
  float jetEta = 2.4;  // to allow b tagging
  
  
  
  //////////////////////////////////////
  ///  Working points for b tagging  ///
  //////////////////////////////////////
  
  /// Updated 04/03/16, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
  
  float CSVv2Loose =  0.460;
  float CSVv2Medium = 0.800;
  float CSVv2Tight = 0.935;
  
  
  
  ////////////////////////////////
  ///  Define TLorentzVectors  ///
  ////////////////////////////////
  
  // Matching
  vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
  TLorentzVector topQuark, antiTopQuark;
  // Transfer functions
  vector<TLorentzVector> partonTLV, jetTLV;
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
  
  double gammaProp, numTopPropagator;
  double denomTopPropagator_gen_matched, topPropagator_gen_matched;
  double denomTopPropagator_reco_matched, topPropagator_reco_matched;
  double denomTopPropagator_reco_unmatched, topPropagator_reco_unmatched;
  double likelihood_gen_matched[sizeListTopMass][sizeListTopWidth] = {{0}};
  double likelihood_reco_matched[sizeListTopMass][sizeListTopWidth] = {{0}};
  double likelihood_reco_unmatched[sizeListTopMass][sizeListTopWidth] = {{0}};
  
  // Temporarily, until calculated from TTbar sample
  float chi2WMass = 80.385;
  float sigmaChi2WMass = 15;
  float chi2TopMass = 172.5;
  float sigmaChi2TopMass = 40;
  
  int nofChi2FreeWithNoBTaggedJets = 0;
  int nofChi2FreeWith1BTaggedJet = 0;
  int nofChi2FreeWith2BTaggedJets = 0;
  int nofChi2FreeWith3BTaggedJets = 0;
  int nofCorrectlyMatched_chi2Free = 0;
  int nofNotCorrectlyMatched_chi2Free = 0;
  
  
  ////////////////////////////////
  ///  Output file top masses  ///
  ////////////////////////////////
  
  ofstream ofTopMassMatlab;
  ofTopMassMatlab.open("dataTopMassMatlab.txt");
  ofTopMassMatlab << "[ " << endl;
  
  ofstream ofTopMassMath;
  ofTopMassMath.open("dataTopMassMathematica.txt");
  ofTopMassMath << "data = {" << endl;
  
  float mu_orig = 172.070;
  float sigma_orig = 1.14112;
  float mu_new = mu_orig;
  //float sigma_new = sigma_orig/2.;
  float sigma_new = sigma_orig*1.5;
  
  //weight = sigma_orig/sigma_new * TMath::Exp( pow( (x-mu_orig)/sigma_orig , 2 )/2. - pow( (x-mu_new)/sigma_new , 2 )/2. );
  
  /// Bins in genTopMass_matched
  string listBins[5] = {"_bin_Infto1p0", "_bin_1p0to0p5", "_bin_0p5", "_bin_0p5to1p0","_bin_1p0toInf"};
  double listBinEdges[4] = {mu_orig - sigma_orig, mu_orig - sigma_orig/2., mu_orig + sigma_orig/2., mu_orig + sigma_orig};  // bin 1: -Inf -> mu - 1s; bin 2: mu - 1s -> mu - 0.5s; bin 3: mu - 0.5s -> mu + 0.5s; bin 4: mu + 0.5s -> mu + 1s; bin 5: mu + 1s -> Inf
  
  
  
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
    nofEventsJetLeptonCleaned = 0;
    nofChi2FreeWithNoBTaggedJets = 0; nofChi2FreeWith1BTaggedJet = 0; nofChi2FreeWith2BTaggedJets = 0; nofChi2FreeWith3BTaggedJets = 0;
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
    
    
    /// book triggers
    if (applyTriggers) { trigger->bookTriggers(isData);}
    
    
    
    ///////////////////////////////////////////
    ///  Initialise Jet Energy Corrections  ///
    ///////////////////////////////////////////
    
    vector<JetCorrectorParameters> vCorrParam;
    string pathCalJEC = "../TopTreeAnalysisBase/Calibrations/JECFiles/";

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
    //for (unsigned int ievt = 0; ievt < 201; ievt++)
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
      if (((int)nEvents[d])%10000 == 0)
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
      
      
      
      /////////////////
      ///  Pile-up  ///
      /////////////////
      
      
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
      
      
      
      ////////////////////////////
      ///  Include trigger set up here when using data
      ////////////////////////////
      
      bool trigged = false;
      bool fileChanged = false;
      
      /// Fill selection table before trigger
      selecTableSemiMu.Fill(d,0,scaleFactor);
      
      
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
      
      
      
      //////////////////////////////////////
      ///  Jet Energy Scale Corrections  ///
      //////////////////////////////////////
      
      if (applyJER && ! isData)
      {
        if (applyJERdown)    jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus", false);
        else if (applyJERup) jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus", false);
        else                 jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal", false);
        
        
        /// Example how to apply JES systematics
        //jetTools->correctJetJESUnc(init_jets_corrected, "minus", 1);
        //jetTools->correctJetJESUnc(init_jets_corrected, "plus", 1);
        //cout << "JER smeared!!! " << endl;
      }
      
      /// Example how to apply JES systematics
      //if (applyJESdown && ! isData)    jetTools->correctJetJESUnc(init_jets_corrected, "minus", 1);
      //else if (applyJESup && ! isData) jetTools->correctJetJESUnc(init_jets_corrected, "plus", 1);
      
      if (applyJEC)
      {
        jetTools->correctJets(init_jets_corrected, event->fixedGridRhoFastjetAll(), isData);
      }
      
      
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      //Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->fixedGridRhoFastjetAll());
      
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
      
      if (applyJetLeptonCleaning)
      {
        if(verbose > 3) cout << "  - Applying jet/lepton cleaning... " << endl; 
        
        vector<TRootPFJet*> selectedJetsBC;
        selectedJetsBC = selectedJets;
        selectedJets.clear();
        int origSizeJets = selectedJetsBC.size();

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
          if ( origSizeJets != selectedJets.size() ) cout << "--> original = " << origSizeJets  << " after cleaning = " << selectedJets.size() << endl;
        }
        nofEventsJetLeptonCleaned++;
        
      }  // end jet cleaning
      
      
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
      hasExactly4Jets = false;
      has1bjet = false;
      has2bjets = false;
      muonSFID = muonSFIso = muonSFTrig = 1.;
      
      
      /// Continue with selection table
      if (isGoodPV)
      {
        selecTableSemiMu.Fill(d,2,scaleFactor);
        if (selectedMuons.size() == 1)
        {
          selecTableSemiMu.Fill(d,3,scaleFactor);
          /// Apply muon scale factor
          if (applyLeptonSF && ! isData)
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
                //selecTableSemiMu.Fill(d,7,scaleFactor);

                if ( selectedBJets.size() > 0 )
                {
                  selecTableSemiMu.Fill(d,7,scaleFactor);
                  //selecTableSemiMu.Fill(d,8,scaleFactor);
                  has1bjet = true;
                  eventSelected = true;

                  if ( selectedBJets.size() > 1 )
                  {
                    selecTableSemiMu.Fill(d,8,scaleFactor);
                    //selecTableSemiMu.Fill(d,9,scaleFactor);
                    has2bjets = true;
                  }  // at least 2 b-tagged jets

                  if (applyBTagSF && ! isData)
                  {
                    double bTagSF = bTagHistoTool_M->getMCEventWeight(selectedJets);
                    scaleFactor = scaleFactor*bTagSF;
                  }

                }  // at least 1 b-tagged jets

                if ( selectedJets.size() == 4 ) hasExactly4Jets = true;
              }  // at least 4 jets
            }  // no veto electrons
          }  // no additional loose muons (tight muon is also loose muon)
        }  // 1 good muon
      }  // good PV
      
      
//       if ( applyTriggers && ! trigged ) { continue;}
//       selecTableSemiMu.Fill(d,9,scaleFactor);
      
      
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
      
      if ( selectedBJets.size() == 1 ) nofEventsWith1BJet++;
      else if ( selectedBJets.size() > 1 ) nofEventsWith2BJets++;
      
            
      /// B-tagging
      if (calculateBTagSF && ! isData)
      {
        bTagHistoTool_M->FillMCEfficiencyHistos(selectedJets);
      }
      
      
      
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
        
        if (all4PartonsMatched && calculateResolutionFunctions)
        {
          
          partonTLV.clear(); jetTLV.clear();  // vector<TLV>
          //genMuTLV.Clear(); selMuTLV.Clear(); genElTLV.Clear(); selElTLV.Clear();
          
          for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
          {
            /// JetPartonPair[i].first  = jet number
            /// JetPartonPair[i].second = parton number
            /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
            
            partonTLV.push_back(mcParticlesTLV[JetPartonPair[iMatch].second]);
            jetTLV.push_back(selectedJetsTLV[JetPartonPair[iMatch].first]);
          }
          
          rf->fillJets(partonTLV, jetTLV);
          
          if (muonmatched) rf->fillMuon((TLorentzVector) *mcParticles[genmuon], (TLorentzVector) *selectedMuons[0]);
          //if (electronmatched) rf->fillElectron(...)
          
        }  // end rf
        
        
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
            gammaProp = sqrt( pow( listTopMass[jMass], 4 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 ) );
            numTopPropagator = ( 2 * sqrt(2) * listTopMass[jMass] * listTopWidth[jWidth] * gammaProp ) / ( TMath::Pi() * sqrt( pow(listTopMass[jMass], 2) + gammaProp ) );
            
            /// Generated mass
            denomTopPropagator_gen_matched = pow( pow(topMassGen_matched, 2) - pow(listTopMass[jMass], 2), 2 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 );
            
            topPropagator_gen_matched = numTopPropagator/denomTopPropagator_gen_matched;
            
            likelihood_gen_matched[jMass][jWidth] += -TMath::Log10(topPropagator_gen_matched);
            
            /// Reconstructed mass
            denomTopPropagator_reco_matched = pow( pow(topMassReco_matched, 2) - pow(listTopMass[jMass], 2), 2 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 );
            
            topPropagator_reco_matched = numTopPropagator/denomTopPropagator_reco_matched;
            
            likelihood_reco_matched[jMass][jWidth] += -TMath::Log10(topPropagator_reco_matched);
            
          }  /// End loop jWidth
        }  /// End loop jMass
        
        
        /// Binned plots
        string thisBin = "";
        if ( topMassGen_matched < listBinEdges[0] ) thisBin = listBins[0];
        else if ( topMassGen_matched < listBinEdges[1] ) thisBin = listBins[1];
        else if ( topMassGen_matched < listBinEdges[2] ) thisBin = listBins[2];
        else if ( topMassGen_matched < listBinEdges[3] ) thisBin = listBins[3];
        else thisBin = listBins[4];
        
        histo1D[("topMass_reco_matched"+thisBin).c_str()]->Fill(topMassReco_matched);
        
        double weight = 1.;
        //weight = sigma_orig/sigma_new * TMath::Exp( pow( (topMassReco_matched-mu_orig)/sigma_orig , 2 )/2. - pow( (topMassReco_matched-mu_new)/sigma_new , 2 )/2. );  // sigma_orig > sigma_new
        //weight = sigma_new/sigma_orig * TMath::Exp( pow( (topMassReco_matched-mu_new)/sigma_new , 2 )/2. - pow( (topMassReco_matched-mu_orig)/sigma_orig , 2 )/2. );  // sigma_orig < sigma_new
        if ( weight > 150.) weight = 1.;
        
        /// Fill plots
        histo1D["WMass_reco_matched"]->Fill(WMassReco_matched, weight);
        histo1D["topMass_reco_matched"]->Fill(topMassReco_matched, weight);
        histo1D["topMass_gen_matched"]->Fill(topMassGen_matched, weight);
        if ( all4JetsMatched_MCdef_ )
        {
          histo1D["WMass_reco_first4matched"]->Fill(WMassReco_matched, weight);
          histo1D["topMass_reco_first4matched"]->Fill(topMassReco_matched, weight);
          histo1D["topMass_gen_first4matched"]->Fill(topMassGen_matched, weight);
        }
        
        /// Fill offiles
        ofTopMassMatlab << topMassGen_matched << "  ";
        ofTopMassMath << topMassGen_matched<< ", ";
        
        
        /// Leptonic top mass
        if (muonmatched)
        {
          float lepTopMass_matched_corr = (*selectedMuons[0] + *selectedJets[MCPermutation[3]]).M();  // lept b
          float lepTopMass_matched_wrong_hadrB = (*selectedMuons[0] + *selectedJets[MCPermutation[2]]).M();  // hadr b
          
          bool jetMC0 = false, jetMC1 = false;
          float lepTopMass_matched_wrong_noB = -1;
          if ( selectedJets[MCPermutation[0]]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium )
            jetMC0 = true;
          if ( selectedJets[MCPermutation[1]]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium )
            jetMC1 = true;
          
          if ( jetMC0 && jetMC1 )
          {
            if ( selectedJets[MCPermutation[0]]->Pt() > selectedJets[MCPermutation[1]]->Pt() ) jetMC1 = false;
            else jetMC0 = false;
          }
          if (jetMC0)
            lepTopMass_matched_wrong_noB = ( *selectedMuons[0] + *selectedJets[MCPermutation[0]] ).M();
          else if (jetMC1)
            lepTopMass_matched_wrong_noB = ( *selectedMuons[0] + *selectedJets[MCPermutation[1]] ).M();
          
          /// Fill plots
          histo1D["lepTop_mass_matched_corr"]->Fill(lepTopMass_matched_corr, weight);
          histo1D["lepTop_mass_matched_wrong"]->Fill(lepTopMass_matched_wrong_hadrB, weight);
          histo1D["lepTop_mass_matched_wrong_hadrB"]->Fill(lepTopMass_matched_wrong_hadrB, weight);
          if ( jetMC0 || jetMC1)
          {
            histo1D["lepTop_mass_matched_wrong"]->Fill(lepTopMass_matched_wrong_noB, weight);
            histo1D["lepTop_mass_matched_wrong_noB"]->Fill(lepTopMass_matched_wrong_noB, weight);
          }
          
          if (hasExactly4Jets)
          {
            histo1D["lepTop_mass_matched_corr_4jets"]->Fill(lepTopMass_matched_corr, weight);
            histo1D["lepTop_mass_matched_wrong_4jets"]->Fill(lepTopMass_matched_wrong_hadrB, weight);
            if (jetMC0 || jetMC1) histo1D["lepTop_mass_matched_wrong_4jets"]->Fill(lepTopMass_matched_wrong_noB, weight);
          }
          
          /// ttbar mass
          float ttbarMass_matched_corr = topMassReco_matched + lepTopMass_matched_corr;
          float ttbarMass_matched_wrong_hadrB = topMassReco_matched + lepTopMass_matched_wrong_hadrB;
          float ttbarMass_matched_wrong_noB = topMassReco_matched + lepTopMass_matched_wrong_noB;
          
          histo1D["ttbar_mass_matched_corr"]->Fill(ttbarMass_matched_corr, weight);
          histo1D["ttbar_mass_matched_wrong"]->Fill(ttbarMass_matched_wrong_hadrB, weight);
          histo1D["ttbar_mass_matched_wrong_hadrB"]->Fill(ttbarMass_matched_wrong_hadrB, weight);
          if ( jetMC0 || jetMC1)
          {
            histo1D["ttbar_mass_matched_wrong"]->Fill(ttbarMass_matched_wrong_noB, weight);
            histo1D["ttbar_mass_matched_wrong_noB"]->Fill(ttbarMass_matched_wrong_noB, weight);
          }
          if (hasExactly4Jets)
          {
            histo1D["ttbar_mass_matched_corr_4jets"]->Fill(ttbarMass_matched_corr, weight);
            histo1D["ttbar_mass_matched_wrong_4jets"]->Fill(ttbarMass_matched_wrong_hadrB, weight);
            if (jetMC0 || jetMC1) histo1D["ttbar_mass_matched_wrong_4jets"]->Fill(ttbarMass_matched_wrong_noB, weight);
          }
          
          /// dR(lep,b)
          double dRLepB_matched = ROOT::Math::VectorUtil::DeltaR((TLorentzVector) *selectedMuons[0], (TLorentzVector) *selectedJets[MCPermutation[3]]);  // lept b
          histo1D["dR_Lep_B"]->Fill(dRLepB_matched, weight);
          
          /// 2D plots: M_lb, M_ttbar, dR(lep,b)
          histo2D["2D_matched_lepTopMass_corr_ttbarMass_corr"]->Fill(lepTopMass_matched_corr, ttbarMass_matched_corr, weight);
          histo2D["2D_matched_lepTopMass_wrong_ttbarMass_wrong_hadrB"]->Fill(lepTopMass_matched_wrong_hadrB, ttbarMass_matched_wrong_hadrB, weight);
          histo2D["2D_matched_lepTopMass_corr_dR_Lep_B"]->Fill(lepTopMass_matched_corr, dRLepB_matched, weight);
          histo2D["2D_matched_lepTopMass_wrong_hadrB_dR_Lep_B"]->Fill(lepTopMass_matched_wrong_hadrB, dRLepB_matched, weight);
          histo2D["2D_matched_ttbarMass_corr_dR_Lep_B"]->Fill(ttbarMass_matched_corr, dRLepB_matched, weight);
          histo2D["2D_matched_ttbarMass_wrong_hadrB_dR_Lep_B"]->Fill(ttbarMass_matched_wrong_hadrB, dRLepB_matched, weight);
          
          histo2D["2D_matched_hadTopMass_lepTopMass_corr"]->Fill(topMassReco_matched, lepTopMass_matched_corr, weight);
          histo2D["2D_matched_hadTopMass_lepTopMass_wrong_hadrB"]->Fill(topMassReco_matched, lepTopMass_matched_wrong_hadrB, weight);
          
          
          
          /// Fill binned histos (binned in genTopMass_matched)
          histo1D[("lepTop_mass_matched_corr"+thisBin).c_str()]->Fill(lepTopMass_matched_corr);
          histo1D[("lepTop_mass_matched_wrong"+thisBin).c_str()]->Fill(lepTopMass_matched_wrong_hadrB);
          histo1D[("ttbar_mass_matched_corr"+thisBin).c_str()]->Fill(ttbarMass_matched_corr);
          histo1D[("ttbar_mass_matched_wrong"+thisBin).c_str()]->Fill(ttbarMass_matched_wrong_hadrB);
          histo1D[("dR_Lep_B"+thisBin).c_str()]->Fill(dRLepB_matched);
          
          histo2D[("2D_matched_lepTopMass_corr_ttbarMass_corr"+thisBin).c_str()]->Fill(lepTopMass_matched_corr, ttbarMass_matched_corr);
          histo2D[("2D_matched_lepTopMass_wrong_ttbarMass_wrong"+thisBin).c_str()]->Fill(lepTopMass_matched_wrong_hadrB, ttbarMass_matched_wrong_hadrB);
          histo2D[("2D_matched_lepTopMass_corr_dR_Lep_B"+thisBin).c_str()]->Fill(lepTopMass_matched_corr, dRLepB_matched);
          histo2D[("2D_matched_lepTopMass_wrong_dR_Lep_B"+thisBin).c_str()]->Fill(lepTopMass_matched_wrong_hadrB, dRLepB_matched);
          histo2D[("2D_matched_ttbarMass_corr_dR_Lep_B"+thisBin).c_str()]->Fill(ttbarMass_matched_corr, dRLepB_matched);
          histo2D[("2D_matched_ttbarMass_wrong_dR_Lep_B"+thisBin).c_str()]->Fill(ttbarMass_matched_wrong_hadrB, dRLepB_matched);
          histo2D[("2D_matched_hadTopMass_lepTopMass_corr"+thisBin).c_str()]->Fill(topMassReco_matched, lepTopMass_matched_corr);
          histo2D[("2D_matched_hadTopMass_lepTopMass_wrong"+thisBin).c_str()]->Fill(topMassReco_matched, lepTopMass_matched_wrong_hadrB);
          
        }  // end muonmatched
        
      }  // end TT && matched
      
      
      
      if (doChi2)
      {
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



        ////////////////
        ///   CHI2   ///
        ////////////////

        int labelsRecoFree[3] = {-9999, -9999, -9999};		// 0 = hadronic b-jet, 1,2 = light jets.
        float recoWMass, recoTopMass, recoTopMass2;

        string listPlotInit_b_chi2[2] = {"1b_", "2b_"};
        string plotInit_b_chi2 = "";
        if ( selectedBJets.size() == 1 ) { plotInit_b_chi2 = listPlotInit_b_chi2[0];}
        else { plotInit_b_chi2 = listPlotInit_b_chi2[1];}


        ///---------------///
        ///   FREE CHI2   ///
        ///---------------///

        float WTerm_chi2Free, topTerm_chi2Free, chi2Free;
        float smallestChi2Free = 999999.;
        for (int ijet = 0; ijet < 4; ijet++)
        {
          for (int jjet = ijet+1; jjet < 4; jjet++)
          {
            for (int kjet = 0; kjet < 4; kjet++)
            {
              if ( ijet != kjet && jjet != kjet )
              {
                recoWMass = (*selectedJets[ijet] + *selectedJets[jjet]).M();
                recoTopMass = (*selectedJets[ijet] + *selectedJets[jjet] + *selectedJets[kjet]).M();

                WTerm_chi2Free = pow( (recoWMass - chi2WMass)/sigmaChi2WMass, 2);
                topTerm_chi2Free = pow( (recoTopMass - chi2TopMass)/sigmaChi2TopMass, 2);

                chi2Free = WTerm_chi2Free + topTerm_chi2Free;


                if (chi2Free < smallestChi2Free)
                {
                  smallestChi2Free = chi2Free;
                  labelsRecoFree[0] = kjet;
                  labelsRecoFree[1] = ijet;
                  labelsRecoFree[2] = jjet;
                }
              }
            }
          }
        }

        int nofBs = 0;
        if ( selectedJets[labelsRecoFree[2]]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium ) { nofBs++;}
        if ( selectedJets[labelsRecoFree[3]]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium ) { nofBs++;;}
        if ( selectedJets[labelsRecoFree[1]]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium ) { nofBs++;;}

        if ( nofBs == 3 ) { nofChi2FreeWith3BTaggedJets++;}
        else if ( nofBs == 2 ) { nofChi2FreeWith2BTaggedJets++;}
        else if ( nofBs == 1 ) { nofChi2FreeWith1BTaggedJet++;}
        else { nofChi2FreeWithNoBTaggedJets++;}


        if (labelsRecoFree[0] != -9999 && labelsRecoFree[1] != -9999 && labelsRecoFree[2] != -9999)
        {
          float Wmass_chi2Free = (*selectedJets[labelsRecoFree[1]] + *selectedJets[labelsRecoFree[2]]).M();
          float hadtopmass_chi2Free = (*selectedJets[labelsRecoFree[0]] + *selectedJets[labelsRecoFree[1]] +  *selectedJets[labelsRecoFree[2]]).M();
          float hadtoppt_chi2Free = (*selectedJets[labelsRecoFree[0]] + *selectedJets[labelsRecoFree[1]] +  *selectedJets[labelsRecoFree[2]]).Pt();
          float hadtopht_chi2Free = selectedJets[labelsRecoFree[0]]->Pt() + selectedJets[labelsRecoFree[1]]->Pt() + selectedJets[labelsRecoFree[2]]->Pt();

          /// Make likelihoods
          if ( dataSetName.find("TT") == 0 )
          {
            for (unsigned int jMass = 0; jMass < sizeListTopMass; jMass++)
            {
              for (unsigned int jWidth = 0; jWidth < sizeListTopWidth; jWidth++)
              {
                gammaProp = sqrt( pow( listTopMass[jMass], 4 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 ) );
                numTopPropagator = ( 2 * sqrt(2) * listTopMass[jMass] * listTopWidth[jWidth] * gammaProp ) / ( TMath::Pi() * sqrt( pow(listTopMass[jMass], 2) + gammaProp ) );

                denomTopPropagator_reco_unmatched = pow( pow(hadtopmass_chi2Free, 2) - pow(listTopMass[jMass], 2), 2 ) + pow( listTopMass[jMass] * listTopWidth[jWidth], 2 );

                topPropagator_reco_unmatched = numTopPropagator/denomTopPropagator_reco_unmatched;

                likelihood_reco_unmatched[jMass][jWidth] += -TMath::Log10(topPropagator_reco_unmatched);

              }  /// End loop jWidth
            }  /// End loop jMass
          }


          //Fill histos
          if ( dataSetName.find("TT") == 0 )
          {
            histo1D[(plotInit_b_chi2+"Chi2Free_W_mass_reco_notMatched").c_str()]->Fill(Wmass_chi2Free);
            histo1D[(plotInit_b_chi2+"Chi2Free_top_mass_reco_notMatched").c_str()]->Fill(hadtopmass_chi2Free);
            if (hasExactly4Jets)
            {
              histo1D[(plotInit_b_chi2+"Chi2Free_top_mass_reco_notMatched_4jets").c_str()]->Fill(hadtopmass_chi2Free);
            }
          }

          /// Leptonic top mass
          vector<TRootPFJet*> bJetsAfterChi2;
          for (int i = 0; i < selectedJets.size(); i++)
          {
            if ( i != labelsRecoFree[0] && i != labelsRecoFree[1] && i != labelsRecoFree[2] 
                && selectedJets[i]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2Medium )
              bJetsAfterChi2.push_back(selectedJets[i]);
          }

          if ( bJetsAfterChi2.size() > 0 )
          {
            double min_M_lb_afterChi2 = 99999., Mlb_temp_afterChi2 = -1.;
            double dRLepB_afterChi2 = -1.;
            for (unsigned int i = 0; i < bJetsAfterChi2.size(); i++)
            {
              Mlb_temp_afterChi2 = (*selectedMuons[0] + *bJetsAfterChi2[i]).M();
              if ( Mlb_temp_afterChi2 < min_M_lb_afterChi2 )
              {
                min_M_lb_afterChi2 = Mlb_temp_afterChi2;
                dRLepB_afterChi2 = ROOT::Math::VectorUtil::DeltaR( (TLorentzVector)*bJetsAfterChi2[i], (TLorentzVector)*selectedMuons[0]);
              }
            }

            double ttbarMass_afterChi2 = min_M_lb_afterChi2 + hadtopmass_chi2Free;


            if ( dataSetName.find("TT") == 0 )
            {
              histo1D[(plotInit_b_chi2+"Chi2Free_lepTop_mass_notMatched").c_str()]->Fill(min_M_lb_afterChi2);
              histo1D[(plotInit_b_chi2+"Chi2Free_ttbar_mass_notMatched").c_str()]->Fill(ttbarMass_afterChi2);
              histo1D[(plotInit_b_chi2+"Chi2Free_dR_lep_b_notMatched").c_str()]->Fill(dRLepB_afterChi2);
              histo2D[(plotInit_b_chi2+"2D_unmatched_hadTopMass_lepTopMass").c_str()]->Fill(hadtopmass_chi2Free, min_M_lb_afterChi2);

              if (hasExactly4Jets)
              {
                histo1D[(plotInit_b_chi2+"Chi2Free_lepTop_mass_notMatched_4jets").c_str()]->Fill(min_M_lb_afterChi2);
                histo1D[(plotInit_b_chi2+"Chi2Free_ttbar_mass_notMatched_4jets").c_str()]->Fill(ttbarMass_afterChi2);
                histo1D[(plotInit_b_chi2+"Chi2Free_dR_lep_b_notMatched_4jets").c_str()]->Fill(dRLepB_afterChi2);
              }

              if ( min_M_lb_afterChi2 < 200 )
              {
                histo1D[(plotInit_b_chi2+"Chi2Free_hadTop_mass_notMatched_mlb_cut").c_str()]->Fill(hadtopmass_chi2Free);
                histo1D[(plotInit_b_chi2+"Chi2Free_ttbar_mass_notMatched_mlb_cut").c_str()]->Fill(ttbarMass_afterChi2);
                histo1D[(plotInit_b_chi2+"Chi2Free_dR_lep_b_notMatched_mlb_cut").c_str()]->Fill(dRLepB_afterChi2);
                if (hasExactly4Jets)
                {
                  histo1D[(plotInit_b_chi2+"Chi2Free_hadTop_mass_notMatched_4jets_mlb_cut").c_str()]->Fill(hadtopmass_chi2Free);
                  histo1D[(plotInit_b_chi2+"Chi2Free_ttbar_mass_notMatched_4jets_mlb_cut").c_str()]->Fill(ttbarMass_afterChi2);
                  histo1D[(plotInit_b_chi2+"Chi2Free_dR_lep_b_notMatched_4jets_mlb_cut").c_str()]->Fill(dRLepB_afterChi2);
                }
              }
            }  // end TT

          }  // end bjetsAfterChi2 > 0

        }  // end labels chi2 filled



        ////////////////////////////////
        // CHECK MATCHED COMBINATION  //
        ////////////////////////////////

        if ( dataSetName.find("TT") == 0 && all4PartonsMatched )
        {
          if ( labelsRecoFree[0] == MCPermutation[2] 
              && ( (labelsRecoFree[1] == MCPermutation[0] && labelsRecoFree[2] == MCPermutation[1]) 
                || (labelsRecoFree[1] == MCPermutation[1] && labelsRecoFree[2] == MCPermutation[0]) ) )
            nofCorrectlyMatched_chi2Free++;
          else
            nofNotCorrectlyMatched_chi2Free++;
        }
        
      }  // end doChi2
      
      
      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  /// Loop on events
    
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofEventsWith1BJet << " events with 1 b tagged jet." << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofEventsWith2BJets << " events with 2 b tagged jets." << endl;
    if (doChi2)
    {
      cout << "Number of matched events: " << nofMatchedEvents << endl;
      cout << "Correctly matched for chi2Free:     " << nofCorrectlyMatched_chi2Free << endl;
      cout << "Not correctly matched for chi2Free: " << nofNotCorrectlyMatched_chi2Free << endl;
      if ( nofCorrectlyMatched_chi2Free != 0 || nofNotCorrectlyMatched_chi2Free != 0 )
        cout << "   ===> This means that " << 100*(float)nofCorrectlyMatched_chi2Free / (float)(nofCorrectlyMatched_chi2Free + nofNotCorrectlyMatched_chi2Free) << "% is correctly matched." << endl;
    }
    
    if (isData)
    {
      weightMuonHLTv2 = ((double) nofEventsHLTv2) / ((double) (nofEventsHLTv2 + nofEventsHLTv3));
      weightMuonHLTv3 = ((double) nofEventsHLTv3) / ((double) (nofEventsHLTv2 + nofEventsHLTv3));
      cout << "The muon trigger scale factors will be scaled by " << weightMuonHLTv2 << " for HLTv2 and " << weightMuonHLTv3 << " for HLTv3." << endl;
    }
    cout << endl;
    
    cout << "Chi2: Jets that give best reconstruction of top quark: " << endl;
    cout << "               " << setw(7) << right << nofChi2FreeWithNoBTaggedJets << " chi2 combinations have no b-tagged jets" << endl;
    cout << "               " << setw(7) << right << nofChi2FreeWith1BTaggedJet << " chi2 combinations have 1 b-tagged jet" << endl;
    cout << "               " << setw(7) << right << nofChi2FreeWith2BTaggedJets << " chi2 combinations have 2 b-tagged jets" << endl;
    cout << "               " << setw(7) << right << nofChi2FreeWith3BTaggedJets << " chi2 combinations have 3 b-tagged jets" << endl;
    cout << endl;
    
    if (applyJetLeptonCleaning)
    {
      cout << "Number of (not necessarily selected) events that had jet cleaning " << nofEventsJetLeptonCleaned << endl;
      cout << endl;
    }
    
    
    /// offiles
    ofTopMassMatlab << " ];" << endl;
    ofTopMassMath << "};" << endl;
    ofTopMassMatlab.close();
    ofTopMassMath.close();
    
    
    /// Fill histogram log likelihood && Transfer functions
    if ( dataSetName.find("TT") == 0 )
    {
      for (unsigned int jMass = 0; jMass < sizeListTopMass; jMass++)
      {
        for (unsigned int jWidth = 0; jWidth < sizeListTopWidth; jWidth++)
        {
          histo2D["logLikeWidthMass_gen_matched"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_gen_matched[jMass][jWidth]);
          histo2D["logLikeWidthMass_gen_matched_zoom"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_gen_matched[jMass][jWidth]);
          histo2D["logLikeWidthMass_reco_matched"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_reco_matched[jMass][jWidth]);
          histo2D["logLikeWidthMass_reco_matched_zoom"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_reco_matched[jMass][jWidth]);
          histo2D["logLikeWidthMass_reco_unmatched"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_reco_unmatched[jMass][jWidth]);
          histo2D["logLikeWidthMass_reco_unmatched_zoom"]->Fill(listTopMass[jMass], listTopWidth[jWidth], likelihood_reco_unmatched[jMass][jWidth]);
          
        }
      }
      
      
      /// Transfer functions
      if (calculateResolutionFunctions)
      {
        string rfFileName = "PlotsForResolutionFunctions.root";
        TFile *foutRF = new TFile(rfFileName.c_str(), "RECREATE");
        foutRF->cd();

        rf->writeHistograms();

        foutRF->Close();

        rf->writeTable(rfFileName);

        delete foutRF;
      }
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
  
  /// To write plots b tagging:
  delete bTagHistoTool_M;
  
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  //string pathPNG = "PlotsOneFourth/";
  mkdir(pathPNG.c_str(),0777);
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile ((pathPNG+rootFileName).c_str(), "RECREATE");
  
  ///Write histograms
  fout->cd();
  
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
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
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
  selecTableSemiMu.Write((pathPNG+selectiontableMu).c_str(), true, true, true, true, true, true, false);
  
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
