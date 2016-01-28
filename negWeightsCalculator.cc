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
  
  
  string rootFileName = "negWeights_output_"+dateString+".root";
  string pathPNG = "Plots_negWeights_"+dateString+"/";
  
  
  
  /////////////////////
  ///  Configuration
  /////////////////////
  
  bool hasNegWeight = false;
  int nofNegWeights = 0;
  int nofPosWeights = 0;
  
  
  /// xml file
  string xmlFileName ="config/topWidth_negWeights.xml";
  
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
  
  
  
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;
  
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets(datasets, xmlfile); cout << "Number of datasets: " << datasets.size() << endl;
  for (unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = 1;
  cout << " - Analysis environment luminosity is "<< Luminosity << endl;
  
  //nof selected events
  //double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  Double_t mc_baseweight = 0;
  
  
  
  ////////////////////////////////////
  ///  Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
	
  map<string,TH1F*> histo1D;
  string histoname = "";
  
  histo1D["sumWeights_ST_t"] = new TH1F("sumWeights_ST_t","Sum of the weights of the ST t-channel amc@nlo sample", 1, -0.5, 0.5);
  histo1D["sumWeights_Wjets"] = new TH1F("sumWeights_Wjets","Sum of the weights of the W+jets amc@nlo sample", 1, -0.5, 0.5);
  histo1D["sumWeights_DYjets"] = new TH1F("sumWeights_DYjets","Sum of the weights of the DY+jets amc@nlo sample", 1, -0.5, 0.5);
  
  
  
  ////////////////////////////////////
  ///  MultiSamplePlot
  ////////////////////////////////////
  
  map<string,MultiSamplePlot*> MSPlot;
  
  MSPlot["nloWeight"] = new MultiSamplePlot(datasets,"nloWeight", 200, -2.0, 2.0, "weights for amc@nlo samples");
  MSPlot["weightIndex"] = new MultiSamplePlot(datasets,"weightIndex", 5, -2.5, 2.5, "0: None; 1: scale_variation 1; 2: Central scale variation 1");
  
  
  
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
  
  for (unsigned int d = 0; d < datasets.size(); d++)
  { 
    nofNegWeights = 0;
    nofPosWeights = 0;
    double nloSF = 1.;
    double sumWeights = 0.;
    int sumWeightsInt = 0.; 
    string dataSetName = datasets[d]->Name();
    
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << "/ title : " << datasets[d]->Title() << endl;
      cout << "      -> Equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
      cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    }
    
    histoname = "";
    if ( dataSetName.find("ST_t") == 0 ) histoname = "sumWeights_ST_t";
    else if ( dataSetName.find("WJets") == 0 ) histoname = "sumWeights_Wjets";
    else if ( dataSetName.find("DY") == 0 ) histoname = "sumWeights_DYjets";
    
    
    //open files and load
    cout << "LoadEvent" << endl;
    treeLoader.LoadDataset(datasets[d], anaEnv);
    cout << "LoadEvent" << endl;
    
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    
    
    if (verbose > 1)
      //cout << "	Loop over events " << endl;
      cout << "	Loop over events  (" << datasets[d]->NofEvtsToRunOver() << ")" << endl;
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    //for (unsigned int ievt = 0; ievt < 1000; ievt++)
    {
      
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
      
      if (ievt%50000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
        
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
      
      datasets[d]->eventTree()->LoadTree(ievt);
      
      
      hasNegWeight = false;
      if ( event->getWeight(1001) != -9999. )
      {
        mc_baseweight = event->getWeight(1001)/abs(event->originalXWGTUP());
        //mc_scaleupweight = event->getWeight(1005)/abs(event->originalXWGTUP());
        //mc_scaledownweight = event->getWeight(1009)/abs(event->originalXWGTUP());
        if ( mc_baseweight >= 0 ) 
        {
          nofPosWeights++;
          MSPlot["weightIndex"]->Fill(1., datasets[d], false, Luminosity);
        }
        else
        {
          hasNegWeight = true;
          nofNegWeights++;
          MSPlot["weightIndex"]->Fill(-1., datasets[d], false, Luminosity);
        }
      }
      if ( event->getWeight(1) != -9999. )
      {
        mc_baseweight = event->getWeight(1)/abs(event->originalXWGTUP());
        //mc_scaleupweight = event->getWeight(5)/abs(event->originalXWGTUP());
        //mc_scaledownweight = event->getWeight(9)/abs(event->originalXWGTUP());
        if ( mc_baseweight >= 0 )
        {
          nofPosWeights++;
          MSPlot["weightIndex"]->Fill(2., datasets[d], false, Luminosity);
        }
        else
        {
          hasNegWeight = true;
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
      sumWeightsInt += (int) mc_baseweight;
      
      if ( histoname != "" )
        histo1D[histoname.c_str()]->Fill(0., mc_baseweight);
      
      
      
      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  /// Loop on events
    
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofPosWeights << " events with positive weights and " << nofNegWeights << " events with negative weights." << endl;
    cout << "                          " << nofPosWeights - nofNegWeights << " net events" << endl;
    cout << "The sum of the weights is " << sumWeights << ", or " << sumWeightsInt << " when the entries are forced to be integers." << endl;
      
    /// Determine scale factor due to negative weights
    nloSF = ((double) (nofPosWeights - nofNegWeights))/((double) (nofPosWeights + nofNegWeights));
    cout << "This corresponds to an event scale factor of " << nloSF << endl;
    
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  /// Loop on datasets
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  mkdir(pathPNG.c_str(),0777);
  
  ///Write histograms
  fout->cd();
  for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    //cout << "MSPlot: " << it->first << endl;
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(name, 0, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
    temp->Write(fout, name, true, pathPNG, "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
  }
  
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
