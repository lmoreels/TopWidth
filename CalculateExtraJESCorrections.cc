#include <TStyle.h>
#include <TPaveText.h>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <array>
#include <TRandom3.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>


// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
//#include "../macros/Style.C"

// user defined
#include "Tools/interface/ResolutionFunctions.h"
#include "Tools/interface/EventReweighting.h"


using namespace std;
using namespace TopTree;


bool test = false;
bool testHistos = false;
bool makePlots = true;


string datasets[] = {"nominal", "tuneup", "tunedown", "isrup", "isrdown", "fsrup", "fsrdown", "hdampup", "hdampdown",  "herwig"};
string dataSetName[] = {"TT_nominal", "TT_tune_up", "TT_tune_down", "TT_isr_up", "TT_isr_down", "TT_fsr_up", "TT_fsr_down", "TT_hdamp_up", "TT_hdamp_down", "TT_herwigpp"};
int nDatasets = sizeof(datasets)/sizeof(datasets[0]);


string ntupleDateMC = "171121";
string ntupleDateSyst = "171123";
int verbose = 2;

string pathNtuples = "";
string pathNtuplesMC = "";
string pathNtuplesSyst = "";
string pathOutput = "";

bool isHerwig = false;

int nofHardSelected = 0;
int nofMETCleaned = 0;
int nofAfterDRmincut = 0;
int nofAfterDRmaxcut = 0;
int nofAfterLastCut = 0;


/// Lumi per data era
double lumi_runBCDEF = 19.67550334113;  // 1/fb
double lumi_runGH = 16.146177597883;  // 1/fb
double Luminosity = (lumi_runBCDEF + lumi_runGH)*1000;  // 1/pb

///  Working points for b tagging  // Updated 13/04/17, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
double CSVv2Loose  = 0.5426;
double CSVv2Medium = 0.8484;
double CSVv2Tight  = 0.9535;


// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TGraphErrors*> graph;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;


/// Function prototypes
struct HighestPt
{
  bool operator()( TLorentzVector j1, TLorentzVector j2 ) const
  {
    return j1.Pt() > j2.Pt() ;
  }
};

struct greater
{
  template<class T>
  bool operator()(T const &a, T const &b) const { return a > b; }
};


string ConvertDoubleToString(double Number);
string DotReplace(double var);
string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();
bool fexists(const char *filename);
void InitSetUp();
void GetMetaData(TTree* tree, bool isData);
void InitTree(TTree* tree, bool isData);
void InitHisto1D();
void InitHisto2D();
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearVars();
void ClearObjects();
long GetNEvents(TTree* fChain, string var, bool isData);
long GetNEvents(TTree* fChain, string var, unsigned int index, bool isData);



// Declaration of leaf types
Int_t           run_num;
Long64_t        evt_num;
Int_t           lumi_num;
Int_t           nvtx;
Int_t           npu;
Int_t           nJets;
Int_t           jet_nConstituents[20];   //[nJets]
Int_t           jet_nChConstituents[20];   //[nJets]
Int_t           jet_charge[20];   //[nJets]
Double_t        jet_pt[20];   //[nJets]
Double_t        jet_phi[20];   //[nJets]
Double_t        jet_eta[20];   //[nJets]
Double_t        jet_E[20];   //[nJets]
Double_t        jet_M[20];   //[nJets]
Double_t        jet_bdiscr[20];   //[nJets]
Double_t        jet_hadronFlavour[20];   //[nJets]
Int_t           nGenJets;
Int_t           genJet_pdgId[33];   //[nGenJets]
Float_t         genJet_charge[33];   //[nGenJets]
Double_t        genJet_pt[33];   //[nGenJets]
Double_t        genJet_phi[33];   //[nGenJets]
Double_t        genJet_eta[33];   //[nGenJets]
Double_t        genJet_E[33];   //[nGenJets]
Double_t        genJet_M[33];   //[nGenJets]

Long64_t        nEvents;
Long64_t        nEventsSel;


// List of branches
TBranch        *b_run_num;   //!
TBranch        *b_evt_num;   //!
TBranch        *b_lumi_num;   //!
TBranch        *b_nJets;   //!
TBranch        *b_jet_nConstituents;   //!
TBranch        *b_jet_nChConstituents;   //!
TBranch        *b_jet_charge;   //!
TBranch        *b_jet_pt;   //!
TBranch        *b_jet_phi;   //!
TBranch        *b_jet_eta;   //!
TBranch        *b_jet_E;   //!
TBranch        *b_jet_M;   //!
TBranch        *b_jet_bdiscr;   //!
TBranch        *b_jet_hadronFlavour;   //!
TBranch        *b_nGenJets;   //!
TBranch        *b_genJet_pdgId;   //!
TBranch        *b_genJet_charge;   //!
TBranch        *b_genJet_pt;   //!
TBranch        *b_genJet_phi;   //!
TBranch        *b_genJet_eta;   //!
TBranch        *b_genJet_E;   //!
TBranch        *b_genJet_M;   //!

TBranch        *b_nEvents;   //!
TBranch        *b_nEventsSel;   //!


long nEventsDataSet;

vector<unsigned int> bJetId;


/// Define TLVs
TLorentzVector jet, genjet;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
vector<TLorentzVector> genJets;

/// Matching
int tempIndex;
double minDR = 9999., tempDR;
vector<pair<int,int>> match;

double barrelCut = 1.3;
double ptBins[] = {30., 50., 70., 100., 140., 200., 250., 300.};
int nPtBins = sizeof(ptBins)/sizeof(ptBins[0]);
string cats[4] = {"_B_barrel", "_B_endcap", "_Light_barrel", "_Light_endcap"};
int nCats = sizeof(cats)/sizeof(cats[0]);

int hBins = 50;
string histoName, graphName, graphNameRef;
string binLowEdge, binUpEdge;

bool isB, isBarrel;
double tempJetPt, ptRatio;

TF1 *gaus;
int thisPoint, nPoints;
double fitmin, fitmax;
double meanPt, meanRatio, errPt, errRatio, SF, errSF;


/// Meta
string strSyst = "";
double eqLumi;

string dateString = "";

int main(int argc, char* argv[])
{
  dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  cout << "* The following corrections are applied:    *" << endl;
  cout << "*   - Jet Energy Corrections                *" << endl;
  cout << "*   - Jet Energy Resolution:                *" << endl;
  cout << "*   - Jet/lepton Cleaning                   *" << endl;
  cout << "*   - MET Cleaning                          *" << endl;
  cout << "*********************************************" << endl;
  
  
  clock_t start = clock();
  
  if ( argc > 1 )
  {
    for (int i = 1; i < argc; i++)
    {
      if ( strcmp(argv[i], "--makePlots") == 0 || strcmp(argv[i], "--plots") == 0 )
      {
        makePlots = true;
      }
    }
  }
  
  
  if (testHistos)
  {
    makePlots = true;
  }
  
  pathOutput = "JESCorrections/";
  mkdir(pathOutput.c_str(),0777);
  if (makePlots)
  {
    if (testHistos)
    {
      pathOutput += "test/";
      mkdir(pathOutput.c_str(),0777);
    }
    // Give timestamp to output path
    pathOutput += dateString+"/";
    mkdir(pathOutput.c_str(),0777);
  }
  
  
  pathNtuplesMC = "NtupleOutput/MergedTuples/mu/"+ntupleDateMC+"/";
  pathNtuplesSyst = "NtupleOutput/MergedTuples/mu/"+ntupleDateSyst+"/";
  cout << "Using Ntuples from " << ntupleDateMC << " for nominal MC and " << ntupleDateSyst << " for systematics." << endl;
  if (testHistos) cout << "Testing histogram consistency..." << endl;
  
  
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  if (makePlots)
  {
    InitHisto1D();
    InitHisto2D();
  }
  
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  //string dataSetName;
  double timePerDataSet[nDatasets];
  
  int nEntries;
  bool skipEvent = false;
  
  
  
  ////////////////////////////////////
  /// Loop over datasets
  ////////////////////////////////////
  
  //for (int d = 0; d < 1; d++)
  for (int d = 0; d < nDatasets; d++)   //Loop through datasets
  {
    clock_t startDataSet = clock();
    
    ClearMetaData();
    
    isHerwig = false;
    if ( dataSetName[d].find("herwig") != std::string::npos || dataSetName[d].find("HERWIG") != std::string::npos || dataSetName[d].find("Herwig") != std::string::npos ) isHerwig = true;
    
    cout << "   Dataset " << d << ": " << dataSetName[d] << endl;
    
    if ( dataSetName[d].find("nominal") != std::string::npos ) pathNtuples = pathNtuplesMC;
    else pathNtuples = pathNtuplesSyst;
    
    string ntupleFileName = pathNtuples+"Ntuples_"+dataSetName[d]+".root";
    
    tFileMap[(dataSetName[d]).c_str()] = new TFile(ntupleFileName.c_str(),"READ"); //create TFile for each dataset
    
    string tTreeName = "tree";
    string tStatsTreeName = "stats";
    
    /// Get meta data
    tStatsTree[(dataSetName[d]).c_str()] = (TTree*)tFileMap[(dataSetName[d]).c_str()]->Get(tStatsTreeName.c_str());
    GetMetaData(tStatsTree[(dataSetName[d]).c_str()], false);
    
    tStatsTree[(dataSetName[d]).c_str()]->GetEntry(0);
    
    
    
    
    
    /// Get data
    tTree[(dataSetName[d]).c_str()] = (TTree*)tFileMap[(dataSetName[d]).c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
    nEntries = (int)tTree[(dataSetName[d]).c_str()]->GetEntries();
    cout << "                nEntries  : " << nEntries << endl;
    
    
    
    // Set branch addresses and branch pointers
    InitTree(tTree[(dataSetName[d]).c_str()], false);
    
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    int endEvent = nEntries;
    if (test || testHistos) endEvent = 10;
    for (int ievt = 0; ievt < endEvent; ievt++)
    {
      ClearObjects();
      
      if (ievt%10000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
      
      
      /// Load event
      tTree[(dataSetName[d]).c_str()]->GetEntry(ievt);
      
      
      
      //////////////////////
      ///  Fill objects  ///
      //////////////////////
      
      for (int iGenJet = 0; iGenJet < nGenJets; iGenJet++)
      {
        genjet.Clear();
        genjet.SetPtEtaPhiE(genJet_pt[iGenJet], genJet_eta[iGenJet], genJet_phi[iGenJet], genJet_E[iGenJet]);
        genJets.push_back(genjet);
      }
      
      for (int iJet = 0; iJet < nJets; iJet++)
      {
        jet.Clear();
        jet.SetPtEtaPhiE(jet_pt[iJet], jet_eta[iJet], jet_phi[iJet], jet_E[iJet]);
        selectedJets.push_back(jet);
      }
      
//       for (int iJet = 0; iJet < selectedJets.size(); iJet++)
//       {
//         if ( jet_bdiscr[iJet] > CSVv2Medium )
//         {
//           selectedBJets.push_back(selectedJets[iJet]);
//           bJetId.push_back(iJet);  /// selectedBJets[j] = selectedJets[bJetId[j]]
//         }
//       }
//       
//       if ( selectedJets.size() > 4 ) continue;
//       nofHardSelected++;
//       
//       skipEvent = false;
//       for (int iJet = 0; iJet < selectedJets.size(); iJet++)
//       {
//         for (int jJet = iJet+1; jJet < selectedJets.size(); jJet++)
//         {
//           tempDR = ROOT::Math::VectorUtil::DeltaR(selectedJets[iJet], selectedJets[jJet]);
//           if ( tempDR < 0.6 ) skipEvent = true;
//         }
//       }
//       if (skipEvent) continue;
//       nofAfterDRmincut++;
      
      
      
      /////////////////////////
      ///  GENJET MATCHING  ///
      /////////////////////////
      
      vector<bool> mLocked(genJets.size(),false);   // when locked, genJet is already matched to a recoJet
      for (int iJet = 0; iJet < selectedJets.size(); iJet++)
      {
        minDR = 9999.;
        tempIndex = -1.;
        for (int iGenJet = 0; iGenJet < genJets.size(); iGenJet++)
        {
          if (mLocked[iGenJet]) continue;
          
          tempDR = ROOT::Math::VectorUtil::DeltaR(selectedJets[iJet], genJets[iGenJet]);
          if( tempDR < 0.4 && tempDR < minDR )
          {
            minDR = tempDR;
            tempIndex = iGenJet;
          }
        }
        if ( tempIndex != -1 )
        {
          mLocked[tempIndex] = true;
          match.push_back(pair<int,int>(iJet,tempIndex));
        }
      }
      
      
      /////////////////////
      ///  FILL HISTOS  ///
      /////////////////////
      
      for (int iMatch = 0; iMatch < match.size(); iMatch++)
      {
        tempJetPt = selectedJets[match[iMatch].first].Pt();
        ptRatio = tempJetPt / genJets[match[iMatch].second].Pt();
        
        histoName = "";
        if ( jet_bdiscr[match[iMatch].first] > CSVv2Medium ) histoName += "_B";
        else histoName += "_Light";
        if ( fabs(jet_eta[match[iMatch].first]) < barrelCut ) histoName += "_barrel";
        else histoName += "_endcap";
        
        for (int iBin = 1; iBin < nPtBins; iBin++)
        {
          if ( tempJetPt > ptBins[iBin-1] && tempJetPt <= ptBins[iBin] )
          {
            histoName += "_Pt"+ConvertDoubleToString(ptBins[iBin-1])+"To"+ConvertDoubleToString(ptBins[iBin]);
            histo1D["Ratio_"+dataSetName[d]+histoName]->Fill(ptRatio);
            continue;
          }
        }
      }
      
      
      
    }  // end loop events
    
    cout << endl;
    
    
    
    ///////////////////
    ///  FIT PLOTS  ///
    ///////////////////
    
    for (int iBin = 1; iBin < nPtBins; iBin++)
    {
      binLowEdge = ConvertDoubleToString(ptBins[iBin-1]);
      binUpEdge = ConvertDoubleToString(ptBins[iBin]);
      
      errPt = (ptBins[iBin] - ptBins[iBin-1])/2.;
      meanPt = ptBins[iBin-1] + errPt;
      
      for (int iCat = 0; iCat < nCats; iCat++)
      {
        graphName = "Ratio_"+dataSetName[d]+cats[iCat];
        histoName = graphName+"_Pt"+binLowEdge+"To"+binUpEdge;
        histo1D[histoName]->Scale(1./histo1D[histoName]->Integral(0,hBins+1));
        fitmin = histo1D[histoName]->GetXaxis()->GetBinLowEdge(histo1D[histoName]->FindFirstBinAbove(0.04));
        if ( fitmin < 0. ) fitmin = 0.;
        fitmax = histo1D[histoName]->GetXaxis()->GetBinUpEdge(histo1D[histoName]->FindLastBinAbove(0.04));
        if ( fitmax > 2. ) fitmax = 2.;
        
        gaus = new TF1("gaus","gaus",fitmin,fitmax);
        
        histo1D[histoName]->Fit("gaus","R");
        meanRatio = gaus->GetParameter(1);
        errRatio = gaus->GetParError(1);
        
        if ( graph[graphName] == NULL ) graph[graphName] = new TGraphErrors();
        thisPoint = graph[graphName]->GetN();
        graph[graphName]->SetPoint(thisPoint, meanPt, meanRatio);
        graph[graphName]->SetPointError(thisPoint, errPt, errRatio);
      }
    }
    
    
    
    
    
    tFileMap[(dataSetName[d]).c_str()]->Close();
    
    timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
    
  }  // end loop datasets
  
  
  cout << endl << "Processing time per dataset: " << endl;
  for (unsigned int d = 0; d < nDatasets; d++)
  {
    cout << dataSetName[d] << ": " << timePerDataSet[d] << " s" << endl;
  }
  
  
  
  //////////////////
  ///  Make SFs  ///
  //////////////////
  
  for (unsigned int d = 1; d < nDatasets; d++)
  {
    cout << "Making scale factors for dataset " << d << " : " << dataSetName[d] << endl;
    
    for (int iCat = 0; iCat < nCats; iCat++)
    {
      graph["SF_"+dataSetName[d]+cats[iCat]] =  new TGraphErrors();
      
      graphName = "Ratio_"+dataSetName[d]+cats[iCat];
      graphNameRef = "Ratio_"+dataSetName[0]+cats[iCat];
      
      double *pts, *ptErrs, *dRatios, *dRatioErrs, *refRatios, *refRatioErrs;
      pts = graph[graphName]->GetX();
      ptErrs = graph[graphName]->GetEX();
      dRatios = graph[graphName]->GetY();
      dRatioErrs = graph[graphName]->GetEY();
      refRatios = graph[graphNameRef]->GetY();
      refRatioErrs = graph[graphNameRef]->GetEY();
      
      nPoints = graph[graphName]->GetN();
      for (int iPoint = 0; iPoint < nPoints; iPoint++)
      {
        meanPt = pts[iPoint];
        errPt  = ptErrs[iPoint];
        
        SF = refRatios[iPoint]/dRatios[iPoint];
        errSF = sqrt( pow(refRatioErrs[iPoint]/dRatios[iPoint], 2) + pow(dRatioErrs[iPoint]*refRatios[iPoint]/(dRatios[iPoint]*dRatios[iPoint]), 2) );
        
        thisPoint = graph["SF_"+dataSetName[d]+cats[iCat]]->GetN();
        graph["SF_"+dataSetName[d]+cats[iCat]]->SetPoint(thisPoint, meanPt, SF);
        graph["SF_"+dataSetName[d]+cats[iCat]]->SetPointError(thisPoint, errPt, errSF);
      }
    }  // end cats
  }
  
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  string rootFileName = "PtRatioPlots.root";

  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile ((pathOutput+rootFileName).c_str(), "RECREATE");
  cout << "   Output file is " << pathOutput+rootFileName << endl;
  
  ///Write histograms
  fout->cd();
  
  for (std::map<std::string,TGraphErrors*>::const_iterator it = graph.begin(); it != graph.end(); it++)
  {
    TGraphErrors *temp = it->second;
    temp->SetName((it->first).c_str());
    temp->SetTitle((it->first).c_str());
    temp->Write();
  }
  
  // 1D
  //TDirectory* th1dir = fout->mkdir("1D_histograms");
  //th1dir->cd();
  gStyle->SetOptStat(1111);
  for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathOutput+it->first+".png").c_str() );
  }
  
//   // 2D
//   TDirectory* th2dir = fout->mkdir("2D_histograms");
//   th2dir->cd();
//   gStyle->SetPalette(55);
//   for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
//   {
//     TH2F *temp = it->second;
//     temp->Write();
//     TCanvas* tempCanvas = TCanvasCreator(temp, it->first, "colz");
//     tempCanvas->SaveAs( (pathOutput+it->first+".png").c_str() );
//   }
  
  
  
  
  fout->Close();
  
  delete fout;
  
  double time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  if ( time >= 60 )
  {
    int mins = time/60;
    double secs = time - mins*60;
    
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
  
}  // end main



/// Functions

string ConvertDoubleToString(double Number)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

string DotReplace(double var)
{
  string str = ConvertDoubleToString(var);
  replace(str.begin(), str.end(), '.', 'p');
  return str;
}

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

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile.good();
}


void GetMetaData(TTree* tree, bool isData)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  
  tree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
  tree->SetBranchAddress("nEventsSel", &nEventsSel, &b_nEventsSel);
}

void InitTree(TTree* tree, bool isData)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  
  tree->SetBranchAddress("run_num", &run_num, &b_run_num);
  tree->SetBranchAddress("evt_num", &evt_num, &b_evt_num);
  tree->SetBranchAddress("lumi_num", &lumi_num, &b_lumi_num);
  tree->SetBranchAddress("nJets", &nJets, &b_nJets);
  tree->SetBranchAddress("jet_nConstituents", jet_nConstituents, &b_jet_nConstituents);
  tree->SetBranchAddress("jet_nChConstituents", jet_nChConstituents, &b_jet_nChConstituents);
  tree->SetBranchAddress("jet_charge", jet_charge, &b_jet_charge);
  tree->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
  tree->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
  tree->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
  tree->SetBranchAddress("jet_E", jet_E, &b_jet_E);
  tree->SetBranchAddress("jet_M", jet_M, &b_jet_M);
  tree->SetBranchAddress("jet_bdiscr", jet_bdiscr, &b_jet_bdiscr);
  tree->SetBranchAddress("jet_hadronFlavour", jet_hadronFlavour, &b_jet_hadronFlavour);
  tree->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets);
  tree->SetBranchAddress("genJet_pdgId", genJet_pdgId, &b_genJet_pdgId);
  tree->SetBranchAddress("genJet_charge", genJet_charge, &b_genJet_charge);
  tree->SetBranchAddress("genJet_pt", genJet_pt, &b_genJet_pt);
  tree->SetBranchAddress("genJet_phi", genJet_phi, &b_genJet_phi);
  tree->SetBranchAddress("genJet_eta", genJet_eta, &b_genJet_eta);
  tree->SetBranchAddress("genJet_E", genJet_E, &b_genJet_E);
  tree->SetBranchAddress("genJet_M", genJet_M, &b_genJet_M);
}

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  string binLowEdge, binUpEdge;
  
  for (int d = 0; d < nDatasets; d++)
  {
    for (int iBin = 1; iBin < nPtBins; iBin++)
    {
      binLowEdge = ConvertDoubleToString(ptBins[iBin-1]);
      binUpEdge = ConvertDoubleToString(ptBins[iBin]);
      for (int iCat = 0; iCat < nCats; iCat++)
      {
        histo1D["Ratio_"+dataSetName[d]+cats[iCat]+"_Pt"+binLowEdge+"To"+binUpEdge] = new TH1F(("Ratio_"+dataSetName[d]+cats[iCat]+"_Pt"+binLowEdge+"To"+binUpEdge).c_str(), "; p_{T,jet} / p_{T,genjet}", hBins, 0., 2.);
      }
    }
  }
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
}

//bool jet_matching(TLorentzVector jet, TLorentzVector genjet, double resolution)
bool jet_matching(TLorentzVector jet, TLorentzVector genjet)
{
  bool is_matched = false;
  float dR = ROOT::Math::VectorUtil::DeltaR(jet,genjet);
  
  //float pt = jet.pt();
  //float dPt = jet.pt() - genjet.pt();
  
  // Matching in dR (0.4 = AK4)
  if ( dR <= 0.4/2. )
  {
    // Matching in pT (vs resolution)
    //if ( dPt <= 3*resolution*pt)
    //{
      is_matched = true;
    //}
  }
  return is_matched;
}

void ClearMetaData()
{
  nEvents = 0;
  nEventsSel = 0;
  
  strSyst = "";
  nEventsDataSet = 0;
  eqLumi = 1.;
  
  nofHardSelected = 0;
  nofAfterDRmincut = 0;
  nofAfterDRmaxcut = 0;
  nofAfterLastCut = 0;
}

void ClearLeaves()
{
  run_num = -1;
  evt_num = -1;
  lumi_num = -1;
  nJets = -1;
  for (Int_t i = 0; i < 20; i++)
  {
    jet_charge[i] = 0;
    jet_nConstituents[i] = -1;
    jet_nChConstituents[i] = -1;
    jet_pt[i] = 0.;
    jet_phi[i] = 0.;
    jet_eta[i] = 0.;
    jet_E[i] = 0.;
    jet_M[i] = 0.;
    jet_bdiscr[i] = -1.;
  }
  nGenJets = 0.;
  for (Int_t i = 0; i < 20; i++)
  {
    genJet_pdgId[i] = 0.;
    genJet_charge[i] = 0.;
    genJet_pt[i] = 0.;
    genJet_phi[i] = 0.;
    genJet_eta[i] = 0.;
    genJet_E[i] = 0.;
    genJet_M[i] = 0.;
  }
}

void ClearTLVs()
{
  jet.Clear();
  genjet.Clear();
  selectedJets.clear();
  selectedBJets.clear();
  genJets.clear();
}

void ClearVars()
{
  bJetId.clear();
  match.clear();
  
  isB = false;
  isBarrel = false;
  
  tempIndex = -1;
  tempDR = 9999.;
  minDR = 9999.;
  
  tempJetPt = 0.;
  ptRatio = -1.;
}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearVars();
}

long GetNEvents(TTree* fChain, string var, bool isData)
{
  GetMetaData(fChain, isData);
  long varNew = 0;
  for (unsigned int iEntry = 0; iEntry < fChain->GetEntries(); iEntry++)
  {
    fChain->GetEntry(iEntry);
    varNew += (fChain->FindLeaf(var.c_str()))->GetValueLong64();
  }
  
  return varNew;
}

long GetNEvents(TTree* fChain, string var, unsigned int index, bool isData)
{
  GetMetaData(fChain, isData);
  long varNew = 0;
  for (unsigned int iEntry = 0; iEntry < fChain->GetEntries(); iEntry++)
  {
    fChain->GetEntry(iEntry);
    varNew += (fChain->FindLeaf(var.c_str()))->GetValueLong64(index);
  }
  
  return varNew;
}

