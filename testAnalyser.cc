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
#include "../TopTreeAnalysisBase/Selection/interface/Selection.h"
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
  int nb_bTags = 0;
  float genTopMass = 173;
  float genTopWidth = 1.4915;
  
  
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
  //map<string,TH2F*> histo2D;
  
  histo1D["muon_pT"] = new TH1F("muon_pT","Transverse momentum of the muon; p_{T} [GeV]", 200, 0, 200);
  histo1D["muon_Eta"] = new TH1F("muon_Eta","Pseudorapidity of the muon; #eta", 60, -3, 3);
  histo1D["leadingJet_pT"] = new TH1F("leadingJet_pT","Transverse momentum of the leading jet; p_{T} [GeV]", 800, 0, 800);
  histo1D["Ht_4leadingJets"] = new TH1F("Ht_4leadingJets","Scalar sum of transverse momenta of the 4 leading jets; H_{T} [GeV]", 120, 0, 1200);
  histo1D["W_Mass_Reco_matched"] = new TH1F("W_Mass_Reco_matched","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 500, 0, 500);
  histo1D["top_Mass_Reco_matched"] = new TH1F("top_Mass_Reco_matched","Reconstructed top mass of matched events; M_{t} [GeV]", 500, 0, 500);
  
  
  
  ////////////////////////////////////
  ///  MultiSamplePlot
  ////////////////////////////////////
  
  map<string,MultiSamplePlot*> MSPlot;
  
  MSPlot["muon_pT"] = new MultiSamplePlot(datasets, "muon_pT", 22, 0, 440, "p_{T} [GeV]");
  MSPlot["muon_Eta"] = new MultiSamplePlot(datasets, "muon_Eta", 60, -3, 3, "Eta");
  MSPlot["leadingJet_pT"] = new MultiSamplePlot(datasets, "leadingJet_pT", 40, 0, 800, "p_{T} [GeV]");
  MSPlot["Ht_4leadingJets"] = new MultiSamplePlot(datasets,"Ht_4leadingJets", 120, 0, 1200, "H_{T} [GeV]");
  
  
  
  ////////////////////////////////////
  ///  Selection table
  ////////////////////////////////////
  
  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  //CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  //CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
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
  
  
  
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
  
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    if ( (datasets[d]->Name()).find("TT") != 0 ) { cout << "Only consider TTbar at the moment" << endl; continue; }
    
    cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
    nofSelectedEvents = 0;
    
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
    
    //for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    for (unsigned int ievt = 0; ievt < 2000; ievt++)
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
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
      
      // the default selection is fine for normal use - if you want a special selection you can use the functions here
      selection.setJetCuts(20,2.5,0.01,1.,0.98,0.3,0.1); //  void setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon);
      //selection.setMuonCuts(1,2.5,0.2,0.5,0.3,1,0.5,5,0);
      selection.setMuonCuts(20,2.5,0.12,0.2,0.3,1,0.5,5,0); // void setMuonCuts(float Pt, float Eta, float RelIso, float d0, float DRJets, int NMatchedStations, float Dz, int NTrackerLayersWithMeas, int NValidPixelHits);
      selection.setElectronCuts(20,2.5,0.1,2.0,0.5,0.4,0); // void setElectronCuts(float Pt, float Eta, float RelIso, float d0, float MVAId, float DRJets, int MaxMissingHits);
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setLooseElectronCuts(20,2.5,0.15,0.);
      
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
      
      vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
      //vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(selectedJets);
      //vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
      //vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(30,2.5,0.2);
      
      //if (selectedJets.size() >= 4)
      //  if (selectedJets[3]->Pt() < 30) selectedJets.clear();
      
      vector<TRootMCParticle*> mcParticles;
      
      if ( dataSetName.find("TT") == 0 )
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      eventSelected = false;
      has1bjet = false;
      has2bjets = false;
      nb_bTags = 0;
      
      selecTableSemiMu.Fill(d,0,scaleFactor);
      /// At the moment do not use trigger
      selecTableSemiMu.Fill(d,1,scaleFactor);
      if (isGoodPV)
      {
        selecTableSemiMu.Fill(d,2,scaleFactor);
        if (selectedMuons.size() == 1)
        {
          selecTableSemiMu.Fill(d,3,scaleFactor);
          //if (vetoMuons.size() == 1) {
          //  selecTableSemiMu.Fill(d,4,scaleFactor);
          //  if (vetoElectronsSemiMu.size() == 0) {
          //    selecTableSemiMu.Fill(d,5,scaleFactor);
              if ( selectedJets.size() >= 4 )
              {
                //selecTableSemiMu.Fill(d,6,scaleFactor);
                selecTableSemiMu.Fill(d,4,scaleFactor);
                eventSelected = true;
                
                for (unsigned int i = 0; i < selectedJets.size(); i++)
                {
                  if (selectedJets[i]->btag_combinedSecondaryVertexBJetTags() > 0.679) nb_bTags++;
                }
                		
                if ( nb_bTags >= 1 )
                {
                  //selecTableSemiMu.Fill(d,7,scaleFactor);
                  selecTableSemiMu.Fill(d,5,scaleFactor);
                  has1bjet = true;
                  
                  if ( nb_bTags >= 2 )
                  {
                    //selecTableSemiMu.Fill(d,8,scaleFactor);
                    selecTableSemiMu.Fill(d,6,scaleFactor);
                    has2bjets = true;
                  }
                }

              }  // at least 4 jets
          //  }  // no loose electrons
          //}  // no loose muons
        }  // 1 good muon
      }  // good PV
      
      
      
      /// Do some stuff with selected events
      
      if (! eventSelected )
      {
        //cout << "Event no. " << ievt << " was not selected. " << endl;
        continue;
      }
      
      nofSelectedEvents++;
      
      if (verbose > 2)
        cout << "  Event " << ievt << " is selected" << endl;
      
      
      
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
        sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)
        
        vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
        TLorentzVector topQuark, antiTopQuark;
        
        bool muPlusFromTop = false, muMinusFromTop = false;
        mcParticlesMatching_.clear();
        
        for (unsigned int i = 0; i < mcParticles.size(); i++)
        {
          if (verbose > 3)
            cout << i << "  Status: " << mcParticles[i]->status() << "  pdgId: " << mcParticles[i]->type() << "  Mother: " << mcParticles[i]->motherType() << "  Granny: " << mcParticles[i]->grannyType() << endl;
          
          //if ( mcParticles[i]->status() != 3) continue;  // Pythia8: 3 is not correct status anymore!!
          
          if ( mcParticles[i]->type() == pdgID_top )
            topQuark = *mcParticles[i];
          else if( mcParticles[i]->type() == -pdgID_top )
            antiTopQuark = *mcParticles[i];
					
          if ( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )		// mu-, W-, tbar
          {
            muMinusFromTop = true;
            genmuon = i;  // FIX ME: More than 1 mcparticle corresponds!
            if ( verbose > 2 && (ievt == 1134 || ievt == 1903) )
            //if (ievt == 406 || ievt == 803)
              cout << i << "  Status: " << mcParticles[i]->status() << "  pdgId: " << mcParticles[i]->type() << "  Mother: " << mcParticles[i]->motherType() << "  Granny: " << mcParticles[i]->grannyType() << "  Pt: " << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
          }
          if ( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )		// mu+, W+, t
          {
            muPlusFromTop = true;
            genmuon = i;  // FIX ME: More than 1 mcparticle corresponds!
            if ( verbose > 2 && (ievt == 1134 || ievt == 1903) )
            //if (ievt == 406 || ievt == 803)
              cout << i << "  Status: " << mcParticles[i]->status() << "  pdgId: " << mcParticles[i]->type() << "  Mother: " << mcParticles[i]->motherType() << "  Granny: " << mcParticles[i]->grannyType() << "  Pt: " << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
	    		}
          
          if ( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 )  //light/b quarks, 6 should stay hardcoded, OR gluon
          {
            mcParticlesTLV.push_back(*mcParticles[i]);
            mcParticlesMatching_.push_back(mcParticles[i]);
            if ( verbose > 2 && (ievt == 1134 || ievt == 1903) )
            //if (ievt == 406 || ievt == 803)
              cout << i << "  Status: " << mcParticles[i]->status() << "  pdgId: " << mcParticles[i]->type() << "  Mother: " << mcParticles[i]->motherType() << "  Granny: " << mcParticles[i]->grannyType() << "  Pt: " << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
          }
          
        }
        
        // take all the selectedJets_ to study the radiation stuff, selectedJets_ are already ordened in decreasing Pt()
        for (unsigned int i = 0; i < selectedJets.size(); i++)
          selectedJetsTLV.push_back(*selectedJets[i]);
        
        if (verbose > 2)
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
                exit(1);
              }
            }
          }
          if ( fabs(mcParticlesMatching_[j]->type()) == 5 )
          {
            if ( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -pdgID_top )
                || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == pdgID_top ) )  // if mu+ (top decay leptonic) and mother is antitop ---> hadronic b
            {
              hadronicBJet_ = JetPartonPair[i];
              MCPermutation[2] = JetPartonPair[i].first;
            }
            else if ( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == pdgID_top )
              || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == -pdgID_top ) )
            {
              leptonicBJet_ = JetPartonPair[i];
              MCPermutation[3] = JetPartonPair[i].first;
            }
          }
        }  /// End loop over Jet Parton Pairs
        
        
        if (hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999 && hadronicBJet_.first != 9999 && leptonicBJet_.first != 9999) {
          
          all4PartonsMatched = true;
          nofMatchedEvents++;
          if (hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4 && leptonicBJet_.first < 4)
            all4JetsMatched_MCdef_ = true;
	  		}
        else { cout << "Size JetPartonPair: " << JetPartonPair.size() << ". Not all partons matched!" << endl; }
        
        if (hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4)
          hadronictopJetsMatched_MCdef_ = true;
        if (genmuon != -9999 && ROOT::Math::VectorUtil::DeltaR( (TLorentzVector)*mcParticles[genmuon], (TLorentzVector)*selectedMuons[0]) < 0.3)  // FIX ME: More than 1 mcparticle corresponds to genmuon!
          muonmatched = true;
        
        
      }  /// End matching
      
      
      
      //////////////////////////////////
      ///  TOP PROPAGATOR (MATCHED)  ///
      //////////////////////////////////
      
      if ( all4PartonsMatched )
      {
        /// MCPermutation = JetPartonPair[i].first  = jet number
        ///                 JetPartonPair[i].second = parton number
        /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
        
        float WMassReco = (*selectedJets[MCPermutation[0]] + *selectedJets[MCPermutation[1]]).M();
        float topMassReco = (*selectedJets[MCPermutation[0]] + *selectedJets[MCPermutation[1]] + *selectedJets[MCPermutation[2]]).M();
        
        /// Fill plots
        histo1D["W_Mass_Reco_matched"]->Fill(WMassReco);
        histo1D["top_Mass_Reco_matched"]->Fill(topMassReco);
      }
      
      
      
      ////////////////////
      ///  FILL PLOTS  ///
      ////////////////////
      
      double HT = selectedJets[0]->Pt()+selectedJets[1]->Pt()+selectedJets[2]->Pt()+selectedJets[3]->Pt();
      
      if ( dataSetName.find("TT") == 0 )
      {
        histo1D["muon_pT"]->Fill(selectedMuons[0]->Pt());
        histo1D["muon_Eta"]->Fill(selectedMuons[0]->Eta());
        histo1D["leadingJet_pT"]->Fill(selectedJets[0]->Pt());
        histo1D["Ht_4leadingJets"]->Fill(HT);
      }
      
      
      /// Fill MSPlots
      
      MSPlot["muon_pT"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["muon_Eta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor); 
      MSPlot["leadingJet_pT"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["Ht_4leadingJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
      
      
      
//       /// Find b jets with highest pT
//       
//       int labelBTag1 = -9999;
//       int labelBTag2 = -9999;
//       float PtBTag1 = -9999.;
//       float PtBTag2 = -9999.;
//       for (unsigned int i = 0; i < selectedJets.size(); i++) {
//         if ( ! has1bjet ) break;
//         if (selectedJets[i]->btag_combinedSecondaryVertexBJetTags() > 0.679) {		// CSVM
//           if ( ! has2bjets ) {
//             labelBTag1 = i;
//             PtBTag1 = selectedJets[i]->Pt();
//             break;
//           }
//           else {
//             if (selectedJets[i]->Pt() > PtBTag1) {
//               // Save previous as second best
//               if(labelBTag1 >= 0){
//                 labelBTag2 = labelBTag1;
//                 PtBTag2 = PtBTag1;
//               }
//               // Keep new one
//               labelBTag1 = i;
//               PtBTag1 = selectedJets[i]->Pt();
//             }
//             else if (selectedJets[i]->Pt() > PtBTag2) {
//               labelBTag2 = i;
//               PtBTag2 = selectedJets[i]->Pt();
//             }
//           }
//         }
//       }
      
      
      
      
      
      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  //loop on events
    
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    cout << "Number of matched events: " << nofMatchedEvents << endl;
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  // Loop on datasets
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  string pathPNG = "Plots/";
  mkdir(pathPNG.c_str(),0777);
  
//   ///Write histograms
//   fout->cd();
//   for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
//   {
//     //cout << "MSPlot: " << it->first << endl;
//     MultiSamplePlot *temp = it->second;
//     string name = it->first;
//     temp->Draw(name, 0, false, false, false, 1);
//     temp->Write(fout, name, true, pathPNG+"MSPlot/");
//   }
//   
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
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " s to run the program" << endl;
    
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
  
}
