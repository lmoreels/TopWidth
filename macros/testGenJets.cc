////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to test genJet properties     /////
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
  
  string rootFileName = "testGenJets_output.root";
  
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
  histo1D["muon_eta"] = new TH1F("muon_Eta","Pseudorapidity of the muon; #eta", 60, -3, 3);
  histo1D["leadingJet_pT"] = new TH1F("leadingJet_pT","Transverse momentum of the leading jet; p_{T} [GeV]", 800, 0, 800);
  histo1D["Ht_4leadingJets"] = new TH1F("Ht_4leadingJets","Scalar sum of transverse momenta of the 4 leading jets; H_{T} [GeV]", 120, 0, 1200);
  
  histo1D["genJets_pT"] = new TH1F("genJets_pT","Transverse momentum of the genJet; p_{T} [GeV]", 800, 0, 800);
  histo1D["genJets_eta"] = new TH1F("genJets_eta","Pseudorapidity of the genJet; #eta", 100, -5, 5);
  histo1D["genJets_nConstit"] = new TH1F("genJets_nConstit","Number of constituents of genJets; # constituents; # genJets", 101, -0.5, 100.5);
  histo1D["genJets_maxDistance"] = new TH1F("genJets_maxDistance","Maximum distance between genJet and constituent ; Distance; # genJets", 100, 0, 10);
  histo1D["genJets_emEnergy"] = new TH1F("genJets_emEnergy","EM energy of genJet; Energy [GeV]; # genJets", 200, 0, 4000);
  histo1D["genJets_hadEnergy"] = new TH1F("genJets_hadEnergy","Hadronic energy of genJet; Energy [GeV]; # genJets", 200, 0, 4000);
  histo1D["genJets_invEnergy"] = new TH1F("genJets_invEnergy","Invisible energy of genJet; Energy [GeV]; # genJets", 200, 0, 4000);
  histo1D["genJets_n90"] = new TH1F("genJets_n90","Number of constituents that carry 90% of energy of genJets; # constituents; # genJets", 101, -0.5, 100.5);
  histo1D["genJets_n60"] = new TH1F("genJets_n60","Number of constituents that carry 60% of energy of genJets; # constituents; # genJets", 101, -0.5, 100.5);
  
  histo1D["diff_jets_pT"] = new TH1F("diff_jets_pT","Vectorial difference between the transverse momentum of jets and their corresponding genJet; p_{T} [GeV]", 200, 0, 50);
  histo1D["diff_jets_scal_pT"] = new TH1F("diff_jets_scal_pT","Scalar difference between the transverse momentum of jets and their corresponding genJet; p_{T} [GeV]", 400, -50, 50);
  histo1D["diff_jets_eta"] = new TH1F("diff_jets_eta","Difference between the pseudorapidity of jets and their corresponding genJet; #eta", 100, -5, 5);
  histo1D["diff_jets_scal_eta"] = new TH1F("diff_jets_scal_eta","Scalar difference between the pseudorapidity of jets and their corresponding genJet; #eta", 100, -5, 5);
  histo1D["diff_jets_phi"] = new TH1F("diff_jets_phi","Difference between the phi angle of jets and their corresponding genJet; #phi", 80, -4, 4);
  histo1D["diff_jets_scal_phi"] = new TH1F("diff_jets_scal_phi","Scalar difference between the phi angle of jets and their corresponding genJet; #phi", 80, -4, 4);
  histo1D["diff_jets_deltaR"] = new TH1F("diff_jets_deltaR","Difference between the deltaR of jets and their corresponding genJet; #Delta R", 60, -0.5, 1);
  
  
  
  ////////////////////////////////////
  ///  Selection table
  ////////////////////////////////////
  
  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  //CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
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
    
    //if ( dataSetName.find("TT") != 0 ) continue;
      
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
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
    //for (unsigned int ievt = 0; ievt < 2000; ievt++)
    {
      
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
      
      if (ievt%1000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
      
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
      
      if (! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
        genjets = treeLoader.LoadGenJet(ievt,false);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      
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
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(25, 2.5, 0.12, "Tight", "Spring15");  // GetSelectedMuons(float PtThr, float etaThr, float relIso, string WorkingPoint, string ProductionCampaign)
      vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(20, 2.5, "Tight", "Spring15_50ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
      vector<TRootMuon*> vetoMuons = selection.GetSelectedMuons(20, 2.5, 0.2, "Loose", "Spring15");  // GetSelectedMuons(float PtThr, float etaThr, float relIso, string WorkingPoint, string ProductionCampaign)
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedElectrons(15, 2.5, "Loose", "Spring15_50ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
      
      
      //if (selectedJets.size() >= 4)
      //  if (selectedJets[3]->Pt() < 30) selectedJets.clear();
      
      //vector<TRootMCParticle*> mcParticles;
      
      //if ( dataSetName.find("TT") == 0 )
      //{
      //  treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles, false);
      //  sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      //}
      
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
            if (vetoElectronsSemiMu.size() == 0) {
              //selecTableSemiMu.Fill(d,5,scaleFactor);
              selecTableSemiMu.Fill(d,4,scaleFactor);
              if ( selectedJets.size() >= 4 )
              {
                //selecTableSemiMu.Fill(d,6,scaleFactor);
                selecTableSemiMu.Fill(d,5,scaleFactor);
                eventSelected = true;
                
                for (unsigned int i = 0; i < selectedJets.size(); i++)
                {
                  if (selectedJets[i]->btag_combinedSecondaryVertexBJetTags() > 0.679) nb_bTags++;
                }
                		
                if ( nb_bTags >= 1 )
                {
                  //selecTableSemiMu.Fill(d,7,scaleFactor);
                  selecTableSemiMu.Fill(d,6,scaleFactor);
                  has1bjet = true;
                  
                  if ( nb_bTags >= 2 )
                  {
                    //selecTableSemiMu.Fill(d,8,scaleFactor);
                    selecTableSemiMu.Fill(d,7,scaleFactor);
                    has2bjets = true;
                  }
                }

              }  // at least 4 jets
            }  // no loose electrons
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
      
      if (verbose > 3)
        cout << endl << "  Event " << ievt << " is selected" << endl;
      if (verbose > 4)
        cout << "Event Id: " << event->eventId() << "  Run Id: " << event->runId() << "  Lumi block Id: " << event->lumiBlockId() << endl;
      
      
      
      ////////////////////////////
      ///  (GEN) JET MATCHING  ///
      ////////////////////////////
      
      vector<TLorentzVector> genJetsTLV, selectedJetsTLV;
      
      // Put the objects in TLorentzVectors, already ordened in decreasing Pt()
      for (unsigned int i = 0; i < selectedJets.size(); i++)
        selectedJetsTLV.push_back(*selectedJets[i]);
      for (unsigned int i = 0; i < genjets.size(); i++)
        genJetsTLV.push_back(*genjets[i]);
      
      JetPartonMatching matching = JetPartonMatching(genJetsTLV, selectedJetsTLV, 2, true, true, 0.3);		// genJets, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
      
      if (matching.getNumberOfAvailableCombinations() != 1)
        cerr << "matching.getNumberOfAvailableCombinations() = " << matching.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
      
      vector< pair<unsigned int, unsigned int> > GenJetPair; // First one is jet number, second one is genJet number
      
      for (unsigned int i = 0; i < genJetsTLV.size(); i++)
      {
        int matchedJetNumber = matching.getMatchForParton(i, 0);  // Get match for genJet
        if (matchedJetNumber > -1)
          GenJetPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
      }
      
      
      
      ////////////////////
      ///  FILL PLOTS  ///
      ////////////////////
      
      double HT = selectedJets[0]->Pt()+selectedJets[1]->Pt()+selectedJets[2]->Pt()+selectedJets[3]->Pt();
      TLorentzVector diffGenJet;
      double diffPT, diffEta, diffPhi, deltaR;
      
      histo1D["muon_pT"]->Fill(selectedMuons[0]->Pt());
      histo1D["muon_eta"]->Fill(selectedMuons[0]->Eta());
      histo1D["leadingJet_pT"]->Fill(selectedJets[0]->Pt());
      histo1D["Ht_4leadingJets"]->Fill(HT);
      for (unsigned int iGenJet = 0; iGenJet < genjets.size(); iGenJet++)
      {
        histo1D["genJets_pT"]->Fill(genjets[iGenJet]->Pt());
        histo1D["genJets_eta"]->Fill(genjets[iGenJet]->Eta());
        histo1D["genJets_nConstit"]->Fill(genjets[iGenJet]->nConstituents());
        histo1D["genJets_maxDistance"]->Fill(genjets[iGenJet]->maxDistance());
        histo1D["genJets_emEnergy"]->Fill(genjets[iGenJet]->emEnergy());
        histo1D["genJets_hadEnergy"]->Fill(genjets[iGenJet]->hadEnergy());
        histo1D["genJets_invEnergy"]->Fill(genjets[iGenJet]->invisibleEnergy());
        histo1D["genJets_n90"]->Fill(genjets[iGenJet]->n90());
        histo1D["genJets_n60"]->Fill(genjets[iGenJet]->n60());
      }

      for (unsigned int iPair = 0; iPair < GenJetPair.size(); iPair++)
      {
        diffGenJet = selectedJetsTLV[GenJetPair[iPair].first] - genJetsTLV[GenJetPair[iPair].second];
        deltaR = selectedJetsTLV[GenJetPair[iPair].first].DeltaR(genJetsTLV[GenJetPair[iPair].second]);
        
        histo1D["diff_jets_pT"]->Fill(diffGenJet.Pt());
        histo1D["diff_jets_eta"]->Fill(diffGenJet.Eta());
        histo1D["diff_jets_phi"]->Fill(diffGenJet.Phi());
        histo1D["diff_jets_deltaR"]->Fill(deltaR);
        
        diffPT = selectedJetsTLV[GenJetPair[iPair].first].Pt() - genJetsTLV[GenJetPair[iPair].second].Pt();
        diffEta = selectedJetsTLV[GenJetPair[iPair].first].Eta() - genJetsTLV[GenJetPair[iPair].second].Eta();
        diffPhi = selectedJetsTLV[GenJetPair[iPair].first].Phi() - genJetsTLV[GenJetPair[iPair].second].Phi();
        if ( diffPhi > TMath::Pi() )
          diffPhi -= TMath::Pi();
        else if ( diffPhi < -TMath::Pi() )
          diffPhi += TMath::Pi();
        
        histo1D["diff_jets_scal_pT"]->Fill(diffPT);
        histo1D["diff_jets_scal_eta"]->Fill(diffEta);
        histo1D["diff_jets_scal_phi"]->Fill(diffPhi);
      }
      
      
      
      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  /// Loop on events
    
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    cout << "Number of matched events: " << nofMatchedEvents << endl;
    
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  /// Loop on datasets
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  string pathPNG = "Plots_genJets/";
  mkdir(pathPNG.c_str(),0777);
  
  ///Write histograms
  fout->cd();
  // 1D
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
  string selectiontableMu = "SelectionTable_testGenJets_SemiMu.tex";
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
