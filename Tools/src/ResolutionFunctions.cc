//
//  ResolutionFunctions.cc
//  
//
//  Created by Lieselotte Moreels on 25/01/16.
//
//

#include "../interface/ResolutionFunctions.h"

//const std::string ResolutionFunctions::histoNames[] = {"Eparton_vs_Eparton-Ebjet", "Etparton_vs_Etparton-Etbjet", "Ptparton_vs_Ptparton-Ptbjet", "Eparton_vs_Thparton-Thbjet", "Eparton_vs_Etaparton-Etabjet", "Eparton_vs_Phiparton-Phibjet", "Eparton_vs_Eparton-Enonbjet", "Etparton_vs_Etparton-Etnonbjet", "Ptparton_vs_Ptparton-Ptnonbjet", "Eparton_vs_Thparton-Thnonbjet", "Eparton_vs_Etaparton-Etanonbjet", "Eparton_vs_Phiparton-Phinonbjet", "EgenMu_vs_EgenMu-ErecMu", "EtgenMu_vs_EtgenMu-EtrecMu", "PtgenMu_vs_PtgenMu-PtrecMu", "PtgenMu_vs_ThgenMu-ThrecMu", "PtgenMu_vs_EtagenMu-EtarecMu", "PtgenMu_vs_PhigenMu-PhirecMu", "EgenEl_vs_EgenEl-ErecEl", "EtgenEl_vs_EtgenEl-EtrecEl", "PtgenEl_vs_PtgenEl-PtrecEl", "genEl_vs_ThgenEl-ThrecEl", "EgenEl_vs_EtagenEl-EtarecEl", "EgenEl_vs_PhigenEl-PhirecEl", "Etparton_vs_Etparton-Etbjet_B", "Etparton_vs_Etparton-Etbjet_O", "Etparton_vs_Etparton-Etbjet_E", "Etparton_vs_Etparton-Etnonbjet_B", "Etparton_vs_Etparton-Etnonbjet_O", "Etparton_vs_Etparton-Etnonbjet_E"};

const std::string ResolutionFunctions::histoNames[] = {"Etbjet_vs_Etparton-Etbjet_B", "Etbjet_vs_Etparton-Etbjet_O", "Etbjet_vs_Etparton-Etbjet_E", "Etbjet_vs_Thparton-Thbjet_B", "Etbjet_vs_Thparton-Thbjet_O", "Etbjet_vs_Thparton-Thbjet_E", "Etbjet_vs_Phiparton-Phibjet_B", "Etbjet_vs_Phiparton-Phibjet_O", "Etbjet_vs_Phiparton-Phibjet_E", "Etnonbjet_vs_Etparton-Etnonbjet_B", "Etnonbjet_vs_Etparton-Etnonbjet_O", "Etnonbjet_vs_Etparton-Etnonbjet_E", "Etnonbjet_vs_Thparton-Thnonbjet_B", "Etnonbjet_vs_Thparton-Thnonbjet_O", "Etnonbjet_vs_Thparton-Thnonbjet_E", "Etnonbjet_vs_Phiparton-Phinonbjet_B", "Etnonbjet_vs_Phiparton-Phinonbjet_O", "Etnonbjet_vs_Phiparton-Phinonbjet_E", "EtrecMu_vs_EtgenMu-EtrecMu_B", "EtrecMu_vs_EtgenMu-EtrecMu_O", "EtrecMu_vs_EtgenMu-EtrecMu_E", "EtrecMu_vs_ThgenMu-ThrecMu_B", "EtrecMu_vs_ThgenMu-ThrecMu_O", "EtrecMu_vs_ThgenMu-ThrecMu_E", "EtrecMu_vs_PhigenMu-PhirecMu_B", "EtrecMu_vs_PhigenMu-PhirecMu_O", "EtrecMu_vs_PhigenMu-PhirecMu_E"};

//const std::string ResolutionFunctions::histoDescription[] = {"b jet energy", "b jet Et", "b jet pt", "b jet theta", "b jet eta", "b jet phi", "non-b jet energy", "non-b jet Et", "non-b jet pt", "non-b jet theta", "non-b jet eta", "non-b jet phi", "muon energy", "muon Et", "muon pt", "muon theta", "muon eta", "muon phi", "electron energy", "electron Et", "electron pt", "electron theta", "electron eta", "electron phi", "b jet Et barrel", "b jet Et overlap area", "b jet Et endcap", "non-b jet Et barrel", "non-b jet Et overlap area", "non-b jet Et endcap"};

const std::string ResolutionFunctions::histoDescription[] = {"b jet Et barrel", "b jet Et overlap area", "b jet Et endcap", "b jet theta barrel", "b jet theta overlap area", "b jet theta endcap", "b jet phi barrel", "b jet phi overlap area", "b jet phi endcap", "non-b jet Et barrel", "non-b jet Et overlap area", "non-b jet Et endcap", "non-b jet theta barrel", "non-b jet theta overlap area", "non-b jet theta endcap", "non-b jet phi barrel", "non-b jet phi overlap area", "non-b jet phi endcap", "muon Et barrel", "muon Et overlap", "muon Et endcap", "muon theta barrel", "muon theta overlap", "muon theta endcap", "muon phi barrel", "muon phi overlap", "muon phi endcap"};

Double_t ResolutionFunctions::sGaus(Double_t *x, Double_t *par)
{
  Double_t narrowGaus = TMath::Exp( - TMath::Power( (x[0] - par[0])/par[1] , 2) /2. )/(par[1]*sqrt(2*TMath::Pi()));
  return par[2] * narrowGaus + par[3];
}

Double_t ResolutionFunctions::dblGaus(Double_t *x, Double_t *par)
{
  Double_t norm = 1./TMath::Sqrt(2.*TMath::Pi()) * par[5]/TMath::Sqrt( pow(par[1], 2) + par[2]*pow(par[4], 2) );
  Double_t narrowGaus = TMath::Exp( - TMath::Power( (x[0] - par[0])/par[1] , 2) /2. );
  Double_t broadGaus = TMath::Exp( - TMath::Power( (x[0] - par[3])/par[4] , 2) /2. );
  return norm * ( narrowGaus + par[2] * broadGaus );
}

Double_t ResolutionFunctions::dblGausParFill(Double_t *x, Double_t *par)
{
  Double_t par0 = par[0] + x[0] * par[1];
  Double_t par1 = par[2] + x[0] * par[3];
  Double_t par2 = par[4] + x[0] * par[5];
  Double_t par3 = par[6] + x[0] * par[7];
  Double_t par4 = par[8] + x[0] * par[9];
  Double_t par5 = par[10] + x[0] * par[11];
  
  Double_t norm = 1./TMath::Sqrt(2.*TMath::Pi()) * par5/TMath::Sqrt( pow(par1, 2) + par2*pow(par4, 2) );
  Double_t narrowGaus = TMath::Exp( - TMath::Power( (x[1] - par0)/par1 , 2) /2. );
  Double_t broadGaus = TMath::Exp( - TMath::Power( (x[1] - par3)/par4 , 2) /2. );
  return norm * ( narrowGaus + par2 * broadGaus );
}

ResolutionFunctions::ResolutionFunctions(bool calculateResolutionFunctions, bool _verbose):
verbose(true), muon(false), electron(false), getHistos(false), histoRes2D(), inputFileName("PlotsForResolutionFunctions.root"), nHistos(sizeof(histoNames)/sizeof(histoNames[0])), nJets(4)
{
  verbose = _verbose;
  if (verbose) std::cout << "ResolutionFunctions::ResolutionFunctions - Initialising..." << std::endl;
  
  if (calculateResolutionFunctions)
  {
    this->bookHistograms();
  }
  else if (verbose)
  {
    std::cout << "Using resolution functions from file..." << std::endl;
  }
}

ResolutionFunctions::~ResolutionFunctions()
{
  /// Clean output files, etc.
}

std::string ResolutionFunctions::toStr(int number)
{
  std::ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << number;
  return convert.str();
}

void ResolutionFunctions::bookHistograms()
{
  if (verbose) std::cout << "ResolutionFunctions::bookHistograms - Initialising..." << std::endl;
  
  /// Energy
  histoRes2D["Eparton_vs_Enonbjet"] = new TH2F("Eparton_vs_Enonbjet","Eparton_vs_Enonbjet", 40, 0, 400, 40, 0, 400);
  histoRes2D["Enonbjet_vs_Eparton-Enonbjet"] = new TH2F("Enonbjet_vs_Eparton-Enonbjet","Enonbjet_vs_Eparton-Enonbjet", 6, 0, 300, 100, -100, 100);
  histoRes2D["Eparton_vs_Ebjet"] = new TH2F("Eparton_vs_Ebjet","Eparton_vs_Ebjet", 40, 0, 400, 40, 0, 400);
  histoRes2D["Ebjet_vs_Eparton-Ebjet"] = new TH2F("Ebjet_vs_Eparton-Ebjet","Ebjet_vs_Eparton-Ebjet", 6, 0, 300, 100, -100, 100);
  histoRes2D["EgenEl_vs_ErecEl"] = new TH2F("EgenEl_vs_ErecEl","EgenEl_vs_ErecEl", 40, 0, 400, 40, 0, 400);
  histoRes2D["ErecEl_vs_EgenEl-ErecEl"] = new TH2F("ErecEl_vs_EgenEl-ErecEl","ErecEl_vs_EgenEl-ErecEl", 5, 0, 400, 100, -50, 50);
  histoRes2D["EgenMu_vs_ErecMu"] = new TH2F("EgenMu_vs_ErecMu","EgenMu_vs_ErecMu", 40, 0, 400, 40, 0, 400);
  histoRes2D["ErecMu_vs_EgenMu-ErecMu"] = new TH2F("ErecMu_vs_EgenMu-ErecMu","ErecMu_vs_EgenMu-ErecMu", 5, 0, 400, 100, -4, 4); //-0.003, 0.003);
  
  /// Et
  histoRes2D["Etparton_vs_Etnonbjet"] = new TH2F("Etparton_vs_Etnonbjet","Etparton_vs_Etnonbjet", 30, 0, 300, 30, 0, 300);
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet"] = new TH2F("Etnonbjet_vs_Etparton-Etnonbjet","Etnonbjet_vs_Etparton-Etnonbjet", 6, 0, 250, 100, -100, 100);
  histoRes2D["Etparton_vs_Etbjet"] = new TH2F("Etparton_vs_Etbjet","Etparton_vs_Etbjet", 30, 0, 300, 30, 0, 300);
  histoRes2D["Etbjet_vs_Etparton-Etbjet"] = new TH2F("Etbjet_vs_Etparton-Etbjet","Etbjet_vs_Etparton-Etbjet", 6, 0, 250, 100, -100, 100);
  histoRes2D["EtgenEl_vs_EtrecEl"] = new TH2F("EtgenEl_vs_EtrecEl","EtgenEl_vs_EtrecEl", 30, 0, 300, 30, 0, 300);
  histoRes2D["EtrecEl_vs_EtgenEl-EtrecEl"] = new TH2F("EtrecEl_vs_EtgenEl-EtrecEl","EtrecEl_vs_EtgenEl-EtrecEl", 6, 0, 250, 100, -50, 50);
  histoRes2D["EtgenMu_vs_EtrecMu"] = new TH2F("EtgenMu_vs_EtrecMu","EtgenMu_vs_EtrecMu", 30, 0, 300, 30, 0, 300);
  histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu"] = new TH2F("EtrecMu_vs_EtgenMu-EtrecMu","EtrecMu_vs_EtgenMu-EtrecMu", 6, 0, 250, 100, -3, 3);
  
  /// Et binned in eta
  histoRes2D["Etparton_vs_Etnonbjet_B"] = new TH2F("Etparton_vs_Etnonbjet_B","Etparton_vs_Etnonbjet_B", 30, 0, 300, 30, 0, 300);
  histoRes2D["Etparton_vs_Etnonbjet_O"] = new TH2F("Etparton_vs_Etnonbjet_O","Etparton_vs_Etnonbjet_O", 30, 0, 300, 30, 0, 300);
  histoRes2D["Etparton_vs_Etnonbjet_E"] = new TH2F("Etparton_vs_Etnonbjet_E","Etparton_vs_Etnonbjet_E", 30, 0, 300, 30, 0, 300);
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_B"] = new TH2F("Etnonbjet_vs_Etparton-Etnonbjet_B","; E_{T} (light jet); E_{T} (parton) - E_{T} (light jet)", 6, 0, 250, 100, -100, 100);
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_O"] = new TH2F("Etnonbjet_vs_Etparton-Etnonbjet_O","; E_{T} (light jet); E_{T} (parton) - E_{T} (light jet)", 6, 0, 250, 100, -100, 100);
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_E"] = new TH2F("Etnonbjet_vs_Etparton-Etnonbjet_E","; E_{T} (light jet); E_{T} (parton) - E_{T} (light jet)", 6, 0, 250, 100, -100, 100);
  histoRes2D["Etparton_vs_Etbjet_B"] = new TH2F("Etparton_vs_Etbjet_B","Etparton_vs_Etbjet_B", 30, 0, 300, 30, 0, 300);
  histoRes2D["Etparton_vs_Etbjet_O"] = new TH2F("Etparton_vs_Etbjet_O","Etparton_vs_Etbjet_O", 30, 0, 300, 30, 0, 300);
  histoRes2D["Etparton_vs_Etbjet_E"] = new TH2F("Etparton_vs_Etbjet_E","Etparton_vs_Etbjet_E", 30, 0, 300, 30, 0, 300);
  histoRes2D["Etbjet_vs_Etparton-Etbjet_B"] = new TH2F("Etbjet_vs_Etparton-Etbjet_B","; E_{T} (b jet); E_{T} (parton) - E_{T} (b jet)", 6, 0, 250, 100, -100, 100);
  histoRes2D["Etbjet_vs_Etparton-Etbjet_O"] = new TH2F("Etbjet_vs_Etparton-Etbjet_O","; E_{T} (b jet); E_{T} (parton) - E_{T} (b jet)", 6, 0, 250, 100, -100, 100);
  histoRes2D["Etbjet_vs_Etparton-Etbjet_E"] = new TH2F("Etbjet_vs_Etparton-Etbjet_E","; E_{T} (b jet); E_{T} (parton) - E_{T} (b jet)", 6, 0, 250, 100, -100, 100);
  histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu_B"] = new TH2F("EtrecMu_vs_EtgenMu-EtrecMu_B","; E_{T} (#mu_{reco}); E_{T} (#mu_{gen}) - E_{T} (#mu_{reco})", 4, 0, 200, 60, -3, 3);
  histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu_O"] = new TH2F("EtrecMu_vs_EtgenMu-EtrecMu_O","; E_{T} (#mu_{reco}); E_{T} (#mu_{gen}) - E_{T} (#mu_{reco})", 4, 0, 200, 60, -3, 3);
  histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu_E"] = new TH2F("EtrecMu_vs_EtgenMu-EtrecMu_E","; E_{T} (#mu_{reco}); E_{T} (#mu_{gen}) - E_{T} (#mu_{reco})", 4, 0, 200, 60, -3, 3);
  
  /// pT
  histoRes2D["Ptparton_vs_Ptnonbjet"] = new TH2F("Ptparton_vs_Ptnonbjet","Ptparton_vs_Ptnonbjet", 30, 0, 300, 30, 0, 300);
  histoRes2D["Ptnonbjet_vs_Ptparton-Ptnonbjet"] = new TH2F("Ptnonbjet_vs_Ptparton-Ptnonbjet","Ptnonbjet_vs_Ptparton-Ptnonbjet", 6, 0, 250, 100, -100, 100);
  histoRes2D["Ptparton_vs_Ptbjet"] = new TH2F("Ptparton_vs_Ptbjet","Ptparton_vs_Ptbjet", 30, 0, 300, 30, 0, 300);
  histoRes2D["Ptbjet_vs_Ptparton-Ptbjet"] = new TH2F("Ptbjet_vs_Ptparton-Ptbjet","Ptbjet_vs_Ptparton-Ptbjet", 6, 0, 250, 100, -100, 100);
  histoRes2D["PtgenEl_vs_PtrecEl"] = new TH2F("PtgenEl_vs_PtrecEl","PtgenEl_vs_PtrecEl", 30, 0, 300, 30, 0, 300);
  histoRes2D["PtrecEl_vs_PtgenEl-PtrecEl"] = new TH2F("PtrecEl_vs_PtgenEl-PtrecEl","PtrecEl_vs_PtgenEl-PtrecEl", 5, 0, 250, 100, -50, 50);
  histoRes2D["PtgenMu_vs_PtrecMu"] = new TH2F("PtgenMu_vs_PtrecMu","PtgenMu_vs_PtrecMu", 30, 0, 300, 30, 0, 300);
  histoRes2D["PtrecMu_vs_PtgenMu-PtrecMu"] = new TH2F("PtrecMu_vs_PtgenMu-PtrecMu","PtrecMu_vs_PtgenMu-PtrecMu", 5, 0, 250, 100, -3, 3);
  
  /// Theta
  histoRes2D["Thparton_vs_Thnonbjet"] = new TH2F("Thparton_vs_Thnonbjet","Thparton_vs_Thnonbjet", 60, 0, 3.15, 60, 0, 3.15);
  histoRes2D["Thnonbjet_vs_Thparton-Thnonbjet"] = new TH2F("Thnonbjet_vs_Thparton-Thnonbjet","Thnonbjet_vs_Thparton-Thnonbjet", 5, 0.1, 3.0, 60, -0.15, 0.15);
  histoRes2D["Thparton_vs_Thbjet"] = new TH2F("Thparton_vs_Thbjet","Thparton_vs_Thbjet", 60, 0, 3.15, 60, 0, 3.15);
  histoRes2D["Thbjet_vs_Thparton-Thbjet"] = new TH2F("Thbjet_vs_Thparton-Thbjet","Thbjet_vs_Thparton-Thbjet", 5, 0.1, 3.0, 60, -0.15, 0.15);
  histoRes2D["ThgenEl_vs_ThrecEl"] = new TH2F("ThgenEl_vs_ThrecEl","ThgenEl_vs_ThrecEl", 60, 0, 3.15, 60, 0, 3.15);
  histoRes2D["ThrecEl_vs_ThgenEl-ThrecEl"] = new TH2F("ThrecEl_vs_ThgenEl-ThrecEl","ThrecEl_vs_ThgenEl-ThrecEl", 5, 0.2, 3.0, 60, -0.05, 0.05);
  histoRes2D["ThgenMu_vs_ThrecMu"] = new TH2F("ThgenMu_vs_ThrecMu","ThgenMu_vs_ThrecMu", 60, 0, 3.15, 60, 0, 3.15);
  histoRes2D["ThrecMu_vs_ThgenMu-ThrecMu"] = new TH2F("ThrecMu_vs_ThgenMu-ThrecMu","ThrecMu_vs_ThgenMu-ThrecMu", 5, 0.2, 3.0, 60, -0.05, 0.05);
  
  /// Energy vs theta
  histoRes2D["Eparton_vs_Thnonbjet"] = new TH2F("Eparton_vs_Thnonbjet","Eparton_vs_Thnonbjet", 40, 0, 400, 60, 0, 3.15);
  histoRes2D["Enonbjet_vs_Thparton-Thnonbjet"] = new TH2F("Enonbjet_vs_Thparton-Thnonbjet","Enonbjet_vs_Thparton-Thnonbjet", 6, 0, 250, 100, -0.15, 0.15);
  histoRes2D["Eparton_vs_Thbjet"] = new TH2F("Eparton_vs_Thbjet","Eparton_vs_Thbjet", 40, 0, 400, 60, 0, 3.15);
  histoRes2D["Ebjet_vs_Thparton-Thbjet"] = new TH2F("Ebjet_vs_Thparton-Thbjet","Ebjet_vs_Thparton-Thbjet", 6, 0, 300, 100, -0.15, 0.15);
  histoRes2D["EgenEl_vs_ThrecEl"] = new TH2F("EgenEl_vs_ThrecEl","EgenEl_vs_ThrecEl", 40, 0, 400, 60, 0, 3.15);
  histoRes2D["ErecEl_vs_ThgenEl-ThrecEl"] = new TH2F("ErecEl_vs_ThgenEl-ThrecEl","ErecEl_vs_ThgenEl-ThrecEl", 5, 0, 400, 100, -0.05, 0.05);
  histoRes2D["PtgenMu_vs_ThrecMu"] = new TH2F("PtgenMu_vs_ThrecMu","PtgenMu_vs_ThrecMu", 30, 0, 300, 60, 0, 3.15);
  histoRes2D["PtrecMu_vs_ThgenMu-ThrecMu"] = new TH2F("PtrecMu_vs_ThgenMu-ThrecMu","PtrecMu_vs_ThgenMu-ThrecMu", 5, 0, 250, 100, -0.002, 0.002);
  
  /// Et vs theta, binned in eta
  histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_B"] = new TH2F("Etnonbjet_vs_Thparton-Thnonbjet_B","; E_{T} (light jet); #theta(parton) - #theta(light jet)", 5, 0, 250, 120, -0.1, 0.1);
  histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_O"] = new TH2F("Etnonbjet_vs_Thparton-Thnonbjet_O","; E_{T} (light jet); #theta(parton) - #theta(light jet)", 5, 0, 250, 100, -0.1, 0.1);
  histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_E"] = new TH2F("Etnonbjet_vs_Thparton-Thnonbjet_E","; E_{T} (light jet); #theta(parton) - #theta(light jet)", 5, 0, 250, 100, -0.1, 0.1);
  histoRes2D["Etbjet_vs_Thparton-Thbjet_B"] = new TH2F("Etbjet_vs_Thparton-Thbjet_B","; E_{T} (b jet); #theta(parton) - #theta(b jet)", 5, 0, 250, 120, -0.1, 0.1);
  histoRes2D["Etbjet_vs_Thparton-Thbjet_O"] = new TH2F("Etbjet_vs_Thparton-Thbjet_O","; E_{T} (b jet); #theta(parton) - #theta(b jet)", 5, 0, 250, 120, -0.1, 0.1);
  histoRes2D["Etbjet_vs_Thparton-Thbjet_E"] = new TH2F("Etbjet_vs_Thparton-Thbjet_E","; E_{T} (b jet); #theta(parton) - #theta(b jet)", 5, 0, 250, 120, -0.1, 0.1);
  histoRes2D["EtrecMu_vs_ThgenMu-ThrecMu_B"] = new TH2F("EtrecMu_vs_ThgenMu-ThrecMu_B","; E_{T} (#mu_{reco}); #theta(#mu_{gen}) - #theta(#mu_{reco})", 4, 0, 200, 80, -0.002, 0.002);
  histoRes2D["EtrecMu_vs_ThgenMu-ThrecMu_O"] = new TH2F("EtrecMu_vs_ThgenMu-ThrecMu_O","; E_{T} (#mu_{reco}); #theta(#mu_{gen}) - #theta(#mu_{reco})", 4, 0, 200, 80, -0.002, 0.002);
  histoRes2D["EtrecMu_vs_ThgenMu-ThrecMu_E"] = new TH2F("EtrecMu_vs_ThgenMu-ThrecMu_E","; E_{T} (#mu_{reco}); #theta(#mu_{gen}) - #theta(#mu_{reco})", 4, 0, 200, 80, -0.002, 0.002);
  
  /// Eta
  histoRes2D["Etaparton_vs_Etanonbjet"] = new TH2F("Etaparton_vs_Etanonbjet","Etaparton_vs_Etanonbjet", 60, -2.5, 2.5, 60, 0, 2.5);
  histoRes2D["Etanonbjet_vs_Etaparton-Etanonbjet"] = new TH2F("Etanonbjet_vs_Etaparton-Etanonbjet","Etanonbjet_vs_Etaparton-Etanonbjet", 5, 0.1, 2.4, 60, -0.15, 0.15);
  histoRes2D["Etaparton_vs_Etabjet"] = new TH2F("Etaparton_vs_Etabjet","Etaparton_vs_Etabjet", 60, -2.5, 2.5, 60, 0, 2.5);
  histoRes2D["Etabjet_vs_Etaparton-Etabjet"] = new TH2F("Etabjet_vs_Etaparton-Etabjet","Etabjet_vs_Etaparton-Etabjet", 5, 0.1, 2.4, 60, -0.15, 0.15);
  histoRes2D["EtagenEl_vs_EtarecEl"] = new TH2F("EtagenEl_vs_EtarecEl","EtagenEl_vs_EtarecEl", 60, -2.5, 2.5, 60, 0, 2.5);
  histoRes2D["EtarecEl_vs_EtagenEl-EtarecEl"] = new TH2F("EtarecEl_vs_EtagenEl-EtarecEl","EtarecEl_vs_EtagenEl-EtarecEl", 5, 0.1, 2.1, 60, -0.05, 0.05);
  histoRes2D["EtagenMu_vs_EtarecMu"] = new TH2F("EtagenMu_vs_EtarecMu","EtagenMu_vs_EtarecMu", 60, -2.5, 2.5, 60, 0, 2.5);
  histoRes2D["EtarecMu_vs_EtagenMu-EtarecMu"] = new TH2F("EtarecMu_vs_EtagenMu-EtarecMu","EtarecMu_vs_EtagenMu-EtarecMu", 5, 0.1, 2.1, 60, -0.05, 0.05);
  
  /// Energy, pT vs eta
  histoRes2D["Eparton_vs_Etanonbjet"] = new TH2F("Eparton_vs_Etanonbjet","Eparton_vs_Etanonbjet", 40, 0, 400, 60, -2.5, 2.5);
  histoRes2D["Enonbjet_vs_Etaparton-Etanonbjet"] = new TH2F("Enonbjet_vs_Etaparton-Etanonbjet","Enonbjet_vs_Etaparton-Etanonbjet", 6, 0, 250, 100, -0.15, 0.15);
  histoRes2D["Eparton_vs_Etabjet"] = new TH2F("Eparton_vs_Etabjet","Eparton_vs_Etabjet", 40, 0, 400, 60, -2.5, 2.5);
  histoRes2D["Ebjet_vs_Etaparton-Etabjet"] = new TH2F("Ebjet_vs_Etaparton-Etabjet","Ebjet_vs_Etaparton-Etabjet", 6, 0, 300, 100, -0.15, 0.15);
  histoRes2D["EgenEl_vs_EtarecEl"] = new TH2F("EgenEl_vs_EtarecEl","EgenEl_vs_EtarecEl", 40, 0, 400, 60, -2.5, 2.5);
  histoRes2D["ErecEl_vs_EtagenEl-EtarecEl"] = new TH2F("ErecEl_vs_EtagenEl-EtarecEl","ErecEl_vs_EtagenEl-EtarecEl", 5, 0, 400, 100, -0.05, 0.05);
  histoRes2D["PtgenMu_vs_EtarecMu"] = new TH2F("PtgenMu_vs_EtarecMu","PtgenMu_vs_EtarecMu", 30, 0, 300, 60, -2.5, 2.5);
  histoRes2D["PtrecMu_vs_EtagenMu-EtarecMu"] = new TH2F("PtrecMu_vs_EtagenMu-EtarecMu","PtrecMu_vs_EtagenMu-EtarecMu", 5, 0, 250, 100, -0.05, 0.05);
  
  /// Phi
  histoRes2D["Phiparton_vs_Phinonbjet"] = new TH2F("Phiparton_vs_Phinonbjet","Phiparton_vs_Phinonbjet", 120, -3.2, 3.2, 120, -3.2, 3.2);
  histoRes2D["Phinonbjet_vs_Phiparton-Phinonbjet"] = new TH2F("Phinonbjet_vs_Phiparton-Phinonbjet","Phinonbjet_vs_Phiparton-Phinonbjet", 5, -3.2, 3.2, 120, -0.3, 0.3);
  histoRes2D["Phiparton_vs_Phibjet"] = new TH2F("Phiparton_vs_Phibjet","Phiparton_vs_Phibjet", 120, -3.2, 3.2, 120, -3.2, 3.2);
  histoRes2D["Phibjet_vs_Phiparton-Phibjet"] = new TH2F("Phibjet_vs_Phiparton-Phibjet","Phibjet_vs_Phiparton-Phibjet", 5, -3.2, 3.2, 120, -0.3, 0.3);
  histoRes2D["PhigenEl_vs_PhirecEl"] = new TH2F("PhigenEl_vs_PhirecEl","PhigenEl_vs_PhirecEl", 120, -3.2, 3.2, 120, -3.2, 3.2);
  histoRes2D["PhirecEl_vs_PhigenEl-PhirecEl"] = new TH2F("PhirecEl_vs_PhigenEl-PhirecEl","PhirecEl_vs_PhigenEl-PhirecEl", 5, -3.2, 3.2, 120, -0.05, 0.05);
  histoRes2D["PhigenMu_vs_PhirecMu"] = new TH2F("PhigenMu_vs_PhirecMu","PhigenMu_vs_PhirecMu", 120, -3.2, 3.2, 120, -3.2, 3.2);
  histoRes2D["PhirecMu_vs_PhigenMu-PhirecMu"] = new TH2F("PhirecMu_vs_PhigenMu-PhirecMu","PhirecMu_vs_PhigenMu-PhirecMu", 5, -3.2, 3.2, 120, -0.05, 0.05);
  
  /// Energy, pT vs phi
  histoRes2D["Eparton_vs_Phinonbjet"] = new TH2F("Eparton_vs_Phinonbjet","Eparton_vs_Phinonbjet", 40, 0, 400, 120, -3.2, 3.2);
  histoRes2D["Enonbjet_vs_Phiparton-Phinonbjet"] = new TH2F("Enonbjet_vs_Phiparton-Phinonbjet","Enonbjet_vs_Phiparton-Phinonbjet", 6, 0, 250, 100, -0.3, 0.3);
  histoRes2D["Eparton_vs_Phibjet"] = new TH2F("Eparton_vs_Phibjet","Eparton_vs_Phibjet", 40, 0, 400, 120, -3.2, 3.2);
  histoRes2D["Ebjet_vs_Phiparton-Phibjet"] = new TH2F("Ebjet_vs_Phiparton-Phibjet","Ebjet_vs_Phiparton-Phibjet", 6, 0, 300, 100, -0.2, 0.2);
  histoRes2D["EgenEl_vs_PhirecEl"] = new TH2F("EgenEl_vs_PhirecEl","EgenEl_vs_PhirecEl", 40, 0, 400, 120, -3.2, 3.2);
  histoRes2D["ErecEl_vs_PhigenEl-PhirecEl"] = new TH2F("ErecEl_vs_PhigenEl-PhirecEl","ErecEl_vs_PhigenEl-PhirecEl", 5, 0, 400, 100, -0.05, 0.05);
  histoRes2D["PtgenMu_vs_PhirecMu"] = new TH2F("PtgenMu_vs_PhirecMu","PtgenMu_vs_PhirecMu", 30, 0, 300, 120, -3.2, 3.2);
  histoRes2D["PtrecMu_vs_PhigenMu-PhirecMu"] = new TH2F("PtrecMu_vs_PhigenMu-PhirecMu","PtrecMu_vs_PhigenMu-PhirecMu", 5, 0, 250, 100, -0.002, 0.002);
  
  /// Et vs phi, binned in eta
  histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_B"] = new TH2F("Etnonbjet_vs_Phiparton-Phinonbjet_B","; E_{T} (light jet); #phi(parton) - #phi(light jet)", 5, 0, 250, 120, -0.18, 0.18);
  histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_O"] = new TH2F("Etnonbjet_vs_Phiparton-Phinonbjet_O","; E_{T} (light jet); #phi(parton) - #phi(light jet)", 5, 0, 250, 100, -0.18, 0.18);
  histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_E"] = new TH2F("Etnonbjet_vs_Phiparton-Phinonbjet_E","; E_{T} (light jet); #phi(parton) - #phi(light jet)", 5, 0, 250, 100, -0.18, 0.18);
  histoRes2D["Etbjet_vs_Phiparton-Phibjet_B"] = new TH2F("Etbjet_vs_Phiparton-Phibjet_B","; E_{T} (b jet); #phi(parton) - #phi(b jet)", 5, 0, 250, 120, -0.18, 0.18);
  histoRes2D["Etbjet_vs_Phiparton-Phibjet_O"] = new TH2F("Etbjet_vs_Phiparton-Phibjet_O","; E_{T} (b jet); #phi(parton) - #phi(b jet)", 5, 0, 250, 100, -0.18, 0.18);
  histoRes2D["Etbjet_vs_Phiparton-Phibjet_E"] = new TH2F("Etbjet_vs_Phiparton-Phibjet_E","; E_{T} (b jet); #phi(parton) - #phi(b jet)", 5, 0, 250, 100, -0.18, 0.18);
  histoRes2D["EtrecMu_vs_PhigenMu-PhirecMu_B"] = new TH2F("EtrecMu_vs_PhigenMu-PhirecMu_B","; E_{T} (#mu_{reco}); #phi(#mu_{gen}) - #phi(#mu_{reco})", 4, 0, 200, 80, -0.002, 0.002);
  histoRes2D["EtrecMu_vs_PhigenMu-PhirecMu_O"] = new TH2F("EtrecMu_vs_PhigenMu-PhirecMu_O","; E_{T} (#mu_{reco}); #phi(#mu_{gen}) - #phi(#mu_{reco})", 4, 0, 200, 80, -0.002, 0.002);
  histoRes2D["EtrecMu_vs_PhigenMu-PhirecMu_E"] = new TH2F("EtrecMu_vs_PhigenMu-PhirecMu_E","; E_{T} (#mu_{reco}); #phi(#mu_{gen}) - #phi(#mu_{reco})", 4, 0, 200, 80, -0.002, 0.002);
  
  if (verbose) std::cout << "                                    - Histograms booked" << std::endl;
}

void ResolutionFunctions::fillJets(std::vector<TLorentzVector> &parton, std::vector<TLorentzVector> &jet)
{
  /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
  nJets = parton.size();
  if ( jet.size() < nJets ) nJets = jet.size();
  
  /// Energy
  histoRes2D["Eparton_vs_Enonbjet"]->Fill(parton[0].E(),jet[0].E());
  histoRes2D["Enonbjet_vs_Eparton-Enonbjet"]->Fill(jet[0].E(),parton[0].E()-jet[0].E());
  histoRes2D["Eparton_vs_Enonbjet"]->Fill(parton[1].E(),jet[1].E());
  histoRes2D["Enonbjet_vs_Eparton-Enonbjet"]->Fill(jet[1].E(),parton[1].E()-jet[1].E());
  histoRes2D["Eparton_vs_Ebjet"]->Fill(parton[2].E(),jet[2].E());
  histoRes2D["Ebjet_vs_Eparton-Ebjet"]->Fill(jet[2].E(),parton[2].E()-jet[2].E());
  if ( nJets > 3 )
  {
    histoRes2D["Eparton_vs_Ebjet"]->Fill(parton[3].E(),jet[3].E());
    histoRes2D["Ebjet_vs_Eparton-Ebjet"]->Fill(jet[3].E(),parton[3].E()-jet[3].E());
  }
  
  /// Et
  histoRes2D["Etparton_vs_Etnonbjet"]->Fill(parton[0].Et(),jet[0].Et());
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet"]->Fill(jet[0].Et(),parton[0].Et()-jet[0].Et());
  histoRes2D["Etparton_vs_Etnonbjet"]->Fill(parton[1].Et(),jet[1].Et());
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet"]->Fill(jet[1].Et(),parton[1].Et()-jet[1].Et());
  histoRes2D["Etparton_vs_Etbjet"]->Fill(parton[2].Et(),jet[2].Et());
  histoRes2D["Etbjet_vs_Etparton-Etbjet"]->Fill(jet[2].Et(),parton[2].Et()-jet[2].Et());
  if ( nJets > 3 )
  {
    histoRes2D["Etparton_vs_Etbjet"]->Fill(parton[3].Et(),jet[3].Et());
    histoRes2D["Etbjet_vs_Etparton-Etbjet"]->Fill(jet[3].Et(),parton[3].Et()-jet[3].Et());
  }
  
  /// pT
  histoRes2D["Ptparton_vs_Ptnonbjet"]->Fill(parton[0].Pt(),jet[0].Pt());
  histoRes2D["Ptnonbjet_vs_Ptparton-Ptnonbjet"]->Fill(jet[0].Pt(),parton[0].Pt()-jet[0].Pt());
  histoRes2D["Ptparton_vs_Ptnonbjet"]->Fill(parton[1].Pt(),jet[1].Pt());
  histoRes2D["Ptnonbjet_vs_Ptparton-Ptnonbjet"]->Fill(jet[1].Pt(),parton[1].Pt()-jet[1].Pt());
  histoRes2D["Ptparton_vs_Ptbjet"]->Fill(parton[2].Pt(),jet[2].Pt());
  histoRes2D["Ptbjet_vs_Ptparton-Ptbjet"]->Fill(jet[2].Pt(),parton[2].Pt()-jet[2].Pt());
  if ( nJets > 3 )
  {
    histoRes2D["Ptparton_vs_Ptbjet"]->Fill(parton[3].Pt(),jet[3].Pt());
    histoRes2D["Ptbjet_vs_Ptparton-Ptbjet"]->Fill(jet[3].Pt(),parton[3].Pt()-jet[3].Pt());
  }
  
  /// Theta
  histoRes2D["Thparton_vs_Thnonbjet"]->Fill(parton[0].Theta(),jet[0].Theta());
  histoRes2D["Thnonbjet_vs_Thparton-Thnonbjet"]->Fill(jet[0].Theta(),parton[0].Theta()-jet[0].Theta());
  histoRes2D["Thparton_vs_Thnonbjet"]->Fill(parton[1].Theta(),jet[1].Theta());
  histoRes2D["Thnonbjet_vs_Thparton-Thnonbjet"]->Fill(jet[1].Theta(),parton[1].Theta()-jet[1].Theta());
  histoRes2D["Thparton_vs_Thbjet"]->Fill(parton[2].Theta(),jet[2].Theta());
  histoRes2D["Thbjet_vs_Thparton-Thbjet"]->Fill(jet[2].Theta(),parton[2].Theta()-jet[2].Theta());
  if ( nJets > 3 )
  {
    histoRes2D["Thparton_vs_Thbjet"]->Fill(parton[3].Theta(),jet[3].Theta());
    histoRes2D["Thbjet_vs_Thparton-Thbjet"]->Fill(jet[3].Theta(),parton[3].Theta()-jet[3].Theta());
  }
  histoRes2D["Eparton_vs_Thnonbjet"]->Fill(parton[0].E(),jet[0].Theta());
  histoRes2D["Enonbjet_vs_Thparton-Thnonbjet"]->Fill(jet[0].E(),parton[0].Theta()-jet[0].Theta());
  histoRes2D["Eparton_vs_Thnonbjet"]->Fill(parton[1].E(),jet[1].Theta());
  histoRes2D["Enonbjet_vs_Thparton-Thnonbjet"]->Fill(jet[1].E(),parton[1].Theta()-jet[1].Theta());
  histoRes2D["Eparton_vs_Thbjet"]->Fill(parton[2].E(),jet[2].Theta());
  histoRes2D["Ebjet_vs_Thparton-Thbjet"]->Fill(jet[2].E(),parton[2].Theta()-jet[2].Theta());
  if ( nJets > 3 )
  {
    histoRes2D["Eparton_vs_Thbjet"]->Fill(parton[3].E(),jet[3].Theta());
    histoRes2D["Ebjet_vs_Thparton-Thbjet"]->Fill(jet[3].E(),parton[3].Theta()-jet[3].Theta());
  }
  
  /// Eta
  histoRes2D["Etaparton_vs_Etanonbjet"]->Fill(parton[0].Eta(),jet[0].Eta());
  histoRes2D["Etanonbjet_vs_Etaparton-Etanonbjet"]->Fill(fabs(jet[0].Eta()),parton[0].Eta()-jet[0].Eta());
  histoRes2D["Etaparton_vs_Etanonbjet"]->Fill(parton[1].Eta(),jet[1].Eta());
  histoRes2D["Etanonbjet_vs_Etaparton-Etanonbjet"]->Fill(fabs(jet[1].Eta()),parton[1].Eta()-jet[1].Eta());
  histoRes2D["Etaparton_vs_Etabjet"]->Fill(parton[2].Eta(),jet[2].Eta());
  histoRes2D["Etabjet_vs_Etaparton-Etabjet"]->Fill(fabs(jet[2].Eta()),parton[2].Eta()-jet[2].Eta());
  if ( nJets > 3 )
  {
    histoRes2D["Etaparton_vs_Etabjet"]->Fill(parton[3].Eta(),jet[3].Eta());
    histoRes2D["Etabjet_vs_Etaparton-Etabjet"]->Fill(fabs(jet[3].Eta()),parton[3].Eta()-jet[3].Eta());
  }
  histoRes2D["Eparton_vs_Etanonbjet"]->Fill(parton[0].E(),jet[0].Eta());
  histoRes2D["Enonbjet_vs_Etaparton-Etanonbjet"]->Fill(jet[0].E(),parton[0].Eta()-jet[0].Eta());
  histoRes2D["Eparton_vs_Etanonbjet"]->Fill(parton[1].E(),jet[1].Eta());
  histoRes2D["Enonbjet_vs_Etaparton-Etanonbjet"]->Fill(jet[1].E(),parton[1].Eta()-jet[1].Eta());
  histoRes2D["Eparton_vs_Etabjet"]->Fill(parton[2].E(),jet[2].Eta());
  histoRes2D["Ebjet_vs_Etaparton-Etabjet"]->Fill(jet[2].E(),parton[2].Eta()-jet[2].Eta());
  if ( nJets > 3 )
  {
    histoRes2D["Eparton_vs_Etabjet"]->Fill(parton[3].E(),jet[3].Eta());
    histoRes2D["Ebjet_vs_Etaparton-Etabjet"]->Fill(jet[3].E(),parton[3].Eta()-jet[3].Eta());
  }

  /// Phi
  histoRes2D["Phiparton_vs_Phinonbjet"]->Fill(parton[0].Phi(),jet[0].Phi());
  histoRes2D["Phiparton_vs_Phinonbjet"]->Fill(parton[1].Phi(),jet[1].Phi());
  histoRes2D["Eparton_vs_Phinonbjet"]->Fill(parton[0].E(),jet[0].Phi());
  histoRes2D["Eparton_vs_Phinonbjet"]->Fill(parton[1].E(),jet[1].Phi());
  float DeltaPhi_nonbjet1 = ROOT::Math::VectorUtil::DeltaPhi(parton[0],jet[0]);
  histoRes2D["Phinonbjet_vs_Phiparton-Phinonbjet"]->Fill(jet[0].Phi(),DeltaPhi_nonbjet1);
  histoRes2D["Enonbjet_vs_Phiparton-Phinonbjet"]->Fill(jet[0].E(),DeltaPhi_nonbjet1);
  float DeltaPhi_nonbjet2 = ROOT::Math::VectorUtil::DeltaPhi(parton[1],jet[1]);
  histoRes2D["Phinonbjet_vs_Phiparton-Phinonbjet"]->Fill(jet[1].Phi(),DeltaPhi_nonbjet2);
  histoRes2D["Enonbjet_vs_Phiparton-Phinonbjet"]->Fill(jet[1].E(),DeltaPhi_nonbjet2);

  histoRes2D["Phiparton_vs_Phibjet"]->Fill(parton[2].Phi(),jet[2].Phi());
  if ( nJets > 3 ) histoRes2D["Phiparton_vs_Phibjet"]->Fill(parton[3].Phi(),jet[3].Phi());
  histoRes2D["Eparton_vs_Phibjet"]->Fill(parton[2].E(),jet[2].Phi());
  if ( nJets > 3 ) histoRes2D["Eparton_vs_Phibjet"]->Fill(parton[3].E(),jet[3].Phi());
  float DeltaPhi_bjet1 = ROOT::Math::VectorUtil::DeltaPhi(parton[2],jet[2]);
  histoRes2D["Phibjet_vs_Phiparton-Phibjet"]->Fill(jet[2].Phi(),DeltaPhi_bjet1);
  histoRes2D["Ebjet_vs_Phiparton-Phibjet"]->Fill(jet[2].E(),DeltaPhi_bjet1);
  float DeltaPhi_bjet2 = 999.;
  if ( nJets > 3 )
  {
    DeltaPhi_bjet2 = ROOT::Math::VectorUtil::DeltaPhi(parton[3],jet[3]);
    histoRes2D["Phibjet_vs_Phiparton-Phibjet"]->Fill(jet[3].Phi(),DeltaPhi_bjet2);
    histoRes2D["Ebjet_vs_Phiparton-Phibjet"]->Fill(jet[3].E(),DeltaPhi_bjet2);
  }
  
  /// Et, theta, phi binned in eta
  float eta0 = fabs(jet[0].Eta());
  if ( eta0 <= 1.3 )
  {
    histoRes2D["Etparton_vs_Etnonbjet_B"]->Fill(parton[0].Et(),jet[0].Et());
    histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_B"]->Fill(jet[0].Et(),parton[0].Et()-jet[0].Et());
    histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_B"]->Fill(jet[0].Et(),parton[0].Theta()-jet[0].Theta());
    histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_B"]->Fill(jet[0].Et(),DeltaPhi_nonbjet1);
  }
  else if ( eta0 > 1.3 && eta0 <= 1.5 )
  {
    histoRes2D["Etparton_vs_Etnonbjet_O"]->Fill(parton[0].Et(),jet[0].Et());
    histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_O"]->Fill(jet[0].Et(),parton[0].Et()-jet[0].Et());
    histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_O"]->Fill(jet[0].Et(),parton[0].Theta()-jet[0].Theta());
    histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_O"]->Fill(jet[0].Et(),DeltaPhi_nonbjet1);
  }
  else
  {
    histoRes2D["Etparton_vs_Etnonbjet_E"]->Fill(parton[0].Et(),jet[0].Et());
    histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_E"]->Fill(jet[0].Et(),parton[0].Et()-jet[0].Et());
    histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_E"]->Fill(jet[0].Et(),parton[0].Theta()-jet[0].Theta());
    histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_E"]->Fill(jet[0].Et(),DeltaPhi_nonbjet1);
  }
  float eta1 = fabs(jet[1].Eta());
  if ( eta1 <= 1.3 )
  {
    histoRes2D["Etparton_vs_Etnonbjet_B"]->Fill(parton[1].Et(),jet[1].Et());
    histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_B"]->Fill(jet[1].Et(),parton[1].Et()-jet[1].Et());
    histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_B"]->Fill(jet[1].Et(),parton[1].Theta()-jet[1].Theta());
    histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_B"]->Fill(jet[1].Et(),DeltaPhi_nonbjet2);
  }
  else if ( eta1 > 1.3 && eta1 <= 1.5 )
  {
    histoRes2D["Etparton_vs_Etnonbjet_O"]->Fill(parton[1].Et(),jet[1].Et());
    histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_O"]->Fill(jet[1].Et(),parton[1].Et()-jet[1].Et());
    histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_O"]->Fill(jet[1].Et(),parton[1].Theta()-jet[1].Theta());
    histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_O"]->Fill(jet[1].Et(),DeltaPhi_nonbjet2);
  }
  else
  {
    histoRes2D["Etparton_vs_Etnonbjet_E"]->Fill(parton[1].Et(),jet[1].Et());
    histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_E"]->Fill(jet[1].Et(),parton[1].Et()-jet[1].Et());
    histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_E"]->Fill(jet[1].Et(),parton[1].Theta()-jet[1].Theta());
    histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_E"]->Fill(jet[1].Et(),DeltaPhi_nonbjet2);
  }
  float eta2 = fabs(jet[2].Eta());
  if ( eta2 <= 1.3 )
  {
    histoRes2D["Etparton_vs_Etbjet_B"]->Fill(parton[2].Et(),jet[2].Et());
    histoRes2D["Etbjet_vs_Etparton-Etbjet_B"]->Fill(jet[2].Et(),parton[2].Et()-jet[2].Et());
    histoRes2D["Etbjet_vs_Thparton-Thbjet_B"]->Fill(jet[2].Et(),parton[2].Theta()-jet[2].Theta());
    histoRes2D["Etbjet_vs_Phiparton-Phibjet_B"]->Fill(jet[2].Et(),DeltaPhi_bjet1);
  }
  else if ( eta2 > 1.3 && eta2 <= 1.5 )
  {
    histoRes2D["Etparton_vs_Etbjet_O"]->Fill(parton[2].Et(),jet[2].Et());
    histoRes2D["Etbjet_vs_Etparton-Etbjet_O"]->Fill(jet[2].Et(),parton[2].Et()-jet[2].Et());
    histoRes2D["Etbjet_vs_Thparton-Thbjet_O"]->Fill(jet[2].Et(),parton[2].Theta()-jet[2].Theta());
    histoRes2D["Etbjet_vs_Phiparton-Phibjet_O"]->Fill(jet[2].Et(),DeltaPhi_bjet1);
  }
  else
  {
    histoRes2D["Etparton_vs_Etbjet_E"]->Fill(parton[2].Et(),jet[2].Et());
    histoRes2D["Etbjet_vs_Etparton-Etbjet_E"]->Fill(jet[2].Et(),parton[2].Et()-jet[2].Et());
    histoRes2D["Etbjet_vs_Thparton-Thbjet_E"]->Fill(jet[2].Et(),parton[2].Theta()-jet[2].Theta());
    histoRes2D["Etbjet_vs_Phiparton-Phibjet_E"]->Fill(jet[2].Et(),DeltaPhi_bjet1);
  }
  if ( nJets > 3 )
  {
    float eta3 = fabs(jet[3].Eta());
    if ( eta3 <= 1.3 )
    {
      histoRes2D["Etparton_vs_Etbjet_B"]->Fill(parton[3].Et(),jet[3].Et());
      histoRes2D["Etbjet_vs_Etparton-Etbjet_B"]->Fill(jet[3].Et(),parton[3].Et()-jet[3].Et());
      histoRes2D["Etbjet_vs_Thparton-Thbjet_B"]->Fill(jet[3].Et(),parton[3].Theta()-jet[3].Theta());
      histoRes2D["Etbjet_vs_Phiparton-Phibjet_B"]->Fill(jet[3].Et(),DeltaPhi_bjet2);
    }
    else if ( eta3 > 1.3 && eta3 <= 1.5 )
    {
      histoRes2D["Etparton_vs_Etbjet_O"]->Fill(parton[3].Et(),jet[3].Et());
      histoRes2D["Etbjet_vs_Etparton-Etbjet_O"]->Fill(jet[3].Et(),parton[3].Et()-jet[3].Et());
      histoRes2D["Etbjet_vs_Thparton-Thbjet_O"]->Fill(jet[3].Et(),parton[3].Theta()-jet[3].Theta());
      histoRes2D["Etbjet_vs_Phiparton-Phibjet_O"]->Fill(jet[3].Et(),DeltaPhi_bjet2);
    }
    else
    {
      histoRes2D["Etparton_vs_Etbjet_E"]->Fill(parton[3].Et(),jet[3].Et());
      histoRes2D["Etbjet_vs_Etparton-Etbjet_E"]->Fill(jet[3].Et(),parton[3].Et()-jet[3].Et());
      histoRes2D["Etbjet_vs_Thparton-Thbjet_E"]->Fill(jet[3].Et(),parton[3].Theta()-jet[3].Theta());
      histoRes2D["Etbjet_vs_Phiparton-Phibjet_E"]->Fill(jet[3].Et(),DeltaPhi_bjet2);
    }
  }
}

void ResolutionFunctions::fillMuon(TLorentzVector genMu, TLorentzVector recMu)
{
  
  muon = true;
  
  histoRes2D["EgenMu_vs_ErecMu"]->Fill(genMu.E(),recMu.E());
  histoRes2D["ErecMu_vs_EgenMu-ErecMu"]->Fill(recMu.E(),genMu.E()-recMu.E());
  
  histoRes2D["EtgenMu_vs_EtrecMu"]->Fill(genMu.Et(),recMu.Et());
  histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu"]->Fill(recMu.Et(),genMu.Et()-recMu.Et());
  
  histoRes2D["PtgenMu_vs_PtrecMu"]->Fill(genMu.Pt(),recMu.Pt());
  histoRes2D["PtrecMu_vs_PtgenMu-PtrecMu"]->Fill(recMu.Pt(),genMu.Pt()-recMu.Pt());
  
  histoRes2D["ThgenMu_vs_ThrecMu"]->Fill(genMu.Theta(),recMu.Theta());
  histoRes2D["ThrecMu_vs_ThgenMu-ThrecMu"]->Fill(recMu.Theta(),genMu.Theta()-recMu.Theta());
  histoRes2D["PtgenMu_vs_ThrecMu"]->Fill(genMu.Pt(),recMu.Theta());
  histoRes2D["PtrecMu_vs_ThgenMu-ThrecMu"]->Fill(recMu.Pt(),genMu.Theta()-recMu.Theta());
  
  histoRes2D["EtagenMu_vs_EtarecMu"]->Fill(genMu.Eta(),recMu.Eta());
  histoRes2D["EtarecMu_vs_EtagenMu-EtarecMu"]->Fill(recMu.Eta(),genMu.Eta()-recMu.Eta());
  histoRes2D["PtgenMu_vs_EtarecMu"]->Fill(genMu.Pt(),recMu.Eta());
  histoRes2D["PtrecMu_vs_EtagenMu-EtarecMu"]->Fill(recMu.Pt(),genMu.Eta()-recMu.Eta());
  
  histoRes2D["PhigenMu_vs_PhirecMu"]->Fill(genMu.Phi(),recMu.Phi());
  histoRes2D["PtgenMu_vs_PhirecMu"]->Fill(genMu.Pt(),recMu.Phi());
  float DeltaPhi = ROOT::Math::VectorUtil::DeltaPhi(genMu,recMu);
  histoRes2D["PhirecMu_vs_PhigenMu-PhirecMu"]->Fill(recMu.Phi(),DeltaPhi);
  histoRes2D["PtrecMu_vs_PhigenMu-PhirecMu"]->Fill(recMu.Pt(),DeltaPhi);
  
  /// Et, theta, phi binned in eta
  float eta = fabs(recMu.Eta());
  if ( eta <= 0.9 )
  {
    histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu_B"]->Fill(recMu.Et(),genMu.Et()-recMu.Et());
    histoRes2D["EtrecMu_vs_ThgenMu-ThrecMu_B"]->Fill(recMu.Et(),genMu.Theta()-recMu.Theta());
    histoRes2D["EtrecMu_vs_PhigenMu-PhirecMu_B"]->Fill(recMu.Et(),DeltaPhi);
  }
  else if ( eta > 0.9 && eta <= 1.2 )
  {
    histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu_O"]->Fill(recMu.Et(),genMu.Et()-recMu.Et());
    histoRes2D["EtrecMu_vs_ThgenMu-ThrecMu_O"]->Fill(recMu.Et(),genMu.Theta()-recMu.Theta());
    histoRes2D["EtrecMu_vs_PhigenMu-PhirecMu_O"]->Fill(recMu.Et(),DeltaPhi);
  }
  else
  {
    histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu_E"]->Fill(recMu.Et(),genMu.Et()-recMu.Et());
    histoRes2D["EtrecMu_vs_ThgenMu-ThrecMu_E"]->Fill(recMu.Et(),genMu.Theta()-recMu.Theta());
    histoRes2D["EtrecMu_vs_PhigenMu-PhirecMu_E"]->Fill(recMu.Et(),DeltaPhi);
  }
  
}

void ResolutionFunctions::fillElectron(TLorentzVector genEl, TLorentzVector recEl)
{
  
  electron = true;
  
  histoRes2D["EgenEl_vs_ErecEl"]->Fill(genEl.E(),recEl.E());
  histoRes2D["ErecEl_vs_EgenEl-ErecEl"]->Fill(recEl.E(),genEl.E()-recEl.E());
  
  histoRes2D["EtgenEl_vs_EtrecEl"]->Fill(genEl.Et(),recEl.Et());
  histoRes2D["EtrecEl_vs_EtgenEl-EtrecEl"]->Fill(recEl.Et(),genEl.Et()-recEl.Et());
  
  histoRes2D["PtgenEl_vs_PtrecEl"]->Fill(genEl.Pt(),recEl.Pt());
  histoRes2D["PtrecEl_vs_PtgenEl-PtrecEl"]->Fill(recEl.Pt(),genEl.Pt()-recEl.Pt());
  
  histoRes2D["ThgenEl_vs_ThrecEl"]->Fill(genEl.Theta(),recEl.Theta());
  histoRes2D["ThrecEl_vs_ThgenEl-ThrecEl"]->Fill(recEl.Theta(),genEl.Theta()-recEl.Theta());
  histoRes2D["EgenEl_vs_ThrecEl"]->Fill(genEl.E(),recEl.Theta());
  histoRes2D["ErecEl_vs_ThgenEl-ThrecEl"]->Fill(recEl.E(),genEl.Theta()-recEl.Theta());
  
  histoRes2D["EtagenEl_vs_EtarecEl"]->Fill(genEl.Eta(),recEl.Eta());
  histoRes2D["EtarecEl_vs_EtagenEl-EtarecEl"]->Fill(recEl.Eta(),genEl.Eta()-recEl.Eta());
  histoRes2D["EgenEl_vs_EtarecEl"]->Fill(genEl.E(),recEl.Eta());
  histoRes2D["ErecEl_vs_EtagenEl-EtarecEl"]->Fill(recEl.E(),genEl.Eta()-recEl.Eta());
  
  histoRes2D["PhigenEl_vs_PhirecEl"]->Fill(genEl.Phi(),recEl.Phi());
  histoRes2D["EgenEl_vs_PhirecEl"]->Fill(genEl.E(),recEl.Phi());
  float DeltaPhi = ROOT::Math::VectorUtil::DeltaPhi(genEl,recEl);
  histoRes2D["PhirecEl_vs_PhigenEl-PhirecEl"]->Fill(recEl.Phi(),DeltaPhi);
  histoRes2D["ErecEl_vs_PhigenEl-PhirecEl"]->Fill(recEl.E(),DeltaPhi);
  
}

void ResolutionFunctions::writeHistograms()
{
  /// Define output file outside of function
  // Give it as argument?
  
  histoRes2D["Eparton_vs_Enonbjet"]->Write();
  histoRes2D["Eparton_vs_Ebjet"]->Write();
  histoRes2D["Etparton_vs_Etnonbjet"]->Write();
  histoRes2D["Etparton_vs_Etbjet"]->Write();
  histoRes2D["Ptparton_vs_Ptnonbjet"]->Write();
  histoRes2D["Ptparton_vs_Ptbjet"]->Write();
  histoRes2D["Thparton_vs_Thnonbjet"]->Write();
  histoRes2D["Thparton_vs_Thbjet"]->Write();
  histoRes2D["Etaparton_vs_Etanonbjet"]->Write();
  histoRes2D["Etaparton_vs_Etabjet"]->Write();
  histoRes2D["Phiparton_vs_Phinonbjet"]->Write();
  histoRes2D["Phiparton_vs_Phibjet"]->Write();
  histoRes2D["Eparton_vs_Thnonbjet"]->Write();
  histoRes2D["Eparton_vs_Thbjet"]->Write();
  histoRes2D["Eparton_vs_Etanonbjet"]->Write();
  histoRes2D["Eparton_vs_Etabjet"]->Write();
  histoRes2D["Eparton_vs_Phinonbjet"]->Write();
  histoRes2D["Eparton_vs_Phibjet"]->Write();
  
  histoRes2D["Enonbjet_vs_Eparton-Enonbjet"]->Write();
  histoRes2D["Ebjet_vs_Eparton-Ebjet"]->Write();
  histoRes2D["Enonbjet_vs_Thparton-Thnonbjet"]->Write();
  histoRes2D["Ebjet_vs_Thparton-Thbjet"]->Write();
  histoRes2D["Enonbjet_vs_Etaparton-Etanonbjet"]->Write();
  histoRes2D["Ebjet_vs_Etaparton-Etabjet"]->Write();
  histoRes2D["Enonbjet_vs_Phiparton-Phinonbjet"]->Write();
  histoRes2D["Ebjet_vs_Phiparton-Phibjet"]->Write();
  
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet"]->Write();
  histoRes2D["Etbjet_vs_Etparton-Etbjet"]->Write();
  histoRes2D["Ptnonbjet_vs_Ptparton-Ptnonbjet"]->Write();
  histoRes2D["Ptbjet_vs_Ptparton-Ptbjet"]->Write();
  histoRes2D["Thnonbjet_vs_Thparton-Thnonbjet"]->Write();
  histoRes2D["Thbjet_vs_Thparton-Thbjet"]->Write();
  histoRes2D["Etanonbjet_vs_Etaparton-Etanonbjet"]->Write();
  histoRes2D["Thbjet_vs_Thparton-Thbjet"]->Write();
  histoRes2D["Phinonbjet_vs_Phiparton-Phinonbjet"]->Write();
  histoRes2D["Phibjet_vs_Phiparton-Phibjet"]->Write();
  
  // binned in eta
  histoRes2D["Etparton_vs_Etnonbjet_B"]->Write();
  histoRes2D["Etparton_vs_Etnonbjet_O"]->Write();
  histoRes2D["Etparton_vs_Etnonbjet_E"]->Write();
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_B"]->Write();
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_O"]->Write();
  histoRes2D["Etnonbjet_vs_Etparton-Etnonbjet_E"]->Write();
  histoRes2D["Etparton_vs_Etbjet_B"]->Write();
  histoRes2D["Etparton_vs_Etbjet_O"]->Write();
  histoRes2D["Etparton_vs_Etbjet_E"]->Write();
  histoRes2D["Etbjet_vs_Etparton-Etbjet_B"]->Write();
  histoRes2D["Etbjet_vs_Etparton-Etbjet_O"]->Write();
  histoRes2D["Etbjet_vs_Etparton-Etbjet_E"]->Write();
  
  histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_B"]->Write();
  histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_O"]->Write();
  histoRes2D["Etnonbjet_vs_Thparton-Thnonbjet_E"]->Write();
  histoRes2D["Etbjet_vs_Thparton-Thbjet_B"]->Write();
  histoRes2D["Etbjet_vs_Thparton-Thbjet_O"]->Write();
  histoRes2D["Etbjet_vs_Thparton-Thbjet_E"]->Write();
  histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_B"]->Write();
  histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_O"]->Write();
  histoRes2D["Etnonbjet_vs_Phiparton-Phinonbjet_E"]->Write();
  histoRes2D["Etbjet_vs_Phiparton-Phibjet_B"]->Write();
  histoRes2D["Etbjet_vs_Phiparton-Phibjet_O"]->Write();
  histoRes2D["Etbjet_vs_Phiparton-Phibjet_E"]->Write();
  
  
  if (muon)
  {
    histoRes2D["EgenMu_vs_ErecMu"]->Write();
    histoRes2D["EtgenMu_vs_EtrecMu"]->Write();
    histoRes2D["PtgenMu_vs_PtrecMu"]->Write();
    histoRes2D["ThgenMu_vs_ThrecMu"]->Write();
    histoRes2D["EtagenMu_vs_EtarecMu"]->Write();
    histoRes2D["PhigenMu_vs_PhirecMu"]->Write();
    histoRes2D["PtgenMu_vs_ThrecMu"]->Write();
    histoRes2D["PtgenMu_vs_EtarecMu"]->Write();
    histoRes2D["PtgenMu_vs_PhirecMu"]->Write();
    
    histoRes2D["PtrecMu_vs_PtgenMu-PtrecMu"]->Write();
    histoRes2D["PtrecMu_vs_ThgenMu-ThrecMu"]->Write();
    histoRes2D["PtrecMu_vs_EtagenMu-EtarecMu"]->Write();
    histoRes2D["PtrecMu_vs_PhigenMu-PhirecMu"]->Write();
    
    histoRes2D["ErecMu_vs_EgenMu-ErecMu"]->Write();
    histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu"]->Write();
    histoRes2D["ThrecMu_vs_ThgenMu-ThrecMu"]->Write();
    histoRes2D["EtarecMu_vs_EtagenMu-EtarecMu"]->Write();
    histoRes2D["PhirecMu_vs_PhigenMu-PhirecMu"]->Write();
    
    histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu_B"]->Write();
    histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu_O"]->Write();
    histoRes2D["EtrecMu_vs_EtgenMu-EtrecMu_E"]->Write();
    histoRes2D["EtrecMu_vs_ThgenMu-ThrecMu_B"]->Write();
    histoRes2D["EtrecMu_vs_ThgenMu-ThrecMu_O"]->Write();
    histoRes2D["EtrecMu_vs_ThgenMu-ThrecMu_E"]->Write();
    histoRes2D["EtrecMu_vs_PhigenMu-PhirecMu_B"]->Write();
    histoRes2D["EtrecMu_vs_PhigenMu-PhirecMu_O"]->Write();
    histoRes2D["EtrecMu_vs_PhigenMu-PhirecMu_E"]->Write();
  }
  
  if (electron)
  {
    histoRes2D["EgenEl_vs_ErecEl"]->Write();
    histoRes2D["EtgenEl_vs_EtrecEl"]->Write();
    histoRes2D["PtgenEl_vs_PtrecEl"]->Write();
    histoRes2D["ThgenEl_vs_ThrecEl"]->Write();
    histoRes2D["EtagenEl_vs_EtarecEl"]->Write();
    histoRes2D["PhigenEl_vs_PhirecEl"]->Write();
    histoRes2D["EgenEl_vs_ThrecEl"]->Write();
    histoRes2D["EgenEl_vs_EtarecEl"]->Write();
    histoRes2D["EgenEl_vs_PhirecEl"]->Write();
    
    histoRes2D["ErecEl_vs_EgenEl-ErecEl"]->Write();
    histoRes2D["ErecEl_vs_ThgenEl-ThrecEl"]->Write();
    histoRes2D["ErecEl_vs_EtagenEl-EtarecEl"]->Write();
    histoRes2D["ErecEl_vs_PhigenEl-PhirecEl"]->Write();
    
    histoRes2D["EtrecEl_vs_EtgenEl-EtrecEl"]->Write();
    histoRes2D["PtrecEl_vs_PtgenEl-PtrecEl"]->Write();
    histoRes2D["ThrecEl_vs_ThgenEl-ThrecEl"]->Write();
    histoRes2D["EtarecEl_vs_EtagenEl-EtarecEl"]->Write();
    histoRes2D["PhirecEl_vs_PhigenEl-PhirecEl"]->Write();
  }
  
  if (verbose) std::cout << "ResolutionFunctions::writeHistograms - Histograms written to file" << std::endl;
  //this->makeFit();
}

void ResolutionFunctions::makeFit()
{
  if (verbose) std::cout << "                             - Starting fit procedure... " << std::endl;
  if (verbose) std::cout << "Single Gaussian? " << useSingleG << std::endl;
  
  for (int f = 0; f < nHistos; f++)
  {
    // 0-5  : b jet, 6-11 : non-b jet, 12-17: muon, 18-23: electron, 24-26: b-jet Et in eta bins, 27-29: non-b jet Et in eta bins
    //if (! muon && (f == 12 || f == 13 || f == 14 || f == 15 || f == 16 || f == 17) ) continue;
    //if (! electron && (f == 18 || f == 19 || f == 20 || f == 21 || f == 22 || f == 23) ) continue;
    
    if (electron) std::cout << "ResolutionFunctions::WriteOutputFiles -- WARNING: Electron fitting not yet implemented" << std::endl;
    //if ( f == 4 || f == 10 || f == 16 || f == 22 ) continue;  // eta later
    
    // now: 0-2: b jet Et (BOE), 3-5: b jet theta (BOE), 6-8: b jet phi (BOE), 9-11: non-b jet Et (BOE), 12-14: non-b jet theta (BOE), 15-17: non-b jet phi (BOE)
    
    bool fitEt = false, fitTheta = false, fitPhi = false, fitEtMu = false, fitThetaMu = false, fitPhiMu = false;
    if      ( f ==  0 || f ==  1 || f ==  2 || f ==  9 || f == 10 || f == 11 ) fitEt = true;
    else if ( f ==  3 || f ==  4 || f ==  5 || f == 12 || f == 13 || f == 14 ) fitTheta = true;
    else if ( f ==  6 || f ==  7 || f ==  8 || f == 15 || f == 16 || f == 17 ) fitPhi = true;
    else if ( f == 18 || f == 19 || f == 20 ) fitEtMu = true;
    else if ( f == 21 || f == 22 || f == 23 ) fitThetaMu = true;
    else if ( f == 24 || f == 25 || f == 26 ) fitPhiMu = true;
    else continue;
    
    bool fitAnglesMuon = false;
    //if ( fitThetaMu || fitPhiMu ) fitAnglesMuon = true;
    
    //if ( f < 9 ) continue;  // TEMPORARILY (disable b jet fits)
    
    std::cout << "  ***Current histogram:  " << histoNames[f] << "***" << std::endl;
    
    TH2F* histo=histoRes2D[(histoNames[f]).c_str()];
    if (getHistos)
      histo=fitHisto2D[(histoNames[f]).c_str()];
    
    
    int nBins = histo->GetXaxis()->GetNbins();
    if (verbose) std::cout << "nbins: " << nBins << std::endl;
    int nPar = 6;
    //int nPar = 5;
    
    /// Create one histogram for each function parameter -> 6 histograms for each 2D plot
    TH1D **hlist = new TH1D*[nPar];
    std::string parnames[6]={"a1","a2","a3","a4","a5","a6"};
    //std::string parnames[5]={"a1","a2","a3","a4","a5"};
    std::string name=""; std::string title="";
    //const TArrayD *bins = histo->GetXaxis()->GetXbins();
    for (int iPar = 0; iPar < nPar; iPar++)
    {
      name = std::string(histo->GetName())+ "_" + parnames[iPar];
      title = std::string(histo->GetName())+ ": Fitted value of " + parnames[iPar] ;
      hlist[iPar] = new TH1D(name.c_str(), title.c_str(), nBins, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
      hlist[iPar]->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    }
    
    /// Loop on all bins in X, generate a projection along Y
    int cut = 0;  // require a minimum number of bins in the slice to be filled
    for (int xBin = 1; xBin < nBins+1; xBin++)
    {
      /// Make projection
      std::string projection_title = std::string(histo->GetName())+"_sliceXbin"+toStr(xBin);
      TH1D *hp = histo->ProjectionY(projection_title.c_str(), xBin, xBin, "e");
      if(xBin == nBins) hp = histo->ProjectionY(projection_title.c_str(), xBin, xBin+1, "e");  //include overflow in last bin
      if(xBin == 1) hp = histo->ProjectionY(projection_title.c_str(), xBin-1, xBin, "e");  //include underflow in first bin
      if (hp == 0) continue;
      float nEntries = float(hp->GetEntries());
      if (nEntries == 0 || nEntries < cut) {delete hp; continue;}
      
      hp->GetXaxis()->SetTitle(histo->GetYaxis()->GetTitle());
      hp->GetYaxis()->SetTitle("");
      
      /// Normalise histogram
      Double_t scale = 1./hp->Integral();
      hp->Scale(scale);
      
      if (verbose) std::cout << "Integral of the histo is " << hp->Integral() << std::endl;
      /// Declare the fit function
      //  ! Its range depends on the jet/lepton energy range (hence, the Y-axis)
      
      double maxX = 0.;
      double fitEdge = 100.;  // largest for Et
      if (fitTheta)
      {
        fitEdge = 0.08;
        if ( xBin > 1 ) fitEdge = 0.06;
        if ( histoNames[f].find("_B") == std::string::npos )
        {
          if ( xBin > 2 ) fitEdge = 0.04;
          if ( xBin > 4 ) fitEdge = 0.03;
          if ( histoNames[f].find("_E") != std::string::npos && xBin > 3 ) fitEdge = 0.02;
        }
        else
        {
          if ( xBin > 4 ) fitEdge = 0.05;
        }
      }
      else if (fitPhi)
      {
        fitEdge = 0.18;
        if ( xBin > 1 ) fitEdge = 0.13;
        if ( xBin > 2 ) fitEdge = 0.08;
        if ( xBin > 3 ) fitEdge = 0.06;
        if ( histoNames[f].find("_B") == std::string::npos && xBin > 4 ) fitEdge = 0.05;
      }
      else if (fitEtMu)
      {
        fitEdge = 3.;
        if ( xBin == 1 ) fitEdge = 2.5;
      }
      else if (fitThetaMu)
      {
        fitEdge = 0.0015;
        if ( histoNames[f].find("_B") == std::string::npos ) fitEdge = 0.001;
      }
      else if (fitPhiMu)
      {
        fitEdge = 0.0015;
      }
      
      if (useSingleG)
      {
        maxX = hp->GetXaxis()->GetBinCenter(hp->GetMaximumBin());
        //maxX = hp->GetMean();
        fitEdge = 0.8*(hp->GetStdDev());
        if ( xBin > 3 ) fitEdge = 1.*(hp->GetStdDev());
        if ( fitEt )
        {
          if (histoNames[f].find("_E") != std::string::npos) fitEdge = 0.9*(hp->GetStdDev());
          if ( xBin > 3 ) fitEdge = 1.5*(hp->GetStdDev());
        }
      }
      if (useSingleG || fitAnglesMuon /*&& ! fitEt*/) nPar = 4;
      
      TF1 *myfit;
      if (useSingleG /*&& ! fitEt*/) myfit = new TF1("myfit", sGaus, maxX - fitEdge, maxX + fitEdge, nPar);
      else if (fitAnglesMuon) myfit = new TF1("myfit", sGaus, -fitEdge, fitEdge, nPar);
      else myfit = new TF1("myfit", dblGaus, -fitEdge, fitEdge, nPar);
      
      //  Give names to the parameters
      myfit->SetParName(0,"a1");  // central value of first, narrow gaussian
      myfit->SetParName(1,"a2");  // sigma value of first, narrow gaussian
      myfit->SetParName(2,"a3");  // relative scale factor gaussians
                                  // sGaus: amplitude
      myfit->SetParName(3,"a4");  // central value of second, broad gaussian
                                  // sGaus: shift up/down (because of second, non-fitted gaus)
      if (! useSingleG && ! fitAnglesMuon /*|| fitEt*/)
      {
        myfit->SetParName(4,"a5");  // sigma value of second, broad gaussian
        myfit->SetParName(5,"a6");  // amplitude
      }
      //  Set initial values
      //if ( f == 0 || f == 1 || f == 2 || f == 6 || f == 7 || f == 8 || f == 24 || f == 25 || f == 26 || f == 27 || f == 28 || f == 29)  // energy, Et, pt of (non-)b jet, with & without eta bins
      if ( fitEt )  // Et of (non-)b jet (BOE)
      {
        // Mean around zero
        //myfit->SetParLimits(0,-10,10);
        myfit->SetParLimits(0, hp->GetMean() - hp->GetRMS()/1.5, hp->GetMean() + hp->GetRMS()/1.5);
        if (! useSingleG) myfit->SetParLimits(3, -50., 35.);
        // Restrict sigma to be positive
        myfit->SetParLimits(1, 2., 28.);
        if (! useSingleG) myfit->SetParLimits(4, 8., 110.);
        
        if ( xBin == 1 )
        {
          if (useSingleG) myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 1., 0.);
          else myfit->SetParameters(hp->GetMean(), hp->GetRMS()/4., 0.3, 2., 10., 1.);
          if (! useSingleG)
          {
            myfit->SetParLimits(1, 2., 8.);
            myfit->SetParLimits(4, 9., 110.);
          }
          else myfit->SetParLimits(4, 7., 110.);
        }
        else if ( xBin == 2 )
        {
          if (useSingleG) myfit->SetParameters(hp->GetMean(), hp->GetRMS()/2., 1., 0.);
          else myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 0.3, 2., 14.4, 1.);
          if (! useSingleG) myfit->SetParLimits(1, 2., 10.);
        }
        else if ( xBin == 3 || xBin == 4 )
        {
          if (useSingleG) myfit->SetParameters(hp->GetMean(), hp->GetRMS()/2., 1., 0.);
          else myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 0.9, 5., 30.0, 1.);
          if (! useSingleG)
          {
            if ( xBin == 3 ) myfit->SetParLimits(1, 8., 11.);
            else myfit->SetParLimits(1, 11., 13.);
            myfit->SetParLimits(3, -25., 25.);
            myfit->SetParLimits(4, 18., 50.);
          }
        }
        else if ( xBin == 5 )
        {
          if (useSingleG) myfit->SetParameters(hp->GetMean(), hp->GetRMS()/2., 1., 0.);
          else myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 2., 6.6, 50., 1.);
          if (! useSingleG)
          {
            myfit->SetParLimits(1, 12.5, 15.);
            myfit->SetParLimits(4, 25., 125.);
          }
          //if ( f > 6 )  // Et, pt of non-b jet
            //fix such that small gaussian has smallest width
        }
        else if ( xBin == 6 )
        {
          if (useSingleG) myfit->SetParameters(hp->GetMean(), hp->GetRMS()/2., 1., 0.);
          else myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 2., 6.6, 50., 1.);
          if (! useSingleG)
          {
            myfit->SetParLimits(1, 13.5, 17.);
            myfit->SetParLimits(4, 28., 135.);
          }
          //if ( f > 6 )  // Et, pt of non-b jet
            //fix such that small gaussian has smallest width
        }
        
      }
      //else if ( f == 3 || f == 5 || f == 9 || f == 11 )  // theta & phi of jet (CHECK: also eta??)
      else if ( fitTheta || fitPhi )  // theta & phi of (non-b) jet (BOE)
      {
        // Mean around zero
        myfit->SetParLimits(0,-0.005,0.005);
        // Restrict sigma to be positive
        myfit->SetParLimits(1,1e-3,15e-3);
        if (fitPhi) myfit->SetParLimits(1,5e-3,20e-3);
        
        if (! useSingleG)
        {
          if (fitTheta) myfit->SetParLimits(2,0.05,0.6);
          else myfit->SetParLimits(2,0.06,0.6);
          myfit->SetParLimits(3,-0.05,0.05);  // mean around zero
          myfit->SetParLimits(4,1e-2,3e-1);   // sigma positive
        }
        
        if (useSingleG) myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 1., 0.);
        else if (fitPhi) myfit->SetParameters(hp->GetMean(), hp->GetRMS()/4., 0.3, hp->GetMean(), hp->GetRMS()*1.5, 0.03);
        else myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 0.3, hp->GetMean(), hp->GetRMS()*2., 0.03);
        
        if ( xBin == 1 )
        {
          if (! useSingleG)
          {
            if (fitPhi)
            {
              myfit->SetParameter(1, 0.023);
              myfit->SetParLimits(1,20e-3,25e-3);
            }
            else myfit->SetParLimits(1,5e-3,18e-3);
            myfit->SetParLimits(4,2e-2,6e-1);
          }
        }
        else if ( xBin > 1 && xBin < 4 )
        {
          if (! useSingleG)
          {
            if ( xBin == 2 )
            {
              if (fitTheta) myfit->SetParLimits(1, 5e-3, 13e-3);
              else if (fitPhi)
              {
                myfit->SetParameter(1, 0.018);
                myfit->SetParLimits(1,13e-3,20e-3);
                myfit->SetParLimits(4,2e-2,6e-1);
                if ( histoNames[f].find("_B") != std::string::npos ) myfit->SetParLimits(1,13e-3,18e-3);
              }
            }
            else if ( xBin == 3 && fitTheta && histoNames[f].find("_E") != std::string::npos )
            {
              myfit->SetParLimits(1,3e-3,5e-3);
            }
            myfit->SetParameter(2, 0.2);
          }
        }
        else
        {
          myfit->SetParLimits(1, 2e-3, 1e-2);
          if (fitPhi) myfit->SetParLimits(1,5e-3,12e-3);
          if (! useSingleG)
          {
            if (fitTheta)
            {
              if ( xBin == 4 ) myfit->SetParLimits(1, 25e-4, 8e-3);
              else myfit->SetParLimits(1, 25e-4, 7e-3);
              if ( histoNames[f].find("_E") != std::string::npos )
              {
                myfit->SetParLimits(1, 2e-3, 35e-4);
                if ( xBin > 4 ) myfit->SetParLimits(1, 1e-3, 3e-3);
              }
              else if ( xBin > 4 && histoNames[f].find("_O") != std::string::npos ) myfit->SetParLimits(1, 25e-4, 6e-3);
            }
            myfit->SetParameter(2, 0.1);
            myfit->SetParLimits(2,0.06,0.45);
            myfit->SetParLimits(4,7e-3,3e-1);
            if (fitPhi) myfit->SetParLimits(4,15e-3,3e-1);
            if ( fitPhi && xBin > 4 && histoNames[f].find("_E") != std::string::npos ) myfit->SetParLimits(1,5e-3,1e-2);
          }
        }
      }
      
      /////////
      /// MUON
      /////////
      
      else if (fitEtMu)   /// muon Et (BOE)
      {
        // Mean around zero
        myfit->SetParLimits(0,-0.5,0.5);
        // Restrict sigma to be positive
        myfit->SetParLimits(1,5e-2,2.);
        if (! useSingleG)
        {
          myfit->SetParLimits(2,0.01,0.6);
          myfit->SetParLimits(3,-1.,1.);
          myfit->SetParLimits(4,8e-1,30.);
        }
        
        myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 0.3, hp->GetMean(), hp->GetRMS()*1.5, 0.03);
//         if ( xBin < 3 && histoNames[f].find("_B") != std::string::npos )
//         {
//           if ( xBin == 1 ) myfit->SetParLimits(1,5e-2,0.42);
//           else myfit->SetParLimits(1,5e-2,0.85);
//         }
        //if ( xBin == 3 /*&& histoNames[f].find("_B") != std::string::npos*/ ) myfit->SetParLimits(1,5e-2,1.6);
        //if ( ! useSingleG && xBin == 4 && histoNames[f].find("_B") != std::string::npos )
        if ( xBin > 2 )
        {
          if ( xBin == 3 ) myfit->SetParLimits(1,0.8,1.6);
          else myfit->SetParLimits(1,0.8,2.);
          if (! useSingleG) myfit->SetParLimits(2,0.01,1.5);
        }
      }
      else if (fitThetaMu || fitPhiMu)
      {
        // Mean around zero
        myfit->SetParLimits(0,-0.001,0.001);
        // Restrict sigma to be positive
        myfit->SetParLimits(1,9e-5,3e-4);
        if (fitPhiMu) myfit->SetParLimits(1,9e-5,25e-5);
        if (! useSingleG && ! fitAnglesMuon)
        {
          if (fitThetaMu) myfit->SetParLimits(2,0.05,0.8);
          //else myfit->SetParLimits(2,0.06,0.8);
          myfit->SetParLimits(4,2e-4,8e-4);
        }
        
        if (fitAnglesMuon) myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 0.1, 0.);
        else myfit->SetParameters(hp->GetMean(), hp->GetRMS()/3., 0.25, hp->GetMean(), hp->GetRMS()*1.5, 0.03);
        
        if ( histoNames[f].find("_B") != std::string::npos )
        {
          myfit->SetParLimits(1,1e-4,3e-4);
          if (fitPhiMu) myfit->SetParLimits(1,7e-5,25e-5);
          if (! useSingleG)
          {
            if (fitThetaMu) myfit->SetParLimits(2,0.2,0.3);
          }
        }
        //else if ( fitPhiMu && histoNames[f].find("_O") != std::string::npos ) myfit->SetParLimits(1,9e-5,25e-5);
        else if ( histoNames[f].find("_E") != std::string::npos )
        {
          myfit->SetParLimits(1,6e-5,3e-4);
          //if (fitPhiMu) myfit->SetParLimits(1,1e-4,25e-5);
          if (! useSingleG)
          {
            if (fitThetaMu) myfit->SetParLimits(2,0.07,0.3);
          }
        }
      }
/*      else if ( f == 12 )  // energy of muon
      {
        //myfit->SetParameter(0, -0.008);  //central value of first, broad gaussian
        //myfit->SetParameter(1, 0.01);    //sigma value of first, broad gaussian
        //myfit->SetParameter(2, 24.0);    //constant value of first, broad gaussian
        //myfit->SetParameter(3, -0.0002); //central value of second, narrow gaussian
        //myfit->SetParameter(4, 0.00002); //sigma value of second, narrow gaussian
        //myfit->SetParameter(5, 19.0);    //constant value of second, narrow gaussian
        myfit->SetParameter(0, -0.008);  //central value of first, broad gaussian
        myfit->SetParameter(1, 0.002);   //sigma value of first, broad gaussian
        myfit->SetParameter(2, 100.0);   //constant value of second, narrow gaussian
        myfit->SetParameter(3, -0.0002); //central value of second, narrow gaussian
        myfit->SetParameter(4, 0.0002);  //sigma value of second, narrow gaussian
      }
      else if ( f == 13 )  // Et of muon
      {
        //test
      }
      else if ( f == 14 )  // pt of muon
      {
        //test
      }
      else if ( f == 15 || f == 17 )  // theta & phi of muon (CHECK: also eta??)
      {
        //myfit->SetParameter(0, 0.0); 		 //central value of first, broad gaussian
        //myfit->SetParameter(1, 0.01);    //sigma value of first, broad gaussian
        //myfit->SetParameter(2, 24.0); 	 //constant value of first, broad gaussian
        //myfit->SetParameter(3, 0.0);     //central value of second, narrow gaussian
        //myfit->SetParameter(4, 0.00001); //sigma value of second, narrow gaussian
        //myfit->SetParameter(5, 4.0);     //constant value of second, narrow gaussian
        myfit->SetParameter(0, 0.0); 		 //central value of first, broad gaussian
        myfit->SetParameter(1, 0.0005);  //sigma value of first, broad gaussian
        myfit->SetParameter(2, 100.0); 	 //constant value of second, narrow gaussian
        myfit->SetParameter(3, 0.0);     //central value of second, narrow gaussian
        myfit->SetParameter(4, 0.00001); //sigma value of second, narrow gaussian
      }
*/      
      ///  Fit
      std::string func_title = std::string(histo->GetName())+"_sliceXbin"+toStr(xBin)+"_Fitted";
      myfit->SetName(func_title.c_str());
      hp->Fit(myfit, "R");
      gStyle->SetOptFit(0111);
      
      if (! useSingleG && ! fitAnglesMuon && myfit->GetParameter(1) > myfit->GetParameter(4) )
      {
        double temp = myfit->GetParameter(1);
        double tempErr = myfit->GetParError(1);
        double addErr = 0.003;
        if (fitEt) addErr = 3.;
        else if (fitEtMu) addErr = 0.2;
        else if (fitThetaMu || fitPhiMu) addErr = 0.00002;
        myfit->SetParameter(1, myfit->GetParameter(4));
        myfit->SetParError(1, myfit->GetParError(4)+addErr);
        myfit->SetParameter(4, temp);
        myfit->SetParError(4, tempErr);
      }
      
      /// Manually increase error on fit parameter when few events and error is small
      if ( hp->GetEntries() < 2000 )
      {
        double tempErr = myfit->GetParError(1);
        if ( fitEt && tempErr < 1. ) myfit->SetParError(1, 3.*tempErr);
        else if ( tempErr/myfit->GetParameter(1) < 0.01 ) myfit->SetParError(1, 3.*tempErr);
      }
      
      if (! useSingleG && fitEt )
      {
        double tempErr = myfit->GetParError(1);
        if ( tempErr < 5e-2 ) myfit->SetParError(1, 1.);
        if ( xBin == 4 && histoNames[f].find("_O") != std::string::npos )
        {
          if ( tempErr < 0.1 ) myfit->SetParError(1, 3.*tempErr);
          else if ( tempErr < 1. ) myfit->SetParError(1, 2.*tempErr);
          if ( tempErr < 5e-2 ) myfit->SetParError(1, 1.);
        }
        if ( xBin < 3 && histoNames[f].find("_B") != std::string::npos )
        {
          if ( tempErr < 0.1 ) myfit->SetParError(1, 2.*tempErr);
          if ( tempErr < 5e-2 ) myfit->SetParError(1, 1.);
        }
      }
      if (! useSingleG && (fitTheta || fitPhi) && xBin < 4 )
      {
        double tempErr = myfit->GetParError(1);
        if (fitTheta)
        {
          //myfit->SetParError(1, 3.*tempErr);
          if ( tempErr < 4e-5 ) myfit->SetParError(1, 0.0005);
        }
        else if (fitPhi)
        {
          if ( myfit->GetParameter(1) > 0.029 ) myfit->SetParError(1, 0.005);
          if ( tempErr < 5e-4 ) myfit->SetParError(1, 0.001);
        }
      }
      if (! useSingleG && fitPhi && xBin == 5 )
      {
        double tempErr = myfit->GetParError(1);
        myfit->SetParError(1, 1.5*tempErr);
      }
      if (! useSingleG && fitEtMu )
      {
        double tempErr = myfit->GetParError(1);
        if ( tempErr < 1e-3 ) myfit->SetParError(1, 0.5);
      }
      
      int npFits = myfit->GetNumberFitPoints();
      if (npFits > nPar && npFits >= cut)
      {
        for (int iPar = 0; iPar < nPar; iPar++)
        {
          //std::cout << "myfit->GetParameter("<<iPar<<") " << myfit->GetParameter(iPar) << std::endl;
          //std::cout << "myfit->GetParError("<<iPar<<") " << myfit->GetParError(iPar) << std::endl;
          hlist[iPar]->SetBinContent(xBin, myfit->GetParameter(iPar)); // fill histogram for parameter i
          //if ( fitTheta && xBin == 1 ) hlist[iPar]->SetBinError(xBin, 3.*myfit->GetParError(iPar));
          //else if ( fitPhi && xBin == 1 ) hlist[iPar]->SetBinError(xBin, 0.0005);
          //else
          hlist[iPar]->SetBinError(xBin,myfit->GetParError(iPar));
        }
        //hchi2->Fill(histo->GetXaxis()->GetBinCenter(binOn),myfit->GetChisquare()/(npfits-npar));
      }
      hp->Write();
      myfit->Write();
      
      if (verbose) std::cout << "Integral of the fit is " << myfit->Integral(-100,100) << std::endl;
      
      delete hp;
      delete myfit;
      
    }  // end loop bins
    
    
    /// Define the fitfunction for all parameters (6):
    //  ai = ai0 + ai1*Ep + ai2*sqrt(Ep)
    //  Its range depends on the parton energy range (hence, the X-axis)
    //TF1 *myfit2 = new TF1("myfit2", "[0]+[1]*x+[2]*sqrt(x)", histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
    double minfit2 = histo->GetXaxis()->GetXmin();
    double maxfit2 = histo->GetXaxis()->GetXmax();
//    if (! useSingleG && fitEt && histoNames[f].find("_B") == std::string::npos )
//      maxfit2 = histo->GetXaxis()->GetBinUpEdge(5) - 1e-5;
//    else if ( fitPhi && histoNames[f].find("_E") != std::string::npos )
//      maxfit2 = histo->GetXaxis()->GetBinUpEdge(4) - 1e-5;
    if (fitEtMu) maxfit2 = histo->GetXaxis()->GetBinUpEdge(3) - 1e-5;
    
    
    TF1 *myfit2;
    if (fitEt || fitEtMu) myfit2 = new TF1("myfit2", "[0]+[1]*x", minfit2, maxfit2);
    else myfit2 = new TF1("myfit2", "[0]+[1]*x+[2]*sqrt(x)", minfit2, maxfit2);
    // Give names to the parameters
    myfit2->SetParName(0,"ai0");
    myfit2->SetParName(1,"ai1");
    if (! fitEt && ! fitEtMu) myfit2->SetParName(2,"ai2");

    for (int iPar = 0; iPar < nPar; iPar++){
      int paramname = iPar+1;
      std::string func_title2 = std::string(histo->GetName())+"_a"+toStr(paramname)+"_Fitted";
      myfit2->SetName(func_title2.c_str());
      hlist[iPar]->Fit(myfit2,"R");
      hlist[iPar]->SetStats(0);
      hlist[iPar]->Write();
      myfit2->Write();
    }
    
    delete [] hlist;
    delete myfit2;
    
  }  /// end loop on histos (f)
  
  
  //writeTable();
}

void ResolutionFunctions::makeFit(std::string inputFileName, std::string outputFileName)
{
  this->makeFit(inputFileName, outputFileName, false);
}

void ResolutionFunctions::makeFit(std::string inputFileName, std::string outputFileName, bool simplify)
{
  if (verbose) std::cout << "ResolutionFunctions::makeFit - Getting histograms from file  " << inputFileName << std::endl;
  getHistos = true;
  useSingleG = simplify;
  
  TFile *fin = new TFile(inputFileName.c_str(),"read");
  fin->cd();
  
  for (int iHisto = 0; iHisto < nHistos; iHisto++)
  {
    if ( fin->GetListOfKeys()->Contains((histoNames[iHisto]).c_str()) )
    {
      fitHisto2D[histoNames[iHisto]] = (TH2F*) fin->Get((histoNames[iHisto]).c_str());
//       if ( iHisto == 12 || iHisto == 13 || iHisto == 14 || iHisto == 15 || iHisto == 16 || iHisto == 17)
//         muon = true;
//       if ( iHisto == 18 || iHisto == 19 || iHisto == 20 || iHisto == 21 || iHisto == 22 || iHisto == 23)
//         electron = true;
    }
  }
  
  
  TFile *foutRF = new TFile(outputFileName.c_str(), "RECREATE");
  foutRF->cd();
  
  for (int f = 0; f < nHistos; f++)
  {
    //if (! muon &&     (f == 12 || f == 13 || f == 14 || f == 15 || f == 16 || f == 17) ) continue;
    //if (! electron && (f == 18 || f == 19 || f == 20 || f == 21 || f == 22 || f == 23) ) continue;
    fitHisto2D[histoNames[f]]->Write();
  }
  
  this->makeFit();
  
  fin->Close();
  foutRF->Close();
  
  delete fin;
  delete foutRF;
}

std::vector<std::array<double, 2> > ResolutionFunctions::getParameters(std::string inputFileName, std::string varName, std::string objName, std::string binName)
{
  int varId = -1, objId = -1, binId = -1;
  if ( varName.std::string::find("Et") != std::string::npos ) varId = 0;
  else if ( varName.std::string::find("theta") != std::string::npos ) varId = 1;
  else if ( varName.std::string::find("phi") != std::string::npos ) varId = 2;
  
  if ( objName.std::string::find("nonbjet") != std::string::npos ) objId = 1;
  else if ( objName.std::string::find("bjet") != std::string::npos ) objId = 0;
  else if ( objName.std::string::find("mu") != std::string::npos ) objId = 2;
  else if ( objName.std::string::find("el") != std::string::npos ) objId = 3;
  
  if ( binName.std::string::find("B") != std::string::npos ) binId = 0;
  else if ( binName.std::string::find("O") != std::string::npos ) binId = 1;
  else if ( binName.std::string::find("E") != std::string::npos ) binId = 2;
  
  if ( varId == -1 || objId == -1 || binId == -1 )
  {
    if ( varId == -1 ) std::cout << "ResolutionFunctions::Variable " << varName << " is unknown. Parameters cannot be found..." << std::endl;
    if ( objId == -1 ) std::cout << "ResolutionFunctions::Object " << objName << " is unknown. Parameters cannot be found..." << std::endl;
    if ( binId == -1 ) std::cout << "ResolutionFunctions::Eta bin " << binName << " is unknown. Parameters cannot be found..." << std::endl;
    return {{0,0}};
  }
  
  int f = 9*objId+3*varId+binId;
  
  if (verbose)
    std::cout << "ResolutionFunctions::Getting resolution function for the " << histoDescription[f] << std::endl;
  
  TF1 *TF_par;
  TFile* rf = new TFile(inputFileName.c_str(),"READ");
  rf->cd();
  
  int nParams = ((TF1*) rf->Get((histoNames[f]+"_sliceXbin1_Fitted").c_str()))->GetNpar();
  if (verbose) std::cout << "ResolutionFunctions::Number of parameters: " << nParams << std::endl;
  
  std::vector<std::array<double, 2> > params;
  params.clear();
  for (int iPar = 0; iPar < nParams; iPar++)
  {
    std::string name = histoNames[f]+"_a"+toStr(iPar+1)+"_Fitted";
    TF_par = (TF1*)rf->Get(name.c_str());
    params.push_back({TF_par->GetParameter(0),TF_par->GetParameter(1)});
  }
  
  rf->Close();
  return params;
  
}

TF2* ResolutionFunctions::getFitFunction2D(std::string inputFileName, std::string varName, std::string objName, std::string binName)
{
  std::vector<std::array<double, 2> > params = this->getParameters(inputFileName, varName, objName, binName);
  
  if ( params.size() < 6 )
  {
    std::cout << "ResolutionFunctions::Using simplified fit. Cannot reconstruct accurate fit function..." << std::endl;
    return NULL;
  }
  
  double funcMin = -80., funcMax = 80.;
  if ( varName.std::string::find("theta") != std::string::npos ) { funcMin = -0.1; funcMax = 0.1; }
  if ( varName.std::string::find("phi") != std::string::npos ) { funcMin = -0.15; funcMax = 0.15; }
  
  TF2 *f2 = new TF2("f2",dblGausParFill,0.,200.,funcMin,funcMax, 12);
  for (int iPar = 0; iPar < 12; iPar++)
  {
    int par = (int) ((double)iPar/2.);
    //if ( iPar%2 != 0 ) par += 1;
    f2->SetParameter(iPar, params[par][iPar%2]);
    //if (verbose) std::cout << "Parameter " << iPar << " set to " << params[par][iPar%2] << std::endl;
  }
  
  return f2;
}

TF1* ResolutionFunctions::getFitFunction1D(std::string inputFileName, std::string varName, std::string objName, std::string binName)
{
  TF2* f2 = this->getFitFunction2D(inputFileName, varName, objName, binName);
  if ( f2 == NULL ) return NULL;
  
  // Make projection on x axis (remove depency on reco/gen difference)
  TF12 *f2x = new TF12("f2x", f2, 0, "x");
  // Convert to independent TF1
  TF1 *f1x = (TF1*) gROOT->GetFunction("f2x");
  
  return f1x; 
}

TF1* ResolutionFunctions::getResolutionFunction(std::string inputFileName, std::string varName, std::string objName, std::string binName)
{
  std::vector<std::array<double, 2> > params = this->getParameters(inputFileName, varName, objName, binName);
  
  TF1 *f = new TF1("f", "[0] + [1]*x", 0., 250.);
  // fill params for narrow gaussian sigma
  f->FixParameter(0, params[1][0]);
  f->FixParameter(1, params[1][1]);
  
  return f; 
}

double ResolutionFunctions::getResolution(std::string inputFileName, std::string varName, std::string objName, double var, std::string binName)
{
  TF1* f = (TF1*) this->getResolutionFunction(inputFileName, varName, objName, binName);
  
  return f->Eval(var);
}

void ResolutionFunctions::writeTable(std::string inputFileName)
{
  if (verbose) std::cout << "ResolutionFunctions::Writing table with parameters..." << std::endl;
  
  std::ofstream myResolutionFunctions;
  std::string myResolutionFunctions_TABLE = "ResolutionFunctions_TABLE.txt";
  myResolutionFunctions.open(myResolutionFunctions_TABLE.c_str());
  
  for(int f = 0; f < nHistos; f++)
  {
    TF1 *TF_par1,*TF_par2,*TF_par3,*TF_par4,*TF_par5,*TF_par6;
    std::string name1 = histoNames[f]+"_a1_Fitted";
    std::string name2 = histoNames[f]+"_a2_Fitted";
    std::string name3 = histoNames[f]+"_a3_Fitted";
    std::string name4 = histoNames[f]+"_a4_Fitted";
    std::string name5 = histoNames[f]+"_a5_Fitted";
    std::string name6 = histoNames[f]+"_a6_Fitted";
    
    TFile* rf = new TFile(inputFileName.c_str(),"READ");
    rf->cd();
    TF_par1 = (TF1*)rf->Get(name1.c_str());
    TF_par2 = (TF1*)rf->Get(name2.c_str());
    TF_par3 = (TF1*)rf->Get(name3.c_str());
    TF_par4 = (TF1*)rf->Get(name4.c_str());
    TF_par5 = (TF1*)rf->Get(name5.c_str());
    TF_par6 = (TF1*)rf->Get(name6.c_str());
    
    /// Write values to table
    if (TF_par1 && TF_par2 && TF_par3 && TF_par4 && TF_par5 && TF_par6)
    {
      myResolutionFunctions<< std::endl;
      myResolutionFunctions<<"\\begin{table}" << std::endl;
      myResolutionFunctions<<"\\caption{Parameters of the resolution function for the " << histoDescription[f]  << "}" << std::endl;
      myResolutionFunctions<<"\\label{tab:}" << std::endl;
      myResolutionFunctions<<"\\centering" << std::endl;
      //myResolutionFunctions<<"\\begin{tabular}{c|ccc}" << std::endl;
      myResolutionFunctions<<"\\begin{tabular}{c|cc}" << std::endl;
      myResolutionFunctions<<"\\hline" << std::endl;
      //myResolutionFunctions << "Type	& $a_{i0}$ & $a_{i1}$ ($\\sqrt{E}$) & $a_{i2}$ ($E$)" << "\\\\" << std::endl;
      myResolutionFunctions << "Type	& $a_{i0}$ & $a_{i1}$ ($E$)" << "\\\\" << std::endl;
      myResolutionFunctions<<"\\hline" << std::endl;
      myResolutionFunctions << "Mean narrow gaussian & $a_{10}$ = " << TF_par1->GetParameter(0) << "$\\pm$" << TF_par1->GetParError(0) << " & $a_{11}$ = " << TF_par1->GetParameter(1) << "$\\pm$" << TF_par1->GetParError(1) /*<< " & $a_{12}$ = " << TF_par1->GetParameter(2) << "$\\pm$" << TF_par1->GetParError(2)*/ << "\\\\" << std::endl;
      myResolutionFunctions << "Width narrow gaussian & $a_{20}$ = " << TF_par2->GetParameter(0) << "$\\pm$" << TF_par2->GetParError(0) << " & $a_{21}$ = " << TF_par2->GetParameter(1) << "$\\pm$" << TF_par2->GetParError(1) /*<< " & $a_{22}$ = " << TF_par2->GetParameter(2) << "$\\pm$" << TF_par2->GetParError(2)*/ << "\\\\" << std::endl;
      myResolutionFunctions << "Scale factor gaussians & $a_{30}$ = " << TF_par3->GetParameter(0) << "$\\pm$" << TF_par3->GetParError(0) << " & $a_{31}$ = " << TF_par3->GetParameter(1) << "$\\pm$" << TF_par3->GetParError(1) /*<< " & $a_{32}$ = " << TF_par3->GetParameter(2) << "$\\pm$" << TF_par3->GetParError(2)*/ << "\\\\" << std::endl;
      myResolutionFunctions << "Mean broad gaussian & $a_{40}$ = " << TF_par4->GetParameter(0) << "$\\pm$" << TF_par4->GetParError(0) << " & $a_{41}$ = " << TF_par4->GetParameter(1) << "$\\pm$" << TF_par4->GetParError(1) /*<< " & $a_{42}$ = " << TF_par4->GetParameter(2) << "$\\pm$" << TF_par4->GetParError(2)*/ << "\\\\" << std::endl;
      myResolutionFunctions << "Width broad gaussian & $a_{50}$ = " << TF_par5->GetParameter(0) << "$\\pm$" << TF_par5->GetParError(0) << " & $a_{51}$ = " << TF_par5->GetParameter(1) << "$\\pm$" << TF_par5->GetParError(1) /*<< " & $a_{52}$ = " << TF_par5->GetParameter(2) << "$\\pm$" << TF_par5->GetParError(2)*/ << "\\\\" << std::endl;
      myResolutionFunctions << "Amplitude & $a_{60}$ = " << TF_par6->GetParameter(0) << "$\\pm$" << TF_par6->GetParError(0) << " & $a_{61}$ = " << TF_par6->GetParameter(1) << "$\\pm$" << TF_par6->GetParError(1) /*<< " & $a_{62}$ = " << TF_par6->GetParameter(2) << "$\\pm$" << TF_par6->GetParError(2)*/ << "\\\\" << std::endl;
      myResolutionFunctions<<"\\hline" << std::endl;
      myResolutionFunctions<<"\\end{tabular}"<< std::endl;
      myResolutionFunctions<<"\\end{table}"<<std::endl;
      myResolutionFunctions<< std::endl;
    }
    rf->Close();
  }
  
  myResolutionFunctions.close();
  
}
