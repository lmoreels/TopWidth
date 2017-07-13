#include <iostream>
#include <string>
#include <cstdlib>
#include <TFile.h>
#include <TObject.h>
#include <TH2.h>
#include <map>

using namespace std;


map<string,TH2F*> histo;


void RenameHistograms()
{
  string sysType = "down";  /// choose between "central", "up" and "down"
  
  string fileNameIn = "PlotsForBTagSFs_nominal_"+sysType;
  string fileNameOut = fileNameIn+"_test";
  
  TFile* fin = new TFile((fileNameIn+".root").c_str(),"READ");
  
  TH2F* AllBJets = (TH2F*)fin->Get("BtaggedJets");
  TH2F* BBJets = (TH2F*)fin->Get("BtaggedBJets");
  TH2F* CBJets = (TH2F*)fin->Get("BtaggedCJets");
  TH2F* LightBJets = (TH2F*)fin->Get("BtaggedLightJets");
  TH2F* TotBJets = (TH2F*)fin->Get("TotalNofBJets");
  TH2F* TotCJets = (TH2F*)fin->Get("TotalNofCJets");
  TH2F* TotLightJets = (TH2F*)fin->Get("TotalNofLightJets");
  
  TFile* fout = new TFile((fileNameOut+".root").c_str(),"RECREATE");
  
  histo["BtaggedJets_"+sysType]       = (TH2F*) AllBJets->Clone(("BtaggedJets_"+sysType).c_str());
  histo["BtaggedBJets_"+sysType]      = (TH2F*) BBJets->Clone(("BtaggedBJets_"+sysType).c_str());
  histo["BtaggedCJets_"+sysType]      = (TH2F*) CBJets->Clone(("BtaggedCJets_"+sysType).c_str());
  histo["BtaggedLightJets_"+sysType]  = (TH2F*) LightBJets->Clone(("BtaggedLightJets_"+sysType).c_str());
  histo["TotalNofBJets_"+sysType]     = (TH2F*) TotBJets->Clone(("TotalNofBJets_"+sysType).c_str());  
  histo["TotalNofCJets_"+sysType]     = (TH2F*) TotCJets->Clone(("TotalNofCJets_"+sysType).c_str());  
  histo["TotalNofLightJets_"+sysType] = (TH2F*) TotLightJets->Clone(("TotalNofLightJets_"+sysType).c_str());
  
  
  for(std::map<std::string,TH2F*>::const_iterator it = histo.begin(); it != histo.end(); it++)
  {
    TH2F *temp = it->second;
    temp->Write();
  }
  
  fout->Close();
  fin->Close();
  
  delete fin;
  delete fout;
}
