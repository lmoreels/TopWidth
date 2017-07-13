#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <TFile.h>
#include <TH2.h>

using namespace std;

map<string,TH2F*> histo2D;

void MakeEfficiencyPlots()
{
  string sysType = "central";
  string histoFileName = "PlotsForBTagSFs_nominal_"+sysType+".root";
  TFile* f_ = new TFile(histoFileName.c_str(),"READ");
  f_->cd();
  histo2D["BtaggedJets"] = (TH2F*) f_->Get(("BtaggedJets_"+sysType).c_str());
  histo2D["BtaggedBJets"] = (TH2F*) f_->Get(("BtaggedBJets_"+sysType).c_str());
  histo2D["BtaggedCJets"] = (TH2F*) f_->Get(("BtaggedCJets_"+sysType).c_str());
  histo2D["BtaggedLightJets"] = (TH2F*) f_->Get(("BtaggedLightJets_"+sysType).c_str());
  histo2D["TotalNofBJets"] = (TH2F*) f_->Get(("TotalNofBJets_"+sysType).c_str());
  histo2D["TotalNofCJets"] = (TH2F*) f_->Get(("TotalNofCJets_"+sysType).c_str());
  histo2D["TotalNofLightJets"] = (TH2F*) f_->Get(("TotalNofLightJets_"+sysType).c_str());
  
  int nBinsX = histo2D["BtaggedJets"]->GetNbinsX();
  int nBinsY = histo2D["BtaggedJets"]->GetNbinsY();
  float minX = histo2D["BtaggedJets"]->GetXaxis()->GetBinLowEdge(1);
  float maxX = histo2D["BtaggedJets"]->GetXaxis()->GetBinUpEdge(nBinsX);
  float minY = histo2D["BtaggedJets"]->GetYaxis()->GetBinLowEdge(1);
  float maxY = histo2D["BtaggedJets"]->GetYaxis()->GetBinUpEdge(nBinsY);
  
  histo2D["effBJets"] = new TH2F("effBJets", "effBJets", nBinsX, minX, maxX, nBinsY, minY, maxY);
  histo2D["effCJets"] = new TH2F("effCJets", "effCJets", nBinsX, minX, maxX, nBinsY, minY, maxY);
  histo2D["effLightJets"] = new TH2F("effLightJets", "effLightJets", nBinsX, minX, maxX, nBinsY, minY, maxY);
  
  float num, denom;
  for (int xBin = 1; xBin < nBinsX+1; xBin++)
  {
    for (int yBin = 1; yBin < nBinsY+1; yBin++)
    {
      /// b jets
      denom = histo2D["TotalNofBJets"]->GetBinContent(xBin, yBin);
      num = histo2D["BtaggedBJets"]->GetBinContent(xBin, yBin);
      if (denom == 0.) histo2D["effBJets"]->SetBinContent(xBin, yBin, 0.);
      else histo2D["effBJets"]->SetBinContent(xBin, yBin, num/denom);
      
      /// c jets
      denom = histo2D["TotalNofCJets"]->GetBinContent(xBin, yBin);
      num = histo2D["BtaggedCJets"]->GetBinContent(xBin, yBin);
      if (denom == 0.) histo2D["effCJets"]->SetBinContent(xBin, yBin, 0.);
      else histo2D["effCJets"]->SetBinContent(xBin, yBin, num/denom);
      
      /// light jets
      denom = histo2D["TotalNofLightJets"]->GetBinContent(xBin, yBin);
      num = histo2D["BtaggedLightJets"]->GetBinContent(xBin, yBin);
      if (denom == 0.) histo2D["effLightJets"]->SetBinContent(xBin, yBin, 0.);
      else histo2D["effLightJets"]->SetBinContent(xBin, yBin, num/denom);
      
    }
  }
  
  TFile* fout = new TFile("BTagEfficiencies.root","RECREATE");
  fout->cd();
  
  histo2D["effBJets"]->SetStats(0);
  histo2D["effCJets"]->SetStats(0);
  histo2D["effLightJets"]->SetStats(0);
  
  histo2D["effBJets"]->Write();
  histo2D["effCJets"]->Write();
  histo2D["effLightJets"]->Write();
  
  fout->Write();
  fout->Close();
  f_->Close();
  
  delete f_;
  delete fout;
}

