#include <iostream>
#include <string>
#include <cstdlib>
#include <TFile.h>
#include <TObject.h>
#include <TH2.h>
#include <map>

using namespace std;


map<string,TH2F*> histo2D;


void CalculateSystBTagSFs()
{
  string sysType[3] = {"central", "up", "down"};
  string mainSys[3] = {"fsr_up", "fsr_down", "herwigpp"};
  
  string fileNameIn, fileNameOut, fileNameNom;
  string dirNom = "BTagHistos/170711/Merged/";
  string dirIn = "BTagHistos/170731/Merged/";
  
  /// Only calculate for central SFs
  int iSys = 0;
  
  fileNameNom = "BTagSFs_TT_nominal_mujets_"+sysType[iSys];
  TFile* fnom = new TFile((dirNom+fileNameNom+".root").c_str(),"READ");
  fnom->cd();
  
  histo2D["BtaggedJets_nominal_"+sysType[iSys]]       = (TH2F*) fnom->Get(("BtaggedJets_"+sysType[iSys]).c_str());
  histo2D["BtaggedBJets_nominal_"+sysType[iSys]]      = (TH2F*) fnom->Get(("BtaggedBJets_"+sysType[iSys]).c_str());
  histo2D["BtaggedCJets_nominal_"+sysType[iSys]]      = (TH2F*) fnom->Get(("BtaggedCJets_"+sysType[iSys]).c_str());
  histo2D["BtaggedLightJets_nominal_"+sysType[iSys]]  = (TH2F*) fnom->Get(("BtaggedLightJets_"+sysType[iSys]).c_str());
  histo2D["TotalNofBJets_nominal_"+sysType[iSys]]     = (TH2F*) fnom->Get(("TotalNofBJets_"+sysType[iSys]).c_str());
  histo2D["TotalNofCJets_nominal_"+sysType[iSys]]     = (TH2F*) fnom->Get(("TotalNofCJets_"+sysType[iSys]).c_str());
  histo2D["TotalNofLightJets_nominal_"+sysType[iSys]] = (TH2F*) fnom->Get(("TotalNofLightJets_"+sysType[iSys]).c_str());
  
  for (int i = 0; i < 3; i++)
  {
    fileNameIn = "BTagSFs_TT_"+mainSys[i]+"_mujets_"+sysType[iSys];
    fileNameOut = fileNameIn+"_diff";
    
    TFile* fin = new TFile((dirIn+fileNameIn+".root").c_str(),"READ");
    fin->cd();
    
    histo2D["BtaggedJets_"+mainSys[i]+"_"+sysType[iSys]]       = (TH2F*) fin->Get(("BtaggedJets_"+sysType[iSys]).c_str());
    histo2D["BtaggedBJets_"+mainSys[i]+"_"+sysType[iSys]]      = (TH2F*) fin->Get(("BtaggedBJets_"+sysType[iSys]).c_str());
    histo2D["BtaggedCJets_"+mainSys[i]+"_"+sysType[iSys]]      = (TH2F*) fin->Get(("BtaggedCJets_"+sysType[iSys]).c_str());
    histo2D["BtaggedLightJets_"+mainSys[i]+"_"+sysType[iSys]]  = (TH2F*) fin->Get(("BtaggedLightJets_"+sysType[iSys]).c_str());
    histo2D["TotalNofBJets_"+mainSys[i]+"_"+sysType[iSys]]     = (TH2F*) fin->Get(("TotalNofBJets_"+sysType[iSys]).c_str());
    histo2D["TotalNofCJets_"+mainSys[i]+"_"+sysType[iSys]]     = (TH2F*) fin->Get(("TotalNofCJets_"+sysType[iSys]).c_str());
    histo2D["TotalNofLightJets_"+mainSys[i]+"_"+sysType[iSys]] = (TH2F*) fin->Get(("TotalNofLightJets_"+sysType[iSys]).c_str());
    
    
    /// Calculate efficiencies
    int nBinsX = histo2D["BtaggedJets_nominal_"+sysType[iSys]]->GetNbinsX();
    int nBinsY = histo2D["BtaggedJets_nominal_"+sysType[iSys]]->GetNbinsY();
    float minX = histo2D["BtaggedJets_nominal_"+sysType[iSys]]->GetXaxis()->GetBinLowEdge(1);
    float maxX = histo2D["BtaggedJets_nominal_"+sysType[iSys]]->GetXaxis()->GetBinUpEdge(nBinsX);
    float minY = histo2D["BtaggedJets_nominal_"+sysType[iSys]]->GetYaxis()->GetBinLowEdge(1);
    float maxY = histo2D["BtaggedJets_nominal_"+sysType[iSys]]->GetYaxis()->GetBinUpEdge(nBinsY);
    
    histo2D["effBJets_nominal_"+sysType[iSys]] = new TH2F(("effBJets_nominal_"+sysType[iSys]).c_str(), ("effBJets_nominal_"+sysType[iSys]+"; p_{T}; |#eta|").c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    histo2D["effCJets_nominal_"+sysType[iSys]] = new TH2F(("effCJets_nominal_"+sysType[iSys]).c_str(), ("effCJets_nominal_"+sysType[iSys]+"; p_{T}; |#eta|").c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    histo2D["effLightJets_nominal_"+sysType[iSys]] = new TH2F(("effLightJets_nominal_"+sysType[iSys]).c_str(), ("effLightJets_nominal_"+sysType[iSys]+"; p_{T}; |#eta|").c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);

    histo2D["effBJets_"+mainSys[i]+"_"+sysType[iSys]] = new TH2F(("effBJets_"+mainSys[i]+"_"+sysType[iSys]).c_str(), ("effBJets_"+mainSys[i]+"_"+sysType[iSys]+"; p_{T}; |#eta|").c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    histo2D["effCJets_"+mainSys[i]+"_"+sysType[iSys]] = new TH2F(("effCJets_"+mainSys[i]+"_"+sysType[iSys]).c_str(), ("effCJets_"+mainSys[i]+"_"+sysType[iSys]+"; p_{T}; |#eta|").c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    histo2D["effLightJets_"+mainSys[i]+"_"+sysType[iSys]] = new TH2F(("effLightJets_"+mainSys[i]+"_"+sysType[iSys]).c_str(), ("effLightJets_"+mainSys[i]+"_"+sysType[iSys]+"; p_{T}; |#eta|").c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    
    histo2D["relEffBJets_"+mainSys[i]+"_"+sysType[iSys]] = new TH2F(("relEffBJets_"+mainSys[i]+"_"+sysType[iSys]).c_str(), ("relEffBJets_"+mainSys[i]+"_"+sysType[iSys]+"; p_{T}; |#eta|").c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    histo2D["relEffCJets_"+mainSys[i]+"_"+sysType[iSys]] = new TH2F(("relEffCJets_"+mainSys[i]+"_"+sysType[iSys]).c_str(), ("relEffCJets_"+mainSys[i]+"_"+sysType[iSys]+"; p_{T}; |#eta|").c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    histo2D["relEffLightJets_"+mainSys[i]+"_"+sysType[iSys]] = new TH2F(("relEffLightJets_"+mainSys[i]+"_"+sysType[iSys]).c_str(), ("relEffLightJets_"+mainSys[i]+"_"+sysType[iSys]+"; p_{T}; |#eta|").c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    
    float num, denom;
    for (int xBin = 1; xBin < nBinsX+1; xBin++)
    {
      for (int yBin = 1; yBin < nBinsY+1; yBin++)
      {
        /// b jets
        denom = histo2D["TotalNofBJets_nominal_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        num = histo2D["BtaggedBJets_nominal_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        if (denom == 0.) histo2D["effBJets_nominal_"+sysType[iSys]]->SetBinContent(xBin, yBin, 0.);
        else histo2D["effBJets_nominal_"+sysType[iSys]]->SetBinContent(xBin, yBin, num/denom);
        denom = histo2D["TotalNofBJets_"+mainSys[i]+"_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        num = histo2D["BtaggedBJets_"+mainSys[i]+"_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        if (denom == 0.) histo2D["effBJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, 0.);
        else histo2D["effBJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, num/denom);
        
        /// c jets
        denom = histo2D["TotalNofCJets_nominal_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        num = histo2D["BtaggedCJets_nominal_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        if (denom == 0.) histo2D["effCJets_nominal_"+sysType[iSys]]->SetBinContent(xBin, yBin, 0.);
        else histo2D["effCJets_nominal_"+sysType[iSys]]->SetBinContent(xBin, yBin, num/denom);
        denom = histo2D["TotalNofCJets_"+mainSys[i]+"_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        num = histo2D["BtaggedCJets_"+mainSys[i]+"_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        if (denom == 0.) histo2D["effCJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, 0.);
        else histo2D["effCJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, num/denom);
        
        /// light jets
        denom = histo2D["TotalNofLightJets_nominal_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        num = histo2D["BtaggedLightJets_nominal_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        if (denom == 0.) histo2D["effLightJets_nominal_"+sysType[iSys]]->SetBinContent(xBin, yBin, 0.);
        else histo2D["effLightJets_nominal_"+sysType[iSys]]->SetBinContent(xBin, yBin, num/denom);
        denom = histo2D["TotalNofLightJets_"+mainSys[i]+"_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        num = histo2D["BtaggedLightJets_"+mainSys[i]+"_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        if (denom == 0.) histo2D["effLightJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, 0.);
        else histo2D["effLightJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, num/denom);
        
        /// Calculate relative efficiency: eff_nom/eff_var
        denom = histo2D["effBJets_"+mainSys[i]+"_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        num = histo2D["effBJets_nominal_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        if (denom == 0.) histo2D["relEffBJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, 0.);
        else histo2D["relEffBJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, num/denom);
        
        denom = histo2D["effCJets_"+mainSys[i]+"_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        num = histo2D["effCJets_nominal_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        if (denom == 0.) histo2D["relEffCJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, 0.);
        else histo2D["relEffCJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, num/denom);
        
        denom = histo2D["effLightJets_"+mainSys[i]+"_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        num = histo2D["effLightJets_nominal_"+sysType[iSys]]->GetBinContent(xBin, yBin);
        if (denom == 0.) histo2D["relEffLightJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, 0.);
        else histo2D["relEffLightJets_"+mainSys[i]+"_"+sysType[iSys]]->SetBinContent(xBin, yBin, num/denom);
      }
    }      
    
    TFile* fout = new TFile((fileNameOut+".root").c_str(),"RECREATE");
    fout->cd();
    
    histo2D["effBJets_nominal_"+sysType[iSys]]->SetStats(0);
    histo2D["effCJets_nominal_"+sysType[iSys]]->SetStats(0);
    histo2D["effLightJets_nominal_"+sysType[iSys]]->SetStats(0);
    
    histo2D["effBJets_"+mainSys[i]+"_"+sysType[iSys]]->SetStats(0);
    histo2D["effCJets_"+mainSys[i]+"_"+sysType[iSys]]->SetStats(0);
    histo2D["effLightJets_"+mainSys[i]+"_"+sysType[iSys]]->SetStats(0);
    
    histo2D["relEffBJets_"+mainSys[i]+"_"+sysType[iSys]]->SetStats(0);
    histo2D["relEffCJets_"+mainSys[i]+"_"+sysType[iSys]]->SetStats(0);
    histo2D["relEffLightJets_"+mainSys[i]+"_"+sysType[iSys]]->SetStats(0);
    
    histo2D["effBJets_nominal_"+sysType[iSys]]->Write();
    histo2D["effCJets_nominal_"+sysType[iSys]]->Write();
    histo2D["effLightJets_nominal_"+sysType[iSys]]->Write();
    
    histo2D["effBJets_"+mainSys[i]+"_"+sysType[iSys]]->Write();
    histo2D["effCJets_"+mainSys[i]+"_"+sysType[iSys]]->Write();
    histo2D["effLightJets_"+mainSys[i]+"_"+sysType[iSys]]->Write();
    
    histo2D["relEffBJets_"+mainSys[i]+"_"+sysType[iSys]]->Write();
    histo2D["relEffCJets_"+mainSys[i]+"_"+sysType[iSys]]->Write();
    histo2D["relEffLightJets_"+mainSys[i]+"_"+sysType[iSys]]->Write();
    
    fout->Close();
    fin->Close();
    
    delete fin;
    delete fout;
  }
  
  fnom->Close();
  
  delete fnom;
}
