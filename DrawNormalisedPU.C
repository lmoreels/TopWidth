//#include <stdio.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TLine.h>
#include <map>
#include <TArray.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TColor.h>


using namespace std;

void DrawNormalisedPU()
{
  map<string,TH1F*> histo;
  
  string pathFiles = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopTreeAnalysisBase/Calibrations/PileUpReweighting/";
  string outputPath = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/";
  string fileNameData = "pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet.root";
  string fileNameMC = "pileup_MC_RunIISpring16MiniAODv2-Asympt.root";
  
  TFile *fileInD = new TFile((pathFiles+fileNameData).c_str(),"read");
  TFile *fileInMC = new TFile((pathFiles+fileNameMC).c_str(),"read");
  TH1F* data = (TH1F*) fileInD->Get("pileup");
  TH1F* MC = (TH1F*) fileInMC->Get("pileup");
  
  TFile *fileOut = new TFile((outputPath+"PU_normalised_plots.root").c_str(),"RECREATE");
  fileOut->cd();
  
  TCanvas *c1 = new TCanvas("c1","Normalised plot");
  c1->cd();
    
  /// Make legend
  TLegend *leg = new TLegend(0.65,0.56,0.9,0.9);
  
  data->SetLineWidth(2);
  data->SetLineColor(kBlack);
  data->SetTitle("Comparison of number of vertices (normalised)");
  data->GetXaxis()->SetTitle("# vertices");
  data->SetStats(0);
  data->DrawNormalized();
  MC->SetLineWidth(2);
  MC->SetLineColor(kRed);
  MC->DrawNormalized("same");
  leg->AddEntry(data,"data","l");
  leg->AddEntry(MC,"MC","l");
  leg->Draw();
  c1->Update();
  c1->Write();
  c1->SaveAs("nVtx_comparison_normalised.png");
  
  data->Scale(1./data->Integral());
  MC->Scale(1./MC->Integral());
  
  TH1F* ratio = (TH1F*)data->Clone("ratio");
  
  int nBinsX = data->GetNbinsX();
  
  float binContentData, binContentMC, binContentNew;
  for (int iBin = 0; iBin < nBinsX+1; iBin++)
  {
    binContentData = data->GetBinContent(iBin);
    binContentMC = MC->GetBinContent(iBin);
    if ( binContentData == 0. ) binContentNew = 0.;
    else binContentNew = binContentData/binContentMC;
    ratio->SetBinContent(iBin, binContentNew);
    ratio->SetEntries(ratio->GetEntries()-1);
  }
  
  ratio->Write();
  
  
  
  fileOut->Close();
  fileInD->Close();
  fileInMC->Close();

}
