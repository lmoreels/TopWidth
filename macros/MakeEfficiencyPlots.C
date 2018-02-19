#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <map>

using namespace std;

map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

string ConvertFloatToString(float Number)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

void MakeEfficiencyPlots()
{
  //string sysType = "central";
  string histoFileName = "input/PlotsForBTagSFs_nominal.root ";
  TFile* f_ = new TFile(histoFileName.c_str(),"READ");
  f_->cd();
  histo2D["BtaggedJets"]       = (TH2F*) f_->Get("BtaggedJets_");
  histo2D["BtaggedBJets"]      = (TH2F*) f_->Get("BtaggedBJets_");
  histo2D["BtaggedCJets"]      = (TH2F*) f_->Get("BtaggedCJets_");
  histo2D["BtaggedLightJets"]  = (TH2F*) f_->Get("BtaggedLightJets_");
  histo2D["TotalNofBJets"]     = (TH2F*) f_->Get("TotalNofBJets_");
  histo2D["TotalNofCJets"]     = (TH2F*) f_->Get("TotalNofCJets_");
  histo2D["TotalNofLightJets"] = (TH2F*) f_->Get("TotalNofLightJets_");
  
  const int nBinsX = histo2D["BtaggedJets"]->GetNbinsX();
  const int nBinsY = histo2D["BtaggedJets"]->GetNbinsY();
  float minX = histo2D["BtaggedJets"]->GetXaxis()->GetBinLowEdge(1);
  float maxX = histo2D["BtaggedJets"]->GetXaxis()->GetBinUpEdge(nBinsX);
  float minY = histo2D["BtaggedJets"]->GetYaxis()->GetBinLowEdge(1);
  float maxY = histo2D["BtaggedJets"]->GetYaxis()->GetBinUpEdge(nBinsY);
  
  float binWidthX = (maxX - minX)/(float)nBinsX;
  binWidthX = floor(binWidthX*1000.+0.5)/1000.;
  float binWidthY = (maxY - minY)/(float)nBinsY;
  binWidthY = floor(binWidthY*1000.+0.5)/1000.;
  
  histo2D["effBJets"] = new TH2F("effBJets", "effBJets; p_{T} (GeV); |#eta|", nBinsX, minX, maxX, nBinsY, minY, maxY);
  histo2D["effCJets"] = new TH2F("effCJets", "effCJets; p_{T} (GeV); |#eta|", nBinsX, minX, maxX, nBinsY, minY, maxY);
  histo2D["effLightJets"] = new TH2F("effLightJets", "effLightJets; p_{T} (GeV); |#eta|", nBinsX, minX, maxX, nBinsY, minY, maxY);
  
  histo1D["eff_pT_B"]      = new TH1F("eff_pT_B",      ("eff_pT_B; p_{T} (GeV); efficiency / "+ConvertFloatToString(binWidthX)+" GeV").c_str(), nBinsX, minX, maxX);
  //histo1D["eff_pT_B_test"]      = new TH1F("eff_pT_B_test",      "eff_pT_B;     p_{T}; #epsilon", nBinsX, minX, maxX);
  histo1D["eff_pT_C"]      = new TH1F("eff_pT_C",      ("eff_pT_C; p_{T} (GeV); efficiency / "+ConvertFloatToString(binWidthX)+" GeV").c_str(), nBinsX, minX, maxX);
  histo1D["eff_pT_Light"]  = new TH1F("eff_pT_Light",  ("eff_pT_Light; p_{T} (GeV); efficiency / "+ConvertFloatToString(binWidthX)+" GeV").c_str(), nBinsX, minX, maxX);
  histo1D["eff_eta_B"]     = new TH1F("eff_eta_B",     ("eff_eta_B; |#eta|; efficiency / "+ConvertFloatToString(binWidthY)+" units").c_str(), nBinsY, minY, maxY);
  histo1D["eff_eta_C"]     = new TH1F("eff_eta_C",     ("eff_eta_C; |#eta|; efficiency / "+ConvertFloatToString(binWidthY)+" units").c_str(), nBinsY, minY, maxY);
  histo1D["eff_eta_Light"] = new TH1F("eff_eta_Light", ("eff_eta_Light; |#eta|; efficiency / "+ConvertFloatToString(binWidthY)+" units").c_str(), nBinsY, minY, maxY);
  
  float num, denom;
  float sumNumPT[3] = {0.}, sumDenomPT[3] = {0.};
  float sumNumEta[3][nBinsY], sumDenomEta[3][nBinsY];
  for (int yBin = 1; yBin < nBinsY+1; yBin++)
  {
    for (int i = 0; i < 3; i++)
    {
      sumNumEta[i][yBin-1] = 0.;
      sumDenomEta[i][yBin-1] = 0.;
    }
  } 
  
  for (int xBin = 1; xBin < nBinsX+1; xBin++)
  {
    for (int i = 0; i < 3; i++)
    {
      sumNumPT[i] = 0.;
      sumDenomPT[i] = 0.;
    }
    for (int yBin = 1; yBin < nBinsY+1; yBin++)
    {
      /// b jets
      denom = histo2D["TotalNofBJets"]->GetBinContent(xBin, yBin);
      num = histo2D["BtaggedBJets"]->GetBinContent(xBin, yBin);
      if (denom == 0.) histo2D["effBJets"]->SetBinContent(xBin, yBin, 0.);
      else histo2D["effBJets"]->SetBinContent(xBin, yBin, num/denom);
      sumNumPT[0] += num;
      sumDenomPT[0] += denom;
      sumNumEta[0][yBin-1] += num;
      sumDenomEta[0][yBin-1] += denom;
      
      /// c jets
      denom = histo2D["TotalNofCJets"]->GetBinContent(xBin, yBin);
      num = histo2D["BtaggedCJets"]->GetBinContent(xBin, yBin);
      if (denom == 0.) histo2D["effCJets"]->SetBinContent(xBin, yBin, 0.);
      else histo2D["effCJets"]->SetBinContent(xBin, yBin, num/denom);
      sumNumPT[1] += num;
      sumDenomPT[1] += denom;
      sumNumEta[1][yBin-1] += num;
      sumDenomEta[1][yBin-1] += denom;
      
      /// light jets
      denom = histo2D["TotalNofLightJets"]->GetBinContent(xBin, yBin);
      num = histo2D["BtaggedLightJets"]->GetBinContent(xBin, yBin);
      if (denom == 0.) histo2D["effLightJets"]->SetBinContent(xBin, yBin, 0.);
      else histo2D["effLightJets"]->SetBinContent(xBin, yBin, num/denom);
      sumNumPT[2] += num;
      sumDenomPT[2] += denom;
      sumNumEta[2][yBin-1] += num;
      sumDenomEta[2][yBin-1] += denom;
      
    }
    histo1D["eff_pT_B"]->SetBinContent(xBin, sumNumPT[0]/sumDenomPT[0]);
    histo1D["eff_pT_C"]->SetBinContent(xBin, sumNumPT[1]/sumDenomPT[1]);
    histo1D["eff_pT_Light"]->SetBinContent(xBin, sumNumPT[2]/sumDenomPT[2]);
    histo1D["eff_pT_B"]->SetBinError(xBin, 0.001);
    histo1D["eff_pT_C"]->SetBinError(xBin, 0.001);
    histo1D["eff_pT_Light"]->SetBinError(xBin, 0.001);
  }
  for (int yBin = 1; yBin < nBinsY+1; yBin++)
  {
    histo1D["eff_eta_B"]->SetBinContent(yBin, sumNumEta[0][yBin-1]/sumDenomEta[0][yBin-1]);
    histo1D["eff_eta_C"]->SetBinContent(yBin, sumNumEta[1][yBin-1]/sumDenomEta[1][yBin-1]);
    histo1D["eff_eta_Light"]->SetBinContent(yBin, sumNumEta[2][yBin-1]/sumDenomEta[2][yBin-1]);
    histo1D["eff_eta_B"]->SetBinError(yBin, 0.001);
    histo1D["eff_eta_C"]->SetBinError(yBin, 0.001);
    histo1D["eff_eta_Light"]->SetBinError(yBin, 0.001);
  }
  
//   float sumPT[3][nBinsX], sumEta[3][nBinsY];
//   for (int xBin = 1; xBin < nBinsX+1; xBin++)
//   {
//     sumPT[0][xBin-1] = 0.;
//     sumPT[1][xBin-1] = 0.;
//     sumPT[2][xBin-1] = 0.;
//   }
//   for (int yBin = 1; yBin < nBinsY+1; yBin++)
//   {
//     sumEta[0][yBin-1] = 0.;
//     sumEta[1][yBin-1] = 0.;
//     sumEta[2][yBin-1] = 0.;
//   } 
//   for (int xBin = 1; xBin < nBinsX+1; xBin++)
//   {
//     for (int yBin = 1; yBin < nBinsY+1; yBin++)
//     {
//       sumPT[0][xBin-1] += histo2D["effBJets"]->GetBinContent(xBin, yBin);
//       sumPT[1][xBin-1] += histo2D["effCJets"]->GetBinContent(xBin, yBin);
//       sumPT[2][xBin-1] += histo2D["effLightJets"]->GetBinContent(xBin, yBin);
//       
//       sumEta[0][yBin-1] += histo2D["effBJets"]->GetBinContent(xBin, yBin);
//       sumEta[1][yBin-1] += histo2D["effCJets"]->GetBinContent(xBin, yBin);
//       sumEta[2][yBin-1] += histo2D["effLightJets"]->GetBinContent(xBin, yBin);
//     }
//   }
  
//   for (int xBin = 1; xBin < nBinsX+1; xBin++)
//   {
//     histo1D["eff_pT_B"]->SetBinContent(xBin, sumPT[0][xBin-1]/nBinsY);
//     histo1D["eff_pT_C"]->SetBinContent(xBin, sumPT[1][xBin-1]/nBinsY);
//     histo1D["eff_pT_Light"]->SetBinContent(xBin, sumPT[2][xBin-1]/nBinsY);
//   }
//   for (int yBin = 1; yBin < nBinsY+1; yBin++)
//   {
//     histo1D["eff_eta_B"]->SetBinContent(yBin, sumEta[0][yBin-1]/nBinsX);
//     histo1D["eff_eta_C"]->SetBinContent(yBin, sumEta[1][yBin-1]/nBinsX);
//     histo1D["eff_eta_Light"]->SetBinContent(yBin, sumEta[2][yBin-1]/nBinsX);
//   }
  
  TFile* fout = new TFile("BTagEfficiencies.root","RECREATE");
  fout->cd();
  
  histo2D["effBJets"]->SetStats(0);
  histo2D["effCJets"]->SetStats(0);
  histo2D["effLightJets"]->SetStats(0);
  
  histo1D["eff_pT_B"]->SetStats(0);
  histo1D["eff_pT_C"]->SetStats(0);
  histo1D["eff_pT_Light"]->SetStats(0);
  histo1D["eff_eta_B"]->SetStats(0);
  histo1D["eff_eta_C"]->SetStats(0);
  histo1D["eff_eta_Light"]->SetStats(0);
  //histo1D["eff_pT_B_test"]->SetStats(0);
  
  histo2D["effBJets"]->Write();
  histo2D["effCJets"]->Write();
  histo2D["effLightJets"]->Write();
  
  histo1D["eff_pT_B"]->Write();
  histo1D["eff_pT_C"]->Write();
  histo1D["eff_pT_Light"]->Write();
  histo1D["eff_eta_B"]->Write();
  histo1D["eff_eta_C"]->Write();
  histo1D["eff_eta_Light"]->Write();
  //histo1D["eff_pT_B_test"]->Write();
  
  double labelsize = 0.05;
  TCanvas *cp = new TCanvas("Efficiency p_{T}","Efficiency p_{T}");
  cp->cd();
  cp->SetTopMargin(0.05);
  cp->SetBottomMargin(0.15);
  cp->SetRightMargin(0.05);
  cp->SetLeftMargin(0.12);
  histo1D["eff_pT_B"]->SetTitle("");
  histo1D["eff_pT_B"]->SetMarkerStyle(24);
  histo1D["eff_pT_B"]->SetMarkerColor(4);
  histo1D["eff_pT_B"]->SetLineColor(4);
  histo1D["eff_pT_B"]->SetLineWidth(2);
  histo1D["eff_pT_B"]->GetYaxis()->SetRangeUser(0.,1.);
  histo1D["eff_pT_B"]->Draw("P");
  histo1D["eff_pT_B"]->GetXaxis()->SetLabelSize(labelsize);
  histo1D["eff_pT_B"]->GetXaxis()->SetTitleSize(labelsize);
  histo1D["eff_pT_B"]->GetXaxis()->SetTitleOffset(1.1);
  histo1D["eff_pT_B"]->GetYaxis()->SetLabelSize(labelsize);
  histo1D["eff_pT_B"]->GetYaxis()->SetTitleSize(labelsize);
  histo1D["eff_pT_B"]->GetYaxis()->SetTitleOffset(1.05);
  histo1D["eff_pT_C"]->SetMarkerStyle(26);
  histo1D["eff_pT_C"]->SetMarkerColor(418);
  histo1D["eff_pT_C"]->SetLineColor(418);
  histo1D["eff_pT_C"]->SetLineWidth(2);
  histo1D["eff_pT_C"]->Draw("P same");
  histo1D["eff_pT_Light"]->SetMarkerStyle(32);
  histo1D["eff_pT_Light"]->SetMarkerColor(2);
  histo1D["eff_pT_Light"]->SetLineColor(2);
  histo1D["eff_pT_Light"]->SetLineWidth(2);
  histo1D["eff_pT_Light"]->Draw("P same");
  cp->Update();
  TLegend *lp = new TLegend(0.70,0.73,0.94,0.93);
  lp->SetBorderSize(0);
  lp->SetFillStyle(0);
  lp->AddEntry(histo1D["eff_pT_B"],"b jets","lp");
  lp->AddEntry(histo1D["eff_pT_C"],"c jets","lp");
  lp->AddEntry(histo1D["eff_pT_Light"],"light-flavour jets","lp");
  lp->Draw();
  cp->Update();
  cp->Write();
  cp->SaveAs("BTagEfficiencies_pT.png");
  
  TCanvas *ce = new TCanvas("Efficiency #eta","Efficiency #eta");
  ce->cd();
  ce->SetTopMargin(0.05);
  ce->SetBottomMargin(0.15);
  ce->SetRightMargin(0.05);
  ce->SetLeftMargin(0.12);
  histo1D["eff_eta_B"]->SetTitle("");
  histo1D["eff_eta_B"]->SetMarkerStyle(24);
  histo1D["eff_eta_B"]->SetMarkerColor(4);
  histo1D["eff_eta_B"]->SetLineColor(4);
  histo1D["eff_eta_B"]->SetLineWidth(2);
  histo1D["eff_eta_B"]->GetYaxis()->SetRangeUser(0.,1.);
  histo1D["eff_eta_B"]->Draw("P");
  histo1D["eff_eta_B"]->GetXaxis()->SetLabelSize(labelsize);
  histo1D["eff_eta_B"]->GetXaxis()->SetTitleSize(labelsize);
  histo1D["eff_eta_B"]->GetXaxis()->SetTitleOffset(1.05);
  histo1D["eff_eta_B"]->GetYaxis()->SetLabelSize(labelsize);
  histo1D["eff_eta_B"]->GetYaxis()->SetTitleSize(labelsize);
  histo1D["eff_eta_B"]->GetYaxis()->SetTitleOffset(1.05);
  histo1D["eff_eta_C"]->SetMarkerStyle(26);
  histo1D["eff_eta_C"]->SetMarkerColor(418);
  histo1D["eff_eta_C"]->SetLineColor(418);
  histo1D["eff_eta_C"]->SetLineWidth(2);
  histo1D["eff_eta_C"]->Draw("P same");
  histo1D["eff_eta_Light"]->SetMarkerStyle(32);
  histo1D["eff_eta_Light"]->SetMarkerColor(2);
  histo1D["eff_eta_Light"]->SetLineColor(2);
  histo1D["eff_eta_Light"]->SetLineWidth(2);
  histo1D["eff_eta_Light"]->Draw("P same");
  ce->Update();
  //TLegend *le = new TLegend(0.65,0.68,0.9,0.88);
  TLegend *le = new TLegend(0.70,0.73,0.94,0.93);
  le->SetBorderSize(0);
  le->SetFillStyle(0);
  le->AddEntry(histo1D["eff_pT_B"],"b jets","lp");
  le->AddEntry(histo1D["eff_pT_C"],"c jets","lp");
  le->AddEntry(histo1D["eff_pT_Light"],"light-flavour jets","lp");
  le->Draw();
  ce->Update();
  ce->Write();
  ce->SaveAs("BTagEfficiencies_eta.png");
  
  fout->Write();
  fout->Close();
  f_->Close();
  
  delete f_;
  delete fout;
}

