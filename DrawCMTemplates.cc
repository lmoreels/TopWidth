//#include <stdio.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
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
#include <TLegendEntry.h>
#include <TAttLine.h>
#include <TAttText.h>
#include <TColor.h>
#include <TLatex.h>
//#include <TMathText.h>
#include <TPaveStats.h>


using namespace std;


string pathOutput = "test/";
string pathTempPlots = "LikelihoodTemplates/";

string ConvertDoubleToString(double Number)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

string PDotReplace(double var)
{
  std::string str = ConvertDoubleToString(var);
  std::replace(str.begin(), str.end(), 'p', '.');
  return str;
}

string PDotReplace(string str)
{
  std::replace(str.begin(), str.end(), 'p', '.');
  return str;
}

double GetBinWidth(TH1F* histo, int precision = 2)
{
  int nBins = histo->GetNbinsX();
  double min = histo->GetXaxis()->GetBinLowEdge(1);
  double max = histo->GetXaxis()->GetBinLowEdge(nBins);
  double binWidth = (max-min)/nBins;
  cout << "Bin width = " << binWidth;
  
  if      ( precision == 0 ) binWidth = floor(binWidth+0.5);  // round number to 0 decimal places
  else if ( precision == 1 ) binWidth = floor(binWidth*10.+0.5)/10.;  // round number to 1 decimal place
  else if ( precision == 2 ) binWidth = floor(binWidth*100.+0.5)/100.;  // round number to 2 decimal places
  else if ( precision == 3 ) binWidth = floor(binWidth*1000.+0.5)/1000.;  // round number to 3 decimal places
  else                      binWidth = floor(binWidth*10000.+0.5)/10000.;  // round number to 4 decimal places
  cout << " ; rounded: " << binWidth << endl;
  
  return binWidth;
}

void DrawCMTemplates()
{
  /// Define input files
  //string pathPlotsTemplates[] = {pathTempPlots+"171223_1022/", pathTempPlots+"171222_2301/", pathTempPlots+"171223_1330/", pathTempPlots+"171223_1609/", pathTempPlots+"171223_1414/", pathTempPlots+"171223_1510/", pathTempPlots+"171223_1438/"};
  //string bounds[] = {"1p6", "1p7", "1p8", "1p85", "1p9", "1p95", "2p0"};
  string pathPlotsTemplates[] = {pathTempPlots+"171222_1600/", pathTempPlots+"171222_2041/", pathTempPlots+"171222_1301/", pathTempPlots+"171222_2223/", pathTempPlots+"171223_1540/", pathTempPlots+"171223_2023/", pathTempPlots+"171223_2100/", pathTempPlots+"171223_2127/", pathTempPlots+"171223_2151/"};
  string bounds[] = {"50", "55", "60", "65", "75", "80", "85", "90", "100"};
  int nBounds = sizeof(bounds)/sizeof(bounds[0]);
  string fileNameTemplates = "TGraphFunctions.root";
  
  
  TFile *fileOut = new TFile((pathOutput+"Templates_CM_overlay.root").c_str(),"RECREATE");
  fileOut->cd();
  
  TCanvas *c1, *c2, *c3, *c4;
  map<string,TF1*> TF;
  map<string,TH1F*> histo;
  map<string,TGraph*> graph;
  string histoName = "";
  string labelRedTopMass = "m_{lb,r}";
  string legendXWidth = "#Gamma_{t,gen} x ";
  
  string binWidthStr;
  
  
  //gStyle->SetOptStat(111111); //8x
  gStyle->SetOptStat(1110);  // display entries, mean & RMS, but no title
  gStyle->SetOptFit(0111);
  gStyle->SetPalette(1,0);
  Color_t listColour[] = {kBlue-2, kRed+1, kGreen-3, kMagenta, kOrange};
  Color_t listColour2[] = {kBlue-2, kGreen-3, kRed+1};
  Color_t listRainbow[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+1, kCyan+1, kBlue+2, kMagenta, kViolet-5, kPink+10};
  Color_t listRainbowShort[] = {kPink-1, kOrange+7, kOrange-4, kGreen+1, kBlue+2, kMagenta};
  cout << "Defined colours" << endl;
  Style_t listMarkers[] = {24, 26, 32, 25, 30, 27, 28, 46, 42};
  
  
  
  string widths[] = {"0p5", "0p75", "1", "1p25", "1p5"};
  //string widths[] = {"2p5", "3", "3p5"};
  const int nWidths = sizeof(widths)/sizeof(widths[0]);
  Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen+2/*, kCyan+1*/, kBlue+2, kMagenta};
  
  
  for (int iBound = 0; iBound < nBounds; iBound++)
  {
    /// Get graphs
    TFile *fileInTemplates = new TFile((pathPlotsTemplates[iBound]+fileNameTemplates).c_str(),"read");
    fileInTemplates->cd();
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      histoName = "CM_widthx"+widths[iWidth];
      graph[histoName+"_bound"+bounds[iBound]] = (TGraph*) fileInTemplates->Get(("g"+histoName).c_str());
    }
    
    /// Draw graphs
    fileOut->cd();
    c1 = new TCanvas(("CM template "+bounds[iBound]).c_str(),("Overlay CM templates "+bounds[iBound]).c_str());
    c1->cd();
    
    TLegend *leg = new TLegend(0.71,0.52,0.89,0.86);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    
    
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      histoName = "CM_widthx"+widths[iWidth]+"_bound"+bounds[iBound];
      //graph[histoName]->SetLineWidth(2);
      graph[histoName]->SetLineColor(colours[iWidth]);
      if (iWidth == 0)
      {
        graph[histoName]->SetTitle("");
        graph[histoName]->Draw("AL");
        graph[histoName]->GetXaxis()->SetTitle(labelRedTopMass.c_str());
        graph[histoName]->GetYaxis()->SetTitle("\\mathscr{L}(\\mathrm{\\Gamma_{t}})");
        //graph[histoName]->GetXaxis()->SetRangeUser(168,177);
        graph[histoName]->Draw("AL");
      }
      else graph[histoName]->Draw("same");
      leg->AddEntry(graph[histoName],(legendXWidth+PDotReplace(widths[iWidth])).c_str(),"l");
    }
    leg->Draw();
    
    fileOut->cd();
    c1->Update();
    c1->Write();
    c1->SaveAs((pathOutput+"Template_CM_overlay_"+bounds[iBound]+".png").c_str());
    
    fileInTemplates->Close();
  }
  
  
  fileOut->Close();

}
  
  
	/*
	 For all histogram types: nbins, xlow, xup
	 bin = 0;       underflow bin
	 bin = 1;       first bin with low-edge xlow INCLUDED
	 bin = nbins;   last bin with upper-edge xup EXCLUDED
	 bin = nbins+1; overflow bin
	 */
