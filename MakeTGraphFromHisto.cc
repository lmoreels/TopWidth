#include <stdio.h>
#include "TStyle.h"
#include <ctime>
#include <cmath>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <map>
#include <array>
#include <vector>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>

using namespace std;


string inputDate = "170425_1455/";
string suffix = "";
map<string,TH1F*> histo;
map<string,TGraph*> graph;

/// Functions
string ConvertIntToString(int Number, int pad);
string MakeTimeStamp();
bool fexists(const char *filename);
void ClearVars();
void DrawGraph(TH1F* h, TGraph* g, string name);
void WriteOutput(int nPoints, double *arrayCentre, double *arrayContent);


int binMin, binMax;
ofstream txtOutput;


int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "********************************" << endl;
  cout << "***  Make TGraph from Histo  ***" << endl;
  cout << "********************************" << endl;
  cout << "* Current time: " << dateString << "    *" << endl;
  cout << "********************************" << endl;
  
  string inputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputPlots/mu/"+inputDate+"NtuplePlots_nominal.root";
  string outputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/TGraphFits/GraphFunctions"+suffix+".root";
  
  TFile *fileIn = new TFile(inputFileName.c_str(),"read");
  TFile *fileOut = new TFile(outputFileName.c_str(),"recreate");
  
  string listCats[] = {"CP", "WP", "UP"};
  const int sizeCats = sizeof(listCats)/sizeof(listCats[0]);
  int totEvents = 0;
  int nEvents[sizeCats] = {0};
  double fracCats[sizeCats] = {0};
  
  int nBins[sizeCats] = {0};
  
  std::array<std::vector<double>, sizeCats> vBinCentres;
  std::array<std::vector<double>, sizeCats> vBinContents;
  
  TH1F *histoTotal;
  
  
  for (int iCat = 0; iCat < sizeCats; iCat++)
  {
    ClearVars();
    
    histo[listCats[iCat]] = (TH1F*) fileIn->Get(("1D_histograms/Red_top_mass_"+listCats[iCat]+"_60b").c_str());
    histo[listCats[iCat]]->Smooth();
    nBins[iCat] = histo[listCats[iCat]]->GetNbinsX();
    if ( iCat > 0 && nBins[iCat] != nBins[iCat-1] )
    {
      cerr << "Histograms should have the same number of bins! Exiting..." << endl;
      exit(1);
    }
    
    binMin = histo[listCats[iCat]]->FindBin(0.6);
    binMax = histo[listCats[iCat]]->FindBin(1.4);
    for (int iBin = binMin; iBin < binMax+1; iBin++) { nEvents[iCat] += histo[listCats[iCat]]->GetBinContent(iBin);}
    totEvents += nEvents[iCat];
    
    cout << listCats[iCat] << ": " << nEvents[iCat] << " events" << endl;
    
    /// Normalise histo on relevant subdomain
    Double_t integral = histo[listCats[iCat]]->Integral(binMin, binMax);
    histo[listCats[iCat]]->Scale(1./integral);
    
    /// Get bin centres & contents
    for (int iBin = binMin; iBin < binMax+1; iBin++)
    {
      (vBinCentres[iCat]).push_back(histo[listCats[iCat]]->GetBinCenter(iBin));
      (vBinContents[iCat]).push_back(histo[listCats[iCat]]->GetBinContent(iBin));
    }
    
  }  /// end cats
  
  cout << "Total: " << totEvents << " events" << endl;
  
  /// Make arrays as input for TGraph
  const int nPoints = vBinCentres[0].size();
  double binCentreArray[nPoints] = {0.}, binContentArray[sizeCats][nPoints] = {0.}, totalBinContentArray[nPoints] = {0.};
  for (int i = 0; i < nPoints; i++)
  {
    binCentreArray[i] = vBinCentres[0][i];
    for (int iCat = 0; iCat < sizeCats; iCat++)
    {
      binContentArray[iCat][i] = vBinContents[iCat][i];
      // Calculate event fractions (one time)
      if ( i == 0 ) { fracCats[iCat] = (double)nEvents[iCat]/(double)totEvents;}
    }
  }
  
  /// Make TGraphs
  for (int iCat = 0; iCat < sizeCats; iCat++)
  {
    graph[listCats[iCat]] = new TGraph(nPoints, binCentreArray, binContentArray[iCat]);
    graph[listCats[iCat]]->Write();
    DrawGraph(histo[listCats[iCat]], graph[listCats[iCat]], "Graph_Red_top_mass_"+listCats[iCat]);
    
    fracCats[iCat] = (double)nEvents[iCat]/(double)totEvents;
    for (int i = 0; i < nPoints; i++) totalBinContentArray[i] += fracCats[iCat]*binContentArray[iCat][i];
    
    histo[listCats[iCat]]->Scale(fracCats[iCat]);
    if ( iCat == 0 )
    {
      histoTotal = (TH1F*) histo[listCats[iCat]]->Clone("SumOfCats");
      histoTotal->SetTitle("#frac{n_{CP}}{n_{tot}} * f_{CP}(x|#Gamma) + #frac{n_{WP}}{n_{tot}} * f_{WP}(x|#Gamma) + #frac{n_{UP}}{n_{tot}} * f_{UP}(x|#Gamma)");
    }
    else histoTotal->Add(histo[listCats[iCat]]);
  }
  
  cout << "The integral of the weighted likelihood histogram is " << histoTotal->Integral(binMin, binMax) << endl;
  
  /// Make likelihood function
  graph["likelihood"] = new TGraph(nPoints, binCentreArray, totalBinContentArray);
  DrawGraph(histoTotal, graph["likelihood"], "Graph_likelihood");
  
  WriteOutput(nPoints, binCentreArray, totalBinContentArray);
  
  
  /// Close files
  fileOut->Close();
  fileIn->Close();
  
  return 0;
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

void ClearVars()
{
  binMin = 0;
  binMax = 100;
}

void DrawGraph(TH1F* h, TGraph* g, string name)
{
  TCanvas *c = new TCanvas(name.c_str(), name.c_str());
  c->cd();
  h->GetXaxis()->SetTitle("m_{t}/<m_{t}> [GeV]");
  h->Draw("C");
  g->SetLineColor(kRed);
  g->Draw("same");
  c->Update();
  c->Write();
  c->SaveAs((name+".png").c_str());
}

void WriteOutput(int nPoints, double *arrayCentre, double *arrayContent)
{
  string outputTxtName = "TGraphFits/output_tgraph"+suffix+".txt";
  txtOutput.open(outputTxtName.c_str());
  
  txtOutput << "BinCentres:  ";
  for (int i = 0; i < nPoints; i++)
  {
    txtOutput << setw(5) << right << arrayCentre[i] << "  ";
  }
  txtOutput << endl;
  txtOutput << "BinContents:  ";
  for (int i = 0; i < nPoints; i++)
  {
    txtOutput << setw(5) << right << arrayContent[i] << "  ";
  }
  txtOutput << endl;
  
  txtOutput << endl;
  txtOutput.close();
}
