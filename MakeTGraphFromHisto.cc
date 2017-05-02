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
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TCanvas.h>

using namespace std;


bool useDifferentWidths = true;

string inputDate = "170426_0952";
string suffix = "_90b";
pair<double,string> input[] = { pair<double,string>(0.2,"170502_1755"), pair<double,string>(0.3,"170502_1756"), pair<double,string>(0.4,"170502_1757"), pair<double,string>(0.5,"170426_0949"), pair<double,string>(0.66,"170426_0950"), pair<double,string>(0.75,"170426_0951"), pair<double,string>(0.8,"170502_1131"), pair<double,string>(0.9,"170502_1132"), pair<double,string>(1.,"170426_0952"), pair<double,string>(1.2,"170502_1254"), pair<double,string>(1.5,"170502_1133"), pair<double,string>(2.,"170426_0956"), pair<double,string>(2.5,"170502_1134"), pair<double,string>(3.,"170426_0957"), pair<double,string>(3.5,"170502_1135"), pair<double,string>(4.,"170426_0958"), pair<double,string>(4.5,"170502_1152"), pair<double,string>(5.,"170502_1153"), pair<double,string>(5.5,"170502_1154"), pair<double,string>(6.,"170502_1155"), pair<double,string>(6.5,"170502_1156"),pair<double,string>(7.,"170502_1251"), pair<double,string>(7.5,"170502_1252"),pair<double,string>(8.,"170502_1253") };

map<string,TH1F*> histo, histoTotal;
map<string,TGraph*> graph;

/// Functions
string ConvertDoubleToString(double Number);
string DotReplace(double var);
string ConvertIntToString(int Number, int pad);
string MakeTimeStamp();
bool fexists(const char *filename);
void ClearVars();
void DrawGraph(TH1F* h, TGraph* g, string name);
void DrawGraph2D(TGraph2D* g, string name);
void DrawLikelihoods(double nWidths, pair<double,string> *input);
void WriteOutput(int nPoints, double *arrayCentre, double *arrayContent, string name);
void WriteToFile(TGraph2D* g, string name);


int binMin, binMax;
string histoName = "";
ofstream txtOutput;


int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "********************************" << endl;
  cout << "***  Make TGraph from Histo  ***" << endl;
  cout << "********************************" << endl;
  cout << "* Current time: " << dateString << "    *" << endl;
  cout << "********************************" << endl;
  
  
  string outputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/TGraphFits/GraphFunctions.root";
  TFile *fileOut = new TFile(outputFileName.c_str(),"recreate");
  
  TGraph2D *gLL2D = new TGraph2D();
    
  string listCats[] = {"CP", "WP", "UP"};
  const int sizeCats = sizeof(listCats)/sizeof(listCats[0]);
  int totEvents = 0;
  int nEvents[sizeCats] = {0};
  double fracCats[sizeCats] = {0};
  
  int nBins[sizeCats] = {0};
  
  std::array<std::vector<double>, sizeCats> vBinCentres;
  std::array<std::vector<double>, sizeCats> vBinContents;
  
  double nWidths = 1;
  if (useDifferentWidths) nWidths = sizeof(input)/sizeof(input[0]);
  cout << "Running over " << nWidths << " widths" << endl;
  
  for (int iWidth = 0; iWidth < nWidths; iWidth++)
  {
    if (useDifferentWidths)
    {
      suffix = "_widthx"+DotReplace(input[iWidth].first)+"_90b";
      inputDate = input[iWidth].second;
    }
    
    string inputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputPlots/mu/"+inputDate+"/NtuplePlots_nominal.root";
    TFile *fileIn = new TFile(inputFileName.c_str(),"read");
    
    fileOut->cd();
    
    totEvents = 0;
    for (int i = 0; i < sizeCats; i++)
    {
      nBins[i] = 0;
      nEvents[i] = 0;
      fracCats[i] = 0.;
      vBinCentres[i].clear();
      vBinContents[i].clear();
    }
    
    for (int iCat = 0; iCat < sizeCats; iCat++)
    {
      ClearVars();
      
      histoName = listCats[iCat]+suffix;
      
      histo[histoName] = (TH1F*) fileIn->Get(("1D_histograms/Red_top_mass_"+listCats[iCat]+"_90b").c_str());
      histo[histoName]->Smooth();
      nBins[iCat] = histo[histoName]->GetNbinsX();
      if ( iCat > 0 && nBins[iCat] != nBins[iCat-1] )
      {
        cerr << "Histograms should have the same number of bins! Exiting..." << endl;
        exit(1);
      }
      
      binMin = histo[histoName]->FindBin(0.6);
      binMax = histo[histoName]->FindBin(1.4)+1;
      for (int iBin = binMin; iBin < binMax+1; iBin++) { nEvents[iCat] += histo[histoName]->GetBinContent(iBin);}
      totEvents += nEvents[iCat];
      
      cout << histoName << ": " << nEvents[iCat] << " events" << endl;
      
      /// Normalise histo on relevant subdomain
      Double_t integral = histo[histoName]->Integral(binMin, binMax);
      histo[histoName]->Scale(1./integral);
      
      /// Get bin centres & contents
      for (int iBin = binMin; iBin < binMax+1; iBin++)
      {
        (vBinCentres[iCat]).push_back(histo[histoName]->GetBinCenter(iBin));
        (vBinContents[iCat]).push_back(histo[histoName]->GetBinContent(iBin));
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
      histoName = listCats[iCat]+suffix;
      
      graph[histoName] = new TGraph(nPoints, binCentreArray, binContentArray[iCat]);
      graph[histoName]->SetName(histoName.c_str());
      graph[histoName]->SetTitle(histoName.c_str());
      graph[histoName]->Write();
      DrawGraph(histo[histoName], graph[histoName], "Graph_Red_top_mass_"+histoName);

      fracCats[iCat] = (double)nEvents[iCat]/(double)totEvents;
      for (int i = 0; i < nPoints; i++) totalBinContentArray[i] += fracCats[iCat]*binContentArray[iCat][i];

      histo[histoName]->Scale(fracCats[iCat]);
      if ( iCat == 0 )
      {
        histoTotal[suffix] = (TH1F*) histo[histoName]->Clone("TotalProbability");
        histoTotal[suffix]->SetTitle("#frac{n_{CP}}{n_{tot}} * f_{CP}(x|#Gamma) + #frac{n_{WP}}{n_{tot}} * f_{WP}(x|#Gamma) + #frac{n_{UP}}{n_{tot}} * f_{UP}(x|#Gamma)");
      }
      else histoTotal[suffix]->Add(histo[histoName]);
    }
    
    cout << "The integral of the weighted probability histogram is " << histoTotal[suffix]->Integral(binMin, binMax) << endl;
    
    double likelihoodArray[nPoints] = {0.};
    for (int i = 0; i < nPoints; i++)
    {
      likelihoodArray[i] = -TMath::Log(totalBinContentArray[i]); 
    }
    
    /// Make likelihood function
    graph["likelihood"+suffix] = new TGraph(nPoints, binCentreArray, likelihoodArray);
    graph["likelihood"+suffix]->SetName(("likelihood"+suffix).c_str());
    graph["likelihood"+suffix]->SetTitle(("likelihood"+suffix).c_str());
    graph["likelihood"+suffix]->Write();
    DrawGraph(histoTotal[suffix], graph["likelihood"+suffix], "Graph_likelihood"+suffix);
    
    WriteOutput(nPoints, binCentreArray, totalBinContentArray, suffix);
    
    // Fill TGraph2D
    if (useDifferentWidths)
    {
      for (int iCentre = 0; iCentre < nPoints; iCentre++)
      {
        gLL2D->SetPoint(gLL2D->GetN(), binCentreArray[iCentre], input[iWidth].first, likelihoodArray[iCentre]);
        //cout << setw(5) << gLL2D->GetN() << "  Width: " << setw(5) << input[iWidth].first << "   Bin Centre: " << setw(10) << binCentreArray[iCentre] << "   Bin Content: " << totalBinContentArray[iCentre] << endl;
      }
    }
    
    /// Close input file
    fileIn->Close();
  }  // end nWidths
  
  if (useDifferentWidths)
  {
    WriteToFile(gLL2D, "LogLikelihoodFunction");
    
    fileOut->cd();
    DrawLikelihoods(nWidths, input);
    gLL2D->SetName("2D_likelihood");
    gLL2D->SetTitle("2D_likelihood");
    gLL2D->Write();
    DrawGraph2D(gLL2D, "2D_likelihood");
  }
  //cout << "TGraph evaluated in (1,1): " << gLL2D->Interpolate(1,1) << " ;  (1,1.01): " << gLL2D->Interpolate(1,1.01) << " ;  (1,1.02): " << gLL2D->Interpolate(1,1.02) << " ;  (1,1.03): " << gLL2D->Interpolate(1,1.03) << " ;  (1,1.05): " << gLL2D->Interpolate(1,1.05) << endl;
  
  /// Close output file
  fileOut->Close();
  
  
  return 0;
}

string ConvertDoubleToString(double Number)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

string DotReplace(double var)
{
  string str = ConvertDoubleToString(var);
  replace(str.begin(), str.end(), '.', 'p');
  return str;
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
  histoName = "";
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
  c->SaveAs(("TGraphFits/"+name+".png").c_str());
  
  delete c;
}

void DrawGraph2D(TGraph2D* g, string name)
{
  TCanvas *c = new TCanvas(name.c_str(), name.c_str());
  c->cd();
  gStyle->SetPalette(1);
  g->SetMaxIter(500000);
  g->SetName(name.c_str());
  g->SetTitle((name+"; m_{t}/<m_{t}> [GeV]; #Gamma/#Gamma_{SM}").c_str());
  //g->GetHistogram("empty")->GetXaxis()->SetTitleOffset(1.5);  // Does not work
  //g->GetHistogram("empty")->GetYaxis()->SetTitleOffset(1.5);  // Does not work
  g->Draw("surf1");
  c->Update();
  c->Write();
  c->SaveAs(("TGraphFits/"+name+"_surf1.png").c_str());
  g->Draw("colz");
  c->Update();
  c->Write();
  c->SaveAs(("TGraphFits/"+name+"_colz.png").c_str());
  g->SetMarkerSize(0.5);
  g->Draw("tri1 p0");
  c->Update();
  c->Write();
  c->SaveAs(("TGraphFits/"+name+"_tri1p0.png").c_str());
  
  delete c;
}

void DrawLikelihoods(double nWidths, pair<double,string> *input)
{
  Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+1, kCyan+1, kBlue+2, kMagenta, kViolet-5, kPink+10};
  TCanvas* c2 = new TCanvas("TotalProbability", "TotalProbability");
  c2->cd();
  graph["likelihood_widthx"+DotReplace(input[0].first)+"_90b"]->SetLineColor(colours[0]);
  graph["likelihood_widthx"+DotReplace(input[0].first)+"_90b"]->Draw();
  c2->Update();
  for (int i = 1; i < nWidths; i++)
  {
    graph["likelihood_widthx"+DotReplace(input[i].first)+"_90b"]->SetLineColor(colours[i%10]);
    graph["likelihood_widthx"+DotReplace(input[i].first)+"_90b"]->Draw("same");
    c2->Update();
  }
  c2->Write();
  c2->SaveAs("TGraphFits/Probabilities_90b.png");
  c2->Close();
  
  delete c2;
}

void WriteOutput(int nPoints, double *arrayCentre, double *arrayContent, string name)
{
  string outputTxtName = "TGraphFits/output_tgraph"+name+".txt";
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

void WriteToFile(TGraph2D* g, string name)
{
  TFile *file2DGraph = new TFile((name+".root").c_str(),"recreate");
  file2DGraph->cd();
  g->Write();
  g->SetName(name.c_str());
  g->SetTitle((name+"; m_{t}/<m_{t}> [GeV]; #Gamma/#Gamma_{SM}").c_str());
  g->Write();
  file2DGraph->Close();
  
  delete file2DGraph;
}
