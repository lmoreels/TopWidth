////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////



#include "TStyle.h"
#include <stdio.h>
#include <ctime>
#include <cmath>
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
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TGraph.h>
#include <TGraph2D.h>




using namespace std;
// using namespace reweight;
// using namespace TopTree;

bool test = false;
bool getGraphFromFile = false;
bool constructGraph = true;
bool constructGraph2D = false;

map<string,TGraph*> graph;

Double_t widthArray[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.};
const int nWidthsLL = sizeof(widthArray)/sizeof(widthArray[0]);


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


int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "*   Current time: " << dateString << "               *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  string inputFileName;
  pair<double, double> input[] = {pair<double, double>(1.,2.), pair<double, double>(1.,0.3), pair<double, double>(1.,0.4), pair<double, double>(1.,0.5), pair<double, double>(1.,0.6), pair<double, double>(1.,1.)};
  
  if (getGraphFromFile)
  {
    TGraph2D* graph = 0;
    TH2 *h = 0;
    
    inputFileName = "LogLikelihoodFunction.root";
    cout << "Reading TGraph2D from file " << inputFileName << endl;
    
    TFile *inputFileLL = new TFile(inputFileName.c_str(), "read");
    inputFileLL->GetObject("LogLikelihoodFunction", graph);
    inputFileLL->GetObject("hLogLikelihoodFunction", h);
    //graph = (TGraph2D*) inputFileLL->Get("LogLikelihoodFunction");
    //h = (TH2D*) inputFileLL->Get("hLogLikelihoodFunction");
    
    h->SetDirectory(0);
    //graph->SetHistogram(h);
    //graph->GetHistogram();
    //graph->SetHistogram(h); // SOMETHING WRONG HERE: SEG FAULT
    //graph->SetMaxIter(500000000);
    
    cout << "N points: " << graph->GetN() << endl;
    for (int i = 0; i < sizeof(input)/sizeof(input[0]); i++)
      cout << "Value in (" << input[i].first << "," << input[i].second << "): " << setw(5) << graph->Interpolate(input[i].first,input[i].second) << endl;
    
    
    TCanvas *c = new TCanvas("c", "c");
    c->cd();
    graph->Draw("colz");
    c->SaveAs("testTGraphs.png");
    delete c;
  }
  
  if (constructGraph2D)
  {
    inputFileName = "TGraphFits/output_tgraph2d_total.txt";
    cout << "Constructing TGraph2D from file " << inputFileName << endl;
    
    TGraph2D *g = new TGraph2D(inputFileName.c_str());
    g->GetHistogram("");
    g->SetMaxIter(500000000);
    
    cout << "N points: " << g->GetN() << endl;
    for (int i = 0; i < sizeof(input)/sizeof(input[0]); i++)
      cout << "Value in (" << input[i].first << "," << input[i].second << "): " << g->Interpolate(input[i].first,input[i].second) << endl;
    
    TCanvas *c = new TCanvas("c", "c");
    c->cd();
    g->Draw("tri1 p0");
    c->SaveAs("testTGraphs.png");
    delete c;
    
  }
  
  if (constructGraph)
  {
    cout << "Constructing TGraphs for each width" << endl;
    string suffix;
    for (int iWidth = 0; iWidth < nWidthsLL; iWidth++)
    {
      suffix = "widthx"+DotReplace(widthArray[iWidth]);
      inputFileName = "TGraphFits/output_tgraph1d_"+suffix+".txt";
      graph[suffix] = new TGraph(inputFileName.c_str());
    }
    
    for (int i = 0; i < sizeof(input)/sizeof(input[0]); i++)
      cout << "Value in (" << input[i].first << "," << input[i].second << "): " << graph["widthx"+DotReplace(input[i].second)]->Eval(input[i].first) << endl;
    
  }
  
  
  
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
