#include <stdio.h>
#include "TStyle.h"
#include <ctime>
#include <cmath>
#include <TMath.h>
#include <iostream>
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
#include <TFile.h>

using namespace std;


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

/// Declare fit functions
Double_t gaussian(Double_t *x,Double_t *par) {
  Double_t arg = 0;
  if ( par[2] != 0 ) arg = (x[0] - par[1])/par[2];
  Double_t fitval = TMath::Exp(-0.5*arg*arg);
  if ( par[2] != 0 ) fitval = fitval/(par[2]*sqrt(2*TMath::Pi()));
  return fitval;
}

Double_t lorentzian(Double_t *x,Double_t *par) {
  Double_t arg1 = par[1]/2.;
  Double_t arg2 = x[0]-par[2];
  Double_t fitval = (par[0]/TMath::Pi())*arg1/(arg2*arg2 + arg1*arg1);
  return fitval;
}

// Sum of gaussian and lorentzian
Double_t fitf(Double_t *x, Double_t *par) {
  return gaussian(x,par) + lorentzian(x,&par[3]);
}



int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "               *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  
  /// Declare histos to be fitted
  const string histoNames[] = {""};
  int sizeHistos = sizeof(histoNames)/sizeof(histoNames[0]);
  
  /// Declare input and output files
  string inputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputPlots/mu/160901_1151/NtuplePlots_nominal.root";
  string outputFileName = "";
  
  TFile *fin = new TFile(inputFileName.c_str(), "READ");
  fin->cd();
  
  TFile *fout = new TFile(outputFileName.c_str(), "RECREATE");
  fout->cd();
  
  
  for (int iHisto = 0; iHisto < sizeHistos; iHisto++)
  {
    if (! fin->GetListOfKeys()->Contains((histoNames[iHisto]).c_str()) )
    {
      cout << " *** Histogram " << histoNames[iHisto] << " does not exist... Proceeding to next histogram..." << endl;
      continue;
    }
    cout << " *** Processing histogram " << histoNames[iHisto] << endl;
    
    /// Get histos
    TH1F* histo = (TH1F*) fin->Get((histoNames[iHisto]).c_str());
    histo->Write();
    
    /// Declare fit function
    int nBins = histo->GetXaxis()->GetNbins();
    const int nPar = 6;
    TH1D **hlist = new TH1D*[nPar];
    string parnames[nPar] = {"a1","a2","a3","a4","a5","a6"};
    string name = "";
    string title = "";
    
    for (int iPar = 0; iPar < nPar; iPar++)
    {
      name = std::string(histo->GetName())+ "_" + parnames[iPar];
      title = std::string(histo->GetName())+ ": Fitted value of " + parnames[iPar] ;
      hlist[iPar] = new TH1D(name.c_str(), title.c_str(), nBins, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
      hlist[iPar]->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    }
    
    // Normalise histogram
    Double_t scale = 1./histo->Integral();
    histo->Scale(scale);
    
    
    /// Declare fit function
    TF1 *myfit = new TF1("fit", fitf, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(), nPar);  // root name, fit function, range, nPar
    
    
    //  Give names to the parameters: 0->a1: gaus constant; 1->a2: gaus mean; 2->a3: gaus sigma; 3->a4: ltz constant; 4->a5: ltz gamma (FWHM); 5->a6: ltz displacement
    myfit->SetParNames("a1","a2","a3","a4","a5","a6");
    
    myfit->SetParameters(0.5, histo->GetMean(), 0.58*histo->GetRMS(), 0.5, 1.36*histo->GetRMS(), 1.);  // sigma = 0.57735027 * RMS; FWHM = 1.359556 * RMS (= 2.354820 * sigma)
    myfit->SetParLimits(2,0,1.e+6);
    myfit->SetParLimits(4,0,1.e+6);
    
    
    ///  Fit
    std::string func_title = std::string(histo->GetName())+"_Fitted";
    myfit->SetName(func_title.c_str());
    histo->Fit(myfit);
    histo->Write();
    myfit->Write();
    
    delete myfit;
    
  }  /// end loop histos
  
  
  
  
  
  
  
  
  fin->Close();
  fout->Close();
  
  delete fin;
  delete fout;
  
  
  ///// end
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