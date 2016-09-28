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
// par[0] = scale; par[1] = mean gaus & lorentz; par[2] = sigma gaus; par[3] = width lorentz
Double_t fitf(Double_t *x, Double_t *par) {
  Double_t arg_g = 0;
  if ( par[2] != 0 ) arg_g = (x[0] - par[1])/par[2];
  Double_t gaus = TMath::Exp(-0.5*arg_g*arg_g);
  if ( par[2] != 0 ) gaus = gaus/(par[2]*sqrt(2*TMath::Pi()));
  
  Double_t arg_l = par[3]/2.;
  Double_t arg_lx = x[0]-par[1];
  Double_t lorentz = (1./TMath::Pi())*arg_l/(arg_lx*arg_lx + arg_l*arg_l);
  
  Double_t fitfunc = par[0]*gaus + (1 - par[0])*lorentz;
  return fitfunc*par[4];
  //return gaussian(x,par) + lorentzian(x,&par[3]);
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
  const string histDir = "1D_histograms/";
  const string histoNames[] = {"mTop_div_aveMTop_TT_matched_reco", "mTop_div_aveMTop_TT_corr_match_chi2_reco"};
  int sizeHistos = sizeof(histoNames)/sizeof(histoNames[0]);
  
  /// Declare input and output files
  //string systInput = "160928_1258/NtuplePlots_nominal.root";  //nominal
  string systInput = "160928_1152/NtuplePlots_JERup.root ";  //JER
  string inputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputPlots/mu/"+systInput;
  string pathOutput = "OutputVoigt/";
  mkdir(pathOutput.c_str(),0777);
  string outputFileName = pathOutput+"VoigtFit_nominal.root";
  
  TFile *fin = new TFile(inputFileName.c_str(), "READ");
  fin->cd();
  
  TFile *fout = new TFile(outputFileName.c_str(), "RECREATE");
  fout->cd();
  
  
  for (int iHisto = 0; iHisto < sizeHistos; iHisto++)
  {
    //if (! fin->GetListOfKeys()->Contains(("1D_histograms/"+histoNames[iHisto]).c_str()) )  // does not work for histograms in dir
    //{
      //cout << " *** Histogram " << histoNames[iHisto] << " does not exist... Proceeding to next histogram..." << endl;
      //continue;
    //}
    cout << " *** Processing histogram " << histoNames[iHisto] << endl;
    
    /// Get histos
    TH1F* histo = (TH1F*) fin->Get((histDir+histoNames[iHisto]).c_str());
    histo->Write();
    
    /// Declare fit function
    int nBins = histo->GetXaxis()->GetNbins();
    const int nPar = 5;
    string parnames[nPar] = {"a1","a2","a3","a4","a5"};
    //string parnames[nPar] = {"a1","a2","a3","a4","a5","a6"};
    
    
    // Normalise histogram
    Double_t scale = 1./histo->Integral();
    histo->Scale(scale);
    
    
    /// Declare fit function
    TF1 *myfit = new TF1("fit", fitf, 0.5, 2., nPar);  // root name, fit function, range, nPar
    
    
    //  Give names to the parameters: // 0->a1 = scale; 1->a2 = mean gaus & lorentz; 2->a3 = sigma gaus; 3->a4 = width lorentz
    myfit->SetParNames("ps","#mu","#sigma","#gamma","norm");
    
    myfit->SetParameters(0.5, 1., 0.58*histo->GetRMS(), 1.36*histo->GetRMS(), 1.);  // sigma = 0.57735027 * RMS; FWHM = 1.359556 * RMS (= 2.354820 * sigma)
    myfit->SetParLimits(0,0.,1.);
    myfit->SetParLimits(1,0.95,1.05);
    myfit->SetParLimits(2,0.001,1.e+3);
    myfit->SetParLimits(3,0.001,1.e+3);
    
    
    ///  Fit
    std::string func_title = std::string(histo->GetName())+"_Fitted";
    myfit->SetName(func_title.c_str());
    histo->Fit(myfit);
    gStyle->SetOptFit(0111);
    histo->SetStats(1);
    histo->Write();
    myfit->Write();
    
    delete myfit;
    
  }  /// end loop histos
  
  
  
  
  
  fout->Write();
  
  
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
 
