#include <stdio.h>
//#include "TSystem.h"
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
#include <TCanvas.h>

using namespace std;


string systStr = "nominal";

string whichDate(string syst)
{
  if ( syst.find("nominal") == 0 ) return "161124_1329/NtuplePlots_nominal.root";
  //else if ( syst.find("JERup") == 0 ) return "161116_1401/NtuplePlots_JERup.root";
  //else if ( syst.find("JERdown") == 0 ) return "161116_1444/NtuplePlots_JERdown.root";
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    return "161124_1329/NtuplePlots_nominal.root";
  }
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
Double_t voigt(Double_t *x, Double_t *par) {
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

// Crystal Ball
Double_t crysBall(Double_t *x, Double_t *par) {
  // params: alpha, n, sigma, mu
  if ( par[2] < 0. ) return 0.;
  Double_t alpha = fabs(par[0]);
  Double_t A = pow( par[1]/alpha , par[1]) * exp(-alpha*alpha/2.);
  Double_t B = par[1]/alpha - alpha;
  //Double_t C = par[1]/alpha * 1/(par[1]-1) * exp(-alpha*alpha/2.);
  //Double_t D = sqrt(TMath::Pi()/2.) * (1 + erf(alpha/sqrt(2)));
  //Double_t N = 1/(par[2]*(C+D));
  
  Double_t ref = (x[0] - par[3])/par[2];  // (x-mean)/sigma
  if ( par[0] < 0 ) ref = -ref;
  //Double_t fitfunc = N;
  Double_t fitfunc = 1.;
  if ( ref > -alpha ) fitfunc = fitfunc * exp(-ref*ref/2.);
  else if (ref <= -alpha ) fitfunc = fitfunc * A * pow ( B - ref , -par[1]);
  return fitfunc*par[4];
}


int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  //gSystem->Load("libMathCore.so");
  
  /// Declare histos to be fitted
  const string histDir = "1D_histograms/";
  pair<const string, int> histoNames[] = { {"mTop_div_aveMTop_TT_matched_reco", 1}, {"mTop_div_aveMTop_TT_corr_match_chi2_reco", 1}, {"mTop_div_aveMTop_TT_wrong_perm_chi2_reco", 1}, {"mTop_div_aveMTop_TT_wrong_match_chi2_reco", 0}, {"mTop_div_aveMTop_TT_wrong_jets_chi2_reco", 0} };
  int sizeHistos = sizeof(histoNames)/sizeof(histoNames[0]);
  
  /// Declare input and output files
  string inputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputPlots/mu/"+whichDate(systStr);
  string pathOutput = "OutputVoigt/";
  mkdir(pathOutput.c_str(),0777);
  string outputFileName = pathOutput+"VoigtFit_"+systStr+".root";
  
  
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
    cout << " *** Processing histogram " << histoNames[iHisto].first << endl;
    
    /// Get histos
    TH1F* histo = (TH1F*) fin->Get((histDir+histoNames[iHisto].first).c_str());
    histo->Write();
    
    // Normalise histogram
    Double_t scale = 1./histo->Integral();
    histo->Scale(scale);
    
    
    /// Declare fit function
    int nBins = histo->GetXaxis()->GetNbins();
    const int nPar = 5;
    
    TF1 *myfit;
    
    if ( histoNames[iHisto].second == 1 )
    {
      myfit = new TF1("fit", voigt, 0.5, 2., nPar);  // root name, fit function, range, nPar
      
      myfit->SetParNames("ps","#mu","#sigma","#gamma","norm");
      myfit->SetParameters(0.5, 1., 0.58*histo->GetRMS(), 1.36*histo->GetRMS(), 1.);  // sigma = 0.57735027 * RMS; FWHM = 1.359556 * RMS (= 2.354820 * sigma)
      
      myfit->SetParLimits(0,0.,1.);
      myfit->SetParLimits(1,0.95,1.05);
      myfit->SetParLimits(2,0.001,1.e+3);
      myfit->SetParLimits(3,0.001,1.e+3);
    }
    else
    {
      myfit = new TF1("fit", crysBall, 0.6, 2., nPar);  // root name, fit function, range, nPar
      
      myfit->SetParNames("alpha","n","#sigma","#mu","norm");
      myfit->SetParameters(-0.59, 17.7, 0.2, 0.8, 1.);   // wrong match
      
      myfit->SetParLimits(2,0.00001,1.e+3);
      myfit->SetParLimits(3, 0.74, 0.86);
    }
    
    
    ///  Fit
    std::string func_title = std::string(histo->GetName())+"_Fitted";
    myfit->SetName(func_title.c_str());
    histo->Fit(myfit, "R");
    gStyle->SetOptFit(0111);
    histo->SetStats(1);
    histo->Write();
    myfit->Write();
    
    if ( histoNames[iHisto].second == 0 )
    {
      TCanvas* c1 = new TCanvas("c1", "Fit function expanded");
      c1->cd();
      histo->SetLineColor(kBlue);
      histo->Draw();
      
      TF1 *fcb = new TF1("fcb", crysBall, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(), 5);
      fcb->SetParameters(myfit->GetParameter(0), myfit->GetParameter(1), myfit->GetParameter(2), myfit->GetParameter(3), myfit->GetParameter(4));
      fcb->SetLineColor(kRed);
      fcb->Draw("same");
      c1->Update();
      c1->Write();
      c1->SaveAs((pathOutput+func_title+".png").c_str());
      c1->Close();
      
      delete fcb;
      delete c1;
    }
    
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
 
