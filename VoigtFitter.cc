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


bool usePredef = true;
bool useTest = false;
bool debug = false;
string systStr = "nominal";
string suffix = "_widthx0p66";

string whichDate(string syst, string suff)
{
  if ( syst.find("nominal") == 0 ) 
  {
    if ( suff.find("widthx0p5") != std::string::npos )       return "170223_1244/NtuplePlots_nominal.root";
    else if ( suff.find("widthx0p66") != std::string::npos ) return "170223_1245/NtuplePlots_nominal.root";
    else if ( suff.find("widthx0p75") != std::string::npos ) return "170223_1246/NtuplePlots_nominal.root";
    else if ( suff.find("widthx1") != std::string::npos )    return "170223_1255/NtuplePlots_nominal.root";
    //else if ( suff.find("widthx1p5") != std::string::npos )  return "170213_1306/NtuplePlots_nominal.root";
    else if ( suff.find("widthx2") != std::string::npos )    return "170223_1256/NtuplePlots_nominal.root";
    else if ( suff.find("widthx3") != std::string::npos )    return "170223_1257/NtuplePlots_nominal.root";
    else if ( suff.find("widthx4") != std::string::npos )    return "170223_1258/NtuplePlots_nominal.root";
    else                                                     return "170223_1255/NtuplePlots_nominal.root";
    
  }
  //else if ( syst.find("JERup") == 0 ) return "161116_1401/NtuplePlots_JERup.root";
  //else if ( syst.find("JERdown") == 0 ) return "161116_1444/NtuplePlots_JERdown.root";
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    suffix = "";
    return "170223_1255/NtuplePlots_nominal.root";
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
  Double_t fitval = (par[0]/TMath::Pi())*arg1/TMath::Max(1.e-10,(arg2*arg2 + arg1*arg1));
  return fitval;
}

// Sum of gaussian and lorentzian
// par[0] = scale; par[1] = mean gaus & lorentz; par[2] = sigma gaus; par[3] = width lorentz
Double_t voigt(Double_t *x, Double_t *par) {
  if ( par[2] <= 0. ) return 0.;
  if ( par[3] <= 0. ) return 0.;
  Double_t arg_g = (x[0] - par[1])/par[2];
  Double_t gaus = TMath::Exp(-0.5*arg_g*arg_g)/(par[2]*sqrt(2*TMath::Pi()));
  
  Double_t arg_l = par[3]/2.;
  Double_t arg_lx = x[0]-par[1];
  Double_t lorentz = (1./TMath::Pi())*arg_l/(arg_lx*arg_lx + arg_l*arg_l);
  
  Double_t fitfunc = par[0]*gaus + (1 - par[0])*lorentz;
  return fitfunc*par[4];
  //return gaussian(x,par) + lorentzian(x,&par[3]);
}

Double_t voigt_predef(Double_t *x, Double_t *par) {
  // Computation of Voigt function (normalised).
  // Voigt is a convolution of
  // gauss(xx) = 1/(sqrt(2*pi)*sigma) * exp(xx*xx/(2*sigma*sigma)
  // and
  // lorentz(xx) = (1/pi) * (lg/2) / (xx*xx + g*g/4)
  // functions.
  //
  // The Voigt function is known to be the real part of Faddeeva function also
  // called complex error function [2].
  //
  // The algoritm was developed by J. Humlicek [1].
  // This code is based on fortran code presented by R. J. Wells [2].
  // Translated and adapted by Miha D. Puc
  //
  // To calculate the Faddeeva function with relative error less than 10^(-r).
  // r can be set by the the user subject to the constraints 2 <= r <= 5.
  
  /// Voigt(x, sigma, lg, r)
  return TMath::Voigt(x[0]-par[0], par[1], par[2], par[3])*par[4];
}

// Crystal Ball
Double_t crysBall(Double_t *x, Double_t *par) {
  // params: alpha, n, sigma, mu
  if ( par[2] <= 0. ) return 0.;
  Double_t alpha = fabs(par[0]);
  Double_t A = pow( par[1]/alpha , par[1]) * exp(-alpha*alpha/2.);
  Double_t B = par[1]/alpha - alpha;
  //Double_t C = par[1]/alpha * 1/(par[1]-1) * exp(-alpha*alpha/2.);
  //Double_t D = sqrt(TMath::Pi()/2.) * (1 + erf(alpha/sqrt(2)));
  //Double_t N = 1/(par[2]*(C+D));
  
  Double_t ref = (x[0] - par[3])/par[2];  // (x-mean)/sigma
  if ( par[0] < 0. ) ref = -ref;
  //Double_t fitfunc = N;
  Double_t fitfunc = 1.;
  if ( ref > -alpha ) fitfunc = fitfunc * exp(-ref*ref/2.);
  else if (ref <= -alpha ) fitfunc = fitfunc * A * pow ( B - ref , -par[1]);
  return fitfunc*par[4];
}

Double_t psCrysBall(Double_t *x, Double_t *par) {
  // mu, sigma, k_L, k_H  // low/high shoulder
  if ( par[1] <= 0. ) return 0.;
  Double_t ref = ( x[0] - par[0] ) / par[1];
  Double_t fitfunc = exp( -0.5*pow(ref , 2) );
  if ( ref <= -par[2] )    fitfunc = exp( -0.5*pow(par[2], 2) + par[2]*ref );
  else if ( ref > par[3] ) fitfunc = exp( -0.5*pow(par[3], 2) - par[3]*ref );
  return fitfunc*par[4];
}

Double_t testFit(Double_t *x, Double_t *par) {
  return crysBall(x,par) + lorentzian(x,&par[3]);
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
  
  if ( argc == 2 ) suffix = "_widthx"+(string)argv[1]; cout << "Suffix is " << suffix << endl;
  
  //gSystem->Load("libMathCore.so");
  
  /// Declare histos to be fitted
  const string histDir = "1D_histograms/";
  pair<const string, int> histoNames[] = { /*{"mTop_div_aveMTop_TT_matched_jets", 1},*/ {"mTop_div_aveMTop_TT_reco_CP", 1}, {"mTop_div_aveMTop_TT_reco_WP", 0}, /*{"mTop_div_aveMTop_TT_reco_WP_WOk", 0}, {"mTop_div_aveMTop_TT_reco_WP_WNotOk", 0},*/ {"mTop_div_aveMTop_TT_reco_UP", 2}, {"mTop_div_aveMTop_TT_reco_WPUP", 3} };
  int sizeHistos = sizeof(histoNames)/sizeof(histoNames[0]);
  
  /// Declare input and output files
  string inputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputPlots/mu/"+whichDate(systStr, suffix);
  string pathOutput = "OutputVoigt/ex4jets/170223/";
  mkdir(pathOutput.c_str(),0777);
  string outputFileName = pathOutput+"VoigtFit_"+systStr+suffix+".root";
  cout << "Output file: " << outputFileName << endl;
  
  ofstream fileOut;
  string outputOFileName = pathOutput+"FitParams"+suffix+".txt";
  fileOut.open(outputOFileName.c_str());
  
  TFile *fin = new TFile(inputFileName.c_str(), "READ");
  fin->cd();
  
  TFile *fout = new TFile(outputFileName.c_str(), "RECREATE");
  fout->cd();
  
  
  for (int iHisto = 0; iHisto < sizeHistos; iHisto++)
  {
    //if (! fin->GetListOfKeys()->Contains(("1D_histograms/"+histoNames[iHisto]).c_str()) )  // does not work for histograms in dir
    //{
      //cout << "*** Histogram " << histoNames[iHisto] << " does not exist... Proceeding to next histogram..." << endl;
      //continue;
    //}
    cout << "*** Processing histogram " << histoNames[iHisto].first << endl;
    
    /// Get histos
    TH1F* histo = (TH1F*) fin->Get((histDir+histoNames[iHisto].first).c_str());
    histo->Write();
    
    // Normalise histogram
    Double_t scale = 1./histo->Integral();  // histo->Integral(0, histoIn->GetNbinsX()+1);
    histo->Scale(scale);
    //TH1F *histo = (TH1F*) histoIn->Clone("histo");
    
    double fitMin = -9., fitMax = -9.;
    double baseline = 18e-5;
    if ( histoNames[iHisto].second != 1 ) baseline = 22e-5;
    for (int iBin = 0; iBin < histo->GetNbinsX(); iBin++)
    {
      if ( fitMax != -9. ) break;
      if ( fitMin != -9. && iBin < histo->GetXaxis()->FindBin(1.2) ) continue;
      double binContent = histo->GetBinContent(iBin) - baseline;
      if ( fitMin == -9. && binContent > 0.) fitMin = histo->GetXaxis()->GetBinCenter(iBin);
      if ( fitMin != -9. && fitMax == -9. && binContent < 0. ) fitMax = histo->GetXaxis()->GetBinCenter(iBin);
    }
    if ( fitMin < 0.2 ) fitMin = 0.2;
    if ( fitMax > 2.2 ) fitMax = 2.2;
    cout << "Fit min = " << fitMin << "; fit Max = " << fitMax << endl;
    
    /// Declare fit function
    int nBins = histo->GetXaxis()->GetNbins();
    const int nPar = 5;
    
    TF1 *myfit;
    
    if ( histoNames[iHisto].second == 1 )  // CP
    {
      if (usePredef)
      {
        myfit = new TF1("fit", voigt_predef, fitMin, fitMax, nPar);  // root name, fit function, range, nPar
        
        myfit->SetParNames("#mu","#sigma","#gamma", "r", "norm");
        myfit->SetParameters(1.007, 0.0729, 1.36*histo->GetRMS(), 2.062, 0.002591);  // sigma = 0.57735027 * RMS; FWHM = 1.359556 * RMS (= 2.354820 * sigma)
        
        //myfit->SetParLimits(0, 0.95, 1.05);
        myfit->FixParameter(0, 1.007);
        //myfit->SetParLimits(1, 1e-4, 1.e+3);
        myfit->FixParameter(1, 0.0749);
        myfit->SetParLimits(2, 1e-4, 1.e+3);
        //myfit->SetParLimits(3, 2., 5.);
        myfit->FixParameter(3, 2.031);
        myfit->FixParameter(4, 0.002578);
      }
      else
      {
        myfit = new TF1("fit", voigt, fitMin, fitMax, nPar);  // root name, fit function, range, nPar
        
        myfit->SetParNames("ps","#mu","#sigma","#gamma","norm");
        myfit->SetParameters(0.5, 1., 0.58*histo->GetRMS(), 1.36*histo->GetRMS(), 1.);  // sigma = 0.57735027 * RMS; FWHM = 1.359556 * RMS (= 2.354820 * sigma)
        
        myfit->SetParLimits(0,0.,1.);
        myfit->SetParLimits(1,0.95,1.05);
        myfit->SetParLimits(2,1e-4,1.e+3);
        myfit->SetParLimits(3,1e-4,1.e+3);
      }
    }
    else  // WP/UP
    {
      if (useTest)
      {
        myfit = new TF1("fit", testFit, fitMin, fitMax, 8);  // root name, fit function, range, nPar
        
        myfit->SetParNames("alpha","n","#sigma","#mu","norm","scaleL", "meanL", "sigmaL");
        myfit->SetParameters(-0.59, 17.7, 0.2, 0.8, 1., 0.05, 0.5, 0.0001);
        
        myfit->SetParLimits(1, 0.001, 50);
        myfit->SetParLimits(2, 1.e-4, 1.e+3);
        myfit->SetParLimits(3, 0.70, 0.92);
        myfit->SetParLimits(5, 1e-4, 1.e+3);
        myfit->SetParLimits(6, 0.45, 0.6);
        myfit->SetParLimits(7, 1e-6, 1.e+2);
      }
      else
      {
        myfit = new TF1("fit", crysBall, fitMin, fitMax, nPar);  // root name, fit function, range, nPar
        
        myfit->SetParNames("alpha","n","#sigma","#mu","norm");
        
        myfit->SetParLimits(1, 0.1, 20);
        myfit->SetParLimits(2, 1.e-4, 1.e+3);
        myfit->SetParLimits(3, 0.60, 0.92);
        
        if ( histoNames[iHisto].second == 0 ) // WP
        {
          myfit->SetParameters(-0.457, 12.7, 0.1835, 0.7747, 0.003768);
          myfit->FixParameter(0, -0.52);
          //myfit->SetParLimits(1, 0.1, 30);
          myfit->FixParameter(1, 9.5);
          myfit->FixParameter(2, 0.1835);
          myfit->FixParameter(3, 0.7747);
          myfit->FixParameter(4, 0.003886);
        }
        else  // UP
        {
          myfit->SetParameters(-0.5966, 12., 0.1835, 0.805, 0.004106);
          myfit->FixParameter(0, -0.595);
          myfit->FixParameter(2, 0.1835);
          myfit->FixParameter(3, 0.8055);
          
          if ( histoNames[iHisto].second == 2 ) // UP
          {
            myfit->FixParameter(1, 11.14);
            myfit->FixParameter(4, 0.004097);
          }
          else if ( histoNames[iHisto].second == 3 ) // WPUP
          {
            myfit->FixParameter(1, 9.434);
            myfit->FixParameter(4, 0.004088);
          }
        }
        
      }
    }
    
    
    ///  Fit
    std::string func_title = std::string(histo->GetName())+"_Fitted";
    myfit->SetName(func_title.c_str());
    histo->Fit(myfit, "R");
    gStyle->SetOptFit(0111);
    histo->SetStats(1);
    histo->Write();
    myfit->Write();
    
    
    /// Write fit params to file
    if ( histoNames[iHisto].second == 1 ) fileOut << "Voigt: ";
    else fileOut << "CrystalBall: ";
    fileOut << histoNames[iHisto].first << endl;
    for (int iPar = 0; iPar < nPar; iPar++)
    {
      fileOut << iPar << "   " << myfit->GetParName(iPar) << "   " << myfit->GetParameter(iPar) << "   " << myfit->GetParError(iPar) << endl;
    }
    fileOut << endl;
    
    /// Plot function with fit params over entire range
    TCanvas* c1 = new TCanvas(func_title.c_str(), func_title.c_str());
    c1->cd();
    histo->SetLineColor(kBlue);
    histo->Draw();
    
    TF1 *func;
    if ( histoNames[iHisto].second == 1 )
    {
      if (usePredef) func = new TF1("func", voigt_predef, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(), 5);
      else func = new TF1("func", voigt, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(), 5);
    }
    else
    {  
      func = new TF1("fcb", crysBall, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(), 5);
    }
    func->SetParameters(myfit->GetParameter(0), myfit->GetParameter(1), myfit->GetParameter(2), myfit->GetParameter(3), myfit->GetParameter(4));
    func->SetLineColor(kRed);
    func->Draw("same");
    c1->Update();
    c1->Write();
    c1->SaveAs((pathOutput+func_title+suffix+".png").c_str());
    c1->Close();
    
    if (debug)
      cout << "Integrate histogram: " << histo->Integral() << "; Histo sum .. : " << histo->GetSumOfWeights() << "; Integral function: " << func->Integral(histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(), 0) << endl;
    
    
    delete func;
    delete c1;
    
    delete myfit;
    
  }  /// end loop histos
  
  
  
  
  
  fout->Write();
  
  
  fin->Close();
  fout->Close();
  fileOut.close();
  
  delete fin;
  delete fout;
  
  
  ///// end
  double time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  if ( time >= 60 )
  {
    int mins = time/60;
    double secs = time - mins*60;
    
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
 
