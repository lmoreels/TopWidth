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


bool debug = false;
bool runLocally = false;
string systStr = "nominal";
string suffix = "";
string widthArray[] = {"0p5", "0p66", "0p75", "1", "2", "3", "4"};
const int nWidths = sizeof(widthArray)/sizeof(widthArray[0]);

string whichDate(string syst, string suff)
{
  if ( syst.find("nominal") == 0 ) 
  {
    if ( suff.find("widthx0p5") != std::string::npos )       return "170331_1259/NtuplePlots_nominal.root";
    else if ( suff.find("widthx0p66") != std::string::npos ) return "170331_1308/NtuplePlots_nominal.root";
    else if ( suff.find("widthx0p75") != std::string::npos ) return "170331_1312/NtuplePlots_nominal.root";
    else if ( suff.find("widthx1") != std::string::npos )    return "170331_1320/NtuplePlots_nominal.root";
    else if ( suff.find("widthx2") != std::string::npos )    return "170331_1324/NtuplePlots_nominal.root";
    else if ( suff.find("widthx3") != std::string::npos )    return "170331_1328/NtuplePlots_nominal.root";
    else if ( suff.find("widthx4") != std::string::npos )    return "170331_1333/NtuplePlots_nominal.root";
    else                                                     return "170331_1320/NtuplePlots_nominal.root";
    
  }
  //else if ( syst.find("JERup") == 0 ) return "161116_1401/NtuplePlots_JERup.root";
  //else if ( syst.find("JERdown") == 0 ) return "161116_1444/NtuplePlots_JERdown.root";
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    suffix = "";
    return "170331_1320/NtuplePlots_nominal.root";
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
  const string histoName = "mTop_div_aveMTop_TT_matched_partons";
  
  string pathOutput = "OutputVoigt/ex4jets/gen/170331/";
  if (runLocally) pathOutput = "/Users/lmoreels/cernbox/TopWidth/TopTrees/tempPlots/Voigt/Ex4jets/gen/170331/";
  mkdir(pathOutput.c_str(),0777);
  
  for (int iWidth = 0; iWidth < nWidths; iWidth++)
  {
    suffix = "_widthx"+widthArray[iWidth];
    cout << "Fitting width " << widthArray[iWidth] << endl;
    
    /// Declare input and output files
    string inputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputPlots/mu/"+whichDate(systStr, suffix);
    if (runLocally) inputFileName = "/Users/lmoreels/cernbox/TopWidth/TopTrees/tempPlots/"+whichDate(systStr, suffix);
    
    string outputFileName = pathOutput+"VoigtFit_"+systStr+suffix+".root";
    cout << "Output file: " << outputFileName << endl;
    
    ofstream fileOut;
    string outputOFileName = pathOutput+"FitParams"+suffix+".txt";
    fileOut.open(outputOFileName.c_str());
    
    TFile *fin = new TFile(inputFileName.c_str(), "READ");
    fin->cd();
    
    TFile *fout = new TFile(outputFileName.c_str(), "RECREATE");
    fout->cd();
    
    
    /// Get histos
    TH1F* histo = (TH1F*) fin->Get((histDir+histoName).c_str());
    histo->Write();
    
    // Normalise histogram
    Double_t scale = 1./histo->Integral();  // histo->Integral(0, histoIn->GetNbinsX()+1);
    histo->Scale(scale);
    //TH1F *histo = (TH1F*) histoIn->Clone("histo");
    
    double fitMin = -9., fitMax = -9.;
    double baseline = 6e-4;
    for (int iBin = 0; iBin < histo->GetNbinsX(); iBin++)
    {
      if ( fitMax != -9. ) break;
      if ( fitMin != -9. && iBin < histo->GetXaxis()->FindBin(1.02) ) continue;
      double binContent = histo->GetBinContent(iBin) - baseline;
      if ( fitMin == -9. && binContent > 0.) fitMin = histo->GetXaxis()->GetBinCenter(iBin);
      if ( fitMin != -9. && fitMax == -9. && binContent < 0. ) fitMax = histo->GetXaxis()->GetBinCenter(iBin-1);
    }
    if ( fitMin < 0.93 ) fitMin = 0.93;
    if ( fitMax > 1.07 ) fitMax = 1.07;
    cout << "Fit min = " << fitMin << "; fit Max = " << fitMax << endl;
    
    /// Declare fit function
    int nBins = histo->GetXaxis()->GetNbins();
    const int nPar = 5;
    
    TF1 *myfit = new TF1("fit", voigt_predef, fitMin, fitMax, nPar);  // root name, fit function, range, nPar
    
    myfit->SetParNames("#mu","#sigma","#gamma", "r", "norm");
    myfit->SetParameters(1.002, 0.0005, 1.36*histo->GetRMS(), 3.5, 0.0004);  // sigma = 0.57735027 * RMS; FWHM = 1.359556 * RMS (= 2.354820 * sigma)
    
    //myfit->SetParLimits(0, 0.99, 1.01);
    myfit->FixParameter(0, 1.002);
    //myfit->SetParLimits(1, 1e-4, 1e-1);
    myfit->FixParameter(1, 0.0005);
    myfit->SetParLimits(2, 1e-5, 1e-1);
    //myfit->SetParLimits(3, 2., 5.);
    myfit->FixParameter(3, 4.4);
    myfit->FixParameter(4, 0.00051);
    
    
    ///  Fit
    std::string func_title = std::string(histo->GetName())+"_Fitted";
    myfit->SetName(func_title.c_str());
    histo->Fit(myfit, "R");
    gStyle->SetOptFit(0111);
    histo->SetStats(1);
    histo->Write();
    myfit->Write();
    
    
    /// Write fit params to file
    fileOut << "Voigt: " << histoName << endl;
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
    
    TF1 *func = new TF1("func", voigt_predef, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(), 5);
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
    
    
    fout->Write();
    
    fin->Close();
    fout->Close();
    fileOut.close();
    
    delete fin;
    delete fout;
    
  }  /// end loop widths
  
  
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
 
