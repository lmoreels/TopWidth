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


bool test = true;
bool runLocally = false;
bool checkNormFunctions = false;
bool printFractions = false;

const int nCP = 354711;  // [0.6, 1.4]  //355056;  // [0.5, 1.5]
const int nWP = 114878;  // [0.6, 1.4]  //123818;  // [0.5, 1.5]
const int nUP = 244934;  // [0.6, 1.4]  //253221;  // [0.5, 1.5] tt only!

const double mu_CP = 1.01, sigma_CP = 0.0665, r_CP = 2.0222, norm_CP = 0.002556;
const double alpha_WP = -0.39, n_WP = 12.95, sigma_WP = 0.185, mu_WP = 0.889, norm_WP = 0.003414;
const double alpha_UP = -1.001, n_UP = 2.4094, sigma_UP = 0.1977, mu_UP = 0.97, norm_UP = 0.004225;
const double alpha_WPUP = -0.8595, n_WPUP = 3.044, sigma_WPUP = 0.2128, mu_WPUP = 0.97, norm_WPUP = 0.003862;
const double gammaConv[2] = {0.0314305, 0.00763251};

double widthArray[] = {0.5, 0.66, 0.75, 1., 2., 3., 4.};
const int nWidths = sizeof(widthArray)/sizeof(widthArray[0]);


/// Define vars
int nTot;
double f_CP, f_WP, f_UP;
TF1 *likelihood;

string inputFileName;
ifstream fileIn;
string var;
double val;
vector<double> widths;
vector<double> LLvalues;

/// Define functions
string ConvertDoubleToString(double Number);
string DotReplace(double var);
string ConvertIntToString(int Number, int pad);
string MakeTimeStamp();
bool fexists(const char *filename);
void ClearVars();
Double_t voigt(Double_t *x, Double_t *par);
Double_t crysBall_WP(Double_t *x, Double_t *par);
Double_t crysBall_UP(Double_t *x, Double_t *par);
Double_t combinedProb(Double_t *x, Double_t *par);
Double_t logLikelihood(Double_t *x, Double_t *par);
Double_t widthToGammaTranslation(Double_t *x);
void DrawFunction(TF1* function, string name, double width, bool writeToFile);
void DrawFunction(TF1* function, string name, std::vector<double> widths, double min, double max, bool writeToFile);
void DrawLikelihood(float width);


int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "********************************" << endl;
  cout << "***  Likelihood calculation  ***" << endl;
  cout << "********************************" << endl;
  cout << "* Current time: " << dateString << "    *" << endl;
  cout << "********************************" << endl;
  
  
  /// Check if functions are normalised
  if (checkNormFunctions)
  {
    double lowerRange = 1e-10, upperRange = 20.;
    double combMin = 1e-10, combMax = 2.;
    
    TF1 *voigt_cp = new TF1("voigt_cp", voigt, combMin, combMax, 1);
    voigt_cp->SetParameter(0,1.5);
    double sf_cp = 1./voigt_cp->Integral(lowerRange, upperRange, 0);
    std::cout << "Integral of the voigt function is " << voigt_cp->Integral(lowerRange, upperRange, 0) << ", so norm scale factor should be " << sf_cp << std::endl;
    TF1 *crysball_wp = new TF1("crysball_wp", crysBall_WP, combMin, combMax, 0);
    double sf_wp = 1./crysball_wp->Integral(lowerRange, upperRange, 0);
    std::cout << "Integral of the crystal ball function (WP) is " << crysball_wp->Integral(lowerRange, upperRange, 0) << ", so norm scale factor should be " << sf_wp << std::endl;
    TF1 *crysball_up = new TF1("crysball_up", crysBall_UP, combMin, combMax, 0);
    double sf_up = 1./crysball_up->Integral(lowerRange, upperRange, 0);
    std::cout << "Integral of the crystal ball function (UP) is " << crysball_up->Integral(lowerRange, upperRange, 0) << ", so norm scale factor should be " << sf_up << std::endl;
    
    TF1 *combi = new TF1("combi", combinedProb, combMin, combMax, 1);
    std::cout << "Integral of the combination is " << combi->Integral(lowerRange, upperRange, 0) << ", so norm scale factor should be " << 1./combi->Integral(lowerRange, upperRange, 0) << std::endl;
    
    /// LIKELIHOOD ///
    // prob = f_CP * Voigt + f_WP * CB + f_UP * CB
    // likelihood = - log (prob)
    
    likelihood = new TF1("likelihood", logLikelihood, combMin, combMax, 1);
    
    // Check if functions are normalised
    std::cout << "Integral of the likelihood is " << likelihood->Integral(lowerRange, upperRange, 0) << std::endl;
    
    
    /// Integrate "by hand"
    double evalin, testwidth = 1.;
    double stepSize = (combMax - combMin)/2000;
    double manualInt_CP = 0, manualInt_WP = 0, manualInt_UP = 0;
    for (int i = 0; i < 2000; i++)
    {
      evalin = i*stepSize;
      manualInt_CP += voigt(&evalin, &testwidth);
      manualInt_WP += crysBall_WP(&evalin, &testwidth);
      manualInt_UP += crysBall_UP(&evalin, &testwidth);
    }
    cout << "Manual integral calculation (2000 points):" << endl;
    cout << "CP: " <<  manualInt_CP << endl;
    cout << "WP: " <<  manualInt_WP << endl;
    cout << "UP: " <<  manualInt_UP << endl;
    
    TFile *foutNorm = new TFile("likelihood/CheckPDFNormalisation.root", "RECREATE");
    foutNorm->cd();
    DrawFunction(voigt_cp, "checkNorm_Voigt", 1., true);
    DrawFunction(crysball_wp, "checkNorm_CB_WP", 1., true);
    DrawFunction(crysball_up, "checkNorm_CB_UP", 1., true);
    DrawFunction(combi, "checkNorm_combi", 1., true);
    DrawFunction(likelihood, "loglikelihood", 0.025, true);
    //vector<double> gammas = {0.030, 0.035, 0.04, 0.045, 0.05, 0.055, 0.060, 0.065};  // C++11 only !
    vector<double> gammas = {0.0349, 0.03654, 0.03739, 0.03899, 0.04764, 0.05472, 0.06279};
    DrawFunction(voigt_cp, "voigt_multipleG_cut", gammas, 0.7, 1.3, true);
    DrawFunction(likelihood, "loglikelihood_multipleG_cut", gammas, combMin, combMax, true);
    vector<double> widths = {0.5, 0.66, 0.75, 1., 2., 3., 4.};
    DrawFunction(voigt_cp, "voigt_multipleW_cut", widths, 0.5, 1.5, true);
    DrawFunction(likelihood, "loglikelihood_multipleW_cut", widths, combMin, combMax, true);
    foutNorm->Close();
    delete foutNorm;
  }
  
  
  
  for (int iWidth = 0; iWidth < nWidths; iWidth++)
  {
    ClearVars();
    
    /// Get loglikelihood values from file
    inputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/output_loglikelihood_widthx"+DotReplace(widthArray[iWidth])+".txt";
    if (! fexists(inputFileName.c_str()) )
    {
      cout << "WARNING: File " << inputFileName << " does not exist." << endl;
      exit(1);
    }
    fileIn.open(inputFileName.c_str());
    cout << "Opening " << inputFileName << "..." << endl;
    
    string line;
    while( getline(fileIn, line) )
    {
      istringstream iss(line);
      if ( line.find("Width") == 0 )
      {
        iss >> var;
        while ( iss >> val )
        {
          widths.push_back(val);
        }
      }
      if ( line.find("LL") == 0 )
      {
        iss >> var;
        while ( iss >> val )
        {
          LLvalues.push_back(val);
        }
      }
    }
    
    fileIn.close();
    
    if (test)
    {
      for (int i = 0; i < widths.size(); i++)
        cout << widths[i] << "  ";
      cout << endl;
    }
    
    /// Make TGraph
    const int nPoints = widths.size();
    double widthArr[nPoints], LLArr[nPoints];
    for (int i = 0; i < nPoints; i++)
    {
      widthArr[i] = widths[i];
      LLArr[i] = LLvalues[i];
    }
    
    TGraph *g1 = new TGraph(widths.size(), widthArr, LLArr);
    
    /// Fit minimum with parabola
    double fitminArray[] = {0.6, 0.6, 0.6, 0.8, 1.0, 1.5, 1.8};
    double fitmaxArray[] = {1.4, 1.5, 1.6, 1.8, 2.5, 3.0, 3.2};
    double fitmin = fitminArray[iWidth], fitmax = fitmaxArray[iWidth];
    TF1 *parabola = new TF1("parabola", "pol2", fitmin, fitmax);
    g1->Fit(parabola,"R");
    
    if (test)
    {
      TCanvas* c1 = new TCanvas("c1", "LLike vs. width");
      c1->cd();
      g1->Draw("AL");
      c1->Update();
      c1->SaveAs(("loglikelihoodVSwidth_"+DotReplace(widthArray[iWidth])+".png").c_str());
      c1->Close();

      delete c1;
    }
    
    double minimum = parabola->GetMinimumX(fitmin, fitmax);
    cout << "For " << widthArray[iWidth] << " times SM width of the top quark, the minimum can be found at " << minimum << endl;
    
    delete parabola;
    delete g1;
    
  }  // end loop widths
  
  
  return 0;
}


/// Functions
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
  var = "";
  val = -1.;
  widths.clear();
  LLvalues.clear();
}

Double_t voigt(Double_t *x, Double_t *par) {
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
  return norm_CP*TMath::Voigt(x[0]-mu_CP, sigma_CP, par[0], r_CP)/*1.268018*/;
}

Double_t crysBall_WP(Double_t *x, Double_t *par) {
  // params: alpha, n, sigma, mu
  if ( sigma_WP <= 0. ) return 0.;
  Double_t alpha = fabs(alpha_WP);
  Double_t A = pow( n_WP/alpha , n_WP) * exp(-alpha*alpha/2.);
  Double_t B = n_WP/alpha - alpha;
  //Double_t C = n_WP/alpha * 1/(n_WP-1) * exp(-alpha*alpha/2.);
  //Double_t D = sqrt(TMath::Pi()/2.) * (1 + erf(alpha/sqrt(2)));
  //Double_t N = 1/(sigma_WP*(C+D));
  
  Double_t ref = (x[0] - mu_WP)/sigma_WP;  // (x-mean)/sigma
  if ( alpha_WP < 0 ) ref = -ref;
  //Double_t fitfunc = N;
  Double_t fitfunc = 1.;
  if ( ref > -alpha ) fitfunc = fitfunc * exp(-ref*ref/2.);
  else if (ref <= -alpha ) fitfunc = fitfunc * A * pow ( B - ref , -n_WP);
  return norm_WP*fitfunc/*1.804797*/;
}

Double_t crysBall_UP(Double_t *x, Double_t *par) {
  // params: alpha, n, sigma, mu
  if ( sigma_UP <= 0. ) return 0.;
  Double_t alpha = fabs(alpha_UP);
  Double_t A = pow( n_UP/alpha , n_UP) * exp(-alpha*alpha/2.);
  Double_t B = n_UP/alpha - alpha;
  //Double_t C = n_UP/alpha * 1/(n_UP-1) * exp(-alpha*alpha/2.);
  //Double_t D = sqrt(TMath::Pi()/2.) * (1 + erf(alpha/sqrt(2)));
  //Double_t N = 1/(sigma_UP*(C+D));
  
  Double_t ref = (x[0] - mu_UP)/sigma_UP;  // (x-mean)/sigma
  if ( alpha_UP < 0 ) ref = -ref;
  //Double_t fitfunc = N;
  Double_t fitfunc = 1.;
  if ( ref > -alpha ) fitfunc = fitfunc * exp(-ref*ref/2.);
  else if (ref <= -alpha ) fitfunc = fitfunc * A * pow ( B - ref , -n_UP);
  return norm_UP*fitfunc/*1.6098745*/;
}

Double_t combinedProb(Double_t *x, Double_t *par) {
  /// Fraction of events that is CP, WP and UP
  nTot = nCP + nWP + nUP;
  f_CP = (double)nCP/(double)nTot;
  f_WP = (double)nWP/(double)nTot;
  f_UP = (double)nUP/(double)nTot;
  if (printFractions)
  {
    cout << "f_CP = " << f_CP << "; f_WP = " << f_WP << "; f_UP = " << f_UP << "; Total = " << f_CP+f_WP+f_UP << endl;
    printFractions = false;
  }
  
  return f_CP*voigt(x, par) + f_WP*crysBall_WP(x, 0) + f_UP*crysBall_UP(x, 0) ;
}

Double_t logLikelihood(Double_t *x, Double_t *par) {
  return -TMath::Log(combinedProb(x, par));
}

Double_t widthToGammaTranslation(Double_t *x) {
  return gammaConv[0] + gammaConv[1] * x[0];
}

void DrawFunction(TF1* function, string name, double width, bool writeToFile)
{
  TCanvas* c2 = new TCanvas(name.c_str(), name.c_str());
  c2->cd();
  function->FixParameter(0,width);
  function->Draw();
  c2->Update();
  if (writeToFile) c2->Write();
  c2->SaveAs(("likelihood/"+name+".png").c_str());
  c2->Close();
  
  delete c2;
}

void DrawFunction(TF1* function, string name, std::vector<double> widths, double min, double max, bool writeToFile)
{
  Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+1, kCyan+1, kBlue+2, kMagenta, kViolet-5, kPink+10};
  TCanvas* c2 = new TCanvas(name.c_str(), name.c_str());
  c2->cd();
  function->FixParameter(0,widths[0]);
  function->SetLineColor(colours[0]);
  function->SetRange(min, max);
  function->DrawCopy();
  c2->Update();
  for (int i = 1; i < widths.size(); i++)
  {
    function->FixParameter(0,widths[i]);
    function->SetLineColor(colours[i%10]);
    function->DrawCopy("same");
    c2->Update();
  }
  if (writeToFile) c2->Write();
  c2->SaveAs(("likelihood/"+name+".png").c_str());
  c2->Close();
  
  delete c2;
}

void DrawLikelihood(float width)
{
  TCanvas* c3 = new TCanvas("c3", "Probability function");
  c3->cd();
  likelihood->FixParameter(0,width);
  //likelihood->SetLineColor(kGreen);
  likelihood->Draw();
  c3->Update();
  //c3->Write();
  c3->SaveAs("loglikelihood.png");
  c3->Close();
  
  delete c3;
}
