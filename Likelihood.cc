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
bool checkNormFunctions = true;
bool printFractions = true;

const int nCP = 424007; //424657;
const int nWP = 177160; //217980;
const int nUP = 773143; //959425;

const double mu_CP = 1.007, sigma_CP = 0.0729, r_CP = 2.062, norm_CP = 0.002591;
const double alpha_WP = -0.457, n_WP = 20., sigma_WP = 0.1835, mu_WP = 0.7747, norm_WP = 0.003768;
const double alpha_UP = -0.5966, n_UP = 12., sigma_UP = 0.1835, mu_UP = 0.805, norm_UP = 0.004106;

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
string ConvertIntToString(int Number, int pad);
string MakeTimeStamp();
bool fexists(const char *filename);
void ClearVars();
Double_t voigt(Double_t *x, Double_t *par);
Double_t crysBall_WP(Double_t *x, Double_t *par);
Double_t crysBall_UP(Double_t *x, Double_t *par);
Double_t combinedProb(Double_t *x, Double_t *par);
Double_t logLikelihood(Double_t *x, Double_t *par);
void DrawFunction(TF1* function, string name, double width, bool writeToFile);
void DrawFunction(TF1* function, string name, bool writeToFile, double width1, double width2, double width3, double width4, double width5, double width6, double width7/*, double width8, double width9, double width10*/);
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
    
    TFile *foutNorm = new TFile("CheckPDFNormalisation.root", "RECREATE");
    foutNorm->cd();
    DrawFunction(voigt_cp, "checkNorm_Voigt", 1., true);
    DrawFunction(crysball_wp, "checkNorm_CB_WP", 1., true);
    DrawFunction(crysball_up, "checkNorm_CB_UP", 1., true);
    DrawFunction(combi, "checkNorm_combi", 1., true);
    DrawFunction(likelihood, "loglikelihood", 0.025, true);
    //DrawFunction(voigt_cp, "voigt_multipleW", true, 0.5, 0.66, 0.75, 1., 2., 3., 4.);
    DrawFunction(voigt_cp, "voigt_multipleG_cut", true, 0.5, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80);
    //DrawFunction(likelihood, "loglikelihood_multipleW", true, 0.5, 0.66, 0.75, 1., 2., 3., 4.);
    DrawFunction(likelihood, "loglikelihood_multipleG_cut", true, 0.5, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80);
    foutNorm->Close();
    delete foutNorm;
  }
  
  
  ClearVars();
  
  /// Get loglikelihood values from file
  inputFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/output_loglikelihood_gamma.txt";
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
  double fitmin = 0.8, fitmax = 1.5;
  TF1 *parabola = new TF1("parabola", "pol2", fitmin, fitmax);
  g1->Fit(parabola,"R");
  
  if (test)
  {
    TCanvas* c1 = new TCanvas("c1", "LLike vs. width");
    c1->cd();
    g1->Draw("AL");
    c1->Update();
    c1->SaveAs("loglikelihoodVSwidth.png");
    c1->Close();
    
    delete c1;
  }
  
  double minimum = parabola->GetMinimumX(fitmin, fitmax);
  cout << "The minimum can be found at " << minimum << " times the SM width of the top quark." << endl;
  
  
  
  return 0;
}


/// Functions
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

void DrawFunction(TF1* function, string name, double width, bool writeToFile)
{
  TCanvas* c2 = new TCanvas(name.c_str(), name.c_str());
  c2->cd();
  function->FixParameter(0,width);
  function->Draw();
  c2->Update();
  if (writeToFile) c2->Write();
  c2->SaveAs((name+".png").c_str());
  c2->Close();
  
  delete c2;
}

void DrawFunction(TF1* function, string name, bool writeToFile, double width1, double width2, double width3, double width4, double width5, double width6, double width7/*, double width8, double width9, double width10*/)
{
  vector<double> widths;
  if ( &width1 != 0 ) widths.push_back(width1);
  if ( &width2 != 0 ) widths.push_back(width2);
  if ( &width3 != 0 ) widths.push_back(width3);
  if ( &width4 != 0 ) widths.push_back(width4);
  if ( &width5 != 0 ) widths.push_back(width5);
  if ( &width6 != 0 ) widths.push_back(width6);
  if ( &width7 != 0 ) widths.push_back(width7);
/*  if ( &width8 != 0 ) widths.push_back(width8);
  if ( &width9 != 0 ) widths.push_back(width9);
  if ( &width10 != 0 ) widths.push_back(width10);*/
  
  Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+2, kCyan+1, kBlue+2, kMagenta};
  TCanvas* c2 = new TCanvas(name.c_str(), name.c_str());
  c2->cd();
  function->FixParameter(0,widths[0]);
  function->SetLineColor(colours[0]);
  function->DrawCopy();
  c2->Update();
  for (int i = 1; i < widths.size(); i++)
  {
    function->FixParameter(0,widths[i]);
    function->SetLineColor(colours[i]);
    function->DrawCopy("same");
    c2->Update();
  }
  if (writeToFile) c2->Write();
  c2->SaveAs((name+".png").c_str());
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
