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
#include <TCanvas.h>

using namespace std;


bool test = true;
bool runLocally = false;

const int nCP = 588412;
const int nWP = 1064848;
const int nUP = 1648250;

const double mu_CP = 0.9984, sigma_CP = 0.08913, r_CP = 2.47, norm_CP = 0.002558;
const double alpha_WP = -0.3614, n_WP = 20, sigma_WP = 0.1278, mu_WP = 0.7436, norm_WP = 0.004463;
const double alpha_UP = -0.3639, n_UP = 20, sigma_UP = 0.1439, mu_UP = 0.7167, norm_UP = 0.003993;

/// Define vars
int nTot;
double f_CP, f_WP, f_UP;
TF1 *likelihood;

/// Define functions
string ConvertIntToString(int Number, int pad);
string MakeTimeStamp();
bool ClearVars();
Double_t voigt(Double_t *x, Double_t *par);
Double_t crysBall_WP(Double_t *x, Double_t *par);
Double_t crysBall_UP(Double_t *x, Double_t *par);
Double_t combinedProb(Double_t *x, Double_t *par);
Double_t logLikelihood(Double_t *x, Double_t *par);
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
  TF1 *voigt_cp = new TF1("voigt_cp", voigt, 1e-10, 4, 1);
  voigt_cp->SetParameter(0,1.5);
  std::cout << "Integral of the voigt function is " << voigt_cp->Integral(1e-10,1.9e+3,0) << ", so norm scale factor should be " << 1./voigt_cp->Integral(1e-10,1.9e+3,0) << std::endl;
  TF1 *crysball_wp = new TF1("crysball_wp", crysBall_WP, 1e-10, 4, 0);
  std::cout << "Integral of the crystal ball function (WP) is " << crysball_wp->Integral(1e-10,1.9e+3,0) << ", so norm scale factor should be " << 1./crysball_wp->Integral(1e-10,1.9e+3,0) << std::endl;
  TF1 *crysball_up = new TF1("crysball_up", crysBall_UP, 1e-10, 4, 0);
  std::cout << "Integral of the crystal ball function (UP) is " << crysball_up->Integral(1e-10,1.9e+3,0) << ", so norm scale factor should be " << 1./crysball_up->Integral(1e-10,1.9e+3,0) << std::endl;
  
  TF1 *combi = new TF1("combi", combinedProb, 1e-10, 4, 1);
  std::cout << "Integral of the combination is " << combi->Integral(1e-10,1.9e+3,0) << ", so norm scale factor should be " << 1./combi->Integral(1e-10,1.9e+3,0) << std::endl;
  
  /// LIKELIHOOD ///
  // prob = f_CP * Voigt + f_WP * CB + f_UP * CB
  // likelihood = - log (prob)
  
  likelihood = new TF1("likelihood", logLikelihood, 1e-10, 2, 1);
  
  // Check if functions are normalised
  std::cout << "Integral of the likelihood is " << likelihood->Integral(1e-10,1.9e+3,0) << std::endl;
  
  // Now sum over all x_i to get proper (-)log likelihood
  
  if (test) DrawLikelihood(1.5);
 
  
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

bool ClearVars()
{
  nTot = -1.;
  f_CP = -1.;
  f_WP = -1.;
  f_UP = -1.;

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
  return TMath::Voigt(x[0]-mu_CP, sigma_CP, par[0], r_CP)*1.260434;
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
  return fitfunc*1.8048;
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
  return fitfunc*1.609878;
}

Double_t combinedProb(Double_t *x, Double_t *par) {
  /// Fraction of events that is CP, WP and UP
  nTot = nCP + nWP + nUP;
  f_CP = (double)nCP/(double)nTot;
  f_WP = (double)nWP/(double)nTot;
  f_UP = (double)nUP/(double)nTot;
  
  return f_CP*voigt(x, par) + f_WP*crysBall_WP(x, 0) + f_UP*crysBall_UP(x, 0);
}

Double_t logLikelihood(Double_t *x, Double_t *par) {
  return -TMath::Log(combinedProb(x, par));
}

void DrawLikelihood(float width) {
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
