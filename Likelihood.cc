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


/// Define functions
string ConvertIntToString(int Number, int pad);
string MakeTimeStamp();
bool ClearVars();
Double_t voigt(Double_t *x, Double_t *par);
Double_t crysBall(Double_t *x, Double_t *par);


int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "********************************" << endl;
  cout << "***  Likelihood calculation  ***" << endl;
  cout << "********************************" << endl;
  cout << "* Current time: " << dateString << "    *" << endl;
  cout << "********************************" << endl;
  
  ClearVars();
  
  /// Fraction of events that is CP, WP and UP
  nTot = nCP + nWP + nUP;
  f_CP = (double)nCP/(double)nTot;
  f_WP = (double)nWP/(double)nTot;
  f_UP = (double)nUP/(double)nTot;
  
  
  /// LIKELIHOOD ///
  // prob = f_CP * Voigt + f_WP * CB + f_UP * CB
  // likelihood = - log (prob)
  
  //TF1 *function = new TF1("function", func, -1, 3, 15);
  //function->SetParameters(mu_CP, sigma_CP, 1., r_CP, norm_CP, alpha_WP, n_WP, sigma_WP, mu_WP, norm_WP, alpha_UP, n_UP, sigma_UP, mu_UP, norm_UP);
  
  //TFormula *function = new TFormula("function", f_CP*Voigt());
  
  TF1 *crysBall_WP = new TF1("crysBall_WP",crysBall);
  crysBall_WP->SetParameters(alpha_WP, n_WP, sigma_WP, mu_WP, norm_WP);
  TF1 *crysBall_UP = new TF1("crysBall_UP",crysBall);
  crysBall_UP->SetParameters(alpha_UP, n_UP, sigma_UP, mu_UP, norm_UP);
  
  TF1 *voigt_CP = new TF1("voigt_CP",voigt);
  voigt_CP->SetParameter(0, mu_CP);
  voigt_CP->SetParameter(1, sigma_CP);
  voigt_CP->SetParameter(3, r_CP);
  voigt_CP->SetParameter(4, norm_CP);
  
  TF1 *likelihood = new TF1("likelihood", "-TMath::Log( f_CP*voigt_CP + f_WP*crysBall_WP + f_UP*crysBall_UP )" );
  
  cout << "Created log likelihood function!" << endl;
  
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
  return TMath::Voigt(x[0]-par[0], par[1], par[2], par[3])*par[4];
}

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
  if ( par[0] < 0 ) ref = -ref;
  //Double_t fitfunc = N;
  Double_t fitfunc = 1.;
  if ( ref > -alpha ) fitfunc = fitfunc * exp(-ref*ref/2.);
  else if (ref <= -alpha ) fitfunc = fitfunc * A * pow ( B - ref , -par[1]);
  return fitfunc*par[4];
}
