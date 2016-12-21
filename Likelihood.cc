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

/// Define vars
float width[] = {0.25, 0.33, 0.5, 1., 2., 3., 4.};
string widthStr[] = {"0p25", "0p33", "0p5", "1", "2", "3", "4"};
int nWidths = sizeof(widthStr)/sizeof(widthStr[0]);
vector<double> muVoigt, sigmaVoigt, gammaVoigt, rVoigt, normVoigt;
vector<double> muErrVoigt, sigmaErrVoigt, gammaErrVoigt, rErrVoigt, normErrVoigt;
vector<double> alphaCB_WP, nCB_WP, sigmaCB_WP, muCB_WP, normCB_WP;
vector<double> alphaErrCB_WP, nErrCB_WP, sigmaErrCB_WP, muErrCB_WP, normErrCB_WP;
vector<double> alphaCB_UP, nCB_UP, sigmaCB_UP, muCB_UP, normCB_UP;
vector<double> alphaErrCB_UP, nErrCB_UP, sigmaErrCB_UP, muErrCB_UP, normErrCB_UP;
int paramNb;
string paramName;
double paramValue, paramError;

char dataLine[1024];
string inputFileName;
ifstream fileIn;
streampos currentPosition;

/// Define functions
string ConvertIntToString(int Number, int pad);
string MakeTimeStamp();
bool fexists(const char *filename);	// check if file exists
bool ClearVectors();
bool ClearVars();
Double_t voigt(Double_t *x, Double_t *par)
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
  ClearVectors();
  
  string pathInput = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputVoigt/";
  if (runLocally) pathInput = "/Users/lmoreels/cernbox/TopWidth/TopTrees/tempPlots/Voigt/161220/";
  
  /// Get fit parameters from files
  for (int iWidth = 0; iWidth < nWidths; iWidth++)
  {
    /// Open input file
    inputFileName = "FitParams_widthx"+widthStr[iWidth]+".txt";
    if (! fexists((pathInput+inputFileName).c_str()) )
    {
      cout << "WARNING: File " << inputFileName << " does not exist." << endl;
      continue;
    }
    fileIn.open((pathInput+inputFileName).c_str());
    cout << "Opening " << inputFileName << "..." << endl;
    fileIn.getline(dataLine,sizeof(dataLine));
    currentPosition = fileIn.tellg();
    
    /// Loop over input file
    do // until file ends
    {
      //cout << dataLine << endl;
      if ( /*string(dataLine).find("Voigt") != std::string::npos && */ string(dataLine).find("CP") != std::string::npos )
      {
        cout << "  - " << dataLine << ": Adding parameter values for width " << widthStr[iWidth] << "." << endl;
        ClearVars();
        
        fileIn.seekg(currentPosition);
        fileIn.getline(dataLine,sizeof(dataLine));
        currentPosition = fileIn.tellg();
        
        do // until empty line
        {
          istringstream iss(dataLine);
          iss >> paramNb >> paramName >> paramValue >> paramError;
          if ( paramNb == 0 ) // mu
          {
            muVoigt   .push_back(paramValue);
            muErrVoigt.push_back(paramError);
          }
          if ( paramNb == 1 ) // sigma
          {
            sigmaVoigt   .push_back(paramValue);
            sigmaErrVoigt.push_back(paramError);
          }
          if ( paramNb == 2 ) // gamma
          {
            gammaVoigt   .push_back(paramValue);
            gammaErrVoigt.push_back(paramError);
          }
          if ( paramNb == 3 ) // r
          {
            rVoigt   .push_back(paramValue);
            rErrVoigt.push_back(paramError);
          }
          if ( paramNb == 4 ) // norm
          {
            normVoigt   .push_back(paramValue);
            normErrVoigt.push_back(paramError);
          }
          
          fileIn.seekg(currentPosition);
          fileIn.getline(dataLine,sizeof(dataLine));
          currentPosition = fileIn.tellg();
          
        } while (! string(dataLine).empty() );
        
      }  // end Voigt CP
      else if ( string(dataLine).find("CrystalBall") != std::string::npos && string(dataLine).find("WP") != std::string::npos && string(dataLine).find("UP") == std::string::npos && string(dataLine).find("Ok") == std::string::npos)
      {
        cout << "  - " << dataLine << ": Adding parameter values for width " << widthStr[iWidth] << "." << endl;
        ClearVars();
        
        fileIn.seekg(currentPosition);
        fileIn.getline(dataLine,sizeof(dataLine));
        currentPosition = fileIn.tellg();
        
        do // until empty line
        {
          istringstream iss(dataLine);
          iss >> paramNb >> paramName >> paramValue >> paramError;
          if ( paramNb == 0 ) // alpha
          {
            alphaCB_WP   .push_back(paramValue);
            alphaErrCB_WP.push_back(paramError);
          }
          if ( paramNb == 1 ) // n
          {
            nCB_WP   .push_back(paramValue);
            nCB_WP.push_back(paramError);
          }
          if ( paramNb == 2 ) // sigma
          {
            sigmaCB_WP   .push_back(paramValue);
            sigmaErrCB_WP.push_back(paramError);
          }
          if ( paramNb == 3 ) // mu
          {
            muCB_WP   .push_back(paramValue);
            muErrCB_WP.push_back(paramError);
          }
          if ( paramNb == 4 ) // norm
          {
            normCB_WP   .push_back(paramValue);
            normErrCB_WP.push_back(paramError);
          }
          
          fileIn.seekg(currentPosition);
          fileIn.getline(dataLine,sizeof(dataLine));
          currentPosition = fileIn.tellg();
          
        } while (! string(dataLine).empty() );
        
      }  // end CB WP
      else if ( string(dataLine).find("CrystalBall") != std::string::npos && string(dataLine).find("UP") != std::string::npos && string(dataLine).find("WP") == std::string::npos  )
      {
        cout << "  - " << dataLine << ": Adding parameter values for width " << widthStr[iWidth] << "." << endl;
        ClearVars();
        
        fileIn.seekg(currentPosition);
        fileIn.getline(dataLine,sizeof(dataLine));
        currentPosition = fileIn.tellg();
        
        do // until empty line
        {
          istringstream iss(dataLine);
          iss >> paramNb >> paramName >> paramValue >> paramError;
          if ( paramNb == 0 ) // alpha
          {
            alphaCB_UP   .push_back(paramValue);
            alphaErrCB_UP.push_back(paramError);
          }
          if ( paramNb == 1 ) // n
          {
            nCB_UP   .push_back(paramValue);
            nCB_UP.push_back(paramError);
          }
          if ( paramNb == 2 ) // sigma
          {
            sigmaCB_UP   .push_back(paramValue);
            sigmaErrCB_UP.push_back(paramError);
          }
          if ( paramNb == 3 ) // mu
          {
            muCB_UP   .push_back(paramValue);
            muErrCB_UP.push_back(paramError);
          }
          if ( paramNb == 4 ) // norm
          {
            normCB_UP   .push_back(paramValue);
            normErrCB_UP.push_back(paramError);
          }
          
          fileIn.seekg(currentPosition);
          fileIn.getline(dataLine,sizeof(dataLine));
          currentPosition = fileIn.tellg();
          
        } while (! string(dataLine).empty() );
        
      }  // end CB UP
      
      fileIn.seekg(currentPosition);
      fileIn.getline(dataLine,sizeof(dataLine));
      currentPosition = fileIn.tellg();
    
    } while ( fileIn.good() );
    
    /// Close input file
    fileIn.close();
    
  }  // end loop widths (to get fit params)
  
  if (test)
  {
    cout << "Widths: " << width[0] << "  " << width[1] << "  " << width[2] << "  " << width[3] << "  " << width[4] << "  " << width[5] << "  " << width[6] << endl;
    cout << "muVoigts: " << muVoigt[0] << "  " << muVoigt[1] << "  " << muVoigt[2] << "  " << muVoigt[3] << "  " << muVoigt[4] << "  " << muVoigt[5] << "  " << muVoigt[6] << endl;
    cout << "muErrVoigts: " << muErrVoigt[0] << "  " << muErrVoigt[1] << "  " << muErrVoigt[2] << "  " << muErrVoigt[3] << "  " << muErrVoigt[4] << "  " << muErrVoigt[5] << "  " << muErrVoigt[6] << endl;
    cout << "muCB_WPs: " << muCB_WP[0] << "  " << muCB_WP[1] << "  " << muCB_WP[2] << "  " << muCB_WP[3] << "  " << muCB_WP[4] << "  " << muCB_WP[5] << "  " << muCB_WP[6] << endl;
    cout << "muErrCB_WPs: " << muErrCB_WP[0] << "  " << muErrCB_WP[1] << "  " << muErrCB_WP[2] << "  " << muErrCB_WP[3] << "  " << muErrCB_WP[4] << "  " << muErrCB_WP[5] << "  " << muErrCB_WP[6] << endl;
    cout << "muCB_UPs: " << muCB_UP[0] << "  " << muCB_UP[1] << "  " << muCB_UP[2] << "  " << muCB_UP[3] << "  " << muCB_UP[4] << "  " << muCB_UP[5] << "  " << muCB_UP[6] << endl;
    cout << "muErrCB_UPs: " << muErrCB_UP[0] << "  " << muErrCB_UP[1] << "  " << muErrCB_UP[2] << "  " << muErrCB_UP[3] << "  " << muErrCB_UP[4] << "  " << muErrCB_UP[5] << "  " << muErrCB_UP[6] << endl;
  }
  
  /// LIKELIHOOD ///
  // prob = f_CP * Voigt + f_WP * CB + f_UP * CB
  // likelihood = - log (prob)
  
  
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
  return ifile;
}

bool ClearVectors()
{
  muVoigt.clear(); muErrVoigt.clear();
  sigmaVoigt.clear(); sigmaErrVoigt.clear();
  gammaVoigt.clear(); gammaErrVoigt.clear();
  rVoigt.clear(); rErrVoigt.clear();
  normVoigt.clear(); normErrVoigt.clear();
  
  alphaCB_WP.clear(); alphaErrCB_WP.clear();
  nCB_WP.clear(); nErrCB_WP.clear();
  sigmaCB_WP.clear(); sigmaErrCB_WP.clear();
  muCB_WP.clear(); muErrCB_WP.clear();
  normCB_WP.clear(); normErrCB_WP.clear();
  
  alphaCB_UP.clear(); alphaErrCB_UP.clear();
  nCB_UP.clear(); nErrCB_UP.clear();
  sigmaCB_UP.clear(); sigmaErrCB_UP.clear();
  muCB_UP.clear(); muErrCB_UP.clear();
  normCB_UP.clear(); normErrCB_UP.clear();
}

bool ClearVars()
{
  paramNb = -1;
  paramName = "";
  paramValue = 9999.;
  paramError = 9999.;

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
