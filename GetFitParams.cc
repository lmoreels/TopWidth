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


bool runLocally = false;


/// Define vars
float width[] = {0.5, 0.66, 0.75, 1., 2., 3., 4.};
string widthStr[] = {"0p5", "0p66", "0p75", "1", "2", "3", "4"};
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
bool PrintParams();


int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "********************************" << endl;
  cout << "***   Get Fit Parameters   ***" << endl;
  cout << "********************************" << endl;
  cout << "* Current time: " << dateString << "    *" << endl;
  cout << "********************************" << endl;
  
  ClearVars();
  ClearVectors();
  
  string pathInput = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputVoigt/ex4jets/";
  if (runLocally) pathInput = "/Users/lmoreels/cernbox/TopWidth/TopTrees/tempPlots/Voigt/161220/";
  
  /// Get fit parameters from files
  for (int iWidth = 0; iWidth < nWidths; iWidth++)
  {
    /// Open input file
    inputFileName = pathInput+"FitParams_widthx"+widthStr[iWidth]+".txt";
    if (! fexists(inputFileName.c_str()) )
    {
      cout << "WARNING: File " << inputFileName << " does not exist." << endl;
      continue;
    }
    fileIn.open(inputFileName.c_str());
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
            nErrCB_WP.push_back(paramError);
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
            nErrCB_UP.push_back(paramError);
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
  
  PrintParams();
  
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

bool PrintParams()
{
  cout << endl << "=== MATHEMATICA ===" << endl;
  cout << endl << "width={" << width[0] << "," << width[1] << "," << width[2] << "," << width[3] << "," << width[4] << "," << width[5] << "," << width[6] << "}; (* times SM width *)" << endl;
  cout << endl << "-- VOIGT CP --" << endl;
  cout << "mu={" << muVoigt[0] << "," << muVoigt[1] << "," << muVoigt[2] << "," << muVoigt[3] << "," << muVoigt[4] << "," << muVoigt[5] << "," << muVoigt[6] << "};" << endl;
  cout << "onzmu={" << muErrVoigt[0] << "," << muErrVoigt[1] << "," << muErrVoigt[2] << "," << muErrVoigt[3] << "," << muErrVoigt[4] << "," << muErrVoigt[5] << "," << muErrVoigt[6] << "};" << endl;
  cout << "sigma={" << sigmaVoigt[0] << "," << sigmaVoigt[1] << "," << sigmaVoigt[2] << "," << sigmaVoigt[3] << "," << sigmaVoigt[4] << "," << sigmaVoigt[5] << "," << sigmaVoigt[6] << "};" << endl;
  cout << "onzsigma={" << sigmaErrVoigt[0] << "," << sigmaErrVoigt[1] << "," << sigmaErrVoigt[2] << "," << sigmaErrVoigt[3] << "," << sigmaErrVoigt[4] << "," << sigmaErrVoigt[5] << "," << sigmaErrVoigt[6] << "};" << endl;
  cout << "gamma={" << gammaVoigt[0] << "," << gammaVoigt[1] << "," << gammaVoigt[2] << "," << gammaVoigt[3] << "," << gammaVoigt[4] << "," << gammaVoigt[5] << "," << gammaVoigt[6] << "};" << endl;
  cout << "onzgamma={" << gammaErrVoigt[0] << "," << gammaErrVoigt[1] << "," << gammaErrVoigt[2] << "," << gammaErrVoigt[3] << "," << gammaErrVoigt[4] << "," << gammaErrVoigt[5] << "," << gammaErrVoigt[6] << "};" << endl;
  cout << "r={" << rVoigt[0] << "," << rVoigt[1] << "," << rVoigt[2] << "," << rVoigt[3] << "," << rVoigt[4] << "," << rVoigt[5] << "," << rVoigt[6] << "};" << endl;
  cout << "onzr={" << rErrVoigt[0] << "," << rErrVoigt[1] << "," << rErrVoigt[2] << "," << rErrVoigt[3] << "," << rErrVoigt[4] << "," << rErrVoigt[5] << "," << rErrVoigt[6] << "};" << endl;
  cout << "norm={" << normVoigt[0] << "," << normVoigt[1] << "," << normVoigt[2] << "," << normVoigt[3] << "," << normVoigt[4] << "," << normVoigt[5] << "," << normVoigt[6] << "};" << endl;
  cout << "onznorm={" << normErrVoigt[0]*1e+6 << "," << normErrVoigt[1]*1e+6 << "," << normErrVoigt[2]*1e+6 << "," << normErrVoigt[3]*1e+6 << "," << normErrVoigt[4]*1e+6 << "," << normErrVoigt[5]*1e+6 << "," << normErrVoigt[6]*1e+6 << "}*10^(-6);" << endl;
  
  cout << endl << "-- CRYSTALBALL WP --" << endl;
  cout << "alpha={" << alphaCB_WP[0] << "," << alphaCB_WP[1] << "," << alphaCB_WP[2] << "," << alphaCB_WP[3] << "," << alphaCB_WP[4] << "," << alphaCB_WP[5] << "," << alphaCB_WP[6] << "};" << endl;
  cout << "onzalpha={" << alphaErrCB_WP[0] << "," << alphaErrCB_WP[1] << "," << alphaErrCB_WP[2] << "," << alphaErrCB_WP[3] << "," << alphaErrCB_WP[4] << "," << alphaErrCB_WP[5] << "," << alphaErrCB_WP[6] << "};" << endl;
  cout << "n={" << nCB_WP[0] << "," << nCB_WP[1] << "," << nCB_WP[2] << "," << nCB_WP[3] << "," << nCB_WP[4] << "," << nCB_WP[5] << "," << nCB_WP[6] << "};" << endl;
  cout << "onzn={" << nErrCB_WP[0] << "," << nErrCB_WP[1] << "," << nErrCB_WP[2] << "," << nErrCB_WP[3] << "," << nErrCB_WP[4] << "," << nErrCB_WP[5] << "," << nErrCB_WP[6] << "};" << endl;
  cout << "sigma={" << sigmaCB_WP[0] << "," << sigmaCB_WP[1] << "," << sigmaCB_WP[2] << "," << sigmaCB_WP[3] << "," << sigmaCB_WP[4] << "," << sigmaCB_WP[5] << "," << sigmaCB_WP[6] << "};" << endl;
  cout << "onzsigma={" << sigmaErrCB_WP[0] << "," << sigmaErrCB_WP[1] << "," << sigmaErrCB_WP[2] << "," << sigmaErrCB_WP[3] << "," << sigmaErrCB_WP[4] << "," << sigmaErrCB_WP[5] << "," << sigmaErrCB_WP[6] << "};" << endl;
  cout << "mu={" << muCB_WP[0] << "," << muCB_WP[1] << "," << muCB_WP[2] << "," << muCB_WP[3] << "," << muCB_WP[4] << "," << muCB_WP[5] << "," << muCB_WP[6] << "};" << endl;
  cout << "onzmu={" << muErrCB_WP[0] << "," << muErrCB_WP[1] << "," << muErrCB_WP[2] << "," << muErrCB_WP[3] << "," << muErrCB_WP[4] << "," << muErrCB_WP[5] << "," << muErrCB_WP[6] << "};" << endl;
  cout << "norm={" << normCB_WP[0] << "," << normCB_WP[1] << "," << normCB_WP[2] << "," << normCB_WP[3] << "," << normCB_WP[4] << "," << normCB_WP[5] << "," << normCB_WP[6] << "};" << endl;
  cout << "onznorm={" << normErrCB_WP[0]*1e+6 << "," << normErrCB_WP[1]*1e+6 << "," << normErrCB_WP[2]*1e+6 << "," << normErrCB_WP[3]*1e+6 << "," << normErrCB_WP[4]*1e+6 << "," << normErrCB_WP[5]*1e+6 << "," << normErrCB_WP[6]*1e+6 << "}*10^(-6);" << endl;
  
  cout << endl << "-- CRYSTALBALL UP --" << endl;
  cout << "alpha={" << alphaCB_UP[0] << "," << alphaCB_UP[1] << "," << alphaCB_UP[2] << "," << alphaCB_UP[3] << "," << alphaCB_UP[4] << "," << alphaCB_UP[5] << "," << alphaCB_UP[6] << "};" << endl;
  cout << "onzalpha={" << alphaErrCB_UP[0] << "," << alphaErrCB_UP[1] << "," << alphaErrCB_UP[2] << "," << alphaErrCB_UP[3] << "," << alphaErrCB_UP[4] << "," << alphaErrCB_UP[5] << "," << alphaErrCB_UP[6] << "};" << endl;
  cout << "n={" << nCB_UP[0] << "," << nCB_UP[1] << "," << nCB_UP[2] << "," << nCB_UP[3] << "," << nCB_UP[4] << "," << nCB_UP[5] << "," << nCB_UP[6] << "};" << endl;
  cout << "onzn={" << nErrCB_UP[0] << "," << nErrCB_UP[1] << "," << nErrCB_UP[2] << "," << nErrCB_UP[3] << "," << nErrCB_UP[4] << "," << nErrCB_UP[5] << "," << nErrCB_UP[6] << "};" << endl;
  cout << "sigma={" << sigmaCB_UP[0] << "," << sigmaCB_UP[1] << "," << sigmaCB_UP[2] << "," << sigmaCB_UP[3] << "," << sigmaCB_UP[4] << "," << sigmaCB_UP[5] << "," << sigmaCB_UP[6] << "};" << endl;
  cout << "onzsigma={" << sigmaErrCB_UP[0] << "," << sigmaErrCB_UP[1] << "," << sigmaErrCB_UP[2] << "," << sigmaErrCB_UP[3] << "," << sigmaErrCB_UP[4] << "," << sigmaErrCB_UP[5] << "," << sigmaErrCB_UP[6] << "};" << endl;
  cout << "mu={" << muCB_UP[0] << "," << muCB_UP[1] << "," << muCB_UP[2] << "," << muCB_UP[3] << "," << muCB_UP[4] << "," << muCB_UP[5] << "," << muCB_UP[6] << "};" << endl;
  cout << "onzmu={" << muErrCB_UP[0] << "," << muErrCB_UP[1] << "," << muErrCB_UP[2] << "," << muErrCB_UP[3] << "," << muErrCB_UP[4] << "," << muErrCB_UP[5] << "," << muErrCB_UP[6] << "};" << endl;
  cout << "norm={" << normCB_UP[0] << "," << normCB_UP[1] << "," << normCB_UP[2] << "," << normCB_UP[3] << "," << normCB_UP[4] << "," << normCB_UP[5] << "," << normCB_UP[6] << "};" << endl;
  cout << "onznorm={" << normErrCB_UP[0]*1e+6 << "," << normErrCB_UP[1]*1e+6 << "," << normErrCB_UP[2]*1e+6 << "," << normErrCB_UP[3]*1e+6 << "," << normErrCB_UP[4]*1e+6 << "," << normErrCB_UP[5]*1e+6 << "," << normErrCB_UP[6]*1e+6 << "}*10^(-6);" << endl;
  
  cout << endl << "=== LIKELIHOOD ===" << endl;
  cout << "const double mu_CP = " << muVoigt[0] << ", sigma_CP = " << sigmaVoigt[0] << ", r_CP = " << rVoigt[0] << ", norm_CP = " << normVoigt[0] << ";" << endl;
  cout << "const double alpha_WP = " << alphaCB_WP[0] << ", n_WP = " << nCB_WP[0] << ", sigma_WP = " << sigmaCB_WP[0] << ", mu_WP = " << muCB_WP[0] << ", norm_WP = " << normCB_WP[0] << ";" << endl;
  cout << "const double alpha_UP = " << alphaCB_UP[0] << ", n_UP = " << nCB_UP[0] << ", sigma_UP = " << sigmaCB_UP[0] << ", mu_UP = " << muCB_UP[0] << ", norm_UP = " << normCB_UP[0] << ";" << endl;
}
