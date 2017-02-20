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
vector<double> alphaCB_WPUP, nCB_WPUP, sigmaCB_WPUP, muCB_WPUP, normCB_WPUP;
vector<double> alphaErrCB_WPUP, nErrCB_WPUP, sigmaErrCB_WPUP, muErrCB_WPUP, normErrCB_WPUP;
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
void writeVoigtMathematicaOutput(std::vector<double> mu, std::vector<double> sigma, std::vector<double> gamma, std::vector<double> r, std::vector<double> norm, std::vector<double> muErr, std::vector<double> sigmaErr, std::vector<double> gammaErr, std::vector<double> rErr, std::vector<double> normErr);
void writeCBMathematicaOutput(std::vector<double> alpha, std::vector<double> n, std::vector<double> sigma, std::vector<double> mu, std::vector<double> norm, std::vector<double> alphaErr, std::vector<double> nErr, std::vector<double> sigmaErr, std::vector<double> muErr, std::vector<double> normErr);

int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "*********************************" << endl;
  cout << "***    Get Fit Parameters     ***" << endl;
  cout << "*********************************" << endl;
  cout << "*   Current time: " << dateString << "   *" << endl;
  cout << "*********************************" << endl;
  
  ClearVars();
  ClearVectors();
  
  string pathInput = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/OutputVoigt/ex4jets/170220/";
  if (runLocally) pathInput = "/Users/lmoreels/cernbox/TopWidth/TopTrees/tempPlots/Voigt/170220/";
  
  /// Get fit parameters from files
  for (int iWidth = 0; iWidth < nWidths; iWidth++)
  {
    /// Open input file
    inputFileName = pathInput+"FitParams_widthx"+widthStr[iWidth]+".txt";
    if (! fexists(inputFileName.c_str()) )
    {
      cerr << "WARNING: File " << inputFileName << " not found..." << endl;
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
      else if ( string(dataLine).find("CrystalBall") != std::string::npos && string(dataLine).find("WPUP") != std::string::npos )
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
            alphaCB_WPUP   .push_back(paramValue);
            alphaErrCB_WPUP.push_back(paramError);
          }
          if ( paramNb == 1 ) // n
          {
            nCB_WPUP   .push_back(paramValue);
            nErrCB_WPUP.push_back(paramError);
          }
          if ( paramNb == 2 ) // sigma
          {
            sigmaCB_WPUP   .push_back(paramValue);
            sigmaErrCB_WPUP.push_back(paramError);
          }
          if ( paramNb == 3 ) // mu
          {
            muCB_WPUP   .push_back(paramValue);
            muErrCB_WPUP.push_back(paramError);
          }
          if ( paramNb == 4 ) // norm
          {
            normCB_WPUP   .push_back(paramValue);
            normErrCB_WPUP.push_back(paramError);
          }
          
          fileIn.seekg(currentPosition);
          fileIn.getline(dataLine,sizeof(dataLine));
          currentPosition = fileIn.tellg();
          
        } while (! string(dataLine).empty() );
        
      }  // end CB WPUP
      
      fileIn.seekg(currentPosition);
      fileIn.getline(dataLine,sizeof(dataLine));
      currentPosition = fileIn.tellg();
    
    } while ( fileIn.good() );
    
    /// Close input file
    fileIn.close();
    
  }  // end loop widths (to get fit params)
  
  if ( muVoigt.size() > 0 )  // if this vector not filled, none filled ==> all files not found
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
  cout << "width={";
  for (int i = 0; i < nWidths; i++)
  {
    cout << width[i];
    if ( i != nWidths-1) cout << ",";
  }
  cout << "}; (* times SM width *)" << endl;
  
  cout << endl << "-- VOIGT CP --" << endl;
  writeVoigtMathematicaOutput(muVoigt, sigmaVoigt, gammaVoigt, rVoigt, normVoigt, muErrVoigt, sigmaErrVoigt, gammaErrVoigt, rErrVoigt, normErrVoigt);
  
  cout << endl << "-- CRYSTALBALL WP --" << endl;
  writeCBMathematicaOutput(alphaCB_WP, nCB_WP, sigmaCB_WP, muCB_WP, normCB_WP, alphaErrCB_WP, nErrCB_WP, sigmaErrCB_WP, muErrCB_WP, normErrCB_WP);
  
  cout << endl << "-- CRYSTALBALL UP --" << endl;
  writeCBMathematicaOutput(alphaCB_UP, nCB_UP, sigmaCB_UP, muCB_UP, normCB_UP, alphaErrCB_UP, nErrCB_UP, sigmaErrCB_UP, muErrCB_UP, normErrCB_UP);
  
  cout << endl << "-- CRYSTALBALL WPUP --" << endl;
  writeCBMathematicaOutput(alphaCB_WPUP, nCB_WPUP, sigmaCB_WPUP, muCB_WPUP, normCB_WPUP, alphaErrCB_WPUP, nErrCB_WPUP, sigmaErrCB_WPUP, muErrCB_WPUP, normErrCB_WPUP);
  
  cout << endl << "=== LIKELIHOOD ===" << endl;
  cout << "const double mu_CP = " << muVoigt[0] << ", sigma_CP = " << sigmaVoigt[0] << ", r_CP = " << rVoigt[0] << ", norm_CP = " << normVoigt[0] << ";" << endl;
  cout << "const double alpha_WP = " << alphaCB_WP[0] << ", n_WP = " << nCB_WP[0] << ", sigma_WP = " << sigmaCB_WP[0] << ", mu_WP = " << muCB_WP[0] << ", norm_WP = " << normCB_WP[0] << ";" << endl;
  cout << "const double alpha_UP = " << alphaCB_UP[0] << ", n_UP = " << nCB_UP[0] << ", sigma_UP = " << sigmaCB_UP[0] << ", mu_UP = " << muCB_UP[0] << ", norm_UP = " << normCB_UP[0] << ";" << endl;
  cout << "const double alpha_WPUP = " << alphaCB_WPUP[0] << ", n_WPUP = " << nCB_WPUP[0] << ", sigma_WPUP = " << sigmaCB_WPUP[0] << ", mu_WPUP = " << muCB_WPUP[0] << ", norm_WPUP = " << normCB_WPUP[0] << ";" << endl;
}

void writeVoigtMathematicaOutput(std::vector<double> mu, std::vector<double> sigma, std::vector<double> gamma, std::vector<double> r, std::vector<double> norm, std::vector<double> muErr, std::vector<double> sigmaErr, std::vector<double> gammaErr, std::vector<double> rErr, std::vector<double> normErr)
{
  cout << "mu={";
  for (int i = 0; i < nWidths; i++) { cout << mu[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onzmu={";
  for (int i = 0; i < nWidths; i++) { cout << muErr[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "sigma={";
  for (int i = 0; i < nWidths; i++) { cout << sigma[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onzsigma={";
  for (int i = 0; i < nWidths; i++) { cout << sigmaErr[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "gamma={";
  for (int i = 0; i < nWidths; i++) { cout << gamma[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onzgamma={";
  for (int i = 0; i < nWidths; i++) { cout << gammaErr[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "r={";
  for (int i = 0; i < nWidths; i++) { cout << r[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onzr={";
  for (int i = 0; i < nWidths; i++) { cout << rErr[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "norm={";
  for (int i = 0; i < nWidths; i++) { cout << norm[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onznorm={";
  for (int i = 0; i < nWidths; i++) { cout << normErr[i]*1e+6; if ( i != nWidths-1) cout << ","; }
  cout << "}*10^(-6);" << endl;
}  


void writeCBMathematicaOutput(std::vector<double> alpha, std::vector<double> n, std::vector<double> sigma, std::vector<double> mu, std::vector<double> norm, std::vector<double> alphaErr, std::vector<double> nErr, std::vector<double> sigmaErr, std::vector<double> muErr, std::vector<double> normErr)
{
  cout << "alpha={";
  for (int i = 0; i < nWidths; i++) { cout << alpha[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onzalpha={";
  for (int i = 0; i < nWidths; i++) { cout << alphaErr[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "n={";
  for (int i = 0; i < nWidths; i++) { cout << n[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onzn={";
  for (int i = 0; i < nWidths; i++) { cout << nErr[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "sigma={";
  for (int i = 0; i < nWidths; i++) { cout << sigma[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onzsigma={";
  for (int i = 0; i < nWidths; i++) { cout << sigmaErr[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "mu={";
  for (int i = 0; i < nWidths; i++) { cout << mu[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onzmu={";
  for (int i = 0; i < nWidths; i++) { cout << muErr[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "norm={";
  for (int i = 0; i < nWidths; i++) { cout << norm[i]; if ( i != nWidths-1) cout << ","; }
  cout << "};" << endl << "onznorm={";
  for (int i = 0; i < nWidths; i++) { cout << normErr[i]*1e+6; if ( i != nWidths-1) cout << ","; }
  cout << "}*10^(-6);" << endl;
}  
