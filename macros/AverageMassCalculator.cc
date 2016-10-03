#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <cmath>
#include "TMath.h"

using namespace std;

/// Define inputs
string inputDate = "160928_1611";
string dataSetNames[] = {"TT", "ST_t_top", "ST_t_antitop", "ST_tW_top", "ST_tW_antitop", "DYJets", "WJets", "data"};
string pathInput = "averageMass/";
string inputFiles[] = {"mass_matched_"+dataSetNames[0]+"_"+inputDate, "mass_chi2_matched_"+dataSetNames[0]+"_"+inputDate, "mass_chi2_"+dataSetNames[0]+"_"+inputDate, "mass_chi2_"+dataSetNames[1]+"_"+inputDate, "mass_chi2_"+dataSetNames[2]+"_"+inputDate, "mass_chi2_"+dataSetNames[3]+"_"+inputDate, "mass_chi2_"+dataSetNames[4]+"_"+inputDate, "mass_chi2_"+dataSetNames[5]+"_"+inputDate, "mass_chi2_"+dataSetNames[6]+"_"+inputDate, "mass_chi2_"+dataSetNames[7]+"_"+inputDate};
int nInputs = sizeof(inputFiles)/sizeof(inputFiles[0]);

string inputFileName;
string thisDataSet;

/// Define functions
bool fexists(const char *filename);	// check if file exists
void ClearVars(bool isNewDataset);
void WriteToFile(std::ofstream &fout, std::string thisDataSet, double meanW, double meanTop);

/// Define vars
char dataLine[1024];
int nEntries, nEntriesChi2, eventId;
double massW, massTop;
double sumW, sumTop, sumWChi2, sumTopChi2;
double meanW, meanTop;

ifstream fileIn;
streampos currentPosition;
ofstream fileOut;

int main()
{
  string outputFileName = pathInput+"averageMass_"+inputDate+".txt";
  fileOut.open(outputFileName.c_str());
  cout << "Creating output file " << outputFileName << "..." << endl;
  fileOut << "# Dataset          meanW    meanTop" << endl;
  
  nEntriesChi2 = 0;
  sumWChi2 = 0.;
  sumTopChi2 = 0.;
  
  for (int iFile = 0; iFile < nInputs; iFile++)
  {
    ClearVars(true);
    
    inputFileName = pathInput+inputFiles[iFile]+".txt";
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
    while ( fileIn.good() )
    {
      ClearVars(false);
      
      fileIn.seekg(currentPosition);
      fileIn.getline(dataLine,sizeof(dataLine));
      istringstream iss(dataLine);
      iss >> eventId >> massW >> massTop;
      currentPosition = fileIn.tellg();
      
      nEntries++;
      sumW += massW;
      sumTop += massTop;
      
      if ( iFile > 1 ) // chi2
      {
        nEntriesChi2++;
        sumWChi2 += massW;
        sumTopChi2 += massTop;
      }
      
    }  // end while
    
    
    /// Calculate mean
    meanW = sumW/((double)nEntries);
    meanTop = sumTop/((double)nEntries);
    
    
    /// Store mean in file
    if ( iFile > 2 )
    {
      thisDataSet = dataSetNames[iFile-2];
    }
    else if ( iFile == 0 )
    {
      thisDataSet = dataSetNames[0]+"_match";
    }
    else if ( iFile == 1 )
    {
      thisDataSet = dataSetNames[0]+"_chi2_match";
    }
    else if ( iFile == 2 )
    {
      thisDataSet = dataSetNames[0]+"_chi2";
    }
    
    WriteToFile(fileOut, thisDataSet, meanW, meanTop);
    
    
    /// Close input file
    fileIn.close();
    
  }  // end loop files
  
  /// Calculate mean all chi2
  ClearVars(true);
  meanW = sumWChi2/((double)nEntriesChi2);
  meanTop = sumTopChi2/((double)nEntriesChi2);
  
  WriteToFile(fileOut, "All Chi2", meanW, meanTop);
  
  
  /// Close output file
  fileOut.close();
  
  cout << endl << " - Goodbye" << endl;
  
}


bool fexists(const char *filename)
{
	ifstream ifile(filename);
	return ifile;
}

void ClearVars(bool isNewDataset)
{
  eventId = -1;
  massW = 0;
  massTop = 0;
  
  if (isNewDataset)
  {
    thisDataSet = "";
    inputFileName = "";
    nEntries = 0;
    sumW = 0;
    sumTop = 0;
    meanW = 0;
    meanTop = 0;
  }
}

void WriteToFile(std::ofstream &fout, std::string thisDataSet, double meanW, double meanTop)
{
  fout << left << setw(16) << thisDataSet;
  cout.setf(ios::fixed,ios::floatfield);  // Add zero to obtain the asked number of digits
  //cout.precision(3);  // setprecision(3) --> set maximum number of meaningful digits to 3. When 'fixed' or 'scientific', only count after decimal point.
  fout << "   " << fixed << showpoint << setprecision(3) << meanW;  // showpoint --> also when decimal part is zero
  fout << "   " << fixed << showpoint << setprecision(3) << meanTop << endl;
}
