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
string inputDate = "170613_1638";  // Should be >= 170511_1314 !
string dataSetNames[] = {"TT", "ST_t_top", "ST_t_antitop", "ST_tW_top", "ST_tW_antitop", "DYJets", "WJets", "data"};
string pathInput = "averageMass/";
string inputFiles[] = {"mass_genp_matched_"+dataSetNames[0]+"_"+inputDate, "mass_genj_matched_"+dataSetNames[0]+"_"+inputDate, "mass_reco_matched_"+dataSetNames[0]+"_"+inputDate, "mass_reco_notCorrectMatch_"+dataSetNames[0]+"_"+inputDate, "mass_reco_notMatched_"+dataSetNames[0]+"_"+inputDate, "mass_reco_wrongPerm_"+dataSetNames[0]+"_"+inputDate, "mass_reco_"+dataSetNames[0]+"_nominal_"+inputDate, "mass_reco_"+dataSetNames[1]+"_"+inputDate, "mass_reco_"+dataSetNames[2]+"_"+inputDate, "mass_reco_"+dataSetNames[3]+"_"+inputDate, "mass_reco_"+dataSetNames[4]+"_"+inputDate, "mass_reco_"+dataSetNames[5]+"_"+inputDate, "mass_reco_"+dataSetNames[6]+"_"+inputDate, "mass_reco_"+dataSetNames[7]+"_"+inputDate};
int nInputs = sizeof(inputFiles)/sizeof(inputFiles[0]);

string inputFileName;
string thisDataSet;

/// Define functions
bool fexists(const char *filename);	// check if file exists
void ClearVars(bool isNewDataset);
void WriteToFile(std::ofstream &fout, std::string thisDataSet, double meanTop);

/// Define vars
char dataLine[1024];
int nEntries, nEntriesAllMC, nEntriesAllSamples, eventId;
double massTop;
double sumTop, sumTopAllMC, sumTopAllSamples;
double meanTop;

ifstream fileIn;
streampos currentPosition;
ofstream fileOut;

int main()
{
  string outputFileName = pathInput+"averageMass_"+inputDate+".txt";
  fileOut.open(outputFileName.c_str());
  cout << "Creating output file " << outputFileName << "..." << endl;
  fileOut << "# Dataset                     meanTop" << endl;
  
  nEntriesAllMC = 0;
  sumTopAllMC = 0.;
  nEntriesAllSamples = 0;
  sumTopAllSamples = 0.;
  
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
      iss >> eventId >> massTop;
      currentPosition = fileIn.tellg();
      
      nEntries++;
      sumTop += massTop;
      
      if ( iFile > 5 && iFile < nInputs-1) // reco
      {
        nEntriesAllMC++;
        sumTopAllMC += massTop;
        
        nEntriesAllSamples++;
        sumTopAllSamples += massTop;
      }
      else if ( iFile == nInputs-1 ) // data
      {
        nEntriesAllSamples++;
        sumTopAllSamples += massTop;
      }
      
    }  // end while
    
    
    /// Calculate mean
    meanTop = sumTop/((double)nEntries);
    
    
    /// Store mean in file
    if ( iFile > 6 )
    {
      thisDataSet = dataSetNames[iFile-6];
    }
    else if ( iFile == 0 )
    {
      thisDataSet = dataSetNames[0]+"_genp_match";
    }
    else if ( iFile == 1 )
    {
      thisDataSet = dataSetNames[0]+"_genj_match";
    }
    else if ( iFile == 2 )
    {
      thisDataSet = dataSetNames[0]+"_reco_match";
    }
    else if ( iFile == 3 )
    {
      thisDataSet = dataSetNames[0]+"_reco_wrongMatch_WP/UP";
    }
    else if ( iFile == 4 )
    {
      thisDataSet = dataSetNames[0]+"_reco_noMatch";
    }
    else if ( iFile == 5 )
    {
      thisDataSet = dataSetNames[0]+"_reco_wrongPerm";
    }
    else if ( iFile == 6 )
    {
      thisDataSet = dataSetNames[0]+"_reco";
    }
    
    WriteToFile(fileOut, thisDataSet, meanTop);
    
    
    /// Close input file
    fileIn.close();
    
  }  // end loop files
  
  /// Calculate mean reco all MC
  ClearVars(true);
  meanTop = sumTopAllMC/((double)nEntriesAllMC);
  
  WriteToFile(fileOut, "Reco All MC", meanTop);
  
  /// Calculate mean reco all samples
  ClearVars(true);
  meanTop = sumTopAllSamples/((double)nEntriesAllSamples);
  
  WriteToFile(fileOut, "Reco All Samples", meanTop);
  
  
  /// Close output file
  fileOut.close();
  
  cout << endl << " - Goodbye" << endl;
  
}


bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile.good();
}

void ClearVars(bool isNewDataset)
{
  eventId = -1;
  massTop = 0.;
  
  if (isNewDataset)
  {
    thisDataSet = "";
    inputFileName = "";
    nEntries = 0;
    sumTop = 0.;
    meanTop = 0.;
  }
}

void WriteToFile(std::ofstream &fout, std::string thisDataSet, double meanTop)
{
  fout << left << setw(27) << thisDataSet;
  cout.setf(ios::fixed,ios::floatfield);  // Add zero to obtain the asked number of digits
  //cout.precision(3);  // setprecision(3) --> set maximum number of meaningful digits to 3. When 'fixed' or 'scientific', only count after decimal point.
  fout << "   " << fixed << showpoint << setprecision(3) << meanTop << endl;  // showpoint --> also when decimal part is zero
}
