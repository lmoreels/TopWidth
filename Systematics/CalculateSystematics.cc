#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
//#include <pair>
#include <map>
#include <math.h>
#include <TMath.h>


using namespace std;



bool verbose_ = true;
string inputFileDir = "OutputLikelihood/170816_0947/";
string listFileName = "list_syst.txt";


string line, lineList;
string thisSyst;
double thisWidthValue, thisWidthUnc;

int nSys = -1;
vector<string> systName;
vector<double> systWidthValue;
vector<double> systWidthUnc;

vector<pair<string,double>> systDifference;

int indexNom;
int indexLepTrkUP, indexLepTrkDOWN, indexLepIdUP, indexLepIdDOWN, indexLepIsoUP, indexLepIsoDOWN, indexLepTrigUP, indexLepTrigDOWN, indexBTagUP, indexBTagDOWN, indexPuUP, indexPuDOWN, indexTopPt, indexLumiUP, indexLumiDOWN, indexRenFac1002, indexRenFac1003, indexRenFac1004, indexRenFac1005, indexRenFac1007, indexRenFac1009, indexIsrUP, indexIsrDOWN, indexFsrUP, indexFsrDOWN, indexHdampUP, indexHdampDOWN, indexTuneUP, indexTuneDOWN, indexCrErd, indexCrQcdErd, indexCrGluonMove, indexCrGluonMoveErd, indexMassUP, indexMassDOWN, indexHerwig;

void ClearIndices();
void ClearVars();
bool foundUP(string name);
bool foundDOWN(string name);
void IndexSystematics();
void WriteTable();

bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile.good();
}

void testfexists(const char *filename)
{
  if (! fexists(filename) )
  {
    std::cerr << "WARNING: File " << filename << " does not exist." << std::endl;
    exit(1);
  }
}

void ReadBashScript()
{
  int pid = fork();
  if ( pid == 0 )
  {
    execl("/bin/sh", "sh", (char*) "./ListFiles.sh", inputFileDir.c_str(), listFileName.c_str(), (char*) NULL);
  }
  int status;
  //waitpid(pid, &status, 0);
  while (-1 == waitpid(pid, &status, 0)) { ;}
  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
  {
    cerr << "Process (pid " << pid << ") failed" << endl;
    exit(1);
  }
}

int main()
{
  if (verbose_) cout << "Beginning of the program" << endl;
  ClearVars();
  
  /// List all files in directory and write to file
  ReadBashScript();
  
  /// Read file with list
  testfexists(listFileName.c_str());
  ifstream fileList(listFileName.c_str());
  while ( getline(fileList, lineList) )
  {
    /// Get input from files
    ifstream fileIn(lineList);
    if (verbose_) cout << "Opening " << lineList << "..." << endl;
    while ( getline(fileIn, line) )
    {
      if ( ! line.empty() )
      {
        istringstream iss(line);
        iss >> thisSyst >> thisWidthValue >> thisWidthUnc;
        
        systName.push_back(thisSyst);
        systWidthValue.push_back(thisWidthValue);
        systWidthUnc.push_back(thisWidthUnc);
      }
    }
    
    fileIn.close();
    
  }
  
  fileList.close();
  
  
  /// Process systematics
  nSys = systName.size();
  if ( nSys == -1 )
  {
    cerr << "Systematics not filled... Exiting..." << endl;
    exit(1);
  }
  if (verbose_) cout << "Number of measurements: " << nSys << endl;
  
//   for (int iSys = 0; iSys < nSys; iSys++)
//   {
//     if ( systName[iSys].find("nominal") != std::string::npos )
//     {
//       indexNom = iSys;
//       break;
//     }
//   }
//   if (verbose_) cout << "Index for nominal value is " << indexNom << endl;
  
  
  double nomVal = systWidthValue[indexNom];
  double diff;
  for (int iSys = 0; iSys < nSys; iSys++)
  {
    diff = systWidthValue[iSys] - nomVal;
    systDifference.push_back(pair<string,double>(systName[iSys],diff));
  }
  
  IndexSystematics();
  
  
  if (verbose_) cout << "End of the test program" << endl;
  
  return 0;
}

void ClearIndices()
{
  indexNom = -1;
  indexLepTrkUP = -1; indexLepTrkDOWN = -1;
  indexLepIdUP = -1; indexLepIdDOWN = -1;
  indexLepIsoUP = -1; indexLepIsoDOWN = -1;
  indexLepTrigUP = -1; indexLepTrigDOWN = -1;
  indexBTagUP = -1; indexBTagDOWN = -1;
  indexPuUP = -1; indexPuDOWN = -1;
  indexTopPt = -1;
  indexLumiUP = -1; indexLumiDOWN = -1;
  indexRenFac1002 = -1;
  indexRenFac1003 = -1;
  indexRenFac1004 = -1;
  indexRenFac1005 = -1;
  indexRenFac1007 = -1;
  indexRenFac1009 = -1;
  indexIsrUP = -1; indexIsrDOWN = -1;
  indexFsrUP = -1; indexFsrDOWN = -1;
  indexHdampUP = -1; indexHdampDOWN = -1;
  indexTuneUP = -1; indexTuneDOWN = -1;
  indexCrErd = -1;
  indexCrQcdErd = -1;
  indexCrGluonMove = -1;
  indexCrGluonMoveErd = -1;
  indexMassUP = -1; indexMassDOWN = -1;
  indexHerwig = -1;
}

void ClearVars()
{
  ClearIndices();
  
  systName.clear();
  systWidthValue.clear();
  systWidthUnc.clear();
  systDifference.clear();
}

bool foundUP(string name)
{
  if ( name.find("up") != std::string::npos || name.find("UP") != std::string::npos || name.find("Up") != std::string::npos ) return true;
  else return false;
}

bool foundDOWN(string name)
{
  if ( name.find("down") != std::string::npos || name.find("DOWN") != std::string::npos || name.find("Down") != std::string::npos ) return true;
  else return false;
}

void IndexSystematics()
{
  for (int iSys = 0; iSys < nSys; iSys++)
  {
    if ( systName[iSys].find("nominal") != std::string::npos ) indexNom = iSys;
    else if ( systName[iSys].find("lepton") != std::string::npos || systName[iSys].find("Lepton") != std::string::npos )
    {
      if ( systName[iSys].find("Id") != std::string::npos || systName[iSys].find("ID") != std::string::npos || systName[iSys].find("id") != std::string::npos )
      {
        if ( foundUP(systName[iSys]) ) indexLepIdUP = iSys;
        else if ( foundDOWN(systName[iSys]) ) indexLepIdDOWN = iSys;
      }
      else if ( systName[iSys].find("Iso") != std::string::npos || systName[iSys].find("ISO") != std::string::npos || systName[iSys].find("iso") != std::string::npos )
      {
        if ( foundUP(systName[iSys]) ) indexLepIsoUP = iSys;
        else if ( foundDOWN(systName[iSys]) ) indexLepIsoDOWN = iSys;
      }
      else if ( systName[iSys].find("Trig") != std::string::npos || systName[iSys].find("TRIG") != std::string::npos || systName[iSys].find("trig") != std::string::npos )
      {
        if ( foundUP(systName[iSys]) ) indexLepTrigUP = iSys;
        else if ( foundDOWN(systName[iSys]) ) indexLepTrigDOWN = iSys;
      }
      else if ( systName[iSys].find("Trk") != std::string::npos || systName[iSys].find("TRK") != std::string::npos || systName[iSys].find("trk") != std::string::npos )
      {
        if ( foundUP(systName[iSys]) ) indexLepTrkUP = iSys;
        else if ( foundDOWN(systName[iSys]) ) indexLepTrkDOWN = iSys;
      }
    }
    else if ( systName[iSys].find("bTag") != std::string::npos || systName[iSys].find("BTag") != std::string::npos || systName[iSys].find("btag") != std::string::npos )
    {
      if ( foundUP(systName[iSys]) ) indexBTagUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexBTagDOWN = iSys;
    }
    else if ( systName[iSys].find("Pu") != std::string::npos || systName[iSys].find("PU") != std::string::npos || systName[iSys].find("pu") != std::string::npos )
    {
      if ( foundUP(systName[iSys]) ) indexPuUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexPuDOWN = iSys;
    }
    else if ( systName[iSys].find("Lumi") != std::string::npos || systName[iSys].find("LUMI") != std::string::npos || systName[iSys].find("lumi") != std::string::npos )
    {
      if ( foundUP(systName[iSys]) ) indexLumiUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexLumiDOWN = iSys;
    }
    else if ( systName[iSys].find("topPt") != std::string::npos || systName[iSys].find("TopPt") != std::string::npos || systName[iSys].find("TOPPT") != std::string::npos || systName[iSys].find("toppt") != std::string::npos ) indexTopPt = iSys;
    else if ( systName[iSys].find("renFac") != std::string::npos || systName[iSys].find("RenFac") != std::string::npos || systName[iSys].find("RENFAC") != std::string::npos || systName[iSys].find("renfac") != std::string::npos )
    {
      if ( systName[iSys].find("1002") != std::string::npos ) indexRenFac1002 = iSys;
      else if ( systName[iSys].find("1003") != std::string::npos ) indexRenFac1003 = iSys;
      else if ( systName[iSys].find("1004") != std::string::npos ) indexRenFac1004 = iSys;
      else if ( systName[iSys].find("1005") != std::string::npos ) indexRenFac1005 = iSys;
      else if ( systName[iSys].find("1007") != std::string::npos ) indexRenFac1007 = iSys;
      else if ( systName[iSys].find("1009") != std::string::npos ) indexRenFac1009 = iSys;
    }
    else if ( systName[iSys].find("isr") != std::string::npos || systName[iSys].find("ISR") != std::string::npos || systName[iSys].find("Isr") != std::string::npos )
    {
      if ( foundUP(systName[iSys]) ) indexIsrUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexIsrDOWN = iSys;
    }
    else if ( systName[iSys].find("fsr") != std::string::npos || systName[iSys].find("FSR") != std::string::npos || systName[iSys].find("Fsr") != std::string::npos )
    {
      if ( foundUP(systName[iSys]) ) indexFsrUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexFsrDOWN = iSys;
    }
    else if ( systName[iSys].find("mass") != std::string::npos || systName[iSys].find("MASS") != std::string::npos || systName[iSys].find("Mass") != std::string::npos )
    {
      if ( systName[iSys].find("169") != std::string::npos ) indexMassUP = iSys;
      //if ( systName[iSys].find("171") != std::string::npos ) indexMassUP = iSys;
      //else if ( systName[iSys].find("173") != std::string::npos ) indexMassDOWN = iSys;
      else if ( systName[iSys].find("175") != std::string::npos ) indexMassDOWN = iSys;
    }
    //else if ( systName[iSys].find("") != std::string::npos )
    
  }
  
  if (verbose_)
  {
    cout << "Test indices: " << endl;
    for (int iSys = 0; iSys < nSys; iSys++)
    {
      cout << systDifference[iSys].first << " " << iSys << "    ";
    }
    cout << endl << endl;
    
    cout << "indexLepTrkUP  " << indexLepTrkUP << endl;
    cout << "indexLepTrkDOWN  " << indexLepTrkDOWN << endl;
    cout << "indexLepIdUP  " << indexLepIdUP << endl;
    cout << "indexLepIdDOWN  " << indexLepIdDOWN << endl;
    cout << "indexLepIsoUP  " << indexLepIsoUP << endl;
    cout << "indexLepIsoDOWN  " << indexLepIsoDOWN << endl;
    cout << "indexPuUP  " << indexPuUP << endl;
    cout << "indexTopPt  " << indexTopPt << endl;
    cout << "indexRenFac1003  " << indexRenFac1003 << endl;
  }
}

void WriteTable()
{
  cout << "Add code to write LaTeX table here" << endl;
}
