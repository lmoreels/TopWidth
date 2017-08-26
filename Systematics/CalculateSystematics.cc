#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <math.h>
#include <TMath.h>


using namespace std;



bool verbose_ = true;
string inputFileDir = "Systematics/temp/"; //"OutputLikelihood/170816_0947/";
string listFileName = "list_syst.txt";

const double calCurvePar_[2] = {0.00986675, 0.963861};
const double calCurveParUnc_[2] = {0.0759725, 0.0301143};


string line, lineList;
string thisSyst;
double thisWidth, thisWidthValue, thisWidthUnc;
pair<double,double> thisCorrWidthValue;
pair<double,double> testData;

int nSys = -1;
vector<string> systName;
vector<double> systWidthValue;
vector<double> systWidthUnc;

vector<double> systDifference, systRelDifference;
vector<double> systDifferenceUnc, systRelDifferenceUnc;

double nomVal, nomValUnc;

int indexNom;
int indexLepTrkUP, indexLepTrkDOWN, indexLepIdUP, indexLepIdDOWN, indexLepIsoUP, indexLepIsoDOWN, indexLepTrigUP, indexLepTrigDOWN, indexBTagUP, indexBTagDOWN, indexPuUP, indexPuDOWN, indexTopPt, indexLumiUP, indexLumiDOWN, indexRenFac1002, indexRenFac1003, indexRenFac1004, indexRenFac1005, indexRenFac1007, indexRenFac1009, indexIsrUP, indexIsrDOWN, indexFsrUP, indexFsrDOWN, indexHdampUP, indexHdampDOWN, indexTuneUP, indexTuneDOWN, indexCrErd, indexCrQcdErd, indexCrGluonMove, indexCrGluonMoveErd, indexMassUP, indexMassDOWN, indexHerwig;

void ClearIndices();
void ClearVars();
bool foundUP(string name);
bool foundDOWN(string name);
void IndexSystematics();
void TestIndexing();
std::pair<double,double> ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma);
void WriteShiftTable();
void WriteRelTable();

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
    /// Get input from files in list
    ifstream fileIn(lineList);
    if (verbose_) cout << "Opening " << lineList << "..." << endl;
    while ( getline(fileIn, line) )
    {
      if ( ! line.empty() )
      {
        istringstream iss(line);
        iss >> thisSyst >> thisWidth >> thisWidthValue >> thisWidthUnc;
        //cout << thisSyst << "  " << thisWidthValue << "  " << thisWidthUnc << endl;
        
        thisCorrWidthValue = ApplyCalibrationCurve(thisWidthValue, thisWidthUnc);
        
        systName.push_back(thisSyst);
        systWidthValue.push_back(thisCorrWidthValue.first);
        systWidthUnc.push_back(thisCorrWidthValue.second);
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
  
  IndexSystematics();
  
  nomVal = systWidthValue[indexNom];
  nomValUnc = systWidthUnc[indexNom];
  double diff, uncDiff, relDiff, uncRelDiff;
  for (int iSys = 0; iSys < nSys; iSys++)
  {
    diff = systWidthValue[iSys] - nomVal;
    systDifference.push_back(diff);
    uncDiff = TMath::Sqrt(systWidthUnc[iSys]*systWidthUnc[iSys] + nomValUnc*nomValUnc);
    systDifferenceUnc.push_back(uncDiff);
    relDiff = 100*diff/nomVal;
    systRelDifference.push_back(relDiff);
    uncRelDiff = TMath::Sqrt( systWidthUnc[iSys]*systWidthUnc[iSys]/(nomVal*nomVal) + nomValUnc*nomValUnc*systWidthValue[iSys]*systWidthValue[iSys]/pow(nomVal,4) );
    systRelDifferenceUnc.push_back(uncRelDiff);
  }
  
  
  WriteShiftTable();
  WriteRelTable();
  
  
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
  
  nomVal = 0.001;
  nomValUnc = 0.001;
  
  systName.clear();
  systWidthValue.clear();
  systWidthUnc.clear();
  systDifference.clear();
  systRelDifference.clear();
  systDifferenceUnc.clear();
  systRelDifferenceUnc.clear();
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
    else if ( systName[iSys].find("hdamp") != std::string::npos )
    {
      if ( foundUP(systName[iSys]) ) indexHdampUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexHdampDOWN = iSys;
    }
    else if ( systName[iSys].find("tune") != std::string::npos )
    {
      if ( foundUP(systName[iSys]) ) indexTuneUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexTuneDOWN = iSys;
    }
    else if ( systName[iSys].find("ERD") != std::string::npos )
    {
      if ( systName[iSys].find("mpi") != std::string::npos ) indexCrErd = iSys;
      else if ( systName[iSys].find("qcd") != std::string::npos ) indexCrQcdErd = iSys;
      else if ( systName[iSys].find("gluon") != std::string::npos ) indexCrGluonMoveErd = iSys;
    }
    else if ( systName[iSys].find("gluon") != std::string::npos ) indexCrGluonMove = iSys;
    //else if ( systName[iSys].find("") != std::string::npos )
    
    else
      cerr << "Systematic " << systName[iSys] << " not found! Please update list..." << endl;
  }
  
  //if (verbose_) TestIndexing();
  
}

void TestIndexing()
{
  cout << endl << "Test indices: " << endl;
  for (int iSys = 0; iSys < nSys; iSys++)
  {
    cout << systName[iSys] << " " << iSys << "    ";
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

std::pair<double,double> ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma)
{
  /// return 'thisInputWidth' and uncertainty
  //  thisOutputWidth = Par_[0] + Par_[1] * thisInputWidth
  
  double thisInputWidth = (thisOutputWidth - calCurvePar_[0])/calCurvePar_[1];
  double thisInputWidthSigma = TMath::Sqrt( thisOutputWidthSigma*thisOutputWidthSigma + calCurveParUnc_[0]*calCurveParUnc_[0] + calCurveParUnc_[1]*calCurveParUnc_[1]*thisInputWidth*thisInputWidth )/calCurvePar_[1];
  
  return std::pair<double,double>(thisInputWidth,thisInputWidthSigma);
  
}

void WriteShiftTable()
{
  string space = "     ";
  string hline = space+"\\hline";
  string interline = space+"&&&&\\\\[-6pt]";
  string headerextra = "[+3pt]";
  string multicolBegin = "\\multicolumn{2}{r}{";
  string multicolEnd = "$\\qquad$}";
  
  
  ofstream fileOut("systematicsShiftFull.tex");
  fileOut << "\\begin{table}[htp]" << endl;
  fileOut << " \\caption{\\fixme{Put caption here}}" << endl;
  fileOut << " \\begin{center}" << endl;
  fileOut << "  \\begin{minipage}[p]{\\textwidth}" << endl;
  fileOut << "   \\begin{center}" << endl;
  fileOut << "    %\\hspace{-2.6cm}" << endl;
  fileOut << "    {\\small" << endl;
  fileOut << "    \\begin{tabular}{l S[table-format=1.4] S[table-format=1.4] S[table-format=1.4] S[table-format=1.4]}" << endl;  /// adapt!!
  
  fileOut << hline << endl;
  fileOut << hline << endl;
  fileOut << interline << endl;
  fileOut << space << "\\multicolumn{3}{l}{Nominal value} & \\multicolumn{2}{l}{$" << fixed << setprecision(4) << nomVal << " \\pm " << nomValUnc << "$} \\\\" << headerextra << endl;
  fileOut << hline << endl;
  fileOut << interline << endl;
  fileOut << space << "Systematic & \\multicolumn{2}{c}{$+1\\sigma$} & \\multicolumn{2}{c}{$-1\\sigma$} \\\\" << endl;
  fileOut << space << " & {Shift} & {   Unc   } & {Shift} & {   Unc   } \\\\" << headerextra << endl;
  fileOut << hline << endl;
  
  fileOut << interline << endl;
  fileOut << space << "\\textit{Lepton SFs} &  &  &  &  \\\\" << endl;
  if ( indexLepIdUP != -1 && indexLepIdDOWN != -1 )
    fileOut << space << "\\tabsp Id & " << fixed << setprecision(4) << systDifference[indexLepIdUP] << " & " << systDifferenceUnc[indexLepIdUP] << " & " << systDifference[indexLepIdDOWN] << " & " << systDifferenceUnc[indexLepIdDOWN] << " \\\\" << endl;
  if ( indexLepIsoUP != -1 && indexLepIsoDOWN != -1 )
    fileOut << space << "\\tabsp Isolation & " << fixed << setprecision(4) << systDifference[indexLepIsoUP] << " & " << systDifferenceUnc[indexLepIsoUP] << " & " << systDifference[indexLepIsoDOWN] << " & " << systDifferenceUnc[indexLepIsoDOWN] << " \\\\" << endl;
  if ( indexLepTrigUP != -1 && indexLepTrigDOWN != -1 )
    fileOut << space << "\\tabsp Trigger & " << fixed << setprecision(4) << systDifference[indexLepTrigUP] << " & " << systDifferenceUnc[indexLepTrigUP] << " & " << systDifference[indexLepTrigDOWN] << " & " << systDifferenceUnc[indexLepTrigDOWN] << " \\\\" << endl;
  if ( indexLepTrkUP != -1 && indexLepTrkDOWN != -1 )
    fileOut << space << "\\tabsp Tracking & " << fixed << setprecision(4) << systDifference[indexLepTrkUP] << " & " << systDifferenceUnc[indexLepTrkUP] << " & " << systDifference[indexLepTrkDOWN] << " & " << systDifferenceUnc[indexLepTrkDOWN] << " \\\\" << endl;
  
  if ( indexBTagUP != -1 && indexBTagDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\bq~tagging SFs & " << fixed << setprecision(4) << systDifferenceUnc[indexBTagUP] << " & " << systDifferenceUnc[indexBTagUP] << " & " << systDifference[indexBTagDOWN] << " & " << systDifferenceUnc[indexBTagDOWN] << " \\\\" << endl;
  }
  
  if ( indexPuUP != -1 && indexPuDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << fixed << setprecision(4) << "\\pileup\\ SFs & " << systDifference[indexPuUP] << " & " << systDifferenceUnc[indexPuUP] << " & " << systDifference[indexPuDOWN] << " & " << systDifferenceUnc[indexPuDOWN] << " \\\\" << endl;
  }
  
  if ( indexLumiUP != -1 && indexLumiDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Luminosity & " << fixed << setprecision(4) << systDifference[indexLumiUP] << " & " << systDifferenceUnc[indexLumiUP] << " & " << systDifference[indexLumiDOWN] << " & " << systDifferenceUnc[indexLumiDOWN] << " \\\\" << endl;
  }
  
  if ( indexMassUP != -1 && indexMassDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Top mass & +" << fixed << setprecision(4) << systDifference[indexMassUP] << " & " << systDifferenceUnc[indexMassUP] << " & " << systDifference[indexMassDOWN] << " & " << systDifferenceUnc[indexMassDOWN] << " \\\\" << endl;
  }
  
  if ( indexTopPt != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Top $\\pT$ reweighting & \\multicolumn{3}{r}{" << fixed << setprecision(4) << systDifference[indexTopPt] << "} & " << systDifferenceUnc[indexTopPt] << " \\\\" << endl;
  }
  
  if ( indexRenFac1002 != -1 || indexRenFac1003 != -1 || indexRenFac1004 != -1 || indexRenFac1005 != -1 || indexRenFac1007 != -1 || indexRenFac1009 != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{Renormalisation/factorisation scale} &  &  &  &  \\\\" << endl;
    if ( indexRenFac1002 != -1 )
      fileOut << space << "\\tabsp 1002 & \\multicolumn{3}{r}{" << systDifference[indexRenFac1002] << "} & " << systDifferenceUnc[indexRenFac1002] << " \\\\" << endl;
    if ( indexRenFac1003 != -1 )
      fileOut << space << "\\tabsp 1003 & \\multicolumn{3}{r}{" << systDifference[indexRenFac1003] << "} & " << systDifferenceUnc[indexRenFac1003] << " \\\\" << endl;
    if ( indexRenFac1004 != -1 )
      fileOut << space << "\\tabsp 1004 & \\multicolumn{3}{r}{" << systDifference[indexRenFac1004] << "} & " << systDifferenceUnc[indexRenFac1004] << " \\\\" << endl;
    if ( indexRenFac1005 != -1 )
      fileOut << space << "\\tabsp 1005 & \\multicolumn{3}{r}{" << systDifference[indexRenFac1005] << "} & " << systDifferenceUnc[indexRenFac1005] << " \\\\" << endl;
    if ( indexRenFac1007 != -1 )
      fileOut << space << "\\tabsp 1007 & \\multicolumn{3}{r}{" << systDifference[indexRenFac1007] << "} & " << systDifferenceUnc[indexRenFac1007] << " \\\\" << endl;
    if ( indexRenFac1009 != -1 )
      fileOut << space << "\\tabsp 1009 & \\multicolumn{3}{r}{" << systDifference[indexRenFac1009] << "} & " << systDifferenceUnc[indexRenFac1009] << " \\\\" << endl;
  }
  
  if ( indexIsrUP != -1 && indexIsrDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "ISR & " << systDifference[indexIsrUP] << " & " << systDifferenceUnc[indexIsrUP] << " & " << systDifference[indexIsrDOWN] << " & " << systDifferenceUnc[indexIsrDOWN] << " \\\\" << endl;
  }
  if ( indexFsrUP != -1 && indexFsrDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "FSR & " << systDifference[indexFsrUP] << " & " << systDifferenceUnc[indexFsrUP] << " & " << systDifference[indexFsrDOWN] << " & " << systDifferenceUnc[indexFsrDOWN] << " \\\\" << endl;
  }
  if ( indexHdampUP != -1 && indexHdampDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\hdamp & " << systDifference[indexHdampUP] << " & " << systDifferenceUnc[indexHdampUP] << " & " << systDifference[indexHdampDOWN] << " & " << systDifferenceUnc[indexHdampDOWN] << " \\\\" << endl;
  }
  if ( indexTuneUP != -1 && indexTuneDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Tune & " << systDifference[indexTuneUP] << " & " << systDifferenceUnc[indexTuneUP] << " & " << systDifference[indexTuneDOWN] << " & " << systDifferenceUnc[indexTuneDOWN] << " \\\\" << endl;
  }
  if ( indexCrErd != -1 || indexCrQcdErd != -1 || indexCrGluonMove != -1 || indexCrGluonMoveErd != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{Colour reconnection} &  &  &  &  \\\\" << endl;
    if ( indexCrErd != -1 )
      fileOut << space << "\\tabsp MPI (ERD) & \\multicolumn{3}{r}{" << systDifference[indexCrErd] << "} & " << systDifferenceUnc[indexCrErd] << " \\\\" << endl;
    if ( indexCrQcdErd != -1 )
      fileOut << space << "\\tabsp QCD-based (ERD) & \\multicolumn{3}{r}{" << systDifference[indexCrQcdErd] << "} & " << systDifferenceUnc[indexCrQcdErd] << " \\\\" << endl;
    if ( indexCrGluonMove != -1 )
      fileOut << space << "\\tabsp Gluon move & \\multicolumn{3}{r}{" << systDifference[indexCrGluonMove] << "} & " << systDifferenceUnc[indexCrGluonMove] << " \\\\" << endl;
    if ( indexCrGluonMoveErd != -1 )
      fileOut << space << "\\tabsp Gluon move (ERD) & \\multicolumn{3}{r}{" << systDifference[indexCrGluonMoveErd] << "} & " << systDifferenceUnc[indexCrGluonMoveErd] << " \\\\" << endl;
  }
  
  if ( indexHerwig != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Herwig & \\multicolumn{3}{r}{" << systDifference[indexHerwig] << "} & " << systDifferenceUnc[indexHerwig] << " \\\\" << endl;
  }
  
  //....
  
  
  fileOut << interline << endl;
  fileOut << hline << endl;
  fileOut << hline << endl;
  
  fileOut << "    \\end{tabular}}" << endl;
  fileOut << "   \\end{center}" << endl;
  fileOut << "  \\end{minipage}" << endl;
  fileOut << " \\end{center}" << endl;
  fileOut << " \\label{tab:systFullShift}" << endl;
  fileOut << "\\end{table}%" << endl;
  
  cout << "Written table with absolute shifts + stat. uncertainties" << endl;
  
}

void WriteRelTable()
{
  string space = "     ";
  string hline = space+"\\hline";
  string interline = space+"&&\\\\[-6pt]";
  string headerextra = "[+3pt]";
  
  
  ofstream fileOut("systematicsRelFull.tex");
  fileOut << "\\begin{table}[htp]" << endl;
  fileOut << " \\caption{\\fixme{Put caption here}}" << endl;
  fileOut << " \\begin{center}" << endl;
  fileOut << "  \\begin{minipage}[p]{\\textwidth}" << endl;
  fileOut << "   \\begin{center}" << endl;
  fileOut << "    %\\hspace{-2.6cm}" << endl;
  fileOut << "    {\\small" << endl;
  fileOut << "    \\begin{tabular}{l S[table-format=1.4] S[table-format=1.4]}" << endl;  /// adapt!!
  
  fileOut << hline << endl;
  fileOut << hline << endl;
  fileOut << interline << endl;
  fileOut << space << "Systematic & {$+1\\sigma$ ($\\%$)} & {$-1\\sigma$ ($\\%$)} \\\\" << headerextra << endl;
  fileOut << hline << endl;
  
  fileOut << interline << endl;
  fileOut << space << "\\textit{Lepton SFs} &  & \\\\" << endl;
  if ( indexLepIdUP != -1 && indexLepIdDOWN != -1 )
    fileOut << space << "\\tabsp Id & +" << fixed << setprecision(4) << systRelDifference[indexLepIdUP] << " & " << fixed << setprecision(4) << systRelDifference[indexLepIdDOWN] << " \\\\" << endl;
  if ( indexLepIsoUP != -1 && indexLepIsoDOWN != -1 )
    fileOut << space << "\\tabsp Isolation & +" << fixed << setprecision(4) << systRelDifference[indexLepIsoUP] << " & " << fixed << setprecision(4) << systRelDifference[indexLepIsoDOWN] << " \\\\" << endl;
  if ( indexLepTrigUP != -1 && indexLepTrigDOWN != -1 )
    fileOut << space << "\\tabsp Trigger & +" << fixed << setprecision(4) << systRelDifference[indexLepTrigUP] << " & " << fixed << setprecision(4) << systRelDifference[indexLepTrigDOWN] << " \\\\" << endl;
  if ( indexLepTrkUP != -1 && indexLepTrkDOWN != -1 )
    fileOut << space << "\\tabsp Tracking & " << fixed << setprecision(4) << systRelDifference[indexLepTrkUP] << " & " << fixed << setprecision(4) << systRelDifference[indexLepTrkDOWN] << " \\\\" << endl;
  
  if ( indexBTagUP != -1 && indexBTagDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\bq~tagging SFs & +" << fixed << setprecision(4) << systRelDifference[indexBTagUP] << " & " << fixed << setprecision(4) << systRelDifference[indexBTagDOWN] << " \\\\" << endl;
  }
  
  if ( indexPuUP != -1 && indexPuDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << fixed << setprecision(4) << "\\pileup\\ SFs & +" << systRelDifference[indexPuUP] << " & " << fixed << setprecision(4) << systRelDifference[indexPuDOWN] << " \\\\" << endl;
  }
  
  if ( indexLumiUP != -1 && indexLumiDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Luminosity & +" << fixed << setprecision(4) << systRelDifference[indexLumiUP] << " & " << fixed << setprecision(4) << systRelDifference[indexLumiDOWN] << " \\\\" << endl;
  }
  
  if ( indexMassUP != -1 && indexMassDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Top mass & +" << systRelDifference[indexMassUP] << " & " << systRelDifference[indexMassDOWN] << " \\\\" << endl;
  }
  
  if ( indexTopPt != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Top $\\pT$ reweighting & \\multicolumn{2}{r}{" << fixed << setprecision(4) << systRelDifference[indexTopPt] << "$\\qquad\\quad$} \\\\" << endl;
  }
  
  if ( indexRenFac1002 != -1 || indexRenFac1003 != -1 || indexRenFac1004 != -1 || indexRenFac1005 != -1 || indexRenFac1007 != -1 || indexRenFac1009 != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{Renormalisation/factorisation scale} &  & \\\\" << endl;
    if ( indexRenFac1002 != -1 )
      fileOut << space << "\\tabsp 1002 & \\multicolumn{2}{r}{" << systRelDifference[indexRenFac1002] << "$\\qquad\\quad$} \\\\" << endl;
    if ( indexRenFac1003 != -1 )
      fileOut << space << "\\tabsp 1003 & \\multicolumn{2}{r}{" << systRelDifference[indexRenFac1003] << "$\\qquad\\quad$} \\\\" << endl;
    if ( indexRenFac1004 != -1 )
      fileOut << space << "\\tabsp 1004 & \\multicolumn{2}{r}{" << systRelDifference[indexRenFac1004] << "$\\qquad\\quad$} \\\\" << endl;
    if ( indexRenFac1005 != -1 )
      fileOut << space << "\\tabsp 1005 & \\multicolumn{2}{r}{" << systRelDifference[indexRenFac1005] << "$\\qquad\\quad$} \\\\" << endl;
    if ( indexRenFac1007 != -1 )
      fileOut << space << "\\tabsp 1007 & \\multicolumn{2}{r}{" << systRelDifference[indexRenFac1007] << "$\\qquad\\quad$} \\\\" << endl;
    if ( indexRenFac1009 != -1 )
      fileOut << space << "\\tabsp 1009 & \\multicolumn{2}{r}{" << systRelDifference[indexRenFac1009] << "$\\qquad\\quad$} \\\\" << endl;
  }
  
  if ( indexIsrUP != -1 && indexIsrDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "ISR & +" << systRelDifference[indexIsrUP] << " & " << systRelDifference[indexIsrDOWN] << " \\\\" << endl;
  }
  if ( indexFsrUP != -1 && indexFsrDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "FSR & +" << systRelDifference[indexFsrUP] << " & " << systRelDifference[indexFsrDOWN] << " \\\\" << endl;
  }
  if ( indexHdampUP != -1 && indexHdampDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\hdamp & +" << systRelDifference[indexHdampUP] << " & " << systRelDifference[indexHdampDOWN] << " \\\\" << endl;
  }
  if ( indexTuneUP != -1 && indexTuneDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Tune & +" << systRelDifference[indexTuneUP] << " & " << systRelDifference[indexTuneDOWN] << " \\\\" << endl;
  }
  if ( indexCrErd != -1 || indexCrQcdErd != -1 || indexCrGluonMove != -1 || indexCrGluonMoveErd != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{Colour reconnection} &  & \\\\" << endl;
    if ( indexCrErd != -1 )
      fileOut << space << "\\tabsp MPI (ERD) & \\multicolumn{2}{r}{" << systRelDifference[indexCrErd] << "$\\qquad\\quad$} \\\\" << endl;
    if ( indexCrQcdErd != -1 )
      fileOut << space << "\\tabsp QCD-based (ERD) & \\multicolumn{2}{r}{" << systRelDifference[indexCrQcdErd] << "$\\qquad\\quad$} \\\\" << endl;
    if ( indexCrGluonMove != -1 )
      fileOut << space << "\\tabsp Gluon move & \\multicolumn{2}{r}{" << systRelDifference[indexCrGluonMove] << "$\\qquad\\quad$} \\\\" << endl;
    if ( indexCrGluonMoveErd != -1 )
      fileOut << space << "\\tabsp Gluon move (ERD) & \\multicolumn{2}{r}{" << systRelDifference[indexCrGluonMoveErd] << "$\\qquad\\quad$} \\\\" << endl;
  }
  
  if ( indexHerwig != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Herwig & \\multicolumn{2}{c}{" << systRelDifference[indexHerwig] << "} \\\\" << endl;
  }
  
  
  //....
  
  fileOut << interline << endl;
  fileOut << hline << endl;
  fileOut << hline << endl;
  
  fileOut << "    \\end{tabular}}" << endl;
  fileOut << "   \\end{center}" << endl;
  fileOut << "  \\end{minipage}" << endl;
  fileOut << " \\end{center}" << endl;
  fileOut << " \\label{tab:systFullRel}" << endl;
  fileOut << "\\end{table}%" << endl;
  
  cout << "Written table with relative uncertainties" << endl;
}
