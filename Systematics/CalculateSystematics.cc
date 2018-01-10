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
#include <algorithm>


using namespace std;



bool verbose_ = true;
bool doMass = false;
bool is2D = false;
string suffix = "";
string inputFileDir = "Systematics/temp"+suffix+"/"; //"OutputLikelihood/170816_0947/";
string listFileName = "list_syst.txt";


/// 1D width
// parameter = comb
//const double calCurvePar_[2] = {-0.249918, 1.08351};   //{-0.2526, 1.084};
//const double calCurveParUnc_[2] = {0.0351629, 0.0139145};   //{0.03321, 0.01345};
// parameter = redMlb
//const double calCurvePar_[2] = {-0.267184, 1.08957};   //{-0.245, 1.083};
//const double calCurveParUnc_[2] = {0.0550971, 0.0248913};   //{0.05185, 0.02447};
// parameter = redTopMass_old
const double calCurvePar_[2] = {-0.152752, 1.09383};   //{-0.15551, 1.09669};        // 180107
const double calCurveParUnc_[2] = {0.0436359, 0.0166395};   //{0.0469114, 0.0175116};  // 180107
/// 1D mass
const double calCurveParMass_[2] = {12.0431, 0.930349};      // 180108
const double calCurveParUncMass_[2] = {2.81821, 0.0163266};  // 180108
/// 2D
const double calCurveParW_[2] = {-0.197866, 1.07643};
const double calCurveParUncW_[2] = {0.0451828, 0.0180112};
const double calCurveParM_[2] = {0., 1.};
const double calCurveParUncM_[2] = {0., 0.};


string line, lineList;
string thisSyst;
double thisWidth, thisMass, thisWidthValue, thisWidthUnc, thisMassValue, thisMassUnc;
pair<double,double> thisCorrWidthValue, thisCorrMassValue;
pair<double,double> testData;

bool isRate[80] = {0};

int nSys = -1;
vector<string> systName;
vector<double> genWidths;
vector<double> genMasses;
vector<double> systWidthValue;
vector<double> systWidthUnc;
vector<double> systMassValue;
vector<double> systMassUnc;

vector<double> systDifference, systRelDifference;
vector<double> systDifferenceUnc, systRelDifferenceUnc;
vector<double> systMassDifference, systRelMassDifference;
vector<double> systMassDifferenceUnc, systRelMassDifferenceUnc;

double nomVal, nomValUnc;
double totalShiftUp, totalShiftDown, maxScale, minScale, maxCR, minCR;
double shiftUp, shiftDown, tmp;


int indexNom;
int indexLepTrkUP, indexLepTrkDOWN, indexLepIdUP, indexLepIdDOWN, indexLepIsoUP, indexLepIsoDOWN, indexLepTrigUP, indexLepTrigDOWN, indexBTagUP, indexBTagDOWN, indexPuUP, indexPuDOWN, indexTopPt, indexLumiUP, indexLumiDOWN, indexRenFac1002, indexRenFac1003, indexRenFac1004, indexRenFac1005, indexRenFac1007, indexRenFac1009, indexFragCentral, indexFragUP, indexFragDOWN, indexFragPeterson, indexFragSemiLepBrUP, indexFragSemiLepBrDOWN, indexIsrUP, indexIsrDOWN, indexFsrUP, indexFsrDOWN, indexHdampUP, indexHdampDOWN, indexTuneUP, indexTuneDOWN, indexCrErd, indexCrQcdErd, indexCrGluonMove, indexCrGluonMoveErd, indexMassUP, indexMassDOWN, indexHerwig, indexJesUP, indexJesDOWN, indexJerUP, indexJerDOWN, indexPdfAlphaSUP, indexPdfAlphaSDOWN, indexPdfVar, indexRateCMUP, indexRateCMDOWN, indexRateSTtUP, indexRateSTtDOWN, indexRateSTtWUP, indexRateSTtWDOWN, indexRateOtherUP, indexRateOtherDOWN;

void ClearIndices();
void ClearVars();
bool foundUP(string name);
bool foundDOWN(string name);
double GetMinimum(int i1, int i2, int i3 = -1, int i4 = -1, int i5 = -1, int i6 = -1, int i7 = -1, int i8 = -1, int i9 = -1, int i10 = -1);
double GetMaximum(int i1, int i2, int i3 = -1, int i4 = -1, int i5 = -1, int i6 = -1, int i7 = -1, int i8 = -1, int i9 = -1, int i10 = -1);
double isPositive(double upVar);
double isNegative(double downVar);
void IndexSystematics();
void TestIndexing();
std::pair<double,double> ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma);
std::pair<double,double> ApplyCalibrationCurveMass(double thisOutputMass, double thisOutputMassSigma);
std::pair<double,double> ApplyCalibrationCurveW(double thisOutputWidth, double thisOutputWidthSigma);
std::pair<double,double> ApplyCalibrationCurveM(double thisOutputMass, double thisOutputMassSigma);
void WriteShiftTable(double scale = 1.);
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

string ConvertDoubleToString(double Number)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

string DotReplace(double var)
{
  string str = ConvertDoubleToString(var);
  replace(str.begin(), str.end(), '.', 'p');
  return str;
}

int main()
{
  if (verbose_) cout << "Beginning of the program" << endl;
  ClearVars();
  
  if (is2D) suffix = "2D";
  else if (doMass) suffix = "Mass";
  else suffix = "";
  inputFileDir = "Systematics/temp"+suffix+"/";
  
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
        iss >> thisSyst >> thisWidth >> thisMass >> thisWidthValue >> thisWidthUnc;
        if (is2D) iss >> thisMassValue >> thisMassUnc;
        //cout << thisSyst << "  " << thisWidthValue << "  " << thisWidthUnc << endl;
        
        if (is2D)
        {
          thisCorrWidthValue = ApplyCalibrationCurveW(thisWidthValue, thisWidthUnc);
          thisCorrMassValue = ApplyCalibrationCurveM(thisMassValue, thisMassUnc);
        }
        else if (doMass) thisCorrWidthValue = ApplyCalibrationCurveMass(thisWidthValue, thisWidthUnc);
        else thisCorrWidthValue = ApplyCalibrationCurve(thisWidthValue, thisWidthUnc);
        
        systName.push_back(thisSyst);
        genWidths.push_back(thisWidth);
        genMasses.push_back(thisMass);
        systWidthValue.push_back(thisCorrWidthValue.first);
        systWidthUnc.push_back(thisCorrWidthValue.second);
        if (is2D)
        {
          systMassValue.push_back(thisCorrMassValue.first);
          systMassUnc.push_back(thisCorrMassValue.second);
        }
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
    if      ( iSys == indexFsrUP )    diff *= 1./sqrt(2.);
    else if ( iSys == indexFsrDOWN )  diff *= 1./sqrt(2.);
//    else if ( iSys == indexMassUP )   diff *= 1./3.;
//    else if ( iSys == indexMassDOWN ) diff *= 1./3.;
    
    systDifference.push_back(diff);
    uncDiff = TMath::Sqrt(systWidthUnc[iSys]*systWidthUnc[iSys] + nomValUnc*nomValUnc);
    systDifferenceUnc.push_back(uncDiff);
    relDiff = 100*diff/nomVal;
    systRelDifference.push_back(relDiff);
    uncRelDiff = TMath::Sqrt( systWidthUnc[iSys]*systWidthUnc[iSys]/(nomVal*nomVal) + nomValUnc*nomValUnc*systWidthValue[iSys]*systWidthValue[iSys]/pow(nomVal,4) );
    systRelDifferenceUnc.push_back(uncRelDiff);
  }
  
  vector<double> scaleQ2 = {systDifference[indexRenFac1002], systDifference[indexRenFac1003], systDifference[indexRenFac1004], systDifference[indexRenFac1005], systDifference[indexRenFac1007], systDifference[indexRenFac1009], systDifference[indexIsrUP], systDifference[indexIsrDOWN], systDifference[indexFsrUP], systDifference[indexFsrDOWN]};
  maxScale = *max_element(scaleQ2.begin(), scaleQ2.end());
  minScale = *min_element(scaleQ2.begin(), scaleQ2.end());
  cout << "Maximum shift for Q^2 scale is " << maxScale << " and minimum " << minScale << endl;
  
  vector<double> CR = {systDifference[indexCrErd], systDifference[indexCrQcdErd], systDifference[indexCrGluonMove], systDifference[indexCrGluonMoveErd]};
  maxCR = *max_element(CR.begin(), CR.end());
  minCR = *min_element(CR.begin(), CR.end());
  cout << "Maximum shift for colour reconnection is " << maxCR << " and minimum " << minCR << endl;
  
  
  // Only add positive shift to positive unc & neg. shifts to neg. unc.
  shiftUp = 0.;
  if ( indexLepIdUP != -1 && indexLepIdDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexLepIdUP, indexLepIdDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexLepIsoUP != -1 && indexLepIsoDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexLepIsoUP, indexLepIsoDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexLepTrigUP != -1 && indexLepTrigDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexLepTrigUP, indexLepTrigDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexLepTrkUP != -1 && indexLepTrkDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexLepTrkUP, indexLepTrkDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexBTagUP != -1 && indexBTagDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexBTagUP, indexBTagDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexPuUP != -1 && indexPuDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexPuUP, indexPuDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexLumiUP != -1 && indexLumiDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexLumiUP, indexLumiDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexMassUP != -1 && indexMassDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexMassUP, indexMassDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexTopPt != -1 )
  {
    tmp = isPositive(systDifference[indexTopPt]);
    shiftUp += tmp*tmp;
  }
  if ( indexRenFac1002 != -1 || indexRenFac1003 != -1 || indexRenFac1004 != -1 || indexRenFac1005 != -1 || indexRenFac1007 != -1 || indexRenFac1009 != -1 || (indexIsrUP != -1 && indexIsrDOWN != -1) || (indexFsrUP != -1 && indexFsrDOWN != -1) )
  {
    tmp = isPositive(GetMaximum(indexRenFac1002, indexRenFac1003, indexRenFac1004, indexRenFac1005, indexRenFac1007, indexRenFac1009, indexIsrUP, indexIsrDOWN, indexFsrUP, indexFsrDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexHdampUP != -1 && indexHdampDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexHdampUP, indexHdampDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexFragCentral != -1 || indexFragUP != -1 || indexFragDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexFragCentral, indexFragUP, indexFragDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexFragPeterson != -1 )
  {
    tmp = isPositive(systDifference[indexFragPeterson]);
    shiftUp += tmp*tmp;
  }
  if ( indexFragSemiLepBrUP != -1 && indexFragSemiLepBrDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexFragSemiLepBrUP, indexFragSemiLepBrDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexPdfAlphaSUP != -1 && indexPdfAlphaSDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexPdfAlphaSUP, indexPdfAlphaSDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexPdfVar != -1 )
  {
    shiftUp += systDifference[indexPdfVar]*systDifference[indexPdfVar];
  }
  if ( indexTuneUP != -1 && indexTuneDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexTuneUP, indexTuneDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexCrErd != -1 || indexCrQcdErd != -1 || indexCrGluonMove != -1 || indexCrGluonMoveErd != -1 )
  {
    tmp = isPositive(GetMaximum(indexCrErd, indexCrQcdErd, indexCrGluonMove, indexCrGluonMoveErd));
    shiftUp += tmp*tmp;
  }
  if ( indexJesUP != -1 && indexJesDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexJesUP, indexJesDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexJerUP != -1 && indexJerDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexJerUP, indexJerDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexRateCMUP != -1 && indexRateCMDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexRateCMUP, indexRateCMDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexRateSTtUP != -1 && indexRateSTtDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexRateSTtUP, indexRateSTtDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexRateSTtWUP != -1 && indexRateSTtWDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexRateSTtWUP, indexRateSTtWDOWN));
    shiftUp += tmp*tmp;
  }
  if ( indexRateOtherUP != -1 && indexRateOtherDOWN != -1 )
  {
    tmp = isPositive(GetMaximum(indexRateOtherUP, indexRateOtherDOWN));
    shiftUp += tmp*tmp;
  }
  shiftUp = sqrt(shiftUp);
  
  
  shiftDown = 0.;
  if ( indexLepIdUP != -1 && indexLepIdDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexLepIdUP, indexLepIdDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexLepIsoUP != -1 && indexLepIsoDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexLepIsoUP, indexLepIsoDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexLepTrigUP != -1 && indexLepTrigDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexLepTrigUP, indexLepTrigDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexLepTrkUP != -1 && indexLepTrkDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexLepTrkUP, indexLepTrkDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexBTagUP != -1 && indexBTagDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexBTagUP, indexBTagDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexPuUP != -1 && indexPuDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexPuUP, indexPuDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexLumiUP != -1 && indexLumiDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexLumiUP, indexLumiDOWN));
    shiftDown += tmp*tmp;
  }
  if ( ! doMass && indexMassUP != -1 && indexMassDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexMassUP, indexMassDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexTopPt != -1 )
  {
    tmp = isNegative(systDifference[indexTopPt]);
    shiftDown += tmp*tmp;
  }
  if ( indexRenFac1002 != -1 || indexRenFac1003 != -1 || indexRenFac1004 != -1 || indexRenFac1005 != -1 || indexRenFac1007 != -1 || indexRenFac1009 != -1 || (indexIsrUP != -1 && indexIsrDOWN != -1) || (indexFsrUP != -1 && indexFsrDOWN != -1) )
  {
    tmp = isNegative(GetMinimum(indexRenFac1002, indexRenFac1003, indexRenFac1004, indexRenFac1005, indexRenFac1007, indexRenFac1009, indexIsrUP, indexIsrDOWN, indexFsrUP, indexFsrDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexHdampUP != -1 && indexHdampDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexHdampUP, indexHdampDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexFragCentral != -1 || indexFragUP != -1 || indexFragDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexFragCentral, indexFragUP, indexFragDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexFragPeterson != -1 )
  {
    tmp = isNegative(systDifference[indexFragPeterson]);
    shiftDown += tmp*tmp;
  }
  if ( indexFragSemiLepBrUP != -1 && indexFragSemiLepBrDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexFragSemiLepBrUP, indexFragSemiLepBrDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexPdfAlphaSUP != -1 && indexPdfAlphaSDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexPdfAlphaSUP, indexPdfAlphaSDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexPdfVar != -1 )
  {
    shiftDown += systDifference[indexPdfVar]*systDifference[indexPdfVar];
  }
  if ( indexTuneUP != -1 && indexTuneDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexTuneUP, indexTuneDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexCrErd != -1 || indexCrQcdErd != -1 || indexCrGluonMove != -1 || indexCrGluonMoveErd != -1 )
  {
    tmp = isNegative(GetMinimum(indexCrErd, indexCrQcdErd, indexCrGluonMove, indexCrGluonMoveErd));
    shiftDown += tmp*tmp;
  }
  if ( indexJesUP != -1 && indexJesDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexJesUP, indexJesDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexJerUP != -1 && indexJerDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexJerUP, indexJerDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexRateCMUP != -1 && indexRateCMDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexRateCMUP, indexRateCMDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexRateSTtUP != -1 && indexRateSTtDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexRateSTtUP, indexRateSTtDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexRateSTtWUP != -1 && indexRateSTtWDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexRateSTtWUP, indexRateSTtWDOWN));
    shiftDown += tmp*tmp;
  }
  if ( indexRateOtherUP != -1 && indexRateOtherDOWN != -1 )
  {
    tmp = isNegative(GetMinimum(indexRateOtherUP, indexRateOtherDOWN));
    shiftDown += tmp*tmp;
  }
  shiftDown = sqrt(shiftDown);
  
  cout << "Total systematic uncertainty: " << endl;
  cout << "Up: " << shiftUp << "; down: " << shiftDown << endl;
  
  // Test on fake data
  testData = ApplyCalibrationCurve(0.680374, 0.137884);  // comb
  //testData = ApplyCalibrationCurve(1.21739, 0.240619);   // mlb
  //testData = ApplyCalibrationCurve(0.638646, 0.131186);  // 180103
  //testData = ApplyCalibrationCurve(0.643469, 0.13191);   // 180107
  if (doMass) testData = ApplyCalibrationCurveMass(172.152, 0.107253);
  cout << "Values before calibration curve: " << endl;
  cout << "   MC:   width = " << nomVal << " +- " << nomValUnc << endl;
  cout << "   data: width = " << testData.first << " +- " << testData.second << " (stat.) +" << shiftUp << " -" << shiftDown << " (syst.)";
  if (doMass) cout << "   = +" << sqrt((testData.second)*(testData.second) + shiftUp*shiftUp) << " -" << sqrt((testData.second)*(testData.second) + shiftDown*shiftDown) << " (total)";
  cout << endl;
  
  //cout << "Combination of statistical and systematic uncertainty gives +" << TMath::Sqrt( (testData.second)* (testData.second) + shiftUp*shiftUp ) << " and -" << TMath::Sqrt( (testData.second)* (testData.second) + shiftDown*shiftDown ) << endl;
  
  if (! doMass) cout << endl << "=====> Data: " << 1.31*testData.first << " +- " << 1.31*testData.second << " (stat.) +" << 1.31*shiftUp << " -" << 1.31*shiftDown << " (syst.)   = +" << 1.31*sqrt((testData.second)*(testData.second) + shiftUp*shiftUp) << " -" << 1.31*sqrt((testData.second)*(testData.second) + shiftDown*shiftDown) << " (total)" << endl;
  
  WriteShiftTable();
  if (! doMass) WriteShiftTable(1.31);
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
  indexFragCentral       = -1;
  indexFragUP            = -1;
  indexFragDOWN          = -1;
  indexFragPeterson      = -1;
  indexFragSemiLepBrUP   = -1;
  indexFragSemiLepBrDOWN = -1;
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
  indexJesUP = -1; indexJesDOWN = -1;
  indexJerUP = -1; indexJerDOWN = -1;
  indexPdfAlphaSUP = -1; indexPdfAlphaSDOWN = -1;
  indexPdfVar = -1;
  indexRateCMUP = -1; indexRateCMDOWN = -1;
  indexRateSTtUP = -1; indexRateSTtDOWN = -1;
  indexRateSTtWUP = -1; indexRateSTtWDOWN = -1;
  indexRateOtherUP = -1; indexRateOtherDOWN = -1;
}

void ClearVars()
{
  ClearIndices();
  
  nomVal = 0.001;
  nomValUnc = 0.001;
  
  systName.clear();
  genWidths.clear();
  genMasses.clear();
  systWidthValue.clear();
  systWidthUnc.clear();
  systMassValue.clear();
  systMassUnc.clear();
  systDifference.clear();
  systRelDifference.clear();
  systDifferenceUnc.clear();
  systRelDifferenceUnc.clear();
  systMassDifference.clear();
  systRelMassDifference.clear();
  systMassDifferenceUnc.clear();
  systRelMassDifferenceUnc.clear();
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

double GetMinimum(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10)
{
  vector<double> v;
  v.clear();
  if ( i1 != -1 )
  {
    if ( isRate[i1] ) v.push_back(systDifference[i1]);
    else
    {
      if ( fabs(systDifference[i1]) >= fabs(systDifferenceUnc[i1]) ) v.push_back(systDifference[i1]);
      else v.push_back(systDifferenceUnc[i1]*systDifference[i1]/fabs(systDifference[i1]));
    }
  }
  if ( i2 != -1 )
  {
    if ( isRate[i2] ) v.push_back(systDifference[i2]);
    else
    {
      if ( fabs(systDifference[i2]) >= fabs(systDifferenceUnc[i2]) ) v.push_back(systDifference[i2]);
      else v.push_back(systDifferenceUnc[i2]*systDifference[i2]/fabs(systDifference[i2]));
    }
  }
  if ( i3 != -1 )
  {
    if ( isRate[i3] ) v.push_back(systDifference[i3]);
    else
    {
      if ( fabs(systDifference[i3]) >= fabs(systDifferenceUnc[i3]) ) v.push_back(systDifference[i3]);
      else v.push_back(systDifferenceUnc[i3]*systDifference[i3]/fabs(systDifference[i3]));
    }
  }
  if ( i4 != -1 )
  {
    if ( isRate[i4] ) v.push_back(systDifference[i4]);
    else
    {
      if ( fabs(systDifference[i4]) >= fabs(systDifferenceUnc[i4]) ) v.push_back(systDifference[i4]);
      else v.push_back(systDifferenceUnc[i4]*systDifference[i4]/fabs(systDifference[i4]));
    }
  }
  if ( i5 != -1 )
  {
    if ( isRate[i5] ) v.push_back(systDifference[i5]);
    else
    {
      if ( fabs(systDifference[i5]) >= fabs(systDifferenceUnc[i5]) ) v.push_back(systDifference[i5]);
      else v.push_back(systDifferenceUnc[i5]*systDifference[i5]/fabs(systDifference[i5]));
    }
  }
  if ( i6 != -1 )
  {
    if ( isRate[i6] ) v.push_back(systDifference[i6]);
    else
    {
      if ( fabs(systDifference[i6]) >= fabs(systDifferenceUnc[i6]) ) v.push_back(systDifference[i6]);
      else v.push_back(systDifferenceUnc[i6]*systDifference[i6]/fabs(systDifference[i6]));
    }
  }
  if ( i7 != -1 )
  {
    if ( isRate[i7] ) v.push_back(systDifference[i7]);
    else
    {
      if ( fabs(systDifference[i7]) >= fabs(systDifferenceUnc[i7]) ) v.push_back(systDifference[i7]);
      else v.push_back(systDifferenceUnc[i7]*systDifference[i7]/fabs(systDifference[i7]));
    }
  }
  if ( i8 != -1 )
  {
    if ( isRate[i8] ) v.push_back(systDifference[i8]);
    else
    {
      if ( fabs(systDifference[i8]) >= fabs(systDifferenceUnc[i8]) ) v.push_back(systDifference[i8]);
      else v.push_back(systDifferenceUnc[i8]*systDifference[i8]/fabs(systDifference[i8]));
    }
  }
  if ( i9 != -1 )
  {
    if ( isRate[i9] ) v.push_back(systDifference[i9]);
    else
    {
      if ( fabs(systDifference[i9]) >= fabs(systDifferenceUnc[i9]) ) v.push_back(systDifference[i9]);
      else v.push_back(systDifferenceUnc[i9]*systDifference[i9]/fabs(systDifference[i9]));
    }
  }
  if ( i10 != -1 )
  {
    if ( isRate[i10] ) v.push_back(systDifference[i10]);
    else
    {
      if ( fabs(systDifference[i10]) >= fabs(systDifferenceUnc[i10]) ) v.push_back(systDifference[i10]);
      else v.push_back(systDifferenceUnc[i10]*systDifference[i10]/fabs(systDifference[i10]));
    }
  }
  
  if ( v.size() == 0 ) return 0.;
  else if ( v.size() == 1 ) return v[0];
  
  double min = *min_element(v.begin(), v.end());
  if ( fabs(min) < 1e-6 ) min = 0.;
  
  return min;
}

double GetMaximum(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10)
{
  vector<double> v;
  v.clear();
  if ( i1 != -1 )
  {
    if ( isRate[i1] ) v.push_back(systDifference[i1]);
    else
    {
      if ( fabs(systDifference[i1]) >= fabs(systDifferenceUnc[i1]) ) v.push_back(systDifference[i1]);
      else v.push_back(systDifferenceUnc[i1]*systDifference[i1]/fabs(systDifference[i1]));
    }
  }
  if ( i2 != -1 )
  {
    if ( isRate[i2] ) v.push_back(systDifference[i2]);
    else
    {
      if ( fabs(systDifference[i2]) >= fabs(systDifferenceUnc[i2]) ) v.push_back(systDifference[i2]);
      else v.push_back(systDifferenceUnc[i2]*systDifference[i2]/fabs(systDifference[i2]));
    }
  }
  if ( i3 != -1 )
  {
    if ( isRate[i3] ) v.push_back(systDifference[i3]);
    else
    {
      if ( fabs(systDifference[i3]) >= fabs(systDifferenceUnc[i3]) ) v.push_back(systDifference[i3]);
      else v.push_back(systDifferenceUnc[i3]*systDifference[i3]/fabs(systDifference[i3]));
    }
  }
  if ( i4 != -1 )
  {
    if ( isRate[i4] ) v.push_back(systDifference[i4]);
    else
    {
      if ( fabs(systDifference[i4]) >= fabs(systDifferenceUnc[i4]) ) v.push_back(systDifference[i4]);
      else v.push_back(systDifferenceUnc[i4]*systDifference[i4]/fabs(systDifference[i4]));
    }
  }
  if ( i5 != -1 )
  {
    if ( isRate[i5] ) v.push_back(systDifference[i5]);
    else
    {
      if ( fabs(systDifference[i5]) >= fabs(systDifferenceUnc[i5]) ) v.push_back(systDifference[i5]);
      else v.push_back(systDifferenceUnc[i5]*systDifference[i5]/fabs(systDifference[i5]));
    }
  }
  if ( i6 != -1 )
  {
    if ( isRate[i6] ) v.push_back(systDifference[i6]);
    else
    {
      if ( fabs(systDifference[i6]) >= fabs(systDifferenceUnc[i6]) ) v.push_back(systDifference[i6]);
      else v.push_back(systDifferenceUnc[i6]*systDifference[i6]/fabs(systDifference[i6]));
    }
  }
  if ( i7 != -1 )
  {
    if ( isRate[i7] ) v.push_back(systDifference[i7]);
    else
    {
      if ( fabs(systDifference[i7]) >= fabs(systDifferenceUnc[i7]) ) v.push_back(systDifference[i7]);
      else v.push_back(systDifferenceUnc[i7]*systDifference[i7]/fabs(systDifference[i7]));
    }
  }
  if ( i8 != -1 )
  {
    if ( isRate[i8] ) v.push_back(systDifference[i8]);
    else
    {
      if ( fabs(systDifference[i8]) >= fabs(systDifferenceUnc[i8]) ) v.push_back(systDifference[i8]);
      else v.push_back(systDifferenceUnc[i8]*systDifference[i8]/fabs(systDifference[i8]));
    }
  }
  if ( i9 != -1 )
  {
    if ( isRate[i9] ) v.push_back(systDifference[i9]);
    else
    {
      if ( fabs(systDifference[i9]) >= fabs(systDifferenceUnc[i9]) ) v.push_back(systDifference[i9]);
      else v.push_back(systDifferenceUnc[i9]*systDifference[i9]/fabs(systDifference[i9]));
    }
  }
  if ( i10 != -1 )
  {
    if ( isRate[i10] ) v.push_back(systDifference[i10]);
    else
    {
      if ( fabs(systDifference[i10]) >= fabs(systDifferenceUnc[i10]) ) v.push_back(systDifference[i10]);
      else v.push_back(systDifferenceUnc[i10]*systDifference[i10]/fabs(systDifference[i10]));
    }
  }
  
  if ( v.size() == 0 ) return 0.;
  else if ( v.size() == 1 ) return v[0];
  
  double max = *max_element(v.begin(), v.end());
  if ( fabs(max) < 1e-6 ) max = 0.;
  
  return max;
}

double isPositive(double upVar)
{
  if ( upVar >= 0. ) return upVar;
  else return 0.;
}

double isNegative(double downVar)
{
  if ( downVar <= 0. ) return downVar;
  else return 0.;
}

void IndexSystematics()
{
  for (int iSys = 0; iSys < nSys; iSys++)
  {
    if ( systName[iSys].find("nominal") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( genWidths[iSys] == 1 && genMasses[iSys] == 172.5 ) indexNom = iSys;
      else if ( genWidths[iSys] == 1 && genMasses[iSys] == 171.5 ) indexMassDOWN = iSys;
      else if ( genWidths[iSys] == 1 && genMasses[iSys] == 173.5 ) indexMassUP = iSys;
    }
    else if ( systName[iSys].find("lepton") != std::string::npos || systName[iSys].find("Lepton") != std::string::npos )
    {
      isRate[iSys] = true;
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
      isRate[iSys] = true;
      if ( foundUP(systName[iSys]) ) indexBTagUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexBTagDOWN = iSys;
    }
    else if ( systName[iSys].find("hdamp") != std::string::npos || systName[iSys].find("Hdamp") != std::string::npos )
    {
      isRate[iSys] = false;
      if ( foundUP(systName[iSys]) ) indexHdampUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexHdampDOWN = iSys;
    }
    else if ( systName[iSys].find("Pu") != std::string::npos || systName[iSys].find("PU") != std::string::npos || systName[iSys].find("pu") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( foundUP(systName[iSys]) ) indexPuUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexPuDOWN = iSys;
    }
    else if ( systName[iSys].find("Lumi") != std::string::npos || systName[iSys].find("LUMI") != std::string::npos || systName[iSys].find("lumi") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( foundUP(systName[iSys]) ) indexLumiUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexLumiDOWN = iSys;
    }
    else if ( systName[iSys].find("topPt") != std::string::npos || systName[iSys].find("TopPt") != std::string::npos || systName[iSys].find("TOPPT") != std::string::npos || systName[iSys].find("toppt") != std::string::npos )
    {
      isRate[iSys] = true;
      indexTopPt = iSys;
    }
    else if ( systName[iSys].find("renFac") != std::string::npos || systName[iSys].find("RenFac") != std::string::npos || systName[iSys].find("RENFAC") != std::string::npos || systName[iSys].find("renfac") != std::string::npos )
    {
      isRate[iSys] = true;
      if      ( systName[iSys].find("1002") != std::string::npos ) indexRenFac1002 = iSys;
      else if ( systName[iSys].find("1003") != std::string::npos ) indexRenFac1003 = iSys;
      else if ( systName[iSys].find("1004") != std::string::npos ) indexRenFac1004 = iSys;
      else if ( systName[iSys].find("1005") != std::string::npos ) indexRenFac1005 = iSys;
      else if ( systName[iSys].find("1007") != std::string::npos ) indexRenFac1007 = iSys;
      else if ( systName[iSys].find("1009") != std::string::npos ) indexRenFac1009 = iSys;
    }
    else if ( systName[iSys].find("frag") != std::string::npos || systName[iSys].find("Frag") != std::string::npos || systName[iSys].find("FRAG") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( systName[iSys].find("semiLepBr") != std::string::npos || systName[iSys].find("SemiLepBr") != std::string::npos )
      {
        if ( foundUP(systName[iSys]) ) indexFragSemiLepBrUP = iSys;
        else if ( foundDOWN(systName[iSys]) ) indexFragSemiLepBrDOWN = iSys;
      }
      else if ( systName[iSys].find("central") != std::string::npos || systName[iSys].find("Central") != std::string::npos ) indexFragCentral = iSys;
      else if ( systName[iSys].find("peterson") != std::string::npos || systName[iSys].find("Peterson") != std::string::npos ) indexFragPeterson = iSys;
      else if ( foundUP(systName[iSys]) )        indexFragUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexFragDOWN = iSys;
    }
    else if ( systName[iSys].find("pdf") != std::string::npos || systName[iSys].find("Pdf") != std::string::npos || systName[iSys].find("PDF") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( systName[iSys].find("alpha") != std::string::npos || systName[iSys].find("Alpha") != std::string::npos )
      {
        if ( foundUP(systName[iSys]) ) indexPdfAlphaSUP = iSys;
        else if ( foundDOWN(systName[iSys]) ) indexPdfAlphaSDOWN = iSys;
      }
      else if ( systName[iSys].find("pdfVar") != std::string::npos ) indexPdfVar = iSys;
    }
    else if ( systName[iSys].find("isr") != std::string::npos || systName[iSys].find("ISR") != std::string::npos || systName[iSys].find("Isr") != std::string::npos )
    {
      isRate[iSys] = false;
      if ( foundUP(systName[iSys]) ) indexIsrUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexIsrDOWN = iSys;
    }
    else if ( systName[iSys].find("fsr") != std::string::npos || systName[iSys].find("FSR") != std::string::npos || systName[iSys].find("Fsr") != std::string::npos )
    {
      isRate[iSys] = false;
      if ( foundUP(systName[iSys]) ) indexFsrUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexFsrDOWN = iSys;
    }
//     else if ( systName[iSys].find("mass") != std::string::npos || systName[iSys].find("MASS") != std::string::npos || systName[iSys].find("Mass") != std::string::npos )
//     {
//       if ( systName[iSys].find("169") != std::string::npos ) indexMassDOWN = iSys;
//       //if ( systName[iSys].find("171") != std::string::npos ) indexMassDOWN = iSys;
//       //else if ( systName[iSys].find("173") != std::string::npos ) indexMassUP = iSys;
//       else if ( systName[iSys].find("175") != std::string::npos ) indexMassUP = iSys;
//     }
    else if ( systName[iSys].find("tune") != std::string::npos )
    {
      isRate[iSys] = false;
      if ( foundUP(systName[iSys]) ) indexTuneUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexTuneDOWN = iSys;
    }
    else if ( systName[iSys].find("ERD") != std::string::npos )
    {
      isRate[iSys] = false;
      if ( systName[iSys].find("mpi") != std::string::npos ) indexCrErd = iSys;
      else if ( systName[iSys].find("qcd") != std::string::npos ) indexCrQcdErd = iSys;
      else if ( systName[iSys].find("gluon") != std::string::npos ) indexCrGluonMoveErd = iSys;
    }
    else if ( systName[iSys].find("gluon") != std::string::npos )
    {
      isRate[iSys] = false;
      indexCrGluonMove = iSys;
    }
    else if ( systName[iSys].find("jes") != std::string::npos || systName[iSys].find("JES") != std::string::npos || systName[iSys].find("Jes") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( foundUP(systName[iSys]) ) indexJesUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexJesDOWN = iSys;
    }
    else if ( systName[iSys].find("jer") != std::string::npos || systName[iSys].find("JER") != std::string::npos || systName[iSys].find("Jer") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( foundUP(systName[iSys]) ) indexJerUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexJerDOWN = iSys;
    }
    else if ( systName[iSys].find("herwig") != std::string::npos || systName[iSys].find("Herwig") != std::string::npos )
    {
      isRate[iSys] = false;
      indexHerwig = iSys;
    }
    else if ( systName[iSys].find("rateGood") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( foundUP(systName[iSys]) ) indexRateCMUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexRateCMDOWN = iSys;
    }
    else if ( systName[iSys].find("rateSTtW") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( foundUP(systName[iSys]) ) indexRateSTtWUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexRateSTtWDOWN = iSys;
    }
    else if ( systName[iSys].find("rateSTt") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( foundUP(systName[iSys]) ) indexRateSTtUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexRateSTtDOWN = iSys;
    }
    else if ( systName[iSys].find("rateOther") != std::string::npos )
    {
      isRate[iSys] = true;
      if ( foundUP(systName[iSys]) ) indexRateOtherUP = iSys;
      else if ( foundDOWN(systName[iSys]) ) indexRateOtherDOWN = iSys;
    }
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
  cout << "indexFragSemiLepBrUP  " << indexFragSemiLepBrUP << endl;
  cout << "indexFragSemiLepBrDOWN  " << indexFragSemiLepBrDOWN << endl;
  cout << "indexFragUP  " << indexFragUP << endl;
  cout << "indexFragDOWN  " << indexFragDOWN << endl;
  cout << "indexHdampUP  " << indexHdampUP << endl;
  cout << "indexHdampDOWN  " << indexHdampDOWN << endl;
}

std::pair<double,double> ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma)
{
  /// return 'thisInputWidth' and uncertainty
  //  thisOutputWidth = Par_[0] + Par_[1] * thisInputWidth
  
  double thisInputWidth = (thisOutputWidth - calCurvePar_[0])/calCurvePar_[1];
  double thisInputWidthSigma = TMath::Sqrt( thisOutputWidthSigma*thisOutputWidthSigma + calCurveParUnc_[0]*calCurveParUnc_[0] + calCurveParUnc_[1]*calCurveParUnc_[1]*thisInputWidth*thisInputWidth )/calCurvePar_[1];
  
  return std::pair<double,double>(thisInputWidth,thisInputWidthSigma);
}

std::pair<double,double> ApplyCalibrationCurveMass(double thisOutputMass, double thisOutputMassSigma)
{
  /// return 'thisInputMass' and uncertainty
  //  thisOutputMass = ParMass_[0] + ParMass_[1] * thisInputMass
  
  double thisInputMass = (thisOutputMass - calCurveParMass_[0])/calCurveParMass_[1];
  //double thisInputMassSigma = TMath::Sqrt( thisOutputMassSigma*thisOutputMassSigma + calCurveParUncMass_[0]*calCurveParUncMass_[0] + calCurveParUncMass_[1]*calCurveParUncMass_[1]*thisInputMass*thisInputMass )/calCurveParMass_[1];
  
  //return std::pair<double,double>(thisInputMass,thisInputMassSigma);
  return std::pair<double,double>(thisInputMass,thisOutputMassSigma);
}

std::pair<double,double> ApplyCalibrationCurveW(double thisOutputWidth, double thisOutputWidthSigma)
{
  /// return 'thisInputWidth' and uncertainty
  //  thisOutputWidth = ParW_[0] + ParW_[1] * thisInputWidth
  
  double thisInputWidth = (thisOutputWidth - calCurveParW_[0])/calCurveParW_[1];
  double thisInputWidthSigma = TMath::Sqrt( thisOutputWidthSigma*thisOutputWidthSigma + calCurveParUncW_[0]*calCurveParUncW_[0] + calCurveParUncW_[1]*calCurveParUncW_[1]*thisInputWidth*thisInputWidth )/calCurveParW_[1];
  
  return std::pair<double,double>(thisInputWidth,thisInputWidthSigma);
}

std::pair<double,double> ApplyCalibrationCurveM(double thisOutputMass, double thisOutputMassSigma)
{
  /// return 'thisInputMass' and uncertainty
  //  thisOutputMass = ParM_[0] + ParM_[1] * thisInputMass
  
  double thisInputMass = (thisOutputMass - calCurveParM_[0])/calCurveParM_[1];
  double thisInputMassSigma = TMath::Sqrt( thisOutputMassSigma*thisOutputMassSigma + calCurveParUncM_[0]*calCurveParUncM_[0] + calCurveParUncM_[1]*calCurveParUncM_[1]*thisInputMass*thisInputMass )/calCurveParM_[1];
  
  return std::pair<double,double>(thisInputMass,thisInputMassSigma);
  
}

void WriteShiftTable(double scale)
{
  string space = "     ";
  string hline = space+"\\hline";
  string interline = space+"&&&&\\\\[-6pt]";
  string headerextra = "[+3pt]";
  string multicolBegin = "\\multicolumn{2}{r}{";
  string multicolEnd = "$\\qquad$}";
  
  string fileName = "systematicsShiftFull";
  if ( scale != 1. ) fileName += "x"+DotReplace(scale);
  ofstream fileOut((fileName+suffix+".tex").c_str());
  fileOut << "\\begin{table}[htp]" << endl;
  fileOut << " \\caption{\\fixme{Systematics compared to $" << scale << "\\times s$, for $m_r = m_t/<m_t>_{\\CM}$, range 0.65 to 1.4}}" << endl;
  fileOut << " \\begin{center}" << endl;
  fileOut << "  \\begin{minipage}[p]{\\textwidth}" << endl;
  fileOut << "   \\begin{center}" << endl;
  fileOut << "    %\\hspace{-2.6cm}" << endl;
  fileOut << "    {\\small" << endl;
  fileOut << "    \\begin{tabular}{l S[table-format=1.4] S[table-format=1.4] S[table-format=2.4] S[table-format=1.4]}" << endl;  /// adapt!!
  
  fileOut << hline << endl;
  fileOut << hline << endl;
  fileOut << interline << endl;
  fileOut << space << "\\multicolumn{3}{l}{Nominal value} & \\multicolumn{2}{l}{$" << fixed << setprecision(4) << scale*nomVal << " \\pm " << scale*nomValUnc << "$} \\\\" << headerextra << endl;
  fileOut << hline << endl;
  fileOut << interline << endl;
  fileOut << space << "Systematic & \\multicolumn{2}{c}{$+1\\sigma$} & \\multicolumn{2}{c}{$-1\\sigma$} \\\\" << endl;
  fileOut << space << " & {Shift} & {   Unc   } & {Shift} & {   Unc   } \\\\" << headerextra << endl;
  fileOut << hline << endl;
  
  fileOut << interline << endl;
  if ( (indexLepIdUP != -1 && indexLepIdDOWN != -1) || (indexLepIsoUP != -1 && indexLepIsoDOWN != -1) || (indexLepTrigUP != -1 && indexLepTrigDOWN != -1) || (indexLepTrkUP != -1 && indexLepTrkDOWN != -1) )
    fileOut << space << "\\textit{Lepton SFs} &  &  &  &  \\\\" << endl;
  if ( indexLepIdUP != -1 && indexLepIdDOWN != -1 )
    fileOut << space << "\\tabsp Id        & " << fixed << setprecision(4) << scale*systDifference[indexLepIdUP] << " &  & " << scale*systDifference[indexLepIdDOWN] << " &  \\\\" << endl;
  if ( indexLepIsoUP != -1 && indexLepIsoDOWN != -1 )
    fileOut << space << "\\tabsp Isolation & " << fixed << setprecision(4) << scale*systDifference[indexLepIsoUP] << " &  & " << scale*systDifference[indexLepIsoDOWN] << " &  \\\\" << endl;
  if ( indexLepTrigUP != -1 && indexLepTrigDOWN != -1 )
    fileOut << space << "\\tabsp Trigger   & " << fixed << setprecision(4) << scale*systDifference[indexLepTrigUP] << " &  & " << scale*systDifference[indexLepTrigDOWN] << " &  \\\\" << endl;
  if ( indexLepTrkUP != -1 && indexLepTrkDOWN != -1 )
    fileOut << space << "\\tabsp Tracking  & " << fixed << setprecision(4) << scale*systDifference[indexLepTrkUP] << " &  & " << scale*systDifference[indexLepTrkDOWN] << " &  \\\\" << endl;
  
  if ( indexBTagUP != -1 && indexBTagDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\bq~tagging SFs & " << fixed << setprecision(4) << scale*systDifference[indexBTagUP] << " &  & " << scale*systDifference[indexBTagDOWN] << " &  \\\\" << endl;
  }
  
  if ( indexPuUP != -1 && indexPuDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << fixed << setprecision(4) << "\\Pileup\\ SFs & " << scale*systDifference[indexPuUP] << " &  & " << scale*systDifference[indexPuDOWN] << " &  \\\\" << endl;
  }
  
  if ( indexLumiUP != -1 && indexLumiDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Luminosity & " << fixed << setprecision(4) << scale*systDifference[indexLumiUP] << " &  & " << scale*systDifference[indexLumiDOWN] << " &  \\\\" << endl;
  }
  
  if ( ! doMass && indexMassUP != -1 && indexMassDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Top mass & " << fixed << setprecision(4) << scale*systDifference[indexMassUP] << " &  & " << scale*systDifference[indexMassDOWN] << " &  \\\\" << endl;
  }
  
  if ( indexTopPt != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Top $\\pT$ reweighting &  &  & " << fixed << setprecision(4) << scale*systDifference[indexTopPt] << " &  \\\\" << endl;
  }
  
  if ( indexRenFac1002 != -1 || indexRenFac1003 != -1 || indexRenFac1004 != -1 || indexRenFac1005 != -1 || indexRenFac1007 != -1 || indexRenFac1009 != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{Renormalisation/factorisation scale} &  &  &  &  \\\\" << endl;
    if ( indexRenFac1002 != -1 )
      fileOut << space << "\\tabsp $\\muR = 1\\,,\\quad \\muF = 2$   &  &  & " << scale*systDifference[indexRenFac1002] << " &  \\\\" << endl;
    if ( indexRenFac1003 != -1 )
      fileOut << space << "\\tabsp $\\muR = 1\\,,\\quad \\muF = 0.5$ &  &  & " << scale*systDifference[indexRenFac1003] << " &  \\\\" << endl;
    if ( indexRenFac1004 != -1 )
      fileOut << space << "\\tabsp $\\muR = 2\\,,\\quad \\muF = 1$   &  &  & " << scale*systDifference[indexRenFac1004] << " &  \\\\" << endl;
    if ( indexRenFac1005 != -1 )
      fileOut << space << "\\tabsp $\\muR = 2\\,,\\quad \\muF = 2$   &  &  & " << scale*systDifference[indexRenFac1005] << " &  \\\\" << endl;
    if ( indexRenFac1007 != -1 )
      fileOut << space << "\\tabsp $\\muR = 0.5\\,,\\, \\muF = 1$    &  &  & " << scale*systDifference[indexRenFac1007] << " &  \\\\" << endl;
    if ( indexRenFac1009 != -1 )
      fileOut << space << "\\tabsp $\\muR = 0.5\\,,\\, \\muF = 0.5$  &  &  & " << scale*systDifference[indexRenFac1009] << " &  \\\\" << endl;
  }
  
  if ( indexIsrUP != -1 && indexIsrDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "ISR & " << scale*systDifference[indexIsrUP] << " & " << scale*systDifferenceUnc[indexIsrUP] << " & " << scale*systDifference[indexIsrDOWN] << " & " << scale*systDifferenceUnc[indexIsrDOWN] << " \\\\" << endl;
  }
  if ( indexFsrUP != -1 && indexFsrDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "FSR & " << scale*systDifference[indexFsrUP] << " & " << scale*systDifferenceUnc[indexFsrUP] << " & " << scale*systDifference[indexFsrDOWN] << " & " << scale*systDifferenceUnc[indexFsrDOWN] << " \\\\" << endl;
  }
  if ( indexHdampUP != -1 && indexHdampDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\hdamp & " << scale*systDifference[indexHdampUP] << " & " << scale*systDifferenceUnc[indexHdampUP] << " & " << scale*systDifference[indexHdampDOWN] << " & " << scale*systDifferenceUnc[indexHdampDOWN] << " \\\\" << endl;
  }
  
  if ( indexFragCentral != -1 || indexFragUP != -1 || indexFragDOWN != -1 || indexFragPeterson != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{\\bq\\ fragmentation} &  &  &  &  \\\\" << endl;
    if ( indexFragUP != -1 && indexFragDOWN != -1 )
      fileOut << space << "\\tabsp Bowler--Lund & " << scale*systDifference[indexFragUP] << " &  & " << scale*systDifference[indexFragDOWN] << " &  \\\\" << endl;
    if ( indexFragCentral != -1 )
      fileOut << space << "\\tabsp central &  &  & " << scale*systDifference[indexFragCentral] << " &  \\\\" << endl;
    if ( indexFragPeterson != -1 )
      fileOut << space << "\\tabsp Peterson &  &  & " << scale*systDifference[indexFragPeterson] << " &  \\\\" << endl;
  }
  
  if ( indexFragSemiLepBrUP != -1 && indexFragSemiLepBrDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\bq\\ \\semilep\\ BR & " << scale*systDifference[indexFragSemiLepBrUP] << " &  & " << scale*systDifference[indexFragSemiLepBrDOWN] << " &  \\\\" << endl;
  }
  
  if ( (indexPdfAlphaSUP != -1 && indexPdfAlphaSDOWN != -1) || indexPdfVar != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{PDF} &  &  &  &  \\\\" << endl;
    if ( indexPdfVar != -1 )
      fileOut << space << "\\tabsp Replicas & " << scale*systDifference[indexPdfVar] << " &  & -" << scale*systDifference[indexPdfVar] << " &  \\\\" << endl;
    if ( indexPdfAlphaSUP != -1 && indexPdfAlphaSDOWN != -1 )
      fileOut << space << "\\tabsp $\\alpha_{s}$ & " << scale*systDifference[indexPdfAlphaSUP] << " &  & " << scale*systDifference[indexPdfAlphaSDOWN] << " &  \\\\" << endl;
  }
  
  if ( indexTuneUP != -1 && indexTuneDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "Tune & " << scale*systDifference[indexTuneUP] << " & " << scale*systDifferenceUnc[indexTuneUP] << " & " << scale*systDifference[indexTuneDOWN] << " & " << scale*systDifferenceUnc[indexTuneDOWN] << " \\\\" << endl;
  }
  if ( indexCrErd != -1 || indexCrQcdErd != -1 || indexCrGluonMove != -1 || indexCrGluonMoveErd != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{Colour reconnection} &  &  &  &  \\\\" << endl;
    if ( indexCrErd != -1 )
      fileOut << space << "\\tabsp MPI (ERD)        &  &  & " << scale*systDifference[indexCrErd] << " & " << scale*systDifferenceUnc[indexCrErd] << " \\\\" << endl;
    if ( indexCrQcdErd != -1 )
      fileOut << space << "\\tabsp QCD-based (ERD)  &  &  & " << scale*systDifference[indexCrQcdErd] << " & " << scale*systDifferenceUnc[indexCrQcdErd] << " \\\\" << endl;
    if ( indexCrGluonMove != -1 )
      fileOut << space << "\\tabsp Gluon move       &  &  & " << scale*systDifference[indexCrGluonMove] << " & " << scale*systDifferenceUnc[indexCrGluonMove] << " \\\\" << endl;
    if ( indexCrGluonMoveErd != -1 )
      fileOut << space << "\\tabsp Gluon move (ERD) &  &  & " << scale*systDifference[indexCrGluonMoveErd] << " & " << scale*systDifferenceUnc[indexCrGluonMoveErd] << " \\\\" << endl;
  }
  
//   if ( indexHerwig != -1 )
//   {
//     fileOut << interline << endl;
//     fileOut << space << "Herwig &  &  & " << scale*systDifference[indexHerwig] << " & " << scale*systDifferenceUnc[indexHerwig] << " \\\\" << endl;
//   }
  
  if ( indexJesUP != -1 && indexJesDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "JES & " << scale*systDifference[indexJesUP] << " &  & " << scale*systDifference[indexJesDOWN] << " &  \\\\" << endl;
    //fileOut << space << "JES & " << scale*systDifference[indexJesUP] << " & " << scale*systDifferenceUnc[indexJesUP] << " & " << scale*systDifference[indexJesDOWN] << " & " << scale*systDifferenceUnc[indexJesDOWN] << " \\\\" << endl;
  }
  if ( indexJerUP != -1 && indexJerDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "JER & " << scale*systDifference[indexJerUP] << " &  & " << scale*systDifference[indexJerDOWN] << " &  \\\\" << endl;
    //fileOut << space << "JER & " << scale*systDifference[indexJerUP] << " & " << scale*systDifferenceUnc[indexJerUP] << " & " << scale*systDifference[indexJerDOWN] << " & " << scale*systDifferenceUnc[indexJerDOWN] << " \\\\" << endl;
  }
  
  if ( ( indexRateCMUP != -1 && indexRateCMDOWN != -1 ) || ( indexRateSTtUP != -1 && indexRateSTtDOWN != -1 ) || ( indexRateSTtWUP != -1 && indexRateSTtWDOWN != -1 ) || ( indexRateOtherUP != -1 && indexRateOtherDOWN != -1 ) )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{Cross section var.} &  &  &  &  \\\\" << endl;
    if ( indexRateCMUP != -1 && indexRateCMDOWN != -1 )
      fileOut << space << "\\tabsp \\# \\CM events & " << scale*systDifference[indexRateCMUP] << " &  & -" << scale*systDifference[indexRateCMDOWN] << " &  \\\\" << endl;
    if ( indexRateSTtUP != -1 && indexRateSTtDOWN != -1 )
      fileOut << space << "\\tabsp ST \\tq channel & " << scale*systDifference[indexRateSTtUP] << " &  & -" << scale*systDifference[indexRateSTtDOWN] << " &  \\\\" << endl;
    if ( indexRateSTtWUP != -1 && indexRateSTtWDOWN != -1 )
      fileOut << space << "\\tabsp ST \\tq\\W channel & " << scale*systDifference[indexRateSTtWUP] << " &  & -" << scale*systDifference[indexRateSTtWDOWN] << " &  \\\\" << endl;
    if ( indexRateOtherUP != -1 && indexRateOtherDOWN != -1 )
      fileOut << space << "\\tabsp Other & " << scale*systDifference[indexRateOtherUP] << " &  & " << scale*systDifference[indexRateOtherDOWN] << " &  \\\\" << endl;
  }
  
  //....
  
  fileOut << interline << endl;
  fileOut << hline << endl;
  fileOut << interline << endl;
  fileOut << space << "Total & \\multicolumn{2}{c}{" << scale*shiftUp << "} & \\multicolumn{2}{c}{-" << scale*shiftDown << "} \\\\" << endl;
  
  fileOut << interline << endl;
  fileOut << hline << endl;
  fileOut << hline << endl;
  
  fileOut << "    \\end{tabular}}" << endl;
  fileOut << "   \\end{center}" << endl;
  fileOut << "  \\end{minipage}" << endl;
  fileOut << " \\end{center}" << endl;
  fileOut << " \\label{tab:systFullShift";
  if ( scale != 1. ) fileOut << "x"+DotReplace(scale);
  fileOut << "}" << endl;
  fileOut << "\\end{table}%" << endl;
  
  cout << "Written table with absolute shifts + stat. uncertainties" << endl;
  
}

void WriteRelTable()
{
  string space = "     ";
  string hline = space+"\\hline";
  string interline = space+"&&\\\\[-6pt]";
  string headerextra = "[+3pt]";
  
  
  ofstream fileOut(("systematicsRelFull"+suffix+".tex").c_str());
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
  
  if ( ! doMass && indexMassUP != -1 && indexMassDOWN != -1 )
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
  
  if ( indexFragCentral != -1 || indexFragUP != -1 || indexFragDOWN != -1 || indexFragPeterson != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\textit{\\bq\\ fragmentation} &  &  \\\\" << endl;
    if ( indexFragUP != -1 && indexFragDOWN != -1 )
      fileOut << space << "\\tabsp up & " << systRelDifference[indexFragUP] << " & " << systRelDifference[indexFragDOWN] << " \\\\" << endl;
    if ( indexFragCentral != -1 )
      fileOut << space << "\\tabsp central & \\multicolumn{2}{r}{" << systRelDifference[indexFragCentral] << " & " << "$\\qquad\\quad$} \\\\" << endl;
    if ( indexFragPeterson != -1 )
      fileOut << space << "\\tabsp Peterson & \\multicolumn{2}{r}{" << systRelDifference[indexFragPeterson] << " & " << "$\\qquad\\quad$} \\\\" << endl;
  }
  
  if ( indexFragSemiLepBrUP != -1 && indexFragSemiLepBrDOWN != -1 )
  {
    fileOut << interline << endl;
    fileOut << space << "\\bq\\ \\semilep\\ BR & " << systRelDifference[indexFragSemiLepBrUP] << " & " << systRelDifference[indexFragSemiLepBrDOWN] << " \\\\" << endl;
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
