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


bool is2D = false;
bool do95 = false;

bool useComb = true;
bool useLep = false;
bool useHad = false;

string inputFileDir = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/Systematics/temp/";
string prefix = "";

/// 1D width
// parameter = comb
const double calCurveParComb[2] = {-0.249918, 1.08351};
const double calCurveParUncComb[2] = {0.0351629, 0.0139145};
// parameter = redMlb
const double calCurveParLep[2] = {-0.267184, 1.08957};
const double calCurveParUncLep[2] = {0.0550971, 0.0248913};
// parameter = redTopMass_old
const double calCurveParHad[2] = {-0.152752, 1.09383};
const double calCurveParUncHad[2] = {0.0436359, 0.0166395};

double calCurvePar_[2] = {0., 1.};
double calCurveParUnc_[2] = {0., 0.};

ifstream fileIn;
ofstream fileOut;

string line, fileName;
string thisSyst;
double thisWidth, thisMass, thisWidthValue, thisWidthUnc, thisMassValue, thisMassUnc;
pair<double,double> nomWidthValue, thisCorrWidthValue, thisCorrMassValue;
double diff, uncDiff;
double inputWidth = 1., inputMass = 172.5;

bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile.good();
}

std::pair<double,double> ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma, string variation = "", int direction = 0)
{
  /// return 'thisInputWidth' and uncertainty
  //  thisOutputWidth = Par_[0] + Par_[1] * thisInputWidth
  double par0 = calCurvePar_[0], par1 = calCurvePar_[1];
  double uncPar0 = calCurveParUnc_[0], uncPar1 = calCurveParUnc_[1];
  double nSigma = 1.;
  if (do95) nSigma = 2.;
  if ( ! variation.empty() && direction != 0 )
  {
    if ( variation.find("const") != std::string::npos )
    {
      if ( direction == 1 )      par0 += nSigma * uncPar0;
      else if (direction == -1 ) par0 -= nSigma * uncPar0;
    }
    else if ( variation.find("slope") != std::string::npos )
    {
      if ( direction == 1 )      par1 += nSigma * uncPar1;
      else if (direction == -1 ) par1 -= nSigma * uncPar1;
    }
  }
  
  double thisInputWidth = (thisOutputWidth - par0)/par1;
  double thisInputWidthSigma = TMath::Sqrt( thisOutputWidthSigma*thisOutputWidthSigma + uncPar0*uncPar0 + uncPar1*uncPar1*thisInputWidth*thisInputWidth )/par1;
  
  return std::pair<double,double>(thisInputWidth,thisInputWidthSigma);
}

int main()
{
  cout << "**** Beginning of the program ****" << endl;
  
  if (do95)
  {
    prefix = "sigma95_";
    inputFileDir = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/Systematics/temp95/";
  }
  
  /// Set up
  if (useComb)
  {
    if (! do95) inputFileDir += "180110_comb/";
    calCurvePar_[0] = calCurveParComb[0];
    calCurvePar_[1] = calCurveParComb[1];
    calCurveParUnc_[0] = calCurveParUncComb[0];
    calCurveParUnc_[1] = calCurveParUncComb[1];
  }
  else if (useLep)
  {
    if (! do95) inputFileDir += "180110_lep/";
    calCurvePar_[0] = calCurveParLep[0];
    calCurvePar_[1] = calCurveParLep[1];
    calCurveParUnc_[0] = calCurveParUncLep[0];
    calCurveParUnc_[1] = calCurveParUncLep[1];
  }
  else if (useHad)
  {
    if (! do95) inputFileDir += "180110_had/";
    calCurvePar_[0] = calCurveParHad[0];
    calCurvePar_[1] = calCurveParHad[1];
    calCurveParUnc_[0] = calCurveParUncHad[0];
    calCurveParUnc_[1] = calCurveParUncHad[1];
  }
  
  
  /// Get nominal value
  fileName = inputFileDir+prefix+"result_minimum_nominal_widthx1_mass172p5.txt";
  if (! fexists(fileName.c_str()))
  {
    cerr << "File " << fileName << " does not exist..." << endl;
    exit(1);
  }
  
  fileIn.open(fileName.c_str());
  cout << "Opening " << fileName << "..." << endl;
  while ( getline(fileIn, line) )
  {
    if ( ! line.empty() )
    {
      istringstream iss(line);
      iss >> thisSyst >> thisWidth >> thisMass >> thisWidthValue >> thisWidthUnc;
      if (is2D) iss >> thisMassValue >> thisMassUnc;
    }
  }

  fileIn.close();
  
  /// Nominal correction calibration curve
  nomWidthValue = thisCorrWidthValue = ApplyCalibrationCurve(thisWidthValue, thisWidthUnc);
  cout << "Nominal   : " << nomWidthValue.first << " +- " << nomWidthValue.second << endl;
  
  /// Vary parameters calibration curve
  thisSyst = "ccConstUp";
  thisCorrWidthValue = ApplyCalibrationCurve(thisWidthValue, thisWidthUnc, "const", 1);
  cout << "Const up  : " << thisCorrWidthValue.first << " +- " << thisCorrWidthValue.second << endl;
  diff = thisCorrWidthValue.first - nomWidthValue.first;
  uncDiff = TMath::Sqrt(thisCorrWidthValue.second*thisCorrWidthValue.second + nomWidthValue.second*nomWidthValue.second);
  cout << "  ===>  Difference = " << diff << " +- " << uncDiff << endl;
  
  fileName = inputFileDir+prefix+"result_minimum_"+thisSyst+"_widthx1_mass172p5.txt";
  fileOut.open(fileName.c_str());
  fileOut << std::setw(18) << std::left << thisSyst << "  " << std::setw(5) << std::left << std::setprecision(5) << thisWidth << "   " << std::setw(5) << std::left << std::setprecision(5) << thisMass << "   " << std::setw(20) << std::right << std::setprecision(20) << thisCorrWidthValue.first << "  " << std::setw(20) << std::right << std::setprecision(20) << thisCorrWidthValue.second << std::endl;
  fileOut.close();
  
  thisSyst = "ccConstDown";
  thisCorrWidthValue = ApplyCalibrationCurve(thisWidthValue, thisWidthUnc, "const", -1);
  cout << "Const down: " << thisCorrWidthValue.first << " +- " << thisCorrWidthValue.second << endl;
  diff = thisCorrWidthValue.first - nomWidthValue.first;
  uncDiff = TMath::Sqrt(thisCorrWidthValue.second*thisCorrWidthValue.second + nomWidthValue.second*nomWidthValue.second);
  cout << "  ===>  Difference = " << diff << " +- " << uncDiff << endl;
  
  fileName = inputFileDir+prefix+"result_minimum_"+thisSyst+"_widthx1_mass172p5.txt";
  fileOut.open(fileName.c_str());
  fileOut << std::setw(18) << std::left << thisSyst << "  " << std::setw(5) << std::left << std::setprecision(5) << thisWidth << "   " << std::setw(5) << std::left << std::setprecision(5) << thisMass << "   " << std::setw(20) << std::right << std::setprecision(20) << thisCorrWidthValue.first << "  " << std::setw(20) << std::right << std::setprecision(20) << thisCorrWidthValue.second << std::endl;
  fileOut.close();
  
  
  thisSyst = "ccSlopeUp";
  thisCorrWidthValue = ApplyCalibrationCurve(thisWidthValue, thisWidthUnc, "slope", 1);
  cout << "Slope up  : " << thisCorrWidthValue.first << " +- " << thisCorrWidthValue.second << endl;
  diff = thisCorrWidthValue.first - nomWidthValue.first;
  uncDiff = TMath::Sqrt(thisCorrWidthValue.second*thisCorrWidthValue.second + nomWidthValue.second*nomWidthValue.second);
  cout << "  ===>  Difference = " << diff << " +- " << uncDiff << endl;
  
  fileName = inputFileDir+prefix+"result_minimum_"+thisSyst+"_widthx1_mass172p5.txt";
  fileOut.open(fileName.c_str());
  fileOut << std::setw(18) << std::left << thisSyst << "  " << std::setw(5) << std::left << std::setprecision(5) << thisWidth << "   " << std::setw(5) << std::left << std::setprecision(5) << thisMass << "   " << std::setw(20) << std::right << std::setprecision(20) << thisCorrWidthValue.first << "  " << std::setw(20) << std::right << std::setprecision(20) << thisCorrWidthValue.second << std::endl;
  fileOut.close();
  
  thisSyst = "ccSlopeDown";
  thisCorrWidthValue = ApplyCalibrationCurve(thisWidthValue, thisWidthUnc, "slope", -1);
  cout << "Slope down: " << thisCorrWidthValue.first << " +- " << thisCorrWidthValue.second << endl;
  diff = thisCorrWidthValue.first - nomWidthValue.first;
  uncDiff = TMath::Sqrt(thisCorrWidthValue.second*thisCorrWidthValue.second + nomWidthValue.second*nomWidthValue.second);
  cout << "  ===>  Difference = " << diff << " +- " << uncDiff << endl;
  
  fileName = inputFileDir+prefix+"result_minimum_"+thisSyst+"_widthx1_mass172p5.txt";
  fileOut.open(fileName.c_str());
  fileOut << std::setw(18) << std::left << thisSyst << "  " << std::setw(5) << std::left << std::setprecision(5) << thisWidth << "   " << std::setw(5) << std::left << std::setprecision(5) << thisMass << "   " << std::setw(20) << std::right << std::setprecision(20) << thisCorrWidthValue.first << "  " << std::setw(20) << std::right << std::setprecision(20) << thisCorrWidthValue.second << std::endl;
  fileOut.close();
  
  
  cout << "**** End of the program ****" << endl;
  cout << " - Goodbye" << endl;
  return 0;
}
