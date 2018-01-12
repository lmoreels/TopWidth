#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <math.h>
#include <TMath.h>
#include <algorithm>


using namespace std;


bool is2D = false;
string inputFileDir = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/Systematics/temp/180110_comb/";
//string inputFileDir = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/Systematics/temp/180110_lep/";
//string inputFileDir = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/Systematics/temp/180110_had/";

/// 1D width
// parameter = comb
const double calCurvePar_[2] = {-0.249918, 1.08351};   //{-0.2526, 1.084};
const double calCurveParUnc_[2] = {0.0351629, 0.0139145};   //{0.03321, 0.01345};
// parameter = redMlb
//const double calCurvePar_[2] = {-0.267184, 1.08957};   //{-0.245, 1.083};
//const double calCurveParUnc_[2] = {0.0550971, 0.0248913};   //{0.05185, 0.02447};
// parameter = redTopMass_old
//const double calCurvePar_[2] = {-0.152752, 1.09383};   //{-0.15551, 1.09669};        // 180107
//const double calCurveParUnc_[2] = {0.0436359, 0.0166395};   //{0.0469114, 0.0175116};  // 180107

ifstream fileIn;
ofstream fileOut;

string line, fileName;
string thisSyst;
double thisWidth, thisMass, thisWidthValue, thisWidthUnc, thisMassValue, thisMassUnc;
double nomVal, nomUnc;
vector<double> pdfWidthValues, pdfWidthUnc;
pair<double,double> thisCorrWidthValue, thisCorrMassValue;


bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile.good();
}

string ConvertIntToString(int Number, int pad = 0)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  if ( pad > 1 ) { convert << std::setw(pad) << std::setfill('0');}
  convert << Number;
  return convert.str();
}

std::pair<double,double> ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma)
{
  /// return 'thisInputWidth' and uncertainty
  //  thisOutputWidth = Par_[0] + Par_[1] * thisInputWidth
  
  double thisInputWidth = (thisOutputWidth - calCurvePar_[0])/calCurvePar_[1];
  double thisInputWidthSigma = TMath::Sqrt( thisOutputWidthSigma*thisOutputWidthSigma + calCurveParUnc_[0]*calCurveParUnc_[0] + calCurveParUnc_[1]*calCurveParUnc_[1]*thisInputWidth*thisInputWidth )/calCurvePar_[1];
  
  return std::pair<double,double>(thisInputWidth,thisInputWidthSigma);
}


int main()
{
  cout << "**** Beginning of the program ****" << endl;
  
  pdfWidthValues.clear();
  pdfWidthUnc.clear();
  
  /// Get results pdf replicas
  for (int ipdf = 0; ipdf < 100; ipdf++)
  {
    fileName = inputFileDir+"pdf/result_minimum_pdfVar"+ConvertIntToString(ipdf)+"_widthx1_mass172p5.txt";
    if (! fexists(fileName.c_str()))
    {
      cerr << "File " << fileName << " does not exist..." << endl;
    }
    
    fileIn.open(fileName.c_str());
    //cout << "Opening " << fileName << "..." << endl;
    while ( getline(fileIn, line) )
    {
      if ( ! line.empty() )
      {
        istringstream iss(line);
        iss >> thisSyst >> thisWidth >> thisMass >> thisWidthValue >> thisWidthUnc;
        if (is2D) iss >> thisMassValue >> thisMassUnc;
        
        thisCorrWidthValue = ApplyCalibrationCurve(thisWidthValue, thisWidthUnc);
        
        pdfWidthValues.push_back(thisCorrWidthValue.first);
        pdfWidthUnc.push_back(thisCorrWidthValue.second);
      }
    }
    
    fileIn.close();
  }
  
  int nPdfs = pdfWidthValues.size();
  if ( nPdfs != 100 )
    cerr << "WARNING: Did not find all 100 PDF replicas... Proceed with caution!" << endl;
  
  
  /// Get nominal value
  fileName = inputFileDir+"result_minimum_nominal_widthx1_mass172p5.txt";
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
      
      thisCorrWidthValue = ApplyCalibrationCurve(thisWidthValue, thisWidthUnc);
      
      nomVal = thisCorrWidthValue.first;
      nomUnc = thisCorrWidthValue.second;
    }
  }

  fileIn.close();
  
  
  /// Calculate RMS
  double tmp;
  double thisRMS = 0.;
  for (int ipdf = 0; ipdf < 100; ipdf++)
  {
    tmp = pdfWidthValues[ipdf] - nomVal;
    thisRMS += tmp*tmp;
  }
  thisRMS *= 1./(nPdfs - 1.);
  thisRMS = sqrt(thisRMS);
  
  double thisRMSUnc = 0.;   // not necessary
  
  cout << "Variation due to PDFs = " << thisRMS << endl;
  
  /// Write output
  thisSyst = "pdfVar";
  fileName = inputFileDir+"result_minimum_"+thisSyst+"_widthx1_mass172p5.txt";
  fileOut.open(fileName.c_str());
  fileOut << std::setw(18) << std::left << thisSyst << "  " << std::setw(5) << std::left << std::setprecision(5) << thisWidth << "   " << std::setw(5) << std::left << std::setprecision(5) << thisMass << "   " << std::setw(20) << std::right << std::setprecision(20) << thisRMS << "  " << std::setw(20) << std::left << std::setprecision(20) << thisRMSUnc << std::endl;
  fileOut.close();
    
  cout << "**** End of the program ****" << endl;
  cout << " - Goodbye" << endl;
  return 0;
}

