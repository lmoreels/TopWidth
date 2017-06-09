////////////////////////////////////////////////////////////////
/////                                                      /////
/////  Preliminary macro to test likelihood calculation    /////
/////                                                      /////
////////////////////////////////////////////////////////////////



#include "TStyle.h"
#include <ctime>
#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <map>
#include <array>
#include <vector>


// user defined
#include "Tools/interface/HelperTools.h"
#include "Tools/interface/Likelihood.h"


using namespace std;
// using namespace reweight;
// using namespace TopTree;

bool test = false;
bool testFit = true;


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


int main (int argc, char *argv[])
{
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "               *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  string outputDir = "test/";
  mkdir(("OutputLikelihood/"+outputDir).c_str(),0777);
  
  string dirNames[] = {"33", "34", "35", "36", "36", "38", "39", "41", "42", "42", "45", "46", "47"};
  double inputWidths[] = {0.2, 0.4, 0.6, 0.8, 1., 2., 3., 4., 5., 6., 7., 8., 9.};
  int nWidths = sizeof(inputWidths)/sizeof(inputWidths[0]);
  
  HelperTools* tls = new HelperTools();
  Likelihood* like = new Likelihood(0.6, 1.4, "LikelihoodTemplates/170608_1926/", outputDir+dateString, false, false, true);  // makeTGraphs, calculateGoodEvtLL, verbose
  
  string pathFile = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/OutputLikelihood/170609_11";
  string llFileName = "";
  
  if (testFit)
  {
    for (int iWidth = 0; iWidth < nWidths; iWidth++)
    {
      llFileName = pathFile+dirNames[iWidth]+"/output_loglikelihood_widthx"+tls->DotReplace(inputWidths[iWidth])+".txt";
      like->GetOutputWidth(llFileName, "", inputWidths[iWidth]);
    }
  }
  
  
  
  
  double time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  if ( time >= 60 )
  {
    int mins = time/60;
    float secs = time - mins*60;
    
    if (mins >= 60 )
    {
      int hours = mins/60;
      mins = mins - hours*60;
      cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
    }
    else
      cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
  }
    
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  cout << " - Goodbye" << endl;

  return 0;
  
}
