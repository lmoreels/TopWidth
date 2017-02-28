////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////



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
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TF1.h>
#include <TF2.h>


// user defined
#include "Tools/interface/ResolutionFunctions.h"


using namespace std;
// using namespace reweight;
// using namespace TopTree;

bool test = true;
bool testFit = true;
bool testRead = false;

std::map<std::string, TF1*> errorFuncMap;

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

void DrawFunction(TF2* f, string varName, string saveName)
{
  TCanvas *c1 = new TCanvas("c1",("Resolution function for "+varName).c_str(),10,200,900,700);
  
  f->SetContour(48);
  f->SetFillColor(45);
  
  // Draw this function with color mesh option
  f->SetLineWidth(1);
  f->SetLineColor(5);
  f->Draw("surf1");
  
  f->GetHistogram()->GetXaxis()->SetTitle(varName.c_str());
  f->GetHistogram()->GetYaxis()->SetTitle((varName+" difference").c_str());
  f->GetHistogram()->GetXaxis()->SetTitleOffset(1.6);
  f->GetHistogram()->GetYaxis()->SetTitleOffset(1.6);
  
  c1->Update();
  c1->SaveAs((saveName+".png").c_str());
  
  delete c1;
}

void DrawFunction(TF1* f, string varName, string saveName)
{
  TCanvas *c1 = new TCanvas("c1",("Projection of resolution function for "+varName).c_str(),10,200,900,700);
  
  f->Draw();
  
  f->GetHistogram()->GetXaxis()->SetTitle(varName.c_str());
  f->GetHistogram()->GetXaxis()->SetTitleOffset(1.1);
  
  c1->Update();
  c1->SaveAs((saveName+".png").c_str());
  
  delete c1;
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
  
  
  
  ResolutionFunctions* rf = new ResolutionFunctions(false, true);
  
  string inputFileName = "PlotsForResolutionFunctions.root";
  
  string rfFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/PlotsForResolutionFunctions_testFit.root";
  //TFile *foutRF = new TFile(rfFileName.c_str(), "RECREATE");
  //foutRF->cd();
  //TDirectory* rootDir = foutRF->mkdir(dateString.c_str());
  //rootDir->cd();
  
  if (testFit)
  {
    rf->makeFit(inputFileName, rfFileName);
  }
  
  if (testRead)
  { 
    TF2* f2 = (TF2*) rf->getFitFunction2D(rfFileName, "Et", "bjet", "B");
    DrawFunction(f2, "bjet Et (barrel)", "RF_bjet_Et_B_fitF");
    
    TF1 *f1 = (TF1*) rf->getFitFunction1D(rfFileName, "Et", "bjet", "B");
    DrawFunction(f1, "bjet Et (barrel)", "RF_bjet_Et_B_fitF_x");
    
    TF1 *f = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "B");
    //DrawFunction(f, "bjet Et (barrel)", "RF_bjet_Et_B");
    DrawFunction(f, "non-bjet Et (barrel)", "RF_nonbjet_Et_B");
    
    if (test)
    {
      cout << "Calculate resolution from function: " << f->Eval(100.) << endl;
      cout << "Resolution of   (5) = " << rf->getResolution(rfFileName, "Et", "nonbjet", 5., "B") << "     (100) = " << rf->getResolution(rfFileName, "Et", "nonbjet", 100., "B") << endl;
      
      errorFuncMap["nonbjetEt_B"] = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "B");
      DrawFunction(errorFuncMap["nonbjetEt_B"], "non-bjet Et (barrel)", "RF_nonbjet_Et_B2");
      cout << errorFuncMap["nonbjetEt_B"]->Eval(100.) << endl;
    }
    
    delete f1;
    delete f2;
    delete f;
  }
  
  //foutRF->Close();
  
  //rf->writeTable(rfFileName);
  
  //delete foutRF;
  
  
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
