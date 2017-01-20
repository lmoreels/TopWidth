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
#include <TF2.h>


// user defined
#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFunctions.h"


using namespace std;
// using namespace reweight;
// using namespace TopTree;

bool test = false;


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

Double_t dblGaus(Double_t *x, Double_t *par)
{
  Double_t par0 = par[0] + x[0] * par[1];
  Double_t par1 = par[2] + x[0] * par[3];
  Double_t par2 = par[4] + x[0] * par[5];
  Double_t par3 = par[6] + x[0] * par[7];
  Double_t par4 = par[8] + x[0] * par[9];
  Double_t par5 = par[10] + x[0] * par[11];
  
  Double_t norm = 1./TMath::Sqrt(2.*TMath::Pi()) * par5/TMath::Sqrt( pow(par1, 2) + par2*pow(par4, 2) );
  Double_t narrowGaus = TMath::Exp( - TMath::Power( (x[1] - par0)/par1 , 2) /2. );
  Double_t broadGaus = TMath::Exp( - TMath::Power( (x[1] - par3)/par4 , 2) /2. );
  return norm * ( narrowGaus + par2 * broadGaus );
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
  
  
  
  ResolutionFunctions* rf = new ResolutionFunctions(false);
  
  string inputFileName = "PlotsForResolutionFunctions.root";
  
  string rfFileName = "/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/PlotsForResolutionFunctions_test.root";
  //TFile *foutRF = new TFile(rfFileName.c_str(), "RECREATE");
  //foutRF->cd();
  //TDirectory* rootDir = foutRF->mkdir(dateString.c_str());
  //rootDir->cd();
  
  //rf->makeFit(inputFileName, rfFileName);
  
  std::vector<std::array<double, 2> > myParams = rf->getParameters(rfFileName, "E", "bjet", true);
  if (test) cout << myParams[0][0] << "  " << myParams[0][1] << endl;
  if (test) cout << myParams[1][0] << "  " << myParams[1][1] << endl;
  
  //foutRF->Close();
  
  //rf->writeTable(rfFileName);
  
  //delete foutRF;
  
  /// Create function
  TF2 *f2 = new TF2("f2",dblGaus,0.,200.,-80.,80., 12);
  for (int iPar = 0; iPar < 12; iPar++)
  {
    int par = (int) ((double)iPar/2.);
    //if ( iPar%2 != 0 ) par += 1;
    f2->SetParameter(iPar, myParams[par][iPar%2]);
    if (test) cout << "Parameter " << iPar << " set to " << myParams[par][iPar%2] << endl;
  }
  //cout << f2->Eval(10.,5.) << endl;
  
  DrawFunction(f2, "bjet energy", "RF_bjet_E");


  delete f2;
  
  
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
