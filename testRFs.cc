////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////



#include <TStyle.h>
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
#include <TLatex.h>
#include <TF1.h>
#include <TF2.h>


// user defined
#include "Tools/interface/ResolutionFunctions.h"


using namespace std;
// using namespace reweight;
// using namespace TopTree;

bool test = false;
bool testFit = false;
bool testRead = true;

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

int DrawFunction(TF2* f, string varName, string saveName)
{
  if ( f == NULL ) return 1;
  
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
  
  return 0;
}

int DrawFunction(TF2* f, string objectName, string varName, string region, string saveName)
{
  if ( f == NULL ) return 1;
  
  string axisTitle = "";
  string varLabel = "";
  string unitLabel = "";
  if ( varName.find("Theta") != std::string::npos || varName.find("theta") != std::string::npos )
    varLabel = "#theta";
  else if ( varName.find("Phi") != std::string::npos || varName.find("phi") != std::string::npos )
    varLabel = "#phi";
  else if ( varName.find("Et") != std::string::npos || varName.find("ET") != std::string::npos || varName.find("et") != std::string::npos )
  {
    varLabel = "E_{T}";
    unitLabel = "GeV";
  }
  
  string objectLabel = "";
  if ( objectName.find("non-b jet") != std::string::npos || objectName.find("nonb jet") != std::string::npos || objectName.find("nonbjet") != std::string::npos || objectName.find("non b jet") != std::string::npos )
    objectLabel = "light jet";
  else if ( objectName.find("b jet") != std::string::npos || objectName.find("bjet") != std::string::npos || objectName.find("B jet") != std::string::npos )
    objectLabel = "b jet";
  
  TLatex txt;
  txt.SetTextColor(1);
  txt.SetTextSize(0.03);
  
  TCanvas *c1 = new TCanvas("c1",("Resolution function for "+varName+" ("+objectName+") "+region).c_str(),10,200,900,700);
  
  f->SetContour(48);
  f->SetFillColor(45);
  
  // Draw this function with color mesh option
  f->SetLineWidth(1);
  f->SetLineColor(5);
  f->Draw("surf1");
  
  f->SetTitle("");
  axisTitle = varLabel+" ("+objectLabel+")";
  if (! unitLabel.empty() ) axisTitle += " ("+unitLabel+")";
  f->GetHistogram()->GetXaxis()->SetTitle(axisTitle.c_str());
  axisTitle = varLabel+" (parton) - "+varLabel+" ("+objectLabel+")";
  if (! unitLabel.empty() ) axisTitle += " ("+unitLabel+")";
  f->GetHistogram()->GetYaxis()->SetTitle(axisTitle.c_str());
  f->GetHistogram()->GetXaxis()->SetTitleOffset(1.8);
  f->GetHistogram()->GetYaxis()->SetTitleOffset(2.0);
  if ( unitLabel.empty() ) f->GetHistogram()->GetYaxis()->SetTitleOffset(2.2);
  
  if ( region.find("B") != std::string::npos || region.find("Barrel") != std::string::npos || region.find("barrel") != std::string::npos )
    txt.DrawLatexNDC(0.75, 0.82, "#font[132]{|#eta| < 1.3}");
  else if ( region.find("O") != std::string::npos || region.find("Overlap") != std::string::npos || region.find("overlap") != std::string::npos )
    txt.DrawLatexNDC(0.75, 0.82, "#font[132]{1.3 < |#eta| < 1.5}");
  else if ( region.find("E") != std::string::npos || region.find("Endcap") != std::string::npos || region.find("endcap") != std::string::npos )
    txt.DrawLatexNDC(0.75, 0.82, "#font[132]{|#eta| > 1.5}");
  
  
  c1->Update();
  c1->SaveAs((saveName+".png").c_str());
  
  delete c1;
  
  return 0;
}

int DrawFunction(TF1* f, string varName, string saveName)
{
  if ( f == NULL ) return 1;
  
  TCanvas *c1 = new TCanvas("c1",("Projection of resolution function for "+varName).c_str(),10,200,900,700);
  
  if ( saveName.std::string::find("_x") == std::string::npos )
    f->SetTitle(("Resolution function for "+varName).c_str());
  f->Draw();
  
  f->GetHistogram()->GetXaxis()->SetTitle(varName.c_str());
  f->GetHistogram()->GetXaxis()->SetTitleOffset(1.1);
  
  c1->Update();
  c1->SaveAs((saveName+".png").c_str());
  
  delete c1;
  
  return 0;
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
  
  string inputFileName = "input/PlotsForResolutionFunctions_171025.root";
  
  string rfFileName = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/PlotsForResolutionFunctions_testFit_171025.root";
  //TFile *foutRF = new TFile(rfFileName.c_str(), "RECREATE");
  //foutRF->cd();
  //TDirectory* rootDir = foutRF->mkdir(dateString.c_str());
  //rootDir->cd();
  
  if (testFit)
  {
    rf->makeFit(inputFileName, rfFileName, false);
  }
  
  if (testRead)
  {
    rfFileName = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/input/PlotsForResolutionFunctions_testFit_171025.root";
    
    TF2* f2 = (TF2*) rf->getFitFunction2D(rfFileName, "Et", "bjet", "B");
    //DrawFunction(f2, "b jet Et (barrel)", "RF_bjet_Et_B_fitF");
    DrawFunction(f2, "b jet", "Et" , "B", "RF_bjet_Et_B_fitF_clean");
    
    TF2* f2nb = (TF2*) rf->getFitFunction2D(rfFileName, "Et", "nonbjet", "B");
    //DrawFunction(f2nb, "non-b jet Et (barrel)", "RF_nonbjet_Et_B_fitF");
    DrawFunction(f2nb, "non-b jet", "Et" , "B", "RF_nonbjet_Et_B_fitF_clean");
    
    TF2* f2nbo = (TF2*) rf->getFitFunction2D(rfFileName, "Et", "nonbjet", "O");
    DrawFunction(f2nbo, "non-b jet", "Et" , "O", "RF_nonbjet_Et_O_fitF_clean");
    
    TF2* f2nbe = (TF2*) rf->getFitFunction2D(rfFileName, "Et", "nonbjet", "E");
    DrawFunction(f2nbe, "non-b jet", "Et" , "E", "RF_nonbjet_Et_E_fitF_clean");
    
    TF1 *f1 = (TF1*) rf->getFitFunction1D(rfFileName, "Et", "bjet", "B");
    DrawFunction(f1, "b jet Et (barrel)", "RF_bjet_Et_B_fitF_x");
    
    TF2* f2th = (TF2*) rf->getFitFunction2D(rfFileName, "theta", "bjet", "B");
    //DrawFunction(f2th, "b jet theta (barrel)", "RF_bjet_theta_B_fitF");
    DrawFunction(f2th, "b jet", "theta" , "B", "RF_bjet_theta_B_fitF_clean");
    
    TF2* f2thnb = (TF2*) rf->getFitFunction2D(rfFileName, "theta", "nonbjet", "B");
    //DrawFunction(f2thnb, "non-b jet theta (barrel)", "RF_nonbjet_theta_B_fitF");
    DrawFunction(f2thnb, "non-b jet", "theta" , "B", "RF_nonbjet_theta_B_fitF_clean");
    
    TF2* f2ph = (TF2*) rf->getFitFunction2D(rfFileName, "phi", "bjet", "B");
    //DrawFunction(f2ph, "b jet phi (barrel)", "RF_bjet_phi_B_fitF");
    DrawFunction(f2ph, "b jet", "phi" , "B", "RF_bjet_phi_B_fitF_clean");
    
    TF2* f2phnb = (TF2*) rf->getFitFunction2D(rfFileName, "phi", "nonbjet", "B");
    //DrawFunction(f2phnb, "non-b jet phi (barrel)", "RF_nonbjet_phi_B_fitF");
    DrawFunction(f2phnb, "non-b jet", "phi" , "B", "RF_nonbjet_phi_B_fitF_clean");
    
    TF1 *f = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "bjet", "B");
    DrawFunction(f, "b jet Et (barrel)", "RF_bjet_Et_B");
    
    TF1 *fnb = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "B");
    DrawFunction(fnb, "non-b jet Et (barrel)", "RF_nonbjet_Et_B");
    
    TF1 *fth = (TF1*) rf->getResolutionFunction(rfFileName, "theta", "bjet", "B");
    DrawFunction(fth, "b jet theta (barrel)", "RF_bjet_theta_B");
    
    TF1 *fthnb = (TF1*) rf->getResolutionFunction(rfFileName, "theta", "bjet", "B");
    DrawFunction(fthnb, "non-b jet theta (barrel)", "RF_nonbjet_theta_B");
    
    TF1 *fph = (TF1*) rf->getResolutionFunction(rfFileName, "phi", "bjet", "B");
    DrawFunction(fph, "b jet phi (barrel)", "RF_bjet_phi_B");
    
    TF1 *fphnb = (TF1*) rf->getResolutionFunction(rfFileName, "phi", "nonbjet", "B");
    DrawFunction(fphnb, "non-b jet phi (barrel)", "RF_nonbjet_phi_B");
    
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
