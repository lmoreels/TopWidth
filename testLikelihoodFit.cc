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
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TFitResult.h>
//#include <Math/Minimizer.h>
#include <Math/MinimizerOptions.h>
//#include <Math/Functor.h>
//#include <Math/GSLMinimizer.h>
//#include <Math/GSLSimAnMinimizer.h>


// user defined
#include "Tools/interface/HelperTools.h"
#include "Tools/interface/Likelihood.h"

using namespace std;


bool verbose_ = true;

Double_t func(Double_t *x, Double_t *p)
{
  return p[0] + p[1]*x[0] + p[2]*x[0]*x[0] + p[3]*x[1] + p[4]*x[1]*x[1] + p[5]*x[0]*x[1]*x[1] + p[6]*x[1]*x[1]*x[1]*x[1];
}

int main (int argc, char *argv[])
{
  HelperTools *tls_ = new HelperTools();
  
  string dateString = tls_->MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "               *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  string outputDir = "test/";
  mkdir(("OutputLikelihood/"+outputDir).c_str(),0777);
  outputDir += dateString;
  
  string inputDir = "LikelihoodTemplates/171011_1932/";
  Likelihood* like = new Likelihood(0.6, 1.4, inputDir, outputDir, false, false, true);  // makeTGraphs, verbose
  
  string pathFile = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/OutputLikelihood/171011_224/";
  string llFileName = "likelihoodValues.txt";
  
//   vector<double> vecWidths, vecMasses, vecLLValues;
//   double thisWidth, thisMass, thisLLValue;
//   string line;
//   
//   /// Open file
//   if (! tls_->fexists((pathFile+llFileName).c_str()) )
//   {
//     std::cerr << "WARNING: File " << pathFile+llFileName << " does not exist." << std::endl;
//     exit(1);
//   }
//   fileIn_.open((pathFile+llFileName).c_str());
//   if (verbose_) std::cout << "Opening " << pathFile+llFileName << "..." << std::endl;
//   
//   vecWidths.clear(); vecMasses.clear(); vecLLValues.clear();
//   
//   while( getline(fileIn_, line) )
//   {
//     std::istringstream iss(line);
//     if ( line.find("Width") != std::string::npos )
//     {
//       iss >> var;
//       while ( iss >> val )
//       {
//         vecWidthFromFile_.push_back(val);
//       }
//     }
//     if ( line.find("LL") != std::string::npos )
//     {
//       iss >> var;
//       while ( iss >> val )
//       {
//         vecLLValsFromFile_.push_back(val);
//       }
//     }
//   }
//   
//   fileIn_.close();
//   
//   std::pair<int,int> locMin = like->LocMinArray(nn, LLvalues);
//   std::cout << "Index of minimum LL value is (" << locMin.first << ", " << locMin.second << ")" << "           " << LLvalues[locMin.first][locMin.second] << std::endl;

  std::pair<double,double> locMin = std::pair<double,double>(0.95,172.5);
  double centreVal = locMin.first;
  
  double interval = 0.15;
  if ( centreVal <= interval ) interval = centreVal - 0.1;
  if ( centreVal > 2.8 ) interval = 0.4;
  double fitmax = centreVal + interval;
  double fitmin = centreVal - interval;
  
  double fitminY = locMin.second - 0.25;
  if ( fitminY < 169.5 ) fitminY = 169.5;
  double fitmaxY = locMin.second + 0.25;
  if ( fitmaxY > 175.5 ) fitmaxY = 175.5;
  
  if (verbose_)
    std::cout << "Look for minimum around (" << centreVal << ", " << locMin.second << ")" << std::endl << "Boundaries are [" << fitmin << ", " << fitmax << "] for width and [" << fitminY << ", " << fitmaxY << "] for mass" << std::endl;
  
  
  string plotName = "loglikelihood_vs_width_vs_mass";
  TFile *filePlots_ = new TFile(("OutputLikelihood/"+outputDir+"/File_"+plotName+".root").c_str(), "RECREATE");
  filePlots_->cd();
  
  TGraph2D *g = new TGraph2D((pathFile+llFileName).c_str());
  g->SetMaxIter(500000);
  
  //TF2 *parabola = new TF2("parabola", "[0]++[1]*x++[2]*x*x++[3]*y++[4]*y*y++[5]*x*x*y++[6]*x*y*y++[7]*x*x*x", fitmin, fitmax, fitminY, fitmaxY);
  TF2 *parabola = new TF2("parabola", "1260541.164*[0] + [1]*x + [2]*x*x + [3]*y + [4]*y*y + [5]*x*y*y + [6]*y*y*y*y", fitmin, fitmax, fitminY, fitmaxY);
  //TF2 *parabola = new TF2("parabola", "1260541.164 + [1]*(x - 172.5)*(x - 172.5) + [2]*y", fitmin, fitmax, fitminY, fitmaxY);
  //TF2 *parabola = new TF2("parabola", "[0] + [1]*x + [2]*x*x + [3]*y + [4]*y*y + [5]*x*y + [6]*x*x*y + [7]*x*y*y", fitmin, fitmax, fitminY, fitmaxY);
  //TF2 *parabola = new TF2("parabola", "x*x*[0] + x*[1] + y*y*[2] + y*[3] + [4] + x*y*[5] + x*x*y*[6] + x*y*y*[7] + x*x*x*[8] + y*y*y*[9]", fitmin, fitmax, fitminY, fitmaxY);
  //TF2 *parabola = new TF2("parabola", "[0]*pow(x-[1],2) + [2]*pow(y-[3], 2) + [4]", fitmin, fitmax, fitminY, fitmaxY);
  //parabola->SetParameters(100., centreVal, 200., massArray_[locMin.second], LLvalues[locMin.first][locMin.second], 0.02);
  //parabola->SetParameter(1, centreVal);
  //parabola->SetParameter(3, massArray_[locMin.second]);
  //parabola->SetParameter(4, LLvalues[locMin.first][locMin.second]);
//  parabola->SetParLimits(1, fitmin, fitmax);
  //parabola->SetParLimits(3, fitminY, fitmaxY);
  //parabola->SetParLimits(4, LLvalues[locMin.first][locMin.second] - 30., LLvalues[locMin.first][locMin.second] + 30.);
  //parabola->SetParameters(10229.4, 5.50558e+06,-40805.9,7.02645e+06,-4.02035e+08,-63944.6,-57.6472,185.66,-75.6738,78.9921);
  //parabola->Write();
  parabola->SetParameter(0, 1.);
  parabola->SetParameter(6, 0.01);
  
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiMin","BFGS2");
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit","Simplex");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-3);
  TFitResultPtr p = g->Fit(parabola,"SWRMEV");
  if (verbose_) p->Print("V");
  
  // Choose method upon creation between:
  // kConjugateFR, kConjugatePR, kVectorBFGS,
  // kVectorBFGS2, kSteepestDescent
//   ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS );
//   min.SetTolerance(1e-3);
//   ROOT::Math::Functor f(&g,2); 
//   double step[2] = {0.001,0.01};
//   double variable[2] = {locMin.first, locMin.second};
//   min.SetFunction(f);
//   // Set the free variables to be minimized!
//   min.SetVariable(0,"x",variable[0], step[0]);
//   min.SetVariable(1,"y",variable[1], step[1]);
//   min.Minimize(); 
//   
//   const double *xs = min.X();
//   cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << endl;
//        << g(xs) << endl;
  
  g->Write();
  TF2 *fitresult = (TF2*) g->FindObject("parabola");
  fitresult->Write();
  g->Write();
  
  
  
  /// Calculate output width & mass
  double outputWidth, outputMass;
  double minimum = fitresult->GetMinimumXY(outputWidth, outputMass);
  
  TF12 *parbx = new TF12("parbx", fitresult, outputMass, "x");
  TF12 *parby = new TF12("parby", fitresult, outputWidth, "y");
  
  double lowerSigma = parbx->GetX(minimum + 0.5, fitmin-interval, outputWidth);
  double upperSigma = parbx->GetX(minimum + 0.5, outputWidth, fitmax+interval);
  double sigma = (upperSigma - lowerSigma)/2.;
  if ( lowerSigma <= fitmin && upperSigma >= fitmax )
    std::cerr << "Likelihood::CalculateOutputWidth: ERROR: Uncertainty calculation limited by fit boundaries. Do not trust..." << std::endl;
  else if ( lowerSigma <= fitmin ) sigma = upperSigma - outputWidth;
  else if ( upperSigma >= fitmax ) sigma = outputWidth - lowerSigma;
  
  //std::cout << "Minimum -log(like) value is " << parabola->Eval(outputWidth) << std::endl;
  std::cout << "Minimum -log(like) value is " << minimum << std::endl;
  if (verbose_) std::cout << "Position of minimum is (" << outputWidth << ", " << outputMass << ")" << std::endl;
  if (verbose_) std::cout << "lowerSigma = " << lowerSigma << " and upperSigma = " << upperSigma << std::endl;
  
  std::cout << "For an input width of 1 the minimum can be found at " << outputWidth << " and the uncertainty is " << sigma << std::endl;
  
  like->DrawOutputLogLikelihood(g, fitresult, fitmin, fitmax, fitminY, fitmaxY, plotName, true);
  like->DrawOutputLogLikelihood(g, fitresult, outputWidth, outputMass, minimum, plotName, true);
  
  filePlots_->Close();
  cout << "File closed" << endl;
  delete filePlots_;
  
  delete parbx;
  delete parby;
  delete fitresult;
  delete parabola;
  delete g;
  
  //like->GetOutputWidth(llFileName, "", inputWidths[iWidth]);
  
  
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

