//#include <stdio.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TLine.h>
#include <map>
#include <TArray.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TColor.h>


using namespace std;


bool runLocally = false;

void DrawTemplatesSameSFs()
{
  map<string,TGraph*> graph;
  
  string pathFiles = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/LikelihoodTemplates/";
  string outputPath = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/test/";
  
  if (runLocally)
  {
    pathFiles = "/Users/lmoreels/cernbox/TopWidth/TopTrees/tempPlots/testSFs/Templates/";
    outputPath = pathFiles;
  }
  
  string fileName[] = {"170719_1800/", "170719_1802/", "170719_1803/", "170719_1804/", "170719_2119/"};
  int sizeFiles = sizeof(fileName)/sizeof(fileName[0]);
  string fileType[] = {"noSFs", "allSFs", "noBTagSFs", "noPUSFs", "allSFs_ttbarOnly"};
  
  string cats[] = {"CM", "WM", "NM"};
  int sizeCats = sizeof(cats)/sizeof(cats[0]);
  
  string widths[] = {"0p5", "0p55", "0p6", "0p65", "0p7", "0p75", "0p8", "0p85", "0p9"};
  int sizeWidths = sizeof(widths)/sizeof(widths[0]);
  
  Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+2, kCyan+1, kBlue+2, kMagenta, kViolet-5, kPink+10};
  
  string plotName = "", graphName = "";
  
  TFile *fileOut = new TFile((outputPath+"Templates_stacked.root").c_str(),"RECREATE");
  fileOut->cd();
  
  for (int iFile = 0; iFile < sizeFiles; iFile++)
  {
    TFile *fileIn = new TFile((pathFiles+fileName[iFile]+"TGraphFunctions.root").c_str(),"read");
    
    for (int iCat = 0; iCat < sizeCats; iCat++)
    {
      fileOut->cd();
      
      plotName = cats[iCat]+"_"+fileType[iFile];
      TCanvas *c1 = new TCanvas(plotName.c_str(),plotName.c_str());
      c1->cd();
      
      TLegend *leg = new TLegend(0.65,0.56,0.9,0.9);
      
      for (int iWidth = 0; iWidth < sizeWidths; iWidth++)
      {
        graphName = "g"+cats[iCat]+"_widthx"+widths[iWidth];
        graph[graphName] = (TGraph*) fileIn->Get(graphName.c_str());
        
        graph[graphName]->SetLineColor(colours[iWidth]);
        if ( iWidth == 0 )
        {
          graph[graphName]->SetTitle(plotName.c_str());
          graph[graphName]->Draw();
        }
        else graph[graphName]->Draw("same");
        leg->AddEntry(graph[graphName], ("width x "+widths[iWidth]).c_str(), "l");
      }
      
      leg->Draw();
      c1->Update();
      fileOut->cd();
      c1->Write();
      c1->SaveAs((outputPath+"Templates_"+plotName+".png").c_str());
      
    }
    
    fileIn->Close();
  }
  
  fileOut->Close();
  
}
