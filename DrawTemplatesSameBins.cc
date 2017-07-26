//#include <stdio.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphSmooth.h>
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
bool combineTemplatesByReading = true;
bool combineTemplatesByAdding = false;


ifstream fileIn;
map<string,TGraph*> graph, graphSm;
map<string,vector<double>> vecBinCentres;
map<string,vector<double>> vecBinContents;

string pathFiles = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/LikelihoodTemplates/";
string outputPath = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/test/";

string fileName[] = {/*"170721_0937/", "170721_1018/", */"170721_1212/", "170720_1321/"/*, "170720_1838/", "170720_2235/"*/};
const int sizeFiles = sizeof(fileName)/sizeof(fileName[0]);
string fileType[] = {/*"20b", "45b", */"60b", "90b"/*, "450b", "900b"*/};

string cats[] = {"CM", "WM", "NM"};
const int sizeCats = sizeof(cats)/sizeof(cats[0]);
double fracs[sizeCats] = {0.479916, 0.148529, 0.371554};  // FOR 0.6 -- 1.4 !!!

string widths[] = {"0p2", "0p3", "0p4", "0p5", "0p75", "1", "2"};
const int sizeWidths = sizeof(widths)/sizeof(widths[0]);

//Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+2, kCyan+1, kBlue+2, kMagenta, kViolet-5, kPink+10};
Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen+2, kCyan+1, kBlue+2, kMagenta};


void DrawGraphs(map<string,TGraph*> graph, string type, string plotName)
{
  string graphName = "";
  string name = plotName+"_"+type;
  TCanvas *c1 = new TCanvas(name.c_str(),name.c_str());
  c1->cd();
  
  TLegend *leg = new TLegend(0.70,0.56,0.9,0.9);
  for (int iWidth = 0; iWidth < sizeWidths; iWidth++)
  {
    graphName = type+"_"+plotName+"_widthx"+widths[iWidth];
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
  c1->Write();
  c1->SaveAs((outputPath+"Templates_"+name+".png").c_str());
  
  delete leg;
  delete c1;
}

void DrawTemplatesSameBins()
{
  if (runLocally)
  {
    pathFiles = "/Users/lmoreels/cernbox/TopWidth/TopTrees/tempPlots/testSFs/Templates/";
    outputPath = pathFiles;
  }
  
  string plotName = "", graphName = "", graphNameSm = "", graphNameOrig = "";
  
  TFile *fileOut = new TFile((outputPath+"Templates_stacked.root").c_str(),"RECREATE");
  fileOut->cd();
  
  for (int iFile = 0; iFile < sizeFiles; iFile++)
  {
    TFile *fileIn = new TFile((pathFiles+fileName[iFile]+"TGraphFunctions.root").c_str(),"read");
    
    for (int iCat = 0; iCat < sizeCats; iCat++)
    {
      fileOut->cd();
      
      plotName = cats[iCat]+"_"+fileType[iFile];
      
      for (int iWidth = 0; iWidth < sizeWidths; iWidth++)
      {
        graphNameOrig = "g"+cats[iCat]+"_widthx"+widths[iWidth];
        graphName = fileType[iFile]+"_"+cats[iCat]+"_widthx"+widths[iWidth];
        graphNameSm = fileType[iFile]+"_"+cats[iCat]+"_Sm_widthx"+widths[iWidth];
        graph[graphName] = (TGraph*) fileIn->Get(graphNameOrig.c_str());
        TGraphSmooth* gs = new TGraphSmooth((cats[iCat]+"_widthx"+widths[iWidth]).c_str());
        graphSm[graphNameSm] = gs->SmoothSuper(graph[graphName],"",0,0);
      }
      
      DrawGraphs(graph, fileType[iFile], cats[iCat]);
      DrawGraphs(graphSm, fileType[iFile], cats[iCat]+"_Sm");
      
    }
    
    if (combineTemplatesByAdding)
    {
      /// To do
      //  Use graph->Eval
    }
    
    fileIn->Close();
  }
  
  if (combineTemplatesByReading)
  {
    double thisCentre, thisContent;
    string pathTxt = "";
    string line;
    
    for (int iFile = 0; iFile < sizeFiles; iFile++)
    {
      fileOut->cd();
      
      plotName = "TotalProb_"+fileType[iFile];
      
      for (int iWidth = 0; iWidth < sizeWidths; iWidth++)
      {
        for (int iCat = 0; iCat < sizeCats; iCat++)
        {
          pathTxt = pathFiles+fileName[iFile]+"OutputTxt/output_func_"+cats[iCat]+"_widthx"+widths[iWidth]+".txt";
          graphName = fileType[iFile]+"_"+cats[iCat]+"_widthx"+widths[iWidth];
          //graph[graphName] = new TGraph(pathTxt.c_str());
          (vecBinCentres[graphName]).clear();
          (vecBinContents[graphName]).clear();
          
          fileIn.open(pathTxt.c_str());
          cout << "Opening " << pathTxt << "..." << endl;
          
          while( getline(fileIn, line) )
          {
            istringstream iss(line);
            iss >> thisCentre >> thisContent;
            (vecBinCentres[graphName]).push_back(thisCentre);
            (vecBinContents[graphName]).push_back(thisContent);
          }
          fileIn.close();
          
        }  // end cats
        
        const int nPoints = (vecBinCentres[fileType[iFile]+"_"+cats[0]+"_widthx"+widths[iWidth]]).size();
        double binCentres[nPoints], binContents[sizeCats][nPoints], totalBinContent[nPoints];
        for (int i = 0; i < nPoints; i++)
        {
          binCentres[i] = (vecBinCentres[fileType[iFile]+"_"+cats[0]+"_widthx"+widths[iWidth]]).at(i);
          totalBinContent[i] = 0.;
          for (int iCat = 0; iCat < sizeCats; iCat++)
          {
            binContents[iCat][i] = (vecBinContents[fileType[iFile]+"_"+cats[iCat]+"_widthx"+widths[iWidth]]).at(i);
            totalBinContent[i] += fracs[iCat] * binContents[iCat][i];
          }
        }
        graphName = fileType[iFile]+"_TotalProb_"+"widthx"+widths[iWidth];
        graphNameSm = fileType[iFile]+"_TotalProb_Sm_"+"widthx"+widths[iWidth];
        graph[graphName] = new TGraph(nPoints, binCentres, totalBinContent);
        TGraphSmooth* gs2 = new TGraphSmooth(("TotalProb_widthx"+widths[iWidth]).c_str());
        graphSm[graphNameSm] = gs2->SmoothSuper(graph[graphName],"",0,0);
        
      }  // end widths
      
      fileOut->cd();
      DrawGraphs(graph, fileType[iFile], "TotalProb");
      DrawGraphs(graphSm, fileType[iFile], "TotalProb_Sm");
      
    }  // end files (/bins)
    
    
  }  // end combineTemplates
  
  
  fileOut->Close();
  
}
