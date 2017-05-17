#ifndef INTERFACE_LIKELIHOOD_H
#define INTERFACE_LIKELIHOOD_H

#include <stdio.h>
#include <TStyle.h>
#include <sys/stat.h>
#include <errno.h>
#include <cmath>
#include <TMath.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TFile.h>

// user defined
#include "EventReweighting.h"


class Likelihood{
  public:
    Likelihood(double min, double max, std::string outputDirName, std::string date, bool verbose);
    ~Likelihood();
    void BookHistograms();
    void FillHistograms(double redMass, double massForWidthSF, bool isTTbar, std::string catSuffix);
    void WriteHistograms(std::string histoFileName);
    void ConstructTGraphsFromHisto(std::string name);
    void ConstructTGraphsFromFile(std::string name);
    int CalculateOutputWidth(double *evalWidths, double *LLvalues, std::string inputType);
      
  private:
    bool verbose_;
    std::string date_;
    std::string outputDirName_;
    std::string inputFileName_;
    std::string suffix_;
    std::string histoName_;
    double minRedMass_;
    double maxRedMass_;
    
    static const double widthArray_[];
    static const std::string listCats_[];
    
    static const int nWidths_;
    static const int nCats_;
    std::vector<std::string> vecWidthStr_;
    
    std::string thisWidth_;
    double thisWidthSF_;
    TFile *file_;
    TFile *fileTGraphs_;
    TGraph2D *gLL2D_;
    std::ofstream txtOutput_;
    EventReweighting *rew;
    
    std::map<std::string,TH1D*> histo_;
    std::map<std::string,TH1D*> histoSm_;
    std::map<std::string,TH1D*> histoTotal_;
    std::map<std::string,TGraph*> graph_;
    std::map<std::string,std::vector<double>> vecBinCentres_;
    std::map<std::string,std::vector<double>> vecBinContents_;
    
    bool fexists(const char *filename);
    std::string ConvertDoubleToString(double Number);
    std::string DotReplace(double var);
    void DrawGraph(TH1D* h, TGraph* g, std::string name);
    void DrawGraph(TGraph2D* g, std::string name);
    void DrawLikelihoods();
    void WriteOutput(int nPoints, double width, double *arrayCentre, double *arrayContent, std::string name, int dim);
    void CombineOutput();
    
};


#endif
