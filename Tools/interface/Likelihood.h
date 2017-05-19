#ifndef INTERFACE_LIKELIHOOD_H
#define INTERFACE_LIKELIHOOD_H

#include <stdio.h>
#include <TStyle.h>
#include <TSystem.h>
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
    Likelihood(double min, double max, std::string outputDirName, bool makeHistograms, bool calculateGoodEvtLL, bool verbose);
    ~Likelihood();
    void BookHistograms();
    void FillHistograms(double redMass, double massForWidthSF, bool isTTbar, bool isData, std::string catSuffix);
    void WriteHistograms(std::string histoFileName);
    void ConstructTGraphsFromHisto(std::string name);
    bool ConstructTGraphsFromFile();
    bool ConstructTGraphsFromFile(std::string name);
    void CalculateLikelihood(double redMass, bool isData);
    void CalculateLikelihood(double redMass, double massForWidthSF, double inputWidth, bool isTTbar, bool isData);
    void CalculateCPLikelihood(double redMass, double massForWidthSF, double inputWidth, bool isTTbar, bool isData);
    void CalculateGenLikelihood(double redMass, double massForWidthSF, double inputWidth, bool isTTbar, bool isData);
    void GetOutputWidth(double inputWidth);
    void GetOutputWidth(double inputWidth, std::string type);
    void GetOutputWidth(std::string inputFileName, double inputWidth);
    void PrintLikelihoodOutput(std::string llFileName);
    void PrintLikelihoodOutputData(std::string llFileName);
    void PrintMtmLikelihoodOutput(std::string llFileName);
      
  private:
    bool verbose_;
    std::string outputDirName_;
    std::string dirNameTGraphTxt_;
    std::string dirNameLLTxt_;
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
    
    bool calculateGoodEvtLL_;
    static double loglike_[];
    static double loglike_data_[];
    static double loglike_per_evt_[];
    static double loglike_good_evts_[];
    static double loglike_good_evts_data_[];
    static double loglike_CP_[];
    static double loglike_CP_per_evt_[];
    static double loglike_CP_good_evts_[];
    static double loglike_gen_[];
    static double loglike_gen_per_evt_[];
    static double loglike_gen_good_evts_[];
    
    bool calledLLCalculation_;
    bool calledCPLLCalculation_;
    bool calledGenLLCalculation_;
    
    std::ifstream fileIn_;
    std::vector<double> vecWidthFromFile_;
    std::vector<double> vecLLValsFromFile_;
    std::vector<double> vecGoodLLValsFromFile_;
    std::pair<double,double> output_;
    
    bool fexists(const char *filename);
    std::string ConvertDoubleToString(double Number);
    std::string DotReplace(double var);
    int LocMinArray(int n, double* array);
    void DrawGraph(TH1D* h, TGraph* g, std::string name);
    void DrawGraph(TGraph2D* g, std::string name);
    void DrawLikelihoods();
    void WriteOutput(int nPoints, double width, double *arrayCentre, double *arrayContent, std::string name, int dim);
    void CombineOutput();
    std::pair<double,double> CalculateOutputWidth(std::string inputFileName, std::string plotName);
    std::pair<double,double> CalculateOutputWidth(int nn, double* LLvalues, std::string plotName);
    std::pair<double,double> CalculateOutputWidth(int nn, double* evalWidths, double* LLvalues, std::string plotName);
    void DrawOutputLogLikelihood(TGraph* g, TF1* f, double maxX, double maxY, std::string name);
    void ReadLLValuesFromFile(std::string inputFileName);
    
};


#endif
