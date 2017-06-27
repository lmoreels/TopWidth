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
#include <TLine.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TFile.h>

// user defined
#include "EventReweighting.h"
#include "HelperTools.h"


class Likelihood{
  public:
    Likelihood(double min, double max, std::string outputDirName, std::string date, bool makeHistograms, bool calculateGoodEvtLL, bool verbose);
    ~Likelihood();
    /// Make TGraphs
    void BookHistograms();
    void FillHistograms(double redMass, double lumiWeight, double hadTopMassForWidthSF, double lepTopMassForWidthSF, bool isTTbar, bool isData, std::string catSuffix);
    void WriteHistograms(std::string histoFileName);
    void ConstructTGraphsFromHisto(std::string name);
    /// Get TGraphs (in order to calculate likelihood)
    bool ConstructTGraphsFromFile();
    bool ConstructTGraphsFromFile(std::string name);
    /// Calculate likelihood
    void CalculateLikelihood(double redMass, double lumiWeight, bool isData);
    void CalculateLikelihood(double redMass, double lumiWeight, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData);
    void CalculateCMLikelihood(double redMass, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData);
    void CalculateTempLikelihood(double redMass, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData);
    void CalculateGenLikelihood(double redMass, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData);
    /// Get output width
    void GetOutputWidth(double inputWidth);
    void GetOutputWidth(double inputWidth, bool writeToFile = false);
    void GetOutputWidth(double inputWidth, std::string type, bool writeToFile = false);
    void GetOutputWidth(std::string inputFileName, double inputWidth, bool writeToFile = false);
    void GetOutputWidth(std::string inputFileName, std::string inputDir, double inputWidth, bool writeToFile = false);
    /// Use pseudo experiments
    int InitPull(int nPsExp);
    void AddPsExp(int thisPsExp, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData);
    void CalculatePull(double inputWidth);
    /// Calibration curve
    std::pair<double,double> ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma);
    /// Write output from likelihood calculation to file
    void PrintLikelihoodOutput(std::string llFileName);
    void PrintLikelihoodOutputData(std::string llFileName);
    void PrintMtmLikelihoodOutput(std::string llFileName);
      
  private:
    bool verbose_;
    std::string outputDirName_;
    std::string dirNameTGraphTxt_;
    std::string dirNameLLTxt_;
    std::string dirNamePull_;
    std::string inputFileName_;
    std::string suffix_;
    std::string histoName_;
    double minRedMass_;
    double maxRedMass_;
    
    HelperTools *tls_;
    EventReweighting *rew_;
    
    static const double widthArray_[];
    static const std::string listCats_[];
    
    static const int nWidths_;
    static const int nCats_;
    
    static std::string stringWidthArray_[];
    
    static double loglike_[];
    static double loglike_data_[];
    static double loglike_per_evt_[];
    static double loglike_good_evts_[];
    static double loglike_good_evts_data_[];
    static double loglike_CM_[];
    static double loglike_CM_per_evt_[];
    static double loglike_CM_good_evts_[];
    static double loglike_temp_[];
    static double loglike_temp_per_evt_[];
    static double loglike_temp_good_evts_[];
    static double loglike_gen_[];
    static double loglike_gen_per_evt_[];
    static double loglike_gen_good_evts_[];
    
    static double loglike_pull_[][1000];
    static double loglike_pull_single_[];
    
    std::string thisWidth_;
    double thisWidthSF_;
    TFile *file_;
    TFile *fileTGraphs_;
    TFile *filePlots_;
    TGraph2D *gLL2D_;
    std::ofstream txtOutput_;
    std::ofstream txtOutputLL_;
    std::ofstream txtOutputPsExp_;
    
    std::map<std::string,TH1D*> histo_;
    std::map<std::string,TH1D*> histoSm_;
    std::map<std::string,TH1D*> histoTotal_;
    std::map<std::string,TGraph*> graph_;
    std::map<std::string,std::vector<double>> vecBinCentres_;
    std::map<std::string,std::vector<double>> vecBinContents_;
    
    bool calculateGoodEvtLL_;
    bool calledLLCalculation_;
    bool calledCMLLCalculation_;
    bool calledTempLLCalculation_;
    bool calledGenLLCalculation_;
    
    std::ifstream fileIn_;
    std::vector<double> vecWidthFromFile_;
    std::vector<double> vecLLValsFromFile_;
    std::vector<double> vecGoodLLValsFromFile_;
    std::pair<double,double> output_;
    
    static const double calCurvePar_[2];
    static const double calCurveParUnc_[2];
    
    int nPsExp_;
    
    int LocMinArray(int n, double* array);
    /// Get output width for pseudo experiments
    std::pair<double,double> GetOutputWidth(double inputWidth, int thisPsExp);
    /// Calculate output width (used in getters)
    std::pair<double,double> CalculateOutputWidth(std::string inputFileName, std::string inputDir, std::string plotName, bool writeToFile = false);
    std::pair<double,double> CalculateOutputWidth(int nn, double* LLvalues, std::string plotName, bool writeToFile = false);
    std::pair<double,double> CalculateOutputWidth(int nn, double* evalWidths, double* LLvalues, std::string plotName);
    std::pair<double,double> CalculateOutputWidth(int nn, double* evalWidths, double* LLvalues, std::string plotName, bool writeToFile = false);
    /// Draw functions
    void DrawGraph(TH1D* h, TGraph* g, std::string name);
    void DrawGraph(TGraph2D* g, std::string name);
    void DrawLikelihoods();
    void DrawOutputLogLikelihood(TGraph* g, TF1* f, double minX, double maxX, double maxY, std::string name, bool writeToFile);
    /// Read file
    void ReadLLValuesFromFile(std::string inputFileName, std::string inputDir);
    /// Write file
    void WritePsExpOutput(std::pair<double,double> *outputWidth, std::pair<double,double> *inputWidth, double genWidth);
    void WriteOutput(int nPoints, double width, double *arrayCentre, double *arrayContent, std::string name, int dim);
    void CombineOutput();
};


#endif
