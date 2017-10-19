#ifndef INTERFACE_LIKELIHOOD2D_H
#define INTERFACE_LIKELIHOOD2D_H

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
#include <TProfile.h>
#include <TH1.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphSmooth.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
#include <TFitResult.h>
#include <Math/MinimizerOptions.h>
//#include <TPaveStats.h>
//#include <TList.h>
//#include <TText.h>

// user defined
#include "EventReweighting.h"
#include "HelperTools.h"


class Likelihood2D{
  public:
    Likelihood2D(double min, double max, std::string outputDirName, std::string date, bool useHadTopOnly, bool makeHistograms, bool verbose);
    ~Likelihood2D();
    void SetMass(double mass);
    void ClearLikelihoods();
    std::vector<double> GetWidths();
    std::vector<double> GetMasses();
    /// Make TGraphs
    void BookHistograms();
    void FillHistograms(double redMass, double relativeSF, double hadTopMassForWidthSF, double lepTopMassForWidthSF, bool isTTbar, bool isData, std::string catSuffix);
    void WriteHistograms(std::string histoFileName);
    void ConstructTGraphsFromHisto(std::string tGraphFileName, std::vector<std::string> datasetNames, std::vector<int> includeDataset);
    /// Get TGraphs (in order to calculate likelihood)
    bool ConstructTGraphsFromFile();
    bool ConstructTGraphsFromFile(std::string name);
    bool ConstructTGraphsFromFile(std::vector<std::string> datasetNames, std::vector<int> includeDataset);
    /// Calculate likelihood
    void CalculateLikelihood(double redMass, double relativeSF, bool isData);
    std::vector<double> CalculateLikelihood(double redMass, double relativeSF, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData);
    void CalculateCMLikelihood(double redMass, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData);
    void CalculateTempLikelihood(double redMass, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData);
    void CalculateGenLikelihood(double redMass, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData);
    /// Get output width
    void GetOutputWidth(double inputWidth, bool writeToFile = false, bool makeNewFile = false);
    void GetOutputWidth(double inputWidth, std::string type, bool writeToFile = false, bool makeNewFile = false);
    void GetOutputWidth(std::string inputFileName, double inputWidth, bool writeToFile = false, bool makeNewFile = false);
    void GetOutputWidth(std::string inputFileName, std::string inputDir, double inputWidth, bool writeToFile = false, bool makeNewFile = false);
    /// Use pseudo experiments
    int InitPull(int nPsExp);
    void AddPsExp(int thisPsExp, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData);
    void CalculatePull(double inputWidth);
    /// Calibration curve
    std::pair<double,double> ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma);
    /// Calculate fractions
    void AddToFraction(int d, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, bool isTTbar, bool isCM, bool isWM, bool isUM);
    void CalculateFractions(std::vector<std::string> datasetNames);
    /// Make 2D likelihood plot
    void Make2DGraph(std::string name, bool makeNewFile);
    void DrawOutputLogLikelihood(TGraph2D* g, TF2* f, double minX, double maxX, double minY, double maxY, std::string name, bool writeToFile);
    void DrawOutputLogLikelihood(TGraph2D* g, TF2* f, double minWidth, double minMass, double minimum, std::string name, bool writeToFile);
    /// Write output from likelihood calculation to file
    void PrintLikelihoodOutput(std::string llFileName);
    void PrintLikelihoodOutputData(std::string llFileName);
    void PrintMtmLikelihoodOutput(std::string llFileName);
      
  private:
    bool verbose_;
    bool rewHadOnly_;
    std::string outputDirName_;
    std::string dirNameTGraphTxt_;
    std::string dirNameNEvents_;
    std::string dirNameLLTxt_;
    std::string dirNamePull_;
    std::string inputFileName_;
    std::string suffix_;
    std::string histoName_;
    std::string histoNameSm_;
    std::string rangeRedMass_;
    double minRedMass_;
    double maxRedMass_;
    
    HelperTools *tls_;
    EventReweighting *rew_;
    
    static const double widthArray_[];
    static const double massArray_[];
    static const std::string listCats_[];
    
    static const int nWidths_;
    static const int nMasses_;
    static const int nCats_;
    
    static std::string stringWidthArray_[];
    static std::string stringMassArray_[];
    static std::string stringSuffix_[][61];
    
    static double loglike_[][61];
    static double loglike_data_[][61];
    static double loglike_per_evt_[][61];
    static double loglike_CM_[][61];
    static double loglike_CM_per_evt_[][61];
    static double loglike_temp_[][61];
    static double loglike_temp_per_evt_[][61];
    static double loglike_gen_[][61];
    static double loglike_gen_per_evt_[][61];
    
    static double loglike_pull_[][61][1000];
    static double loglike_pull_single_[][61];
    
    std::vector<double> vecWidths_;
    std::vector<double> vecMasses_;
    std::vector<double> vecLogLike_;
    
    std::string thisWidth_;
    std::string thisMass_;
    double thisWidthSF_;
    double thisMassSF_;
    TFile *file_;
    TFile *fileTGraphs_;
    TFile *filePlots_;
    TGraph2D *gLL2D_;
    std::ofstream txtOutput_;
    std::ofstream txtOutputLL_;
    std::ofstream txtOutputPsExp_;
    std::ofstream txtOutputFractions_;
    
    std::map<std::string,TH1D*> histo_;
    std::map<std::string,TH1D*> histoSm_;
    std::map<std::string,TH1D*> histoTotal_;
    std::map<std::string,TGraph*> graph_;
    std::map<std::string,std::vector<double>> vecBinCentres_;
    std::map<std::string,std::vector<double>> vecBinContents_;
    
    bool calledLLCalculation_;
    bool calledCMLLCalculation_;
    bool calledTempLLCalculation_;
    bool calledGenLLCalculation_;
    
    std::ifstream fileIn_;
    std::vector<double> vecWidthFromFile_;
    std::vector<double> vecLLValsFromFile_;
    std::pair<double,double> output_;
    
    static const double calCurvePar_[2];
    static const double calCurveParUnc_[2];
    
    int nPsExp_;
    static double nEventsCMFractions_[][25];
    static double nEventsWMFractions_[][25];
    static double nEventsUMFractions_[][25];
    
    int LocMinArray(int n, double* array);
    std::pair<int,int> LocMinArray(int n, double (*array)[61]);
    void ClearArray(int size, int* array);
    void ClearArray(int size, double* array);
    void ClearArray2D(int size, double (*array)[61]);
    void MakeTable(double* array, int n, double min, double max);
    Double_t DDparabola(Double_t *x, Double_t *par);
    Double_t Rosenbrock(Double_t *x, Double_t *par);
    void GetHistogram(int iCat, int iMass);
    void MakeGraph(int iCat, int iMass, int nPoints, double* centres, double* contents, std::string name, bool drawGraph = false);
    void MakeGraphSmooth(int iCat, int iMass, int nPoints, double* centres, double* contents, std::string name, bool drawGraph = false);
    bool ReadInput(std::string name);
    /// Get output width for pseudo experiments
    std::pair<double,double> GetOutputWidth(double inputWidth, int thisPsExp);
    /// Calculate output width (used in getters)
    std::pair<double,double> CalculateOutputWidth(std::string inputFileName, std::string inputDir, std::string plotName, bool writeToFile = false, bool makeNewFile = false);
    std::pair<double,double> CalculateOutputWidth(int nn, double (*LLvalues)[61], std::string plotName, bool writeToFile = false, bool makeNewFile = false);
    std::pair<double,double> CalculateOutputWidth(int nn, double* evalWidths, double (*LLvalues)[61], std::string plotName, bool writeToFile = false, bool makeNewFile = false);
    /// Calculate & get fractions for TGraphs
    void GetFractions(double *fractions, int nCats, std::vector<std::string> datasetNames, std::vector<int> includeDataset);
    void GetFractions(double *fractions, int nCats, std::string mass, std::vector<std::string> datasetNames, std::vector<int> includeDataset);
    /// Draw functions
    void DrawGraph(TH1D* h, TGraph* g, std::string name);
    void DrawGraph(TGraph2D* g, std::string name);
    void DrawLikelihoods();
    void DrawOutputLogLikelihood(TGraph* g, TF1* f, double minX, double maxX, double maxY, std::string name, bool writeToFile);
//     void DrawOutputLogLikelihood(TGraph2D* g, TF2* f, double minX, double maxX, double minY, double maxY, std::string name, bool writeToFile);
//     void DrawOutputLogLikelihood(TGraph2D* g, TF2* f, double minWidth, double minMass, double minimum, std::string name, bool writeToFile);
    /// Read file
    void ReadLLValuesFromFile(std::string inputFileName, std::string inputDir);
    /// Write file
    void WritePsExpOutput(std::pair<double,double> *outputWidth, std::pair<double,double> *inputWidth, double genWidth);
    void WriteFuncOutput(int nPoints, double *arrayCentre, double *arrayContent, std::string name);
    void WriteOutput(int nPoints, double width, double *arrayCentre, double *arrayContent, std::string name, int dim);
    void CombineOutput(std::string mass);
};


#endif
