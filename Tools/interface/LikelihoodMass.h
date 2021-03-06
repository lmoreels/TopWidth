#ifndef INTERFACE_LIKELIHOODMASS_H
#define INTERFACE_LIKELIHOODMASS_H

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
#include <TGraphSmooth.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
//#include <TPaveStats.h>
//#include <TList.h>
//#include <TText.h>

// user defined
#include "EventReweighting.h"
#include "HelperTools.h"


class LikelihoodMass{
  public:
    LikelihoodMass(double min, double max, std::string outputDirName, std::string date, bool useHadTopOnly, bool makeHistograms, bool verbose);
    ~LikelihoodMass();
    void ClearLikelihoods();
    std::vector<double> GetMasses();
    /// Make TGraphs
    void BookHistograms();
    void FillHistograms(double redMass, double relativeSF, double hadTopMassForMassSF, double lepTopMassForMassSF, bool isTTbar, bool isData, std::string catSuffix);
    void WriteHistograms(std::string histoFileName);
    void ConstructTGraphsFromHisto(std::string tGraphFileName, std::vector<std::string> datasetNames, std::vector<int> includeDataset);
    /// Get TGraphs (in order to calculate likelihood)
    bool ConstructTGraphsFromFile();
    bool ConstructTGraphsFromFile(std::string name);
    bool ConstructTGraphsFromFile(std::vector<std::string> datasetNames, std::vector<int> includeDataset);
    /// Calculate likelihood
    void CalculateLikelihood(double redMass, double relativeSF, bool isData);
    std::vector<double> CalculateLikelihood(double redMass, double relativeSF, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData);
    std::vector<double> CalculateLikelihood(double redMass, double relativeSF, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputWidth, double inputMass, bool isTTbar, bool isData);
    void CalculateCMLikelihood(double redMass, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData);
    void CalculateCMLikelihood(double redMass, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputWidth, double inputMass, bool isTTbar, bool isData);
    void CalculateTempLikelihood(double redMass, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData);
    void CalculateTempLikelihood(double redMass, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputWidth, double inputMass, bool isTTbar, bool isData);
    void CalculateGenLikelihood(double redMass, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData);
    /// Get output mass
    void GetOutputMass(double inputWidth, double inputMass, bool writeToFile = false, bool makeNewFile = false);
    void GetOutputMass(double inputWidth, double inputMass, std::string type, bool writeToFile = false, bool makeNewFile = false);
    void GetOutputMass(std::string inputFileName, double inputWidth, double inputMass, bool writeToFile = false, bool makeNewFile = false);
    void GetOutputMass(std::string inputFileName, std::string inputDir, double inputWidth, double inputMass, bool writeToFile = false, bool makeNewFile = false);
    /// Use pseudo experiments
    int InitPull(int nPsExp);
    void AddPsExp(int thisPsExp, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData);
    void AddPsExp(int thisPsExp, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputWidth, double inputMass, bool isTTbar, bool isData);
    void CalculatePull(double inputMass);
    void CalculatePull(double inputWidth, double inputMass);
    /// Calibration curve
    std::pair<double,double> ApplyCalibrationCurve(double thisOutputMass, double thisOutputMassSigma);
    /// Calculate fractions
    void AddToFraction(int d, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, bool isTTbar, bool isCM, bool isWM, bool isUM);
    void CalculateFractions(std::vector<std::string> datasetNames);
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
    
    static const double massArray_[];
    static const std::string listCats_[];
    
    static const int nMasses_;
    static const int nCats_;
    
    static std::string stringMassArray_[];
    static std::string stringSuffix_[];
    
    static double loglike_[];
    static double loglike_data_[];
    static double loglike_per_evt_[];
    static double loglike_CM_[];
    static double loglike_CM_per_evt_[];
    static double loglike_temp_[];
    static double loglike_temp_per_evt_[];
    static double loglike_gen_[];
    static double loglike_gen_per_evt_[];
    
    static double loglike_pull_[][1000];
    static double loglike_pull_single_[];
    
    std::vector<double> vecMasses_;
    std::vector<double> vecLogLike_;
    
    std::string thisMass_;
    double thisMassSF_;
    TFile *file_;
    TFile *fileTGraphs_;
    TFile *filePlots_;
    //TGraph2D *gLL2D_;
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
    std::vector<double> vecMassFromFile_;
    std::vector<double> vecLLValsFromFile_;
    std::pair<double,double> output_;
    
    static const double calCurvePar_[2];
    static const double calCurveParUnc_[2];
    
    int nPsExp_;
    static double nEventsCMFractions_[][25];
    static double nEventsWMFractions_[][25];
    static double nEventsUMFractions_[][25];
    
    int LocMinArray(int n, double* array);
    void ClearArray(int size, int* array);
    void ClearArray(int size, double* array);
    void ClearArray2D(int size, double (*array)[3]);
    void MakeTable(double* array, int n, double min, double max);
    Double_t DDparabola(Double_t *x, Double_t *par);
    void GetHistogram(int iCat);
    void MakeGraph(int iCat, int nPoints, double* centres, double* contents, std::string name, bool drawGraph = false);
    void MakeGraphSmooth(int iCat, int nPoints, double* centres, double* contents, std::string name, bool drawGraph = false);
    bool ReadInput(std::string name);
    /// Get output mass for pseudo experiments
    std::pair<double,double> GetOutputMass(double inputWidth, double inputMass, int thisPsExp);
    /// Calculate output mass (used in getters)
    std::pair<double,double> CalculateOutputMass(std::string inputFileName, std::string inputDir, std::string plotName, bool writeToFile = false, bool makeNewFile = false);
    std::pair<double,double> CalculateOutputMass(int nn, double* LLvalues, std::string plotName, bool writeToFile = false, bool makeNewFile = false);
    std::pair<double,double> CalculateOutputMass(int nn, double* evalMasses, double* LLvalues, std::string plotName, bool writeToFile = false, bool makeNewFile = false);
    /// Calculate & get fractions for TGraphs
    void GetFractions(double *fractions, int nCats, std::vector<std::string> datasetNames, std::vector<int> includeDataset);
    void GetFractions(double *fractions, int nCats, std::string mass, std::vector<std::string> datasetNames, std::vector<int> includeDataset);
    /// Draw functions
    void DrawGraph(TH1D* h, TGraph* g, std::string name);
    void DrawGraph(TGraph2D* g, std::string name);
    void DrawLikelihoods();
    void DrawOutputLogLikelihood(TGraph* g, TF1* f, double minX, double maxX, double maxY, std::string name, bool writeToFile);
    /// Read file
    void ReadLLValuesFromFile(std::string inputFileName, std::string inputDir);
    /// Write file
    void WritePsExpOutput(std::pair<double,double> *outputMass, std::pair<double,double> *inputMass, double genMass);
    void WriteFuncOutput(int nPoints, double *arrayCentre, double *arrayContent, std::string name);
    void WriteOutput(int nPoints, double mass, double *arrayCentre, double *arrayContent, std::string name, int dim);
    void CombineOutput();
};


#endif
