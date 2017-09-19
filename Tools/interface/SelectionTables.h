#ifndef INTERFACE_SELECTIONTABLES_H
#define INTERFACE_SELECTIONTABLES_H

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <TMath.h>
#include <TFile.h>
#include <boost/algorithm/string/replace.hpp>

#include "TopTreeAnalysisBase/Content/interface/Dataset.h"


class SelectionTables{
  public:
    SelectionTables(std::vector<Dataset*> datasets);
    ~SelectionTables();
    void SetPrecision(int i);
    void SetLumi(double lumi);
    void SetEqLumi(int d, double eqLumi);
    void AddCutStep(std::string cutStepName);
    void SetUpTable();
    void Fill(int d, int cutStep, double value);
    void Fill(int d, std::vector<double> values);
    void CalculateTable(bool scaleEvents = true);
    void Write(std::string filename, bool writeError, bool writeMerged, bool writeLandscape, bool writeVertical = false);
    void Write(ofstream& fout, bool writeError, bool writeMerged, bool writeLandscape, bool writeVertical = false);
    
  private:
    int precision_; // nb of digits to display after the decimal point
    
    int numberOfCuts_;
    int numberOfDatasets_, numberOfDatasetsMerged_;
    
    std::vector<std::string> listOfCuts_;
    std::vector<Dataset*> listOfDatasets_, listOfDatasetsMerged_;
    
    double lumi_;
    double* eqLumi_;
    
    double** nofEventsRaw_;
    double** nofEventsRawError_;
    double** nofEventsRawVar_;
    double** nofEvents_;
    double** nofEventsError_;
    double** nofEventsVar_;

    double** nofEventsMerged_;
    double** nofEventsErrorMerged_;
    double** nofEventsVarMerged_;
    
    double* totalEvents_;
    double* totalEventsError_;
    double* totalEventsVar_;
    
    /* Calculate the binomial error on the efficiency p */
    double BinomialError(double p, double n) const { return (sqrt( (p*(1-p)/n)*p )); };
    /* Calculate the error on the number of selected events, rescaled by factor*/
    double ErrorCalculator( double number, double p, double factor) const { return (sqrt(number*(1-p))*factor); };
    
    void MergeDatasets();
    void WriteTable(ofstream& fout, double** listTable_,double** listTableError_, bool writeError, bool writeMerged, bool writeLandscape);
    void WriteTableVertical(ofstream& fout, double** listTable_,double** listTableError_, bool writeError, bool writeMerged);
  
};


#endif
