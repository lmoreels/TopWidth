//
//  ResolutionFunctions.h
//  
//
//  Created by Lieselotte Moreels on 25/01/16.
//
//

#ifndef INTERFACE_RESOLUTIONFUNCTIONS_H
#define INTERFACE_RESOLUTIONFUNCTIONS_H

#include <stdio.h>
#include <TROOT.h>
#include <TStyle.h>
#include <cmath>
#include <TMath.h>
#include <Math/GenVector/VectorUtil.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <errno.h>
#include <map>
#include <array>
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF12.h>
#include <TF2.h>
#include <TFile.h>
#include <TLorentzVector.h>


class ResolutionFunctions{
  public:
    ResolutionFunctions(bool calculateResolutionFunctions, bool _verbose);
    ~ResolutionFunctions();
    std::string toStr(int number);
    void bookHistograms();
    void fillJets(std::vector<TLorentzVector> &parton, std::vector<TLorentzVector> &jet);
    void fillMuon(TLorentzVector genMu, TLorentzVector recMu);
    void fillElectron(TLorentzVector genEl, TLorentzVector recEl);
    void writeHistograms();
    void makeFit();
    void makeFit(std::string inputFileName, std::string outputFileName);
    void makeFit(std::string inputFileName, std::string outputFileName, bool simplify);
    TF2* getFitFunction2D(std::string inputFileName, std::string varName, std::string objName, std::string binName);
    TF1* getFitFunction1D(std::string inputFileName, std::string varName, std::string objName, std::string binName);
    TF1* getResolutionFunction(std::string inputFileName, std::string varName, std::string objName, std::string binName);
    double getResolution(std::string inputFileName, std::string varName, std::string objName, double var, std::string binName);
    void writeTable(std::string inputFileName);
  
  
  private:
    bool verbose;
    bool muon;
    bool electron;
    bool getHistos;
    bool useSingleG;
    std::map<std::string,TH2F*> histoRes2D;
    std::map<std::string,TH2F*> fitHisto2D;
    //std::map<std::string,std::vector<std::array<double, 2> > > fitParams;
    std::string inputFileName;
    int nHistos;
    unsigned int nJets;
    static const std::string histoNames[];
    static const std::string histoDescription[];
    static Double_t sGaus(Double_t *x, Double_t *par);
    static Double_t dblGaus(Double_t *x, Double_t *par);
    static Double_t dblGausParFill(Double_t *x, Double_t *par);
    std::vector<std::array<double, 2> > getParameters(std::string inputFileName, std::string varName, std::string objName, std::string binName);
};



#endif
