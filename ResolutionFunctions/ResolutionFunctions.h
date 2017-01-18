//
//  ResolutionFunctions.h
//  
//
//  Created by Lieselotte Moreels on 25/01/16.
//
//

#ifndef ResolutionFunctions_h
#define ResolutionFunctions_h

#include <stdio.h>
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
#include <TFile.h>
#include <TLorentzVector.h>


class ResolutionFunctions{
  public:
    ResolutionFunctions(bool calculateResolutionFunctions);
    ~ResolutionFunctions();
    std::string toStr(int number);
    void bookHistograms();
    void fillJets(std::vector<TLorentzVector> &parton, std::vector<TLorentzVector> &jet);
    void fillMuon(TLorentzVector genMu, TLorentzVector recMu);
    void fillElectron(TLorentzVector genEl, TLorentzVector recEl);
    void writeHistograms();
    void makeFit();
    void makeFit(std::string inputFileName, std::string outputFileName);
    void writeTable(std::string inputFileName);

    
  private:
    bool muon;
    bool electron;
    bool getHistos;
    std::map<std::string,TH2F*> histoRes2D;
    std::map<std::string,TH2F*> fitHisto2D;
    std::string inputFileName;
    int nHistos;
    static const std::string histoNames[];
    static const std::string histoDescription[];
    static Double_t dblGaus(Double_t *x, Double_t *par);

};



#endif
