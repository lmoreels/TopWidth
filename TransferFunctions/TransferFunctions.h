//
//  TransferFunctions.h
//  
//
//  Created by Lieselotte Moreels on 25/01/16.
//
//

#ifndef TransferFunctions_h
#define TransferFunctions_h

#include <stdio.h>
#include <cmath>
#include <TMath.h>
#include <Math/GenVector/VectorUtil.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
//#include <ofstream>
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


class TransferFunctions{
  public:
    TransferFunctions(bool calculateTransferFunctions);
    ~TransferFunctions();
    std::string toStr(int number);
    void bookHistograms();
    void fillJets(std::vector<TLorentzVector> &parton, std::vector<TLorentzVector> &jet);
    void fillMuon(TLorentzVector *genMu, TLorentzVector *recMu);
    void fillElectron(TLorentzVector *genEl, TLorentzVector *recEl);
    void writeHistograms();
    void makeFit();
    void writeTable(std::string inputFileName);

    
  private:
    bool muon;
    bool electron;
    std::map<std::string,TH2F*> histoTr2D;
    std::string inputFileName;
    int nHistos;
    static const std::string histoNames[];
    static const std::string histoDescription[];

};



#endif
