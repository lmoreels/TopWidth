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
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <errno.h>
#include <map>
#include <vector>
#include <TH2.h>
#include <TLorentzVector.h>


class TransferFunctions{
  public:
    TransferFunctions(bool calculateTransferFunctions);
    ~TransferFunctions();
    void bookHistograms();
    void fillJets(std::vector<TLorentzVector*> parton, std::vector<TLorentzVector*> jet);
    void fillMuon(TLorentzVector *genMu, TLorentzVector *recMu);
    void fillElectron(TLorentzVector *genEl, TLorentzVector *recEl);
    void writeOutputFiles();
    void writeTable();

    
  private:
    bool muon;
    bool electron;
    std::map<std::string,TH2F*> histoTr2D;
    std::string filename;

};



#endif
