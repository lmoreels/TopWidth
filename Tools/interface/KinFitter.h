#ifndef INTERFACE_KINFITTER_H
#define INTERFACE_KINFITTER_H

#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>
#include <cmath>
#include <algorithm>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <TF1.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>

// load TopTreeAnalysisBase
#include "TopTreeAnalysisBase/KinFitter/interface/TKinFitter.h"
#include "TopTreeAnalysisBase/KinFitter/interface/TAbsFitParticle.h"
#include "TopTreeAnalysisBase/KinFitter/interface/TFitConstraintM.h"
#include "TopTreeAnalysisBase/KinFitter/interface/TFitParticleEtThetaPhi.h"

// load user defined
#include "ResolutionFunctions.h"


class KinFitter{
  public:
    KinFitter(bool addWConstr, bool addEqMassConstr);
    KinFitter(std::string _rfFileName, bool addWConstr, bool addEqMassConstr);
    ~KinFitter();
    TKinFitter* doFit(TLorentzVector jet1, TLorentzVector jet2, int _verbosity);
    TKinFitter* doFit(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector jet3, TLorentzVector jet4, TLorentzVector lepton, TLorentzVector neutrino, int _verbosity);
    std::vector<TLorentzVector> getCorrectedJets();
      
  private:
    bool addWConstr_;
    bool addEqMassConstr_;
    TKinFitter *fitter_;
    TAbsFitParticle *jet1_;
    TAbsFitParticle *jet2_;
    TAbsFitParticle *jet3_;
    TAbsFitParticle *jet4_;
    TAbsFitParticle *lepton_;
    TAbsFitParticle *neutrino_;
    TFitConstraintM *consMW_;
    TFitConstraintM *consMNu_;
    TFitConstraintM *consEqM_;
    TMatrixD mErr = TMatrixD(3,3);
    TMatrixD mErrJet1 = TMatrixD(3,3);
    TMatrixD mErrJet2 = TMatrixD(3,3);
    TMatrixD mErrJet3 = TMatrixD(3,3);
    TMatrixD mErrJet4 = TMatrixD(3,3);
    TMatrixD mErrLepton = TMatrixD(3,3);
    TMatrixD mErrNeutrino = TMatrixD(3,3);
    std::string rfFileName;
    std::map<std::string, TF1*> errorFuncMap;
    //std::map<std::string, TMatrixD(3,3)> mErr;
    TLorentzVector tempJet;
    std::vector<TLorentzVector> corrJets;
    
    void GetErrorFunctions();
    TMatrixD SetErrors(TLorentzVector part, std::string type);
    void SetupFitter();
};


#endif
