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
    KinFitter();
    KinFitter(std::string _rfFileName);
    ~KinFitter();
    TKinFitter* doFit(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector jet3);
      
  private:
    TKinFitter *fitter_;
    TAbsFitParticle *jet1_;
    TAbsFitParticle *jet2_;
    TAbsFitParticle *jet3_;
    TFitConstraintM *consMW_;
    TMatrixD mErrJet1 = TMatrixD(3,3);
    TMatrixD mErrJet2 = TMatrixD(3,3);
    TMatrixD mErrJet3 = TMatrixD(3,3);
    std::string rfFileName;
    map<std::string, TF1*> errorFuncMap;
    
    void GetErrorFunctions();
    void SetErrors(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector jet3);
};


#endif
