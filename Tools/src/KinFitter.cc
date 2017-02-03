#include "../interface/KinFitter.h"

KinFitter::KinFitter():
rfFileName("/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/PlotsForResolutionFunctions_testFit.root")
{
  // Initialise error matrices
  mErrJet1.Zero();
  mErrJet2.Zero();
  mErrJet3.Zero();
  
  this->GetErrorFunctions();
}

KinFitter::KinFitter(std::string _rfFileName)
{
  rfFileName = _rfFileName;
  
  KinFitter();
}

KinFitter::~KinFitter()
{
  // Destructor
}

void KinFitter::GetErrorFunctions()
{
  ResolutionFunctions *rf = new ResolutionFunctions(false, false);
  
  // At the moment only errors for Et
  errorFuncMap["bjetEt_B"]    = rf->getResolutionFunction(rfFileName, "Et", "bjet", "B");
  errorFuncMap["bjetEt_O"]    = rf->getResolutionFunction(rfFileName, "Et", "bjet", "O");
  errorFuncMap["bjetEt_E"]    = rf->getResolutionFunction(rfFileName, "Et", "bjet", "E");
  errorFuncMap["nonbjetEt_B"] = rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "B");
  errorFuncMap["nonbjetEt_O"] = rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "O");
  errorFuncMap["nonbjetEt_E"] = rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "E");
}

void KinFitter::SetErrors(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector jet3)
{
  Double_t Et1 = jet1.Et();
  Double_t Et2 = jet2.Et();
  Double_t Et3 = jet3.Et();
  Double_t Eta1 = fabs(jet1.Eta());
  Double_t Eta2 = fabs(jet2.Eta());
  Double_t Eta3 = fabs(jet3.Eta());
  if ( Et1 > 250 ) Et1 = 250;
  if ( Et2 > 250 ) Et2 = 250;
  if ( Et3 > 250 ) Et3 = 250;
  
  /// jet 1 & jet 2 "from W boson" ==> light jets; jet 3 ==> b jet
  if ( Eta1 <= 1.3 )                    mErrJet1(0,0) = pow( errorFuncMap["nonbjetEt_B"]->Eval(Et1), 2);
  else if ( Eta1 > 1.3 && Eta1 <= 1.5 ) mErrJet1(0,0) = pow( errorFuncMap["nonbjetEt_O"]->Eval(Et1), 2);
  else                                  mErrJet1(0,0) = pow( errorFuncMap["nonbjetEt_E"]->Eval(Et1), 2);
  
  if ( Eta2 <= 1.3 )                    mErrJet2(0,0) = pow( errorFuncMap["nonbjetEt_B"]->Eval(Et2), 2);
  else if ( Eta2 > 1.3 && Eta2 <= 1.5 ) mErrJet2(0,0) = pow( errorFuncMap["nonbjetEt_O"]->Eval(Et2), 2);
  else                                  mErrJet2(0,0) = pow( errorFuncMap["nonbjetEt_E"]->Eval(Et2), 2);
  
  if ( Eta3 <= 1.3 )                    mErrJet3(0,0) = pow( errorFuncMap["bjetEt_B"]->Eval(Et3), 2);
  else if ( Eta3 > 1.3 && Eta3 <= 1.5 ) mErrJet3(0,0) = pow( errorFuncMap["bjetEt_O"]->Eval(Et3), 2);
  else                                  mErrJet3(0,0) = pow( errorFuncMap["bjetEt_E"]->Eval(Et3), 2);
  
  /// No errors on theta & phi atm
  mErrJet1(1,1) = 1.; mErrJet1(2,2) = 1.;
  mErrJet2(1,1) = 1.; mErrJet2(2,2) = 1.;
  mErrJet3(1,1) = 1.; mErrJet3(2,2) = 1.;
}

TKinFitter* KinFitter::doFit(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector jet3)
{
  this->SetErrors(jet1, jet2, jet3);
  
  /// Define particles
  jet1_ = new TFitParticleEtThetaPhi( "jet1", "jet1", &jet1, &mErrJet1);
  jet2_ = new TFitParticleEtThetaPhi( "jet2", "jet2", &jet2, &mErrJet2);
  jet3_ = new TFitParticleEtThetaPhi( "jet3", "jet3", &jet3, &mErrJet3);
  
  /// Define constraints
  consMW_ = new TFitConstraintM( "WMassConstraint", "WMassConstraint", 0, 0, 80.385);  // pdg2014
  consMW_->addParticles1(jet1_, jet2_);
  
  /// Set up fitter
  //  Add constraints & particles
  fitter_ = new TKinFitter("WMassFit", "WMassFit");
  fitter_->addConstraint( consMW_ );
  fitter_->addMeasParticles( jet1_, jet2_, jet3_ );
  
  /// Set convergence criteria (?)  --> functions do not exist in our kinfitter?
  fitter_->setMaxNbIter(30);
  fitter_->setMaxDeltaS(5e-5);
  fitter_->setMaxF(1e-4);
  fitter_->setVerbosity(0);
  
  /// Perform fit
  fitter_->fit();
  
  return fitter_;
}