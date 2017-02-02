#include "../interface/KinFitter.h"

KinFitter::KinFitter(std::string _rfFileName):rfFileName("/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/PlotsForResolutionFunctions_testFit.root")
{
  rfFileName = _rfFileName;
  
  // Initialise error matrices
  mErrJet1.Zero();
  mErrJet2.Zero();
  mErrJet3.Zero();
}

KinFitter::~KinFitter()
{
  // Destructor
}

void KinFitter::GetErrorFunctions()
{
  ResolutionFunctions *rf = new ResolutionFunctions(false);
  
  // At the moment only errors for Et
  errorFuncMap["bjetEt_B"]    = rf->getResolutionFunction1D(rfFileName, "Et", "bjet", "B", true);
  errorFuncMap["bjetEt_O"]    = rf->getResolutionFunction1D(rfFileName, "Et", "bjet", "O", true);
  errorFuncMap["bjetEt_E"]    = rf->getResolutionFunction1D(rfFileName, "Et", "bjet", "E", true);
  errorFuncMap["nonbjetEt_B"] = rf->getResolutionFunction1D(rfFileName, "Et", "nonbjet", "B", true);
  errorFuncMap["nonbjetEt_O"] = rf->getResolutionFunction1D(rfFileName, "Et", "nonbjet", "O", true);
  errorFuncMap["nonbjetEt_E"] = rf->getResolutionFunction1D(rfFileName, "Et", "nonbjet", "E", true);
}

void KinFitter::SetErrors(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector jet3)
{
  this->GetErrorFunctions();
  
  Double_t Et1 = jet1.Et();
  Double_t Et2 = jet2.Et();
  Double_t Et3 = jet3.Et();
  Double_t Eta1 = jet1.Eta();
  Double_t Eta2 = jet2.Eta();
  Double_t Eta3 = jet3.Eta();
  if ( Et1 > 200 ) Et1 = 200;
  if ( Et2 > 200 ) Et2 = 200;
  if ( Et3 > 200 ) Et3 = 200;
  
  /// jet 1 & jet 2 "from W boson" ==> light jets; jet 3 ==> b jet
  if ( Eta1 <= 1.3 )                    mErrJet1(0,0) = errorFuncMap["nonbjetE_B"]->Eval(Et1);
  else if ( Eta1 > 1.3 && Eta1 <= 1.5 ) mErrJet1(0,0) = errorFuncMap["nonbjetE_O"]->Eval(Et1);
  else                                  mErrJet1(0,0) = errorFuncMap["nonbjetE_E"]->Eval(Et1);
  
  if ( Eta2 <= 1.3 )                    mErrJet2(0,0) = errorFuncMap["nonbjetE_B"]->Eval(Et2);
  else if ( Eta2 > 1.3 && Eta2 <= 1.5 ) mErrJet2(0,0) = errorFuncMap["nonbjetE_O"]->Eval(Et2);
  else                                  mErrJet2(0,0) = errorFuncMap["nonbjetE_E"]->Eval(Et2);
  
  if ( Eta3 <= 1.3 )                    mErrJet3(0,0) = errorFuncMap["bjetE_B"]->Eval(Et3);
  else if ( Eta3 > 1.3 && Eta3 <= 1.5 ) mErrJet3(0,0) = errorFuncMap["bjetE_O"]->Eval(Et3);
  else                                  mErrJet3(0,0) = errorFuncMap["bjetE_E"]->Eval(Et3);  
}

void KinFitter::doFit(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector jet3)
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
  
  /// Perform fit
  fitter_->fit();
  
}
