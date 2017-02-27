#include "../interface/KinFitter.h"

KinFitter::KinFitter(bool addWConstr, bool addEqMassConstr):
rfFileName("/user/lmoreels/CMSSW_7_6_5/src/TopBrussels/TopWidth/PlotsForResolutionFunctions_testFit.root"), errorFuncMap()
{
  std::cout << "KinFitter::KinFitter - Initialising..." << std::endl;
  
  addWConstr_ = addWConstr;
  addEqMassConstr_ = addEqMassConstr;
  
  // Initialise error matrices
  mErrJet1.Zero();
  mErrJet2.Zero();
  mErrJet3.Zero();
  mErrJet4.Zero();
  mErrLepton.Zero();
  mErrNeutrino.Zero();
  
  this->GetErrorFunctions();
  this->SetupFitter();
}

KinFitter::KinFitter(std::string _rfFileName, bool addWConstr, bool addEqMassConstr):errorFuncMap()
{
  std::cout << "KinFitter::KinFitter - Initialising..." << std::endl;
  
  rfFileName = _rfFileName;
  addWConstr_ = addWConstr;
  addEqMassConstr_ = addEqMassConstr;
  
  // Initialise error matrices
  mErrJet1.Zero();
  mErrJet2.Zero();
  mErrJet3.Zero();
  mErrJet4.Zero();
  mErrLepton.Zero();
  mErrNeutrino.Zero();
  
  this->GetErrorFunctions();
  this->SetupFitter();
}

KinFitter::~KinFitter()
{
  // Destructor
}

void KinFitter::GetErrorFunctions()
{
  ResolutionFunctions *rf = new ResolutionFunctions(false, false);
  
  // At the moment only errors for Et
  errorFuncMap["bjetEt_B"]    = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "bjet", "B");
  errorFuncMap["bjetEt_O"]    = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "bjet", "O");
  errorFuncMap["bjetEt_E"]    = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "bjet", "E");
  errorFuncMap["nonbjetEt_B"] = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "B");
  errorFuncMap["nonbjetEt_O"] = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "O");
  errorFuncMap["nonbjetEt_E"] = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "E");
  
  delete rf;
}

TMatrixD KinFitter::SetErrors(TLorentzVector part, std::string type)
{
  mErr.Zero();
    
  if ( type.find("lepton") != std::string::npos || type.find("muon") != std::string::npos || type.find("electron") != std::string::npos || type.find("neutrino") != std::string::npos )
  {
    mErr(0,0) = 1.;
  }
  else
  {
    Double_t Et = part.Et();
    Double_t Eta = fabs(part.Eta());
    if ( Et > 250. ) Et = 250.;
    
    /// Errors on Et
    if ( Eta <= 1.3 )                   mErr(0,0) = pow( errorFuncMap[type+"Et_B"]->Eval(Et), 2);
    else if ( Eta > 1.3 && Eta <= 1.5 ) mErr(0,0) = pow( errorFuncMap[type+"Et_O"]->Eval(Et), 2);
    else                                mErr(0,0) = pow( errorFuncMap[type+"Et_E"]->Eval(Et), 2);
  }
  
  /// No errors on theta & phi atm
  mErr(1,1) = 1.; mErr(2,2) = 1.;
  
  return mErr;
}

void KinFitter::SetupFitter()
{
  /// Define particles
  jet1_ = new TFitParticleEtThetaPhiEMomFix( "jet1", "jet1", 0, &mErrJet1);
  jet2_ = new TFitParticleEtThetaPhiEMomFix( "jet2", "jet2", 0, &mErrJet2);
  jet3_ = new TFitParticleEtThetaPhiEMomFix( "jet3", "jet3", 0, &mErrJet3);
  jet4_ = new TFitParticleEtThetaPhiEMomFix( "jet4", "jet4", 0, &mErrJet4);
  lepton_ = new TFitParticleEtThetaPhiEMomFix( "lepton", "lepton", 0, &mErrLepton);
  neutrino_ = new TFitParticleEtThetaPhiEMomFix( "neutrino", "neutrino", 0, &mErrNeutrino);  // should be unmeasured part!!
  
  /// Define constraints
  consMW_ = new TFitConstraintM( "WMassConstraint", "WMassConstraint", 0, 0, 80.385);  // pdg2014
  consMW_->addParticles1(jet1_, jet2_);
  consMNu_ = new TFitConstraintM( "NeutrinoMassConstraint", "NeutrinoMassConstraint", 0, 0, 0.);
  consMNu_->addParticle1(neutrino_);
  consEqM_ = new TFitConstraintM( "EqualMassConstraint", "EqualMassConstraint", 0, 0, 0.);
  consEqM_->addParticles1(jet1_, jet2_, jet3_);
  consEqM_->addParticles2(jet4_, lepton_, neutrino_);
  
  /// Set up fitter
  //  Add constraints & particles
  fitter_ = new TKinFitter("WMassFit", "WMassFit");
  if (addWConstr_)
  {
    fitter_->addConstraint( consMW_ );
    if (! addEqMassConstr_) fitter_->addMeasParticles( jet1_, jet2_);
  }
  if (addEqMassConstr_)
  {
    fitter_->addConstraint( consEqM_ );
    fitter_->addConstraint( consMNu_ );
    fitter_->addMeasParticles( jet1_, jet2_, jet3_, jet4_, lepton_ );
    fitter_->addUnmeasParticle( neutrino_ );
  }
  
  /// Set convergence criteria
  fitter_->setMaxNbIter(30);
  fitter_->setMaxDeltaS(5e-5);
  fitter_->setMaxF(1e-4);
  fitter_->setVerbosity(0);
}

TKinFitter* KinFitter::doFit(TLorentzVector jet1, TLorentzVector jet2, int _verbosity)
{
  if (addEqMassConstr_)
    std::cerr << "KinFitter::doFit: Too few particles added to apply equal mass constraint!" << endl;
  
  /// Get errors
  mErrJet1 = SetErrors(jet1, "nonbjet");
  mErrJet2 = SetErrors(jet2, "nonbjet");
    
  /// Fill particles
  jet1_->setIni4Vec(&jet1); jet1_->setCovMatrix(&mErrJet1);
  jet2_->setIni4Vec(&jet2); jet2_->setCovMatrix(&mErrJet2);
  
  if ( _verbosity != 0 )
    fitter_->setVerbosity(_verbosity);
  
  /// Perform fit
  fitter_->fit();
  
  return fitter_;
}

TKinFitter* KinFitter::doFit(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector jet3, TLorentzVector jet4, TLorentzVector lepton, TLorentzVector neutrino, int _verbosity)
{
  /// Get errors
  mErrJet1 = SetErrors(jet1, "nonbjet");
  mErrJet2 = SetErrors(jet2, "nonbjet");
  if (addEqMassConstr_)
  {
    mErrJet3 = SetErrors(jet3, "nonbjet");
    mErrJet4 = SetErrors(jet4, "nonbjet");
    mErrLepton = SetErrors(lepton, "muon");
    mErrNeutrino = SetErrors(neutrino, "neutrino");
  }
  
  /// Fill particles
  jet1_->setIni4Vec(&jet1); jet1_->setCovMatrix(&mErrJet1);
  jet2_->setIni4Vec(&jet2); jet2_->setCovMatrix(&mErrJet2);
  if (addEqMassConstr_)
  {
    jet3_->setIni4Vec(&jet3); jet3_->setCovMatrix(&mErrJet3);
    jet4_->setIni4Vec(&jet4); jet4_->setCovMatrix(&mErrJet4);
    lepton_->setIni4Vec(&lepton); lepton_->setCovMatrix(&mErrLepton);
    //neutrino??
  }
  
  if ( _verbosity != 0 )
    fitter_->setVerbosity(_verbosity);
  
  /// Perform fit
  fitter_->fit();
  
  return fitter_;
}

vector<TLorentzVector> KinFitter::getCorrectedJets()
{
  tempJet.Clear();
  corrJets.clear();
  
  if ( fitter_->getStatus() == 0 )
  {
    tempJet.SetXYZM( jet1_->getCurr4Vec()->X(), jet1_->getCurr4Vec()->Y(), jet1_->getCurr4Vec()->Z(), jet1_->getCurr4Vec()->M() );
    corrJets.push_back(tempJet);
    tempJet.SetXYZM( jet2_->getCurr4Vec()->X(), jet2_->getCurr4Vec()->Y(), jet2_->getCurr4Vec()->Z(), jet2_->getCurr4Vec()->M() );
    corrJets.push_back(tempJet);
    if (addEqMassConstr_)
    {
      tempJet.SetXYZM( jet3_->getCurr4Vec()->X(), jet3_->getCurr4Vec()->Y(), jet3_->getCurr4Vec()->Z(), jet3_->getCurr4Vec()->M() );
      corrJets.push_back(tempJet);
      tempJet.SetXYZM( jet4_->getCurr4Vec()->X(), jet4_->getCurr4Vec()->Y(), jet4_->getCurr4Vec()->Z(), jet4_->getCurr4Vec()->M() );
      corrJets.push_back(tempJet);
      tempJet.SetXYZM( lepton_->getCurr4Vec()->X(), lepton_->getCurr4Vec()->Y(), lepton_->getCurr4Vec()->Z(), lepton_->getCurr4Vec()->M() );
      corrJets.push_back(tempJet);
      // neutrino?
    }
  }
  else std::cerr << "KinFitter::getCorrectedJets: Fit did not converge! Returning zero..." << std::endl;
  
  return corrJets;
}
