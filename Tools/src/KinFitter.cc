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
  
  errorFuncMap["bjetEt_B"]    = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "bjet", "B");
  errorFuncMap["bjetEt_O"]    = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "bjet", "O");
  errorFuncMap["bjetEt_E"]    = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "bjet", "E");
  errorFuncMap["nonbjetEt_B"] = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "B");
  errorFuncMap["nonbjetEt_O"] = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "O");
  errorFuncMap["nonbjetEt_E"] = (TF1*) rf->getResolutionFunction(rfFileName, "Et", "nonbjet", "E");
  
  errorFuncMap["bjetTheta_B"]    = (TF1*) rf->getResolutionFunction(rfFileName, "theta", "bjet", "B");
  errorFuncMap["bjetTheta_O"]    = (TF1*) rf->getResolutionFunction(rfFileName, "theta", "bjet", "O");
  errorFuncMap["bjetTheta_E"]    = (TF1*) rf->getResolutionFunction(rfFileName, "theta", "bjet", "E");
  errorFuncMap["nonbjetTheta_B"]    = (TF1*) rf->getResolutionFunction(rfFileName, "theta", "nonbjet", "B");
  errorFuncMap["nonbjetTheta_O"]    = (TF1*) rf->getResolutionFunction(rfFileName, "theta", "nonbjet", "O");
  errorFuncMap["nonbjetTheta_E"]    = (TF1*) rf->getResolutionFunction(rfFileName, "theta", "nonbjet", "E");
  
  errorFuncMap["bjetPhi_B"]    = (TF1*) rf->getResolutionFunction(rfFileName, "phi", "bjet", "B");
  errorFuncMap["bjetPhi_O"]    = (TF1*) rf->getResolutionFunction(rfFileName, "phi", "bjet", "O");
  errorFuncMap["bjetPhi_E"]    = (TF1*) rf->getResolutionFunction(rfFileName, "phi", "bjet", "E");
  errorFuncMap["nonbjetPhi_B"]    = (TF1*) rf->getResolutionFunction(rfFileName, "phi", "nonbjet", "B");
  errorFuncMap["nonbjetPhi_O"]    = (TF1*) rf->getResolutionFunction(rfFileName, "phi", "nonbjet", "O");
  errorFuncMap["nonbjetPhi_E"]    = (TF1*) rf->getResolutionFunction(rfFileName, "phi", "nonbjet", "E");
  
  delete rf;
}

TMatrixD KinFitter::SetErrors(TLorentzVector part, std::string type)
{
  mErr.Zero();
  mErr(1,1) = 1e-6;
  mErr(2,2) = 1e-6;
    
  if ( type.find("lepton") != std::string::npos || type.find("muon") != std::string::npos || type.find("electron") != std::string::npos )
  {
    //if ( part.Et() < 100 ) mErr(0,0) = 0.5;
    //else 
      mErr(0,0) = 1e-3;
  }
  else if ( type.find("neutrino") != std::string::npos )
  {
    mErr(0,0) = 200.;
    mErr(1,1) = 9999.;
    mErr(2,2) = 1e-4;
  }
  else
  {
    Double_t Et = part.Et();
    Double_t Eta = fabs(part.Eta());
    if ( Et > 250. ) Et = 250.;
    
    /// Get errors from functions
    if ( Eta <= 1.3 )
    {
      mErr(0,0) = pow( errorFuncMap[type+"Et_B"]->Eval(Et)   , 2);
      mErr(1,1) = pow( errorFuncMap[type+"Theta_B"]->Eval(Et), 2);
      mErr(2,2) = pow( errorFuncMap[type+"Phi_B"]->Eval(Et)  , 2);
    }
    else if ( Eta > 1.3 && Eta <= 1.5 )
    {
      mErr(0,0) = pow( errorFuncMap[type+"Et_O"]->Eval(Et)   , 2);
      mErr(1,1) = pow( errorFuncMap[type+"Theta_O"]->Eval(Et), 2);
      mErr(2,2) = pow( errorFuncMap[type+"Phi_O"]->Eval(Et)  , 2);
    }
    else
    {
      mErr(0,0) = pow( errorFuncMap[type+"Et_E"]->Eval(Et)   , 2);
      mErr(1,1) = pow( errorFuncMap[type+"Theta_E"]->Eval(Et), 2);
      mErr(2,2) = pow( errorFuncMap[type+"Phi_E"]->Eval(Et)  , 2);
    }
  }
  
  return mErr;
}

void KinFitter::SetupFitter()
{
  /// Define particles
  jet1_ = new TFitParticleEtThetaPhi( "jet1", "jet1", 0, &mErrJet1);
  jet2_ = new TFitParticleEtThetaPhi( "jet2", "jet2", 0, &mErrJet2);
  jet3_ = new TFitParticleEtThetaPhi( "jet3", "jet3", 0, &mErrJet3);
  jet4_ = new TFitParticleEtThetaPhi( "jet4", "jet4", 0, &mErrJet4);
  lepton_ = new TFitParticleEtThetaPhiEMomFix( "lepton", "lepton", 0, &mErrLepton);
  neutrino_ = new TFitParticleEtThetaPhiEMomFix( "neutrino", "neutrino", 0, &mErrNeutrino);  // should be unmeasured part!!
  
  /// Define constraints
  consMWhad_ = new TFitConstraintM( "WMassConstraint", "WMassConstraint", 0, 0, 80.385);  // pdg2017
  consMWhad_->addParticles1(jet1_, jet2_);
  consMWlep_ = new TFitConstraintM( "WMassConstraint_lep", "WMassConstraint_lep", 0, 0, 80.385);  // pdg2017
  consMWlep_->addParticles1(lepton_, neutrino_);
  consMNu_ = new TFitConstraintM( "NeutrinoMassConstraint", "NeutrinoMassConstraint", 0, 0, 0.);
  consMNu_->addParticle1(neutrino_);
  consEqM_ = new TFitConstraintM( "EqualMassConstraint", "EqualMassConstraint", 0, 0, 0.);
  consEqM_->addParticles1(jet1_, jet2_, jet3_);
  consEqM_->addParticles2(jet4_, lepton_, neutrino_);
  consSumPx_ = new TFitConstraintEp("sumPx", "sumPx", 0, TFitConstraintEp::pX, 0.);
  consSumPx_->addParticles(jet1_, jet2_, jet3_, jet4_, lepton_, neutrino_);
  consSumPy_ = new TFitConstraintEp("sumPy", "sumPy", 0, TFitConstraintEp::pY, 0.);
  consSumPy_->addParticles(jet1_, jet2_, jet3_, jet4_, lepton_, neutrino_);
  
  /// Set up fitter
  //  Add constraints & particles
  fitter_ = new TKinFitter("WMassFit", "WMassFit");
  if (addEqMassConstr_)
  {
    fitter_->addConstraint( consMWhad_ );
    fitter_->addConstraint( consMWlep_);
    fitter_->addConstraint( consEqM_ );
    //fitter_->addConstraint( consMNu_ );
    fitter_->addConstraint( consSumPx_ );
    fitter_->addConstraint( consSumPy_ );
    fitter_->addMeasParticles( jet1_, jet2_, jet3_, jet4_, lepton_, neutrino_ );
    //fitter_->addUnmeasParticle( neutrino_ );
    //fitter_->addMeasParticle( neutrino_ );
  }
  else if (addWConstr_)
  {
    fitter_->addConstraint( consMWhad_ );
    fitter_->addMeasParticles( jet1_, jet2_);
  }
  
  /// Set convergence criteria
  fitter_->setMaxNbIter(30);
  fitter_->setMaxDeltaS(1e-5);
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
    mErrJet3 = SetErrors(jet3, "bjet");
    mErrJet4 = SetErrors(jet4, "bjet");
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
    neutrino_->setIni4Vec(&neutrino); neutrino_->setCovMatrix(&mErrNeutrino);
  }
  
  /// Set Px & Py constraints
  if (addEqMassConstr_)
  {
    consSumPx_->setConstraint( jet1.Px() + jet2.Px() + jet3.Px() + jet4.Px() + lepton.Px() + neutrino.Px() );
    consSumPy_->setConstraint( jet1.Py() + jet2.Py() + jet3.Py() + jet4.Py() + lepton.Py() + neutrino.Py() );
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
      tempJet.SetXYZM( neutrino_->getCurr4Vec()->X(), neutrino_->getCurr4Vec()->Y(), neutrino_->getCurr4Vec()->Z(), neutrino_->getCurr4Vec()->M() );
      corrJets.push_back(tempJet);
    }
  }
  else std::cerr << "KinFitter::getCorrectedJets: Fit did not converge! Returning zero..." << std::endl;
  
  return corrJets;
}
