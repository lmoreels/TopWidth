#include "../interface/EventReweighting.h"

EventReweighting::EventReweighting(bool evtCorr):
useEvtCorr_(evtCorr), corr0_(0.), corr1_(1.), genTopWidth_(1.31), genTopMass_(172.5)
{
  if (useEvtCorr_)  // from 76X !!
  {
    corr0_ = 0.0080432;
    corr1_ = 0.99195679;
  }
}

EventReweighting::~EventReweighting()
{
  
}

double EventReweighting::BreitWigner(double topMass, double scale)
{
  double BWmass = genTopMass_;
  double BWgamma = scale*genTopWidth_;
  double gammaterm = sqrt( pow(BWmass, 4) + pow(BWmass*BWgamma, 2) );
  double numerator = 2.*sqrt(2.)*BWmass*BWgamma*gammaterm/( TMath::Pi()*sqrt( pow(BWmass, 2) + gammaterm ) );
  double denominator = pow(pow(topMass, 2) - pow(BWmass, 2), 2) + pow(topMass, 4)*pow(BWgamma/BWmass, 2);
  
  return numerator/denominator;
}

double EventReweighting::BreitWignerNonRel(double topMass, double scale)
{
  double BWmass = genTopMass_;
  double BWgamma = scale*genTopWidth_/2.;
  double bw = BWgamma/( pow(topMass - BWmass, 2) + pow(BWgamma, 2) );
  
  return bw/(TMath::Pi());
}

double EventReweighting::MassBreitWigner(double topMass, double scale)
{
  double BWmass = scale;
  double BWgamma = genTopWidth_;
  double gammaterm = sqrt( pow(BWmass, 4) + pow(BWmass*BWgamma, 2) );
  double numerator = 2.*sqrt(2.)*BWmass*BWgamma*gammaterm/( TMath::Pi()*sqrt( pow(BWmass, 2) + gammaterm ) );
  double denominator = pow(pow(topMass, 2) - pow(BWmass, 2), 2) + pow(topMass, 4)*pow(BWgamma/BWmass, 2);
  
  return numerator/denominator;
}

double EventReweighting::MassBreitWignerNonRel(double topMass, double scale)
{
  double BWmass = scale;
  double BWgamma = genTopWidth_/2.;
  double bw = BWgamma/( pow(topMass - BWmass, 2) + pow(BWgamma, 2) );
  
  return bw/(TMath::Pi());
}

double EventReweighting::BBreitWigner(double topMass, double mass, double scaleW)
{
  double BWmass = mass;
  double BWgamma = scaleW*genTopWidth_;
  double gammaterm = sqrt( pow(BWmass, 4) + pow(BWmass*BWgamma, 2) );
  double numerator = 2.*sqrt(2.)*BWmass*BWgamma*gammaterm/( TMath::Pi()*sqrt( pow(BWmass, 2) + gammaterm ) );
  double denominator = pow(pow(topMass, 2) - pow(BWmass, 2), 2) + pow(topMass, 4)*pow(BWgamma/BWmass, 2);
  
  return numerator/denominator;
}

double EventReweighting::BBreitWignerNonRel(double topMass, double mass, double scaleW)
{
  double BWmass = mass;
  double BWgamma = scaleW*genTopWidth_/2.;
  double bw = BWgamma/( pow(topMass - BWmass, 2) + pow(BWgamma, 2) );
  
  return bw/(TMath::Pi());
}

double EventReweighting::EventWeightCalculator(double topMass, double scale)
{
  double corrNEvts = 1./(corr0_*scale+corr1_);
  return corrNEvts * BreitWigner(topMass, scale)/BreitWigner(topMass, 1.);
}

double EventReweighting::EventWeightCalculatorNonRel(double topMass, double scale)
{
  double corrNEvts = 1./(corr0_*scale+corr1_);
  return corrNEvts * BreitWignerNonRel(topMass, scale)/BreitWignerNonRel(topMass, 1.);
}

double EventReweighting::MassEventWeightCalculator(double topMass, double scale)
{
  double corrNEvts = 1./(corr0_*scale+corr1_);
  return corrNEvts * MassBreitWigner(topMass, scale)/MassBreitWigner(topMass, genTopMass_);
}

double EventReweighting::MassEventWeightCalculatorNonRel(double topMass, double scale)
{
  double corrNEvts = 1./(corr0_*scale+corr1_);
  return corrNEvts * MassBreitWignerNonRel(topMass, scale)/MassBreitWignerNonRel(topMass, genTopMass_);
}

double EventReweighting::BEventWeightCalculator(double topMass, double mass, double scaleW)
{
  double corrNEvts = 1./(corr0_*scaleW+corr1_);
  return corrNEvts * BBreitWigner(topMass, mass, scaleW)/BBreitWigner(topMass, genTopMass_, 1.);
}

double EventReweighting::BEventWeightCalculator(double topMass, double mass1, double mass2, double width1, double width2)
{
  /// Simultaneously reweight events from mass1 to mass2
  //                            and from width1 to width2
  double corrNEvts = 1./(corr0_*width2+corr1_);
  return corrNEvts * BBreitWigner(topMass, mass2, width2)/BBreitWigner(topMass, mass1, width1);
}

double EventReweighting::BEventWeightCalculatorNonRel(double topMass, double mass, double scaleW)
{
  double corrNEvts = 1./(corr0_*scaleW+corr1_);
  return corrNEvts * BBreitWignerNonRel(topMass, mass, scaleW)/BBreitWignerNonRel(topMass, genTopMass_, 1.);
}

double EventReweighting::BEventWeightCalculatorNonRel(double topMass, double mass1, double mass2, double width1, double width2)
{
  /// Simultaneously reweight events from mass1 to mass2
  //                            and from width1 to width2
  double corrNEvts = 1./(corr0_*width2+corr1_);
  return corrNEvts * BBreitWignerNonRel(topMass, mass2, width2)/BBreitWignerNonRel(topMass, mass1, width1);
}
