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

double EventReweighting::BreitWignerNonRel(double topMass, double scale)
{
  double BWmass = genTopMass_;
  double BWgamma = scale*genTopWidth_/2.;
  double bw = BWgamma/( pow(topMass - BWmass, 2) + pow(BWgamma, 2) );
  
  return bw/(TMath::Pi());
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

double EventReweighting::EventWeightCalculator(double topMass, double scale)
{
  double corrNEvts = 1./(corr0_*scale+corr1_);
  return corrNEvts * BreitWigner(topMass, scale)/BreitWigner(topMass, 1.);
}
