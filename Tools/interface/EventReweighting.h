#ifndef INTERFACE_EVENTREWEIGHTING_H
#define INTERFACE_EVENTREWEIGHTING_H

#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>
#include <cmath>
#include <TMath.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>


class EventReweighting{
  public:
    EventReweighting(bool evtCorr);
    ~EventReweighting();
    double EventWeightCalculator(double topMass, double scale);
    double EventWeightCalculatorNonRel(double topMass, double scale);
    double MassEventWeightCalculator(double topMass, double scale);
    double MassEventWeightCalculatorNonRel(double topMass, double scale);
    double BEventWeightCalculator(double topMass, double mass, double scaleW);
    double BEventWeightCalculator(double topMass, double mass1, double mass2, double width1, double width2);
    double BEventWeightCalculatorNonRel(double topMass, double mass, double scaleW);
    double BEventWeightCalculatorNonRel(double topMass, double mass1, double mass2, double width1, double width2);
  
  private:
    bool useEvtCorr_;
    double corr0_;
    double corr1_;
    double genTopWidth_;
    double genTopMass_;
    double BreitWigner(double topMass, double scale);
    double BreitWignerNonRel(double topMass, double scale);
    double MassBreitWigner(double topMass, double scale);
    double MassBreitWignerNonRel(double topMass, double scale);
    double BBreitWigner(double topMass, double mass, double scaleW);
    double BBreitWignerNonRel(double topMass, double mass, double scaleW);
};


#endif
