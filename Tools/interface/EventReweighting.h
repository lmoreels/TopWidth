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
  
  private:
    bool useEvtCorr_;
    double corr0_;
    double corr1_;
    double genTopWidth_;
    double genTopMass_;
    double BreitWigner(double topMass, double scale);
    double BreitWignerNonRel(double topMass, double scale);
};


#endif
