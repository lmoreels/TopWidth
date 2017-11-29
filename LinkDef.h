#ifdef __CINT__
#include "Tools/interface/EventReweighting.h"
#include "Tools/interface/HelperTools.h"
#include "Tools/interface/KinFitter.h"
#include "Tools/interface/Likelihood.h"
#include "Tools/interface/Likelihood2D.h"
#include "Tools/interface/LikelihoodMass.h"
#include "Tools/interface/ResolutionFunctions.h"
#include "Tools/interface/SelectionTables.h"
#include "Tools/interface/Trigger.h"
#else
#include "Tools/interface/EventReweighting.h"
#include "Tools/interface/HelperTools.h"
#include "Tools/interface/KinFitter.h"
#include "Tools/interface/Likelihood.h"
#include "Tools/interface/Likelihood2D.h"
#include "Tools/interface/LikelihoodMass.h"
#include "Tools/interface/ResolutionFunctions.h"
#include "Tools/interface/SelectionTables.h"
#include "Tools/interface/Trigger.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class EventReweighting+;
#pragma link C++ class HelperTools+;
#pragma link C++ class KinFitter+;
#pragma link C++ class Likelihood+;
#pragma link C++ class Likelihood2D+;
#pragma link C++ class LikelihoodMass+;
#pragma link C++ class ResolutionFunctions+;
#pragma link C++ class SelectionTables+;
#pragma link C++ class Trigger+;

#endif
