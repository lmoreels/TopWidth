#ifndef INTERFACE_TRIGGER_H
#define INTERFACE_TRIGGER_H

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <errno.h>
#include <iostream>
#include <map>

//user code
#include "/user/lmoreels/CMSSW_7_4_15_patch1/src/TopBrussels/TopTreeProducer/interface/TRootRun.h"
#include "/user/lmoreels/CMSSW_7_4_15_patch1/src/TopBrussels/TopTreeProducer/interface/TRootEvent.h"

#include "/user/lmoreels/CMSSW_7_4_15_patch1/src/TopBrussels/TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "/user/lmoreels/CMSSW_7_4_15_patch1/src/TopBrussels/TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "/user/lmoreels/CMSSW_7_4_15_patch1/src/TopBrussels/TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "/user/lmoreels/CMSSW_7_4_15_patch1/src/TopBrussels/TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "/user/lmoreels/CMSSW_7_4_15_patch1/src/TopBrussels/TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "/user/lmoreels/CMSSW_7_4_15_patch1/src/TopBrussels/TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "/user/lmoreels/CMSSW_7_4_15_patch1/src/TopBrussels/TopTreeAnalysisBase/Selection/interface/Run2Selection.h"


class Trigger{
	public:
		Trigger(bool isMuon, bool isElectron);
		~Trigger();
		void bookTriggers(bool isData);
		void checkAvail(int currentRunTrig, vector<Dataset*> datasets, unsigned int d, TTreeLoader* treeLoader, TRootEvent* event, bool verbose);
		int checkIfFired();


	private:
		bool muon;
		bool electron;
		bool trigged;
		bool redotrigmap;
	  std::vector<std::string> triggerList;
	  int currentRunTrig;
	  int previousRunTrig;
		string currentFilenameTrig;
		string previousFilenameTrig;
		int iFileTrig;
	  std::map<std::string,std::pair<int,bool> > triggermap;
};

#endif
