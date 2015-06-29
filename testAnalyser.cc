////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////



#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>

//used TopTreeAnalysis classes
#include "../TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/Selection.h"
#include "../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"


using namespace std;
using namespace reweight;
using namespace TopTree;


int main (int argc, char *argv[])
{
  
  clock_t start = clock();
  
  
  /////////////////////
  // Configuration
  /////////////////////

  string xmlFileName ="config/test.xml";

  if (argc > 1)
      xmlFileName = (string)argv[1];

  const char *xmlfile = xmlFileName.c_str();
  
  cout << " - Using config file " << xmlfile << endl;
  
  
  AnalysisEnvironment anaEnv;
  cout << " - Loading environment ..." << endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  
  cout << anaEnv.JetCollection << " " <<  anaEnv.METCollection << " "
    << anaEnv.ElectronCollection << " " << anaEnv.MuonCollection << " "
    << anaEnv.PrimaryVertexCollection << " " << anaEnv.GenJetCollection << " "
    << endl;
  
  float oldLuminosity = anaEnv.Luminosity;  // in 1/pb
  cout << "Luminosity from analysis environment is " << oldLuminosity << endl;
  
  
  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets(datasets, xmlfile);
  
  
  ////////////
  
  
  
  
  
  
  
  
  
  
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " s to run the program" << endl;
    
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
  
}
