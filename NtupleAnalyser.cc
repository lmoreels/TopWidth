#include <TStyle.h>
#include <TPaveText.h>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <array>
#include <TRandom3.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TGraph.h>
#include <TGraph2D.h>

// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
//#include "../macros/Style.C"

#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"

// user defined
#include "Tools/interface/ResolutionFunctions.h"
#include "Tools/interface/KinFitter.h"
#include "Tools/interface/EventReweighting.h"
#include "Tools/interface/Likelihood.h"
#include "Tools/interface/Likelihood2D.h"
#include "Tools/interface/LikelihoodMass.h"
#include "Tools/interface/SelectionTables.h"


using namespace std;
using namespace TopTree;


bool test = false;
bool testHistos = false;
bool testTTbarOnly = false;
bool skipData = false;
bool unblind = false;
bool doGenOnly = false;
bool makePlots = false;
bool makeControlPlots = false;
bool makeLikelihoodPlots = false;
bool calculateAverageMassAllMC = false;
bool calculateResolutionFunctions = false;
bool calculateAverageMass = false;
bool calculateFractions = false;
bool makeTGraphs = false;
bool useTTTemplates = false;
bool calculateLikelihood = true;
bool doPseudoExps = false;
bool useNewVar = true;

bool doLikeW = true;
bool doLikeM = false;
bool doLike2D = false;

bool doKinFit = true;
bool applyKinFitCut = true;
double kinFitCutValue = 15.;
double kinFitMinCutValue = 0.;

bool doMETCleaning = true;
bool applyLeptonSF = true;
bool applyPU = true;
bool applyBTagSF = true;

bool runTTbar = true;
bool runSTtW = true;
bool runSTt = true;
bool runOther = true;


bool rewHadTopOnly = false;
bool applyWidthSF = false;
double scaleWidth = 0.6;


bool runListWidths = false;
double listWidths[] = {0.2, 0.4, 0.5, 0.6, 0.8, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.};
int nWidths = sizeof(listWidths)/sizeof(listWidths[0]);
double thisWidth, thisMass;
double origWidth = 1., origMass = 172.5;

bool runSystematics = false;
bool runRateSystematics = false;
bool runSampleSystematics = false;
bool useOtherXml = false;
bool makeMassTemplates = false;
int nSystematics;
string thisSystematic;
string thisDataSetName = "";

string listRateSyst[] = {"nominal", "leptonIdSFup", "leptonIdSFdown", "leptonIsoSFup", "leptonIsoSFdown", "leptonTrigSFup", "leptonTrigSFdown", "leptonTrkSFup", "leptonTrkSFdown", "puSFup", "puSFdown", "btagSFup", "btagSFdown", "topPtReweighting", "lumiup", "lumidown", "renFac1002", "renFac1003", "renFac1004", "renFac1005", "renFac1007", "renFac1009","fragUp", "fragCentral", "fragDown", "fragPeterson", "fragSemiLepBrUp", "fragSemiLepBrDown"};
int nRateSystematics = sizeof(listRateSyst)/sizeof(listRateSyst[0]);
string listSampleSyst[] = {"tuneup", "tunedown", "isrup", "isrdown", "fsrup", "fsrdown", "hdampup", "hdampdown", "mpiERD", "qcdERD", "gluonMove", "gluonMoveERD", "mass169p5", "mass171p5", "mass173p5", "mass175p5", "herwig"};
string dataSetNameSyst[] = {"TT_tune_up", "TT_tune_down", "TT_isr_up", "TT_isr_down", "TT_fsr_up", "TT_fsr_down", "TT_hdamp_up", "TT_hdamp_down", "TT_erdOn", "TT_QCD_erdOn", "TT_gluon_move_erdOff", "TT_gluon_move_erdOn", "TT_mass169p5", "TT_mass171p5", "TT_mass173p5", "TT_mass175p5", "TT_herwigpp"};
int nSampleSystematics = sizeof(listSampleSyst)/sizeof(listSampleSyst[0]);
const char *xmlSyst = "config/topWidth_extra.xml";


TFile *fileWidths;


bool newTrees = true;
string systStr = "nominal";
pair<string,string> whichDate(string syst)
{
  if ( syst.find("nominal") != std::string::npos )
  {
    if (newTrees)
      return pair<string,string>("171121","171021");
    else
      return pair<string,string>("171015","171021");
  }
  else if ( syst.find("JESup") != std::string::npos )
  {
    newTrees = false;
    return pair<string,string>("171022","171021");
  }
  else if ( syst.find("JESdown") != std::string::npos )
  {
    newTrees = false;
    return pair<string,string>("171023","171021");
  }
  else if ( syst.find("JERup") != std::string::npos ) return pair<string,string>("171128","171126");
  else
  {
    cout << "WARNING: No valid systematic given! Will use nominal sample..." << endl;
    return pair<string,string>("171121","171021");
  }
}
pair<string,string> ntupleDate = whichDate(systStr);
string ntupleSystDate = "171123";
string ntupleSystDateExtraJES = "171124";
int verbose = 2;

string pathNtuples = "";
string pathNtuplesMC = "";
string pathNtuplesData = "";
string pathNtuplesSyst = "";
string pathOutput = "";
string outputDirLL = "LikelihoodTemplates/";
string inputDirLL = "";
string inputDateLL = "171128_2250/";  // 1D likelihood; no WM; redMlbMass, 30 < jet pT < 250, 0.2 -> 2.1, chi2 < 15, mt_alt > 240, m_lb < 240 + cuts
//string inputDateLL = "171128_1923/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171128_1721/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.6 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171128_1351/";  // 1D likelihood; no WM; new redTopMass, 30 < jet pT < 250, 0.6 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171127_1814/";  // 1D likelihood; no WM; new redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171127_1552/";  // 1D likelihood; no WM; old redTopMass, 35 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171127_1446/";  // 1D likelihood; no WM; old redTopMass, 40 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171126_1216/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171125_1529/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171125_1325/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171124_2050/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 200, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171124_1839/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240, Ht > 220 + cuts
//string inputDateLL = "171124_1449/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240, jet pT_aKF > 25 + cuts
//string inputDateLL = "171124_1154/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171123_1611/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171123_1432/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171122_2020/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171121_2216/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171121_1556/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171117_2002/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240, Ht > 220 + cuts
//string inputDateLL = "171114_1828/";  // 2D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171117_1235/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240, jets pT_aKF > 30 + cuts
//string inputDateLL = "171117_1623/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240, Ht > 220 + cuts
//string inputDateLL = "171116_1551/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240, Ht > 200 + cuts
//string inputDateLL = "171116_1033/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 250 + cuts
//string inputDateLL = "171112_1220/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171112_1125/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171111_2042/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171111_2012/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts aKF
//string inputDateLL = "171111_1937/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts aKF
//string inputDateLL = "171111_1844/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts aKF
//string inputDateLL = "171111_1751/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts bKF
//string inputDateLL = "171111_1705/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240 + cuts
//string inputDateLL = "171111_1542/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 15, mt_alt > 240
//string inputDateLL = "171111_1235/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 10, mt_alt > 240
//string inputDateLL = "171111_1141/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 5, mt_alt > 240
//string inputDateLL = "171109_1257/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 2, mt_alt > 240
//string inputDateLL = "171108_2124/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, 0.65 -> 1.4, chi2 < 2
//string inputDateLL = "171106_1136/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, ACDE, 0.65 -> 1.4, no btag SFs
//string inputDateLL = "171105_1314/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, ACDE, 0.65 -> 1.4, 4-5 jets
//string inputDateLL = "171105_2016/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, no min dR, ACDE, 0.65 -> 1.4
//string inputDateLL = "171104_2048/";  // 1D likelihood; no WM; old redTopMass, 30 < jet pT < 250, ACDE, 0.65 -> 1.4
//string inputDateLL = "171104_1749/";  // 1D likelihood; no WM; new redTopMass, 30 < jet pT < 250, ACDE
//string inputDateLL = "171104_1537/";  // 1D likelihood; no WM; new redTopMass, 30 < jet pT < 250, 1.65 < m_3/2 < 2.45
//string inputDateLL = "171104_1113/";  // 1D likelihood; no WM; 30 < jet pT < 250; 1.5 < NEW VAR < 2.6
//string inputDateLL = "171102_1053/";  // 2D likelihood; no WM; old definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF, ACDE
//string inputDateLL = "171031_1739/";  // 1D likelihood; no WM; old definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF, ACDE
//string inputDateLL = "171031_1519/";  // 1D likelihood; no WM; old definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF, C4
//string inputDateLL = "171031_1440/";  // 1D likelihood; no WM; old definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF, C3
//string inputDateLL = "171031_1100/";  // 1D likelihood; no WM; old definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF, C2
//string inputDateLL = "171031_0932/";  // 1D likelihood; no WM; old definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF, C1
//string inputDateLL = "171030_1850/";  // 1D likelihood; no WM; old definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF, dR cuts bKF
//string inputDateLL = "171030_1049/";  // 1D likelihood; no WM; old definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF
//string inputDateLL = "171030_1141/";  // 1D likelihood; no WM; new definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF, no chi2 cut
//string inputDateLL = "171028_1015/";  // 1D likelihood; no WM; new definition redTopMass, 30 < jet pT < 250, aKF, no dR cuts, no mass cuts before KF
//string inputDateLL = "171026_1943/";  // 1D likelihood; no WM; new definition redTopMass, 30 < jet pT < 250, aKF, no chi2 cut
//string inputDateLL = "171027_1033/";  // 1D likelihood; no WM; new definition redTopMass, 30 < jet pT < 250, aKF, no mass cuts before KF
//string inputDateLL = "171026_2244/";  // 1D likelihood; no WM; new definition redTopMass, 30 < jet pT < 250, aKF, no mass cuts
//string inputDateLL = "171026_1704/";  // 1D likelihood; no WM; new definition redTopMass, 30 < jet pT < 250, bKF
//string inputDateLL = "171026_1404/";  // 1D likelihood; no WM; new definition redTopMass, 30 < jet pT < 250
//string inputDateLL = "171009_0931/";  // 1D likelihood, new variable: m_3/2
//string inputDateLL = "171011_1932/";  // 2D likelihood, new definition redTopMass
//string inputDateLL = "171018_1034/";  // 1D likelihood, mt/<mt>_CM + (<mt>_CM - <mt>_all)/<mt>_CM
//string inputDateLL = "171016_1639/";  // 1D likelihood, new definition redTopMass, 30 < jet pT < 200
//string inputDateLL = "171013_1418/";  // 1D likelihood, new definition redTopMass, 30 < jet pT < 250
//string inputDateLL = "171012_1630/";  // 1D likelihood, new definition redTopMass
//string inputDateLL = "171002_1139/";  // reweight mass, range [171.5, 173.5] GeV
//string inputDateLL = "170929_1018/";  // reweight mass
//string inputDateLL = "170922_1045/";  // dR > 0.8
//string inputDateLL = "170918_1300/";  // diff masses
//string inputDateLL = "170829_1010/";  // use relativeSF instead of lumiWeight
//string inputDateLL = "170814_1625/";  // all SFs, all samples; widthSF had&lep; redtopmass [0.6, 1.4]; other widths
//string inputDateLL = "170802_1317/";  // all SFs, ttbar-only; widthSF had&lep; redtopmass [0.6, 1.4]
//string inputDateLL = "170802_1222/";  // all SFs, all samples; widthSF had&lep; redtopmass [0.6, 1.4]
//string inputDateLL = "170801_1451/";  // all SFs, ttbar-only; mlb cut 160; widthSF had-only; redtopmass [0.7, 1.3]
//string inputDateLL = "170801_1252/";  // all SFs, ttbar-only; mlb cut 160; widthSF had & lep; redtopmass [0.7, 1.3]
//string inputDateLL = "170731_1824/";  // all SFs, ttbar-only; no smoothing; mlb cut 160; widthSF had-only
//string inputDateLL = "170731_1823/";  // all SFs, ttbar-only; no smoothing; mlb cut 160; widthSF had & lep
//string inputDateLL = "170728_1418/";  // w/o W+1/2jets, all SFs, ttbar-only; no smoothing; widthSF had & lep
//string inputDateLL = "170728_1324/";  // w/o W+1/2jets, all SFs, ttbar-only; no smoothing; widthSF had-only
string whichTemplates()
{
  //if (rewHadTopOnly) return "170728_1324/";
  //else return "170802_1317/";
  if (useNewVar) return "171009_1346/";
  else return "171006_1130/";
}

bool isData = false;
bool isTTbar = false;
bool isTTsemilep = false;
bool isTTother = false;
bool isST = false;
bool isSTtW = false;
bool isSTt = false;
bool isOther = false;
bool isHerwig = false;

int nofHardSelected = 0;
int nofMETCleaned = 0;
int nofAfterDRmincut = 0;
int nofAfterDRmaxcut = 0;
int nofAfterLastCut = 0;
int nofTTsemilep = 0;
int nofTTdilep = 0;
int nofTThadr = 0;
int nofMatchedEvents = 0;
int nofHadrMatchedEvents = 0;
int nofHadrMatchedEventsAKF = 0;
int nofCorrectlyMatched = 0;
int nofNotCorrectlyMatched = 0;
int nofUnmatched = 0;
int nofUnmatchedTTsemilep = 0;
int nofUnmatchedTTother = 0;
int nofCorrectlyMatchedAKF = 0;
int nofNotCorrectlyMatchedAKF = 0;
int nofUnmatchedAKF = 0;
int nofUnmatchedTTsemilepAKF = 0;
int nofUnmatchedTTotherAKF = 0;
int nofCorrectlyMatchedAKFNoCut = 0;
int nofNotCorrectlyMatchedAKFNoCut = 0;
int nofUnmatchedAKFNoCut = 0;
int nofUnmatchedTTsemilepAKFNoCut = 0;
int nofUnmatchedTTotherAKFNoCut = 0;
int nofCM = 0, nofWM = 0, nofUM = 0;
int nofCM_TT = 0, nofWM_TT = 0, nofUM_TT = 0;
double nofCMl = 0., nofWMl = 0., nofUMl = 0.;
double nofCM_weighted = 0., nofWM_weighted = 0., nofUM_weighted = 0., nofUM_TTsemilep_weighted = 0., nofUM_TTother_weighted = 0., nofUM_other_weighted = 0.;
double nofCMout_weighted = 0., nofWMout_weighted = 0., nofUMout_weighted = 0., nofUMout_TTsemilep_weighted = 0., nofUMout_TTother_weighted = 0., nofUMout_other_weighted = 0.;
double nofCM_gate1 = 0., nofWM_gate1 = 0., nofUM_gate1 = 0.;
double nofCM_gate2 = 0., nofWM_gate2 = 0., nofUM_gate2 = 0.;
double nofBKF_weighted[15], nofAKF_weighted[15];

/// Lumi per data era
double lumi_runBCDEF = 19.67550334113;  // 1/fb
double lumi_runGH = 16.146177597883;  // 1/fb
double Luminosity = (lumi_runBCDEF + lumi_runGH)*1000;  // 1/pb

///  Working points for b tagging  // Updated 13/04/17, https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
double CSVv2Loose  = 0.5426;
double CSVv2Medium = 0.8484;
double CSVv2Tight  = 0.9535;

/// Average top mass
// TT_genp_match, TT_genj_match, TT_reco_match, TT_reco_wrongMatch_WM/UM, TT_reco_noMatch, TT_reco_wrongPerm, TT_reco, ST_t_top, ST_t_antitop, ST_tW_top, ST_tW_antitop, DYJets, WJets, data, Reco, All, MC, Reco, All, Samples
// also background in CM/WM/UM cats (unlike name suggests)
const int nofAveMasses = 16;
//  KF chi2 < 15
std::array<double, 14> aveTopMass = {171.908, 170.684, 169.412, 168.525, 168.482, 170.976, 169.012, 166.927, 168.547, 170.504, 171.287, 168.810, 169.042, 168.928};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240, 0.6->1.4 + cuts 28_1707
//std::array<double, 14> aveTopMass = {172.028, 173.264, 170.350, 170.994, 170.944, 174.429, 170.502, 168.937, 170.578, 173.522, 173.469, 170.274, 170.556, 170.419};   // Res 171025, 35 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 27_1543
//std::array<double, 14> aveTopMass = {172.139, 175.837, 171.244, 173.307, 173.252, 177.886, 171.866, 168.288, 172.296, 175.581, 176.613, 171.728, 171.937, 171.836};   // Res 171025, 40 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 27_1433
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.412, 168.525, 168.482, 170.976, 169.012, 166.927, 168.547, 170.504, 171.287, 168.810, 169.042, 168.928};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 25_1515
//std::array<double, 14> aveTopMass = {171.909, 170.808, 169.251, 169.240, 169.215, 170.358, 169.162, 167.877, 169.679, 171.486, 171.689, 168.899, 169.208, 169.057};   // Res 171025, 30 < jet pT, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 25_1103
//std::array<double, 14> aveTopMass = {171.897, 170.451, 168.912, 168.779, 168.725, 172.158, 168.784, 167.303, 169.138, 170.777, 171.562, 168.571, 168.828, 168.702};   // Res 171025, 30 < jet pT < 200, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 24_2038
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.497, 170.062, 170.026, 171.955, 169.628, 168.757, 171.222, 172.540, 173.200, 169.496, 169.693, 169.598};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 24_1825
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.110, 169.035, 168.993, 171.557, 169.004, 167.699, 169.615, 171.119, 171.758, 168.753, 169.051, 168.905};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 22_2217
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.433, 171.782, 171.718, 176.760, 170.159, 171.264, 174.114, 173.200, 173.957, 170.032, 170.239, 170.138};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 22_2009
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.110, 169.035, 168.993, 171.574, 169.004, 167.690, 169.612, 171.141, 171.754, 168.753, 169.051, 168.905};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 21_1545
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.497, 170.062, 170.026, 171.974, 169.628, 168.746, 171.226, 172.564, 173.196, 169.496, 169.694, 169.598};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 17_1612
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.591, 172.721, 172.658, 177.768, 170.545, 172.230, 175.428, 173.910, 174.732, 170.482, 170.634, 170.560};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 17_1209
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.302, 169.554, 169.516, 171.756, 169.318, 168.513, 170.241, 171.810, 172.579, 169.117, 169.376, 169.249};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 16_1541
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.191, 169.293, 169.253, 172.063, 169.155, 167.466, 170.560, 171.449, 172.168, 168.952, 169.207, 169.082};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 250 + cuts 16_1021
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.110, 169.035, 168.993, 171.574, 169.004, 167.690, 169.612, 171.141, 171.754, 168.753, 169.051, 168.905};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 12_1201
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.456, 169.446, 169.410, 171.704, 169.377, 167.763, 170.071, 171.770, 173.175, 169.193, 169.435, 169.318};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 12_1115
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.600, 171.573, 171.502, 177.151, 170.148, 171.977, 173.367, 173.608, 174.462, 170.056, 170.236, 170.148};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 11_2034
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.603, 171.724, 171.644, 177.944, 170.186, 172.027, 173.617, 173.874, 174.423, 170.093, 170.278, 170.188};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 11_2004
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.316, 172.304, 172.199, 180.262, 170.298, 173.642, 175.893, 173.529, 174.191, 170.200, 170.396, 170.300};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 11_1925
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.469, 171.696, 171.613, 178.093, 170.113, 172.493, 174.230, 173.833, 173.942, 170.032, 170.206, 170.121};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 11_1828
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.170, 172.252, 172.198, 176.559, 170.325, 174.415, 175.673, 173.721, 174.571, 170.237, 170.426, 170.334};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 11_1735
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.471, 171.522, 171.458, 176.973, 170.062, 172.493, 173.836, 173.869, 174.530, 169.935, 170.151, 170.046};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240 + cuts 11_1653
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.170, 172.466, 172.379, 178.599, 170.419, 174.362, 175.983, 173.429, 174.293, 170.389, 170.523, 170.457};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240
//  KF chi2 < 10
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.267, 172.409, 172.328, 178.314, 170.327, 173.790, 174.784, 173.574, 174.104, 170.237, 170.424, 170.333};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240
//  KF chi2 < 5
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.439, 172.192, 172.114, 177.828, 170.204, 173.751, 174.407, 174.451, 173.826, 170.153, 170.305, 170.231};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240
//  KF chi2 < 2
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.607, 171.837, 171.761, 176.490, 170.117, 171.209, 172.600, 173.741, 174.812, 170.194, 170.210, 170.202};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245, mt_alt > 240
//std::array<double, 14> aveTopMass = {171.908, 170.684, 168.512, 169.145, 168.780, 169.851, 168.782, 169.103, 169.446, 170.999, 171.637, 168.373, 168.843, 168.609};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 100 < mt < 245
//  KF chi2 < 5
//std::array<double, 14> aveTopMass = {171.900, 170.981, 169.553, 173.018, 172.479, 178.290, 170.653, 173.835, 175.415, 175.097, 174.673, 170.580, 170.760, 170.672};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, ACDE, 4-5 jets
//std::array<double, 14> aveTopMass = {171.828, 170.636, 169.365, 170.133, 170.047, 177.062, 169.491, 171.657, 172.082, 173.465, 173.326, 169.442, 169.584, 169.515};   // Res 171025, 30 < jet pT < 250, ACDE
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.439, 172.192, 172.114, 177.828, 170.204, 173.751, 174.407, 174.451, 173.826, 170.153, 170.305, 170.231};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, ACDE
//std::array<double, 14> aveTopMass = {171.908, 170.684, 168.100, 167.805, 168.344, 166.350, 167.960, 169.410, 168.439, 168.042, 167.928, 167.667, 167.965, 167.817};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 1.5 < m_3/2 < 2.6
//std::array<double, 14> aveTopMass = {171.908, 170.684, 168.058, 166.747, 167.211, 165.554, 167.427, 167.523, 167.474, 167.345, 167.913, 167.113, 167.431, 167.273};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, 1.5 < NEW VAR < 2.6
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.439, 172.192, 172.114, 177.828, 170.204, 173.751, 174.407, 174.451, 173.826, 170.153, 170.305, 170.231};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, mass cuts aKF, ACDE
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.437, 172.243, 172.171, 177.223, 170.210, 173.759, 175.264, 174.091, 173.749, 170.232, 170.309, 170.271};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, mass cuts aKF, C4
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.493, 172.102, 172.067, 174.055, 170.209, 174.537, 173.635, 174.719, 174.755, 170.185, 170.306, 170.247};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, mass cuts aKF, C3
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.217, 169.463, 169.441, 168.786, 169.268, 165.849, 164.705, 170.667, 169.883, 169.001, 169.273, 169.142};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, mass cuts aKF, C2
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.435, 171.977, 171.942, 173.420, 170.122, 174.510, 174.045, 173.970, 174.572, 170.224, 170.213, 170.218};   // Res 171025, 30 < jet pT < 250, dR jets > 0.6, mass cuts aKF, C1
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.058, 168.896, 168.902, 157.660, 168.996, 164.330, 164.642, 169.500, 169.494, 168.958, 168.989, 168.974};  // Res 171025, 30 < jet pT < 250, dR jets > 0.6, mass cuts aKF, dR cuts bKF
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.485, 172.310, 172.235, 177.941, 170.265, 173.727, 174.783, 174.617, 173.917, 170.186, 170.369, 170.279};  // Res 171025, 30 < jet pT < 250, dR jets > 0.6, mass cuts aKF
//std::array<double, 14> aveTopMass = {171.908, 170.684, 169.073, 168.860, 168.867, 156.244, 168.995, 164.122, 164.603, 169.416, 169.499, 169.031, 168.987, 169.008};  // Res 171025, 30 < jet pT < 250, dR jets > 0.6
//std::array<double, 14> aveTopMass = {172.009, 170.380, 167.449, 171.916, 171.340, 173.257, 169.846, 172.896, 173.663, 172.385, 173.262, 169.387, 169.933, 169.658};  // 30 < jet pT < 200
//std::array<double, 14> aveTopMass = {172.019, 170.603, 167.778, 202.264, 201.573, 203.810, 185.976, 256.266, 255.607, 220.936, 221.428, 188.305, 188.393, 188.378};  // Res 170915, dR(jets) > 0.8
//std::array<double, 14> aveTopMass = {171.828, 170.662, 167.968, 198.481, 198.418, 198.640, 183.102, 252.884, 250.668, 230.103, 227.942, 185.725, 185.985, 185.942};  // Res 170915
//std::array<double, 14> aveTopMass = {171.828, 170.662, 167.841, 197.327, 196.972, 198.197, 182.357, 248.361, 247.185, 227.921, 226.553, 184.765, 185.085, 185.033};  // Res 170912
//std::array<double, 14> aveTopMass = {171.833, 169.809, 167.636, 197.975, 197.718, 198.582, 182.317, 252.174, 249.964, 229.383, 227.814, 184.794, 185.096, 185.046};  // with SFs
//std::array<double, 14> aveTopMass = {171.826, 169.746, 167.511, 197.053, 196.687, 197.911, 181.895, 249.468, 247.437, 227.530, 226.099, 184.794, 184.594, 184.624};  // Res 170608 Single Gaus
//std::array<double, 14> aveTopMass = {171.826, 169.746, 167.556, 197.087, 196.662, 198.143, 182.150, 249.229, 246.893, 226.933, 225.681, 185.024, 184.880, 184.902};  // no DYJets, no WJets // Res 170608

/// # events after kin fitter
//  KF chi2 < 5
int nEventsAKF[] = {100949, 547108, 8927, 8847, 7051, 4283, 3, 16, 55, 276, 1, 18, 123, 1415};

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlotCP;
map<string,TH1D*> histo1DLike;
map<string,TGraph*> graph;

map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;

//map<string,TNtuple*> ntuple;

vector < Dataset* > datasets, datasetsMSP, datasetsTemp, datasetsSyst;
vector<string> dataSetNames;
vector<int> includeDataSets;

ofstream txtMassGenPMatched, txtMassGenJMatched, txtMassRecoCM, txtMassRecoWMUM, txtMassRecoUM, txtMassRecoWM, txtMassReco;
double sumTopMass, sumMlb, sumEvents, sumCMTopMass, sumCMMlb, sumCMEvents;
Double_t aveMlbMass = 99.9389, aveMlbMassCM = 97.9464;

/// Function prototypes
struct HighestPt
{
  bool operator()( TLorentzVector j1, TLorentzVector j2 ) const
  {
    return j1.Pt() > j2.Pt() ;
  }
};

struct greater
{
  template<class T>
  bool operator()(T const &a, T const &b) const { return a > b; }
};


string ConvertDoubleToString(double Number);
string DotReplace(double var);
string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();
bool fexists(const char *filename);
void InitSetUp();
void GetMetaData(TTree* tree, bool isData);
void InitTree(TTree* tree, bool isData);
void InitMSPlots();
void InitHisto1D();
void InitHisto2D();
void InitHisto1DMatch();
void InitHisto2DMatch();
void InitLikelihoodPlots();
void TruthMatching(vector<TLorentzVector> partons, vector<TLorentzVector> selectedJets, pair<unsigned int, unsigned int> *MCPermutation);
void ClearMetaData();
void ClearLeaves();
void ClearTLVs();
void ClearMatching();
void ClearVars();
void ClearObjects();
void FillControlPlots(vector<Dataset *> datasets, int d);
void FillControlPlots(vector<Dataset *> datasets, int d, string suffix);
void FillMatchingPlots();
void FillKinFitPlots(bool doneKinFit);
void FillCatsPlots(string catSuffix);
void FillMSPlots(int d, bool doneKinFit);
void FillLikelihoodPlots();
void WriteLikelihoodPlots();
long GetNEvents(TTree* fChain, string var, bool isData);
long GetNEvents(TTree* fChain, string var, unsigned int index, bool isData);
void GetEraFraction(double* fractions);
void SetUpSelectionTable();
void FillSelectionTable(int d, string dataSetName);
void FillSelectionTable2(int d, string dataSetName);
void WriteSelectionTable();
void CheckSystematics(vector<int> vJER, vector<int> vJES, vector<int> vPU);
void PrintKFDebug(int ievt);



// Declaration of leaf types
Int_t           run_num;
Long64_t        evt_num;
Int_t           lumi_num;
Int_t           nvtx;
Int_t           npu;
Double_t        rho;
Bool_t          isTrigged;
Bool_t          hasExactly4Jets;
Bool_t          hasJetLeptonCleaning;
Bool_t          hasErasedBadOrCloneMuon;
Bool_t          passedMETFilter = true;
Bool_t          isDataRunB;
Bool_t          isDataRunC;
Bool_t          isDataRunD;
Bool_t          isDataRunE;
Bool_t          isDataRunF;
Bool_t          isDataRunG;
Bool_t          isDataRunH;
Int_t           nMuons;
Double_t        muon_charge[1];   //[nMuons]
Double_t        muon_pt[1];   //[nMuons]
Double_t        muon_phi[1];   //[nMuons]
Double_t        muon_eta[1];   //[nMuons]
Double_t        muon_E[1];   //[nMuons]
Double_t        muon_M[1];   //[nMuons]
Double_t        muon_d0[1];   //[nMuons]
Double_t        muon_chargedHadronIso[1];   //[nMuons]
Double_t        muon_neutralHadronIso[1];   //[nMuons]
Double_t        muon_photonIso[1];   //[nMuons]
Double_t        muon_puChargedHadronIso[1];   //[nMuons]
Double_t        muon_relIso[1];   //[nMuons]
Double_t        muon_pfIso[1];   //[nMuons]
Int_t           nJets;
Int_t           jet_nConstituents[20];   //[nJets]
Int_t           jet_nChConstituents[20];   //[nJets]
Double_t        jet_charge[20];   //[nJets]
Double_t        jet_pt[20];   //[nJets]
Double_t        jet_phi[20];   //[nJets]
Double_t        jet_eta[20];   //[nJets]
Double_t        jet_E[20];   //[nJets]
Double_t        jet_M[20];   //[nJets]
Double_t        jet_bdiscr[20];   //[nJets]
Double_t        jet_hadronFlavour[20];   //[nJets]
Double_t        met_px;
Double_t        met_py;
Double_t        met_pt;
Double_t        met_phi;
Double_t        met_eta;
Double_t        met_Et;
Double_t        met_E;
Double_t        met_corr_px;
Double_t        met_corr_py;
Double_t        met_corr_pt;
Double_t        met_corr_phi;
Double_t        met_corr_eta;
Double_t        met_corr_Et;
Double_t        met_corr_E;
Int_t           nGenJets;
Int_t           genJet_pdgId[40];   //[nGenJets]
Double_t        genJet_charge[40];   //[nGenJets]
Double_t        genJet_pt[40];   //[nGenJets]
Double_t        genJet_phi[40];   //[nGenJets]
Double_t        genJet_eta[40];   //[nGenJets]
Double_t        genJet_E[40];   //[nGenJets]
Double_t        genJet_M[40];   //[nGenJets]
Int_t           nMCParticles;
Int_t           mc_status[200];   //[nMCParticles]
Int_t           mc_pdgId[200];   //[nMCParticles]
Int_t           mc_mother[200];   //[nMCParticles]
Int_t           mc_granny[200];   //[nMCParticles]
Double_t        mc_pt[200];   //[nMCParticles]
Double_t        mc_phi[200];   //[nMCParticles]
Double_t        mc_eta[200];   //[nMCParticles]
Double_t        mc_E[200];   //[nMCParticles]
Double_t        mc_M[200];   //[nMCParticles]
Bool_t          mc_isLastCopy[200];   //[nMCParticles]
Bool_t          mc_isPromptFinalState[200];   //[nMCParticles]
Bool_t          mc_isHardProcess[200];   //[nMCParticles]
Bool_t          mc_fromHardProcessFinalState[200];   //[nMCParticles]
Bool_t          hasGenTop;
Bool_t          hasGenAntiTop;
Double_t        baseweight;
Double_t        weight1001;
Double_t        weight1002;
Double_t        weight1003;
Double_t        weight1004;
Double_t        weight1005;
Double_t        weight1007;
Double_t        weight1009;
Double_t        upFragWeight;
Double_t        centralFragWeight;
Double_t        downFragWeight;
Double_t        petersonFragWeight;
Double_t        semilepbrUp;
Double_t        semilepbrDown;
Double_t        pdfWeights[100];
Double_t        pdfAlphaSUp;
Double_t        pdfAlphaSDown;
Double_t        btagSF;
Double_t        btagSF_up;
Double_t        btagSF_down;
Double_t        puSF;
Double_t        puSF_up;
Double_t        puSF_down;
Double_t        muonIdSF_BCDEF[1];   //[nMuons]
Double_t        muonIdSF_GH[1];   //[nMuons]
Double_t        muonIdSF_up_BCDEF[1];   //[nMuons]
Double_t        muonIdSF_up_GH[1];   //[nMuons]
Double_t        muonIdSF_down_BCDEF[1];   //[nMuons]
Double_t        muonIdSF_down_GH[1];   //[nMuons]
Double_t        muonIsoSF_BCDEF[1];   //[nMuons]
Double_t        muonIsoSF_GH[1];   //[nMuons]
Double_t        muonIsoSF_up_BCDEF[1];   //[nMuons]
Double_t        muonIsoSF_up_GH[1];   //[nMuons]
Double_t        muonIsoSF_down_BCDEF[1];   //[nMuons]
Double_t        muonIsoSF_down_GH[1];   //[nMuons]
Double_t        muonTrigSF_BCDEF[1];   //[nMuons]
Double_t        muonTrigSF_GH[1];   //[nMuons]
Double_t        muonTrigSF_up_BCDEF[1];   //[nMuons]
Double_t        muonTrigSF_up_GH[1];   //[nMuons]
Double_t        muonTrigSF_down_BCDEF[1];   //[nMuons]
Double_t        muonTrigSF_down_GH[1];   //[nMuons]
Double_t        muonTrackSF_eta[1];   //[nMuons]
Double_t        muonTrackSF_aeta[1];   //[nMuons]
Double_t        muonTrackSF_nPV[1];   //[nMuons]

Long64_t        nEvents;
Long64_t        nEventsSel;
Int_t           cutFlow[10];
Int_t           cutFlow2[10];
Double_t        cutFlowWeighted[10];
Double_t        cutFlow2Weighted[10];
Int_t           appliedJER;
Int_t           appliedJES;
Int_t           appliedPU;
Long64_t        nofEventsRunB;
Long64_t        nofEventsRunCD;
Long64_t        nofEventsRunEF;
Long64_t        nofEventsRunG;
Long64_t        nofEventsRunH;
Long64_t        nofSelEventsRunB;
Long64_t        nofSelEventsRunCD;
Long64_t        nofSelEventsRunEF;
Long64_t        nofSelEventsRunG;
Long64_t        nofSelEventsRunH;
Long64_t        nofEventsWithGenTop;
Long64_t        nofEventsWithGenAntiTop;
Long64_t        nofTTEventsWithoutAGenTop;
Double_t        sumWeight1001;
Double_t        sumWeight1002;
Double_t        sumWeight1003;
Double_t        sumWeight1004;
Double_t        sumWeight1005;
Double_t        sumWeight1007;
Double_t        sumWeight1009;
Double_t        sumUpFragWeight;
Double_t        sumCentralFragWeight;
Double_t        sumDownFragWeight;
Double_t        sumPetersonFragWeight;
Double_t        sumSemilepbrUp;
Double_t        sumSemilepbrDown;
Double_t        sumPdfWeights[100];
Double_t        sumPdfAlphaSUp;
Double_t        sumPdfAlphaSDown;

// List of branches
TBranch        *b_run_num;   //!
TBranch        *b_evt_num;   //!
TBranch        *b_lumi_num;   //!
TBranch        *b_nvtx;   //!
TBranch        *b_npu;   //!
TBranch        *b_rho;   //!
TBranch        *b_isTrigged;   //!
TBranch        *b_hasExactly4Jets;   //!
TBranch        *b_hasJetLeptonCleaning;   //!
TBranch        *b_hasErasedBadOrCloneMuon;   //!
//TBranch        *b_passedMETFilter;   //!
TBranch        *b_isDataRunB;   //!
TBranch        *b_isDataRunC;   //!
TBranch        *b_isDataRunD;   //!
TBranch        *b_isDataRunE;   //!
TBranch        *b_isDataRunF;   //!
TBranch        *b_isDataRunG;   //!
TBranch        *b_isDataRunH;   //!
TBranch        *b_nMuons;   //!
TBranch        *b_muon_charge;   //!
TBranch        *b_muon_pt;   //!
TBranch        *b_muon_phi;   //!
TBranch        *b_muon_eta;   //!
TBranch        *b_muon_E;   //!
TBranch        *b_muon_M;   //!
TBranch        *b_muon_d0;   //!
TBranch        *b_muon_chargedHadronIso;   //!
TBranch        *b_muon_neutralHadronIso;   //!
TBranch        *b_muon_photonIso;   //!
TBranch        *b_muon_puChargedHadronIso;   //!
TBranch        *b_muon_relIso;   //!
TBranch        *b_muon_pfIso;   //!
TBranch        *b_nJets;   //!
TBranch        *b_jet_nConstituents;   //!
TBranch        *b_jet_nChConstituents;   //!
TBranch        *b_jet_charge;   //!
TBranch        *b_jet_pt;   //!
TBranch        *b_jet_phi;   //!
TBranch        *b_jet_eta;   //!
TBranch        *b_jet_E;   //!
TBranch        *b_jet_M;   //!
TBranch        *b_jet_bdiscr;   //!
TBranch        *b_jet_hadronFlavour;   //!
TBranch        *b_met_px;   //!
TBranch        *b_met_py;   //!
TBranch        *b_met_pt;   //!
TBranch        *b_met_phi;   //!
TBranch        *b_met_eta;   //!
TBranch        *b_met_Et;   //!
TBranch        *b_met_E;   //!
TBranch        *b_met_corr_px;   //!
TBranch        *b_met_corr_py;   //!
TBranch        *b_met_corr_pt;   //!
TBranch        *b_met_corr_phi;   //!
TBranch        *b_met_corr_eta;   //!
TBranch        *b_met_corr_Et;   //!
TBranch        *b_met_corr_E;   //!
TBranch        *b_nGenJets;   //!
TBranch        *b_genJet_pdgId;   //!
TBranch        *b_genJet_charge;   //!
TBranch        *b_genJet_pt;   //!
TBranch        *b_genJet_phi;   //!
TBranch        *b_genJet_eta;   //!
TBranch        *b_genJet_E;   //!
TBranch        *b_genJet_M;   //!
TBranch        *b_nMCParticles;   //!
TBranch        *b_mc_status;   //!
TBranch        *b_mc_pdgId;   //!
TBranch        *b_mc_mother;   //!
TBranch        *b_mc_granny;   //!
TBranch        *b_mc_pt;   //!
TBranch        *b_mc_phi;   //!
TBranch        *b_mc_eta;   //!
TBranch        *b_mc_E;   //!
TBranch        *b_mc_M;   //!
TBranch        *b_mc_isLastCopy;   //!
TBranch        *b_mc_isPromptFinalState;   //!
TBranch        *b_mc_isHardProcess;   //!
TBranch        *b_mc_fromHardProcessFinalState;   //!
TBranch        *b_hasGenTop;   //!
TBranch        *b_hasGenAntiTop;   //!
TBranch        *b_weight1001;   //!
TBranch        *b_weight1002;   //!
TBranch        *b_weight1003;   //!
TBranch        *b_weight1004;   //!
TBranch        *b_weight1005;   //!
TBranch        *b_weight1007;   //!
TBranch        *b_weight1009;   //!
TBranch        *b_upFragWeight;   //!
TBranch        *b_centralFragWeight;   //!
TBranch        *b_downFragWeight;   //!
TBranch        *b_petersonFragWeight;   //!
TBranch        *b_semilepbrUp;   //!
TBranch        *b_semilepbrDown;   //!
TBranch        *b_pdfWeights;   //!
TBranch        *b_pdfAlphaSUp;   //!
TBranch        *b_pdfAlphaSDown;   //!
TBranch        *b_btagSF;   //!
TBranch        *b_btagSF_up;   //!
TBranch        *b_btagSF_down;   //!
TBranch        *b_puSF;   //!
TBranch        *b_puSF_up;   //!
TBranch        *b_puSF_down;   //!
TBranch        *b_muonIdSF_BCDEF;   //!
TBranch        *b_muonIdSF_GH;   //!
TBranch        *b_muonIdSF_up_BCDEF;   //!
TBranch        *b_muonIdSF_up_GH;   //!
TBranch        *b_muonIdSF_down_BCDEF;   //!
TBranch        *b_muonIdSF_down_GH;   //!
TBranch        *b_muonIsoSF_BCDEF;   //!
TBranch        *b_muonIsoSF_GH;   //!
TBranch        *b_muonIsoSF_up_BCDEF;   //!
TBranch        *b_muonIsoSF_up_GH;   //!
TBranch        *b_muonIsoSF_down_BCDEF;   //!
TBranch        *b_muonIsoSF_down_GH;   //!
TBranch        *b_muonTrigSF_BCDEF;   //!
TBranch        *b_muonTrigSF_GH;   //!
TBranch        *b_muonTrigSF_up_BCDEF;   //!
TBranch        *b_muonTrigSF_up_GH;   //!
TBranch        *b_muonTrigSF_down_BCDEF;   //!
TBranch        *b_muonTrigSF_down_GH;   //!
TBranch        *b_muonTrackSF_eta;   //!
TBranch        *b_muonTrackSF_aeta;   //!
TBranch        *b_muonTrackSF_nPV;   //!

TBranch        *b_nEvents;   //!
TBranch        *b_nEventsSel;   //!
TBranch        *b_cutFlow;   //!
TBranch        *b_cutFlow2;   //!
TBranch        *b_cutFlowWeighted;   //!
TBranch        *b_cutFlow2Weighted;   //!
TBranch        *b_appliedJER;   //!
TBranch        *b_appliedJES;   //!
TBranch        *b_appliedPU;   //!
TBranch        *b_nofEventsRunB;   //!
TBranch        *b_nofEventsRunCD;   //!
TBranch        *b_nofEventsRunEF;   //!
TBranch        *b_nofEventsRunG;   //!
TBranch        *b_nofEventsRunH;   //!
TBranch        *b_nofSelEventsRunB;   //!
TBranch        *b_nofSelEventsRunCD;   //!
TBranch        *b_nofSelEventsRunEF;   //!
TBranch        *b_nofSelEventsRunG;   //!
TBranch        *b_nofSelEventsRunH;   //!
TBranch        *b_nofEventsWithGenTop;   //!
TBranch        *b_nofEventsWithGenAntiTop;   //!
TBranch        *b_nofTTEventsWithoutAGenTop;   //!
TBranch        *b_sumWeight1001;   //!
TBranch        *b_sumWeight1002;   //!
TBranch        *b_sumWeight1003;   //!
TBranch        *b_sumWeight1004;   //!
TBranch        *b_sumWeight1005;   //!
TBranch        *b_sumWeight1007;   //!
TBranch        *b_sumWeight1009;   //!
TBranch        *b_sumUpFragWeight;   //!
TBranch        *b_sumCentralFragWeight;   //!
TBranch        *b_sumDownFragWeight;   //!
TBranch        *b_sumPetersonFragWeight;   //!
TBranch        *b_sumSemilepbrUp;   //!
TBranch        *b_sumSemilepbrDown;   //!
TBranch        *b_sumPdfWeights;   //!
TBranch        *b_sumPdfAlphaSUp;   //!
TBranch        *b_sumPdfAlphaSDown;   //!


long nEventsDataSet;
double xSection;
double lumiWeight, scaleFactor, widthSF, relativeSF, eqLumi_TT;
double thisLeptonSF, thisLeptonIdSF, thisLeptonIsoSF, thisLeptonTrigSF;
double renFacSumNom, renFacSum1002, renFacSum1003, renFacSum1004, renFacSum1005, renFacSum1007, renFacSum1009;
double fragUpSum, fragCentralSum, fragDownSum, fragPetersonSum, semilepbrUpSum, semilepbrDownSum;
double topPtRewSF, topPtSF, antiTopPtSF;
bool foundTop22, foundAntiTop22;
bool foundTop62, foundAntiTop62;
bool foundLastCopyTop, foundLastCopyAntitop;
vector<unsigned int> bJetId;
double bdiscrTop, bdiscrTop2, tempbdiscr;
int labelB1, labelB2;
int labelsReco[4];
double massHadTopQ, massLepTopQ;

string catSuffix = "";
string catSuffixList[] = {"_CM", "_WM", "_UM"};
bool isCM, isWM, isUM, isUM_TTsemilep, isUM_TTother;
unsigned int dMSPmax;


/// Define TLVs
TLorentzVector muon, jet, mcpart;
TLorentzVector WCandidate;
vector<TLorentzVector> selectedLepton;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> selectedBJets;
vector<TLorentzVector> selectedJetsAKF;
vector<TLorentzVector> selectedJetsKFcorrected;
vector<TLorentzVector> mcParticles;
vector<TLorentzVector> partons;
vector<TLorentzVector> partonsMatched;
vector<TLorentzVector> jetsMatched;
vector<TLorentzVector*> selJets;

/// Matching
int pdgID_top = 6; //top quark

bool doMatching = true;
bool all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
bool all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
bool hadronicTopJetsMatched = false;
bool hadronicTopJetsMatched_MCdef_ = false;
pair<unsigned int, unsigned int> MCPermutation[4] = {pair<unsigned int,unsigned int>(9999,9999)};
int topQuark = -9999, antiTopQuark = -9999;
int genmuon = -9999;
bool muonmatched = false;
bool foundMuPlus = false, foundMuMinus = false;
bool muPlusFromTop = false, muMinusFromTop = false;
vector<unsigned int> partonId;
vector<unsigned int> foundB;
vector<unsigned int> foundLight;
unsigned int nPartons;

/// KinFitter
bool doneKinFit = false;
KinFitter *kf;
TKinFitter* kFitter;
bool addWMassKF = true;
bool addEqMassKF = false;
int kFitVerbosity = 0;
double kFitChi2 = 99.;
int nofAcceptedKFit = 0;
double nofAcceptedKFitWeighted = 0.;


/// Likelihood
//Double_t aveTopMassLL = aveTopMass[2];  // only CM events
Double_t aveTopMassLL = 169.073; //175.984; //168.058; //167.647; // 100 < mt < 245  // 168.879; // all MC events with 100 < mt < 250  //220.651;  // all MC events
Double_t aveTopMassCM = aveTopMass[2];
Double_t maxRedTopMass = 0., minRedTopMass = 9999.;
Double_t minTopMass = 100., maxTopMass = 245.;
Double_t minCutRedTopMass = 0.65, maxCutRedTopMass = 1.4;
int nWidthsLike = 0, nMassesLike = 0;
vector<double> widthsLike, massesLike;
vector<double> loglike_per_evt;
double redTopMassArray[] = {0.65, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1., 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.4};
int nLikeMasses = sizeof(redTopMassArray)/sizeof(redTopMassArray[0]) - 1;
Likelihood *likeW;  // like1D
LikelihoodMass *likeM;  // likeMass
Likelihood2D *like2D;  // like2D


ofstream txtDebugTopMass;

/// Pseudo experiments
const int nPseudoExps = 400;
int nPsExps = nPseudoExps;
TRandom3 random3;
double toyValues[nPseudoExps]; // = random3.Uniform(0,1);
double toyMax;
int nEvtsInPseudoExp[nPseudoExps][15] = {0}, nDataEvts;
double nEvtsInPseudoExpW[nPseudoExps][15] = {0.};


/// Variables
double M3, Ht, min_Mlb, dRLepB;
double M3_aKF, Ht_aKF;
double reco_W_mass_bKF, reco_top_mass_bKF, reco_top_mass_alt_bKF, reco_top_mass_rescaled_bKF, reco_W_pt_bKF, reco_top_pt_bKF, reco_top_pt_alt_bKF, reco_mlb_bKF, reco_dRLepB_lep_bKF, reco_dRLepB_had_bKF, reco_dRlight_bKF, reco_dRblight_qsum_bKF, reco_dRblight_min_bKF, reco_dRblight_max_bKF, reco_dRbW_bKF, reco_dRbb_bKF, reco_dPhi_bb_bKF, reco_dPhi_light_bKF, reco_dPhi_bW_bKF, reco_ttbar_mass_bKF, redTopMass_old_bKF, redTopMass_new_bKF, redTopMass_bKF, redMlbMass_bKF, reco_mbjj_div_mjj_bKF, reco_mTW_bKF;
double reco_W_mass_aKF, reco_top_mass_aKF, reco_top_mass_alt_aKF, reco_W_pt_aKF, reco_top_pt_aKF, reco_top_pt_alt_aKF, reco_mlb_aKF, reco_dRLepB_lep_aKF, reco_dRLepB_had_aKF, reco_dRlight_aKF, reco_dRblight_qsum_aKF, reco_dRblight_min_aKF, reco_dRblight_max_aKF, reco_dRbW_aKF, reco_dRbb_aKF, reco_dPhi_bb_aKF, reco_dPhi_light_aKF, reco_dPhi_bW_aKF, reco_ttbar_mass_aKF, redTopMass_old, redTopMass_new, redTopMass, redMlbMass, reco_mbjj_div_mjj, reco_new_var, reco_mTW_aKF;
double tempDR;

double matched_W_mass_q, matched_top_mass_q;
double matched_W_mass_j, matched_top_mass_j, matched_top_mass_j_akF;
double matched_mlb_corr, matched_ttbarMass_corr, matched_dR_lep_b_corr;
double matched_mlb_wrong, matched_ttbarMass_wrong, matched_dR_lep_b_wrong;


// JEC
string pathCalJEC = "../TopTreeAnalysisBase/Calibrations/JECFiles/";
vector<JetCorrectorParameters> vCorrParam;
JetCorrectionUncertainty *jecUnc;
void InitJEC(bool isData, string dataSetName);


/// Meta
string strSyst = "";
double eqLumi;
vector<int> vJER, vJES, vPU;
SelectionTables *selTab, *selTabKF;

bool CharSearch( char str[], char substr[] )
{
  char *output = NULL;
  output = strstr(str, substr);
  if (output) return true;
  else return false;
}

string dateString = "";

int main(int argc, char* argv[])
{
  dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  cout << "* The following corrections are applied:    *" << endl;
  cout << "*   - Jet Energy Corrections                *" << endl;
  if ( systStr.find("JECdown") != std::string::npos )
    cout << "*   - Jet Energy Scale: scale down          *" << endl;
  else if ( systStr.find("JECup") != std::string::npos )
    cout << "*   - Jet Energy Scale: scale up            *" << endl;
  cout << "*   - Jet Energy Resolution: ";
  if ( systStr.find("JERdown") != std::string::npos )    cout << "scale down     *" << endl;
  else if ( systStr.find("JERup") != std::string::npos ) cout << "scale up       *" << endl;
  else                                                   cout << "nominal        *" << endl;
  cout << "*   - Jet/lepton Cleaning                   *" << endl;
  if (doMETCleaning) cout << "*   - MET Cleaning                          *" << endl;
  cout << "*********************************************" << endl;
  cout << "* The following scale factors are applied:  *" << endl;
  if (applyLeptonSF) cout << "*   - Lepton scale factors: nominal         *" << endl;
  if (applyPU)       cout << "*   - Pile up: nominal                      *" << endl;
  if (applyBTagSF)   cout << "*   - B tag scale factors: nominal          *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  
  string channel = "mu";
  //if ( argc == 1 ) channel = "mu";
  //else if ( CharSearch(argv[1], "mu") || CharSearch(argv[1], "Mu") || CharSearch(argv[1], "MU") || CharSearch(argv[1], "m") ) channel = "mu";
  //else if ( CharSearch(argv[1], "el") || CharSearch(argv[1], "El") || CharSearch(argv[1], "EL") || CharSearch(argv[1], "e") ) channel = "el";
  //else if ( (argv[1]).find("all") != std::string::npos || (argv[1]).find("All") != std::string::npos || (argv[1]).find("ALL") != std::string::npos ) channel = "all";
  
  if ( argc > 1 )
  {
    for (int i = 1; i < argc; i++)
    {
      if ( strcmp(argv[i], "--runWidth") == 0 )
      {
        if ( i+1 >= argc )
        {
          cerr << "Cannot find argument. Exiting..." << endl;
          exit(1);
        }
        double argi = strtod(argv[i+1], NULL);
        applyWidthSF = true;
        scaleWidth = argi;
        
        calculateFractions = false;
        makeTGraphs = false;
        calculateLikelihood = true;
        makePlots = false;
      }
      else if ( strcmp(argv[i], "--runSyst") == 0 || strcmp(argv[i], "--runSystematics") == 0 )
      {
        if ( i+2 >= argc )
        {
          cerr << "Cannot find argument. Exiting..." << endl;
          exit(1);
        }
        int argi = strtol(argv[i+2], NULL, 10);
        
        if ( strcmp(argv[i+1], "rate") == 0 || strcmp(argv[i+1], "Rate") == 0 )
        {
          if ( argi > nRateSystematics )
          {
            cerr << "Argument " << argi << " is not a valid rate systematic. Exiting" << endl;
            exit(1);
          }
          runSystematics = true;
          nSystematics = 1;
          thisSystematic = listRateSyst[argi];
        }
        else if ( strcmp(argv[i+1], "sample") == 0 || strcmp(argv[i+1], "Sample") == 0 || strcmp(argv[i+1], "shape") == 0 || strcmp(argv[i+1], "Shape") == 0 )
        {
          if ( argi > nSampleSystematics )
          {
            cerr << "Argument " << argi << " is not a valid sample systematic. Exiting" << endl;
            exit(1);
          }
          runSystematics = true;
          useOtherXml = true;
          nSystematics = 1;
          thisSystematic = listSampleSyst[argi];
          thisDataSetName = dataSetNameSyst[argi];
          if ( thisDataSetName.find("fsr") != std::string::npos || thisDataSetName.find("isr") != std::string::npos || thisDataSetName.find("herwig") != std::string::npos ) ntupleSystDate = ntupleSystDateExtraJES;
        }
        else
        {
          cerr << "Did not find systematic " << argv[i] << " " << argv[i+1] << " " << argv[i+2] << ". Exiting..." << endl;
          exit(1);
        }
        cout << "Running systematic " << thisSystematic << endl;
        
//         calculateFractions = false;
//         makeTGraphs = false;
//         calculateLikelihood = true;
//         makePlots = false;
      }
      
      if ( strcmp(argv[i], "--makePlots") == 0 || strcmp(argv[i], "--plots") == 0 )
      {
        makePlots = true;
      }
    }
  }
  
  
  
  if (calculateResolutionFunctions)
  {
    testTTbarOnly = true;
    doGenOnly = true;
    makePlots = false;
    calculateAverageMass = false;
    calculateAverageMassAllMC = false;
    calculateLikelihood = false;
    makeTGraphs = false;
    doPseudoExps = false;
    doKinFit = false;
    runListWidths = false;
    runSystematics = false;
  }
  if (calculateAverageMass)
  {
    calculateAverageMassAllMC = false;
    calculateFractions = false;
    makeTGraphs = false;
    calculateLikelihood = false;
    doPseudoExps = false;
    makePlots = false;
    doGenOnly = false;
    runListWidths = false;
    runSystematics = false;
  }
  if (test) makePlots = false;
  if (testHistos)
  {
    makePlots = true;
    makeControlPlots = true;
    doGenOnly = false;
    makeTGraphs = false;
    runListWidths = false;
    runSystematics = false;
  }
  if (doGenOnly)
  {
    doPseudoExps = false;
    runSystematics = false;
  }
  if (calculateFractions)
  {
    calculateLikelihood = false;
    doPseudoExps = false;
    runListWidths = false;
    runSystematics = false;
  }
  if (makeTGraphs)
  {
    calculateFractions = true;
    runListWidths = false;
    calculateLikelihood = false;
  }
  if (calculateLikelihood) makeTGraphs = false;
  else
  {
    makeLikelihoodPlots = false;
    doPseudoExps = false;
  }
  if (doPseudoExps)
  {
    makePlots = false;
    runListWidths = false;
  }
  if (runListWidths)
  {
    makePlots = false;
    unblind = false;
  }
  if (runSystematics)
  {
    runListWidths = false;
    unblind = false;
    if (! makeMassTemplates && ! makePlots)
    {
      makeTGraphs = false;
      calculateLikelihood = true;
    }
    doPseudoExps = false;
  }
  else
  {
    runRateSystematics = false;
    runSampleSystematics = false;
  }
  //if (! runRateSystematics && ! runSampleSystematics) runSystematics = false;
  
  if (! makePlots) makeControlPlots = false;
  
  pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  if (makePlots || makeLikelihoodPlots)
  {
    // Add channel to output path
    pathOutput += channel+"/";
    mkdir(pathOutput.c_str(),0777);
    if (doPseudoExps)
    {
      pathOutput+= "PseudoExp/";
      mkdir(pathOutput.c_str(),0777);
    }
    if (testHistos)
    {
      pathOutput += "test/";
      mkdir(pathOutput.c_str(),0777);
    }
    // Give timestamp to output path
    pathOutput += dateString+"/";
    mkdir(pathOutput.c_str(),0777);
  }
  
  if (makeTGraphs)
  {
    mkdir(outputDirLL.c_str(),0777);
    outputDirLL += dateString+"/";
    mkdir(outputDirLL.c_str(),0777);
  }
  else if (calculateLikelihood)
  {
    //if (useNewVar) inputDateLL = whichTemplates();
    inputDirLL = outputDirLL+inputDateLL;
  }
  
  if (useNewVar)
  {
    /// m_bjj/m_jj
    minCutRedTopMass = 1.5;
    maxCutRedTopMass = 2.6;
    /// red mlb
    minCutRedTopMass = 0.2;
    maxCutRedTopMass = 2.1;
  }
  
  pathNtuplesMC = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.first+"/";
  pathNtuplesData = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleDate.second+"/";
  cout << "Using Ntuples from " << ntupleDate.first << " for MC and " << ntupleDate.second << " for data. This corresponds to systematics: " << systStr << endl;
  pathNtuplesSyst = "NtupleOutput/MergedTuples/"+channel+"/"+ntupleSystDate+"/";
  if (calculateAverageMass) cout << "Calculating average mass values..." << endl;
  if (calculateLikelihood) cout << "Calculating -loglikelihood values using templates from " << inputDateLL << endl;
  if (doPseudoExps)       cout << "              for pseudo experiments" << endl;
  if (testHistos) cout << "Testing histogram consistency..." << endl;
  if (doGenOnly) cout << "Running only matching..." << endl;
  if (applyWidthSF) cout << "TTbar sample width will be scaled by a factor " << scaleWidth << endl;
  else scaleWidth = 1.;
  
  bool saveSetUp[] = {makePlots, makeControlPlots, makeLikelihoodPlots, calculateFractions, makeTGraphs, useTTTemplates, calculateLikelihood, doPseudoExps};
  
  if (calculateAverageMassAllMC)
  {
    makeTGraphs = false;
    calculateLikelihood = false;
    doPseudoExps = false;
    makePlots = false;
    doGenOnly = false;
  }
  
  bool runAll = false;
  if (runTTbar && runSTtW && runSTt && runOther) runAll = true;
  
  if (testTTbarOnly || (runTTbar && ! runSTtW && ! runSTt && ! runOther) )
    cout << "Only running ttbar events..." << endl;
  else if (runAll) cout << "Running all datasets" << endl;
  else if (runTTbar) cout << "Running ttbar events..." << endl;
  else if (runSTtW || runSTt)
  {
    if (runSTtW && runSTt && ! runTTbar && ! runOther) cout << "Only running ST events..." << endl;
    else if (runSTtW && ! runSTt && ! runTTbar && ! runOther) cout << "Only running ST tW events..." << endl;
    else if (! runSTtW && runSTt && ! runTTbar && ! runOther) cout << "Only running ST t events..." << endl;
    else if (runSTtW && runSTt) cout << "Running ST events..." << endl;
    else if (runSTtW) cout << "Running ST tW events..." << endl;
    else if (runSTt) cout << "Running ST t events..." << endl;
  }
  else if (runOther)
  {
    if (! runTTbar && ! runSTtW && ! runSTt) cout << "Only running DY+jets and W+jets events... " << endl;
    else cout << "Running DY+jets and W+jets events... " << endl;
  }
  
  
  /// xml file
  string xmlFileName ="config/topWidth.xml";
  
  const char *xmlfile = xmlFileName.c_str();
  
  cout << " - Using config file " << xmlfile << endl;
  
  
  
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  datasets.clear();
  TTreeLoader treeLoader; 
  treeLoader.LoadDatasets(datasets, xmlfile);
  treeLoader.LoadDatasets(datasetsMSP, xmlfile);
  
  int dTT = -1;
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    string dataSetName = datasets[d]->Name();
    
    if ( dataSetName.find("TT") != std::string::npos )
      dataSetNames.push_back("TT_nominal");   // Temporarily ! Change to "TT" when making new (final?) templates
    else
      dataSetNames.push_back(dataSetName);
    
    
    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA") != std::string::npos )
      includeDataSets.push_back(0);
//      Luminosity = datasets[d]->EquivalentLumi();
    else if ( dataSetName.find("TT") != std::string::npos )
      includeDataSets.push_back(1);
    else if ( dataSetName.find("ST") != std::string::npos )
      includeDataSets.push_back(1);
    else
      includeDataSets.push_back(1);
    
    if ( dataSetName.find("QCD") != std::string::npos )
    {
      datasets[d]->SetColor(kYellow);
      datasetsMSP[d]->SetColor(kYellow);
    }
    if ( dataSetName.find("TT") != std::string::npos )
    {
      datasets[d]->SetTitle("t#bar{t}");
      datasets[d]->SetColor(kRed+1);
      datasetsMSP[dTT+2]->SetName("TT_other");
      datasetsMSP[dTT+2]->SetTitle("t#bar{t} unm. other");
      datasetsMSP[dTT+2]->SetColor(kRed-10);
      
      dTT = d;
    }
    //if ( dataSetName.find("TTbarJets_Other") != std::string::npos ) datasets[d]->SetColor(kRed-7);
    if ( dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      //datasets[d]->SetColor(kGreen-3);
      datasets[d]->SetColor(kBlue-2);
      datasetsMSP[d]->SetTitle("W#rightarrowl#nu");
      //datasetsMSP[d]->SetColor(kGreen-3);
      datasetsMSP[d]->SetColor(kBlue-2);
    }
    if ( dataSetName.find("ZJets") != std::string::npos || dataSetName.find("DY") != std::string::npos )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{#font[122]{\55}}");
      //datasets[d]->SetColor(kAzure-2);
      //datasets[d]->SetColor(kMagenta);
      datasets[d]->SetColor(kAzure+6);
      datasetsMSP[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{#font[122]{\55}}");
      //datasetsMSP[d]->SetColor(kAzure-2);
      //datasetsMSP[d]->SetColor(kMagenta);
      datasetsMSP[d]->SetColor(kAzure+6);
    }
    if ( dataSetName.find("ST") != std::string::npos || dataSetName.find("SingleTop") != std::string::npos )
    {
      datasets[d]->SetTitle("ST");
      //datasets[d]->SetColor(kBlue-2);
      datasets[d]->SetColor(kOrange-4);  // 595, 615, 800
      datasetsMSP[d]->SetTitle("ST");
      //datasetsMSP[d]->SetColor(kBlue-2);
      datasetsMSP[d]->SetColor(kOrange-4);
      //if ( dataSetName.find("tW") != std::string::npos )
      //{
      //  datasets[d]->SetTitle("ST tW");
      //  datasets[d]->SetColor(kBlue-4);
      //}
      //else
      //  datasets[d]->SetTitle("ST t");
    }
  }
  
  treeLoader.LoadDatasets(datasetsTemp, xmlfile);
  datasetsTemp[dTT]->SetName("TT_UM");
  datasetsTemp[dTT]->SetTitle("t#bar{t} unm. signal");
  datasetsTemp[dTT]->SetColor(kRed-9);
  datasetsMSP.insert(datasetsMSP.begin()+dTT, 1, datasetsTemp[dTT]);
  datasetsTemp.clear();
  treeLoader.LoadDatasets(datasetsTemp, xmlfile);
  datasetsTemp[dTT]->SetName("TT_WM");
  datasetsTemp[dTT]->SetTitle("t#bar{t} wrong");
  datasetsTemp[dTT]->SetColor(kRed-7);
  datasetsMSP.insert(datasetsMSP.begin()+dTT, 1, datasetsTemp[dTT]);
  datasetsTemp.clear();
  treeLoader.LoadDatasets(datasetsTemp, xmlfile);
  datasetsTemp[dTT]->SetName("TT_CM");
  datasetsTemp[dTT]->SetTitle("t#bar{t} correct");
  datasetsTemp[dTT]->SetColor(kRed-3);
  datasetsMSP.insert(datasetsMSP.begin()+dTT, 1, datasetsTemp[dTT]);
  
  dMSPmax = datasetsMSP.size();
  
//   cout << "Content datasets  " << endl;
//   for (int d = 0; d < datasets.size(); d++)
//     cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
//   cout << "Content datasetsMSP  " << endl;
//   for (int d = 0; d < datasetsMSP.size(); d++)
//     cout << "   Dataset " << d << ": " << datasetsMSP[d]->Name() << " / title : " << datasetsMSP[d]->Title() << endl;
  
  
  /// Load systematic samples
  if (runSampleSystematics)
  {
    datasetsSyst.clear();
    treeLoader.LoadDatasets(datasetsSyst, xmlSyst);
    cout << "Found " << datasetsSyst.size() << " extra systematics samples..." << endl;
  }
  
  /// Selection table
  if (! runListWidths && ! runSystematics && systStr.find("nominal") != std::string::npos )
  {
    //SelectionTables *selTab = new SelectionTables(datasets);
    SetUpSelectionTable();
  }
  
  
  ////////////////////////
  ///  Initialise ...  ///
  ////////////////////////
  
  EventReweighting *rew = new EventReweighting(false);  // no correction for number of events
  ResolutionFunctions* rf = new ResolutionFunctions(calculateResolutionFunctions, true);
//  KinFitter *kf;
//  Likelihood *like;
  
  if (! calculateResolutionFunctions)
  {
    kf = new KinFitter("input/PlotsForResolutionFunctions_testFit_171025.root", addWMassKF, addEqMassKF);
  }
  
  
  InitSetUp();
  
  vJER.clear(); vJES.clear(); vPU.clear();
  
  
  if (calculateAverageMass)
  {
    mkdir("averageMass/",0777);
    txtMassGenPMatched.open(("averageMass/mass_genp_matched_TT_"+dateString+".txt").c_str());
    txtMassGenJMatched.open(("averageMass/mass_genj_matched_TT_"+dateString+".txt").c_str());
    txtMassRecoCM.open(("averageMass/mass_reco_matched_TT_"+dateString+".txt").c_str());
    txtMassRecoWMUM.open(("averageMass/mass_reco_notCorrectMatch_TT_"+dateString+".txt").c_str());
    txtMassRecoUM.open(("averageMass/mass_reco_notMatched_TT_"+dateString+".txt").c_str());
    txtMassRecoWM.open(("averageMass/mass_reco_wrongPerm_TT_"+dateString+".txt").c_str());
  }
  
  txtDebugTopMass.open("debug_missing_topQ.txt");
  txtDebugTopMass << "## Events where top quark(s) not found in genParticle collection" << endl;
  txtDebugTopMass << "#  If lepton charge > 0 : leptonically decaying top, hadronically decaying antitop" << endl;
  txtDebugTopMass << "#  If lepton charge < 0 : hadronically decaying top, leptonically decaying antitop" << endl;
  
  
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName;
  double timePerDataSet[datasets.size()];
  
  int nEntries;
  double fracDataEras[2] = {-1.};
  //GetEraFraction(fracDataEras);
  fracDataEras[0] = lumi_runBCDEF/(lumi_runBCDEF + lumi_runGH);
  fracDataEras[1] = lumi_runGH/(lumi_runBCDEF + lumi_runGH);
  if ( fracDataEras[0] == -1. || fracDataEras[1] == -1. )
  {
    cerr << "Something went wrong with the fraction calculation for muon SFs!" << endl;
    exit(1);
  }
  cout << "The muon scale factors will be scaled by " << fracDataEras[0] << " for eras B-F and " << fracDataEras[1] << " for eras G-H." << endl;
  
  int dMSP;
  bool hasFoundTTbar = false;
  bool doReweighting = false;
  
  /// Loop over systematics or widths
  if (runRateSystematics) nSystematics = nRateSystematics;
  else if (runSampleSystematics) nSystematics = nSampleSystematics;
  int endSys = nSystematics;
  if (runListWidths) endSys = nWidths;
  if (! runSystematics && ! runListWidths) endSys = 1;
  for (int iSys = 0; iSys < endSys; iSys++)
  {
    if (runListWidths)
    {
      thisWidth = listWidths[iSys];
      cout << endl << "Running over widths... Now at width " << thisWidth << " x SM width... " << endl;
    }
    else
    {
      thisWidth = scaleWidth;
      cout << endl << "Width for ttbar sample is " << thisWidth << " x SM width..." << endl;
    }
    
    thisMass = 172.5;
    
    if (runRateSystematics)
    {
      thisSystematic = listRateSyst[iSys];
      cout << endl << "Running over systematics...  " << thisSystematic << endl;
      if (! applyLeptonSF && thisSystematic.find("lepton") != std::string::npos ) continue;
      else if (! applyBTagSF && thisSystematic.find("btag") != std::string::npos ) continue;
      else if (! applyPU && thisSystematic.find("pu") != std::string::npos ) continue;
      else if ( thisSystematic.find("lumiup") != std::string::npos ) Luminosity *= 1.025;
      else if ( thisSystematic.find("lumidown") != std::string::npos ) Luminosity *= 0.975;
    }
    else if (runSampleSystematics)
    {
      thisSystematic = listSampleSyst[iSys];
      cout << endl << "Running over systematics...  " << thisSystematic << endl;
      
      useOtherXml = true;
      thisDataSetName = dataSetNameSyst[iSys];
      
      if ( thisSystematic.find("mass") != std::string::npos )
      {
        if ( thisSystematic.find("169") != std::string::npos )
        {
          origMass = 169.5;
          thisMass = 169.5;
          origWidth = 0.9424;  // width +-2%/GeV --> in sample = 1.23
          thisWidth = 1.;
        }
        else if ( thisSystematic.find("171") != std::string::npos )
        {
          origMass = 171.5;
          thisMass = 171.5;
          origWidth = 0.9776;  // width +-2%/GeV --> in sample = 1.28
          thisWidth = 1.;
        }
        else if ( thisSystematic.find("173") != std::string::npos )
        {
          origMass = 173.5;
          thisMass = 173.5;
          origWidth = 1.0234;  // width +-2%/GeV --> in sample = 1.34
          thisWidth = 1.;
        }
        else if ( thisSystematic.find("175") != std::string::npos )
        {
          origMass = 175.5;
          thisMass = 175.5;
          origWidth = 1.065;  // width +-2%/GeV --> in sample = 1.39
          thisWidth = 1.;
        }
        cout << "Rescaling mass sample from (" << origWidth << ", " << origMass << ") to (" << thisWidth << thisMass << ") ..." << endl;
      }
    }
    
    
    /// Clear counters and likelihood
    nofCM = 0; nofWM = 0; nofUM = 0;
    nofCM_TT = 0; nofWM_TT = 0; nofUM_TT = 0;
    nofCMl = 0; nofWMl = 0; nofUMl = 0;
    nofCM_weighted = 0.; nofWM_weighted = 0.; nofUM_weighted = 0.; nofUM_TTsemilep_weighted = 0.; nofUM_TTother_weighted = 0.; nofUM_other_weighted = 0.;
    nofCMout_weighted = 0.; nofWMout_weighted = 0.; nofUMout_weighted = 0.; nofUMout_TTsemilep_weighted = 0.; nofUMout_TTother_weighted = 0.; nofUMout_other_weighted = 0.;
    nofCM_gate1 = 0.; nofWM_gate1 = 0.; nofUM_gate1 = 0.; nofCM_gate2 = 0.; nofWM_gate2 = 0.; nofUM_gate2 = 0.;
    for (unsigned int i = 0; i < dMSPmax; i++)
    {
      nofBKF_weighted[i] = 0.;
      nofAKF_weighted[i] = 0.;
    }
    sumEvents = 0.; sumCMEvents = 0.;
    sumTopMass = 0.; sumCMTopMass = 0.;
    sumMlb = 0.; sumCMMlb = 0.;
    
    
    if (calculateLikelihood)
    {
      if (doLikeW)       likeW->ClearLikelihoods();
      else if (doLikeM)  likeM->ClearLikelihoods();
      else if (doLike2D) like2D->ClearLikelihoods();
    }
    
    hasFoundTTbar = false;
    doReweighting = false;
    bool skipEvent = false;
    
    eqLumi_TT = -1.;
    
    
    
    ////////////////////////////////////
    /// Loop over datasets
    ////////////////////////////////////
    
    //for (int d = 0; d < 1; d++)
    for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
    {
      clock_t startDataSet = clock();
      
      ClearMetaData();
      
      dataSetName = datasets[d]->Name();
      isData = false; isTTbar = false; isST = false; isSTtW = false; isSTt = false; isOther = false; isHerwig = false;
      if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data") != std::string::npos || dataSetName.find("DATA") != std::string::npos )
      {
        isData = true;
        pathNtuples = pathNtuplesData;
      }
      else if ( dataSetName.find("TT") != std::string::npos )
      {
        isTTbar = true;
        hasFoundTTbar = true;
        if ( dataSetName.find("width") != std::string::npos || dataSetName.find("Width") != std::string::npos )
        {
          if ( dataSetName.find("x0p2") != std::string::npos ) thisWidth = 0.2;
          else if ( dataSetName.find("x0p5") != std::string::npos ) thisWidth = 0.5;
          else if ( dataSetName.find("x4") != std::string::npos ) thisWidth = 4.;
          else if ( dataSetName.find("x8") != std::string::npos ) thisWidth = 8.;
          applyWidthSF = false;
          if (runListWidths)
          {
            runListWidths = false;
            endSys = 1;
          }
          else if (runSystematics)
          {
            cerr << "ERROR: Trying to run systematics on ttbar sample with non-nominal width... Exiting..." << endl;
            exit(1);
          }
        }
        if (useOtherXml) dataSetName = thisDataSetName;
        if ( dataSetName.find("herwig") != std::string::npos || dataSetName.find("HERWIG") != std::string::npos || dataSetName.find("Herwig") != std::string::npos ) isHerwig = true;
        if ( dataSetName.find("isr") != std::string::npos || dataSetName.find("fsr") != std::string::npos || dataSetName.find("hdamp") != std::string::npos || dataSetName.find("erdOn") != std::string::npos || dataSetName.find("ERD") != std::string::npos || dataSetName.find("gluon") != std::string::npos || dataSetName.find("tune") != std::string::npos || dataSetName.find("mass") != std::string::npos || dataSetName.find("mtop") != std::string::npos || isHerwig )
        {
          applyWidthSF = false;
          runListWidths = false;
        }
      }
      else if ( dataSetName.find("ST") != std::string::npos )
      {
        isST = true;
        if ( dataSetName.find("tW") != std::string::npos ) isSTtW = true;
        else isSTt = true;
      }
      else
      {
        isOther = true;
      }
      
      if (verbose > 1)
      {
        cout << "   Dataset " << d << ": " << dataSetName << " / title : " << datasets[d]->Title() << endl;
      }
      
      if (isData && skipData)
      {
        cout << "Skipping data..." << endl;
        continue;
      }
      
      if (isData && calculateAverageMassAllMC)
      {
        cout << "Skipping data when calculating average mass" << endl;
        continue;
      }
      
      if (isData && (runListWidths || runSystematics) && ! makePlots)
      {
        cout << "Skipping data";
        if (runSystematics) cout << " when running systematics";
        cout << "..." << endl;
        continue;
      }
      
      if (! isData && useTTTemplates && includeDataSets[d] == 0 )
      {
        cout << "Skipping dataset, because not included in likelihood calculation" << endl;
        continue;
      }
      
      doReweighting = false;
      if ( isTTbar && (applyWidthSF || runListWidths) ) doReweighting = true;
      
      if (! isData)
      {
        pathNtuples = pathNtuplesMC;
      }
      
      if (testTTbarOnly && ! isTTbar)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      
      if (! runTTbar && isTTbar)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      if (! runSTtW && isSTtW)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      if (! runSTt && isSTt)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      if (! runOther && isOther)
      {
        cout << "Skipping dataset..." << endl;
        continue;
      }
      
//       if (doPseudoExps && isOther)
//       {
//         cout << "Skipping DY+jets & W+jets datasets when doing pseudo experiments." << endl;
//         continue;
//       }
      
      
      ///////////////////////////////////////////
      ///  Initialise Jet Energy Corrections  ///
      ///////////////////////////////////////////

      vCorrParam.clear();
      InitJEC(isData, dataSetName);
      JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); //true means redo also L1
      
      
      
      
      string ntupleFileName = pathNtuples+"Ntuples_"+dataSetName+".root";
      if (useOtherXml && isTTbar) ntupleFileName = pathNtuplesSyst+"Ntuples_"+dataSetName+".root";
      
      /// Change name of ttbar dataset to TT (whether it is nominal or widthxX)
      if (isTTbar) dataSetName = "TT";
      
      tFileMap[dataSetName.c_str()] = new TFile(ntupleFileName.c_str(),"READ"); //create TFile for each dataset
      
      string tTreeName = "tree";
      string tStatsTreeName = "stats";
      
      /// Get meta data
      tStatsTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tStatsTreeName.c_str());
      GetMetaData(tStatsTree[dataSetName.c_str()], isData);
      
      tStatsTree[(dataSetName).c_str()]->GetEntry(0);
      vJER.push_back(appliedJER);
      vJES.push_back(appliedJES);
      vPU.push_back(appliedPU);
      
      /// eqLumi calculation
      if (isData)
      {
        lumiWeight = 1.;
        if (! runListWidths && ! runSystematics && systStr.find("nominal") != std::string::npos )
          selTab->SetEqLumi(d, Luminosity);
      }
      else
      {
        nEventsDataSet = GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData);
        if (isTTbar && ! isHerwig)
        {
          nEventsDataSet -= GetNEvents(tStatsTree[(dataSetName).c_str()], "nofTTEventsWithoutAGenTop", isData);
        }
        xSection = datasets[d]->Xsection();  // pb
        eqLumi = (double)nEventsDataSet/xSection;  // 1/pb
        if (! runListWidths && ! runSystematics && systStr.find("nominal") != std::string::npos )
          selTab->SetEqLumi(d, eqLumi);
        if (isTTbar) eqLumi_TT = eqLumi;
        
        lumiWeight = Luminosity/eqLumi;
        if ( eqLumi_TT != -1. )
          relativeSF = eqLumi_TT/eqLumi;
        else if ( makeTGraphs || calculateLikelihood )
        {
          cerr << "ERROR: Cannot find eqLumi of ttbar sample for relative scale factor... " << endl;
          exit(1);
        }
      }
      
      if (doPseudoExps)
      {
        lumiWeight = 1.;  // select same amount as in data, so no reweighting necessary
        
        if (! isData) toyMax = Luminosity/eqLumi;  // CHECK: better nDataEvts/nEvtsPassKinFit ??
                                                   //else nDataEvts = GetNEvents(tStatsTree[(dataSetName).c_str()], "nEventsSel", 1);
        else nDataEvts = nEventsAKF[d];
        cout << "PseudoExperiments::Number of selected data events: " << nDataEvts << endl;
        
        if (test)
          cout << "      Lumi : " << Luminosity << "/pb; eqLumi: " << eqLumi << "/pb." << endl;
        cout << "PseudoExperiments::Lumi/eqLumi = " << toyMax;
//         if (! isData)
//         {
//           toyMax *= 0.908;  // small overshoot in MC --> scale down --> fixed now
//           cout << " x 0.908 = " << toyMax;
//         }
        cout << endl;
      }
      
      if (calculateAverageMass)
      {
        txtMassReco.open(("averageMass/mass_reco_"+dataSetName+"_"+dateString+".txt").c_str());
      }
      
      
      /// Get data
      tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
      nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
      cout << "                nEntries  : " << nEntries << endl;
      if (isData) cout << "                Lumi    : " << Luminosity << "/pb" << endl;
      else
      {
        cout << "                eqLumi    : " << eqLumi << "/pb = " << nEventsDataSet << " / " << xSection << " pb" << endl;
        if (isTTbar && ! isHerwig)
          cout << "                              ( " << nEventsDataSet << " = " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData) << " - " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofTTEventsWithoutAGenTop", isData) << " )" << endl;
        cout << "                lumiWeight: " << lumiWeight << endl;
      }
      
      if (test && isData)
      {
        cout << "NEvents data:   " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nEvents", isData) << endl;
        cout << "NEvents Run B:  " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunB", isData) << endl;
        cout << "NEvents Run CD: " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunCD", isData) << endl;
        cout << "NEvents Run EF: " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunEF", isData) << endl;
        cout << "NEvents Run G:  " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunG", isData) << endl;
        cout << "NEvents Run H:  " << GetNEvents(tStatsTree[(dataSetName).c_str()], "nofEventsRunH", isData) << endl;
      }
      
      if (runSystematics && isTTbar )
      {
        if ( thisSystematic.find("renFac") != std::string::npos )
        {
          renFacSumNom = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1001", isData);
          renFacSum1002 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1002", isData);
          renFacSum1003 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1003", isData);
          renFacSum1004 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1004", isData);
          renFacSum1005 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1005", isData);
          renFacSum1007 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1007", isData);
          renFacSum1009 = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1009", isData);
          cout << endl << "            Sum weight1001: " << (long)renFacSumNom << endl;
        }
        else if ( thisSystematic.find("frag") != std::string::npos )
        {
          renFacSumNom = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumWeight1001", isData);
          fragUpSum       = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumUpFragWeight", isData);
          fragCentralSum  = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumCentralFragWeight", isData);
          fragDownSum     = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumDownFragWeight", isData);
          fragPetersonSum = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumPetersonFragWeight", isData);
          semilepbrUpSum   = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumSemilepbrUp", isData);
          semilepbrDownSum = GetNEvents(tStatsTree[(dataSetName).c_str()], "sumSemilepbrDown", isData);
        }
      }
      
      // Set branch addresses and branch pointers
      InitTree(tTree[dataSetName.c_str()], isData);
      
      
      
      ////////////////////////////////////
      ///  Loop on events
      ////////////////////////////////////
      
      int endEvent = nEntries;
      if (test || testHistos) endEvent = 1000;
      for (int ievt = 0; ievt < endEvent; ievt++)
      {
        ClearObjects();
        
        if (ievt%10000 == 0)
          std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
        
        //if (! isTTbar) continue;
        //if (ievt > 588000) break;
        //if (ievt != 2062 && ievt != 75831 && ievt != 113603 && ievt != 115687 && ievt != 155732 && ievt != 163161 && ievt != 186900 && ievt != 215759 && ievt != 233634 && ievt != 238021 && ievt != 243052 && ievt != 243674 && ievt != 266399 && ievt != 317190 && ievt != 317752 && ievt != 325854 && ievt != 330813 && ievt != 333620 && ievt != 347247 && ievt != 439571 && ievt != 450329 && ievt != 491328 && ievt != 510024 && ievt != 514196 && ievt != 538345 && ievt != 570225 && ievt != 576194 && ievt != 577278 && ievt != 587570) continue;
        
        
        /// Load event
        tTree[(dataSetName).c_str()]->GetEntry(ievt);
        
        
        /// Scale factors
        if (! isData)
        {
          thisLeptonIdSF = fracDataEras[0]*muonIdSF_BCDEF[0] + fracDataEras[1]*muonIdSF_GH[0];
          thisLeptonIsoSF = fracDataEras[0]*muonIsoSF_BCDEF[0] + fracDataEras[1]*muonIsoSF_GH[0];
          thisLeptonTrigSF = fracDataEras[0]*muonTrigSF_BCDEF[0] + fracDataEras[1]*muonTrigSF_GH[0];
          thisLeptonSF = muonTrackSF_eta[0] * thisLeptonIdSF * thisLeptonIsoSF * thisLeptonTrigSF;
          if (runSystematics)
          {
            if (runSampleSystematics)
            {
              if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
              if (applyBTagSF) { scaleFactor *= btagSF;}
              if (applyPU) { scaleFactor *= puSF;}
            }
            //else if (runRateSystematics)
            else
            {
              if ( thisSystematic.find("lepton") != std::string::npos )  // apply some lepton SF systematic
              {
                if (applyBTagSF) { scaleFactor *= btagSF;}
                if (applyPU) { scaleFactor *= puSF;}
                if ( thisSystematic.find("IdSFup") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * (fracDataEras[0]*muonIdSF_up_BCDEF[0] + fracDataEras[1]*muonIdSF_up_GH[0]) * thisLeptonIsoSF * thisLeptonTrigSF;}
                else if ( thisSystematic.find("IdSFdown") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * (fracDataEras[0]*muonIdSF_down_BCDEF[0] + fracDataEras[1]*muonIdSF_down_GH[0]) * thisLeptonIsoSF * thisLeptonTrigSF;}
                else if ( thisSystematic.find("IsoSFup") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * thisLeptonIdSF * (fracDataEras[0]*muonIsoSF_up_BCDEF[0] + fracDataEras[1]*muonIsoSF_up_GH[0]) * thisLeptonTrigSF;}
                else if ( thisSystematic.find("IsoSFdown") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * thisLeptonIdSF * (fracDataEras[0]*muonIsoSF_down_BCDEF[0] + fracDataEras[1]*muonIsoSF_down_GH[0]) * thisLeptonTrigSF;}
                else if ( thisSystematic.find("TrigSFup") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * thisLeptonIdSF * thisLeptonIsoSF * (fracDataEras[0]*muonTrigSF_up_BCDEF[0] + fracDataEras[1]*muonTrigSF_up_GH[0]);}
                else if ( thisSystematic.find("TrigSFdown") != std::string::npos ) { scaleFactor *= muonTrackSF_eta[0] * thisLeptonIdSF * thisLeptonIsoSF * (fracDataEras[0]*muonTrigSF_down_BCDEF[0] + fracDataEras[1]*muonTrigSF_down_GH[0]);}
                else if ( thisSystematic.find("TrkSFup") != std::string::npos ) { scaleFactor *= 1.01*thisLeptonSF;}
                else if ( thisSystematic.find("TrkSFdown") != std::string::npos ) { scaleFactor *= 0.99*thisLeptonSF;}
              }
              else if ( thisSystematic.find("btag") != std::string::npos )
              {
                if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
                if (applyPU) { scaleFactor *= puSF;}
                if ( thisSystematic.find("up") != std::string::npos ) { scaleFactor *= btagSF_up;}
                else if ( thisSystematic.find("down") != std::string::npos ) { scaleFactor *= btagSF_down;}
              }
              else if ( thisSystematic.find("pu") != std::string::npos )
              {
                if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
                if (applyBTagSF) { scaleFactor *= btagSF;}
                if ( thisSystematic.find("up") != std::string::npos ) { scaleFactor *= puSF_up;}
                else if ( thisSystematic.find("down") != std::string::npos ) { scaleFactor *= puSF_down;}
              }
              else if ( thisSystematic.find("renFac") != std::string::npos && isTTbar )
              {
                if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
                if (applyBTagSF) { scaleFactor *= btagSF;}
                if (applyPU) { scaleFactor *= puSF;}

                if      ( thisSystematic.find("1002") != std::string::npos ) scaleFactor *= weight1002 * renFacSumNom/renFacSum1002;
                else if ( thisSystematic.find("1003") != std::string::npos ) scaleFactor *= weight1003 * renFacSumNom/renFacSum1003;
                else if ( thisSystematic.find("1004") != std::string::npos ) scaleFactor *= weight1004 * renFacSumNom/renFacSum1004;
                else if ( thisSystematic.find("1005") != std::string::npos ) scaleFactor *= weight1005 * renFacSumNom/renFacSum1005;
                else if ( thisSystematic.find("1007") != std::string::npos ) scaleFactor *= weight1007 * renFacSumNom/renFacSum1007;
                else if ( thisSystematic.find("1009") != std::string::npos ) scaleFactor *= weight1009 * renFacSumNom/renFacSum1009;
              }
              else if ( thisSystematic.find("frag") != std::string::npos && isTTbar )
              {
                if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
                if (applyBTagSF) { scaleFactor *= btagSF;}
                if (applyPU) { scaleFactor *= puSF;}
                
                if      ( thisSystematic.find("fragUp") != std::string::npos )       scaleFactor *= upFragWeight * renFacSumNom/fragUpSum;
                else if ( thisSystematic.find("fragCentral") != std::string::npos )  scaleFactor *= centralFragWeight * renFacSumNom/fragCentralSum;
                else if ( thisSystematic.find("fragDown") != std::string::npos )     scaleFactor *= downFragWeight * renFacSumNom/fragDownSum;
                else if ( thisSystematic.find("fragPeterson") != std::string::npos ) scaleFactor *= petersonFragWeight * renFacSumNom/fragPetersonSum;
                else if ( thisSystematic.find("fragSemiLepBrUp") != std::string::npos )   scaleFactor *= semilepbrUp * renFacSumNom/semilepbrUpSum;
                else if ( thisSystematic.find("fragSemiLepBrDown") != std::string::npos ) scaleFactor *= semilepbrDown * renFacSumNom/semilepbrDownSum;
              }
              else
              {
                if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
                if (applyBTagSF) { scaleFactor *= btagSF;}
                if (applyPU) { scaleFactor *= puSF;}
              }
            }
          }
          else
          {
            if (makePlots && passedMETFilter)
            {
              MSPlot["nPVs_beforePU_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF);
              MSPlot["nPVs_afterPU_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF*puSF);
              MSPlot["nPVs_afterPU_up_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF*puSF_up);
              MSPlot["nPVs_afterPU_down_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF*puSF_down);
              MSPlot["nPVs_afterPU_begin16_"]->Fill(nvtx, datasets[d], true, lumi_runBCDEF*1000/eqLumi*thisLeptonSF*btagSF*puSF);
              MSPlot["nPVs_afterPU_end16_"]->Fill(nvtx, datasets[d], true, lumi_runGH*1000/eqLumi*thisLeptonSF*btagSF*puSF);
            }
            
            if (applyLeptonSF) { scaleFactor *= thisLeptonSF;}
            if (applyBTagSF) { scaleFactor *= btagSF;}
            if (applyPU) { scaleFactor *= puSF;}
            
            if (makePlots)
            {
              MSPlot["btag_SF_"]->Fill(btagSF, datasets[d], true, lumiWeight*scaleFactor);
            }
          }
        }
        else if (makePlots && passedMETFilter)
        {
          MSPlot["nPVs_beforePU_"]->Fill(nvtx, datasets[d], false, lumiWeight);
          MSPlot["nPVs_afterPU_"]->Fill(nvtx, datasets[d], false, lumiWeight);
          MSPlot["nPVs_afterPU_up_"]->Fill(nvtx, datasets[d], false, lumiWeight);
          MSPlot["nPVs_afterPU_down_"]->Fill(nvtx, datasets[d], false, lumiWeight);
          if ( isDataRunB || isDataRunC || isDataRunD || isDataRunE || isDataRunF )
            MSPlot["nPVs_afterPU_begin16_"]->Fill(nvtx, datasets[d], false, lumiWeight);
          else if ( isDataRunG || isDataRunH )
            MSPlot["nPVs_afterPU_end16_"]->Fill(nvtx, datasets[d], false, lumiWeight);
        }
        
        
        
        //////////////////////
        ///  Fill objects  ///
        //////////////////////
        
        muon.SetPtEtaPhiE(muon_pt[0], muon_eta[0], muon_phi[0], muon_E[0]);
        selectedLepton.push_back(muon);
        
        for (int iJet = 0; iJet < nJets; iJet++)
        {
          jet.Clear();
          jet.SetPtEtaPhiE(jet_pt[iJet], jet_eta[iJet], jet_phi[iJet], jet_E[iJet]);
          if ( jet_pt[iJet] < 250. ) selectedJets.push_back(jet);
        }
        
        if ( selectedJets.size() < 4 ) continue;
        
        for (int iJet = 0; iJet < selectedJets.size(); iJet++)
        {
          if ( jet_bdiscr[iJet] > CSVv2Medium )
          {
            selectedBJets.push_back(selectedJets[iJet]);
            bJetId.push_back(iJet);  /// selectedBJets[j] = selectedJets[bJetId[j]]
          }
        }
        
        if ( selectedBJets.size() < 2 ) continue;
        
        if (makePlots && passedMETFilter)
        {
          MSPlot["nJets_"]->Fill(selectedJets.size(), datasets[d], true, lumiWeight*scaleFactor);
          MSPlot["rho_"]->Fill(rho, datasets[d], true, lumiWeight*scaleFactor);
        }
        
        if ( ! calculateResolutionFunctions && selectedJets.size() > 4 ) continue;
        nofHardSelected++;
        
        if (! passedMETFilter) continue;
        nofMETCleaned++;
        
        skipEvent = false;
        for (int iJet = 0; iJet < selectedJets.size(); iJet++)
        {
          for (int jJet = iJet+1; jJet < selectedJets.size(); jJet++)
          {
            tempDR = ROOT::Math::VectorUtil::DeltaR(selectedJets[iJet], selectedJets[jJet]);
            if ( tempDR < 0.6 ) skipEvent = true;
          }
        }
        if (skipEvent) continue;
        nofAfterDRmincut++;
        
//         for (int iJet = 0; iJet < selectedJets.size(); iJet++)
//         {
//           if ( jet_bdiscr[iJet] > CSVv2Medium )
//           {
//             selectedBJets.push_back(selectedJets[iJet]);
//             bJetId.push_back(iJet);  /// selectedBJets[j] = selectedJets[bJetId[j]]
//           }
//         }
        //std::sort(selectedBJets.begin(),selectedBJets.end(),HighestPt());  // already the case
        
        /// label jets with highest b discr
        for (int iJet = 0; iJet < selectedBJets.size(); iJet++)
        {
          tempbdiscr = jet_bdiscr[bJetId[iJet]];
          if ( tempbdiscr > bdiscrTop )
          {
            bdiscrTop2 = bdiscrTop;
            bdiscrTop = tempbdiscr;
            
            labelB2 = labelB1;
            labelB1 = iJet;
          }
          else if ( tempbdiscr > bdiscrTop2 )
          {
            bdiscrTop2 = tempbdiscr;
            labelB2 = iJet;
          }
        }
        
        
        if (! doReweighting ) widthSF = 1.;
        //else if ( applyWidthSF && ! isTTbar ) widthSF = 1.;  // also for data
        
        if (makePlots)
        {
          if (isData)
          {
            MSPlot["nPVs_beforePU_aSel_"]->Fill(nvtx, datasets[d], false, lumiWeight*thisLeptonSF*btagSF);
            MSPlot["nPVs_afterPU_aSel_"]->Fill(nvtx, datasets[d], false, lumiWeight*thisLeptonSF*btagSF*puSF);
          }
          else
          {
            MSPlot["nPVs_beforePU_aSel_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF);
            MSPlot["nPVs_afterPU_aSel_"]->Fill(nvtx, datasets[d], true, lumiWeight*thisLeptonSF*btagSF*puSF);
          }
          MSPlot["leadingJet_pT_bKF_noTopCut_"]->Fill(selectedJets[0].Pt(), datasets[d], true, lumiWeight*scaleFactor);
        }
        
        
        
        /////////////////////////////
        ///  JET PARTON MATCHING  ///
        /////////////////////////////
        
        //if ( isTTbar || dataSetName.find("ST") != std::string::npos )
        if (! isData && ! isHerwig)
        {
          for (int iMC = 0; iMC < nMCParticles; iMC++)
          {
            mcpart.Clear();
            mcpart.SetPtEtaPhiE(mc_pt[iMC], mc_eta[iMC], mc_phi[iMC], mc_E[iMC]);
            mcParticles.push_back(mcpart);
          }
          
          foundTop22 = false; foundAntiTop22 = false; foundTop62 = false; foundAntiTop62 = false;
          foundLastCopyTop = false; foundLastCopyAntitop = false;
          for (unsigned int i = 0; i < mcParticles.size(); i++)
          {
            if ( test && verbose > 4 )
              cout << setw(3) << right << i << "  Status: " << setw(2) << mc_status[i] << "  pdgId: " << setw(3) << mc_pdgId[i] << "  Mother: " << setw(4) << mc_mother[i] << "  Granny: " << setw(4) << mc_granny[i] << "  Pt: " << setw(7) << left << mc_pt[i] << "  Eta: " << mc_eta[i] << endl;
            
            /// Find tops
            if ( mc_pdgId[i] == pdgID_top )  // isLastCopy() == status 62
            {
              if ( mc_status[i] == 22 )
              {
                topQuark = i;
                foundTop22 = true;
              }
              if ( topQuark == -9999 && mc_status[i] == 62 ) topQuark = i;
              
              if (makePlots && isTTbar)
              {
                histo1D["genTop_status"]->Fill(mc_status[i]);
                if ( mc_status[i] == 22 )
                  histo1D["genTop_status22_pT"]->Fill(mc_pt[i]);
                else if ( mc_status[i] == 62 )
                  histo1D["genTop_status62_pT"]->Fill(mc_pt[i]);
              }
              if (isTTbar && mc_isLastCopy[i])
              {
                foundLastCopyTop = true;
                if (makePlots) histo1D["genTop_isLastCopy_status"]->Fill(mc_status[i]);
              }
            }
            else if ( mc_pdgId[i] == -pdgID_top )
            {
              if ( mc_status[i] == 22 )
              {
                antiTopQuark = i;
                foundAntiTop22 = true;
              }
              if ( antiTopQuark == -9999 && mc_status[i] == 62 ) antiTopQuark = i;
              
              if (makePlots && isTTbar)
              {
                histo1D["genAntitop_status"]->Fill(mc_status[i]);
                if ( mc_status[i] == 22 )
                  histo1D["genAntitop_status22_pT"]->Fill(mc_pt[i]);
                else if ( mc_status[i] == 62 )
                  histo1D["genAntitop_status62_pT"]->Fill(mc_pt[i]);
              }
              if (isTTbar && mc_isLastCopy[i])
              {
                foundLastCopyAntitop = true;
                if (makePlots) histo1D["genAntitop_isLastCopy_status"]->Fill(mc_status[i]);
              }
            }
            
            if ( mc_pdgId[i] == pdgID_top && mc_status[i] == 62 )  // top
            {
              foundTop62 = true;
              if ( mcParticles[i].Pt() < 800. ) topPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
              else topPtSF = TMath::Exp(0.0615-0.0005*800.);
            }
            else if ( mc_pdgId[i] == -pdgID_top && mc_status[i] == 62 )  // antitop
            {
              foundAntiTop62 = true;
              if ( mcParticles[i].Pt() < 800. ) antiTopPtSF = TMath::Exp(0.0615-0.0005*mcParticles[i].Pt());
              else antiTopPtSF = TMath::Exp(0.0615-0.0005*800.);
            }
            
            
            /// Status restriction: Final state particle or particle from hardest process
            if ( (mc_status[i] > 1 && mc_status[i] <= 20) || mc_status[i] >= 30 ) continue;
            
            /// Muons   // from top or from W (in ST tW channel)
            if ( mc_pdgId[i] == 13 && mc_mother[i] == -24 )		// mu-, W-
            {
              foundMuMinus = true;
              if ( mc_granny[i] == -pdgID_top ) muMinusFromTop = true;  // t~
              
              if ( mc_status[i] == 23 ) genmuon = i;
              else if ( genmuon == -9999 ) genmuon = i;
            }
            if ( mc_pdgId[i] == -13 && mc_mother[i] == 24 )		// mu+, W+
            {
              foundMuPlus = true;
              if ( mc_granny[i] == pdgID_top ) muPlusFromTop = true;  // t
              
              if ( mc_status[i] == 23 ) genmuon = i;
              else if ( genmuon == -9999 ) genmuon = i;
            }
            
            /// Partons/gluons
            if ( abs(mc_pdgId[i]) < 6 /*|| abs(mc_pdgId[i]) == 21 */)  //light/b quarks, 6 should stay hardcoded, OR gluon
            {
              partons.push_back(mcParticles[i]);
              partonId.push_back(i);  /// partons[j] = mcParticles[partonId[j]]
              
              if ( abs(mc_pdgId[i]) == 5 && abs(mc_mother[i]) == pdgID_top )
                foundB.push_back(i);
              else if ( abs(mc_pdgId[i]) < 5 && abs(mc_mother[i]) == 24 && abs(mc_granny[i]) == pdgID_top )
                foundLight.push_back(i);
            }
            
          }  // end loop mcParticles
          
          if (makePlots && isTTbar)
          {
            if (foundLastCopyTop) histo1D["genTop_hasLastCopy"]->Fill(1.);
            else histo1D["genTop_hasLastCopy"]->Fill(0.);
            if (foundLastCopyAntitop) histo1D["genAntitop_hasLastCopy"]->Fill(1.);
            else histo1D["genAntitop_hasLastCopy"]->Fill(0.);
            if (foundTop22) histo1D["genTop_hasStatus22"]->Fill(1.);
            else histo1D["genTop_hasStatus22"]->Fill(0.);
            if (foundAntiTop22) histo1D["genAntitop_hasStatus22"]->Fill(1.);
            else histo1D["genAntitop_hasStatus22"]->Fill(0.);
          }
          
          if (verbose > 3)
          {
            cout << "Size mcParticles:   " << mcParticles.size() << endl;
            cout << "Size partons:       " << partons.size() << endl;
            cout << "Size selectedJets:  " << selectedJets.size() << endl;
          }
          
          if ( muMinusFromTop && muPlusFromTop )
          {
            if (test) cout << "Both tops decay leptonically... Event " << ievt << " will not be matched." << endl;
            doMatching = false;
          }
          else if ( foundMuMinus && foundMuPlus )
          {
            if (test) cout << "Found fully leptonic decay of ST tW... Event " << ievt << " will not be matched." << endl;
            doMatching = false;
          }
          
          if ( isTTbar && (topQuark == -9999 || antiTopQuark == -9999) )
          {
            if (makePlots) histo1D["leadingJet_pT_reco_no22or62"]->Fill(selectedJets[0].Pt());
            if ( thisWidth == 1 )
            {
              txtDebugTopMass << "Event " << ievt << ";  Event nb. " << evt_num << "; Run nb. " << run_num << endl;
              txtDebugTopMass << "Top mass id: " << topQuark << "; antiTop mass id: " << antiTopQuark << "; Lepton charge: " << muon_charge[0] << endl;
            }
            continue;
          }
//           else if ( dataSetName.find("ST_tW") != std::string::npos && ( (foundMuMinus && topQuark == -9999) || (foundMuPlus && antiTopQuark == -9999) ) )
//           {
//             if ( thisWidth == 1 )
//             {
//               txtDebugTopMass << "Event " << ievt << ";  Event nb. " << evt_num << "; Run nb. " << run_num << endl;
//               txtDebugTopMass << "Top mass id: " << topQuark << "; antiTop mass id: " << antiTopQuark << "; Lepton charge: " << muon_charge[0] << endl;
//             }
//             continue;
//           }
          
          nPartons = foundB.size() + foundLight.size();
          if (isTTbar)
          {
            if ( muMinusFromTop && muPlusFromTop )
            {
              nofTTdilep++;
              isTTother = true;
              if (makePlots) histo1D["nPartons_dilep"]->Fill(nPartons);
            }
            else if ( muMinusFromTop || muPlusFromTop )
            {
              nofTTsemilep++;
              isTTsemilep = true;
              if (makePlots) histo1D["nPartons_semilep"]->Fill(nPartons);
            }
            else
            {
              nofTThadr++;
              isTTother = true;
              if (makePlots) histo1D["nPartons_allhad"]->Fill(nPartons);
              if ( foundMuMinus && foundMuPlus )
                cout << "Event: " << setw(6) << right << ievt << ": Found mu from W boson, but not from top" << endl;
            }
          }
          
          if (runSystematics)
          {
//            if (! makeMassTemplates) doMatching = false;
            if ( thisSystematic.find("topPtRew") != std::string::npos )
            {
              if (isTTbar) topPtRewSF = TMath::Sqrt(topPtSF*antiTopPtSF);
              else topPtRewSF = 1.;
              scaleFactor *= topPtRewSF;
            }
          }
          
          
          
          /////////////////////////////////////////
          ///  Scale factor ttbar sample width  ///
          /////////////////////////////////////////
          
          if ( muon_charge[0] > 0. )
          {
            massHadTopQ = (mcParticles[antiTopQuark]).M();
            massLepTopQ = (mcParticles[topQuark]).M();
          }
          else if ( muon_charge[0] < 0. )
          {
            massHadTopQ = (mcParticles[topQuark]).M();
            massLepTopQ =  (mcParticles[antiTopQuark]).M();
          }
          
          if ( doReweighting || (runSystematics && thisSystematic.find("mass") != std::string::npos) )
          {
            if (rewHadTopOnly) widthSF = rew->BEventWeightCalculatorNonRel(massHadTopQ, origMass, thisMass, origWidth, thisWidth);
            else widthSF = rew->BEventWeightCalculatorNonRel(massHadTopQ, origMass, thisMass, origWidth, thisWidth) * rew->BEventWeightCalculatorNonRel(massLepTopQ, origMass, thisMass, origWidth, thisWidth);
            
            if ( widthSF != widthSF )  // widthSF = NaN
            {
              continue;
            }
            
            if (makePlots) histo1D["width_SF"]->Fill(widthSF);
          }  // end applyWidthSF
          
          
          
          //////////////////
          ///  Matching  ///
          //////////////////
          
          if (doMatching)
          {
            TruthMatching(partons, selectedJets, MCPermutation);
            
            if (test && verbose > 3 && all4PartonsMatched)
            {
              for (unsigned int iMatch = 0; iMatch < 4; iMatch++)
              {
                cout << "Event  " << right << setw(4) << ievt << ";  Matched parton " << iMatch << "  Pt: " << setw(7) << left << mc_pt[partonId[MCPermutation[iMatch].second]] << "  Eta: " << mc_eta[partonId[MCPermutation[iMatch].second]] << "  Phi: " << mc_phi[partonId[MCPermutation[iMatch].second]] << endl;
                cout << "Event  " << right << setw(4) << ievt << ";  Matched jet    " << iMatch << "  Pt: " << setw(7) << left << jet_pt[MCPermutation[iMatch].first] << "  Eta: " << jet_eta[MCPermutation[iMatch].first] << "  Phi: " << jet_phi[MCPermutation[iMatch].first] << endl;
              }
            }
            
            if (hadronicTopJetsMatched)
            {
              
              for (unsigned int iMatch = 0; iMatch < 3; iMatch++)
              {
                /// MCPermutation[i].first  = jet number
                /// MCPermutation[i].second = parton number
                /// 0,1: light jets from W; 2: hadronic b jet; 3: leptonic b jet
                
                partonsMatched.push_back(partons[MCPermutation[iMatch].second]);
                jetsMatched.push_back(selectedJets[MCPermutation[iMatch].first]);
              }
              if (all4PartonsMatched)
              {
                partonsMatched.push_back(partons[MCPermutation[3].second]);
                jetsMatched.push_back(selectedJets[MCPermutation[3].first]);
              }
              
              matched_W_mass_j = (jetsMatched[0] + jetsMatched[1]).M();
              matched_W_mass_q = (partonsMatched[0] + partonsMatched[1]).M();
              matched_top_mass_j = (jetsMatched[0] + jetsMatched[1] + jetsMatched[2]).M();
              matched_top_mass_q = (partonsMatched[0] + partonsMatched[1] + partonsMatched[2]).M();
              
              if (calculateAverageMass)
              {
                txtMassGenPMatched << ievt << "  " << matched_top_mass_q << "  " << scaleFactor << "  " << lumiWeight << endl;
                txtMassGenJMatched << ievt << "  " << matched_top_mass_j << "  " << scaleFactor << "  " << lumiWeight << endl;
              }
              
              if (isTTbar && makePlots)
              {
                FillMatchingPlots();
              }
              
            }  // end hadronicTopJetsMatched
            
            
            
            ///////////////////
            ///  Resolution functions
            ///////////////////
            
            if (isTTbar && calculateResolutionFunctions)
            {
              if (hadronicTopJetsMatched) rf->fillJets(partonsMatched, jetsMatched);
              
              if (muonmatched) rf->fillMuon(mcParticles[genmuon], selectedLepton[0]);
              //if (electronmatched) rf->fillElectron(...)
              
            }  // end rf
            
          }  // end doMatching
          
          
        }  // end if not data
        
        
        
        if (doGenOnly) continue;
        
        
        
        ////////////////////////////////
        ///  Jet Energy Corrections  ///
        ////////////////////////////////
        
        //cout << "Original jet  p_x: " << selectedJets[0].Px() << "  p_y: " << selectedJets[0].Py() << "  p_z: " << selectedJets[0].Pz() << "  E: " << selectedJets[0].E() << "  Eta: " << selectedJets[0].Eta() << "  Phi: " << selectedJets[0].Phi() << endl;
//         selJets.clear();
//         for (int ijet = 0; ijet < selectedJets.size(); ijet++)
//         {
//           selJets.push_back(&selectedJets[ijet]);
//         }
//         if (! isData)
//         {
//           //jetTools->correctJetJESUnc(selJets, "minus", 1); if (nofAfterDRmincut == 1) cout << "Applying ad-hoc JESdown corrections" << endl;
//           jetTools->correctJetJESUnc(selJets, "plus", 1); if (nofAfterDRmincut == 1) cout << "Applying ad-hoc JESup corrections" << endl;
//         }
//         for (int ijet = 0; ijet < selectedJets.size(); ijet++)
//         {
//           selectedJets[ijet] = *selJets[ijet];
//         }
        //cout << "Corrected jet p_x: " << selectedJets[0].Px() << "  p_y: " << selectedJets[0].Py() << "  p_z: " << selectedJets[0].Pz() << "  E: " << selectedJets[0].E() << "  Eta: " << selectedJets[0].Eta() << "  Phi: " << selectedJets[0].Phi() << endl;
        
//         skipEvent = false;
//         for (int ijet = 0; ijet < selectedJets.size(); ijet++)
//         {
//           if ( selectedJets[ijet].Pt() < 30. || selectedJets[ijet].Pt() > 250. )
//           {
//             skipEvent = true;
//             break;
//           }
//         }
//         if (skipEvent) continue;
        
        
        
        /////////////////////////////////////
        ///  Reconstruction of Top quark  ///
        /////////////////////////////////////
        
        //int labelsReco[4] = {-9999, -9999, -9999, -9999};  // 0,1: light jets; 2: hadronic b; 3: leptonic b
        double deltaR;
        double minDeltaR = 9999.;
        
        for (int kjet = 0; kjet < selectedBJets.size(); kjet++)
        {
          if ( selectedBJets.size() > 2 && kjet != labelB1 && kjet != labelB2 ) continue;
          
          deltaR = ROOT::Math::VectorUtil::DeltaR(selectedBJets[kjet], selectedLepton[0]);
          
          if (deltaR < minDeltaR)
          {
            minDeltaR = deltaR;
            labelsReco[3] = bJetId[kjet];
          }
        }
        
        if ( labelsReco[3] == -9999 ) continue;
        
        if ( labelsReco[3] == bJetId[labelB1] ) labelsReco[2] = bJetId[labelB2];
        else if ( labelsReco[3] == bJetId[labelB2] ) labelsReco[2] = bJetId[labelB1];
        else cerr << endl << "Seems like something went wrong with the b jets..." << endl;
        
        for (int ijet = 0; ijet < selectedJets.size(); ijet++)
        {
          if ( ijet == labelsReco[2] || ijet == labelsReco[3] ) continue;
          
          if ( labelsReco[0] == -9999 ) labelsReco[0] = ijet;
          else if ( labelsReco[1] == -9999 ) labelsReco[1] = ijet;
          //else cerr << endl << "Seems like there are too many jets..." << endl;
        }
        
        if ( labelsReco[0] == -9999 || labelsReco[1] == -9999 || labelsReco[2] == -9999 ) continue;
        
        
        
        ///////////////////////////////////
        ///  CHECK MATCHED COMBINATION  ///
        ///////////////////////////////////
        
        ///
        // 3 possibilities:
        // - correct top match: 3 jets selected with reco method correspond to the 3 matched jets (n.b. this is also true when the jets originating from the W boson and the b jet do not exactly correspond to the matched jets, because we are only interested in the reconstructed top quark.)
        // - wrong permutation: the correct jet combination exists in the selected jets, but is not chosen by the reco method.
        // - wrong (no) match:  the correct jet combination does not exist in the selected jets (e.g. when one jet is not selected.)
        
        
        if (! isData)
        {
          if (hadronicTopJetsMatched)
          {
            /// Correct match
            if ( ( labelsReco[0] == MCPermutation[0].first || labelsReco[0] == MCPermutation[1].first || labelsReco[0] == MCPermutation[2].first ) && ( labelsReco[1] == MCPermutation[0].first || labelsReco[1] == MCPermutation[1].first || labelsReco[1] == MCPermutation[2].first ) && ( labelsReco[2] == MCPermutation[0].first || labelsReco[2] == MCPermutation[1].first || labelsReco[2] == MCPermutation[2].first ) )  // correct jets for top quark
            {
              isCM = true;
              nofCorrectlyMatched++;
            }
            else  // wrong permutation
            {
              isWM = true;
              nofNotCorrectlyMatched++;
            }
          }  // end hadrTopMatch
          else  // no match
          {
            isUM = true;
            nofUnmatched++;
            if (isTTsemilep)
            {
              isUM_TTsemilep = true;
              nofUnmatchedTTsemilep++;
            }
            else if (isTTother)
            {
              isUM_TTother = true;
              nofUnmatchedTTother++;
            }
          }
          
          
          if ( (! isCM && ! isWM && ! isUM) || (isCM && isWM) || (isCM && isUM) || (isWM && isUM) )
            cerr << "Something wrong with trigger logic CM/WM/UM !! " << endl;
          
        }  // not Data
        
        
        if (isCM) catSuffix = catSuffixList[0];
        else if (isWM) catSuffix = catSuffixList[1];
        else if (isUM) catSuffix = catSuffixList[2];
        
        dMSP = d;
        if (hasFoundTTbar && ! isTTbar) dMSP = d+3;
        else if (isTTbar && isWM) dMSP = d+1;
        else if (isTTbar && isUM_TTsemilep) dMSP = d+2;
        else if (isTTbar && isUM_TTother) dMSP = d+3;
        
        nofBKF_weighted[dMSP] += lumiWeight*scaleFactor*widthSF;
        
        /// Fill variables before performing kinFit
        WCandidate.SetPxPyPzE((selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).Px(), (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).Py(), (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).Pz(), (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]]).E());
        reco_W_mass_bKF = WCandidate.M();
        reco_W_pt_bKF = WCandidate.Pt();
        reco_mTW_bKF = sqrt( pow(WCandidate.M(),2) + pow(WCandidate.Px(),2) + pow(WCandidate.Py(),2) );
        reco_top_mass_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).M();
        reco_top_mass_alt_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[3]]).M();
        reco_top_mass_rescaled_bKF = 80.385/reco_W_mass_bKF * reco_top_mass_bKF;
        reco_top_pt_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[2]]).Pt();
        reco_top_pt_alt_bKF = (selectedJets[labelsReco[0]] + selectedJets[labelsReco[1]] + selectedJets[labelsReco[3]]).Pt();
        reco_mlb_bKF = (selectedLepton[0] + selectedJets[labelsReco[3]]).M();
        reco_dRLepB_lep_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[3]], selectedLepton[0] );  // deltaR between lepton and leptonic b jet
        reco_dRLepB_had_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
        reco_dRlight_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[0]], selectedJets[labelsReco[1]] );  // deltaR between light jets
        reco_dRblight_min_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedJets[labelsReco[0]] );
        reco_dRblight_max_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedJets[labelsReco[1]] );
        if ( reco_dRblight_min_bKF > reco_dRblight_max_bKF )
        {
          tempDR = reco_dRblight_min_bKF;
          reco_dRblight_min_bKF = reco_dRblight_max_bKF;
          reco_dRblight_max_bKF = tempDR;
        }
        reco_dRbW_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], (selectedJets[labelsReco[0]]+selectedJets[labelsReco[1]]) );
        reco_dRblight_qsum_bKF = sqrt( pow ( ROOT::Math::VectorUtil::DeltaR(selectedJets[labelsReco[2]], selectedJets[labelsReco[0]]), 2.) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJets[labelsReco[2]], selectedJets[labelsReco[1]]), 2.) ) ;  // quadratic sum of deltaR between b and light jets
        reco_dRbb_bKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedJets[labelsReco[3]]);
        reco_dPhi_bb_bKF = ROOT::Math::VectorUtil::DeltaPhi(selectedJets[labelsReco[2]], selectedJets[labelsReco[3]]);
        reco_dPhi_light_bKF = ROOT::Math::VectorUtil::DeltaPhi(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]]);
        reco_dPhi_bW_bKF = ROOT::Math::VectorUtil::DeltaPhi(selectedJets[labelsReco[2]], (selectedJets[labelsReco[0]]+selectedJets[labelsReco[1]]) );
        reco_ttbar_mass_bKF = reco_mlb_bKF + reco_top_mass_bKF;
        
        redTopMass_old_bKF = reco_top_mass_bKF/aveTopMassCM;
        redTopMass_new_bKF = reco_top_mass_bKF/aveTopMassLL;
        redMlbMass_bKF = reco_mlb_bKF/aveMlbMassCM;
        redTopMass_bKF = redTopMass_old_bKF;
        
        reco_mbjj_div_mjj_bKF = reco_top_mass_bKF/reco_W_mass_bKF;
        
        if (isCM)      nofCM_gate1 += lumiWeight*scaleFactor*widthSF;
        else if (isWM) nofWM_gate1 += lumiWeight*scaleFactor*widthSF;
        else if (isUM) nofUM_gate1 += lumiWeight*scaleFactor*widthSF;
        
        /// Max dR cut
//        if ( reco_dRlight_bKF > 4. ) continue;
//        if ( reco_dRblight_min_bKF > 3.5 ) continue;
        
        
//         if ( reco_dRlight_bKF > 2.4 ) continue;
//         if ( reco_dRblight_min_bKF > 2.2 ) continue;
//         if ( reco_dRblight_max_bKF > 3.2 ) continue;
//         if ( reco_dRbW_bKF > 2.5 || reco_dRbW_bKF < 0.4 ) continue;
//         if ( reco_dRbb_bKF < 1.6 || reco_dRbb_bKF > 4.) continue;
//         if ( reco_dRLepB_lep_bKF > 2.2 ) continue;
//         if ( reco_dRLepB_had_bKF > 3.6 || reco_dRLepB_had_bKF < 1.5 ) continue;
        nofAfterDRmaxcut++;
        
        
        if ( reco_top_mass_bKF < minTopMass || reco_top_mass_bKF > maxTopMass ) continue;
        //if ( fabs( reco_top_mass_bKF - reco_top_mass_alt_bKF ) < 75. ) continue;
        
        //if ( reco_mbjj_div_mjj_bKF < 1.65 || reco_mbjj_div_mjj_bKF > 2.45 ) continue;
        
//         if ( selectedJets[3].Pt() > 140. ) continue;
//        if ( reco_ttbar_mass_bKF > 500. ) continue;
//         if ( reco_mlb_bKF < 40. || reco_mlb_bKF > 200. ) continue;
//        if ( reco_top_mass_bKF - reco_top_mass_alt_bKF > 75. ) continue;
        //Ht = selectedJets[0].Pt() + selectedJets[1].Pt() + selectedJets[2].Pt() + selectedJets[3].Pt();
        //if ( Ht < 220. ) continue;
        
        if (isCM)      nofCM_gate2 += lumiWeight*scaleFactor*widthSF;
        else if (isWM) nofWM_gate2 += lumiWeight*scaleFactor*widthSF;
        else if (isUM) nofUM_gate2 += lumiWeight*scaleFactor*widthSF;
        
        
        
        if (makePlots)
        {
          if (isTTbar && doKinFit) FillKinFitPlots(doneKinFit);
          
          FillMSPlots(dMSP, doneKinFit);
          //FillControlPlots(datasetsMSP, dMSP);
          
          if (! isData) histo2D["qqb1_vs_qqb2"+catSuffix]->Fill(reco_top_mass_bKF, reco_top_mass_alt_bKF, widthSF);
          
          MSPlot["leadingJet_pT_"]->Fill(selectedJets[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
          for (int iJet = 0; iJet < selectedJets.size(); iJet++)
          {
            MSPlot["jet_pT_allJets_"]->Fill(selectedJets[iJet].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
          }
          MSPlot["W_mass_"]->Fill(reco_W_mass_bKF, datasets[d], true, lumiWeight*scaleFactor*widthSF);
          MSPlot["top_mass_"]->Fill(reco_top_mass_bKF, datasets[d], true, lumiWeight*scaleFactor*widthSF);
        }
        
        
        
        ////////////////////////////
        ///   Kinematic Fit      ///
        ////////////////////////////
        
        if (doKinFit)
        {
          /*if (addEqMassKF)
           kFitter = kf->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]], selectedJets[labelsReco[2]], selectedJets[labelsReco[3]], selectedLepton[0], TLV NEUTRINO, kFitVerbosity);
           else if (addWMassKF)*/
          kFitter = kf->doFit(selectedJets[labelsReco[0]], selectedJets[labelsReco[1]], kFitVerbosity);
          
          if ( kFitter->getStatus() != 0 )  // did not converge
          {
            if (test && verbose > 2) cout << "Event " << ievt << ": Fit did not converge..." << endl;
            continue;
          }
          
          kFitChi2 = kFitter->getS();
          if (test && verbose > 4) cout << "Fit converged: Chi2 = " << kFitChi2 << endl;
          
          doneKinFit = true;
          if (makePlots)
          {
            MSPlot["KF_Chi2_wide"]->Fill(kFitChi2, datasetsMSP[dMSP], true, lumiWeight*scaleFactor*widthSF);
            if (! isData)
            {
              selectedJetsKFcorrected.clear();
              selectedJetsKFcorrected = kf->getCorrectedJets();
              reco_top_mass_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJets[labelsReco[2]]).M();
              WCandidate.SetPxPyPzE((selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Px(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Py(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Pz(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).E());
              reco_mTW_aKF = sqrt( pow(WCandidate.M(),2) + pow(WCandidate.Px(),2) + pow(WCandidate.Py(),2) );
              histo1D["KF_Chi2"+catSuffix+"_wide"]->Fill(kFitChi2);
              histo2D["KF_top_mass_orig_vs_corr_noCut"+catSuffix]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
              histo2D["KF_top_mass_chi2_noCut_logY"+catSuffix]->Fill(reco_top_mass_aKF, TMath::Log10(kFitChi2), widthSF);
              histo2D["KF_W_mass_T_orig_vs_corr_noCut"+catSuffix]->Fill(reco_mTW_bKF, reco_mTW_aKF, widthSF);
              /// Combine DY & W+jets
              if ( dataSetName.find("DY") != std::string::npos )
              {
                histo2D["KF_top_mass_orig_vs_corr_noCut_DYJets"]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
              }
              else if ( dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos )
              {
                histo2D["KF_top_mass_orig_vs_corr_noCut_WJets"]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
              }
              else
              {
                histo2D["KF_top_mass_orig_vs_corr_noCut_"+dataSetName]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
              }
            }
            if (isTTbar) histo1D["KF_Chi2_TT"]->Fill(kFitChi2);
          }
          
          if (isCM) nofCorrectlyMatchedAKFNoCut++;
          else if (isWM) nofNotCorrectlyMatchedAKFNoCut++;
          else if (isUM)
          {
            nofUnmatchedAKFNoCut++;
            if (isTTsemilep) nofUnmatchedTTsemilepAKFNoCut++;
            else if (isTTother) nofUnmatchedTTotherAKFNoCut++;
          }
          
          
          if ( applyKinFitCut && kFitChi2 > kinFitCutValue ) continue;
          if ( applyKinFitCut && kFitChi2 < kinFitMinCutValue ) continue;
          nofAcceptedKFit++;
          nofAcceptedKFitWeighted += scaleFactor;
          if (hadronicTopJetsMatched) nofHadrMatchedEventsAKF++;
          if (isCM) nofCorrectlyMatchedAKF++;
          else if (isWM) nofNotCorrectlyMatchedAKF++;
          else if (isUM)
          {
            nofUnmatchedAKF++;
            if (isTTsemilep) nofUnmatchedTTsemilepAKF++;
            else if (isTTother) nofUnmatchedTTotherAKF++;
          }
          
          selectedJetsKFcorrected.clear();
          selectedJetsKFcorrected = kf->getCorrectedJets();
          
        }
        
        nofAKF_weighted[dMSP] += lumiWeight*scaleFactor*widthSF;
        
        
        /// Reconstruct event
        //  Hadronic variables  // OBS: only W mass constraint ! Jet3 = selectedJets[labelsReco[2]] !
        if (! doKinFit)
        {
          selectedJetsKFcorrected.clear();
          selectedJetsKFcorrected.push_back(selectedJets[labelsReco[0]]);
          selectedJetsKFcorrected.push_back(selectedJets[labelsReco[1]]);
          selectedJetsKFcorrected.push_back(selectedJets[labelsReco[2]]);
        }
        else if ( selectedJetsKFcorrected.size() == 2 ) selectedJetsKFcorrected.push_back(selectedJets[labelsReco[2]]);
        
        /// Make pT ordered jet collection after KF
        selectedJetsAKF = selectedJetsKFcorrected;
        selectedJetsAKF.push_back(selectedJets[labelsReco[3]]);
        std::sort(selectedJetsAKF.begin(),selectedJetsAKF.end(),HighestPt());
        
        
        ////////////////////////////////
        ///  Jet Energy Corrections  ///
        ////////////////////////////////
        
//         //cout << "Original jet  p_x: " << selectedJets[0].Px() << "  p_y: " << selectedJets[0].Py() << "  p_z: " << selectedJets[0].Pz() << "  E: " << selectedJets[0].E() << endl;
//         selJets.clear();
//         for (int ijet = 0; ijet < selectedJetsKFcorrected.size(); ijet++)
//         {
//           selJets.push_back(&selectedJetsKFcorrected[ijet]);
//         }
//         if (! isData)
//         {
//           //jetTools->correctJetJESUnc(selJets, "minus", 1); if (nofAfterDRmincut == 1) cout << "Applying ad-hoc JESdown corrections" << endl;
//           jetTools->correctJetJESUnc(selJets, "plus", 1); if (nofAfterDRmincut == 1) cout << "Applying ad-hoc JESup corrections" << endl;
//         }
//         for (int ijet = 0; ijet < selectedJetsKFcorrected.size(); ijet++)
//         {
//           selectedJetsKFcorrected[ijet] = *selJets[ijet];
//         }
        
        
        
        /// Define variables
        WCandidate.SetPxPyPzE((selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Px(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Py(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Pz(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).E());
        reco_W_mass_aKF = WCandidate.M();
        reco_W_pt_aKF = WCandidate.Pt();
        reco_mTW_aKF = sqrt( pow(WCandidate.M(),2) + pow(WCandidate.Px(),2) + pow(WCandidate.Py(),2) );
        
        reco_top_mass_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).M();
        reco_top_mass_alt_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJets[labelsReco[3]]).M();
        reco_top_pt_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).Pt();
        reco_top_pt_alt_aKF = (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1] + selectedJets[labelsReco[3]]).Pt();
        reco_mlb_aKF = (selectedLepton[0] + selectedJets[labelsReco[3]]).M();
        reco_dRLepB_lep_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[3]], selectedLepton[0] );  // deltaR between lepton and leptonic b jet
        reco_dRLepB_had_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJets[labelsReco[2]], selectedLepton[0] );  // deltaR between lepton and hadronic b jet
        reco_dRlight_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJetsKFcorrected[0], selectedJetsKFcorrected[1] );  // deltaR between light jets
        reco_dRblight_min_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJetsKFcorrected[2], selectedJetsKFcorrected[0] );
        reco_dRblight_max_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJetsKFcorrected[2], selectedJetsKFcorrected[1] );
        if ( reco_dRblight_min_aKF > reco_dRblight_max_aKF )
        {
          tempDR = reco_dRblight_min_aKF;
          reco_dRblight_min_aKF = reco_dRblight_max_aKF;
          reco_dRblight_max_aKF = tempDR;
        }
        reco_dRbW_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJetsKFcorrected[2], (selectedJetsKFcorrected[0]+selectedJetsKFcorrected[1]) );
        reco_dRblight_qsum_aKF = sqrt( pow ( ROOT::Math::VectorUtil::DeltaR(selectedJetsKFcorrected[2], selectedJetsKFcorrected[0]), 2.) + pow( ROOT::Math::VectorUtil::DeltaR(selectedJetsKFcorrected[2], selectedJetsKFcorrected[1]), 2.) );  // quadratic sum of deltaR between b and light jets
        reco_dRbb_aKF = ROOT::Math::VectorUtil::DeltaR( selectedJetsKFcorrected[2], selectedJets[labelsReco[3]]);
        reco_dPhi_bb_aKF = ROOT::Math::VectorUtil::DeltaPhi(selectedJetsKFcorrected[2], selectedJets[labelsReco[3]]);
        reco_dPhi_light_aKF = ROOT::Math::VectorUtil::DeltaPhi(selectedJetsKFcorrected[0], selectedJetsKFcorrected[1]);
        reco_dPhi_bW_aKF = ROOT::Math::VectorUtil::DeltaPhi(selectedJetsKFcorrected[2], (selectedJetsKFcorrected[0]+selectedJetsKFcorrected[1]) );
        reco_ttbar_mass_aKF = reco_mlb_aKF + reco_top_mass_aKF;
        
        reco_mbjj_div_mjj = reco_top_mass_bKF/reco_W_mass_bKF;
        
        /// Temporarily !
        //reco_top_mass_aKF += 1.5; // Manually shift top mass distribution by X GeV
        
        redTopMass_old = reco_top_mass_aKF/aveTopMassCM;
        redTopMass_new = reco_top_mass_aKF/aveTopMassLL;
        redMlbMass = reco_mlb_aKF/aveMlbMassCM;
        redTopMass = redTopMass_old;
        //redTopMass = reco_top_mass_aKF/aveTopMassLL;
        //redTopMass = reco_top_mass_bKF/aveTopMassLL;
        //redTopMass = reco_top_mass_aKF/aveTopMassCM + (aveTopMassCM - aveTopMassLL)/aveTopMassCM;
        //reco_new_var = reco_mbjj_div_mjj;
        reco_new_var = redMlbMass;
        
//         if (isCM)      nofCM_gate1 += lumiWeight*scaleFactor*widthSF;
//         else if (isWM) nofWM_gate1 += lumiWeight*scaleFactor*widthSF;
//         else if (isUM) nofUM_gate1 += lumiWeight*scaleFactor*widthSF;
        
        if ( reco_top_mass_aKF > maxTopMass || reco_top_mass_aKF < minTopMass ) continue;
        //if ( reco_top_mass_aKF - reco_top_mass_bKF > 25 || reco_top_mass_aKF - reco_top_mass_bKF < -25 ) continue;
        //if ( reco_top_pt_aKF < 100. ) continue;
        //if ( selectedJetsKFcorrected[0].Pt() < 25. || selectedJetsKFcorrected[1].Pt() < 25. ) continue;
        
//        if ( fabs(reco_top_mass_aKF - reco_top_mass_alt_aKF) < 20. ) continue;
        
        if ( reco_top_mass_alt_aKF < maxTopMass-5 ) continue;
        if ( reco_mlb_aKF > maxTopMass-5 ) continue;
        
        //if ( reco_mbjj_div_mjj_bKF < 1.5 || reco_mbjj_div_mjj_bKF > 2.6 ) continue;
        
//         if (isCM)      nofCM_gate2 += lumiWeight*scaleFactor*widthSF;
//         else if (isWM) nofWM_gate2 += lumiWeight*scaleFactor*widthSF;
//         else if (isUM) nofUM_gate2 += lumiWeight*scaleFactor*widthSF;
        
        nofAfterLastCut++;
        
        if (calculateAverageMassAllMC && ! isData)
        {
          sumEvents  += lumiWeight*scaleFactor*widthSF;
          sumTopMass += reco_top_mass_aKF*lumiWeight*scaleFactor*widthSF;
          //sumTopMass += reco_top_mass_bKF*lumiWeight*scaleFactor*widthSF;
          sumMlb += reco_mlb_aKF*lumiWeight*scaleFactor*widthSF;
          if (isCM)
          {
            sumCMEvents += lumiWeight*scaleFactor*widthSF;
            sumCMTopMass += reco_top_mass_aKF*lumiWeight*scaleFactor*widthSF;
            sumCMMlb += reco_mlb_aKF*lumiWeight*scaleFactor*widthSF;
          }
        }
        
        if ( reco_top_mass_aKF < 0. )
          PrintKFDebug(ievt);
        
        if (calculateAverageMass) txtMassReco << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
        
        if ( doKinFit && makePlots )
        {
          if (isTTbar) FillKinFitPlots(doneKinFit);
          if (! isData)
          {
            histo1D["allSim_top_mass"]->Fill(reco_top_mass_aKF, lumiWeight*scaleFactor*widthSF);
            histo1D["allSim_red_top_mass_old"]->Fill(redTopMass_old, lumiWeight*scaleFactor*widthSF);
            histo1D["allSim_red_top_mass_new"]->Fill(redTopMass_new, lumiWeight*scaleFactor*widthSF);
            histo1D["allSim_red_top_mass"]->Fill(redTopMass, lumiWeight*scaleFactor*widthSF);
            histo1D["allSim_mass_bjj_div_m_jj"]->Fill(reco_mbjj_div_mjj, lumiWeight*scaleFactor*widthSF);
          }
        }
        
        
        
        ////////////////////
        ///  Likelihood  ///
        ////////////////////
        
        if ( isWM && (makeTGraphs || calculateLikelihood || calculateFractions) )
        {
          isWM = false;
          isUM = true;
          catSuffix = catSuffixList[2];
        }
        
        if (makeTGraphs && ! isData)
        {
          if (useNewVar)
          {
            if (doLikeW)       likeW->FillHistograms(reco_new_var, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isData, catSuffix);
            else if (doLikeM)  likeM->FillHistograms(reco_new_var, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isData, catSuffix);
            else if (doLike2D) like2D->FillHistograms(reco_new_var, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isData, catSuffix);
          }
          else
          {
            if (doLikeW)       likeW->FillHistograms(redTopMass, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isData, catSuffix);
            else if (doLikeM)  likeM->FillHistograms(redTopMass, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isData, catSuffix);
            else if (doLike2D) like2D->FillHistograms(redTopMass, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isData, catSuffix);
          }
        }
        
        if (calculateLikelihood)
        {
          if (useNewVar)
          {
            if (doLikeW)       loglike_per_evt = likeW->CalculateLikelihood(reco_new_var, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
            else if (doLikeM)  loglike_per_evt = likeM->CalculateLikelihood(reco_new_var, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
            else if (doLike2D) loglike_per_evt = like2D->CalculateLikelihood(reco_new_var, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
          }
          else
          {
            if (doLikeW)       loglike_per_evt = likeW->CalculateLikelihood(redTopMass, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
            else if (doLikeM)  loglike_per_evt = likeM->CalculateLikelihood(redTopMass, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
            else if (doLike2D) loglike_per_evt = like2D->CalculateLikelihood(redTopMass, relativeSF*scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
          }
          if ( ! doPseudoExps && ! useTTTemplates && ! runSystematics && isCM )  // isCM ensures ! isData
          {
            if (useNewVar)
            {
              if (doLikeW)       likeW->CalculateCMLikelihood(reco_new_var, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
              else if (doLikeM)  likeM->CalculateCMLikelihood(reco_new_var, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
              else if (doLike2D) like2D->CalculateCMLikelihood(reco_new_var, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
            }
            else
            {
              if (doLikeW)       likeW->CalculateCMLikelihood(redTopMass, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
              else if (doLikeM)  likeM->CalculateCMLikelihood(redTopMass, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
              else if (doLike2D) like2D->CalculateCMLikelihood(redTopMass, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
            }
          }
//            if ( ! doPseudoExps && ! useTTTemplates && ! runSystematics && ( isCM || isWM ) )
//            {
//              if (useNewVar) like->CalculateCMLikelihood(reco_new_var, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
//              else like->CalculateTempLikelihood(redTopMass, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
//            }
        }
        
        /// Calculate fraction of events in category outside interval
        if (! isData)
        {
          if (isCM) nofCMout_weighted += lumiWeight*scaleFactor*widthSF;
          else if (isWM) nofWMout_weighted += lumiWeight*scaleFactor*widthSF;
          else if (isUM)
          {
            nofUMout_weighted += lumiWeight*scaleFactor*widthSF;
            if (isTTsemilep) nofUMout_TTsemilep_weighted += lumiWeight*scaleFactor*widthSF;
            else if (isTTother) nofUMout_TTother_weighted += lumiWeight*scaleFactor*widthSF;
            else nofUMout_other_weighted += lumiWeight*scaleFactor*widthSF;
          }
        }
        
        if ( redTopMass > maxRedTopMass ) maxRedTopMass = redTopMass;
        if ( redTopMass < minRedTopMass ) minRedTopMass = redTopMass;
        if (! isData)
        {
          if ( (useNewVar && reco_new_var > minCutRedTopMass && reco_new_var < maxCutRedTopMass) 
              || (! useNewVar && redTopMass > minCutRedTopMass && redTopMass < maxCutRedTopMass) )
          {
            if (isCM)
            {
              nofCM++;
              nofCMl += lumiWeight*scaleFactor;
              nofCM_weighted += lumiWeight*scaleFactor*widthSF;
              if (isTTbar) nofCM_TT++;
            }
            else if (isWM)
            {
              nofWM++;
              nofWMl += lumiWeight*scaleFactor;
              nofWM_weighted += lumiWeight*scaleFactor*widthSF;
              if (isTTbar) nofWM_TT++;
            }
            else if (isUM)
            {
              nofUM++;
              nofUMl += lumiWeight*scaleFactor;
              nofUM_weighted += lumiWeight*scaleFactor*widthSF;
              if (isTTbar) nofUM_TT++;
              if (isTTsemilep) nofUM_TTsemilep_weighted += lumiWeight*scaleFactor*widthSF;
              else if (isTTother) nofUM_TTother_weighted += lumiWeight*scaleFactor*widthSF;
              else nofUM_other_weighted += lumiWeight*scaleFactor*widthSF;
            }
            
            if (calculateFractions)
            {
              if (doLikeW)       likeW->AddToFraction(d, lumiWeight*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isCM, isWM, isUM);
              else if (doLikeM)  likeM->AddToFraction(d, lumiWeight*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isCM, isWM, isUM);
              else if (doLike2D) like2D->AddToFraction(d, lumiWeight*scaleFactor, massHadTopQ, massLepTopQ, isTTbar, isCM, isWM, isUM);
            }
            
//             if (isTTbar && makeLikelihoodPlots)
//             {
//               FillLikelihoodPlots();
//             }
          }
        }
        
        if (calculateAverageMass && ! isData)
        {
          if (isCM) txtMassRecoCM << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
          else
          {
            txtMassRecoWMUM << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
            if (isWM)
              txtMassRecoWM << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
            else if (isUM)
              txtMassRecoUM << ievt << "  " << reco_top_mass_aKF << "  " << scaleFactor << "  " << lumiWeight << endl;
          }
        }  // end aveMassCalc
        
        
        /// Pseudo experiments
        if (doPseudoExps && ! isData)
        {
          random3.RndmArray(nPsExps, toyValues);
          for (int iPsExp = 0; iPsExp < nPsExps; iPsExp++)
          {
            //toyValue = random3.Rndm();
            //if ( toyValue > toyMax ) continue;
            if ( toyValues[iPsExp] > toyMax ) continue;
            (nEvtsInPseudoExp[iPsExp][d])++;
            nEvtsInPseudoExpW[iPsExp][d] += scaleFactor;
            if ( (useNewVar && reco_new_var > minCutRedTopMass && reco_new_var < maxCutRedTopMass) 
              || (! useNewVar && redTopMass > minCutRedTopMass && redTopMass < maxCutRedTopMass) )
            {
              if (doLikeW)       likeW->AddPsExp(iPsExp, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
              else if (doLikeM)  likeM->AddPsExp(iPsExp, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
              else if (doLike2D) like2D->AddPsExp(iPsExp, scaleFactor, massHadTopQ, massLepTopQ, thisWidth, thisMass, doReweighting, isData);
            }
            
            /// Fill plots only for first pseudo experiment
            if ( makePlots && iPsExp == 0 )
            {
              /// Combine DY & W+jets
              if ( dataSetName.find("DY") != std::string::npos )
              {
                histo1D["red_top_mass_DYJets"]->Fill(redTopMass);
                histo2D["KF_top_mass_orig_vs_corr_DYJets"]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
              }
              else if ( dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos )
              {
                histo1D["red_top_mass_WJets"]->Fill(redTopMass);
                histo2D["KF_top_mass_orig_vs_corr_WJets"]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
              }
              else
              {
                histo1D[("red_top_mass_"+dataSetName).c_str()]->Fill(redTopMass);
                histo2D["KF_top_mass_orig_vs_corr_"+dataSetName]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
              }
              
              FillCatsPlots(catSuffix);
              
//               int dMSP = d;
//               if (hasFoundTTbar && ! isTTbar) dMSP = d+3;
//               else if (isTTbar && isWM) dMSP = d+1;
//               else if (isTTbar && isUM_TTsemilep) dMSP = d+2;
//               else if (isTTbar && isUM_TTother) dMSP = d+3;
              
              FillMSPlots(dMSP, doneKinFit);
            }
          }
        }
        
        //Fill histos
        if ( makePlots && (! doPseudoExps || isData) )
        {
          /// Combine DY & W+jets
          if ( dataSetName.find("DY") != std::string::npos )
          {
            histo1D["red_top_mass_DYJets"]->Fill(redTopMass);
            histo2D["KF_top_mass_orig_vs_corr_DYJets"]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
          }
          else if ( dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos )
          {
            histo1D["red_top_mass_WJets"]->Fill(redTopMass);
            histo2D["KF_top_mass_orig_vs_corr_WJets"]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
          }
          else
          {
            histo1D[("red_top_mass_"+dataSetName).c_str()]->Fill(redTopMass);
            histo2D["KF_top_mass_orig_vs_corr_"+dataSetName]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
          }
          
          FillCatsPlots(catSuffix);
          
          FillMSPlots(dMSP, doneKinFit);
          
          MSPlot["leadingJet_pT_aKF_"]->Fill(selectedJetsAKF[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
          for (int iJet = 0; iJet < selectedJetsAKF.size(); iJet++)
          {
            MSPlot["jet_pT_allJets_aKF_"]->Fill(selectedJetsAKF[iJet].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
          }
          
        }  // end makePlots
        
        
        
        
      }  // end loop events
      
      
      if (calculateAverageMassAllMC)
      {
        if ( d == datasets.size()-1 )
        {
          aveTopMassLL = sumTopMass/sumEvents;
          cout << endl << "Average top mass for all MC is    " << aveTopMassLL << endl;
          //aveTopMassCM = sumCMTopMass/sumCMEvents;
          cout << endl << "Average top mass for CM events is " << aveTopMassCM << endl;
          
          aveMlbMass = sumMlb/sumEvents;
          aveMlbMassCM = sumCMMlb/sumCMEvents;
          cout << endl << "Average M_{lb} mass is " << aveMlbMass << " and for CM events " << aveMlbMassCM << endl;
          
          /// Rerun with new average top mass
          if ( saveSetUp[0] || saveSetUp[1] || saveSetUp[2] || saveSetUp[3] ||saveSetUp[4] || saveSetUp[5] || saveSetUp[6] || saveSetUp[7] )
          {
            d = -1;
            calculateAverageMassAllMC = false;
            //  Reset booleans
            makePlots = saveSetUp[0];
            makeControlPlots = saveSetUp[1];
            makeLikelihoodPlots = saveSetUp[2];
            calculateFractions = saveSetUp[3];
            makeTGraphs = saveSetUp[4];
            useTTTemplates = saveSetUp[5];
            calculateLikelihood = saveSetUp[6];
            doPseudoExps = saveSetUp[7];
            //  Set up environment with new booleans
            InitSetUp();
            
            nofCM = 0; nofWM = 0; nofUM = 0;
            nofCM_TT = 0; nofWM_TT = 0; nofUM_TT = 0;
            nofCMl = 0; nofWMl = 0; nofUMl = 0;
            nofCM_weighted = 0.; nofWM_weighted = 0.; nofUM_weighted = 0.; nofUM_TTsemilep_weighted = 0.; nofUM_TTother_weighted = 0.; nofUM_other_weighted = 0.;
            nofCMout_weighted = 0.; nofWMout_weighted = 0.; nofUMout_weighted = 0.; nofUMout_TTsemilep_weighted = 0.; nofUMout_TTother_weighted = 0.; nofUMout_other_weighted = 0.;
            nofCM_gate1 = 0.; nofWM_gate1 = 0.; nofUM_gate1 = 0.; nofCM_gate2 = 0.; nofWM_gate2 = 0.; nofUM_gate2 = 0.;
            
            hasFoundTTbar = false;
            doReweighting = false;
            
            continue;
          }
        }
        else continue;
      }
      
      cout << endl;  /// Stronger selection in this analyser compared to Ntuples ==> endEvent --> nofHardSelected
      cout << "Number of events with exactly 4 jets with pT > 30 GeV: " << nofHardSelected << " (" << 100*((float)nofHardSelected/(float)endEvent) << "%)" << endl;
      cout << "Number of events with clean MET: " << nofMETCleaned << " (" << 100*((float)nofMETCleaned/(float)nofHardSelected) << "%)" << endl;
      cout << "Number of events after min dR cut: " << nofAfterDRmincut << " (" << 100*((float)nofAfterDRmincut/(float)nofMETCleaned) << "%)" << endl;
      cout << "Number of events after max dR cut: " << nofAfterDRmaxcut << " (" << 100*((float)nofAfterDRmaxcut/(float)nofAfterDRmincut) << "%)" << endl;
      if (doKinFit) cout << "Number of events accepted by kinFitter: " << nofAcceptedKFit << " (" << 100*((float)nofAcceptedKFit/(float)nofAfterDRmaxcut) << "%)" << endl;
      if (doKinFit) cout << "Number of events selected after last cut: " << nofAfterLastCut << " (" << 100*((float)nofAfterLastCut/(float)nofAcceptedKFit) << "%)" << endl;
      else cout << "Number of events selected after last cut: " << nofAfterLastCut << " (" << 100*((float)nofAfterLastCut/(float)nofAfterDRmaxcut) << "%)" << endl;
      
      if (isTTbar)
      {
        cout << "Number of semileptonic events: " << setw(8) << right << nofTTsemilep << " (" << 100*((float)nofTTsemilep/(float)(nofTTsemilep+nofTTdilep+nofTThadr)) << "%)" << endl;
        cout << "Number of dileptonic events:   " << setw(8) << right << nofTTdilep << " (" << 100*((float)nofTTdilep/(float)(nofTTsemilep+nofTTdilep+nofTThadr)) << "%)" << endl;
        cout << "Number of all-hadronic events: " << setw(8) << right << nofTThadr << " (" << 100*((float)nofTThadr/(float)(nofTTsemilep+nofTTdilep+nofTThadr)) << "%)" << endl;
      }
      
      //if ( isTTbar || dataSetName.find("ST") != std::string::npos )
      if (! isData && nofHadrMatchedEvents > 0 )
      {
        cout << "Number of matched events: " << setw(8) << right << nofMatchedEvents << endl;
        cout << "Number of events with hadronic top matched (before KF): " << setw(8) << right << nofHadrMatchedEvents << " (" << 100*((float)nofHadrMatchedEvents/(float)nofAfterDRmaxcut) << "%)" << endl;
        if (doKinFit) cout << "Number of events with hadronic top matched (after KF):  " << setw(8) << right << nofHadrMatchedEventsAKF << " (" << 100*((float)nofHadrMatchedEventsAKF/(float)nofAfterLastCut) << "%)" << endl;
        if (! doGenOnly)
        {
          cout << "Correctly matched reconstructed events:     " << setw(8) << right << nofCorrectlyMatched << endl;
          cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatched << endl;
          cout << "Unmatched reconstructed events            : " << setw(8) << right << nofUnmatched << endl;
          cout << "                      whereof TT semilep  : " << setw(8) << right << nofUnmatchedTTsemilep << endl;
          cout << "                              TT other    : " << setw(8) << right << nofUnmatchedTTother << endl;
          if ( nofCorrectlyMatched != 0 || nofNotCorrectlyMatched != 0 )
            cout << "   ===> This means that " << 100*(float)nofCorrectlyMatched / (float)(nofCorrectlyMatched + nofNotCorrectlyMatched) << "% of matched events is correctly matched." << endl;
          
          if (doKinFit)
          {
            cout << "                        " << 100*(float)nofCorrectlyMatched / (float)nofAfterDRmaxcut << "% of all events is correctly matched before kinfitter." << endl;
            cout << " --- Kinematic fit --- Before chi2 cut --- " << endl;
            cout << "Correctly matched reconstructed events    : " << setw(8) << right << nofCorrectlyMatchedAKFNoCut << endl;
            cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatchedAKFNoCut << endl;
            cout << "Unmatched reconstructed events            : " << setw(8) << right << nofUnmatchedAKFNoCut << endl;
            cout << "                      whereof TT semilep  : " << setw(8) << right << nofUnmatchedTTsemilepAKFNoCut << endl;
            cout << "                              TT other    : " << setw(8) << right << nofUnmatchedTTotherAKFNoCut << endl;
            if ( nofCorrectlyMatchedAKFNoCut != 0 || nofNotCorrectlyMatchedAKFNoCut != 0 )
              cout << "   ===> This means that " << 100*(float)nofCorrectlyMatchedAKFNoCut / (float)(nofCorrectlyMatchedAKFNoCut + nofNotCorrectlyMatchedAKFNoCut) << "% of matched events is correctly matched after KF." << endl;
            
            cout << "                        " << 100*(float)nofCorrectlyMatchedAKFNoCut / (float)(nofCorrectlyMatchedAKFNoCut+nofNotCorrectlyMatchedAKFNoCut+nofUnmatchedAKFNoCut) << "% of all events accepted by kinfitter is correctly matched." << endl;
            
            cout << " --- Kinematic fit --- After chi2 cut --- " << endl;
            cout << "Correctly matched reconstructed events (after KF): " << setw(8) << right << nofCorrectlyMatchedAKF << endl;
            cout << "Not correctly matched reconstructed events: " << setw(8) << right << nofNotCorrectlyMatchedAKF << endl;
            cout << "Unmatched reconstructed events            : " << setw(8) << right << nofUnmatchedAKF << endl;
            cout << "                      whereof TT semilep  : " << setw(8) << right << nofUnmatchedTTsemilepAKF << endl;
            cout << "                              TT other    : " << setw(8) << right << nofUnmatchedTTotherAKF << endl;
            if ( nofCorrectlyMatchedAKF != 0 || nofNotCorrectlyMatchedAKF != 0 )
              cout << "   ===> This means that " << 100*(float)nofCorrectlyMatchedAKF / (float)(nofCorrectlyMatchedAKF + nofNotCorrectlyMatchedAKF) << "% of matched events is correctly matched after KF." << endl;
            
            cout << "                        " << 100*(float)nofCorrectlyMatchedAKF / (float)nofAcceptedKFit << "% of all events accepted by kinfitter is correctly matched." << endl;
          }
          else cout << "                        " << 100*(float)nofCorrectlyMatched / (float)nofAfterLastCut << "% of all events is correctly matched." << endl;
        }
        
        
        /// Resolution functions
        if (isTTbar && calculateResolutionFunctions)
        {
          string rfFileName = "PlotsForResolutionFunctions.root";
          string rfFitFileName = "PlotsForResolutionFunctions_Fitted.root";
          TFile *foutRF = new TFile(rfFileName.c_str(), "RECREATE");
          foutRF->cd();
          
          rf->writeHistograms();
          
          foutRF->Close();
          
          rf->makeFit(rfFileName, rfFitFileName);
          rf->writeTable(rfFitFileName);
          
          delete foutRF;
        }
        
      }  // end ! isData
      
      
      /// Make selection table
      if (! runListWidths && ! runSystematics && ! calculateResolutionFunctions && ! calculateAverageMass && systStr.find("nominal") != std::string::npos )
        FillSelectionTable(d, dataSetName);
      
      
      if (calculateAverageMass) txtMassReco.close();
      
      tFileMap[dataSetName.c_str()]->Close();
      
      timePerDataSet[d] = ((double)clock() - startDataSet) / CLOCKS_PER_SEC;
      
    }  // end loop datasets
    
    if (! runListWidths && ! runSystematics && systStr.find("nominal") != std::string::npos ) WriteSelectionTable();
    
    if (! doGenOnly && ! testTTbarOnly)
    {
      cout << "Number of events with " << minCutRedTopMass << " < mt/<mt> < " << maxCutRedTopMass << " : CM: " << nofCM << " (" << 100*(double)nofCM/((double)(nofCM+nofWM+nofUM)) << "%)   WM: " << nofWM << " (" << 100*(double)nofWM/((double)(nofCM+nofWM+nofUM)) << "%)   UM: " << nofUM << " (" << 100*(double)nofUM/((double)(nofCM+nofWM+nofUM)) << "%)   Total: " << nofCM+nofWM+nofUM << endl;
      cout << "Number of events with " << minCutRedTopMass << " < mt/<mt> < " << maxCutRedTopMass << " : CM: " << nofCMl << " (" << 100*nofCMl/(nofCMl+nofWMl+nofUMl) << "%)   WM: " << nofWMl << " (" << 100*nofWMl/(nofCMl+nofWMl+nofUMl) << "%)   UM: " << nofUMl << " (" << 100*nofUMl/(nofCMl+nofWMl+nofUMl) << "%)   Total: " << nofCMl+nofWMl+nofUMl << endl;
      cout << "                                  weighted: CM: " << nofCM_weighted << " (" << 100*nofCM_weighted/(nofCM_weighted+nofWM_weighted+nofUM_weighted) << "%)   WM: " << nofWM_weighted << " (" << 100*nofWM_weighted/(nofCM_weighted+nofWM_weighted+nofUM_weighted) << "%)   UM: " << nofUM_weighted << " (" << 100*nofUM_weighted/(nofCM_weighted+nofWM_weighted+nofUM_weighted) << "%),  whereof TT semilep UM: " << nofUM_TTsemilep_weighted << " (" << 100*nofUM_TTsemilep_weighted/(nofCM_weighted+nofWM_weighted+nofUM_weighted) << "%),  and TT other UM: " << nofUM_TTother_weighted << " (" << 100*nofUM_TTother_weighted/(nofCM_weighted+nofWM_weighted+nofUM_weighted) << "%),  and other UM: " << nofUM_other_weighted << " (" << 100*nofUM_other_weighted/(nofCM_weighted+nofWM_weighted+nofUM_weighted) << "%)   Total: " << (int)(nofCM_weighted+nofWM_weighted+nofUM_weighted) << endl;
      cout << "                               (TTbar only) CM: " << nofCM_TT << "               WM: " << nofWM_TT << "               UM: " << nofUM_TT << endl;
      cout << endl << "Number of events outside interval:          CM: " << nofCMout_weighted << " (" << 100*nofCMout_weighted/(nofCMout_weighted+nofWMout_weighted+nofUMout_weighted) << "%)   WM: " << nofWMout_weighted << " (" << 100*nofWMout_weighted/(nofCMout_weighted+nofWMout_weighted+nofUMout_weighted) << "%)   UM: " << nofUMout_weighted << " (" << 100*nofUMout_weighted/(nofCMout_weighted+nofWMout_weighted+nofUMout_weighted) << "%),  whereof TT semilep UM: " << nofUMout_TTsemilep_weighted << " (" << 100*nofUMout_TTsemilep_weighted/(nofCMout_weighted+nofWMout_weighted+nofUMout_weighted) << "%),  and TT other UM: " << nofUMout_TTother_weighted << " (" << 100*nofUMout_TTother_weighted/(nofCMout_weighted+nofWMout_weighted+nofUMout_weighted) 
          << "%),  and other UM: " << nofUMout_other_weighted << " (" << 100*nofUMout_other_weighted/(nofCMout_weighted+nofWMout_weighted+nofUMout_weighted) << "%)   Total: " << (int)(nofCMout_weighted+nofWMout_weighted+nofUMout_weighted) << endl;
    }
    
    if (! doGenOnly)
    {
      cout << endl;
      cout << "Number of events before cut:  CM: " << nofCM_gate1 << " (" << 100*nofCM_gate1/(nofCM_gate1+nofWM_gate1+nofUM_gate1) << "%)   WM: " << nofWM_gate1 << " (" << 100*nofWM_gate1/(nofCM_gate1+nofWM_gate1+nofUM_gate1) << "%)   UM: " << nofUM_gate1 << " (" << 100*nofUM_gate1/(nofCM_gate1+nofWM_gate1+nofUM_gate1) << "%)   Total: " << nofCM_gate1+nofWM_gate1+nofUM_gate1 << endl;
      cout << "Number of events after cut :  CM: " << nofCM_gate2 << " (" << 100*nofCM_gate2/(nofCM_gate2+nofWM_gate2+nofUM_gate2) << "%)   WM: " << nofWM_gate2 << " (" << 100*nofWM_gate2/(nofCM_gate2+nofWM_gate2+nofUM_gate2) << "%)   UM: " << nofUM_gate2 << " (" << 100*nofUM_gate2/(nofCM_gate2+nofWM_gate2+nofUM_gate2) << "%)   Total: " << nofCM_gate2+nofWM_gate2+nofUM_gate2 << endl;
    }
    
    
    if (calculateFractions)
    {
      if (doLikeW)       likeW->CalculateFractions(dataSetNames);
      else if (doLikeM)  likeM->CalculateFractions(dataSetNames);
      else if (doLike2D) like2D->CalculateFractions(dataSetNames);
    }
    if (makeTGraphs)
    {
      if (doLikeW)
      {
        likeW->WriteHistograms("ReducedTopMassPlots.root");
        likeW->ConstructTGraphsFromHisto("TGraphFunctions.root", dataSetNames, includeDataSets);
      }
      else if (doLikeM)
      {
        likeM->WriteHistograms("ReducedTopMassPlots.root");
        likeM->ConstructTGraphsFromHisto("TGraphFunctions.root", dataSetNames, includeDataSets);
      }
      else if (doLike2D)
      {
        like2D->WriteHistograms("ReducedTopMassPlots.root");
        like2D->ConstructTGraphsFromHisto("TGraphFunctions.root", dataSetNames, includeDataSets);
      }
    }
    
    if (calculateLikelihood)
    {
      cout << "Minimum reduced top mass: " << minRedTopMass << endl;
      cout << "Maximum reduced top mass: " << maxRedTopMass << endl;
      
      /// Print output to file  /// TEMPORARILY !
//       string llFileName = "output_loglikelihood_widthx"+DotReplace(thisWidth);
//       if (doGenOnly) llFileName = "output_loglikelihood_parton_widthx"+DotReplace(thisWidth);
//       //if (useToys) llFileName = "output_loglikelihood_toys";
//       if (runSystematics) llFileName = "output_loglikelihood_"+thisSystematic;
//       like->PrintLikelihoodOutput(llFileName+".txt");
//       if (unblind) like->PrintLikelihoodOutputData(llFileName+"_data.txt");
//       like->PrintMtmLikelihoodOutput(llFileName+"_Mtm.txt");
      
      /// Calculate output width
      if ( runSystematics || runListWidths)
      {
        if (runSystematics) cout << "Output width for " << thisSystematic << ": " << endl;
        else cout << endl << "Standard output width: " << endl;
        if (runSystematics)
        {
          if (runRateSystematics || runSampleSystematics)
          {
            fileWidths->cd();
            if (doLikeW)       likeW->GetOutputWidth(thisWidth, thisMass, thisSystematic, true, false);
            else if (doLikeM)  likeM->GetOutputMass(thisWidth, thisMass, thisSystematic, true, false);
            else if (doLike2D) like2D->GetOutputWidth(thisWidth, thisMass, thisSystematic, true, false);
          }
          else
          {
            if (doLikeW)       likeW->GetOutputWidth(thisWidth, thisMass, thisSystematic, true, true);
            else if (doLikeM)  likeM->GetOutputMass(thisWidth, thisMass, thisSystematic, true, true);
            else if (doLike2D) like2D->GetOutputWidth(thisWidth, thisMass, thisSystematic, true, true);
          }
        }
        else
        {
          fileWidths->cd();
          if (doLikeW)       likeW->GetOutputWidth(thisWidth, thisMass, true, false);
          else if (doLikeM)  likeM->GetOutputMass(thisWidth, thisMass, true, false);
          else if (doLike2D) like2D->GetOutputWidth(thisWidth, thisMass, true, false);
        }
//         if (! useTTTemplates && ! runSystematics)
//         {
//           cout << "Output width for correctly matched events (using likelihood with only CM template): " << endl;
//           like->GetOutputWidth(thisWidth, thisMass, "CM", true, false);
//           cout << "Output width for correctly & wrongly matched events (using likelihood with only CM & WM templates): " << endl;
//           like->GetOutputWidth(thisWidth, thisMass, "matched", true, false);
//         }
      }
      else
      {
        cout << endl << "Standard output width: " << endl;
        if (doLikeW)
        {
          likeW->GetOutputWidth(thisWidth, thisMass, systStr, true, true);
          if (unblind) likeW->GetOutputWidth(thisWidth, thisMass, "data", true, true);
        }
        else if (doLikeM)
        {
          likeM->GetOutputMass(thisWidth, thisMass, systStr, true, true);
          if (unblind) likeM->GetOutputMass(thisWidth, thisMass, "data", true, true);
        }
        else if (doLike2D)
        {
          like2D->Make2DGraph("2D_"+systStr, true);  // fileName, makeNewFile
          like2D->GetOutputWidth(thisWidth, thisMass, systStr, true, true);
          if (unblind) like2D->GetOutputWidth(thisWidth, thisMass, "data", true, true);
        }
        if (! doPseudoExps && ! useTTTemplates)
        {
          cout << endl << "Output width for correctly matched events (using likelihood with only CM template): " << endl;
          if (doLikeW)       likeW->GetOutputWidth(thisWidth, thisMass, "CM", true, true);
          else if (doLikeM)  likeM->GetOutputMass(thisWidth, thisMass, "CM", true, true);
          else if (doLike2D) like2D->GetOutputWidth(thisWidth, thisMass, "CM", true, true);
//          cout << "Output width for correctly & wrongly matched events (using likelihood with only CM & WM templates): " << endl;
//          like->GetOutputWidth(thisWidth, "matched", true, true);
//           //cout << "Output width for generated events (using likelihood with only CM template): " << endl;
//           //like->GetOutputWidth(scaleWidth, "gen", true);
        }
//         //cout << "Output width from file (standard calculation): " << endl;
//         //like->GetOutputWidth(llFileName+".txt", scaleWidth, true);
      }
    }
    
    if (doPseudoExps)
    {
      if (doLikeW)       likeW->CalculatePull(thisWidth, thisMass);
      else if (doLikeM)  likeM->CalculatePull(thisWidth, thisMass);
      else if (doLike2D) like2D->CalculatePull(thisWidth, thisMass);
    }
    
    
    cout << endl << "Processing time per dataset: " << endl;
    for (unsigned int d = 0; d < datasets.size(); d++)
    {
      cout << datasets[d]->Name() << ": " << timePerDataSet[d] << " s" << endl;
    }
    
    if (doPseudoExps)
    {
      cout << endl << "PseudoExperiments::Number of selected data events: " << nDataEvts << endl;
      double totMCPsExp, totMCPsExpW;
      for (int i = 0; i < 3; i++)
      {
        totMCPsExp = 0;
        totMCPsExpW = 0.;
        cout << "PseudoExperiment " << std::setw(3) << std::right << i << endl;
        for (unsigned int d = 1; d < datasets.size(); d++)
        {
          totMCPsExp += nEvtsInPseudoExp[i][d];
          totMCPsExpW += nEvtsInPseudoExpW[i][d];
          cout << "                      " << datasets[d]->Name() << ": " << nEvtsInPseudoExp[i][d] << endl;
        }
        cout << "                                            Total MC: " << totMCPsExp << endl;
        cout << "                                           (weighted: " << totMCPsExpW << ")" << endl;
      }
    }
    
    
    ///  Check Shape Changing Systematics
    if (! testTTbarOnly && runAll && ! runListWidths && ! runSystematics && ! useTTTemplates && ! doPseudoExps && systStr.find("nominal") != std::string::npos)
      CheckSystematics(vJER, vJES, vPU);
    
    
  }  // end loop systematics/widths
  
  if ( calculateLikelihood && (runRateSystematics || runSampleSystematics || runListWidths) )
  {
    fileWidths->Close();
    delete fileWidths;
  }
  
  if (calculateAverageMass)
  {
    txtMassGenPMatched.close();
    txtMassGenJMatched.close();
    txtMassRecoCM.close();
    txtMassRecoWMUM.close();
    txtMassRecoUM.close();
    txtMassRecoWM.close();
  }
  
  if (applyWidthSF) txtDebugTopMass.close();
  
  
  
  if (test)
  {
    cout << "Exiting because of test..." << endl;
    exit(1);
  }
  else if (calculateAverageMass)
  {
    cout << "Average mass calculated. Exiting..." << endl;
    exit(1);
  }
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  if (makePlots)
  {
    string rootFileName = "NtuplePlots_"+systStr+".root";
    if (! doGenOnly)
    {
      mkdir((pathOutput+"MSPlot/").c_str(),0777);
      if (makeControlPlots) mkdir((pathOutput+"MSPlot/ControlPlots/").c_str(),0777);
    }
    
    cout << " - Recreate output file ..." << endl;
    TFile *fout = new TFile ((pathOutput+rootFileName).c_str(), "RECREATE");
    cout << "   Output file is " << pathOutput+rootFileName << endl;
    
    ///Write histograms
    fout->cd();
    
    if (! doGenOnly)
    {
      for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
      {
        //cout << "MSPlot: " << it->first << endl;
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        temp->Draw(name, 1, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
        temp->Write(fout, name, true, pathOutput+"MSPlot", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
      }
      
      for (map<string,MultiSamplePlot*>::const_iterator it = MSPlotCP.begin(); it != MSPlotCP.end(); it++)
      {
        //cout << "MSPlot: " << it->first << endl;
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        temp->Draw(name, 1, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
        temp->Write(fout, name, true, pathOutput+"MSPlot/ControlPlots", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
      }
    }
    
    // 1D
    TDirectory* th1dir = fout->mkdir("1D_histograms");
    th1dir->cd();
    gStyle->SetOptStat(1111);
    for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
      TH1F *temp = it->second;
      int N = temp->GetNbinsX();
      temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
      temp->SetBinContent(N+1,0);
      temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
      temp->Write();
      TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
      tempCanvas->SaveAs( (pathOutput+it->first+".png").c_str() );
    }
    
    // 2D
    TDirectory* th2dir = fout->mkdir("2D_histograms");
    th2dir->cd();
    gStyle->SetPalette(55);
    for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {
      TH2F *temp = it->second;
      temp->Write();
      TCanvas* tempCanvas = TCanvasCreator(temp, it->first, "colz");
      tempCanvas->SaveAs( (pathOutput+it->first+".png").c_str() );
    }
    
    fout->Close();
    
    delete fout;
  }
  if (makeLikelihoodPlots) WriteLikelihoodPlots();
  
  
  double time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  if ( time >= 60 )
  {
    int mins = time/60;
    double secs = time - mins*60;
    
    if (mins >= 60 )
    {
      int hours = mins/60;
      mins = mins - hours*60;
      cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
    }
    else
      cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
  }
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  cout << " - Goodbye" << endl;
  
  return 0;
  
}  // end main



/// Functions

string ConvertDoubleToString(double Number)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

string DotReplace(double var)
{
  string str = ConvertDoubleToString(var);
  replace(str.begin(), str.end(), '.', 'p');
  return str;
}

string ConvertIntToString(int Number, int pad)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  if ( pad > 1 ) { convert << std::setw(pad) << std::setfill('0');}
  convert << Number;
  return convert.str();
}

string MakeTimeStamp()
{
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  int year = now->tm_year - 100;  /// + 1900 to get current year
  int month = now->tm_mon + 1;
  int day = now->tm_mday;
  int hour = now->tm_hour;
  int min = now->tm_min;
  //int sec = now->tm_sec;
  
  string year_str = ConvertIntToString(year, 2);
  string month_str = ConvertIntToString(month, 2);
  string day_str = ConvertIntToString(day, 2);
  string hour_str = ConvertIntToString(hour, 2);
  string min_str = ConvertIntToString(min, 2);
  //string sec_str = ConvertIntToString(sec, 2);
  
  string date_str = year_str + month_str + day_str + "_" + hour_str + min_str;
  return date_str;
}

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile.good();
}

void InitSetUp()
{
  if (makeTGraphs || calculateFractions)
  {
    if (doLikeW)       likeW = new Likelihood(minCutRedTopMass, maxCutRedTopMass, outputDirLL, dateString, rewHadTopOnly, makeTGraphs, true);  // verbose
    else if (doLikeM)  likeM = new LikelihoodMass(minCutRedTopMass, maxCutRedTopMass, outputDirLL, dateString, rewHadTopOnly, makeTGraphs, true);  // verbose
    else if (doLike2D) like2D = new Likelihood2D(minCutRedTopMass, maxCutRedTopMass, outputDirLL, dateString, rewHadTopOnly, makeTGraphs, true);  // verbose
  }
  if (calculateLikelihood)
  {
    if (doLikeW)       likeW = new Likelihood(minCutRedTopMass, maxCutRedTopMass, inputDirLL, dateString, rewHadTopOnly, makeTGraphs, true);  // verbose
    else if (doLikeM)  likeM = new LikelihoodMass(minCutRedTopMass, maxCutRedTopMass, inputDirLL, dateString, rewHadTopOnly, makeTGraphs, true);  // verbose
    else if (doLike2D) like2D = new Likelihood2D(minCutRedTopMass, maxCutRedTopMass, inputDirLL, dateString, rewHadTopOnly, makeTGraphs, true);  // verbose
    
    if (useTTTemplates)
    {
      if (doLikeW)       calculateLikelihood = likeW->ConstructTGraphsFromFile(dataSetNames, includeDataSets);
      else if (doLikeM)  calculateLikelihood = likeM->ConstructTGraphsFromFile(dataSetNames, includeDataSets);
      else if (doLike2D) calculateLikelihood = like2D->ConstructTGraphsFromFile(dataSetNames, includeDataSets);
    }
    else
    {
      if (doLikeW)
      {
        calculateLikelihood = likeW->ConstructTGraphsFromFile();
        if (! runSystematics) calculateLikelihood = likeW->ConstructTGraphsFromFile("CorrectMatchLikelihood_");
//        if (! runSystematics) calculateLikelihood = likeW->ConstructTGraphsFromFile("MatchLikelihood_");
      }
      else if (doLikeM)
      {
        calculateLikelihood = likeM->ConstructTGraphsFromFile();
        if (! runSystematics) calculateLikelihood = likeM->ConstructTGraphsFromFile("CorrectMatchLikelihood_");
//        if (! runSystematics) calculateLikelihood = likeM->ConstructTGraphsFromFile("MatchLikelihood_");
      }
      else if (doLike2D)
      {
        calculateLikelihood = like2D->ConstructTGraphsFromFile();
        if (! runSystematics) calculateLikelihood = like2D->ConstructTGraphsFromFile("CorrectMatchLikelihood_");
      }
    }
    widthsLike.clear();
    massesLike.clear();
    if (doLikeW)       widthsLike = likeW->GetWidths();
    else if (doLikeM)  massesLike = likeM->GetMasses();
    else if (doLike2D)
    {
      widthsLike = like2D->GetWidths();
      massesLike = like2D->GetMasses();
    }
    nWidthsLike = widthsLike.size();
    nMassesLike = massesLike.size();
    
    if (runRateSystematics || runSampleSystematics) fileWidths = new TFile(("OutputLikelihood/"+dateString+"/OutputWidths_syst.root").c_str(), "RECREATE");
    else if (runListWidths) fileWidths = new TFile(("OutputLikelihood/"+dateString+"/OutputWidths.root").c_str(), "RECREATE");
  }
  
  if (doPseudoExps)
  {
    if (doLikeW)       nPsExps = likeW->InitPull(nPseudoExps);
    else if (doLikeM)  nPsExps = likeM->InitPull(nPseudoExps);
    else if (doLike2D) nPsExps = like2D->InitPull(nPseudoExps);
  }
  
  if (makePlots)
  {
    if (! doGenOnly) InitMSPlots();
    InitHisto1D();
    InitHisto2D();
  }
  if (makeLikelihoodPlots)
  {
    InitLikelihoodPlots();
  }
}

void GetMetaData(TTree* tree, bool isData)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  
  tree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
  tree->SetBranchAddress("nEventsSel", &nEventsSel, &b_nEventsSel);
  tree->SetBranchAddress("cutFlow", cutFlow, &b_cutFlow);
  tree->SetBranchAddress("cutFlow2", cutFlow2, &b_cutFlow2);
  tree->SetBranchAddress("cutFlowWeighted", cutFlowWeighted, &b_cutFlowWeighted);
  tree->SetBranchAddress("cutFlow2Weighted", cutFlow2Weighted, &b_cutFlow2Weighted);
  tree->SetBranchAddress("appliedJER", &appliedJER, &b_appliedJER);
  tree->SetBranchAddress("appliedJES", &appliedJES, &b_appliedJES);
  tree->SetBranchAddress("appliedPU", &appliedPU, &b_appliedPU);
  if (isData)
  {
    tree->SetBranchAddress("nofEventsRunB", &nofEventsRunB, &b_nofEventsRunB);
    tree->SetBranchAddress("nofEventsRunCD", &nofEventsRunCD, &b_nofEventsRunCD);
    tree->SetBranchAddress("nofEventsRunEF", &nofEventsRunEF, &b_nofEventsRunEF);
    tree->SetBranchAddress("nofEventsRunG", &nofEventsRunG, &b_nofEventsRunG);
    tree->SetBranchAddress("nofEventsRunH", &nofEventsRunH, &b_nofEventsRunH);
    tree->SetBranchAddress("nofSelEventsRunB", &nofSelEventsRunB, &b_nofEventsRunB);
    tree->SetBranchAddress("nofSelEventsRunCD", &nofSelEventsRunCD, &b_nofSelEventsRunCD);
    tree->SetBranchAddress("nofSelEventsRunEF", &nofSelEventsRunEF, &b_nofSelEventsRunEF);
    tree->SetBranchAddress("nofSelEventsRunG", &nofSelEventsRunG, &b_nofSelEventsRunG);
    tree->SetBranchAddress("nofSelEventsRunH", &nofSelEventsRunH, &b_nofSelEventsRunH);
  }
  else if (isTTbar)
  {
    if (newTrees)
    {
      tree->SetBranchAddress("nofEventsWithGenTop", &nofEventsWithGenTop, &b_nofEventsWithGenTop);
      tree->SetBranchAddress("nofEventsWithGenAntiTop", &nofEventsWithGenAntiTop, &b_nofEventsWithGenAntiTop);
    }
    tree->SetBranchAddress("nofTTEventsWithoutAGenTop", &nofTTEventsWithoutAGenTop, &b_nofTTEventsWithoutAGenTop);
    tree->SetBranchAddress("sumWeight1001", &sumWeight1001, &b_sumWeight1001);
    tree->SetBranchAddress("sumWeight1002", &sumWeight1002, &b_sumWeight1002);
    tree->SetBranchAddress("sumWeight1003", &sumWeight1003, &b_sumWeight1003);
    tree->SetBranchAddress("sumWeight1004", &sumWeight1004, &b_sumWeight1004);
    tree->SetBranchAddress("sumWeight1005", &sumWeight1005, &b_sumWeight1005);
    tree->SetBranchAddress("sumWeight1007", &sumWeight1007, &b_sumWeight1007);
    tree->SetBranchAddress("sumWeight1009", &sumWeight1009, &b_sumWeight1009);
    if (newTrees)
    {
      tree->SetBranchAddress("sumUpFragWeight", &sumUpFragWeight, &b_sumUpFragWeight);
      tree->SetBranchAddress("sumCentralFragWeight", &sumCentralFragWeight, &b_sumCentralFragWeight);
      tree->SetBranchAddress("sumDownFragWeight", &sumDownFragWeight, &b_sumDownFragWeight);
      tree->SetBranchAddress("sumPetersonFragWeight", &sumPetersonFragWeight, &b_sumPetersonFragWeight);
      tree->SetBranchAddress("sumSemilepbrUp", &sumSemilepbrUp, &b_sumSemilepbrUp);
      tree->SetBranchAddress("sumSemilepbrDown", &sumSemilepbrDown, &b_sumSemilepbrDown);
      tree->SetBranchAddress("sumPdfWeights", sumPdfWeights, &b_sumPdfWeights);
      tree->SetBranchAddress("sumPdfAlphaSUp", &sumPdfAlphaSUp, &b_sumPdfAlphaSUp);
      tree->SetBranchAddress("sumPdfAlphaSDown", &sumPdfAlphaSDown, &b_sumPdfAlphaSDown);
    }
  }
}

void InitTree(TTree* tree, bool isData)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  
  tree->SetBranchAddress("run_num", &run_num, &b_run_num);
  tree->SetBranchAddress("evt_num", &evt_num, &b_evt_num);
  tree->SetBranchAddress("lumi_num", &lumi_num, &b_lumi_num);
  tree->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  tree->SetBranchAddress("npu", &npu, &b_npu);
  tree->SetBranchAddress("rho", &rho, &b_rho);
  tree->SetBranchAddress("isTrigged", &isTrigged, &b_isTrigged);
  tree->SetBranchAddress("hasExactly4Jets", &hasExactly4Jets, &b_hasExactly4Jets);
  tree->SetBranchAddress("hasJetLeptonCleaning", &hasJetLeptonCleaning, &b_hasJetLeptonCleaning);
  tree->SetBranchAddress("hasErasedBadOrCloneMuon", &hasErasedBadOrCloneMuon, &b_hasErasedBadOrCloneMuon);
//  tree->SetBranchAddress("passedMETFilter", &passedMETFilter, &b_passedMETFilter);
  tree->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
  tree->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
  tree->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
  tree->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
  tree->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
  tree->SetBranchAddress("muon_E", muon_E, &b_muon_E);
  tree->SetBranchAddress("muon_M", muon_M, &b_muon_M);
  tree->SetBranchAddress("muon_d0", muon_d0, &b_muon_d0);
  tree->SetBranchAddress("muon_chargedHadronIso", muon_chargedHadronIso, &b_muon_chargedHadronIso);
  tree->SetBranchAddress("muon_neutralHadronIso", muon_neutralHadronIso, &b_muon_neutralHadronIso);
  tree->SetBranchAddress("muon_photonIso", muon_photonIso, &b_muon_photonIso);
  tree->SetBranchAddress("muon_puChargedHadronIso", muon_puChargedHadronIso, &b_muon_puChargedHadronIso);
  tree->SetBranchAddress("muon_relIso", muon_relIso, &b_muon_relIso);
  tree->SetBranchAddress("muon_pfIso", muon_pfIso, &b_muon_pfIso);
  tree->SetBranchAddress("nJets", &nJets, &b_nJets);
  tree->SetBranchAddress("jet_nConstituents", jet_nConstituents, &b_jet_nConstituents);
  tree->SetBranchAddress("jet_nChConstituents", jet_nChConstituents, &b_jet_nChConstituents);
  tree->SetBranchAddress("jet_charge", jet_charge, &b_jet_charge);
  tree->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
  tree->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
  tree->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
  tree->SetBranchAddress("jet_E", jet_E, &b_jet_E);
  tree->SetBranchAddress("jet_M", jet_M, &b_jet_M);
  tree->SetBranchAddress("jet_bdiscr", jet_bdiscr, &b_jet_bdiscr);
  if (newTrees && ! isData) tree->SetBranchAddress("jet_hadronFlavour", jet_hadronFlavour, &b_jet_hadronFlavour);   // TEMPORARILY !
  tree->SetBranchAddress("met_px", &met_px, &b_met_px);
  tree->SetBranchAddress("met_py", &met_py, &b_met_py);
  tree->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
  tree->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
  tree->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
  tree->SetBranchAddress("met_Et", &met_Et, &b_met_Et);
  tree->SetBranchAddress("met_E", &met_E, &b_met_E);
  tree->SetBranchAddress("met_corr_px", &met_corr_px, &b_met_corr_px);
  tree->SetBranchAddress("met_corr_py", &met_corr_py, &b_met_corr_py);
  tree->SetBranchAddress("met_corr_pt", &met_corr_pt, &b_met_corr_pt);
  tree->SetBranchAddress("met_corr_phi", &met_corr_phi, &b_met_corr_phi);
  tree->SetBranchAddress("met_corr_eta", &met_corr_eta, &b_met_corr_eta);
  tree->SetBranchAddress("met_corr_Et", &met_corr_Et, &b_met_corr_Et);
  tree->SetBranchAddress("met_corr_E", &met_corr_E, &b_met_corr_E);
  if (! isData)
  {
    if (newTrees)
    {
      tree->SetBranchAddress("nGenJets", &nGenJets, &b_nGenJets);
      tree->SetBranchAddress("genJet_pdgId", genJet_pdgId, &b_genJet_pdgId);
      tree->SetBranchAddress("genJet_charge", genJet_charge, &b_genJet_charge);
      tree->SetBranchAddress("genJet_pt", genJet_pt, &b_genJet_pt);
      tree->SetBranchAddress("genJet_phi", genJet_phi, &b_genJet_phi);
      tree->SetBranchAddress("genJet_eta", genJet_eta, &b_genJet_eta);
      tree->SetBranchAddress("genJet_E", genJet_E, &b_genJet_E);
      tree->SetBranchAddress("genJet_M", genJet_M, &b_genJet_M);
    }
    if (! isHerwig)
    {
      tree->SetBranchAddress("nMCParticles", &nMCParticles, &b_nMCParticles);
      tree->SetBranchAddress("mc_status", mc_status, &b_mc_status);
      tree->SetBranchAddress("mc_pdgId", mc_pdgId, &b_mc_pdgId);
      tree->SetBranchAddress("mc_mother", mc_mother, &b_mc_mother);
      tree->SetBranchAddress("mc_granny", mc_granny, &b_mc_granny);
      tree->SetBranchAddress("mc_pt", mc_pt, &b_mc_pt);
      tree->SetBranchAddress("mc_phi", mc_phi, &b_mc_phi);
      tree->SetBranchAddress("mc_eta", mc_eta, &b_mc_eta);
      tree->SetBranchAddress("mc_E", mc_E, &b_mc_E);
      tree->SetBranchAddress("mc_M", mc_M, &b_mc_M);
      tree->SetBranchAddress("mc_isLastCopy", mc_isLastCopy, &b_mc_isLastCopy);
      tree->SetBranchAddress("mc_isPromptFinalState", mc_isPromptFinalState, &b_mc_isPromptFinalState);
      tree->SetBranchAddress("mc_isHardProcess", mc_isHardProcess, &b_mc_isHardProcess);
      tree->SetBranchAddress("mc_fromHardProcessFinalState", mc_fromHardProcessFinalState, &b_mc_fromHardProcessFinalState);
      tree->SetBranchAddress("hasGenTop", &hasGenTop, &b_hasGenTop);
      tree->SetBranchAddress("hasGenAntiTop", &hasGenAntiTop, &b_hasGenAntiTop);
    }
    if (isTTbar)
    {
      tree->SetBranchAddress("weight1001", &weight1001, &b_weight1001);
      tree->SetBranchAddress("weight1002", &weight1002, &b_weight1002);
      tree->SetBranchAddress("weight1003", &weight1003, &b_weight1003);
      tree->SetBranchAddress("weight1004", &weight1004, &b_weight1004);
      tree->SetBranchAddress("weight1005", &weight1005, &b_weight1005);
      tree->SetBranchAddress("weight1007", &weight1007, &b_weight1007);
      tree->SetBranchAddress("weight1009", &weight1009, &b_weight1009);
      tree->SetBranchAddress("upFragWeight", &upFragWeight, &b_upFragWeight);
      tree->SetBranchAddress("centralFragWeight", &centralFragWeight, &b_centralFragWeight);
      tree->SetBranchAddress("downFragWeight", &downFragWeight, &b_downFragWeight);
      tree->SetBranchAddress("petersonFragWeight", &petersonFragWeight, &b_petersonFragWeight);
      tree->SetBranchAddress("semilepbrUp", &semilepbrUp, &b_semilepbrUp);
      tree->SetBranchAddress("semilepbrDown", &semilepbrDown, &b_semilepbrDown);
      if (newTrees)
      {
        tree->SetBranchAddress("pdfWeights", pdfWeights, &b_pdfWeights);
        tree->SetBranchAddress("pdfAlphaSUp", &pdfAlphaSUp, &b_pdfAlphaSUp);
        tree->SetBranchAddress("pdfAlphaSDown", &pdfAlphaSDown, &b_pdfAlphaSDown);
      }
    }
    tree->SetBranchAddress("btagSF", &btagSF, &b_btagSF);
    tree->SetBranchAddress("btagSF_up", &btagSF_up, &b_btagSF_up);
    tree->SetBranchAddress("btagSF_down", &btagSF_down, &b_btagSF_down);
    tree->SetBranchAddress("puSF", &puSF, &b_puSF);
    tree->SetBranchAddress("puSF_up", &puSF_up, &b_puSF_up);
    tree->SetBranchAddress("puSF_down", &puSF_down, &b_puSF_down);
    tree->SetBranchAddress("muonIdSF_BCDEF", muonIdSF_BCDEF, &b_muonIdSF_BCDEF);
    tree->SetBranchAddress("muonIdSF_GH", muonIdSF_GH, &b_muonIdSF_GH);
    tree->SetBranchAddress("muonIdSF_up_BCDEF", muonIdSF_up_BCDEF, &b_muonIdSF_up_BCDEF);
    tree->SetBranchAddress("muonIdSF_up_GH", muonIdSF_up_GH, &b_muonIdSF_up_GH);
    tree->SetBranchAddress("muonIdSF_down_BCDEF", muonIdSF_down_BCDEF, &b_muonIdSF_down_BCDEF);
    tree->SetBranchAddress("muonIdSF_down_GH", muonIdSF_down_GH, &b_muonIdSF_down_GH);
    tree->SetBranchAddress("muonIsoSF_BCDEF", muonIsoSF_BCDEF, &b_muonIsoSF_BCDEF);
    tree->SetBranchAddress("muonIsoSF_GH", muonIsoSF_GH, &b_muonIsoSF_GH);
    tree->SetBranchAddress("muonIsoSF_up_BCDEF", muonIsoSF_up_BCDEF, &b_muonIsoSF_up_BCDEF);
    tree->SetBranchAddress("muonIsoSF_up_GH", muonIsoSF_up_GH, &b_muonIsoSF_up_GH);
    tree->SetBranchAddress("muonIsoSF_down_BCDEF", muonIsoSF_down_BCDEF, &b_muonIsoSF_down_BCDEF);
    tree->SetBranchAddress("muonIsoSF_down_GH", muonIsoSF_down_GH, &b_muonIsoSF_down_GH);
    tree->SetBranchAddress("muonTrigSF_BCDEF", muonTrigSF_BCDEF, &b_muonTrigSF_BCDEF);
    tree->SetBranchAddress("muonTrigSF_GH", muonTrigSF_GH, &b_muonTrigSF_GH);
    tree->SetBranchAddress("muonTrigSF_up_BCDEF", muonTrigSF_up_BCDEF, &b_muonTrigSF_up_BCDEF);
    tree->SetBranchAddress("muonTrigSF_up_GH", muonTrigSF_up_GH, &b_muonTrigSF_up_GH);
    tree->SetBranchAddress("muonTrigSF_down_BCDEF", muonTrigSF_down_BCDEF, &b_muonTrigSF_down_BCDEF);
    tree->SetBranchAddress("muonTrigSF_down_GH", muonTrigSF_down_GH, &b_muonTrigSF_down_GH);
    tree->SetBranchAddress("muonTrackSF_eta", muonTrackSF_eta, &b_muonTrackSF_eta);
    tree->SetBranchAddress("muonTrackSF_aeta", muonTrackSF_aeta, &b_muonTrackSF_aeta);
    tree->SetBranchAddress("muonTrackSF_nPV", muonTrackSF_nPV, &b_muonTrackSF_nPV);
  }
  else
  {
    tree->SetBranchAddress("isDataRunB", &isDataRunB, &b_isDataRunB);
    tree->SetBranchAddress("isDataRunC", &isDataRunC, &b_isDataRunC);
    tree->SetBranchAddress("isDataRunD", &isDataRunD, &b_isDataRunD);
    tree->SetBranchAddress("isDataRunE", &isDataRunE, &b_isDataRunE);
    tree->SetBranchAddress("isDataRunF", &isDataRunF, &b_isDataRunF);
    tree->SetBranchAddress("isDataRunG", &isDataRunG, &b_isDataRunG);
    tree->SetBranchAddress("isDataRunH", &isDataRunH, &b_isDataRunH);
  }
}

void InitMSPlots()
{
  /// Do not split ttbar dataset
  MSPlot["nPVs_beforePU_"] = new MultiSamplePlot(datasets, "nPVs_beforePU_", 46, -0.5, 45.5, "# PVs");
  MSPlot["nPVs_afterPU_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_", 46, -0.5, 45.5, "# PVs");
  MSPlot["nPVs_afterPU_up_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_up_", 46, -0.5, 45.5, "# PVs");
  MSPlot["nPVs_afterPU_down_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_down_", 46, -0.5, 45.5, "# PVs");
  MSPlot["nPVs_beforePU_aSel_"] = new MultiSamplePlot(datasets, "nPVs_beforePU_aSel_", 46, -0.5, 45.5, "# PVs");
  MSPlot["nPVs_afterPU_aSel_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_aSel_", 46, -0.5, 45.5, "# PVs");
  MSPlot["nPVs_afterPU_begin16_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_begin16_", 46, -0.5, 45.5, "# PVs");
  MSPlot["nPVs_afterPU_end16_"] = new MultiSamplePlot(datasets, "nPVs_afterPU_end16_", 46, -0.5, 45.5, "# PVs");
  MSPlot["rho_"] = new MultiSamplePlot(datasets, "#rho", 41, -0.5, 40.5, "#rho");
  MSPlot["nJets_"] = new MultiSamplePlot(datasets, "nJets_", 13, -0.5, 12.5, "# jets");
  MSPlot["leadingJet_pT_"] = new MultiSamplePlot(datasets, "leadingJet_pT_", 40, 0, 400, "p_{T}", "GeV");
  MSPlot["leadingJet_pT_bKF_noTopCut_"] = new MultiSamplePlot(datasets, "leadingJet_pT_bKF_noTopCut_", 40, 0, 400, "p_{T}", "GeV");
  MSPlot["jet_pT_allJets_"] = new MultiSamplePlot(datasets, "jet_pT_allJets_", 40, 0, 400, "p_{T}", "GeV");
  MSPlot["leadingJet_pT_aKF_"] = new MultiSamplePlot(datasets, "leadingJet_pT_aKF_", 40, 0, 400, "p_{T}", "GeV");
  MSPlot["jet_pT_allJets_aKF_"] = new MultiSamplePlot(datasets, "jet_pT_allJets_aKF_", 40, 0, 400, "p_{T}", "GeV");
  MSPlot["btag_SF_"] = new MultiSamplePlot(datasets, "btag_SF_", 80, 0., 2., "btag SF");
  MSPlot["W_mass_"] = new MultiSamplePlot(datasets, "W mass before kinFitter_", 60, 0, 300, "m_{W}", "GeV");
  MSPlot["top_mass_"] = new MultiSamplePlot(datasets, "Top mass before kinFitter_", 50, 0, 500, "m_{t}", "GeV");
  
  /// Control plots
  if (makeControlPlots)
  {
    MSPlotCP["muon_pT"] = new MultiSamplePlot(datasetsMSP, "muon_pT", 22, 0, 440, "p_{T}", "GeV");
    MSPlotCP["muon_eta"] = new MultiSamplePlot(datasetsMSP, "muon_eta", 30, -3, 3, "Eta");
    MSPlotCP["muon_phi"] = new MultiSamplePlot(datasetsMSP, "muon_phi", 32, -3.2, 3.2, "Phi");
    MSPlotCP["muon_relIso"] = new MultiSamplePlot(datasetsMSP, "muon_relIso", 20, 0, 0.2, "relIso");
    MSPlotCP["muon_d0"] = new MultiSamplePlot(datasetsMSP, "muon_d0", 60, 0, 0.003, "d_{0}");
    MSPlotCP["leadingJet_pT"] = new MultiSamplePlot(datasetsMSP, "leadingJet_pT", 60, 0, 600, "p_{T}", "GeV");
    MSPlotCP["jet2_pT"] = new MultiSamplePlot(datasetsMSP, "jet2_pT", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet3_pT"] = new MultiSamplePlot(datasetsMSP, "jet3_pT", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet4_pT"] = new MultiSamplePlot(datasetsMSP, "jet4_pT", 40, 0, 200, "p_{T}", "GeV");
    MSPlotCP["jet_pT_allJets"] = new MultiSamplePlot(datasetsMSP, "jet_pT_allJets", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["Ht_4leadingJets"] = new MultiSamplePlot(datasetsMSP,"Ht_4leadingJets", 60, 0, 1200, "H_{T}", "GeV");
    MSPlotCP["met_pT"] = new MultiSamplePlot(datasetsMSP, "met_pT", 40, 0, 400, "E_{T}^{miss} p_{T}", "GeV");
    MSPlotCP["met_eta"] = new MultiSamplePlot(datasetsMSP, "met_eta", 30, -3, 3, "E_{T}^{miss} Eta");
    MSPlotCP["met_phi"] = new MultiSamplePlot(datasetsMSP, "met_phi", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
    MSPlotCP["met_corr_pT"] = new MultiSamplePlot(datasetsMSP, "met_corr_pT", 40, 0, 400, "E_{T}^{miss} p_{T}", "GeV");
    MSPlotCP["met_corr_eta"] = new MultiSamplePlot(datasetsMSP, "met_corr_eta", 30, -3, 3, "E_{T}^{miss} Eta");
    MSPlotCP["met_corr_phi"] = new MultiSamplePlot(datasetsMSP, "met_corr_phi", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
    
    MSPlotCP["muon_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_pT_aKF", 22, 0, 440, "p_{T}", "GeV");
    MSPlotCP["muon_eta_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_eta_aKF", 30, -3, 3, "Eta");
    MSPlotCP["muon_phi_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_phi_aKF", 32, -3.2, 3.2, "Phi");
    MSPlotCP["muon_relIso_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_relIso_aKF", 20, 0, 0.2, "relIso");
    MSPlotCP["muon_d0_aKF"] = new MultiSamplePlot(datasetsMSP, "muon_d0_aKF", 60, 0, 0.003, "d_{0}");
    MSPlotCP["leadingJet_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "leadingJet_pT_aKF", 60, 0, 600, "p_{T}", "GeV");
    MSPlotCP["jet2_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "jet2_pT_aKF", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet3_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "jet3_pT_aKF", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["jet4_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "jet4_pT_aKF", 40, 0, 200, "p_{T}", "GeV");
    MSPlotCP["jet_pT_allJets_aKF"] = new MultiSamplePlot(datasetsMSP, "jet_pT_allJets_aKF", 40, 0, 400, "p_{T}", "GeV");
    MSPlotCP["Ht_4leadingJets_aKF"] = new MultiSamplePlot(datasetsMSP,"Ht_4leadingJets_aKF", 60, 0, 1200, "H_{T}", "GeV");
    MSPlotCP["met_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "met_pT_aKF", 40, 0, 400, "E_{T}^{miss} p_{T}", "GeV");
    MSPlotCP["met_eta_aKF"] = new MultiSamplePlot(datasetsMSP, "met_eta_aKF", 30, -3, 3, "E_{T}^{miss} Eta");
    MSPlotCP["met_phi_aKF"] = new MultiSamplePlot(datasetsMSP, "met_phi_aKF", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
    MSPlotCP["met_corr_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "met_corr_pT_aKF", 40, 0, 400, "E_{T}^{miss} p_{T}", "GeV");
    MSPlotCP["met_corr_eta_aKF"] = new MultiSamplePlot(datasetsMSP, "met_corr_eta_aKF", 30, -3, 3, "E_{T}^{miss} Eta");
    MSPlotCP["met_corr_phi_aKF"] = new MultiSamplePlot(datasetsMSP, "met_corr_phi_aKF", 32, -3.2, 3.2, "E_{T}^{miss} Phi");
    
    
    MSPlotCP["nJets"] = new MultiSamplePlot(datasetsMSP, "nJets", 13, -0.5, 12.5, "# jets");
    MSPlotCP["nBJets"] = new MultiSamplePlot(datasetsMSP, "nBJets", 9, -0.5, 8.5, "# b jets");
    MSPlotCP["CSVv2Discr_allJets"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_allJets", 48, 0.0, 1.2, "CSVv2 discriminant value");
    MSPlotCP["CSVv2Discr_leadingJet"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_leadingJet", 48, 0.0, 1.2, "CSVv2 discriminant value of leading jet");
    MSPlotCP["CSVv2Discr_jet2"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet2", 48, 0.0, 1.2, "CSVv2 discriminant value of jet2");
    MSPlotCP["CSVv2Discr_jet3"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet3", 48, 0.0, 1.2, "CSVv2 discriminant value of jet3");
    MSPlotCP["CSVv2Discr_jet4"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet4", 48, 0.0, 1.2, "CSVv2 discriminant value of jet4");
    MSPlotCP["CSVv2Discr_highest"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_highest", 48, 0.0, 1.2, "Highest CSVv2 discriminant value");
    MSPlotCP["CSVv2Discr_jetNb"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jetNb", 8, -0.5, 7.5, "Jet number (in order of decreasing p_{T}) with highest CSVv2 discriminant value");
    
    MSPlotCP["nJets_aKF"] = new MultiSamplePlot(datasetsMSP, "nJets_aKF", 13, -0.5, 12.5, "# jets");
    MSPlotCP["nBJets_aKF"] = new MultiSamplePlot(datasetsMSP, "nBJets_aKF", 9, -0.5, 8.5, "# b jets");
    MSPlotCP["CSVv2Discr_allJets_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_allJets_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value");
    MSPlotCP["CSVv2Discr_leadingJet_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_leadingJet_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value of leading jet");
    MSPlotCP["CSVv2Discr_jet2_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet2_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value of jet2");
    MSPlotCP["CSVv2Discr_jet3_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet3_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value of jet3");
    MSPlotCP["CSVv2Discr_jet4_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jet4_aKF", 48, 0.0, 1.2, "CSVv2 discriminant value of jet4");
    MSPlotCP["CSVv2Discr_highest_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_highest_aKF", 48, 0.0, 1.2, "Highest CSVv2 discriminant value");
    MSPlotCP["CSVv2Discr_jetNb_aKF"] = new MultiSamplePlot(datasetsMSP, "CSVv2Discr_jetNb_aKF", 8, -0.5, 7.5, "Jet number (in order of decreasing p_{T}) with highest CSVv2 discriminant value");
    
    MSPlotCP["M3"] = new MultiSamplePlot(datasetsMSP, "M3", 40, 60, 460, "M_{3}", "GeV");
    MSPlotCP["min_Mlb"] = new MultiSamplePlot(datasetsMSP, "min_Mlb", 40, 0, 400, "M_{lb}", "GeV");
    MSPlotCP["dR_Lep_B"] = new MultiSamplePlot(datasetsMSP, "dR_Lep_B", 25, 0, 5, "#Delta R(l,b)");
    
    MSPlotCP["M3_aKF"] = new MultiSamplePlot(datasetsMSP, "M3_aKF", 40, 60, 460, "M_{3}", "GeV");
    MSPlotCP["min_Mlb_aKF"] = new MultiSamplePlot(datasetsMSP, "min_Mlb_aKF", 40, 0, 400, "M_{lb}", "GeV");
    MSPlotCP["dR_Lep_B_aKF"] = new MultiSamplePlot(datasetsMSP, "dR_Lep_B_aKF", 25, 0, 5, "#Delta R(l,b)");
  }
  
  /// SFs
  MSPlot["scaleFactor"] = new MultiSamplePlot(datasetsMSP, "scaleFactor", 80, 0., 2., "SF");
  MSPlot["btag_SF"] = new MultiSamplePlot(datasetsMSP, "btag_SF", 80, 0., 2., "btag SF");
  MSPlot["pu_SF"] = new MultiSamplePlot(datasetsMSP, "pu_SF", 80, 0., 2., "pu SF");
  
  /// Reco
  MSPlot["W_mass"] = new MultiSamplePlot(datasetsMSP, "W mass before kinFitter", 60, 0, 300, "m_{W}", "GeV");
  MSPlot["W_mass_T"] = new MultiSamplePlot(datasetsMSP, "Transverse W mass before kinFitter", 70, 0, 350, "m_{T,W}", "GeV");
  MSPlot["top_mass"] = new MultiSamplePlot(datasetsMSP, "Top mass before kinFitter", 50, 0, 500, "m_{t}", "GeV");
  MSPlot["top_mass_400"] = new MultiSamplePlot(datasetsMSP, "Top mass before kinFitter (400)", 40, 0, 400, "m_{t}", "GeV");
  MSPlot["top_mass_zoom"] = new MultiSamplePlot(datasetsMSP, "Top mass before kinFitter (zoomed)", 40, 110, 230, "m_{t}", "GeV");
  MSPlot["top_mass_alt"] = new MultiSamplePlot(datasetsMSP, "Alternative top mass before kinFitter", 50, 100, 600, "m_{t,alt}", "GeV");
  MSPlot["top_mass_rescaled"] = new MultiSamplePlot(datasetsMSP, "Top mass rescaled by 80.385/m_{W} before kinFitter", 50, 0, 500, "m_{t}", "GeV");
  MSPlot["red_top_mass_manyBins"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass (many bins)", 640, 0.4, 2., "m_{r}");
  MSPlot["red_top_mass"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass", 32, 0.4, 2., "m_{r}");
  MSPlot["red_top_mass_old"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass (old definition)", 32, 0.4, 2., "m_{r, old}");
  MSPlot["red_top_mass_new"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass (new definition)", 32, 0.4, 2., "m_{r, new}");
  MSPlot["W_pT"] = new MultiSamplePlot(datasetsMSP, "W p_{T} before kinFitter", 60, 0, 300, "p_{T}", "GeV");
  MSPlot["top_pT"] = new MultiSamplePlot(datasetsMSP, "Top p_{T} before kinFitter", 80, 0, 400, "p_{T}", "GeV");
  
  MSPlot["light_jets_pT"] = new MultiSamplePlot(datasetsMSP, "p_{T} of the light jets before kinFitter", 80, 0, 400, "p_{T}", "GeV");
  MSPlot["hadr_b_jet_pT"] = new MultiSamplePlot(datasetsMSP, "p_{T} of the hadronic b jet before KF cuts", 80, 0, 400, "p_{T}", "GeV");
  MSPlot["lept_b_jet_pT"] = new MultiSamplePlot(datasetsMSP, "p_{T} of the leptonic b jet before KF cuts", 80, 0, 400, "p_{T}", "GeV");
  
  MSPlot["top_mass_cand_diff"] = new MultiSamplePlot(datasetsMSP, "Mass difference between top quark candidates before kinFitter", 60, -150, 150, "m_{t} - m_{t'}", "GeV");
    MSPlot["top_pt_cand_diff"] = new MultiSamplePlot(datasetsMSP, "P_{T} difference between top quark candidates before kinFitter", 60, -150, 150, "p_{T,t} - p_{T,t'}", "GeV");
  
  MSPlot["mlb"] = new MultiSamplePlot(datasetsMSP, "mlb before kinFitter", 80, 0, 800, "m_{lb}", "GeV");
  MSPlot["red_mlb_mass"] = new MultiSamplePlot(datasetsMSP, "Reduced m_{lb} mass before kinFitter", 100, 0., 2.5, "m_{lb,r}");
  
  MSPlot["ttbar_mass"] = new MultiSamplePlot(datasetsMSP, "ttbar mass before kinFitter", 50, 0, 1000, "m_{t#bar{t}}", "GeV");
  MSPlot["mass_bjj_div_m_jj"] = new MultiSamplePlot(datasetsMSP, "Top quark mass divided by W boson mass before kinfit cut", 50, 1., 3.5, "m_{bjj} /m_{jj}");
  MSPlot["mass_bjj_div_m_jj_2"] = new MultiSamplePlot(datasetsMSP, "Top quark mass divided by W boson mass before kinfit cut 2", 40, 0.5, 2., "m_{bjj} /m_{jj}");
  MSPlot["mass_bjj_div_m_lb"] = new MultiSamplePlot(datasetsMSP, "Reconstructed top quark mass divided by mass of lepton and leptonic b jet before kinfit cut", 100, 0.5, 5.5, "m_{bjj} /m_{lb}");
  
  MSPlot["dR_lep_b_min"] = new MultiSamplePlot(datasetsMSP, "Minimum dR(lep,b) before kinFitter", 25, 0, 5, "#Delta R(l,b_{l})");
  MSPlot["dR_lep_b_max"] = new MultiSamplePlot(datasetsMSP, "Maximum dR(lep,b) before kinFitter", 25, 0, 5, "#Delta R(l,b_{h})");
  MSPlot["dR_light_jets"] = new MultiSamplePlot(datasetsMSP, "dR between light jets before kinFitter", 25, 0, 5, "#Delta R(j_{1},j_{2})");
  MSPlot["dR_b_light_min"] = new MultiSamplePlot(datasetsMSP, "Minimum dR between hadronic b and one of the light jets before kinFitter", 25, 0, 5, "#Delta R(b,j_{1})");
  MSPlot["dR_b_light_max"] = new MultiSamplePlot(datasetsMSP, "Maximum dR between hadronic b and one of the light jets before kinFitter", 25, 0, 5, "#Delta R(b,j_{2})");
  MSPlot["dR_b_W"] = new MultiSamplePlot(datasetsMSP, "Minimum dR between hadronic b jet and the reconstructed W boson before kinFitter", 25, 0, 5, "#Delta R(b,W)");
  MSPlot["dR_b_light_sum"] = new MultiSamplePlot(datasetsMSP, "Quadratic sum of dR between hadronic b and light jets before kinFitter", 50, 0, 10, "#sqrt{(#Delta R(b_{h},j_{1}))^{2} + (#Delta R(b_{h},j_{2}))^{2}}");
  MSPlot["dR_bb"] = new MultiSamplePlot(datasetsMSP, "dR between hadronic and leptonic b before kinFitter", 25, 0, 5, "#Delta R(b_{h},b_{l})");
  MSPlot["dPhi_light"] = new MultiSamplePlot(datasetsMSP, "d#phi between light jets before kinFitter", 32, -3.2, 3.2, "#Delta#phi(j_{1},j_{2})");
  MSPlot["dPhi_bb"] = new MultiSamplePlot(datasetsMSP, "d#phi between hadronic and leptonic b before kinFitter", 32, -3.2, 3.2, "#Delta#phi(b_{h},b_{l})");
  MSPlot["dPhi_bW"] = new MultiSamplePlot(datasetsMSP, "d#phi between hadronic b jet and the reconstructed W boson before kinFitter", 32, -3.2, 3.2, "#Delta#phi(b,W)");
  
  if (doKinFit)
  {
    MSPlot["scaleFactor_aKF"] = new MultiSamplePlot(datasetsMSP, "scaleFactor_aKF", 80, 0., 2., "SF");
    MSPlot["btag_SF_aKF"] = new MultiSamplePlot(datasetsMSP, "btag_SF_aKF", 80, 0., 2., "btag SF");
    MSPlot["pu_SF_aKF"] = new MultiSamplePlot(datasetsMSP, "pu_SF_aKF", 80, 0., 2., "pu SF");
    MSPlot["W_mass_aKF"] = new MultiSamplePlot(datasetsMSP, "W mass after kinFitter", 40, 0, 200, "m_{W}", "GeV");
    MSPlot["W_mass_T_aKF"] = new MultiSamplePlot(datasetsMSP, "Transverse W mass after kinFitter", 70, 0, 350, "m_{T,W}", "GeV");
    MSPlot["top_mass_aKF"] = new MultiSamplePlot(datasetsMSP, "Top mass after kinFitter", 40, 0, 400, "m_{t}", "GeV");
    MSPlot["top_mass_aKF_zoom"] = new MultiSamplePlot(datasetsMSP, "Top mass after kinFitter (zoomed)", 40, 110, 230, "m_{t}", "GeV");
    MSPlot["top_mass_alt_aKF"] = new MultiSamplePlot(datasetsMSP, "Alternative top mass after kinFitter", 50, 100, 600, "m_{t,alt}", "GeV");
    MSPlot["red_top_mass_aKF_manyBins"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass after KF (many bins)", 640, 0.4, 2., "m_{r}");
    MSPlot["red_top_mass_aKF"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass after KF", 32, 0.4, 2., "m_{r}");
    MSPlot["red_top_mass_old_aKF"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass after KF (old definition)", 32, 0.4, 2., "m_{r, old}");
    MSPlot["red_top_mass_new_aKF"] = new MultiSamplePlot(datasetsMSP, "Reduced top quark mass after KF (new definition)", 32, 0.4, 2., "m_{r, new}");
    MSPlot["W_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "W p_{T} after kinFitter", 60, 0, 300, "p_{T}", "GeV");
    MSPlot["top_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "Top p_{T} after kinFitter", 80, 0, 400, "p_{T}", "GeV");
    
    MSPlot["light_jets_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "p_{T} of the light jets after kinFitter", 80, 0, 400, "p_{T}", "GeV");
    MSPlot["hadr_b_jet_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "p_{T} of the hadronic b jet after KF cuts", 80, 0, 400, "p_{T}", "GeV");
    MSPlot["lept_b_jet_pT_aKF"] = new MultiSamplePlot(datasetsMSP, "p_{T} of the leptonic b jet after KF cuts", 80, 0, 400, "p_{T}", "GeV");
    
    MSPlot["mlb_aKF"] = new MultiSamplePlot(datasetsMSP, "mlb after kinFitter", 80, 0, 800, "m_{lb}", "GeV");
    MSPlot["red_mlb_mass_aKF"] = new MultiSamplePlot(datasetsMSP, "Reduced m_{lb} mass after kinFitter", 100, 0., 2.5, "m_{lb,r}");
    
    MSPlot["ttbar_mass_aKF"] = new MultiSamplePlot(datasetsMSP, "ttbar mass after kinFitter", 50, 0, 1000, "m_{t#bar{t}}", "GeV");
    MSPlot["mass_bjj_div_m_jj_aKF"] = new MultiSamplePlot(datasetsMSP, "Top quark mass divided by W boson mass", 50, 1., 3.5, "m_{bjj} /m_{jj}");
    MSPlot["mass_bjj_div_m_jj_2_aKF"] = new MultiSamplePlot(datasetsMSP, "Top quark mass divided by W boson mass 2", 40, 0.5, 2., "m_{bjj} /m_{jj}");
    MSPlot["mass_bjj_div_m_lb_aKF"] = new MultiSamplePlot(datasetsMSP, "Reconstructed top quark mass divided by mass of lepton and leptonic b jet after kinfit cut", 100, 0.5, 5.5, "m_{bjj} /m_{lb}");
    
    MSPlot["dR_lep_b_min_aKF"] = new MultiSamplePlot(datasetsMSP, "Minimum dR(lep,b) after kinFitter", 25, 0, 5, "#Delta R(l,b_{l})");
    MSPlot["dR_lep_b_max_aKF"] = new MultiSamplePlot(datasetsMSP, "Maximum dR(lep,b) after kinFitter", 25, 0, 5, "#Delta R(l,b_{h})");
    MSPlot["dR_light_jets_aKF"] = new MultiSamplePlot(datasetsMSP, "dR between light jets after kinFitter", 25, 0, 5, "#Delta R(j_{1},j_{2})");
    MSPlot["dR_b_light_min_aKF"] = new MultiSamplePlot(datasetsMSP, "Minimum dR between hadronic b and one of the light jets after kinFitter", 25, 0, 5, "#Delta R(b,j_{1})");
    MSPlot["dR_b_light_max_aKF"] = new MultiSamplePlot(datasetsMSP, "Maximum dR between hadronic b and one of the light jets after kinFitter", 25, 0, 5, "#Delta R(b,j_{2})");
    MSPlot["dR_b_W_aKF"] = new MultiSamplePlot(datasetsMSP, "Minimum dR between hadronic b jet and the reconstructed W boson after kinFitter", 25, 0, 5, "#Delta R(b,W)");
    MSPlot["dR_b_light_sum_aKF"] = new MultiSamplePlot(datasetsMSP, "Quadratic sum of dR between hadronic b and light jets after kinFitter", 50, 0, 10, "#sqrt{(#Delta R(b_{h},j_{1}))^{2} + (#Delta R(b_{h},j_{2}))^{2}}");
    MSPlot["dR_bb_aKF"] = new MultiSamplePlot(datasetsMSP, "dR between hadronic and leptonic b after kinFitter", 25, 0, 5, "#Delta R(b_{h},b_{l})");
    MSPlot["dPhi_light_aKF"] = new MultiSamplePlot(datasetsMSP, "d#phi between light jets after kinFitter", 32, -3.2, 3.2, "#Delta#phi(j_{1},j_{2})");
    MSPlot["dPhi_bb_aKF"] = new MultiSamplePlot(datasetsMSP, "d#phi between hadronic and leptonic b after kinFitter", 32, -3.2, 3.2, "#Delta#phi(b_{h},b_{l})");
    MSPlot["dPhi_bW_aKF"] = new MultiSamplePlot(datasetsMSP, "d#phi between hadronic b jet and the reconstructed W boson after kinFitter", 32, -3.2, 3.2, "#Delta#phi(b,W)");
    
    MSPlot["KF_Chi2"] = new MultiSamplePlot(datasetsMSP, "Chi2 value of kinFitter", 50, 0, 5, "#chi^{2}");
    MSPlot["KF_Chi2_narrow"] = new MultiSamplePlot(datasetsMSP, "Chi2 value of kinFitter (zoomed)", 50, 0, 2, "#chi^{2}");
    MSPlot["KF_Chi2_wide"] = new MultiSamplePlot(datasetsMSP, "Chi2 value of kinFitter (large)", 50, 0, 20, "#chi^{2}");
    
    MSPlot["W_mass_diff"] = new MultiSamplePlot(datasetsMSP, "Difference between W mass before and after kinFitter", 100, -50, 50, "Delta m_{W}", "GeV");
    MSPlot["W_mass_T_diff"] = new MultiSamplePlot(datasetsMSP, "Difference between transverse W mass before and after kinFitter", 150, -75, 75, "Delta m_{W}", "GeV");
    MSPlot["W_pT_diff"] = new MultiSamplePlot(datasetsMSP, "Difference between W p_{T} before and after kinFitter", 100, -50, 50, "#Delta p_{T}", "GeV");
    MSPlot["top_mass_diff"] = new MultiSamplePlot(datasetsMSP, "Difference between top mass before and after kinFitter", 150, -75, 75, "Delta m_{t}", "GeV");
    MSPlot["top_pT_diff"] = new MultiSamplePlot(datasetsMSP, "Difference between top p_{T} before and after kinFitter", 100, -50, 50, "#Delta p_{T}", "GeV");
    //MSPlot["top_px_diff"] = new MultiSamplePlot(datasetsMSP, "Difference between top p_{x} before and after kinFitter", 100, -50, 50, "#Delta p_{x}", "GeV");
    //MSPlot["top_py_diff"] = new MultiSamplePlot(datasetsMSP, "Difference between top p_{y} before and after kinFitter", 100, -50, 50, "#Delta p_{y}", "GeV");
    //MSPlot["top_pz_diff"] = new MultiSamplePlot(datasetsMSP, "Difference between top p_{z} before and after kinFitter", 100, -50, 50, "#Delta p_{z}", "GeV");
    
    MSPlot["top_mass_cand_diff_aKF"] = new MultiSamplePlot(datasetsMSP, "Mass difference between top quark candidates after kinFitter", 60, -150, 150, "m_{t} - m_{t'}", "GeV");
    MSPlot["top_pt_cand_diff_aKF"] = new MultiSamplePlot(datasetsMSP, "P_{T} difference between top quark candidates after kinFitter", 60, -150, 150, "p_{T,t} - p_{T,t'}", "GeV");
  }
}

void InitHisto1D()
{
  TH1::SetDefaultSumw2();
  
  InitHisto1DMatch();
  
  histo1D["genTop_status"] = new TH1F("genTop_status", "Status code of the generated top quark; status", 100, -0.5, 99.5);
  histo1D["genTop_isLastCopy_status"] = new TH1F("genTop_isLastCopy_status", "Status code of the last copy of the generated top quark; status", 100, -0.5, 99.5);
  histo1D["genTop_hasLastCopy"] = new TH1F("genTop_hasLastCopy", "Generated top quark has isLastCopy tag; has isLastCopy", 2, -0.5, 1.5);
  histo1D["genTop_hasStatus22"] = new TH1F("genTop_hasStatus22", "Generated top quark has status 22; has status 22", 2, -0.5, 1.5);
  histo1D["genAntitop_status"] = new TH1F("genAntitop_status", "Status code of the generated antitop quark; status", 100, -0.5, 99.5);
  histo1D["genAntitop_isLastCopy_status"] = new TH1F("genAntitop_isLastCopy_status", "Status code of the last copy of the generated antitop quark; status", 100, -0.5, 99.5);
  histo1D["genAntitop_hasLastCopy"] = new TH1F("genAntitop_hasLastCopy", "Generated antitop quark has isLastCopy tag; has isLastCopy", 2, -0.5, 1.5);
  histo1D["genAntitop_hasStatus22"] = new TH1F("genAntitop_hasStatus22", "Generated antitop quark has status 22; has status 22", 2, -0.5, 1.5);
  
  histo1D["genTop_status22_pT"] = new TH1F("genTop_status22_pT", "pT of the generated top quark with status 22; p_{T} (GeV)", 250, 0., 250.);
  histo1D["genTop_status62_pT"] = new TH1F("genTop_status62_pT", "pT of the generated top quark with status 62; p_{T} (GeV)", 250, 0., 250.);
  histo1D["genAntitop_status22_pT"] = new TH1F("genAntitop_status22_pT", "pT of the generated antitop quark with status 22; p_{T} (GeV)", 250, 0., 250.);
  histo1D["genAntitop_status62_pT"] = new TH1F("genAntitop_status62_pT", "pT of the generated antitop quark with status 62; p_{T} (GeV)", 250, 0., 250.);
  
  histo1D["leadingJet_pT_reco_no22or62"] = new TH1F("leadingJet_pT_reco_no22or62", "pT of the leading reconstructed jet for events that have no top quark with status 22 or 62; p_{T} (GeV)", 250, 0., 500.);
  
  histo1D["nPartons_dilep"] = new TH1F("nPartons_dilep", "nPartons_dilep", 10, -0.5, 9.5);
  histo1D["nPartons_semilep"] = new TH1F("nPartons_semilep", "nPartons_semilep", 10, -0.5, 9.5);
  histo1D["nPartons_allhad"] = new TH1F("nPartons_allhad", "nPartons_allhad", 10, -0.5, 9.5);
  
  /// SFs
  histo1D["width_SF"] = new TH1F("width_SF", "Scale factor to change the ttbar distribution width; width SF", 500, 0, 5);
  //histo1D["btag_SF"] = new TH1F("btag_SF", "b tag scale factor; btag SF", 80, 0, 2);
  
  if (makeControlPlots)
  {
    histo1D["dR_jets_min"]  = new TH1F("dR_jets_min","Minimal delta R between two jets; #Delta R(j_{1},j_{2})", 35, 0, 3.5);
//    histo1D["dR_jets_sum"]  = new TH1F("dR_jets_sum","Sum of delta R between all jets; #Delta R(j_{1},j_{2})", 60, 0, 12);
  }
  
  
  /// Systematic comparison
  histo1D["allSim_top_mass"] = new TH1F("allSim_top_mass","Reconstructed top mass for all simulated samples; m_{t}", 32, 130, 210);
  histo1D["allSim_red_top_mass"] = new TH1F("allSim_red_top_mass","Reduced top mass for all simulated samples; m_{r}", 64, 0.4, 2.);
  histo1D["allSim_red_top_mass_old"] = new TH1F("allSim_red_top_mass_old","Reduced top mass for all simulated samples (old definition); m_{r, old}", 64, 0.4, 2.);
  histo1D["allSim_red_top_mass_new"] = new TH1F("allSim_red_top_mass_new","Reduced top mass for all simulated samples (new definition); m_{r, new}", 64, 0.4, 2.);
  histo1D["allSim_mass_bjj_div_m_jj"] = new TH1F("allSim_mass_bjj_div_m_jj","Reconstructed top quark mass divided by reconstructed W boson mass for all simulated samples; m_{bjj} /m_{jj}",  50, 1., 3.5);
  
  /// m_3/2
  histo1D["mass_bjj_div_m_jj_CM"] = new TH1F("mass_bjj_div_m_jj_CM","Reconstructed top quark mass divided by reconstructed W boson mass (CM); m_{bjj} /m_{jj}",  50, 1., 3.5);
  histo1D["mass_bjj_div_m_jj_WM"] = new TH1F("mass_bjj_div_m_jj_WM","Reconstructed top quark mass divided by reconstructed W boson mass (WM); m_{bjj} /m_{jj}",  50, 1., 3.5);
  histo1D["mass_bjj_div_m_jj_UM"] = new TH1F("mass_bjj_div_m_jj_UM","Reconstructed top quark mass divided by reconstructed W boson mass (UM); m_{bjj} /m_{jj}",  50, 1., 3.5);
  
  histo1D["mass_bjj_div_m_jj_2_CM"] = new TH1F("mass_bjj_div_m_jj_2_CM","Reconstructed top quark mass divided by reconstructed W boson mass (CM); m_{bjj} /m_{jj}",  50, 0.5, 2.);
  histo1D["mass_bjj_div_m_jj_2_WM"] = new TH1F("mass_bjj_div_m_jj_2_WM","Reconstructed top quark mass divided by reconstructed W boson mass (WM); m_{bjj} /m_{jj}",  50, 0.5, 2.);
  histo1D["mass_bjj_div_m_jj_2_UM"] = new TH1F("mass_bjj_div_m_jj_2_UM","Reconstructed top quark mass divided by reconstructed W boson mass (UM); m_{bjj} /m_{jj}",  50, 0.5, 2.);
  
  /// m_t/<m_t>
  histo1D["red_top_mass_TT_CM"] = new TH1F("red_top_mass_TT_CM","Reduced top mass for matched TT sample (reco, correct top match); m_{r}", 64, 0.4, 2.);
  histo1D["red_top_mass_TT_WMUM"] = new TH1F("red_top_mass_TT_WMUM","Reduced top mass for unmatched TT sample (reco, no & wrong top match); m_{r}", 64, 0.4, 2.);
  histo1D["red_top_mass_TT_UM"] = new TH1F("red_top_mass_TT_UM","Reduced top mass for unmatched TT sample (reco, no top match); m_{r}", 64, 0.4, 2.);
  histo1D["red_top_mass_TT_WM"] = new TH1F("red_top_mass_TT_WM","Reduced top mass for matched TT sample (reco, wrong top match: wrong permutation); m_{r}", 64, 0.4, 2.);
  
  histo1D["red_top_mass_bkgd"] = new TH1F("red_top_mass_bkgd","Reduced top mass for background samples; m_{r}", 64, 0.4, 2.);
  histo1D["red_top_mass_bkgd_manyBins"] = new TH1F("red_top_mass_bkgd_manyBins","Reduced top mass for background samples; m_{r}", 640, 0.4, 2.);
  histo1D["red_top_mass_TT"] = new TH1F("red_top_mass_TT","Reduced top mass for TT sample; m_{r}", 640, 0.4, 2.);
  histo1D["red_top_mass_ST_tW_top"] = new TH1F("red_top_mass_ST_tW_top","Reduced top mass for ST tW top sample; m_{r}", 640, 0.4, 2.);
  histo1D["red_top_mass_ST_tW_antitop"] = new TH1F("red_top_mass_ST_tW_antitop","Reduced top mass for ST tW antitop sample; m_{r}", 640, 0.4, 2.);
  histo1D["red_top_mass_ST_t_top"] = new TH1F("red_top_mass_ST_t_top","Reduced top mass for ST t top sample; m_{r}", 640, 0.4, 2.);
  histo1D["red_top_mass_ST_t_antitop"] = new TH1F("red_top_mass_ST_t_antitop","Reduced top mass for ST t antitop sample; m_{r}", 640, 0.4, 2.);
  histo1D["red_top_mass_DYJets"] = new TH1F("red_top_mass_DYJets","Reduced top mass for DY+Jets sample; m_{r}", 640, 0.4, 2.);
  histo1D["red_top_mass_WJets"] = new TH1F("red_top_mass_WJets","Reduced top mass for W+Jets sample; m_{r}", 640, 0.4, 2.);
  histo1D["red_top_mass_data"] = new TH1F("red_top_mass_data","Reduced top mass for data sample; m_{r}", 640, 0.4, 2.);
  
  /// mlb
  histo1D["mlb_CM"]  = new TH1F("minMlb_CM","m_{lb} mass for events that have a correct hadronic top match (CM); m_{lb} [GeV]", 400, 0, 800);
  histo1D["mlb_WM"]  = new TH1F("minMlb_WM","m_{lb} for events that have a wrong hadronic top match (WM); m_{lb} [GeV]", 400, 0, 800);
  histo1D["mlb_UM"]  = new TH1F("minMlb_UM","m_{lb} for events that have no hadronic top match (UM); m_{lb} [GeV]", 400, 0, 800);
  histo1D["red_mlb_mass_TT_CM"] = new TH1F("red_mlb_mass_TT_CM","Reduced m_{lb} mass for matched TT sample (reco, correct top match); m_{lb,r}", 100, 0., 2.5);
  histo1D["red_mlb_mass_TT_WMUM"] = new TH1F("red_mlb_mass_TT_WMUM","Reduced m_{lb} mass for unmatched TT sample (reco, no correct top match); m_{lb,r}", 100, 0., 2.5);
  histo1D["red_mlb_mass_TT_UM"] = new TH1F("red_mlb_mass_TT_UM","Reduced m_{lb} mass for unmatched TT sample (reco, no top match); m_{lb,r}", 100, 0., 2.5);
  histo1D["red_mlb_mass_TT_WM"] = new TH1F("red_mlb_mass_TT_WM","Reduced m_{lb} mass for matched TT sample (reco, wrong top match); m_{lb,r}", 100, 0., 2.5);
  
  /// dR
  //  lepton, b(lep)
  histo1D["dR_lep_b_lep_CM"]  = new TH1F("dR_lep_b_lep_CM","Minimal delta R between the lepton and the leptonic b jet (reco, correct match); #Delta R(l,b_{l})", 25, 0, 5);
  histo1D["dR_lep_b_lep_WM"]  = new TH1F("dR_lep_b_lep_WM","Minimal delta R between the lepton and the leptonic b jet (reco, wrong permutation); #Delta R(l,b_{l})", 25, 0, 5);
  histo1D["dR_lep_b_lep_UM"]  = new TH1F("dR_lep_b_lep_UM","Minimal delta R between the lepton and the leptonic b jet (reco, no match); #Delta R(l,b_{l})", 25, 0, 5);
  //  lepton, b(hadr)
  histo1D["dR_lep_b_had_CM"]  = new TH1F("dR_lep_b_had_CM","Minimal delta R between the lepton and the hadronic b jet (reco, correct match); #Delta R(l,b_{h})", 25, 0, 5);
  histo1D["dR_lep_b_had_WM"]  = new TH1F("dR_lep_b_had_WM","Minimal delta R between the lepton and the hadronic b jet (reco, wrong permutation); #Delta R(l,b_{h})", 25, 0, 5);
  histo1D["dR_lep_b_had_UM"]  = new TH1F("dR_lep_b_had_UM","Minimal delta R between the lepton and the hadronic b jet (reco, no match); #Delta R(l,b_{h})", 25, 0, 5);
  
  /// ttbar mass
  histo1D["ttbar_mass_CM"] = new TH1F("ttbar_mass_CM","Reconstructed mass of the top quark pair (reco, correct match); m_{t#bar{t}} [GeV]", 500, 0, 1000);
  histo1D["ttbar_mass_WM"] = new TH1F("ttbar_mass_WM","Reconstructed mass of the top quark pair (reco, wrong permutation); m_{t#bar{t}} [GeV]", 500, 0, 1000);
  histo1D["ttbar_mass_UM"] = new TH1F("ttbar_mass_UM","Reconstructed mass of the top quark pair (reco, no match); m_{t#bar{t}} [GeV]", 500, 0, 1000);
  
  
  /// KinFitter
  if (doKinFit)
  {
    histo1D["KF_Chi2_TT"] = new TH1F("KF_Chi2_TT", "Chi2 value of kinFitter (TT); #chi^{2}", 200, 0, 20);
    histo1D["KF_W_mass_orig_TT"] = new TH1F("KF_W_mass_orig_TT", "W mass before kinFitter (TT); m_{W} [GeV]", 250, 0, 500);
    histo1D["KF_top_mass_orig_TT"] = new TH1F("KF_top_mass_orig_TT", "Top mass before kinFitter (TT); m_{t} [GeV]", 400, 0, 800);
    histo1D["KF_top_pt_orig_TT"] = new TH1F("KF_top_pt_orig_TT", "Top pt before kinFitter (TT); p_{T} [GeV]", 400, 0, 800);
    histo1D["KF_dR_lep_b_orig_TT"] = new TH1F("KF_dR_lep_b_orig_TT","Minimal delta R between the lepton and a b jet before kinFitter (for reconstructed TT events); #Delta R(l,b)", 25, 0, 5);
    histo1D["KF_W_mass_corr_TT"] = new TH1F("KF_W_mass_corr_TT", "W mass after kinFitter (TT); m_{W,kf} [GeV]", 250, 0, 500);
    histo1D["KF_top_mass_corr_TT"] = new TH1F("KF_top_mass_corr_TT", "Top mass after kinFitter (TT); m_{t,kf} [GeV]", 400, 0, 800);
    histo1D["KF_top_pt_corr_TT"] = new TH1F("KF_top_pt_corr_TT", "Top pt after kinFitter (TT); p_{T,kf} [GeV]", 400, 0, 800);
    histo1D["KF_dR_lep_b_corr_TT"] = new TH1F("KF_dR_lep_b_corr_TT","Minimal delta R between the lepton and a b jet after kinFitter (for reconstructed TT events); #Delta R(l,b)", 25, 0, 5);
    
    histo1D["KF_top_mass_orig_ex4j_chi2cut5_TT"] = new TH1F("KF_top_mass_orig_ex4j_chi2cut5_TT", "Top mass before kinFitter (TT) (exactly 4 jets - KF chi2 < 5); m_{t} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_ex4j_chi2cut5_TT"] = new TH1F("KF_top_mass_corr_ex4j_chi2cut5_TT", "Top mass after kinFitter (TT) (exactly 4 jets - KF chi2 < 5); m_{t,kf} [GeV]", 400, 0, 800);
    
    histo1D["KF_top_mass_orig_ex4j_chi2cut2_TT"] = new TH1F("KF_top_mass_orig_ex4j_chi2cut2_TT", "Top mass before kinFitter (TT) (exactly 4 jets - KF chi2 < 2); m_{t} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_ex4j_chi2cut2_TT"] = new TH1F("KF_top_mass_corr_ex4j_chi2cut2_TT", "Top mass after kinFitter (TT) (exactly 4 jets - KF chi2 < 2); m_{t,kf} [GeV]", 400, 0, 800);
    
    histo1D["KF_top_mass_orig_ex4j_chi2cut1p5_TT"] = new TH1F("KF_top_mass_orig_ex4j_chi2cut1p5_TT", "Top mass before kinFitter (TT) (exactly 4 jets - KF chi2 < 1.5); m_{t} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_ex4j_chi2cut1p5_TT"] = new TH1F("KF_top_mass_corr_ex4j_chi2cut1p5_TT", "Top mass after kinFitter (TT) (exactly 4 jets - KF chi2 < 1.5); m_{t,kf} [GeV]", 400, 0, 800);
    
    histo1D["KF_jet0_Et_diff_TT"] = new TH1F("KF_jet0_Et_diff_TT", "Et difference after kinFitter for jet0 (TT); E_{T,kf} - E_{T,orig} [GeV]", 400, -50, 50);
    histo1D["KF_jet1_Et_diff_TT"] = new TH1F("KF_jet1_Et_diff_TT", "Et difference after kinFitter for jet1 (TT); E_{T,kf} - E_{T,orig} [GeV]", 400, -50, 50);
    
    histo1D["KF_top_mass_corr_CM"] = new TH1F("KF_top_mass_corr_CM", "Top mass after kinFitter for correct match (CM); m_{t,kf} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_WM"] = new TH1F("KF_top_mass_corr_WM", "Top mass after kinFitter for wrong permutations (WM); m_{t,kf} [GeV]", 400, 0, 800);
    histo1D["KF_top_mass_corr_UM"] = new TH1F("KF_top_mass_corr_UM", "Top mass after kinFitter for no match (UM); m_{t,kf} [GeV]", 400, 0, 800);
    
    histo1D["KF_Chi2_CM"] = new TH1F("KF_Chi2_CM", "Chi2 value of kinFitter (CM); #chi^{2}", 200, 0, 5);
    histo1D["KF_Chi2_WM"] = new TH1F("KF_Chi2_WM", "Chi2 value of kinFitter (WM); #chi^{2}", 200, 0, 5);
    histo1D["KF_Chi2_UM"] = new TH1F("KF_Chi2_UM", "Chi2 value of kinFitter (UM); #chi^{2}", 200, 0, 5);
    histo1D["KF_Chi2_CM_wide"] = new TH1F("KF_Chi2_CM_wide", "Chi2 value of kinFitter (CM); #chi^{2}", 200, 0, 20);
    histo1D["KF_Chi2_WM_wide"] = new TH1F("KF_Chi2_WM_wide", "Chi2 value of kinFitter (WM); #chi^{2}", 200, 0, 20);
    histo1D["KF_Chi2_UM_wide"] = new TH1F("KF_Chi2_UM_wide", "Chi2 value of kinFitter (UM); #chi^{2}", 200, 0, 20);
  }
}

void InitHisto2D()
{
  TH2::SetDefaultSumw2();
  
  InitHisto2DMatch();
  
  histo2D["puSF_vs_nTruePU"] = new TH2F("puSF_vs_nTruePU","nTruePU vs. pu SF (TT); pu SF; nTruePU", 80, 0, 2, 60, 0, 60);
  histo2D["puSF_vs_nVtx"] = new TH2F("puSF_vs_nVtx","nVtx vs. pu SF (TT); pu SF; nVtx", 80, 0, 2, 60, 0, 60);
  histo2D["nVtx_vs_nTruePU"] = new TH2F("nVtx_vs_nTruePU","nTruePU vs. nVtx (TT); nVtx; nTruePU", 60, 0, 60, 60, 0, 60);
  
  /// Reco
  histo2D["dR_lep_b_lep_vs_had_CM"] = new TH2F("dR_lep_b_lep_vs_had_CM","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, correct match); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_lep_b_lep_vs_had_WM"] = new TH2F("dR_lep_b_lep_vs_had_WM","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, wrong match); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_lep_b_lep_vs_had_UM"] = new TH2F("dR_lep_b_lep_vs_had_UM","#DeltaR(l,b_{had}) vs. #Delta R(l,b_{l}) (reco, unmatched); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  
  histo2D["dR_b_light_min_max_CM"]  = new TH2F("dR_b_light_min_max_CM","Maximum delta R vs. minimum delta R between the hadronic b jet and the light jet (aKF, correct match); min #Delta R(b,j); max #Delta R(b,j)", 25, 0, 5, 25, 0, 5);
  histo2D["dR_b_light_min_max_WM"]  = new TH2F("dR_b_light_min_max_WM","Maximum delta R vs. minimum delta R between the hadronic b jet and the light jet (aKF, wrong match); min #Delta R(b,j); max #Delta R(b,j)", 25, 0, 5, 25, 0, 5);
  histo2D["dR_b_light_min_max_UM"]  = new TH2F("dR_b_light_min_max_UM","Maximum delta R vs. minimum delta R between the hadronic b jet and the light jet (aKF, unmatched); min #Delta R(b,j); max #Delta R(b,j)", 25, 0, 5, 25, 0, 5);
  histo2D["dR_b_light_min_dR_light_CM"]  = new TH2F("dR_b_light_min_dR_light_CM","delta R between light jets vs. minimum delta R between the hadronic b jet and the light jet (aKF, correct match); min #Delta R(b,j); #Delta R(j_{1},j_{2})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_b_light_min_dR_light_WM"]  = new TH2F("dR_b_light_min_dR_light_WM","delta R between light jets vs. minimum delta R between the hadronic b jet and the light jet (aKF, wrong match); min #Delta R(b,j); #Delta R(j_{1},j_{2})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_b_light_min_dR_light_UM"]  = new TH2F("dR_b_light_min_dR_light_UM","delta R between light jets vs. minimum delta R between the hadronic b jet and the light jet (aKF, unmatched); min #Delta R(b,j); #Delta R(j_{1},j_{2})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_bW_dR_light_CM"]  = new TH2F("dR_bW_dR_light_CM","delta R between light jets vs. delta R between the hadronic b jet and the W candidate (aKF, correct match); #Delta R(b,W); #Delta R(j_{1},j_{2})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_bW_dR_light_WM"]  = new TH2F("dR_bW_dR_light_WM","delta R between light jets vs. delta R between the hadronic b jet and the W candidate (aKF, wrong match); #Delta R(b,W); #Delta R(j_{1},j_{2})", 25, 0, 5, 25, 0, 5);
  histo2D["dR_bW_dR_light_UM"]  = new TH2F("dR_bW_dR_light_UM","delta R between light jets vs. delta R between the hadronic b jet and the W candidate (aKF, unmatched); #Delta R(b,W); #Delta R(j_{1},j_{2})", 25, 0, 5, 25, 0, 5);
  
  histo2D["ttbar_mass_vs_minMlb_CM"] = new TH2F("ttbar_mass_vs_minMlb_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_WM"] = new TH2F("ttbar_mass_vs_minMlb_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  histo2D["ttbar_mass_vs_minMlb_UM"] = new TH2F("ttbar_mass_vs_minMlb_UM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, unmatched); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  
  histo2D["qqb1_vs_qqb2_CM"] = new TH2F("qqb1_vs_qqb2_CM","Reconstructed top quark mass using the hadronic vs. leptonic b jet (reco bKF, correct match); m_{qqb_{1}} [GeV]; m_{qqb_{2}} [GeV]", 500, 0, 500, 500, 0, 500);
  histo2D["qqb1_vs_qqb2_WM"] = new TH2F("qqb1_vs_qqb2_WM","Reconstructed top quark mass using the hadronic vs. leptonic b jet (reco bKF, wrong match); m_{qqb_{1}} [GeV]; m_{qqb_{2}} [GeV]", 500, 0, 500, 500, 0, 500);
  histo2D["qqb1_vs_qqb2_UM"] = new TH2F("qqb1_vs_qqb2_UM","Reconstructed top quark mass using the hadronic vs. leptonic b jet (reco bKF, unmatched); m_{qqb_{1}} [GeV]; m_{qqb_{2}} [GeV]", 500, 0, 500, 500, 0, 500);
  histo2D["qqb1_vs_qqb2_aKF_CM"] = new TH2F("qqb1_vs_qqb2_aKF_CM","Reconstructed top quark mass using the hadronic vs. leptonic b jet (reco, correct match); m_{qqb_{1}} [GeV]; m_{qqb_{2}} [GeV]", 500, 0, 500, 500, 0, 500);
  histo2D["qqb1_vs_qqb2_aKF_WM"] = new TH2F("qqb1_vs_qqb2_aKF_WM","Reconstructed top quark mass using the hadronic vs. leptonic b jet (reco, wrong match); m_{qqb_{1}} [GeV]; m_{qqb_{2}} [GeV]", 500, 0, 500, 500, 0, 500);
  histo2D["qqb1_vs_qqb2_aKF_UM"] = new TH2F("qqb1_vs_qqb2_aKF_UM","Reconstructed top quark mass using the hadronic vs. leptonic b jet (reco, unmatched); m_{qqb_{1}} [GeV]; m_{qqb_{2}} [GeV]", 500, 0, 500, 500, 0, 500);
  
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCut_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCut_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCut_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCut_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCut_UM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCut_UM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCut_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCut_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCut_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCut_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCut_UM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCut_UM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCutHard_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCutHard_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRlepCutHard_UM"] = new TH2F("ttbar_mass_vs_minMlb_dRlepCutHard_UM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCutHard_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCutHard_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRhadCutHard_UM"] = new TH2F("ttbar_mass_vs_minMlb_dRhadCutHard_UM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCuts_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCuts_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCuts_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCuts_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCuts_UM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCuts_UM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_CM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCutsHard_CM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, correct match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_WM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCutsHard_WM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, wrong permutations); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  //  histo2D["ttbar_mass_vs_minMlb_dRBothCutsHard_UM"] = new TH2F("ttbar_mass_vs_minMlb_dRBothCutsHard_UM","Reconstructed mass of the top quark pair vs. minimal delta R between the lepton and the (supposed to be) leptonic b jet (reco, no match); min(M_{lb}) [GeV]; M_{t#bar{t}} [GeV]", 400, 0, 800, 500, 0, 1000);
  
  /// KinFitter
  if (doKinFit)
  {
    histo2D["KF_W_mass_orig_vs_corr_TT"] = new TH2F("KF_W_mass_orig_vs_corr_TT", "W mass made with KF corrected jets vs. original jets (TT); m_{W,orig} [GeV]; m_{W,kf} [GeV]", 250, 0, 500, 250, 0, 500);
    //histo2D["KF_top_mass_orig_vs_corr_TT"] = new TH2F("KF_top_mass_orig_vs_corr_TT", "Top mass made with KF corrected jets vs. original jets (TT); m_{t,orig} [GeV]; m_{t,kf} [GeV]", 400, 0, 800, 400, 0, 800);
    
    histo2D["KF_top_mass_orig_vs_corr_CM"] = new TH2F("KF_top_mass_orig_vs_corr_CM", "Top mass made with KF corrected jets vs. original jets (CM); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_WM"] = new TH2F("KF_top_mass_orig_vs_corr_WM", "Top mass made with KF corrected jets vs. original jets (WM); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_UM"] = new TH2F("KF_top_mass_orig_vs_corr_UM", "Top mass made with KF corrected jets vs. original jets (UM); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_CM"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_CM", "Top mass made with KF corrected jets vs. original jets (CM); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_WM"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_WM", "Top mass made with KF corrected jets vs. original jets (WM); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_UM"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_UM", "Top mass made with KF corrected jets vs. original jets (UM); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    
    histo2D["KF_top_mass_orig_vs_corr_TT"] = new TH2F("KF_top_mass_orig_vs_corr_TT", "Top mass made with KF corrected jets vs. original jets (TT); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_ST_tW_top"] = new TH2F("KF_top_mass_orig_vs_corr_ST_tW_top", "Top mass made with KF corrected jets vs. original jets (ST_tW_top); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_ST_tW_antitop"] = new TH2F("KF_top_mass_orig_vs_corr_ST_tW_antitop", "Top mass made with KF corrected jets vs. original jets (ST_tW_antitop); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_ST_t_top"] = new TH2F("KF_top_mass_orig_vs_corr_ST_t_top", "Top mass made with KF corrected jets vs. original jets (ST_t_top); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_ST_t_antitop"] = new TH2F("KF_top_mass_orig_vs_corr_ST_t_antitop", "Top mass made with KF corrected jets vs. original jets (ST_t_antitop); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_DYJets"] = new TH2F("KF_top_mass_orig_vs_corr_DYJets", "Top mass made with KF corrected jets vs. original jets (DYJets); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_WJets"] = new TH2F("KF_top_mass_orig_vs_corr_WJets", "Top mass made with KF corrected jets vs. original jets (WJets); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_data"] = new TH2F("KF_top_mass_orig_vs_corr_data", "Top mass made with KF corrected jets vs. original jets (data); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    
    histo2D["KF_top_mass_orig_vs_corr_noCut_TT"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_TT", "Top mass made with KF corrected jets vs. original jets (TT); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_ST_tW_top"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_ST_tW_top", "Top mass made with KF corrected jets vs. original jets (ST_tW_top); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_ST_tW_antitop"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_ST_tW_antitop", "Top mass made with KF corrected jets vs. original jets (ST_tW_antitop); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_ST_t_top"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_ST_t_top", "Top mass made with KF corrected jets vs. original jets (ST_t_top); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_ST_t_antitop"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_ST_t_antitop", "Top mass made with KF corrected jets vs. original jets (ST_t_antitop); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_DYJets"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_DYJets", "Top mass made with KF corrected jets vs. original jets (DYJets); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_WJets"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_WJets", "Top mass made with KF corrected jets vs. original jets (WJets); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_top_mass_orig_vs_corr_noCut_data"] = new TH2F("KF_top_mass_orig_vs_corr_noCut_data", "Top mass made with KF corrected jets vs. original jets (data); m_{t,orig} (GeV); m_{t} (GeV)", 250, 0, 500, 250, 0, 500);
    
    histo2D["KF_W_mass_T_orig_vs_corr_CM"] = new TH2F("KF_W_mass_T_orig_vs_corr_CM", "Transverse W mass made with KF corrected jets vs. original jets (CM); m_{T,W,orig} (GeV); m_{T,W} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_W_mass_T_orig_vs_corr_WM"] = new TH2F("KF_W_mass_T_orig_vs_corr_WM", "Transverse W mass made with KF corrected jets vs. original jets (WM); m_{T,W,orig} (GeV); m_{T,W} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_W_mass_T_orig_vs_corr_UM"] = new TH2F("KF_W_mass_T_orig_vs_corr_UM", "Transverse W mass made with KF corrected jets vs. original jets (UM); m_{T,W,orig} (GeV); m_{T,W} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_W_mass_T_orig_vs_corr_noCut_CM"] = new TH2F("KF_W_mass_T_orig_vs_corr_noCut_CM", "Transverse W mass made with KF corrected jets vs. original jets (CM); m_{T,W,orig} (GeV); m_{T,W} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_W_mass_T_orig_vs_corr_noCut_WM"] = new TH2F("KF_W_mass_T_orig_vs_corr_noCut_WM", "Transverse W mass made with KF corrected jets vs. original jets (WM); m_{T,W,orig} (GeV); m_{T,W} (GeV)", 250, 0, 500, 250, 0, 500);
    histo2D["KF_W_mass_T_orig_vs_corr_noCut_UM"] = new TH2F("KF_W_mass_T_orig_vs_corr_noCut_UM", "Transverse W mass made with KF corrected jets vs. original jets (UM); m_{T,W,orig} (GeV); m_{T,W} (GeV)", 250, 0, 500, 250, 0, 500);
    
    histo2D["KF_jets_Et_diff_TT"] = new TH2F("KF_jets_Et_diff_TT", "Et difference after kinFitter for jet1 in function of jet0 (TT); E_{T,0,kf} - E_{T,0,orig} [GeV]; E_{T,1,kf} - E_{T,1,orig} [GeV]", 400, -50, 50, 400, -50, 50);
    histo2D["KF_W_px_orig_vs_corr_TT"] = new TH2F("KF_W_px_orig_vs_corr_TT", "W p_{x} made with KF corrected jets vs. original jets (TT); p_{x,orig} [GeV]; p_{x,kf} [GeV]", 400, 0, 800, 400, 0, 800);
    histo2D["KF_W_py_orig_vs_corr_TT"] = new TH2F("KF_W_py_orig_vs_corr_TT", "W p_{y} made with KF corrected jets vs. original jets (TT); p_{y,orig} [GeV]; p_{y,kf} [GeV]", 400, 0, 800, 400, 0, 800);
    histo2D["KF_W_pz_orig_vs_corr_TT"] = new TH2F("KF_W_pz_orig_vs_corr_TT", "W p_{z} made with KF corrected jets vs. original jets (TT); p_{z,orig} [GeV]; p_{z,kf} [GeV]", 400, 0, 800, 400, 0, 800);
    
    histo2D["KF_chi2_dR_b_light_min_CM"] = new TH2F("KF_chi2_dR_b_light_min_CM","Minimum delta R between the hadronic b jet and the light jet vs. KF #chi^{2} (correct match); #chi^{2}; min #Delta R(b,j)", 50, 0, 5, 25, 0, 5);
    histo2D["KF_chi2_dR_b_light_min_WM"] = new TH2F("KF_chi2_dR_b_light_min_WM","Minimum delta R between the hadronic b jet and the light jet vs. KF #chi^{2} (wrong match); #chi^{2}; min #Delta R(b,j)", 50, 0, 5, 25, 0, 5);
    histo2D["KF_chi2_dR_b_light_min_UM"] = new TH2F("KF_chi2_dR_b_light_min_UM","Minimum delta R between the hadronic b jet and the light jet vs. KF #chi^{2} (unmatched); #chi^{2}; min #Delta R(b,j)", 50, 0, 5, 25, 0, 5);
    histo2D["KF_chi2_dR_light_CM"] = new TH2F("KF_chi2_dR_light_CM","Delta R between the light jets vs. KF #chi^{2} (correct match); #chi^{2}; #Delta R(j_{1},j_{2})", 50, 0, 5, 25, 0, 5);
    histo2D["KF_chi2_dR_light_WM"] = new TH2F("KF_chi2_dR_light_WM","Delta R between the light jets vs. KF #chi^{2} (wrong match); #chi^{2}; #Delta R(j_{1},j_{2})", 50, 0, 5, 25, 0, 5);
    histo2D["KF_chi2_dR_light_UM"] = new TH2F("KF_chi2_dR_light_UM","Delta R between the light jets vs. KF #chi^{2} (unmatched); #chi^{2}; #Delta R(j_{1},j_{2})", 50, 0, 5, 25, 0, 5);
    histo2D["KF_chi2_dR_bW_CM"] = new TH2F("KF_chi2_dR_bW_CM","Delta R between the hadronic b jet and the W candidate vs. KF #chi^{2} (correct match); #chi^{2}; #Delta R(b,W)", 50, 0, 5, 25, 0, 5);
    histo2D["KF_chi2_dR_bW_WM"] = new TH2F("KF_chi2_dR_bW_WM","Delta R between the hadronic b jet and the W candidate vs. KF #chi^{2} (wrong match); #chi^{2}; #Delta R(b,W)", 50, 0, 5, 25, 0, 5);
    histo2D["KF_chi2_dR_bW_UM"] = new TH2F("KF_chi2_dR_bW_UM","Delta R between the hadronic b jet and the W candidate vs. KF #chi^{2} (unmatched); #chi^{2}; #Delta R(b,W)", 50, 0, 5, 25, 0, 5);
    
    histo2D["KF_chi2_top_mass_CM"] = new TH2F("KF_chi2_top_mass_CM","Top quark mass vs. KF #chi^{2} (correct match); #chi^{2}; m_{t} (GeV)", 150, 0, 15, 300, 0, 600);
    histo2D["KF_chi2_top_mass_WM"] = new TH2F("KF_chi2_top_mass_WM","Top quark mass vs. KF #chi^{2} (wrong match); #chi^{2}; m_{t} (GeV)", 150, 0, 15, 300, 0, 600);
    histo2D["KF_chi2_top_mass_UM"] = new TH2F("KF_chi2_top_mass_UM","Top quark mass vs. KF #chi^{2} (no match); #chi^{2}; m_{t} (GeV)", 150, 0, 15, 300, 0, 600);
    histo2D["KF_top_mass_chi2_CM"] = new TH2F("KF_top_mass_chi2_CM","KF #chi^{2} vs. top quark mass (correct match); m_{t} (GeV); #chi^{2}", 200, 0, 400, 150, 0, 15);
    histo2D["KF_top_mass_chi2_WM"] = new TH2F("KF_top_mass_chi2_WM","KF #chi^{2} vs. top quark mass (wrong match); m_{t} (GeV); #chi^{2}", 200, 0, 400, 150, 0, 15);
    histo2D["KF_top_mass_chi2_UM"] = new TH2F("KF_top_mass_chi2_UM","KF #chi^{2} vs. top quark mass (no match); m_{t} (GeV); #chi^{2}", 200, 0, 400, 150, 0, 15);
    histo2D["KF_top_mass_chi2_logY_CM"]  = new TH2F("KF_top_mass_chi2_logY_CM","KF #chi^{2} vs. top quark mass (correct match); m_{t} (GeV); Log_{10} #chi^{2}", 200, 0, 400, 50, 0.000001, 1);
    histo2D["KF_top_mass_chi2_logY_WM"]  = new TH2F("KF_top_mass_chi2_logY_WM","KF #chi^{2} vs. top quark mass (wrong match); m_{t} (GeV); Log_{10} #chi^{2}", 200, 0, 400, 50, 0.000001, 1);
    histo2D["KF_top_mass_chi2_logY_UM"]  = new TH2F("KF_top_mass_chi2_logY_UM","KF #chi^{2} vs. top quark mass (no match); m_{t} (GeV); Log_{10} #chi^{2}", 200, 0, 400, 50, 0.000001, 1);
    
    histo2D["KF_top_mass_chi2_noCut_logY_CM"]  = new TH2F("KF_top_mass_chi2_noCut_logY_CM","KF #chi^{2} vs. top quark mass (correct match); m_{t} (GeV); Log_{10} #chi^{2}", 200, 0, 400, 50, 0.000001, 1);
    histo2D["KF_top_mass_chi2_noCut_logY_WM"]  = new TH2F("KF_top_mass_chi2_noCut_logY_WM","KF #chi^{2} vs. top quark mass (wrong match); m_{t} (GeV); Log_{10} #chi^{2}", 200, 0, 400, 50, 0.000001, 1);
    histo2D["KF_top_mass_chi2_noCut_logY_UM"]  = new TH2F("KF_top_mass_chi2_noCut_logY_UM","KF #chi^{2} vs. top quark mass (no match); m_{t} (GeV); Log_{10} #chi^{2} ", 200, 0, 400, 50, 0.000001, 1);
  }
}

void InitHisto1DMatch()
{
  TH1::SetDefaultSumw2();
  
  histo1D["matched_W_mass_reco"] = new TH1F("matched_W_mass_reco","Reconstructed hadronic W mass of matched events; M_{W} [GeV]", 125, 0, 250);
  histo1D["matched_top_mass_reco"] = new TH1F("matched_top_mass_reco","Reconstructed top mass of matched events; M_{t} [GeV]", 175, 50, 400);
  histo1D["matched_top_mass_gen"] = new TH1F("matched_top_mass_gen","Generated top mass of matched events; M_{t} [GeV]", 1200, 150, 190);
  
  histo1D["matched_red_top_mass_TT_partons"] = new TH1F("matched_red_top_mass_TT_partons","Reduced top mass for matched TT sample (using matched partons); m_{r}", 800, 0.8, 1.2);
  histo1D["matched_red_top_mass_TT_jets"] = new TH1F("matched_red_top_mass_TT_jets","Reduced top mass for matched TT sample (using jets from matched partons); m_{r}", 640, 0.4, 2.);
  
  histo1D["matched_mlb_corr"]  = new TH1F("matched_mlb_corr","Reconstructed leptonic top mass using correctly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["matched_mlb_wrong"] = new TH1F("matched_mlb_wrong","Reconstructed leptonic top mass using wrongly matched events; M_{lb} [GeV]", 80, 0, 800);
  histo1D["matched_ttbar_mass_corr"] = new TH1F("matched_ttbar_mass_corr","Reconstructed mass of the top quark pair using correctly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["matched_ttbar_mass_wrong"] = new TH1F("matched_ttbar_mass_wrong","Reconstructed mass of the top quark pair using wrongly matched events; M_{t#bar{t}} [GeV]", 100, 0, 1000);
  histo1D["matched_dR_lep_b_corr"] = new TH1F("matched_dR_lep_b_corr","Delta R between the lepton and the leptonic b jet for matched events; #Delta R(l,b)", 25, 0, 5);
  histo1D["matched_dR_lep_b_wrong"] = new TH1F("matched_dR_lep_b_wrong","Delta R between the lepton and the hadronic b jet for matched events; #Delta R(l,b_{had})", 25, 0, 5);
  
  if (doKinFit) histo1D["matched_top_mass_reco_aKF"] = new TH1F("matched_top_mass_reco_aKF", "Top mass after kinFitter for matched events; m_{t,kf} [GeV]", 80, 0, 400);
}

void InitHisto2DMatch()
{
  TH2::SetDefaultSumw2();
  
  /// Matched events
  histo2D["matched_mlb_corr_mlb_wrong"] = new TH2F("matched_mlb_corr_mlb_wrong","Wrongly constructed M_{lb} vs. correctly constructed M_{lb}; M_{lb_{lep}}; M_{lb_{had}}", 80, 0, 800, 80, 0, 800);
  histo2D["matched_dR_lep_b_corr_dR_lep_b_wrong"] = new TH2F("matched_dR_lep_b_corr_dR_lep_b_wrong","Wrongly constructed dR(l,b) vs. correctly constructed dR(l,b); #Delta R(l,b_{lep}); #Delta R(l,b_{had})", 25, 0, 5, 25, 0, 5);
  histo2D["matched_mlb_dR_lep_b_corr"] = new TH2F("matched_mlb_dR_lep_b_corr","dR(l,b) vs. M_{lb}; M_{lb_{lep}}; #Delta R(l,b_{lep})", 80, 0, 800, 25, 0, 5);
  histo2D["matched_mlb_dR_lep_b_wrong"] = new TH2F("matched_mlb_dR_lep_b_wrong","dR(l,b) vs. M_{lb}, both wrongly matched; M_{lb_{had}}; #Delta R(l,b_{had})", 80, 0, 800, 25, 0, 5);
}

void InitLikelihoodPlots()
{
  for (int iMass = 0; iMass < nLikeMasses; iMass++)
  {
    histo1DLike["loglike_vs_width_m"+DotReplace(redTopMassArray[iMass])] = new TH1D(("loglike_vs_width_m"+DotReplace(redTopMassArray[iMass])).c_str(), ("loglike_vs_width_m"+DotReplace(redTopMassArray[iMass])+"; width (a.u.); -log(likelihood)").c_str(), nWidthsLike, -0.5, nWidthsLike-0.5);
//     histo1DLike["loglike_vs_width_CM_m"+DotReplace(redTopMassArray[iMass])] = new TH1D(("loglike_vs_width_CM_m"+DotReplace(redTopMassArray[iMass])).c_str(), ("loglike_vs_width_CM_m"+DotReplace(redTopMassArray[iMass])+"; width (a.u.); -log(likelihood)").c_str(), nWidthsLike, -0.5, nWidthsLike-0.5);
//     histo1DLike["loglike_vs_width_WM_m"+DotReplace(redTopMassArray[iMass])] = new TH1D(("loglike_vs_width_WM_m"+DotReplace(redTopMassArray[iMass])).c_str(), ("loglike_vs_width_WM_m"+DotReplace(redTopMassArray[iMass])+"; width (a.u.); -log(likelihood)").c_str(), nWidthsLike, -0.5, nWidthsLike-0.5);
//     histo1DLike["loglike_vs_width_UM_m"+DotReplace(redTopMassArray[iMass])] = new TH1D(("loglike_vs_width_UM_m"+DotReplace(redTopMassArray[iMass])).c_str(), ("loglike_vs_width_UM_m"+DotReplace(redTopMassArray[iMass])+"; width (a.u.); -log(likelihood)").c_str(), nWidthsLike, -0.5, nWidthsLike-0.5);
  }
}

void TruthMatching(vector<TLorentzVector> partons, vector<TLorentzVector> selectedJets, pair<unsigned int, unsigned int> *MCPermutation)  /// MCPermutation: 0,1 hadronic W jet; 2 hadronic b jet; 3 leptonic b jet
{
  JetPartonMatching matching = JetPartonMatching(partons, selectedJets, 2, true, true, 0.3);  // partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
  
  if (matching.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matching.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matching.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
  }
  
  if (verbose > 2)
    cout << "Matching done" << endl;
  
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int j = JetPartonPair[i].second;
    
    if ( fabs(mc_pdgId[partonId[j]]) < 6 )  //light/b quarks, 6 should stay hardcoded
    {
      if (isTTbar)
      {
        if ( ( muPlusFromTop && mc_mother[partonId[j]] == -24 /*&& mc_granny[partonId[j]] == -pdgID_top*/ )
          || ( muMinusFromTop && mc_mother[partonId[j]] == 24 /*&& mc_granny[partonId[j]] == pdgID_top*/ ) )  // if mu+, check if mother of particle is W- and granny tbar --> then it is a quark from W- decay
        {
          if (verbose > 3)
            cout << "Light jet: " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;

          if (MCPermutation[0].first == 9999)
          {
            MCPermutation[0] = JetPartonPair[i];
          }
          else if (MCPermutation[1].first == 9999)
          {
            MCPermutation[1] = JetPartonPair[i];
          }
          else
          {
            cerr << "Found a third jet coming from a W boson which comes from a top quark..." << endl;
            cerr << " -- muMinusFromMtop: " << muMinusFromTop << " muPlusFromMtop: " << muPlusFromTop << endl;
            cerr << " -- pdgId: " << mc_pdgId[partonId[j]] << " mother: " << mc_mother[partonId[j]] << " granny: " << mc_granny[partonId[j]] << " Pt: " << mc_pt[partonId[j]] << endl;
  //          cerr << " -- ievt: " << ievt << endl;
            exit(1);
          }
        }
      }
//       else   /// normally only ST tW should give results
//       {
//         if ( ( foundMuPlus && mc_mother[partonId[j]] == -24 && mc_granny[partonId[j]] == -pdgID_top )
//           || ( foundMuMinus && mc_mother[partonId[j]] == 24 && mc_granny[partonId[j]] == pdgID_top ) )  // if mu+, check if mother of particle is W- and granny tbar --> then it is a quark from W- decay
//         {
//           if (MCPermutation[0].first == 9999)
//           {
//             MCPermutation[0] = JetPartonPair[i];
//           }
//           else if (MCPermutation[1].first == 9999)
//           {
//             MCPermutation[1] = JetPartonPair[i];
//           }
//           else
//           {
//             cerr << "Found a third jet coming from a W boson which comes from a top quark..." << endl;
//             cerr << " -- muMinusFromMtop: " << muMinusFromTop << " muPlusFromMtop: " << muPlusFromTop << endl;
//             cerr << " -- pdgId: " << mc_pdgId[partonId[j]] << " mother: " << mc_mother[partonId[j]] << " granny: " << mc_granny[partonId[j]] << " Pt: " << mc_pt[partonId[j]] << endl;
//   //          cerr << " -- ievt: " << ievt << endl;
//             exit(1);
//           }
//         }
//       }
    }
    if ( fabs(mc_pdgId[partonId[j]]) == 5 )
    {
      if ( ( (muPlusFromTop || foundMuPlus) && mc_mother[partonId[j]] == -pdgID_top )
          || ( (muMinusFromTop || foundMuMinus) && mc_mother[partonId[j]] == pdgID_top ) )  // if mu+ (top decay leptonic) and mother is antitop ---> hadronic b
      {
        if (verbose > 3)
          cout << "b jet:     " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;
        
        MCPermutation[2] = JetPartonPair[i];
      }
      else if ( ( (muPlusFromTop || foundMuPlus) && mc_mother[partonId[j]] == pdgID_top )
               || ( (muMinusFromTop || foundMuMinus) && mc_mother[partonId[j]] == -pdgID_top ) )  // if mu+ (top decay leptonic) and mother is top ---> leptonic b
      {
        if (verbose > 3)
          cout << "b jet:     " << j << "  Status: " << mc_status[partonId[j]] << "  pdgId: " << mc_pdgId[partonId[j]] << "  Mother: " << mc_mother[partonId[j]] << "  Granny: " << mc_granny[partonId[j]] << "  Pt: " << mc_pt[partonId[j]] << "  Eta: " << mc_eta[partonId[j]] << "  Phi: " << mc_phi[partonId[j]] << "  Mass: " << mc_M[partonId[j]] << endl;
        
        MCPermutation[3] = JetPartonPair[i];
      }
    }
  }  /// End loop over Jet Parton Pairs
  
  
  if ( MCPermutation[0].first != 9999 && MCPermutation[1].first != 9999 && MCPermutation[2].first != 9999 && MCPermutation[3].first != 9999 )
  {
    all4PartonsMatched = true;
    nofMatchedEvents++;
    if (MCPermutation[0].first < 4 && MCPermutation[1].first < 4 && MCPermutation[2].first < 4 && MCPermutation[3].first < 4)
      all4JetsMatched_MCdef_ = true;
  }
  else if (verbose > 3) cout << "Size JetPartonPair: " << JetPartonPair.size() << ". Not all partons matched!" << endl;
  
  if ( MCPermutation[0].first != 9999 && MCPermutation[1].first != 9999 && MCPermutation[2].first != 9999 )
  {
    hadronicTopJetsMatched = true;
    nofHadrMatchedEvents++;
    if ( MCPermutation[0].first < 4 && MCPermutation[1].first < 4 && MCPermutation[2].first < 4 )
      hadronicTopJetsMatched_MCdef_ = true;
  }
  
  if ( genmuon != -9999 && ROOT::Math::VectorUtil::DeltaR(mcParticles[genmuon], selectedLepton[0]) < 0.1 )
    muonmatched = true;
}

void ClearMetaData()
{
  nEvents = 0;
  nEventsSel = 0;
  for (int i = 0; i < 10; i++)
  {
    cutFlow[i] = 0;
    cutFlow2[i] = 0;
    cutFlowWeighted[i] = 0;
    cutFlow2Weighted[i] = 0;
  }
  appliedJER = 999;
  appliedJES = 999;
  appliedPU = 999;
  if (isData)
  {
    nofEventsRunB = 0;
    nofEventsRunCD = 0;
    nofEventsRunEF = 0;
    nofEventsRunG = 0;
    nofEventsRunH = 0;
    nofSelEventsRunB = 0;
    nofSelEventsRunCD = 0;
    nofSelEventsRunEF = 0;
    nofSelEventsRunG = 0;
    nofSelEventsRunH = 0;
  }
  else
  {
    nofTTEventsWithoutAGenTop = 0;
    sumWeight1001 = 0.;
    sumWeight1002 = 0.;
    sumWeight1003 = 0.;
    sumWeight1004 = 0.;
    sumWeight1005 = 0.;
    sumWeight1007 = 0.;
    sumWeight1009 = 0.;
    sumUpFragWeight = 0.;
    sumCentralFragWeight = 0.;
    sumDownFragWeight = 0.;
    sumPetersonFragWeight = 0.;
    sumSemilepbrUp = 0.;
    sumSemilepbrDown = 0.;
    for (Int_t i = 0; i < 100; i++)
    {
      sumPdfWeights[i] = 0.;
    }
    sumPdfAlphaSUp = 0.;
    sumPdfAlphaSDown = 0.;
  }
  
  strSyst = "";
  nEventsDataSet = 0;
  xSection = 1.;
  eqLumi = 1.;
  lumiWeight = 1.;
  relativeSF = 1.;
  
  nofHardSelected = 0;
  nofMETCleaned = 0;
  nofAfterDRmincut = 0;
  nofAfterDRmaxcut = 0;
  nofAfterLastCut = 0;
  nofTTsemilep = 0;
  nofTTdilep = 0;
  nofTThadr = 0;
  nofMatchedEvents = 0;
  nofHadrMatchedEvents = 0;
  nofHadrMatchedEventsAKF = 0;
  nofCorrectlyMatched = 0;
  nofNotCorrectlyMatched = 0;
  nofUnmatched = 0;
  nofCorrectlyMatchedAKF = 0;
  nofNotCorrectlyMatchedAKF = 0;
  nofUnmatchedAKF = 0;
  nofCorrectlyMatchedAKFNoCut = 0;
  nofNotCorrectlyMatchedAKFNoCut = 0;
  nofUnmatchedAKFNoCut = 0;
  nofAcceptedKFit = 0;
  nofAcceptedKFitWeighted = 0.;
  
  toyMax = 1.;
}

void ClearLeaves()
{
  run_num = -1;
  evt_num = -1;
  lumi_num = -1;
  nvtx = -1;
  npu = -1;
  rho = -999.;
  isTrigged = false;
  hasExactly4Jets = false;
  hasJetLeptonCleaning = false;
  nMuons = -1;
  muon_charge[0] = 0.;
  muon_pt[0] = 0.;
  muon_phi[0] = 0.;
  muon_eta[0] = 0.;
  muon_E[0] = 0.;
  muon_M[0] = 0.;
  muon_d0[0] = 999.;
  muon_chargedHadronIso[0] = 999.;
  muon_neutralHadronIso[0] = 999.;
  muon_photonIso[0] = 999.;
  muon_puChargedHadronIso[0] = 999.;
  muon_relIso[0] = 999.;
  muon_pfIso[0] = 999.;
  nJets = -1;
  for (Int_t i = 0; i < 20; i++)
  {
    jet_charge[i] = 0.;
    jet_nConstituents[i] = -1;
    jet_nChConstituents[i] = -1;
    jet_pt[i] = 0.;
    jet_phi[i] = 0.;
    jet_eta[i] = 0.;
    jet_E[i] = 0.;
    jet_M[i] = 0.;
    jet_bdiscr[i] = -1.;
    jet_hadronFlavour[i] = -1.;
  }
  met_px = 0.;
  met_py = 0.;
  met_pt = 0.;
  met_phi = 0.;
  met_eta = 0.;
  met_Et = 0.;
  met_E = 0.;
  met_corr_px = 0.;
  met_corr_py = 0.;
  met_corr_pt = 0.;
  met_corr_phi = 0.;
  met_corr_eta = 0.;
  met_corr_Et = 0.;
  met_corr_E = 0.;
  
  if (! isData)
  {
    nGenJets = -1;
    for (Int_t i = 0; i < 40; i++)
    {
      genJet_pdgId[i] = 0;
      genJet_charge[i] = 0.;
      genJet_pt[i] = 0.;
      genJet_phi[i] = 0.;
      genJet_eta[i] = 0.;
      genJet_E[i] = 0.;
      genJet_M[i] = 0.;
    }
    nMCParticles = -1;
    for (Int_t i = 0; i < 200; i++)
    {
      mc_status[i] = -1;
      mc_pdgId[i] = 0;
      mc_mother[i] = 0;
      mc_granny[i] = 0;
      mc_pt[i] = 0.;
      mc_phi[i] = 0.;
      mc_eta[i] = 0.;
      mc_E[i] = 0.;
      mc_M[i] = 0.;
    }
    hasGenTop = false;
    hasGenAntiTop = false;
    
    weight1001 = 1.;
    weight1002 = 1.;
    weight1003 = 1.;
    weight1004 = 1.;
    weight1005 = 1.;
    weight1007 = 1.;
    weight1009 = 1.;
    upFragWeight = 1.;
    centralFragWeight = 1.;
    downFragWeight = 1.;
    petersonFragWeight = 1.;
    semilepbrUp = 1.;
    semilepbrDown = 1.;
    for (Int_t i = 0; i < 100; i++)
    {
      pdfWeights[i] = 1.;
    }
    pdfAlphaSUp = 1.;
    pdfAlphaSDown = 1.;
    btagSF = 1.;
    btagSF_up = 1.;
    btagSF_down = 1.;
    puSF = 1.;
    puSF_up = 1.;
    puSF_down = 1.;
    muonIdSF_BCDEF[0] = 1.;
    muonIdSF_GH[0] = 1.;
    muonIdSF_up_BCDEF[0] = 1.;
    muonIdSF_up_GH[0] = 1.;
    muonIdSF_down_BCDEF[0] = 1.;
    muonIdSF_down_GH[0] = 1.;
    muonIsoSF_BCDEF[0] = 1.;
    muonIsoSF_GH[0] = 1.;
    muonIsoSF_up_BCDEF[0] = 1.;
    muonIsoSF_up_GH[0] = 1.;
    muonIsoSF_down_BCDEF[0] = 1.;
    muonIsoSF_down_GH[0] = 1.;
    muonTrigSF_BCDEF[0] = 1.;
    muonTrigSF_GH[0] = 1.;
    muonTrigSF_up_BCDEF[0] = 1.;
    muonTrigSF_up_GH[0] = 1.;
    muonTrigSF_down_BCDEF[0] = 1.;
    muonTrigSF_down_GH[0] = 1.;
    muonTrackSF_eta[0] = 1.;
    muonTrackSF_aeta[0] = 1.;
    muonTrackSF_nPV[0] = 1.;
  }
  else
  {
    isDataRunB = false;
    isDataRunC = false;
    isDataRunD = false;
    isDataRunE = false;
    isDataRunF = false;
    isDataRunG = false;
    isDataRunH = false;
  }
}

void ClearTLVs()
{
  muon.Clear();
  jet.Clear();
  mcpart.Clear();
  WCandidate.Clear();
  selectedLepton.clear();
  selectedJets.clear();
  selectedBJets.clear();
  selectedJetsAKF.clear();
  selectedJetsKFcorrected.clear();
  mcParticles.clear();
  partons.clear();
  partonsMatched.clear();
  jetsMatched.clear();
}

void ClearMatching()
{
  partons.clear();
  partonsMatched.clear();
  jetsMatched.clear();
  
  doMatching = true;
  all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
  all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
  hadronicTopJetsMatched = false;
  hadronicTopJetsMatched_MCdef_ = false;
  for (int i = 0; i < 4; i++)
  {
    MCPermutation[i] = pair<unsigned int,unsigned int>(9999,9999);
  }
  topQuark = -9999;
  antiTopQuark = -9999;
  genmuon = -9999;
  muonmatched = false;
  foundMuPlus = false;
  foundMuMinus = false;
  muPlusFromTop = false;
  muMinusFromTop = false;
  partonId.clear();
  foundB.clear();
  foundLight.clear();
  nPartons = -1;
  
  matched_W_mass_q = -1.;
  matched_top_mass_q = -1.;
  matched_W_mass_j = -1.;
  matched_top_mass_j = -1.;
  matched_mlb_corr = -1.;
  matched_ttbarMass_corr = -1.;
  matched_dR_lep_b_corr = 999.;
  matched_mlb_wrong = -1.;
  matched_ttbarMass_wrong = -1.;
  matched_dR_lep_b_wrong = 999.;
}

void ClearVars()
{
  ClearMatching();
  
  scaleFactor = 1.;
  widthSF = 1.;
  thisLeptonSF = 1.;
  thisLeptonIdSF = 1.;
  thisLeptonIsoSF = 1.;
  thisLeptonTrigSF = 1.;
  topPtRewSF = 1.;
  topPtSF = 1.;
  antiTopPtSF = 1.;
  bJetId.clear();
  bdiscrTop = -99.;
  bdiscrTop2 = -99.;
  tempbdiscr = -99.;
  labelB1 = -9999;
  labelB2 = -9999;
  for (int i = 0; i < 4; i++)
  {
    labelsReco[i] = -9999;
  }
  massHadTopQ = 0.01;
  massLepTopQ = 0.01;
  catSuffix = "";
  isCM = false;
  isWM = false;
  isUM = false;
  isUM_TTsemilep = false;
  isUM_TTother = false;
  isTTsemilep = false;
  isTTother = false;
  doneKinFit = false;
  kFitVerbosity = false;
  kFitChi2 = 99.;
  //toyValue = -1.;
  for (int i = 0; i < nPsExps; i++)
  {
    toyValues[i] = 1.;
  }
  
  M3 = -1.;
  M3_aKF = -1.;
  Ht = -1.;
  Ht_aKF = -1;
  min_Mlb = 9999.;
  dRLepB = -1.;
  reco_W_mass_bKF = -1.;
  reco_top_mass_bKF = -1.;
  reco_top_mass_alt_bKF = -1.;
  reco_W_pt_bKF = -1.;
  reco_top_pt_bKF = -1.;
  reco_top_pt_alt_bKF = -1.;
  reco_mlb_bKF = -1.;
  reco_dRLepB_lep_bKF = -1.;
  reco_dRLepB_had_bKF = -1.;
  reco_dRlight_bKF = -1.;
  reco_dRblight_min_bKF = -1.;
  reco_dRblight_max_bKF = -1.;
  reco_dRbW_bKF = -1.;
  reco_dRblight_qsum_bKF = -1.;
  reco_dRbb_bKF = -1.;
  reco_dPhi_bb_bKF = -1.;
  reco_dPhi_light_bKF = -1.;
  reco_dPhi_bW_bKF = -1.;
  reco_ttbar_mass_bKF = -1.;
  redTopMass_old_bKF = -1.;
  redTopMass_new_bKF = -1.;
  redTopMass_bKF = -1.;
  redMlbMass_bKF = -1.;
  reco_mbjj_div_mjj_bKF = -1.;
  reco_mTW_bKF = -1.;
  reco_W_mass_aKF = -1.;
  reco_top_mass_aKF = -1.;
  reco_top_mass_alt_aKF = -1.;
  reco_W_pt_aKF = -1.;
  reco_top_pt_aKF = -1.;
  reco_top_pt_alt_aKF = -1.;
  reco_mlb_aKF = -1.;
  reco_dRLepB_lep_aKF = -1.;
  reco_dRLepB_had_aKF = -1.;
  reco_dRlight_aKF = -1.;
  reco_dRblight_min_aKF = -1.;
  reco_dRblight_max_aKF = -1.;
  reco_dRbW_aKF = -1.;
  reco_dRblight_qsum_aKF = -1.;
  reco_dRbb_aKF = -1.;
  reco_dPhi_bb_aKF = -1.;
  reco_dPhi_light_aKF = -1.;
  reco_dPhi_bW_aKF = -1.;
  reco_ttbar_mass_aKF = -1.;
  redTopMass_old = -1.;
  redTopMass_new = -1.;
  redTopMass = -1.;
  redMlbMass = -1.;
  reco_mbjj_div_mjj = -1.;
  reco_new_var = -1.;
  reco_mTW_aKF = -1.;
  loglike_per_evt.clear();
}

void ClearObjects()
{
  ClearLeaves();
  ClearTLVs();
  ClearVars();
}

void FillControlPlots(vector<Dataset *> datasets, int d)
{
  FillControlPlots(datasets, d, "");
}

void FillControlPlots(vector<Dataset *> datasets, int d, string suffix)
{
  M3 = (selectedJets[0] + selectedJets[1] + selectedJets[2]).M();
  Ht = selectedJets[0].Pt() + selectedJets[1].Pt() + selectedJets[2].Pt() + selectedJets[3].Pt();
  
  double Mlb_temp = -1;
  for (unsigned int i = 0; i < selectedBJets.size(); i++)
  {
    Mlb_temp = (selectedLepton[0] + selectedBJets[i]).M();
    if ( Mlb_temp < min_Mlb )
    {
      min_Mlb = Mlb_temp;
      dRLepB = ROOT::Math::VectorUtil::DeltaR( selectedBJets[i], selectedLepton[0]);
    }
  }
  
  MSPlotCP["muon_pT"+suffix]->Fill(selectedLepton[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["muon_eta"+suffix]->Fill(selectedLepton[0].Eta(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["muon_phi"+suffix]->Fill(selectedLepton[0].Phi(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["muon_relIso"+suffix]->Fill(muon_relIso[0], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["muon_d0"+suffix]->Fill(muon_d0[0], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  if (suffix.find("aKF") != std::string::npos )
  {
    MSPlotCP["leadingJet_pT"+suffix]->Fill(selectedJetsAKF[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet2_pT"+suffix]->Fill(selectedJetsAKF[1].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet3_pT"+suffix]->Fill(selectedJetsAKF[2].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet4_pT"+suffix]->Fill(selectedJetsAKF[3].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    Ht_aKF = selectedJetsAKF[0].Pt() + selectedJetsAKF[1].Pt() + selectedJetsAKF[2].Pt() + selectedJetsAKF[3].Pt();
    MSPlotCP["Ht_4leadingJets"+suffix]->Fill(Ht_aKF, datasets[d], true, lumiWeight*scaleFactor*widthSF);
    M3_aKF = (selectedJetsAKF[0] + selectedJetsAKF[1] + selectedJetsAKF[2]).M();
    MSPlotCP["M3"+suffix]->Fill(M3_aKF, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  }
  else
  {
    MSPlotCP["leadingJet_pT"+suffix]->Fill(selectedJets[0].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet2_pT"+suffix]->Fill(selectedJets[1].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet3_pT"+suffix]->Fill(selectedJets[2].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["jet4_pT"+suffix]->Fill(selectedJets[3].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["Ht_4leadingJets"+suffix]->Fill(Ht, datasets[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlotCP["M3"+suffix]->Fill(M3, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  }
  MSPlotCP["met_pT"+suffix]->Fill(met_pt, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_eta"+suffix]->Fill(met_eta, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_phi"+suffix]->Fill(met_phi, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_corr_pT"+suffix]->Fill(met_corr_pt, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_corr_eta"+suffix]->Fill(met_corr_eta, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["met_corr_phi"+suffix]->Fill(met_corr_phi, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  MSPlotCP["min_Mlb"+suffix]->Fill(min_Mlb, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["dR_Lep_B"+suffix]->Fill(dRLepB, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  MSPlotCP["nJets"+suffix]->Fill(selectedJets.size(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["nBJets"+suffix]->Fill(selectedBJets.size(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  MSPlotCP["CSVv2Discr_leadingJet"+suffix]->Fill(jet_bdiscr[0], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["CSVv2Discr_jet2"+suffix]->Fill(jet_bdiscr[1], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["CSVv2Discr_jet3"+suffix]->Fill(jet_bdiscr[2], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["CSVv2Discr_jet4"+suffix]->Fill(jet_bdiscr[3], datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  int labelB = -1;
  double highestBDiscr = -999.;
  for (int iJet = 0; iJet < selectedJets.size(); iJet++)
  {
    if ( suffix.find("aKF") != std::string::npos )
      MSPlotCP["jet_pT_allJets"+suffix]->Fill(selectedJetsAKF[iJet].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    else
      MSPlotCP["jet_pT_allJets"+suffix]->Fill(selectedJets[iJet].Pt(), datasets[d], true, lumiWeight*scaleFactor*widthSF);
    
    MSPlotCP["CSVv2Discr_allJets"+suffix]->Fill(jet_bdiscr[iJet], datasets[d], true, lumiWeight*scaleFactor*widthSF);
    if ( jet_bdiscr[iJet] > highestBDiscr )
    {
      highestBDiscr = jet_bdiscr[iJet];
      labelB = iJet;
    }
  }
  MSPlotCP["CSVv2Discr_highest"+suffix]->Fill(highestBDiscr, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  MSPlotCP["CSVv2Discr_jetNb"+suffix]->Fill(labelB, datasets[d], true, lumiWeight*scaleFactor*widthSF);
  
  if (isTTbar)
  {
    pair<unsigned int,unsigned int> hasMinDR = pair<unsigned int,unsigned int>(9999,9999);
    double tempMin = 999.;
    double tempDR;
    for (unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
    {
      for (unsigned int jJet = iJet+1; jJet < selectedJets.size(); jJet++)
      {
        tempDR = ROOT::Math::VectorUtil::DeltaR( selectedJets[iJet], selectedJets[jJet]);
        if ( tempDR < tempMin)
        {
          tempMin = tempDR;
          hasMinDR = pair<unsigned int,unsigned int>(iJet,jJet);
        }
      }
    }
    tempDR = ROOT::Math::VectorUtil::DeltaR( selectedJets[hasMinDR.first], selectedJets[hasMinDR.second]);
    histo1D["dR_jets_min"]->Fill(tempDR);
  }
}

void FillMatchingPlots()
{
  if (hadronicTopJetsMatched)
  {  
    histo1D["matched_W_mass_reco"]->Fill(matched_W_mass_j, widthSF);
    histo1D["matched_top_mass_reco"]->Fill(matched_top_mass_j, widthSF);
    histo1D["matched_top_mass_gen"]->Fill(matched_top_mass_q, widthSF);
    
    histo1D["matched_red_top_mass_TT_partons"]->Fill(matched_top_mass_q/aveTopMass[0], widthSF);
    histo1D["matched_red_top_mass_TT_jets"]->Fill(matched_top_mass_j/aveTopMass[1], widthSF);
    
    if (doKinFit)
    {
      histo1D["matched_top_mass_reco_aKF"]->Fill(matched_top_mass_j_akF, widthSF);
    }
    
    if (all4PartonsMatched && muonmatched)
    {
      matched_mlb_corr = (selectedLepton[0] + jetsMatched[3]).M();  // lept b
      matched_mlb_wrong = (selectedLepton[0] + jetsMatched[2]).M();  // hadr b
      matched_ttbarMass_corr = matched_top_mass_j + matched_mlb_corr;
      matched_ttbarMass_wrong = matched_top_mass_j + matched_mlb_wrong;
      matched_dR_lep_b_corr = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[3]);  // lept b
      matched_dR_lep_b_wrong = ROOT::Math::VectorUtil::DeltaR(selectedLepton[0], jetsMatched[2]);  // hadr b
      
      histo1D["matched_mlb_corr"]->Fill(matched_mlb_corr);
      histo1D["matched_mlb_wrong"]->Fill(matched_mlb_wrong);
      histo1D["matched_ttbar_mass_corr"]->Fill(matched_ttbarMass_corr);
      histo1D["matched_ttbar_mass_wrong"]->Fill(matched_ttbarMass_wrong);
      histo1D["matched_dR_lep_b_corr"]->Fill(matched_dR_lep_b_corr);
      histo1D["matched_dR_lep_b_wrong"]->Fill(matched_dR_lep_b_wrong);
      
      histo2D["matched_mlb_corr_mlb_wrong"]->Fill(matched_mlb_corr, matched_mlb_wrong);
      histo2D["matched_dR_lep_b_corr_dR_lep_b_wrong"]->Fill(matched_dR_lep_b_corr, matched_dR_lep_b_wrong);
      histo2D["matched_mlb_dR_lep_b_corr"]->Fill(matched_mlb_corr, matched_dR_lep_b_corr);
      histo2D["matched_mlb_dR_lep_b_wrong"]->Fill(matched_mlb_wrong, matched_dR_lep_b_wrong);
    }  // end all4PartonsMatched && muonMatched
  }  // end hadronicTopJetsMatched
}

void FillKinFitPlots(bool doneKinFit)
{
  if (isTTbar)
  {
    if (! doneKinFit)
    {
      histo1D["KF_W_mass_orig_TT"]->Fill(reco_W_mass_bKF);
      histo1D["KF_top_mass_orig_TT"]->Fill(reco_top_mass_bKF);
      histo1D["KF_top_pt_orig_TT"]->Fill(reco_top_pt_bKF);
      histo1D["KF_dR_lep_b_orig_TT"]->Fill(reco_dRLepB_lep_bKF);
    }
    else
    {
      histo1D["KF_W_mass_corr_TT"]->Fill(reco_W_mass_aKF);
      histo2D["KF_W_mass_orig_vs_corr_TT"]->Fill(reco_W_mass_bKF, reco_W_mass_aKF);
      histo2D["KF_W_px_orig_vs_corr_TT"]->Fill( (selectedJets[0] + selectedJets[1]).Px(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Px() );
      histo2D["KF_W_py_orig_vs_corr_TT"]->Fill( (selectedJets[0] + selectedJets[1]).Py(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Py() );
      histo2D["KF_W_pz_orig_vs_corr_TT"]->Fill( (selectedJets[0] + selectedJets[1]).Pz(), (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).Pz() );
      
      histo1D["KF_top_mass_corr_TT"]->Fill(reco_top_mass_aKF);
      histo2D["KF_top_mass_orig_vs_corr_TT"]->Fill(reco_top_mass_bKF, reco_top_mass_aKF);
      histo1D["KF_top_pt_corr_TT"]->Fill(reco_top_pt_aKF);
      histo1D["KF_dR_lep_b_corr_TT"]->Fill(reco_dRLepB_lep_aKF);
      
      histo1D["KF_jet0_Et_diff_TT"]->Fill((selectedJetsKFcorrected[0] - selectedJets[labelsReco[0]]).Et());
      histo1D["KF_jet1_Et_diff_TT"]->Fill((selectedJetsKFcorrected[1] - selectedJets[labelsReco[1]]).Et());
      histo2D["KF_jets_Et_diff_TT"]->Fill((selectedJetsKFcorrected[0] - selectedJets[labelsReco[0]]).Et(), (selectedJetsKFcorrected[1] - selectedJets[labelsReco[1]]).Et());
    }
    
    if ( kFitChi2 < 5. )
    {
      if (! doneKinFit) histo1D["KF_top_mass_orig_ex4j_chi2cut5_TT"]->Fill(reco_top_mass_bKF);
      else histo1D["KF_top_mass_corr_ex4j_chi2cut5_TT"]->Fill(reco_top_mass_aKF);
      
      if ( kFitChi2 < 2. )
      {
        if (! doneKinFit) histo1D["KF_top_mass_orig_ex4j_chi2cut2_TT"]->Fill(reco_top_mass_bKF);
        else histo1D["KF_top_mass_corr_ex4j_chi2cut2_TT"]->Fill(reco_top_mass_aKF);
        
        if ( kFitChi2 < 1.5 )
        {
          if (! doneKinFit) histo1D["KF_top_mass_orig_ex4j_chi2cut1p5_TT"]->Fill(reco_top_mass_bKF);
          else histo1D["KF_top_mass_corr_ex4j_chi2cut1p5_TT"]->Fill(reco_top_mass_aKF);
        }  // 1.5
      }  // 2
    }  // 5
  }
}

void FillCatsPlots(string catSuffix)
{
  if (isTTbar)
  {
    histo2D["puSF_vs_nTruePU"]->Fill(puSF,npu);
    histo2D["puSF_vs_nVtx"]->Fill(puSF,nvtx);
    histo2D["nVtx_vs_nTruePU"]->Fill(nvtx,npu);
  }
  
  if (! isData)
  {
    histo1D["red_top_mass_bkgd"]->Fill(redTopMass, widthSF);
    histo1D["red_top_mass_bkgd_manyBins"]->Fill(redTopMass, widthSF);
    if ( isWM || isUM )
    {
      histo1D["red_top_mass_TT_WMUM"]->Fill(redTopMass, widthSF);
      histo1D["red_mlb_mass_TT_WMUM"]->Fill(redMlbMass, widthSF);
    }
    histo1D["red_top_mass_TT"+catSuffix]->Fill(redTopMass, widthSF);
    histo1D["red_mlb_mass_TT"+catSuffix]->Fill(redMlbMass, widthSF);
    histo1D["mass_bjj_div_m_jj"+catSuffix]->Fill(reco_mbjj_div_mjj, widthSF);
    histo1D["mass_bjj_div_m_jj_2"+catSuffix]->Fill(80.385/172.5*reco_mbjj_div_mjj, widthSF);
    histo1D["mlb"+catSuffix]->Fill(reco_mlb_aKF, widthSF);
    histo1D["dR_lep_b_lep"+catSuffix]->Fill(reco_dRLepB_lep_aKF, widthSF);
    histo1D["dR_lep_b_had"+catSuffix]->Fill(reco_dRLepB_had_aKF, widthSF);
    histo1D["ttbar_mass"+catSuffix]->Fill(reco_ttbar_mass_aKF, widthSF);
    histo2D["dR_lep_b_lep_vs_had"+catSuffix]->Fill(reco_dRLepB_lep_aKF, reco_dRLepB_had_aKF, widthSF);
    histo2D["dR_b_light_min_max"+catSuffix]->Fill(reco_dRblight_min_aKF, reco_dRblight_max_aKF, widthSF);
    histo2D["dR_b_light_min_dR_light"+catSuffix]->Fill(reco_dRblight_min_aKF, reco_dRlight_aKF, widthSF);
    histo2D["dR_bW_dR_light"+catSuffix]->Fill(reco_dRbW_aKF, reco_dRlight_aKF, widthSF);
    histo2D["ttbar_mass_vs_minMlb"+catSuffix]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
    histo2D["qqb1_vs_qqb2_aKF"+catSuffix]->Fill(reco_top_mass_aKF, reco_top_mass_alt_aKF, widthSF);
//    if ( reco_dRLepB_lep < 3. )
//    {
//      histo2D[("ttbar_mass_vs_minMlb_dRlepCut"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//      if ( reco_dRLepB_had > 1.2 )
//        histo2D[("ttbar_mass_vs_minMlb_dRBothCuts"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//      
//      if ( reco_dRLepB_lep < 2. )
//      {
//        histo2D[("ttbar_mass_vs_minMlb_dRlepCutHard"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//        if ( reco_dRLepB_had > 2. )
//          histo2D[("ttbar_mass_vs_minMlb_dRBothCutsHard"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//      }
//    }
//    if ( reco_dRLepB_had > 1.2 )
//    {
//      histo2D[("ttbar_mass_vs_minMlb_dRhadCut"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//      if ( reco_dRLepB_had > 2. )
//        histo2D[("ttbar_mass_vs_minMlb_dRhadCutHard"+catSuffix).c_str()]->Fill(reco_mlb_aKF, reco_ttbar_mass_aKF, widthSF);
//    }
    if (doKinFit)
    {
      histo1D["KF_Chi2"+catSuffix]->Fill(kFitChi2);
      histo1D["KF_top_mass_corr"+catSuffix]->Fill(reco_top_mass_aKF, widthSF);
      histo2D["KF_top_mass_orig_vs_corr"+catSuffix]->Fill(reco_top_mass_bKF, reco_top_mass_aKF, widthSF);
      
      histo2D["KF_W_mass_T_orig_vs_corr"+catSuffix]->Fill(reco_mTW_bKF, reco_mTW_aKF, widthSF);
      
      histo2D["KF_chi2_dR_b_light_min"+catSuffix]->Fill(kFitChi2, reco_dRblight_min_aKF, widthSF);
      histo2D["KF_chi2_dR_light"+catSuffix]->Fill(kFitChi2, reco_dRlight_aKF, widthSF);
      histo2D["KF_chi2_dR_bW"+catSuffix]->Fill(kFitChi2, reco_dRbW_aKF, widthSF);
      
      histo2D["KF_chi2_top_mass"+catSuffix]->Fill(kFitChi2, reco_top_mass_aKF, widthSF);
      histo2D["KF_top_mass_chi2"+catSuffix]->Fill(reco_top_mass_aKF, kFitChi2, widthSF);
      histo2D["KF_top_mass_chi2_logY"+catSuffix]->Fill(reco_top_mass_aKF, TMath::Log10(kFitChi2), widthSF);
    }
  }
}

void FillMSPlots(int d, bool doneKinFit)
{
  string suffix = "";
  if (doneKinFit) suffix = "_aKF";
  
  if (makeControlPlots) FillControlPlots(datasetsMSP, d, suffix);
  
  if (! isData)
  {
    MSPlot["scaleFactor"+suffix]->Fill(scaleFactor, datasetsMSP[d], false, lumiWeight*scaleFactor*widthSF);
    MSPlot["btag_SF"+suffix]->Fill(btagSF, datasetsMSP[d], false, lumiWeight*scaleFactor*widthSF);
    MSPlot["pu_SF"+suffix]->Fill(puSF, datasetsMSP[d], false, lumiWeight*scaleFactor*widthSF);
  }
  
  if (! doneKinFit)
  {
    MSPlot["W_mass"]->Fill(reco_W_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["W_pT"]->Fill(reco_W_pt_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["W_mass_T"]->Fill(reco_mTW_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_mass"]->Fill(reco_top_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_mass_400"]->Fill(reco_top_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    if ( reco_top_mass_bKF < 230 && reco_top_mass_bKF > 110 )
    {
      MSPlot["top_mass_zoom"]->Fill(reco_top_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    MSPlot["top_mass_alt"]->Fill(reco_top_mass_alt_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_mass_rescaled"]->Fill(reco_top_mass_rescaled_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    if ( redTopMass_bKF < 2.4 )
    {
      MSPlot["red_top_mass_manyBins"]->Fill(redTopMass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      MSPlot["red_top_mass"]->Fill(redTopMass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      MSPlot["red_top_mass_old"]->Fill(redTopMass_old_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      MSPlot["red_top_mass_new"]->Fill(redTopMass_new_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    MSPlot["top_pT"]->Fill(reco_top_pt_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    
    MSPlot["light_jets_pT"]->Fill(selectedJets[labelsReco[0]].Pt(), datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["light_jets_pT"]->Fill(selectedJets[labelsReco[1]].Pt(), datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["hadr_b_jet_pT"]->Fill(selectedJets[labelsReco[2]].Pt(), datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["lept_b_jet_pT"]->Fill(selectedJets[labelsReco[3]].Pt(), datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    
    MSPlot["mlb"]->Fill(reco_mlb_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["red_mlb_mass"]->Fill(redMlbMass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["ttbar_mass"]->Fill(reco_ttbar_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["mass_bjj_div_m_jj"]->Fill(reco_mbjj_div_mjj_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["mass_bjj_div_m_jj_2"]->Fill(80.385/172.5*reco_mbjj_div_mjj_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["mass_bjj_div_m_lb"]->Fill(reco_top_mass_bKF/reco_mlb_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_lep_b_min"]->Fill(reco_dRLepB_lep_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_lep_b_max"]->Fill(reco_dRLepB_had_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_light_jets"]->Fill(reco_dRlight_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_b_light_min"]->Fill(reco_dRblight_min_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_b_light_max"]->Fill(reco_dRblight_max_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_b_W"]->Fill(reco_dRbW_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_b_light_sum"]->Fill(reco_dRblight_qsum_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_bb"]->Fill(reco_dRbb_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dPhi_light"]->Fill(reco_dPhi_light_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dPhi_bb"]->Fill(reco_dPhi_bb_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dPhi_bW"]->Fill(reco_dPhi_bW_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    
    MSPlot["top_mass_cand_diff"]->Fill(reco_top_mass_bKF-reco_top_mass_alt_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_pt_cand_diff"]->Fill(reco_top_pt_bKF-reco_top_pt_alt_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
  }
  else if (doKinFit)
  {
    if ( kFitChi2 < 5 )
    {
      MSPlot["KF_Chi2"]->Fill(kFitChi2, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      if ( kFitChi2 < 2 )
        MSPlot["KF_Chi2_narrow"]->Fill(kFitChi2, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    
    MSPlot["W_mass"+suffix]->Fill(reco_W_mass_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["W_pT"+suffix]->Fill(reco_W_pt_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["W_mass_T"+suffix]->Fill(reco_mTW_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_mass"+suffix]->Fill(reco_top_mass_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    if ( reco_top_mass_aKF < 230 && reco_top_mass_aKF > 110 )
    {
      MSPlot["top_mass"+suffix+"_zoom"]->Fill(reco_top_mass_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    MSPlot["top_mass_alt"+suffix]->Fill(reco_top_mass_alt_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    if ( redTopMass < 2.4 )
    {
      MSPlot["red_top_mass"+suffix+"_manyBins"]->Fill(redTopMass, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      MSPlot["red_top_mass"+suffix]->Fill(redTopMass, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      MSPlot["red_top_mass_old"+suffix]->Fill(redTopMass_old, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
      MSPlot["red_top_mass_new"+suffix]->Fill(redTopMass_new, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    }
    MSPlot["top_pT"+suffix]->Fill(reco_top_pt_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    
    MSPlot["light_jets_pT"+suffix]->Fill(selectedJetsKFcorrected[0].Pt(), datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["light_jets_pT"+suffix]->Fill(selectedJetsKFcorrected[1].Pt(), datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["hadr_b_jet_pT"+suffix]->Fill(selectedJets[labelsReco[2]].Pt(), datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["lept_b_jet_pT"+suffix]->Fill(selectedJets[labelsReco[3]].Pt(), datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    
    MSPlot["mlb"+suffix]->Fill(reco_mlb_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["red_mlb_mass"+suffix]->Fill(redMlbMass, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["ttbar_mass"+suffix]->Fill(reco_ttbar_mass_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["mass_bjj_div_m_jj"+suffix]->Fill(reco_mbjj_div_mjj, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["mass_bjj_div_m_jj_2"+suffix]->Fill(80.385/172.5*reco_mbjj_div_mjj, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["mass_bjj_div_m_lb"+suffix]->Fill(reco_top_mass_aKF/reco_mlb_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_lep_b_min"+suffix]->Fill(reco_dRLepB_lep_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_lep_b_max"+suffix]->Fill(reco_dRLepB_had_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_light_jets"+suffix]->Fill(reco_dRlight_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_b_light_min"+suffix]->Fill(reco_dRblight_min_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_b_light_max"+suffix]->Fill(reco_dRblight_max_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_b_W"+suffix]->Fill(reco_dRbW_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_b_light_sum"+suffix]->Fill(reco_dRblight_qsum_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dR_bb"+suffix]->Fill(reco_dRbb_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dPhi_light"+suffix]->Fill(reco_dPhi_light_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dPhi_bb"+suffix]->Fill(reco_dPhi_bb_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["dPhi_bW"+suffix]->Fill(reco_dPhi_bW_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    
    MSPlot["W_mass_diff"]->Fill(reco_W_mass_aKF-reco_W_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["W_mass_T_diff"]->Fill(reco_mTW_aKF-reco_mTW_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["W_pT_diff"]->Fill(reco_W_pt_aKF-reco_W_pt_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_mass_diff"]->Fill(reco_top_mass_aKF-reco_top_mass_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_pT_diff"]->Fill(reco_top_pt_aKF-reco_top_pt_bKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    //MSPlot["top_px_diff"]
    //MSPlot["top_py_diff"]
    //MSPlot["top_pz_diff"]
    
    MSPlot["top_mass_cand_diff"+suffix]->Fill(reco_top_mass_aKF-reco_top_mass_alt_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
    MSPlot["top_pt_cand_diff"+suffix]->Fill(reco_top_pt_aKF-reco_top_pt_alt_aKF, datasetsMSP[d], true, lumiWeight*scaleFactor*widthSF);
  }
}

void FillLikelihoodPlots()
{
  for (int iMass = 0; iMass < nLikeMasses; iMass++)
  {
    loglike_per_evt.clear();
    if ( redTopMassArray[iMass] > minCutRedTopMass && maxCutRedTopMass > redTopMassArray[iMass])
    {
      if (doLikeW)       loglike_per_evt = likeW->CalculateLikelihood(redTopMassArray[iMass], 1., 172.5, 172.5, 1., 172.5, false, false);
      else if (doLikeM)  loglike_per_evt = likeM->CalculateLikelihood(redTopMassArray[iMass], 1., 172.5, 172.5, 1., 172.5, false, false);
      else if (doLike2D) loglike_per_evt = like2D->CalculateLikelihood(redTopMassArray[iMass], 1., 172.5, 172.5, 1., 172.5, false, false);
      for (int iWidth = 0; iWidth < nWidthsLike; iWidth++)  // Deliberately filling iWidth and not actual width
      {
//        if (redTopMass > redTopMassArray[iMass] && redTopMass < redTopMassArray[iMass+1])
//        {
          histo1DLike["loglike_vs_width_m"+DotReplace(redTopMassArray[iMass])]->SetBinContent(iWidth, loglike_per_evt[iWidth]);
//          if (isCM) histo1DLike["loglike_vs_width_CM_m"+DotReplace(redTopMassArray[iMass])]->SetBinContent(iWidth, loglike_per_evt[iWidth]);
//          if (isWM) histo1DLike["loglike_vs_width_WM_m"+DotReplace(redTopMassArray[iMass])]->SetBinContent(iWidth, loglike_per_evt[iWidth]);
//          if (isUM) histo1DLike["loglike_vs_width_UM_m"+DotReplace(redTopMassArray[iMass])]->SetBinContent(iWidth, loglike_per_evt[iWidth]);
//        }
      }  // end width

//      break;  // only fill once per event (because only one redTopMass)
    }
  }  // end masses
}

void WriteLikelihoodPlots()
{
  FillLikelihoodPlots();
  
  mkdir((pathOutput+"Likelihood/").c_str(),0777);
  string fileName = pathOutput+"Likelihood/LikelihoodPlots.root";
  cout << " - Recreate output file ..." << endl;
  TFile *foutLike = new TFile (fileName.c_str(), "RECREATE");
  cout << "   Output file is " << fileName << endl;
  foutLike->cd();
  
  //gStyle->SetOptStat(1111);
  for (std::map<std::string,TH1D*>::const_iterator it = histo1DLike.begin(); it != histo1DLike.end(); it++)
  {
    TH1D *temp = it->second;
    for (int iBin = 0; iBin < temp->GetNbinsX(); iBin++)
    {
      temp->GetXaxis()->SetBinLabel(iBin+1, (DotReplace(widthsLike[iBin])).c_str());
    }
    temp->SetStats(0);
    temp->Write();
//     TCanvas* tempCanvas = new TCanvas((it->first).c_str(), (it->first).c_str());
//     tempCanvas->cd();
//     temp->Draw();
//     tempCanvas->SaveAs( (pathOutput+"Likelihood/"+it->first+".png").c_str() );
  }
  
  foutLike->Close();
  delete foutLike;
}

long GetNEvents(TTree* fChain, string var, bool isData)
{
  GetMetaData(fChain, isData);
  long varNew = 0;
  for (unsigned int iEntry = 0; iEntry < fChain->GetEntries(); iEntry++)
  {
    fChain->GetEntry(iEntry);
    varNew += (fChain->FindLeaf(var.c_str()))->GetValueLong64();
  }
  
  return varNew;
}

long GetNEvents(TTree* fChain, string var, unsigned int index, bool isData)
{
  GetMetaData(fChain, isData);
  long varNew = 0;
  for (unsigned int iEntry = 0; iEntry < fChain->GetEntries(); iEntry++)
  {
    fChain->GetEntry(iEntry);
    varNew += (fChain->FindLeaf(var.c_str()))->GetValueLong64(index);
  }
  
  return varNew;
}

void GetEraFraction(double* fractions)
{
  TFile* file = new TFile((pathNtuplesData+"Ntuples_data.root").c_str(),"READ");
  TTree* fChain = (TTree*) file->Get("stats");
  
  //fChain->SetBranchAddress("nofEventsHLTv2", &nofEventsHLTv2, &b_nofEventsHLTv2);
  //fChain->SetBranchAddress("nofEventsHLTv3", &nofEventsHLTv3, &b_nofEventsHLTv3);
  
  long nBCDEF = GetNEvents(fChain, "nofEventsRunB", 1) + GetNEvents(fChain, "nofEventsRunCD", 1) + GetNEvents(fChain, "nofEventsRunEF", 1);
  long nGH = GetNEvents(fChain, "nofEventsRunG", 1) + GetNEvents(fChain, "nofEventsRunH", 1);
  
  fractions[0] = ((double)nBCDEF) / ((double)nBCDEF + (double)nGH);
  fractions[1] = ((double)nGH) / ((double)nBCDEF + (double)nGH);
  
  file->Close();
  ClearMetaData();
}

void SetUpSelectionTable()
{
  selTab = new SelectionTables(datasets);
  selTab->SetLumi(Luminosity);
  selTab->AddCutStep("preselected");
  selTab->AddCutStep("triggered");
  selTab->AddCutStep("PV \\& filters");
  selTab->AddCutStep("1 muon");
  selTab->AddCutStep("veto muon");
  selTab->AddCutStep("veto elec");
  selTab->AddCutStep("$\\geq 4$ jets");
  
//   // Cutflow 1
//   selTab->AddCutStep("$\\geq 1$ \\bq~jet");
//   selTab->AddCutStep("$\\geq 2$ \\bq~jets");
//   selTab->AddCutStep("$= 4$ jets");
  
  // Cutflow 2
  selTab->AddCutStep("$= 4$ jets");
  selTab->AddCutStep("$\\geq 1$ \\bq~jet");
  selTab->AddCutStep("$\\geq 2$ \\bq~jets");
  
  selTab->AddCutStep("$\\forall j,j': \\Delta R(j,j') > 0.6$");
  selTab->AddCutStep("KF $\\chi^{2} < 5$");
  selTab->SetUpTable();
  
  selTabKF = new SelectionTables(datasetsMSP);
  selTabKF->SetLumi(Luminosity);
  selTabKF->AddCutStep("Before KF");
  selTabKF->AddCutStep("KF $\\chi^2 < 5$");
  selTabKF->SetUpTable();
}

void FillSelectionTable(int d, string dataSetName)
{
  vector<double> cutFlowValues;
  cutFlowValues.clear();
  for (unsigned int i = 0; i < 10; i++)
  {
//     // Cutflow 1
//     cutFlowValues.push_back(GetNEvents(tStatsTree[(dataSetName).c_str()], "cutFlowWeighted", i, isData));
    // Cutflow 2
    cutFlowValues.push_back(GetNEvents(tStatsTree[(dataSetName).c_str()], "cutFlow2Weighted", i, isData));
  }
  
  selTab->Fill(d, cutFlowValues);
  selTab->Fill(d, 10, nofAfterDRmincut);
  selTab->Fill(d, 11, nofAcceptedKFitWeighted);
}

void WriteSelectionTable()
{
  selTab->CalculateTable();
  selTab->Write("SelectionTable_semiMu_notMerged.tex", true, false, true);
  selTab->Write("SelectionTable_semiMu.tex", true, true, true);
  
  for (unsigned int d = 0; d < dMSPmax; d++)
  {
    selTabKF->Fill(d, 0, nofBKF_weighted[d]);
    selTabKF->Fill(d, 1, nofAKF_weighted[d]);
  }
  selTabKF->CalculateTable(false);
  selTabKF->Write("SelectionTableKF_semiMu.tex", false, true, false, true);
}

void CheckSystematics(vector<int> vJER, vector<int> vJES, vector<int> vPU)
{
  int sumJER = 0, sumJES = 0, sumPU = 0;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    sumJER += vJER[d];
    sumJES += vJES[d];
    sumPU  += vPU[d];
  }
  
  if ( sumJER == 0 && sumJES == 0 && sumPU == 0 )
  {
    strSyst = "nominal";
  }
  else if ( abs(sumJER) == datasets.size() && sumJES == 0 && sumPU == 0 )
  {
    if ( sumJER > 0 ) { strSyst = "JERup";}
    else if ( sumJER < 0 ) { strSyst = "JERdown";}
  }
  else if ( sumJER == 0 && abs(sumJES) == datasets.size() && sumPU == 0 )
  {
    if ( sumJES > 0 ) { strSyst = "JESup";}
    else if ( sumJES < 0 ) { strSyst = "JESdown";}
  }
  else if ( sumJER == 0 && sumJES == 0 && abs(sumPU) == datasets.size() )
  {
    if ( sumPU > 0 ) { strSyst = "PUup";}
    else if ( sumPU < 0 ) { strSyst = "PUdown";}
  }
  else
  {
    cerr << "Shape changing systematics not consistent accross datasets or multiple applied at once" << endl;
    cerr << "Exiting...." << endl;
    exit(1);
  }
  cout << " - Systematics confirmed to be " << strSyst << endl;
}

void PrintKFDebug(int ievt)
{
  cout << endl <<"Event " << setw(5) << right << ievt << "   ";
  cout << "Top mass after kinFit is negative! I.e. " << reco_top_mass_aKF << " Before kinFit: " << reco_top_mass_bKF << "  ievt " << ievt << endl;
  cout << "Mass jet 1 & jet 2: " << (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[1]).M() << "; Mass jet 1 & jet 3: " << (selectedJetsKFcorrected[0] + selectedJetsKFcorrected[2]).M() << "; Mass jet 2 & jet 3: " << (selectedJetsKFcorrected[1] + selectedJetsKFcorrected[2]).M() << endl;
  //if (test)
  //{
  cout << "Original:   Jet 1: pT " << selectedJets[labelsReco[0]].Pt() << "; Jet 2: pT " << selectedJets[labelsReco[1]].Pt() << "; Jet 3: pT " << selectedJets[labelsReco[2]].Pt() << endl;
  cout << "Corrected:  Jet 1: pT " << selectedJetsKFcorrected[0].Pt() << "; Jet 2: pT " << selectedJetsKFcorrected[1].Pt() << "; Jet 3: pT " << selectedJetsKFcorrected[2].Pt() << endl;
  cout << "Original:   Jet 1: px " << selectedJets[labelsReco[0]].Px() << "; py " << selectedJets[labelsReco[0]].Py() << "; pz " << selectedJets[labelsReco[0]].Pz() << "; E " << selectedJets[labelsReco[0]].E() << endl;
  cout << "Corrected:  Jet 1: px " << selectedJetsKFcorrected[0].Px() << "; py " << selectedJetsKFcorrected[0].Py() << "; pz " << selectedJetsKFcorrected[0].Pz() << "; E " << selectedJetsKFcorrected[0].E() << endl;
  cout << "Original:   Jet 2: px " << selectedJets[labelsReco[1]].Px() << "; py " << selectedJets[labelsReco[1]].Py() << "; pz " << selectedJets[labelsReco[1]].Pz() << "; E " << selectedJets[labelsReco[1]].E() << endl;
  cout << "Corrected:  Jet 2: px " << selectedJetsKFcorrected[1].Px() << "; py " << selectedJetsKFcorrected[1].Py() << "; pz " << selectedJetsKFcorrected[1].Pz() << "; E " << selectedJetsKFcorrected[1].E() << endl;
  cout << "Original:   Jet 1: pt " << selectedJets[labelsReco[0]].Pt() << "; eta " << selectedJets[labelsReco[0]].Eta() << "; phi " << selectedJets[labelsReco[0]].Phi() << "; M " << selectedJets[labelsReco[0]].M() << endl;
  cout << "Corrected:  Jet 1: pt " << selectedJetsKFcorrected[0].Pt() << "; eta " << selectedJetsKFcorrected[0].Eta() << "; phi " << selectedJetsKFcorrected[0].Phi() << "; M " << selectedJetsKFcorrected[0].M() << endl;
  cout << "Original:   Jet 2: pt " << selectedJets[labelsReco[1]].Pt() << "; eta " << selectedJets[labelsReco[1]].Eta() << "; phi " << selectedJets[labelsReco[1]].Phi() << "; M " << selectedJets[labelsReco[1]].M() << endl;
  cout << "Corrected:  Jet 2: pt " << selectedJetsKFcorrected[1].Pt() << "; eta " << selectedJetsKFcorrected[1].Eta() << "; phi " << selectedJetsKFcorrected[1].Phi() << "; M " << selectedJetsKFcorrected[1].M() << endl;
  cout << "Original:   Jet 1: x " << selectedJets[labelsReco[0]].X() << "; y " << selectedJets[labelsReco[0]].Y() << "; z " << selectedJets[labelsReco[0]].Z() << "; t " << selectedJets[labelsReco[0]].T() << endl;
  cout << "Corrected:  Jet 1: x " << selectedJetsKFcorrected[0].X() << "; y " << selectedJetsKFcorrected[0].Y() << "; z " << selectedJetsKFcorrected[0].Z() << "; t " << selectedJetsKFcorrected[0].T() << endl;
  cout << "Original:   Jet 2: x " << selectedJets[labelsReco[1]].X() << "; y " << selectedJets[labelsReco[1]].Y() << "; z " << selectedJets[labelsReco[1]].Z() << "; t " << selectedJets[labelsReco[1]].T() << endl;
  cout << "Corrected:  Jet 2: x " << selectedJetsKFcorrected[1].X() << "; y " << selectedJetsKFcorrected[1].Y() << "; z " << selectedJetsKFcorrected[1].Z() << "; t " << selectedJetsKFcorrected[1].T() << endl;
  //}
}

void InitJEC(bool isData, string dataSetName)
{
  if (! isData)
  {
    JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt");
    vCorrParam.push_back(*L1JetCorPar);
    JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt");
    vCorrParam.push_back(*L2JetCorPar);
    JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt");
    vCorrParam.push_back(*L3JetCorPar);
    jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
    //jecUnc = new JetCorrectionUncertainty("input/JEC_Uncertainty_temp.txt");
    /// Do not use total jes uncertainties
    //JetCorrectorParameters *MCUncCorPar = new JetCorrectorParameters(pathCalJEC+"/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt","TotalNoFlavor");
    //jecUnc = new JetCorrectionUncertainty(*MCUncCorPar);
  }
  else if ( dataSetName.find("Run2016B") != std::string::npos || dataSetName.find("Run2016C") != std::string::npos || dataSetName.find("Run2016D") != std::string::npos )
  {
    JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt");
    vCorrParam.push_back(*L1JetCorPar);
    JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt");
    vCorrParam.push_back(*L2JetCorPar);
    JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt");
    vCorrParam.push_back(*L3JetCorPar);
    JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt");
    vCorrParam.push_back(*L2L3ResJetCorPar);
    jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt");
  }
  else if ( dataSetName.find("Run2016E") != std::string::npos || dataSetName.find("Run2016F") != std::string::npos )
  {
    JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt");
    vCorrParam.push_back(*L1JetCorPar);
    JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt");
    vCorrParam.push_back(*L2JetCorPar);
    JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt");
    vCorrParam.push_back(*L3JetCorPar);
    JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt");
    vCorrParam.push_back(*L2L3ResJetCorPar);
    jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFchs.txt");
  }
  else if ( dataSetName.find("Run2016G") != std::string::npos )
  {
    JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt");
    vCorrParam.push_back(*L1JetCorPar);
    JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt");
    vCorrParam.push_back(*L2JetCorPar);
    JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt");
    vCorrParam.push_back(*L3JetCorPar);
    JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt");
    vCorrParam.push_back(*L2L3ResJetCorPar);
    jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFchs.txt");
  }
  else if ( dataSetName.find("Run2016H") != std::string::npos )
  {
    JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt");
    vCorrParam.push_back(*L1JetCorPar);
    JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt");
    vCorrParam.push_back(*L2JetCorPar);
    JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt");
    vCorrParam.push_back(*L3JetCorPar);
    JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt");
    vCorrParam.push_back(*L2L3ResJetCorPar);
    jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt");
  }
}

