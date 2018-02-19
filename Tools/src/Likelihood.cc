#include "../interface/Likelihood.h"

const double Likelihood::widthArray_[] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5., 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6., 6.1, 6.2, 6.3, 6.4, 6.5};

//const std::string Likelihood::listCats_[] = {"CM", "WM", "UM"};
const std::string Likelihood::listCats_[] = {"CM", "UM"};

const int Likelihood::nWidths_ = sizeof(widthArray_)/sizeof(widthArray_[0]);
const int Likelihood::nCats_ = sizeof(listCats_)/sizeof(listCats_[0]);

std::string Likelihood::stringWidthArray_[nWidths_] = {""};
std::string Likelihood::stringSuffix_[nCats_] = {""};

double Likelihood::loglike_[nWidths_] = {0.};
double Likelihood::loglike_data_[nWidths_] = {0.};
double Likelihood::loglike_per_evt_[nWidths_] = {0.};
double Likelihood::loglike_CM_[nWidths_] = {0.};
double Likelihood::loglike_CM_per_evt_[nWidths_] = {0.};
double Likelihood::loglike_temp_[nWidths_] = {0.};
double Likelihood::loglike_temp_per_evt_[nWidths_] = {0.};
double Likelihood::loglike_gen_[nWidths_] = {0.};
double Likelihood::loglike_gen_per_evt_[nWidths_] = {0.};

double Likelihood::loglike_pull_[nWidths_][1000] = {{0.}};
double Likelihood::loglike_pull_single_[nWidths_] = {0.};

//const double Likelihood::calCurvePar_[2] = {0., 1.};  // at the moment no output calibration
//const double Likelihood::calCurveParUnc_[2] = {0., 0.};  // at the moment no output calibration
//const double Likelihood::calCurvePar_[2] = {-0.152752, 1.09383};  // HAD
//const double Likelihood::calCurveParUnc_[2] = {0.0436359, 0.0166395};  // HAD
//const double Likelihood::calCurvePar_[2] = {-0.267184, 1.08957};  // LEP
//const double Likelihood::calCurveParUnc_[2] = {0.0550971, 0.0248913};  // LEP
const double Likelihood::calCurvePar_[2] = {-0.249918, 1.08351};  // COMB
const double Likelihood::calCurveParUnc_[2] = {0.0351629, 0.0139145};  // COMB

double Likelihood::nEventsCMFractions_[nWidths_][25] = {{0.}};
double Likelihood::nEventsWMFractions_[nWidths_][25] = {{0.}};
double Likelihood::nEventsUMFractions_[nWidths_][25] = {{0.}};


int Likelihood::LocMinArray(int n, double* array)
{
  if ( n == 0 ) return -1;
  
  int locmin = 0;
  double min = array[0];
  for (int i = 0; i < n; i++)
  {
    if ( array[i] < min )
    {
      min = array[i];
      locmin = i;
    }
  }
  return locmin;
}

void Likelihood::ClearArray(int size, int* array)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = 0;
  }
}

void Likelihood::ClearArray(int size, double* array)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = 0.;
  }
}

void Likelihood::ClearArray2D(int size, double (*array)[3])
{
  for (int i = 0; i < size; i++)
  {
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      array[i][iCat] = 0.;
    }
  }
}

void Likelihood::MakeTable(double* array, int n, double min, double max)
{
  double dist = (max - min)/((double)(n-1));
  for (int i = 0; i < n; i++)
  {
    array[i] = min + i * dist;
  }
}

Double_t Likelihood::fitFunc(Double_t *x, Double_t *par)
{
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

Likelihood::Likelihood(double min, double max, std::string outputDirName, std::string date, bool useHadTopOnly, bool useNewVar, bool makeHistograms, bool verbose):
verbose_(verbose), rewHadOnly_(useHadTopOnly), useHadVar_(true), outputDirName_(outputDirName), dirNameTGraphTxt_("OutputTxt/"), dirNameNEvents_("OutputNEvents/"), dirNameLLTxt_("OutputLikelihood/"+date+"/"), dirNamePull_("PseudoExp/"), inputFileName_(""), suffix_(""), histoName_(""), minRedMass_(min), maxRedMass_(max), histo_(), histoSm_(), histoTotal_(), graph_(), vecBinCentres_(), vecBinContents_(), calledLLCalculation_(false), calledCMLLCalculation_(false), calledGenLLCalculation_(false), vecWidthFromFile_(), vecLLValsFromFile_()
{
  tls_ = new HelperTools();
  rew_ = new EventReweighting(false);  // no correction for number of events
  
  if (useNewVar) useHadVar_ = false;
  
  rangeRedMass_ = tls_->DotReplace(minRedMass_)+"To"+tls_->DotReplace(maxRedMass_);
  dirNameNEvents_ += rangeRedMass_+"/";
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    stringWidthArray_[iWidth] = tls_->DotReplace(widthArray_[iWidth]);
  }
  
  this->ClearLikelihoods();
  
  if (makeHistograms)
  {
    mkdir(outputDirName_.c_str(),0777);
    mkdir((outputDirName_+dirNameTGraphTxt_).c_str(),0777);
    this->BookHistograms();
  }
  else mkdir(dirNameLLTxt_.c_str(),0777);
}

Likelihood::Likelihood(double minhad, double maxhad, double minlep, double maxlep, std::string outputDirName, std::string outputDirName2, std::string date, bool useHadTopOnly, bool makeHistograms, bool verbose):
verbose_(verbose), rewHadOnly_(useHadTopOnly), useHadVar_(true), outputDirName_(outputDirName), outputDirName2_(outputDirName2), dirNameTGraphTxt_("OutputTxt/"), dirNameNEvents_("OutputNEvents/"), dirNameLLTxt_("OutputLikelihood/"+date+"/"), dirNamePull_("PseudoExp/"), inputFileName_(""), suffix_(""), histoName_(""), minRedMass_(minhad), maxRedMass_(maxhad), minRedMass2_(minlep), maxRedMass2_(maxlep), histo_(), histoSm_(), histoTotal_(), graph_(), vecBinCentres_(), vecBinContents_(), calledLLCalculation_(false), calledCMLLCalculation_(false), calledGenLLCalculation_(false), vecWidthFromFile_(), vecLLValsFromFile_()
{
  if (makeHistograms)
  {
    std::cerr << "Likelihood::Likelihood: Cannot make histograms for both variables simulataneously... Exiting..." << std::endl;
    exit(1);
  }
  
  tls_ = new HelperTools();
  rew_ = new EventReweighting(false);  // no correction for number of events
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    stringWidthArray_[iWidth] = tls_->DotReplace(widthArray_[iWidth]);
  }
  
  this->ClearLikelihoods();
  
  bool tmpbool;
  tmpbool = this->ConstructTGraphsFromFile(outputDirName_, "had_");
  tmpbool = this->ConstructTGraphsFromFile(outputDirName2_, "lep_");
  if (! tmpbool)
  {
    std::cerr << "Likelihood::Likelihood: Something went wrong when loading graphs... Exiting..." << std::endl;
    exit(1);
  }
  
  mkdir(dirNameLLTxt_.c_str(),0777);
}

Likelihood::~Likelihood()
{
  delete file_;
  delete rew_;
  delete tls_;
  
  if ( fileTGraphs_ != NULL )
  {
    fileTGraphs_->Close();
    delete fileTGraphs_;
  }
}

void Likelihood::ClearLikelihoods()
{
  for (int i = 0; i < nWidths_; i++)
  {
    loglike_[i] = 0.;
    loglike_data_[i] = 0.;
    loglike_per_evt_[i] = 0.;
    loglike_CM_[i] = 0.;
    loglike_temp_[i] = 0.;
    loglike_gen_[i] = 0.;
  }
}

std::vector<double> Likelihood::GetWidths()
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    vecWidths_.push_back(widthArray_[iWidth]);
  }
  return vecWidths_;
}

void Likelihood::BookHistograms()
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    thisWidth_ = stringWidthArray_[iWidth];
    
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      if (useHadVar_) histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", "+listCats_[iCat]+"; m_{r}").c_str(), 90, 0.5, 2.0);  // 90 bins in range [0.5, 2.0]
      /// red mlb
      else histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", "+listCats_[iCat]+"; m_{lb,r}").c_str(), 80, 0., 2.5);  // 75 bins
//       histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_79b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_79b").c_str(),("Reduced top mass for width "+thisWidth_+", "+listCats_[iCat]+"; m_{lb,r}").c_str(), 79, 0., 2.5);
//       histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_82b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_82b").c_str(),("Reduced top mass for width "+thisWidth_+", "+listCats_[iCat]+"; m_{lb,r}").c_str(), 82, 0., 2.5);
//       histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_86b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_86b").c_str(),("Reduced top mass for width "+thisWidth_+", "+listCats_[iCat]+"; m_{lb,r}").c_str(), 86, 0., 2.5);
       /// dR(lep,b)
       //histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", "+listCats_[iCat]+"; dR(l,b_{l})").c_str(), 80, 0.4, 4.4);  // 75,90 bins
      /// New var
      //histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", "+listCats_[iCat]+"; m_{3/2}").c_str(), 50, 1., 3.5);
    }
  }
}

void Likelihood::FillHistograms(double redMass, double relativeSF, double hadTopMassForWidthSF, double lepTopMassForWidthSF, bool isTTbar, bool isData, std::string catSuffix)
{
  if ( ! isData )
  {
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      thisWidth_ = stringWidthArray_[iWidth];
      if (isTTbar)
      {
        if (rewHadOnly_) thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., widthArray_[iWidth]);
        else thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., widthArray_[iWidth]) * rew_->BEventWeightCalculatorNonRel(lepTopMassForWidthSF, 172.5, 172.5, 1., widthArray_[iWidth]);
      }
      else thisWidthSF_ = 1.;
      
      histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_90b").c_str()]->Fill(redMass, relativeSF*thisWidthSF_);
      // mlb
//      histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_79b").c_str()]->Fill(redMass, relativeSF*thisWidthSF_);
//      histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_82b").c_str()]->Fill(redMass, relativeSF*thisWidthSF_);
//      histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_86b").c_str()]->Fill(redMass, relativeSF*thisWidthSF_);
    }
  }
}

void Likelihood::WriteHistograms(std::string histoFileName)
{
  file_ = new TFile((outputDirName_+histoFileName).c_str(), "RECREATE");
  file_->cd();
  gStyle->SetOptStat(1111);
  for (std::map<std::string,TH1D*>::const_iterator it = histo_.begin(); it != histo_.end(); it++)
  {
    TH1D *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
  }
  file_->Close();
  delete file_;
}

void Likelihood::GetHistogram(int iCat)
{
  int binMin = 0;
  int binMax = 1000;

  /// Get histo to smooth
  histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat];
  histoSm_[histoName_] = (TH1D*) histo_["Red_top_mass_"+histoName_+"_90b"]->Clone(histoName_.c_str());
  
  if (! useHadVar_)
  {
    int nBins = histoSm_[histoName_]->GetNbinsX();
    // Remove high spikes in flanks
    for (int i = 2; i < nBins-1; i++)
    {
      if ( histoSm_[histoName_]->GetBinCenter(i) < 0.9 && histoSm_[histoName_]->GetBinContent(i) > histoSm_[histoName_]->GetBinContent(i+1) && histoSm_[histoName_]->GetBinContent(i) > histoSm_[histoName_]->GetBinContent(i+2) )
        histoSm_[histoName_]->SetBinContent(i, (histoSm_[histoName_]->GetBinContent(i-1)+histoSm_[histoName_]->GetBinContent(i+1))/2. );
      //else if ( histoSm_[histoName_]->GetBinCenter(i) > 1.2 && histoSm_[histoName_]->GetBinContent(i) > histoSm_[histoName_]->GetBinContent(i-1) && histoSm_[histoName_]->GetBinContent(i) > histoSm_[histoName_]->GetBinContent(i-2) )
      else if ( histoSm_[histoName_]->GetBinCenter(i) > 1.5 && histoSm_[histoName_]->GetBinContent(i) > histoSm_[histoName_]->GetBinContent(i-1) && histoSm_[histoName_]->GetBinContent(i) > histoSm_[histoName_]->GetBinContent(i-2) )
        histoSm_[histoName_]->SetBinContent(i, (histoSm_[histoName_]->GetBinContent(i-1)+histoSm_[histoName_]->GetBinContent(i+1))/2. );
    }
    // Remove small dips due to stats
//    for (int i = 2; i < nBins-1; i++)
//    {
//      if ( histoSm_[histoName_]->GetBinContent(i) < histoSm_[histoName_]->GetBinContent(i+1) && histoSm_[histoName_]->GetBinContent(i) < histoSm_[histoName_]->GetBinContent(i-1) )
//        histoSm_[histoName_]->SetBinContent(i, (histoSm_[histoName_]->GetBinContent(i-1)+histoSm_[histoName_]->GetBinContent(i+1))/2. );
//    }
//    for (int i = 2; i < nBins-1; i++) // second iteration
//    {
//      if ( histoSm_[histoName_]->GetBinContent(i) < histoSm_[histoName_]->GetBinContent(i+1) && histoSm_[histoName_]->GetBinContent(i) < histoSm_[histoName_]->GetBinContent(i-1) )
//        histoSm_[histoName_]->SetBinContent(i, (histoSm_[histoName_]->GetBinContent(i-1)+histoSm_[histoName_]->GetBinContent(i+1))/2. );
//    }
  }
  
  if ( iCat != 0 ) histoSm_[histoName_]->Smooth(3);
  else if (! useHadVar_) histoSm_[histoName_]->Smooth(1);
  histoSm_[histoName_]->Write();

  //nBins[iCat] = histoSm_[histoName_]->GetNbinsX();
  binMin = histoSm_[histoName_]->FindBin(minRedMass_);
  binMax = histoSm_[histoName_]->FindBin(maxRedMass_)+1;
  
  /// Normalise histo on relevant subdomain
  Double_t integral = histoSm_[histoName_]->Integral(binMin, binMax);
  histoSm_[histoName_]->Scale(1./integral);

  /// Get bin centres & contents
  for (int iBin = binMin; iBin < binMax+1; iBin++)
  {
    (vecBinCentres_[histoName_]).push_back(histoSm_[histoName_]->GetBinCenter(iBin));
    (vecBinContents_[histoName_]).push_back(histoSm_[histoName_]->GetBinContent(iBin));
  }
}

void Likelihood::MakeGraph(int iCat, int nPoints, double* centres, double* contents, std::string name, bool drawGraph)
{
  histoName_ = name+stringSuffix_[iCat];
  graph_[histoName_] = new TGraph(nPoints, centres, contents);
  graph_[histoName_]->SetName(("g"+histoName_).c_str());
  graph_[histoName_]->SetTitle(histoName_.c_str());
  graph_[histoName_]->Write();
  if (drawGraph) this->DrawGraph(histoSm_[histoName_], graph_[histoName_], "Graph_Red_top_mass_"+histoName_);
}

void Likelihood::MakeGraphSmooth(int iCat, int nPoints, double* centres, double* contents, std::string name, bool drawGraph)
{
  this->MakeGraph(iCat, nPoints, centres, contents, name, drawGraph);
  
  histoName_ = name+stringSuffix_[iCat];
  histoNameSm_ = listCats_[iCat]+"_Sm_"+stringSuffix_[iCat];
  TGraphSmooth* gs = new TGraphSmooth(histoName_.c_str());
  //graph_[histoNameSm_] = gs->SmoothSuper(graph_[histoName_],"",0,0);
  graph_[histoNameSm_] = gs->SmoothSuper(graph_[histoName_],"",3);
  graph_[histoNameSm_]->SetName(("g"+histoNameSm_).c_str());
  graph_[histoNameSm_]->SetTitle(histoNameSm_.c_str());
  graph_[histoNameSm_]->Write();
  if (drawGraph) this->DrawGraph(histoSm_[histoName_], graph_[histoNameSm_], "Graph_Red_top_mass_"+histoNameSm_);
}

void Likelihood::ConstructTGraphsFromHisto(std::string tGraphFileName, std::vector<std::string> datasetNames, std::vector<int> includeDataset)
{
  TGraph2D *gLL2D_ = new TGraph2D();
  
  /// Define vars for likelihood calculation
  const int nEval = 50;
  double evalPoints[nEval], outputValues[nEval], outputValuesTemp[nEval], likelihoodValues[nEval], likelihoodValuesTemp[nEval];
  this->ClearArray(nEval, evalPoints);
  this->MakeTable(evalPoints, nEval, minRedMass_, maxRedMass_);
  
  /// Determine fractions based on number of events within reduced top mass range
  //  (is not dependent on number of bins, so the same for all widths)
  double fracCats[nCats_] = {0.}, fracCatsTemp[nCats_-1] = {0.};
  this->GetFractions(fracCats, nCats_, datasetNames, includeDataset);
  this->GetFractions(fracCatsTemp, nCats_-1, datasetNames, includeDataset);
  if (verbose_)
  {
    for (int iCat = 0; iCat < nCats_; iCat++)
      std::cout << "   # " << listCats_[iCat] << ": " << fracCats[iCat] << "*100%";
    std::cout << std::endl;
  }
  
  /// Make output file
  fileTGraphs_ = new TFile((outputDirName_+tGraphFileName).c_str(), "RECREATE");
  fileTGraphs_->cd();
  
  /// WM & UM distribution are independent of the width
  for (int iCat = 1; iCat < nCats_; iCat++)
  {
    stringSuffix_[iCat] = "widthx1";
    this->GetHistogram(iCat);
  }
  
  /// Make arrays as input for TGraph
  const int nPoints = (vecBinCentres_[listCats_[1]+"_"+stringSuffix_[1]]).size();
  double binCentreArray[nPoints], binContentArray[nCats_][nPoints];
  this->ClearArray(nPoints, binCentreArray);
  for (int iCat = 0; iCat < nCats_; iCat++) this->ClearArray(nPoints, binContentArray[iCat]);
  
  for (int i = 0; i < nPoints; i++)
  {
    binCentreArray[i] = (vecBinCentres_[listCats_[1]+"_"+stringSuffix_[1]]).at(i);
    for (int iCat = 1; iCat < nCats_; iCat++)
    {
      binContentArray[iCat][i] = (vecBinContents_[listCats_[iCat]+"_"+stringSuffix_[iCat]]).at(i);
    }
  }
  
  for (int iCat = 1; iCat < nCats_; iCat++)
  {
    this->WriteFuncOutput(nPoints, binCentreArray, binContentArray[iCat], listCats_[iCat]+"_"+stringSuffix_[iCat]);
    this->MakeGraphSmooth(iCat, nPoints, binCentreArray, binContentArray[iCat], listCats_[iCat]+"_", true);
    histoSm_[histoName_]->Scale(fracCats[iCat]);
  }
  
  
  /// Loop over widths
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    stringSuffix_[0] = "widthx"+stringWidthArray_[iWidth];
    this->GetHistogram(0);
    
    this->ClearArray(nPoints, binContentArray[0]);
    for (int i = 0; i < nPoints; i++)
      binContentArray[0][i] = (vecBinContents_[listCats_[0]+"_"+stringSuffix_[0]]).at(i);
    
    this->WriteFuncOutput(nPoints, binCentreArray, binContentArray[0], listCats_[0]+"_"+stringSuffix_[0]);
    this->MakeGraphSmooth(0, nPoints, binCentreArray, binContentArray[0], listCats_[0]+"_", true);
    
    /// Make likelihood functions
    this->ClearArray(nEval, outputValues);
    this->ClearArray(nEval, outputValuesTemp);
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat];
      histoNameSm_ = listCats_[iCat]+"_Sm_"+stringSuffix_[iCat];
      
      if ( iCat == 0 ) histoSm_[histoName_]->Scale(fracCats[iCat]);
      
      this->ClearArray(nEval, likelihoodValues);
      for (int i = 0; i < nEval; i++)
      {
        outputValues[i] += fracCats[iCat] * graph_[histoName_]->Eval(evalPoints[i]);
        if ( iCat == 0 ) likelihoodValues[i] = -TMath::Log(graph_[histoName_]->Eval(evalPoints[i]));
        else likelihoodValues[i] = -TMath::Log(outputValues[i]);
        if ( iCat != nCats_-1)
        {
          outputValuesTemp[i] += fracCatsTemp[iCat] * graph_[histoName_]->Eval(evalPoints[i]);
          likelihoodValuesTemp[i] = -TMath::Log(outputValuesTemp[i]);
        }
      }
      
      if ( iCat == 0 )  // outputValues = CM
      {
        histoTotal_[stringSuffix_[0]] = (TH1D*) histoSm_[histoName_]->Clone(("TotalProbability_"+stringSuffix_[0]).c_str());
        histoTotal_[stringSuffix_[0]]->SetTitle("#frac{n_{CM}}{n_{tot}} * f_{CM}(x|#Gamma) + #frac{n_{WM}}{n_{tot}} * f_{WM}(x) + #frac{n_{UM}}{n_{tot}} * f_{UM}(x)");
        
        this->MakeGraph(0, nEval, evalPoints, likelihoodValues, "likelihood_CM_", false);
        this->WriteOutput(nEval, widthArray_[iWidth], evalPoints, likelihoodValues, "CorrectMatchLikelihood_"+stringSuffix_[0], 1);
      }
      else if ( listCats_[iCat].find("WM") != std::string::npos && iCat != nCats_-1 )  // outputValues = CM + WM
      {
        histoTotal_[stringSuffix_[0]]->Add(histoSm_[histoName_]);
        
        this->MakeGraph(0, nEval, evalPoints, likelihoodValuesTemp, "likelihood_CMWM_", false);
        this->WriteOutput(nEval, widthArray_[iWidth], evalPoints, likelihoodValuesTemp, "MatchLikelihood_"+stringSuffix_[0], 1);
      }
      else if ( iCat == nCats_-1 ) // outputValues = CM + WM + UM
      {
        histoTotal_[stringSuffix_[0]]->Add(histoSm_[histoName_]);
        histoTotal_[stringSuffix_[0]]->Write();
        if (verbose_) std::cout << "Likelihood::ConstructTGraphs: The integral of the weighted probability histogram is " << histoTotal_[stringSuffix_[0]]->Integral(histoTotal_[stringSuffix_[0]]->FindBin(minRedMass_), histoTotal_[stringSuffix_[0]]->FindBin(maxRedMass_)+1) << std::endl;
        this->MakeGraph(0, nEval, evalPoints, outputValues, "TotalProbability_", false);
        this->DrawGraph(histoTotal_[stringSuffix_[0]], graph_["TotalProbability_"+stringSuffix_[0]], "Graph_totalProbability_"+stringSuffix_[0]);
//         if (verbose_)
//         {
//           graph_["TotalProbability_"+stringSuffix_[0]+"_test"] = graph_["TotalProbability_"+stringSuffix_[0]];
//           graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->SetPoint(graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->GetN(), maxRedMass_, 0.);
//           graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->SetPoint(graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->GetN(), minRedMass_, 0.);
//           std::cout << "Likelihood::ConstructTGraphs: The integral of the weighted probability graph is " << graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->Integral() << std::endl;
//         }
        
        this->MakeGraph(0, nEval, evalPoints, likelihoodValues, "likelihood_", false);
        this->WriteOutput(nEval, widthArray_[iWidth], evalPoints, likelihoodValues, stringSuffix_[0], 1);  // for TGraph
        this->WriteOutput(nEval, widthArray_[iWidth], evalPoints, likelihoodValues, stringSuffix_[0], 2);  // for TGraph2D
      }
      
    }  // end cats
    
    /// Fill TGraph2D
    for (int iCentre = 0; iCentre < nEval; iCentre++)
    {
      //if ( likelihoodValues[iCentre] < 1. )
      //if ( widthArray_[iWidth] > 0.35 && widthArray_[iWidth] < 0.65 && evalPoints[iCentre] > 0.7 && evalPoints[iCentre] < 0.78 )
      //  std::cout << "Point " << gLL2D_->GetN() << ":  " << evalPoints[iCentre] << "  " << widthArray_[iWidth] << "   " << likelihoodValues[iCentre] << std::endl;
      gLL2D_->SetPoint(gLL2D_->GetN(), evalPoints[iCentre], widthArray_[iWidth], likelihoodValues[iCentre]);
    }
    
  }  // end widths
  
  this->CombineOutput();
  
  this->DrawLikelihoods();
  gLL2D_->SetName("2D_likelihood");
  gLL2D_->SetTitle("2D_likelihood");
  gLL2D_->Write();
  this->DrawGraph(gLL2D_, "2D_likelihood");
  
  
  fileTGraphs_->Close();
  delete fileTGraphs_;
  delete gLL2D_;
}

bool Likelihood::ConstructTGraphsFromFile()
{
  return this->ConstructTGraphsFromFile("", outputDirName_, "");
}

bool Likelihood::ConstructTGraphsFromFile(std::string name)
{
  return this->ConstructTGraphsFromFile(name, outputDirName_, "");
}

bool Likelihood::ConstructTGraphsFromFile(std::string dir, std::string var)
{
  return this->ConstructTGraphsFromFile("", dir, var);
}

bool Likelihood::ConstructTGraphsFromFile(std::string name, std::string dir, std::string var)
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    suffix_ = name+"widthx"+stringWidthArray_[iWidth];
    inputFileName_ = dir+dirNameTGraphTxt_+"output_tgraph1d_"+suffix_+".txt";
    if ( ! tls_->fexists(inputFileName_.c_str()) )
    {
      std::cerr << "Likelihood::ConstructTGraphs: File " << inputFileName_ << " not found!!" << std::endl;
      std::cerr << "                              Aborting the likelihood calculation..." << std::endl;
      return false;
    }
    graph_[var+suffix_] = new TGraph(inputFileName_.c_str());
  }
  std::cout << "Likelihood::ConstructTGraphs: Constructed TGraphs for likelihood measurements (using " << nWidths_ << " widths)" << std::endl;
  return true;
}

bool Likelihood::ReadInput(std::string name)
{
  std::string line;
  double thisCentre, thisContent;
  (vecBinCentres_[name]).clear();
  (vecBinContents_[name]).clear();
  
  inputFileName_ = outputDirName_+dirNameTGraphTxt_+"output_func_"+name+".txt";
  if ( ! tls_->fexists(inputFileName_.c_str()) )
  {
    std::cerr << "Likelihood::ConstructTGraphs: File " << inputFileName_ << " not found!!" << std::endl;
    //std::cerr << "                          Aborting the likelihood calculation..." << std::endl;
    return false;
  }
  
  fileIn_.open(inputFileName_.c_str());
  if (verbose_) std::cout << "Opening " << inputFileName_ << "..." << std::endl;
  while( getline(fileIn_, line) )
  {
    std::istringstream iss(line);
    iss >> thisCentre >> thisContent;
    (vecBinCentres_[name]).push_back(thisCentre);
    (vecBinContents_[name]).push_back(thisContent);
  }
  fileIn_.close();
  return true;
}

bool Likelihood::ConstructTGraphsFromFile(std::vector<std::string> datasetNames, std::vector<int> includeDataset)
{
  const int nEval = 50;
  double evalPoints[nEval], outputValues[nEval], likelihoodValues[nEval];
  this->ClearArray(nEval, evalPoints);
  this->MakeTable(evalPoints, nEval, minRedMass_, maxRedMass_);
  
  /// Get fractions
  double fracs[nCats_];
  this->ClearArray(nCats_, fracs);
  this->GetFractions(fracs, nCats_, datasetNames, includeDataset);
  
  fileTGraphs_ = new TFile((dirNameLLTxt_+"TGraphs.root").c_str(), "RECREATE");
  std::cout << "Likelihood::ConstructTGraphs: Creating output file " << dirNameLLTxt_+"TGraphs.root" << std::endl;
  fileTGraphs_->cd();
  
  /// WM & UM distribution are independent of the width
  for (int iCat = 1; iCat < nCats_; iCat++)
  {
    stringSuffix_[iCat] = "widthx1";
    histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat];
    
    if (! this->ReadInput(histoName_)) return false;
  }
  
  /// Make arrays as input for TGraph
  const int nPoints = (vecBinCentres_[listCats_[1]+"_"+stringSuffix_[1]]).size();
  double binCentreArray[nPoints], binContentArray[nCats_][nPoints];
  this->ClearArray(nPoints, binCentreArray);
  for (int iCat = 0; iCat < nCats_; iCat++) this->ClearArray(nPoints, binContentArray[iCat]);
  
  for (int i = 0; i < nPoints; i++)
  {
    binCentreArray[i] = (vecBinCentres_[listCats_[1]+"_"+stringSuffix_[1]]).at(i);
    for (int iCat = 1; iCat < nCats_; iCat++)
    {
      binContentArray[iCat][i] = (vecBinContents_[listCats_[iCat]+"_"+stringSuffix_[iCat]]).at(i);
    }
  }
  
  /// Make TGraphs
  for (int iCat = 1; iCat < nCats_; iCat++)
    this->MakeGraphSmooth(iCat, nPoints, binCentreArray, binContentArray[iCat], listCats_[iCat]+"_");
  
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    stringSuffix_[0] = "widthx"+stringWidthArray_[iWidth];
    histoName_ = listCats_[0]+"_"+stringSuffix_[0];
    
    if (! this->ReadInput(histoName_)) return false;
    
    this->ClearArray(nPoints, binContentArray[0]);
    for (int i = 0; i < nPoints; i++)
      binContentArray[0][i] = (vecBinContents_[histoName_]).at(i);
    
    this->MakeGraphSmooth(0, nPoints, binCentreArray, binContentArray[0], listCats_[0]+"_");
    
    /// Make likelihood functions
    this->ClearArray(nEval, outputValues);
    this->ClearArray(nEval, likelihoodValues);
    for (int i = 0; i < nEval; i++)
    {
      for (int iCat = 0; iCat < nCats_; iCat++)
      {
        histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat];
        histoNameSm_ = listCats_[iCat]+"_Sm_"+stringSuffix_[iCat];
        outputValues[i] += fracs[iCat] * graph_[histoName_]->Eval(evalPoints[i]);
        
      }
      likelihoodValues[i] = -TMath::Log(outputValues[i]);
    }
    
    this->MakeGraph(0, nEval, evalPoints, likelihoodValues, "");
    
  }  // end widths
  
  std::cout << "Likelihood::ConstructTGraphs: Constructed TGraphs for likelihood measurements (using " << nWidths_ << " widths)" << std::endl;
  return true;
}

void Likelihood::CalculateLikelihood(double redMass, double relativeSF, bool isData)
{
  this->CalculateLikelihood(redMass, relativeSF, 1., 1., 1., false, isData);  // isTTbar = false ==> thisWidthSF_ = 1.;
}

std::vector<double> Likelihood::CalculateLikelihood(double redMass, double relativeSF, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData)
{
  return this->CalculateLikelihood(redMass, relativeSF, hadTopMassForWidthSF, lepTopMassForWidthSF, inputWidth, 172.5, isTTbar, isData);
}

std::vector<double> Likelihood::CalculateLikelihood(double redMass, double relativeSF, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  if (! isData && ! calledLLCalculation_) calledLLCalculation_ = true;
  
  vecLogLike_.clear();
  
  if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (isTTbar)
    {
      if (rewHadOnly_) thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, inputMass, 1., inputWidth);
      else thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, inputMass, 1., inputWidth) * rew_->BEventWeightCalculatorNonRel(lepTopMassForWidthSF, 172.5, inputMass, 1., inputWidth);
    }
    else thisWidthSF_ = 1.;
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      histoName_ = "widthx"+stringWidthArray_[iWidth];
      
      loglike_per_evt_[iWidth] = graph_[histoName_]->Eval(redMass);
      vecLogLike_.push_back(loglike_per_evt_[iWidth]);
      if (! isData) loglike_[iWidth] += loglike_per_evt_[iWidth]*relativeSF*thisWidthSF_;
      else loglike_data_[iWidth] += loglike_per_evt_[iWidth];
    }
  }
  
  return vecLogLike_;
}

void Likelihood::CalculateLikelihood(double redMass, double lepMass, double relativeSF, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  if (! isData && ! calledLLCalculation_) calledLLCalculation_ = true;
  
  if ( redMass > minRedMass_ && redMass < maxRedMass_ && lepMass > minRedMass2_ && lepMass < maxRedMass2_ )
  {
    if (isTTbar)
    {
      if (rewHadOnly_) thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, inputMass, 1., inputWidth);
      else thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, inputMass, 1., inputWidth) * rew_->BEventWeightCalculatorNonRel(lepTopMassForWidthSF, 172.5, inputMass, 1., inputWidth);
    }
    else thisWidthSF_ = 1.;
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      histoName_ = "widthx"+stringWidthArray_[iWidth];
      
      loglike_per_evt_[iWidth] = graph_["had_"+histoName_]->Eval(redMass)+graph_["lep_"+histoName_]->Eval(lepMass);
      if (! isData) loglike_[iWidth] += loglike_per_evt_[iWidth]*relativeSF*thisWidthSF_;
      else loglike_data_[iWidth] += loglike_per_evt_[iWidth];
    }
  }
}

void Likelihood::CalculateCMLikelihood(double redMass, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  this->CalculateCMLikelihood(redMass, scaleFactor, hadTopMassForWidthSF, lepTopMassForWidthSF, inputWidth, isTTbar, isData);
}

void Likelihood::CalculateCMLikelihood(double redMass, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "Likelihood::Cannot calculate loglikelihood for matched events when running over data" << std::endl;
    std::cerr << "            Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledCMLLCalculation_) calledCMLLCalculation_ = true;
    
    if (isTTbar)
    {
      if (rewHadOnly_) thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., inputWidth);
      else thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., inputWidth) * rew_->BEventWeightCalculatorNonRel(lepTopMassForWidthSF, 172.5, 172.5, 1., inputWidth);
    }
    else thisWidthSF_ = 1.;
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      thisWidth_ = stringWidthArray_[iWidth];
      
      loglike_CM_per_evt_[iWidth] = graph_["CorrectMatchLikelihood_widthx"+thisWidth_]->Eval(redMass);
      loglike_CM_[iWidth] += loglike_CM_per_evt_[iWidth]*scaleFactor*thisWidthSF_;
    }
  }
}

void Likelihood::CalculateTempLikelihood(double redMass, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  this->CalculateTempLikelihood(redMass, scaleFactor, hadTopMassForWidthSF, lepTopMassForWidthSF, inputWidth, isTTbar, isData);
}

void Likelihood::CalculateTempLikelihood(double redMass, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "Likelihood::Cannot calculate loglikelihood for correctly/wrondly matched events when running over data" << std::endl;
    std::cerr << "            Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledTempLLCalculation_) calledTempLLCalculation_ = true;
    
    if (isTTbar)
    {
      if (rewHadOnly_) thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., inputWidth);
      else thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., inputWidth) * rew_->BEventWeightCalculatorNonRel(lepTopMassForWidthSF, 172.5, 172.5, 1., inputWidth);
    }
    else thisWidthSF_ = 1.;
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      thisWidth_ = stringWidthArray_[iWidth];
      
      loglike_temp_per_evt_[iWidth] = graph_["MatchLikelihood_widthx"+thisWidth_]->Eval(redMass);
      loglike_temp_[iWidth] += loglike_temp_per_evt_[iWidth]*scaleFactor*thisWidthSF_;
    }
  }
}

void Likelihood::CalculateGenLikelihood(double redMass, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "Likelihood::Cannot calculate loglikelihood for generated events when running over data" << std::endl;
    std::cerr << "            Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledGenLLCalculation_) calledGenLLCalculation_ = true;
    
    if (isTTbar)
    {
      if (rewHadOnly_) thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., inputWidth);
      else thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., inputWidth) * rew_->BEventWeightCalculatorNonRel(lepTopMassForWidthSF, 172.5, 172.5, 1., inputWidth);
    }
    else thisWidthSF_ = 1.;
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      thisWidth_ = stringWidthArray_[iWidth];
      
      loglike_gen_per_evt_[iWidth] = graph_["CorrectMatchLikelihood_widthx"+thisWidth_]->Eval(redMass);
      loglike_gen_[iWidth] += loglike_gen_per_evt_[iWidth]*thisWidthSF_;
    }
  }
}

void Likelihood::GetOutputWidth(double inputWidth, double inputMass, bool writeToFile, bool makeNewFile)
{
  this->GetOutputWidth(inputWidth, inputMass, "", writeToFile, makeNewFile);
}

void Likelihood::GetOutputWidth(double inputWidth, double inputMass, std::string type, bool writeToFile, bool makeNewFile)
{
  std::string loglikePlotName = "loglikelihood_vs_width_";
  if (! type.empty() ) loglikePlotName += type+"_";
  loglikePlotName += "widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass);
  
  if ( type.find("CM") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_CM_, loglikePlotName, writeToFile, makeNewFile, type);
  else if ( type.find("match") != std::string::npos || type.find("Match") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_temp_, loglikePlotName, writeToFile, makeNewFile, type);
  else if ( type.find("gen") != std::string::npos || type.find("Gen") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_gen_, loglikePlotName, writeToFile, makeNewFile, type);
  else if ( type.find("data") != std::string::npos || type.find("Data") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_data_, loglikePlotName, writeToFile, makeNewFile, type);
  else
    output_ = this->CalculateOutputWidth(nWidths_, loglike_, loglikePlotName, writeToFile, makeNewFile, type);
  
  std::cout << "For an input width of " << inputWidth << " the minimum can be found at " << output_.first << " and the uncertainty is " << output_.second << std::endl;
  
  std::string fileName = dirNameLLTxt_+"result_minimum_";
  if (! type.empty() ) fileName += type+"_";
  fileName += "widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass)+".txt";
  txtOutputLL_.open(fileName.c_str());
  if (! type.empty() ) txtOutputLL_ << std::setw(18) << std::left << type << "  ";
  else txtOutputLL_ << std::setw(18) << std::left << "nominal  ";
  txtOutputLL_ << std::setw(5) << std::left << std::setprecision(5) << inputWidth << "   " << std::setw(5) << std::left << std::setprecision(5) << inputMass << "   " << std::setw(20) << std::right << std::setprecision(20) << output_.first << "  " << std::setw(20) << std::right << std::setprecision(20) << output_.second << std::endl;
  txtOutputLL_.close();
}

void Likelihood::GetOutputWidth(std::string inputFileName, double inputWidth, double inputMass, bool writeToFile, bool makeNewFile)
{
  this->GetOutputWidth(inputFileName, dirNameLLTxt_, inputWidth, inputMass, writeToFile, makeNewFile);
}

void Likelihood::GetOutputWidth(std::string inputFileName, std::string inputDir, double inputWidth, double inputMass, bool writeToFile, bool makeNewFile)
{
  if (verbose_) std::cout << "Using LogLikelihood values from file" << std::endl;
  std::string loglikePlotName = "loglikelihood_vs_width_";
  //if (! type.empty() ) loglikePlotName += type+"_";
  loglikePlotName += "ff_widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass);
  
  output_ = this->CalculateOutputWidth(inputFileName, inputDir, loglikePlotName, writeToFile, makeNewFile);
  
  std::cout << "For an input width of " << inputWidth << " the minimum can be found at " << output_.first << " and the uncertainty is " << output_.second << std::endl;
}

std::pair<double,double> Likelihood::GetOutputWidth(double inputWidth, double inputMass, int thisPsExp)
{
  std::string loglikePlotName = dirNamePull_+"loglikelihood_vs_width_psExp_"+tls_->ConvertIntToString(thisPsExp)+"_widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass);
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    loglike_pull_single_[iWidth] = loglike_pull_[iWidth][thisPsExp];
  }
  
  return this->CalculateOutputWidth(nWidths_, loglike_pull_single_, loglikePlotName, true, false);
}

std::pair<double,double> Likelihood::CalculateOutputWidth(std::string inputFileName, std::string inputDir, std::string plotName, bool writeToFile, bool makeNewFile)
{
  this->ReadLLValuesFromFile(inputFileName, inputDir);
  const int nn = vecWidthFromFile_.size();
  if ( nn == 0 )
    return std::pair<double,double>(-1.,-1.);
  double arrWidth[nn], arrLLVals[nn];
  for (int i = 0; i < nn; i++)
  {
    arrWidth[i] = vecWidthFromFile_.at(i);
    arrLLVals[i] = vecLLValsFromFile_.at(i);
  }
  return this->CalculateOutputWidth(nn, arrWidth, arrLLVals, plotName, writeToFile, makeNewFile);
}

std::pair<double,double> Likelihood::CalculateOutputWidth(int nn, double* LLvalues, std::string plotName, bool writeToFile, bool makeNewFile, std::string type)
{
  return this->CalculateOutputWidth(nn, (double*)widthArray_, LLvalues, plotName, writeToFile, makeNewFile, type);
}

std::pair<double,double> Likelihood::CalculateOutputWidth(int nn, double* evalWidths, double* LLvalues, std::string plotName, bool writeToFile, bool makeNewFile, std::string type)
{
  TGraph *g = new TGraph(nn, evalWidths, LLvalues);
  
  int locMin = LocMinArray(nn, LLvalues);
  //std::cout << "Index of minimum LL value is " << locMin << std::endl;
  double centreVal = evalWidths[locMin];
  if ( centreVal <= 0.3 )
  {
    double tempArray[nn-3];
    for (int i = 0; i < nn-3; i++) { tempArray[i] = LLvalues[i+3]; }
    int tempMin = LocMinArray(nn-3, tempArray); //std::cout << tempMin << "  " << evalWidths[tempMin] << std::endl;
    if ( LLvalues[locMin+2] > tempArray[tempMin+2] ) centreVal = evalWidths[tempMin+3];
  }
  
  
  //double interval = 0.4;
  double interval = 0.3;  // mlb
  if ( centreVal <= interval ) interval = centreVal - 0.1;
  //if ( centreVal > 3.8 ) interval = 0.8;
  if ( centreVal > 3.8 ) interval = 0.6;  // mlb
  double fitmax = centreVal + interval;
  double fitmin = centreVal - interval;
  
  if ( centreVal > 0.15 && fitmin < 0.15 ) fitmin = 0.15;
  if ( centreVal > 0.2 && fitmin < 0.2 ) fitmin = 0.2;
  if ( centreVal > 0.35 && fitmin < 0.3 ) fitmin = 0.3;
  //if ( centreVal > 0.4 && fitmin < 0.3 ) fitmin = 0.3;
  if ( centreVal > 0.5 && fitmin < 0.35 ) fitmin = 0.35;
  if ( centreVal > 0.55 && fitmin < 0.4 ) fitmin = 0.4;
  if ( centreVal > 0.65 && fitmin < 0.45 ) fitmin = 0.45;
  if ( centreVal > 0.75 && fitmin < 0.5 ) fitmin = 0.5;
  if ( centreVal > 1.1 && fitmin < 0.8 ) fitmin = 0.8;
  
  // mlb
  if ( centreVal > 0.3 && centreVal < 1.6 ) fitmin += 0.05;
  //if ( centreVal < 1.1 )  fitmin += 0.05;  // --> += 0.1
  if ( fitmin >= centreVal ) fitmin = centreVal - 0.05;
  
  
  if ( centreVal < 0.35 && fitmax > 0.4 ) fitmax = 0.4;
  if ( centreVal > 0.35 && centreVal < 1.2 ) fitmax = centreVal + (centreVal - fitmin);
  
  
  if (verbose_) std::cout << "Likelihood::CalculateOutputWidth: Look for minimum in interval [" << fitmin << "," << fitmax << "] around " << centreVal << std::endl;
  
  TF1 *func = new TF1("func", "[0]+[1]*x+[2]*x*x", fitmin, fitmax);
  
  g->Fit(func,"RWME");
  
  double outputWidth = func->GetMinimumX(fitmin, fitmax);
  
  double LLreduced[nn];
  for (int i = 0; i < nn; i++)
    LLreduced[i] = LLvalues[i] - func->Eval(outputWidth);
  TGraph *g2 = new TGraph(nn, evalWidths, LLreduced);
  g2->Fit(func,"RWME");
  
  double funcMin = fitmin - 3*interval;
  if ( funcMin < 0. ) funcMin = 0.;
  double funcMax = fitmax + 3*interval;
  TF1 *parabola = new TF1("parabola", fitFunc, funcMin, funcMax, 3);
  parabola->SetParameters(func->GetParameter(0), func->GetParameter(1), func->GetParameter(2));
  
  outputWidth = parabola->GetMinimumX(fitmin, fitmax);
  double lowerSigma = parabola->GetX(parabola->Eval(outputWidth) + 0.5, fitmin-interval, outputWidth);
  double upperSigma = parabola->GetX(parabola->Eval(outputWidth) + 0.5, outputWidth, fitmax+interval);
  double sigma = (upperSigma - lowerSigma)/2.;
  if ( lowerSigma <= fitmin && upperSigma >= fitmax )
    std::cerr << "Likelihood::CalculateOutputWidth: ERROR: Uncertainty calculation limited by fit boundaries. Do not trust..." << std::endl;
  else if ( lowerSigma <= fitmin ) sigma = upperSigma - outputWidth;
  else if ( upperSigma >= fitmax ) sigma = outputWidth - lowerSigma;
  
  std::cout << "Minimum -log(like) value is " << parabola->Eval(outputWidth) << std::endl;
  std::cout << "===> Minimum at " << outputWidth << " +- " << sigma << std::endl;
  
  double lower95 = parabola->GetX(parabola->Eval(outputWidth) + 2., funcMin, outputWidth);
  double upper95 = parabola->GetX(parabola->Eval(outputWidth) + 2., outputWidth, funcMax);
  double sigma95 = (upper95 - lower95)/2.;
  if ( lower95 <= funcMin && upper95 >= funcMax )
    std::cerr << "Likelihood::CalculateOutputWidth: ERROR: 95% uncertainty calculation limited by function boundaries. Do not trust..." << std::endl;
  else if ( lower95 <= funcMin ) sigma95 = upper95 - outputWidth;
  else if ( upper95 >= funcMax ) sigma95 = outputWidth - lower95;
  
  std::cout << "===> 95% CL :   " << sigma95 << std::endl;
  
  // Write 95% output
  if ( plotName.find("widthx1_mass172p5") != std::string::npos )
  {
    std::string fileName = dirNameLLTxt_+"sigma95_result_minimum_";
    if (! type.empty() ) fileName += type+"_";
    fileName += "widthx1_mass172p5.txt";
    txtOutputLL_.open(fileName.c_str());
    if (! type.empty() ) txtOutputLL_ << std::setw(18) << std::left << type << "  ";
    else txtOutputLL_ << std::setw(18) << std::left << "nominal  ";
    txtOutputLL_ << "1   172.5   " << std::setw(20) << std::right << std::setprecision(20) << outputWidth << "  " << std::setw(20) << std::right << std::setprecision(20) << sigma95 << std::endl;
    txtOutputLL_.close();
  }
  
  
  if (makeNewFile)
  {
    filePlots_ = new TFile((dirNameLLTxt_+"File_"+plotName+".root").c_str(), "RECREATE");
    filePlots_->cd();
  }
  
  this->DrawOutputLogLikelihood(g2, parabola, 0, 15, 1200., plotName+"_full", writeToFile);
  this->DrawOutputLogLikelihood(g2, parabola, outputWidth-2., outputWidth+2., 1.5*g2->Eval(outputWidth+2.), plotName, writeToFile);
  this->DrawOutputLogLikelihood(g2, parabola, outputWidth-0.5, outputWidth+0.5, std::max(g2->Eval(outputWidth-3.*sigma),g2->Eval(outputWidth+0.5)), plotName+"_zoom", writeToFile);
  
  if (makeNewFile)
  {
    filePlots_->Close();
    delete filePlots_;
  }
  
  delete g2;
  delete g;
  delete parabola;
  
  return std::pair<double,double>(outputWidth,sigma);
}

std::pair<double,double> Likelihood::ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma)
{
  /// return 'thisInputWidth' and uncertainty
  //  thisOutputWidth = Par_[0] + Par_[1] * thisInputWidth

  double thisInputWidth = (thisOutputWidth - calCurvePar_[0])/calCurvePar_[1];
  double thisInputWidthSigma = TMath::Sqrt( thisOutputWidthSigma*thisOutputWidthSigma + calCurveParUnc_[0]*calCurveParUnc_[0] + calCurveParUnc_[1]*calCurveParUnc_[1]*thisInputWidth*thisInputWidth )/calCurvePar_[1];
  
  return std::pair<double,double>(thisInputWidth,thisInputWidthSigma);
  
}

int Likelihood::InitPull(int nPsExp)
{
  mkdir((dirNameLLTxt_+dirNamePull_).c_str(),0777);
  
  nPsExp_ = nPsExp;
  if ( nPsExp > 1000 )
  {
    std::cout << "Likelihood::Pull: Warning: Only 1000 pseudo experiments will be performed" << std::endl;
    return 1000;
  }
  else if (verbose_) std::cout << "Likelihood::Pull: Performing " << nPsExp << " pseudo experiments" << std::endl;
  
  return nPsExp;
}

void Likelihood::AddPsExp(int thisPsExp, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  this->AddPsExp(thisPsExp, scaleFactor, hadTopMassForWidthSF, lepTopMassForWidthSF, inputWidth, isTTbar, isData);
}

void Likelihood::AddPsExp(int thisPsExp, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, bool isTTbar, bool isData)
{
  if (isData) std::cerr << "Likelihood::Pull: Will not use data for pseudo experiments..." << std::endl;
  else if (calledLLCalculation_)
  {
    if (isTTbar)
    {
      if (rewHadOnly_) thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., inputWidth);
      else thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., inputWidth) * rew_->BEventWeightCalculatorNonRel(lepTopMassForWidthSF, 172.5, 172.5, 1., inputWidth);
    }
    else thisWidthSF_ = 1.;
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      loglike_pull_[iWidth][thisPsExp] += loglike_per_evt_[iWidth]*scaleFactor*thisWidthSF_;
    }
  }
  else std::cerr << "Likelihood::Pull: Did not calculate likelihoods! Cannot get input for pseudo experiments..." << thisPsExp << std::endl;
}

void Likelihood::CalculatePull(double inputWidth)
{
  this->CalculatePull(inputWidth, 172.5);
}

void Likelihood::CalculatePull(double inputWidth, double inputMass)
{
  std::string fileName = dirNameLLTxt_+dirNamePull_+"PseudoExperiments_widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass)+".root";
  TFile *filePull = new TFile(fileName.c_str(),"RECREATE");
  filePull->cd();
  TH1D *hPull = new TH1D("hPull", "; (#Gamma_{j} - <#Gamma>)/#sigma_{#Gamma_{j}}", 32, -4., 4.);
  TH1D *hAveWidthOut = new TH1D("hAveWidthOut", "; #Gamma_{out}", 40, 0., 2.);
  TH1D *hAveWidthIn = new TH1D("hAveWidthIn", "; #Gamma_{in}", 40, 0., 2.);
  TH1D *hUncAveWidthOut = new TH1D("hUncAveWidthOut", "; #sigma_{#Gamma_{out}}", 30, 0., 0.3);
  TH1D *hUncAveWidthIn = new TH1D("hUncAveWidthIn", "; #sigma_{#Gamma_{in}}", 30, 0., 0.3);
  
  /// Calculate output width for all pseudo experiments
  //  Transform to input width via calibration curve & calculate average input width
  std::pair<double,double> thisOutputWidth[nPsExp_], thisInputWidth[nPsExp_];
  double aveInputWidth = 0.;
  
  for (int iPsExp = 0; iPsExp < nPsExp_; iPsExp++)
  {
    // Clear
    thisOutputWidth[iPsExp] = std::pair<double,double>(-1.,-1.);
    thisInputWidth[iPsExp] = std::pair<double,double>(-1.,-1.);
    
    // Fill
    thisOutputWidth[iPsExp] = this->GetOutputWidth(inputWidth, inputMass, iPsExp);
    thisInputWidth[iPsExp] = this->ApplyCalibrationCurve(thisOutputWidth[iPsExp].first, thisOutputWidth[iPsExp].second);
    if ( thisInputWidth[iPsExp].first != -1. ) aveInputWidth += thisInputWidth[iPsExp].first;
    else std::cerr << "Likelihood::CalculatePull: Input width for pseudo experiment " << iPsExp << " is equal to -1! Ignoring this pseudo experiment... " << std::endl;
    hAveWidthOut->Fill(thisOutputWidth[iPsExp].first);
    hAveWidthIn->Fill(thisInputWidth[iPsExp].first);
    hUncAveWidthOut->Fill(thisOutputWidth[iPsExp].second);
    hUncAveWidthIn->Fill(thisInputWidth[iPsExp].second);
  }
  aveInputWidth = aveInputWidth/nPsExp_;
  //std::cout << "Output width : " << thisOutputWidth[0].first << " ;   sigma : " << thisOutputWidth[0].second << std::endl;
  //std::cout << "Input width  : " << thisInputWidth[0].first << " ;   sigma : " << thisInputWidth[0].second << std::endl;
  std::cout << "Average width for pseudo experiments is " << aveInputWidth << std::endl;
  
  this->WritePsExpOutput(thisOutputWidth, thisInputWidth, inputWidth);
  gStyle->SetOptStat(1110);  // display entries, mean & RMS, but no title
  hAveWidthOut->Write();
  hAveWidthIn->Write();
  hUncAveWidthOut->Write();
  hUncAveWidthIn->Write();
  gStyle->SetOptStat(10);  // display entries, but no title, mean & RMS
  TF1 *wOut = new TF1("wOut", "gaus", hAveWidthOut->GetMean()-0.3, hAveWidthOut->GetMean()+0.3);
  hAveWidthOut->Fit(wOut,"R");
  gStyle->SetOptFit(0111);
  hAveWidthOut->Write();
  wOut->Write();
  TF1 *wIn = new TF1("wIn", "gaus", hAveWidthIn->GetMean()-0.3, hAveWidthIn->GetMean()+0.3);
  hAveWidthIn->Fit(wIn,"R");
  gStyle->SetOptFit(0111);
  hAveWidthIn->Write();
  wIn->Write();
  
  /// Fill histogram with (Gamma_j - <Gamma>)/sigma_j
  gStyle->SetOptStat(1110);  // display entries, mean & RMS, but no title
  double fillValue;
  
  for (int iPsExp = 0; iPsExp < nPsExp_; iPsExp++)
  {
    if ( thisInputWidth[iPsExp].first != -1. )
    {
      fillValue = (thisInputWidth[iPsExp].first - aveInputWidth)/thisInputWidth[iPsExp].second;
      hPull->Fill(fillValue);
    }
  }
  hPull->Write();
  
  /// Fit distribution with Gaussian function
  //  Should have mean = 0 and sigma = 1
  //double fitMin = hPull->GetXaxis()->GetXmin();
  //double fitMax = hPull->GetXaxis()->GetXmax();
  
  gStyle->SetOptStat(10);  // display entries, but no title, mean & RMS
  TF1 *gaus = new TF1("gaus", "gaus", -3., 3.);
  hPull->Fit(gaus,"R");
  gStyle->SetOptFit(0111);
  hPull->Write();
  gaus->Write();
  
  std::cout << "The pseudo experiment distribution is fitted with a Gaussian function with parameters ";
  for (int iPar = 0; iPar < gaus->GetNpar(); iPar++)
  {
    std::cout << gaus->GetParameter(iPar) << "+-" << gaus->GetParError(iPar);
    if ( iPar != gaus->GetNpar()-1) std::cout << " ,  ";
  }
  std::cout << std::endl;
  
  // Close file
  filePull->Close();
  
  delete gaus;
  delete filePull;
}

void Likelihood::AddToFraction(int d, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, bool isTTbar, bool isCM, bool isWM, bool isUM)
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    /// Calculate widthSF
    if (isTTbar)
    {
      if (rewHadOnly_) thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., widthArray_[iWidth]);
      else thisWidthSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForWidthSF, 172.5, 172.5, 1., widthArray_[iWidth]) * rew_->BEventWeightCalculatorNonRel(lepTopMassForWidthSF, 172.5, 172.5, 1., widthArray_[iWidth]);
    }
    else thisWidthSF_ = 1.;
    
    if      (isCM) nEventsCMFractions_[iWidth][d] += scaleFactor * thisWidthSF_;
    else if (isWM) nEventsWMFractions_[iWidth][d] += scaleFactor * thisWidthSF_;
    else if (isUM) nEventsUMFractions_[iWidth][d] += scaleFactor * thisWidthSF_;
  }
}

void Likelihood::CalculateFractions(std::vector<std::string> datasetNames)
{
  std::string fileName = "";
  int nDatasets = datasetNames.size();
  mkdir(dirNameNEvents_.c_str(),0777);
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    for (int d = 0; d < nDatasets; d++)
    {
      fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_widthx"+stringWidthArray_[iWidth]+"_CM.txt";
      txtOutputFractions_.open(fileName.c_str());
      txtOutputFractions_ << nEventsCMFractions_[iWidth][d] << std::endl;
      txtOutputFractions_.close();
      fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_widthx"+stringWidthArray_[iWidth]+"_WM.txt";
      txtOutputFractions_.open(fileName.c_str());
      txtOutputFractions_ << nEventsWMFractions_[iWidth][d] << std::endl;
      txtOutputFractions_.close();
      fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_widthx"+stringWidthArray_[iWidth]+"_UM.txt";
      txtOutputFractions_.open(fileName.c_str());
      txtOutputFractions_ << nEventsUMFractions_[iWidth][d] << std::endl;
      txtOutputFractions_.close();
    }  // end datasets
  }  // end widths
}

void Likelihood::GetFractions(double *fractions, int nCats, std::vector<std::string> datasetNames, std::vector<int> includeDataset)
{
  Likelihood::GetFractions(fractions, nCats, "1", datasetNames, includeDataset);
}

void Likelihood::GetFractions(double *fractions, int nCats, std::string width, std::vector<std::string> datasetNames, std::vector<int> includeDataset)
{
  std::string line;
  std::string fileName;
  double val;
  const int nDatasets = datasetNames.size();
  std::vector<double> nEvents[nDatasets];
  double totEvents[nCats_];
  this->ClearArray(nCats_, totEvents);
  double totalNbEvents = 0.;
  
  for (int d = 0; d < nDatasets; d++)
  {
    nEvents[d].clear();
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_widthx"+width+"_"+listCats_[iCat]+".txt";
      if (! tls_->fexists(fileName.c_str()) )
      {
        std::cerr << "WARNING: File " << fileName << " does not exist." << std::endl;
        exit(1);
      }
      fileIn_.open(fileName.c_str());
      if (verbose_) std::cout << "Opening " << fileName << "..." << std::endl;
      while( getline(fileIn_, line) )
      {
        std::istringstream iss(line);
        while ( iss >> val )
        {
          //iss >> val;
          (nEvents[d]).push_back(val);
        }
      }
      fileIn_.close();
    }  // end cats
    if ( (nEvents[d]).size() != nCats_ )
    {
      std::cerr << "Something wrong with number of categories..." << std::endl;
      exit(1);
    }
    
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      totEvents[iCat] += (double)includeDataset[d] * (nEvents[d]).at(iCat);
    }
  }  // end datsets
  
  for (int iCat = 0; iCat < nCats; iCat++)
  {
    totalNbEvents += totEvents[iCat];
  }
  for (int iCat = 0; iCat < nCats; iCat++)
  {
    fractions[iCat] = totEvents[iCat] / totalNbEvents;
    if ( verbose_ && iCat == nCats-1 )
    {
      std::cout << "Width: " << width;
      for (int iCat = 0; iCat < nCats_; iCat++)
        std::cout << "   # " << listCats_[iCat] << ": " << 100.*fractions[iCat] << "%";
      std::cout << std::endl;
    }
  }
}

void Likelihood::DrawGraph(TH1D* h, TGraph* g, std::string name)
{
  TCanvas *c = new TCanvas(name.c_str(), name.c_str());
  c->cd();
  h->GetXaxis()->SetTitle("m_{r}");
  h->Draw("C");
  g->SetLineColor(kRed);
  g->Draw("same");
  c->Update();
  c->Write();
  c->SaveAs((outputDirName_+name+".png").c_str());
  
  gSystem->ProcessEvents();
  delete c;
}

void Likelihood::DrawGraph(TGraph2D* g, std::string name)
{
  TCanvas *c = new TCanvas(name.c_str(), name.c_str());
  c->cd();
  gStyle->SetPalette(1);
  g->SetMaxIter(500000);
  g->Draw("AP");
  g->SetName(name.c_str());
  g->SetTitle((name+"; m_{r}; #Gamma/#Gamma_{t}").c_str());
  g->GetXaxis()->SetTitleOffset(1.5);
  g->GetYaxis()->SetTitleOffset(1.5);
  g->Draw("surf1");
  c->Update();
  c->Write();
  c->SaveAs((outputDirName_+name+"_surf1.png").c_str());
  g->Draw("colz");
  c->Update();
  c->Write();
  c->SaveAs((outputDirName_+name+"_colz.png").c_str());
  g->SetMarkerSize(0.5);
  g->Draw("tri1 p0");
  c->Update();
  c->Write();
  c->SaveAs((outputDirName_+name+"_tri1p0.png").c_str());
  
  delete c;
}

void Likelihood::DrawLikelihoods()
{
  Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+1, kCyan+1, kBlue+2, kMagenta, kViolet-5, kPink+10};
  TCanvas* c2 = new TCanvas("-log(likelihood)", "-log(likelihood)");
  c2->cd();
  graph_["likelihood_widthx"+stringWidthArray_[0]]->SetLineColor(colours[0]);
  graph_["likelihood_widthx"+stringWidthArray_[0]]->Draw();
  std::string temp = graph_["likelihood_widthx"+stringWidthArray_[0]]->GetTitle();
  graph_["likelihood_widthx"+stringWidthArray_[0]]->SetTitle("-LogLikelihood");
  c2->Update();
  for (int i = 1; i < nWidths_; i++)
  {
    graph_["likelihood_widthx"+stringWidthArray_[i]]->SetLineColor(colours[i%10]);
    graph_["likelihood_widthx"+stringWidthArray_[i]]->Draw("same");
    c2->Update();
  }
  c2->Write();
  c2->SaveAs((outputDirName_+"LogLikelihood_multi.png").c_str());
  c2->Close();
  
  graph_["likelihood_widthx"+stringWidthArray_[0]]->SetTitle(temp.c_str());
  
  delete c2;
}

void Likelihood::DrawOutputLogLikelihood(TGraph* g, TF1* f, double minX, double maxX, double maxY, std::string name, bool writeToFile)
{
  std::string outputFileName = dirNameLLTxt_+name+".png";
  TCanvas* c1 = new TCanvas(name.c_str(), "-LogLikelihood vs. width");
  c1->cd();
  g->Draw("AP");
  if ( minX < 0. ) minX = 0.;
  g->SetTitle("");
  g->GetXaxis()->SetRangeUser(minX,maxX);
  g->GetYaxis()->SetRangeUser(-0.05*maxY,maxY);
  g->GetXaxis()->SetTitle("#Gamma_{t} /#Gamma_{t,gen}");
  g->GetYaxis()->SetTitle("\\Delta\\mathscr{L}(\\mathrm{\\Gamma_{t}})");
  g->SetMarkerStyle(2);  //kPlus
  g->Draw("AP");
  //if (writeToFile) g->Write();
  c1->Update();
  f->Draw("same");
  c1->Modified();
  c1->Update();
  TLine l;
  l.SetLineColor(921);
  l.DrawLine(minX, 0.5, maxX, 0.5);
  c1->Update();
  
//   TPaveStats* sb = (TPaveStats*) g->GetListOfFunctions()->FindObject("stats");
//   TList *listOfLines = sb->GetListOfLines();
//   TText *tconst = sb->GetLineWith("const");
//   listOfLines->Remove(tconst);
//   if ( name.find("zoom") != std::string::npos)
//   {
//     sb->SetBorderSize(0);
//     sb->SetX1NDC(0.30);
//     sb->SetX2NDC(0.70);
//     sb->SetY2NDC(0.89);
//   }
//   c1->Modified();
//   c1->Update();
  
  if (writeToFile) c1->Write();
  c1->SaveAs(outputFileName.c_str());
  c1->Close();
  
  delete c1;
}

void Likelihood::ReadLLValuesFromFile(std::string inputFileName, std::string inputDir)
{
  if (! tls_->fexists((inputDir+inputFileName).c_str()) )
  {
    std::cerr << "WARNING: File " << inputDir+inputFileName << " does not exist." << std::endl;
    exit(1);
  }
  fileIn_.open((inputDir+inputFileName).c_str());
  if (verbose_) std::cout << "Opening " << inputDir+inputFileName << "..." << std::endl;
  
  vecWidthFromFile_.clear();
  vecLLValsFromFile_.clear();
  
  std::string line;
  std::string var = "";
  double val = -1.;
  while( getline(fileIn_, line) )
  {
    std::istringstream iss(line);
    if ( line.find("Width") != std::string::npos )
    {
      iss >> var;
      while ( iss >> val )
      {
        vecWidthFromFile_.push_back(val);
      }
    }
    if ( line.find("LL") != std::string::npos )
    {
      iss >> var;
      while ( iss >> val )
      {
        vecLLValsFromFile_.push_back(val);
      }
    }
  }
  
  fileIn_.close();
  
  if ( vecWidthFromFile_.size() != vecLLValsFromFile_.size() )
  {
    std::cerr << "Likelihood::ReadLLValuesFromFile : Input file " << inputFileName << " cannot be read:" << std::endl;
    std::cerr << "                                   # widths != # likelihood values !!" << std::endl;
    vecWidthFromFile_.clear();
    vecLLValsFromFile_.clear();
  }
}

void Likelihood::WritePsExpOutput(std::pair<double,double> *outputWidth, std::pair<double,double> *inputWidth, double genWidth)
{
  std::string outputTxtName = dirNameLLTxt_+dirNamePull_+"PseudoExperiments_output_widthx"+tls_->DotReplace(genWidth)+".txt";
  txtOutputPsExp_.open(outputTxtName.c_str());
  txtOutputPsExp_ << "iPsExp   output width        unc           input width        unc" << std::endl;
  for (int iPsExp = 0; iPsExp < nPsExp_; iPsExp++)
  {
    txtOutputPsExp_ << std::setw(4) << std::right << iPsExp << "   " << std::setw(8) << std::left << outputWidth[iPsExp].first << "   " << std::setw(8) << std::left << outputWidth[iPsExp].second << "   " << std::setw(8) << std::left << inputWidth[iPsExp].first << std::setw(8) << std::left << inputWidth[iPsExp].second << std::endl;
  }
  txtOutputPsExp_.close();
}

void Likelihood::WriteFuncOutput(int nPoints, double *arrayCentre, double *arrayContent, std::string name)
{
  std::string outputTxtName = outputDirName_+dirNameTGraphTxt_+"output_func_"+name+".txt";
  txtOutput_.open(outputTxtName.c_str());
  
  for (int i = 0; i < nPoints; i++)
  {
    txtOutput_ << std::setw(8) << std::right << arrayCentre[i];
    txtOutput_ << "  " << std::setw(8) << std::right << arrayContent[i] << std::endl;
  }
  txtOutput_.close();
}

void Likelihood::WriteOutput(int nPoints, double width, double *arrayCentre, double *arrayContent, std::string name, int dim)
{
  std::string dimension = "";
  if ( dim == 1 ) dimension = "1d";
  else if ( dim == 2 ) dimension = "2d";
  std::string outputTxtName = outputDirName_+dirNameTGraphTxt_+"output_tgraph"+dimension+"_"+name+".txt";
  txtOutput_.open(outputTxtName.c_str());
  
  for (int i = 0; i < nPoints; i++)
  {
    txtOutput_ << std::setw(8) << std::right << arrayCentre[i];
    if ( dim == 2 ) txtOutput_ << "   " << std::setw(5) << std::right << width;
    txtOutput_ << "  " << std::setw(8) << std::right << arrayContent[i] << std::endl;
  }
  txtOutput_.close();
}

void Likelihood::CombineOutput()
{
  std::ofstream combFile((outputDirName_+dirNameTGraphTxt_+"output_tgraph2d_total"+".txt").c_str());
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    std::string inFileName = outputDirName_+dirNameTGraphTxt_+"output_tgraph2d_widthx"+stringWidthArray_[iWidth]+".txt";
    std::ifstream infile(inFileName.c_str());
    combFile << infile.rdbuf();
    infile.close();
  }
  combFile.close();
}

void Likelihood::PrintLikelihoodOutput(std::string llFileName)
{
  std::ofstream txtOutputLogLike((dirNameLLTxt_+llFileName).c_str());
  txtOutputLogLike << "Widths:      ";
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    txtOutputLogLike << std::setw(5) << std::right << widthArray_[iWidth] << "  ";
  }
  txtOutputLogLike << std::endl << "LLikelihood: ";
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    txtOutputLogLike << std::setw(12) << std::right << loglike_[iWidth] << "  ";
  }
  txtOutputLogLike << std::endl;
  txtOutputLogLike.close();
}

void Likelihood::PrintLikelihoodOutputData(std::string llFileName)
{
  std::ofstream txtOutputLogLike((dirNameLLTxt_+llFileName).c_str());
  txtOutputLogLike << "Widths:      ";
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    txtOutputLogLike << std::setw(5) << std::right << widthArray_[iWidth] << "  ";
  }
  txtOutputLogLike << std::endl << "LLikelihood: ";
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    txtOutputLogLike << std::setw(17) << std::right << std::setprecision(15) << loglike_data_[iWidth] << "  ";
  }
  txtOutputLogLike << std::endl;
  txtOutputLogLike.close();
}

void Likelihood::PrintMtmLikelihoodOutput(std::string llFileName)
{
  std::ofstream txtOutputLogLike((dirNameLLTxt_+llFileName).c_str());
  txtOutputLogLike << "likelihood values (all MC) : {";
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    txtOutputLogLike << std::fixed << std::setprecision(15) << loglike_[iWidth];
    if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
  }
  txtOutputLogLike << "};" << std::endl;
//   txtOutputLogLike << "likelihood values (data-only) : {";
//   for (int iWidth = 0; iWidth < nWidths_; iWidth++)
//   {
//     txtOutputLogLike << std::fixed << std::setprecision(15) << loglike_data_[iWidth][1];
//     if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
//   }
//   txtOutputLogLike << "};" << std::endl;
  if (calledCMLLCalculation_)
  {
    txtOutputLogLike << "CM likelihood values : {";
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      txtOutputLogLike << std::fixed << std::setprecision(15) << loglike_CM_[iWidth];
      if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
    }
    txtOutputLogLike << "};";
  }
  txtOutputLogLike << std::endl << "widths: ";
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    txtOutputLogLike << widthArray_[iWidth];
    if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
  }
  if (calledGenLLCalculation_)
  {
    txtOutputLogLike << std::endl << "PARTON likelihood values (all matched MC) : {";
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      txtOutputLogLike << std::fixed << std::setprecision(15) << loglike_gen_[iWidth];
      if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
    }
    txtOutputLogLike << "};";
  }
  txtOutputLogLike << std::endl << std::endl;
}
