#include "../interface/LikelihoodMass.h"

const double LikelihoodMass::massArray_[] = {171., 171.05, 171.1, 171.15, 171.2, 171.25, 171.3, 171.35, 171.4, 171.45, 171.5, 171.55, 171.6, 171.65, 171.7, 171.75, 171.8, 171.85, 171.9, 171.95, 172., 172.05, 172.1, 172.15, 172.2, 172.25, 172.3, 172.35, 172.4, 172.45, 172.5, 172.55, 172.6, 172.65, 172.7, 172.75, 172.8, 172.85, 172.9, 172.95, 173., 173.05, 173.1, 173.15, 173.2, 173.25, 173.3, 173.35, 173.4, 173.45, 173.5, 173.55, 173.6, 173.65, 173.7, 173.75, 173.8, 173.85, 173.9, 173.95, 174.};

//const std::string LikelihoodMass::listCats_[] = {"CM", "WM", "UM"};
const std::string LikelihoodMass::listCats_[] = {"CM", "UM"};

const int LikelihoodMass::nMasses_ = sizeof(massArray_)/sizeof(massArray_[0]);
const int LikelihoodMass::nCats_ = sizeof(listCats_)/sizeof(listCats_[0]);

std::string LikelihoodMass::stringMassArray_[nMasses_] = {""};
std::string LikelihoodMass::stringSuffix_[nCats_] = {""};

double LikelihoodMass::loglike_[nMasses_] = {0.};
double LikelihoodMass::loglike_data_[nMasses_] = {0.};
double LikelihoodMass::loglike_per_evt_[nMasses_] = {0.};
double LikelihoodMass::loglike_CM_[nMasses_] = {0.};
double LikelihoodMass::loglike_CM_per_evt_[nMasses_] = {0.};
double LikelihoodMass::loglike_temp_[nMasses_] = {0.};
double LikelihoodMass::loglike_temp_per_evt_[nMasses_] = {0.};
double LikelihoodMass::loglike_gen_[nMasses_] = {0.};
double LikelihoodMass::loglike_gen_per_evt_[nMasses_] = {0.};

double LikelihoodMass::loglike_pull_[nMasses_][1000] = {{0.}};
double LikelihoodMass::loglike_pull_single_[nMasses_] = {0.};

const double LikelihoodMass::calCurvePar_[2] = {0., 1.};  // at the moment no output calibration
const double LikelihoodMass::calCurveParUnc_[2] = {0., 0.};  // at the moment no output calibration

double LikelihoodMass::nEventsCMFractions_[nMasses_][25] = {{0.}};
double LikelihoodMass::nEventsWMFractions_[nMasses_][25] = {{0.}};
double LikelihoodMass::nEventsUMFractions_[nMasses_][25] = {{0.}};


int LikelihoodMass::LocMinArray(int n, double* array)
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

void LikelihoodMass::ClearArray(int size, int* array)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = 0;
  }
}

void LikelihoodMass::ClearArray(int size, double* array)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = 0.;
  }
}

void LikelihoodMass::ClearArray2D(int size, double (*array)[3])
{
  for (int i = 0; i < size; i++)
  {
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      array[i][iCat] = 0.;
    }
  }
}

void LikelihoodMass::MakeTable(double* array, int n, double min, double max)
{
  double dist = (max - min)/((double)(n-1));
  for (int i = 0; i < n; i++)
  {
    array[i] = min + i * dist;
  }
}

LikelihoodMass::LikelihoodMass(double min, double max, std::string outputDirName, std::string date, bool useHadTopOnly, bool makeHistograms, bool verbose):
verbose_(verbose), rewHadOnly_(useHadTopOnly), outputDirName_(outputDirName), dirNameTGraphTxt_("OutputTxt/"), dirNameNEvents_("OutputNEvents/"), dirNameLLTxt_("OutputLikelihood/"+date+"/"), dirNamePull_("PseudoExp/"), inputFileName_(""), suffix_(""), histoName_(""), minRedMass_(min), maxRedMass_(max), histo_(), histoSm_(), histoTotal_(), graph_(), vecBinCentres_(), vecBinContents_(), calledLLCalculation_(false), calledCMLLCalculation_(false), calledGenLLCalculation_(false), vecMassFromFile_(), vecLLValsFromFile_()
{
  tls_ = new HelperTools();
  rew_ = new EventReweighting(false);  // no correction for number of events
  
  rangeRedMass_ = tls_->DotReplace(minRedMass_)+"To"+tls_->DotReplace(maxRedMass_);
  dirNameNEvents_ += rangeRedMass_+"/";
  
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    stringMassArray_[iMass] = tls_->DotReplace(massArray_[iMass]);
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

LikelihoodMass::~LikelihoodMass()
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

void LikelihoodMass::ClearLikelihoods()
{
  for (int i = 0; i < nMasses_; i++)
  {
    loglike_[i] = 0.;
    loglike_data_[i] = 0.;
    loglike_per_evt_[i] = 0.;
    loglike_CM_[i] = 0.;
    loglike_temp_[i] = 0.;
    loglike_gen_[i] = 0.;
  }
}

std::vector<double> LikelihoodMass::GetMasses()
{
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    vecMasses_.push_back(massArray_[iMass]);
  }
  return vecMasses_;
}

void LikelihoodMass::BookHistograms()
{
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    thisMass_ = stringMassArray_[iMass];
    
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      histo_[("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_60b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_60b").c_str(),("Reduced top mass for mass "+thisMass_+", "+listCats_[iCat]+"; m_{r}").c_str(), 60, 0.5, 2.0);
      histo_[("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_75b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_75b").c_str(),("Reduced top mass for mass "+thisMass_+", "+listCats_[iCat]+"; m_{r}").c_str(), 75, 0.5, 2.0);
      histo_[("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_90b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_90b").c_str(),("Reduced top mass for mass "+thisMass_+", "+listCats_[iCat]+"; m_{r}").c_str(), 90, 0.5, 2.0);  // 90 bins in range [0.5, 2.0]
      //histo_[("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_90b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_90b").c_str(),("Reduced top mass for mass "+thisMass_+", "+listCats_[iCat]+"; m_{r}").c_str(), 125, 0., 2.5);  // 90 bins in range [0.5, 2.0] --> red mlb
      //histo_[("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_90b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_90b").c_str(),("Reduced top mass for mass "+thisMass_+", "+listCats_[iCat]+"; m_{3/2}").c_str(), 50, 1., 3.5);  /// ONLY FOR NEW VAR !!!
//      histo_[("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_100b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_mass"+thisMass_+"_100b").c_str(),("Reduced top mass for mass "+thisMass_+", "+listCats_[iCat]+"; m_{r}").c_str(), 100, 0.5, 2.0);
    }
  }
}

void LikelihoodMass::FillHistograms(double redMass, double relativeSF, double hadTopMassForMassSF, double lepTopMassForMassSF, bool isTTbar, bool isData, std::string catSuffix)
{
  if ( ! isData )
  {
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      thisMass_ = stringMassArray_[iMass];
      if (isTTbar)
      {
        if (rewHadOnly_) thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, massArray_[iMass], 1., 1.);
        else thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, massArray_[iMass], 1., 1.) * rew_->BEventWeightCalculatorNonRel(lepTopMassForMassSF, 172.5, massArray_[iMass], 1., 1.);
      }
      else thisMassSF_ = 1.;
      
      histo_[("Red_top_mass"+catSuffix+"_mass"+thisMass_+"_60b").c_str()]->Fill(redMass, relativeSF*thisMassSF_);
      histo_[("Red_top_mass"+catSuffix+"_mass"+thisMass_+"_75b").c_str()]->Fill(redMass, relativeSF*thisMassSF_);
      histo_[("Red_top_mass"+catSuffix+"_mass"+thisMass_+"_90b").c_str()]->Fill(redMass, relativeSF*thisMassSF_);
//      histo_[("Red_top_mass"+catSuffix+"_mass"+thisMass_+"_100b").c_str()]->Fill(redMass, relativeSF*thisMassSF_);
    }
  }
}

void LikelihoodMass::WriteHistograms(std::string histoFileName)
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

void LikelihoodMass::GetHistogram(int iCat)
{
  int binMin = 0;
  int binMax = 1000;

  /// Get histo to smooth
  histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat];
  histoSm_[histoName_] = (TH1D*) histo_["Red_top_mass_"+histoName_+"_90b"]->Clone(histoName_.c_str());
  if ( iCat != 0 ) histoSm_[histoName_]->Smooth(3);
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

void LikelihoodMass::MakeGraph(int iCat, int nPoints, double* centres, double* contents, std::string name, bool drawGraph)
{
  histoName_ = name+stringSuffix_[iCat];
  graph_[histoName_] = new TGraph(nPoints, centres, contents);
  graph_[histoName_]->SetName(("g"+histoName_).c_str());
  graph_[histoName_]->SetTitle(histoName_.c_str());
  graph_[histoName_]->Write();
  if (drawGraph) this->DrawGraph(histoSm_[histoName_], graph_[histoName_], "Graph_Red_top_mass_"+histoName_);
}

void LikelihoodMass::MakeGraphSmooth(int iCat, int nPoints, double* centres, double* contents, std::string name, bool drawGraph)
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

void LikelihoodMass::ConstructTGraphsFromHisto(std::string tGraphFileName, std::vector<std::string> datasetNames, std::vector<int> includeDataset)
{
  TGraph2D *gLL2D_ = new TGraph2D();
  
  /// Define vars for likelihood calculation
  const int nEval = 50;
  double evalPoints[nEval], outputValues[nEval], outputValuesTemp[nEval], likelihoodValues[nEval], likelihoodValuesTemp[nEval];
  this->ClearArray(nEval, evalPoints);
  this->MakeTable(evalPoints, nEval, minRedMass_, maxRedMass_);
  
  /// Determine fractions based on number of events within reduced top mass range
  //  (is not dependent on number of bins, so the same for all masses)
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
  
  
  
  /// Loop over masses
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      stringSuffix_[iCat] = "mass"+stringMassArray_[iMass];
      this->GetHistogram(iCat);
    }
    
    /// Make arrays as input for TGraph
    const int nPoints = (vecBinCentres_[listCats_[0]+"_"+stringSuffix_[0]]).size();
    double binCentreArray[nPoints], binContentArray[nCats_][nPoints];
    this->ClearArray(nPoints, binCentreArray);
    
    /// Make likelihood functions
    this->ClearArray(nEval, outputValues);
    this->ClearArray(nEval, outputValuesTemp);
    
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      this->ClearArray(nPoints, binContentArray[iCat]);
      for (int i = 0; i < nPoints; i++)
      {
        if ( iCat == 0 ) binCentreArray[i] = (vecBinCentres_[listCats_[iCat]+"_"+stringSuffix_[iCat]]).at(i);
        binContentArray[iCat][i] = (vecBinContents_[listCats_[iCat]+"_"+stringSuffix_[iCat]]).at(i);
      }
      
      this->WriteFuncOutput(nPoints, binCentreArray, binContentArray[iCat], listCats_[iCat]+"_"+stringSuffix_[iCat]);
      this->MakeGraph(iCat, nPoints, binCentreArray, binContentArray[iCat], listCats_[iCat]+"_", true);
      
      histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat];
      histoNameSm_ = listCats_[iCat]+"_Sm_"+stringSuffix_[iCat];
      
      histoSm_[histoName_]->Scale(fracCats[iCat]);
      
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
        histoTotal_[stringSuffix_[0]]->SetTitle("#frac{n_{CM}}{n_{tot}} * f_{CM}(x|m_{t}) + #frac{n_{WM}}{n_{tot}} * f_{WM}(x) + #frac{n_{UM}}{n_{tot}} * f_{UM}(x)");
        
        this->MakeGraph(0, nEval, evalPoints, likelihoodValues, "likelihood_CM_", false);
        this->WriteOutput(nEval, massArray_[iMass], evalPoints, likelihoodValues, "CorrectMatchLikelihood_"+stringSuffix_[0], 1);
      }
      else histoTotal_[stringSuffix_[0]]->Add(histoSm_[histoName_]);
      
      if ( listCats_[iCat].find("WM") != std::string::npos && iCat != nCats_-1 )  // outputValues = CM + WM
      {
        this->MakeGraph(0, nEval, evalPoints, likelihoodValuesTemp, "likelihood_CMWM_", false);
        this->WriteOutput(nEval, massArray_[iMass], evalPoints, likelihoodValuesTemp, "MatchLikelihood_"+stringSuffix_[0], 1);
      }
      else if ( iCat == nCats_-1 ) // outputValues = CM + WM + UM
      {
        histoTotal_[stringSuffix_[0]]->Write();
        if (verbose_) std::cout << "LikelihoodMass::ConstructTGraphs: The integral of the weighted probability histogram is " << histoTotal_[stringSuffix_[0]]->Integral(histoTotal_[stringSuffix_[0]]->FindBin(minRedMass_), histoTotal_[stringSuffix_[0]]->FindBin(maxRedMass_)+1) << std::endl;
        this->MakeGraph(0, nEval, evalPoints, outputValues, "TotalProbability_", false);
        this->DrawGraph(histoTotal_[stringSuffix_[0]], graph_["TotalProbability_"+stringSuffix_[0]], "Graph_totalProbability_"+stringSuffix_[0]);
//         if (verbose_)
//         {
//           graph_["TotalProbability_"+stringSuffix_[0]+"_test"] = graph_["TotalProbability_"+stringSuffix_[0]];
//           graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->SetPoint(graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->GetN(), maxRedMass_, 0.);
//           graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->SetPoint(graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->GetN(), minRedMass_, 0.);
//           std::cout << "LikelihoodMass::ConstructTGraphs: The integral of the weighted probability graph is " << graph_["TotalProbability_"+stringSuffix_[0]+"_test"]->Integral() << std::endl;
//         }
        
        this->MakeGraph(0, nEval, evalPoints, likelihoodValues, "likelihood_", false);
        this->WriteOutput(nEval, massArray_[iMass], evalPoints, likelihoodValues, stringSuffix_[0], 1);  // for TGraph
        this->WriteOutput(nEval, massArray_[iMass], evalPoints, likelihoodValues, stringSuffix_[0], 2);  // for TGraph2D
      }
      
    }  // end cats
    
    /// Fill TGraph2D
    for (int iCentre = 0; iCentre < nEval; iCentre++)
    {
      gLL2D_->SetPoint(gLL2D_->GetN(), evalPoints[iCentre], massArray_[iMass], likelihoodValues[iCentre]);
    }
    
  }  // end masses
  
  CombineOutput();
  
  DrawLikelihoods();
  gLL2D_->SetName("2D_likelihood");
  gLL2D_->SetTitle("2D_likelihood");
  gLL2D_->Write();
  DrawGraph(gLL2D_, "2D_likelihood");
  
  
  fileTGraphs_->Close();
  delete fileTGraphs_;
  delete gLL2D_;
}

bool LikelihoodMass::ConstructTGraphsFromFile()
{
  return this->ConstructTGraphsFromFile("");
}

bool LikelihoodMass::ConstructTGraphsFromFile(std::string name)
{
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    suffix_ = name+"mass"+stringMassArray_[iMass];
    inputFileName_ = outputDirName_+dirNameTGraphTxt_+"output_tgraph1d_"+suffix_+".txt";
    if ( ! tls_->fexists(inputFileName_.c_str()) )
    {
      std::cerr << "LikelihoodMass::ConstructTGraphs: File " << inputFileName_ << " not found!!" << std::endl;
      std::cerr << "                              Aborting the likelihood calculation..." << std::endl;
      continue;
    }
    graph_[suffix_] = new TGraph(inputFileName_.c_str());
  }
  std::cout << "LikelihoodMass::ConstructTGraphs: Constructed TGraphs for likelihood measurements (using " << nMasses_ << " masses)" << std::endl;
  return true;
}

bool LikelihoodMass::ReadInput(std::string name)
{
  std::string line;
  double thisCentre, thisContent;
  (vecBinCentres_[name]).clear();
  (vecBinContents_[name]).clear();
  
  inputFileName_ = outputDirName_+dirNameTGraphTxt_+"output_func_"+name+".txt";
  if ( ! tls_->fexists(inputFileName_.c_str()) )
  {
    std::cerr << "LikelihoodMass::ConstructTGraphs: File " << inputFileName_ << " not found!!" << std::endl;
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

bool LikelihoodMass::ConstructTGraphsFromFile(std::vector<std::string> datasetNames, std::vector<int> includeDataset)
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
  std::cout << "LikelihoodMass::ConstructTGraphs: Creating output file " << dirNameLLTxt_+"TGraphs.root" << std::endl;
  fileTGraphs_->cd();
  
  
  /// Make arrays as input for TGraph
  const int nPoints = (vecBinCentres_[listCats_[1]+"_"+stringSuffix_[1]]).size();
  double binCentreArray[nPoints], binContentArray[nCats_][nPoints];
  
  
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    this->ClearArray(nPoints, binCentreArray);
    
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      stringSuffix_[iCat] = "mass"+stringMassArray_[iCat];
      histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat];
      //histoNameSm_ = listCats_[iCat]+"_Sm_"+stringSuffix_[iCat];
      
      if (! this->ReadInput(histoName_)) return false;
      
      this->ClearArray(nPoints, binContentArray[iCat]);
      for (int i = 0; i < nPoints; i++)
      {
        if ( iCat == 0 ) binCentreArray[i] = (vecBinCentres_[listCats_[iCat]+"_"+stringSuffix_[iCat]]).at(i);
        binContentArray[iCat][i] = (vecBinContents_[histoName_]).at(i);
      }
      
      /// Make TGraphs
      this->MakeGraph(iCat, nPoints, binCentreArray, binContentArray[iCat], listCats_[iCat]+"_");
      
      /// Make likelihood functions
      this->ClearArray(nEval, outputValues);
      this->ClearArray(nEval, likelihoodValues);
      for (int i = 0; i < nEval; i++)
      {
        outputValues[i] += fracs[iCat] * graph_[histoName_]->Eval(evalPoints[i]);
        likelihoodValues[i] = -TMath::Log(outputValues[i]);
      }
    }  // end cats
    
    this->MakeGraph(0, nEval, evalPoints, likelihoodValues, "");
    
  }  // end masses
  
  std::cout << "LikelihoodMass::ConstructTGraphs: Constructed TGraphs for likelihood measurements (using " << nMasses_ << " masses)" << std::endl;
  return true;
}

void LikelihoodMass::CalculateLikelihood(double redMass, double relativeSF, bool isData)
{
  this->CalculateLikelihood(redMass, relativeSF, 1., 1., 1., false, isData);  // isTTbar = false ==> thisMassSF_ = 1.;
}

std::vector<double> LikelihoodMass::CalculateLikelihood(double redMass, double relativeSF, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  return this->CalculateLikelihood(redMass, relativeSF, hadTopMassForMassSF, lepTopMassForMassSF, inputMass, isTTbar, isData);
}

std::vector<double> LikelihoodMass::CalculateLikelihood(double redMass, double relativeSF, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData)
{
  if (! isData && ! calledLLCalculation_) calledLLCalculation_ = true;
  
  vecLogLike_.clear();
  
  if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (isTTbar)
    {
      if (rewHadOnly_) thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.);
      else thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.) * rew_->BEventWeightCalculatorNonRel(lepTopMassForMassSF, 172.5, inputMass, 1., 1.);
    }
    else thisMassSF_ = 1.;
    
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      histoName_ = "mass"+stringMassArray_[iMass];
      
      loglike_per_evt_[iMass] = graph_[histoName_]->Eval(redMass);
      vecLogLike_.push_back(loglike_per_evt_[iMass]);
      if (! isData) loglike_[iMass] += loglike_per_evt_[iMass]*relativeSF*thisMassSF_;
      else loglike_data_[iMass] += loglike_per_evt_[iMass];
    }
  }
  
  return vecLogLike_;
}

void LikelihoodMass::CalculateCMLikelihood(double redMass, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  this->CalculateCMLikelihood(redMass, scaleFactor, hadTopMassForMassSF, lepTopMassForMassSF, inputMass, isTTbar, isData);
}

void LikelihoodMass::CalculateCMLikelihood(double redMass, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "LikelihoodMass::Cannot calculate loglikelihood for matched events when running over data" << std::endl;
    std::cerr << "            Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledCMLLCalculation_) calledCMLLCalculation_ = true;
    
    if (isTTbar)
    {
      if (rewHadOnly_) thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.);
      else thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.) * rew_->BEventWeightCalculatorNonRel(lepTopMassForMassSF, 172.5, inputMass, 1., 1.);
    }
    else thisMassSF_ = 1.;
    
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      thisMass_ = stringMassArray_[iMass];
      
      loglike_CM_per_evt_[iMass] = graph_["CorrectMatchLikelihood_mass"+thisMass_]->Eval(redMass);
      loglike_CM_[iMass] += loglike_CM_per_evt_[iMass]*scaleFactor*thisMassSF_;
    }
  }
}

void LikelihoodMass::CalculateTempLikelihood(double redMass, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  this->CalculateTempLikelihood(redMass, scaleFactor, hadTopMassForMassSF, lepTopMassForMassSF, inputMass, isTTbar, isData);
}

void LikelihoodMass::CalculateTempLikelihood(double redMass, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "LikelihoodMass::Cannot calculate loglikelihood for correctly/wrondly matched events when running over data" << std::endl;
    std::cerr << "            Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledTempLLCalculation_) calledTempLLCalculation_ = true;
    
    if (isTTbar)
    {
      if (rewHadOnly_) thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.);
      else thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.) * rew_->BEventWeightCalculatorNonRel(lepTopMassForMassSF, 172.5, inputMass, 1., 1.);
    }
    else thisMassSF_ = 1.;
    
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      thisMass_ = stringMassArray_[iMass];
      
      loglike_temp_per_evt_[iMass] = graph_["MatchLikelihood_mass"+thisMass_]->Eval(redMass);
      loglike_temp_[iMass] += loglike_temp_per_evt_[iMass]*scaleFactor*thisMassSF_;
    }
  }
}

void LikelihoodMass::CalculateGenLikelihood(double redMass, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "LikelihoodMass::Cannot calculate loglikelihood for generated events when running over data" << std::endl;
    std::cerr << "            Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledGenLLCalculation_) calledGenLLCalculation_ = true;
    
    if (isTTbar)
    {
      if (rewHadOnly_) thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.);
      else thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.) * rew_->BEventWeightCalculatorNonRel(lepTopMassForMassSF, 172.5, inputMass, 1., 1.);
    }
    else thisMassSF_ = 1.;
    
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      thisMass_ = stringMassArray_[iMass];
      
      loglike_gen_per_evt_[iMass] = graph_["CorrectMatchLikelihood_mass"+thisMass_]->Eval(redMass);
      loglike_gen_[iMass] += loglike_gen_per_evt_[iMass]*thisMassSF_;
    }
  }
}

void LikelihoodMass::GetOutputMass(double inputWidth, double inputMass, bool writeToFile, bool makeNewFile)
{
  this->GetOutputMass(inputWidth, inputMass, "", writeToFile, makeNewFile);
}

void LikelihoodMass::GetOutputMass(double inputWidth, double inputMass, std::string type, bool writeToFile, bool makeNewFile)
{
  std::string loglikePlotName = "loglikelihood_vs_mass_";
  if (! type.empty() ) loglikePlotName += type+"_";
  loglikePlotName += "widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass);
  
  if ( type.find("CM") != std::string::npos )
    output_ = this->CalculateOutputMass(nMasses_, loglike_CM_, loglikePlotName, writeToFile, makeNewFile);
  else if ( type.find("match") != std::string::npos || type.find("Match") != std::string::npos )
    output_ = this->CalculateOutputMass(nMasses_, loglike_temp_, loglikePlotName, writeToFile, makeNewFile);
  else if ( type.find("gen") != std::string::npos || type.find("Gen") != std::string::npos )
    output_ = this->CalculateOutputMass(nMasses_, loglike_gen_, loglikePlotName, writeToFile, makeNewFile);
  else if ( type.find("data") != std::string::npos || type.find("Data") != std::string::npos )
    output_ = this->CalculateOutputMass(nMasses_, loglike_data_, loglikePlotName, writeToFile, makeNewFile);
  else
    output_ = this->CalculateOutputMass(nMasses_, loglike_, loglikePlotName, writeToFile, makeNewFile);
  
  std::cout << "For an input mass of " << inputMass << " the minimum can be found at " << output_.first << " and the uncertainty is " << output_.second << std::endl;
  
  std::string fileName = dirNameLLTxt_+"result_minimum_";
  if (! type.empty() ) fileName += type+"_";
  fileName += "widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass)+".txt";
  txtOutputLL_.open(fileName.c_str());
  if (! type.empty() ) txtOutputLL_ << std::setw(18) << std::left << type << "  ";
  else txtOutputLL_ << std::setw(18) << std::left << "nominal  ";
  txtOutputLL_ << std::setw(5) << std::left << std::setprecision(5) << inputWidth << "   " << std::setw(5) << std::left << std::setprecision(5) << inputMass << "   " << std::setw(20) << std::right << std::setprecision(20) << output_.first << "  " << std::setw(20) << std::right << std::setprecision(20) << output_.second << std::endl;
  txtOutputLL_.close();
}

void LikelihoodMass::GetOutputMass(std::string inputFileName, double inputWidth, double inputMass, bool writeToFile, bool makeNewFile)
{
  this->GetOutputMass(inputFileName, dirNameLLTxt_, inputWidth, inputMass, writeToFile, makeNewFile);
}

void LikelihoodMass::GetOutputMass(std::string inputFileName, std::string inputDir, double inputWidth, double inputMass, bool writeToFile, bool makeNewFile)
{
  if (verbose_) std::cout << "Using LogLikelihood values from file" << std::endl;
  std::string loglikePlotName = "loglikelihood_vs_mass_";
  //if (! type.empty() ) loglikePlotName += type+"_";
  loglikePlotName += "ff_widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass);
  
  output_ = this->CalculateOutputMass(inputFileName, inputDir, loglikePlotName, writeToFile, makeNewFile);
  
  std::cout << "For an input mass of " << inputMass << " the minimum can be found at " << output_.first << " and the uncertainty is " << output_.second << std::endl;
}

std::pair<double,double> LikelihoodMass::GetOutputMass(double inputWidth, double inputMass, int thisPsExp)
{
  std::string loglikePlotName = dirNamePull_+"loglikelihood_vs_mass_psExp_"+tls_->ConvertIntToString(thisPsExp)+"_widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass);
  
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    loglike_pull_single_[iMass] = loglike_pull_[iMass][thisPsExp];
  }
  
  return this->CalculateOutputMass(nMasses_, loglike_pull_single_, loglikePlotName, true, false);
}

std::pair<double,double> LikelihoodMass::CalculateOutputMass(std::string inputFileName, std::string inputDir, std::string plotName, bool writeToFile, bool makeNewFile)
{
  this->ReadLLValuesFromFile(inputFileName, inputDir);
  const int nn = vecMassFromFile_.size();
  if ( nn == 0 )
    return std::pair<double,double>(-1.,-1.);
  double arrMass[nn], arrLLVals[nn];
  for (int i = 0; i < nn; i++)
  {
    arrMass[i] = vecMassFromFile_.at(i);
    arrLLVals[i] = vecLLValsFromFile_.at(i);
  }
  return this->CalculateOutputMass(nn, arrMass, arrLLVals, plotName, writeToFile, makeNewFile);
}

std::pair<double,double> LikelihoodMass::CalculateOutputMass(int nn, double* LLvalues, std::string plotName, bool writeToFile, bool makeNewFile)
{
  return this->CalculateOutputMass(nn, (double*)massArray_, LLvalues, plotName, writeToFile, makeNewFile);
}

std::pair<double,double> LikelihoodMass::CalculateOutputMass(int nn, double* evalMasses, double* LLvalues, std::string plotName, bool writeToFile, bool makeNewFile)
{
  int locMin = LocMinArray(nn, LLvalues);
  //std::cout << "Index of minimum LL value is " << locMin << std::endl;
  if ( evalMasses[locMin] <= 171.5 )
  {
    double tempArray[nn-10];
    for (int i = 0; i < nn-10; i++) { tempArray[i] = LLvalues[i+10]; }
    int tempMin = LocMinArray(nn-10, tempArray); //std::cout << tempMin << "  " << evalMasses[tempMin] << std::endl;
    if ( LLvalues[locMin+2] > tempArray[tempMin+2] ) locMin = tempMin+10;
  }
  else if ( evalMasses[locMin] > 173.5 )
  {
    double tempArray[nn-10];
    for (int i = 0; i < nn-10; i++) { tempArray[i] = LLvalues[i]; }
    int tempMin = LocMinArray(nn-10, tempArray); //std::cout << tempMin << "  " << evalMasses[tempMin] << std::endl;
    if ( LLvalues[locMin+2] > tempArray[tempMin+2] ) locMin = tempMin;
  }
  double centreVal = evalMasses[locMin];
  
  double temp;
  TGraph *g = new TGraph();
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    temp = LLvalues[iMass] - LLvalues[locMin];
    g->SetPoint(g->GetN(), massArray_[iMass], temp);
  }
  
  double interval = 0.15;
//   if ( centreVal <= interval ) interval = centreVal - 0.1;
//   if ( centreVal > 3.8 ) interval = 0.8;
  double fitmax = centreVal + interval;
  double fitmin = centreVal - interval;
//   //if ( centreVal < 1.6 ) fitmin += 0.1;
//   if ( centreVal > 0.15 && fitmin < 0.15 ) fitmin = 0.15;
//   if ( centreVal > 0.2 && fitmin < 0.2 ) fitmin = 0.2;
//   if ( centreVal > 0.35 && fitmin < 0.3 ) fitmin = 0.3;
//   //if ( centreVal > 0.4 && fitmin < 0.3 ) fitmin = 0.3;
//   //if ( centreVal > 0.5 && fitmin < 0.3 ) fitmin = 0.3;
//   if ( centreVal > 0.55 && fitmin < 0.4 ) fitmin = 0.4;
//   if ( centreVal > 0.65 && fitmin < 0.45 ) fitmin = 0.45;
//   if ( centreVal > 0.75 && fitmin < 0.5 ) fitmin = 0.5;
//   if ( centreVal > 1.1 && fitmin < 0.8 ) fitmin = 0.8;
//   
//   if ( centreVal < 0.35 && fitmax > 0.4 ) fitmax = 0.4;
//   if ( centreVal > 0.35 && centreVal < 1.2 ) fitmax = centreVal + (centreVal - fitmin);
  
  
  if (verbose_) std::cout << "LikelihoodMass::CalculateOutputMass: Look for minimum in interval [" << fitmin << "," << fitmax << "] around " << centreVal << std::endl;
  
  //TF1 *parabola = new TF1("parabola", "pol2", fitmin, fitmax);
  TF1 *parabola = new TF1("parabola", "[0]*(x - [1])*(x - [1]) + [2]", fitmin, fitmax);
  parabola->SetParameter(0, 10.);
  parabola->SetParameter(1, centreVal);
  
  g->Fit(parabola,"RWME");
  
  double outputMass = parabola->GetMinimumX(fitmin, fitmax);
  double lowerSigma = parabola->GetX(parabola->Eval(outputMass) + 0.5, fitmin-interval, outputMass);
  double upperSigma = parabola->GetX(parabola->Eval(outputMass) + 0.5, outputMass, fitmax+interval);
  double sigma = (upperSigma - lowerSigma)/2.;
  if ( lowerSigma <= fitmin && upperSigma >= fitmax )
    std::cerr << "LikelihoodMass::CalculateOutputMass: ERROR: Uncertainty calculation limited by fit boundaries. Do not trust..." << std::endl;
  else if ( lowerSigma <= fitmin ) sigma = upperSigma - outputMass;
  else if ( upperSigma >= fitmax ) sigma = outputMass - lowerSigma;
  
  std::cout << "Minimum -log(like) value is " << parabola->Eval(outputMass) << std::endl;
  std::cout << "===> Minimum at " << outputMass << " +- " << sigma << std::endl;
  
  double LLreduced[nn];
  for (int i = 0; i < nn; i++)
    LLreduced[i] = LLvalues[i] - LLvalues[locMin] - parabola->Eval(outputMass);
  TGraph *g2 = new TGraph(nn, evalMasses, LLreduced);
  //parabola->FixParameter(1, outputMass);
  g2->Fit(parabola,"RWME");
  
  if (makeNewFile)
  {
    filePlots_ = new TFile((dirNameLLTxt_+"File_"+plotName+".root").c_str(), "RECREATE");
    filePlots_->cd();
  }
  
  this->DrawOutputLogLikelihood(g2, parabola, 0, 15, 1200., plotName+"_full", writeToFile);
  this->DrawOutputLogLikelihood(g2, parabola, outputMass-2., outputMass+2., 1.5*g2->Eval(outputMass+2.), plotName, writeToFile);
  this->DrawOutputLogLikelihood(g2, parabola, outputMass-0.5, outputMass+0.5, std::max(g2->Eval(outputMass-3.*sigma),g2->Eval(outputMass+0.5)), plotName+"_zoom", writeToFile);
  
  if (makeNewFile)
  {
    filePlots_->Close();
    delete filePlots_;
  }
  
  delete g2;
  delete g;
  delete parabola;
  
  return std::pair<double,double>(outputMass,sigma);
}

std::pair<double,double> LikelihoodMass::ApplyCalibrationCurve(double thisOutputMass, double thisOutputMassSigma)
{
  /// return 'thisInputMass' and uncertainty
  //  thisOutputMass = Par_[0] + Par_[1] * thisInputMass

  double thisInputMass = (thisOutputMass - calCurvePar_[0])/calCurvePar_[1];
  double thisInputMassSigma = TMath::Sqrt( thisOutputMassSigma*thisOutputMassSigma + calCurveParUnc_[0]*calCurveParUnc_[0] + calCurveParUnc_[1]*calCurveParUnc_[1]*thisInputMass*thisInputMass )/calCurvePar_[1];
  
  return std::pair<double,double>(thisInputMass,thisInputMassSigma);
  
}

int LikelihoodMass::InitPull(int nPsExp)
{
  mkdir((dirNameLLTxt_+dirNamePull_).c_str(),0777);
  
  nPsExp_ = nPsExp;
  if ( nPsExp > 1000 )
  {
    std::cout << "LikelihoodMass::Pull: Warning: Only 1000 pseudo experiments will be performed" << std::endl;
    return 1000;
  }
  else if (verbose_) std::cout << "LikelihoodMass::Pull: Performing " << nPsExp << " pseudo experiments" << std::endl;
  
  return nPsExp;
}

void LikelihoodMass::AddPsExp(int thisPsExp, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  this->AddPsExp(thisPsExp, scaleFactor, hadTopMassForMassSF, lepTopMassForMassSF, inputMass, isTTbar, isData);
}

void LikelihoodMass::AddPsExp(int thisPsExp, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, double inputMass, bool isTTbar, bool isData)
{
  if (isData) std::cerr << "LikelihoodMass::Pull: Will not use data for pseudo experiments..." << std::endl;
  else if (calledLLCalculation_)
  {
    if (isTTbar)
    {
      if (rewHadOnly_) thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.);
      else thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, inputMass, 1., 1.) * rew_->BEventWeightCalculatorNonRel(lepTopMassForMassSF, 172.5, inputMass, 1., 1.);
    }
    else thisMassSF_ = 1.;
    
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      loglike_pull_[iMass][thisPsExp] += loglike_per_evt_[iMass]*scaleFactor*thisMassSF_;
    }
  }
  else std::cerr << "LikelihoodMass::Pull: Did not calculate likelihoods! Cannot get input for pseudo experiments..." << thisPsExp << std::endl;
}

void LikelihoodMass::CalculatePull(double inputMass)
{
  this->CalculatePull(1., inputMass);
}

void LikelihoodMass::CalculatePull(double inputWidth, double inputMass)
{
  std::string fileName = dirNameLLTxt_+dirNamePull_+"PseudoExperiments_widthx"+tls_->DotReplace(inputWidth)+"_mass"+tls_->DotReplace(inputMass)+".root";
  TFile *filePull = new TFile(fileName.c_str(),"RECREATE");
  filePull->cd();
  TH1D *hPull = new TH1D("hPull", "; (m_{t,j} - <m_{t}>)/#sigma_{m_{t,j}}", 32, -4., 4.);
  TH1D *hAveMassOut = new TH1D("hAveMassOut", "; m_{t,out}", 60, 171., 174.);
  TH1D *hAveMassIn = new TH1D("hAveMassIn", "; m_{t,in}", 60, 171., 174.);
  TH1D *hUncAveMassOut = new TH1D("hUncAveMassOut", "; #sigma_{m_{t,out}}", 30, 0., 0.3);
  TH1D *hUncAveMassIn = new TH1D("hUncAveMassIn", "; #sigma_{m_{t,in}}", 30, 0., 0.3);
  
  /// Calculate output mass for all pseudo experiments
  //  Transform to input mass via calibration curve & calculate average input mass
  std::pair<double,double> thisOutputMass[nPsExp_], thisInputMass[nPsExp_];
  double aveInputMass = 0.;
  
  for (int iPsExp = 0; iPsExp < nPsExp_; iPsExp++)
  {
    // Clear
    thisOutputMass[iPsExp] = std::pair<double,double>(-1.,-1.);
    thisInputMass[iPsExp] = std::pair<double,double>(-1.,-1.);
    
    // Fill
    thisOutputMass[iPsExp] = this->GetOutputMass(inputWidth, inputMass, iPsExp);
    thisInputMass[iPsExp] = this->ApplyCalibrationCurve(thisOutputMass[iPsExp].first, thisOutputMass[iPsExp].second);
    if ( thisInputMass[iPsExp].first != -1. ) aveInputMass += thisInputMass[iPsExp].first;
    else std::cerr << "LikelihoodMass::CalculatePull: Input mass for pseudo experiment " << iPsExp << " is equal to -1! Ignoring this pseudo experiment... " << std::endl;
    hAveMassOut->Fill(thisOutputMass[iPsExp].first);
    hAveMassIn->Fill(thisInputMass[iPsExp].first);
    hUncAveMassOut->Fill(thisOutputMass[iPsExp].second);
    hUncAveMassIn->Fill(thisInputMass[iPsExp].second);
  }
  aveInputMass = aveInputMass/nPsExp_;
  //std::cout << "Output mass : " << thisOutputMass[0].first << " ;   sigma : " << thisOutputMass[0].second << std::endl;
  //std::cout << "Input mass  : " << thisInputMass[0].first << " ;   sigma : " << thisInputMass[0].second << std::endl;
  std::cout << "Average mass for pseudo experiments is " << aveInputMass << std::endl;
  
  this->WritePsExpOutput(thisOutputMass, thisInputMass, inputMass);
  gStyle->SetOptStat(1110);  // display entries, mean & RMS, but no title
  hAveMassOut->Write();
  hAveMassIn->Write();
  hUncAveMassOut->Write();
  hUncAveMassIn->Write();
  gStyle->SetOptStat(10);  // display entries, but no title, mean & RMS
  TF1 *wOut = new TF1("wOut", "gaus", hAveMassOut->GetMean()-0.3, hAveMassOut->GetMean()+0.3);
  hAveMassOut->Fit(wOut,"R");
  gStyle->SetOptFit(0111);
  hAveMassOut->Write();
  wOut->Write();
  TF1 *wIn = new TF1("wIn", "gaus", hAveMassIn->GetMean()-0.3, hAveMassIn->GetMean()+0.3);
  hAveMassIn->Fit(wIn,"R");
  gStyle->SetOptFit(0111);
  hAveMassIn->Write();
  wIn->Write();
  
  /// Fill histogram with (Gamma_j - <Gamma>)/sigma_j
  gStyle->SetOptStat(1110);  // display entries, mean & RMS, but no title
  double fillValue;
  
  for (int iPsExp = 0; iPsExp < nPsExp_; iPsExp++)
  {
    if ( thisInputMass[iPsExp].first != -1. )
    {
      fillValue = (thisInputMass[iPsExp].first - aveInputMass)/thisInputMass[iPsExp].second;
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

void LikelihoodMass::AddToFraction(int d, double scaleFactor, double hadTopMassForMassSF, double lepTopMassForMassSF, bool isTTbar, bool isCM, bool isWM, bool isUM)
{
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    /// Calculate massSF
    if (isTTbar)
    {
      if (rewHadOnly_) thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, massArray_[iMass], 1., 1.);
      else thisMassSF_ = rew_->BEventWeightCalculatorNonRel(hadTopMassForMassSF, 172.5, massArray_[iMass], 1., 1.) * rew_->BEventWeightCalculatorNonRel(lepTopMassForMassSF, 172.5, massArray_[iMass], 1., 1.);
    }
    else thisMassSF_ = 1.;
    
    if      (isCM) nEventsCMFractions_[iMass][d] += scaleFactor * thisMassSF_;
    else if (isWM) nEventsWMFractions_[iMass][d] += scaleFactor * thisMassSF_;
    else if (isUM) nEventsUMFractions_[iMass][d] += scaleFactor * thisMassSF_;
  }
}

void LikelihoodMass::CalculateFractions(std::vector<std::string> datasetNames)
{
  std::string fileName = "";
  int nDatasets = datasetNames.size();
  mkdir(dirNameNEvents_.c_str(),0777);
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    for (int d = 0; d < nDatasets; d++)
    {
      fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_mass"+stringMassArray_[iMass]+"_CM.txt";
      txtOutputFractions_.open(fileName.c_str());
      txtOutputFractions_ << nEventsCMFractions_[iMass][d] << std::endl;
      txtOutputFractions_.close();
      fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_mass"+stringMassArray_[iMass]+"_WM.txt";
      txtOutputFractions_.open(fileName.c_str());
      txtOutputFractions_ << nEventsWMFractions_[iMass][d] << std::endl;
      txtOutputFractions_.close();
      fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_mass"+stringMassArray_[iMass]+"_UM.txt";
      txtOutputFractions_.open(fileName.c_str());
      txtOutputFractions_ << nEventsUMFractions_[iMass][d] << std::endl;
      txtOutputFractions_.close();
    }  // end datasets
  }  // end masses
}

void LikelihoodMass::GetFractions(double *fractions, int nCats, std::vector<std::string> datasetNames, std::vector<int> includeDataset)
{
  LikelihoodMass::GetFractions(fractions, nCats, "172p5", datasetNames, includeDataset);
}

void LikelihoodMass::GetFractions(double *fractions, int nCats, std::string mass, std::vector<std::string> datasetNames, std::vector<int> includeDataset)
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
      fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_mass"+mass+"_"+listCats_[iCat]+".txt";
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
      std::cout << "Mass: " << mass;
      for (int iCat = 0; iCat < nCats_; iCat++)
        std::cout << "   # " << listCats_[iCat] << ": " << 100.*fractions[iCat] << "%";
      std::cout << std::endl;
    }
  }
}

void LikelihoodMass::DrawGraph(TH1D* h, TGraph* g, std::string name)
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

void LikelihoodMass::DrawGraph(TGraph2D* g, std::string name)
{
  TCanvas *c = new TCanvas(name.c_str(), name.c_str());
  c->cd();
  gStyle->SetPalette(1);
  g->SetMaxIter(500000);
  g->Draw("AP");
  g->SetName(name.c_str());
  g->SetTitle((name+"; m_{r}; M_{t}").c_str());
  g->GetXaxis()->SetTitleOffset(1.5);
  g->GetYaxis()->SetTitleOffset(1.5);
  c->SetName((name+"_surf1").c_str());
  g->Draw("surf1");
  c->Update();
  c->Write();
  c->SaveAs((outputDirName_+name+"_surf1.png").c_str());
  c->SetName((name+"_contz2").c_str());
  g->Draw("contz2");
  c->Update();
  c->Write();
  c->SaveAs((outputDirName_+name+"_contz2.png").c_str());
  c->SetName((name+"_colz").c_str());
  g->Draw("colz");
  c->Update();
  c->Write();
  c->SaveAs((outputDirName_+name+"_colz.png").c_str());
  c->SetName((name+"_tri1p0").c_str());
  g->SetMarkerSize(0.5);
  g->Draw("tri1 p0");
  c->Update();
  c->Write();
  c->SaveAs((outputDirName_+name+"_tri1p0.png").c_str());
  
  delete c;
}

void LikelihoodMass::DrawLikelihoods()
{
  Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+1, kCyan+1, kBlue+2, kMagenta, kViolet-5, kPink+10};
  TCanvas* c2 = new TCanvas("-log(likelihood)", "-log(likelihood)");
  c2->cd();
  graph_["likelihood_mass"+stringMassArray_[0]]->SetLineColor(colours[0]);
  graph_["likelihood_mass"+stringMassArray_[0]]->Draw();
  std::string temp = graph_["likelihood_mass"+stringMassArray_[0]]->GetTitle();
  graph_["likelihood_mass"+stringMassArray_[0]]->SetTitle("-LogLikelihood");
  c2->Update();
  for (int i = 1; i < nMasses_; i++)
  {
    graph_["likelihood_mass"+stringMassArray_[i]]->SetLineColor(colours[i%10]);
    graph_["likelihood_mass"+stringMassArray_[i]]->Draw("same");
    c2->Update();
  }
  c2->Write();
  c2->SaveAs((outputDirName_+"LogLikelihood_multi.png").c_str());
  
  graph_["likelihood_mass"+stringMassArray_[0]]->SetTitle(temp.c_str());
  
  /// CM likelihood
  c2->Clear();
  c2->SetName("-log(likelihood) CM");
  graph_["likelihood_CM_mass"+stringMassArray_[0]]->SetLineColor(colours[0]);
  temp = graph_["likelihood_CM_mass"+stringMassArray_[0]]->GetTitle();
  graph_["likelihood_CM_mass"+stringMassArray_[0]]->SetTitle("-LogLikelihood CM");
  graph_["likelihood_CM_mass"+stringMassArray_[0]]->Draw();
  c2->Update();
  for (int i = 1; i < nMasses_; i++)
  {
    graph_["likelihood_CM_mass"+stringMassArray_[i]]->SetLineColor(colours[i%10]);
    graph_["likelihood_CM_mass"+stringMassArray_[i]]->Draw("same");
    c2->Update();
  }
  c2->Write();
  c2->SaveAs((outputDirName_+"LogLikelihood_CM_multi.png").c_str());
  c2->Close();
  
  graph_["likelihood_CM_mass"+stringMassArray_[0]]->SetTitle(temp.c_str());
  
  delete c2;
}

void LikelihoodMass::DrawOutputLogLikelihood(TGraph* g, TF1* f, double minX, double maxX, double maxY, std::string name, bool writeToFile)
{
  std::string outputFileName = dirNameLLTxt_+name+".png";
  TCanvas* c1 = new TCanvas(name.c_str(), "-LogLikelihood vs. mass");
  c1->cd();
  g->Draw("AP");
  if ( minX < 0. ) minX = 0.;
  g->SetTitle("");
  g->GetXaxis()->SetRangeUser(minX,maxX);
  g->GetYaxis()->SetRangeUser(-0.05*maxY,maxY);
  g->GetXaxis()->SetTitle("M_{t}");
  g->GetYaxis()->SetTitle("\\Delta\\mathscr{L}(\\mathrm{M_{t}})");
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

void LikelihoodMass::ReadLLValuesFromFile(std::string inputFileName, std::string inputDir)
{
  if (! tls_->fexists((inputDir+inputFileName).c_str()) )
  {
    std::cerr << "WARNING: File " << inputDir+inputFileName << " does not exist." << std::endl;
    exit(1);
  }
  fileIn_.open((inputDir+inputFileName).c_str());
  if (verbose_) std::cout << "Opening " << inputDir+inputFileName << "..." << std::endl;
  
  vecMassFromFile_.clear();
  vecLLValsFromFile_.clear();
  
  std::string line;
  std::string var = "";
  double val = -1.;
  while( getline(fileIn_, line) )
  {
    std::istringstream iss(line);
    if ( line.find("Mass") != std::string::npos )
    {
      iss >> var;
      while ( iss >> val )
      {
        vecMassFromFile_.push_back(val);
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
  
  if ( vecMassFromFile_.size() != vecLLValsFromFile_.size() )
  {
    std::cerr << "LikelihoodMass::ReadLLValuesFromFile : Input file " << inputFileName << " cannot be read:" << std::endl;
    std::cerr << "                                   # masses != # likelihood values !!" << std::endl;
    vecMassFromFile_.clear();
    vecLLValsFromFile_.clear();
  }
}

void LikelihoodMass::WritePsExpOutput(std::pair<double,double> *outputMass, std::pair<double,double> *inputMass, double genMass)
{
  std::string outputTxtName = dirNameLLTxt_+dirNamePull_+"PseudoExperiments_output_mass"+tls_->DotReplace(genMass)+".txt";
  txtOutputPsExp_.open(outputTxtName.c_str());
  txtOutputPsExp_ << "iPsExp   output mass        unc           input mass        unc" << std::endl;
  for (int iPsExp = 0; iPsExp < nPsExp_; iPsExp++)
  {
    txtOutputPsExp_ << std::setw(4) << std::right << iPsExp << "   " << std::setw(8) << std::left << outputMass[iPsExp].first << "   " << std::setw(8) << std::left << outputMass[iPsExp].second << "   " << std::setw(8) << std::left << inputMass[iPsExp].first << std::setw(8) << std::left << inputMass[iPsExp].second << std::endl;
  }
  txtOutputPsExp_.close();
}

void LikelihoodMass::WriteFuncOutput(int nPoints, double *arrayCentre, double *arrayContent, std::string name)
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

void LikelihoodMass::WriteOutput(int nPoints, double mass, double *arrayCentre, double *arrayContent, std::string name, int dim)
{
  std::string dimension = "";
  if ( dim == 1 ) dimension = "1d";
  else if ( dim == 2 ) dimension = "2d";
  std::string outputTxtName = outputDirName_+dirNameTGraphTxt_+"output_tgraph"+dimension+"_"+name+".txt";
  txtOutput_.open(outputTxtName.c_str());
  
  for (int i = 0; i < nPoints; i++)
  {
    txtOutput_ << std::setw(8) << std::right << arrayCentre[i];
    if ( dim == 2 ) txtOutput_ << "   " << std::setw(5) << std::right << mass;
    txtOutput_ << "  " << std::setw(8) << std::right << arrayContent[i] << std::endl;
  }
  txtOutput_.close();
}

void LikelihoodMass::CombineOutput()
{
  std::ofstream combFile((outputDirName_+dirNameTGraphTxt_+"output_tgraph2d_total"+".txt").c_str());
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    std::string inFileName = outputDirName_+dirNameTGraphTxt_+"output_tgraph2d_mass"+stringMassArray_[iMass]+".txt";
    std::ifstream infile(inFileName.c_str());
    combFile << infile.rdbuf();
    infile.close();
  }
  combFile.close();
}

void LikelihoodMass::PrintLikelihoodOutput(std::string llFileName)
{
  std::ofstream txtOutputLogLike((dirNameLLTxt_+llFileName).c_str());
  txtOutputLogLike << "Masses:      ";
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    txtOutputLogLike << std::setw(5) << std::right << massArray_[iMass] << "  ";
  }
  txtOutputLogLike << std::endl << "LLikelihood: ";
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    txtOutputLogLike << std::setw(12) << std::right << loglike_[iMass] << "  ";
  }
  txtOutputLogLike << std::endl;
  txtOutputLogLike.close();
}

void LikelihoodMass::PrintLikelihoodOutputData(std::string llFileName)
{
  std::ofstream txtOutputLogLike((dirNameLLTxt_+llFileName).c_str());
  txtOutputLogLike << "Masses:      ";
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    txtOutputLogLike << std::setw(5) << std::right << massArray_[iMass] << "  ";
  }
  txtOutputLogLike << std::endl << "LLikelihood: ";
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    txtOutputLogLike << std::setw(17) << std::right << std::setprecision(15) << loglike_data_[iMass] << "  ";
  }
  txtOutputLogLike << std::endl;
  txtOutputLogLike.close();
}

void LikelihoodMass::PrintMtmLikelihoodOutput(std::string llFileName)
{
  std::ofstream txtOutputLogLike((dirNameLLTxt_+llFileName).c_str());
  txtOutputLogLike << "likelihood values (all MC) : {";
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    txtOutputLogLike << std::fixed << std::setprecision(15) << loglike_[iMass];
    if ( iMass != nMasses_-1 ) txtOutputLogLike << ", ";
  }
  txtOutputLogLike << "};" << std::endl;
//   txtOutputLogLike << "likelihood values (data-only) : {";
//   for (int iMass = 0; iMass < nMasses_; iMass++)
//   {
//     txtOutputLogLike << std::fixed << std::setprecision(15) << loglike_data_[iMass][1];
//     if ( iMass != nMasses_-1 ) txtOutputLogLike << ", ";
//   }
//   txtOutputLogLike << "};" << std::endl;
  if (calledCMLLCalculation_)
  {
    txtOutputLogLike << "CM likelihood values : {";
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      txtOutputLogLike << std::fixed << std::setprecision(15) << loglike_CM_[iMass];
      if ( iMass != nMasses_-1 ) txtOutputLogLike << ", ";
    }
    txtOutputLogLike << "};";
  }
  txtOutputLogLike << std::endl << "masses: ";
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    txtOutputLogLike << massArray_[iMass];
    if ( iMass != nMasses_-1 ) txtOutputLogLike << ", ";
  }
  if (calledGenLLCalculation_)
  {
    txtOutputLogLike << std::endl << "PARTON likelihood values (all matched MC) : {";
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      txtOutputLogLike << std::fixed << std::setprecision(15) << loglike_gen_[iMass];
      if ( iMass != nMasses_-1 ) txtOutputLogLike << ", ";
    }
    txtOutputLogLike << "};";
  }
  txtOutputLogLike << std::endl << std::endl;
}
