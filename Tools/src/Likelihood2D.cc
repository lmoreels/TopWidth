#include "../interface/Likelihood2D.h"

const double Likelihood2D::widthArray_[] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5., 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.};

//const double Likelihood2D::massArray_[] = {169.5, 170., 170.5, 171., 171.25, 171.5, 171.75, 172., 172.1, 172.2, 172.3, 172.4, 172.5, 172.6, 172.7, 172.8, 172.9, 173., 173.25, 173.5, 173.75, 174., 174.5, 175., 175.5};
const double Likelihood2D::massArray_[] = {171., 171.05, 171.1, 171.15, 171.2, 171.25, 171.3, 171.35, 171.4, 171.45, 171.5, 171.55, 171.6, 171.65, 171.7, 171.75, 171.8, 171.85, 171.9, 171.95, 172., 172.05, 172.1, 172.15, 172.2, 172.25, 172.3, 172.35, 172.4, 172.45, 172.5, 172.55, 172.6, 172.65, 172.7, 172.75, 172.8, 172.85, 172.9, 172.95, 173., 173.05, 173.1, 173.15, 173.2, 173.25, 173.3, 173.35, 173.4, 173.45, 173.5, 173.55, 173.6, 173.65, 173.7, 173.75, 173.8, 173.85, 173.9, 173.95, 174.};

const std::string Likelihood2D::listCats_[] = {"CM", "WM", "UM"};

const int Likelihood2D::nWidths_ = sizeof(widthArray_)/sizeof(widthArray_[0]);
const int Likelihood2D::nCats_ = sizeof(listCats_)/sizeof(listCats_[0]);
//const int Likelihood2D::nMasses_ = 25;
const int Likelihood2D::nMasses_ = 61;

std::string Likelihood2D::stringWidthArray_[nWidths_] = {""};
std::string Likelihood2D::stringMassArray_[nMasses_] = {""};
std::string Likelihood2D::stringSuffix_[nCats_][nMasses_] = {{""}};

double Likelihood2D::loglike_[nWidths_][nMasses_] = {{0.}};
double Likelihood2D::loglike_data_[nWidths_][nMasses_] = {{0.}};
double Likelihood2D::loglike_per_evt_[nWidths_][nMasses_] = {{0.}};
double Likelihood2D::loglike_CM_[nWidths_][nMasses_] = {{0.}};
double Likelihood2D::loglike_CM_per_evt_[nWidths_][nMasses_] = {{0.}};
double Likelihood2D::loglike_temp_[nWidths_][nMasses_] = {{0.}};
double Likelihood2D::loglike_temp_per_evt_[nWidths_][nMasses_] = {{0.}};
double Likelihood2D::loglike_gen_[nWidths_][nMasses_] = {{0.}};
double Likelihood2D::loglike_gen_per_evt_[nWidths_][nMasses_] = {{0.}};

double Likelihood2D::loglike_pull_[nWidths_][nMasses_][1000] = {{{0.}}};
double Likelihood2D::loglike_pull_single_[nWidths_][nMasses_] = {{0.}};

//const double Likelihood2D::calCurvePar_[2] = {0., 1.};  // at the moment no output calibration
//const double Likelihood2D::calCurveParUnc_[2] = {0., 0.};  // at the moment no output calibration
const double Likelihood2D::calCurvePar_[2] = {0.0110942, 0.962622};
const double Likelihood2D::calCurveParUnc_[2] = {0.03277, 0.0131771};

double Likelihood2D::nEventsCMFractions_[nMasses_][25] = {{0.}};
double Likelihood2D::nEventsWMFractions_[nMasses_][25] = {{0.}};
double Likelihood2D::nEventsUMFractions_[nMasses_][25] = {{0.}};


int Likelihood2D::LocMinArray(int n, double* array)
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

std::pair<int,int> Likelihood2D::LocMinArray(int n, double (*array)[nMasses_])
{
  if ( n == 0 ) return std::pair<int,int>(-1, -1);
  
  std::pair<int,int> locmin = std::pair<int,int>(0,0);
  double min = array[0][0];
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < nMasses_; j++)
    {
      if ( array[i][j] < min )
      {
        min = array[i][j];
        locmin = std::pair<int,int>(i,j);
      }
    }
  }
  return locmin;
}

void Likelihood2D::ClearArray(int size, int* array)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = 0;
  }
}

void Likelihood2D::ClearArray(int size, double* array)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = 0.;
  }
}

void Likelihood2D::ClearArray2D(int size, double (*array)[nMasses_])
{
  for (int i = 0; i < size; i++)
  {
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      array[i][iCat] = 0.;
    }
  }
}

void Likelihood2D::MakeTable(double* array, int n, double min, double max)
{
  double dist = (max - min)/((double)(n-1));
  for (int i = 0; i < n; i++)
  {
    array[i] = min + i * dist;
  }
}

Double_t Likelihood2D::DDparabola(Double_t *x, Double_t *par)
{
  Double_t xx = x[0] - par[0];
  Double_t yy = x[1] - par[1];
  
  return par[2] + par[3]*xx*xx + par[4]*yy*yy;
}

Double_t Likelihood2D::Rosenbrock(Double_t *x, Double_t *par)
{
  // f(x,y) = (a - x)^2 + b(y - x^2)^2
  Double_t term1 = par[0]*x[0] - x[1]*x[1];
  Double_t term2 = par[1] - par[2]*x[1];
  
  return par[3]*term1*term1 + par[4]*term2*term2 + par[5]*x[0] + par[6];
}

Likelihood2D::Likelihood2D(double min, double max, std::string outputDirName, std::string date, bool useHadTopOnly, bool makeHistograms, bool verbose):
verbose_(verbose), rewHadOnly_(useHadTopOnly), outputDirName_(outputDirName), dirNameTGraphTxt_("OutputTxt/"), dirNameNEvents_("OutputNEvents/"), dirNameLLTxt_("OutputLikelihood/"+date+"/"), dirNamePull_("PseudoExp/"), inputFileName_(""), suffix_(""), histoName_(""), minRedMass_(min), maxRedMass_(max), histo_(), histoSm_(), histoTotal_(), graph_(), vecBinCentres_(), vecBinContents_(), calledLLCalculation_(false), calledCMLLCalculation_(false), calledGenLLCalculation_(false), vecWidthFromFile_(), vecLLValsFromFile_()
{
  tls_ = new HelperTools();
  rew_ = new EventReweighting(false);  // no correction for number of events
  
  rangeRedMass_ = tls_->DotReplace(minRedMass_)+"To"+tls_->DotReplace(maxRedMass_);
  dirNameNEvents_ += rangeRedMass_+"/";
  //thisMass_ = "172p5";
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    stringWidthArray_[iWidth] = tls_->DotReplace(widthArray_[iWidth]);
  }
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    stringMassArray_[iMass] = tls_->DotReplace(massArray_[iMass]);
  }
  
  ClearLikelihoods();
  
  if (makeHistograms)
  {
    mkdir(outputDirName_.c_str(),0777);
    mkdir((outputDirName_+dirNameTGraphTxt_).c_str(),0777);
    this->BookHistograms();
  }
  else mkdir(dirNameLLTxt_.c_str(),0777);
}

Likelihood2D::~Likelihood2D()
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

void Likelihood2D::SetMass(double mass)
{
  //thisMass_ = tls_->DotReplace(mass);
  std::cout << "Likelihood2D::SetMass: Function obsolete now" << std::endl;
}

void Likelihood2D::ClearLikelihoods()
{
  for (int i = 0; i < nWidths_; i++)
  {
    for (int j = 0; j < nMasses_; j++)
    {
      loglike_[i][j] = 0.;
      loglike_data_[i][j] = 0.;
      loglike_per_evt_[i][j] = 0.;
      loglike_CM_[i][j] = 0.;
      loglike_temp_[i][j] = 0.;
      loglike_gen_[i][j] = 0.;
    }
  }
}

std::vector<double> Likelihood2D::GetWidths()
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    vecWidths_.push_back(widthArray_[iWidth]);
  }
  return vecWidths_;
}

std::vector<double> Likelihood2D::GetMasses()
{
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    vecMasses_.push_back(massArray_[iMass]);
  }
  return vecMasses_;
}

void Likelihood2D::BookHistograms()
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    thisWidth_ = stringWidthArray_[iWidth];
    
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      thisMass_ = stringMassArray_[iMass];
      
      for (int iCat = 0; iCat < nCats_; iCat++)
      {
//        histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_mass"+thisMass_+"_75b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_mass"+thisMass_+"_75b").c_str(),("Reduced top mass for width "+thisWidth_+" and mass "+thisMass_+", "+listCats_[iCat]+"; m_{r}").c_str(), 75, 0.5, 2.0);
        histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_mass"+thisMass_+"_90b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_mass"+thisMass_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+" and mass "+thisMass_+", "+listCats_[iCat]+"; m_{r}").c_str(), 90, 0.5, 2.0);
//        histo_[("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_mass"+thisMass_+"_100b").c_str()] = new TH1D(("Red_top_mass_"+listCats_[iCat]+"_widthx"+thisWidth_+"_mass"+thisMass_+"_100b").c_str(),("Reduced top mass for width "+thisWidth_+" and mass "+thisMass_+", "+listCats_[iCat]+"; m_{r}").c_str(), 100, 0.5, 2.0);
      }
    }
  }
}

void Likelihood2D::FillHistograms(double redMass, double relativeSF, double hadTopMassForWidthSF, double lepTopMassForWidthSF, bool isTTbar, bool isData, std::string catSuffix)
{
  if ( ! isData )
  {
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      thisWidth_ = stringWidthArray_[iWidth];
      if (isTTbar)
      {
        if (rewHadOnly_) thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, widthArray_[iWidth]);
        else thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, widthArray_[iWidth]) * rew_->EventWeightCalculatorNonRel(lepTopMassForWidthSF, widthArray_[iWidth]);
      }
      else thisWidthSF_ = 1.;
      
      for (int iMass = 0; iMass < nMasses_; iMass++)
      {
        thisMass_ = stringMassArray_[iMass];
        
        if (isTTbar)
        {
          if (rewHadOnly_) thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, massArray_[iMass]);
          else thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, massArray_[iMass]) * rew_->MassEventWeightCalculatorNonRel(lepTopMassForWidthSF, massArray_[iMass]);
        }
        else thisMassSF_ = 1.;
        
//        histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_mass"+thisMass_+"_75b").c_str()]->Fill(redMass, relativeSF*thisWidthSF_*thisMassSF_);
        histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_mass"+thisMass_+"_90b").c_str()]->Fill(redMass, relativeSF*thisWidthSF_*thisMassSF_);
//        histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_mass"+thisMass_+"_100b").c_str()]->Fill(redMass, relativeSF*thisWidthSF_*thisMassSF_);
      }
    }
  }
}

void Likelihood2D::WriteHistograms(std::string histoFileName)
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

void Likelihood2D::GetHistogram(int iCat, int iMass)
{
  int binMin = 0;
  int binMax = 1000;

  /// Get histo to smooth
  histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat][iMass];
  histoSm_[histoName_] = (TH1D*) histo_["Red_top_mass_"+histoName_+"_90b"]->Clone(histoName_.c_str());
  if ( iCat != 0 )
  {
    histoSm_[histoName_]->Smooth(3);
    if ( fabs(massArray_[iMass] - 172.5) > 2 ) histoSm_[histoName_]->Smooth(2);
  }
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

void Likelihood2D::MakeGraph(int iCat, int iMass, int nPoints, double* centres, double* contents, std::string name, bool drawGraph)
{
  histoName_ = name+stringSuffix_[iCat][iMass];
  graph_[histoName_] = new TGraph(nPoints, centres, contents);
  graph_[histoName_]->SetName(("g"+histoName_).c_str());
  graph_[histoName_]->SetTitle(histoName_.c_str());
  graph_[histoName_]->Write();
  if (drawGraph) this->DrawGraph(histoSm_[histoName_], graph_[histoName_], "Graph_Red_top_mass_"+histoName_);
}

void Likelihood2D::MakeGraphSmooth(int iCat, int iMass, int nPoints, double* centres, double* contents, std::string name, bool drawGraph)
{
  this->MakeGraph(iCat, iMass, nPoints, centres, contents, name, drawGraph);
  
  histoName_ = name+stringSuffix_[iCat][iMass];
  histoNameSm_ = name+"Sm_"+stringSuffix_[iCat][iMass];
  TGraphSmooth* gs = new TGraphSmooth(histoName_.c_str());
  //graph_[histoNameSm_] = gs->SmoothSuper(graph_[histoName_],"",0,0);
  graph_[histoNameSm_] = gs->SmoothSuper(graph_[histoName_],"",3);
  graph_[histoNameSm_]->SetName(("g"+histoNameSm_).c_str());
  graph_[histoNameSm_]->SetTitle(histoNameSm_.c_str());
  graph_[histoNameSm_]->Write();
  if (drawGraph) this->DrawGraph(histoSm_[histoName_], graph_[histoNameSm_], "Graph_Red_top_mass_"+histoNameSm_);
}

void Likelihood2D::ConstructTGraphsFromHisto(std::string tGraphFileName, std::vector<std::string> datasetNames, std::vector<int> includeDataset)
{
  if ( gLL2D_ == NULL ) gLL2D_ = new TGraph2D();
  
  /// Define vars for likelihood calculation
  const int nEval = 50;
  double evalPoints[nEval], outputValues[nEval], outputValuesTemp[nEval], likelihoodValues[nEval], likelihoodValuesTemp[nEval];
  this->ClearArray(nEval, evalPoints);
  this->MakeTable(evalPoints, nEval, minRedMass_, maxRedMass_);
  
  /// Determine fractions based on number of events within reduced top mass range
  //  (is not dependent on number of bins, so the same for all widths)
  double fracCats[nMasses_][nCats_] = {0.}, fracCatsTemp[nMasses_][nCats_-1] = {0.};
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    //this->GetFractions(fracCats[iMass], nCats_, stringMassArray_[iMass], datasetNames, includeDataset);
    //this->GetFractions(fracCatsTemp[iMass], nCats_-1, stringMassArray_[iMass], datasetNames, includeDataset);
    this->GetFractions(fracCats[iMass], nCats_, datasetNames, includeDataset);
    this->GetFractions(fracCatsTemp[iMass], nCats_-1, datasetNames, includeDataset);
    //if (verbose_) std::cout << "Mass : " << setw(6) << massArray_[iMass] << "   # CM: " << fracCats[iMass][0] << "*100%   # WM: " << fracCats[iMass][1] << "*100%   # UM: " << fracCats[iMass][2] << "*100%  " << std::endl;
  }
  
  /// Make output file
  fileTGraphs_ = new TFile((outputDirName_+tGraphFileName).c_str(), "RECREATE");
  fileTGraphs_->cd();
  
  /// WM & UM distribution are independent of the width
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    for (int iCat = 1; iCat < nCats_; iCat++)
    {
      stringSuffix_[iCat][iMass] = "widthx1_mass"+stringMassArray_[iMass];
      this->GetHistogram(iCat, iMass);
    }
  }
  
  /// Make arrays as input for TGraph
  const int nPoints = (vecBinCentres_[listCats_[1]+"_"+stringSuffix_[1][6]]).size();
  double binCentreArray[nPoints], binContentArray[nCats_][nPoints];
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    this->ClearArray(nPoints, binCentreArray);
    for (int iCat = 0; iCat < nCats_; iCat++) ClearArray(nPoints, binContentArray[iCat]);
    
    for (int i = 0; i < nPoints; i++)
    {
      binCentreArray[i] = (vecBinCentres_[listCats_[1]+"_"+stringSuffix_[1][iMass]]).at(i);
      for (int iCat = 1; iCat < nCats_; iCat++)
      {
        binContentArray[iCat][i] = (vecBinContents_[listCats_[iCat]+"_"+stringSuffix_[iCat][iMass]]).at(i);
      }
    }
    
    for (int iCat = 1; iCat < nCats_; iCat++)
    {
      this->WriteFuncOutput(nPoints, binCentreArray, binContentArray[iCat], listCats_[iCat]+"_"+stringSuffix_[iCat][iMass]);
      this->MakeGraphSmooth(iCat, iMass, nPoints, binCentreArray, binContentArray[iCat], listCats_[iCat]+"_", true);
      histoSm_[listCats_[iCat]+"_"+stringSuffix_[iCat][iMass]]->Scale(fracCats[iMass][iCat]);
    }
  }
  
  /// Loop over widths
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      stringSuffix_[0][iMass] = "widthx"+stringWidthArray_[iWidth]+"_mass"+stringMassArray_[iMass];
      this->GetHistogram(0, iMass);
      
      this->ClearArray(nPoints, binContentArray[0]);
      for (int i = 0; i < nPoints; i++)
        binContentArray[0][i] = (vecBinContents_[listCats_[0]+"_"+stringSuffix_[0][iMass]]).at(i);
      
      this->WriteFuncOutput(nPoints, binCentreArray, binContentArray[0], listCats_[0]+"_"+stringSuffix_[0][iMass]);
      this->MakeGraphSmooth(0, iMass, nPoints, binCentreArray, binContentArray[0], listCats_[0]+"_", true);
      histoSm_[listCats_[0]+"_"+stringSuffix_[0][iMass]]->Scale(fracCats[iMass][0]);
      
      /// Make likelihood functions
      this->ClearArray(nEval, outputValues);
      this->ClearArray(nEval, outputValuesTemp);
      for (int iCat = 0; iCat < nCats_; iCat++)
      {
//         if ( iCat == 0 )
//           stringSuffix_[iCat][iMass] = "widthx"+stringWidthArray_[iWidth]+"_mass"+stringMassArray_[iMass];
//         else stringSuffix_[iCat][iMass] = "widthx1_mass"+stringMassArray_[iMass];
        
        histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat][iMass];
        histoNameSm_ = listCats_[iCat]+"_Sm_"+stringSuffix_[iCat][iMass];
        
        this->ClearArray(nEval, likelihoodValues);
        for (int i = 0; i < nEval; i++)
        {
          outputValues[i] += fracCats[iMass][iCat] * graph_[histoName_]->Eval(evalPoints[i]);
          if ( iCat == 0 ) likelihoodValues[i] = -TMath::Log(graph_[histoName_]->Eval(evalPoints[i]));
          else likelihoodValues[i] = -TMath::Log(outputValues[i]);
          if ( iCat < nCats_-1)
          {
            outputValuesTemp[i] += fracCatsTemp[iMass][iCat] * graph_[histoName_]->Eval(evalPoints[i]);
            likelihoodValuesTemp[i] = -TMath::Log(outputValuesTemp[i]);
          }
        }
        
        if ( iCat == 0 )  // outputValues = CM
        {
          histoTotal_[stringSuffix_[0][iMass]] = (TH1D*) histoSm_[histoName_]->Clone(("TotalProbability_"+stringSuffix_[0][iMass]).c_str());
          histoTotal_[stringSuffix_[0][iMass]]->SetTitle("#frac{n_{CM}}{n_{tot}} * f_{CM}(x|#Gamma) + #frac{n_{WM}}{n_{tot}} * f_{WM}(x) + #frac{n_{UM}}{n_{tot}} * f_{UM}(x)");

          this->MakeGraph(0, iMass, nEval, evalPoints, likelihoodValues, "likelihood_CM_", false);
          this->WriteOutput(nEval, iWidth, evalPoints, likelihoodValues, "CorrectMatchLikelihood_"+stringSuffix_[0][iMass], 1);
        }
        else if ( iCat == 1 )  // outputValues = CM + WM
        {
          histoTotal_[stringSuffix_[0][iMass]]->Add(histoSm_[histoName_]);

          this->MakeGraph(0, iMass, nEval, evalPoints, likelihoodValuesTemp, "likelihood_CMWM_", false);
          this->WriteOutput(nEval, iWidth, evalPoints, likelihoodValuesTemp, "MatchLikelihood_"+stringSuffix_[0][iMass], 1);
        }
        else  // outputValues = CM + WM + UM
        {
          histoTotal_[stringSuffix_[0][iMass]]->Add(histoSm_[histoName_]);
          histoTotal_[stringSuffix_[0][iMass]]->Write();
          if (verbose_) std::cout << "Likelihood2D::ConstructTGraphs: The integral of the weighted probability histogram is " << histoTotal_[stringSuffix_[0][iMass]]->Integral(histoTotal_[stringSuffix_[0][iMass]]->FindBin(minRedMass_), histoTotal_[stringSuffix_[0][iMass]]->FindBin(maxRedMass_)+1) << std::endl;
          this->MakeGraph(0, iMass, nEval, evalPoints, outputValues, "TotalProbability_", false);
          this->DrawGraph(histoTotal_[stringSuffix_[0][iMass]], graph_["TotalProbability_"+stringSuffix_[0][iMass]], "Graph_totalProbability_"+stringSuffix_[0][iMass]);
//           if (verbose_)
//           {
//             graph_["TotalProbability_"+stringSuffix_[0][iMass]+"_test"] = graph_["TotalProbability_"+stringSuffix_[0][iMass]];
//             graph_["TotalProbability_"+stringSuffix_[0][iMass]+"_test"]->SetPoint(graph_["TotalProbability_"+stringSuffix_[0][iMass]+"_test"]->GetN(), maxRedMass_, 0.);
//             graph_["TotalProbability_"+stringSuffix_[0][iMass]+"_test"]->SetPoint(graph_["TotalProbability_"+stringSuffix_[0][iMass]+"_test"]->GetN(), minRedMass_, 0.);
//             std::cout << "Likelihood2D::ConstructTGraphs: The integral of the weighted probability graph is " << graph_["TotalProbability_"+stringSuffix_[0][iMass]+"_test"]->Integral() << std::endl;
//           }
          
          this->MakeGraph(0, iMass, nEval, evalPoints, likelihoodValues, "likelihood_", false);
          this->WriteOutput(nEval, iWidth, evalPoints, likelihoodValues, stringSuffix_[0][iMass], 1);  // for TGraph
          this->WriteOutput(nEval, iWidth, evalPoints, likelihoodValues, stringSuffix_[0][iMass], 2);  // for TGraph2D
        }
        
      }  // end cats
      
      /// Fill TGraph2D
      if ( massArray_[iMass] == 172.5 )
      {
        for (int iCentre = 0; iCentre < nEval; iCentre++)
        {
          gLL2D_->SetPoint(gLL2D_->GetN(), evalPoints[iCentre], widthArray_[iWidth], likelihoodValues[iCentre]);
        }
      }
     
    }  // end masses
  }  // end widths
  
  //CombineOutput();
  
  DrawLikelihoods();
  gLL2D_->SetName("2D_likelihood");
  gLL2D_->SetTitle("2D_likelihood");
  gLL2D_->Write();
  DrawGraph(gLL2D_, "2D_likelihood");
  
  
  fileTGraphs_->Close();
  delete fileTGraphs_;
}

bool Likelihood2D::ConstructTGraphsFromFile()
{
  return this->ConstructTGraphsFromFile("");
}

bool Likelihood2D::ConstructTGraphsFromFile(std::string name)
{
  unsigned int tmp = 0;
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      tmp = 0;
      suffix_ = name+"widthx"+stringWidthArray_[iWidth]+"_mass"+stringMassArray_[iMass];
      inputFileName_ = outputDirName_+dirNameTGraphTxt_+"output_tgraph1d_"+suffix_+".txt";
      if ( ! tls_->fexists(inputFileName_.c_str()) )
      {
        std::cerr << "Likelihood2D::ConstructTGraphs: File " << inputFileName_ << " not found!!" << std::endl;
        //std::cerr << "                          Aborting the likelihood calculation..." << std::endl;
        tmp++;
        if ( tmp == nMasses_ ) return false;
        else continue;
      }
      graph_[suffix_] = new TGraph(inputFileName_.c_str());
    }
  }
  std::cout << "Likelihood2D::ConstructTGraphs: Constructed TGraphs for likelihood measurements (using " << nWidths_ << " widths)" << std::endl;
  return true;
}

bool Likelihood2D::ReadInput(std::string name)
{
  std::string line;
  double thisCentre, thisContent;
  (vecBinCentres_[name]).clear();
  (vecBinContents_[name]).clear();
  
  inputFileName_ = outputDirName_+dirNameTGraphTxt_+"output_func_"+name+".txt";
  if ( ! tls_->fexists(inputFileName_.c_str()) )
  {
    std::cerr << "Likelihood2D::ConstructTGraphs: File " << inputFileName_ << " not found!!" << std::endl;
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

bool Likelihood2D::ConstructTGraphsFromFile(std::vector<std::string> datasetNames, std::vector<int> includeDataset)
{
  const int nEval = 50;
  double evalPoints[nEval], outputValues[nEval], likelihoodValues[nEval];
  this->ClearArray(nEval, evalPoints);
  this->MakeTable(evalPoints, nEval, minRedMass_, maxRedMass_);
  
  /// Get fractions
  double fracs[nMasses_][nCats_];
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    this->ClearArray(nCats_, fracs[iMass]);
    //this->GetFractions(fracs[iMass], nCats_, stringMassArray_[iMass], datasetNames, includeDataset);
    this->GetFractions(fracs[iMass], nCats_, datasetNames, includeDataset);
  }
  
  fileTGraphs_ = new TFile((dirNameLLTxt_+"TGraphs.root").c_str(), "RECREATE");
  std::cout << "Likelihood2D::ConstructTGraphs: Creating output file " << dirNameLLTxt_+"TGraphs.root" << std::endl;
  fileTGraphs_->cd();
  
  /// WM & UM distribution are independent of the width
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    for (int iCat = 1; iCat < nCats_; iCat++)
    {
      stringSuffix_[iCat][iMass] = "widthx1_mass"+stringMassArray_[iMass];
      histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat][iMass];
      
      if (! this->ReadInput(histoName_)) return false;
    }
  }
  
  /// Make arrays as input for TGraph
  const int nPoints = (vecBinCentres_[listCats_[1]+"_"+stringSuffix_[1][6]]).size();
  double binCentreArray[nPoints], binContentArray[nCats_][nPoints];
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    this->ClearArray(nPoints, binCentreArray);
    for (int iCat = 0; iCat < nCats_; iCat++) this->ClearArray(nPoints, binContentArray[iCat]);
    
    for (int i = 0; i < nPoints; i++)
    {
      binCentreArray[i] = (vecBinCentres_[listCats_[1]+"_"+stringSuffix_[1][iMass]]).at(i);
      for (int iCat = 1; iCat < nCats_; iCat++)
      {
        binContentArray[iCat][i] = (vecBinContents_[listCats_[iCat]+"_"+stringSuffix_[iCat][iMass]]).at(i);
      }
    }
    
    /// Make TGraphs
    for (int iCat = 1; iCat < nCats_; iCat++)
      this->MakeGraphSmooth(iCat, iMass, nPoints, binCentreArray, binContentArray[iCat], listCats_[iCat]+"_");
  }
  
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      stringSuffix_[0][iMass] = "widthx"+stringWidthArray_[iWidth]+"_mass"+stringMassArray_[iMass];
      histoName_ = listCats_[0]+"_"+stringSuffix_[0][iMass];
      
      if (! this->ReadInput(histoName_)) return false;
      
      ClearArray(nPoints, binContentArray[0]);
      for (int i = 0; i < nPoints; i++)
        binContentArray[0][i] = (vecBinContents_[histoName_]).at(i);
      
      this->MakeGraphSmooth(0, iMass, nPoints, binCentreArray, binContentArray[0], listCats_[0]+"_");
      
      /// Make likelihood functions
      this->ClearArray(nEval, outputValues);
      this->ClearArray(nEval, likelihoodValues);
      for (int i = 0; i < nEval; i++)
      {
        for (int iCat = 0; iCat < nCats_; iCat++)
        {
          histoName_ = listCats_[iCat]+"_"+stringSuffix_[iCat][iMass];
          histoNameSm_ = listCats_[iCat]+"_Sm_"+stringSuffix_[iCat][iMass];
          outputValues[i] += fracs[iMass][iCat] * graph_[histoName_]->Eval(evalPoints[i]);
        }
        likelihoodValues[i] = -TMath::Log(outputValues[i]);
      }
      
      this->MakeGraph(0, iMass, nEval, evalPoints, likelihoodValues, "");
      
    }  // end masses
    
  }  // end widths
  
  std::cout << "Likelihood2D::ConstructTGraphs: Constructed TGraphs for likelihood measurements (using " << nWidths_ << " widths)" << std::endl;
  return true;
}

void Likelihood2D::CalculateLikelihood(double redMass, double relativeSF, bool isData)
{
  this->CalculateLikelihood(redMass, relativeSF, 1., 1., 1., 172.5, false, isData);  // isTTbar = false ==> thisWidthSF_ = 1.;
}

std::vector<double> Likelihood2D::CalculateLikelihood(double redMass, double relativeSF, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  if (! isData && ! calledLLCalculation_) calledLLCalculation_ = true;
  
  vecLogLike_.clear();
  
  if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (isTTbar)
    {
      if (rewHadOnly_)
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass);
      }
      else
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth) * rew_->EventWeightCalculatorNonRel(lepTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass) * rew_->MassEventWeightCalculatorNonRel(lepTopMassForWidthSF, inputMass);
      }
    }
    else
    {
      thisWidthSF_ = 1.;
      thisMassSF_ = 1.;
    }
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      for (int iMass = 0; iMass < nMasses_; iMass++)
      {
        histoName_ = "widthx"+stringWidthArray_[iWidth]+"_mass"+stringMassArray_[iMass];
        
        loglike_per_evt_[iWidth][iMass] = graph_[histoName_]->Eval(redMass);
        vecLogLike_.push_back(loglike_per_evt_[iWidth][iMass]);
        if (! isData) loglike_[iWidth][iMass] += loglike_per_evt_[iWidth][iMass]*relativeSF*thisWidthSF_*thisMassSF_;
        else loglike_data_[iWidth][iMass] += loglike_per_evt_[iWidth][iMass];
      }
    }
  }
  
  return vecLogLike_;
}

void Likelihood2D::CalculateCMLikelihood(double redMass, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "Likelihood2D::Cannot calculate loglikelihood for matched events when running over data" << std::endl;
    std::cerr << "              Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledCMLLCalculation_) calledCMLLCalculation_ = true;
    
    if (isTTbar)
    {
      if (rewHadOnly_)
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass);
      }
      else
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth) * rew_->EventWeightCalculatorNonRel(lepTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass) * rew_->MassEventWeightCalculatorNonRel(lepTopMassForWidthSF, inputMass);
      }
    }
    else
    {
      thisWidthSF_ = 1.;
      thisMassSF_ = 1.;
    }
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      for (int iMass = 0; iMass < nMasses_; iMass++)
      {
        histoName_ = "widthx"+stringWidthArray_[iWidth]+"_mass"+stringMassArray_[iMass];
        
        loglike_CM_per_evt_[iWidth][iMass] = graph_["CorrectMatchLikelihood_"+histoName_]->Eval(redMass);
        loglike_CM_[iWidth][iMass] += loglike_CM_per_evt_[iWidth][iMass]*scaleFactor*thisWidthSF_*thisMassSF_;
      }
    }
  }
}

void Likelihood2D::CalculateTempLikelihood(double redMass, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "Likelihood2D::Cannot calculate loglikelihood for correctly/wrondly matched events when running over data" << std::endl;
    std::cerr << "              Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledTempLLCalculation_) calledTempLLCalculation_ = true;
    
    if (isTTbar)
    {
      if (rewHadOnly_)
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass);
      }
      else
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth) * rew_->EventWeightCalculatorNonRel(lepTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass) * rew_->MassEventWeightCalculatorNonRel(lepTopMassForWidthSF, inputMass);
      }
    }
    else
    {
      thisWidthSF_ = 1.;
      thisMassSF_ = 1.;
    }
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      for (int iMass = 0; iMass < nMasses_; iMass++)
      {
        histoName_ = "widthx"+stringWidthArray_[iWidth]+"_mass"+stringMassArray_[iMass];
        
        loglike_temp_per_evt_[iWidth][iMass] = graph_["MatchLikelihood_"+histoName_]->Eval(redMass);
        loglike_temp_[iWidth][iMass] += loglike_temp_per_evt_[iWidth][iMass]*scaleFactor*thisWidthSF_*thisMassSF_;
      }
    }
  }
}

void Likelihood2D::CalculateGenLikelihood(double redMass, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "Likelihood2D::Cannot calculate loglikelihood for generated events when running over data" << std::endl;
    std::cerr << "              Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledGenLLCalculation_) calledGenLLCalculation_ = true;
    
    if (isTTbar)
    {
      if (rewHadOnly_)
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass);
      }
      else
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth) * rew_->EventWeightCalculatorNonRel(lepTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass) * rew_->MassEventWeightCalculatorNonRel(lepTopMassForWidthSF, inputMass);
      }
    }
    else
    {
      thisWidthSF_ = 1.;
      thisMassSF_ = 1.;
    }
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      for (int iMass = 0; iMass < nMasses_; iMass++)
      {
        histoName_ = "widthx"+stringWidthArray_[iWidth]+"_mass"+stringMassArray_[iMass];
        
        loglike_gen_per_evt_[iWidth][iMass] = graph_["CorrectMatchLikelihood_"+histoName_]->Eval(redMass);
        loglike_gen_[iWidth][iMass] += loglike_gen_per_evt_[iWidth][iMass]*thisWidthSF_*thisMassSF_;
      }
    }
  }
}

void Likelihood2D::GetOutputWidth(double inputWidth, bool writeToFile, bool makeNewFile)
{
  this->GetOutputWidth(inputWidth, "", writeToFile, makeNewFile);
}

void Likelihood2D::GetOutputWidth(double inputWidth, std::string type, bool writeToFile, bool makeNewFile)
{
  std::string loglikePlotName = "loglikelihood_vs_width_";
  if (! type.empty() ) loglikePlotName += type+"_";
  loglikePlotName += "widthx"+tls_->DotReplace(inputWidth);
  
  if ( type.find("CM") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_CM_, loglikePlotName, writeToFile, makeNewFile);
  else if ( type.find("match") != std::string::npos || type.find("Match") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_temp_, loglikePlotName, writeToFile, makeNewFile);
  else if ( type.find("gen") != std::string::npos || type.find("Gen") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_gen_, loglikePlotName, writeToFile, makeNewFile);
  else if ( type.find("data") != std::string::npos || type.find("Data") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_data_, loglikePlotName, writeToFile, makeNewFile);
  else
    output_ = this->CalculateOutputWidth(nWidths_, loglike_, loglikePlotName, writeToFile, makeNewFile);
  
  std::cout << "For an input width of " << inputWidth << " the minimum can be found at " << output_.first << " and the uncertainty is " << output_.second << std::endl;
  
  std::string fileName = dirNameLLTxt_+"result_minimum_";
  if (! type.empty() ) fileName += type+"_";
  fileName += "widthx"+tls_->DotReplace(inputWidth)+".txt";
  txtOutputLL_.open(fileName.c_str());
  if (! type.empty() ) txtOutputLL_ << std::setw(18) << std::left << type << "  ";
  else txtOutputLL_ << std::setw(18) << std::left << "nominal  ";
  txtOutputLL_ << std::setw(5) << std::left << std::setprecision(5) << inputWidth << "   " << std::setw(20) << std::right << std::setprecision(20) << output_.first << "  " << std::setw(20) << std::right << std::setprecision(20) << output_.second << std::endl;
  txtOutputLL_.close();
}

void Likelihood2D::GetOutputWidth(std::string inputFileName, double inputWidth, bool writeToFile, bool makeNewFile)
{
  this->GetOutputWidth(inputFileName, dirNameLLTxt_, inputWidth, writeToFile, makeNewFile);
}

void Likelihood2D::GetOutputWidth(std::string inputFileName, std::string inputDir, double inputWidth, bool writeToFile, bool makeNewFile)
{
  if (verbose_) std::cout << "Using LogLikelihood values from file" << std::endl;
  std::string loglikePlotName = "loglikelihood_vs_width_";
  //if (! type.empty() ) loglikePlotName += type+"_";
  loglikePlotName += "ff_widthx"+tls_->DotReplace(inputWidth);
  
  output_ = this->CalculateOutputWidth(inputFileName, inputDir, loglikePlotName, writeToFile, makeNewFile);
  
  std::cout << "For an input width of " << inputWidth << " the minimum can be found at " << output_.first << " and the uncertainty is " << output_.second << std::endl;
}

std::pair<double,double> Likelihood2D::GetOutputWidth(double inputWidth, int thisPsExp)
{
  std::string loglikePlotName = dirNamePull_+"loglikelihood_vs_width_psExp_"+tls_->ConvertIntToString(thisPsExp)+"_widthx"+tls_->DotReplace(inputWidth);
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      loglike_pull_single_[iWidth][iMass] = loglike_pull_[iWidth][iMass][thisPsExp];
    }
  }
  
  return this->CalculateOutputWidth(nWidths_, loglike_pull_single_, loglikePlotName, true, false);
}

std::pair<double,double> Likelihood2D::CalculateOutputWidth(std::string inputFileName, std::string inputDir, std::string plotName, bool writeToFile, bool makeNewFile)
{
  this->ReadLLValuesFromFile(inputFileName, inputDir);
  const int nn = vecWidthFromFile_.size();
  if ( nn == 0 )
    return std::pair<double,double>(-1.,-1.);
  double arrWidth[nn], arrLLVals[nn][nMasses_];
  for (int i = 0; i < nn; i++)
  {
    arrWidth[i] = vecWidthFromFile_.at(i);
    arrLLVals[i][6] = vecLLValsFromFile_.at(i);
  }
  return this->CalculateOutputWidth(nn, (double*)arrWidth, arrLLVals, plotName, writeToFile, makeNewFile);
}

std::pair<double,double> Likelihood2D::CalculateOutputWidth(int nn, double (*LLvalues)[nMasses_], std::string plotName, bool writeToFile, bool makeNewFile)
{
  return this->CalculateOutputWidth(nn, (double*)widthArray_, LLvalues, plotName, writeToFile, makeNewFile);
}

std::pair<double,double> Likelihood2D::CalculateOutputWidth(int nn, double* evalWidths, double (*LLvalues)[nMasses_], std::string plotName, bool writeToFile, bool makeNewFile)
{
  std::pair<int,int> locMin = this->LocMinArray(nn, LLvalues);
  std::cout << "Index of minimum LL value is (" << locMin.first << ", " << locMin.second << ")" << "           " << LLvalues[locMin.first][locMin.second] << std::endl;
  
  double temp;
  //TGraph *g = new TGraph(nn, evalWidths, LLvalues);
  TGraph2D *g = new TGraph2D();
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      temp = LLvalues[iWidth][iMass] - LLvalues[locMin.first][locMin.second];
      if ( temp < 50. )
      {
        g->SetPoint(g->GetN(), evalWidths[iWidth], massArray_[iMass], temp);
        if (verbose_) std::cout << evalWidths[iWidth] << "   " << massArray_[iMass] << "   " << std::setprecision(10) << temp << std::endl;
      }
    }
  }
  g->SetMaxIter(500000);
  
  double centreVal = evalWidths[locMin.first];
  if ( centreVal <= 0.3 )
  {
    double tempArray[nn-5][nMasses_];
    for (int i = 0; i < nn-5; i++)
    {
      for (int iMass = 0; iMass < nMasses_; iMass++)
      {
        tempArray[i][iMass] = LLvalues[i+5][iMass];
      }
    }
    std::pair<int,int> tempMin = LocMinArray(nn-5, tempArray); //std::cout << tempMin << "  " << evalWidths[tempMin] << std::endl;
    if ( LLvalues[locMin.first+2][locMin.second] > tempArray[tempMin.first+2][tempMin.second] )
      centreVal = evalWidths[tempMin.first+5];
  }
  
  
//   double interval = 0.4;
  double interval = 0.1;
  if ( centreVal <= interval ) interval = centreVal - 0.1;
//   if ( centreVal > 3.8 ) interval = 0.8;
  if ( centreVal > 2.8 ) interval = 0.4;
  double fitmax = centreVal + interval;
  double fitmin = centreVal - interval;
  //if ( centreVal < 1.6 ) fitmin += 0.1;
//   if ( centreVal > 0.15 && fitmin < 0.15 ) fitmin = 0.15;
//   if ( centreVal > 0.2 && fitmin < 0.2 ) fitmin = 0.2;
//   if ( centreVal > 0.35 && fitmin < 0.3 ) fitmin = 0.3;
//   //if ( centreVal > 0.4 && fitmin < 0.3 ) fitmin = 0.3;
//   //if ( centreVal > 0.5 && fitmin < 0.3 ) fitmin = 0.3;
//   if ( centreVal > 0.55 && fitmin < 0.4 ) fitmin = 0.4;
//   if ( centreVal > 0.75 && fitmin < 0.5 ) fitmin = 0.55;
//   if ( centreVal > 1.1 && fitmin < 0.8 ) fitmin = 0.8;
//   
//   if ( centreVal < 0.35 && fitmax > 0.4 ) fitmax = 0.4;
//   if ( centreVal > 0.35 && centreVal < 1.2 ) fitmax = centreVal + (centreVal - fitmin);
  
  double fitminY = massArray_[locMin.second] - 0.2;
  if ( fitminY < 169.5 ) fitminY = 169.5;
  double fitmaxY = massArray_[locMin.second] + 0.2;
  if ( fitmaxY > 175.5 ) fitmaxY = 175.5;
  
  if (verbose_) std::cout << "Likelihood2D::CalculateOutputWidth: Look for minimum around (" << centreVal << ", " << massArray_[locMin.second] << ")" << std::endl;
  if (verbose_) std::cout << "Boundaries are [" << fitmin << ", " << fitmax << "] for width and [" << fitminY << ", " << fitmaxY << "] for mass" << std::endl;
  
  //ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-3);
  //ROOT::Math::MinimizerOptions::SetDefaultErrorDef(2);  //TVirtualFitter::SetErrorDef(tol);
  //TVirtualFitter::SetPrecision(tol);
  
  if (makeNewFile)
  {
    filePlots_ = new TFile((dirNameLLTxt_+"File_"+plotName+".root").c_str(), "RECREATE");
    filePlots_->cd();
  }
  
  //TF1 *parabola = new TF1("parabola", "pol2", fitmin, fitmax);
  TF2 *parabola = new TF2("parabola", "[0] + [1]*x + [2]*x*x + [3]*y + [4]*y*y + [5]*x*y*y + [6]*y*y*y*y", fitmin, fitmax, fitminY, fitmaxY);
  //TF2 *parabola = new TF2("parabola", "x*x*[0] + x*[1] + y*y*[2] + y*[3] + [4] + x*y*[5] + x*x*y*[6] + x*y*y*[7] + x*x*x*[8] + y*y*y*[9]", fitmin, fitmax, fitminY, fitmaxY);
  //TF2 *parabola = new TF2("parabola", "[0]*pow(x-[1],2) + [2]*pow(y-[3], 2) + [4]", fitmin, fitmax, fitminY, fitmaxY);
  //parabola->SetParameters(100., centreVal, 200., massArray_[locMin.second], LLvalues[locMin.first][locMin.second], 0.02);
  //parabola->SetParameter(1, centreVal);
  //parabola->SetParameter(3, massArray_[locMin.second]);
  //parabola->SetParameter(4, LLvalues[locMin.first][locMin.second]);
//  parabola->SetParLimits(1, fitmin, fitmax);
  //parabola->SetParLimits(3, fitminY, fitmaxY);
  //parabola->SetParLimits(4, LLvalues[locMin.first][locMin.second] - 30., LLvalues[locMin.first][locMin.second] + 30.);
  parabola->SetParameters(0., 200., 50., 100., 400., -0.1, 1.);
  
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-3);
  TFitResultPtr p = g->Fit(parabola,"SWRME");
  if (verbose_) p->Print("V");
  
  g->Write();
  TF2 *fitresult = (TF2*) g->FindObject("parabola");
  fitresult->Write();
  //parabola->Write();
  
  double outputWidth, outputMass;
  double minimum = fitresult->GetMinimumXY(outputWidth, outputMass);
  
  TF12 *parbx = new TF12("parbx", fitresult, outputMass, "x");
  TF12 *parby = new TF12("parby", fitresult, outputWidth, "y");
  
  double lowerSigma = parbx->GetX(minimum + 0.5, fitmin-interval, outputWidth);
  double upperSigma = parbx->GetX(minimum + 0.5, outputWidth, fitmax+interval);
  double sigma = (upperSigma - lowerSigma)/2.;
  if ( lowerSigma <= fitmin && upperSigma >= fitmax )
    std::cerr << "Likelihood2D::CalculateOutputWidth: ERROR: Uncertainty calculation limited by fit boundaries. Do not trust..." << std::endl;
  else if ( lowerSigma <= fitmin ) sigma = upperSigma - outputWidth;
  else if ( upperSigma >= fitmax ) sigma = outputWidth - lowerSigma;
  
  //std::cout << "Minimum -log(like) value is " << parabola->Eval(outputWidth) << std::endl;
  std::cout << "Minimum -log(like) value is " << minimum << std::endl;
  if (verbose_) std::cout << "Position of minimum is (" << outputWidth << ", " << outputMass << ")" << std::endl;
  if (verbose_) std::cout << "lowerSigma = " << lowerSigma << " and upperSigma = " << upperSigma << std::endl;
  
  //double LLreduced[nn][nMasses_];
  //for (int i = 0; i < nn; i++)
    //LLreduced[i] = LLvalues[i] - minimum;
    //LLreduced[i] = LLvalues[i] - parabola->Eval(outputWidth);
//   TGraph2D *g2 = new TGraph2D;
//   for (int iWidth = 0; iWidth < nWidths_; iWidth++)
//   {
//     for (int iMass = 0; iMass < nMasses_; iMass++)
//     {
//       temp = LLvalues[iWidth][iMass] - minimum;
//       //if ( temp < 20. )
//         g2->SetPoint(g2->GetN(), evalWidths[iWidth], massArray_[iMass], temp);
//     }
//   }
//   g2->Fit(parabola,"R");
//   
  
  this->DrawOutputLogLikelihood(g, fitresult, fitmin, fitmax, fitminY, fitmaxY, plotName, writeToFile);
  this->DrawOutputLogLikelihood(g, fitresult, outputWidth, outputMass, minimum, plotName, writeToFile);
//   this->DrawOutputLogLikelihood(g2, parabola, 0, 15, 1200., plotName+"_full", writeToFile);
//   this->DrawOutputLogLikelihood(g2, parabola, outputWidth-2., outputWidth+2., 1.5*g2->Eval(outputWidth+2.), plotName, writeToFile);
//   this->DrawOutputLogLikelihood(g2, parabola, outputWidth-0.5, outputWidth+0.5, std::max(g2->Eval(outputWidth-3.*sigma),g2->Eval(outputWidth+0.5)), plotName+"_zoom", writeToFile);
  
  if (makeNewFile)
  {
    filePlots_->Close();
    delete filePlots_;
  }
  
//  delete g2;
  delete parbx;
  delete parby;
  delete fitresult;
  delete parabola;
//  delete g;
  
  return std::pair<double,double>(outputWidth,sigma);
}

std::pair<double,double> Likelihood2D::ApplyCalibrationCurve(double thisOutputWidth, double thisOutputWidthSigma)
{
  /// return 'thisInputWidth' and uncertainty
  //  thisOutputWidth = Par_[0] + Par_[1] * thisInputWidth

  double thisInputWidth = (thisOutputWidth - calCurvePar_[0])/calCurvePar_[1];
  double thisInputWidthSigma = TMath::Sqrt( thisOutputWidthSigma*thisOutputWidthSigma + calCurveParUnc_[0]*calCurveParUnc_[0] + calCurveParUnc_[1]*calCurveParUnc_[1]*thisInputWidth*thisInputWidth )/calCurvePar_[1];
  
  return std::pair<double,double>(thisInputWidth,thisInputWidthSigma);
  
}

int Likelihood2D::InitPull(int nPsExp)
{
  mkdir((dirNameLLTxt_+dirNamePull_).c_str(),0777);
  
  nPsExp_ = nPsExp;
  if ( nPsExp > 1000 )
  {
    std::cout << "Likelihood2D::Pull: Warning: Only 1000 pseudo experiments will be performed" << std::endl;
    return 1000;
  }
  else if (verbose_) std::cout << "Likelihood2D::Pull: Performing " << nPsExp << " pseudo experiments" << std::endl;
  
  return nPsExp;
}

void Likelihood2D::AddPsExp(int thisPsExp, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, double inputWidth, double inputMass, bool isTTbar, bool isData)
{
  if (isData) std::cerr << "Likelihood2D::Pull: Will not use data for pseudo experiments..." << std::endl;
  else if (calledLLCalculation_)
  {
    if (isTTbar)
    {
      if (rewHadOnly_)
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass);
      }
      else
      {
        thisWidthSF_ = rew_->EventWeightCalculatorNonRel(hadTopMassForWidthSF, inputWidth) * rew_->EventWeightCalculatorNonRel(lepTopMassForWidthSF, inputWidth);
        thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, inputMass) * rew_->MassEventWeightCalculatorNonRel(lepTopMassForWidthSF, inputMass);
      }
    }
    else
    {
      thisWidthSF_ = 1.;
      thisMassSF_ = 1.;
    }
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      for (int iMass = 0; iMass < nMasses_; iMass++)
      {
        loglike_pull_[iWidth][iMass][thisPsExp] += loglike_per_evt_[iWidth][iMass]*scaleFactor*thisWidthSF_*thisMassSF_;
      }
    }
  }
  else std::cerr << "Likelihood2D::Pull: Did not calculate likelihoods! Cannot get input for pseudo experiments..." << thisPsExp << std::endl;
}

void Likelihood2D::CalculatePull(double inputWidth)
{
  std::string fileName = dirNameLLTxt_+dirNamePull_+"PseudoExperiments_widthx"+tls_->DotReplace(inputWidth)+".root";
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
    thisOutputWidth[iPsExp] = this->GetOutputWidth(inputWidth, iPsExp);
    thisInputWidth[iPsExp] = this->ApplyCalibrationCurve(thisOutputWidth[iPsExp].first, thisOutputWidth[iPsExp].second);
    if ( thisInputWidth[iPsExp].first != -1. ) aveInputWidth += thisInputWidth[iPsExp].first;
    else std::cerr << "Likelihood2D::CalculatePull: Input width for pseudo experiment " << iPsExp << " is equal to -1! Ignoring this pseudo experiment... " << std::endl;
    hAveWidthOut->Fill(thisOutputWidth[iPsExp].first);
    hAveWidthIn->Fill(thisInputWidth[iPsExp].first);
    hUncAveWidthOut->Fill(thisOutputWidth[iPsExp].second);
    hUncAveWidthIn->Fill(thisInputWidth[iPsExp].second);
  }
  aveInputWidth = aveInputWidth/nPsExp_;
  //std::cout << "Output width : " << thisOutputWidth[0].first << " ;   sigma : " << thisOutputWidth[0].second << std::endl;
  //std::cout << "Input width  : " << thisInputWidth[0].first << " ;   sigma : " << thisInputWidth[0].second << std::endl;
  std::cout << "Average width for pseudo experiments is " << aveInputWidth << std::endl;
  
  WritePsExpOutput(thisOutputWidth, thisInputWidth, inputWidth);
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

void Likelihood2D::AddToFraction(int d, double scaleFactor, double hadTopMassForWidthSF, double lepTopMassForWidthSF, bool isTTbar, bool isCM, bool isWM, bool isUM)
{
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    /// Calculate massSF
    if (isTTbar)
    {
      if (rewHadOnly_) thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, massArray_[iMass]);
      else thisMassSF_ = rew_->MassEventWeightCalculatorNonRel(hadTopMassForWidthSF, massArray_[iMass]) * rew_->MassEventWeightCalculatorNonRel(lepTopMassForWidthSF, massArray_[iMass]);
    }
    else thisMassSF_ = 1.;
    
    if      (isCM) nEventsCMFractions_[iMass][d] += scaleFactor * thisMassSF_;
    else if (isWM) nEventsWMFractions_[iMass][d] += scaleFactor * thisMassSF_;
    else if (isUM) nEventsUMFractions_[iMass][d] += scaleFactor * thisMassSF_;
  }
}

void Likelihood2D::CalculateFractions(std::vector<std::string> datasetNames)
{
  std::string fileName = "";
  int nDatasets = datasetNames.size();
  for (int iMass = 0; iMass < nMasses_; iMass++)
  {
    for (int d = 0; d < nDatasets; d++)
    {
      fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_mass"+stringMassArray_[iMass]+".txt";
      txtOutputFractions_.open(fileName.c_str());
      txtOutputFractions_ << nEventsCMFractions_[iMass][d] << "   " << nEventsWMFractions_[iMass][d] << "   " << nEventsUMFractions_[iMass][d] << std::endl;
      txtOutputFractions_.close();
    }  // end datasets
  }  // end masses
}

void Likelihood2D::GetFractions(double *fractions, int nCats, std::vector<std::string> datasetNames, std::vector<int> includeDataset)
{
  Likelihood2D::GetFractions(fractions, nCats, "172p5", datasetNames, includeDataset);
}

void Likelihood2D::GetFractions(double *fractions, int nCats, std::string mass, std::vector<std::string> datasetNames, std::vector<int> includeDataset)
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
    fileName = dirNameNEvents_+"nEvents_"+datasetNames[d]+"_mass"+mass+".txt";
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
    if ( (nEvents[d]).size() != nCats_ ) std::cerr << "Something wrong with number of categories..." << std::endl;
    
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
      std::cout << "Mass: " << mass << "   # CM: " << 100.*fractions[0] << "%   # WM: " << 100.*fractions[1];
      if ( nCats > 2 ) std::cout << "%   # UM: " << 100.*fractions[2];
      std::cout << "%  " << std::endl;
    }
  }
}

void Likelihood2D::Make2DGraph(std::string name, bool makeNewFile)
{
  TGraph2D *g = new TGraph2D;
  TGraph2D *g2 = new TGraph2D;
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    for (int iMass = 0; iMass < nMasses_; iMass++)
    {
      g->SetPoint(g->GetN(), widthArray_[iWidth], massArray_[iMass], loglike_[iWidth][iMass]);
      if ( widthArray_[iWidth] > 0.2 && widthArray_[iWidth] < 1.9 )
        g2->SetPoint(g2->GetN(), widthArray_[iWidth], massArray_[iMass], loglike_[iWidth][iMass]);
    }
  }
  
  
  if (makeNewFile)
  {
    filePlots_ = new TFile((dirNameLLTxt_+"File_"+name+".root").c_str(), "RECREATE");
    filePlots_->cd();
  }
  
  TCanvas *c = new TCanvas(name.c_str(), name.c_str());
  c->cd();
  gStyle->SetPalette(1);
  g->SetMaxIter(500000);
  g->Draw("AP");
  g->SetName(name.c_str());
  g->SetTitle((name+";#Gamma_{t}/#Gamma_{t,gen}; M_{t}; L(m_{r}|#Gamma_{t},M_{t})").c_str());
  g->GetXaxis()->SetTitleOffset(1.7);
  g->GetYaxis()->SetTitleOffset(1.7);
  g->Draw("surf3");
  c->Update();
  c->Write();
  c->SaveAs((dirNameLLTxt_+name+"_surf3.png").c_str());
  c->SetName((name+"_zoom").c_str());
  g->GetXaxis()->SetRangeUser(0.1,2.);
  c->Update();
  c->Write();
  c->SaveAs((dirNameLLTxt_+name+"_surf3_zoom.png").c_str());
  c->SetName((name+"_colz_zoom").c_str());
  g->Draw("colz");
  g->GetXaxis()->SetLimits(0.1,2.);
  c->Update();
  c->Write();
  c->SaveAs((dirNameLLTxt_+name+"_colz_zoom.png").c_str());
  c->SetName((name+"_colz_hardzoom").c_str());
  g2->Draw("colz");
  c->Update();
  c->Write();
  c->SaveAs((dirNameLLTxt_+name+"_colz_hardzoom.png").c_str());
  c->SetName((name+"_cont1z_hardzoom").c_str());
  g2->Draw("cont1z");
  c->Update();
  c->Write();
  c->SaveAs((dirNameLLTxt_+name+"_cont1z_hardzoom.png").c_str());
  
  if (makeNewFile)
  {
    filePlots_->Close();
    delete filePlots_;
  }
  
  delete c;
  delete g;
}

void Likelihood2D::DrawGraph(TH1D* h, TGraph* g, std::string name)
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

void Likelihood2D::DrawGraph(TGraph2D* g, std::string name)
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

void Likelihood2D::DrawLikelihoods()
{
  Color_t colours[] = {kRed, kOrange-3, kYellow-7, kGreen-7, kGreen+1, kCyan+1, kBlue+2, kMagenta, kViolet-5, kPink+10};
  TCanvas* c2 = new TCanvas("-log(likelihood)", "-log(likelihood)");
  c2->cd();
  graph_["likelihood_widthx"+stringWidthArray_[0]+"_mass172p5"]->SetLineColor(colours[0]);
  graph_["likelihood_widthx"+stringWidthArray_[0]+"_mass172p5"]->Draw();
  std::string temp = graph_["likelihood_widthx"+stringWidthArray_[0]+"_mass172p5"]->GetTitle();
  graph_["likelihood_widthx"+stringWidthArray_[0]+"_mass172p5"]->SetTitle("-LogLikelihood");
  c2->Update();
  for (int i = 1; i < nWidths_; i++)
  {
    graph_["likelihood_widthx"+stringWidthArray_[i]+"_mass172p5"]->SetLineColor(colours[i%10]);
    graph_["likelihood_widthx"+stringWidthArray_[i]+"_mass172p5"]->Draw("same");
    c2->Update();
  }
  c2->Write();
  c2->SaveAs((outputDirName_+"LogLikelihood_multi.png").c_str());
  c2->Close();
  
  graph_["likelihood_widthx"+stringWidthArray_[0]+"_mass172p5"]->SetTitle(temp.c_str());
  
  delete c2;
}

void Likelihood2D::DrawOutputLogLikelihood(TGraph* g, TF1* f, double minX, double maxX, double maxY, std::string name, bool writeToFile)
{
  std::string outputFileName = dirNameLLTxt_+name+".png";
  TCanvas* c1 = new TCanvas(name.c_str(), "-LogLikelihood vs. width");
  c1->cd();
  g->Draw("AP");
  if ( minX < 0. ) minX = 0.;
  g->SetTitle("");
  g->GetXaxis()->SetLimits(minX,maxX);
  g->GetYaxis()->SetLimits(-0.05*maxY,maxY);
  g->GetXaxis()->SetTitle("#Gamma_{t} /#Gamma_{t,gen}");
  g->GetYaxis()->SetTitle("\\Delta\\mathscr{L}(\\mathrm{\\Gamma_{t}})");
  g->SetMarkerStyle(2);  //kPlus
  g->Draw("P");
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

void Likelihood2D::DrawOutputLogLikelihood(TGraph2D* g, TF2* f, double minX, double maxX, double minY, double maxY, std::string name, bool writeToFile)
{
  std::string outputFileName = dirNameLLTxt_+name+".png";
  TCanvas* c1 = new TCanvas(name.c_str(), "-LogLikelihood vs. width");
  c1->cd();
  gStyle->SetPalette(1);
  g->SetMaxIter(500000);
  g->Draw("AP");
  if ( minX < 0. ) minX = 0.;
  if ( minY < 169.5 ) minY = 169.5;
  if ( maxY > 175.5 ) maxY = 175.5;
  g->SetTitle("");
  g->GetXaxis()->SetLimits(minX,maxX);
  g->GetYaxis()->SetLimits(minY,maxY);
  g->GetXaxis()->SetTitle("#Gamma_{t} /#Gamma_{t,gen}");
  g->GetYaxis()->SetTitle("m_{t} (GeV)");
//  g->GetZaxis()->SetTitle("\\Delta\\mathscr{L}(\\mathrm{\\Gamma_{t}})");
  //g->SetMarkerStyle(2);  //kPlus
  g->SetMarkerStyle(20);
  g->Draw("PCOLZ");
  //if (writeToFile) g->Write();
  c1->Update();
  //f->Draw("same");
  //c1->Modified();
  //c1->Update();
//   TLine l;
//   l.SetLineColor(921);
//   l.DrawLine(minX, 0.5, maxX, 0.5);
//   c1->Update();
  
  if (writeToFile) c1->Write();
  c1->SaveAs(outputFileName.c_str());
  c1->Close();
  
  delete c1;
}

void Likelihood2D::DrawOutputLogLikelihood(TGraph2D* g, TF2* f, double minWidth, double minMass, double minimum, std::string name, bool writeToFile)
{
  std::string outputFileName = dirNameLLTxt_+name;
  gStyle->SetPalette(1);
  g->SetMaxIter(500000);
  TH1D *hx = (TH1D*) g->Project("x");
  TH1D *hy = (TH1D*) g->Project("y");
  TF12 *fx = new TF12("fx", f, minMass, "x");
  TF12 *fy = new TF12("fy", f, minWidth, "y");
  double minX, maxX;
  TCanvas* c1 = new TCanvas((name+"_X").c_str(), "-LogLikelihood vs. width (projection)");
  c1->cd();
  hx->SetTitle("");
  minX = minWidth-0.5;
  if ( minX < 0. ) minX = 0.;
  maxX = minWidth+1.;
  hx->GetXaxis()->SetRangeUser(minX,maxX);
  //hx->GetYaxis()->SetRangeUser(minY,maxY);
  hx->GetXaxis()->SetTitle("#Gamma_{t} /#Gamma_{t,gen}");
  //hy->GetXaxis()->SetTitle("m_{t} (GeV)");
  hx->GetYaxis()->SetTitle("\\Delta\\mathscr{L}(\\mathrm{\\Gamma_{t}})");
  hx->SetMarkerStyle(2);  //kPlus
  hx->Draw("P");
  //if (writeToFile) g->Write();
  c1->Update();
  fx->SetTitle("");
  fx->Draw("same");
  c1->Modified();
  c1->Update();
  TLine l;
  l.SetLineColor(921);
  l.DrawLine(minX, minimum+0.5, maxX, minimum+0.5);
  c1->Update();
  if (writeToFile) c1->Write();
  c1->SaveAs((outputFileName+"_projX.png").c_str());
  
  c1->SetName((name+"_Y").c_str());
  c1->SetTitle("-LogLikelihood vs. mass (projection)");
  hy->SetTitle("");
  minX = minMass-1.5;
  if ( minX < 0. ) minX = 0.;
  maxX = minMass+1.5;
  hy->GetXaxis()->SetRangeUser(minX,maxX);
  //hy->GetYaxis()->SetRangeUser(minY,maxY);
  hy->GetXaxis()->SetTitle("m_{t} (GeV)");
  hy->GetYaxis()->SetTitle("\\Delta\\mathscr{L}(\\mathrm{\\Gamma_{t}})");
  hy->SetMarkerStyle(2);  //kPlus
  hy->Draw("P");
  c1->Update();
  fy->SetTitle("");
  fy->Draw("same");
  c1->Modified();
  c1->Update();
  l.SetLineColor(921);
  l.DrawLine(minX, minimum+0.5, maxX, minimum+0.5);
  c1->Update();
  if (writeToFile) c1->Write();
  c1->SaveAs((outputFileName+"_projY.png").c_str());
  c1->Close();
  
  delete c1;
}

void Likelihood2D::ReadLLValuesFromFile(std::string inputFileName, std::string inputDir)
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
    std::cerr << "Likelihood2D::ReadLLValuesFromFile : Input file " << inputFileName << " cannot be read:" << std::endl;
    std::cerr << "                                   # widths != # likelihood values !!" << std::endl;
    vecWidthFromFile_.clear();
    vecLLValsFromFile_.clear();
  }
}

void Likelihood2D::WritePsExpOutput(std::pair<double,double> *outputWidth, std::pair<double,double> *inputWidth, double genWidth)
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

void Likelihood2D::WriteFuncOutput(int nPoints, double *arrayCentre, double *arrayContent, std::string name)
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

void Likelihood2D::WriteOutput(int nPoints, double width, double *arrayCentre, double *arrayContent, std::string name, int dim)
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

void Likelihood2D::CombineOutput(std::string mass)
{
  std::ofstream combFile((outputDirName_+dirNameTGraphTxt_+"output_tgraph2d_total"+"_mass"+mass+".txt").c_str());
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    std::string inFileName = outputDirName_+dirNameTGraphTxt_+"output_tgraph2d_widthx"+stringWidthArray_[iWidth]+"_mass"+mass+".txt";
    std::ifstream infile(inFileName.c_str());
    combFile << infile.rdbuf();
    infile.close();
  }
  combFile.close();
}

void Likelihood2D::PrintLikelihoodOutput(std::string llFileName)
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
    txtOutputLogLike << std::setw(12) << std::right << loglike_[iWidth][1] << "  ";
  }
  txtOutputLogLike << std::endl;
  txtOutputLogLike.close();
}

void Likelihood2D::PrintLikelihoodOutputData(std::string llFileName)
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

void Likelihood2D::PrintMtmLikelihoodOutput(std::string llFileName)
{
  std::ofstream txtOutputLogLike((dirNameLLTxt_+llFileName).c_str());
  txtOutputLogLike << "likelihood values (all MC) : {";
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    txtOutputLogLike << std::fixed << std::setprecision(15) << loglike_[iWidth][1];
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
