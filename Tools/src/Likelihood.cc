#include "../interface/Likelihood.h"

const double Likelihood::widthArray_[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10.};

const std::string Likelihood::listCats_[] = {"CM", "WM", "NM"};

const int Likelihood::nWidths_ = sizeof(widthArray_)/sizeof(widthArray_[0]);
const int Likelihood::nCats_ = sizeof(listCats_)/sizeof(listCats_[0]);

std::string Likelihood::stringWidthArray_[nWidths_] = {""};

double Likelihood::loglike_[nWidths_] = {0};
double Likelihood::loglike_data_[nWidths_] = {0};
double Likelihood::loglike_per_evt_[nWidths_] = {0};
double Likelihood::loglike_good_evts_[nWidths_] = {0};
double Likelihood::loglike_good_evts_data_[nWidths_] = {0};
double Likelihood::loglike_CM_[nWidths_] = {0};
double Likelihood::loglike_CM_per_evt_[nWidths_] = {0};
double Likelihood::loglike_CM_good_evts_[nWidths_] = {0};
double Likelihood::loglike_gen_[nWidths_] = {0};
double Likelihood::loglike_gen_per_evt_[nWidths_] = {0};
double Likelihood::loglike_gen_good_evts_[nWidths_] = {0};

double Likelihood::loglike_pull_[nWidths_][1000] = {0};
double Likelihood::loglike_pull_single_[nWidths_] = {0};

//const double Likelihood::calCurvePar_[2] = {0., 1.};  // at the moment no output calibration
//const double Likelihood::calCurveParUnc_[2] = {0., 0.};  // at the moment no output calibration
const double Likelihood::calCurvePar_[2] = {0.028244, 0.974992};
const double Likelihood::calCurveParUnc_[2] = {0.00620311, 0.00200199};


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

Likelihood::Likelihood(double min, double max, std::string outputDirName, std::string date, bool makeHistograms, bool calculateGoodEvtLL, bool verbose):
verbose_(verbose), outputDirName_(outputDirName), dirNameTGraphTxt_("OutputTxt/"), dirNameLLTxt_("OutputLikelihood/"+date+"/"), dirNamePull_("PseudoExp/"), inputFileName_(""), suffix_(""), histoName_(""), minRedMass_(min), maxRedMass_(max), histo_(), histoSm_(), histoTotal_(), graph_(), vecBinCentres_(), vecBinContents_(), calculateGoodEvtLL_(calculateGoodEvtLL), calledLLCalculation_(false), calledCMLLCalculation_(false), calledGenLLCalculation_(false), vecWidthFromFile_(), vecLLValsFromFile_(), vecGoodLLValsFromFile_()
{
  tls_ = new HelperTools();
  rew_ = new EventReweighting(false);  // no correction for number of events
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    stringWidthArray_[iWidth] = tls_->DotReplace(widthArray_[iWidth]);
  }
  
  for (int i = 0; i < nWidths_; i++)
  {
    loglike_[i] = 0.;
    loglike_data_[i] = 0.;
    loglike_good_evts_[i] = 0.;
    loglike_good_evts_data_[i] = 0.;
    loglike_CM_[i] = 0.;
    loglike_CM_good_evts_[i] = 0.;
    loglike_gen_[i] = 0.;
    loglike_gen_good_evts_[i] = 0.;
  }
  
  if (makeHistograms)
  {
    mkdir(outputDirName_.c_str(),0777);
    mkdir((outputDirName_+dirNameTGraphTxt_).c_str(),0777);
    this->BookHistograms();
  }
  else mkdir(dirNameLLTxt_.c_str(),0777);
}

Likelihood::~Likelihood()
{
  delete file_;
  delete rew_;
  delete tls_;
}

void Likelihood::BookHistograms()
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    thisWidth_ = stringWidthArray_[iWidth];
    
    histo_[("Red_top_mass_CM_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_CM_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", CM; M_{t}/<M_{t}>").c_str(), 90, 0.5, 2.0);
    histo_[("Red_top_mass_CM_widthx"+thisWidth_+"_100b").c_str()] = new TH1D(("Red_top_mass_CM_widthx"+thisWidth_+"_100b").c_str(),("Reduced top mass for width "+thisWidth_+", CM; M_{t}/<M_{t}>").c_str(), 100, 0.5, 2.0);
    histo_[("Red_top_mass_CM_widthx"+thisWidth_+"_120b").c_str()] = new TH1D(("Red_top_mass_CM_widthx"+thisWidth_+"_120b").c_str(),("Reduced top mass for width "+thisWidth_+", CM; M_{t}/<M_{t}>").c_str(), 120, 0.5, 2.0);
    
    histo_[("Red_top_mass_WM_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_WM_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", WM; M_{t}/<M_{t}>").c_str(), 90, 0.5, 2.0);
    histo_[("Red_top_mass_WM_widthx"+thisWidth_+"_100b").c_str()] = new TH1D(("Red_top_mass_WM_widthx"+thisWidth_+"_100b").c_str(),("Reduced top mass for width "+thisWidth_+", WM; M_{t}/<M_{t}>").c_str(), 100, 0.5, 2.0);
    histo_[("Red_top_mass_WM_widthx"+thisWidth_+"_120b").c_str()] = new TH1D(("Red_top_mass_WM_widthx"+thisWidth_+"_120b").c_str(),("Reduced top mass for width "+thisWidth_+", WM; M_{t}/<M_{t}>").c_str(), 120, 0.5, 2.0);
    
    histo_[("Red_top_mass_NM_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_NM_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", NM; M_{t}/<M_{t}>").c_str(), 90, 0.5, 2.0);
    histo_[("Red_top_mass_NM_widthx"+thisWidth_+"_100b").c_str()] = new TH1D(("Red_top_mass_NM_widthx"+thisWidth_+"_100b").c_str(),("Reduced top mass for width "+thisWidth_+", NM; M_{t}/<M_{t}>").c_str(), 100, 0.5, 2.0);
    histo_[("Red_top_mass_NM_widthx"+thisWidth_+"_120b").c_str()] = new TH1D(("Red_top_mass_NM_widthx"+thisWidth_+"_120b").c_str(),("Reduced top mass for width "+thisWidth_+", NM; M_{t}/<M_{t}>").c_str(), 120, 0.5, 2.0);
  }
}

void Likelihood::FillHistograms(double redMass, double lumiWeight, double massForWidthSF, bool isTTbar, bool isData, std::string catSuffix)
{
  if ( ! isData && redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      thisWidth_ = stringWidthArray_[iWidth];
      if (isTTbar) thisWidthSF_ = rew_->EventWeightCalculator(massForWidthSF, widthArray_[iWidth]);
      else thisWidthSF_ = 1.;
      
      histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_90b").c_str()]->Fill(redMass, lumiWeight*thisWidthSF_);
      histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_100b").c_str()]->Fill(redMass, lumiWeight*thisWidthSF_);
      histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_120b").c_str()]->Fill(redMass, lumiWeight*thisWidthSF_);
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

void Likelihood::ConstructTGraphsFromHisto(std::string tGraphFileName)
{
  int binMin, binMax, totEvents;
  int nBins[nCats_] = {0};
  int nEvents[nCats_] = {0};
  double fracCats[nCats_] = {0};
  
  if ( gLL2D_ == NULL ) gLL2D_ = new TGraph2D();
  
  fileTGraphs_ = new TFile((outputDirName_+tGraphFileName).c_str(), "RECREATE");
  fileTGraphs_->cd();
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    /// Clear vars
    totEvents = 0;
    for (int i = 0; i < nCats_; i++)
    {
      nBins[i] = 0;
      nEvents[i] = 0;
      fracCats[i] = 0.;
      //vecBinCentres_[i].clear();
      //vecBinContents_[i].clear();
    }
    
    suffix_ = "widthx"+stringWidthArray_[iWidth];
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      /// Clear vars
      binMin = 0;
      binMax = 100;
      
      /// Get histo to smooth
      histoName_ = listCats_[iCat]+"_"+suffix_;
      histoSm_[histoName_] = (TH1D*) histo_["Red_top_mass_"+histoName_+"_90b"]->Clone(histoName_.c_str());
      histoSm_[histoName_]->Smooth();
      
      nBins[iCat] = histoSm_[histoName_]->GetNbinsX();
      if ( iCat > 0 && nBins[iCat] != nBins[iCat-1] )
      {
        std::cerr << "Likelihood::ConstructTGraphs: Histograms should have the same number of bins! Exiting..." << std::endl;
        exit(1);
      }
      
      binMin = histoSm_[histoName_]->FindBin(0.6);
      binMax = histoSm_[histoName_]->FindBin(1.4)+1;
      for (int iBin = binMin; iBin < binMax+1; iBin++) { nEvents[iCat] += histoSm_[histoName_]->GetBinContent(iBin);}
      totEvents += nEvents[iCat];
      
      if (verbose_) std::cout << "Likelihood::ConstructTGraphs: " << histoName_ << ": " << nEvents[iCat] << " events" << std::endl;
      
      /// Normalise histo on relevant subdomain
      Double_t integral = histoSm_[histoName_]->Integral(binMin, binMax);
      histoSm_[histoName_]->Scale(1./integral);
      
      /// Get bin centres & contents
      for (int iBin = binMin; iBin < binMax+1; iBin++)
      {
        (vecBinCentres_[histoName_]).push_back(histoSm_[histoName_]->GetBinCenter(iBin));
        (vecBinContents_[histoName_]).push_back(histoSm_[histoName_]->GetBinContent(iBin));
      }
      
    }  // end cats
    
    if (verbose_) std::cout << "Likelihood::ConstructTGraphs: Total: " << totEvents << " events" << std::endl;
    
    /// Make arrays as input for TGraph
    const int nPoints = (vecBinCentres_[listCats_[0]+"_"+suffix_]).size();
    double binCentreArray[nPoints] = {0.}, binContentArray[nCats_][nPoints] = {0.}, likelihoodCMArray[nPoints] = {0.}, totalBinContentArray[nPoints] = {0.}, likelihoodArray[nPoints] = {0.};
    for (int i = 0; i < nPoints; i++)
    {
      binCentreArray[i] = (vecBinCentres_[listCats_[0]+"_"+suffix_]).at(i);
      for (int iCat = 0; iCat < nCats_; iCat++)
      {
        binContentArray[iCat][i] = (vecBinContents_[listCats_[iCat]+"_"+suffix_]).at(i);
        // For CM events
        if ( iCat == 0 ) { likelihoodCMArray[i] = -TMath::Log(binContentArray[iCat][i]);}
        // Calculate event fractions (one time)
        if ( i == 0 ) { fracCats[iCat] = (double)nEvents[iCat]/(double)totEvents;}
      }
    }
    
    /// Make TGraphs
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      histoName_ = listCats_[iCat]+"_"+suffix_;
      
      graph_[histoName_] = new TGraph(nPoints, binCentreArray, binContentArray[iCat]);
      graph_[histoName_]->SetName(("g"+histoName_).c_str());
      graph_[histoName_]->SetTitle(histoName_.c_str());
      graph_[histoName_]->Write();
      DrawGraph(histoSm_[histoName_], graph_[histoName_], "Graph_Red_top_mass_"+histoName_);
      if ( iCat == 0 )  //CM
      {
        graph_["likelihood_CM_"+suffix_] = new TGraph(nPoints, binCentreArray, likelihoodCMArray);
        graph_["likelihood_CM_"+suffix_]->SetName(("likelihood_CM_"+suffix_).c_str());
        graph_["likelihood_CM_"+suffix_]->SetTitle(("likelihood_CM_"+suffix_).c_str());
        graph_["likelihood_CM_"+suffix_]->Write();
        
        WriteOutput(nPoints, iWidth, binCentreArray, likelihoodCMArray, "CorrectMatchLikelihood_"+suffix_, 1);
      }
      
      
      /// Calculate total probability
      fracCats[iCat] = (double)nEvents[iCat]/(double)totEvents;
      for (int i = 0; i < nPoints; i++) { totalBinContentArray[i] += fracCats[iCat]*binContentArray[iCat][i];}
      
      histoSm_[histoName_]->Scale(fracCats[iCat]);
      if ( iCat == 0 )
      {
        histoTotal_[suffix_] = (TH1D*) histoSm_[histoName_]->Clone("TotalProbability");
        histoTotal_[suffix_]->SetTitle("#frac{n_{CM}}{n_{tot}} * f_{CM}(x|#Gamma) + #frac{n_{WM}}{n_{tot}} * f_{WM}(x|#Gamma) + #frac{n_{NM}}{n_{tot}} * f_{NM}(x|#Gamma)");
      }
      else histoTotal_[suffix_]->Add(histoSm_[histoName_]);
    }
    
    if (verbose_) std::cout << "Likelihood::ConstructTGraphs: The integral of the weighted probability histogram is " << histoTotal_[suffix_]->Integral(binMin, binMax) << std::endl;
    
    graph_["TotalProbability_"+suffix_] = new TGraph(nPoints, binCentreArray, totalBinContentArray);
    graph_["TotalProbability_"+suffix_]->SetName(("TotalProbability_"+suffix_).c_str());
    graph_["TotalProbability_"+suffix_]->SetTitle(("TotalProbability_"+suffix_).c_str());
    graph_["TotalProbability_"+suffix_]->Write();
    DrawGraph(histoTotal_[suffix_], graph_["TotalProbability_"+suffix_], "Graph_totalProbability_"+suffix_);
    
    
    /// Negative log(likelihood)
    for (int i = 0; i < nPoints; i++)
    {
      likelihoodArray[i] = -TMath::Log(totalBinContentArray[i]); 
    }
    
    /// Make likelihood function
    graph_["likelihood_"+suffix_] = new TGraph(nPoints, binCentreArray, likelihoodArray);
    graph_["likelihood_"+suffix_]->SetName(("likelihood_"+suffix_).c_str());
    graph_["likelihood_"+suffix_]->SetTitle(("likelihood_"+suffix_).c_str());
    graph_["likelihood_"+suffix_]->Write();
    
    WriteOutput(nPoints, widthArray_[iWidth], binCentreArray, likelihoodArray, suffix_, 1);  // for TGraph
    WriteOutput(nPoints, widthArray_[iWidth], binCentreArray, likelihoodArray, suffix_, 2);  // for TGraph2D
    
    /// Fill TGraph2D
    for (int iCentre = 0; iCentre < nPoints; iCentre++)
    {
      gLL2D_->SetPoint(gLL2D_->GetN(), binCentreArray[iCentre], widthArray_[iWidth], likelihoodArray[iCentre]);
    }
  }  // end widths
  
  CombineOutput();
  
  DrawLikelihoods();
  gLL2D_->SetName("2D_likelihood");
  gLL2D_->SetTitle("2D_likelihood");
  gLL2D_->Write();
  DrawGraph(gLL2D_, "2D_likelihood");
  
  
  fileTGraphs_->Close();
  delete fileTGraphs_;
}

bool Likelihood::ConstructTGraphsFromFile()
{
  return this->ConstructTGraphsFromFile("");
}

bool Likelihood::ConstructTGraphsFromFile(std::string name)
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    suffix_ = name+"widthx"+stringWidthArray_[iWidth];
    inputFileName_ = outputDirName_+dirNameTGraphTxt_+"output_tgraph1d_"+suffix_+".txt";
    if ( ! tls_->fexists(inputFileName_.c_str()) )
    {
      std::cerr << "Likelihood::ConstructTGraphs: File " << inputFileName_ << " not found!!" << std::endl;
      //std::cerr << "                          Aborting the likelihood calculation..." << std::endl;
      return false;
    }
    graph_[suffix_] = new TGraph(inputFileName_.c_str());
  }
  std::cout << "Likelihood::ConstructTGraphs: Constructed TGraphs for likelihood measurements (using " << nWidths_ << " widths)" << std::endl;
  return true;
}

void Likelihood::CalculateLikelihood(double redMass, double lumiWeight, bool isData)
{
  this->CalculateLikelihood(redMass, lumiWeight, 1., 1., false, isData);  // isTTbar = false ==> thisWidthSF_ = 1.;
}

void Likelihood::CalculateLikelihood(double redMass, double lumiWeight, double massForWidthSF, double inputWidth, bool isTTbar, bool isData)
{
  if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! isData && ! calledLLCalculation_) calledLLCalculation_ = true;
    
    if (isTTbar) thisWidthSF_ = rew_->EventWeightCalculator(massForWidthSF, inputWidth);
    else thisWidthSF_ = 1.;
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      thisWidth_ = stringWidthArray_[iWidth];
      
      loglike_per_evt_[iWidth] = graph_["widthx"+thisWidth_]->Eval(redMass);
      if (! isData) loglike_[iWidth] += loglike_per_evt_[iWidth]*lumiWeight*thisWidthSF_;
      else loglike_data_[iWidth] += loglike_per_evt_[iWidth];
    }
    
    if (calculateGoodEvtLL_)
    {
      bool isGoodLL = false;
      for (int iWidth = 1; iWidth < nWidths_-1; iWidth++)
      {
        if ( loglike_per_evt_[0] > loglike_per_evt_[iWidth] && loglike_per_evt_[iWidth] < loglike_per_evt_[nWidths_-1] )
        {
          isGoodLL = true;
          break;
        }
      }
      if (isGoodLL)
      {
//        nofGoodEvtsLL[d]++;
        for (int iWidth = 0; iWidth < nWidths_; iWidth++)
        {
          if (! isData) loglike_good_evts_[iWidth] += loglike_per_evt_[iWidth]*lumiWeight*thisWidthSF_;
          else loglike_good_evts_data_[iWidth] += loglike_per_evt_[iWidth];
        }
      }
//       else
//       {
//         nofBadEvtsLL[d]++;
//       }
      
    }
  }
}

void Likelihood::CalculateCMLikelihood(double redMass, double massForWidthSF, double inputWidth, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "Likelihood::Cannot calculate loglikelihood for matched events when running over data" << std::endl;
    std::cerr << "            Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledCMLLCalculation_) calledCMLLCalculation_ = true;
    
    if (isTTbar) thisWidthSF_ = rew_->EventWeightCalculator(massForWidthSF, inputWidth);
    else thisWidthSF_ = 1.;
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      thisWidth_ = stringWidthArray_[iWidth];
      
      loglike_CM_per_evt_[iWidth] = graph_["CorrectMatchLikelihood_widthx"+thisWidth_]->Eval(redMass);
      loglike_CM_[iWidth] += loglike_CM_per_evt_[iWidth]*thisWidthSF_;
    }
    
    if (calculateGoodEvtLL_)
    {
      bool isGoodLL = false;
      for (int iWidth = 1; iWidth < nWidths_-1; iWidth++)
      {
        if ( loglike_CM_per_evt_[0] > loglike_CM_per_evt_[iWidth] && loglike_CM_per_evt_[iWidth] < loglike_CM_per_evt_[nWidths_-1] )
        {
          isGoodLL = true;
          break;
        }
      }
      if (isGoodLL)
      {
//        nofGoodEvtsLL[d]++;
        for (int iWidth = 0; iWidth < nWidths_; iWidth++)
        {
          loglike_CM_good_evts_[iWidth] += loglike_CM_per_evt_[iWidth]*thisWidthSF_;
        }
      }
//       else
//       {
//         nofBadEvtsLL[d]++;
//       }
      
    }
  }
}

void Likelihood::CalculateGenLikelihood(double redMass, double massForWidthSF, double inputWidth, bool isTTbar, bool isData)
{
  if (isData)
  {
    std::cerr << "Likelihood::Cannot calculate loglikelihood for generated events when running over data" << std::endl;
    std::cerr << "            Something went wrong here... Check..." << std::endl;
  }
  else if ( redMass > minRedMass_ && redMass < maxRedMass_ )
  {
    if (! calledGenLLCalculation_) calledGenLLCalculation_ = true;
    
    if (isTTbar) thisWidthSF_ = rew_->EventWeightCalculator(massForWidthSF, inputWidth);
    else thisWidthSF_ = 1.;
    
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      thisWidth_ = stringWidthArray_[iWidth];
      
      loglike_gen_per_evt_[iWidth] = graph_["CorrectMatchLikelihood_widthx"+thisWidth_]->Eval(redMass);
      loglike_gen_[iWidth] += loglike_gen_per_evt_[iWidth]*thisWidthSF_;
    }
    
    if (calculateGoodEvtLL_)
    {
      bool isGoodLL = false;
      for (int iWidth = 1; iWidth < nWidths_-1; iWidth++)
      {
        if ( loglike_gen_per_evt_[0] > loglike_gen_per_evt_[iWidth] && loglike_gen_per_evt_[iWidth] < loglike_gen_per_evt_[nWidths_-1] )
        {
          isGoodLL = true;
          break;
        }
      }
      if (isGoodLL)
      {
//        nofGoodEvtsLL[d]++;
        for (int iWidth = 0; iWidth < nWidths_; iWidth++)
        {
          loglike_gen_good_evts_[iWidth] += loglike_gen_per_evt_[iWidth]*thisWidthSF_;
        }
      }
//       else
//       {
//         nofBadEvtsLL[d]++;
//       }
      
    }
  }
}

void Likelihood::GetOutputWidth(double inputWidth)
{
  this->GetOutputWidth(inputWidth, "");
}

void Likelihood::GetOutputWidth(double inputWidth, std::string type)
{
  std::string loglikePlotName = "loglikelihood_vs_width_";
  if (! type.empty() ) loglikePlotName += type+"_";
  loglikePlotName += "widthx"+tls_->DotReplace(inputWidth);
  
  if ( type.find("CM") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_CM_, loglikePlotName);
  else if ( type.find("gen") != std::string::npos || type.find("Gen") != std::string::npos )
    output_ = this->CalculateOutputWidth(nWidths_, loglike_gen_, loglikePlotName);
  else
    output_ = this->CalculateOutputWidth(nWidths_, loglike_, loglikePlotName);
  
  std::cout << "For an input width of " << inputWidth << " the minimum can be found at " << output_.first << " and the uncertainty is " << output_.second << std::endl;
}

void Likelihood::GetOutputWidth(std::string inputFileName, double inputWidth)
{
  if (verbose_) std::cout << "Using LogLikelihood values from file" << std::endl;
  std::string loglikePlotName = "loglikelihood_vs_width_";
  //if (! type.empty() ) loglikePlotName += type+"_";
  loglikePlotName += "ff_widthx"+tls_->DotReplace(inputWidth);
  
  output_ = this->CalculateOutputWidth(inputFileName, loglikePlotName);
  
  std::cout << "For an input width of " << inputWidth << " the minimum can be found at " << output_.first << " and the uncertainty is " << output_.second << std::endl;
}

std::pair<double,double> Likelihood::GetOutputWidth(double inputWidth, int thisPsExp)
{
  std::string loglikePlotName = dirNamePull_+"loglikelihood_vs_width_psExp_"+tls_->ConvertIntToString(thisPsExp)+"_widthx"+tls_->DotReplace(inputWidth);
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    loglike_pull_single_[iWidth] = loglike_pull_[iWidth][thisPsExp];
  }
  
  return this->CalculateOutputWidth(nWidths_, loglike_pull_single_, loglikePlotName, true);
}

std::pair<double,double> Likelihood::CalculateOutputWidth(std::string inputFileName, std::string plotName)
{
  this->ReadLLValuesFromFile(inputFileName);
  const int nn = vecWidthFromFile_.size();
  if ( nn == 0 )
    return std::pair<double,double>(-1.,-1.);
  double arrWidth[nn], arrLLVals[nn];
  for (int i = 0; i < nn; i++)
  {
    arrWidth[i] = vecWidthFromFile_.at(i);
    arrLLVals[i] = vecLLValsFromFile_.at(i);
  }
  return this->CalculateOutputWidth(nn, arrWidth, arrLLVals, plotName, false);
}

std::pair<double,double> Likelihood::CalculateOutputWidth(int nn, double* LLvalues, std::string plotName)
{
  return this->CalculateOutputWidth(nn, (double*)widthArray_, LLvalues, plotName, false);
}

std::pair<double,double> Likelihood::CalculateOutputWidth(int nn, double* LLvalues, std::string plotName, bool writeToFile)
{
  return this->CalculateOutputWidth(nn, (double*)widthArray_, LLvalues, plotName, writeToFile);
}

std::pair<double,double> Likelihood::CalculateOutputWidth(int nn, double* evalWidths, double* LLvalues, std::string plotName)
{
  return this->CalculateOutputWidth(nn, evalWidths, LLvalues, plotName, false);
}

std::pair<double,double> Likelihood::CalculateOutputWidth(int nn, double* evalWidths, double* LLvalues, std::string plotName, bool writeToFile)
{
  TGraph *g = new TGraph(nn, evalWidths, LLvalues);
  
  int locMin = LocMinArray(nn, LLvalues);
  //std::cout << "Index of minimum LL value is " << locMin << std::endl;
  double centreVal = evalWidths[locMin];
  double fitmax = centreVal + 0.8;
  double fitmin = centreVal - 0.8;
  if ( fitmin < 0. ) fitmin = 0.;
  
  if (verbose_) std::cout << "Likelihood::CalculateOutputWidth: Look for minimum around " << centreVal << std::endl;
  
  TF1 *parabola = new TF1("parabola", "pol2", fitmin, fitmax);
  
  g->Fit(parabola,"R");
  
  double outputWidth = parabola->GetMinimumX(fitmin, fitmax);
  double lowerSigma = parabola->GetX(parabola->Eval(outputWidth) + 0.5, fitmin, outputWidth);
  double upperSigma = parabola->GetX(parabola->Eval(outputWidth) + 0.5, outputWidth, fitmax);
  double sigma = (upperSigma - lowerSigma)/2.;
  
  double LLreduced[nn];
  for (int i = 0; i < nn; i++)
    LLreduced[i] = LLvalues[i] - parabola->Eval(outputWidth);
  TGraph *g2 = new TGraph(nn, evalWidths, LLreduced);
  g2->Fit(parabola,"R");
  
  
  this->DrawOutputLogLikelihood(g2, parabola, outputWidth-3., outputWidth+3., 1.5*g2->Eval(outputWidth+3.), plotName, writeToFile);
  this->DrawOutputLogLikelihood(g2, parabola, outputWidth-3.*sigma, outputWidth+3.*sigma, std::max(g2->Eval(outputWidth-3.*sigma),g2->Eval(outputWidth+3.*sigma)), plotName+"_zoom", writeToFile);
  
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

void Likelihood::AddPsExp(int thisPsExp, bool isData)
{
  if (isData) std::cerr << "Likelihood::Pull: Will not use data for pseudo experiments..." << std::endl;
  else if (calledLLCalculation_)
  {
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      loglike_pull_[iWidth][thisPsExp] += loglike_per_evt_[iWidth]*thisWidthSF_;
    }
  }
  else std::cerr << "Likelihood::Pull: Did not calculate likelihoods! Cannot get input for pseudo experiments..." << std::endl;
}

void Likelihood::CalculatePull(double inputWidth)
{
  std::string fileName = dirNameLLTxt_+dirNamePull_+"PseudoExperiments_widthx"+tls_->DotReplace(inputWidth)+".root";
  TFile *filePull = new TFile(fileName.c_str(),"RECREATE");
  
  /// Calculate output width for all pseudo experiments
  //  Transform to input width via calibration curve & calculate average input width
  std::pair<double,double> thisOutputWidth[nPsExp_] = {std::pair<double,double>(-1.,-1.)};
  std::pair<double,double> thisInputWidth[nPsExp_] = {std::pair<double,double>(-1.,-1.)};
  double aveInputWidth = 0.;
  
  for (int iPsExp = 0; iPsExp < nPsExp_; iPsExp++)
  {
    thisOutputWidth[iPsExp] = this->GetOutputWidth(inputWidth, iPsExp);
    thisInputWidth[iPsExp] = this->ApplyCalibrationCurve(thisOutputWidth[iPsExp].first, thisOutputWidth[iPsExp].second);
    aveInputWidth += thisInputWidth[iPsExp].first;
  }
  aveInputWidth = aveInputWidth/nPsExp_;
  std::cout << "Output width : " << thisOutputWidth[0].first << " ;   sigma : " << thisOutputWidth[0].second << std::endl;
  std::cout << "Input width  : " << thisInputWidth[0].first << " ;   sigma : " << thisInputWidth[0].second << std::endl;
  std::cout << "Average width for pseudo experiments is " << aveInputWidth << std::endl;
  
  /// Fill histogram with (Gamma_j - <Gamma>)/sigma_j
  double fillValue;
  filePull->cd();
  TH1D *hPull = new TH1D("hPull", "; (#Gamma_{j} - <#Gamma>)/#sigma_{#Gamma_{j}}", 50, -5., 5.);
  
  for (int iPsExp = 0; iPsExp < nPsExp_; iPsExp++)
  {
    fillValue = (thisInputWidth[iPsExp].first - aveInputWidth)/thisInputWidth[iPsExp].second;
    hPull->Fill(fillValue);
  }
  hPull->Write();
  
  /// Fit distribution with Gaussian function
  //  Should have mean = 0 and sigma = 1
  double fitMin = hPull->GetXaxis()->GetXmin();
  double fitMax = hPull->GetXaxis()->GetXmax();
  
  TF1 *gaus = new TF1("gaus", "gaus", fitMin, fitMax);
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

void Likelihood::DrawGraph(TH1D* h, TGraph* g, std::string name)
{
  TCanvas *c = new TCanvas(name.c_str(), name.c_str());
  c->cd();
  h->GetXaxis()->SetTitle("m_{t}/<m_{t}> [GeV]");
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
  g->Draw();
  g->SetName(name.c_str());
  g->SetTitle((name+"; m_{t}/<m_{t}> [GeV]; #Gamma/#Gamma_{SM}").c_str());
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
  g->GetXaxis()->SetRangeUser(minX,maxX);
  g->GetYaxis()->SetRangeUser(-0.05*maxY,maxY);
  g->GetXaxis()->SetTitle("#Gamma/#Gamma_{SM}");
  g->GetYaxis()->SetTitle("-Log(likelihood)");
  g->SetMarkerStyle(2);  //kPlus
  g->Draw("AP");
  //if (writeToFile) g->Write();
  c1->Update();
  f->Draw("same");
  c1->Modified();
  c1->Update();
  if (writeToFile) c1->Write();
  c1->SaveAs(outputFileName.c_str());
  c1->Close();
  
  delete c1;
}

void Likelihood::ReadLLValuesFromFile(std::string inputFileName)
{
  if (! tls_->fexists((dirNameLLTxt_+inputFileName).c_str()) )
  {
    std::cerr << "WARNING: File " << dirNameLLTxt_+inputFileName << " does not exist." << std::endl;
    exit(1);
  }
  fileIn_.open((dirNameLLTxt_+inputFileName).c_str());
  if (verbose_) std::cout << "Opening " << dirNameLLTxt_+inputFileName << "..." << std::endl;
  
  vecWidthFromFile_.clear();
  vecLLValsFromFile_.clear();
  vecGoodLLValsFromFile_.clear();
  
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
    if ( line.find("GoodLL") != std::string::npos )
    {
      iss >> var;
      while ( iss >> val )
      {
        vecGoodLLValsFromFile_.push_back(val);
      }
    }
    else if ( line.find("LL") != std::string::npos )
    {
      iss >> var;
      while ( iss >> val )
      {
        vecLLValsFromFile_.push_back(val);
      }
    }
  }
  
  fileIn_.close();
  
  if ( vecWidthFromFile_.size() != vecLLValsFromFile_.size() || vecWidthFromFile_.size() != vecGoodLLValsFromFile_.size() )
  {
    std::cerr << "Likelihood::ReadLLValuesFromFile : Input file " << inputFileName << " cannot be read:" << std::endl;
    std::cerr << "                                   # widths != # likelihood values !!" << std::endl;
    vecWidthFromFile_.clear();
    vecLLValsFromFile_.clear();
    vecGoodLLValsFromFile_.clear();
  }
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
  std::ofstream combFile((outputDirName_+dirNameTGraphTxt_+"output_tgraph2d_total.txt").c_str());
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
  if (calculateGoodEvtLL_)
  {
    txtOutputLogLike << std::endl << "GoodLLikelihood: ";
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      txtOutputLogLike << std::setw(12) << std::right << loglike_good_evts_[iWidth] << "  ";
    }
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
    txtOutputLogLike << std::setw(12) << std::right << loglike_data_[iWidth] << "  ";
  }
  if (calculateGoodEvtLL_)
  {
    txtOutputLogLike << std::endl << "GoodLLikelihood: ";
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      txtOutputLogLike << std::setw(12) << std::right << loglike_good_evts_data_[iWidth] << "  ";
    }
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
    txtOutputLogLike << loglike_[iWidth]/(1e+2);
    if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
  }
  txtOutputLogLike << "} *10^2;" << std::endl;
//   txtOutputLogLike << "likelihood values (data-only) : {";
//   for (int iWidth = 0; iWidth < nWidths_; iWidth++)
//   {
//     txtOutputLogLike << loglike_data_[iWidth];
//     if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
//   }
//   txtOutputLogLike << "};" << std::endl;
  if (calculateGoodEvtLL_)
  {
    txtOutputLogLike << "likelihood values (only good MC events) : {";
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      txtOutputLogLike << loglike_good_evts_[iWidth]/(1e+2);
      if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
    }
  }
  txtOutputLogLike << "} *10^2;" << std::endl;
//   txtOutputLogLike << "likelihood values (only good data events) : {";
//   for (int iWidth = 0; iWidth < nWidths_; iWidth++)
//   {
//     txtOutputLogLike << loglike_good_evts_data_[iWidth];
//     if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
//   }
//   txtOutputLogLike << "};" << std::endl;
  if (calledCMLLCalculation_)
  {
    txtOutputLogLike << "CM likelihood values : {";
    for (int iWidth = 0; iWidth < nWidths_; iWidth++)
    {
      txtOutputLogLike << loglike_CM_[iWidth]/(1e+2);
      if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
    }
    txtOutputLogLike << "} *10^2;";
    if (calculateGoodEvtLL_)
    {
      txtOutputLogLike << std::endl << "CM likelihood values (only good events) : {";
      for (int iWidth = 0; iWidth < nWidths_; iWidth++)
      {
        txtOutputLogLike << loglike_CM_good_evts_[iWidth]/(1e+2);
        if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
      }
    }
    txtOutputLogLike << "} *10^2;";
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
      txtOutputLogLike << loglike_gen_[iWidth]/(1e+2);
      if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
    }
    if (calculateGoodEvtLL_)
    {
      txtOutputLogLike << "} *10^2;" << std::endl << "PARTON likelihood values (only good matched events) : {";
      for (int iWidth = 0; iWidth < nWidths_; iWidth++)
      {
        txtOutputLogLike << loglike_gen_good_evts_[iWidth]/(1e+2);
        if ( iWidth != nWidths_-1 ) txtOutputLogLike << ", ";
      }
    }
    txtOutputLogLike << "} *10^2;";
  }
  txtOutputLogLike << std::endl << std::endl;
}
