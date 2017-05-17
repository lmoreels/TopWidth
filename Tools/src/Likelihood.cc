#include "../interface/Likelihood.h"

const double Likelihood::widthArray_[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.};

const std::string Likelihood::listCats_[] = {"CP", "WP", "UP"};

const int Likelihood::nWidths_ = sizeof(widthArray_)/sizeof(widthArray_[0]);
const int Likelihood::nCats_ = sizeof(listCats_)/sizeof(listCats_[0]);

bool Likelihood::fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile.good();
}

std::string Likelihood::ConvertDoubleToString(double Number)
{
  std::ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

std::string Likelihood::DotReplace(double var)
{
  std::string str = ConvertDoubleToString(var);
  std::replace(str.begin(), str.end(), '.', 'p');
  return str;
}

Likelihood::Likelihood(double min, double max, std::string outputDirName, std::string date, bool verbose):
verbose_(verbose),  date_(date), outputDirName_(outputDirName), inputFileName_(""), suffix_(""), histoName_(""), minRedMass_(min), maxRedMass_(max), vecWidthStr_(), histo_(), histoSm_(), histoTotal_(), graph_(), vecBinCentres_(), vecBinContents_()
{
  rew = new EventReweighting(false);  // no correction for number of events
  
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    vecWidthStr_.push_back(DotReplace(widthArray_[iWidth]));
  }
  
  this->BookHistograms();
}

Likelihood::~Likelihood()
{
  delete file_;
  delete rew;
}

void Likelihood::BookHistograms()
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    thisWidth_ = vecWidthStr_[iWidth];
    
    histo_[("Red_top_mass_CP_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_CP_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", CP; M_{t}/<M_{t}>").c_str(), 90, 0.5, 2.0);
    histo_[("Red_top_mass_CP_widthx"+thisWidth_+"_100b").c_str()] = new TH1D(("Red_top_mass_CP_widthx"+thisWidth_+"_100b").c_str(),("Reduced top mass for width "+thisWidth_+", CP; M_{t}/<M_{t}>").c_str(), 100, 0.5, 2.0);
    histo_[("Red_top_mass_CP_widthx"+thisWidth_+"_120b").c_str()] = new TH1D(("Red_top_mass_CP_widthx"+thisWidth_+"_120b").c_str(),("Reduced top mass for width "+thisWidth_+", CP; M_{t}/<M_{t}>").c_str(), 120, 0.5, 2.0);
    
    histo_[("Red_top_mass_WP_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_WP_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", WP; M_{t}/<M_{t}>").c_str(), 90, 0.5, 2.0);
    histo_[("Red_top_mass_WP_widthx"+thisWidth_+"_100b").c_str()] = new TH1D(("Red_top_mass_WP_widthx"+thisWidth_+"_100b").c_str(),("Reduced top mass for width "+thisWidth_+", WP; M_{t}/<M_{t}>").c_str(), 100, 0.5, 2.0);
    histo_[("Red_top_mass_WP_widthx"+thisWidth_+"_120b").c_str()] = new TH1D(("Red_top_mass_WP_widthx"+thisWidth_+"_120b").c_str(),("Reduced top mass for width "+thisWidth_+", WP; M_{t}/<M_{t}>").c_str(), 120, 0.5, 2.0);
    
    histo_[("Red_top_mass_UP_widthx"+thisWidth_+"_90b").c_str()] = new TH1D(("Red_top_mass_UP_widthx"+thisWidth_+"_90b").c_str(),("Reduced top mass for width "+thisWidth_+", UP; M_{t}/<M_{t}>").c_str(), 90, 0.5, 2.0);
    histo_[("Red_top_mass_UP_widthx"+thisWidth_+"_100b").c_str()] = new TH1D(("Red_top_mass_UP_widthx"+thisWidth_+"_100b").c_str(),("Reduced top mass for width "+thisWidth_+", UP; M_{t}/<M_{t}>").c_str(), 100, 0.5, 2.0);
    histo_[("Red_top_mass_UP_widthx"+thisWidth_+"_120b").c_str()] = new TH1D(("Red_top_mass_UP_widthx"+thisWidth_+"_120b").c_str(),("Reduced top mass for width "+thisWidth_+", UP; M_{t}/<M_{t}>").c_str(), 120, 0.5, 2.0);
  }
}

void Likelihood::FillHistograms(double redMass, double massForWidthSF, bool isTTbar, std::string catSuffix)
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    thisWidth_ = vecWidthStr_[iWidth];
    if (isTTbar) thisWidthSF_ = rew->EventWeightCalculator(massForWidthSF, widthArray_[iWidth]);
    else thisWidthSF_ = 1.;
    
    histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_90b").c_str()]->Fill(redMass, thisWidthSF_);
    histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_100b").c_str()]->Fill(redMass, thisWidthSF_);
    histo_[("Red_top_mass"+catSuffix+"_widthx"+thisWidth_+"_120b").c_str()]->Fill(redMass, thisWidthSF_);
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
    
    suffix_ = "widthx"+DotReplace(widthArray_[iWidth]);
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
        std::cerr << "Histograms should have the same number of bins! Exiting..." << std::endl;
        exit(1);
      }
      
      binMin = histoSm_[histoName_]->FindBin(0.6);
      binMax = histoSm_[histoName_]->FindBin(1.4)+1;
      for (int iBin = binMin; iBin < binMax+1; iBin++) { nEvents[iCat] += histoSm_[histoName_]->GetBinContent(iBin);}
      totEvents += nEvents[iCat];
      
      if (verbose_) std::cout << histoName_ << ": " << nEvents[iCat] << " events" << std::endl;
      
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
    
    if (verbose_) std::cout << "Total: " << totEvents << " events" << std::endl;
    
    /// Make arrays as input for TGraph
    const int nPoints = (vecBinCentres_[listCats_[0]+"_"+suffix_]).size();
    double binCentreArray[nPoints] = {0.}, binContentArray[nCats_][nPoints] = {0.}, totalBinContentArray[nPoints] = {0.};
    for (int i = 0; i < nPoints; i++)
    {
      binCentreArray[i] = (vecBinCentres_[listCats_[0]+"_"+suffix_]).at(i);
      for (int iCat = 0; iCat < nCats_; iCat++)
      {
        binContentArray[iCat][i] = (vecBinContents_[listCats_[iCat]+"_"+suffix_]).at(i);
        // Calculate event fractions (one time)
        if ( i == 0 ) { fracCats[iCat] = (double)nEvents[iCat]/(double)totEvents;}
      }
    }
    
    /// Make TGraphs
    for (int iCat = 0; iCat < nCats_; iCat++)
    {
      histoName_ = listCats_[iCat]+suffix_;
      
      graph_[histoName_] = new TGraph(nPoints, binCentreArray, binContentArray[iCat]);
      graph_[histoName_]->SetName(histoName_.c_str());
      graph_[histoName_]->SetTitle(histoName_.c_str());
      graph_[histoName_]->Write();
      DrawGraph(histoSm_[histoName_], graph_[histoName_], "Graph_Red_top_mass_"+histoName_);
      if ( iCat == 0 )  //CP
        WriteOutput(nPoints, iWidth, binCentreArray, binContentArray[iCat], "CorrectMatchLikelihood", 1);
      
      
      /// Calculate total probability
      fracCats[iCat] = (double)nEvents[iCat]/(double)totEvents;
      for (int i = 0; i < nPoints; i++) totalBinContentArray[i] += fracCats[iCat]*binContentArray[iCat][i];
      
      histoSm_[histoName_]->Scale(fracCats[iCat]);
      if ( iCat == 0 )
      {
        histoTotal_[suffix_] = (TH1D*) histoSm_[histoName_]->Clone("TotalProbability");
        histoTotal_[suffix_]->SetTitle("#frac{n_{CP}}{n_{tot}} * f_{CP}(x|#Gamma) + #frac{n_{WP}}{n_{tot}} * f_{WP}(x|#Gamma) + #frac{n_{UP}}{n_{tot}} * f_{UP}(x|#Gamma)");
      }
      else histoTotal_[suffix_]->Add(histoSm_[histoName_]);
    }
    
    if (verbose_) std::cout << "The integral of the weighted probability histogram is " << histoTotal_[suffix_]->Integral(binMin, binMax) << std::endl;
    
    graph_["TotalProbability_"+suffix_] = new TGraph(nPoints, binCentreArray, totalBinContentArray);
    graph_["TotalProbability_"+suffix_]->SetName(("TotalProbability_"+suffix_).c_str());
    graph_["TotalProbability_"+suffix_]->SetTitle(("TotalProbability_"+suffix_).c_str());
    graph_["TotalProbability_"+suffix_]->Write();
    DrawGraph(histoTotal_[suffix_], graph_["TotalProbability_"+suffix_], "Graph_totalProbability_"+suffix_);
    
    
    /// Negative log(likelihood)
    double likelihoodArray[nPoints] = {0.};
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

void Likelihood::ConstructTGraphsFromFile(std::string name)
{
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    suffix_ = name+"widthx"+DotReplace(widthArray_[iWidth]);
    inputFileName_ = "TGraphTemplates/output_tgraph1d_"+suffix_+".txt";
    if ( ! fexists(inputFileName_.c_str()) )
    {
      std::cerr << "ERROR::ConstructTGraphs : File " << inputFileName_ << " not found!!" << std::endl;
      //std::cerr << "                          Aborting the likelihood calculation..." << std::endl;
      //calculateLikelihood = false;
      break;
    }
    graph_[suffix_] = new TGraph(inputFileName_.c_str());
  }
  //if (calculateLikelihood)
    std::cout << "Constructed TGraphs for likelihood measurements (using " << nWidths_ << " widths)" << std::endl;
}

int Likelihood::CalculateOutputWidth(double *evalWidths, double *LLvalues, std::string inputType)
{
  if ( sizeof(evalWidths)/sizeof(evalWidths[0]) != sizeof(LLvalues)/sizeof(LLvalues[0]) )
  {
    std::cerr << "Likelihood array has not the same size as widths array. Will not calculate output width..." << std::endl;
    return 1;
  }
  int N = sizeof(evalWidths)/sizeof(evalWidths[0]);
  if ( N == 0 )
  {
    std::cerr << "Arrays are empty. Will not calculate output width..." << std::endl;
    return 1;
  }
  TGraph *g = new TGraph(N, evalWidths, LLvalues);
  
  double centreVal = evalWidths[TMath::LocMin(N, LLvalues)];  // CHECK !
  double fitmin = centreVal - 0.8;
  double fitmax = centreVal + 0.8;
  TF1 *parabola = new TF1("parabola", "pol2", fitmin, fitmax);
  
  g->Fit(parabola,"R");
  
  double outputWidth = parabola->GetMinimumX(fitmin, fitmax);
  double lowerSigma = parabola->GetX(parabola->Eval(outputWidth) + 0.5, fitmin, outputWidth);
  double upperSigma = parabola->GetX(parabola->Eval(outputWidth) + 0.5, outputWidth, fitmax);
  
  return 0;
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
  
  delete c;
}

void Likelihood::DrawGraph(TGraph2D* g, std::string name)
{
  TCanvas *c = new TCanvas(name.c_str(), name.c_str());
  c->cd();
  gStyle->SetPalette(1);
  g->SetMaxIter(500000);
  g->SetName(name.c_str());
  g->SetTitle((name+"; m_{t}/<m_{t}> [GeV]; #Gamma/#Gamma_{SM}").c_str());
  //g->GetHistogram("empty")->GetXaxis()->SetTitleOffset(1.5);  // Does not work
  //g->GetHistogram("empty")->GetYaxis()->SetTitleOffset(1.5);  // Does not work
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
  graph_["likelihood_widthx"+vecWidthStr_[0]]->SetLineColor(colours[0]);
  graph_["likelihood_widthx"+vecWidthStr_[0]]->Draw();
  std::string temp = graph_["likelihood_widthx"+vecWidthStr_[0]]->GetTitle();
  graph_["likelihood_widthx"+vecWidthStr_[0]]->SetTitle("-LogLikelihood");
  c2->Update();
  for (int i = 1; i < nWidths_; i++)
  {
    graph_["likelihood_widthx"+vecWidthStr_[i]]->SetLineColor(colours[i%10]);
    graph_["likelihood_widthx"+vecWidthStr_[i]]->Draw("same");
    c2->Update();
  }
  c2->Write();
  c2->SaveAs((outputDirName_+"LogLikelihood_multi.png").c_str());
  c2->Close();
  
  graph_["likelihood_widthx"+vecWidthStr_[0]]->SetTitle(temp.c_str());
  
  delete c2;
}

void Likelihood::WriteOutput(int nPoints, double width, double *arrayCentre, double *arrayContent, std::string name, int dim)
{
  std::string dimension = "";
  if ( dim == 1 ) dimension = "1d";
  else if ( dim == 2 ) dimension = "2d";
  std::string outputTxtName = outputDirName_+"output_tgraph"+dimension+"_"+name+".txt";
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
  std::ofstream combFile((outputDirName_+"output_tgraph2d_total.txt").c_str());
  for (int iWidth = 0; iWidth < nWidths_; iWidth++)
  {
    std::string inFileName = outputDirName_+"output_tgraph2d_widthx"+vecWidthStr_[iWidth]+".txt";
    std::ifstream infile(inFileName.c_str());
    combFile << infile.rdbuf();
    infile.close();
  }
  combFile.close();
}
