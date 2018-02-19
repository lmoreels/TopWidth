//#include <stdio.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TLine.h>
#include <map>
#include <TArray.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TColor.h>


using namespace std;


string ConvertDoubleToString(double Number)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  convert << Number;
  return convert.str();
}

string PDotReplace(double var)
{
  std::string str = ConvertDoubleToString(var);
  std::replace(str.begin(), str.end(), 'p', '.');
  return str;
}

string PDotReplace(string str)
{
  std::replace(str.begin(), str.end(), 'p', '.');
  return str;
}

double GetBinWidth(TH1F* histo, int precision = 2)
{
  int nBins = histo->GetNbinsX();
  double min = histo->GetXaxis()->GetBinLowEdge(1);
  double max = histo->GetXaxis()->GetBinLowEdge(nBins);
  double binWidth = (max-min)/nBins;
  cout << "Bin width = " << binWidth;
  
  if      ( precision == 0 ) binWidth = floor(binWidth+0.5);  // round number to 0 decimal places
  else if ( precision == 1 ) binWidth = floor(binWidth*10.+0.5)/10.;  // round number to 1 decimal place
  else if ( precision == 2 ) binWidth = floor(binWidth*100.+0.5)/100.;  // round number to 2 decimal places
  else if ( precision == 3 ) binWidth = floor(binWidth*1000.+0.5)/1000.;  // round number to 3 decimal places
  else                      binWidth = floor(binWidth*10000.+0.5)/10000.;  // round number to 4 decimal places
  cout << " ; rounded: " << binWidth << endl;
  
  return binWidth;
}

int main()
{
  map<string,TH1F*> histo;
  Color_t listRainbow[] = {kRed, kOrange-3, kYellow-7, /*kGreen-7,*/ kGreen+1, kCyan+1,/* kBlue+2,*/ kMagenta/*, kViolet-5, kPink+10*/};
  Color_t listColour[] = {kBlue-2, kRed+1, kOrange-3, kGreen-7, kAzure+1, kViolet-9, kMagenta};
  Color_t listColour2[] = {kBlue-2, kGreen-3, kRed+1};
  Style_t listMarkers[] = {24, 26, 32, 25, 30, 27, 28, 46, 42};
  
  string suffixName = "";
  
//   string sysNames[] = {"nominal", "tune_up", "tune_down", "ISR_up", "ISR_down", "FSR_up", "FSR_down", "hdamp_up", "hdamp_down", "CR_ERD", "CR_QCD_ERD", "CR_gluon_move", "CR_GM_ERD", "herwig", "JES_up", "JES_down", "JER_up", "JER_down"};
//   int drawSys[] = {0, 3, 4, 5, 6, 7, 8}; suffixName = "_shower";
//   //int drawSys[] = {0, 1, 2, 9, 10, 11, 12}; suffixName = "_CR";
//   //int drawSys[] = {0, 13}; suffixName = "_herwig";
//   int nSys = sizeof(drawSys)/sizeof(drawSys[0]);
//   //int nSys = sizeof(sysNames)/sizeof(sysNames[0])
  
  
  std::tuple<std::string,std::string> nominal_ = std::make_tuple("nominal", "180121_2012");
  std::tuple<std::string,std::string> input_single_[] = {std::make_tuple("topPt", "180121_2022"),
       std::make_tuple("fragP", "180121_2049")};
  int nSingleSys = sizeof(input_single_)/sizeof(input_single_[0]);
  std::tuple<std::string,std::string,std::string> input_ud_[] = {std::make_tuple("JER", "180121_2155", "180121_2153"),
       std::make_tuple("JES", "180121_2159", "180121_2157"), std::make_tuple("leptonId", "180121_2014", "180121_2038"),
       std::make_tuple("leptonIso", "180121_2059", "180121_2117"),
       std::make_tuple("leptonTrig", "180121_2119", "180121_2121"),
       std::make_tuple("leptonTrk", "180121_2123", "180121_2125"),
       std::make_tuple("btag", "180121_2018", "180121_2020"), std::make_tuple("pileup", "180121_2127", "180121_2016"),
       std::make_tuple("lumi", "180121_2026", "180121_2028"), std::make_tuple("fragBL", "180121_2043", "180121_2047"),
       std::make_tuple("BRsemilep", "180121_2051", "180121_2053"),
       std::make_tuple("pdfAlphaS", "180121_2055", "180121_2222"),
       std::make_tuple("rateCM", "180121_2224", "180121_2103"),
       std::make_tuple("rateSTt", "180121_2105", "180121_2107"),
       std::make_tuple("rateSTtW", "180121_2109", "180121_2111"),
       std::make_tuple("rateOther", "180121_2113", "180121_2115"),
       std::make_tuple("UE", "180121_2129", "180121_2131"), std::make_tuple("ISR", "180121_2137", "180121_2139"),
       std::make_tuple("FSR", "180121_2141", "180121_2143"), std::make_tuple("matching", "180121_2145", "180121_2147"),
       std::make_tuple("mass", "180122_1159", "180122_1014")};
  int nUDSys = sizeof(input_ud_)/sizeof(input_ud_[0]);
  std::tuple<std::string,std::string> input_cr_[] = {std::make_tuple("gluonMove", "180121_2226"), std::make_tuple("gluonMoveERD", "180121_2135"), std::make_tuple("ERD", "180121_2149"), std::make_tuple("QCDERD", "180121_2228")};
  int nCRSys = sizeof(input_cr_)/sizeof(input_cr_[0]);
  std::tuple<std::string,std::string> input_renfac_[] = {std::make_tuple("1002", "180121_2030"), std::make_tuple("1003", "180121_2215"), std::make_tuple("1004", "180121_2217"), std::make_tuple("1005", "180121_2220"), std::make_tuple("1007", "180121_2039"), std::make_tuple("1009", "180121_2041")};
  int nRenfacSys = sizeof(input_renfac_)/sizeof(input_renfac_[0]);
  
  //std::get<0>(input_);
  
  //string pathFiles = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/OutputPlots/mu/"+date;
  string pathFiles = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/OutputPlots/mu/";
  string pathOutput = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/Systematics/Plots/";
  
  TFile *fileOut = new TFile((pathOutput+"Systematics_plots"+suffixName+".root").c_str(),"RECREATE");
  
  string fileName = "NtuplePlots_nominal.root";
  
  string varNames[] = {"red_top_mass", "red_mlb_mass"};
  string axisLabels[] = {"m_{r}", "m_{lb,r}"};
  int nVars = sizeof(varNames)/sizeof(varNames[0]);
  
  string prefixes[] = {"allSim_", "allData_"};
  int nPrefixes = sizeof(prefixes)/sizeof(prefixes[0]);
  
  
  string histoName = "", legendName = "";
  string sysName = "";
  string binWidthStr = "";
  
  /// Get nominal histograms
  TFile *fileNom = new TFile((pathFiles+std::get<1>(nominal_)+"/"+fileName).c_str(),"read");
  fileNom->cd();
  for (int iVar = 0; iVar < nVars; iVar++)
  {
    histo[varNames[iVar]+"_data"] =  (TH1F*) fileNom->Get(("1D_histograms/allData_"+varNames[iVar]).c_str())->Clone((varNames[iVar]+"_data").c_str());
    histo[varNames[iVar]+"_data"]->SetDirectory(0);
    histo[varNames[iVar]+"_nominal_MC"] =  (TH1F*) fileNom->Get(("1D_histograms/allSim_"+varNames[iVar]).c_str())->Clone((varNames[iVar]+"_nominal_MC").c_str());
    histo[varNames[iVar]+"_nominal_MC"]->SetDirectory(0);
    /// TEMPORARILY
    histo[varNames[iVar]+"_nominal_MC"]->Scale(0.962);
  }
  //fileNom->Close();
  
  
  /// Get histograms for up/down systematics
  TFile *fileInUp, *fileInDown;
  for (int iSys = 0; iSys < nUDSys; iSys++)
  {
    sysName = std::get<0>(input_ud_[iSys]);
//     if ( sysNames[drawSys[iSys]].find("JES_up") != std::string::npos ) fileName = "NtuplePlots_JESup.root";
//     else if ( sysNames[drawSys[iSys]].find("JES_down") != std::string::npos ) fileName = "NtuplePlots_JESdown.root";
//     else if ( sysNames[drawSys[iSys]].find("JER_up") != std::string::npos ) fileName = "NtuplePlots_JERup.root";
//     else if ( sysNames[drawSys[iSys]].find("JER_down") != std::string::npos ) fileName = "NtuplePlots_JERdown.root";
//     else fileName = "NtuplePlots_nominal.root";
    if ( sysName.find("JES") != std::string::npos ) fileName = "NtuplePlots_JESup.root";
    else if ( sysName.find("JER") != std::string::npos ) fileName = "NtuplePlots_JERup.root";
    else fileName = "NtuplePlots_nominal.root";
    
    
    fileInUp = new TFile((pathFiles+std::get<1>(input_ud_[iSys])+"/"+fileName).c_str(),"read");
    fileInUp->cd();
    for (int iVar = 0; iVar < nVars; iVar++)
    {
      histoName = varNames[iVar]+"_"+sysName+"_up";
      histo[histoName] = (TH1F*) fileInUp->Get(("1D_histograms/allSim_"+varNames[iVar]).c_str())->Clone(histoName.c_str());
      histo[histoName]->SetDirectory(0);  // keep histogram when file is closed
      /// TEMPORARILY
      histo[histoName]->Scale(0.962);
    }
    //fileInUp->Close();
    
    if ( sysName.find("JES") != std::string::npos ) fileName = "NtuplePlots_JESdown.root";
    else if ( sysName.find("JER") != std::string::npos ) fileName = "NtuplePlots_JERdown.root";
    else fileName = "NtuplePlots_nominal.root";
    
    fileInDown = new TFile((pathFiles+std::get<2>(input_ud_[iSys])+"/"+fileName).c_str(),"read");
    fileInDown->cd();
    for (int iVar = 0; iVar < nVars; iVar++)
    {
      histoName = varNames[iVar]+"_"+sysName+"_down";
      histo[histoName] = (TH1F*) fileInDown->Get(("1D_histograms/allSim_"+varNames[iVar]).c_str())->Clone(histoName.c_str());
      histo[histoName]->SetDirectory(0);  // keep histogram when file is closed
      /// TEMPORARILY
      histo[histoName]->Scale(0.962);
    }
    //fileInDown->Close();
    
    
    string suffixes[] = {"_nominal_MC", "_"+sysName+"_up", "_"+sysName+"_down", "_data"};
    string suffixNames[] = {"nominal", "up", "down", "data"};
    int nSuffixes = sizeof(suffixes)/sizeof(suffixes[0]);
    
    /// Make plots
    fileOut->cd();
    
    for (int iVar = 0; iVar < nVars; iVar++)
    {
      // Fix scale of FSR systematic
      if ( sysName.find("FSR") != std::string::npos )
      {
        double factor = 1./sqrt(2.);
        double nBins = histo[varNames[iVar]+suffixes[0]]->GetNbinsX();
        double binNom, binSys, binNew;
        for (int iSys = 1; iSys < 3; iSys++)
        {
          for (int iBin = 1; iBin < nBins+1; iBin++)
          {
            binNom = histo[varNames[iVar]+suffixes[0]]->GetBinContent(iBin);
            binSys = histo[varNames[iVar]+suffixes[iSys]]->GetBinContent(iBin);
            binNew = binNom - (binNom - binSys)*factor;
            histo[varNames[iVar]+suffixes[iSys]]->SetBinContent(iBin, binNew);
          }
        }
        
      }
      
      // Draw plots
      for (int iType = 0; iType < 2; iType++)  // Draw full/normalised
      {
        TCanvas *c1;
        if ( iType == 0 ) c1 = new TCanvas((varNames[iVar]+"_"+sysName+"_overlay").c_str(),(varNames[iVar]+"_"+sysName+"_overlay").c_str());
        else c1 = new TCanvas((varNames[iVar]+"_"+sysName+"_norm_overlay").c_str(),(varNames[iVar]+"_"+sysName+"_norm_overlay").c_str());
        c1->cd();
        
        
        /// Make legend
        TLegend *leg = new TLegend(0.65,0.55,0.89,0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        
        for (int iSuffix = 0; iSuffix < nSuffixes-1; iSuffix++)
        {
          histoName = varNames[iVar]+suffixes[iSuffix];
          histo[histoName]->SetLineWidth(2);
          histo[histoName]->SetLineColor(listColour2[iSuffix]);
          //histo[histoName]->SetMarkerColor(listColour2[iSuffix]);
          //histo[histoName]->SetMarkerStyle(listMarkers[iSuffix%8]);
          
          
          if ( iSuffix == 0 )
          {
            binWidthStr = ConvertDoubleToString(GetBinWidth(histo[histoName], 2));
            legendName = suffixNames[iSuffix];
            
            histo[histoName]->SetTitle("");
            histo[histoName]->GetXaxis()->SetTitle((axisLabels[iVar]).c_str());
            //histo[histoName]->GetXaxis()->SetRangeUser(0.,2.1.);
            if ( iType == 0 ) histo[histoName]->GetYaxis()->SetTitle(("Events / "+binWidthStr+" units").c_str());
            else histo[histoName]->GetYaxis()->SetTitle(("Normalised events / "+binWidthStr+" units").c_str());
            histo[histoName]->GetYaxis()->SetTitleOffset(1.4);
            if ( sysName.find("FSR") != std::string::npos )
              histo[histoName]->SetMaximum(1.05*histo[varNames[iVar]+suffixes[2]]->GetMaximum());
            histo[histoName]->SetStats(0);
            if ( iType == 0 ) histo[histoName]->Draw("e hist");
            else histo[histoName]->DrawNormalized("e hist");
          }
          else
          {
            legendName = sysName+" "+suffixNames[iSuffix];
            if ( iType == 0 ) histo[histoName]->Draw("e hist same");
            else histo[histoName]->DrawNormalized("e hist same");
          }
          leg->AddEntry(histo[histoName],(legendName).c_str(),"lp");
          c1->Update();
          
        }  // end suffixes
        
        // Draw nominal on top, no legend entry
        histoName = varNames[iVar]+suffixes[0];
        histo[histoName]->SetLineWidth(2);
        histo[histoName]->SetLineColor(listColour2[0]);
        histo[histoName]->SetMarkerColor(listColour2[0]);
        //histo[histoName]->SetMarkerStyle(listMarkers[0]);
        if ( iType == 0 ) histo[histoName]->Draw("e hist same");
        else histo[histoName]->DrawNormalized("e hist same");
        
        // Draw data with legend entry
        histoName = varNames[iVar]+suffixes[nSuffixes-1];
        legendName = suffixNames[nSuffixes-1];
        histo[histoName]->SetLineColor(kBlack);
        histo[histoName]->SetMarkerColor(kBlack);
        histo[histoName]->SetMarkerStyle(20);
        if ( iType == 0 ) histo[histoName]->Draw("E same");
        else histo[histoName]->DrawNormalized("E same");
        leg->AddEntry(histo[histoName],(legendName).c_str(),"lp");
        c1->Update();
        
        leg->Draw();
        c1->Update();
        c1->Write();
        if ( iType == 0 ) c1->SaveAs((pathOutput+"Syst_overlay_"+varNames[iVar]+"_"+sysName+".png").c_str());
        else c1->SaveAs((pathOutput+"Syst_normalised_overlay_"+varNames[iVar]+"_"+sysName+".png").c_str());
        
      }  // end full/normalised
      
    }  // end vars
    
    fileInUp->Close();
    fileInDown->Close();
    
  }  // end up/down sys
  
  
  /// Get histograms for single systematics
  TFile *fileIn;
  for (int iSys = 0; iSys < nSingleSys; iSys++)
  {
    sysName = std::get<0>(input_single_[iSys]);
    
    fileIn = new TFile((pathFiles+std::get<1>(input_single_[iSys])+"/"+fileName).c_str(),"read");
    fileIn->cd();
    for (int iVar = 0; iVar < nVars; iVar++)
    {
      histoName = varNames[iVar]+"_"+sysName;
      histo[histoName] = (TH1F*) fileIn->Get(("1D_histograms/allSim_"+varNames[iVar]).c_str())->Clone(histoName.c_str());
      histo[histoName]->SetDirectory(0);  // keep histogram when file is closed
      /// TEMPORARILY
      histo[histoName]->Scale(0.962);
    }
    
    
    string suffixes[] = {"_nominal_MC", "_"+sysName, "_data"};
    string suffixNames[] = {"nominal", sysName, "data"};
    int nSuffixes = sizeof(suffixes)/sizeof(suffixes[0]);
    
    /// Make plots
    fileOut->cd();
    
    for (int iVar = 0; iVar < nVars; iVar++)
    {
      for (int iType = 0; iType < 2; iType++)  // Draw full/normalised
      {
        TCanvas *c1;
        if ( iType == 0 ) c1 = new TCanvas((varNames[iVar]+"_"+sysName+"_overlay").c_str(),(varNames[iVar]+"_"+sysName+"_overlay").c_str());
        else c1 = new TCanvas((varNames[iVar]+"_"+sysName+"_norm_overlay").c_str(),(varNames[iVar]+"_"+sysName+"_norm_overlay").c_str());
        c1->cd();
        
        
        /// Make legend
        TLegend *leg = new TLegend(0.65,0.55,0.89,0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        
        for (int iSuffix = 0; iSuffix < nSuffixes-1; iSuffix++)
        {
          histoName = varNames[iVar]+suffixes[iSuffix];
          histo[histoName]->SetLineWidth(2);
          histo[histoName]->SetLineColor(listColour2[iSuffix]);
          histo[histoName]->SetMarkerColor(listColour2[iSuffix]);
          //histo[histoName]->SetMarkerStyle(listMarkers[iSuffix%8]);
          
          
          if ( iSuffix == 0 )
          {
            binWidthStr = ConvertDoubleToString(GetBinWidth(histo[histoName], 2));

            histo[histoName]->SetTitle("");
            histo[histoName]->GetXaxis()->SetTitle((axisLabels[iVar]).c_str());
            //histo[histoName]->GetXaxis()->SetRangeUser(0.,2.1.);
            if ( iType == 0 ) histo[histoName]->GetYaxis()->SetTitle(("Events / "+binWidthStr+" units").c_str());
            else histo[histoName]->GetYaxis()->SetTitle(("Normalised events / "+binWidthStr+" units").c_str());
            histo[histoName]->GetYaxis()->SetTitleOffset(1.4);
            histo[histoName]->SetStats(0);
            if ( iType == 0 ) histo[histoName]->Draw("e hist");
            else histo[histoName]->DrawNormalized("e hist");
          }
          else
          {
            legendName = sysName+" "+suffixNames[iSuffix];
            if ( iType == 0 ) histo[histoName]->Draw("e hist same");
            else histo[histoName]->DrawNormalized("e hist same");
          }
          leg->AddEntry(histo[histoName],(suffixNames[iSuffix]).c_str(),"lp");
          c1->Update();
          
        }  // end suffixes
        
        // Draw normalised on top, no legend entry
        histoName = varNames[iVar]+suffixes[0];
        histo[histoName]->SetLineWidth(2);
        histo[histoName]->SetLineColor(listColour2[0]);
        histo[histoName]->SetMarkerColor(listColour2[0]);
        //histo[histoName]->SetMarkerStyle(listMarkers[0]);
        if ( iType == 0 ) histo[histoName]->Draw("e hist same");
        else histo[histoName]->DrawNormalized("e hist same");
        
        // Draw data with legend entry
        histoName = varNames[iVar]+suffixes[nSuffixes-1];
        legendName = suffixNames[nSuffixes-1];
        histo[histoName]->SetLineColor(kBlack);
        histo[histoName]->SetMarkerColor(kBlack);
        histo[histoName]->SetMarkerStyle(20);
        if ( iType == 0 ) histo[histoName]->Draw("E same");
        else histo[histoName]->DrawNormalized("E same");
        leg->AddEntry(histo[histoName],(suffixNames[nSuffixes-1]).c_str(),"lp");
        c1->Update();
        
        leg->Draw();
        c1->Update();
        c1->Write();
        if ( iType == 0 ) c1->SaveAs((pathOutput+"Syst_overlay_"+varNames[iVar]+"_"+sysName+".png").c_str());
        else c1->SaveAs((pathOutput+"Syst_normalised_overlay_"+varNames[iVar]+"_"+sysName+".png").c_str());
        
      }  // end full/normalised
      
    }  // end vars
    
    fileIn->Close();
    
  }  // end single sys
  
  
  /// Get histograms for colour reconnection
  for (int iSys = 0; iSys < nCRSys; iSys++)
  {
    sysName = std::get<0>(input_cr_[iSys]);
    //sysName = "cr";
    
    fileIn = new TFile((pathFiles+std::get<1>(input_cr_[iSys])+"/"+fileName).c_str(),"read");
    fileIn->cd();
    for (int iVar = 0; iVar < nVars; iVar++)
    {
      histoName = varNames[iVar]+"_"+sysName;
      histo[histoName] = (TH1F*) fileIn->Get(("1D_histograms/allSim_"+varNames[iVar]).c_str())->Clone(histoName.c_str());
      histo[histoName]->SetDirectory(0);  // keep histogram when file is closed
      /// TEMPORARILY
      histo[histoName]->Scale(0.962);
    }
    
    fileIn->Close();
  }
  
  
  string suffixes[] = {"_nominal_MC", "_data"};
  string suffixNames[] = {"nominal", "data"};
  int nSuffixes = sizeof(suffixes)/sizeof(suffixes[0]);
  
  /// Make plots
  fileOut->cd();
  
  for (int iVar = 0; iVar < nVars; iVar++)
  {
    for (int iType = 0; iType < 2; iType++)  // Draw full/normalised
    {
      TCanvas *c1;
      if ( iType == 0 ) c1 = new TCanvas((varNames[iVar]+"_cr_overlay").c_str(),(varNames[iVar]+"_cr_overlay").c_str());
      else c1 = new TCanvas((varNames[iVar]+"_cr_norm_overlay").c_str(),(varNames[iVar]+"_cr_norm_overlay").c_str());
      c1->cd();
      
      
      /// Make legend
      TLegend *leg = new TLegend(0.65,0.55,0.89,0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.04);
      
      // Draw nominal
      histoName = varNames[iVar]+suffixes[0];
      histo[histoName]->SetLineWidth(2);
      histo[histoName]->SetLineColor(listColour[0]);
      histo[histoName]->SetMarkerColor(listColour[0]);
      
      binWidthStr = ConvertDoubleToString(GetBinWidth(histo[histoName], 2));
      histo[histoName]->SetTitle("");
      histo[histoName]->GetXaxis()->SetTitle((axisLabels[iVar]).c_str());
      if ( iType == 0 ) histo[histoName]->GetYaxis()->SetTitle(("Events / "+binWidthStr+" units").c_str());
      else histo[histoName]->GetYaxis()->SetTitle(("Normalised events / "+binWidthStr+" units").c_str());
      histo[histoName]->GetYaxis()->SetTitleOffset(1.4);
      histo[histoName]->SetStats(0);
      if ( iType == 0 ) histo[histoName]->Draw("e hist");
      else histo[histoName]->DrawNormalized("e hist");
      leg->AddEntry(histo[histoName],(suffixNames[0]).c_str(),"lp");
      c1->Update();
      
      // Draw CRs
      for (int iSys = 0; iSys < nCRSys; iSys++)
      {
        sysName = std::get<0>(input_cr_[iSys]);
        histoName = varNames[iVar]+"_"+sysName;
        histo[histoName]->SetLineWidth(2);
        histo[histoName]->SetLineColor(listColour[iSys+1]);
        histo[histoName]->SetMarkerColor(listColour[iSys+1]);
        //histo[histoName]->SetMarkerStyle(listMarkers[(iSys+1)%8]);
        if ( iType == 0 ) histo[histoName]->Draw("e hist same");
        else histo[histoName]->DrawNormalized("e hist same");
        leg->AddEntry(histo[histoName],(sysName).c_str(),"lp");
        c1->Update();
      }  // end sys
      
      // Draw normalised on top, no legend entry
      histoName = varNames[iVar]+suffixes[0];
      histo[histoName]->SetLineWidth(2);
      histo[histoName]->SetLineColor(listColour[0]);
      histo[histoName]->SetMarkerColor(listColour[0]);
      //histo[histoName]->SetMarkerStyle(listMarkers[0]);
      if ( iType == 0 ) histo[histoName]->Draw("e hist same");
      else histo[histoName]->DrawNormalized("e hist same");
      
      // Draw data with legend entry
      histoName = varNames[iVar]+suffixes[nSuffixes-1];
      legendName = suffixNames[nSuffixes-1];
      histo[histoName]->SetLineColor(kBlack);
      histo[histoName]->SetMarkerColor(kBlack);
      histo[histoName]->SetMarkerStyle(20);
      if ( iType == 0 ) histo[histoName]->Draw("E same");
      else histo[histoName]->DrawNormalized("E same");
      leg->AddEntry(histo[histoName],(suffixNames[nSuffixes-1]).c_str(),"lp");
      c1->Update();
      
      leg->Draw();
      c1->Update();
      c1->Write();
      if ( iType == 0 ) c1->SaveAs((pathOutput+"Syst_overlay_"+varNames[iVar]+"_cr.png").c_str());
      else c1->SaveAs((pathOutput+"Syst_normalised_overlay_"+varNames[iVar]+"_cr.png").c_str());
      
    }  // end full/normalised
      
  }  // end vars
  
  // end cr
  
  
  /// Get histograms for renormalisation/factorisation scale
  for (int iSys = 0; iSys < nRenfacSys; iSys++)
  {
    sysName = std::get<0>(input_renfac_[iSys]);
    
    fileIn = new TFile((pathFiles+std::get<1>(input_renfac_[iSys])+"/"+fileName).c_str(),"read");
    fileIn->cd();
    for (int iVar = 0; iVar < nVars; iVar++)
    {
      histoName = varNames[iVar]+"_"+sysName;
      histo[histoName] = (TH1F*) fileIn->Get(("1D_histograms/allSim_"+varNames[iVar]).c_str())->Clone(histoName.c_str());
      histo[histoName]->SetDirectory(0);  // keep histogram when file is closed
      /// TEMPORARILY
      histo[histoName]->Scale(0.962);
    }
    
    fileIn->Close();
  }
  
  /// Make plots
  fileOut->cd();
  
  for (int iVar = 0; iVar < nVars; iVar++)
  {
    for (int iType = 0; iType < 2; iType++)  // Draw full/normalised
    {
      TCanvas *c1;
      if ( iType == 0 ) c1 = new TCanvas((varNames[iVar]+"_renfac_overlay").c_str(),(varNames[iVar]+"_renfac_overlay").c_str());
      else c1 = new TCanvas((varNames[iVar]+"_renfac_norm_overlay").c_str(),(varNames[iVar]+"_renfac_norm_overlay").c_str());
      c1->cd();
      
      
      /// Make legend
      TLegend *leg = new TLegend(0.65,0.55,0.89,0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.04);
      
      // Draw nominal
      histoName = varNames[iVar]+suffixes[0];
      histo[histoName]->SetLineWidth(2);
      histo[histoName]->SetLineColor(listColour[0]);
      histo[histoName]->SetMarkerColor(listColour[0]);
      
      binWidthStr = ConvertDoubleToString(GetBinWidth(histo[histoName], 2));
      histo[histoName]->SetTitle("");
      histo[histoName]->GetXaxis()->SetTitle((axisLabels[iVar]).c_str());
      if ( iType == 0 ) histo[histoName]->GetYaxis()->SetTitle(("Events / "+binWidthStr+" units").c_str());
      else histo[histoName]->GetYaxis()->SetTitle(("Normalised events / "+binWidthStr+" units").c_str());
      histo[histoName]->GetYaxis()->SetTitleOffset(1.4);
      histo[histoName]->SetStats(0);
      if ( iType == 0 ) histo[histoName]->Draw("e hist");
      else histo[histoName]->DrawNormalized("e hist");
      leg->AddEntry(histo[histoName],(suffixNames[0]).c_str(),"lp");
      c1->Update();
      
      // Draw CRs
      for (int iSys = 0; iSys < nRenfacSys; iSys++)
      {
        sysName = std::get<0>(input_renfac_[iSys]);
        histoName = varNames[iVar]+"_"+sysName;
        histo[histoName]->SetLineWidth(2);
        histo[histoName]->SetLineColor(listColour[iSys+1]);
        histo[histoName]->SetMarkerColor(listColour[iSys+1]);
        //histo[histoName]->SetMarkerStyle(listMarkers[(iSys+1)%8]);
        if ( iType == 0 ) histo[histoName]->Draw("e hist same");
        else histo[histoName]->DrawNormalized("e hist same");
        leg->AddEntry(histo[histoName],(sysName).c_str(),"lp");
        c1->Update();
      }  // end sys
      
      // Draw normalised on top, no legend entry
      histoName = varNames[iVar]+suffixes[0];
      histo[histoName]->SetLineWidth(2);
      histo[histoName]->SetLineColor(listColour[0]);
      histo[histoName]->SetMarkerColor(listColour[0]);
      //histo[histoName]->SetMarkerStyle(listMarkers[0]);
      if ( iType == 0 ) histo[histoName]->Draw("e hist same");
      else histo[histoName]->DrawNormalized("e hist same");
      
      // Draw data with legend entry
      histoName = varNames[iVar]+suffixes[nSuffixes-1];
      legendName = suffixNames[nSuffixes-1];
      histo[histoName]->SetLineColor(kBlack);
      histo[histoName]->SetMarkerColor(kBlack);
      histo[histoName]->SetMarkerStyle(20);
      if ( iType == 0 ) histo[histoName]->Draw("E same");
      else histo[histoName]->DrawNormalized("E same");
      leg->AddEntry(histo[histoName],(suffixNames[nSuffixes-1]).c_str(),"lp");
      c1->Update();
      
      leg->Draw();
      c1->Update();
      c1->Write();
      if ( iType == 0 ) c1->SaveAs((pathOutput+"Syst_overlay_"+varNames[iVar]+"_renfac.png").c_str());
      else c1->SaveAs((pathOutput+"Syst_normalised_overlay_"+varNames[iVar]+"_renfac.png").c_str());
      
    }  // end full/normalised
    
  }  // end vars
  
  // end ren/fac
  
  
  fileOut->Close();
  fileNom->Close();
  
  return 0;
}
