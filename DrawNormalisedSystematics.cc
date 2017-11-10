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

int main()
{
  map<string,TH1F*> histo;
  Color_t listRainbow[] = {kRed, kOrange-3, kYellow-7, /*kGreen-7,*/ kGreen+1, kCyan+1,/* kBlue+2,*/ kMagenta/*, kViolet-5, kPink+10*/};
  Style_t listMarkers[] = {24, 26, 32, 25, 30, 27, 28, 46, 42};
  
  string sysNames[] = {"nominal", "tune_up", "tune_down", "ISR_up", "ISR_down", "FSR_up", "FSR_down", "hdamp_up", "hdamp_down", "CR_ERD", "CR_QCD_ERD", "CR_gluon_move", "CR_GM_ERD", "JES_up", "JES_down"};
  int nSys = sizeof(sysNames)/sizeof(sysNames[0]);
  
  /// chi2 < 5
  //pair<string,string> input[] = {pair<string,string>("nominal","1341"), pair<string,string>("tune up","1322"), pair<string,string>("tune down","1334"), pair<string,string>("ISR up","1335"), pair<string,string>("ISR down","1337"),pair<string,string>("FSR up","1338"), pair<string,string>("FSR down","1339"), pair<string,string>("hdamp up","1342"), pair<string,string>("hdamp down","1345"), pair<string,string>("CR ERD","1346"), pair<string,string>("CR QCD ERD","1347"), pair<string,string>("CR gluon move","1351"), pair<string,string>("CR gluon move ERD","1353")};
  /// chi2 < 2
  //pair<string,string> input[] = {pair<string,string>("nominal","1646"), pair<string,string>("tune up","1648"), pair<string,string>("tune down","1649"), pair<string,string>("ISR up","1650"), pair<string,string>("ISR down","1651"),pair<string,string>("FSR up","1652"), pair<string,string>("FSR down","1653"), pair<string,string>("hdamp up","1654"), pair<string,string>("hdamp down","1655"), pair<string,string>("CR ERD","1656"), pair<string,string>("CR QCD ERD","1657"), pair<string,string>("CR gluon move","1658"), pair<string,string>("CR gluon move ERD","1659")};
  /// chi2 < 3
  pair<string,string> input[] = {pair<string,string>("nominal","1846"), pair<string,string>("tune up","1902"), pair<string,string>("tune down","1903"), pair<string,string>("ISR up","1917"), pair<string,string>("ISR down","1920"),pair<string,string>("FSR up","1922"), pair<string,string>("FSR down","1924"), pair<string,string>("hdamp up","1926"), pair<string,string>("hdamp down","1928"), pair<string,string>("CR ERD","1930"), pair<string,string>("CR QCD ERD","1932"), pair<string,string>("CR gluon move","1906"), pair<string,string>("CR gluon move ERD","1907"), pair<string,string>("JES up","2207"), pair<string,string>("JES down","2208")};
  /// 1e-4 < chi2 < 3
  //pair<string,string> input[] = {pair<string,string>("nominal","2010"), pair<string,string>("tune up","2013"), pair<string,string>("tune down","2016"), pair<string,string>("ISR up","2021"), pair<string,string>("ISR down","2023"),pair<string,string>("FSR up","2025"), pair<string,string>("FSR down","2028"), pair<string,string>("hdamp up","2030"), pair<string,string>("hdamp down","2032"), pair<string,string>("CR ERD","2034"), pair<string,string>("CR QCD ERD","2036"), pair<string,string>("CR gluon move","2017"), pair<string,string>("CR gluon move ERD","2019")};
  /// 1e-2 < chi2 < 3
  //pair<string,string> input[] = {pair<string,string>("nominal","2109"), pair<string,string>("tune up","2110"), pair<string,string>("tune down","2112"), pair<string,string>("ISR up","2118"), pair<string,string>("ISR down","2120"),pair<string,string>("FSR up","2122"), pair<string,string>("FSR down","2124"), pair<string,string>("hdamp up","2126"), pair<string,string>("hdamp down","2128"), pair<string,string>("CR ERD","2130"), pair<string,string>("CR QCD ERD","2132"), pair<string,string>("CR gluon move","2114"), pair<string,string>("CR gluon move ERD","2116")};
  
  string pathFiles = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/OutputPlots/mu/171107_";
  string pathOutput = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/Systematics/Plots/";
  
  TFile *fileOut = new TFile((pathOutput+"Systematics_normalised_plots.root").c_str(),"RECREATE");
  
  string fileName = "NtuplePlots_nominal.root";
  
  string varNames[] = {"jet4_pT", "Ht_4leadingJets", "top_mass", "top_mass_aKF", "top_mass_diff", "top_mass_cand_diff", "mlb", "ttbar_mass", "dR_b_W", "dR_b_light_max", "dR_b_light_min", "dR_bb", "dR_lep_b_max", "dR_lep_b_min"};
  string axisLabels[] = {"jet 4 p_{T} (GeV)", "H_{T} (GeV)", "m_{t} (GeV)", "m_{t} (GeV)", "m_{t,aKF} - m_{t,bKF} (GeV)", "m_{t} - m_{t'} (GeV)", "m_{lb} (GeV)", "m_{t} + m_{lb} (GeV)", "#Delta R(b,W)", "max #Delta R(b,light)", "min #Delta R(b,light)", "#Delta R(b,b)", "#Delta R(l,b_{h})", "#Delta R(l,b_{l})"};
  int nVars = sizeof(varNames)/sizeof(varNames[0]);
  
  string suffixes[] = {"", "_CM"};
  int nSuffixes = sizeof(suffixes)/sizeof(suffixes[0]);
  
  
  /// Get histograms
  for (int iSys = 0; iSys < nSys; iSys++)
  {
    if ( iSys != 0 && iSys != 5 && iSys != 6 && iSys != 13 && iSys != 14 ) continue;
    
    if ( iSys == 13 ) fileName = "NtuplePlots_JESup.root";
    else if ( iSys == 14 ) fileName = "NtuplePlots_JESdown.root";
    else fileName = "NtuplePlots_nominal.root";
    
    TFile *fileIn = new TFile((pathFiles+input[iSys].second+"/"+fileName).c_str(),"read");
    fileIn->cd();
    
    for (int iVar = 0; iVar < nVars; iVar++)
    {
      histo[varNames[iVar]+"_"+sysNames[iSys]+"_CM"] = (TH1F*) fileIn->Get(("MultiSamplePlot_"+varNames[iVar]+"/"+varNames[iVar]+"_TT_CM_").c_str());
      histo[varNames[iVar]+"_"+sysNames[iSys]+"_CM"]->SetDirectory(0);  // keep histogram when file is closed
      histo[varNames[iVar]+"_"+sysNames[iSys]+"_WM"] = (TH1F*) fileIn->Get(("MultiSamplePlot_"+varNames[iVar]+"/"+varNames[iVar]+"_TT_WM_").c_str());
      histo[varNames[iVar]+"_"+sysNames[iSys]+"_WM"]->SetDirectory(0);  // keep histogram when file is closed
      histo[varNames[iVar]+"_"+sysNames[iSys]+"_UM"] = (TH1F*) fileIn->Get(("MultiSamplePlot_"+varNames[iVar]+"/"+varNames[iVar]+"_TT_UM_").c_str());
      histo[varNames[iVar]+"_"+sysNames[iSys]+"_UM"]->SetDirectory(0);  // keep histogram when file is closed
      histo[varNames[iVar]+"_"+sysNames[iSys]+"_other"] = (TH1F*) fileIn->Get(("MultiSamplePlot_"+varNames[iVar]+"/"+varNames[iVar]+"_TT_other_").c_str());
      histo[varNames[iVar]+"_"+sysNames[iSys]+"_other"]->SetDirectory(0);  // keep histogram when file is closed
      
      histo[varNames[iVar]+"_"+sysNames[iSys]] = (TH1F*) histo[varNames[iVar]+"_"+sysNames[iSys]+"_CM"]->Clone((varNames[iVar]+"_"+sysNames[iSys]).c_str());
      histo[varNames[iVar]+"_"+sysNames[iSys]]->Add(histo[varNames[iVar]+"_"+sysNames[iSys]+"_WM"]);
      histo[varNames[iVar]+"_"+sysNames[iSys]]->Add(histo[varNames[iVar]+"_"+sysNames[iSys]+"_UM"]);
      histo[varNames[iVar]+"_"+sysNames[iSys]]->Add(histo[varNames[iVar]+"_"+sysNames[iSys]+"_other"]);
      histo[varNames[iVar]+"_"+sysNames[iSys]]->SetDirectory(0);
    }
    
    fileIn->Close();
  }
  
  
  fileOut->cd();
  
  for (int iVar = 0; iVar < nVars; iVar++)
  {
    
    for (int iSuffix = 0; iSuffix < nSuffixes; iSuffix++)
    {
      TCanvas *c1 = new TCanvas((varNames[iVar]+suffixes[iSuffix]+"_overlay_norm").c_str(),(varNames[iVar]+suffixes[iSuffix]+"_overlay_norm").c_str());
      c1->cd();

      /// Make legend
      TLegend *leg = new TLegend(0.65,0.55,0.89,0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);

      for (int iSys = 0; iSys < nSys; iSys++)
      {
        if ( iSys != 0 && iSys != 5 && iSys != 6 && iSys != 13 && iSys != 14 ) continue;
        
        histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->SetLineWidth(2);
        histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->SetMarkerStyle(listMarkers[iSys%8]);

        if ( iSys == 0 )
        {
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->SetLineColor(kBlack);
        histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->SetMarkerColor(kBlack);
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->SetTitle("");
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->GetXaxis()->SetTitle((axisLabels[iVar]).c_str());
          //histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->GetXaxis()->SetRangeUser(-60.,50.);
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->GetYaxis()->SetTitle("# normalised events");
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->GetYaxis()->SetTitleOffset(1.3);
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->SetStats(0);
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->DrawNormalized("P");
        }
        else
        {
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->SetLineColor(listRainbow[iSys%6]);
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->SetMarkerColor(listRainbow[iSys%6]);
          histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]]->DrawNormalized("P same");
        }
        leg->AddEntry(histo[varNames[iVar]+"_"+sysNames[iSys]+suffixes[iSuffix]],(sysNames[iSys]).c_str(),"lp");
        c1->Update();
      }
      leg->Draw();

      c1->Update();
      c1->Write();
      c1->SaveAs((pathOutput+"Syst_overlay_"+varNames[iVar]+suffixes[iSuffix]+"_normalised.png").c_str());
      c1->SetName((varNames[iVar]+suffixes[iSuffix]+"_overlay_norm_logY").c_str());
      c1->SetLogy(1);
      c1->Write();
      c1->SaveAs((pathOutput+"Syst_overlay_"+varNames[iVar]+suffixes[iSuffix]+"_normalised_logY.png").c_str());
      c1->SetLogy(0);
    
    }  // end suffixes
    
  }  // end vars
  
  
  fileOut->Close();
  
  return 0;
}
