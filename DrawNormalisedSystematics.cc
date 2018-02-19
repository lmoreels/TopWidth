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
  Color_t listRainbow[] = {kGreen+1, kRed};
  //Color_t listRainbow[] = {kRed, kOrange-3, kYellow-7, /*kGreen-7,*/ kGreen+1, kCyan+1,/* kBlue+2,*/ kMagenta/*, kViolet-5, kPink+10*/};
  Style_t listMarkers[] = {24, 26, 32, 25, 30, 27, 28, 46, 42};
  
  string suffixName = "";
  
  string sysNames[] = {"nominal", "tune_up", "tune_down", "ISR_up", "ISR_down", "FSR_up", "FSR_down", "hdamp_up", "hdamp_down", "CR_ERD", "CR_QCD_ERD", "CR_gluon_move", "CR_GM_ERD", "herwig", "JES_up", "JES_down", "JER_up", "JER_down"};
  int drawSys[] = {0, 5, 6}; suffixName = "_fsr";
  //int drawSys[] = {0, 3, 4, 5, 6, 7, 8}; suffixName = "_shower";
  //int drawSys[] = {0, 1, 2, 9, 10, 11, 12}; suffixName = "_CR";
  //int drawSys[] = {0, 13}; suffixName = "_herwig";
  int nSys = sizeof(drawSys)/sizeof(drawSys[0]);
  //int nSys = sizeof(sysNames)/sizeof(sysNames[0])
  
  /// chi2 < 15
  //pair<string,string> input[] = {pair<string,string>("nominal", "171116_1324"), pair<string,string>("tune up", "171113_1712"), pair<string,string>("tune down", "171113_1714"), pair<string,string>("ISR up", "171113_1722"), pair<string,string>("ISR down", "171113_1724"),pair<string,string>("FSR up", "171113_1726"), pair<string,string>("FSR down", "171113_1728"), pair<string,string>("hdamp up", "171113_1730"), pair<string,string>("hdamp down", "171113_1732"), pair<string,string>("CR ERD", "171113_1734"), pair<string,string>("CR QCD ERD", "171113_1736"), pair<string,string>("CR gluon move", "171113_1717"), pair<string,string>("CR gluon move ERD", "171113_1718"), pair<string,string>("herwig", "171113_1720")};
  // mlb < 220 GeV
  pair<string,string> input[] = {pair<string,string>("nominal", /*"171204_1701"*/"180121_2012"), pair<string,string>("tune up", "171204_1800"), pair<string,string>("tune down", ""), pair<string,string>("ISR up", "171204_1805"), pair<string,string>("ISR down", "171204_1807"),pair<string,string>("FSR up", /*"171204_1808"*/"180121_2141"), pair<string,string>("FSR down", /*"171205_1009"*/"180121_2143"), pair<string,string>("hdamp up", "171204_1813"), pair<string,string>("hdamp down", "171204_1815"), pair<string,string>("CR ERD", "171204_1817"), pair<string,string>("CR QCD ERD", "171204_1819"), pair<string,string>("CR gluon move", "171204_1803"), pair<string,string>("CR gluon move ERD", "171204_1804"), pair<string,string>("herwig", ""), pair<string,string>("JES up", "171204_1827"), pair<string,string>("JES down", "171204_1825"), pair<string,string>("JER up", "171204_1823"), pair<string,string>("JES down", "171205_1011")};
  
  string date = "171107_";
  /// chi2 < 5
  //pair<string,string> input[] = {pair<string,string>("nominal","1341"), pair<string,string>("tune up","1322"), pair<string,string>("tune down","1334"), pair<string,string>("ISR up","1335"), pair<string,string>("ISR down","1337"),pair<string,string>("FSR up","1338"), pair<string,string>("FSR down","1339"), pair<string,string>("hdamp up","1342"), pair<string,string>("hdamp down","1345"), pair<string,string>("CR ERD","1346"), pair<string,string>("CR QCD ERD","1347"), pair<string,string>("CR gluon move","1351"), pair<string,string>("CR gluon move ERD","1353")};
  /// chi2 < 2
  //pair<string,string> input[] = {pair<string,string>("nominal","1646"), pair<string,string>("tune up","1648"), pair<string,string>("tune down","1649"), pair<string,string>("ISR up","1650"), pair<string,string>("ISR down","1651"),pair<string,string>("FSR up","1652"), pair<string,string>("FSR down","1653"), pair<string,string>("hdamp up","1654"), pair<string,string>("hdamp down","1655"), pair<string,string>("CR ERD","1656"), pair<string,string>("CR QCD ERD","1657"), pair<string,string>("CR gluon move","1658"), pair<string,string>("CR gluon move ERD","1659")};
  /// chi2 < 3
  //pair<string,string> input[] = {pair<string,string>("nominal","1846"), pair<string,string>("tune up","1902"), pair<string,string>("tune down","1903"), pair<string,string>("ISR up","1917"), pair<string,string>("ISR down","1920"),pair<string,string>("FSR up","1922"), pair<string,string>("FSR down","1924"), pair<string,string>("hdamp up","1926"), pair<string,string>("hdamp down","1928"), pair<string,string>("CR ERD","1930"), pair<string,string>("CR QCD ERD","1932"), pair<string,string>("CR gluon move","1906"), pair<string,string>("CR gluon move ERD","1907"), pair<string,string>("JES up","2207"), pair<string,string>("JES down","2208")};
  /// 1e-4 < chi2 < 3
  //pair<string,string> input[] = {pair<string,string>("nominal","2010"), pair<string,string>("tune up","2013"), pair<string,string>("tune down","2016"), pair<string,string>("ISR up","2021"), pair<string,string>("ISR down","2023"),pair<string,string>("FSR up","2025"), pair<string,string>("FSR down","2028"), pair<string,string>("hdamp up","2030"), pair<string,string>("hdamp down","2032"), pair<string,string>("CR ERD","2034"), pair<string,string>("CR QCD ERD","2036"), pair<string,string>("CR gluon move","2017"), pair<string,string>("CR gluon move ERD","2019")};
  /// 1e-2 < chi2 < 3
  //pair<string,string> input[] = {pair<string,string>("nominal","2109"), pair<string,string>("tune up","2110"), pair<string,string>("tune down","2112"), pair<string,string>("ISR up","2118"), pair<string,string>("ISR down","2120"),pair<string,string>("FSR up","2122"), pair<string,string>("FSR down","2124"), pair<string,string>("hdamp up","2126"), pair<string,string>("hdamp down","2128"), pair<string,string>("CR ERD","2130"), pair<string,string>("CR QCD ERD","2132"), pair<string,string>("CR gluon move","2114"), pair<string,string>("CR gluon move ERD","2116")};
  
  //string pathFiles = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/OutputPlots/mu/"+date;
  //string pathFiles = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidthArchive/OutputPlots/mu/";
  string pathFiles = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/OutputPlots/mu/";
  string pathOutput = "/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/Systematics/Plots/";
  
  TFile *fileOut = new TFile((pathOutput+"Systematics_normalised_plots"+suffixName+"_new.root").c_str(),"RECREATE");
  
  string fileName = "NtuplePlots_nominal.root";
  
  string varNames[] = {"jet3_pT", "jet3_pT_aKF", "jet4_pT", "jet4_pT_aKF", "Ht_4leadingJets", "Ht_4leadingJets_aKF", "top_pT_aKF", "top_mass", "top_mass_aKF", "top_mass_alt", "top_mass_alt_aKF", "top_mass_diff", "top_mass_cand_diff", "top_mass_cand_diff_aKF", "mlb", "mlb_aKF", "ttbar_mass", "ttbar_mass_aKF", "W_mass_diff", "W_mass_T_diff", "dR_b_W_aKF", "dR_b_light_max_aKF", "dR_b_light_min_aKF", "dR_bb_aKF", "dR_lep_b_max_aKF", "dR_lep_b_min_aKF"};
  string axisLabels[] = {"jet 3 p_{T} (GeV)", "jet 3 p_{T} (GeV)", "jet 4 p_{T} (GeV)", "jet 4 p_{T} (GeV)", "H_{T} (GeV)", "H_{T} (GeV)", "top p_{T} (GeV)", "m_{bjj} (GeV)", "m_{bjj} (GeV)", "m_{b'jj} (GeV)", "m_{b'jj} (GeV)", "m_{bjj,aKF} - m_{bjj,bKF} (GeV)", "m_{bjj} - m_{b'jj} (GeV)", "m_{bjj} - m_{b'jj} (GeV)", "m_{lb} (GeV)", "m_{lb} (GeV)", "m_{bjj} + m_{lb} (GeV)", "m_{bjj} + m_{lb} (GeV)", "m_{jj,aKF} - m_{jj,bKF} (GeV)", "m_{T,W,aKF} - m_{T,W,bKF} (GeV)", "#Delta R(b,W)", "max #Delta R(b,light)", "min #Delta R(b,light)", "#Delta R(b,b)", "#Delta R(l,b_{h})", "#Delta R(l,b_{l})"};
  int nVars = sizeof(varNames)/sizeof(varNames[0]);
  
  string suffixes[] = {"", "_CM"};
  int nSuffixes = sizeof(suffixes)/sizeof(suffixes[0]);
  
  
  /// Get histograms
  for (int iSys = 0; iSys < nSys; iSys++)
  {
    if ( sysNames[drawSys[iSys]].find("JES_up") != std::string::npos ) fileName = "NtuplePlots_JESup.root";
    else if ( sysNames[drawSys[iSys]].find("JES_down") != std::string::npos ) fileName = "NtuplePlots_JESdown.root";
    else if ( sysNames[drawSys[iSys]].find("JER_up") != std::string::npos ) fileName = "NtuplePlots_JERup.root";
    else if ( sysNames[drawSys[iSys]].find("JER_down") != std::string::npos ) fileName = "NtuplePlots_JERdown.root";
    else fileName = "NtuplePlots_nominal.root";
    
    TFile *fileIn = new TFile((pathFiles+input[drawSys[iSys]].second+"/"+fileName).c_str(),"read");
    fileIn->cd();
    
    for (int iVar = 0; iVar < nVars; iVar++)
    {
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_CM"] = (TH1F*) fileIn->Get(("MultiSamplePlot_"+varNames[iVar]+"/"+varNames[iVar]+"_TT_CM_").c_str());
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_CM"]->SetDirectory(0);  // keep histogram when file is closed
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_WM"] = (TH1F*) fileIn->Get(("MultiSamplePlot_"+varNames[iVar]+"/"+varNames[iVar]+"_TT_WM_").c_str());
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_WM"]->SetDirectory(0);  // keep histogram when file is closed
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_UM"] = (TH1F*) fileIn->Get(("MultiSamplePlot_"+varNames[iVar]+"/"+varNames[iVar]+"_TT_UM_").c_str());
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_UM"]->SetDirectory(0);  // keep histogram when file is closed
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_other"] = (TH1F*) fileIn->Get(("MultiSamplePlot_"+varNames[iVar]+"/"+varNames[iVar]+"_TT_other_").c_str());
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_other"]->SetDirectory(0);  // keep histogram when file is closed
      
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]] = (TH1F*) histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_CM"]->Clone((varNames[iVar]+"_"+sysNames[drawSys[iSys]]).c_str());
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]]->Add(histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_WM"]);
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]]->Add(histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_UM"]);
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]]->Add(histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+"_other"]);
      histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]]->SetDirectory(0);
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
      c1->SetTickx(1);
      c1->SetTicky(1);

      /// Make legend
      TLegend *leg = new TLegend(0.65,0.55,0.89,0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);

      for (int iSys = 0; iSys < nSys; iSys++)
      {
        histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->SetLineWidth(2);
        histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->SetMarkerStyle(listMarkers[iSys%8]);

        if ( iSys == 0 )
        {
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->SetLineColor(kBlack);
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->SetMarkerColor(kBlack);
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->SetTitle("");
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->GetXaxis()->SetTitle((axisLabels[iVar]).c_str());
          //histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->GetXaxis()->SetRangeUser(-60.,50.);
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->GetYaxis()->SetTitle("# normalised events");
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->GetYaxis()->SetTitleOffset(1.3);
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->SetStats(0);
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->DrawNormalized("P");
        }
        else
        {
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->SetLineColor(listRainbow[(iSys-1)%6]);
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->SetMarkerColor(listRainbow[(iSys-1)%6]);
          histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]]->DrawNormalized("P same");
        }
        leg->AddEntry(histo[varNames[iVar]+"_"+sysNames[drawSys[iSys]]+suffixes[iSuffix]],(sysNames[drawSys[iSys]]).c_str(),"lp");
        c1->Update();
      }
      leg->Draw();
      
      histo[varNames[iVar]+"_"+sysNames[drawSys[0]]+suffixes[iSuffix]]->DrawNormalized("P same");  // Draw nominal on top
      c1->Update();
      c1->Write();
      c1->SaveAs((pathOutput+"Syst_overlay_"+varNames[iVar]+suffixes[iSuffix]+"_normalised"+suffixName+"_new.png").c_str());
      c1->SetName((varNames[iVar]+suffixes[iSuffix]+"_overlay_norm_logY").c_str());
      c1->SetLogy(1);
      c1->Write();
      c1->SaveAs((pathOutput+"Syst_overlay_"+varNames[iVar]+suffixes[iSuffix]+"_normalised"+suffixName+"_logY_new.png").c_str());
      c1->SetLogy(0);
    
    }  // end suffixes
    
  }  // end vars
  
  
  fileOut->Close();
  
  return 0;
}
