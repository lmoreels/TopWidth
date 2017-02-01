#include "../interface/Trigger.h"

Trigger::Trigger(bool isMuon, bool isElectron, bool trigSingleLep, bool trigDoubleLep):
muon(false), electron(false), singleLep(false), doubleLep(false), fullHadr(false), trigged(false), redotrigmap(false), triggerList(), currentRunTrig(0), previousRunTrig(-1), currentFilenameTrig(""), previousFilenameTrig(""), iFileTrig(-1), treeNumberTrig(-1), triggermap()
{
  if (isMuon)
  {
    muon = true;
  }
  if (isElectron)
  {
    electron = true;
  }
  
  if (trigSingleLep)
  {
    singleLep = true;
  }
  if (trigDoubleLep)
  {
    doubleLep = true;
  }
  if (! trigSingleLep && ! trigDoubleLep)
  {
    fullHadr = true;
    cout << "TRIGGER::TRIGGER - Do you really want no lepton triggers?" << endl;
  }
  else if ( trigSingleLep && trigDoubleLep )
  {
    cout << "TRIGGER::TRIGGER - Doing logical OR between double and single lepton triggers..." << endl;
  }
}

Trigger::~Trigger()
{
  
}

void Trigger::bookTriggers(bool isData)
{
  /// This function is called in the dataset loop
  //  Reset all quantities here
  triggerList.clear();
  currentRunTrig = 0;
  previousRunTrig = -1;
  currentFilenameTrig = previousFilenameTrig = "";
  iFileTrig = -1;
  treeNumberTrig = -1;
  triggermap.clear();
  
  /// Add relevant triggers to triggerlist
  //  Recommended triggers for TOP analyses
  //  Last updated: 30 March 2016. https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger
  if (singleLep)
  {
    if (muon)   // Triggers for data & MC are the same
    {
//      triggerList.push_back("HLT_IsoMu18_v*");  // Not available for Run C: 254231 - 254907
      triggerList.push_back("HLT_IsoMu20_v*");
      triggerList.push_back("HLT_IsoTkMu20_v*");
    }
    if (electron)
    {
      triggerList.push_back("HLT_Ele23_WPLoose_Gsf_v*"); 	
    }
  }
  
  if (doubleLep)   // Triggers for data & MC are the same
  {
    if (muon && ! electron)
    {
      triggerList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
    }
    if (! muon && electron)
    {
      triggerList.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    }
    if (muon && electron)
    {
      triggerList.push_back("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
      triggerList.push_back("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*");
    }
  }
  
  if (fullHadr)
  {
    triggerList.push_back("HLT_PFHT450_SixJet40_PFBTagCSV0p72_v*");
    triggerList.push_back("HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v*");
  }
  
  for(UInt_t iTrig = 0; iTrig < triggerList.size(); iTrig++)
  {
    triggermap[triggerList[iTrig]] = std::pair<int,bool> (-999,false);
  }
  
}

void Trigger::checkAvail(int currentRunTrig, vector < Dataset* > datasets, unsigned int d, TTreeLoader *treeLoader, TRootEvent* event, bool verbose)
{
  redotrigmap = false;
  treeNumberTrig = datasets[d]->eventTree()->GetTreeNumber();
  currentFilenameTrig = datasets[d]->eventTree()->GetFile()->GetName();
  if (previousFilenameTrig != currentFilenameTrig)
  {
    previousFilenameTrig = currentFilenameTrig;
    iFileTrig++;
    redotrigmap = true;
    cout << endl << "*****File changed!!! => iFile = " << iFileTrig << " new file is " << currentFilenameTrig << " in sample " << datasets[d]->Name() << " *****" << endl;
  }
  if (previousRunTrig != currentRunTrig)
  {
    previousRunTrig = currentRunTrig;
    cout << "*****Run changed!!! => new run = " << previousRunTrig << " *****" << endl;
    redotrigmap=true;
  }
  
  if (verbose && redotrigmap)
  {
    treeLoader->ListTriggers(currentRunTrig, treeNumberTrig);
  }
  
  
  // get trigger info:
  for(std::map<std::string,std::pair<int,bool> >::iterator iter = triggermap.begin(); iter != triggermap.end(); iter++)
  {
    if (redotrigmap)
    {
      Int_t loc = treeLoader->iTrigger(iter->first, currentRunTrig, treeNumberTrig);
      string trigname = iter->first;
      cout << "trigname: " << trigname << "  location: " << loc << /*"  " << currentFilenameTrig << "  " << currentRunTrig <<*/ endl;
      
      iter->second.first = loc;
    }
    // and check if it exists and if it fired:
    if (iter->second.first >= 0 && iter->second.first != 9999) // trigger exists
      iter->second.second = event->trigHLT(iter->second.first);
    else
      iter->second.second = false;
  }   
}


int Trigger::checkIfFired()
{
  // now check if the appropriate triggers fired for each analysis:
  trigged = 0;
  
  for(UInt_t itrig = 0; itrig < triggerList.size() && trigged == 0; itrig++)
  {
    if (triggermap[triggerList[itrig]].second)
      trigged = 1;
  }
  
  return trigged;
}