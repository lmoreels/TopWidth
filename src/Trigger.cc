#include "../interface/Trigger.h"

Trigger::Trigger(bool isMuon, bool isElectron):
muon(false), electron(false), trigged(false), redotrigmap(false), triggerList(), currentRun(0), previousRun(-1), currentFilename(""),
previousFilename(""), iFile(-1), triggermap()
{
  if (isMuon)
  {
    muon = true;
  }
  else if (isElectron)
  {
    electron = true;
  }
  else
  {
    cout << "TRIGGER::TRIGGER - No selected lepton..." << endl;
  }
}

Trigger::~Trigger()
{
  
}

void Trigger::bookTriggers(bool isData)
{
  triggerList.clear();
  if (muon)
  {
    if (isData)
    {
      triggerList.push_back("HLT_IsoMu18_v*");  // Data
    }
    else
    {
      triggerList.push_back("HLT_IsoMu17_eta2p1_v*");  // MC
    }
  }
  
  if (electron)
  {
    if (isData)
     {
       triggerList.push_back("HLT_Ele23_WPLoose_Gsf_v*");  // Data, restricted to eta < 2.1
     }
     else
     {
       triggerList.push_back("HLT_Ele22_eta2p1_WP75_Gsf_v*");  // MC
     }   	
   }
  
  for(UInt_t iTrig = 0; iTrig < triggerList.size(); iTrig++)
  {
    triggermap[triggerList[iTrig]] = std::pair<int,bool> (-999,false);
  }
  
  // for(std::map<std::string,std::pair<int,bool> >::iterator trigiter = triggermap.begin(); trigiter != triggermap.end(); trigiter++){
  //     std::pair<int,bool> bla = trigiter->second;
  //     std::string bla2 = trigiter->first; 
  // }    
}

void Trigger::checkAvail(int currentRun, vector < Dataset* > datasets, unsigned int d, TTreeLoader *treeLoader, TRootEvent* event)
{
  redotrigmap = false;
  currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
  if (previousFilename != currentFilename)
  {
    previousFilename = currentFilename;
    iFile++;
    redotrigmap = true;
    cout << "File changed!!! => iFile = " << iFile << " new file is " << datasets[d]->eventTree()->GetFile()->GetName() << " in sample " << datasets[d]->Name() << endl;
  }
  if (previousRun != currentRun)
  {
    previousRun = currentRun;
    //cout << "*****!!!!new run!!!! => new run = " << previousRun << " *****" << endl;
    redotrigmap=true;
  }
  
  if (redotrigmap)
  {
    //treeLoader->ListTriggers(currentRun, iFile);
  }
  
  
  // get trigger info:
  for(std::map<std::string,std::pair<int,bool> >::iterator iter = triggermap.begin(); iter != triggermap.end(); iter++)
  {
    if (redotrigmap)
    {
      Int_t loc = treeLoader->iTrigger(iter->first, currentRun, iFile);
      string trigname = iter->first;
      //cout << "trigname: " << trigname << "  location: " << loc << endl;
      
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