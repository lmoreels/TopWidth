#include "../interface/SelectionTables.h"

SelectionTables::SelectionTables(vector<Dataset*> datasets):
precision_(0), listOfDatasets_(datasets), lumi_(1.)
{
  numberOfDatasets_ = listOfDatasets_.size();
  numberOfDatasetsMerged_ = numberOfDatasets_;
  
  eqLumi_ = new double[numberOfDatasets_];
  for (int d = 0; d < numberOfDatasets_; d++)
  {
    eqLumi_[d] = 1.;
  }
}

SelectionTables::~SelectionTables()
{
  
}

void SelectionTables::SetPrecision(int i)
{
  precision_ = i;
}

void SelectionTables::SetLumi(double lumi)
{
  lumi_ = lumi;
}

void SelectionTables::SetEqLumi(int d, double eqLumi)
{
  eqLumi_[d] = eqLumi;
}

void SelectionTables::AddCutStep(string cutStepName)
{
  listOfCuts_.push_back(cutStepName);
}

void SelectionTables::SetUpTable()
{
  numberOfCuts_ = listOfCuts_.size();
  
  nofEventsRaw_      = new double*[numberOfCuts_];
  nofEventsRawError_ = new double*[numberOfCuts_];
  nofEventsRawVar_   = new double*[numberOfCuts_];
  nofEvents_         = new double*[numberOfCuts_];
  nofEventsError_    = new double*[numberOfCuts_];
  nofEventsVar_      = new double*[numberOfCuts_];
  totalEvents_       = new double[numberOfCuts_];
  totalEventsError_  = new double[numberOfCuts_];
  totalEventsVar_    = new double[numberOfCuts_];
  for (int iCut = 0; iCut < numberOfCuts_; iCut++)
  {
    nofEventsRaw_[iCut]      = new double[numberOfDatasets_];
    nofEventsRawError_[iCut] = new double[numberOfDatasets_];
    nofEventsRawVar_[iCut]   = new double[numberOfDatasets_];
    nofEvents_[iCut]         = new double[numberOfDatasets_];
    nofEventsError_[iCut]    = new double[numberOfDatasets_];
    nofEventsVar_[iCut]      = new double[numberOfDatasets_];
    totalEvents_[iCut]       = 0.;
    totalEventsError_[iCut]  = 0.;
    totalEventsVar_[iCut]    = 0.;
    for (int d = 0; d < numberOfDatasets_; d++)
    {
      nofEventsRaw_[iCut][d]      = 0.;
      nofEventsRawError_[iCut][d] = 0.;
      nofEventsRawVar_[iCut][d]   = 0.;
      nofEvents_[iCut][d]         = 0.;
      nofEventsError_[iCut][d]    = 0.;
      nofEventsVar_[iCut][d]      = 0.;
    }
  }
}

void SelectionTables::Fill(int d, int cutStep, double value)
{
  if ( cutStep >= numberOfCuts_ )
  {
    std::cerr << "SelectionTables::Fill: Cut step number beyond scope..." << std::endl;
  }
  else if ( d >= numberOfDatasets_ )
  {
    std::cerr << "SelectionTables::Fill: Dataset number beyond scope..." << std::endl;
  }
  else
  {
    nofEventsRaw_[cutStep][d] = value;
  }
}

void SelectionTables::Fill(int d, vector<double> values)
{
  if ( values.size() > (unsigned int)numberOfCuts_ )
  {
    std::cerr << "SelectionTables::Fill: Too many values for initialised number of cuts. Only filling table up to cut number " << numberOfCuts_ << std::endl;
  }
  if ( d > numberOfDatasets_ )
  {
    std::cerr << "SelectionTables::Fill: Dataset number beyond scope..." << std::endl;
  }
  else
  {
    for (int iCut = 0; iCut < numberOfCuts_; iCut++)
    {
      nofEventsRaw_[iCut][d] = values[iCut];
    }
  }
}

void SelectionTables::MergeDatasets()
{
  //dataSetTitle_ = datasets[d]->Title();
  string dataSetName, dataSetNameMerged;
  bool addedST = false, addedDY = false, addedW = false;
  for (int d = 0; d < numberOfDatasets_; d++)
  {
    dataSetName = listOfDatasets_[d]->Name();
    
    if ( dataSetName.find("ST") != std::string::npos )
    {
      if (! addedST)
      {
        listOfDatasetsMerged_.push_back(new Dataset("ST_Merged", "\\text{single top}", true, listOfDatasets_[d]->Color(), listOfDatasets_[d]->LineStyle(), listOfDatasets_[d]->LineWidth(), 1., 1.));
        addedST = true;
      }
    }
    else if ( dataSetName.find("Jets") != std::string::npos )
    {
      if( dataSetName.find("DY") != std::string::npos )
      {
        if (! addedDY)
        {
          listOfDatasetsMerged_.push_back(new Dataset("DYJets_Merged", "\\text{DY${}+{}$jets}", true, listOfDatasets_[d]->Color(), listOfDatasets_[d]->LineStyle(), listOfDatasets_[d]->LineWidth(), 1., 1.));
          addedDY = true;
        }
      }
      else if ( dataSetName.find("W") != std::string::npos )
      {
        if (! addedW)
        {
          listOfDatasetsMerged_.push_back(new Dataset("WJets_Merged", "\\text{W${}+{}$jets}", true, listOfDatasets_[d]->Color(), listOfDatasets_[d]->LineStyle(), listOfDatasets_[d]->LineWidth(), 1., 1.));
          addedW = true;
        }
      }
    }
    else
    {
      if ( dataSetName.find("TT") != std::string::npos ) listOfDatasets_[d]->SetTitle("\\ttbar");
      listOfDatasetsMerged_.push_back(listOfDatasets_[d]);
    }
  }
  
  numberOfDatasetsMerged_ = listOfDatasetsMerged_.size();
  
  nofEventsMerged_      = new double*[numberOfCuts_];
  nofEventsErrorMerged_ = new double*[numberOfCuts_];
  nofEventsVarMerged_   = new double*[numberOfCuts_];
  for (int iCut = 0; iCut < numberOfCuts_; iCut++)
  {
    nofEventsMerged_[iCut]      = new double[numberOfDatasetsMerged_];
    nofEventsErrorMerged_[iCut] = new double[numberOfDatasetsMerged_];
    nofEventsVarMerged_[iCut]   = new double[numberOfDatasetsMerged_];
    for (int dM = 0; dM < numberOfDatasetsMerged_; dM++)
    {
      nofEventsMerged_[iCut][dM]      = 0.;
      nofEventsErrorMerged_[iCut][dM] = 0.;
      nofEventsVarMerged_[iCut][dM]   = 0.;
      
      dataSetNameMerged = listOfDatasetsMerged_[dM]->Name();
      for (int d = 0; d < numberOfDatasets_; d++)
      {
        dataSetName = listOfDatasets_[d]->Name();
        if ( (dataSetNameMerged.find("TT") != std::string::npos && dataSetName.find("TT") != std::string::npos) || 
            (dataSetNameMerged.find("ST") != std::string::npos && dataSetName.find("ST") != std::string::npos) || 
            (dataSetNameMerged.find("WJets") != std::string::npos && dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos) || 
            (dataSetNameMerged.find("DYJets") != std::string::npos && dataSetName.find("DY") != std::string::npos && dataSetName.find("Jets") != std::string::npos) ||
            (dataSetNameMerged.find("data") != std::string::npos && dataSetName.find("data") != std::string::npos) )
        {
          //nofEventsRawMerged_[iCut][dM]      += nofEventsRaw_[iCut][d];   // not allowed unless all constituents have the same eqLumi
          nofEventsMerged_[iCut][dM]    += nofEvents_[iCut][d];
          nofEventsVarMerged_[iCut][dM] += nofEventsVar_[iCut][d];
        }
        
      }  // end d
      
      nofEventsErrorMerged_[iCut][dM] = sqrt(nofEventsVarMerged_[iCut][dM]);
      
    }  // end dM
  }  // end cut
  
}

void SelectionTables::CalculateTable()
{
  string dataSetName;
  double lumiWeight;
  
  for (int d = 0; d < numberOfDatasets_; d++)
  {
    dataSetName = listOfDatasets_[d]->Name();
    if ( dataSetName.find("data") != std::string::npos ) lumiWeight = 1.;
    else lumiWeight = lumi_/eqLumi_[d];
    if ( eqLumi_[d] == 1. ) std::cout << "SelectionTables::CalculateTable: WARNING: eqLumi = 1 for dataset " << dataSetName << ". Please check if this is correct..." << endl;
      
    for (int iCut = 0; iCut < numberOfCuts_; iCut++)
    {
      if ( nofEventsRaw_[0][d] == 0. ) nofEventsRawError_[iCut][d] = 0.;
      else nofEventsRawError_[iCut][d] = ErrorCalculator(nofEventsRaw_[iCut][d], nofEventsRaw_[iCut][d]/nofEventsRaw_[0][d], 1.);
      nofEventsRawVar_[iCut][d]   = pow(nofEventsRawError_[iCut][d], 2.);
      
      nofEvents_[iCut][d]      = nofEventsRaw_[iCut][d] * lumiWeight;
      if ( nofEvents_[0][d] == 0. ) nofEventsError_[iCut][d] = 0.;
      else nofEventsError_[iCut][d] = ErrorCalculator(nofEvents_[iCut][d], nofEvents_[iCut][d]/nofEvents_[0][d], 1.);
      nofEventsVar_[iCut][d]   = pow(nofEventsError_[iCut][d], 2.);
      
      totalEvents_[iCut]    += nofEvents_[iCut][d];
      totalEventsVar_[iCut] += nofEventsVar_[iCut][d];
    }  // end cut
  }  // end d
  
  for (int iCut = 0; iCut < numberOfCuts_; iCut++)
  {
    totalEventsError_[iCut] = sqrt(totalEventsVar_[iCut]);
  }
  
  this->MergeDatasets();
}

void SelectionTables::WriteTable(ofstream& fout, double** listTable_,double** listTableError_, bool writeError, bool writeMerged, bool writeLandscape)
{
  int nDatasets = ( writeMerged ? numberOfDatasetsMerged_ : numberOfDatasets_ );
  if( precision_ >= 0 ) fout << fixed << setprecision(precision_);
  
  string interline = "";
  for (int d = 0; d < nDatasets; d++) interline += "&";
  interline += "&\\\\[-6pt]";
  string headerextra = "[+3pt]";
  
  if(writeLandscape) fout << "\\begin{landscape}" << std::endl;
  fout << "\\begin{table}" << std::endl;
  fout << "\\caption{Effect of each selection requirement on the number of selected events. ($" << lumi_/1000. << "\\fbinv$ of int. lumi.)}" << std::endl;
  fout << "\\label{tab:selectionTable}" << std::endl;
  fout << "\\begin{center}" << std::endl;
  fout << "\\begin{tabular}{r|";
  for (int d = 1; d < nDatasets; d++) fout << "c";
  fout << "|c|c}" << std::endl;
  fout << "\\hline" << std::endl;
  fout << "\\hline" << std::endl;
  fout << interline << std::endl;
  fout << " & ";
  for (int d = 1; d < nDatasets; d++)
  {
    if (writeMerged) fout << "$" << listOfDatasetsMerged_[d]->Title() << "$";
    else fout << "$" << listOfDatasets_[d]->Title() << "$";
    
    if ( d < nDatasets-1 ) fout << " & ";
    else fout << " & Total exp. & Observed \\\\" << headerextra << std::endl;
  }
  fout << "\\hline" << std::endl;
  fout << interline << std::endl;
  for (int iCut = 0; iCut < numberOfCuts_; iCut++)
  {
    fout << listOfCuts_[iCut] << " & ";
    for (int d = 1; d < nDatasets; d++)
    {
      fout << listTable_[iCut][d];
      if (writeError) fout << " $\\pm$ " << listTableError_[iCut][d];
      
      if ( d < nDatasets-1 ) fout << " & ";
      else
      {
        fout << " & " << totalEvents_[iCut];
        if (writeError) fout << " $\\pm$ " << totalEventsError_[iCut];
        fout << " & " << listTable_[iCut][0];
        if (writeError) fout << " $\\pm$ " << listTableError_[iCut][0];
        fout << " \\\\" << std::endl;
      }
    }
  }
  fout << interline << std::endl;
  fout << "\\hline" << std::endl;
  fout << "\\hline" << std::endl;
  fout << "\\end{tabular}" << std::endl;
  fout << "\\end{center}" << std::endl;
  fout << "\\end{table}" << std::endl;
  if (writeLandscape) fout << "\\end{landscape}%" << std::endl;
}

void SelectionTables::Write(string filename, bool writeError, bool writeMerged, bool writeLandscape)
{
  ofstream fout(filename.c_str());
  this->Write(fout, writeError, writeMerged, writeLandscape);
}

void SelectionTables::Write(ofstream& fout, bool writeError, bool writeMerged, bool writeLandscape)
{
  if (writeMerged) this->WriteTable(fout, nofEventsMerged_, nofEventsErrorMerged_, writeError, writeMerged, writeLandscape);
  else this->WriteTable(fout, nofEvents_, nofEventsError_, writeError, writeMerged, writeLandscape);
}
