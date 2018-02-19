#include "../interface/SelectionTables.h"

SelectionTables::SelectionTables(std::vector<Dataset*> datasets):
precision_(0), listOfDatasets_(datasets), lumi_(1.), extremeMerge_(false)
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
  slumi_.precision(3); slumi_ << lumi_/1000.;
}

void SelectionTables::SetEqLumi(int d, double eqLumi)
{
  eqLumi_[d] = eqLumi;
}

void SelectionTables::SetExtremeMerge(bool merge)
{
  extremeMerge_ = merge;
}

void SelectionTables::AddCutStep(std::string cutStepName)
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

void SelectionTables::Fill(int d, std::vector<double> values)
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
  std::string dataSetName, dataSetNameMerged;
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
      if (extremeMerge_)
      {
        if (! addedDY && ! addedW)
        {
          listOfDatasetsMerged_.push_back(new Dataset("Other_Merged", "\\text{other}", true, listOfDatasets_[d]->Color(), listOfDatasets_[d]->LineStyle(), listOfDatasets_[d]->LineWidth(), 1., 1.));
          addedDY = true;
          addedW = true;
        }
      }
      else if ( dataSetName.find("DY") != std::string::npos )
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
        if ( /*(dataSetNameMerged.find("TT") != std::string::npos && dataSetName.find("TT") != std::string::npos) || */
            (dataSetNameMerged.find("ST") != std::string::npos && dataSetName.find("ST") != std::string::npos) || 
            (extremeMerge_ && dataSetNameMerged.find("Other") != std::string::npos && dataSetName.find("Jets") != std::string::npos) || 
            (! extremeMerge_ && dataSetNameMerged.find("WJets") != std::string::npos && dataSetName.find("W") != std::string::npos && dataSetName.find("Jets") != std::string::npos) || 
            (! extremeMerge_ && dataSetNameMerged.find("DYJets") != std::string::npos && dataSetName.find("DY") != std::string::npos && dataSetName.find("Jets") != std::string::npos) ||
            (dataSetNameMerged.find("data") != std::string::npos && dataSetName.find("data") != std::string::npos) )
        {
          //nofEventsRawMerged_[iCut][dM]      += nofEventsRaw_[iCut][d];   // not allowed unless all constituents have the same eqLumi
          nofEventsMerged_[iCut][dM]    += nofEvents_[iCut][d];
          nofEventsVarMerged_[iCut][dM] += nofEventsVar_[iCut][d];
        }
        else if ( dataSetNameMerged.find("TT") != std::string::npos )
        {
          nofEventsMerged_[iCut][dM]    = nofEvents_[iCut][dM];
          nofEventsVarMerged_[iCut][dM] = nofEventsVar_[iCut][dM];
        }
        
      }  // end d
      
      nofEventsErrorMerged_[iCut][dM] = sqrt(nofEventsVarMerged_[iCut][dM]);
      
    }  // end dM
  }  // end cut
  
}

void SelectionTables::CalculateTable(bool scaleEvents)
{
  std::string dataSetName;
  double lumiWeight;
  
  for (int d = 0; d < numberOfDatasets_; d++)
  {
    dataSetName = listOfDatasets_[d]->Name();
    if ( dataSetName.find("data") != std::string::npos ) lumiWeight = 1.;
    else lumiWeight = lumi_/eqLumi_[d];
    if ( scaleEvents && eqLumi_[d] == 1. ) std::cout << "SelectionTables::CalculateTable: WARNING: eqLumi = 1 for dataset " << dataSetName << ". Please check if this is correct..." << endl;
      
    for (int iCut = 0; iCut < numberOfCuts_; iCut++)
    {
      if ( nofEventsRaw_[0][d] == 0. ) nofEventsRawError_[iCut][d] = 0.;
      else nofEventsRawError_[iCut][d] = ErrorCalculator(nofEventsRaw_[iCut][d], nofEventsRaw_[iCut][d]/nofEventsRaw_[0][d], 1.);
      nofEventsRawVar_[iCut][d]   = pow(nofEventsRawError_[iCut][d], 2.);
      
      if (scaleEvents) nofEvents_[iCut][d] = nofEventsRaw_[iCut][d] * lumiWeight;
      else             nofEvents_[iCut][d] = nofEventsRaw_[iCut][d];
      
      if ( nofEvents_[0][d] == 0. ) nofEventsError_[iCut][d] = 0.;
      else nofEventsError_[iCut][d] = ErrorCalculator(nofEvents_[iCut][d], nofEvents_[iCut][d]/nofEvents_[0][d], 1.);
      nofEventsVar_[iCut][d] = pow(nofEventsError_[iCut][d], 2.);
      
      if ( dataSetName.find("data") == std::string::npos )
      {
        totalEvents_[iCut]    += nofEvents_[iCut][d];
        totalEventsVar_[iCut] += nofEventsVar_[iCut][d];
      }
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
  std::string datasettitle = "";
  
  std::string interline = "";
  for (int d = 0; d < nDatasets; d++) interline += "&";
  interline += "&\\\\[-6pt]";
  std::string headerextra = "[+3pt]";
  
  if(writeLandscape)
  {
    fout << "\\afterpage{" << std::endl;
    fout << "\\begin{landscape}" << std::endl;
  }
  fout << "\\begin{table}" << std::endl;
  fout << "\\caption{Effect of each selection requirement on the number of selected events. ($" << slumi_.str() << "\\fbinv$ of int. lumi.)}" << std::endl;
  fout << "\\label{tab:selectionTable}" << std::endl;
  fout << "\\begin{center}" << std::endl;
  fout << "\\begin{tabular}{r|";
  for (int d = 1; d < nDatasets; d++)
  {
    fout << "c";
  }
  fout << "|c|c}" << std::endl;
  fout << "\\hline" << std::endl;
  fout << "\\hline" << std::endl;
  for (int d = 0; d < nDatasets; d++) fout << "&";
  fout << "&\\\\[-1pt]" << std::endl;
  fout << " & ";
  for (int d = 1; d < nDatasets; d++)
  {
    if (writeMerged)
    {
      datasettitle = listOfDatasetsMerged_[d]->Title();
    }
    else
    {
      datasettitle = listOfDatasets_[d]->Title();
    }
    
    if ( datasettitle.find("t#bar{t}") != std::string::npos )
    {
      boost::replace_all(datasettitle, "t#bar{t}", "$\\ttbar$");
      boost::replace_all(datasettitle, ". ", ".\\ ");
      fout << datasettitle;
    }
    else fout << "$" << datasettitle << "$";
    
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
      if ( writeError && listTableError_[iCut][d] > 0. ) fout << " $\\pm$ " << listTableError_[iCut][d];
      
      if ( d < nDatasets-1 ) fout << " & ";
      else
      {
        fout << " & " << totalEvents_[iCut];
        if ( writeError && totalEventsError_[iCut] > 0. ) fout << " $\\pm$ " << totalEventsError_[iCut];
        fout << " & " << listTable_[iCut][0];
        //if (writeError) fout << " $\\pm$ " << listTableError_[iCut][0];  // no errors on observed no. of events
        fout << " \\\\" << std::endl;
        if ( iCut < numberOfCuts_-1 ) fout << "[+1pt]" << std::endl;
      }
    }
  }
  for (int d = 0; d < nDatasets; d++) fout << "&";
  fout << "&\\\\[-8pt]" << std::endl;
  fout << "\\multicolumn{" << nDatasets+2 << "}{c}{ }\\\\[-6pt]" << endl;
  fout << "\\hline" << std::endl;
  fout << "\\hline" << std::endl;
  fout << "\\end{tabular}" << std::endl;
  fout << "\\end{center}" << std::endl;
  fout << "\\end{table}" << std::endl;
  if (writeLandscape) fout << "\\end{landscape}}%" << std::endl;
}

void SelectionTables::WriteTableVertical(ofstream& fout, double** listTable_,double** listTableError_, bool writeError, bool writeMerged)
{
  int nDatasets = ( writeMerged ? numberOfDatasetsMerged_ : numberOfDatasets_ );
  if( precision_ >= 0 ) fout << fixed << setprecision(precision_);
  std::string datasettitle = "";
  
  std::string interline = "";
  for (int iCut = 0; iCut < numberOfCuts_; iCut++) interline += "&";
  interline += "\\\\[-6pt]";
  std::string headerextra = "[+3pt]";
  
  fout << "\\begin{table}" << std::endl;
  fout << "\\caption{Number of events before and after applying a kinematic fit and requiring $\\chi^2 < 15$ ($" << slumi_.str() << "\\fbinv$ of int. lumi.)}" << std::endl;
  fout << "\\label{tab:KFtable}" << std::endl;
  fout << "\\begin{center}" << std::endl;
  fout << "\\begin{tabular}{l";
  for (int iCut = 0; iCut < numberOfCuts_; iCut++)
  {
    if (! writeError) fout << "r";
    else fout << "c";
  }
  fout << "}" << std::endl;
  fout << "\\hline" << std::endl;
  fout << "\\hline" << std::endl;
  fout << interline << std::endl;
  for (int iCut = 0; iCut < numberOfCuts_; iCut++) fout << " & " << listOfCuts_[iCut];
  fout << " \\\\" << headerextra << std::endl;
  fout << "\\hline" << std::endl;
  fout << interline << std::endl;
  for (int d = 1; d < nDatasets; d++)
  {
    if (writeMerged)
    {
      datasettitle = listOfDatasetsMerged_[d]->Title();
    }
    else
    {
      datasettitle = listOfDatasets_[d]->Title();
    }
    
    if ( datasettitle.find("t#bar{t}") != std::string::npos )
    {
      boost::replace_all(datasettitle, "t#bar{t}", "$\\ttbar$");
      boost::replace_all(datasettitle, ". ", ".\\ ");
      fout << datasettitle;
    }
    else fout << "$" << datasettitle << "$";
    
    for (int iCut = 0; iCut < numberOfCuts_; iCut++)
    {
      fout << " & " << listTable_[iCut][d];
      if ( writeError && listTableError_[iCut][d] > 0. ) fout << " $\\pm$ " << listTableError_[iCut][d];
      else fout << " \\qquad";
      if ( iCut == numberOfCuts_-1 ) fout << " \\\\" << std::endl;
    }
    if ( d < nDatasets-1 ) fout << "[+1pt]" << std::endl;
  }
  fout << interline << std::endl;
  fout << "\\hline" << std::endl;
  fout << interline << std::endl;
  fout << "Total exp.";
  for (int iCut = 0; iCut < numberOfCuts_; iCut++)
  {
    fout << " & " << totalEvents_[iCut];
    if ( writeError && totalEventsError_[iCut] > 0. ) fout << " $\\pm$ " << totalEventsError_[iCut];
    else fout << " \\qquad";
  }
  fout << " \\\\" << headerextra << std::endl;
  fout << "\\hline" << std::endl;
  fout << interline << std::endl;
  fout << "Observed";
  for (int iCut = 0; iCut < numberOfCuts_; iCut++)
  {
    fout << " & " << listTable_[iCut][0];
    //if (writeError) fout << " $\\pm$ " << listTableError_[iCut][0];  // no errors on observed no. of events
    //else
    fout << " \\qquad";
  }
  fout << " \\\\" << headerextra << std::endl;
  fout << "\\hline" << std::endl;
  fout << "\\hline" << std::endl;
  fout << "\\end{tabular}" << std::endl;
  fout << "\\end{center}" << std::endl;
  fout << "\\end{table}%" << std::endl;
}

void SelectionTables::Write(std::string filename, bool writeError, bool writeMerged, bool writeLandscape, bool writeVertical)
{
  ofstream fout(filename.c_str());
  this->Write(fout, writeError, writeMerged, writeLandscape, writeVertical);
}

void SelectionTables::Write(ofstream& fout, bool writeError, bool writeMerged, bool writeLandscape, bool writeVertical)
{
  if (writeVertical)
  {
    if (writeMerged) this->WriteTableVertical(fout, nofEventsMerged_, nofEventsErrorMerged_, writeError, writeMerged);
    else this->WriteTableVertical(fout, nofEvents_, nofEventsError_, writeError, writeMerged);
  }
  else
  {
    if (writeMerged) this->WriteTable(fout, nofEventsMerged_, nofEventsErrorMerged_, writeError, writeMerged, writeLandscape);
    else this->WriteTable(fout, nofEvents_, nofEventsError_, writeError, writeMerged, writeLandscape);
  }
}
