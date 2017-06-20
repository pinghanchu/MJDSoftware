#include "MJDGat.hh"
#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run]" << endl;
    return 1;
  }

  Int_t fRun = atoi(argv[1]);
  GATAutoCal ds(fRun,fRun);
  string fDataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  TTree *mjdtree=ds.GetMJDTree();
  Double_t starttime = ds.GetStartTime();
  Double_t endtime = ds.GetStopTime();
  Double_t time = (endtime-starttime)/1e8;
  Int_t count = mjdtree->GetEntries();
  cout.precision(15);
  cout << fRun << " " << count << " "<< starttime << " " << endtime << " " << time << endl;
}
