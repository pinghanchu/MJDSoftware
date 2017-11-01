#include "MJDGat.hh"
#include "GATPulserTag.hh"
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
  MJDGat ds(fRun);

  //string fDataSet = ds.GetDataSet();
  //vector<Int_t> fChannel = ds.GetChannel();
  //vector<Int_t> fCryo = ds.GetCryo();
  //vector<Int_t> fStr = ds.GetString();
  //vector<Int_t> fDetpos = ds.GetDetPosition();
  //TChain* fMjdTree = ds.GetMJDTree();
  //ds.SaveHisto();
  ds.SearchEvents();
}
