#include "MJDGat.hh"
#include "GATAutoCal.hh"
#include "GATDataSet.hh"
#include "MJTRun.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <bitset>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run]" << endl;
    return 1;
  }
  Int_t fRun = atoi(argv[1]);
  GATDataSet ds(fRun);

  TChain *c = ds.GetBuiltChain();
  MJTRun *runInfo = new MJTRun();
  c->SetBranchAddress("run",&runInfo);
  c->GetEntry(0);

  bitset<32> event_type = runInfo->GetRunBits();
  //cout << "event type = " << event_type << endl;

  if (event_type.test(12) || event_type.test(13) ){
    cout << fRun << " " << 1 << endl;
    return 1;
  }
  else{
    cout << fRun << " " << 0 << endl;
    return 0;
  }
}
