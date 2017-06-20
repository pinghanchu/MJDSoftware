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
  if(argc != 3 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run] [pulser path]" << endl;
    return 1;
  }

  Int_t fRun = atoi(argv[1]);
  string fPulserPath = argv[2];
  GATAutoCal ds(fRun,fRun);
  string fDataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();

  TFile *pulserfile = new TFile(Form("%spulser_%d.root", fPulserPath.c_str(), fRun), "read");
  TTree *pulsertree = (TTree*)pulserfile->Get("pulsertree");
  Int_t Pulser;
  pulsertree->SetBranchAddress("Pulser", &Pulser);
  Int_t count = 0;
  for(Int_t i=0;i<(Int_t)pulsertree->GetEntries();i++){
    pulsertree->GetEntry(i);
    if(Pulser == 1){
      count++;
    }
  }
  cout << fRun << " " << count << endl;

}
