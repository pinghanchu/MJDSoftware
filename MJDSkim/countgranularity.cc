#include "MJDSkim.hh"
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
  if(argc != 2) {
    cout << "Usage: " << argv[0] << " [dataset]" << endl;
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  
  MJDSkim ds(fDataSet);

  TTree *skimtree = ds.GetSkimTree();
  Int_t mHClean;
  skimtree->SetBranchAddress("mHClean", &mHClean);
  Int_t count = 0;
  for(Int_t i=0;i<(Int_t)skimtree->GetEntries();i++){
    skimtree->GetEntry(i);
    if(mHClean > 1){
      count++;
    }
  }
  cout << count << endl;
}
