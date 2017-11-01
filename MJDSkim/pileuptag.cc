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
  if(argc != 4) {
    cout << "Usage: " << argv[0] << " [data set] [subset] [iscal]" << endl;
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  Int_t fSubSet = atoi(argv[2]);
  Int_t fIsCal = atoi(argv[3]);
  const char* fPileUpPath = "";
  MJDSkim ds(fDataSet,fSubSet,fSubSet,fIsCal);
  ds.PileUpTree(fPileUpPath);
  
  return 1;
}
