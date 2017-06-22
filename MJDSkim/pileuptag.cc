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
  if(argc != 3 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [data set] [subset]" << endl;
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  Int_t fSubSet = atoi(argv[2]);
  const char* fPileUpPath = "";
  MJDSkim ds(fDataSet,fSubSet);
  ds.PileUpTree(fPileUpPath);
}