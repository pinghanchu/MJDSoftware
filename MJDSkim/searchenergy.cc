#include "MJDSkim.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 6 ) {
    cout << "Usage: " << argv[0] << " [dataset] [subset] [iscal] [energy] [window]" << endl;
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  Int_t fSubSet = atoi(argv[2]);
  Int_t fIsCal = atoi(argv[3]);

  Double_t fEnergy = atof(argv[4]);
  Double_t fWindow = atof(argv[5]);

  MJDSkim ds(fDataSet,fSubSet,fIsCal);
  string fOutputFile = Form("data_%d_%d.txt",fDataSet,fSubSet);
  ds.SearchEnergyEvent(fEnergy,fWindow,fOutputFile);
  cout << "The scan is done!" << endl;
}
