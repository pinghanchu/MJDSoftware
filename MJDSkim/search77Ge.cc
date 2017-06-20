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
  if(argc != 4 ) {
    cout << "Usage: " << argv[0] << " [dataset] [energy] [delay time]" << endl;
    return 1;
  }
  Int_t fDataSet = atoi(argv[1]);
  Double_t fEnergy = atof(argv[2]);
  Double_t fTime = atof(argv[3]);
  MJDSkim ds(fDataSet);
  Double_t fQ = 10000;
  Int_t fTimems = (Int_t) 114;
  string fOutputFile = Form("data_%d_%d_%d.txt",fDataSet,(Int_t)fEnergy,fTimems);
  cout << "Searching candidates..." << "energy is " << fEnergy << "; the delayed time is " <<fTime << endl;
  ds.SearchDelayedEvent(fEnergy,fQ,fTime,fOutputFile);
  cout << "Searching candidates is done. "<< endl;
}
