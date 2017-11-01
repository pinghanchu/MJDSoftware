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
  if(argc != 3 ) {
    cout << "Usage: " << argv[0] << " [dataset] [subset]" << endl;
    return 1;
  }
  Int_t fDataSet = atoi(argv[1]);
  Int_t fSubSet = atoi(argv[2]);
  Double_t fEnergy = 215;
  Double_t fTime = 52.9;
  MJDSkim ds(fDataSet,fSubSet,fSubSet,0);
  Int_t fTimems = 53;
  string fOutputFile = Form("77mGe_MuonCoin_data_%d_%d_%d.txt",fDataSet,(Int_t)fEnergy,fTimems);
  ds.SearchMuonCoinEvent(fEnergy,fTime,fOutputFile);

  fOutputFile = Form("77mGe_HighEnr_data_%d_%d_%d_%d.txt",fDataSet,fSubSet,(Int_t)fEnergy,fTimems);
  ds.SearchDelayedHighEnergyEvent(1000,1000,fTime,fOutputFile);

  fOutputFile = Form("77mGe_Multiplicity_data_%d_%d_%d_%d.txt",fDataSet,fSubSet,(Int_t)fEnergy,fTimems);
  ds.SearchDelayedMultiplicityEvent(fEnergy,10000,fTime,fOutputFile);
  
  cout << "Searching candidates is done. "<< endl;
  

  
}
