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
    cout << "Usage: " << argv[0] << " [dataset] [subset] [iscal]" << endl;
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  Int_t fSubSet = atoi(argv[2]);
  Int_t fIsCal = atoi(argv[3]);
  TCanvas *c1 = new TCanvas("c1");
  MJDSkim ds(fDataSet,fSubSet,fIsCal);
  TChain* skimTree = ds.GetSkimTree();
  TFile fhist(Form("hist_%d_%d.root",fDataSet,fSubSet),"update");
  TH1D* h1 = ds.FillHisto(skimTree, "trapENFCal","trapENFCal00","channel%2==0 && isGood==1 && isNat == 1",30000,0,3000);
  h1->Write();

  TH1D* h2 = ds.FillHisto(skimTree, "trapENFCalC","trapENFCalC00","channel%2== 0 && isGood==1 && isNat == 1",30000,0,3000);
  h2->Write();

  TH1D* h3 = ds.FillHisto(skimTree, "trapENFCal","trapENFCal10","channel%2==0 && isGood==1 && isEnr == 1",30000,0,3000);
  h3->Write();

  TH1D* h4 = ds.FillHisto(skimTree, "trapENFCalC","trapENFCalC10","channel%2== 0 && isGood==1 && isEnr == 1",30000,0,3000);
  h4->Write();

  TH1D* h5 = ds.FillHisto(skimTree, "trapENFCal","trapENFCalHG","channel%2==0 && isGood==1",30000,0,3000);
  h5->Write();

  TH1D* h6 = ds.FillHisto(skimTree, "trapENFCalC","trapENFCalCHG","channel%2== 0 && isGood==1",30000,0,3000);
  h6->Write();


  //h1->Draw();
  //c1->Print("Energy.pdf");
}
