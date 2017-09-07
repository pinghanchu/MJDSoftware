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
  if(argc != 5) {
    cout << "Usage: " << argv[0] << " [dataset] [subset1] [subset2] [iscal]" << endl;
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  Int_t fSubSet = atoi(argv[2]);
  Int_t fSubSet2 = atoi(argv[3]);
  Int_t fIsCal = atoi(argv[4]);
  TCanvas *c1 = new TCanvas("c1");
  MJDSkim ds(fDataSet,fSubSet,fSubSet2,fIsCal);
  TChain* skimTree = ds.GetSkimTree();
  string FileName;
  if(fSubSet == fSubSet2){
    FileName = Form("hist_%d_%d.root",fDataSet,fSubSet);
  }else{
    FileName = Form("hist_%d_%d_%d.root",fDataSet,fSubSet,fSubSet2);
  }

  TFile fhist(Form("%s",FileName.c_str()),"update");
  string fEName = "trapENFCalC";

  
  TH1D* h2 = ds.FillHisto(skimTree, Form("%s",fEName.c_str()),Form("%s00",fEName.c_str()),"channel%2== 0 && isGood==1 && isNat == 1",30000,0,3000);
  h2->Write();
 
  TH1D* h1 = ds.FillHisto(skimTree, Form("%s",fEName.c_str()),Form("%s01",fEName.c_str()),"channel%2== 1 && isGood==1 && isNat == 1",30000,0,3000);
  h1->Write();

  TH1D* h4 = ds.FillHisto(skimTree, Form("%s",fEName.c_str()),Form("%s10",fEName.c_str()),"channel%2== 0 && isGood==1 && isEnr == 1",30000,0,3000);
  h4->Write();


  TH1D* h3 = ds.FillHisto(skimTree, Form("%s",fEName.c_str()),Form("%s11",fEName.c_str()),"channel%2== 1 && isGood==1 && isEnr == 1",30000,0,3000);
  h3->Write();



  TH1D* h6 = ds.FillHisto(skimTree, Form("%s",fEName.c_str()),Form("%sHG",fEName.c_str()),"channel%2== 0 && isGood==1",30000,0,3000);
  h6->Write();

  TH1D* h5 = ds.FillHisto(skimTree, Form("%s",fEName.c_str()),Form("%sLG",fEName.c_str()),"channel%2== 1 && isGood==1",30000,0,3000);
  h5->Write();



  vector<Int_t> fChannel = ds.GetChannel();
  cout << fChannel.size() << endl;
  for(size_t i=0;i<fChannel.size();i++){
    cout << fChannel.at(i) << endl;
    TH1D* h = ds.FillHisto(skimTree,"trapENFCalC",Form("trapENFCalC%d",fChannel.at(i)),Form("channel == %d",fChannel.at(i)),300000,0,3000);
    h->Write();
    delete h;
  }

}
