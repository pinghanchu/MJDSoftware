#include "MJDGat.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 3) {
    //cout << "Usage: " << argv[0] << " [dataset] [subset1] [subset2] [iscal]" << endl;
    
    cout << "Usage: " << argv[0] << " [dataset] [run]" << endl;
  }

  Int_t fDataSet = atoi(argv[1]);
  Int_t fRun = atoi(argv[2]);
  //Int_t fSubSet = atoi(argv[2]);
  //Int_t fSubSet2 = atoi(argv[3]);
  //Int_t fIsCal = atoi(argv[4]);
  Int_t fIsCal = 1;
  TCanvas *c1 = new TCanvas("c1");

  MJDGat ds(fRun);
  string DataSetName = ds.GetDataSet();
  //TChain* skimTree = ds.GetSkimTree();
  TChain* mjdTree = new TChain("mjdTree");
  mjdTree->Add(Form("/global/project/projectdirs/majorana/data/mjd/surfmjd/data/gatified/%s/mjd_run%d.root",DataSetName.c_str(),fRun));
  //skimTree->Add(Form("./cal-skim/skimDS%d_run%d_low.root",fDataSet,fSubSet));
  string FileName;
  FileName = Form("hist_%d_%d.root",fDataSet,fRun);
  

  TFile fhist(Form("%s",FileName.c_str()),"update");
  string fEName = "trapECal";
  
  TH1D* h2 = ds.FillHisto(mjdTree, Form("%s",fEName.c_str()),Form("%s00",fEName.c_str()),"channel%2== 0 && isNat == 1",300000,0,3000);
  h2->Write();
 
  TH1D* h1 = ds.FillHisto(mjdTree, Form("%s",fEName.c_str()),Form("%s01",fEName.c_str()),"channel%2== 1 && isNat == 1",300000,0,3000);
  h1->Write();

  TH1D* h4 = ds.FillHisto(mjdTree, Form("%s",fEName.c_str()),Form("%s10",fEName.c_str()),"channel%2== 0 && isEnr == 1",300000,0,3000);
  h4->Write();


  TH1D* h3 = ds.FillHisto(mjdTree, Form("%s",fEName.c_str()),Form("%s11",fEName.c_str()),"channel%2== 1 && isEnr == 1",300000,0,3000);
  h3->Write();



  TH1D* h6 = ds.FillHisto(mjdTree, Form("%s",fEName.c_str()),Form("%sHG",fEName.c_str()),"channel%2== 0 ",300000,0,3000);
  h6->Write();

  TH1D* h5 = ds.FillHisto(mjdTree, Form("%s",fEName.c_str()),Form("%sLG",fEName.c_str()),"channel%2== 1 ",300000,0,3000);
  h5->Write();



  vector<Int_t> fChannel = ds.GetChannel();
  //cout << fChannel.size() << endl;
  for(size_t i=0;i<fChannel.size();i++){
    //cout << fChannel.at(i) << endl;
    TH1D* h = ds.FillHisto(mjdTree,"trapECal",Form("trapECal%d",fChannel.at(i)),Form("channel == %d",fChannel.at(i)),300000,0,3000);
    //TH1D* h = ds.FillHisto(skimTree,"trapECal",Form("trapECal%d",fChannel.at(i)),Form("channel == %d",fChannel.at(i)),300000,0,3000);

    h->Write();
    delete h;
  }

}
