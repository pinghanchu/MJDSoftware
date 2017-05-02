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
    cout << "Usage: " << argv[0] << " [startrun] [endrun]" << endl;
    return 1;
  }

  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  GATAutoCal ds(fStartRun,fEndRun);
  //ds.SetEnergyName(fEnergyName);
  string fDataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();

  Double_t fEnr1 = 66.7;
  Double_t fEnr2 = 10000;
  Double_t fTime = 0.5;
  string fInputFile = Form("./data/wf_%d_%d.txt",fStartRun,fEndRun);

  ifstream fin(Form("%s",fInputFile.c_str()));

  Int_t run,entry,channel;
  Double_t enr,ratio,deltaT,maxFFT,aovere;
  vector<Double_t> Ratio;
  vector<Double_t> DeltaT;
  vector<Double_t> Energy;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run >> entry >> channel >> enr >> ratio >> deltaT >> maxFFT >> aovere ;
      if(enr<72 && enr>62 && ratio > 3.7 && ratio <4.5){
	Ratio.push_back(ratio);
	DeltaT.push_back(deltaT);
	Energy.push_back(enr);
	if(deltaT>2000){
	  cout<<run<<" " <<entry << " "<< channel << endl;
	}
      }
    }
  }

  TCanvas *c1 = new TCanvas("c1");
  if(Ratio.size()>0){
    TH1D *h1 = new TH1D("h1","",10000,0,10000);
    for(size_t i = 0;i<DeltaT.size();i++){
      h1->Fill(DeltaT.at(i));
    }
    h1->SetTitle(";#Delta T(ns)");
    h1->GetYaxis()->SetTitleOffset(1.3);
    h1->Draw();
    c1->Print(Form("deltaT_%d_%d.pdf",fStartRun,fEndRun));
  }
}
