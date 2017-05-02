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
  if(argc != 3 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << "[energy] [window]" << endl;
    return 1;
  }

  Double_t fEnr = atof(argv[1]);
  Double_t fWindow =atof(argv[2]);
  string fInputFile = Form("./data/wf_%d_%d.txt",(Int_t)fEnr,(Int_t)fWindow);

  ifstream fin(Form("%s",fInputFile.c_str()));

  Int_t run,entry,pos,channel;
  Double_t enr,ratio,deltaT,maxFFT,aovere;
  vector<Double_t> Ratio;
  vector<Double_t> DeltaT;
  vector<Double_t> Energy;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run >> entry >> channel >> enr >> ratio >> deltaT >> maxFFT >> aovere ;
      if(ratio > 3.7 && ratio <4.5){
	Ratio.push_back(ratio);
	DeltaT.push_back(deltaT);
	Energy.push_back(enr);
	cout<<run<<" " <<entry << " "<< channel << endl;
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
    c1->Print(Form("deltaT_%d_%d.pdf",(Int_t)fEnr,(Int_t)fWindow));
  }
}
