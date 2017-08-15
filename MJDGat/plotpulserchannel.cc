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
  if(argc != 2 ) {
    cout << "Usage: " << argv[0] << " [data set]" << endl;
    return 1;
  }
  Int_t fDataSet = atoi(argv[1]);
  Int_t posM1[29] = {111,112,113,114,121,122,123,124,131,132,133,134,141,142,143,144,145,151,152,153,154,161,162,163,164,171,172,173,174};
  Int_t posM2[29] = {211,212,213,214,221,222,223,224,225,231,232,233,241,242,243,244,245,251,252,253,254,261,262,263,264,271,272,273,274};
  vector<string> Det;
  vector<Int_t> Pos;
  if(fDataSet == 0 || fDataSet == 1 || fDataSet == 2 || fDataSet ==3){
    for(int ip = 0;ip<29;ip++){
      Pos.push_back(posM1[ip]*10+0);
      Pos.push_back(posM1[ip]*10+1);
    }
  }else if(fDataSet == 4){
    for(int ip = 0;ip<29;ip++){
      Pos.push_back(posM2[ip]*10+0);
      Pos.push_back(posM2[ip]*10+1);
    }
  }else if(fDataSet >= 5){
    for(int ip = 0;ip<29;ip++){
      Pos.push_back(posM1[ip]*10+0);
      Pos.push_back(posM1[ip]*10+1);
    }
    for(int ip = 0;ip<29;ip++){
      Pos.push_back(posM2[ip]*10+0);
      Pos.push_back(posM2[ip]*10+1);
    }
  }
  const Int_t nPos = Pos.size();
  vector<Double_t> Value[nPos];
  vector<Int_t> Run;
  for(size_t i = 0;i<Pos.size();i++){
    ifstream fin(Form("pulser_%d_%d.txt",fDataSet,Pos.at(i)));
    Int_t run,pos,channel;
    Double_t count;
    if(fin.is_open()){
      while(!fin.eof()){
	fin >> run >> pos >> channel >> count ;
	Value[i].push_back(count);
	if(i ==0){
	  Run.push_back(run);
	}
      }
    }
  }
  Run.pop_back();
  for(size_t i=0;i<Pos.size();i++){
    Value[i].pop_back();
  }
  Int_t nRun = Run.size();
  int RunRange = Run.at(nRun-1)-Run.at(0);
  vector<Double_t> ValueNew[nPos];
  for(size_t i=0;i<Pos.size();i++){
    for(int irun = Run.at(0);irun< Run.at(nRun-1);irun++){
      Int_t match=-1;
      for(size_t j=0;j<Run.size();j++){
	if(irun == Run.at(j)){
	  match = j;
	}
      }
      if(match != -1){
	ValueNew[i].push_back(Value[i].at(match));      
      }else{
	ValueNew[i].push_back(0);
      }
    }
  }

  TH2D *h1 = new TH2D("h1","",RunRange,Run.at(0), Run.at(nRun-1),nPos,0,nPos);

  for(size_t i = 0;i<Pos.size();i++){
    for(Int_t ij = 0;ij<RunRange;ij++){
      h1->SetBinContent(ij,i,ValueNew[i].at(ij));

    }
  }
  h1->SetTitle(";Run;Channel(Detector ID)");
  h1->GetYaxis()->SetTitleOffset(1.3);
  TCanvas *c1 = new TCanvas("c1");
  gStyle->SetOptStat(0);

  h1->Draw("COLZ");
  c1->Print(Form("pulserchannel_%d.pdf",fDataSet));
}
