#include "MJDGat.hh"
#include "GATPulserTag.hh"
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
  if(argc != 2) {
    cout << "Usage: " << argv[0] << " [run]" << endl;
    return 1;
  }
  
  Int_t fRun = atoi(argv[1]);
  
  MJDGat ds(fRun);
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDet = ds.GetDetPosition();
  string FileName;
  //TChain *fMjdTree = new TChain("mjdTree");
  TCanvas *c1 = new TCanvas("c1");
  c1->SetLogy();
  ofstream fout(Form("pulser_%d.txt",fRun),ios::app);
  TChain *fMjdTree = ds.GetMJDTree();
  vector<Double_t> fPulser;

  for(size_t i = 0;i<fChannel.size();i++){
    Int_t fPos = fCryo.at(i)*100+fStr.at(i)*10+fDet.at(i);
    TH1D *h = new TH1D(Form("%s","trapE"),"", 8000,0,8000);
    fMjdTree->Draw(Form("trapE>>%s","trapE"),Form("channel == %d",fChannel.at(i)));
    h->SetTitle(";Energy(ADC);Counts");
    h->Rebin(10);
    h->GetXaxis()->SetRangeUser(50,6000);
    h->Draw();
    c1->Print(Form("h_%d_%d.pdf",fPos,fChannel.at(i)));
    Double_t maxbin = h->GetMaximumBin();
    Double_t ymax = h->GetMaximum();
    Double_t xmax = h->GetBinCenter(maxbin);
    if(ymax>5){
      fout << fRun << " " << fPos << " " << fChannel.at(i) << " " << xmax << " " << ymax  <<endl;
      fPulser.push_back(xmax);
    }else{
      fout << fRun << " " << fPos << " " << fChannel.at(i) << " " << 0 << " " << 0 <<endl;
      fPulser.push_back(0);
    }
    delete h;
  }
  for(size_t i=0;i<fPulser.size();i++){
    cout << fPulser.at(i) << ", ";
  }
}
