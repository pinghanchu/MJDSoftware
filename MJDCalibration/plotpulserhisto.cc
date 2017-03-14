// 2016.2.18
// written by Pinghan Chu
// Following the logic in GATAutoCal.cc
// 
#include "GATAutoCal.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TSpectrum.h"
#include <string>
#include <stdlib.h>

Double_t Gaus0(Double_t *v, Double_t *par){
  Double_t arg = 0;
  if (par[2] != 0) arg = (v[0] - par[1])/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}

Double_t FinMaxX(vector<Double_t>* fPx, vector<Double_t>* fPy){
  Double_t temp = 0;
  Double_t maxX = 0;
  for(size_t i = 0;i<fPx->size();i++){
    if(fPy->at(i) > temp){
      maxX = fPx->at(i);
      temp = fPy->at(i);
    }
  }
  return maxX;
}

void FindPeak(TH1D *h, vector<Double_t>* fPx, vector<Double_t>* fPy){
  TSpectrum *s = new TSpectrum(1000,1);
  Int_t nfound = s->Search(h,10,"new",0.0001);
  Double_t *xpeaks = s->GetPositionX();
  Double_t *ypeaks = s->GetPositionY();
  for(Int_t i = 0;i<nfound; i++){
    if(ypeaks[i]>3){
      fPx->push_back(xpeaks[i]);
      fPy->push_back(ypeaks[i]);
    }
  }

  const Int_t npeak = fPx->size();
  Double_t Px[npeak];
  Double_t Py[npeak];
  Int_t n[npeak];
  for(Int_t i = 0;i<npeak;i++){
    Px[i] = fPx->at(i);
    Py[i] = fPy->at(i);
  }
  TMath::Sort(npeak,Px,n,0);
  fPx->clear();
  fPy->clear();

  for(Int_t i = 0;i<npeak;i++){
    Int_t ii = n[i];
    fPx->push_back(Px[ii]);
    fPy->push_back(Py[ii]);
  }
}


TH1D* GetHisto(string FileName, string HistoName){
  TFile *fhist = TFile::Open(Form("%s",FileName.c_str()),"read");
  TH1D *h = (TH1D*)fhist->Get(Form("%s",HistoName.c_str()));
  return h;
}


int main(int argc, char** argv)
{
  if(argc != 4 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [energy name]" << endl;
    return 1;
  }
  int fStartRun = atoi(argv[1]);
  int fEndRun = atoi(argv[2]);
  const char* EnergyName = argv[3];
  string fEName(EnergyName);

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  TCanvas *c1 = new TCanvas("c1");
  c1->SetLogy();
  string FileName = Form("./Hist/pulser/pulser_%d_%d.root",fStartRun,fEndRun);
  string HistoName = "pulserdeltaHG";
  TH1D *h1 = GetHisto(FileName, HistoName);
  h1->SetTitle(";#Delta E(ADC);Counts");
  h1->Draw();
  c1->Print(Form("./Plot/%s_%d_%d.pdf",HistoName.c_str(),fStartRun,fEndRun));
  h1->GetXaxis()->SetRangeUser(-5,5);
  //c1->SetLogy(0);
  TF1 *f1 = new TF1("f1",Gaus0,-5,5,3);
  f1->SetParameter(1,h1->GetMean());
  f1->SetParameter(2,h1->GetRMS());
  h1->Fit(f1);
  h1->Draw();
  c1->Print(Form("./Plot/%szoom_%d_%d.pdf",HistoName.c_str(),fStartRun,fEndRun));
  ofstream fout("./List/pulsertable.txt",ios::app);
  fout << fStartRun << " " << fEndRun << " " 
       << h1->GetMean() << " " << h1->GetRMS() << " " 
       << h1->GetEffectiveEntries() << " " << h1->GetRMS()/sqrt(h1->GetEffectiveEntries()) << " "
       << f1->GetParameter(1) << " " << f1->GetParError(1) << " "  
       << f1->GetParameter(2) <<  " " << f1->GetParError(2) << endl;

  HistoName = "pulserdeltaLG";
  TH1D *h2 = GetHisto(FileName, HistoName);
  h2->SetTitle(";#Delta E(ADC);Counts");
  h2->Draw();
  c1->Print(Form("./Plot/%s_%d_%d.pdf",HistoName.c_str(),fStartRun,fEndRun));
  h2->GetXaxis()->SetRangeUser(-5,5);
  //c1->SetLogy(0);
  //TF1 *f1 = new TF1("f1",Gaus0,-5,5,3);
  f1->SetParameter(1,h2->GetMean());
  f1->SetParameter(2,h2->GetRMS());
  h2->Fit(f1);
  h2->Draw();
  c1->Print(Form("./Plot/%szoom_%d_%d.pdf",HistoName.c_str(),fStartRun,fEndRun));
  //ofstream fout("./List/pulsertable.txt",ios::app);
  fout << fStartRun << " " << fEndRun << " "
       << h2->GetMean() << " " << h2->GetRMS() << " "
       << h2->GetEffectiveEntries() << " " << h2->GetRMS()/sqrt(h2->GetEffectiveEntries()) << " "
       << f1->GetParameter(1) << " " << f1->GetParError(1) << " "
       << f1->GetParameter(2) <<  " " << f1->GetParError(2) << endl;



}


