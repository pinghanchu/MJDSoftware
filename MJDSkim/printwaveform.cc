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
  if(argc != 5 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run] [entry] [channel] [energy]" << endl;
    return 1;
  }

  Int_t fRun = atoi(argv[1]);
  cout << "Run:" << fRun << endl;
  Int_t fEntry = atoi(argv[2]);
  Int_t fChan = atoi(argv[3]);
  Double_t fEnr = atof(argv[4]);
  TCanvas *c1 = new TCanvas("c1");
  MJDSkim ds(0,0,0);
  TH1D* h = ds.GetWaveform(fRun,fEntry,fChan,fEnr);
  TH1D* h2 = ds.GetHistoDerivative(h,10);
  TH1D* h4 = ds.GetHistoDerivative(h2,10);
  TH1*  hFFT = ds.GetHistoFFT(h);
  h->Draw();
  c1->Print(Form("waveform_%d_%d_%d.pdf",fRun,fEntry,fChan));
  h2->Draw();
  c1->Print(Form("waveform1_%d_%d_%d.pdf",fRun,fEntry,fChan));
  h4->Draw();
  c1->Print(Form("waveform2_%d_%d_%d.pdf",fRun,fEntry,fChan));
  c1->SetLogy();
  hFFT->GetXaxis()->SetRangeUser(0,1000);
  hFFT->Draw();
  c1->Print(Form("waveformFFT_%d_%d_%d.pdf",fRun,fEntry,fChan));
  c1->SetLogy(0);
  TH1D *htrap =  ds.TrapezoidalFilter(h, 400,250,7000);
  htrap->Draw();
  c1->Print(Form("waveformTrap_%d_%d_%d.pdf",fRun,fEntry,fChan));

  /*

  //h3->SetName(Form("waveform1_%d_%d_%d",fRun,fEntry,fChan));
  //h5->SetName(Form("waveform2_%d_%d_%d",fRun,fEntry,fChan));
  //hFFT->SetName(Form("waveformFFT_%d_%d_%d",fRun,fEntry,fChan));
  hFFT->GetXaxis()->SetRangeUser(0,1000);

  Double_t maxFFT = ds.GetMax(hFFT, 50.,150.);
  //Double_t maxY0 = ds.GetMax(h, 100,19000);
  h->GetXaxis()->SetRangeUser(100,19000);
  Int_t maxBin0 = h->GetMaximumBin();
  Double_t maxX0 = h->GetBinCenter(maxBin0);

  Double_t Xmin = 4000;
  Double_t Xmax = 19000;
  if(fRun>=14503 && fRun<=15892){
    Xmax = 38000;
  }
  Double_t fResolution = 1;
  Double_t fThreshold = 0.01;
  Double_t fSigma = 5;
  Double_t A = ds.GetMax(h3, Xmin,Xmax);
  Double_t AoverE = A/fEnr;

  vector<Double_t> xp;
  vector<Double_t> yp;
  vector<Double_t> xp1;
  vector<Double_t> yp1;
  vector<Double_t> xp2;
  vector<Double_t> yp2;
  xp.clear();
  yp.clear();
  xp1.clear();
  yp1.clear();
  xp2.clear();
  yp2.clear();

  Int_t nPeak1 = ds.FindPeaks(h3, Xmin, Xmax,fResolution,fSigma, fThreshold, &xp, &yp);
  cout << fRun << " " << fEntry << " " << fChan << " " << fEnr << " " << nPeak1 << " " << maxX0 << endl;
  ofstream fout(Form("wf_%d.txt",fRun),ios::app);
  if(nPeak1>1){
    for(Int_t ip = 0;ip<(Int_t)xp.size();ip++){
      Double_t y = ds.GetYValue(h5, xp.at(ip));
      cout << ip << " " << xp.at(ip) << " "<<yp.at(ip)<< " "<<y <<endl;
      if(abs(y)<5e-3 && xp.at(ip)< maxX0 ){
	xp1.push_back(xp.at(ip));
	yp1.push_back(yp.at(ip));
      }
    }
    vector<Int_t> Index2 = ds.Sort(yp1);
    for(size_t ip1 = Index2.size()-2;ip1<Index2.size();ip1++){
      Int_t ii = Index2.at(ip1);
      xp2.push_back(xp1.at(ii));
      yp2.push_back(yp1.at(ii));
      cout <<ii << " "<< yp1.at(ii) << " " << yp1.at(ii) << endl;
    }
    if(xp2.size()>0){
      fout << fRun << " " << fEntry << " " << fChan << " "<< fEnr << " " 
	   << 1 << " " << yp2.at(1)/yp2.at(0) << " " << xp2.at(0)-xp2.at(1) << " " 
	   << maxFFT << " " << AoverE << endl;
    }
    
    h1->Draw();
    c1->Print(Form("waveform_%d_%d_%d.pdf",fRun,fEntry,fChan));
    h3->Draw();
    c1->Print(Form("waveform1_%d_%d_%d.pdf",fRun,fEntry,fChan));
    h5->Draw();
    c1->Print(Form("waveform2_%d_%d_%d.pdf",fRun,fEntry,fChan));
    c1->SetLogy();
    hFFT->GetXaxis()->SetRangeUser(0,1000);
    hFFT->Draw();
    c1->Print(Form("waveformFFT_%d_%d_%d.pdf",fRun,fEntry,fChan));
  }
  */
  /*
  string fOutputWaveform = Form("waveform_%d.root",fRun);
  TFile fhist(Form("%s",fOutputWaveform.c_str()),"recreate");

  h3->SetName(Form("waveform1_%d_%d_%d",fRun,fEntry,fChan));
  h5->SetName(Form("waveform2_%d_%d_%d",fRun,fEntry,fChan));
  hFFT->SetName(Form("waveformFFT_%d_%d_%d",fRun,fEntry,fChan));
  h->Write();
  //h1->Write();
  h3->Write();
  h5->Write();
  hFFT->Write();
  */
}
