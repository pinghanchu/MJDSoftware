// 2016.11.21 Pinghan Chu
// Fill histograms from the files on PDSF using GATDataSet();
#include "TStyle.h"
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
Double_t Gaus(Double_t *v, Double_t *par){
  Double_t arg = 0;
  arg = (v[0] - (par[1]))/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg)/par[2]/TMath::Pi();
  return fitval;
}
Double_t Slope(Double_t *v, Double_t *par){
  Double_t fitval = par[0]+par[1]*v[0];
  return fitval;
}

Double_t GausSlope(Double_t *v, Double_t *par){
  return Gaus(v,par)+Slope(v,&par[3]);
}

Double_t TwoGausSlope(Double_t *v, Double_t *par){
  return Gaus(v,par)+Gaus(v,&par[3])+Slope(v,&par[6]);
}


Double_t Special(Double_t *v, Double_t *par){
  Double_t arg1 = 0;
  arg1 = (v[0] - (par[1]))/par[2];

  Double_t fitval1 = par[0]*TMath::Exp(-0.5*arg1*arg1)/par[2]/TMath::Pi();

  Double_t arg2 = 0;
  arg2 = (v[0] - (par[1]*79.6139/80.9971))/par[2];

  Double_t fitval2 = par[0]*(2.62/34.06)*TMath::Exp(-0.5*arg1*arg1)/par[2]/TMath::Pi();

  Double_t slope = par[3]+par[4]*v[0];

  Double_t fitval = fitval1+fitval2+slope;

  return fitval;


  return Gaus(v,par)+Gaus(v,&par[3])+Slope(v,&par[6]);
}

Double_t SkewGaus(Double_t *v, Double_t *par){
  Double_t arg = 0;
  if (v[0]>=(par[1])){
    arg = (v[0] - (par[1]))/par[2];
  }
  if (v[0]<(par[1])){
    arg = (v[0] - (par[1]))/par[3];
  }

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg)/(0.5*par[2]+0.5*par[3])/TMath::Pi();
  return fitval;
}

Double_t SkewGausSlope(Double_t *v, Double_t *par){
  return SkewGaus(v,par)+Slope(v,&par[4]);
}


void Fit()
{

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1111);

  ifstream fin("./inputdata.txt");     
  ofstream fout("./output.txt");
  Double_t r1 = 0;
  TCanvas *c1 = new TCanvas("c1");
  vector<Double_t> R1;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> r1;
      R1.push_back(r1);
    }
  }
  R1.pop_back();
  Int_t entries = R1.size();
  TH1D *h1 = new TH1D("h1","",entries,0,entries);
  for(size_t i=0;i<R1.size();i++){
    h1->SetBinContent(i,R1.at(i));
  }
  h1->SetTitle(";Energy(ADC);Counts");
  h1->Draw();
  c1->Print("data.pdf");
  TH1D *h2 = (TH1D*)h1->Clone();
  TH1D *h3 = (TH1D*)h1->Clone();

  h2->GetXaxis()->SetRangeUser(500,1500);
  Int_t maxbin = h2->GetMaximumBin();
  Double_t xmax = h2->GetBinCenter(maxbin);
  Double_t ymax = h2->GetMaximum();
  Double_t ybase = h2->GetBinContent(maxbin-100);
  h2->GetXaxis()->SetRangeUser(xmax-100,xmax+100);

  TF1 *fP1 = new TF1("fP1",SkewGausSlope,xmax-100,xmax+100,6);
  Double_t y2max = ymax*2.62/34.06;
  Double_t x2max = 79.6139/80.9971*xmax;
  h2->GetXaxis()->SetRangeUser(xmax-100,xmax+100);
  fP1->SetParameter(0,ymax);
  fP1->SetParameter(1,xmax);
  fP1->SetParameter(2,10);
  fP1->SetParameter(3,10);
  fP1->SetParameter(4,ybase);
  fP1->SetParameter(5,1);

  fP1->SetParLimits(0,0,10000000);
  fP1->SetParLimits(2,0,100);
  fP1->SetParLimits(3,0,100);

  /*
  TF1 *fP1 = new TF1("fP1",Special,xmax-100,xmax+100,5);
  fP1->SetParameter(0,ymax);
  fP1->SetParameter(1,xmax);
  fP1->SetParameter(2,10);
  fP1->SetParameter(3,ybase);
  fP1->SetParameter(4,1);
  fP1->SetParLimits(0,0,10000000);
  fP1->SetParLimits(2,0,100);
*/

  /*
  TF1 *fP1 = new TF1("fP1",TwoGausSlope,xmax-100,xmax+100,8);
  Double_t y2max = ymax*2.62/34.06;
  Double_t x2max = 79.6139/80.9971*xmax;
  h2->GetXaxis()->SetRangeUser(xmax-100,xmax+100);
  fP1->SetParameter(0,ymax);
  fP1->SetParameter(1,xmax);
  fP1->SetParameter(2,10);
  fP1->SetParameter(3,y2max);
  fP1->SetParameter(4,x2max);
  fP1->SetParameter(5,10);
  fP1->SetParameter(6,ybase);
  fP1->SetParameter(7,1);
  fP1->SetParLimits(0,0,10000000);
  fP1->SetParLimits(3,0,10000000);
  fP1->SetParLimits(2,0,100);
  fP1->SetParLimits(4,0,100);
*/
  h2->Fit(fP1);
  Double_t y1 = fP1->GetParameter(0);
  Double_t y1err = fP1->GetParError(0);
  //  h1->Draw();
  c1->Print("data1.pdf");

  ////////////////////
  h3->GetXaxis()->SetRangeUser(3400,3800);
  maxbin = h3->GetMaximumBin();
  xmax = h3->GetBinCenter(maxbin);
  ymax = h3->GetMaximum();
  ybase = h3->GetBinContent(maxbin-100);
  TF1 *fP2 = new TF1("fP2",SkewGausSlope,xmax-100,xmax+100,6);

  h3->GetXaxis()->SetRangeUser(xmax-100,xmax+100);
  fP2->SetParameter(0,ymax);
  fP2->SetParameter(1,xmax);
  fP2->SetParameter(2,10);
  fP2->SetParameter(3,10);
  fP2->SetParameter(4,ybase);
  fP2->SetParameter(5,1);
  h3->Fit(fP2);
  Double_t y2 = fP2->GetParameter(0);
  Double_t y2err = fP2->GetParError(0);
  Double_t ratio = y1/y2;
  Double_t ratioerr = ratio*sqrt(pow(y1err/y1,2)+pow(y2err/y2,2));

  c1->Print("data2.pdf");
  fout << y1 << " "<< y1err << " "<< y2 <<  " " << y2err << " " << ratio << " "<< ratioerr << endl;
}


