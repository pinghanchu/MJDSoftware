// 2016.2.18
// written by Pinghan Chu
// Following the logic in GATAutoCal.cc
// 
#include "MJDSkim.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include <string>
#include <stdlib.h>
/*
Double_t SkewGaus1(Double_t *v, Double_t *par){
  Double_t arg = 0;
  if (v[0]>=(par[1])){
    arg = (v[0] - (par[1]))/par[2];
  }
  if (v[0]<(par[1])){
    arg = (v[0] - (par[1]))/par[3];
  }

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}
*/
/*
Double_t Gaus1(Double_t *v, Double_t *par){
  Double_t arg = 0;
  arg = (v[0] - (par[1]))/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}

Double_t TwoGaus1(Double_t *v, Double_t *par){
  return Gaus1(v,par)+Gaus1(v,&par[3]);
}


Double_t IsNan(string Input){
  Double_t output = 0;
  if(Input == "nan" || Input == "-nan" || Input == "inf" || Input == "-inf"){
    output = 0;
  }else{
    output = atof(Input.c_str());
  }
  return output;
}
*/

TGraphErrors* GetGraph(vector<Double_t> Px,vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr){
  const Int_t nPy = Py.size();

  Double_t A[nPy];
  Double_t AErr[nPy];
  Double_t B[nPy];
  Double_t BErr[nPy];

  for(Int_t j=0;j<nPy;j++){
    A[j]=Px.at(j);
    AErr[j] = PxErr.at(j);
    B[j]=Py.at(j);
    BErr[j] =PyErr.at(j);
  }
  TGraphErrors* fGraph = new TGraphErrors(nPy,A,B,AErr,BErr);
  fGraph->SetFillStyle(0);
  fGraph->SetFillColor(0);
  return fGraph;
}


int main(int argc, char** argv)
{
  if(argc != 2) {
    cout << "Usage: " << argv[0] << "[run list]" << endl;
    return 1;
  }
  //Int_t fDataSet = atoi(argv[1]);
  string fInputName = argv[1];

  gROOT->ProcessLine(".x $GATDIR/MJDCalibration/MJDTalkPlotStyle.C");
  //cout << fInputName.c_str() << endl;
  ifstream fin(Form("%s",fInputName.c_str()));
  vector<Double_t> Run;
  vector<Double_t> RunErr;
  vector<Double_t> DeltaMu;
  vector<Double_t> DeltaErr;

  Double_t r1,r2,r3,r4,r5,r6,r7,r8;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> r1 >> r2 >> r3 >> r4 >> r5 >> r6 >> r7 >> r8;
      Run.push_back(r1+5);
      RunErr.push_back(5);
      DeltaMu.push_back(r7);
      DeltaErr.push_back(sqrt(r5)/2);
    }
  }
  Run.pop_back();
  RunErr.pop_back();
  DeltaMu.pop_back();
  DeltaErr.pop_back();

  TCanvas *c1 = new TCanvas("c1");
  
  //TGraphErrors *gr1 = GetGraph(StartTime,StartTimeErr,DeltaMu,DeltaErr);

  TGraphErrors *gr3 = GetGraph(Run,RunErr,DeltaMu,DeltaErr);
  gr3->SetTitle(";Energy (keV); Nonlinearity Shift (keV)");
  gr3->GetYaxis()->SetTitleOffset(1.3);
  gr3->SetMarkerStyle(8);
  gr3->SetMarkerSize(1);
  gr3->SetLineColor(4);
  gr3->SetMarkerColor(4);
  gr3->SetMaximum(0.5);
  gr3->SetMinimum(-0.5);
  //gr3->GetXaxis()->SetTimeDisplay(1);
  //gr3->GetXaxis()->SetNdivisions(-503);
  //gr3->GetXaxis()->SetTimeFormat("%Y-%m-%d");
  //gr3->GetXaxis()->SetTimeOffset(0,"gmt");
  gr3->Draw("ap");
  c1->Print("delta_NL.pdf");
  gr3->Draw("ap");
  c1->Print("delta_NL.C");

}


