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

void GetPar(string fInputName, vector<Double_t>* X, vector<Double_t>* XUnc, vector<Double_t>* Y, vector<Double_t>* YUnc){
  ifstream fin(Form("%s",fInputName.c_str()));
  Double_t r1,r2,r3,r4,r5,r6,r7,r8;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> r1 >> r2 >> r3 >> r4 >> r5 >> r6 >> r7 >> r8;
      X->push_back(r1+5);
      XUnc->push_back(5);
      Y->push_back(r7);
      YUnc->push_back(sqrt(r5)/2);
    }
  }
  X->pop_back();
  XUnc->pop_back();
  Y->pop_back();
  YUnc->pop_back();
}

int main(int argc, char** argv)
{
  if(argc != 3) {
    cout << "Usage: " << argv[0] << "[run list] [run list2]" << endl;
    return 1;
  }
  //Int_t fDataSet = atoi(argv[1]);
  string fInputName = argv[1];
  string fInputName2 = argv[2];
  gROOT->ProcessLine(".x $GATDIR/MJDCalibration/MJDTalkPlotStyle.C");
  //cout << fInputName.c_str() << endl;
  vector<Double_t> X1;
  vector<Double_t> X1Unc;
  vector<Double_t> Y1;
  vector<Double_t> Y1Unc;
  vector<Double_t> X2;
  vector<Double_t> X2Unc;
  vector<Double_t> Y2;
  vector<Double_t> Y2Unc;
  GetPar(fInputName,&X1,&X1Unc,&Y1,&Y1Unc);
  GetPar(fInputName2,&X2,&X2Unc,&Y2,&Y2Unc);

  TCanvas *c1 = new TCanvas("c1");
  
  //TGraphErrors *gr1 = GetGraph(StartTime,StartTimeErr,DeltaMu,DeltaErr);

  TGraphErrors *gr1 = GetGraph(X1,X1Unc,Y1,Y1Unc);
  gr1->SetTitle(";Energy (keV); Nonlinearity Shift (keV)");
  gr1->GetYaxis()->SetTitleOffset(1.3);
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerSize(1);
  gr1->SetLineColor(4);
  gr1->SetMarkerColor(4);
  gr1->SetMaximum(0.6);
  gr1->SetMinimum(-0.5);
  TGraphErrors *gr2 = GetGraph(X2,X2Unc,Y2,Y2Unc);
  gr2->SetTitle(";Energy (keV); Nonlinearity Shift (keV)");
  gr2->GetYaxis()->SetTitleOffset(1.3);
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerSize(1);
  gr2->SetLineColor(2);
  gr2->SetMarkerColor(2);
  gr2->SetMaximum(0.6);
  gr2->SetMinimum(-0.5);

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Draw("ap");
  mg->SetTitle(";Energy (keV); Nonlinearity Shift (keV)");

  TLegend* legend = new TLegend(0.6,0.7,0.85,0.85);
  legend->AddEntry(gr1,"Without NL corr.","lp");
  legend->AddEntry(gr2,"With NL corr.","lp");
  legend->Draw();
  gStyle->SetLegendTextSize(0.04);
  legend->SetTextFont(132);
  legend->Draw();

  c1->Print("delta_NL.pdf");
  c1->Print("delta_NL.C");

}


