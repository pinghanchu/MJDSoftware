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
#include "TDatime.h"
#include <string>
#include <stdlib.h>

TDatime GetTime(Int_t date, Double_t time){
  Int_t starthour = (Int_t) time;
  Int_t startmin = (Int_t)((time-(Double_t)starthour)*60);
  Int_t startsec = (Int_t)((time-(Double_t)starthour)*60-(Double_t)startmin)*60;
  Int_t startyear = (Int_t)((Double_t)date/10000);
  Int_t startmonth = (Int_t)(((Double_t)(date-startyear*10000))/100);
  Int_t startdate = date-startyear*10000-startmonth*100;
  TDatime da(startyear,startmonth,startdate,starthour,startmin,startsec);
  return da;
}


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
    cout << "Usage: " << argv[0] << " [data set]" << endl;
    return 1;
  }
  string fDataSet = argv[1];
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  TCanvas *c1 = new TCanvas("c1");
  ifstream fin1("./List/cal/ezfit_HG.txt");

  cout.precision(15);
  vector<Int_t> fStartRun;
  vector<Int_t> fEndRun;
  vector<Double_t> fStartTime;
  vector<Double_t> fEndTime;
  vector<Double_t> fPeak;
  vector<Double_t> fPeakErr;
  vector<Double_t> fFWHM;
  vector<Double_t> fFWHMErr;
  vector<Double_t> fTimeErr;
  Int_t run1,run2,index;
  Double_t starttime,endtime,x1,x2,x3,x4;
  string dataname,ename;
  if(fin1.is_open()){
    while(!fin1.eof()){
      fin1 >> dataname >> run1>>run2 >> starttime >> endtime >> ename >> index >> x1 >> x2 >>x3 >>x4;
      fStartRun.push_back(run1);
      fEndRun.push_back(run2);
      fStartTime.push_back(starttime);
      fEndTime.push_back(endtime);
      fPeak.push_back(abs(x1));
      fPeakErr.push_back(abs(x2));
      fFWHM.push_back(abs(x3));
      fFWHMErr.push_back(abs(x4));
      fTimeErr.push_back(0);
      cout << run1 << endl;
    }
  }
  fStartRun.pop_back();
  fEndRun.pop_back();
  fStartTime.pop_back();
  fEndTime.pop_back();
  fPeak.pop_back();
  fPeakErr.pop_back();
  fFWHM.pop_back();
  fFWHMErr.pop_back();
  fTimeErr.pop_back();

  vector<Int_t> Color;
  Color.push_back(2);
  Color.push_back(3);
  Color.push_back(4);
  Color.push_back(6);
  Color.push_back(7);
  Color.push_back(8);
  Color.push_back(9);
  Color.push_back(12);
  Color.push_back(28);
  Color.push_back(32);
  Color.push_back(37);
  Color.push_back(38);
  Color.push_back(41);
  Color.push_back(44);
  Color.push_back(46);
  Color.push_back(49);

  vector<Int_t> Style;
  Style.push_back(20);
  Style.push_back(24);
  Style.push_back(21);
  Style.push_back(25);
  Style.push_back(22);
  Style.push_back(26);
  Style.push_back(23);
  Style.push_back(32);
  Style.push_back(33);
  Style.push_back(27);
  Style.push_back(34);
  Style.push_back(28);
  Style.push_back(29);
  Style.push_back(30);
  TGraphErrors *gr1 = GetGraph(fStartTime,fTimeErr,fPeak,fPeakErr);
  gr1->SetTitle(";Time;Peak(keV)");
  gr1->GetYaxis()->SetTitleOffset(1.3);
  gr1->SetMarkerStyle(22);
  gr1->SetMarkerSize(1);
  gr1->SetLineColor(Color.at(0));
  gr1->SetMarkerColor(Color.at(0));
  gr1->GetXaxis()->SetTimeDisplay(1);
  gr1->GetXaxis()->SetNdivisions(-503);
  gr1->GetXaxis()->SetTimeFormat("%Y-%m-%d");
  gr1->GetXaxis()->SetTimeOffset(0,"gmt");
  gr1->Draw("AP");
  //gr1->Draw("AP");
  c1->Print("./Plot/peaktime.pdf");

  TH1D *h1 = new TH1D("h1","", 100,2614.5-1.5,2614.5+1.5);
  for(size_t ii = 0;ii<fPeak.size();ii++){
    h1->Fill(fPeak.at(ii));
  }
  h1->SetTitle(";Peak(keV);Subsets");
  h1->Fit("gaus");
  c1->Print("./Plot/peakhist.pdf");

  TGraphErrors *gr2 = GetGraph(fStartTime,fTimeErr,fFWHM,fFWHMErr);
  gr2->SetTitle(";Time;FWHM(keV)");
  gr2->GetYaxis()->SetTitleOffset(1.3);
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerSize(1);
  gr2->SetLineColor(Color.at(0));
  gr2->SetMarkerColor(Color.at(0));
  gr2->GetXaxis()->SetTimeDisplay(1);
  gr2->GetXaxis()->SetNdivisions(-503);
  gr2->GetXaxis()->SetTimeFormat("%Y-%m-%d");
  gr2->GetXaxis()->SetTimeOffset(0,"gmt");
  gr2->Draw("AP");
  c1->Print("./Plot/fwhmtime.pdf");

  TH1D *h2 = new TH1D("h2","", 100,2.5,3.5);
  for(size_t ii = 0;ii<fFWHM.size();ii++){
    h2->Fill(fFWHM.at(ii));
  }
  h2->SetTitle(";FWHM(keV);Subsets");
  h2->Fit("gaus");
  c1->Print("./Plot/fwhmhist.pdf");
}


