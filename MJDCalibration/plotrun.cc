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
  ifstream fin1(Form("./List/runlist/cal.DS%s.txt",fDataSet.c_str()));
  ifstream fin2(Form("./List/runlist/bk.DS%s.txt",fDataSet.c_str()));

  cout.precision(15);
  vector<Int_t> startrun1;
  vector<Int_t> endrun1;
  Int_t r1,r2,r3,r4;
  if(fin1.is_open()){
    while(!fin1.eof()){
      fin1 >> r1 >> r2 >>r3 >>r4;
      startrun1.push_back(r1);
      endrun1.push_back(r2);
    }
  }
  startrun1.pop_back();
  endrun1.pop_back();

  vector<Int_t> startrun2;
  vector<Int_t> endrun2;

  if(fin2.is_open()){
    while(!fin2.eof()){
      fin2 >> r1 >> r2;
      startrun2.push_back(r1);
      endrun2.push_back(r2);
    }
  }
  startrun2.pop_back();
  endrun2.pop_back();

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
  Double_t startrun;
  Double_t endrun;
  Double_t R1,R1err;  
  Double_t startdate = 0;
  Double_t enddate = 0;
  Double_t R2,R2err;
  Double_t starttimeMT=0;
  Double_t endtimeMT =0;
  Double_t startdateMT =0;
  Double_t enddateMT = 0;
  vector<Double_t> X;
  vector<Double_t> XErr;
  vector<Double_t> Y;
  vector<Double_t> YErr;
  
  for(size_t i = 0;i<startrun1.size();i++){
    startrun = startrun1.at(i);
    endrun = endrun1.at(i);
    cout << startrun << " " << endrun << endl;
    //R1 = (endrun+startrun)/2;
    //R1err = (endrun-startrun)/2;
    R1 = startrun;
    R1err = 0;
    GATAutoCal ds(startrun,endrun);
    startdate = ds.GetStartTime();
    enddate =  ds.GetStopTime();
    starttimeMT = ds.GetStartTimeMT();
    endtimeMT =  ds.GetStopTimeMT();
    startdateMT = (Double_t) ds.GetStartDateMT();
    enddateMT =  (Double_t) ds.GetStopDateMT();
    TDatime da1 = GetTime(startdateMT, starttimeMT);
    TDatime da2 = GetTime(enddateMT,endtimeMT);
    TDatime da3(startdate);
    TDatime da4(enddate);
    Double_t t1 = da1.Convert();
    Double_t t2 = da2.Convert();

    //long int temp = (long int) (enddate+startdate)/2;
    //R2 = (Double_t) temp;
    R2 = startdate;
    R2err = 0;
    X.push_back(R2);
    XErr.push_back(R2err);
    Y.push_back(R1);
    YErr.push_back(R1err);
    //da1.Print();
    //da2.Print();
    // da3.Print();
    //da4.Print();
    //cout << startrun << " " << endrun << " " << startdate << " " << enddate << " "<< R1 << " " << R2 << " " 
    //<< starttimeMT <<  " " << startdateMT << " " << endtimeMT << " " << enddateMT << " " 
    //<< t1 << " " << t2 << endl;
  }

  TGraphErrors *gr1 = GetGraph(X,XErr,Y,YErr);
  gr1->SetTitle(";Time;Run");
  gr1->GetYaxis()->SetTitleOffset(1.3);
  gr1->SetMarkerStyle(22);
  gr1->SetMarkerSize(1);
  gr1->SetLineColor(Color.at(0));
  gr1->SetMarkerColor(Color.at(0));
  gr1->GetXaxis()->SetTimeDisplay(1);
  gr1->GetXaxis()->SetNdivisions(-503);
  gr1->GetXaxis()->SetTimeFormat("%Y-%m-%d");
  gr1->GetXaxis()->SetTimeOffset(0,"gmt");

  X.clear();
  Y.clear();
  XErr.clear();
  YErr.clear();
  for(size_t i = 0;i<startrun2.size();i++){
    
    startrun = startrun2.at(i);
    endrun = endrun2.at(i);
    cout << startrun << " " << endrun << endl;
    //R1 = (endrun+startrun)/2;
    //R1err = (endrun-startrun)/2;
    R1 = startrun;
    R1err = 0;
    GATAutoCal ds(startrun,endrun);
    startdate = ds.GetStartTime();
    enddate = ds.GetStopTime();
    starttimeMT = ds.GetStartTimeMT();
    endtimeMT =  ds.GetStopTimeMT();
    startdateMT = (Double_t) ds.GetStartDateMT();
    enddateMT =  (Double_t) ds.GetStopDateMT();
    TDatime da1 = GetTime(startdateMT, starttimeMT);
    TDatime da2 = GetTime(enddateMT,endtimeMT);
    TDatime da3(startdate);
    TDatime da4(enddate);
    Double_t t1 = da1.Convert();
    Double_t t2 = da2.Convert();
    //long int temp = (long int) (enddate+startdate)/2;
    //R2 = (Double_t) temp;
    R2 = startdate;
    R2err = 0;
    X.push_back(R2);
    XErr.push_back(R2err);
    Y.push_back(R1);
    YErr.push_back(R1err);
    // da1.Print();
    //da2.Print();
    //da3.Print();
    //da4.Print();
    //cout << startrun << " " << endrun << " " << startdate << " " << enddate << " "<< R1 << " " << R2 << " " 
    //<< starttimeMT <<  " " << startdateMT << " " << endtimeMT << " " << enddateMT << " " 
    //<< t1 << " " << t2 << endl;
  }

  TGraphErrors *gr2 = GetGraph(X,XErr,Y,YErr);
  gr2->SetTitle(";Time;Run");
  gr2->GetYaxis()->SetTitleOffset(1.3);
  gr2->SetMarkerStyle(26);
  gr2->SetMarkerSize(1);
  gr2->SetLineColor(Color.at(1));
  gr2->SetMarkerColor(Color.at(1));
  gr2->GetXaxis()->SetTimeDisplay(1);
  gr2->GetXaxis()->SetNdivisions(-503);
  gr2->GetXaxis()->SetTimeFormat("%Y-%m-%d");
  gr2->GetXaxis()->SetTimeOffset(0,"gmt");


  TMultiGraph *mg1 = new TMultiGraph();
  mg1->SetTitle(";Time;Run");
  TLegend *leg1 = new TLegend(0.2,0.65,0.4,0.85);
  gStyle->SetLegendBorderSize(0);
  mg1->Add(gr1,"ep");
  mg1->Add(gr2,"ep");
  //mg1->GetXaxis()->SetTimeDisplay(1);
  //mg1->GetXaxis()->SetNdivisions(-503);
  //mg1->GetXaxis()->SetTimeFormat("%Y-%m-%d");
  //mg1->GetXaxis()->SetTimeOffset(0,"gmt");

  leg1->AddEntry(gr1,"Calibration Run");
  leg1->AddEntry(gr2,"Data Run");
  mg1->Draw("AP");
  mg1->GetXaxis()->SetTimeDisplay(1);
  mg1->GetXaxis()->SetNdivisions(-503);
  mg1->GetXaxis()->SetTimeFormat("%Y-%m-%d");
  mg1->GetXaxis()->SetTimeOffset(0,"gmt");

  mg1->Draw("AP");
  //gr1->Draw("AP");
  leg1->Draw();
  c1->Print(Form("./Plot/run_DS%s.pdf",fDataSet.c_str()));

}


