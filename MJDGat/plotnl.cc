// 2016.2.18
// This macro uploads the calibration parameters.
// 
#include "GATAutoCal.hh"
//#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

void GetPar(string Input, vector<Double_t>* Par, vector<Double_t>*  ParErr){
  ifstream fin(Form("%s",Input.c_str()));

  Int_t index;
  string name;
  Double_t par, parerr;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> index >> name >> par >> parerr;
      Par->push_back(par);
      ParErr->push_back(parerr);
    }
  }
  Par->pop_back();
  ParErr->pop_back();
}

void GetCal(string Input, Int_t Channel, vector<Double_t>* Par, vector<Double_t>*  ParErr){
  ifstream fin(Form("%s",Input.c_str()));

  Int_t run1,run2,pos,chan;
  string name;
  Double_t offset,offseterr,slope,slopeerr,roi,roierr;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run1 >> run2 >> name >> pos >> chan >> offset >> offseterr >> slope >> slopeerr >> roi >> roierr;
      if(chan == Channel){
	Par->push_back(offset);
	Par->push_back(slope);
	ParErr->push_back(offseterr);
	ParErr->push_back(slopeerr);
      }
    }
  }
  //Par->pop_back();
  //ParErr->pop_back();
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
  fGraph->SetMarkerStyle(26);
  fGraph->SetMarkerSize(2);
  return fGraph;
}




int main(int argc, char** argv)
{
  if(argc != 7 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [position] [channel] [energy 1] [energy 2]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  Int_t fPos = atoi(argv[3]);
  Int_t fChannel = atoi(argv[4]);
  string fEnr1 = argv[5];
  string fEnr2 = argv[6];

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");

  string FileName1 = Form("parameters_%d_%d_%s_%d_%d.txt",fStartRun,fEndRun,fEnr1.c_str(),fPos,fChannel);
  string FileName2 = Form("parameters_%d_%d_%s_%d_%d.txt",fStartRun,fEndRun,fEnr2.c_str(),fPos,fChannel);
  string FileName3 = Form("calibration_%d_%d_%s_%d.txt",fStartRun,fEndRun,fEnr1.c_str(),fPos);
  string FileName4 = Form("calibration_%d_%d_%s_%d.txt",fStartRun,fEndRun,fEnr2.c_str(),fPos);

  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> Px;
  vector<Double_t> PxErr;
  vector<Double_t> Py;
  vector<Double_t> PyErr;
  Double_t slope;
  Double_t offset;
  Par.clear();
  ParErr.clear();
  Px.clear();
  Py.clear();
  PxErr.clear();
  PyErr.clear();


  Px.push_back(238.632);
  Px.push_back(240.986);
  Px.push_back(277.371);
  Px.push_back(300.087);
  Px.push_back(583.191);
  Px.push_back(727.330);
  Px.push_back(785.37);
  Px.push_back(860.557);
  Px.push_back(2614.533);

  Int_t calpeaks = Px.size();
  
  GetCal(FileName3, fChannel, &Par, &ParErr);
  Int_t pars = Par.size();
  slope  = Par.at(pars-1);
  offset = Par.at(pars-2);
  Par.clear();
  ParErr.clear();
  GetPar(FileName1, &Par,&ParErr);
  for(Int_t i = calpeaks;i<calpeaks*2;i++){
    Int_t ii = i-calpeaks;
    PxErr.push_back(0);
    Py.push_back(Par.at(i)*slope+offset-Px.at(ii));
    PyErr.push_back(ParErr.at(i)*slope);
  }

  TGraphErrors* gr1 = GetGraph(Px,PxErr,Py,PyErr);
  /////////////////////
  Par.clear();
  ParErr.clear();
  Py.clear();
  PxErr.clear();
  PyErr.clear();
  GetCal(FileName4, fChannel, &Par, &ParErr);
  pars = Par.size();
  slope  = Par.at(pars-1);
  offset = Par.at(pars-2);
  Par.clear();
  ParErr.clear();
  GetPar(FileName2, &Par,&ParErr);
  for(Int_t i = calpeaks;i<calpeaks*2;i++){
    Int_t ii = i-calpeaks;
    PxErr.push_back(0);
    Py.push_back(Par.at(i)*slope+offset-Px.at(ii));
    PyErr.push_back(ParErr.at(i)*slope);
  }
  TGraphErrors* gr2 = GetGraph(Px,PxErr,Py,PyErr);
  gr1->SetMarkerColor(2);
  gr2->SetMarkerColor(4);
  gr1->SetLineColor(2);
  gr2->SetLineColor(4);

  TMultiGraph *mDeltaPeak = new TMultiGraph();
  TLegend *legDeltaPeak = new TLegend(0.65,0.65,0.85,0.85);
  legDeltaPeak->SetFillStyle(0);
  //legDeltaPeak->SetNColumns(6);
  legDeltaPeak->SetBorderSize(0);

  mDeltaPeak->Add(gr1,"ep");
  legDeltaPeak->AddEntry(gr1,"W/O NL Corr");
  mDeltaPeak->Add(gr2,"ep");
  legDeltaPeak->AddEntry(gr2,"W NL Corr");
  TCanvas *c1 = new TCanvas("c1");
  
  mDeltaPeak->SetTitle(Form(";Energy (keV);#Delta E (keV)"));
  //mDeltaPeak->SetMaximum(7);
  //mDeltaPeak->SetMinimum(-3);
  mDeltaPeak->Draw("AP");
  legDeltaPeak->Draw();
  c1->Print(Form("deltapeak_%d_%d_%s_%s_%d_%d.pdf",fStartRun,fEndRun,fEnr1.c_str(),fEnr2.c_str(),fPos,fChannel));

}


