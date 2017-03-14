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
#include <string>
#include <stdlib.h>

Double_t Flat(Double_t *v, Double_t *par){
  Double_t fitval = par[0];
  return fitval;
}


void GetCal(Int_t fStartRun, Int_t fEndRun, Int_t fChannel, string fEName, vector<Double_t>* Par){
  ifstream fin(Form("./List/cal/cov_%d_%d.txt",fStartRun,fEndRun));

  Int_t cha,pos;
  string R[4];
  Double_t r[4];
  Int_t run1,run2;
  string ename;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> run1 >> run2 >> ename >> R[0] >> R[1] >> R[2] >>R[3];
      for(Int_t  i = 0;i<4;i++){
        if(R[i]=="nan" || R[i] =="-nan" || R[i] == "inf" || R[i] == "-inf"){
          r[i] = 0;
        }else{
          r[i] = atof(R[i].c_str());
        }
      }
      if(cha == fChannel && ename == fEName.c_str()){
	Par->push_back( r[0] );  // offset
	Par->push_back( r[1] );
	Par->push_back( r[2] ); // slope
	Par->push_back( r[3] );
      }
    }
  }
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
  if(argc != 4 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [energy name]" << endl;
    return 1;
  }
  int fStartRun = atoi(argv[1]);
  int fEndRun = atoi(argv[2]);
  const char* EnergyName = argv[3];
  string fEName(EnergyName);
  //gStyle->SetOptStat(1100);
  //gStyle->SetOptFit(1111);
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  //TCanvas *c1 = new TCanvas("c1");
  ofstream fmiss(Form("./List/miss.cal_%d_%d.txt",fStartRun,fEndRun));
  ofstream falign(Form("./List/align.cal_%d_%d.txt",fStartRun,fEndRun));
  ofstream fbad(Form("./List/bad.cal_%d_%d.txt",fStartRun,fEndRun));
  ifstream fin1("./List/runlist/cal.all.txt");
  vector<Int_t> startrun;
  vector<Int_t> endrun;
  Int_t r1,r2,r3,r4;
  if(fin1.is_open()){
    while(!fin1.eof()){
      fin1 >> r1 >> r2 >> r3 >> r4 ;
      if(r1>= fStartRun && r2<= fEndRun){
	cout << r1 << " " << r2 << endl;
	startrun.push_back(r1);
	endrun.push_back(r2);
      }
    }
  }


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


  GATAutoCal ds(fStartRun,fStartRun);
  ds.SetEnergyName(EnergyName);
  //  ds.SetParameters();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos= ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  const Int_t channels = fChannel.size();
  string Pos;
  string TitleName;
  string PlotName;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> offset;
  vector<Double_t> offseterr;
  vector<Double_t> slope;
  vector<Double_t> slopeerr;
  vector<Double_t> run;
  vector<Double_t> runerr;
  vector<Int_t> GoodBadIndex;

  TGraphErrors *fOffset[channels];
  TGraphErrors *fSlope[channels];
  for(Int_t i=0;i<channels;i++){
    Pos = Form("%d%d%d",fCryo.at(i),fStr.at(i),fDetpos.at(i));
    cout << Pos.c_str() << " " << fChannel.at(i) << endl;
    if(fGoodBad.at(i)>0){
      //      c1->Update();
      run.clear();
      runerr.clear();
      offset.clear();
      offseterr.clear();
      slope.clear();
      slopeerr.clear();
      for(size_t j = 0 ; j<startrun.size();j++){
	Par.clear();        
	GetCal(startrun.at(j),endrun.at(j),fChannel.at(i),fEName,&Par);
	Int_t pars = Par.size();
	Double_t r1=(Double_t) (startrun.at(j)+endrun.at(j))/2.;
	Double_t r2=(Double_t) (-startrun.at(j)+endrun.at(j))/2.;

	run.push_back(r1);
	runerr.push_back(r2);

	if(Par.size()>0){	  
	  offset.push_back(Par.at(pars-4));
	  offseterr.push_back(Par.at(pars-3));
	  slope.push_back(Par.at(pars-2));
	  slopeerr.push_back(Par.at(pars-1));
	}
      }

      for(size_t j = 0;j<slope.size();j++){
	ofstream fnew(Form("./List/cal/cov1_%d_%d.txt",startrun.at(j),endrun.at(j)),ios::app);
	fnew << Pos.c_str() << " " << fChannel.at(i) << " " << startrun.at(j) << " " << endrun.at(j) << " " << fEName.c_str() << " " 
	     << offset.at(j) << " "<< offseterr.at(j) << " "<< slope.at(j) << " " << slopeerr.at(j) << endl;
	fnew.close();
      }
    }
  }


}


