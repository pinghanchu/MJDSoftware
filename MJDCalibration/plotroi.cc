#include "GATAutoCal.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>

void GetROIPar(Int_t fStartRun, Int_t fEndRun, string fChannel, string fEName, Double_t fRef, vector<Double_t>* Par, vector<Double_t>* ParErr){

  ifstream fin(Form("./List/cal/roi_%d_%d.txt",fStartRun,fEndRun));

  string cha, pos;
  Double_t r1,r2,r3,r4,r0;
  string ename;
  Int_t startrun,endrun;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> startrun >> endrun >> ename >> r0 >> r1 >> r2 >>r3 >>r4;
      Int_t tempchan = atoi(cha.c_str());
      if(cha == fChannel && ename == fEName && r0 == fRef){
        Par->push_back(r1);
        ParErr->push_back(r2);
        Par->push_back(r3*2.355);
        ParErr->push_back(r4*2.355);
      }
    }
  }
}



void GetROIPar(Int_t fStartRun, Int_t fEndRun, string fPos, Int_t fGain, string fEName, Double_t fRef, vector<Double_t>* Par, vector<Double_t>* ParErr){

  ifstream fin(Form("./List/cal/roi_%d_%d.txt",fStartRun,fEndRun));

  string cha;
  string pos;
  Double_t r1,r2,r3,r4,r0;
  string ename;        
  Int_t startrun,endrun;
  Int_t gain;
  if(fin.is_open()){   
    while(!fin.eof()){ 
      fin >> pos >> cha >> startrun >> endrun >> ename >> r0 >> r1 >> r2 >>r3 >>r4;
      Int_t tempchan = atoi(cha.c_str());
      gain = tempchan%2;
      if(pos == fPos && gain == fGain && ename == fEName && r0 == fRef){
        Par->push_back(r1);                                  
        ParErr->push_back(r2);                               
        Par->push_back(r3*2.355);
        ParErr->push_back(r4*2.355);
      }                                                         
    }                                                           
  } 
}

TGraphErrors* GetGraph(vector<Double_t> Py, vector<Double_t> PyErr){
  const Int_t nPy = Py.size();

  Double_t A[nPy];
  Double_t AErr[nPy];
  Double_t B[nPy];
  Double_t BErr[nPy];

  for(Int_t j=0;j<nPy;j++){
    A[j]=(Double_t)j;
    AErr[j] =0;
    B[j]=Py.at(j);
    BErr[j] =PyErr.at(j);
  }
  TGraphErrors* fGraph = new TGraphErrors(nPy,A,B,AErr,BErr);
  fGraph->SetMarkerStyle(8);
  fGraph->SetMarkerSize(2);
  return fGraph;
}


int main(int argc, char** argv)
{
  if(argc != 4 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun] [endrun] [energy name]" << endl;
    return 1;
  }
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  string fEName = argv[3];
  ofstream froi(Form("./List/cal/roi1_%d_%d.txt",fStartRun,fEndRun),ios::app);
  froi.precision(4);
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


  vector<Int_t> LineStyle;
  LineStyle.push_back(1);
  LineStyle.push_back(2);
  LineStyle.push_back(3);
  LineStyle.push_back(4);
  LineStyle.push_back(5);
  LineStyle.push_back(6);
  LineStyle.push_back(7);
  LineStyle.push_back(8);
  LineStyle.push_back(9);
  LineStyle.push_back(10);

  GATAutoCal ds(fStartRun,fStartRun+1);
  string HistPathName = "./Hist/";
  string PathName = "./";
  string FileName = Form("roichan_%d_%d.pdf",fStartRun,fEndRun);
  string TitleName;

  vector<Int_t> fChannel= ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fString = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad= ds.GetGoodBad();
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  string Pos;
  TCanvas *c1 = new TCanvas("c1");
  //const Int_t nChannel = fChannel.size();
  //TF1 *fun[nChannel];
  //vector<Int_t> index;

  vector<Double_t> fRef;
  fRef.push_back(60);
  fRef.push_back(2039);
  fRef.push_back(2614.5);
  Int_t fGain;
  TGraphErrors *gr[3];
  vector<Double_t> fwhm;
  vector<Double_t> fwhmerr;
  vector<Int_t> Px;
  vector<Int_t> PxErr;

  for(size_t ir=0;ir<fRef.size();ir++){    
    fwhm.clear();
    fwhmerr.clear();
    //Px.clear();
    //PxErr.clear();
    for(size_t i =0;i<fChannel.size();i++){
      Par.clear();
      ParErr.clear();
      Pos = Form("%d%d%d",fCryo.at(i),fString.at(i),fDetpos.at(i));
      fGain = fChannel.at(i)%2;
      if(fChannel.at(i)%2 == 0 && fGoodBad.at(i)==1){
	GetROIPar(fStartRun,fEndRun,Pos,fGain,fEName, fRef.at(ir), &Par, &ParErr);
	Int_t pars = Par.size();
	if(Par.size()>0){
	  froi << Pos.c_str() << " "<< fChannel.at(i) << " " <<fStartRun << " " << fEndRun << " "
	       << fEName.c_str() << " " << fRef.at(ir) << " " 
	       << Par.at(pars-2)<< "\\pm"<<ParErr.at(pars-2) <<" "
	       << Par.at(pars-1)<< "\\pm"<<ParErr.at(pars-1) <<endl;	
	  fwhm.push_back(Par.at(pars-1));
	  fwhmerr.push_back(ParErr.at(pars-1));
	  if(ir==0){
	    Px.push_back(i);
	    PxErr.push_back(0);
	  }
	}
      }
    }
    gr[ir] = GetGraph(fwhm,fwhmerr);
  }
  TMultiGraph *mg1 = new TMultiGraph();
  TLegend *leg1 = new TLegend(0.2,0.7,0.85,0.85);
  leg1->SetNColumns(3);
  for(size_t ir=0;ir<fRef.size();ir++){
    gr[ir]->SetMarkerStyle(Style.at(ir*2));
    gr[ir]->SetMarkerColor(Color.at(ir));
    gr[ir]->SetLineColor(Color.at(ir));
    mg1->Add(gr[ir],"ep");
    leg1->AddEntry(gr[ir],Form("@%.01f(keV)",fRef.at(ir)),"lp");
  }
  vector<Double_t> ave;
  for(size_t ir=0;ir<fRef.size();ir++){
    Par.clear();
    ParErr.clear();
    GetROIPar(fStartRun,fEndRun,"HG",fEName, fRef.at(ir), &Par, &ParErr);
    ave.push_back(Par.at(Par.size()-1));
  }
  ds.PlotGrid(mg1,leg1,Px,ave,0.,5.5, PathName, FileName);

}


