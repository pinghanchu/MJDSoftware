#include "GATAutoCal.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>

void count(Int_t StartRun,Int_t EndRun,Int_t Channel, string EName, vector<Double_t>* FWHM, vector<Double_t>* FWHMErr){

  ifstream fin(Form("./List/ezfit_%d_%d_0.txt",StartRun,EndRun));
  Int_t cha,pos,index;
  Double_t peak, width, peakerr, widtherr, fwhm, fwhmerr;
  string ename;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >>  pos >> cha >> ename >> peak>> peakerr >> width>>widtherr >> fwhm >> fwhmerr >> index;
      //cout << pos << " " << cha << " " << ename.c_str() << " " << peak << " " << fwhm << endl;
      if(cha == Channel && ename == EName){
	//cout << cha << " " << ename.c_str() << " " << peak << " " << fwhm << endl;
	FWHM->push_back(peak);
	FWHMErr->push_back(peakerr);
	FWHM->push_back(abs(fwhm));
	FWHMErr->push_back(fwhmerr);
      }
    }
  }
  if(FWHM->size()==0){
    FWHM->push_back(0);
    FWHMErr->push_back(0);
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
  if(argc != 5 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun] [endrun] [startrun2] [endrun2]" << endl;
    return 1;
  }
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  Int_t fStartRun2 = atoi(argv[3]);
  Int_t fEndRun2 = atoi(argv[4]);
  GATAutoCal ds(fStartRun,fStartRun);
  GATAutoCal ds2(fStartRun2,fStartRun2);

  string HistPathName = "./Hist/";
  string PathName = "./Plot/";
  string FileName = Form("fwhm_%d_%d.pdf",fStartRun,fEndRun);
  string TitleName;

  vector<Int_t> channel= ds.GetChannel();
  vector<Int_t> str = ds.GetString();
  vector<Int_t> detpos = ds.GetDetPosition();
  vector<Int_t> goodbad= ds.GetGoodBad();
  vector<Int_t> channel2= ds2.GetChannel();
  vector<Int_t> str2 = ds2.GetString();
  vector<Int_t> detpos2 = ds2.GetDetPosition();
  vector<Int_t> goodbad2= ds2.GetGoodBad();

  //const Int_t nchannel = channel.size();
  TGraphErrors *gr1;
  TGraphErrors *gr2;
  vector<Int_t> Px;
  vector<Int_t> PxErr;
  vector<Double_t> Py1;
  vector<Double_t> PyErr1;
  vector<Double_t> Py2;
  vector<Double_t> PyErr2;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  Double_t Pos;
  for(size_t i =0;i<channel.size();i++){
    //cout << i << " " << channel.at(i) << endl;
    Pos = 100+str.at(i)*10+detpos.at(i);
    Par.clear();
    ParErr.clear();
    count(fStartRun,fEndRun,channel.at(i),"trapE", &Par,&ParErr);
    if(channel.at(i)%2 == 0 && goodbad.at(i)==1){
      Px.push_back(Pos);
      PxErr.push_back(0);
      Py1.push_back(Par.at(1)*2614.5/Par.at(0));
      //PyErr1.push_back(ParErr.at(1)*2614.5/Par.at(0));
      PyErr1.push_back(0);
    }
  }

  for(size_t i =0;i<channel.size();i++){
    //cout << i << " " << channel.at(i) << endl;
    Par.clear();
    ParErr.clear();
    count(fStartRun,fEndRun,channel.at(i),"trapENFDBSGCal", &Par,&ParErr);
    if(channel.at(i)%2 == 0 && goodbad.at(i)==1){
      Py2.push_back(Par.at(1));
      PyErr2.push_back(0);
    }
  }

  for(size_t i =0;i<channel2.size();i++){
    //cout << i << " " << channel2.at(i) << endl;
    Par.clear();
    ParErr.clear();
    Pos = 200+str2.at(i)*10+detpos2.at(i);
    count(fStartRun2,fEndRun2,channel2.at(i),"trapE", &Par,&ParErr);
    if(channel2.at(i)%2 == 0 && goodbad2.at(i)==1){
      Px.push_back(Pos);
      PxErr.push_back(0);
      Py1.push_back(Par.at(1)*2614.5/Par.at(0));
      PyErr1.push_back(0);
    }
  }

  for(size_t i =0;i<channel2.size();i++){
    //cout << i << " " << channel2.at(i) << endl;
    Par.clear();
    ParErr.clear();
    count(fStartRun2,fEndRun2,channel2.at(i),"trapENFDBSGCal", &Par,&ParErr);
    if(channel2.at(i)%2 == 0 && goodbad2.at(i)==1){
      Py2.push_back(Par.at(1));
      PyErr2.push_back(0);
    }
  }

  //  cout << Px.size() << endl;

  
  cout << Px.size() << " " << Py1.size() << " " << Py2.size() << endl;
  gr1 = GetGraph(Py1,PyErr1);
  gr2 = GetGraph(Py2,PyErr2);

  cout << " Test..." << endl;
  TMultiGraph *mg1 = new TMultiGraph();

  TLegend *leg1 = new TLegend(0.65,0.7,0.9,0.85);
  gStyle->SetLegendBorderSize(1);  
  gr1->SetMarkerStyle(20);
  gr2->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);
  gr2->SetMarkerSize(1);
  gr1->SetMarkerColor(4);
  gr2->SetMarkerColor(2);
  mg1->Add(gr1,"eP");
  mg1->Add(gr2,"ep");
  mg1->SetTitle(";;FWHM@2614keV");
  leg1->SetFillStyle(0);
  leg1->AddEntry(gr1,"No CT correction","lp");
  leg1->AddEntry(gr2,"With CT correction","lp");
  
  vector<Double_t> ave;
  //ave.push_back(0.);
  ave.push_back(2.905);
  ds.PlotGrid(mg1,leg1,Px,ave,1.5,9., PathName, FileName);

}



