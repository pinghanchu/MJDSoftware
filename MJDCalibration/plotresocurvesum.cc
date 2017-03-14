#include "GATAutoCal.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>

Double_t Reso1(Double_t *v, Double_t *par){
  Double_t fitval = 2.355*TMath::Sqrt(par[0]*par[0]+par[1]*par[1]*v[0]+par[2]*par[2]*v[0]*v[0]);
  return fitval;
}


void GetResoPar(Int_t StartRun, Int_t EndRun, string Channel, string EName, vector<Double_t>* Par, vector<Double_t>* ParErr){

  ifstream fin(Form("./List/cal/reso_%d_%d.txt",StartRun,EndRun));

  string cha, pos;
  Double_t r1,r2,r3,r4,r5,r6;
  string ename;        
  Int_t startrun,endrun;
  if(fin.is_open()){   
    while(!fin.eof()){ 
      fin >> pos >> cha >> startrun >> endrun >> ename >> r1 >> r2 >>r3 >>r4 >> r5 >> r6;
      if(cha == Channel && ename == EName){                     
        Par->push_back(r1);                                  
        ParErr->push_back(r2);                               
        Par->push_back(r3);
        ParErr->push_back(r4);
        Par->push_back(r5);
        ParErr->push_back(r6);
      }                                                         
    }                                                           
  } 
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
  ofstream freso(Form("./List/cal/reso1_%d_%d.txt",fStartRun,fEndRun),ios::app);
  freso.precision(3);

  GATAutoCal ds(fStartRun,fStartRun+1);
  string HistPathName = "./Hist/";
  string PathName = "./Plot/";
  string FileName = Form("./reso_%d_%d.pdf",fStartRun,fEndRun);

  vector<string> DataName;
  DataName.push_back("00"); // Natural HG
  DataName.push_back("01"); // Natural LG
  DataName.push_back("10"); // Enriched HG
  DataName.push_back("11"); // Enriched LG
  DataName.push_back("HG");
  DataName.push_back("LG");
  vector<string> TitleName;
  TitleName.push_back("Natural HG");
  TitleName.push_back("Natural LG");
  TitleName.push_back("Enriched HG");
  TitleName.push_back("Enrichd LG");
  TitleName.push_back("High Gain");
  TitleName.push_back("Low Gain");

  vector<Double_t> Par;
  vector<Double_t> ParErr;
  string Pos;
  TCanvas *c1 = new TCanvas("c1");
  const Int_t nDataName = DataName.size();
  TF1 *fun[nDataName];

  for(size_t i =0;i<DataName.size();i++){
    Par.clear();
    ParErr.clear();
    GetResoPar(fStartRun,fEndRun,DataName.at(i), fEName, &Par, &ParErr);
    Int_t pars = Par.size();
    if(Par.size()>0){
      fun[i] = new TF1("f1",Reso1,0,3000,3);
      fun[i]->SetTitle(";Energy(keV); FWHM(keV)");
      fun[i]->SetLineStyle(1);
      fun[i]->SetLineWidth(1);
      fun[i]->SetParameter(0,Par.at(pars-3));
      fun[i]->SetParameter(1,Par.at(pars-2));
      fun[i]->SetParameter(2,Par.at(pars-1));
      fun[i]->Draw();
      c1->Print(Form("reso_%s_%d_%d_%s.pdf",fEName.c_str(), fStartRun,fEndRun, DataName.at(i).c_str()));
      freso << "000 " << DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun<<" " 
	    << Par.at(pars-3) << "\\pm" << ParErr.at(pars-3) << " " 
	    << Par.at(pars-2) << "\\pm" << ParErr.at(pars-2) << " "
	    << Par.at(pars-1) << "\\pm" << ParErr.at(pars-1) << endl;
    }
  }

  TLegend *leg1 = new TLegend(0.2,0.35,0.4,0.85);
  gStyle->SetLegendBorderSize(0);
  leg1->SetFillStyle(0);


  Int_t j = 0;
  for(size_t i=0;i<DataName.size();i++){
    fun[i]->SetLineColor(Color.at(i));
    leg1->AddEntry(fun[i],Form("%s",TitleName.at(i).c_str()));
  }
  for(size_t i=0;i<DataName.size();i++){
    fun[i]->Draw("same");
  }
  leg1->Draw();

  c1->Print(Form("reso_%s_%d_%d.pdf",fEName.c_str(),fStartRun,fEndRun));

}


