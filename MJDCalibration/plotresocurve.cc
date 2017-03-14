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


void GetResoPar(Int_t StartRun, Int_t EndRun, Int_t FChannel, string EName, vector<Double_t>* Par, vector<Double_t>* ParErr){

  ifstream fin(Form("./List/cal/reso_%d_%d.txt",StartRun,EndRun));

  string cha, pos;
  Double_t r1,r2,r3,r4,r5,r6;
  string ename;        
  Int_t startrun,endrun;
  if(fin.is_open()){   
    while(!fin.eof()){ 
      fin >> pos >> cha >> startrun >> endrun >> ename >> r1 >> r2 >>r3 >>r4 >> r5 >> r6;
      Int_t tempchan = atoi(cha.c_str());
      if(tempchan == FChannel && ename == EName){                     
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
  ofstream freso(Form("./List/cal/reso1_%d_%d.txt",fStartRun,fEndRun),ios::app);
  freso.precision(3);
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
  string PathName = "./Plot/";
  string FileName = Form("reso_%d_%d.pdf",fStartRun,fEndRun);
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
  const Int_t nChannel = fChannel.size();
  TF1 *fun[nChannel];
  vector<Int_t> index;


  for(size_t i =0;i<fChannel.size();i++){
    Par.clear();
    ParErr.clear();
    Pos = Form("%d%d%d",fCryo.at(i),fString.at(i),fDetpos.at(i));
    if(fChannel.at(i)%2 >= 0 && fGoodBad.at(i)==1){
      GetResoPar(fStartRun,fEndRun,fChannel.at(i),fEName, &Par, &ParErr);
      Int_t pars = Par.size();
      if(Par.size()>0){
	index.push_back(i);
	fun[i] = new TF1("f1",Reso1,0,3000,3);
	fun[i]->SetTitle(";Energy(keV); FWHM(keV)");
	fun[i]->SetLineStyle(1);
	fun[i]->SetLineWidth(1);
	fun[i]->SetParameter(0,Par.at(pars-3));
	fun[i]->SetParameter(1,Par.at(pars-2));
	fun[i]->SetParameter(2,Par.at(pars-1));
	fun[i]->Draw();
	c1->Print(Form("reso_%s_%d_%d_%s_%d.pdf",fEName.c_str(),fStartRun,fEndRun,Pos.c_str(),fChannel.at(i)));
	
	freso << Pos.c_str() << " "<< fChannel.at(i) << " " <<fStartRun << " " << fEndRun << " "<< fEName.c_str() << " " 
	      << Par.at(pars-3)<< "\\pm"<<ParErr.at(pars-3) <<" "
	      << Par.at(pars-2)<< "\\pm"<<ParErr.at(pars-2) <<" "
	      << Par.at(pars-1)<< "\\pm"<<ParErr.at(pars-1) <<endl;
	
      }
    }
  }
  cout << index.size() << endl;
  TLegend *leg1 = new TLegend(0.2,0.35,0.8,0.85);
  gStyle->SetLegendBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetNColumns(4);
  gStyle->SetLegendTextSize(0.02);
 
  Int_t ii = 0; 
  Int_t ic,is,id, gain;
  for(size_t i=0;i<index.size();i++){
    ii=index.at(i);
    ic=fCryo.at(ii);
    is=fString.at(ii);
    id=fDetpos.at(ii);
    gain = fChannel.at(ii)%2;
    TitleName = Form("C%dP%dD%d, %d", fCryo.at(ii), fString.at(ii),fDetpos.at(ii),fChannel.at(ii));
    //TitleName = Form("#splitline{C%dP%dD%d}{%d}",fCryo.at(ii),fString.at(ii),fDetpos.at(ii),fChannel.at(ii));
    cout << ii << " " << TitleName.c_str() << endl;
    fun[ii]->SetLineColor(Color.at(is));
    fun[ii]->SetLineStyle(LineStyle.at(id));
    if(gain==0){
      fun[ii]->SetLineWidth(3);
    }else{
      fun[ii]->SetLineWidth(1);
    }
      //Pos = Form("%d%d%d",fCryo.at(ii),fString.at(ii),fDetpos.at(ii));
    fun[ii]->SetMaximum(8);
    fun[ii]->SetMinimum(0);
    leg1->AddEntry(fun[ii],TitleName.c_str());
  }

  for(size_t i=0;i<index.size();i++){
    ii=index.at(i);
    fun[ii]->Draw("same");
  }
  leg1->Draw();
  c1->Print(Form("resochan_%s_%d_%d.pdf",fEName.c_str(),fStartRun,fEndRun));

}


