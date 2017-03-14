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


  GATAutoCal ds(fStartRun,fStartRun+1);
  string HistPathName = "./Hist/";
  string PathName = "./Plot/";
  string FileName = Form("roi_%d_%d.pdf",fStartRun,fEndRun);

  vector<Int_t> fChannel= ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fString = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad= ds.GetGoodBad();
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  string Pos;
  TCanvas *c1 = new TCanvas("c1");

  vector<Double_t> fRef;
  fRef.push_back(60);
  fRef.push_back(2039);
  fRef.push_back(2614.5);
  for(size_t ir=0;ir<fRef.size();ir++){
    for(size_t i =0;i<DataName.size();i++){
      Par.clear();
      ParErr.clear();
      GetROIPar(fStartRun,fEndRun,DataName.at(i),fEName, fRef.at(ir), &Par, &ParErr);
      Int_t pars = Par.size();
      if(Par.size()>0){
	froi << "000" << " "<< DataName.at(i).c_str() << " " <<fStartRun << " " << fEndRun << " "
	     << fEName.c_str() << " " << fRef.at(ir) << " " 
	     << Par.at(pars-2)<< "\\pm"<<ParErr.at(pars-2) <<" "
	     << Par.at(pars-1)<< "\\pm"<<ParErr.at(pars-1) <<endl;	
      }
    }
  }
}


