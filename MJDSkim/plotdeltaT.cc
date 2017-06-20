#include "MJDSkim.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 4 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << "[energy] [window] [input file]" << endl;
    return 1;
  }

  Double_t fEnr = atof(argv[1]);
  Double_t fWindow =atof(argv[2]);
  string fInputFile = argv[3];
  //string fInputFile = Form("./data/wf_%d_%d.txt",(Int_t)fEnr,(Int_t)fWindow);
  cout << fInputFile.c_str() << endl;
  ifstream fin(Form("%s",fInputFile.c_str()));

  Int_t run,entry,channel,nx;
  Double_t enr,ratio,deltaT,maxFFT,aovere,Y1,Y2;
  vector<Double_t> Ratio;
  vector<Double_t> DeltaT;
  vector<Double_t> Energy;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run >> entry >> channel >> enr >> nx >> ratio >> deltaT >> maxFFT >> aovere; 
      //fin >> run >> entry >> channel >> enr >> ratio >> deltaT >> maxFFT >> aovere ;
      //if(ratio > 3.7 && ratio <4.5 && aovere>0.004){
	Ratio.push_back(ratio);
	DeltaT.push_back(deltaT);
	Energy.push_back(enr);
	cout<<run<<" " <<entry << " "<< channel << " " << enr << " "<< ratio << " " <<deltaT << endl;
	//}
    }
  }

  TCanvas *c1 = new TCanvas("c1");
  if(Ratio.size()>0){
    TH1D *h1 = new TH1D("h1","",10000,0,10000);
    TH2D *h2 = new TH2D("h2","",1000,0,10000,100,0,20);
    for(size_t i = 0;i<DeltaT.size();i++){
      h1->Fill(DeltaT.at(i));
      h2->Fill(DeltaT.at(i),Ratio.at(i));
    }
    h1->SetTitle(";#Delta T(ns)");
    h1->GetYaxis()->SetTitleOffset(1.3);
    h1->Draw();
    c1->Print(Form("deltaT_%d_%d.pdf",(Int_t)fEnr,(Int_t)fWindow));
    h2->SetTitle(";#Delta T(ns);Ratio");
    h2->GetYaxis()->SetTitleOffset(1.3);
    h2->Draw("COLZ");
    c1->Print(Form("deltaTratio_%d_%d.pdf",(Int_t)fEnr,(Int_t)fWindow));


  }
}
