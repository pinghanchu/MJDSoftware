#include "MJDSkim.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

Double_t Exp(Double_t *v, Double_t *par){
  Double_t fitval = 0;
  fitval = par[0]*TMath::Exp(-v[0]/par[1]);
  return fitval;
}

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
  Int_t run1,run2,entry1,entry2;
  Double_t time1,time2,enr1,enr2,deltaT1,deltaT2;
  string channel1,channel2;
  string eventid, dataset;
  string isenr,isgood;
  if(fin.is_open()){
    while(!fin.eof()){
      //fin >> run >> entry >> channel >> enr >> nx >> ratio >> deltaT >> maxFFT >> aovere; 
      //fin >> run >> entry >> channel >> enr >> ratio >> deltaT >> maxFFT >> aovere ;
      fin >> eventid >> dataset >> run1 >> entry1 >> channel1 >> time1 >> enr1 >> run2 >> entry2 >> channel2 >> time2 >> enr2 >> deltaT1 >> deltaT2 >> isenr >> isgood; 
      cout << eventid << " "<< isgood << " "<< deltaT2 << endl;
      //if(ratio > 3.7 && ratio <4.5 && aovere>0.004){
      //Ratio.push_back(ratio);
      DeltaT.push_back(deltaT1);
      Energy.push_back(enr2);
	//cout<<run<<" " <<entry << " "<< channel << " " << enr << " "<< ratio << " " <<deltaT << endl;
	//}
    }
  }
  TF1 *fun = new TF1("fun", Exp,0,10,2);
  fun->FixParameter(1,0.5);
  fun->SetParameter(0,1);
  TCanvas *c1 = new TCanvas("c1");
  TH1D *h1 = new TH1D("h1","",200,0,20);
  for(size_t i = 0;i<DeltaT.size();i++){
    h1->Fill(DeltaT.at(i));
  }
  h1->Fit("fun");
  h1->SetTitle(";#Delta T(s)");
  h1->GetYaxis()->SetTitleOffset(1.3);
  //c1->SetLogy();
  h1->Draw();
  //fun->Draw("same");
  c1->Print(Form("deltaT_%d_%d.pdf",(Int_t)fEnr,(Int_t)fWindow));

  /*
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
  */
}
