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
  fitval = par[0]*TMath::Exp(-TMath::Log(2)*v[0]/par[1]);
  return fitval;
}

int main(int argc, char** argv)
{
  if(argc != 2 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << "[input file]" << endl;
    return 1;
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  string fInputFile = argv[1];
  //string fInputFile = Form("./data/wf_%d_%d.txt",(Int_t)fEnr,(Int_t)fWindow);
  cout << fInputFile.c_str() << endl;
  ifstream fin(Form("%s",fInputFile.c_str()));

  //Int_t run,entry,channel,nx;
  //Double_t enr,ratio,deltaT,maxFFT,aovere,Y1,Y2;
  //vector<Double_t> Ratio;
  vector<Double_t> DeltaT;
  vector<Double_t> Energy1;
  vector<Double_t> Energy2;
  Int_t dataset,subset,list1,run1,event1,chan1,list2,run2,event2,chan2;
  Double_t time1,enr1,dcr1,avse1,trapetailmin1,mu1,time2,enr2,dcr2,avse2,trapetailmin2,mu2,deltaT;

  if(fin.is_open()){
    while(!fin.eof()){
      fin >> dataset >> subset >> 
	list1 >> run1 >> event1 >> time1 >> 
	chan1 >> enr1 >> dcr1 >> avse1 >> trapetailmin1 >> mu1 >> 
	list2 >> run2 >> event2 >> time2 >> 
	chan2 >> enr2 >> dcr2 >> avse2 >> trapetailmin2 >> mu2 >>
	deltaT;
      DeltaT.push_back(deltaT);
      Energy1.push_back(enr1);
      Energy2.push_back(enr2);
    }
  }
  DeltaT.pop_back();
  Energy1.pop_back();
  Energy2.pop_back();

  TCanvas *c1 = new TCanvas("c1");
  TH2D *h1 = new TH2D("h1","",400,1000,5000,400,1000,5000);
  TH2D *h2 = new TH2D("h2","",400,1000,5000,100,0,1000);
  TH1D *h3 = new TH1D("h3","",4000,1000,5000);
  TH1D *h5 = new TH1D("h5","",1000,0,1000);
  TH1D *h4 = new TH1D("h4","",4000,1000,5000); 
  for(size_t i = 0;i<Energy1.size();i++){
    h1->Fill(Energy1.at(i),Energy2.at(i));
    h2->Fill(Energy1.at(i),DeltaT.at(i));
    h3->Fill(Energy1.at(i));
    h4->Fill(Energy2.at(i));
    h5->Fill(DeltaT.at(i));
  }
  
  // c1->SetLogx();
  //c1->SetLogy();
  h1->SetTitle(";E_{1}(keV);E_{2}(keV)");
  h1->GetYaxis()->SetTitleOffset(1.3);
  h1->Draw("COLZ");
  c1->Print("77mGe_HighEnr_HighEnr.pdf");

  h2->SetTitle(";E_{1}(keV);#Delta T(s)");
  h2->GetYaxis()->SetTitleOffset(1.3);
  h2->Draw("COLZ");
  c1->Print("77mGe_HighEnr_DeltaT.pdf");

  h3->SetTitle(";E_{1}(keV)");
  h3->GetYaxis()->SetTitleOffset(1.3);
  h3->Draw();
  c1->Print("77mGe_HighEnr1.pdf");
  c1->SetLogy();
  h4->SetTitle(";E_{2}(keV)");
  h4->GetYaxis()->SetTitleOffset(1.3);
  h4->Draw();
  c1->Print("77mGe_HighEnr2.pdf");
  c1->SetLogy(0);
  h5->SetTitle(";#Delta T(s)");
  h5->GetYaxis()->SetTitleOffset(1.3);
  h5->Draw();
  c1->Print("77mGe_DeltaT.pdf");


}
