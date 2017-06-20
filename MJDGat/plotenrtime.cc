#include "MJDGat.hh"
#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2 ) {
    cout << "Usage: " << argv[0] << " [input file]" << endl;
    return 1;
  }
  string fInputFile = argv[1];
  ifstream fin(Form("%s",fInputFile.c_str()));
  Int_t run,entry,channel;
  Double_t enr,deltaT;
  vector<Double_t> DeltaT;
  vector<Double_t> Energy;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run >> entry >> channel >> enr >> deltaT ;
      DeltaT.push_back(deltaT);
      Energy.push_back(enr);
    }
  }
  if(Energy.size()>0){
    DeltaT.pop_back();
    Energy.pop_back();
  }
  TCanvas *c1 = new TCanvas("c1");
  if(Energy.size()>0){
    TH1D *h1 = new TH1D("h1","",10,210,220);

    for(size_t i = 0;i<Energy.size();i++){
      h1->Fill(Energy.at(i));
    }
    h1->SetTitle(";Energy (keV);");
    h1->GetYaxis()->SetTitleOffset(1.3);
    h1->Draw();
    c1->Print("energy.pdf");

    TH1D *h2 = new TH1D("h2","",100,0,1000);

    for(size_t i = 0;i<Energy.size();i++){
      h2->Fill(DeltaT.at(i));
    }
    h2->SetTitle(";#Delta T(s);");
    h2->GetYaxis()->SetTitleOffset(1.3);
    h2->Draw();
    c1->Print("deltaT.pdf");

    TH2D *h3 = new TH2D("h3","",10,210,220,100,0,1000);

    for(size_t i = 0;i<Energy.size();i++){
      h3->Fill(Energy.at(i),DeltaT.at(i));
    }
    h3->SetTitle(";Energy (keV);#Delta T(s)");
    h3->GetYaxis()->SetTitleOffset(1.3);
    h3->Draw("COLZ");
    c1->Print("energydeltaT.pdf");
  }
}
