// 2016.2.18
// This macro uploads the calibration parameters.
// 
#include "GATAutoCal.hh"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

int main(int argc, char** argv)
{
  if(argc != 7 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [position] [channel] [energy 1] [energy 2]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  Int_t fPos = atoi(argv[3]);
  Int_t fChannel = atoi(argv[4]);
  string fEnr1 = argv[5];
  string fEnr2 = argv[6];
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");

  
  TChain *mjdTree = new TChain("mjdTree");
  for(Int_t i = fStartRun;i<=fEndRun;i++){
    mjdTree->Add(Form("./mjd_run%d.root",i));
  }

  TCanvas *c1 = new TCanvas("c1");
  TH2D *h1 = new TH2D("h1","",3000,0,3000,1000,-5,5);
  TH2D *h2 = new TH2D("h2","",3000,0,3000,1000,-5,5);

  mjdTree->Draw(Form("%s[j]-%s[i]:%s[i]>>h1",fEnr1.c_str(),fEnr1.c_str(),fEnr2.c_str()),Form("channel[i] == %d && channel[j] == %d",fChannel,fChannel+1));
  mjdTree->Draw(Form("%s[j]-%s[i]:%s[i]>>h2",fEnr2.c_str(),fEnr2.c_str(),fEnr2.c_str()),Form("channel[i] == %d && channel[j] == %d",fChannel,fChannel+1));

  h1->SetMarkerColor(2);
  h2->SetMarkerColor(4);
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h1->SetTitle(";Energy (keV); E_{low-gain}-E_{high-gain}");
  h2->SetTitle(";Energy (keV); E_{low-gain}-E_{high-gain}");

  TLegend *legDeltaPeak = new TLegend(0.65,0.65,0.85,0.85);
  legDeltaPeak->SetFillStyle(0);
  legDeltaPeak->SetBorderSize(0);
  legDeltaPeak->AddEntry(h1,"W/O NL Corr");
  legDeltaPeak->AddEntry(h2,"W NL Corr");

  h1->Draw();
  h2->Draw("same");
  legDeltaPeak->Draw();
  c1->Print(Form("highlow_%d_%d_%s_%s_%d_%d.pdf",fStartRun,fEndRun,fEnr1.c_str(),fEnr2.c_str(),fPos,fChannel));  

}


