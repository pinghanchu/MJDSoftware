#include "MJDGat.hh"
#include "GATPulserTag.hh"
#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TROOT.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 4) {
    cout << "Usage: " << argv[0] << " [InputFile1] [InputFile2] [Histogram]" << endl;
    return 1;
  }

  string fInput1 = argv[1];
  string fInput2 = argv[2];
  string fHistName = argv[3];
  
  TFile fin1(Form("%s",fInput1.c_str()),"read");
  TFile fin2(Form("%s",fInput2.c_str()),"read");
  TH1D* h1 = (TH1D*)fin1.Get(Form("%s",fHistName.c_str()));
  TH1D* h2 = (TH1D*)fin2.Get(Form("%s",fHistName.c_str()));
  h1->SetLineColor(4);
  h2->SetLineColor(2);
  TCanvas *c1 = new TCanvas("c1");

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  gStyle->SetOptTitle(0);
  //TCanvas *c1 = new TCanvas("c1");
  c1->SetLogy();
  h1->Rebin(10);
  h2->Rebin(10);
  h1->Draw("Hist");
  h2->Draw("Hist Same");
  c1->Print("hist.pdf");

}
