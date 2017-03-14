// 2016.2.18
// written by Pinghan Chu
// Following the logic in GATAutoCal.cc
// 
#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include <string>
#include <stdlib.h>


int main(int argc, char** argv)
{
  if(argc != 3 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number]" << endl;
    return 1;
  }
  int fStartRun = atoi(argv[1]);
  int fEndRun = atoi(argv[2]);
  //  int Channel = atoi(argv[3]);
  //int EnergyPeak = atoi(argv[3]);
  //int EnergyWindow = atoi(argv[4]);
  gStyle->SetOptStat(1100);
  gStyle->SetOptFit(1111);
  TCanvas *c1 = new TCanvas("c1");

  Double_t hour;
  Double_t x1,x2,x3,x4;
  vector<Double_t> Peak;
  vector<Double_t> PeakErr;
  vector<Double_t> Run;
  vector<Double_t> RunErr;
  ifstream fin1(Form("./List/profiletime_%d_%d_2614_10.txt",fStartRun,fEndRun));
  if(fin1.is_open()){
    while(!fin1.eof()){
      fin1 >> hour >> x1 >> x2 >>x3 >>x4;      
      cout << hour << " " << x1 << endl;
      if(x1!=0){
	Peak.push_back(x1);
	PeakErr.push_back(x2);
	Run.push_back(hour);
	RunErr.push_back(0.5);
      }else{
	//	cout << pos << " " << cha << " " << run << " " << x1 << endl;
      }
    }
  }
  Peak.pop_back();
  PeakErr.pop_back();
  Run.pop_back();
  RunErr.pop_back();

  TGraphErrors fPeak;
  fPeak.Set(0);
  for(Int_t j =0;j<(int)Run.size();j++){
    fPeak.SetPoint(j,Run.at(j),Peak.at(j));
    fPeak.SetPointError(j,RunErr.at(j),PeakErr.at(j));  
  }  

  fPeak.SetTitle(Form("run%d-%d;Time(hour);Energy(keV)",fStartRun,fEndRun));
  fPeak.GetYaxis()->SetTitleOffset(5);
  fPeak.SetMarkerStyle(2);
  fPeak.SetMarkerSize(1);
  c1->Update();
  fPeak.Draw("AP");
  c1->Print(Form("./Plot/profilegraph_%d_%d.pdf", fStartRun,fEndRun));

  /*
  TH1D *hPeak = new TH1D("hPeak",Form("run%d-%d",StartRun,EndRun),100,EnergyPeak-5,EnergyPeak+5);
  for(Int_t j =0;j<(int)Run.size();j++){
    hPeak->Fill(Peak.at(j));
  }

  hPeak->SetTitle(Form("run%d-%d;peak(keV);counts",StartRun,EndRun));
  hPeak->GetYaxis()->SetTitleOffset(1.3);
  hPeak->SetMarkerStyle(2);
  hPeak->SetMarkerSize(1);
  c1->Update();
  hPeak->Draw();
  c1->Print(Form("./Plot/profilehisto_%d_%d_%d_%d.png", StartRun,EndRun,EnergyPeak,EnergyWindow));

  */
}


