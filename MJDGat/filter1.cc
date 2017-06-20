#include "MJDGat.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 5 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run] [entry] [channel] [energy] " << endl;
    return 1;
  }

  Int_t fRun = atoi(argv[1]);
  cout << "Run:" << fRun << endl;
  Int_t fEntry = atoi(argv[2]);
  Int_t fChan = atoi(argv[3]);
  Double_t fEnergy = atof(argv[4]);


  MJDGat ds(fRun);
  //ds.SetEnergyName(fEnergyName);
  string fDataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
/*
  Double_t fstarttime = ds.GetStartTime();
  Double_t fendtime = ds.GetStopTime();

  ofstream fout("test.txt",ios::app);
  fout.precision(15);
  fout << fRun << " " << fstarttime << " " << fendtime << " "<<fendtime-fstarttime << endl;
*/
  
  cout << "Read Waveform" << endl;
  string fOutputWaveform = Form("waveform_%d.root",fRun);
  TFile fhist(Form("./data/%s",fOutputWaveform.c_str()),"read");
  ofstream fout(Form("wf_%d.txt",fRun),ios::app);
    
  vector<Double_t> xp;
  vector<Double_t> yp;
  vector<Double_t> xp1;
  vector<Double_t> yp1;
  vector<Double_t> xp2;
  vector<Double_t> yp2;
  
  Double_t Xmin = 4000;
  Double_t Xmax = 19000;
  if(fRun>=14503 && fRun<=15892){
    Xmax = 38000;
  }
  Double_t fResolution = 1;
  Double_t fThreshold = 0.01;
  Double_t fSigma = 5;
  
  xp.clear();
  yp.clear();
  xp1.clear();
  yp1.clear();
  xp2.clear();
  yp2.clear();
  
  TH1D* h1 = (TH1D*)fhist.Get(Form("waveform_%d_%d_%d",fRun,fEntry,fChan));
  TH1D* h3 = (TH1D*)fhist.Get(Form("waveform1_%d_%d_%d",fRun,fEntry,fChan));
  TH1D* h5 = (TH1D*)fhist.Get(Form("waveform2_%d_%d_%d",fRun,fEntry,fChan));
  TH1*  hFFT = (TH1D*)fhist.Get(Form("waveformFFT_%d_%d_%d",fRun,fEntry,fChan));
  Double_t maxFFT = ds.GetMax(hFFT, 50.,150.); 
  Double_t maxY0 = ds.GetMax(h1, 100,19000);
  Int_t maxBin0 = h1->GetMaximumBin();
  Double_t maxX0 = h1->GetBinCenter(maxBin0);
  //cout << max << endl;
  // Cut2 
  Double_t A = ds.GetMax(h3, Xmin,Xmax);
  Double_t AoverE = A/fEnergy;
  
  if(maxFFT<2000 && AoverE>0.004){
    Int_t nPeak1 = ds.FindPeaks(h3, Xmin, Xmax,fResolution,fSigma, fThreshold, &xp, &yp);
    if(nPeak1>1){
      for(Int_t ip = 0;ip<(Int_t)xp.size();ip++){
	
	Double_t y = ds.GetYValue(h5, xp.at(ip));
	if(abs(y)<2e-4){
	  xp1.push_back(xp.at(ip));
	  yp1.push_back(yp.at(ip));
	}
      }
      
      if(xp1.size()>1){
	vector<Int_t> Index2 = ds.Sort(yp1);
	for(size_t ip1 = Index2.size()-2;ip1<Index2.size();ip1++){
	  Int_t ii = Index2.at(ip1);
	  xp2.push_back(xp1.at(ii));
	  yp2.push_back(yp1.at(ii));
	}
	if(xp2.at(0)<maxX0 && xp2.at(1)<maxX0){
	  fout << fRun << " " << fEntry << " " << fChan << " "<< fEnergy << " " << 0 << " " << yp2.at(1)/yp2.at(0) << " " << xp2.at(0)-xp2.at(1) << " " << maxFFT << " " << AoverE << endl;

	    //" "<< maxX0 << " "<< maxY0 << " " << xp2.at(0) << " " << ds.GetYValue(h1,xp2.at(0)) << " " << xp2.at(1) << " "<< ds.GetYValue(h1,xp2.at(1)) << endl;	      
	}
      }
    }
  }	
}
