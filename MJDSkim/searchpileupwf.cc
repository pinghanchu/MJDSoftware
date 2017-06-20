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
  if(argc != 4 ) {
    cout << "Usage: " << argv[0] << " [dataset] [energy] [window]" << endl;
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  Double_t fEnergy = atof(argv[2]);
  Double_t fWindow = atof(argv[3]);

  MJDSkim ds(fDataSet);
  //ds.SetEnergyName(fEnergyName);
  //Double_t fEnr = fEnergy;
  //Double_t fTime = 0.5;
  string fOutputFile = Form("data_%d_%d_%d.txt",fDataSet,(Int_t)fEnergy,(Int_t)fWindow);
  ds.SearchEnergyEvent(fEnergy,fWindow,fOutputFile);

  ifstream fin(Form("%s",fOutputFile.c_str()));

  Int_t run,entry,pos,channel,nX;
  Double_t enr;
  vector<Int_t> Run;
  vector<Int_t> Entry;
  vector<Int_t> Channel;
  vector<Double_t> Enr;
  vector<Int_t> NX;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run >> entry >> pos >> channel >> enr >> nX;
      Run.push_back(run);
      Entry.push_back(entry);
      Channel.push_back(channel);
      Enr.push_back(enr);
      NX.push_back(nX);
    }
  }
  if(Entry.size()>0){
    Run.pop_back();
    Entry.pop_back();
    Channel.pop_back();
    Enr.pop_back();
    NX.pop_back();
  }

  if(Run.size()>0){
    
    string fOutputWaveform = Form("waveform_%d_%d_%d.root",fDataSet,(Int_t)fEnergy,(Int_t)fWindow);
    TFile fhist(Form("%s",fOutputWaveform.c_str()),"recreate");
    ofstream fout(Form("wf_%d_%d_%d.txt",fDataSet,(Int_t)fEnergy,(Int_t)fWindow));
   
    vector<Double_t> xp;
    vector<Double_t> yp;
    vector<Double_t> xp1;
    vector<Double_t> yp1;
    vector<Double_t> xp2;
    vector<Double_t> yp2;
    
    Double_t Xmin = 4000;
    Double_t Xmax = 19000;
    if(fDataSet == 2){
      Xmax = 38000;
    }
    Double_t fResolution = 1;
    Double_t fThreshold = 0.01;
    Double_t fSigma = 5;
    
    for(size_t i=0;i<Run.size();i++){
      cout << Run.at(i) << " " << Entry.at(i) << " " << Channel.at(i) << " " << Enr.at(i) << " " << NX.at(i) << endl;
      //if(NX.at(i)>1 && NX.at(i)<10){
      if(NX.at(i)>0){
	xp.clear();
	yp.clear();
	xp1.clear();
	yp1.clear();
	xp2.clear();
	yp2.clear();
	
	Int_t run = Run.at(i);
	Int_t event = Entry.at(i);
	Int_t chan = Channel.at(i);
	Double_t enr = Enr.at(i);
	TH1D* h = ds.GetWaveform(run,event,chan,enr);
	h->Write();
	TH1D* h1 = ds.GetHistoSmooth(h,10);
	TH1D* h2 = ds.GetHistoDerivative(h1,10);
	TH1D* h3 = ds.GetHistoSmooth(h2,10);
	TH1D* h4 = ds.GetHistoDerivative(h3,10);
	TH1D* h5 = ds.GetHistoSmooth(h4,10);
	TH1*  hFFT = ds.GetHistoFFT(h);
	// Cut1
	Double_t maxFFT = ds.GetMax(hFFT, 50.,150.); 
	Double_t maxY0 = ds.GetMax(h1, 100,19000);
	Int_t maxBin0 = h1->GetMaximumBin();
	Double_t maxX0 = h1->GetBinCenter(maxBin0);

	//cout << max << endl;
	// Cut2 
	Double_t A = ds.GetMax(h3, Xmin,Xmax);
	Double_t AoverE = A/enr;
	
	if(maxFFT<4000 && AoverE>0.001){
	  Int_t nPeak1 = ds.FindPeaks(h3, Xmin, Xmax,fResolution,fSigma, fThreshold, &xp, &yp);
	  if(nPeak1>1){
	    for(Int_t ip = 0;ip<(Int_t)xp.size();ip++){
	      
	      Double_t y = ds.GetYValue(h5, xp.at(ip));
	      if(abs(y)<2e-4 && xp.at(ip)< maxX0){
		xp1.push_back(xp.at(ip));
		yp1.push_back(yp.at(ip));
	      }
	    }
	    
	    if(xp1.size()>1){
	      h3->SetName(Form("waveform1_%d_%d_%d",run,event,chan));
	      h5->SetName(Form("waveform2_%d_%d_%d",run,event,chan));
	      hFFT->SetName(Form("waveformFFT_%d_%d_%d",run,event,chan));
	      h1->Write();
	      h3->GetXaxis()->SetRangeUser(0,20000);
	      h3->Write();
	      h5->Write();
	      hFFT->GetXaxis()->SetRangeUser(0,1000);
	      hFFT->Write();
	      
	      vector<Int_t> Index2 = ds.Sort(yp1);
	      for(size_t ip1 = Index2.size()-2;ip1<Index2.size();ip1++){
		Int_t ii = Index2.at(ip1);
		xp2.push_back(xp1.at(ii));
		yp2.push_back(yp1.at(ii));
	      }
	      fout << run << " " << event << " " << chan << " " << enr << " " << yp2.at(1)/yp2.at(0) << " " << xp2.at(0)-xp2.at(1) << " " << maxFFT << " " << AoverE << endl;
	      
	    }
	  }	
	  delete h;
	  delete h1;
	  delete h2;
	  delete h3;
	  delete h4;
	  delete h5;
	  delete hFFT;
	}
      }
    }
    fhist.Close();
  }

}
