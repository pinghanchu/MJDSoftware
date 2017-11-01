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
  if(argc != 4 ) {
    cout << "Usage: " << argv[0] << " [run] [energy] [delay time]" << endl;
    return 1;
  }
  Int_t fRun = atoi(argv[1]);
  Double_t fEnergy = atof(argv[2]);
  Double_t fTime = atof(argv[3]);
  MJDGat ds(fRun);
  Double_t fQ = 10000;
  Int_t fTimems = (Int_t) 500;
  string fOutputFile = Form("data_73Ga_%d",fRun);
  cout << "Searching candidates..." << "energy is " << fEnergy << "; the delayed time is " <<fTimems << endl;
  ds.SearchDelayedEvent(fEnergy,fQ,fTime,fOutputFile);
  cout << "Searching candidates is done. "<< endl;

  ifstream fin(Form("%s.csv",fOutputFile.c_str()));

  string run1,list1,entry1,channel1;
  string enr1,time1;
  string list2,entry2,channel2;
  string enr2,time2;

  vector<Int_t> Run1;
  vector<Int_t> Entry1;
  vector<Int_t> Channel1;
  vector<Double_t> Enr1;
  vector<Double_t> Time1;
  vector<Int_t> Entry2;
  vector<Int_t> Channel2;
  vector<Double_t> Enr2;
  vector<Double_t> Time2;
  while(getline(fin,run1,',')){
    getline(fin,entry1,',');
    getline(fin,channel1,',');
    getline(fin,enr1,',');
    getline(fin,time1,',');
    getline(fin,entry2,',');
    getline(fin,channel2,',');
    getline(fin,enr2,',');
    getline(fin,time2,'\n');
    Run1.push_back(stoi(run1));
    Entry1.push_back(stoi(entry1));
    Channel1.push_back(stoi(channel1));
    Enr1.push_back(stof(enr1));
    Time1.push_back(stof(time1));
    Entry2.push_back(stoi(entry2));
    Channel2.push_back(stoi(channel2));
    Enr2.push_back(stof(enr2));
    Time2.push_back(stof(time2));
  }

  vector<Int_t> Index;
  if(Run1.size()>0){
    Index.push_back(0);
  }
  if(Run1.size()>0){
    for(size_t i=1;i<Run1.size();i++){
      if( Run1.at(i)!=Run1.at(i-1) || Entry1.at(i)!=Entry1.at(i-1) || Channel1.at(i)!=Channel1.at(i-1)){
	Index.push_back(i);
      }
    }
  }

  string fOutputFile1 = Form("wf_73Ga_%d",fRun);
  ofstream fout(Form("%s.csv",fOutputFile1.c_str()));
  TCanvas *c1 = new TCanvas("c1");
  if(Index.size()>0){    
    vector<Double_t> xp;
    vector<Double_t> yp;
    vector<Double_t> xp1;
    vector<Double_t> yp1;
    vector<Double_t> xp2;
    vector<Double_t> yp2;
    
    Double_t Xmin = 4000;
    Double_t Xmax = 19000;
    if(fRun >= 25672){
      Xmax = 38000;
    }
    Double_t fResolution = 1;
    Double_t fThreshold = 0.01;
    Double_t fSigma = 5;
    
    for(size_t i=0;i<Index.size();i++){
      Int_t ii = Index.at(i);
      xp.clear();
      yp.clear();
      xp1.clear();
      yp1.clear();
      xp2.clear();
      yp2.clear();
      
      Int_t run = Run1.at(ii);
      Int_t event = Entry1.at(ii);
      Int_t chan = Channel1.at(ii);
      Double_t enr = Enr1.at(ii);
      Double_t time = Time1.at(ii);

      TH1D* h = ds.GetWaveform(run,event,chan,enr);
      TH1D* h2 = ds.GetHistoDerivative(h,10);
      TH1D* h4 = ds.GetHistoDerivative(h2,10);
      TH1* hFFT = ds.GetHistoFFT(h);
      string name1 = Form("waveform_%d_%d_%d",run,event,chan);
      string name2 = Form("waveform1_%d_%d_%d",run,event,chan);
      string name3 = Form("waveform2_%d_%d_%d",run,event,chan);

      h2->SetName(Form("%s",name2.c_str()));
      h4->SetName(Form("%s",name3.c_str()));
      h->Draw();
      c1->Print(Form("%s.pdf",name1.c_str()));
      h2->Draw();
      c1->Print(Form("%s.pdf",name2.c_str()));
      h4->Draw();
      c1->Print(Form("%s.pdf",name3.c_str()));

      //h2->GetXaxis()->SetRangeUser(0,20000);
      //h3->Write();
      //h5->Write();
      //hFFT->GetXaxis()->SetRangeUser(0,1000);
      //hFFT->Write();


      // Cut1
      Double_t maxFFT = ds.GetMax(hFFT, 50.,150.); 
      //Double_t maxY0 = ds.GetMax(h1, 100,Xmax);
      Int_t maxBin0 = h->GetMaximumBin();
      Double_t maxX0 = h->GetBinCenter(maxBin0);
      
      //cout << max << endl;
      // Cut2 
      Double_t A = ds.GetMax(h2, Xmin,Xmax);
      Double_t AoverE = A/enr;
      
      Int_t nPeak1 = ds.FindPeaks(h2, Xmin, Xmax,fResolution,fSigma, fThreshold, &xp, &yp);
      if(nPeak1>1 && AoverE>0.002){
	for(Int_t ip = 0;ip<(Int_t)xp.size();ip++){    
	  cout << xp.at(ip) << " " << yp.at(ip) << endl;
	  Double_t y = ds.GetYValue(h4, xp.at(ip));
	  if(abs(y)< 3e-4 && xp.at(ip)< maxX0){
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
	    
	  fout << run << "," << event << "," << chan << "," << enr << "," << time << ","
	       << yp2.at(1)/yp2.at(0) << "," << xp2.at(0)-xp2.at(1) << "," << AoverE << ","<< maxFFT << endl;    
	}
	
	delete h;
	delete h2;
	delete h4;
      }
      //      fhist.Close();

    }
  }

  cout << "Selecting candidates is done. "<< endl;


}
