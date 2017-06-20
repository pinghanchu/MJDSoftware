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
  if(argc != 4 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run] [energy] [window" << endl;
    return 1;
  }

  Int_t fRun = atoi(argv[1]);
  cout << "Run:" << fRun << endl;
  Double_t fEnergy = atof(argv[2]);
  Double_t fWindow = atof(argv[3]);

  MJDGat ds(fRun);
  //ds.SetEnergyName(fEnergyName);
  string fDataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();

  cout << "Energy Screening..." << endl;
  string fOutputFile = Form("data_%d.txt",fRun);
  ds.SearchEnergyEvent(fEnergy,fWindow,fOutputFile);
  /*  
  cout << "Waveform filtering..." << endl;
  ifstream fin(Form("%s",fOutputFile.c_str()));

  Int_t run,entry,channel;
  Double_t enr,nx;
  vector<Int_t> Entry;
  vector<Int_t> Channel;
  vector<Double_t> Enr;
  vector<Double_t> nX;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run >> entry >> channel >> enr >> nx;
      Entry.push_back(entry);
      Channel.push_back(channel);
      Enr.push_back(enr);
      nX.push_back(nx);
    }
  }
  if(Entry.size()>0){
    Entry.pop_back();
    Channel.pop_back();
    Enr.pop_back();
    nX.pop_back();
  }
  
  string fOutputWaveform = Form("waveform_%d.root",fRun);
  TFile fhist(Form("%s",fOutputWaveform.c_str()),"recreate");
      
  //  for(size_t i=0;i<Entry.size();i++){
  for(size_t i=0;i<1000;i++){
    cout << fRun << " " << Entry.at(i) << " " << Channel.at(i) << " " << Enr.at(i) << " "<< nX.at(i) << endl;
    TH1D* h = ds.GetWaveform(fRun,Entry.at(i),Channel.at(i),Enr.at(i));
    TH1D* h1 = ds.GetHistoSmooth(h,10);
    TH1D* h2 = ds.GetHistoDerivative(h1,10);
    TH1D* h3 = ds.GetHistoSmooth(h2,10);
    TH1D* h4 = ds.GetHistoDerivative(h3,10);
    TH1D* h5 = ds.GetHistoSmooth(h4,10);
    TH1*  hFFT = ds.GetHistoFFT(h);
    Double_t maxFFT = ds.GetMax(hFFT, 50.,150.);
    Double_t maxY0 = ds.GetMax(h1, 100,19000);
    Int_t maxBin0 = h1->GetMaximumBin();
    Double_t maxX0 = h1->GetBinCenter(maxBin0);
    Double_t A = ds.GetMax(h3, Xmin,Xmax);
    Double_t AoverE = A/Enr.at(ii);
	
    //if(maxFFT<2000 && AoverE>0.004){
      h3->SetName(Form("waveform1_%d_%d_%d",fRun,Entry.at(i),Channel.at(i)));
      h5->SetName(Form("waveform2_%d_%d_%d",fRun,Entry.at(i),Channel.at(i)));
      hFFT->SetName(Form("waveformFFT_%d_%d_%d",fRun,Entry.at(i),Channel.at(i)));
      
      h1->Write();
      h3->GetXaxis()->SetRangeUser(0,20000);
      h3->Write();
      h5->Write();
      hFFT->GetXaxis()->SetRangeUser(0,1000);
      hFFT->Write();      
      // }
    fhist.Close();
  }
  */
}
