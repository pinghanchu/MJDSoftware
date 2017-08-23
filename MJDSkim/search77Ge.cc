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
  if(argc != 6 ) {
    cout << "Usage: " << argv[0] << " [dataset][subset][isCal] [energy] [delay time]" << endl;
    return 1;
  }
  Int_t fDataSet = atoi(argv[1]);
  Int_t fSubSet = atoi(argv[2]);
  Int_t fIsCal = atoi(argv[3]);
  Double_t fEnergy = atof(argv[4]);
  Double_t fTime = atof(argv[5]);
  MJDSkim ds(fDataSet, fSubSet,fSubSet,fIsCal);
  Double_t fQ = 10000;
  Int_t fTimems = 114;
  string fOutputFile = Form("data_%d_%d_%d_%d.txt",fDataSet,fSubSet,(Int_t)fEnergy,fTimems);
  cout << "Searching candidates..." << "energy is " << fEnergy << "; the delayed time is " <<fTime << endl;
  ds.SearchDelayedEvent(fEnergy,fQ,fTime,fOutputFile);
  cout << "Searching candidates is done. "<< endl;
  /*
  ifstream fin(Form("%s",fOutputFile.c_str()));

  Int_t run1,list1,entry1,channel1;
  Double_t enr1,time1,mu_s1;
  Int_t run2,list2,entry2,channel2;
  Double_t enr2,time2,mu_s2;
  Double_t dcr1,dcr2;
  Double_t difftime;
  Double_t avse1,avse2;
  Double_t trapetailmin1,trapetailmin2;
  string M1, M2;
  vector<Int_t> Run1;
  vector<Int_t> Entry1;
  vector<Int_t> Channel1;
  vector<Double_t> Enr1;
  vector<Double_t> Time1;
  vector<Double_t> Mu_s1;
  vector<Double_t> DCR1;
  vector<Int_t> Run2;
  vector<Int_t> Entry2;
  vector<Int_t> Channel2;
  vector<Double_t> Enr2;
  vector<Double_t> Time2;
  vector<Double_t> Mu_s2;
  vector<Double_t> DCR2;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run1 >> list1 >> entry1 >> channel1 >> enr1 >> time1 >> M1 >> dcr1 >> avse1 >> trapetailmin1 >> 
	run2 >> list2 >> entry2 >> channel2 >> enr2 >> time2 >> M2 >> dcr2 >> avse2 >> trapetailmin2 >> difftime;      
      if(M1 == "nan" || M1 == "-nan" || M1 == "inf" || M1 == "-inf"){
	mu_s1 = 0;
      }else{
	mu_s1 = atof(M1.c_str());
      }
      if(M2 == "nan" || M2 == "-nan" || M2 == "inf" || M2 == "-inf"){
        mu_s2 = 0;
      }else{
        mu_s2 = atof(M2.c_str());
      }
      //cout << count << " "<< run1 << " "<< list1 << " " << entry1 << endl;
      //if(enr1>53-5 && enr1<67+5 && !((enr1>49 && enr1<51) && (dcr1>0.015 && dcr1<0.018))  && !( enr1>49 && enr1<52 && dcr1>0.004 && dcr1<0.006) && !(enr1>48&&enr1<58 && dcr1<-0.004 && dcr1>-0.006)){
      Run1.push_back(run1);
      Entry1.push_back(entry1);
      Channel1.push_back(channel1);
      Enr1.push_back(enr1);
      Time1.push_back(time1);
      Mu_s1.push_back(mu_s1);
      DCR1.push_back(dcr1);
      Run2.push_back(run2);
      Entry2.push_back(entry2);
      Channel2.push_back(channel2);
      Enr2.push_back(enr2);
      Time2.push_back(time2);
      Mu_s2.push_back(mu_s2);
      DCR2.push_back(dcr2);
      
    }
  } 
  vector<Int_t> Index;
  if(Run1.size()>0){
    Index.push_back(0);    
    for(size_t i=1;i<Run1.size();i++){
      if(Run1.at(i)!=Run1.at(i-1) || Entry1.at(i)!=Entry1.at(i-1) || Channel1.at(i)!=Channel1.at(i-1)){
	Index.push_back(i);
      }else{
      }
    }
  }
  cout << Index.size() << endl;
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
    if(fDataSet == 2 || fDataSet == 6){
      Xmax = 38000;
    }
    Double_t fResolution = 1;
    Double_t fThreshold = 0.01;
    Double_t fSigma = 5;
    
    for(size_t i=0;i<Index.size();i++){
      Int_t ii = Index.at(i);
      //cout << Run1.at(ii) << " " << Entry1.at(ii) << " " << Channel1.at(ii) << " " << Enr1.at(ii) << " " << Mu_s1.at(ii) << " " << DCR1.at(ii) << endl;
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
    }
  }
  */
  cout << "Scaning is "<< fEnergy << "done. " << fDataSet << " "<< fSubSet << endl;

}
