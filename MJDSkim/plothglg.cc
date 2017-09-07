#include "MJDSkim.hh"
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
  if(argc != 2) {
    cout << "Usage: " << argv[0] << " [dataset]" << endl;
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  string fInputFile = Form("./List/runlist/DS%d.callist.txt",fDataSet);

  //TCanvas *c1 = new TCanvas("c1");
  ifstream fin(Form("%s",fInputFile.c_str()));

  Int_t startrun,endrun;
  vector<Int_t> StartRun;
  vector<Int_t> EndRun;

  if(fin.is_open()){
    while(!fin.eof()){
      fin>> startrun >> endrun;
      StartRun.push_back(startrun);
      EndRun.push_back(endrun);
    }
  }
  StartRun.pop_back();
  EndRun.pop_back();


  //GATAutoCal ps(StartRun.at(0),StartRun.at(0));
  MJDSkim ps(fDataSet,0,0,1);
  vector<Int_t> channel = ps.GetChannel();
  //vector<Int_t> goodbad = ps.GetGoodBad();
  TChain *skimTree = new TChain("skimTree");
  string path = "GAT-v01-06-125-gd9332b6";

  for(size_t i = 0;i<StartRun.size();i++){
    for(Int_t j = StartRun.at(i);j <= EndRun.at(i);j++){
      skimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%dcal/%s/skimDS%d_run%d_small.root",fDataSet,path.c_str(),fDataSet,j));
    }
  }
  skimTree->SetBranchStatus("*",1);
  vector<Int_t>* fChannel = NULL;
  vector<Double_t>* fEnergy = NULL;
  skimTree->SetBranchAddress("channel",&fChannel);
  skimTree->SetBranchAddress("trapENFCalC",&fEnergy);

  //TFile *newfile = new TFile("small.root","recreate");
  TTree *newtree = skimTree->CloneTree(0);
  for (Int_t i=0;i<skimTree->GetEntries(); i++) {
    skimTree->GetEntry(i);
    for(size_t j = 0;j<fChannel->size();j++){
      if (fEnergy->at(j) > 2037 && fEnergy->at(j) <2041){
	newtree->Fill();
      }
    }
  }
  //newtree->Print();
  // newtree->AutoSave();

  newtree->SetBranchStatus("*",1);
  vector<Int_t>* fNChannel = NULL;
  vector<Double_t>* fNEnergy = NULL;
  newtree->SetBranchAddress("channel",&fNChannel);
  newtree->SetBranchAddress("trapENFCalC",&fNEnergy);


  Int_t fEntries = newtree->GetEntries();
  cout << fEntries << endl;
  vector<Double_t> delta;
  vector<Int_t> hgchan;
  vector<Double_t> hgenr;
  vector<Int_t> lgchan;
  vector<Double_t> lgenr;
  for(Int_t i1=0;i1<fEntries;i1++){
    newtree->GetEntry(i1);
    hgchan.clear();
    lgchan.clear();
    hgenr.clear();
    lgenr.clear();
    for(size_t j=0;j<fNChannel->size();j++){
      Int_t chan = fNChannel->at(j);
      Double_t enr = fNEnergy->at(j);
      
      if(chan%2 ==0){
	hgchan.push_back(chan);
	hgenr.push_back(enr);
      }else{
	lgchan.push_back(chan);
	lgenr.push_back(enr);
      }
      for(size_t k1=0;k1<hgchan.size();k1++){
	for(size_t k2=0;k2<lgchan.size();k2++){
	  if((lgchan.at(k2)-hgchan.at(k1))==1){
	    if(hgenr.at(k1)>2037 && hgenr.at(k1)<2041){
	      Double_t del = lgenr.at(k2)-hgenr.at(k1);
	      if(abs(del)<5){
		delta.push_back(lgenr.at(k2)-hgenr.at(k1));
	      }else{
		cout << i1 << " " <<  lgenr.at(k2)-hgenr.at(k1) << endl;
	      }
	    }
	  }
	}
      }
    }
  }
  Double_t mean = 0;
  Int_t n = delta.size();
  for(size_t i=0;i<delta.size();i++){
    mean = mean + delta.at(i);
  }
  mean = mean/n;
  Double_t sigma = 0;
  Double_t mu = 0;
  for(size_t i=0;i<delta.size();i++){
    sigma = sigma + pow(delta.at(i)-mean,2);
    mu = mu + pow(delta.at(i)-mean,4);
  }
  sigma = sigma/n;
  mu = mu/n;
  Double_t sig1 = mean/2;
  Double_t sig2 = sqrt((mu-(n-3)/(n-1)*pow(sigma,2))/n);
  ofstream fout("NL.txt",ios::app);
  cout << "DS" << fDataSet << " " << mean << " " << sigma << " " << mu << " " << n << " "<<sig1 << " " << sig2 << endl;
  fout << "DS" << fDataSet << " " << mean << " " << sigma << " " << mu << " " << n <<  " "<< sig1 << " " << sig2 << endl;
  
  /*
  vector<Double_t> Mean;
  vector<Double_t> Sigma;
  vector<Double_t> Mu;
  Int_t totN = 0;
  for(Int_t ii =0;ii<channels/2;ii++){
    if(goodbad.at(ii*2)==1){
      TH1D *h1 = new TH1D("h1","",1000, -5,5);
      Int_t chan1 = channel.at(ii*2);
      Int_t chan2 = channel.at(ii*2+2);
      skimTree->Draw("trapENFCalC[i]C-trapENFCalC[j]>>h1",Form("trapENFCalC[i]>2037 && trapENFCalC[j]<2041 && channel[i] == %d && channel[j] == %d",chan1,chan2));
      Double_t mean = h1->GetMean();
      //Double_t rms = h1->GetRMS();
      Double_t entries = h1->GetEntries();
      Mean.push_back(mean);
      skimTree->Draw("pow((trapENFCalC[i]C-trapENFCalC[j]-mean),2)>>h1",Form("trapENFCalC[i]>2037 && trapENFCalC[j]<2041 && channel[i] == %d && channel[j] == %d",chan1,chan2));
      Double_t sigma = h1->GetMean();
      Sigma.push_back(sigma);
      skimTree->Draw("pow((trapENFCalC[i]C-trapENFCalC[j]-mean),4)>>h1",Form("trapENFCalC[i]>2037 && trapENFCalC[j]<2041 && channel[i] == %d && channel[j] == %d",chan1,chan2));
      Double_t mu = h1->GetMean();
      Mu.push_back(sigma);
      totN = totN + entries;
    }
  }
  Double_t s1=0;
  Double_t s2=0;
  Double_t s3=0;
  for(size_t i=0;i<Mean.size();i++){
    s1 = s1+Mean.at(i);
    s2 = s2+Sigma.at(i);
    s3 = s3+Mu.at(i);
  }
  Int_t n = Mean.size();
  s1 = s1/n;
  s2 = s2/n;
  s3 = s3/n;
  Double_t sig = s3-(totN-3)/(totN-1)*s2*s2)
*/
}

