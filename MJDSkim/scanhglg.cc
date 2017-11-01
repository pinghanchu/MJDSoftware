#include "MJDSkim.hh"
#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;
void GetDeltaNL(Int_t DataSet, Double_t EnergyLow, Double_t EnergyUp, vector<Double_t>* Delta){

  string fInputFile = Form("/global/projecta/projectdirs/majorana/users/pchu/git/MJDSoftware/MJDSkim/List/runlist/DS%d.callist.txt",DataSet);
  //cout << fInputFile << endl;
  ifstream fin(Form("%s",fInputFile.c_str()));

  Int_t startrun,endrun;
  vector<Int_t> StartRun;
  vector<Int_t> EndRun;

  if(fin.is_open()){
    while(!fin.eof()){
      fin>> startrun >> endrun;
      //cout << startrun << " " << endrun << endl;
      StartRun.push_back(startrun);
      EndRun.push_back(endrun);
    }
  }
  StartRun.pop_back();
  EndRun.pop_back();

  MJDSkim ps(DataSet,StartRun.at(0),StartRun.at(0),1);
  vector<Int_t> channel = ps.GetChannel();
  TChain *skimTree = new TChain("skimTree");
  string path = "GAT-v01-07";
  
  for(size_t i = 0;i<StartRun.size();i++){
    for(Int_t j = StartRun.at(i);j <= EndRun.at(i);j++){
      skimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%dcal/%s/skimDS%d_run%d_small.root",DataSet,path.c_str(),DataSet,j));
      //cout << j << endl;
    }
  }

  skimTree->SetBranchStatus("*",1);
  vector<Int_t>* fChannel = NULL;
  vector<Double_t>* fEnergy = NULL;
  skimTree->SetBranchAddress("channel",&fChannel);
  skimTree->SetBranchAddress("trapENFCalC",&fEnergy);

  //TFile *newfile = new TFile("small.root","recreate");
  TTree *newtree = skimTree->CloneTree(0);
  cout << skimTree->GetEntries() << endl;
  for (Int_t i=0;i<skimTree->GetEntries(); i++) {
    //cout << i << endl;
    skimTree->GetEntry(i);
    for(size_t j = 0;j<fChannel->size();j++){
      if (fEnergy->at(j) > EnergyLow && fEnergy->at(j) <EnergyUp){
	newtree->Fill();
      }
    }
  }

  newtree->SetBranchStatus("*",1);
  vector<Int_t>* fNChannel = NULL;
  vector<Double_t>* fNEnergy = NULL;
  newtree->SetBranchAddress("channel",&fNChannel);
  newtree->SetBranchAddress("trapENFCalC",&fNEnergy);

  Int_t fEntries = newtree->GetEntries();
  cout << fEntries << endl;
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
	    if(hgenr.at(k1)>EnergyLow && hgenr.at(k1)<EnergyUp){
	      Double_t del = lgenr.at(k2)-hgenr.at(k1);
	      cout << i1 << " " << del << endl;
	      if(abs(del)<5){
		Delta->push_back(lgenr.at(k2)-hgenr.at(k1));
	      }else{
	      }
	    }
	  }
	}
      }
    }
  }
}

int main(int argc, char** argv)
{
  if(argc != 3) {
    cout << "Usage: " << argv[0] << " [low] [up]" << endl;
    return 1;
  }
  Int_t EnergyLow10 = stoi(argv[1]);
  Int_t EnergyUp10 = stoi(argv[2]);
  Double_t EnergyLow = (Double_t)EnergyLow10/10;
  Double_t EnergyUp = (Double_t)EnergyUp10/10;

  vector<Int_t> fDataSet;
  for(int id = 0;id<=5;id++){
    fDataSet.push_back(id);
  }

  vector<Double_t> delta;
  vector<Double_t> DeltaNL;
  DeltaNL.clear();
  for(size_t i = 0;i<fDataSet.size();i++){
    delta.clear();
    cout << fDataSet.at(i) << " " << EnergyLow << " " << EnergyUp << endl;
    GetDeltaNL(fDataSet.at(i), EnergyLow, EnergyUp,&delta);
    cout << delta.size() << endl;
    DeltaNL.insert(DeltaNL.end(), delta.begin(), delta.end());
  }
  Double_t mean = 0;
  Int_t n = DeltaNL.size();
  for(size_t i=0;i<DeltaNL.size();i++){
    mean = mean + DeltaNL.at(i);
  }
  mean = mean/n;
  Double_t sigma = 0;
  Double_t mu = 0;
  for(size_t i=0;i<DeltaNL.size();i++){
    sigma = sigma + pow(DeltaNL.at(i)-mean,2);
    mu = mu + pow(DeltaNL.at(i)-mean,4);
  }
  sigma = sigma/n;
  mu = mu/n;
  Double_t sig1 = mean/2;
  Double_t sig2 = sqrt((mu-(n-3)/(n-1)*pow(sigma,2))/n);
  ofstream fout(Form("NL_%d_%d.txt",EnergyLow10,EnergyUp10),ios::app);
  cout << "DS" << EnergyLow << " " << EnergyUp << " " << mean << " " << sigma << " " << mu << " " << n << " "<<sig1 << " " << sig2 << endl;
  fout << "DS" << EnergyLow << " " << EnergyUp << " " << mean << " " << sigma << " " << mu << " " << n <<  " "<< sig1 << " " << sig2 << endl;  
}

