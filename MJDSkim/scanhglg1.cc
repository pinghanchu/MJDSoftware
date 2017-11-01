#include "MJDSkim.hh"
#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;
void GetDeltaNL(Int_t DataSet, Int_t Run, Double_t EnergyLow, Double_t EnergyUp, string fOutputFile){


  string path = "GAT-v01-07";

  vector<Double_t> Delta;
  //MJDSkim ps(DataSet,,1);
  //vector<Int_t> channel = ps.GetChannel();
  ofstream fout(Form("%s",fOutputFile.c_str()),ios::app);
  Delta.clear();
  TChain *skimTree = new TChain("skimTree");
  //skimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%dcal/%s/skimDS%d_run%d_small.root",DataSet,path.c_str(),DataSet,Run));
  skimTree->Add(Form("/global/projecta/projectdirs/majorana/users/bxyzhu/cal-skim/skimDS%d_run%d_low.root",DataSet,Run));
  skimTree->SetBranchStatus("*",1);
  vector<Int_t>* fChannel = NULL;
  vector<Double_t>* fEnergy = NULL;
  skimTree->SetBranchAddress("channel",&fChannel);
  //skimTree->SetBranchAddress("trapENFCalC",&fEnergy);
  skimTree->SetBranchAddress("trapECal",&fEnergy);
		
  vector<Int_t> hgchan;
  vector<Double_t> hgenr;
  vector<Int_t> lgchan;
  vector<Double_t> lgenr;
  
  for (Int_t i=0;i<skimTree->GetEntries(); i++) {
    skimTree->GetEntry(i);
    hgchan.clear();
    lgchan.clear();
    hgenr.clear();
    lgenr.clear();
    for(size_t j = 0;j<fChannel->size();j++){
      Int_t chan = fChannel->at(j);
      Double_t enr = fEnergy->at(j);
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
	      
	      if(abs(del)<5){
		Delta.push_back(lgenr.at(k2)-hgenr.at(k1));
	      }
	    }
	  }
	}
      }		
    }
  }
  Int_t n = Delta.size();
  if(n>0){
    Double_t mean = 0;
    for(size_t i=0;i<Delta.size();i++){
      mean = mean + Delta.at(i);
    }
    mean = mean/n;
    fout << DataSet << " " << Run << " " << mean << " " << n << endl;
  }else{
    fout << DataSet << " " << Run << " " << 0 << " " << n << endl;
  }
}

int main(int argc, char** argv)
{
  if(argc != 5) {
    cout << "Usage: " << argv[0] << "[dataset] [run] [low] [up]" << endl;
    return 1;
  }
  Int_t DataSet = atoi(argv[1]);
  Int_t Run = atoi(argv[2]);
  Int_t EnergyLow10 = atoi(argv[3]);
  Int_t EnergyUp10 = atoi(argv[4]);
  Double_t EnergyLow = (Double_t)EnergyLow10/1;
  Double_t EnergyUp = (Double_t)EnergyUp10/1;

  string OutputFile = Form("NL_%d_%d_%d.txt",DataSet,EnergyLow10,EnergyUp10);
  GetDeltaNL(DataSet, Run,EnergyLow, EnergyUp,OutputFile);
  cout << DataSet << " " << Run << " " << EnergyLow10 << " " << EnergyUp10 << " is done! " << endl;
}

