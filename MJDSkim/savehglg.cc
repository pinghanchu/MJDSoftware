#include "MJDSkim.hh"
#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;
void GetDeltaNL(Int_t DataSet, Double_t EnergyLow, Double_t EnergyUp){

  string fInputFile = Form("/global/projecta/projectdirs/majorana/users/pchu/git/MJDSoftware/MJDSkim/List/runlist/DS%dcal.list.txt",DataSet);
  //cout << fInputFile << endl;
  ifstream fin(Form("%s",fInputFile.c_str()));

  string runfile;
  vector<string> RunFile;

  if(fin.is_open()){
    while(!fin.eof()){
      fin>> runfile;
      RunFile.push_back(runfile);
    }
  }
  RunFile.pop_back();
  string path = "GAT-v01-07";

  for(size_t i = 0;i<RunFile.size();i++){
    string FileName = Form("$MJDDATADIR/surfmjd/analysis/skim/DS%dcal/%s/%s",DataSet,path.c_str(),RunFile.at(i).c_str());
    cout << FileName << endl;
    TChain *skimTree = new TChain("skimTree");
    skimTree->Add(Form("%s",FileName.c_str()));
    skimTree->SetBranchStatus("*",0);
    vector<Int_t>* fChannel = NULL;
    vector<Double_t>* fEnergy = NULL;
    skimTree->SetBranchStatus("channel",1);
    skimTree->SetBranchStatus("trapENFCalC",1);

    skimTree->SetBranchAddress("channel",&fChannel);
    skimTree->SetBranchAddress("trapENFCalC",&fEnergy);
    
    vector<Int_t>* Channel = NULL;
    vector<Double_t>* Energy = NULL;


    TFile *newfile = new TFile(Form("%s",RunFile.at(i).c_str()),"recreate");
    TTree *newtree = new TTree("newtree","");
    newtree->Branch("Channel",&Channel);
    newtree->Branch("Energy",&Energy);
    
    for (Int_t i=0;i<skimTree->GetEntries(); i++) {
      skimTree->GetEntry(i);
      for(size_t j = 0;j<fChannel->size();j++){
	if (fEnergy->at(j) > EnergyLow && fEnergy->at(j) < EnergyUp){
	  //newtree->Fill();
	  Channel->push_back(fChannel->at(j));
	  Energy->push_back(fEnergy->at(j));
	}
      }
      newtree->Fill();
    }
    newtree->Write();
    delete newfile;    
  }
}

int main(int argc, char** argv)
{
  if(argc != 2) {
    cout << "Usage: " << argv[0] << " [enter p]" << endl;
    return 1;
  }
  Int_t DataSet = atoi(argv[1]);
  GetDeltaNL(DataSet,2037,2041);
  
}

