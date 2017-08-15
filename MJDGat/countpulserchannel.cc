#include "MJDGat.hh"
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
  if(argc != 3) {
    cout << "Usage: " << argv[0] << "[dataset] [run]" << endl;
    return 1;
  }
  Int_t fDataSet =atoi(argv[1]);
  Int_t fRun = atoi(argv[2]);

  GATAutoCal ds(fRun,fRun);
  MJDGat ga(fRun);
  //string fDataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> PulserChannel;
  Int_t PulserNum = ga.PulserCount(&PulserChannel);

  ofstream fout(Form("pulser_%d.txt",fDataSet),ios::app);
  fout << fRun << " " << PulserNum << endl;
  fout.close();

  for(size_t i=0;i<fChannel.size();i++){
    Int_t pos = fCryo.at(i)*1000+fStr.at(i)*100+fDetpos.at(i)*10+(fChannel.at(i)%2);
    ofstream fout(Form("pulser_%d_%d.txt",fDataSet,pos),ios::app);
    fout << fRun << " " << pos << " " << fChannel.at(i) << " " << PulserChannel.at(i) << endl;
    fout.close();
  }

  /*
  TFile *pulserfile = new TFile(Form("%spulser_%d.root", fPulserPath.c_str(), fRun), "read");
  TTree *pulsertree = (TTree*)pulserfile->Get("pulsertree");
  Int_t Pulser;
  pulsertree->SetBranchAddress("Pulser", &Pulser);
  Int_t count = 0;
  for(Int_t i=0;i<(Int_t)pulsertree->GetEntries();i++){
    pulsertree->GetEntry(i);
    if(Pulser == 1){
      count++;
    }
  }
  cout << fRun << " " << count << endl;
  */

}
