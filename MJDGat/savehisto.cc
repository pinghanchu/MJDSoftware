#include "MJDGat.hh"
#include "GATPulserTag.hh"
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
  if(argc != 5) {
    cout << "Usage: " << argv[0] << " [InputFile Path] [InputFile Name] [Cut] [Input Run File]" << endl;
    return 1;
  }

  string fInputPath = argv[1];
  string fInputName = argv[2];
  string fCut = argv[3];
  string fInputRun = argv[4];
  cout << fInputPath << " " <<fCut <<endl;  
  ifstream fin(Form("%s",fInputRun.c_str()));

  Int_t run1,run2;
  vector<Int_t> Run1;
  vector<Int_t> Run2;
  vector<Int_t> Run;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run1 >> run2;
      Run1.push_back(run1);
      Run2.push_back(run2);
    }
  }
  Run1.pop_back();
  Run2.pop_back();
  for(size_t i =0;i<Run1.size();i++){
    for(Int_t ir = Run1.at(i);ir<=Run2.at(i);ir++){
      Run.push_back(ir);
    }
  }

  MJDGat ds(Run.at(0));
  string FileName;
  //TChain *fMjdTree = new TChain("mjdTree");
  //TCanvas *c1 = new TCanvas("c1");
  //c1->SetLogy();
  //ofstream fout("entries.txt");
  //  for(size_t i=0;i<Run.size();i++){
  TChain *fMjdTree = new TChain("mjdTree");
  for(size_t i=0;i<Run.size();i++){
    cout << Run.at(i) << endl;
    FileName = Form("%s%s%d.root",fInputPath.c_str(),fInputName.c_str(),Run.at(i));
    fMjdTree->Add(Form("%s",FileName.c_str()));
  }
  TH1D *h = new TH1D(Form("%s","trapENFCal"),"", 30000,0,3000);
  fMjdTree->Draw(Form("trapENFCal>>%s","trapENFCal"),Form("%s",fCut.c_str()));
  h->SetTitle(";Energy(keV);Counts");
  //h->Rebin(400);
  //h->GetXaxis()->SetRangeUser(300,400);
  //if(h->GetEffectiveEntries()>10){
  //  fout << Run.at(i) << " " << h->GetEffectiveEntries() <<endl;
  // }
    //h->Draw("Hist");
    //c1->Print(Form("hist_%d_1.pdf",Run.at(i)));
  TFile fhist("hist.root","recreate");
  h->Write();
  fhist.Close();
  //}
}
