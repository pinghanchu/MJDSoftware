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
  TCanvas *c1 = new TCanvas("c1");
  c1->SetLogy();
  ofstream fout("pulser.txt",ios::app);
  for(size_t i=0;i<Run.size();i++){
    TChain *fMjdTree = new TChain("mjdTree");
    FileName = Form("%s%s%d.root",fInputPath.c_str(),fInputName.c_str(),Run.at(i));
    fMjdTree->Add(Form("%s",FileName.c_str()));
    TH1D *h = new TH1D(Form("%s","trapE"),"", 80000,0,8000);
    fMjdTree->Draw(Form("trapE>>%s","trapE"),Form("%s",fCut.c_str()));
    h->SetTitle(";Energy(ADC);Counts");
    h->Rebin(10);
    h->GetXaxis()->SetRangeUser(50,6000);
    Double_t maxbin = h->GetMaximumBin();
    Double_t ymax = h->GetMaximum();
    Double_t xmax = h->GetBinCenter(maxbin);
    if(ymax>5){
      fout << Run.at(i) << " " << xmax << " " << ymax << " " << xmax << " " << ymax <<endl;
    }else{
      fout << Run.at(i) << " " << 0 << " " << 0 << " " << xmax << " " << ymax <<endl;
    }
  }
}
