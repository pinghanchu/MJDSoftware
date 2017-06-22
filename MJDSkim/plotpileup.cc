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
  if(argc != 4 ) {
    cout << "Usage: " << argv[0] << " [data set] [start subset] [end subset]" << endl;    
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  Int_t fStartSubSet = atoi(argv[2]);
  Int_t fEndSubSet = atoi(argv[3]);

  TChain* fSkimTree = new TChain("skimTree");
  /*
  for(Int_t i=fStartSubSet;i<=fEndSubSet;i++){
    fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%d/GAT-v01-06/skimDS%d_%d.root",fDataSet,fDataSet,i));
  }
  */

  TChain* fPileUpTree = new TChain("pileupTree");
  /*
  for(Int_t i=fStartSubSet;i<=fEndSubSet;i++){
    fPileUpTree->Add(Form("./pileup_%d_%d.root",fDataSet,i));
  }
  */
  string fInputFile = "runlist.txt";
  ifstream fin(Form("%s",fInputFile.c_str()));

  Int_t dataset,subset;
  vector<Int_t> DataSet;
  vector<Int_t> SubSet;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> dataset >> subset;
      DataSet.push_back(dataset);
      SubSet.push_back(subset);
    }
  }
  DataSet.pop_back();
  SubSet.pop_back();

  for(size_t i=0;i<DataSet.size();i++){
    fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%d/GAT-v01-06/skimDS%d_%d.root",DataSet.at(i),DataSet.at(i),SubSet.at(i)));
    fPileUpTree->Add(Form("./pileup_%d_%d.root",DataSet.at(i),SubSet.at(i)));    
  }
  fSkimTree->AddFriend(fPileUpTree);
  /*
  TH2D *h1 = new TH2D("h1","",300,0,3000,160,-8000,8000);
  TH2D *h2 = new TH2D("h2","",300,0,3000,30,0,30);
  TH2D *h3 = new TH2D("h3","",160,-8000,8000,30,0,30);
  TH1D *h4 = new TH1D("h4","",1000,0,0.2);
  TH1D *h5 = new TH1D("h5","",300,0,30);
  TH1D *h6 = new TH1D("h6","",1600,-8000,8000);

  //trapENFCal vs PileUpDeltaT
  TCanvas *c1 = new TCanvas("c1");
  fSkimTree->Draw("PileUpDeltaT:trapENFCal >> h1","IsPileUp>0 && PileUpAE>0.004 && trapENFCal>40 && trapENFCal<12000");
  //c1->Print(Form("trapENFCal_PileUpDeltaT_%d_%d_%d.png",fDataSet,fStartSubSet,fEndSubSet));
  h1->SetTitle(";trapENFCal(keV);#Delta T(ns)");
  h1->Draw("COLZ");
  c1->Print("trapENFCal_PileUpDeltaT.pdf");
  c1->Update();

  //trapENFCal vs PileUpRatio
  fSkimTree->Draw("PileUpRatio:trapENFCal >> h2","IsPileUp>0 && PileUpAE>0.004 && trapENFCal>40 && trapENFCal<12000");
  //c1->Print(Form("trapENFCal_PileUpRatio_%d_%d_%d.png",fDataSet,fStartSubSet,fEndSubSet));
  h2->SetTitle(";trapENFCal(keV);Ratio");
  h2->Draw("COLZ");
  c1->Print("trapENFCal_PileUpRatio.pdf");
  c1->Update();

  //PileUpRatio  vs PileUpDeltaT
  fSkimTree->Draw("PileUpRatio:PileUpDeltaT >> h3","IsPileUp>0 && PileUpAE>0.004 && trapENFCal>40 && trapENFCal<12000");
  //c1->Print(Form("PileUpRatio_PileUpDeltaT_%d_%d_%d.png",fDataSet,fStartSubSet,fEndSubSet));
  h3->SetTitle(";#Delta T(ns);Ratio");
  h3->Draw("COLZ");
  c1->Print("PileUpRatio_PileUpDeltaT.pdf");
  c1->Update();

  //PileUpAE
  c1->SetLogy();
  fSkimTree->Draw("PileUpAE>>h4","trapENFCal>40 && trapENFCal<12000");
  //c1->Print(Form("PileUpAE_%d_%d_%d.png",fDataSet,fStartSubSet,fEndSubSet));
  h4->SetTitle(";A/E;");
  h4->Draw();
  c1->Print("PileUpAE.pdf");
  c1->Update();
  //PileUpRatio
  fSkimTree->Draw("PileUpRatio>>h5","IsPileUp>0 && PileUpAE>0.004 && trapENFCal>40 && trapENFCal<12000");
  //c1->Print(Form("PileUpRatio_%d_%d_%d.png",fDataSet,fStartSubSet,fEndSubSet));
  h5->SetTitle(";Ratio;");
  h5->Draw();
  c1->Print("PileUpRatio.pdf");
  c1->Update();
  //PileUpDeltaT
  fSkimTree->Draw("PileUpDeltaT>>h6","IsPileUp>0 && PileUpAE>0.004 && trapENFCal>40 && trapENFCal<12000");
  //c1->Print(Form("PileUpDeltaT_%d_%d_%d.png",fDataSet,fStartSubSet,fEndSubSet));
  h6->SetTitle(";#Delta T(ns);");
  h6->Draw();
  c1->Print("PileUpDeltaT.pdf");
  c1->Update();

*/

}
