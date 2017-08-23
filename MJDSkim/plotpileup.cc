#include "MJDSkim.hh"
#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TFile.h"
#include "TROOT.h"
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
  TCanvas *c1 = new TCanvas("c1");
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);  
  
  TChain* fSkimTree = new TChain("skimTree");
  vector<Int_t> DataSet;
  vector<Int_t> SubSet;
  for(int i = 0;i<6;i++){
    DataSet.push_back(i);
  }
  SubSet.push_back(75);
  SubSet.push_back(51);
  SubSet.push_back(7);
  SubSet.push_back(24);
  SubSet.push_back(18);
  SubSet.push_back(112);

  string path = "GAT-v01-06-125-gd9332b6";


  if(fStartSubSet != fEndSubSet && fEndSubSet!=0){
    for(Int_t i=fStartSubSet;i<=fEndSubSet;i++){
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%dcal/%s/skimDS%d_run%d_small.root",fDataSet,path.c_str(),fDataSet,i));
    }
  }else{
    for(size_t i=0;i<DataSet.size();i++){
      for(Int_t j=0;j<=SubSet.at(i);j++){
	fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%d/%s/skimDS%d_%d.root",DataSet.at(i),path.c_str(),DataSet.at(i),j));
      }
    }
  }
  /*
  c1->SetLogy();
  TH1D *t1 = new TH1D("t1","",1000,-10,10);
  fSkimTree->Draw("trapETailMin>>t1");
  t1->SetTitle(";trapETailMin;");
  t1->Draw();
  c1->Print("trapETailMin.pdf");

  TH1D *t2 = new TH1D("t2","",600,-0.3,0.3);
  fSkimTree->Draw("trapETailMin/trapENFCal>>t2");
  t2->SetTitle(";trapETailMin/Energy;");
  t2->Draw();
  c1->Print("trapETailMin_Enr.pdf");

  TH1D *t21 = new TH1D("t21","",600,-0.3,0.3);
  fSkimTree->Draw("trapETailMin/trapENFCal>>t21","trapENFCal>48 && trapENFCal<72");
  t21->SetTitle(";trapETailMin/Energy;");
  t21->Draw();
  c1->Print("trapETailMin_Enr_Zoom_In_1.pdf");

  TH1D *t22 = new TH1D("t22","",600,-0.3,0.3);
  fSkimTree->Draw("trapETailMin/trapENFCal>>t22","trapENFCal>200 && trapENFCal<500");
  t22->SetTitle(";trapETailMin/Energy;");
  t22->Draw();
  c1->Print("trapETailMin_Enr_Zoom_In_2.pdf");
  */
  c1->SetLogy(0);
  
  
  c1->SetLogz();
  /*
  TH2D *t3 = new TH2D("t3","",300,0,3000,1200,-2,10);
  fSkimTree->Draw("trapETailMin:trapENFCal>>t3");
  t3->SetTitle(";Energy(keV);trapETailMin");
  t3->Draw();
  c1->Print("Energy_trapETailMin.pdf");
  TH2D *t31 = new TH2D("t31","",200,180,200,400,0,4);
  fSkimTree->Draw("trapETailMin:trapENFCal>>t31");
  t31->SetTitle(";Energy(keV);trapETailMin");
  t31->Draw();
  c1->Print("Energy_trapETailMin_Zoom_In.pdf");
  */
  /*
  TH2D *t4 = new TH2D("t4","",300,0,3000,600,-0.3,0.3);
  fSkimTree->Draw("trapETailMin/trapENFCal:trapENFCal>>t4");
  t4->SetTitle(";Energy(keV);trapETailMin/Enr");
  t4->Draw();
  c1->Print("Energy_trapETailMin_Enr.pdf");
  */
  /*
  TH1D *h1 = new TH1D("h1","",30000,0,30000);
  fSkimTree->Draw("run>>h1","dcr99<0.006 && dcr99>0.004 && trapENFCal< 60 && trapENFCal>40");
  h1->Draw();
  c1->Print("h1.pdf");
  */
  /*
  //////AvsE
  c1->SetLogy();
  TH1D *a1 = new TH1D("a1","AvsE",60000,-3000,3000);
  fSkimTree->Draw("avse>>a1");
  a1->SetTitle(";AvsE;");
  a1->Draw();
  c1->Print("AvsE.pdf");
  c1->Update();
  c1->SetLogy(0);
  c1->SetLogz();

  TH2D *a2 = new TH2D("a2","AvsE vs Enr",300,0,3000,6000,-3000,3000);
  fSkimTree->Draw("avse:trapENFCal>>a2");
  a2->SetTitle(";Energy(keV);AvsE");
  a2->Draw("COLZ");
  c1->Print("Energy_AvsE.pdf");
  c1->Update();
  */
  /*
  TH2D *a3 = new TH2D("a3","AvsE vs Enr",1000,0,100,3500,-500,3000);
  fSkimTree->Draw("avse:trapENFCal>>a3");
  a3->SetTitle(";Energy(keV);AvsE");
  a3->Draw("COLZ");
  c1->Print("Energy_AvsE_Zoom_In.pdf");
  c1->Update();
  */
  /*
  c1->SetLogz();
  
  TH2D *ad1 = new TH2D("ad1","DCR vs AvsE",5000,-2000,3000,400,-0.2,0.2);
  fSkimTree->Draw("dcr99:avse>>ad1");
  ad1->SetTitle(";AvsE;DCR");
  ad1->Draw("COLZ");
  c1->Print("AvsE_DCR.pdf");
  c1->Update();

  TH2D *ad2 = new TH2D("ad2","DCR vs AvsE",1100,-80,30,300,-0.01,0.02);
  //fSkimTree->Draw("dcr99:avse>>ad2","trapENFCal>48 && trapENFCal<72");
  fSkimTree->Draw("dcr99:avse>>ad2");

  ad2->SetTitle(";AvsE;DCR");
  ad2->Draw("COLZ");
  c1->Print("AvsE_DCR_Zoom_In.pdf");
  c1->Update();
  */
  TH2D *ad2 = new TH2D("ad2","DCR vs AvsE",4500,500,5000,600,-0.001,0.005);
  //fSkimTree->Draw("dcr99:avse>>ad2","trapENFCal>48 && trapENFCal<72");
  fSkimTree->Draw("dcr99:avse>>ad2");
  ad2->SetTitle(";AvsE;DCR");
  ad2->Draw("COLZ");
  c1->Print("AvsE_DCR_High.pdf");



  /*
  ////DCR
  c1->SetLogz();
  //TH2D *d2 = new TH2D("d2","DCR vs Enr",240,48,72,300,-0.01,0.02);
  TH2D *d2 = new TH2D("d2","DCR vs Enr",280,200,3000,300,-0.01,0.02);

  fSkimTree->Draw("dcr99:trapENFCal>>d2");
  d2->SetTitle(";Energy(keV);DCR");
  d2->Draw("COLZ");
  c1->Print("Energy_DCR_ROI.pdf");
  c1->Update();

  c1->SetLogz();
  TH2D *d3 = new TH2D("d3","DCR vs Enr",280,200,3000,1000,-0.05,0.05);
  fSkimTree->Draw("dcr99:trapENFCal>>d3");
  d3->SetTitle(";Energy(keV);DCR");
  d3->Draw("COLZ");
  c1->Print("Energy_DCR.pdf");
  c1->Update();

  c1->SetLogy();
  TH1D *d1 = new TH1D("d1","DCR",10000,-0.5,0.5);
  fSkimTree->Draw("dcr99>>d1");
  d1->SetTitle(";DCR;");
  d1->Draw();
  c1->Print("DCR.pdf");
  c1->Update();

  c1->SetLogy();
  TH1D *d4 = new TH1D("d4","DCR",10000,-0.01,0.01);
  fSkimTree->Draw("dcr99>>d4","trapENFCal>48 && trapENFCal<72");
  d4->SetTitle(";DCR;");
  d4->Draw();
  c1->Print("DCR_EnergyCut_ROI.pdf");
  c1->Update();

  c1->SetLogy();
  TH1D *d5 = new TH1D("d5","DCR",10000,-0.01,0.01);
  fSkimTree->Draw("dcr99>>d5","trapENFCal>200 && trapENFCal<500");
  d5->SetTitle(";DCR;");
  d5->Draw();
  c1->Print("DCR_EnergyCut_ROI_2.pdf");
  c1->Update();



  c1->SetLogy(0);
  c1->SetLogz(0);

  */

  ////////////////////////
  /*
  TChain* fPileUpTree = new TChain("pileupTree");
  //TChain* fTimeDiffTree=new TChain("timediffTree");
  if(fStartSubSet != fEndSubSet && fEndSubSet!=0){
    for(Int_t i=fStartSubSet;i<=fEndSubSet;i++){
      fPileUpTree->Add(Form("./pileup_%d_%d.root",fDataSet,i));
      //fTimeDiffTree->Add(Form("./timediff_%d_%d.root",fDataSet,i));
    }
  }else{
  
    for(size_t i=0;i<DataSet.size();i++){
      for(Int_t j=0;j<=SubSet.at(i);j++){
	fPileUpTree->Add(Form("./pileup_%d_%d.root",DataSet.at(i),j));
	//fTimeDiffTree->Add(Form("./timediff_%d_%d.root",DataSet.at(i),j));
      }
    }
  }

  
  fPileUpTree->SetBranchStatus("*",1);

  Int_t fRun;
  Int_t fEvent;
  vector<Int_t>* fChannel = NULL;
  vector<Double_t>* fEnr = NULL;
  vector<Double_t>* fAE = NULL;
  vector<Double_t>* fA = NULL;
  vector<Double_t>* fRatio = NULL;
  vector<Double_t>* fDeltaT = NULL;
  vector<Double_t>* fDCR = NULL;
  //vector<Double_t>* fTimeDiff = NULL;
  fPileUpTree->SetBranchAddress("Run", &fRun);
  fPileUpTree->SetBranchAddress("Event", &fEvent);
  fPileUpTree->SetBranchAddress("Channel",&fChannel);
  fPileUpTree->SetBranchAddress("Energy",&fEnr);
  fPileUpTree->SetBranchAddress("AE",&fAE);
  fPileUpTree->SetBranchAddress("A",&fA);
  fPileUpTree->SetBranchAddress("DCR",&fDCR);
  fPileUpTree->SetBranchAddress("PileUpRatio",&fRatio);
  fPileUpTree->SetBranchAddress("PileUpDeltaT",&fDeltaT);
  //fPileUpTree->SetBranchAddress("TimeDiff",&fTimeDiff);
  ofstream fout("wf.txt");
  for(Int_t i = 0;i<fPileUpTree->GetEntries();i++){
    fPileUpTree->GetEntry(i);
    for(size_t j=0;j<fChannel->size();j++){
      Double_t AE = fAE->at(j);
      Double_t Ratio = fRatio->at(j);
      Double_t DeltaT = fDeltaT->at(j);
      Double_t DCR = fDCR->at(j);
      Double_t DeltaT1 = fTimeDiff->at(j);
      if(abs(fEnr->at(j)-67)<20 && AE>0.004 && Ratio>3 && Ratio<5 && DeltaT>0 && abs(DCR)<0.001){
	fout << fRun << " " << fEvent << " " << fChannel->at(j) << " "<< fEnr->at(j) << " "
	     << Ratio << " " << DeltaT << " " << DeltaT1 << " " << AE << " " << DCR << endl;
	MJDSkim ds(fDataSet,fStartSubSet,fStartSubSet,0);
	TH1D* h = ds.GetWaveform(fRun,fEvent,fChannel->at(j),fEnr->at(j));
        TH1D* h2 = ds.GetHistoDerivative(h,10);
        h->Draw();
        c1->Print(Form("wf_%d_%d_%d.pdf",fRun,fEvent,fChannel->at(j)));
	h2->Draw();
	c1->Print(Form("wf1_%d_%d_%d.pdf",fRun,fEvent,fChannel->at(j)));
	delete h;
	delete h2;
      }
    }
  }

  */
  /*  
  TH2D *h1 = new TH2D("h1","Energy vs Ratio",300,0,3000,30,0,30);
  TH2D *h2 = new TH2D("h2","Energy vs DeltaT",300,0,3000,160,-8000,8000);
  TH2D *h3 = new TH2D("h3","Energy vs DCR",300,0,3000,1000,-0.5,0.5);
  TH2D *h4 = new TH2D("h4","Energy vs AvsE",300,0,3000,6000,-3000,3000);

  TH2D *h5 = new TH2D("h5","DeltaT vs Ratio",160,-8000,8000,30,0,30);
  TH2D *h6 = new TH2D("h6","DeltaT vs Ratio",800,0,8000,4,2,6);

  TH1D *g1 = new TH1D("g1","Ratio",300,0,30);
  TH1D *g2 = new TH1D("g2","DeltaT",1600,-8000,8000);
  TH1D *g3 = new TH1D("g3","DCR",4000,-0.002,0.002);  
  TH1D *g4 = new TH1D("g4","AvsE",6000,-400,200);
*/
  /*
 //Energy vs PileUpRatio
  fPileUpTree->Draw("PileUpRatio:Energy >> h1","Energy>40 && Energy<12000 && Channel%2==0");
  h1->SetTitle(";Energy(keV);Ratio");
  h1->Draw("COLZ");
  c1->Print("Energy_PileUpRatio.pdf");
  c1->Update();

  //Energy vs PileUpDeltaT
  fPileUpTree->Draw("PileUpDeltaT:Energy >> h2","Energy>40 && Energy<12000 && Channel%2==0");
  //c1->Print(Form("Energy_PileUpDeltaT_%d_%d_%d.png",fDataSet,fStartSubSet,fEndSubSet));
  h2->SetTitle(";Energy(keV);#Delta T(ns)");
  h2->Draw("COLZ");
  c1->Print("Energy_PileUpDeltaT.pdf");
  c1->Update();

  //Energy vs DCR
  fPileUpTree->Draw("DCR:Energy >> h3","Energy>40 && Energy<12000 && Channel%2==0");
  h3->SetTitle(";Energy(keV);DCR");
  h3->Draw("COLZ");
  c1->Print("Energy_PileUpDCR.pdf");
  c1->Update();

  //Energy vs AvsE
  fPileUpTree->Draw("AvsE:Energy >> h4","Energy>40 && Energy<12000 && Channel%2==0");
  h4->SetTitle(";Energy(keV);AvsE");
  h4->Draw("COLZ");
  c1->Print("Energy_PileUpAvsE.pdf");
  c1->Update();
  
  //PileUpRatio  vs PileUpDeltaT
  fPileUpTree->Draw("PileUpRatio:PileUpDeltaT >> h5","Energy>40 && Energy<12000 && Channel%2==0");
  h5->SetTitle(";#Delta T(ns);Ratio");
  h5->Draw("COLZ");
  c1->Print("PileUpRatio_PileUpDeltaT.pdf");
  c1->Update();

  //PileUpRatio  vs PileUpDeltaT
  c1->SetLogz(0);
  c1->SetLogy(0);
  fPileUpTree->Draw("PileUpRatio:PileUpDeltaT >> h6","abs(Energy-67)<20 && Channel%2==0 && PileUpRatio<5 && PileUpRatio>3 && abs(DCR)<0.02 && PileUpDeltaT>250");
  h6->SetTitle(";#Delta T(ns);Ratio");
  h6->Draw("COLZ");
  c1->Print("PileUpRatio_PileUpDeltaT_Cut.pdf");
  c1->Update();
*/
  /*
  c1->SetLogy();
  //PileUpRatio
  fPileUpTree->Draw("PileUpRatio>>g1","Energy>40 && Energy<12000 && Channel%2==0");
  g1->SetTitle(";Ratio;");
  g1->Draw();
  c1->Print("PileUpRatio.pdf");
  c1->Update();
  //PileUpDeltaT
  fPileUpTree->Draw("PileUpDeltaT>>g2","Energy>40 && Energy<12000 && Channel%2==0");
  g2->SetTitle(";#Delta T(ns);");
  g2->Draw();
  c1->Print("PileUpDeltaT.pdf");
  c1->Update();

  //DCR
  fPileUpTree->Draw("DCR>>g3","Energy>40 && Energy<12000 && Channel%2==0");
  g3->SetTitle(";DCR;");
  g3->Draw();
  c1->Print("PileUpDCR.pdf");
  c1->Update();


  //AvsE
  fPileUpTree->Draw("AvsE>>g4","Energy>40 && Energy<12000 && Channel%2==0");
  g4->SetTitle(";AvsE;");
  g4->Draw();
  c1->Print("PileUpAvsE.pdf");
  c1->Update();
*/


  /*
  fTimeDiffTree->Draw("TimeDiff>>h9");
  c1->SetLogy();
  h9->SetTitle(";#Delta T(ns);");
  h9->Draw("");
  c1->Print("TimeDiff.pdf");
  */
}
