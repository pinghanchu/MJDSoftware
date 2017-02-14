#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TLine.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TLatex.h"
#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
//#include <unistd.h>
#include <getopt.h>
using namespace std;


int ctob(int c = 3){
  //cout << c <<endl;
  int remain = 0;    
  int b = 0;         
  int count = 0;     
  while(c>0){        
    remain = c%2;    
    b = b+remain*pow(10,count);
    c=c/2;                     
    count++;                   
  }                            
  return b;                    
}    



Double_t fitf(Double_t *v, Double_t *par){
  Double_t arg = 0;                       
  if (par[2] != 0) arg = (v[0] - par[1])/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg) + par[3];
  return fitval;                                             
}                                                            

Double_t linear(Double_t *v, Double_t *par){

  Double_t fitval = -v[0]+par[0] ;
  return fitval;
}
Double_t inverse(Double_t *v, Double_t *par){

  Double_t fitval = par[0]/v[0] ;
  return fitval;
}



Double_t horizontal(Double_t *v, Double_t *par){

  Double_t fitval=par[0];
  return fitval;
}

int main(int argc, char** argv)
{
  if(argc != 2) {
    cout << "Usage: " << argv[0] << " [random ]" << endl;
    return 1;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetFrameBorderMode(0);
  TCanvas *c1 = new TCanvas("c1", "Human-Readable Name");
  c1->SetFillColor(0);  

  const int input =1;
  ifstream fin[input];
  fin[0].open("./data/coupling.12.txt");
  //  fin[1].open("./data/Long.12.txt");
  //fin[2].open("./data/Chu.12.txt");
  //fin[3].open("./data/coupling.12.opt.txt");

  //fin[3].open("./data/NIST.12.txt");

  double r1, r2;
  vector<double> Lambda[input];
  vector<double> Coupling[input];
  for(int i=0;i<input;i++){
    if(fin[i].is_open()){
      while(!fin[i].eof()){
	fin[i] >> r1 >> r2;
	Lambda[i].push_back(r1);
	Coupling[i].push_back(r2);
      }
    }
    Lambda[i].pop_back();
    Coupling[i].pop_back();
  }
  //for(int i=0;i<4;i++){
  //  for(int j=0;j<Lambda[i].size();j++){
  //    cout << Lambda[i].at(j) << " " << Coupling[i].at(j) << endl;
  //  }
  //}

  const int n1 = Lambda[0].size();
  /*
  const int n2 = Lambda[1].size();
  const int n3 = Lambda[2].size();
  const int n4 = Lambda[3].size();
  */
  double X1[n1],Y1[n1];
  /*
  double X2[n2],Y2[n2];
  double X3[n3],Y3[n3];
  double X4[n4],Y4[n4];
  */
  for(int i = 0;i<n1;i++){
    X1[i] = Lambda[0].at(i);
    Y1[i] = Coupling[0].at(i);
  }
  /*
  for(int i = 0;i<n2;i++){
    X2[i] = Lambda[1].at(i);
    Y2[i] = Coupling[1].at(i);
  }
  for(int i = 0;i<n3;i++){
    X3[i] = Lambda[2].at(i);
    Y3[i] = Coupling[2].at(i);
  }

  for(int i = 0;i<n4;i++){
    X4[i] = Lambda[3].at(i);
    Y4[i] = Coupling[3].at(i);
  }
  */
  TGraph *g[input];
  g[0] = new TGraph(n1,X1,Y1);
  //  g[1] = new TGraph(n2,X2,Y2);
  //g[2] = new TGraph(n3,X3,Y3);
  //g[3] = new TGraph(n4,X4,Y4);


  c1->SetLogx();
  TMultiGraph *mg = new TMultiGraph();
  for(int i=0;i<input;i++){
    g[i]->GetXaxis()->SetTitle("lambda(m)");
    g[i]->GetYaxis()->SetTitle("Log(gg)");
    g[i]->GetYaxis()->SetTitleOffset(1.3);
    g[i]->SetMarkerStyle(7);
    g[i]->SetMarkerSize(1);
    g[i]->SetLineWidth(3);
    mg->Add(g[i],"l");
  }
  
  g[0]->SetLineStyle(7);
  g[0]->SetLineWidth(4);
  g[0]->SetLineColor(2);
  /*
  g[1]->SetLineStyle(3);
  g[1]->SetLineColor(3);

  g[2]->SetLineStyle(4);
  g[2]->SetLineColor(4);

  g[3]->SetLineStyle(7);
  g[3]->SetLineColor(46);
  */
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTopMargin(0.05);
  mg->SetTitle(";#lambda(m);Log(f_{12+13})");
  mg->Draw("a");
  mg->GetXaxis()->SetLimits(0.000001,0.2);
  mg->SetMinimum(-33.);
  mg->SetMaximum(-10.);
  mg->GetXaxis()->SetLabelSize(0.05);
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleOffset(1.1);
  mg->GetYaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleOffset(1.1);
 
  TF1 *f1 = new TF1("f1",inverse,1e-6,2e-1,1);
  f1->SetParameter(0,2e-7);
  TGaxis *ax1 = new TGaxis(2e-1,-10,1e-6,-10,"f1",510,"G+");
  ax1->SetTitle("m(eV)");
  ax1->SetLabelSize(0.);
  ax1->SetLabelOffset(-0.04);
  ax1->SetLabelFont(40);
  ax1->SetTitleFont(42);
  ax1->SetTitleSize(0.);
  ax1->SetTitleOffset(1.2);
  ax1->Draw();

  TLatex latex;
  latex.SetTextSize(0.05);
  latex.DrawLatex(1.4e-4,-25,"case 1");

  /*

  //  TLegend *leg = new TLegend(0.5,0.8,0.92,0.95);
  //  TLegend *leg = new TLegend(0.14,0.15,0.6,0.23);
  TLegend *leg = new TLegend(0.15,0.15,0.7,0.23);
  leg->AddEntry(g[0],"10 fT LANL SERF, gap = 1 cm; case 1","l");
  //  leg->AddEntry(g[3],"1 fT/#sqrt{Hz}, gap = 0.1 mm; case 1","l");
  //leg->AddEntry(g[1],"Lesile, et al, PRD89,114022","l");
  //leg->AddEntry(g[2],"Chu, et al, PRD91,102006","l");
  //  leg->AddEntry(g[3],"Yan, et al, PRL110,082003(g_{A}^{n}g_{V}^{N})","l");
  leg->SetTextFont(132);
  leg->SetBorderSize(0);
  leg->Draw();
  */
  c1->Print("sensitivity_12.pdf");
  return 1;
}


