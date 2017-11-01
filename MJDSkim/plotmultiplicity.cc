#include "MJDSkim.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << "[input file]" << endl;
    return 1;
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  string fInputFile = argv[1];
  cout << fInputFile.c_str() << endl;
  ifstream fin(Form("%s",fInputFile.c_str()));

  vector<Double_t> DeltaT;
  vector<Double_t> Energy1;
  vector<Double_t> Energy2;
  Int_t dataset,subset,list1,run1,event1,chan1,list2,run2,event2,chan2;
  Double_t time1,enr1,dcr1,avse1,trapetailmin1,mu1,time2,enr2,dcr2,avse2,trapetailmin2,deltaT;
  string mu2;
  vector<Int_t> List;
  vector<Int_t> DataSet;
  vector<Int_t> SubSet;
  vector<Int_t> Channel1;
  vector<Int_t> Channel2;
  vector<Int_t> Channel3;
  vector<Int_t> Energy3;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> dataset >> subset >> 
	list1 >> run1 >> event1 >> time1 >> 
	chan1 >> enr1 >> dcr1 >> avse1 >> trapetailmin1 >> mu1 >> 
	list2 >> run2 >> event2 >> time2 >> 
	chan2 >> enr2 >> dcr2 >> avse2 >> trapetailmin2 >> mu2 >>
	deltaT;
      //if(deltaT>800 && deltaT<1000 && enr1>207 && enr1 < 210){
      if(enr2>5){
	cout << dataset << " " << subset << " " << list1 << " " << chan1 << " " << enr1 << " " << chan2 << " "<< enr2 << " "<< deltaT << endl;
	DeltaT.push_back(deltaT);
	Energy1.push_back(enr1);
	Energy2.push_back(enr2);
	Channel1.push_back(chan1);
	Channel2.push_back(chan2);
	DataSet.push_back(dataset);
	SubSet.push_back(subset);
	List.push_back(list1);
      }
    }
  }

  DeltaT.pop_back();
  Energy1.pop_back();
  Energy2.pop_back();
  Channel1.pop_back();
  Channel2.pop_back();
  DataSet.pop_back();
  SubSet.pop_back();
  List.pop_back();
  vector<Int_t> Channel;
  vector<Double_t> Energy;
  
  TCanvas *c1 = new TCanvas("c1");
  TH2D *h1 = new TH2D("h1","",20,205,225,5000,0,5000); // E1 vs E2
  TH2D *h2 = new TH2D("h2","",20,205,225,100,0,1000); // E1 vs DeltaT
  TH1D *h3 = new TH1D("h3","",20,205,225); // E1
  TH1D *h5 = new TH1D("h5","",1000,0,1000); // Delta T
  TH1D *h4 = new TH1D("h4","",5000,0,5000); // E2
  TH2D *h6 = new TH2D("h6","",20,205,225,5000,0,5000);//E1 vs E3
  TH1D *h7 = new TH1D("h7","",5000,0,5000);//E3

  for(size_t i = 0;i<Energy1.size();i++){
    h1->Fill(Energy1.at(i),Energy2.at(i));
    h2->Fill(Energy1.at(i),DeltaT.at(i));
    h3->Fill(Energy1.at(i));
    h4->Fill(Energy2.at(i));
    h5->Fill(DeltaT.at(i));
  }

  // c1->SetLogx();
  //c1->SetLogy();
  h1->SetTitle(";E_{1}(keV);E_{2}(keV)");
  h1->GetYaxis()->SetTitleOffset(1.3);
  h1->Draw("COLZ");
  c1->Print("77mGe_HighEnr1_HighEnr2.pdf");

  h2->SetTitle(";E_{1}(keV);#Delta T(s)");
  h2->GetYaxis()->SetTitleOffset(1.3);
  h2->Draw("COLZ");
  c1->Print("77mGe_HighEnr1_DeltaT.pdf");


  h3->SetTitle(";E_{1}(keV)");
  h3->GetYaxis()->SetTitleOffset(1.3);
  h3->Draw();
  c1->Print("77mGe_HighEnr1.pdf");
  c1->SetLogy();
  c1->SetLogx();
  h4->SetTitle(";E_{2}(keV)");
  h4->GetYaxis()->SetTitleOffset(1.3);
  h4->Draw();
  c1->Print("77mGe_HighEnr2.pdf");
  c1->SetLogy(0);
  c1->SetLogx(0);
  h5->SetTitle(";#Delta T(s)");
  h5->GetYaxis()->SetTitleOffset(1.3);
  h5->Draw();
  c1->Print("77mGe_DeltaT.pdf");

  Int_t listtemp = -1;
  Int_t subsettemp = -1;
  Int_t datasettemp =-1;

  for(size_t i =0;i<List.size();i++){
    if(DataSet.at(i) == datasettemp && SubSet.at(i) == subsettemp && List.at(i) == listtemp){
    }else{
      cout << DataSet.at(i) << " " << SubSet.at(i) << " " << List.at(i) << endl;
      MJDSkim ds(DataSet.at(i),SubSet.at(i),SubSet.at(i),0);
      Channel.clear();
      Energy.clear();
      ds.FindEvent(List.at(i),&Channel, &Energy);
      for(size_t j=0;j<Channel.size();j++){
        if(Channel.at(j)!=Channel1.at(i) && Channel.at(j)%2==0){
          h6->Fill(Energy1.at(i),Energy.at(j));
	  h7->Fill(Energy.at(j));
        }
      }
    }
    listtemp = List.at(i);
    subsettemp = SubSet.at(i);
    datasettemp = DataSet.at(i);
  }


  h6->SetTitle(";E_{1}(keV);E_{3}(keV)");
  h6->GetYaxis()->SetTitleOffset(1.3);
  h6->Draw("COLZ");
  c1->Print("77mGe_HighEnr1_HighEnr3.pdf");
  h7->SetTitle(";E_{3}(keV);");
  h7->GetYaxis()->SetTitleOffset(1.3);
  h7->Draw();
  c1->Print("77mGe_HighEnr3.pdf");




}
