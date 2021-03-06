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
  if(argc != 4 ) {
    cout << "Usage: " << argv[0] << " [dataset] [energy] [window]" << endl;
    return 1;
  }

  Int_t fDataSet = atoi(argv[1]);
  Int_t fStartSubSet = atoi(argv[2]);
  Int_t fEndSubSet = atoi(argv[3]);

  TChain* fSkimTree = new TChain("skimTree");
  TChain* fPileUpTree = new TChain("pileupTree");
  
  for(Int_t i=fStartSubSet;i<=fEndSubSet;i++){
    fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%d/GAT-v01-06/skimDS%d_%d.root",fDataSet,fDataSet,i));
    fPileUpTree->Add(Form("./pileup_%d_%d.root",fDataSet,i));
  }
  
  /*
  vector<Int_t> DataSet;
  vector<Int_t> SubSet;
  for(int i = 0;i<6;i++){
    DataSet.push_back(i);
  }
  SubSet.push_back(77);
  SubSet.push_back(52);
  SubSet.push_back(8);
  SubSet.push_back(25);
  SubSet.push_back(23);
  SubSet.push_back(113);
  for(size_t i=0;i<DataSet.size();i++){
    for(size_t j=0;j<SubSet.at(i);j++){
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%d/GAT-v01-06/skimDS%d_%d.root",DataSet.at(i),DataSet.at(i),j));
      fPileUpTree->Add(Form("./pileup_%d_%d.root",DataSet.at(i),j));
    }
  }
  */
  //fSkimTree->AddFriend(fPileUpTree);
  fSkimTree->SetBranchStatus("*",1);
  fPileUpTree->SetBranchStatus("*",1);

  Double_t fEnergy = 67.;

  vector<Int_t>* IsPileUp = NULL;
  vector<Double_t>* PileUpRatio = NULL;
  vector<Double_t>* PileUpDeltaT =NULL;
  vector<Double_t>* PileUpAE = NULL;
  fPileUpTree->SetBranchAddress("IsPileUp",&IsPileUp);
  fPileUpTree->SetBranchAddress("PileUpRatio", &PileUpRatio);
  fPileUpTree->SetBranchAddress("PileUpDeltaT", &PileUpDeltaT);
  fPileUpTree->SetBranchAddress("PileUpAE", &PileUpAE);
  Int_t fRun;
  Int_t fEvent;
  vector<Int_t>* fChannel = NULL;
  vector<Int_t>* fP = NULL;
  vector<Int_t>* fD = NULL;
  vector<Int_t>* fC = NULL;
  vector<bool>* fIsEnr = NULL;
  vector<bool>* fIsNat = NULL;
  vector<bool>* fIsGood = NULL;
  vector<Double_t>* fTrapENFCal = NULL;
  fSkimTree->SetBranchAddress("run", &fRun);
  fSkimTree->SetBranchAddress("iEvent", &fEvent);
  fSkimTree->SetBranchAddress("channel",&fChannel);
  fSkimTree->SetBranchAddress("P",&fP);
  fSkimTree->SetBranchAddress("D",&fD);
  fSkimTree->SetBranchAddress("C",&fC);
  fSkimTree->SetBranchAddress("isEnr",&fIsEnr);
  fSkimTree->SetBranchAddress("isNat",&fIsNat);
  fSkimTree->SetBranchAddress("isGood",&fIsGood);
  fSkimTree->SetBranchAddress("trapENFCal", &fTrapENFCal);

  ofstream fout("pileup.txt",ios::app);
  //fSkimTree->GetEntry(0);
  //cout << fSkimTree->GetEntries() << " "<< fPileUpTree->GetEntries() << endl;
  for(Int_t ie=0;ie<fSkimTree->GetEntries();ie++){
    fSkimTree->GetEntry(ie);
    fPileUpTree->GetEntry(ie);
    //cout << ie << " " << fRun << " " << fEvent << " " << fChannel->size() << endl;
    //cout << PileUpAE->size() << " " << fChannel->size() << endl;
    for(size_t j =0;j<fChannel->size();j++){
      if( IsPileUp->at(j)>0 && abs(fTrapENFCal->at(j)-fEnergy)<20. && PileUpRatio->at(j)>2 && PileUpRatio->at(j)<10 ){
	//&& PileUpAE->at(j)>0.002 && PileUpRatio->at(j)>2 && PileUpRatio->at(j)<5){
    	//cout << fRun << " "<< fEvent<< " " << fChannel->at(j) << " " << fTrapENFCal->at(j) << " " << PileUpAE->at(j) << " " <<  PileUpRatio->at(j) << " " << PileUpDeltaT->at(j) << endl;
	fout << fRun << " " << fEvent << " "<< fChannel->at(j) << " " << fTrapENFCal->at(j) << " " << PileUpRatio->at(j) << " "<< PileUpDeltaT->at(j) <<  endl;
      }
    }
  }

  /*
  Int_t fDataSet = atoi(argv[1]);
  Double_t fEnergy = atof(argv[2]);
  //Double_t fTime = atof(argv[3]);
  MJDSkim ds(fDataSet);
  Double_t fQ = 10000;
  Int_t fTimems = (Int_t) 500;
  Double_t fWindow = atof(argv[3]);
  string fOutputFile = Form("data_%d_%d_%d.txt",fDataSet,(Int_t)fEnergy,(Int_t)fWindow);
  cout << "Searching candidates..." << "energy is " << fEnergy << "; the window is " <<fWindow << endl;
  ds.SearchEnergyEvent(fEnergy,fWindow,fOutputFile);
  cout << "Searching candidates is done. "<< endl;

  ifstream fin(Form("%s",fOutputFile.c_str()));

  Int_t run1,list1,entry1,channel1,pos1;
  Double_t enr1,time1,mu_s1;
  vector<Int_t> Run1;
  vector<Int_t> Entry1;
  vector<Int_t> Channel1;
  vector<Double_t> Enr1;
  vector<Double_t> Time1;
  vector<Double_t> Mu_s1;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run1 >> list1 >> entry1 >> pos1 >> channel1 >> enr1 >> time1 >> mu_s1;
      if(enr1> 51){
	Run1.push_back(run1);
	Entry1.push_back(entry1);
	Channel1.push_back(channel1);
	Enr1.push_back(enr1);
	Time1.push_back(time1);
	Mu_s1.push_back(mu_s1);
      }
    }
  }
  if(Run1.size()>0){
    Run1.pop_back();
    Entry1.pop_back();
    Channel1.pop_back();
    Enr1.pop_back();
    Time1.pop_back();
    Mu_s1.pop_back();
  }
  */
  /*
  vector<Int_t> Index;
  Index.push_back(0);
  
  for(size_t i=1;i<Run1.size();i++){
    if(Run1.at(i)!=Run1.at(i-1) || Entry1.at(i)!=Entry1.at(i-1) || Channel1.at(i)!=Channel1.at(i-1)){
      Index.push_back(i);
    }
  }
  */
  /*
  if(Run1.size()>0){    
    //string fOutputWaveform = Form("waveform_%d_%d_%d.root",fDataSet,(Int_t)fEnergy,fTimems);
    //TFile fhist(Form("%s",fOutputWaveform.c_str()),"update");
    ofstream fout(Form("wf_%d_%d_%d.txt",fDataSet,(Int_t)fEnergy,(Int_t)fWindow));
    fout.precision(3);

    vector<Double_t> xp;
    vector<Double_t> yp;
    vector<Double_t> xp1;
    vector<Double_t> yp1;
    vector<Double_t> xp2;
    vector<Double_t> yp2;
    
    Double_t Xmin = 4000;
    Double_t Xmax = 19000;
    if(fDataSet == 2 || fDataSet == 6){
      Xmax = 38000;
    }
    Double_t fResolution = 1;
    Double_t fThreshold = 0.01;
    Double_t fSigma = 5;
    
    for(size_t i=0;i<Run1.size();i++){
      Int_t ii = i;
      cout << Run1.at(ii) << " " << Entry1.at(ii) << " " << Channel1.at(ii) << " " << Enr1.at(ii) << " " << Mu_s1.at(ii) << endl;
      xp.clear();
      yp.clear();
      xp1.clear();
      yp1.clear();
      xp2.clear();
      yp2.clear();
      
      Int_t run = Run1.at(ii);
      Int_t event = Entry1.at(ii);
      Int_t chan = Channel1.at(ii);
      Double_t enr = Enr1.at(ii);
      Double_t time = Time1.at(ii);

      TH1D* h = ds.GetWaveform(run,event,chan,enr);
      if(h!=NULL){
	//h->Write();
	TH1D* h1 = ds.GetHistoSmooth(h,10);
	TH1D* h2 = ds.GetHistoDerivative(h1,10);
	TH1D* h3 = ds.GetHistoSmooth(h2,10);
	TH1D* h4 = ds.GetHistoDerivative(h3,10);
	TH1D* h5 = ds.GetHistoSmooth(h4,10);
	TH1*  hFFT = ds.GetHistoFFT(h);
	
	h3->SetName(Form("waveform1_%d_%d_%d",run,event,chan));
	h5->SetName(Form("waveform2_%d_%d_%d",run,event,chan));
	hFFT->SetName(Form("waveformFFT_%d_%d_%d",run,event,chan));
	//h1->Write();
	h3->GetXaxis()->SetRangeUser(0,20000);
	//h3->Write();
	//h5->Write();
	hFFT->GetXaxis()->SetRangeUser(0,1000);
	//hFFT->Write();
	
	
	// Cut1
	Double_t maxFFT = ds.GetMax(hFFT, 50.,150.); 
	//Double_t maxY0 = ds.GetMax(h1, 100,Xmax);
	Int_t maxBin0 = h1->GetMaximumBin();
	Double_t maxX0 = h1->GetBinCenter(maxBin0);
	
	//cout << max << endl;
	// Cut2 
	Double_t A = ds.GetMax(h3, Xmin,Xmax);
	Double_t AoverE = A/enr;
	
	Int_t nPeak1 = ds.FindPeaks(h3, Xmin, Xmax,fResolution,fSigma, fThreshold, &xp, &yp);
	if(nPeak1>1 && AoverE>0.002){
	  for(Int_t ip = 0;ip<(Int_t)xp.size();ip++){	    
	    Double_t y = ds.GetYValue(h5, xp.at(ip));
	    if(abs(y)< 3e-4 && xp.at(ip)< maxX0){
	      xp1.push_back(xp.at(ip));
	      yp1.push_back(yp.at(ip));
	    }
	  }
	  
	  if(xp1.size()>1){
	    vector<Int_t> Index2 = ds.Sort(yp1);
	    for(size_t ip1 = Index2.size()-2;ip1<Index2.size();ip1++){
	      Int_t ii = Index2.at(ip1);
	      xp2.push_back(xp1.at(ii));
	      yp2.push_back(yp1.at(ii));
	    }
	    
	    fout << run << " " << event << " " << chan << " " <<  enr << " " << time << " " 
		 << yp2.at(1)/yp2.at(0) << " " << xp2.at(0)-xp2.at(1) << " " << AoverE << " "<< maxFFT << endl;	    
	  }
	  
	  delete h;
	  delete h1;
	  delete h2;
	  delete h3;
	  delete h4;
	  delete h5;
	  delete hFFT;
	}
      }
      //      fhist.Close();
    }
  }
*/
}
