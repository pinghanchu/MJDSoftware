// 2016.2.18
// written by Pinghan Chu
// Following the logic in GATAutoCal.cc
// 
#include "GATAutoCal.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <string>
#include <stdlib.h>

Double_t Flat(Double_t *v, Double_t *par){
  Double_t fitval = par[0];
  return fitval;
}

Double_t IsNan(string Input){
  Double_t output = 0;
  if(Input == "nan" || Input == "-nan" || Input == "inf" || Input == "-inf"){
    output = 0;
  }else{
    output = atof(Input.c_str());
  }
  return output;
}

void GetCal(Int_t fStartRun, Int_t fEndRun, string fChannel, string fEName, Int_t fIndex,vector<Double_t>* Par, vector<Double_t>* ParErr){

  //cout << "Read: "<<fStartRun << " " << fEndRun << " "<< fChannel << " " << fEName.c_str() << " " << fIndex << endl;
  ifstream fin2(Form("./List/cal/ezfit_%d_%d.txt",fStartRun,fEndRun));
  string x0,x1,x2,x3,x4;
  Double_t r1,r2,r3,r4;
  Int_t run1,run2;
  string cha;
  Int_t pos,index;
  string ename;
  if(fin2.is_open()){
    while(!fin2.eof()){
      fin2 >> pos>> cha >> run1>>run2>>ename >> index >> x0 >> x1 >> x2 >>x3 >>x4;      
      r1 = IsNan(x1);
      r2 = IsNan(x2);
      r3 = IsNan(x3);
      r4 = IsNan(x4);
      //cout << "Load: " << cha<< " " << ename << " " << index << " " << r1 << " " << r2 << " " << r3 <<  " " << r4 << endl;
      if(cha == fChannel && ename == fEName && index==fIndex){	
      //if(cha == fChannel && index == fIndex){
	//cout << "Match: " << cha<< " " << ename << " " << index << " " << r1 << " " << r2 << " " << r3 <<  " " << r4 << endl;  
	Par->push_back( (fStartRun+fEndRun)/2); // Run
	Par->push_back( r1 ); //peak
	Par->push_back( r3 ); //reso
	ParErr->push_back( (fEndRun-fStartRun)/2 );
	ParErr->push_back( r2 );
	ParErr->push_back( r4 );	
      }
    }
  }
}



TGraphErrors* GetGraph(vector<Double_t> Px,vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr){
  const Int_t nPy = Py.size();

  Double_t A[nPy];
  Double_t AErr[nPy];
  Double_t B[nPy];
  Double_t BErr[nPy];

  for(Int_t j=0;j<nPy;j++){
    A[j]=Px.at(j);
    AErr[j] = PxErr.at(j);
    B[j]=Py.at(j);
    BErr[j] =PyErr.at(j);
  }
  TGraphErrors* fGraph = new TGraphErrors(nPy,A,B,AErr,BErr);
  fGraph->SetFillStyle(0);
  fGraph->SetFillColor(0);


  return fGraph;
}


int main(int argc, char** argv)
{
  if(argc != 5 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [energy name] [peak index]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  const char* EnergyName = argv[3];
  Int_t fIndex = atoi(argv[4]);
  string fEName(EnergyName);
  //gStyle->SetOptStat(1100);
  //gStyle->SetOptFit(1111);
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  ofstream fmiss(Form("./List/miss.ezfit_%d_%d.txt",fStartRun,fEndRun));
  ofstream falign(Form("./List/align.ezfit_%d_%d.txt",fStartRun,fEndRun));
  ifstream fin1("./List/runlist/cal.all.txt");
  vector<Int_t> startrun;
  vector<Int_t> endrun;
  Int_t r1,r2,r3,r4;
  if(fin1.is_open()){
    while(!fin1.eof()){
      fin1 >> r1 >> r2 >>r3 >>r4;
      if(r1>=fStartRun && r2 <=fEndRun){
	//cout << r1 << " " << r2 << endl;
	startrun.push_back(r1);
	endrun.push_back(r2);
      }
    }
  }
  if(startrun.at(startrun.size()-1)==startrun.at(startrun.size()-2)){
    startrun.pop_back();
    endrun.pop_back();
  }
  vector<Int_t> Color;
  Color.push_back(2);
  Color.push_back(3);
  Color.push_back(4);
  Color.push_back(6);
  Color.push_back(7);
  Color.push_back(8);
  Color.push_back(9);
  Color.push_back(12);
  Color.push_back(28);
  Color.push_back(32);
  Color.push_back(37);
  Color.push_back(38);
  Color.push_back(41);
  Color.push_back(44);
  Color.push_back(46);
  Color.push_back(49);

  vector<Int_t> Style;
  Style.push_back(20);
  Style.push_back(24);
  Style.push_back(21);
  Style.push_back(25);
  Style.push_back(22);
  Style.push_back(26);
  Style.push_back(23);
  Style.push_back(32);
  Style.push_back(33);
  Style.push_back(27);
  Style.push_back(34);
  Style.push_back(28);
  Style.push_back(29);
  Style.push_back(30);

  vector<string> DataName;
  //DataName.push_back("00"); // Natural HG
  // DataName.push_back("01"); // Natural LG
  //DataName.push_back("10"); // Enriched HG
  //DataName.push_back("11"); // Enriched LG
  DataName.push_back("HG");
  //DataName.push_back("LG");
  vector<string> TitleName;
  //TitleName.push_back("Natural HG");
  //TitleName.push_back("Natural LG");
  //TitleName.push_back("Enriched HG");
  //TitleName.push_back("Enrichd LG");
  TitleName.push_back("High Gain");
  //TitleName.push_back("Low Gain");


  ///////////////////////////////////////////
  GATAutoCal ds(fStartRun,fStartRun+1);
  ds.SetEnergyName(EnergyName);
  ds.SetParameters();
  const Int_t channels = DataName.size();
  Int_t startdate = 0;
  Int_t enddate = 0;

  string PlotName;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> peak;
  vector<Double_t> peakerr;
  vector<Double_t> fwhm;
  vector<Double_t> fwhmerr;
  vector<Double_t> run;
  vector<Double_t> runerr;
  vector<Int_t> GoodBadIndex;
  TGraphErrors *fPeak[channels];
  TGraphErrors *fPeakDelta[channels];
  TGraphErrors *fFWHM[channels];
  for(Int_t i=0;i<channels;i++){
    ofstream fnew(Form("./List/cal/ezfit_%s.txt",DataName.at(i).c_str()),ios::app);
    fnew.precision(12);
    run.clear();
    runerr.clear();
    peak.clear();
    peakerr.clear();
    fwhm.clear();
    fwhmerr.clear();
    for(Int_t j=0;j<(int)startrun.size();j++){
      Int_t run1 = startrun.at(j);
      Int_t run2 = endrun.at(j);
      Par.clear();
      ParErr.clear();
      GetCal(run1,run2,DataName.at(i),fEName,fIndex, &Par,&ParErr);
      Int_t pars = Par.size();
      if(Par.size()>0){
	run.push_back(Par.at(pars-3));
	peak.push_back(Par.at(pars-2));
	runerr.push_back(ParErr.at(pars-3));
	peakerr.push_back(ParErr.at(pars-2));
	fwhm.push_back(Par.at(pars-1));
	fwhmerr.push_back(ParErr.at(pars-1));
      }else{
	fmiss << run1 << " " << run2 << " " << fEName.c_str() << " " << DataName.at(i).c_str() << endl;
      }
    }
    for(size_t ir = 0;ir<run.size();ir++){      
      Int_t r1 = run.at(ir)-runerr.at(ir);
      Int_t r2 = run.at(ir)+runerr.at(ir);
      GATAutoCal ds1(r1,r2);
      startdate = ds1.GetStartTime();
      enddate =  ds1.GetStopTime();

      fnew << DataName.at(i).c_str() << " " << r1 << " " << r2 << " " << startdate << " "<< enddate << " "
	   << fEName.c_str() << " " << fIndex << " " 
	   << peak.at(ir) << " " << peakerr.at(ir) << " " << fwhm.at(ir) << " "<< fwhmerr.at(ir) << endl;
    }
    vector<Double_t> run1;
    vector<Double_t> runerr1;
    vector<Double_t> peak1;
    vector<Double_t> peakerr1;
    vector<Double_t> fwhm1;
    vector<Double_t> fwhmerr1;
    Double_t peakthresholdHG = 6000;
    Double_t peakthresholdLG = 1700;
    run1.clear();
    runerr1.clear();
    peak1.clear();
    peakerr1.clear();
    fwhm1.clear();
    fwhmerr1.clear();
    string Unit;
    for(size_t irun=0;irun<run.size();irun++){
      Int_t r1 = run.at(irun)-runerr.at(irun);
      Int_t r2 = run.at(irun)+runerr.at(irun);
      
      if(fEName  == "trapENFCal" || fEName == "trapECal" || fEName == "trapENMCal"){
	Unit = "KeV";
	run1.push_back(run.at(irun));
	runerr1.push_back(runerr.at(irun));
	peak1.push_back(peak.at(irun));
	peakerr1.push_back(peakerr.at(irun));	    
	fwhm1.push_back(fwhm.at(irun));
	if(fwhmerr.at(irun)>1){
	  fwhmerr1.push_back(sqrt(fwhm.at(irun)));
	  }else{
	  fwhmerr1.push_back(fwhmerr.at(irun));
	}
	if(peak.at(irun)< 2610 || peak.at(irun)>2620){
	  falign << r1 << " " << r2 << " " << fEName.c_str() << " " << DataName.at(i).c_str() << " " << peak.at(irun)<< endl;	    
	}
      }
    }
    
    if(run1.size()>0){
      TF1 *f = new TF1("f",Flat,fStartRun,fEndRun,1);
      GoodBadIndex.push_back(i);
      fPeak[i] = GetGraph(run1,runerr1,peak1,peakerr1);	
      fPeak[i]->Fit(f);
      Double_t mean = f->GetParameter(0);
      fPeak[i]->SetTitle(Form("%s;Run;Peak(%s)",TitleName.at(i).c_str(),Unit.c_str()));
      fPeak[i]->GetYaxis()->SetTitleOffset(1.5);
      fPeak[i]->SetMarkerStyle(2);
      fPeak[i]->SetMarkerSize(1);
      PlotName = Form("./Plot/peak_%s_%s_%d_%d_%d.pdf", fEName.c_str(),DataName.at(i).c_str(),fStartRun,fEndRun,fIndex);
      ds.PlotGraph(fPeak[i], PlotName);
      
      for(size_t ip = 0;ip<peak1.size();ip++){
	peak1.at(ip) = peak1.at(ip)-mean;
	peakerr1.at(ip)=peakerr1.at(ip);
      }
      fPeakDelta[i] = GetGraph(run1,runerr1,peak1,peakerr1);
      fPeakDelta[i]->SetTitle(Form("%s;Run;#Delta Peak(%s)",TitleName.at(i).c_str(),Unit.c_str()));
      fPeakDelta[i]->GetYaxis()->SetTitleOffset(1.5);
      fPeakDelta[i]->SetMarkerStyle(2);
      fPeakDelta[i]->SetMarkerSize(1);
      PlotName = Form("./Plot/peakdelta_%s_%s_%d_%d_%d.pdf", fEName.c_str(),DataName.at(i).c_str(), fStartRun,fEndRun,fIndex);
      ds.PlotGraph(fPeakDelta[i], PlotName);
      
      fFWHM[i] = GetGraph(run1,runerr1,fwhm1,fwhmerr1);
      
      fFWHM[i]->SetTitle(Form("%s;Run;FWHM(%s)",DataName.at(i).c_str(),Unit.c_str()));
      fFWHM[i]->GetYaxis()->SetTitleOffset(1.5);
      fFWHM[i]->SetMarkerStyle(2);
      fFWHM[i]->SetMarkerSize(1);
      PlotName = Form("./Plot/fwhm_%s_%s_%d_%d_%d.pdf", fEName.c_str(),DataName.at(i).c_str(),fStartRun,fEndRun,fIndex);
      ds.PlotGraph(fFWHM[i], PlotName);
    }
    
  }

}


