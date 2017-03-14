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
#include "TFile.h"
#include "TSpectrum.h"
#include <string>
#include <stdlib.h>

Double_t Gaus0(Double_t *v, Double_t *par){
  Double_t arg = 0;
  if (par[2] != 0) arg = (v[0] - par[1])/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}

Double_t Flat(Double_t *v, Double_t *par){
  Double_t fitval = par[0];
  return fitval;
}

Double_t FinMaxX(vector<Double_t>* fPx, vector<Double_t>* fPy){
  Double_t temp = 0;
  Double_t maxX = 0;
  for(size_t i = 0;i<fPx->size();i++){
    if(fPy->at(i) > temp){
      maxX = fPx->at(i);
      temp = fPy->at(i);
    }
  }
  return maxX;
}

void FindPeak(TH1D *h, vector<Double_t>* fPx, vector<Double_t>* fPy){
  TSpectrum *s = new TSpectrum(1000,1);
  Int_t nfound = s->Search(h,10,"new",0.0001);
  Double_t *xpeaks = s->GetPositionX();
  Double_t *ypeaks = s->GetPositionY();
  for(Int_t i = 0;i<nfound; i++){
    if(ypeaks[i]>3){
      fPx->push_back(xpeaks[i]);
      fPy->push_back(ypeaks[i]);
    }
  }

  const Int_t npeak = fPx->size();
  Double_t Px[npeak];
  Double_t Py[npeak];
  Int_t n[npeak];
  for(Int_t i = 0;i<npeak;i++){
    Px[i] = fPx->at(i);
    Py[i] = fPy->at(i);
  }
  TMath::Sort(npeak,Px,n,0);
  fPx->clear();
  fPy->clear();

  for(Int_t i = 0;i<npeak;i++){
    Int_t ii = n[i];
    fPx->push_back(Px[ii]);
    fPy->push_back(Py[ii]);
  }
}


TH1D* GetHisto(string FileName, string HistoName){
  TFile *fhist = TFile::Open(Form("%s",FileName.c_str()),"read");
  TH1D *h = (TH1D*)fhist->Get(Form("%s",HistoName.c_str()));
  return h;
}

void FitPulser(Int_t fRun, string fPos, Int_t fChannel, string fEName, vector<Double_t>*Par){

  GATAutoCal ds(fRun,fRun);
  ofstream fout(Form("./List/pulser/pulser_%d.txt",fRun),ios::app);
  string FileName = Form("./hist_%d.root",fRun);
  string HistoName = Form("%s%s%d",fEName.c_str(),fPos.c_str(),fChannel);
  TH1D *h1 = GetHisto(FileName,HistoName);
  if(fRun<60000000){
    h1->GetXaxis()->SetRangeUser(20,3000);
  }else{
    h1->GetXaxis()->SetRangeUser(110,3000);
  }
  Int_t maxbin = h1->GetMaximumBin();
  Double_t maxX = h1->GetBinCenter(maxbin);
  h1->GetXaxis()->SetRangeUser(maxX-5,maxX+5);
  Int_t entries = h1->GetEffectiveEntries();
  vector<Double_t> par;
  vector<Double_t> parerr;
  if(entries>30){
    string PlotName = Form("%s_%d_%s_%d",fEName.c_str(),fRun,fPos.c_str(),fChannel);
    //TitleName = Form("C%dP%dD%d, channel %d",fCryo.at(ii),fStr.at(ii), fDetpos.at(ii), fChannel.at(ii));
    par.clear();
    parerr.clear();
    ds.SkewGaussFit(h1,maxX,5,&par,&parerr,"",PlotName,"","");
    if(parerr.at(1)<1){
      Par->push_back(par.at(0));
      Par->push_back(parerr.at(0));
      Par->push_back(par.at(1));
      Par->push_back(parerr.at(1));
      Par->push_back(par.at(4));
      Par->push_back(parerr.at(4));
    }else{
      Par->push_back(0);
      Par->push_back(0);
      Par->push_back(0);
      Par->push_back(0);
      Par->push_back(0);
      Par->push_back(0);
    }
  }else{
    Par->push_back(0);
    Par->push_back(0);
    Par->push_back(0);
    Par->push_back(0);
    Par->push_back(0);
    Par->push_back(0);
  }
}

void GetPulser(Int_t fRun, Int_t fChannel, string fEName, vector<Double_t>* Par){
  ifstream fin(Form("./List/pulser/pulser_%d.txt",fRun));

  Int_t cha,pos;
  //string R[2];
  Double_t r[2];
  Int_t run;
  string ename;
  if(fin.is_open()){
    while(!fin.eof()){
      //      fin >> pos >> cha >> run >> ename >> R[0] >> R[1] >> R[2] >>R[3] >>R[4] >> R[5];
      fin >> pos >> cha >> run >> ename >> r[0] >> r[1] ;
      if(cha == fChannel && ename ==fEName){
	Par->push_back( (Double_t)run); // Run
	Par->push_back( r[0] ); //Pulser Mu
	Par->push_back( r[1] );
      }
      /*
	fin >> pos >> cha >> run >> ename >> R[0] >> R[1] >> R[2] >>R[3] >>R[4] >> R[5];
      for(Int_t  i = 0;i<2;i++){
        if(R[i]=="nan" || R[i] =="-nan" || R[i] == "inf" || R[i] == "-inf"){
          r[i] = 0;
        }else{
          r[i] = atof(R[i].c_str());
        }
      }
      
      if(cha == fChannel && ename == fEName.c_str()){
	if(r[3] == 0){
	  r[2] = 0;
	}
	if(r[2]>8000){
	  r[2] =0;
	  r[3] =0;
	}
	Par->push_back( (Double_t)run); // Run
	Par->push_back( r[2] ); //Pulser Mu
	Par->push_back( r[3] );
      }
      */

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
  if(argc != 4 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [energy name]" << endl;
    return 1;
  }
  int fStartRun = atoi(argv[1]);
  int fEndRun = atoi(argv[2]);
  const char* EnergyName = argv[3];
  string fEName(EnergyName);
  //gStyle->SetOptStat(1100);
  //gStyle->SetOptFit(1111);
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  TCanvas *c1 = new TCanvas("c1");
  ofstream fmiss(Form("./List/miss.pulser_%d_%d.txt",fStartRun,fEndRun));
  ofstream falign(Form("./List/align.pulser_%d_%d.txt",fStartRun,fEndRun));
  ifstream fin1("./List/runlist/bk.all.txt");
  vector<Int_t> startrun;
  vector<Int_t> endrun;
  Int_t r1,r2;
  if(fin1.is_open()){
    while(!fin1.eof()){
      fin1 >> r1 >> r2 ;
      if(r1>= fStartRun && r2<= fEndRun){
	startrun.push_back(r1);
	endrun.push_back(r2);
	//	cout << r1 <<" " << r2 << endl;
      }
    }
  }
  if(startrun.at(startrun.size()-1) == startrun.at(startrun.size()-2)){
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


  GATAutoCal ds(fStartRun,fStartRun);
  ds.SetEnergyName(EnergyName);
  //  ds.SetParameters();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos= ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  const Int_t channels = fChannel.size();
  Double_t maxX;
  string Pos;
  string TitleName;
  string PlotName;
  vector<Double_t> fPx;
  vector<Double_t> fPy;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> peak;
  vector<Double_t> peakerr;
  vector<Double_t> run;
  vector<Double_t> runerr;
  vector<Int_t> GoodBadIndex;
  TH1D *fPulserHisto[channels];
  TH1D *fPulserHistoDelta[channels];
  TGraphErrors *fPeak[channels];
  TGraphErrors *fPeakDelta[channels];
  for(Int_t i=0;i<channels;i++){
    Pos = Form("%d%d%d",fCryo.at(i),fStr.at(i),fDetpos.at(i));
    //cout << Pos.c_str() << " " << fChannel.at(i) << endl;
    fPulserHisto[i] = new TH1D(Form("pulser%s%d",Pos.c_str(),fChannel.at(i)),"",200000,0,2000);
    fPulserHistoDelta[i] = new TH1D(Form("pulserdelta%s%d",Pos.c_str(),fChannel.at(i)),"",400000,-2000,2000);
    if(fGoodBad.at(i)>0){
      run.clear();
      runerr.clear();
      peak.clear();
      peakerr.clear();

      for(Int_t j=0;j<(int)startrun.size();j++){
	for(Int_t k=startrun.at(j);k<=endrun.at(j);k++){
	  Int_t irun = k;
	  Par.clear();        
	  GetPulser(irun,fChannel.at(i),fEName,&Par);
	  //cout << k << " " <<fChannel.at(i) << " "<<  Par.size() << endl;
	  Int_t pars = Par.size();
	  
	  if(Par.size()>0){
	    if(Par.at(pars-2)>0){
	      run.push_back(Par.at(pars-3));
	      peak.push_back(Par.at(pars-2));
	      if(Par.at(pars-1)<1){	      
		peakerr.push_back(Par.at(pars-1));
	      }else{
		peakerr.push_back(Par.at(pars-1));
		//falign << k << " " << Pos.c_str() << " " << fChannel.at(i) << " " << Par.at(pars-2) << " "<< Par.at(pars-1) << endl;
	      } 
	      runerr.push_back(0);
	    }else{
	      fmiss << k << " " << Pos.c_str() << " " << fChannel.at(i) << endl;
	    }
	  }else{
	    fmiss << k << " " << Pos.c_str() << " " << fChannel.at(i) << endl;
	  }	  
	}
      }


      for(size_t irun=0;irun<run.size();irun++){	
	Double_t mu = peak.at(irun);
	Double_t muerr = peakerr.at(irun);
	Double_t mu1 = mu;
	Double_t mu2 = mu;
	fPulserHisto[i]->Fill(mu);
	if(irun > 0){	    
	    mu1 = peak.at(irun-1);
	}
	if(irun<run.size()-1){
	  mu2 = peak.at(irun+1);
	}
	Par.clear();
	if(((abs(mu-mu1)>1 || abs(mu-mu2)>1) && mu1 >0 && mu2>0&& mu>0) || muerr>1 ){
	  falign << Pos.c_str() << " " << fChannel.at(i) << " " << (Int_t)run.at(irun) << " " << mu << " " << mu1 << " " << mu2 << " " << muerr << endl;
	  }else if( mu == -99999){
	  falign << Pos.c_str() << " " << fChannel.at(i) << " " << (Int_t)run.at(irun) << " " << mu << " " << mu1 << " " << mu2 << " " << muerr << endl;
	}
      }

      
      //cout << run.size() << " "<< peak.size()<< " "<< peakerr.size() <<endl;
      if(run.size()>0){
	GoodBadIndex.push_back(i);	//TF1 *fun = new TF1("fun",Flat,fStartRun,fEndRun,1);
	fPx.clear();
	fPy.clear();
	FindPeak(fPulserHisto[i], &fPx, &fPy);
	maxX = FinMaxX(&fPx,&fPx);		
	//TF1 *f1 = new TF1("f1",Gaus0,maxX-10,maxX+10,3);       
	//f1->SetParameter(1,maxX);
	fPulserHisto[i]->GetXaxis()->SetRangeUser(maxX-3,maxX+3);
	//fPulserHisto[i]->Fit(f1);
	//maxX = f1->GetParameter(1);
	maxX = fPulserHisto[i]->GetMean();
	cout << maxX << endl;
	TF1 *f2 = new TF1("f2",Flat,fStartRun,fEndRun,1);	
	fPeak[i] = GetGraph(run,runerr,peak,peakerr);
	f2->SetParameter(0,maxX);
	fPeak[i]->Fit(f2);
	//Double_t mean = f->GetParameter(0);
	Int_t ic = fCryo.at(i);
	Int_t is = fStr.at(i);
	Int_t id = fDetpos.at(i);
	//cout << ic << is <<id << " "<< fChannel.at(i) << endl;
	fPeak[i]->SetTitle(Form("C%dP%dD%d, Channel %d;Run;Pulaser(ADC)",ic,is,id,fChannel.at(i)));
	fPeak[i]->GetYaxis()->SetTitleOffset(1.5);
	fPeak[i]->SetMarkerStyle(2);
	fPeak[i]->SetMarkerSize(1);
	PlotName = Form("./Plot/pulser_%s_%s_%d_%d_%d.pdf", fEName.c_str(),Pos.c_str(),fChannel.at(i), fStartRun,fEndRun);
	ds.PlotGraph(fPeak[i], PlotName);
	for(size_t ip = 0;ip<peak.size();ip++){
          peak.at(ip) = peak.at(ip)-maxX;
	  fPulserHistoDelta[i]->Fill(peak.at(ip));
        }
        fPeakDelta[i] = GetGraph(run,runerr,peak,peakerr);
        fPeakDelta[i]->SetTitle(Form("C%dP%dD%d, channel %d;run;peak(ADC)",ic,is,id,fChannel.at(i)));
        fPeakDelta[i]->GetYaxis()->SetTitleOffset(1.5);
        fPeakDelta[i]->SetMarkerStyle(2);
        fPeakDelta[i]->SetMarkerSize(1);
        PlotName = Form("./Plot/pulserdelta_%s_%s_%d_%d_%d.pdf", fEName.c_str(),Pos.c_str(),fChannel.at(i), fStartRun,fEndRun);
        ds.PlotGraph(fPeakDelta[i], PlotName);
      }
    }
  }

  TMultiGraph *mg1 = new TMultiGraph();
  TMultiGraph *mg2 = new TMultiGraph();

  TLegend *leg1 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg2 = new TLegend(0.15,0.55,0.85,0.85);

  gStyle->SetLegendBorderSize(0);

  TFile fhist(Form("./pulser_%d_%d.root",fStartRun,fEndRun),"recreate");
  TH1D *pulserdeltaHG = new TH1D("pulserdeltaHG","",400000,-2000,2000);
  TH1D *pulserdeltaLG = new TH1D("pulserdeltaLG","",400000,-2000,2000);
  for(size_t i=0;i<GoodBadIndex.size();i++){
    Int_t j = GoodBadIndex.at(i);
    fPulserHisto[j]->Write();
    fPulserHistoDelta[j]->Write();

    Int_t ic = fCryo.at(j);
    Int_t is = fStr.at(j);
    Int_t id = fDetpos.at(j);
    Int_t hg = id*2;
    Int_t lg = hg+1;
    Int_t chan = fChannel.at(j);
    if(chan%2==0){
      pulserdeltaHG->Add(fPulserHistoDelta[j]);
    }else{
      pulserdeltaLG->Add(fPulserHistoDelta[j]);
    }
    fPeak[j]->SetMarkerSize(1);
    fPeak[j]->SetMarkerColor(Color.at(is));
    fPeak[j]->SetLineColor(Color.at(is));

    if(chan%2==0){
      fPeak[j]->SetMarkerStyle(Style.at(hg));
      fPeak[j]->SetLineStyle(0);
    }else{
      fPeak[j]->SetMarkerStyle(Style.at(lg));
      fPeak[j]->SetLineStyle(2);
    }

    //Titlename = Form("C%dP%dD%d, channel %d",ic,is,id,chan);
    TitleName = Form("#splitline{C%dP%dD%d}{%d}",ic,is,id,chan);
    //cout << j<<""<<TitleName.c_str() << endl;
    if(chan%2 == 0){
      mg1->Add(fPeak[j],"ep");
      leg1->AddEntry(fPeak[j],Form("%s",TitleName.c_str()));
    }else{
      mg2->Add(fPeak[j],"ep");
      leg2->AddEntry(fPeak[j],Form("%s",TitleName.c_str()));
    }
  }
  pulserdeltaHG->Write();
  pulserdeltaLG->Write();
  fhist.Close();
  gStyle->SetLegendTextSize(0.02);
  leg1->SetFillStyle(0);
  leg1->SetNColumns(6);
  leg1->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetNColumns(6);
  leg2->SetBorderSize(0);

  /*
  mg1->SetMaximum(2000);
  mg2->SetMaximum(2000);
  mg1->SetMinimum(0);
  mg2->SetMinimum(0);
  */
  PlotName = Form("./Plot/pulserhg_%s_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun);
  ds.PlotMultiGraph(mg1, leg1, "", PlotName,"");

  PlotName = Form("./Plot/pulserlg_%s_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun);
  ds.PlotMultiGraph(mg2, leg2, "", PlotName,"");

  TMultiGraph *mg3 = new TMultiGraph();
  TMultiGraph *mg4 = new TMultiGraph();

  TLegend *leg3 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg4 = new TLegend(0.15,0.55,0.85,0.85);

  gStyle->SetLegendBorderSize(0);

  for(size_t i=0;i<GoodBadIndex.size();i++){
    Int_t j = GoodBadIndex.at(i);
    Int_t ic = fCryo.at(j);
    Int_t is = fStr.at(j);
    Int_t id = fDetpos.at(j);
    Int_t hg = id*2;
    Int_t lg = hg+1;
    Int_t chan = fChannel.at(j);
    fPeakDelta[j]->SetMarkerSize(1);
    fPeakDelta[j]->SetMarkerColor(Color.at(is));
    fPeakDelta[j]->SetLineColor(Color.at(is));

    if(chan%2==0){
      fPeakDelta[j]->SetMarkerStyle(Style.at(hg));
      fPeakDelta[j]->SetLineStyle(0);
    }else{
      fPeakDelta[j]->SetMarkerStyle(Style.at(lg));
      fPeakDelta[j]->SetLineStyle(2);
    }

    //Titlename = Form("C%dP%dD%d, channel %d",ic,is,id,chan);
    TitleName = Form("#splitline{C%dP%dD%d}{%d}",ic,is,id,chan);
    //cout << j<<""<<TitleName.c_str() << endl;
    if(chan%2 == 0){
      mg3->Add(fPeakDelta[j],"ep");
      leg3->AddEntry(fPeakDelta[j],Form("%s",TitleName.c_str()));
    }else{
      mg4->Add(fPeakDelta[j],"ep");
      leg4->AddEntry(fPeakDelta[j],Form("%s",TitleName.c_str()));
    }
  }
  gStyle->SetLegendTextSize(0.02);
  leg3->SetFillStyle(0);
  leg3->SetNColumns(6);
  leg3->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetNColumns(6);
  leg4->SetBorderSize(0);


  mg3->SetMaximum(1);
  mg4->SetMaximum(1);
  mg3->SetMinimum(-0.5);
  mg4->SetMinimum(-0.2);


  PlotName = Form("./Plot/pulserdeltahg_%s_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun);
  ds.PlotMultiGraph(mg3, leg3, "", PlotName,"");

  PlotName = Form("./Plot/pulserdeltalg_%s_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun);
  ds.PlotMultiGraph(mg4, leg4, "", PlotName,"");

}


