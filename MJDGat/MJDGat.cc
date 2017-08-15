#include "MJDGat.hh"
#include "GATAutoCal.hh"
#include "MJTRun.hh"
#include "MJAnalysisDoc.hh"
#include "MJAnalysisParameters.hh"
#include "GATWaveformBrowser.hh"
#include "MGTWaveform.hh"
#include "GATPeakShape.hh"
#include "GATMultiPeakFitter.hh"
#include "GATHybridMonteCarlo.hh"
#include "MJProvenance.hh"
#include "MJDatabase.hh"
#include "MJDBUtilities.hh"
#include "KTreeFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include "TMatrixDSym.h"
#include "TPaveStats.h"
#include "Fit/FitResult.h"
#include "Fit/FitConfig.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TVirtualFFT.h"
#include <iostream>
#include <bitset>

//ClassImp(MJDGat)

using namespace std;
using namespace katrin;
using namespace MJDB;
//using namespace GATPeakShapeFunction;

////////////////////////////////////////////////////////

//Set ChannelMap, mjdTree (from GATDataSet), CalibrationPeak, Initial Parameters
MJDGat::MJDGat(Int_t Run) : fDS(Run,Run)
{
  //fMjdTree = fDS.GetMJDTree();
  if(Run>45000000 && Run <50000000){
    fMjdTree = new TChain("mjdTree");
    fMjdTree->Add(Form("/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/GatData/mjd_run%d.root",Run));
  }else{
    fMjdTree = fDS.GetMJDTree();
  }
  fMap = fDS.GetMap();
  //fEntries = fDS.GetEntries();
  fRun = Run;
  //fGATRev = fDS.GetGATRev();
  fIsRadio = fDS.GetIsRadio();
  fMTStartTime = fDS.GetStartTime();
  fMTStopTime = fDS.GetStopTime();
  fDataSet = fDS.GetDataSet();
  fChannel = fDS.GetChannel();;
  fPulserTagChannel = fDS.GetPulserChannel();
  fCryo = fDS.GetCryo();
  fString = fDS.GetString();
  fDetector = fDS.GetDetPosition();
  fGoodBad = fDS.GetGoodBad();
  fEnriched = fDS.GetEnriched();
  fDetectorName = fDS.GetDetectorName();

  //fMjdTree->SetBranchStatus("*",1);
  fMjdTree->SetBranchStatus("*",0);
  //fMjdTree->SetBranchStatus("run",1);
  //fMjdTree->SetBranchStatus("C",1);
  //fMjdTree->SetBranchStatus("P",1);
  //fMjdTree->SetBranchStatus("D",1);
  fMjdTree->SetBranchStatus("channel", 1);
  fMjdTree->SetBranchStatus("trapE", 1);
  fMjdTree->SetBranchStatus("trapENFCal",1);
  fMjdTree->SetBranchStatus("timeinfo", 1);
  fMjdTree->SetBranchStatus("globalTime", 1);
  fMjdTree->SetBranchStatus("localTime", 1);
  fMjdTree->SetBranchStatus("clockTime", 1);
  fMjdTree->SetBranchStatus("tOffset", 1);
  fMjdTree->SetBranchStatus("mH",1);
  //fMjdTree->SetBranchStatus("mL",1);
  if(Run<45000000 || Run >50000000){
    fMjdTree->SetBranchStatus("EventDC1Bits",1);
    //fMjdTree->SetBranchStatus("fastTrapNLCWFsnRisingX",1);
  }
  fMTChannel = NULL;
  fMTTrapE = NULL;
  fMTTrapENFCal = NULL;
  fMTTimeInfo = NULL;
  //fMTfastTrapNLCWFsnRisingX = NULL;

  fMjdTree->SetBranchAddress("channel", &fMTChannel);
  fMjdTree->SetBranchAddress("trapE", &fMTTrapE);
  fMjdTree->SetBranchAddress("trapENFCal", &fMTTrapENFCal);
  fMjdTree->SetBranchAddress("timeinfo",&fMTTimeInfo);
  fMjdTree->SetBranchAddress("mH",&fMTmH);
  if(Run<45000000 || Run >50000000){
    fMjdTree->SetBranchAddress("EventDC1Bits", &fMTEventDC1Bits);
    //fMjdTree->SetBranchAddress("fastTrapNLCWFsnRisingX", &fMTfastTrapNLCWFsnRisingX);
  }

  fMjdTree->GetEntry(0);
  fEntries = fMjdTree->GetEntries();
}


void MJDGat::SeachPulserChannel(string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the second gamma energy
  //fEnr2 : the first Q value
  //////////////////////////////////////
  //map<int,int> detid;
  //for(size_t i=0; i<fChannel.size(); i++) detid[fChannel.at(i)]= i;

  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout.precision(15);
  Int_t pulserhit = 0;
  Int_t pulserevent = 0;
  Int_t hitnum = 0;
  Double_t time;
  Int_t IsPulser = 0;
  Int_t channel = 0;
  for(size_t i=0;i<fEntries;i++){
    fMjdTree->GetEntry(i);
    IsPulser = fMTEventDC1Bits;    
    if(IsPulser>0){
      //time1 = fMTTimeInfo->localTime;
      time = fMTTimeInfo->clockTime/10;
      //hitnum = fMTChannel->size();
      hitnum = 0;
      for(size_t j = 0;j<fMTChannel->size();j++){
	channel = fMTChannel->at(j);
	for(size_t k = 0;k < fChannel.size(); k++){
	  if(channel == fChannel.at(k)){
	    hitnum = hitnum + 1;
	  }
	}	  
      }
      if(hitnum >0){
	fout << fRun << " " <<  time << " " << hitnum << endl;
      }
    }
    //pulserhit = pulserhit + hitnum;
    //pulserevent = pulserevent+1;
  }
  //fout << fRun << " " << pulserevent << " " << pulserhit << endl;
}

void MJDGat::SearchDelayedEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the second gamma energy
  //fEnr2 : the first Q value
  //////////////////////////////////////
  map<int,int> detid;
  for(size_t i=0; i<fChannel.size(); i++) detid[fChannel.at(i)]= i;
  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout.precision(15);
  vector<Int_t> List1;
  vector<Int_t> Chan1;
  vector<Double_t> Enr1;
  for(size_t i=0;i<fEntries;i++){
    fMjdTree->GetEntry(i);
    //if((fRun < 45000000 || fRun>50000000) && fMTEventDC1Bits == 0){
    for(size_t j=0;j<fMTChannel->size();j++){	
      Int_t chan1 = fMTChannel->at(j);
      Double_t enr1 = fMTTrapENFCal->at(j);
      Int_t index1 = detid[chan1];
      //cout << i << " " << chan1 << " " << enr1 << endl;
      if(abs(enr1-fEnr1)<20 && chan1%2==0 && fGoodBad.at(index1) == 1){
	List1.push_back(i);
	Chan1.push_back(chan1);
	Enr1.push_back(enr1);
	//cout << i << " " << chan1 << " " << enr1 << endl;
      }
    }
    //}
  }

  vector<Int_t> List2;
  vector<Int_t> List3;
  vector<Int_t> Chan2;
  vector<Int_t> Chan3;
  vector<Double_t> Enr2;
  vector<Double_t> Enr3;
  vector<Double_t> Time2;
  vector<Double_t> Time3;
  if((Int_t)List1.size()>0){
    for(size_t i=0;i<List1.size();i++){
      fMjdTree->GetEntry(List1.at(i));
      Double_t time1 = fMTTimeInfo->globalTime/1e8;
      Double_t time2 = time1;
      Int_t ii = List1.at(i)-1;
      while(abs(time1-time2)<fTime*20 && ii>0){
	fMjdTree->GetEntry(ii);
	time2 = fMTTimeInfo->globalTime/1e8;
	//if((fRun < 45000000 || fRun>50000000) && fMTEventDC1Bits == 0){	  
	  for(size_t j=0;j<fMTChannel->size();j++){
	    Int_t chan2 = fMTChannel->at(j);
	    Double_t enr2 = fMTTrapENFCal->at(j);	    
	    Int_t index = detid[chan2];
	    if(enr2<fEnr2 && enr2>5 && fGoodBad.at(index)==1 && chan2%2==0){
	      List2.push_back(List1.at(i));
	      Chan2.push_back(Chan1.at(i));
	      Enr2.push_back(Enr1.at(i));
	      Time2.push_back(time1);
	      List3.push_back(ii);
	      Chan3.push_back(chan2);
	      Enr3.push_back(enr2);			     
	      Time3.push_back(time2);
	    }
	  }
	  //}
	ii--;
      }
    }
  }

  for(size_t i=0;i<List2.size();i++){
    fout << fRun << " " << List2.at(i) << " " << Chan2.at(i) << " " << Enr2.at(i) << " " <<Time2.at(i) << " "
	 << List3.at(i) << " " << Chan3.at(i) << " " << Enr3.at(i) << " " << Time3.at(i)<<endl;
  }
}


void MJDGat::SearchEnergyEvent(Double_t fEnr, Double_t fEnrWindow, string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the gamma energy
  //////////////////////////////////////
  map<int,int> detid;
  for(size_t i=0; i<fChannel.size(); i++) {    
    detid[fChannel.at(i)]= i;
  }

  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout.precision(15);
  //vector<Int_t> List1;
  //vector<Int_t> Chan1;
  //vector<Double_t> Enr1;
  for(size_t i=0;i<fEntries;i++){
    fMjdTree->GetEntry(i);
    
    if((fRun < 45000000 || fRun>5000000)){
      for(size_t j=0;j<fMTChannel->size();j++){	
	Int_t chan = fMTChannel->at(j);
	Double_t enr = fMTTrapENFCal->at(j);
	Int_t index = detid[chan];
	/*
	Double_t nX = 1;
	if((fRun < 45000000 || fRun>50000000)){
	  nX = fMTfastTrapNLCWFsnRisingX->at(j);
	}else{
	  nX = 1;
	  }*/
	if(abs(enr-fEnr)< fEnrWindow && chan%2==0 && fGoodBad.at(index) == 1){
	  fout << fRun << " " << i << " " << chan << " " << enr << " " << fMTEventDC1Bits << endl;
	}
      }
    }else{
      for(size_t j=0;j<fMTChannel->size();j++){
        Int_t chan = fMTChannel->at(j);
        Double_t enr = fMTTrapENFCal->at(j);
        Int_t index = detid[chan];
	/*
        Double_t nX = 1;
        if((fRun < 45000000 || fRun>50000000)){
          nX = fMTfastTrapNLCWFsnRisingX->at(j);
        }else{
          nX = 1;
	  }*/

        if(abs(enr-fEnr)< fEnrWindow && chan%2==0 && fGoodBad.at(index) == 1){
          fout << fRun << " " << i << " " << chan << " " << enr << endl;
        }
      }
    }
  }

}


TH1D* MJDGat::GetWaveform(Int_t fR,Int_t fEntry, Int_t fChan,Double_t fEnr){
  GATDataSet ds(fR);
  Int_t Entry = fEntry+1;
  TCut cut1 = Form("channel == %d", fChan);
  GATWaveformBrowser wb;
  wb.LoadWaveforms(ds.GetBuiltChain(), cut1, "fWaveforms", Entry);
  size_t i=wb.GetNWaveforms();
  if(i>0){
    MGTWaveform *w1 = wb.GetWaveform(i-1).get();
    TH1D* h= (TH1D*)w1->GimmeHist();
    h->SetTitle(Form("Run=%d,Entry=%d, Channel=%d, Energy=%f(keV);Time(ns);ADC",fR,fEntry,fChan,fEnr));
    h->SetName(Form("waveform_%d_%d_%d",fR,fEntry,fChan));
    return h;
  }else{
    TH1D* h =NULL;
    return h;
    cout << " No this waveform" << endl;
  }
}



TH1D* MJDGat::GetHistoSmooth(TH1D* hist, Int_t DeltaBin){
  TH1D *h = (TH1D*)hist->Clone();
  h->Reset(0);
  Int_t entries = hist->GetEntries();
  Double_t biny=0;
  Double_t biny1=0;
  Double_t biny2=0;
  Double_t ave1=0;
  Double_t ave2=0;
  for(Int_t i1 = DeltaBin;i1<entries-DeltaBin;i1++){
    Double_t fDummySum=0;
    Double_t fDummyAve=0;
    Int_t count = 0;
    for(Int_t i2 = i1-DeltaBin;i2<i1+DeltaBin;i2++){
      biny = hist->GetBinContent(i2);
      biny1 = hist->GetBinContent(i2-1);
      biny2 = hist->GetBinContent(i2+1);
      
      if( ((abs(biny-biny1)>10 || abs(biny-biny2)>10)&& (i2<800|| i2>1200)) || (biny> 1e7&&i2>1900)){
	fDummySum = fDummySum;	
      }else{
	fDummySum = fDummySum+biny;
	count++;
      }
    }
    fDummyAve = fDummySum/(count);
    h->SetBinContent(i1,fDummyAve);
    if(i1==DeltaBin){
      ave1 = fDummyAve;
    }
    if(i1==entries-DeltaBin-1){
      ave2 = fDummyAve;
    }
  }
  for(Int_t i = 0;i<DeltaBin;i++){
    h->SetBinContent(i,ave1);
  }
  for(Int_t i = entries-DeltaBin;i<entries;i++){
    h->SetBinContent(i,ave2);
  }
  return h;
}


TH1D* MJDGat::GetHistoDerivative(TH1D* hist, Int_t DeltaBin){
  TH1D *h = (TH1D*)hist->Clone();
  h->Reset(0);
  Int_t entries = hist->GetEntries();
  Double_t binx= hist->GetBinWidth(1);
  Double_t biny1=0;
  Double_t biny2=0;
  Double_t slope=0;
  Double_t slope1=0;
  Double_t slope2=0;
  for(Int_t i = DeltaBin;i<entries-DeltaBin;i++){
    biny1 = hist->GetBinContent(i-DeltaBin);
    biny2 = hist->GetBinContent(i+DeltaBin);
    slope = (biny2-biny1)/(binx*2*DeltaBin);
    h->SetBinContent(i,slope);
    if(i==DeltaBin){
      slope1 = slope;
    }
    if(i==entries-DeltaBin-1){
      slope2 = slope;
    }
  }
  for(Int_t i = 0;i<DeltaBin;i++){
    h->SetBinContent(i,slope1);
  }
  for(Int_t i = entries-DeltaBin;i<entries;i++){
    h->SetBinContent(i,slope2);
  }
  return h;
}


TH1* MJDGat::GetHistoFFT(TH1D* hist){
  TH1 *h = NULL;
  TVirtualFFT::SetTransform(0);
  h = hist->FFT(h,"MAG");
  return h;
}


vector<Int_t> MJDGat::Sort(vector<Double_t> X){
  Int_t n = X.size();
  const Int_t n1 = n;
  Double_t x1[n1];
  Int_t index[n1];
  for(Int_t i =0;i<n;i++){
    x1[i] = X.at(i);
  }
  TMath::Sort(n,x1,index,0);
  vector<Int_t> Index;
  for(Int_t i=0;i<n;i++){
    Index.push_back(index[i]);
  }
  return Index;
}

vector<Int_t> MJDGat::Clean(vector<Double_t> X){

  vector<Int_t> Index;
  Index.push_back(0);

  for(size_t i = 1 ; i< X.size();i++){
    Int_t count = 0;
    for(size_t j=0;j<i;j++){
      if(X.at(i) == X.at(j)){
	count++;
      }
    }      
    if(count==0){
      Index.push_back(i);
    }
  }
  return Index;
}

Double_t MJDGat::GetYValue(TH1D* hist, Double_t X){
  Double_t binwidth= hist->GetBinWidth(1);
  Double_t binx=0;
  Double_t biny=0;
  Double_t Y=0;
  for(size_t i = 0;i<hist->GetEntries();i++){
    binx = hist->GetBinCenter(i);
    biny = hist->GetBinContent(i);
    if(abs(X-binx)<binwidth){
      Y = biny;
    }
  }
  return Y;
}
Int_t MJDGat::FindPeaks(TH1D* hist, Double_t Low, Double_t Up, Double_t Resolution, Double_t Sigma, Double_t Threshold, vector<Double_t>* fPositionX, vector<Double_t>* fPositionY){
  Double_t baseline = 0;
  Double_t blrms = 0;
  Int_t count =0;
  for(size_t i = 0; i< hist->GetEntries(); i++){
    if(hist->GetBinCenter(i)>1000 && hist->GetBinCenter(i)<4000){
      baseline = baseline+hist->GetBinContent(i);
      count++;
    }
  }
  baseline = baseline/(Double_t)count;
  blrms = 0;
  count = 0;
  for(size_t i = 0; i< hist->GetEntries(); i++){
    if(hist->GetBinCenter(i)>1000 && hist->GetBinCenter(i)<4000){
      blrms = blrms + pow((hist->GetBinContent(i)-baseline),2);
      count++;
    }
  }
  blrms = sqrt(blrms/(Double_t)count);
  Double_t Ymin = baseline+2*blrms;
  Double_t fResolution = Resolution;
  Double_t fSigma = Sigma;
  Double_t fThreshold = Threshold;
  TSpectrum *s = new TSpectrum(20,fResolution);
  Int_t npeaks = s->Search(hist,fSigma,"new",fThreshold);

  Double_t *xpeaks = s->GetPositionX();
  Double_t *ypeaks = s->GetPositionY();

  for(Int_t i=0;i<npeaks;i++){
    if(xpeaks[i]>Low && xpeaks[i]<Up && ypeaks[i]>Ymin){
      fPositionX->push_back(xpeaks[i]);
      fPositionY->push_back(ypeaks[i]);
    }
  }
  npeaks = fPositionX->size();
  delete s;
  return npeaks;
}

Double_t MJDGat::GetMax(TH1* hist, Double_t Low, Double_t Up){
  hist->GetXaxis()->SetRangeUser(Low,Up);
  Double_t xmax = hist->GetMaximum();
  return xmax;
}

void MJDGat::SaveSubTree(string FileName){
  cout << "Save histogram " << FileName.c_str() << endl;

  TFile *newfile = new TFile(Form("%s.root",FileName.c_str()),"recreate");
  TTree *newtree = fMjdTree->CloneTree(0);

  for(size_t i=0;i<fEntries;i++){
    //cout << i << " " << fEntries << endl;
    fMjdTree->GetEntry(i);
    newtree->Fill();    
  }
  //newtree->Print();
  newtree->AutoSave();
  delete newfile;
}

void MJDGat::IsPileUpTag(Int_t fEvent, vector<Int_t>* IsPileUp, vector<Double_t>* Ratio, vector<Double_t>* DeltaT, vector<Double_t>* AE){

  Double_t Xmin = 4000;
  Double_t Xmax = 19000;
  if((fRun >=14503 && fRun<=15892) || fRun>=25675){
    Xmax = 38000;
  }
  Double_t fResolution = 1;
  Double_t fThreshold = 0.01;
  Double_t fSigma = 5;

  vector<Double_t> xp;
  vector<Double_t> yp;
  vector<Double_t> xp1;
  vector<Double_t> yp1;
  vector<Double_t> xp2;
  vector<Double_t> yp2;

  fMjdTree->GetEntry(fEvent);
  if(fMTEventDC1Bits==0){
    for(size_t j=0;j<fMTChannel->size();j++){
      xp.clear();
      yp.clear();
      xp1.clear();
      yp1.clear();
      xp2.clear();
      yp2.clear();
      
      Int_t fChan = fMTChannel->at(j);
      Double_t fEnr = fMTTrapENFCal->at(j);
      if(fEnr>20 && fEnr<12000){
	TH1D* h = MJDGat::GetWaveform(fRun,fEvent,fChan,fEnr);
	TH1D* h1 = MJDGat::GetHistoSmooth(h,10);
	TH1D* h2 = MJDGat::GetHistoDerivative(h1,10);
	TH1D* h3 = MJDGat::GetHistoSmooth(h2,10);
	TH1D* h4 = MJDGat::GetHistoDerivative(h3,10);
	TH1D* h5 = MJDGat::GetHistoSmooth(h4,10);
	Int_t maxBin0 = h1->GetMaximumBin();
	Double_t maxX0 = h1->GetBinCenter(maxBin0);    
	Double_t A = MJDGat::GetMax(h3, Xmin,Xmax);
	Double_t AoverE = A/fEnr;
	AE->push_back(AoverE);
	
	Int_t nPeak1 = MJDGat::FindPeaks(h3, Xmin, Xmax,fResolution,fSigma, fThreshold, &xp, &yp);
	if(nPeak1>1 && AoverE>0.002){
	  for(Int_t ip = 0;ip<(Int_t)xp.size();ip++){
	    Double_t y = MJDGat::GetYValue(h5, xp.at(ip));
	    if(abs(y)< 3e-4 && xp.at(ip)< maxX0){
	      xp1.push_back(xp.at(ip));
	      yp1.push_back(yp.at(ip));
	    }
	  }
	  if(xp1.size()>1){
	    IsPileUp->push_back(1);
	    vector<Int_t> Index2 = MJDGat::Sort(yp1);
	    for(size_t ip1 = Index2.size()-2;ip1<Index2.size();ip1++){
	      Int_t ii = Index2.at(ip1);
	      xp2.push_back(xp1.at(ii));
	      yp2.push_back(yp1.at(ii));
	    }
	    Double_t ratio = yp2.at(1)/yp2.at(0);
	    Double_t deltaT = xp2.at(0)-xp2.at(1);
	    Ratio->push_back(ratio);
	    DeltaT->push_back(deltaT);	
	  }else{
	    IsPileUp->push_back(0);
	    Ratio->push_back(0);
	    DeltaT->push_back(0);
	  }
	}else{
	  IsPileUp->push_back(0);
	  Ratio->push_back(0);
	  DeltaT->push_back(0);
	}

	delete h;
	delete h1;
	delete h2;
	delete h3;  
	delete h4;
	delete h5;
      }else{
	IsPileUp->push_back(0);
	Ratio->push_back(0);
	DeltaT->push_back(0);
      }
    }
  }else{
    for(size_t j=0;j<fMTChannel->size();j++){
      IsPileUp->push_back(0);
      Ratio->push_back(0);
      DeltaT->push_back(0);
    }
  }
}

void MJDGat::PileUpTree(const char* pathName){
  TFile *newfile = new TFile(Form("%spileup_%d.root", pathName, fRun), "recreate");
  TTree *newtree = new TTree("pileupTree", "pile-up tag");
  vector<Int_t>* IsPileUp = NULL;
  vector<Double_t>* PileUpRatio = NULL;
  vector<Double_t>* PileUpDeltaT =NULL;
  vector<Double_t>* PileUpAE = NULL;
  newtree->Branch("IsPileUp",&IsPileUp);
  newtree->Branch("PileUpRatio", &PileUpRatio);
  newtree->Branch("PileUpDeltaT", &PileUpDeltaT);
  newtree->Branch("PileUpAE", &PileUpAE);
  for(size_t i = 0;i<fEntries;i++){
    IsPileUp->clear();
    PileUpRatio->clear();
    PileUpDeltaT->clear();
    PileUpAE->clear();
    MJDGat::IsPileUpTag(i,IsPileUp,PileUpRatio,PileUpDeltaT,PileUpAE);
    //cout << i << " " << IsPileUp->at(0) << " "<<  fMTChannel->at(0) << " "<< fMTTrapENFCal->at(0) << endl;
    newtree->Fill();
  }
  newtree->Write();
  cout << "Pile-up file is generated...." <<endl;
  delete newfile;
}


Int_t MJDGat::PulserCount(vector<Int_t>* PulserChannel){
  Int_t pulsernum = 0;
  const Int_t nChannel = fChannel.size();

  Int_t pulsercount[nChannel];
  for(Int_t k=0;k<nChannel;k++){
    pulsercount[k] = 0;
  }
  for(size_t i=0;i<fEntries;i++){
    fMjdTree->GetEntry(i);
    if(fMTEventDC1Bits>0){
      pulsernum++;
      for(size_t j=0;j<fMTChannel->size();j++){
	Int_t chan = fMTChannel->at(j);
	for(Int_t k=0;k<nChannel;k++){
	  if(chan == fChannel.at(k)){
	    pulsercount[k]++;
	  }
	}
      }
    }
  }

  for(Int_t k=0;k<nChannel;k++){
    PulserChannel->push_back(pulsercount[k]);
  }

  return pulsernum;
}
