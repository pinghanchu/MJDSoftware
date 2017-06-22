#include "MJDSkim.hh"
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

//ClassImp(MJDSkim)

using namespace std;
using namespace katrin;
using namespace MJDB;
//using namespace GATPeakShapeFunction;

////////////////////////////////////////////////////////

//Set ChannelMap, mjdTree (from GATDataSet), CalibrationPeak, Initial Parameters
MJDSkim::MJDSkim(Int_t DataSet, Int_t SubSet) 
{
  fDataSet = DataSet;
  fSubSet = SubSet;
  fSkimTree = new TChain("skimTree");
  fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%d/GAT-v01-06/skimDS%d_%d.root",fDataSet,fDataSet,fSubSet));

  /*
  if(fDataSet == 0){
    fSkimTree = new TChain("skimTree");
    for(Int_t i =0;i<77;i++){
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS0/GAT-v01-06/skimDS0_%d.root",i));
    }
  }else if(fDataSet == 1){    
    fSkimTree = new TChain("skimTree");
    for(Int_t i =0;i<52;i++){
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS1/GAT-v01-06/skimDS1_%d.root",i));
    }
  }else if(fDataSet == 2){
    fSkimTree = new TChain("skimTree");
    for(Int_t i =0;i<8;i++){
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS2/GAT-v01-06/skimDS2_%d.root",i));
    }
  }else if(fDataSet == 3){
    fSkimTree = new TChain("skimTree");
    for(Int_t i =0;i<25;i++){
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS3/GAT-v01-06/skimDS3_%d.root",i));
    }
  }else if(fDataSet == 4){
    fSkimTree = new TChain("skimTree");
    for(Int_t i =0;i<23;i++){
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS4/GAT-v01-06/skimDS4_%d.root",i));
    }
  }else if(fDataSet == 5){
    fSkimTree = new TChain("skimTree");
    for(Int_t i =0;i<113;i++){
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS5/GAT-v01-06/skimDS5_%d.root",i));
    }
  }else{
    cout << "No this data set!" << endl;
  }
  */
  fSkimTree->SetBranchStatus("*",1);
  
  fChannel = NULL;
  fP = NULL;
  fD = NULL;
  fC = NULL;
  fGain = NULL;
  fIsEnr = NULL;
  fIsNat = NULL;
  fIsGood = NULL;
  fTrapENFCal = NULL;
  ftOffset = NULL;
  fnX = NULL;
  fdtmu_s = NULL;

  fSkimTree->SetBranchAddress("run", &fRun);
  fSkimTree->SetBranchAddress("iEvent", &fEvent);
  fSkimTree->SetBranchAddress("channel",&fChannel);
  fSkimTree->SetBranchAddress("P",&fP);
  fSkimTree->SetBranchAddress("D",&fD);
  fSkimTree->SetBranchAddress("C",&fC);
  fSkimTree->SetBranchAddress("gain",&fGain);
  fSkimTree->SetBranchAddress("isEnr",&fIsEnr);
  fSkimTree->SetBranchAddress("isNat",&fIsNat);
  fSkimTree->SetBranchAddress("isGood",&fIsGood);
  fSkimTree->SetBranchAddress("trapENFCal", &fTrapENFCal);
  fSkimTree->SetBranchAddress("localTime_s",&flocalTime_s);
  fSkimTree->SetBranchAddress("tOffset",&ftOffset);
  fSkimTree->SetBranchAddress("mHClean",&fmH);
  fSkimTree->SetBranchAddress("nX",&fnX);
  fSkimTree->SetBranchAddress("muVeto",&fmuVeto);
  fSkimTree->SetBranchAddress("dtmu_s",&fdtmu_s);
  fSkimTree->SetBranchAddress("muTUnc",&fmuTUnc);
  fSkimTree->GetEntry(0);
  fEntries = fSkimTree->GetEntries();

}

void MJDSkim::SearchDelayedEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the second gamma energy
  //fEnr2 : the first Q value
  //////////////////////////////////////
  //cout << fTime << endl;
  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout.precision(15);
  vector<Int_t> Run1;
  vector<Int_t> List1;
  vector<Int_t> Chan1;
  vector<Double_t> Enr1;
  for(size_t i=0;i<fEntries;i++){
    fSkimTree->GetEntry(i);
    Int_t run1 = fRun;
    for(size_t j=0;j<fChannel->size();j++){	      
      Int_t chan1 = fChannel->at(j);
      Double_t enr1 = fTrapENFCal->at(j);
      if(abs(enr1-fEnr1)<20 && chan1%2==0 && fIsGood->at(j) == 1){
	Run1.push_back(run1);
	List1.push_back(i);
	Chan1.push_back(chan1);
	Enr1.push_back(enr1);
      }
    }
  }

  vector<Int_t> Run2;
  vector<Int_t> Run3;
  vector<Int_t> List2;
  vector<Int_t> List3;
  vector<Int_t> Chan2;
  vector<Int_t> Chan3;
  vector<Double_t> Enr2;
  vector<Double_t> Enr3;
  vector<Double_t> Time2;
  vector<Double_t> Time3;
  vector<Int_t> Event2;
  vector<Int_t> Event3;
  vector<Double_t> Mu_s2;
  vector<Double_t> Mu_s3;
  cout.precision(15);
  if((Int_t)List1.size()>0){
    for(size_t i=0;i<List1.size();i++){
      fSkimTree->GetEntry(List1.at(i));
      Int_t run1 = Run1.at(i);
      Double_t time1 = flocalTime_s;
      Double_t time2 = time1;
      Double_t dtmu_s1 = fdtmu_s->at(0);

      Int_t event1 = fEvent;
      Int_t ii = List1.at(i)-1;
      fSkimTree->GetEntry(ii);
      Int_t run2 = fRun;
      time2 = flocalTime_s;
      Int_t event2 = fEvent;
      
      while( (time1-time2)<fTime*20 && ii>0){
	//	cout << ii << " "<< fTime*20 << " "<< time1-time2 << " " << time1 << " " << time2 << endl;
	//fSkimTree->GetEntry(ii);
	//Int_t run2 = fRun;
	//time2 = fTimestamp->at(0);
	//Int_t event2 = fEvent;
	//cout << run1 << " " << event1 << " " << event2 << " " << ii << " "<< time1-time2 << " " << time1 << " " << time2 << endl;

	for(size_t j=0;j<fChannel->size();j++){
	  Int_t chan2 = fChannel->at(j);
	  Double_t enr2 = fTrapENFCal->at(j);	    
	  Double_t dtmu_s2 = fdtmu_s->at(j);

	  if(enr2<fEnr2 && enr2>5 && fIsGood->at(j)==1 && chan2%2==0){
	    Run2.push_back(run1);
	    List2.push_back(List1.at(i));
	    Chan2.push_back(Chan1.at(i));
	    Enr2.push_back(Enr1.at(i));
	    Time2.push_back(time1);
	    Event2.push_back(event1);
	    Mu_s2.push_back(dtmu_s1);
	    Run3.push_back(run2);
	    List3.push_back(ii);
	    Chan3.push_back(chan2);
	    Enr3.push_back(enr2);			     
	    Time3.push_back(time2);
            Event3.push_back(event2);
            Mu_s3.push_back(dtmu_s2);
	    //cout << time1 << " " << time2 << " "<< time1-time2 << endl;
	  }
	}     
	ii--;
	fSkimTree->GetEntry(ii);
        run2 = fRun;
        time2 = flocalTime_s;
        event2 = fEvent;

      }
    }
  }

  for(size_t i=0;i<List2.size();i++){
    fout << Run2.at(i) << " " << List2.at(i) << " " << Event2.at(i) << " " << Chan2.at(i) << " " << Enr2.at(i) << " " << Time2.at(i) << " " << Mu_s2.at(i) << " "
	 << Run3.at(i) << " " << List3.at(i) << " " << Event3.at(i) << " " << Chan3.at(i) << " " << Enr3.at(i) << " " << Time3.at(i) << " " << Mu_s3.at(i) << endl;
  }

}

void MJDSkim::SearchMuonCoinEvent(Double_t fEnr, Double_t fTime, string fOutputFile){
  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout.precision(15);
  vector<Int_t> Run1;
  vector<Int_t> Event1;
  vector<Int_t> Chan1;
  vector<Double_t> Enr1;
  vector<Double_t> DtMu1;
  for(size_t i=0;i<fEntries;i++){
    fSkimTree->GetEntry(i);
    Int_t run1 = fRun;
    Int_t event1 = fEvent;
    for(size_t j=0;j<fChannel->size();j++){	      
      Int_t chan1 = fChannel->at(j);
      Double_t enr1 = fTrapENFCal->at(j);
      Double_t dtmu1 = fdtmu_s->at(j);
      if(abs(enr1-fEnr)<5 && chan1%2==0 && fIsGood->at(j) == 1 && dtmu1<fTime*20){
	Run1.push_back(run1);
	Event1.push_back(event1);
	Chan1.push_back(chan1);
	Enr1.push_back(enr1);
	DtMu1.push_back(dtmu1);
      }
    }
  }
  for(size_t i=0;i<Run1.size();i++){
    fout << Run1.at(i) << " " << Event1.at(i)<< " " << Chan1.at(i) << " " << Enr1.at(i) << " " << DtMu1.at(i) << endl;
  }
}


void MJDSkim::SearchEnergyEvent(Double_t fEnr, Double_t fEnrWindow, string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the gamma energy
  //////////////////////////////////////
  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout.precision(15);
  for(size_t i=0;i<fEntries;i++){
    fSkimTree->GetEntry(i);
    Int_t irun = fRun;
    Int_t ievent = fEvent;
    Double_t time1 = flocalTime_s;
    for(size_t j=0;j<fChannel->size();j++){
      Int_t ichan = fChannel->at(j);
      Int_t ic = fC->at(j);
      Int_t ip = fP-> at(j);
      Int_t id = fD->at(j);
      Int_t isGood = (Int_t)fIsGood->at(j);
      Double_t ienr = fTrapENFCal->at(j);
      Int_t inX = fnX->at(j);
      string pos = Form("%d%d%d",ic,ip,id);
      Double_t dtmu1 = fdtmu_s->at(j);

      if(abs(ienr-fEnr)< fEnrWindow && ichan%2==0 && isGood == 1){
	fout << irun << " " << i << " " <<  ievent << " " << pos.c_str() << " "<< ichan << " "
	     << ienr << " "<< time1 << " " << dtmu1 << endl;	    
      }      
    }
  }
}


TH1D* MJDSkim::GetWaveform(Int_t fR,Int_t fEntry, Int_t fChan,Double_t fEnr){
  GATDataSet ds(fR);
  Int_t Entry = fEntry+1;
  TCut cut1 = Form("channel == %d", fChan);
  GATWaveformBrowser wb;
  wb.LoadWaveforms(ds.GetBuiltChain(), cut1, "fWaveforms", Entry);
  size_t i=wb.GetNWaveforms();
  if(i>0){
    //cout << i << endl;
    MGTWaveform *w1 = wb.GetWaveform(i-1).get();
    TH1D* h= (TH1D*)w1->GimmeHist();
    h->SetTitle(Form("Run=%d,Entry=%d, Channel=%d, Energy=%f(keV);Time(ns);ADC",fR,fEntry,fChan,fEnr));
    h->SetName(Form("waveform_%d_%d_%d",fR,fEntry,fChan));
    return h;
  }else{
    TH1D* h = NULL;
    return h;
    cout << "No this waveform!" << endl;
  }
}



TH1D* MJDSkim::GetHistoSmooth(TH1D* hist, Int_t DeltaBin){
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
      
      if(abs(biny-biny1)>10 || abs(biny-biny2)>10){
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


TH1D* MJDSkim::GetHistoDerivative(TH1D* hist, Int_t DeltaBin){
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


TH1* MJDSkim::GetHistoFFT(TH1D* hist){
  TH1 *h = NULL;
  TVirtualFFT::SetTransform(0);
  h = hist->FFT(h,"MAG");
  return h;
}


vector<Int_t> MJDSkim::Sort(vector<Double_t> X){
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

vector<Int_t> MJDSkim::Clean(vector<Double_t> X){

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

Double_t MJDSkim::GetYValue(TH1D* hist, Double_t X){
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
Int_t MJDSkim::FindPeaks(TH1D* hist, Double_t Low, Double_t Up, Double_t Resolution, Double_t Sigma, Double_t Threshold, vector<Double_t>* fPositionX, vector<Double_t>* fPositionY){
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

Double_t MJDSkim::GetMax(TH1* hist, Double_t Low, Double_t Up){
  hist->GetXaxis()->SetRangeUser(Low,Up);
  Double_t xmax = hist->GetMaximum();
  return xmax;
}


void MJDSkim::IsPileUpTag(Int_t fList, vector<Int_t>* IsPileUp, vector<Double_t>* Ratio, vector<Double_t>* DeltaT, vector<Double_t>* AE){

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

  fSkimTree->GetEntry(fList);
  
  for(size_t j=0;j<fChannel->size();j++){
    xp.clear();
    yp.clear();
    xp1.clear();
    yp1.clear();
    xp2.clear();
    yp2.clear();
    
    Int_t fChan = fChannel->at(j);
    Double_t fEnr = fTrapENFCal->at(j);
    if(fEnr>40 && fEnr<12000){
      TH1D* h = MJDSkim::GetWaveform(fRun,fEvent,fChan,fEnr);
      TH1D* h1 = MJDSkim::GetHistoSmooth(h,10);
      TH1D* h2 = MJDSkim::GetHistoDerivative(h1,10);
      TH1D* h3 = MJDSkim::GetHistoSmooth(h2,10);
      TH1D* h4 = MJDSkim::GetHistoDerivative(h3,10);
      TH1D* h5 = MJDSkim::GetHistoSmooth(h4,10);
      Int_t maxBin0 = h1->GetMaximumBin();
      Double_t maxX0 = h1->GetBinCenter(maxBin0);    
      Double_t A = MJDSkim::GetMax(h3, Xmin,Xmax);
      Double_t AoverE = A/fEnr;
      AE->push_back(AoverE);
      
      Int_t nPeak1 = MJDSkim::FindPeaks(h3, Xmin, Xmax,fResolution,fSigma, fThreshold, &xp, &yp);
      if(nPeak1>1 && AoverE>0.002){
	for(Int_t ip = 0;ip<(Int_t)xp.size();ip++){
	  Double_t y = MJDSkim::GetYValue(h5, xp.at(ip));
	  if(abs(y)< 3e-4 && xp.at(ip)< maxX0){
	    xp1.push_back(xp.at(ip));
	    yp1.push_back(yp.at(ip));
	  }
	}
	if(xp1.size()>1){
	  IsPileUp->push_back(1);
	  vector<Int_t> Index2 = MJDSkim::Sort(yp1);
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
}


void MJDSkim::PileUpTree(const char* pathName){
  TFile *newfile = new TFile(Form("%spileup_%d_%d.root", pathName, fDataSet,fSubSet), "recreate");
  TTree *newtree = new TTree("pileupTree", "pile-up tag");
  vector<Int_t>* IsPileUp = NULL;
  vector<Double_t>* PileUpRatio = NULL;
  vector<Double_t>* PileUpDeltaT =NULL;
  vector<Double_t>* PileUpAE = NULL;
  vector<Double_t>* Energy = NULL;
  newtree->Branch("IsPileUp",&IsPileUp);
  newtree->Branch("PileUpRatio", &PileUpRatio);
  newtree->Branch("PileUpDeltaT", &PileUpDeltaT);
  newtree->Branch("PileUpAE", &PileUpAE);
  newtree->Branch("Energy",&Energy);
  for(size_t i = 0;i<fEntries;i++){
    IsPileUp->clear();
    PileUpRatio->clear();
    PileUpDeltaT->clear();
    PileUpAE->clear();
    MJDSkim::IsPileUpTag(i,IsPileUp,PileUpRatio,PileUpDeltaT,PileUpAE);
    cout << fRun << " " << fEvent<<  " " << i << " " << IsPileUp->at(0) << " "<<  fChannel->at(0) << " "<< fTrapENFCal->at(0) << endl;
    newtree->Fill();
  }
  newtree->Write();
  cout << "Pile-up file is generated...." <<endl;
  delete newfile;
}
