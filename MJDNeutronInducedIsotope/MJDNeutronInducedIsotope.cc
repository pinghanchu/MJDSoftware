#include "MJDNeutronInducedIsotope.hh"
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

//ClassImp(MJDNeutronInducedIsotope)

using namespace std;
using namespace katrin;
using namespace MJDB;
//using namespace GATPeakShapeFunction;

////////////////////////////////////////////////////////

//Set ChannelMap, mjdTree (from GATDataSet), CalibrationPeak, Initial Parameters
MJDNeutronInducedIsotope::MJDNeutronInducedIsotope(Int_t Run) : fDS(Run,Run)
{
  fMjdTree = fDS.GetMJDTree();
  fMap = fDS.GetMap();
  fEntries = fDS.GetEntries();
  fRun = Run;
  fGATRev = fDS.GetGATRev();
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


  fMjdTree->SetBranchStatus("*",1);

  fMTChannel = NULL;
  fMTTrapENFCal = NULL;
  fMTTimestamp = NULL;
  fMTfastTrapNLCWFsnRisingX = NULL;

  fMjdTree->SetBranchAddress("channel", &fMTChannel);
  fMjdTree->SetBranchAddress("trapENFCal", &fMTTrapENFCal);
  fMjdTree->SetBranchAddress("timestamp",&fMTTimestamp);
  fMjdTree->SetBranchAddress("mH",&fMTmH);
  fMjdTree->SetBranchAddress("EventDC1Bits", &fMTEventDC1Bits);
  fMjdTree->SetBranchAddress("fastTrapNLCWFsnRisingX", &fMTfastTrapNLCWFsnRisingX);
}

void MJDNeutronInducedIsotope::SearchDelayedEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile){
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
    if(fMTEventDC1Bits == 0){
      for(size_t j=0;j<fMTChannel->size();j++){	
	Int_t chan1 = fMTChannel->at(j);
	Double_t enr1 = fMTTrapENFCal->at(j);
	if(abs(enr1-fEnr1)<5 && chan1%2==0){
	  List1.push_back(i);
	  Chan1.push_back(chan1);
	  Enr1.push_back(enr1);
	}
      }
    }
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
      Double_t time1 = fMTTimestamp->at(0)/1e8;
      Double_t time2 = time1;
      Int_t ii = List1.at(i)-1;
      while(abs(time1-time2)<fTime*200 && ii>0){
	fMjdTree->GetEntry(ii);
	time2 = fMTTimestamp->at(0)/1e8;
	if(fMTEventDC1Bits == 0){	  
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
	}
	ii--;
      }
    }
  }

  for(size_t i=0;i<List2.size();i++){
    fout << fRun << " " << List2.at(i) << " " << Chan2.at(i) << " " << Enr2.at(i) << " " <<Time2.at(i) << " "
	 << List3.at(i) << " " << Chan3.at(i) << " " << Enr3.at(i) << " " << Time3.at(i)<<endl;
  }
}


void MJDNeutronInducedIsotope::SearchEnergyEvent(Double_t fEnr1, string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the gamma energy
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
    if(fMTEventDC1Bits == 0){
      for(size_t j=0;j<fMTChannel->size();j++){	
	Int_t chan1 = fMTChannel->at(j);
	Double_t enr1 = fMTTrapENFCal->at(j);
	if(abs(enr1-fEnr1)<5 && chan1%2==0){
	  List1.push_back(i);
	  Chan1.push_back(chan1);
	  Enr1.push_back(enr1);
	}
      }
    }
  }

  for(size_t i=0;i<List1.size();i++){
    fout << fRun << " " << List1.at(i) << " " << Chan1.at(i) << " " << Enr1.at(i) << endl;
  }
}


TH1D* MJDNeutronInducedIsotope::GetWaveform(Int_t fR,Int_t fEntry, Int_t fChan,Double_t fEnr){
  GATDataSet ds(fR);
  Int_t Entry = fEntry+1;
  TCut cut1 = Form("channel == %d", fChan);
  GATWaveformBrowser wb;
  wb.LoadWaveforms(ds.GetBuiltChain(), cut1, "fWaveforms", Entry);
  size_t i=wb.GetNWaveforms();
  MGTWaveform *w1 = wb.GetWaveform(i-1).get();
  TH1D* h= (TH1D*)w1->GimmeHist();
  h->SetTitle(Form("Run=%d,Entry=%d, Channel=%d, Energy=%f(keV);Time(ns);ADC",fR,fEntry,fChan,fEnr));
  h->SetName(Form("waveform_%d_%d_%d",fR,fEntry,fChan));
  return h;
}



TH1D* MJDNeutronInducedIsotope::GetHistoSmooth(TH1D* hist, Int_t DeltaBin){
  TH1D *h = (TH1D*)hist->Clone();
  h->Reset(0);
  Int_t entries = hist->GetEntries();
  Double_t biny=0;
  Double_t ave1=0;
  Double_t ave2=0;
  for(Int_t i1 = DeltaBin;i1<entries-DeltaBin;i1++){
    Double_t fDummySum=0;
    Double_t fDummyAve=0;
    for(Int_t i2 = i1-DeltaBin;i2<i1+DeltaBin;i2++){
      biny = hist->GetBinContent(i2);
      fDummySum = fDummySum+biny;
    }
    fDummyAve = fDummySum/(2*DeltaBin);
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


TH1D* MJDNeutronInducedIsotope::GetHistoDerivative(TH1D* hist, Int_t DeltaBin){
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


TH1* MJDNeutronInducedIsotope::GetHistoFFT(TH1D* hist){
  TH1 *h = NULL;
  TVirtualFFT::SetTransform(0);
  h = hist->FFT(h,"MAG");
  return h;
}


vector<Int_t> MJDNeutronInducedIsotope::Sort(vector<Double_t> X){
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

vector<Int_t> MJDNeutronInducedIsotope::Clean(vector<Double_t> X){

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

Double_t MJDNeutronInducedIsotope::GetYValue(TH1D* hist, Double_t X){
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
Int_t MJDNeutronInducedIsotope::FindPeaks(TH1D* hist, Double_t Low, Double_t Up, Double_t Resolution, Double_t Sigma, Double_t Threshold, vector<Double_t>* fPositionX, vector<Double_t>* fPositionY){
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

  delete s;
  return npeaks;
}
