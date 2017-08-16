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
#include <iomanip> 
//ClassImp(MJDSkim)

using namespace std;
using namespace katrin;
using namespace MJDB;
//using namespace GATPeakShapeFunction;

////////////////////////////////////////////////////////

//Set ChannelMap, mjdTree (from GATDataSet), CalibrationPeak, Initial Parameters
MJDSkim::MJDSkim(Int_t DataSet, Int_t SubSet, Int_t IsCal) 
{
  fDataSet = DataSet;
  fSubSet = SubSet;
  fSkimTree = new TChain("skimTree");
  fIsCal = IsCal;
  string path = "GAT-v01-06-125-gd9332b6";

  if(IsCal == 0){
    fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%d/%s/skimDS%d_%d.root",fDataSet,path.c_str(),fDataSet,fSubSet));
  }else{
    fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%dcal/%s/skimDS%d_run%d_small.root",fDataSet,path.c_str(),fDataSet,fSubSet));
  }

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
  fDCR = NULL;
  fAvsE = NULL;
  fnX = NULL;
  fdtmu_s = NULL;
  fglobalTime = NULL;
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
  fSkimTree->SetBranchAddress("trapENFCalC", &fTrapENFCal);
  fSkimTree->SetBranchAddress("globalTime",&fglobalTime);
  fSkimTree->SetBranchAddress("localTime_s",&flocalTime_s);
  fSkimTree->SetBranchAddress("clockTime_s",&fclockTime_s);
  fSkimTree->SetBranchAddress("tOffset",&ftOffset);
  fSkimTree->SetBranchAddress("dcr99",&fDCR);
  fSkimTree->SetBranchAddress("avse",&fAvsE);
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
  fout << fixed << setprecision(5);

  vector<Int_t> Run1;
  vector<Int_t> List1;
  vector<Int_t> Chan1;
  vector<Double_t> Enr1;
  vector<Double_t> DCR1;
  vector<Double_t> Time1;
  for(size_t i=0;i<fEntries;i++){
    fSkimTree->GetEntry(i);
    Int_t run1 = fRun;
    for(size_t j=0;j<fChannel->size();j++){	      
      Int_t chan1 = fChannel->at(j);
      Double_t enr1 = fTrapENFCal->at(j);
      Double_t dcr1 = fDCR->at(j);
      Double_t time1 = fglobalTime->AsDouble()+ ftOffset->at(j)/1e9;
      //cout << fglobalTime->AsDouble() << " "<< ftOffset->at(j) << endl;
      if(abs(enr1-fEnr1)<20 && chan1%2==0 && fIsGood->at(j) == 1 ){
      //if(abs(enr1-fEnr1)<20){
	Run1.push_back(run1);
	List1.push_back(i);
	Chan1.push_back(chan1);
	Enr1.push_back(enr1);
	DCR1.push_back(dcr1);
	Time1.push_back(time1);
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
  vector<Double_t> DCR2;
  vector<Double_t> DCR3;

  cout.precision(15);
  if((Int_t)List1.size()>0){
    for(size_t i=0;i<List1.size();i++){
      fSkimTree->GetEntry(List1.at(i));
      Int_t run1 = Run1.at(i);
      //Double_t time1 = flocalTime_s;
      Double_t time1 = Time1.at(i);
      Double_t time2 = time1;
      Double_t dtmu_s1 = fdtmu_s->at(0);

      Int_t event1 = fEvent;
      Int_t ii = List1.at(i)-1;
      fSkimTree->GetEntry(ii);
      Int_t run2 = fRun;
      time2 = fglobalTime->AsDouble();
      Int_t event2 = fEvent;
      Double_t dcr2 = DCR1.at(i);
      while( (time1-time2)<fTime*20 && (time1-time2) >0 && ii>0){
	for(size_t j=0;j<fChannel->size();j++){
	  Int_t chan2 = fChannel->at(j);
	  Double_t enr2 = fTrapENFCal->at(j);	    
	  Double_t dtmu_s2 = fdtmu_s->at(j);
	  Double_t dt = ftOffset->at(j)/1e9;
	  Double_t dcr3 = fDCR->at(j);
	  cout << event1 << " " << time1 << " "<< event2 << " "<< time2 << " " << dt << " " << chan2 << endl;
	  if(enr2<fEnr2 && enr2>5 && fIsGood->at(j)==1 && chan2%2==0 && (time1 - (time2+dt))>0){
	    Run2.push_back(run1);
	    List2.push_back(List1.at(i));
	    Chan2.push_back(Chan1.at(i));
	    Enr2.push_back(Enr1.at(i));
	    Time2.push_back(time1);
	    Event2.push_back(event1);
	    Mu_s2.push_back(dtmu_s1);
	    DCR2.push_back(dcr2);
	    Run3.push_back(run2);
	    List3.push_back(ii);
	    Chan3.push_back(chan2);
	    Enr3.push_back(enr2);			     
	    Time3.push_back(time2+dt);
            Event3.push_back(event2);
            Mu_s3.push_back(dtmu_s2);
	    DCR3.push_back(dcr3);
	    //cout << time1 << " " << time2 << " "<< time1-time2 << endl;
	  }
	}     
	ii--;
	fSkimTree->GetEntry(ii);
        run2 = fRun;
        time2 = fglobalTime->AsDouble();
        event2 = fEvent;
      }
    }
  }

  for(size_t i=0;i<List2.size();i++){
    fout << Run2.at(i) << " " << List2.at(i) << " " << Event2.at(i) << " " << Chan2.at(i) << " " << Enr2.at(i) << " " << Time2.at(i) << " " << Mu_s2.at(i) << " " << DCR2.at(i) << " "
	 << Run3.at(i) << " " << List3.at(i) << " " << Event3.at(i) << " " << Chan3.at(i) << " " << Enr3.at(i) << " " << Time3.at(i) << " " << Mu_s3.at(i) << " " << DCR3.at(i) << " " 
	 << Time2.at(i) - Time3.at(i) << endl;

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
      if(abs(enr1-fEnr)<20 && chan1%2==0 && fIsGood->at(j) == 1 && dtmu1<fTime*20){
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
      //Int_t inX = fnX->at(j);
      string pos = Form("%d%d%d",ic,ip,id);
      Double_t dtmu1 = fdtmu_s->at(j);
      Double_t dcr = fDCR->at(j);
      Double_t avse = fAvsE->at(j);
      //cout << irun << " " << ievent <<" " << ienr << endl;
      //if(abs(ienr-fEnr)< fEnrWindow && ichan%2==0 && isGood == 1 && dcr < 0.006 && dcr>0.004){	
      //if(abs(ienr-fEnr)< fEnrWindow && ichan%2==0 &&  dcr>0.004 && dcr<0.006){
      if(abs(ienr-fEnr)< fEnrWindow){

	fout << irun << " " << i << " " <<  ievent << " " << pos.c_str() << " "<< ichan << " "
	     << ienr << " "<< dcr << " "<< avse << endl;	    
      }      
    }
  }
  cout << "SearchEnergy file is generated...." <<endl;

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

TH1D* MJDSkim::TrapezoidalFilter(TH1D* hist, double RampTime, double FlatTime , double DecayTime )
{

  TH1D *h = (TH1D*)hist->Clone();
  h->Reset(0);
  Int_t entries = hist->GetEntries();
  vector<Double_t> anInput;
  for(Int_t i1 = 0;i1<entries;i1++){
    cout << i1 << " "<< hist->GetBinContent(i1) << endl;
    anInput.push_back(hist->GetBinContent(i1));
  }

  // Arbitrarily choosing 1000 as decay constant
  double decayConstant = 0.0;
  if(DecayTime != 0) decayConstant = 1./(exp(1./DecayTime) - 1);
  // double decayConstant = 0;
  double rampStep = RampTime;
  double flatStep = FlatTime;
  double baseline = 0; // No baseline for now
  double norm = rampStep;
  if(decayConstant != 0)norm *= decayConstant;

  std::vector<double> fVector;
  std::vector<double> anOutput;
  if(fVector.size() != anInput.size()) 
    {
      fVector.resize(anInput.size());
      anOutput.resize(anInput.size());
    }

  fVector[0] = anInput[0] - baseline;
  anOutput[0] = (decayConstant+1.)*(anInput[0] - baseline);
  double scratch = 0.0;
  for(size_t i = 1; i < anInput.size(); i++)
    {
      // This is a little tricky with all the ternary operators, but it's faster
      // this way.  We must check the bounds.
      scratch = anInput[i]  - ((i>=rampStep) ? anInput[i-rampStep] : baseline)
	- ((i>=flatStep+rampStep) ? anInput[i-flatStep-rampStep] : baseline)
	+ ((i>=flatStep+2*rampStep) ? anInput[i-flatStep-2*rampStep] : baseline);  
      if(decayConstant != 0.0) 
	{
	  fVector[i] = fVector[i-1] + scratch; 
	  anOutput[i] = (anOutput[i-1] + fVector[i] + decayConstant*scratch);
	  // anOutput[i] /= norm;
	} 
      else anOutput[i] = anOutput[i-1] + scratch;
    }
  for(size_t i = 2*rampStep+flatStep; i < anInput.size(); i++) anOutput[i-(2*rampStep+flatStep)] = anOutput[i];
  anOutput.resize(anOutput.size()-(2*rampStep+flatStep));
  // Rescale event by normalization factor
  for(size_t i = 1; i < anOutput.size(); i++)anOutput[i] = anOutput[i]/norm;
  cout << anOutput.size() << " "<< entries << endl;
  for(size_t i2 = 0;i2<anOutput.size();i2++){
    h->SetBinContent(i2,anOutput.at(i2));
  }
  return h;
  //return anOutput;
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
  h->Draw("wf.pdf");
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
  TSpectrum *s = new TSpectrum(100,fResolution);
  Int_t npeaks = s->Search(hist,fSigma,"new",fThreshold);

  Double_t *xpeaks = s->GetPositionX();
  Double_t *ypeaks = s->GetPositionY();

  for(Int_t i=0;i<npeaks;i++){
    if(xpeaks[i]>Low && xpeaks[i]<Up && ypeaks[i]>Ymin){
      hist->GetXaxis()->SetRangeUser(xpeaks[i]-20,xpeaks[i]+20);
      double ymax = hist->GetMaximum();
      int xbin = hist->GetMaximumBin();
      double xmax = hist->GetBinCenter(xbin);
      fPositionX->push_back(xmax);
      fPositionY->push_back(ymax);
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


void MJDSkim::IsPileUpTag(Int_t fList, vector<Int_t>* IsPileUp, vector<Double_t>* Ratio, vector<Double_t>* DeltaT, vector<Double_t>* AE,vector<Double_t>* Cur){

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
  ofstream fout(Form("pileup_%d.txt",fRun),ios::app);
  for(size_t j=0;j<fChannel->size();j++){
    xp.clear();
    yp.clear();
    xp1.clear();
    yp1.clear();
    xp2.clear();
    yp2.clear();
    
    Int_t fChan = fChannel->at(j);
    Double_t fEnr = fTrapENFCal->at(j);
    Double_t DCR = fDCR->at(j);
    if(fEnr>40 && fEnr<12000 && abs(DCR)<0.002){
      //cout << fEnr << endl;

    //if(fEnr>53-5 && fEnr<67+5){

      TH1D* h = MJDSkim::GetWaveform(fRun,fEvent,fChan,fEnr);
      //cout << fRun << " " << fEvent << " " << fChan << " " << fEnr << endl;
      //TH1D* h1 = MJDSkim::GetHistoSmooth(h,10);
      Double_t maxY = h->GetMaximum();
      if(maxY<10000){
	TH1D* h2 = MJDSkim::GetHistoDerivative(h,10);
	//TH1D* h3 = MJDSkim::GetHistoSmooth(h2,10);
	TH1D* h4 = MJDSkim::GetHistoDerivative(h2,10);
	//TH1D* h5 = MJDSkim::GetHistoSmooth(h4,10);
	
	Int_t maxBin0 = h->GetMaximumBin();
	Double_t maxX0 = h->GetBinCenter(maxBin0);    
	Double_t A = MJDSkim::GetMax(h2, Xmin,Xmax);
	Double_t AoverE = A/fEnr;
	AE->push_back(AoverE);
	Cur->push_back(A);
	Int_t nPeak1 = MJDSkim::FindPeaks(h2, Xmin, Xmax,fResolution,fSigma, fThreshold, &xp, &yp);
	//cout << nPeak1 << endl;
	if(nPeak1>1){
	  for(Int_t ip = 0;ip<(Int_t)xp.size();ip++){
	    Double_t yabsmin =  10;
	    for(Int_t is = 0;is<40;is++){
	      Double_t y = abs(MJDSkim::GetYValue(h4, xp.at(ip)-20+is));
	      if(y<yabsmin){
		yabsmin = y;
	      }
	    }
	    //cout << yabsmin << " "<< xp.at(ip) << " "<< yp.at(ip) << endl;
	    if(abs(yabsmin)< 3e-4 && xp.at(ip)< maxX0){
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
	    fout << fRun << " " << fEvent << " " << fChan << " " << fEnr << " " << ratio << " " << deltaT << " " << AoverE << endl;
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

	delete h2;
	//delete h3;
	delete h4;
	//delete h5;
      }else{
	IsPileUp->push_back(0);
	Ratio->push_back(0);
	DeltaT->push_back(0);
        AE->push_back(0);
	Cur->push_back(0);
      }
      delete h;
      //delete h1;
    }else{
      IsPileUp->push_back(0);
      Ratio->push_back(0);
      DeltaT->push_back(0);
      AE->push_back(0);
      Cur->push_back(0);
    }
  }
}


void MJDSkim::PileUpTree(const char* pathName){
  TFile *newfile = new TFile(Form("%spileup_%d_%d.root", pathName, fDataSet,fSubSet), "recreate");
  TTree *newtree = new TTree("pileupTree", "pile-up tag");

  vector<Int_t>* IsPileUp = NULL;
  vector<Double_t>* PileUpRatio = NULL;
  vector<Double_t>* PileUpDeltaT =NULL;
  vector<Double_t>* AE = NULL;
  vector<Double_t>* A = NULL;

  vector<Int_t> IsPileUp1;
  vector<Double_t> PileUpRatio1;
  vector<Double_t> PileUpDeltaT1;
  vector<Double_t> AE1;
  vector<Double_t> A1;

  newtree->Branch("IsPileUp",&IsPileUp);
  newtree->Branch("PileUpRatio", &PileUpRatio);
  newtree->Branch("PileUpDeltaT", &PileUpDeltaT);
  newtree->Branch("AE", &AE);
  newtree->Branch("A", &A);

  Int_t Run;
  Int_t Event;
  vector<Int_t>* Channel = NULL;
  vector<Int_t>* C = NULL;
  vector<Int_t>* P = NULL;
  vector<Int_t>* D = NULL;
  vector<bool>* IsEnr = NULL;
  vector<bool>* IsNat = NULL;
  vector<bool>* IsGood = NULL;
  vector<Double_t>* TOffset = NULL;
  vector<Double_t>* Dtmu_s = NULL;
  vector<Double_t>* Energy = NULL;
  vector<Double_t>* DCR = NULL;
  vector<Double_t>* AvsE = NULL;
  vector<Double_t>* ToverE = NULL;
  //vector<Double_t>* TimeDiff = NULL;
  Double_t LocalTime;
  Double_t GlobalTime;
  Int_t mHClean;
  bool MuVeto;
  Double_t MuTUnc;

  newtree->Branch("Run", &Run);
  newtree->Branch("Event", &Event);
  newtree->Branch("Channel",&Channel);
  newtree->Branch("Energy",&Energy);
  newtree->Branch("P",&P);
  newtree->Branch("D",&D);
  newtree->Branch("C",&C);
  newtree->Branch("IsEnr",&IsEnr);
  newtree->Branch("IsNat",&IsNat);
  newtree->Branch("IsGood",&IsGood);
  newtree->Branch("LocalTime",&LocalTime);
  newtree->Branch("GlobalTime",&GlobalTime);
  newtree->Branch("TOffset",&TOffset);
  newtree->Branch("mHClean",&mHClean);
  newtree->Branch("MuVeto",&MuVeto);
  newtree->Branch("Dtmu_s",&Dtmu_s);
  newtree->Branch("MuTUnc",&MuTUnc);
  newtree->Branch("DCR",&DCR);
  newtree->Branch("AvsE",&AvsE);
  newtree->Branch("ToverE",&ToverE);
  //newtree->Branch("TimeDiff", &TimeDiff);
  Int_t entries = 0;
  if(fIsCal == 0){
    entries = fSkimTree->GetEntries();
  }else{
    entries = 1000;
  }
  for(Int_t i=0;i<entries;i++){
    Double_t count = 0;
    //cout << i << endl;
    fSkimTree->GetEntry(i);
    IsPileUp1.clear();
    PileUpRatio1.clear();
    PileUpDeltaT1.clear();
    AE1.clear();
    A1.clear();
    IsPileUp->clear();
    PileUpRatio->clear();
    PileUpDeltaT->clear();
    AE->clear();
    A->clear();
    
    MJDSkim::IsPileUpTag(i,&IsPileUp1,&PileUpRatio1,&PileUpDeltaT1,&AE1,&A1);
    
    Channel->clear();
    Energy->clear();
    P->clear();
    D->clear();
    C->clear();
    IsEnr->clear();
    IsNat->clear();
    IsGood->clear();
    TOffset->clear();
    Dtmu_s->clear();
    DCR->clear();
    AvsE->clear();
    
    Run = fRun;
    Event = fEvent;
    LocalTime = flocalTime_s;
    GlobalTime = fglobalTime->AsDouble();
    mHClean = fmH;
    MuVeto = fmuVeto;
    MuTUnc = fmuTUnc;
    
    //cout << IsPileUp1.size() << " "<< fChannel->size() << endl;
    for(size_t j = 0;j<fChannel->size();j++){
      if(IsPileUp1.at(j)>0){
	
	IsPileUp->push_back(IsPileUp1.at(j));
	PileUpRatio->push_back(PileUpRatio1.at(j));
	PileUpDeltaT->push_back(PileUpDeltaT1.at(j));
	AE->push_back(AE1.at(j));
	A->push_back(A1.at(j));
	
	Channel->push_back(fChannel->at(j));
	Energy->push_back(fTrapENFCal->at(j));
	P->push_back(fP->at(j));
	D->push_back(fD->at(j));
	C->push_back(fC->at(j));
	IsEnr->push_back(fIsEnr->at(j));
	IsNat->push_back(fIsNat->at(j));
	IsGood->push_back(fIsGood->at(j));
	TOffset->push_back(ftOffset->at(j)/1e9);
	Dtmu_s->push_back(fdtmu_s->at(j));	       
	DCR->push_back(fDCR->at(j));
	AvsE->push_back(fAvsE->at(j));
      }
    }
    newtree->Fill();      
  }
  newtree->Write();
  cout << "Pile-up file is generated...." <<endl;
  delete newfile;
}

void MJDSkim::TimeDiffTree(const char* pathName){
  TFile *newfile = new TFile(Form("%stimediff_%d_%d.root", pathName, fDataSet,fSubSet), "recreate");
  TTree *newtree = new TTree("timediffTree", "");
  vector<Double_t>* TimeDiff = NULL;
  newtree->Branch("TimeDiff",&TimeDiff);
  Int_t entries = 0;
  if(fIsCal == 0){
    entries = fSkimTree->GetEntries();
  }else{
    entries = 1000;
  }
  for(Int_t i=0;i<entries;i++){
    Double_t count = 0;
    fSkimTree->GetEntry(i);
    TimeDiff->clear();
    if(fmH>1){
      for(size_t j = 0;j<fChannel->size();j++){
	
	vector<Double_t> tOff;
	tOff.clear();
	for(size_t j = 0;j<fChannel->size();j++){
	  if(fChannel->at(j)%2==0){
	    tOff.push_back(ftOffset->at(j)/1e9);
	  }
	}
	vector<Int_t> tInd = MJDSkim::Sort(tOff);	
	vector<Double_t> tDiff;
	tDiff.clear();
   
	if(tOff.size()>0){
	  Double_t t0 = tOff.at(0);
	  for(size_t j=0;j<tOff.size();j++){
	    Double_t t1 = tOff.at(j);
	    TimeDiff->push_back(t1-t0);
	    t0 = t1;
	  }
	}
      }
    }
    newtree->Fill();      
  }
  newtree->Write();
  cout << "TimeDiff file is generated...." <<endl;
  delete newfile;
}

TH1D* MJDSkim::FillHisto(TChain* mTree, string InputParaName, string OutputParaName,string CutName, Int_t Bin, Double_t Low, Double_t Up){
  TH1D *h = new TH1D(Form("%s", OutputParaName.c_str()),"", Bin, Low, Up);
  TCut cut1 = Form("%s",CutName.c_str());
  //TCanvas *c1 = new TCanvas("c1");
  mTree->Draw(Form("%s>>%s",InputParaName.c_str(),OutputParaName.c_str()),cut1);
  return h;
  cout << "histogram is generated" << endl;
}
