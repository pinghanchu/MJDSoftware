#ifndef _MJDGat_hh
#define _MJDGat_hh

#include "GATDataSet.hh"
#include "GATAutoCal.hh"
#include "GATTimeInfo.hh"
#include "MJAnalysisDoc.hh"
#include "MJProvenance.hh"
#include "TObject.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TMatrixDSym.h"
#include "TLegend.h"
#include <vector>
#include <string>

using namespace std;


class MJDGat : public TObject  
{
public:   
  MJDGat(Int_t Run = 1);
  virtual ~MJDGat() {}
  virtual inline Int_t GetRun() { return fRun; }
  virtual inline Int_t GetEntries() { return fEntries; }
  virtual inline Int_t GetIsRadio() { return fIsRadio; }
  virtual inline string GetDataSet(){ return fDataSet; }
  virtual inline vector<Int_t> GetPulserChannel(){ return fPulserTagChannel;}
  virtual inline vector<Int_t> GetChannel(){ return fChannel; }
  virtual inline vector<Int_t> GetCryo(){ return fCryo;}
  virtual inline vector<Int_t> GetString(){ return fString; }
  virtual inline vector<Int_t> GetDetPosition(){ return fDetector; }
  virtual inline vector<Int_t> GetGoodBad(){ return fGoodBad; }
  virtual inline vector<Int_t> GetEnriched(){ return fEnriched;}
  virtual inline vector<string> GetDetectorName(){ return fDetectorName;}
  virtual inline TChain* GetMJDTree(){ return fMjdTree; }
  virtual inline MJTChannelMap* GetMap(){ return fMap; }
  virtual inline Double_t GetStartTime(){return fMTStartTime;}
  virtual inline Double_t GetStopTime(){return fMTStopTime;}
  virtual void SeachPulserChannel(string fOutputFile);
  virtual void SearchDelayedEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile);
  virtual void SearchEnergyEvent(Double_t fEnr1, Double_t fEnrWindow, string fOutputFile);
  virtual TH1D* GetWaveform(Int_t fR,Int_t fEntry, Int_t fChan,Double_t fEnr);
  virtual TH1D* GetHistoSmooth(TH1D* hist, Int_t DeltaBin);
  virtual TH1D* GetHistoDerivative(TH1D* hist, Int_t DeltaBin);
  virtual TH1* GetHistoFFT(TH1D* hist);
  virtual Int_t FindPeaks(TH1D* hist, Double_t Low, Double_t Up, Double_t Resolution, Double_t Sigma, Double_t Threshold, vector<Double_t>* fPositionX, vector<Double_t>* fPositionY);
  virtual Double_t GetYValue(TH1D* hist, Double_t X);
  virtual Double_t GetMax(TH1* hist, Double_t Low, Double_t Up);
  virtual vector<Int_t> Sort(vector<Double_t> X);
  virtual vector<Int_t> Clean(vector<Double_t> X); 
  virtual void SaveSubTree(string FileName);
  virtual void IsPileUpTag(Int_t fEvent, vector<Int_t>* IsPileUp,vector<Double_t>* Ratio, vector<Double_t>* DeltaT, vector<Double_t>* AE);
  virtual void PileUpTree(const char* pathName);

protected:
  GATAutoCal fDS;
  TChain* fMjdTree;
  MJTChannelMap* fMap;

  size_t fEntries;
  Int_t fRun;
  //UInt_t fGATRev;
  Int_t fIsRadio;
  Double_t fMTStartTime;
  Double_t fMTStopTime;
  string fDataSet;
  vector<Int_t> fChannel;
  vector<Int_t> fPulserTagChannel;
  vector<Int_t> fCryo;
  vector<Int_t> fString;
  vector<Int_t> fDetector;
  vector<Int_t> fGoodBad;
  vector<Int_t> fEnriched;
  vector<string> fDetectorName;

  vector<Double_t>* fMTChannel;
  vector<Double_t>* fMTTrapE;
  vector<Double_t>* fMTTrapENFCal;
  //vector<Double_t>* fMTTimestamp;
  GATTimeInfo* fMTTimeInfo;
  Int_t fMTmH;
  Int_t fMTEventDC1Bits;
  vector<Double_t>* fMTfastTrapNLCWFsnRisingX;


private:


};

#endif
