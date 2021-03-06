#ifndef _MJDSkim_hh
#define _MJDSkim_hh

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
#include "TTimeStamp.h"
#include "TLegend.h"
#include <vector>
#include <string>

using namespace std;


class MJDSkim : public TObject  
{
public:   
  MJDSkim(Int_t DataSet = 1, Int_t SubSet = 0, Int_t SubSet2=0, Int_t IsCal = 0);
  virtual ~MJDSkim() {}
  virtual inline TChain* GetSkimTree(){ return fSkimTree; }
  virtual void SetChannel();
  virtual inline vector<Int_t> GetChannel(){return fMTChannel;};
  virtual inline vector<Int_t> GetCryo(){ return fCryo;}
  virtual inline vector<Int_t> GetString(){ return fString; }
  virtual inline vector<Int_t> GetDetPosition(){ return fDetector; }
  virtual inline vector<string> GetDetectorName(){ return fDetectorName;}
  virtual void FindEvent(Int_t List, vector<Int_t>* Channel, vector<Double_t>* Energy);
  virtual void FindDelayedEvent(Int_t List, Int_t Chan, Double_t fTime, string fOutputFile);   
  virtual void SearchDelayedEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile);
  virtual void SearchDelayedHighEnergyEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile);
  virtual void SearchDelayedMultiplicityEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile);

  virtual void SearchMuonCoinEvent(Double_t fEnr, Double_t fTime, string fOutputFile);
  virtual void SearchEnergyEvent(Double_t fEnr1, Double_t fEnrWindow, string fOutputFile);
  virtual TH1D* GetWaveform(Int_t fR,Int_t fEntry, Int_t fChan,Double_t fEnr);
  virtual TH1D* GetHistoSmooth(TH1D* hist, Int_t DeltaBin);
  virtual TH1D* GetHistoDerivative(TH1D* hist, Int_t DeltaBin);
  virtual TH1D* TrapezoidalFilter(TH1D* hist, double RampTime, double FlatTime , double DecayTime);
  virtual TH1* GetHistoFFT(TH1D* hist);
  virtual Int_t FindPeaks(TH1D* hist, Double_t Low, Double_t Up, Double_t Resolution, Double_t Sigma, Double_t Threshold, vector<Double_t>* fPositionX, vector<Double_t>* fPositionY);
  virtual Double_t GetYValue(TH1D* hist, Double_t X);
  virtual Double_t GetMax(TH1* hist, Double_t Low, Double_t Up);
  virtual vector<Int_t> Sort(vector<Double_t> X);
  virtual vector<Int_t> Clean(vector<Double_t> X); 
  virtual void IsPileUpTag(Int_t fList, vector<Int_t>* IsPileUp, vector<Double_t>* Ratio, 
			   vector<Double_t>* DeltaT);
  virtual void PileUpTree(const char* pathName);
  virtual void TimeDiffTree(const char* pathName);
  virtual TH1D* FillHisto(TChain* mTree, string InputParaName, string OutputParaName,string CutName, Int_t Bin, Double_t Low, Double_t Up);

protected:
  Int_t fDataSet;
  Int_t fSubSet;
  Int_t fSubSet2;
  Int_t fIsCal;
  TChain* fSkimTree;
  size_t fEntries;
  Int_t fRun;
  Int_t fEvent;
  Int_t fStartRun;
  MJTChannelMap* fMap;
  vector<Int_t> fMTChannel;
  vector<Int_t> fCryo;
  vector<Int_t> fString;
  vector<Int_t> fDetector;
  vector<string> fDetectorName;
  vector<Int_t>* fChannel;  
  vector<Int_t>* fP;
  vector<Int_t>* fD;
  vector<Int_t>* fC;
  vector<Int_t>* fGain;
  vector<bool>* fIsEnr;
  vector<bool>* fIsNat;
  vector<bool>* fIsGood;
  vector<Double_t>* fTrapENFCal;
  vector<Double_t>* ftOffset;
  vector<Double_t>* fDCR;
  vector<Double_t>* fAvsE;
  vector<Double_t>* fTrapETailMin;
  Double_t flocalTime_s;
  Double_t fclockTime_s;
  TTimeStamp* fglobalTime;
  Int_t fmH ;
  vector<Int_t>* fnX;
  bool fmuVeto;
  vector<Double_t>* fdtmu_s;
  Double_t fmuTUnc;

private:


};

#endif
