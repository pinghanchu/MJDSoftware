#ifndef GATAutoCal_hh
#define GATAutoCal_hh

#include "GATDataSet.hh"
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
/*
Double_t Reso(Double_t *v, Double_t *par){
  Double_t fitval = TMath::Sqrt(par[0]*par[0]+par[1]*par[1]*v[0]+par[2]*par[2]*v[0]*v[0]);
  return fitval;
}
*/

class GATAutoCal : public TObject  
{
public:   
  GATAutoCal(Int_t StartRun = 1, Int_t EndRun = 1);
  virtual ~GATAutoCal() {}

  virtual inline Int_t GetStartRun() { return fStartRun; }
  virtual inline Int_t GetEndRun() { return fEndRun; }
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
  virtual inline vector<Double_t> GetCalibrationPeak(){ return fCalibrationPeak; }
  virtual inline vector<string> GetDetectorName(){ return fDetectorName;}
  virtual inline vector<Double_t> GetPulserCal(){ return fPulserCal;}
  virtual inline vector<Double_t> GetPulser(){ return fPulser;}
  virtual inline void SetEnergyName(const char* energyname) { fEnergyName = energyname; }
  virtual inline Double_t GetTotalTime() { return fDS.GetRunTime(); }
  virtual inline TChain* GetMJDTree(){ return fMjdTree; }
  virtual inline MJTChannelMap* GetMap(){ return fMap; }
  virtual inline TH1D* GetHisto(Int_t index){ return TrapE[index]; }

  virtual void SetMjdTree(TChain *mjdTree);
  virtual void SetRunRange(Int_t StartRun, Int_t EndRun);
  virtual void SetCalibrationPeak();
  virtual void SetCalibrationPeak(vector<Double_t> gammakeV);
  virtual void SetDataSet();
  virtual void SetParameters();
  virtual void SetChannel();
  virtual Int_t GetGATRev();
  virtual Double_t GetStartTime();
  virtual Double_t GetStartTimeMT();
  virtual Int_t GetStartDateMT();
  virtual Double_t GetStopTime();
  virtual Double_t GetStopTimeMT();
  virtual Int_t GetStopDateMT();
 
  virtual void FillHisto(TChain* mTree,Int_t Bin1, Double_t Low1, Double_t Up1, Int_t Bin2, Double_t Low2, Double_t Up2);
  virtual void FillHistoCal(TChain* mTree,vector<Double_t> Offset, vector<Double_t> Slope,Int_t Bin1, Double_t Low1, Double_t Up1, Int_t Bin2, Double_t Low2, Double_t Up2);
  virtual void FillHistoDC(TChain* mTree, vector<Int_t> Channel, vector<string> EnergyName, Int_t Bin1, Double_t Low1, Double_t Up1, Int_t Bin2, Double_t Low2, Double_t Up2);
  virtual void FillHistoDCCal(TChain* mTree,vector<string> EnergyName,vector<Double_t> Offset, vector<Double_t> Slope,Int_t Bin1, Double_t Low1, Double_t Up1, Int_t Bin2, Double_t Low2, Double_t Up2);
  virtual void SaveFile(string PathName,string FileName);
  virtual void LoadFile(string PathName,string FileName);

  virtual Double_t GetMaximumPeak(TH1D *Hist,Double_t Low, Double_t Up);
  virtual vector<Double_t> GetRefPeak(Double_t LastPeak, Double_t FirstPeak);

  virtual void GaussFit(TH1D *Hist, Double_t Mean, Double_t Window, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName,string FileName, string TitleName);
  virtual void SkewGaussFit(TH1D *Hist, Double_t Mean, Double_t Window, vector<Double_t>* Par, vector<Double_t>* ParErr,string PathName, string FileName, string TitleName, string opt);
  virtual void GaussSlopeFit(TH1D *Hist, Double_t Mean, Double_t Window,vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName,string FileName, string TitleName);
  virtual void EzFit(TH1D *Hist, Double_t Mean, Double_t Window,Int_t ReBin, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName,string FileName, string TitleName);
  virtual void NoBGFit(TH1D *Hist, Double_t Mean, Double_t Window,Int_t ReBin, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName,string FileName, string TitleName);
  virtual void MultiPeakFitterFull(TH1F *Hist, Double_t scaleenergy, Double_t scaleamps,string FileName, vector<string>* ParName, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov);
  virtual Int_t MultiPeakFitter(TH1F *Hist,Double_t scaleenergy, Double_t scaleamps,string FileName, vector<string>* ParName, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov);
  virtual Int_t MultiPeakFitter(TH1F *Hist,Double_t scaleenergy, Double_t scaleamps,string FileName, vector<string>* ParName, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov, TH1F* HistoE0);
  virtual void LinearFit(vector<Double_t> Px, vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, vector<Double_t>* Par, vector<Double_t>* ParErr,string PathName,string FileName,string TitleName, Double_t ROIEnergy);

  virtual void ResolutionFit(vector<Double_t> Px, vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName, string FileName, string TitleName, Double_t ROIEnergy);

  virtual TGraphErrors* QuadFit(vector<Double_t> Px,vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName, string FileName, string TitleName);

  virtual Double_t RampTimeFit(vector<Double_t> Px,vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName, string FileName, string TitleName);

  virtual void PlotMultiGraph(TMultiGraph *MG, TLegend *Leg,string PathName, string FileName, string TitleName);

  virtual void PlotGrid(vector<Int_t> Px,vector<Double_t> Py, vector<Double_t> PyErr, Double_t Ave, string PathName,string FileName,string TitleName);
  virtual void PlotGrid(TMultiGraph *MG, TLegend *Leg,vector<Int_t> Px, vector<Double_t> Ave, Double_t Low, Double_t Up, string PathName,string FileName);

  virtual void PlotSpectrum(TH1D *Hist, std::string FileName);
  virtual void PlotSpectrum2(TH1D *H1,TH1D *H2,TLegend *Leg, std::string FileName);
  virtual void PlotGraph(TGraphErrors *Graph, string FileName);
  //virtual void PlotMultiGraph(TMultiGraph *Graph, TLegend *Leg, string FileName);

  virtual void SetUpProvenance(std::string title, std::string yourname,
		       Int_t gatrev,
		       int startrun, int endrun,
		       int coverstartrun, int coverendrun,
		       int startdate, int enddate,
		       int coverstartdate, int coverenddate);
    
  virtual void PutECalMJDB(int channel,std::string detectorid,
		   double scale, double scalerr,
		   double offset, double offerr,
		   std::vector<double> covariance);


  virtual TProfile* PeakProfileTime(TChain* fChain,string fCut,Double_t PeakChannel,Double_t PeakChannelWindow, string PathName,string FileName);

  virtual void PutPeakLocationMJDB(Int_t channel,std::string detectorid,
			   vector<Double_t> &calpeak,
			   vector<Double_t> &calpeakerr,
			   vector<Double_t> &calwidth,
			   vector<Double_t> &calwidtherr);
  //   int &calibrationpeaks);
  virtual void GetPeakLocationMJDB(Int_t channel,
			   vector<Double_t> &calpeak,
			   vector<Double_t> &calpeakerr,
			   vector<Double_t> &calwidth,
			   vector<Double_t> &calwidtherr,
			   int &calibrationpeaks);


protected:
  GATDataSet fDS;
  TChain* fMjdTree;
  TChain* fMjdTree1;
  MJTChannelMap* fMap;
  TH1D *TrapE[125];

  size_t fEntries;
  Int_t fStartRun;
  Int_t fEndRun;
  UInt_t fGATRev;
  Int_t fIsRadio;
  Double_t fMTStartTime;
  Double_t fMTStopTime;
  string fEnergyName;
  string fDataSet;

  vector<Int_t> fChannel;
  vector<Int_t> fPulserTagChannel;
  vector<Int_t> fCryo;
  vector<Int_t> fString;
  vector<Int_t> fDetector;
  vector<Int_t> fGoodBad;
  vector<Int_t> fEnriched;
  vector<Double_t> fPulserCal;
  vector<Double_t> fPulser;
  vector<Double_t> fCalibrationPeak;
  vector<string> fDetectorName;
  vector<Double_t>* fMTTime;
  vector<Int_t>* fMTDate;
  //std::vector<Double_t>* fMTChannel;
  //std::vector<Double_t>* fMTTrapENFCal;
  //std::vector<Double_t>* fMTTrapECal;
  //std::vector<Double_t>* fMTTimestamp;
  //std::vector<Int_t>* fMTwfDCBits;

  Int_t fChannels;
  Int_t fStrings;
  Int_t fCalibrationPeaks;
  Int_t fPulserTagChannels;
  Int_t fTotalChannels;
private:
  MJDB::MJProvenance fACProvenance;
  // use this one copy of the analysisdb
  MJDB::MJAnalysisDoc fACAnalysisDB;
  //ClassDef(GATAutoCal,1)
};

#endif
