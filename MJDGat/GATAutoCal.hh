#ifndef _GATAutoCal_hh
#define _GATAutoCal_hh

#include "GATDataSet.hh"
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

class GATAutoCal : public TObject  
{
public:   
  GATAutoCal(Int_t StartRun = 1, Int_t EndRun = 1);
  virtual ~GATAutoCal();

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

  virtual Double_t IsNan(string Input);
  virtual void SetMjdTree(TChain *mjdTree);
  virtual void SetRunRange(Int_t StartRun, Int_t EndRun);
  virtual void SetCalibrationPeak();
  virtual void SetCalibrationPeak(vector<Double_t> gammakeV);
  virtual void SetDataSet();
  virtual void SetParameters();
  virtual void SetChannel();
  virtual Int_t GetGATRev();
  virtual Double_t GetStartTime();
  virtual Double_t GetStopTime();
 
  virtual TH1D* FillHisto(TChain* mTree, string EnergyName, Int_t Channel, Int_t Bin, Double_t Low, Double_t Up);
  virtual TH1D* FillPulser(TChain* mTree, string EnergyName, Int_t Channel, Int_t Bin, Double_t Low, Double_t Up);
  virtual void SaveHisto(TH1D* Hist, string FileName, string Option);
  virtual TH1D* LoadHisto(string FileName, string HistoName);

  virtual Int_t MultiPeakFit(TH1F *Hist,Double_t scaleenergy, Double_t scaleamps,string FitName, vector<string>* ParName, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov, vector<Double_t>* CalPeak);
  virtual Int_t LinearFit(vector<Double_t> Px, vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, string FitName, string TitleName, Double_t EnergyROI, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov, vector<Int_t>* FitIndex);
  virtual void PlotResolution(vector<Double_t> Par, string FileName);
  virtual void PlotGrid(TMultiGraph *MG, TLegend *Leg,vector<string> Px, vector<Double_t> Ave, Double_t Low, Double_t Up, string FileName);
  virtual void SetUpProvenance(std::string Title, std::string YourName,
			       Int_t gatrev,
			       int startrun, int endrun,
			       int coverstartrun, int coverendrun,
			       int startdate, int enddate,
			       int coverstartdate, int coverenddate);

  virtual void PutECalMJDB(int Channel,std::string DetectorName,
			   double scale, double scalerr,
			   double offset, double offerr,
			   std::vector<double> covariance);
  

  virtual void PutPulserMJDB(std::string Title, std::string YourName,
			     Int_t gatrev,Int_t Run, 
			     Int_t startdate, Int_t enddate,
			     Int_t Channel, std::string DetectorName,
			     Double_t Pulser, Double_t PulserErr);
  /*
  virtual void PutMultiPeakFitMJDB(std::string Title, std::string YourName,
			      Int_t gatrev,Int_t StartRun, Int_t EndRun, 
			      Int_t startdate, Int_t enddate,
			      Int_t Channel, std::string DetectorName,
			      vector<string> ParName,vector<Double_t> PeakFit, vector<Double_t> PeakFitErr);
  */
protected:
  GATDataSet fDS;
  TChain* fMjdTree;
  MJTChannelMap* fMap;
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
  GATTimeInfo* fTimeInfo;
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
