#ifndef _GATPulserTag_hh
#define _GATPulserTag_hh

#include "GATDataSet.hh"
#include "GATTimeInfo.hh"
#include "TObject.h"
#include "TChain.h"
#include <vector>
#include <string>
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

class GATPulserTag : public TObject  
{
public:   
  GATPulserTag(Int_t run = 1);
  virtual ~GATPulserTag() {}

  virtual inline Int_t GetRun() { return fRun; }
  virtual inline Int_t GetEntries() { return fEntries; }
  virtual inline TChain* GetMJDTree(){ return fMjdTree; }

  virtual Int_t GetIsRadio() { return fIsRadio; }
  virtual inline Double_t GetTotalTime() { return fDS.GetRunTime()/CLHEP::second; }

  virtual void SetParameters();
  virtual inline void SetPulserWindow(Double_t pulserwindow) { fPulserWindow = pulserwindow; }
  virtual inline void SetPulserTimeWindow(Double_t pulsertimewindow) { fPulserTimeWindow = pulsertimewindow; }
  virtual inline void SetEnergyName(const char* energyname) { fEnergyName = energyname; }
  virtual void SetDataSet();
  virtual void SetChannel();
  virtual void SetPulser();
  virtual vector<Int_t> GetChannel(){ return fChannel;}
  virtual string GetDataSet(){ return fDataSet;}
  virtual void GetPulserPeakFromDB(Int_t LrunNumber, Int_t UrunNumber, Int_t chanNumber, Double_t* par, Double_t* parerr);

  // first pulser filter, energy, time window
  virtual void Tag1(const std::vector<Int_t>& inputEntry, std::vector<Int_t>& outputEntry); 

  // second pulser filter, pulser monitor channel
  virtual void Tag2(const std::vector<Int_t>& inputEntry, std::vector<Int_t>& outputEntry); 

  // multiplicity 
  virtual void Tag3(const std::vector<Int_t>& inputEntry, std::vector<Int_t>& outputEntry); 

  // get physical event list (without correction)
  virtual void Tag(const std::vector<Int_t>& inputEntry, std::vector<Int_t>& outputEntry); 

  // add addition list
  virtual void TagCorrect(std::vector<Int_t>& inputList, const std::vector<Int_t>& correctList); 

  virtual void GetTagList(const std::vector<Int_t>& inputEntry, std::vector<Int_t>& outputEntry);

  virtual void GetPulserList(const std::vector<Int_t>& inputEntry, std::vector<Int_t>& PList);

  virtual void GetPulserChannelList(const std::vector<Int_t>& inputEntry, 
                                    Int_t index, 
                                    std::vector<Int_t>& PList);

  virtual void GetPulserTime(const std::vector<Int_t>& inputEntry, 
                             Int_t index, 
                             std::vector<Double_t>& pTime);

  virtual Double_t GetPulserDiffTime(std::vector<Double_t> PulserTime);

  virtual void GetDeviateEvent(const std::vector<Int_t>& inputEntry, 
                               const std::vector<Double_t>& PulserTime,
                               Double_t PulserDiffTime,
                               std::vector<Int_t>& DeviateList);

  virtual void PhysicalEventList(std::vector<Int_t>& peList); 

  virtual void PulserTree(const char* pathName);
  virtual void SaveTree(TChain* fChain, Int_t IsPulser,string PathName, string FileName, string PulserPath);

protected:
  GATDataSet fDS;
  Int_t fRun;
  Int_t fMTRun;
  Int_t fmH;
  Int_t fmL;
  size_t fEntries;
  Int_t fChannels;
  std::string fEnergyName;
  std::string fDataSet;
  Double_t fPulserWindow;
  Double_t fPulserTimeWindow;
  std::vector<Int_t> fChannel;
  std::vector<Int_t> fPulserTagChannel;
  std::vector<Int_t> fGoodBad;
  std::vector<Double_t> fPulser;
  //std::vector<Double_t> fPulserCal;
  Int_t fIsRadio;
  Int_t fStrings;
  MJTChannelMap* fMap;
  TChain* fMjdTree;
  std::vector<Int_t>* fMTC;
  std::vector<Int_t>* fMTP;
  std::vector<Int_t>* fMTD;
  std::vector<double>* fMTChannel;
  std::vector<double>* fMTTrapE;
  GATTimeInfo* fTimeInfo;

  //ClassDef(GATPulserTag,1)
};

#endif
