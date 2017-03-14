// 2016.11.15 Pinghan Chu
// 1.This code automatically calibrates "trapE" or "trapENFxxx";not work on the onboard energy "energy". 
// 2.Before applying this code, need to login the database in order to upload the calibration parameters in the end. Ask Robert Varner(varnerrl@ornl.gov) for the password
// 
#include "GATAutoCal.hh"
#include "TFile.h"


TH1D* GetHisto(string FileName, string HistoName){
  TFile *fhist = TFile::Open(Form("%s",FileName.c_str()),"read");
  TH1D *h = (TH1D*)fhist->Get(Form("%s",HistoName.c_str()));
  return h;
}


void GetPeak(TH1D *Input, Double_t Low, Double_t Up, vector<Double_t>* Peak){
  //Locate the 2614 keV peak in the raw spectrum.

  TH1D *h = (TH1D*)Input->Clone();
  h->GetXaxis()->SetRangeUser(Low,Up);
  Int_t maxbin = h->GetMaximumBin();
  Double_t maxX = h->GetBinCenter(maxbin);
  cout << "The maximum peak is " << maxX << endl;
  Int_t nPeak = Peak->size();
  for(size_t i = 0;i<Peak->size();i++){
    Peak->at(i) = Peak->at(i)*(maxX/Peak->at(nPeak-1));
  }
  delete h;
}

int main(int argc, char** argv)
{
  if(argc != 5 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [start run] [endrun] [energy name] [FileName]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  string fEName = argv[3];
  const char* EnergyName = Form("%s",fEName.c_str());
  string FileName = argv[4];

  ofstream ferror(Form("./List/cal/error_%d_%d.txt",fStartRun,fEndRun),ios::app); // calibration constants
  GATAutoCal ds(fStartRun,fStartRun);
  ds.SetEnergyName(EnergyName);
  string fDataSet = ds.GetDataSet();  
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr  = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  vector<Double_t> fCalibrationPeak;
  vector<Double_t> fCalibrationPeakErr;
  vector<Double_t> fPeak = ds.GetCalibrationPeak();
  vector<Int_t> index;
  Int_t pos = 0;
  for(size_t i = 0;i<fChannel.size();i++){
    if(fGoodBad.at(i)==1){
      index.push_back(i);
    }
  }
  Double_t Up = 2700;
  Double_t Low = 2500;
  Int_t rebin = 20;
  Double_t ROIEnergy = 2039;
  string Pos;
  string HistoName;
  string PlotName;
  string TitleName;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> Mu;
  vector<Double_t> MuErr;
  vector<Double_t> FWHM;
  vector<Double_t> FWHMErr;
  vector<Double_t> cov;

  for(size_t i = 0;i<index.size();i++){
    Int_t ii = index.at(i);
    Pos = Form("%d%d%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii));
    HistoName = Form("%s%s%d",fEName.c_str(),Pos.c_str(),fChannel.at(ii));
    TH1D *h1 = GetHisto(FileName,HistoName);

    fPeak.clear();
    fPeak = ds.GetCalibrationPeak();
    h1->Rebin(rebin);
    GetPeak(h1,Low, Up,&fPeak);
    Int_t pars = fPeak.size();
    if(abs(fPeak.at(pars-1)-2614.5)>3){
      ferror << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " 
	     <<fEName.c_str()<< " " << fPeak.at(pars-1) << " " << abs(fPeak.at(pars-1)-2614.5)<<endl;
    }
    delete h1;
  }
}


