// 2016.11.15 Pinghan Chu
// 1.This code automatically calibrates "trapE" or "trapENFxxx";not work on the onboard energy "energy". 
// 2.Before applying this code, need to login the database in order to upload the calibration parameters in the end. Ask Robert Varner(varnerrl@ornl.gov) for the password
// 
#include "GATAutoCal.hh"
#include "TFile.h"

void GetEZero(Int_t fStartRun, Int_t fEndRun, Int_t fChannel, string fEName, vector<Double_t> *Par){
  ifstream fin(Form("./List/ezero_%d_%d.txt",fStartRun,fEndRun));

  Int_t pos,cha,index;
  Double_t ref,peak,peakerr,fwhm,fwhmerr;
  string ename;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> ename >> index >> ref >> peak >> peakerr >> fwhm >> fwhmerr;
      if(cha == fChannel && ename == fEName){
        Par->push_back(peak);
        Par->push_back(peakerr);
        Par->push_back(fwhm);
        Par->push_back(fwhmerr);
      }
    }
    if(Par->size()==0){
      Par->push_back(0);
      Par->push_back(0);
      Par->push_back(0);
      Par->push_back(0);
    }
  }else{
    Par->push_back(0);
    Par->push_back(0);
    Par->push_back(0);
    Par->push_back(0);
  }
}


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
  //Double_t maxX = 1890;
  cout << "The maximum peak is " << maxX << endl;
  Int_t nPeak = Peak->size();
  for(size_t i = 0;i<Peak->size();i++){
    Peak->at(i) = Peak->at(i)*(maxX/Peak->at(nPeak-1));
  }
  delete h;
}

void IsExist(Int_t fStartRun,Int_t fEndRun,Int_t fChannel,string fEName,vector<Double_t>* Par){

  ifstream fin(Form("./List/cal/calibration_%d_%d.txt",fStartRun,fEndRun));

  Int_t cha,pos;
  Double_t r1,r2,r3,r4;
  Int_t run1,run2;
  string ename;
  
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> run1 >> run2 >> ename >> r1 >> r2 >>r3 >>r4;
      if(cha == fChannel && ename == fEName){
        Par->push_back(r1);
        Par->push_back(r2);
        Par->push_back(r3);
        Par->push_back(r4);
      }
    }
  }
}

TH1D* CleanHisto(TH1D* Input, Double_t Mean,Double_t Window){
  Int_t nbin = Input->GetNbinsX();
  TH1D *h = (TH1D*)Input->Clone();
  for(Int_t ib=0;ib<nbin;ib++){
    if((Input->GetBinCenter(ib)>Mean*(1-0.02) && Input->GetBinCenter(ib)<Mean*(1-0.007))){
      //if((Input->GetBinCenter(ib)> Mean-Window && Input->GetBinCenter(ib)<Mean-10*Window)){
      h->SetBinContent(ib,Input->GetBinContent(ib-1000));
      //h->SetBinContent(ib,0);
    }
    //if((Input->GetBinCenter(ib)> Mean+Window && Input->GetBinCenter(ib)<Mean+10*Window)){
    if((Input->GetBinCenter(ib)>Mean*(1+0.003) && Input->GetBinCenter(ib)<Mean*(1+0.02))){
      h->SetBinContent(ib,Input->GetBinContent(ib-10000));
      //h->SetBinContent(ib,0);
    }
  }
  return h;
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

  ofstream fcal(Form("./List/cal/calibration_%d_%d.txt",fStartRun,fEndRun),ios::app); // calibration constants
  ofstream fcov(Form("./List/cal/cov_%d_%d.txt",fStartRun,fEndRun),ios::app); // covariant
  ofstream fpeak(Form("./List/cal/ezfit_%d_%d.txt",fStartRun,fEndRun),ios::app); // peak information
  ofstream freso(Form("./List/cal/reso_%d_%d.txt",fStartRun,fEndRun),ios::app); // ezero peak
  fcal.precision(12);
  fcov.precision(12);
  fpeak.precision(12);
  freso.precision(12);
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

  vector<string> DataName;
  //DataName.push_back("00"); // Natural HG
  //DataName.push_back("01"); // Natural LG
  //DataName.push_back("10"); // Enriched HG
  //DataName.push_back("11"); // Enriched LG
  DataName.push_back("HG");
  //DataName.push_back("LG");
  vector<string> TitleName;
  //TitleName.push_back("Natural HG");
  //TitleName.push_back("Natural LG");
  //TitleName.push_back("Enriched HG");
  //TitleName.push_back("Enrichd LG");
  TitleName.push_back("High Gain");
  //TitleName.push_back("Low Gain");



  Double_t Up = 2700;
  Double_t Low = 2500;
  Int_t rebin = 1;
  Double_t ROIEnergy = 2039;
  string Pos;
  string HistoName;
  string PlotName;
  //string TitleName;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> Mu;
  vector<Double_t> MuErr;
  vector<Double_t> FWHM;
  vector<Double_t> FWHMErr;
  vector<Double_t> cov;

  for(size_t i = 0;i<DataName.size();i++){
   
    Low = 2600;
    Up = 2630;
    rebin = 2;
	
    HistoName = Form("%s%s",fEName.c_str(),DataName.at(i).c_str());
	
    TH1D *h1 = GetHisto(FileName,HistoName);
    TH1D *h2 = GetHisto(FileName,HistoName);
    Mu.clear();
    MuErr.clear();
    FWHM.clear();
    FWHMErr.clear();
    fPeak.clear();
    fPeak = ds.GetCalibrationPeak();
    //fPeak.erase(fPeak.begin()+2);
    h1->Rebin(rebin);
    GetPeak(h1,Low, Up,&fPeak);
    cout << fPeak.at(fPeak.size()-1) << endl;
    fCalibrationPeak.clear();
    fCalibrationPeak = ds.GetCalibrationPeak();  
    //fCalibrationPeak.erase(fCalibrationPeak.begin()+2);
    fCalibrationPeakErr.clear();
    
    for(size_t ip = 0;ip<fPeak.size();ip++){
      if(ip>1){
	rebin = rebin*2;
      }
      
      cout << ip << " " << fPeak.at(ip) << endl;
      
      PlotName = Form("%s_%d_%d_%s_%d",fEName.c_str(),fStartRun,fEndRun,DataName.at(i).c_str(),(Int_t)ip);
      
      ds.EzFit(h2, fPeak.at(ip), 10, rebin, &Par,&ParErr, "", PlotName,TitleName.at(i));
      fpeak << "000" << " " << DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "
	    << ip << " " << fCalibrationPeak.at(ip) << " "
	    << Par.at(1) << " " << ParErr.at(1) << " " << Par.at(12) << " " << ParErr.at(12) << endl;
      Mu.push_back(Par.at(1));
      MuErr.push_back(ParErr.at(1));
      FWHM.push_back(Par.at(12));
      FWHMErr.push_back(ParErr.at(12));
      fCalibrationPeakErr.push_back(0);
    }
    PlotName = Form("%s_%d_%d_%s",fEName.c_str(),fStartRun,fEndRun,DataName.at(i).c_str());
    ROIEnergy = 2039;
    Par.clear();
    ParErr.clear();
    ds.LinearFit(Mu,MuErr,fCalibrationPeak,fCalibrationPeakErr,&Par,&ParErr,"",PlotName,TitleName.at(i),ROIEnergy);
    fcal << "000" << " " << DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "
	 << Par.at(0) << " " << ParErr.at(0) << " " << Par.at(1) << " " << ParErr.at(1) << endl;
    ROIEnergy = Par.at(2);
    cov.clear();
    for(Int_t k = 0;k<2;k++){
      for(Int_t j = 0;j<2;j++){
	Int_t l = k*2+j+3;
	cov.push_back(ParErr.at(l));
      }
    }
    fcov <<"000" << " " << DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "
	 << cov.at(0) << " " << cov.at(1) << " " << cov.at(2) << " " << cov.at(3) << endl;
    
    Par.clear();
    ParErr.clear();
    ds.ResolutionFit(Mu, MuErr,FWHM, FWHMErr,&Par,&ParErr,"",PlotName,TitleName.at(i),ROIEnergy);
    freso <<"000" << " " << DataName.at(i).c_str() << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "
	  << Par.at(0) << " " << ParErr.at(0) << " " << Par.at(1) << " " << ParErr.at(1) << " " << Par.at(2) << " " << ParErr.at(2) << " "
	  << Par.at(3) << " " << ParErr.at(3) << endl;
  }
}


