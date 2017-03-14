// 2016.11.15 Pinghan Chu
// 1.This code automatically calibrates "trapE" or "trapENFxxx";not work on the onboard energy "energy". 
// 2.Before applying this code, need to login the database in order to upload the calibration parameters in the end. Ask Robert Varner(varnerrl@ornl.gov) for the password
// 
#include "GATAutoCal.hh"
#include "TFile.h"

void GetEZero(Int_t fStartRun, Int_t fEndRun, Int_t fChannel, string fEName, vector<Double_t> *Par){
  ifstream fin(Form("./List/cal/ezero_%d_%d.txt",fStartRun,fEndRun));

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
  Double_t minX = 536;
  /*
  for(size_t i = 0;i<Peak->size();i++){
    Peak->at(i) = (Peak->at(i)-238)*((maxX-minX)/(2614.5-238))+minX;
    cout << Peak->at(i) << endl;
  }
  */
  Peak->clear();
  Peak->push_back(minX);
  Peak->push_back(1312);
  Peak->push_back(1638);
  Peak->push_back(5645);
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
void GetParameters(string fEName, string fDataSet, Int_t fStartRun, Int_t fEndRun, Int_t fChannel, vector<Double_t>* Par){
  Double_t Low;
  Double_t Up;
  Double_t rebin;
  if(fEName == "trapENFCal" || fEName == "trapECal" || fEName == "trapENMCal"){
    Low = 2000;
    Up = 3000;
    rebin = 20;
  }else{
    if((fDataSet == "P3END" && fChannel == 145) || (fStartRun >= 45007013 && fEndRun<=45499999 && fChannel == 149)){
      Low = 700;
      Up = 1500;
      rebin = 5;
    }else if(fDataSet == "P3JDY" && fChannel == 663){
      Low = 700;
      Up = 1000;
      rebin = 5;
    }else if(fDataSet == "P3KJR" && fChannel == 694 && fStartRun < 16797){
      Low = 1500;
      Up = 1800;
      rebin = 2;
    }else if(fDataSet == "P3KJR" && fChannel == 695 && fStartRun< 16797){
      Low = 400;
      Up = 600;
      rebin = 1;
    }else if(fDataSet == "P3LQK" && fChannel == 616){
      Low = 700;
      Up = 1000;
      rebin = 5;
    }else if(fDataSet == "P3LQK" && fChannel == 617){
      Low = 200;
      Up = 250;
      rebin = 2;
    }else if(fChannel%2==0){
      Low = 5000;
      Up = 6000;
      rebin = 50;
    }else if(fChannel%2==1){
      Low = 1500;
      Up = 2500;
      rebin = 20;
    }
  }
  Par->push_back(Low);
  Par->push_back(Up);
  Par->push_back(rebin);
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
    if((Input->GetBinCenter(ib)>Mean*(1+0.007) && Input->GetBinCenter(ib)<Mean*(1+0.02))){
      h->SetBinContent(ib,Input->GetBinContent(ib-10000));
      //h->SetBinContent(ib,0);
    }
  }
  return h;
}

int main(int argc, char** argv)
{
  if(argc != 8 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [start run] [endrun] [cover start run] [cover end run] [energy name] [channel] [FileName]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  Int_t fCoverStartRun = atoi(argv[3]);
  Int_t fCoverEndRun = atoi(argv[4]);
  string fEName = argv[5];
  const char* EnergyName = Form("%s",fEName.c_str());
  Int_t fChan = atoi(argv[6]);
  string FileName = argv[7];

  Int_t startrun1=4171;
  Int_t endrun1 = 4201;

  if(fStartRun >=2361 && fEndRun<=8183){
    startrun1 = 4171;
    endrun1 = 4201;
  }else if(fStartRun >=9034 && fEndRun<=13070){
    startrun1 = 9913;
    endrun1 = 9926;
  }else if(fStartRun>=13071 && fEndRun<=18589){
    startrun1 = 13071;
    endrun1 = 13074;
  }else if(fStartRun>=60000000 && fEndRun<65000000){
    startrun1 = 60002368;
    endrun1 = 60002371;
  }else if(fStartRun>=18623 && fEndRun<50000000){
    startrun1 = 20000;
    endrun1 = 20000;
  }


  ofstream fcal(Form("./List/cal/calibration_%d_%d.txt",fStartRun,fEndRun),ios::app); // calibration constants
  ofstream fcal1("./List/calibration.txt",ios::app); // calibration constants
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
  for(size_t i = 0;i<fChannel.size();i++){

    if(fChannel.at(i) == fChan){
      index.push_back(i);
    }
  }
  Double_t Up = 2700;
  Double_t Low = 2500;
  Int_t rebin = 1;
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
    Par.clear();
    //IsExist(fStartRun,fEndRun,fChannel.at(ii),fEName,&Par);
    //if(Par.size()>=0){
      if(fGoodBad.at(ii)==1){
	cout << Pos.c_str() << " " << fChannel.at(ii) << " "<< fGoodBad.at(ii) << endl;
	Par.clear();
	GetParameters(fEName,fDataSet,fStartRun,fEndRun,fChannel.at(ii),&Par);
	Low = Par.at(0);
	Up = Par.at(1);
	rebin = (Int_t)Par.at(2);
	
	HistoName = Form("%s%s%d",fEName.c_str(),Pos.c_str(),fChannel.at(ii));
	
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
	  
	  PlotName = Form("%s_%d_%d_%s_%d_%d",fEName.c_str(),fStartRun,fEndRun,Pos.c_str(),fChannel.at(ii),(Int_t)ip);
	  TitleName = Form("C%dP%dD%d, channel %d",fCryo.at(ii),fStr.at(ii), fDetpos.at(ii), fChannel.at(ii));      
	  
	  //TH1D *h3 = CleanHisto(h2, fPeak.at(ip),5);
	  Par.clear();
	  ParErr.clear();	  
	  ds.EzFit(h2, fPeak.at(ip),30, rebin, &Par,&ParErr, "", PlotName,TitleName);
	  fpeak << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "
		<< ip << " " << fCalibrationPeak.at(ip) << " "
		<< Par.at(1) << " " << ParErr.at(1) << " " << Par.at(23) << " " << ParErr.at(12) << endl;
	  Mu.push_back(Par.at(1));
	  MuErr.push_back(ParErr.at(1));
	  FWHM.push_back(Par.at(23));
	  FWHMErr.push_back(ParErr.at(12));
	  fCalibrationPeakErr.push_back(0);
	  //delete h3;
	}
	
	if(fEName == "trapENF"){
	  Pos = Form("%d%d%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii));
	  HistoName = Form("%sBL%s%d",fEName.c_str(),Pos.c_str(),fChannel.at(ii));      
	  PlotName = Form("%s_%d_%d_%s_%d_%d",fEName.c_str(),fStartRun,fEndRun,Pos.c_str(),fChannel.at(ii),(Int_t)fPeak.size());
	  TitleName = Form("C%dP%dD%d, channel %d",fCryo.at(ii),fStr.at(ii), fDetpos.at(ii), fChannel.at(ii));
	  TH1D *h3 = GetHisto(FileName,HistoName);  
	  //h3->Rebin(1);
          //Int_t maxbin = h3->GetMaximumBin();
	  //Double_t maxX = h3->GetBinCenter(maxbin);
	  //cout << h3->GetRMS() << " "<< h3->GetMaximum() << " " << maxX << endl;
	  
	  if(h3->GetRMS()<1.2 && h3->GetEffectiveEntries()>10){
	    Par.clear();
	    ParErr.clear();
	    //TH1D *h4 = CleanHisto(h3,0,0.5);
	    //h4->Rebin(10);

	    ds.SkewGaussFit(h3,0,2,&Par,&ParErr,"",PlotName,TitleName,"");
	    if(abs(Par.at(1))<0.1){
	      Mu.push_back(Par.at(1));
	      MuErr.push_back(ParErr.at(1));
	      FWHM.push_back(Par.at(4));
	      FWHMErr.push_back(ParErr.at(4));
	      fpeak << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "
		    << fPeak.size() << " " << 0 << " "
		    << Par.at(1) << " " << ParErr.at(1) << " " << Par.at(4) << " " << ParErr.at(4) << endl;
	    }else{
	      Mu.push_back(0);
	      MuErr.push_back(0.02);
	      FWHM.push_back(0);
	      FWHMErr.push_back(0.02);
	    }
	    fCalibrationPeak.push_back(0);
	    fCalibrationPeakErr.push_back(0);
	  }else{
	    Mu.push_back(0);
	    MuErr.push_back(0.005);
	    FWHM.push_back(0);
            FWHMErr.push_back(0.1);
            fCalibrationPeak.push_back(0);
            fCalibrationPeakErr.push_back(0);
	  }
	  delete h3;
	}else if(fEName == "trapENM" || fEName == "trapE"){
	  /*
	  Par.clear();
          GetEZero(startrun1,endrun1,fChannel.at(ii),fEName, &Par);
          if(Par.at(0)!=0){
            Mu.push_back(Par.at(0));
            MuErr.push_back(Par.at(1));
            FWHM.push_back(Par.at(2));
            FWHMErr.push_back(Par.at(3));
            fCalibrationPeak.push_back(0);
            fCalibrationPeakErr.push_back(0);
          }
	  */
	}
	
	PlotName = Form("%s_%d_%d_%s_%d",fEName.c_str(),fStartRun,fEndRun,Pos.c_str(),fChannel.at(ii));
	TitleName = Form("C%dP%dD%d, channel %d",fCryo.at(ii),fStr.at(ii), fDetpos.at(ii), fChannel.at(ii));
	ROIEnergy = 2039;
	Par.clear();
	ParErr.clear();
	ds.LinearFit(Mu,MuErr,fCalibrationPeak,fCalibrationPeakErr,&Par,&ParErr,"",PlotName,TitleName,ROIEnergy);
	fcal << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "
	     << Par.at(0) << " " << ParErr.at(0) << " " << Par.at(1) << " " << ParErr.at(1) << endl;      
        fcal1 << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "
             << Par.at(0) << " " << ParErr.at(0) << " " << Par.at(1) << " " << ParErr.at(1) << endl;
	ROIEnergy = Par.at(2);
	cov.clear();
	for(Int_t k = 0;k<2;k++){
	  for(Int_t j = 0;j<2;j++){
	    Int_t l = k*2+j+3;
	    cov.push_back(ParErr.at(l));
	  }
	}
	fcov <<Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " 
	     << cov.at(0) << " " << cov.at(1) << " " << cov.at(2) << " " << cov.at(3) << endl;
	
	Par.clear();
	ParErr.clear();
	ds.ResolutionFit(Mu, MuErr,FWHM, FWHMErr,&Par,&ParErr,"",PlotName,TitleName,ROIEnergy);
	freso << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "
	      << Par.at(0) << " " << ParErr.at(0) << " " << Par.at(1) << " " << ParErr.at(1) << " " << Par.at(2) << " " << ParErr.at(2) << " "
	      << Par.at(3) << " " << ParErr.at(3) << endl;
	
	delete h1;
	delete h2;
      }else{
	Par.clear();
	fcal << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " 
	     << fEName.c_str() << " 0 0 999999 0 "
	     << endl;
        fcal1 << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " 
	      << fEName.c_str() << " 0 0 999999 9 " <<endl;
             
	fcov << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " 0 0 0 0 "
	     << endl;
	freso << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " 0 0 0 0 0 0 0 0 "
	      << endl;
      }
      //    }
  }
}


