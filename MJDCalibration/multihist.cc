// 2016.11.21 Pinghan Chu
// calibrate sum spectrum of all high gain, all natural detectors etc.
#include "GATAutoCal.hh"
#include "TMatrixD.h"
#include "TFile.h"


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
    }else if(fDataSet == "P3LQG" && fChannel == 1136){
      Low = 400;
      Up = 1000;
      rebin = 5;
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
      Up = 8000;
      rebin = 20;
    }else if(fChannel%2==1){
      Low = 1500;
      Up = 2500;
      rebin = 10;
    }
  }
  Par->push_back(Low);
  Par->push_back(Up);
  Par->push_back(rebin);
}


int main(int argc, char** argv)
{
  if(argc != 5 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [start run] [endrun] [energy name] [position Mod*10+String]" << endl;
    return 1;
  }
  
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  string fEName = argv[3];
  Int_t fPos = atoi(argv[4]);
  const char* EnergyName = Form("%s",fEName.c_str());


  ofstream fout(Form("./List/cal/calibrationlog_%d_%d.txt",fStartRun,fEndRun),ios::app);
  ofstream fcal(Form("./List/cal/calibration_%d_%d.txt",fStartRun,fEndRun),ios::app);
  ofstream fcov(Form("./List/cal/cov_%d_%d.txt",fStartRun,fEndRun),ios::app);
  ofstream freso(Form("./List/cal/reso_%d_%d.txt",fStartRun,fEndRun),ios::app);
  ofstream froi(Form("./List/cal/roi_%d_%d.txt",fStartRun,fEndRun),ios::app);
 
  fcal.precision(15);
  fcov.precision(15);
  freso.precision(15);
  froi.precision(15);
  GATAutoCal ds(fStartRun,fStartRun);
  string fDataSet = ds.GetDataSet();
  ds.SetEnergyName(EnergyName);
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr  = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  vector<Int_t> index;
  const Int_t nChannel = fChannel.size();
  string Pos;
  Int_t pos;
  for(size_t i = 0;i<fChannel.size();i++){
    pos = fCryo.at(i)*10+fStr.at(i);
    //if(fChannel.at(i)==fPos){
    if(pos == fPos){
      index.push_back(i);
    }
  }

  TH1D *Energy[nChannel];
  Double_t Low, Up;
  Int_t rebin = 1;
  string FileName;

  Double_t fScale = 1;
  Double_t fOffset = 0;
  Double_t fScaleErr = 0;
  Double_t fOffsetErr = 0;
  Double_t P0=0;
  Double_t P1=0;
  Double_t P2=0;
  Double_t P0Err=0;
  Double_t P1Err=0;
  Double_t P2Err=0;
  Double_t errsum = 0;
  Int_t Index = 0;
  vector<string> ParName;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> Cov;
  vector<Double_t> cov;

  TFile fhist(Form("./Hist/hist_%d_%d_%d.root",fStartRun,fEndRun,fPos),"read");
 
  for(size_t i = 0;i<index.size();i++){
    Int_t ii = index.at(i);
    Pos = Form("%d%d%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii));
    Par.clear();
    GetParameters(fEName,fDataSet,fStartRun,fEndRun,fChannel.at(ii),&Par);
    Low = Par.at(0);
    Up = Par.at(1);
    rebin = (Int_t)Par.at(2);
    
    if(fGoodBad.at(ii)==1){
      Energy[i] = (TH1D*)fhist.Get(Form("%s%s%d",fEName.c_str(),Pos.c_str(),fChannel.at(ii)));

      FileName = Form("%s_%d_%d_%s_%d",fEName.c_str(),fStartRun,fEndRun,Pos.c_str(),fChannel.at(ii));
      ParName.clear();
      Par.clear();
      ParErr.clear();
      Cov.clear();
      TH1F *h1 = (TH1F*)Energy[i]->Clone();
      TH1F *h2 = (TH1F*)Energy[i]->Clone();
      h1->Rebin(rebin);
      h2->Rebin(rebin);
      h1->GetXaxis()->SetRangeUser(Low,Up);
      Int_t maxbin = h1->GetMaximumBin();
      Double_t scaleenergy = h1->GetBinCenter(maxbin)/2614.5333;
      h1->GetXaxis()->SetRangeUser(2500*scaleenergy,2700*scaleenergy);
      Double_t scaleamps = h1->Integral()*h1->GetBinWidth(maxbin)/35000;
      fout << Low << " " << Up << " " << maxbin << " " <<  scaleenergy << " "<< scaleamps << endl;
      h1->GetXaxis()->SetRange();
      
      ds.MultiPeakFitter(h2,scaleenergy, scaleamps, FileName,&ParName,&Par,&ParErr, &Cov);
      for(size_t ip = 0;ip<ParName.size();ip++){
	fout << ip << " " << ParName.at(ip).c_str() << " " << Par.at(ip) << " " << ParErr.at(ip) << endl;
      }
      
      fOffset = -Par.at(9)/Par.at(10);
      fScale = 1./Par.at(10);
      fOffsetErr = ParErr.at(9)/Par.at(10);
      fScaleErr = ParErr.at(10)/pow(Par.at(10),2);
      P0 = Par.at(11);
      P1 = Par.at(12);
      P2 = Par.at(13);
      P0Err = ParErr.at(11);
      P1Err = ParErr.at(12);
      P2Err = ParErr.at(13);
      
      fcal << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " 
	   << fOffset << " " << fOffsetErr << " " << fScale << " " << fScaleErr << endl;
      freso << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " 
	    << P0 << " " << P0Err << " " 
	    << P1 << " " << P1Err << " "	
	    << P2 << " " << P2Err <<endl;
      TMatrixD mat(2,2);
      for(Int_t i1 = 0;i1<2;i1++){
	for(Int_t i2 = 0;i2<2;i2++){
	  Index = (i1+9)*Par.size()+(i2+9);
	  mat(i1,i2) = Cov.at(Index);
	  cov.push_back(Cov.at(Index));
	}
      }
      fcov << Pos.c_str() << " " << fChannel.at(ii) << " " << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " 
	   << cov.at(0) << " " << cov.at(1) << " " << cov.at(2) << " " << cov.at(3) << endl;
      
      delete h1;
      delete h2;
    }
  }
}


