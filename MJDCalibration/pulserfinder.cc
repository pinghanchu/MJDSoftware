// 2016.11.15 Pinghan Chu
// 1.This code automatically calibrates "trapE" or "trapENFxxx";not work on the onboard energy "energy". 
// 2.Before applying this code, need to login the database in order to upload the calibration parameters in the end. Ask Robert Varner(varnerrl@ornl.gov) for the password
// 
#include "GATAutoCal.hh"
#include "TFile.h"
#include "TSpectrum.h"
#include "TMath.h"

TH1D* GetHisto(string fFileName, string fHistoName,string fCutName){
  TChain *mjdTree = new TChain("mjdTree");
  mjdTree->Add(Form("%s",fFileName.c_str()));
  //TFile *fhist = TFile::Open(Form("%s",FileName.c_str()),"read");
  //TH1D *h = (TH1D*)fhist->Get(Form("%s",HistoName.c_str()));
  TH1D *h = new TH1D("h","",800000,10, 8010);
  mjdTree->Draw(Form("%s>>h",fHistoName.c_str()),Form("%s",fCutName.c_str()));
  return h;
}

void FindPeak(TH1D *h, vector<Double_t>* fPx, vector<Double_t>* fPy){
  TSpectrum *s = new TSpectrum(1000,1);
  Int_t nfound = s->Search(h,10,"new",0.0001);
  Double_t *xpeaks = s->GetPositionX();
  Double_t *ypeaks = s->GetPositionY();
  vector<Double_t> Px;
  vector<Double_t> Py;
  for(Int_t i = 0;i<nfound; i++){
    if(ypeaks[i]>3){
      fPx->push_back(xpeaks[i]);
      fPy->push_back(ypeaks[i]);
    }
  }
}

void Sort(vector<Double_t> *fPx, vector<Double_t> *fPy){
  const Int_t npeak = fPx->size();
  Double_t Px[npeak];
  Double_t Py[npeak];
  Int_t n[npeak];
  for(Int_t i = 0;i<npeak;i++){
    Px[i] = fPx->at(i);
    Py[i] = fPy->at(i);
  }
  TMath::Sort(npeak,Px,n,0);
  fPx->clear();
  fPy->clear();
  
  for(Int_t i = 0;i<npeak;i++){
    Int_t ii = n[i];
    fPx->push_back(Px[ii]);
    fPy->push_back(Py[ii]);
  }
}
/*
void IsExist(Int_t fRun, Int_t fChannel,string fEName,vector<Double_t>* Par){

  ifstream fin(Form("./List/pulser/pulser_%d.txt",fRun));
  Int_t cha,pos;
  string R[6];
  Double_t r[6];
  Int_t run;
  string ename;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> run >> ename >> R[0] >> R[1] >> R[2] >>R[3] >>R[4] >> R[5];
      for(Int_t  i = 0;i<6;i++){
	if(R[i]=="nan" || R[i] =="-nan" || R[i] == "inf" || R[i] == "-inf"){
	  r[i] = 0;
	}else{
	  r[i] = atof(R[i].c_str());
	}
      }
      if(cha == fChannel && ename == fEName){
	for(Int_t  i = 0;i<6;i++){
	  Par->push_back(r[i]);
	}
      }
    }
  }
}


void GetPeak(TH1D *Input, Double_t Low, Double_t Up, vector<Double_t>* Peak){
  //Locate the 2614 keV peak in the raw spectrum.

  TH1D *h = (TH1D*)Input->Clone();
  h->GetXaxis()->SetRangeUser(Low,Up);
  Int_t maxbin = h->GetMaximumBin();
  Double_t maxX = h->GetBinCenter(maxbin);
  Int_t nPeak = Peak->size();
  for(size_t i = 0;i<Peak->size();i++){
    Peak->at(i) = Peak->at(i)*(maxX/Peak->at(nPeak-1));
  }
  delete h;
}
*/
int main(int argc, char** argv)
{
  if(argc != 5 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run] [energy name] [position = module*10+string] [FileName]" << endl;
    return 1;
  }
  Int_t fRun = atoi(argv[1]);
  string fEName = argv[2];
  const char* EnergyName = Form("%s",fEName.c_str());
  Int_t fPos = atoi(argv[3]);
  string fFileName = argv[4];

  ofstream fpulser(Form("./List/pulser/pulser_%d.txt",fRun),ios::app); // calibration constants
  fpulser.precision(12);

  GATAutoCal ds(fRun,fRun);
  ds.SetEnergyName(EnergyName);
  string fDataSet = ds.GetDataSet();  
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr  = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  vector<Double_t> fPeak = ds.GetCalibrationPeak();
  vector<Int_t> index;
  Int_t pos = 0;
  for(size_t i = 0;i<fChannel.size();i++){
    pos = fCryo.at(i)*10+fStr.at(i);
    if(pos == fPos){
      index.push_back(i);
    }
  }

  string Pos;
  vector<Double_t> fPx;
  vector<Double_t> fPy;

  for(size_t i = 0;i<index.size();i++){
    Int_t ii = index.at(i);
    Pos = Form("%d%d%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii));
    
    if(fGoodBad.at(ii)==1 ){
      string fCut = Form("channel == %d && EventDC1Bits>0",fChannel.at(ii));
      TH1D *h1 = GetHisto(fFileName,fEName,fCut);
      fPx.clear();
      fPy.clear();
      FindPeak(h1, &fPx, &fPy);
      if(fPx.size()>1){
	Sort(&fPx,&fPy);      
      }
      Int_t pars = fPx.size();

      if(pars>0){
	h1->GetXaxis()->SetRangeUser(fPx.at(pars-1)-10,fPx.at(pars-1)+10);
	Double_t entries = h1->GetEffectiveEntries();
	Double_t mean = h1->GetMean();
	Double_t rms = h1->GetRMS();
	fpulser << Pos.c_str() << " " << fChannel.at(ii) << " " << fRun << " " << fEName.c_str() << " "
		<< mean << " " << rms/sqrt(entries) << endl;
      }else{
	fpulser << Pos.c_str() << " " << fChannel.at(ii) << " " << fRun << " " << fEName.c_str() << " "
                << -99999 << " " << 0 << endl;
      }
      delete h1;
      /*
      HistoName = Form("%s%s%d",fEName.c_str(),Pos.c_str(),fChannel.at(ii));
      
      TH1D *h1 = GetHisto(FileName,HistoName);
      h1->Rebin(rebin);
      h1->GetXaxis()->SetRangeUser(Low,Up);
      Int_t maxbin = h1->GetMaximumBin();
      Double_t maxX = h1->GetBinCenter(maxbin);
      h1->GetXaxis()->SetRangeUser(maxX-5,maxX+5);
      Int_t entries = h1->GetEffectiveEntries();
      if(entries>30){
	PlotName = Form("%s_%d_%s_%d",fEName.c_str(),fRun,Pos.c_str(),fChannel.at(ii));
	TitleName = Form("C%dP%dD%d, channel %d",fCryo.at(ii),fStr.at(ii), fDetpos.at(ii), fChannel.at(ii));
	Par.clear();
	ParErr.clear();
	ds.SkewGaussFit(h1,maxX,2,&Par,&ParErr,"",PlotName,TitleName,"");
	fpulser << Pos.c_str() << " " << fChannel.at(ii) << " " << fRun << " " << fEName.c_str() << " "
		<< Par.at(0) << " " << ParErr.at(0) << " " << Par.at(1) << " " << ParErr.at(1) << " " << Par.at(4) << " " << ParErr.at(4) << endl;
      }else{
	cout << Pos.c_str() << " " << fChannel.at(ii) << " " << fRun << " No Pulser! "<< "Max peak : " << maxX << " Total entries : " << entries <<endl;
	fpulser << Pos.c_str() << " " << fChannel.at(ii) << " " << fRun << " " << fEName.c_str() << " 0 0 0 0 0 0" << endl;
      }
      delete h1;
    }else{
      fpulser << Pos.c_str() << " " << fChannel.at(ii) << " " << fRun << " " << fEName.c_str() << " 0 0 0 0 0 0" << endl;
      */
    }
  }
}


