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


int main(int argc, char** argv)
{
  if(argc != 9 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [start run] [endrun] [cover startrun] [cover endrun] [energy name] [position Mod*100+String*10+Detpos] [Input File] [Output File Path]" << endl;
    return 1;
  }
  
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  Int_t fCoverStartRun = atoi(argv[3]);
  Int_t fCoverEndRun = atoi(argv[4]);
  string fEName = argv[5];
  Int_t fPos = atoi(argv[6]);
  const char* fEnergyName = Form("%s",fEName.c_str());
  string fInputFile = argv[7];
  string fOutputPath = argv[8];
  
  ofstream fcal(Form("%scalibration_%d_%d_%s_%d.txt", fOutputPath.c_str(),fStartRun,fEndRun,fEName.c_str(),fPos),ios::app);
  ofstream fcov(Form("%scov_%d_%d_%s_%d.txt", fOutputPath.c_str(),fStartRun,fEndRun,fEName.c_str(),fPos),ios::app);
  ofstream freso(Form("%sreso_%d_%d_%s_%d.txt", fOutputPath.c_str(),fStartRun,fEndRun,fEName.c_str(),fPos),ios::app);
 
  fcal.precision(15);
  fcov.precision(15);
  freso.precision(15);

  GATAutoCal ds(fStartRun,fStartRun);
  string fDataSet = ds.GetDataSet();
  Int_t gatrev=ds.GetGATRev();
  ds.SetEnergyName(fEnergyName);
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr  = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  vector<Int_t> index;

  /*
  Int_t startdate = 0;
  Int_t enddate = 0;
  Int_t coverstartdate = 0;
  Int_t coverenddate = 0;
  startdate = (Int_t)ds.GetStartTime();
  enddate = (Int_t) ds.GetStopTime();

  GATAutoCal ds1(coverstartrun,coverstartrun+1);
  coverstartdate = (Int_t) ds1.GetStartTime();

  if( !(coverendrun == 64999999 || coverendrun == 4999999)){
    GATAutoCal ds2(coverendrun-1,coverendrun);
    coverenddate = (Int_t) ds2.GetStopTime();
  }
  cout << "Gat Rev:" <<gatrev << "; [start run] = " << startrun << "; [end run] = " << endrun << "; [start time] = " << startdate << "; [end time] = " << enddate << 
    "; [cover start run] = " << coverstartrun << "; [cover end run] = " << coverendrun << "; [cover start time] = " << coverstartdate << "; [cover end time] = " << coverenddate << endl;
  ds.SetUpProvenance("GATAutoCal Energy Calibration", "P.Chu",gatrev, fStartRun, fEndRun, fCoverStartRun, fCoverEndRun,startdate, enddate,coverstartdate,coverenddate);
*/

  for(size_t i = 0;i<fChannel.size();i++){
    Int_t pos = fCryo.at(i)*100+fStr.at(i)*10+fDetpos.at(i);
    if(pos == fPos){
      index.push_back(i);
    }
  }

  Double_t Low, Up;
  Int_t rebin = 1;
  Double_t EnergyROI = 2039;
  vector<string> ParName;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> CalPeak;
  vector<Double_t> Cov;
  vector<Double_t> CovLinear;
  vector<Double_t> Px;
  vector<Double_t> PxErr;
  vector<Double_t> Py;
  vector<Double_t> PyErr;

  string Pos;
  string FitName;
  string TitleName;
  for(size_t i = 0;i<index.size();i++){
    Int_t ii = index.at(i);
    Pos = Form("%d%d%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii));
    Par.clear();
    GetParameters(fEName,fDataSet,fStartRun,fEndRun,fChannel.at(ii),&Par);
    Low = Par.at(0);
    Up = Par.at(1);
    rebin = (Int_t)Par.at(2);
    FitName = Form("%d_%d_%s_%s_%d",fStartRun,fEndRun,fEName.c_str(), Pos.c_str(),fChannel.at(ii));
    TitleName = Form("C%dP%dD%d, Channel = %d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii),fChannel.at(ii));
    string fHistoName = Form("%s%d",fEName.c_str(),fChannel.at(ii));    

    if(fGoodBad.at(ii)==1){
      TH1D* h = ds.LoadHisto(fInputFile, fHistoName);
      if(h->GetEntries()>0){
	ParName.clear();
	Par.clear();
	ParErr.clear();
	Cov.clear();
	CovLinear.clear();
	TH1F *h1 = (TH1F*)h->Clone();
	TH1F *h2 = (TH1F*)h->Clone();
	h1->Rebin(rebin);
	h2->Rebin(rebin);
	h1->GetXaxis()->SetRangeUser(Low,Up);
	Int_t maxbin = h1->GetMaximumBin();
	Double_t scaleenergy = h1->GetBinCenter(maxbin)/2614.5333;
	h1->GetXaxis()->SetRangeUser(2500*scaleenergy,2700*scaleenergy);
	Double_t scaleamps = h1->Integral()*h1->GetBinWidth(maxbin)/35000;
	cout << "Lower bound:" << Low << "; Upper bound:" << Up << "; Maximum Bin:" << maxbin << "; Scale factor:" <<  scaleenergy << "; Scale amplitude:"<< scaleamps << endl;
	
	h1->GetXaxis()->SetRange();
	Int_t index = 18;
	Int_t isValid = ds.MultiPeakFit(h2,scaleenergy, scaleamps, FitName,&ParName,&Par,&ParErr, &Cov, &CalPeak);
	freso << fStartRun << " " << fEndRun << " " << fEName.c_str() << " "<< fPos << " " << fChannel.at(ii) << " "<< Par.at(index) << " " << ParErr.at(index) << " " 
	      << Par.at(index+1) << " " << ParErr.at(index+1) <<  " " << Par.at(index+2) << " " << ParErr.at(index+2) <<endl;
	
	for(size_t ip = 0;ip<ParName.size();ip++){
	  cout << ip << " " << ParName.at(ip).c_str() << " " << Par.at(ip) << " " << ParErr.at(ip) << endl;
	}
	for(size_t ip = 0;ip<CalPeak.size();ip++){
	  cout << ip << " " << CalPeak.at(ip) << endl;
	}
	cout << "The fitting is valid : " << isValid << endl;
	
	Px.clear();
	PxErr.clear();
	Py.clear();
	PyErr.clear();
	for(Int_t i = 0;i<index/2;i++){
	  Int_t ii = i+index/2;
	  Px.push_back(Par.at(ii));
	  PxErr.push_back(ParErr.at(ii));
	  Py.push_back(CalPeak.at(i));
	  PyErr.push_back(0);
	}
	
	Par.clear();
	ParErr.clear();
	Cov.clear();
	ds.LinearFit(Px,PxErr,Py,PyErr,FitName,TitleName,EnergyROI, &Par,&ParErr,&Cov); 
	fcal << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " << fPos << " " << fChannel.at(ii) << " "<< Par.at(0) << " " << ParErr.at(0) << " " 
	     << Par.at(1) << " " << ParErr.at(1) << " " << Par.at(2) << " " << ParErr.at(2) <<  endl;
	fcov << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " << fPos << " " << fChannel.at(ii) << " "<< Cov.at(0) << " " << Cov.at(1) << " " << Cov.at(2) << " " << Cov.at(3) << endl;
	
	delete h1;
	delete h2;
      }else{
	cout << "This is an empty channel!" <<endl;
      }
    }
  }
}


