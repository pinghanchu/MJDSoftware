// 2016.2.18
// This macro uploads the calibration parameters.
// 
#include "GATAutoCal.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

void linearfit(Int_t StartRun,Int_t EndRun,Int_t Channel,string EName, vector<Double_t> *offset,vector<Double_t> *offseterr,vector<Double_t> *slope,vector<Double_t> *slopeerr){
  ifstream fin(Form("./List/%d_%d/calibration_%d_%d.txt",StartRun, EndRun, StartRun,EndRun));
  Int_t cha,pos;
  Double_t r1,r2,r3,r4,r5,r6;
  Int_t run1,run2;
  string ename;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run1 >> run2 >> ename >> pos >> cha >> r1 >> r2 >>r3 >>r4 >> r5 >> r6;
      if(cha == Channel && ename == EName){
        offset->push_back(r1);
        offseterr->push_back(r2);
	slope->push_back(r3);
	slopeerr->push_back(r4);
      }
    }
  }
}



void covariant(Int_t StartRun,Int_t EndRun,Int_t Channel,string EName, vector<Double_t> *cov){

  ifstream fin(Form("./List/%d_%d/cov_%d_%d.txt",StartRun,EndRun,StartRun,EndRun));

  Int_t cha,pos;
  Double_t r1,r2,r3,r4;
  Int_t run1,run2;
  string ename;
  
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run1 >> run2 >> ename >> pos >> cha >> r1 >> r2 >>r3 >>r4;
      if(cha == Channel && ename == EName){
        cov->push_back(r1);
        cov->push_back(r2);
        cov->push_back(r3);
        cov->push_back(r4);
      }
    }
  }
}


int main(int argc, char** argv)
{
  if(argc != 6 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [cover startrun] [cover endrun] [energy name]" << endl;
    return 1;
  }
  int fStartRun = atoi(argv[1]);
  int fEndRun = atoi(argv[2]);
  int fCoverStartRun = atoi(argv[3]);
  int fCoverEndRun = atoi(argv[4]);
  string fEName = argv[5];

  GATAutoCal ds(fStartRun,fStartRun);
  string fDataSet = ds.GetDataSet();
  Int_t gatrev=ds.GetGATRev();

  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fPulserChannel = ds.GetPulserChannel();
  vector<Int_t> fCryo=ds.GetCryo();
  vector<Int_t> fStr= ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  vector<string> fDetName = ds.GetDetectorName();
  
  const Int_t channels = fChannel.size();
  Int_t startdate = 0;
  Int_t enddate = 0;
  Int_t coverstartdate = 0;
  Int_t coverenddate = 0;
  startdate = (Int_t)ds.GetStartTime();
  enddate = (Int_t) ds.GetStopTime();

  GATAutoCal ds1(fCoverStartRun,fCoverStartRun);
  coverstartdate = (Int_t) ds1.GetStartTime();
   
  if( !(fCoverEndRun == 64999999 || fCoverEndRun == 4999999)){
    GATAutoCal ds2(fCoverEndRun-1,fCoverEndRun);
    coverenddate = (Int_t) ds2.GetStopTime();
  }
  cout << "Gat Rev:" <<gatrev << "; [start run] = " << fStartRun << "; [end run] = " << fEndRun << "; [start time] = " << startdate << "; [end time] = " << enddate <<
    "; [cover start run] = " << fCoverStartRun << "; [cover end run] = " << fCoverEndRun << "; [cover start time] = " << coverstartdate << "; [cover end time] = " << coverenddate << endl;
  ds.SetUpProvenance("GATAutoCal Energy Calibration", "P.Chu",gatrev, fStartRun, fEndRun, fCoverStartRun, fCoverEndRun,startdate, enddate,coverstartdate,coverenddate);

  vector<Double_t> offset;
  vector<Double_t> offseterr;
  vector<Double_t> slope;
  vector<Double_t> slopeerr;
  vector<Double_t> cov;
  vector<Double_t> covcorr;

  Double_t slo, sloer, off, offer;
  Int_t fPos = 0; 
  for(Int_t i = 0;i<channels;i++){
    fPos = fCryo.at(i)*100+fStr.at(i)*10+fDetpos.at(i);
    offset.clear();
    offseterr.clear();
    slope.clear();
    slopeerr.clear();
    cov.clear();
    covcorr.clear();
    cov.push_back(0);
    cov.push_back(0);
    cov.push_back(0);
    cov.push_back(0);

    linearfit(fStartRun,fEndRun,fChannel.at(i),fEName,&offset,&offseterr,&slope,&slopeerr);
    covariant(fStartRun,fEndRun,fChannel.at(i),fEName,&cov);

    Int_t pars = slope.size();
    
    if( fGoodBad.at(i) == 1 ){
      if( pars==0){
	cout << "Error : Good channel " << fChannel.at(i) << " of detector " << fPos << " has no calibration parameters!" << endl;
	slo = 0;
	sloer = 0;
	off = 999999;
	offer = 0;
	cout << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " << fPos << " " << fChannel.at(i) << " " << slo << " " << sloer << " " << off << " " << offer << endl;
	ds.PutECalMJDB(fChannel.at(i), fDetName.at(i), slo,sloer,off,offer, covcorr);
      }else{
	slo = slope.at(pars-1);
	sloer = slopeerr.at(pars-1);
	off = offset.at(pars-1);
	offer = offseterr.at(pars-1);
	for(size_t ic=cov.size()-4;ic<cov.size();ic++){
	  covcorr.push_back(cov.at(ic));
	}
	cout << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " << fPos << " " << fChannel.at(i) << " " << slo << " " << sloer << " " << off << " " << offer << endl;
	ds.PutECalMJDB(fChannel.at(i), fDetName.at(i), slo,sloer,off,offer, covcorr);
      }
    }else{
      cout << "Bad channel " << fChannel.at(i) << " of detector " << fPos << " has no calibration parameters!" << endl;
      slo = 0;
      sloer = 0;
      off = 999999;
      offer = 0;
      cout << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " << fPos << " " << fChannel.at(i) << " " << slo << " " << sloer << " " << off << " " << offer << endl;
      ds.PutECalMJDB(fChannel.at(i), fDetName.at(i), slo,sloer,off,offer, cov);
    }
  }

  cov.clear();
  cov.push_back(0);
  cov.push_back(0);
  cov.push_back(0);
  cov.push_back(0);

  
  for(size_t i=0;i<fPulserChannel.size();i++){
    slo = 0;
    sloer = 0;
    off = 999999;
    offer = 0;
    cout << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " << "000" << " " << fPulserChannel.at(i) << " " << slo << " " << sloer << " " << off << " " << offer << endl;
    ds.PutECalMJDB(fPulserChannel.at(i), "Pulser", slo,sloer,off,offer, cov);
  }


  vector<Double_t> UnknownChannel;
  

  //Unknown channel for DS1.
  //UnknownChannel.push_back(674);
  //UnknownChannel.push_back(675);
  //UnknownChannel.push_back(676);
  // UnknownChannel.push_back(677);
  //Pulser channel for DS5;
  UnknownChannel.push_back(644);
  //Pulser channel for DS4;
  //UnknownChannel.push_back(1334);

  for(size_t i = 0;i<UnknownChannel.size();i++){
    slo = 0;
    sloer = 0;
    off = 999999;
    offer = 0;
    cout << fStartRun << " " << fEndRun << " " << fEName.c_str() << " " << "000" << " " << fPulserChannel.at(i) << " " << slo << " " << sloer << " " << off << " " << offer << endl;
    ds.PutECalMJDB(UnknownChannel.at(i), "Unknown", slo,sloer,off,offer, cov);
  }


}


