// 2016.11.15 Pinghan Chu
// 1.This code automatically calibrates "trapE" or "trapENFxxx";not work on the onboard energy "energy". 
// 2.Before applying this code, need to login the database in order to upload the calibration parameters in the end. Ask Robert Varner(varnerrl@ornl.gov) for the password
// 
#include "GATAutoCal.hh"
#include "TFile.h"


Double_t IsNan(string Input){
  Double_t output = 0;
  if(Input == "nan" || Input == "-nan" || Input == "inf" || Input == "-inf"){
    output = 0;
  }else{
    output = atof(Input.c_str());
  }
  return output;
}


void GetCal(Int_t fStartRun, Int_t fEndRun, Int_t fChannel, string fEName, vector<Double_t>* Par){
  ifstream fin(Form("./List/cal/calibration_%d_%d.txt",fStartRun,fEndRun));
  //ifstream fin(Form("./List/cal/cal1_%d_%d.txt",fStartRun,fEndRun));
  //ifstream fin("./List/calibration.txt");

  Int_t cha,pos;
  string x1,x2,x3,x4;
  Double_t r1,r2,r3,r4;
  Int_t run1,run2;
  string ename;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> run1 >> run2 >> ename >> x1 >> x2 >> x3 >> x4;
      r1 = IsNan(x1);
      r2 = IsNan(x2);
      r3 = IsNan(x3);
      r4 = IsNan(x4);
      if(cha == fChannel && ename == fEName && run1 == fStartRun && run2 == fEndRun){
        Par->push_back( r1);  // offset
        Par->push_back( r2 );
        Par->push_back( r3 ); // slope
        Par->push_back( r4 );
      }
    }
  }
}



void GetEzfit(Int_t fStartRun, Int_t fEndRun, Int_t fChannel, string fEName, Int_t fIndex,vector<Double_t>* Par, vector<Double_t>* ParErr){

  //cout << fStartRun << " " << fEndRun << " "<< fChannel << " " << fEName.c_str() << endl;
  ifstream fin2(Form("./List/cal/ezfit_%d_%d.txt",fStartRun,fEndRun));
  string x0,x1,x2,x3,x4;
  Double_t r1,r2,r3,r4;
  Int_t run1,run2;
  Int_t cha,pos,index;
  string ename;
  if(fin2.is_open()){
    while(!fin2.eof()){
      fin2 >> pos>> cha >> run1>>run2>>ename >> index >> x0 >>x1 >> x2 >>x3 >>x4;     
      r1 = IsNan(x1);
      r2 = IsNan(x2);
      r3 = IsNan(x3);
      r4 = IsNan(x4);
      if(cha == fChannel && ename == fEName && index==fIndex){
        Par->push_back( (fStartRun+fEndRun)/2); // Run
        Par->push_back( r1 ); //peak
        Par->push_back( r3 ); //reso
        ParErr->push_back( (fEndRun-fStartRun)/2 );
        ParErr->push_back( r2 );
        ParErr->push_back( r4 );
      }
    }
  }
}



int main(int argc, char** argv)
{
  if(argc != 5 ) {
    cout << "Usage: " << argv[0] << " [start run] [endrun] [energy name] [channel]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  string fEName = argv[3];
  const char* EnergyName = Form("%s",fEName.c_str());
  Int_t fChan = atoi(argv[4]);

  ofstream ferror(Form("./List/cal/error_%d_%d.txt",fStartRun,fEndRun),ios::app); // calibration constants
  ferror.precision(12);
  ofstream ferrorchan(Form("./List/cal/errorchan_%d_%d.txt",fStartRun,fEndRun),ios::app); // calibration constants
  ferrorchan.precision(12);
  ofstream fcal("./List/calcorr.txt",ios::app);
  fcal.precision(12);
  GATAutoCal ds(fStartRun,fStartRun);
  ds.SetEnergyName(EnergyName);
  string fDataSet = ds.GetDataSet();  
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr  = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  vector<Double_t> fCalibrationPeak = ds.GetCalibrationPeak();
  vector<Int_t> index;
  string Pos;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  Int_t pars;
  Double_t offset = 999999;
  Double_t slope = 0;
  Double_t calpeak = 1;
  Double_t peak = 1;
  Int_t parerrs;
  //cout << fStartRun << " " << fEndRun << endl;
  for(size_t i = 0;i<fChannel.size();i++){
    Pos = Form("%d%d%d",fCryo.at(i),fStr.at(i),fDetpos.at(i));
    cout << Pos.c_str() << " " << fChannel.at(i) << endl;
    if(fChan > 0){
      if(fGoodBad.at(i)==1 && fChannel.at(i) == fChan){
	Par.clear();
	ParErr.clear();
	
	GetCal(fStartRun, fEndRun, fChannel.at(i), fEName, &Par);
	if(Par.size()>0){
	  pars = Par.size();
	  offset = Par.at(pars-4);
	  slope = Par.at(pars-2);
	  fcal << Pos.c_str() << " " << fChannel.at(i) << " " << fStartRun << " " << fEndRun << " "
	       << fEName.c_str() << " "
	       << Par.at(pars-4) << " " << Par.at(pars-3) << " " << Par.at(pars-2) << " " << Par.at(pars-1) << endl;
	  cout << Par.at(pars-4) << " " << Par.at(pars-3) << " " << Par.at(pars-2) << " " << Par.at(pars-1) << endl;
	  
	  //Int_t count = 0;
	  GetEzfit(fStartRun, fEndRun, fChannel.at(i), fEName, 3, &Par, &ParErr);
	  pars = Par.size();
	  parerrs = ParErr.size();
	  peak = Par.at(pars-2);
	  peak = peak*slope+offset;
	  calpeak = fCalibrationPeak.at(3);

	  if(abs(peak-calpeak)>2 && peak < 9000){
	    ferror << fStartRun << " " << fEndRun<< " "<< Pos.c_str() << " " << fChannel.at(i) << " "
		   << fEName.c_str() << " " << 3 << " " << peak << " " <<  abs(peak-calpeak) << " " << endl;
	  }
	  /*
	  if( ParErr.at(parerrs-1) == 0 || ParErr.at(parerrs-2) == 0 ){
	    ferror << fStartRun << " " << fEndRun<< " "<< Pos.c_str() << " " << fChannel.at(i) << " "
		   << fEName.c_str() << " " << 3 << " " << ParErr.at(parerrs-1) << " " <<  ParErr.at(parerrs-2) << " " << endl;
	  }
	  */
	  cout << peak << " " <<  abs(peak-calpeak) << " " << endl;
	}else{
          ferror << fStartRun << " " << fEndRun<< " "<< Pos.c_str() << " " << fChannel.at(i) << " "
                 << fEName.c_str() << " " << 3 << " " << 0 << " " <<  0 << " " << endl;
	}	  
      }
    }else{
      if(fGoodBad.at(i)==1){
        Par.clear();
        ParErr.clear();
        GetCal(fStartRun, fEndRun, fChannel.at(i), fEName, &Par);
        pars = Par.size();
        offset = Par.at(pars-4);
        slope = Par.at(pars-2);
        //Int_t count = 0;
        GetEzfit(fStartRun, fEndRun, fChannel.at(i), fEName, 3, &Par, &ParErr);
        pars = Par.size();
        peak = Par.at(pars-2);
        peak = peak*slope+offset;
        calpeak = fCalibrationPeak.at(3);
        if(abs(peak-calpeak)>2 && peak < 9000){
          ferror << fStartRun << " " << fEndRun<< " "<< Pos.c_str() << " " << fChannel.at(i) << " "
                 << fEName.c_str() << " " << 3 << " " << peak << " " <<  abs(peak-calpeak) << " " << endl;
        }
      }
    }
  }
}


