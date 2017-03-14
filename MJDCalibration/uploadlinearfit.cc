// 2016.2.18
// written by Pinghan Chu
// Following the logic in GATAutoCal.cc
// 
#include "GATAutoCal.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
/*
void multifit(Int_t StartRun,Int_t EndRun,Int_t Channel,string EName, vector<Double_t> *peak,vector<Double_t> *peakerr,vector<Double_t> *width,vector<Double_t> *widtherr){

  ifstream fin(Form("./List/multifit_%d_%d.txt",StartRun,EndRun));

  Int_t pos;
  Int_t cha;
  Int_t index;
  Int_t goodfit;
  Double_t r1,r2,r3,r4;
  string r5,r6,r7,r8;
  string ename;
  Double_t calpeak;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> ename >> index >> calpeak>> r1 >> r2 >>r3 >>r4 >> r5>>r6>>r7>>r8 >> goodfit;
      if(ename == EName && cha == Channel){
        peak->push_back(r1);
        peakerr->push_back(r2);
        width->push_back(r3*2.355);
        widtherr->push_back(r4*2.355);
      }
    }
  }

  while(peak->size()>4){
    peak->pop_back();
    peakerr->pop_back();
    width->pop_back();
    widtherr->pop_back();
  }
  if(peak->size()==0){
    peak->push_back(-1);
    peakerr->push_back(-1);
    width->push_back(0);
    widtherr->push_back(0);
  }

}

*/

void linearfit(Int_t StartRun,Int_t EndRun,Int_t Channel,string EName, vector<Double_t> *offset,vector<Double_t> *offseterr,vector<Double_t> *slope,vector<Double_t> *slopeerr){
  ifstream fin(Form("./List/cal/cal1_%d_%d.txt",StartRun,EndRun));
  Int_t cha,pos;
  Double_t r1,r2,r3,r4;
  Int_t run1,run2;
  string ename;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> run1 >> run2 >> ename >> r1 >> r2 >>r3 >>r4;
      //cout << cha << endl;
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

  ifstream fin(Form("./List/cal/cov1_%d_%d.txt",StartRun,EndRun));

  //ifstream fin(Form("./List/cov1_%d_%d.txt",StartRun,EndRun));
  Int_t cha;
  Double_t r1,r2,r3,r4;
  Int_t run1,run2;
  string ename;
  
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> cha >> run1 >> run2 >> ename >> r1 >> r2 >>r3 >>r4;
      if(cha == Channel && ename == EName){
        cov->push_back(r1);
        cov->push_back(r2);
        cov->push_back(r3);
        cov->push_back(r4);
      }
    }
  }

  while(cov->size()>4){
    cov->pop_back();
  }
}

/*
void GetDecayConstant(Int_t DataSet, vector<Int_t> Channel, vector<Int_t>* DecayConstant){
  ifstream fin(Form("./List/optCT_DS%d.txt",DataSet));
  Int_t pos,cha,dc;
  vector<Int_t> Cha;
  vector<Int_t> DC;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> dc;
      
      Cha.push_back(cha);
      DC.push_back(dc);
    }
  }
  Cha.pop_back();
  DC.pop_back();

  for(size_t i=0;i<Channel.size();i++){
    Int_t count=0;
    Int_t index=0;
    for(size_t j = 0;j<Cha.size();j++){
      if(Channel.at(i) == Cha.at(j)){
	count++;
	index=j;
      }
    }
    if(count>0){
      DecayConstant->push_back(DC.at(index));
    }else{
en     DecayConstant->push_back(70);
    }
  }
}
*/

int main(int argc, char** argv)
{
  if(argc != 6 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [cover startrun] [cover endrun] [energy name]" << endl;
    return 1;
  }
  int startrun = atoi(argv[1]);
  int endrun = atoi(argv[2]);
  int coverstartrun = atoi(argv[3]);
  int coverendrun = atoi(argv[4]);
  string EName = argv[5];
  string EName1;
  if(EName == "trapENFDB"){
    EName1 = Form("%s","trapENF");
  }else{
    EName1 = Form("%s",EName.c_str());
  }

  GATAutoCal ds(startrun,startrun+1);
  string fDataSet = ds.GetDataSet();
  string dataset = ds.GetDataSet();
  Int_t gatrev=ds.GetGATRev();
  //Int_t gatrev = 0;
  vector<Int_t> channel = ds.GetChannel();
  vector<Int_t> pulserchannel = ds.GetPulserChannel();
  vector<Int_t> cryo=ds.GetCryo();
  vector<Int_t> str= ds.GetString();
  vector<Int_t> detpos = ds.GetDetPosition();
  vector<Int_t> goodbad = ds.GetGoodBad();
  vector<string> detname = ds.GetDetectorName();
  
  const Int_t channels = channel.size();
  Int_t startdate = 0;
  Int_t enddate = 0;
  Int_t coverstartdate = 0;
  Int_t coverenddate = 0;
  startdate = (Int_t)ds.GetStartTime();
  enddate = (Int_t) ds.GetStopTime();

  //  GATAutoCal ds1(coverstartrun,coverendrun);
  GATAutoCal ds1(coverstartrun,coverstartrun+1);
  coverstartdate = (Int_t) ds1.GetStartTime();
   
  if( !(coverendrun == 64999999 || coverendrun == 4999999)){
    GATAutoCal ds2(coverendrun-1,coverendrun);
    coverenddate = (Int_t) ds2.GetStopTime();
  }
  cout << gatrev << " " << startrun << " " << endrun << " " << startdate << " " << enddate << " " << coverstartrun << " " << coverendrun << " " << coverstartdate << " " << coverenddate << endl;
  ds.SetUpProvenance("GATAutoCal Energy Calibration", "P.Chu",gatrev, startrun, endrun, coverstartrun, coverendrun,startdate, enddate,coverstartdate,coverenddate);

  vector<Double_t> offset;
  vector<Double_t> offseterr;
  vector<Double_t> slope;
  vector<Double_t> slopeerr;
  vector<Double_t> cov;
  vector<Double_t> cov1;
  vector<Double_t> peak;
  vector<Double_t> peakerr;
  vector<Double_t> width;
  vector<Double_t> widtherr;

  Int_t is,id,pos;
  const char* EnergyName;
  if(EName == "trapENFDB"){
    EnergyName = Form("%s","trapENF");
  }else{
    EnergyName = Form("%s",EName.c_str());
  }
  ds.SetEnergyName(EnergyName);

  for(Int_t i = 0;i<channels;i++){
    offset.clear();
    offseterr.clear();
    slope.clear();
    slopeerr.clear();

    cov.clear();
    peak.clear();
    peakerr.clear();
    width.clear();
    widtherr.clear();

    cov.push_back(0);
    cov.push_back(0);
    cov.push_back(0);
    cov.push_back(0);

    is = str.at(i);
    id = detpos.at(i);
    pos = cryo.at(i)*100+is*10+id;
 
    linearfit(startrun,endrun,channel.at(i),EName,&offset,&offseterr,&slope,&slopeerr);
    Int_t pars = slope.size();
    EnergyName = Form("%s",EName.c_str());

    if(channel.at(i) >0 && (goodbad.at(i)==1 || slope.size()>0)){      
      //linearfit(startrun,endrun,channel.at(i),EName,&offset,&offseterr,&slope,&slopeerr);
      //covariant(startrun,endrun,channel.at(i),EName, &cov);
      Double_t slo = 0;
      Double_t sloer = 0;
      Double_t off = 999999;
      Double_t offer = 0;
      cov1.clear();
      if(slope.size()>0){
	slo = slope.at(pars-1);
	sloer = slopeerr.at(pars-1);
	off = offset.at(pars-1);
	offer = offseterr.at(pars-1);
	for(size_t ic=cov.size()-4;ic<cov.size();ic++){
	  cov1.push_back(cov.at(ic));
	}
      }else{
	slo = 0;
        sloer = 0;
        off = 999999;
        offer = 0;
	cov1.push_back(0);
	cov1.push_back(0);
	cov1.push_back(0);
	cov1.push_back(0);
      }
      //cout << pos << " "<<  channel.at(i) << " " << startrun << " " << endrun  << " " << EnergyName << " " << offset.at(0) << " " << offseterr.at(0) << " " << slope.at(0) << " " << slopeerr.at(0) << endl;
      //cout << pos << " " << channel.at(i) << " " << startrun << " " << endrun << " " << cov.at(0) << " " << cov.at(1) << " " << cov.at(2) << " " << cov.at(3) << endl;
      cout << pos << " "<< channel.at(i)<< " " << detname.at(i).c_str() << " " << startrun << " " << endrun << " " << EnergyName << " " << slo << " " << sloer <<" " << off << " "<< offer << endl;
      //cout << pos << " " << channel.at(i) << " " << startrun << " " << endrun << " " << cov1.at(0) << " " << cov1.at(1) << " " << cov1.at(2) << " " << cov1.at(3) << endl;
      //ds.PutECalMJDB(channel.at(i), detname.at(i), slo,sloer,off,offer, cov1);
      //cout <<  pos << " "<< channel.at(i)<< " " << detname.at(i).c_str() << " " << startrun << " " << endrun << " " << EnergyName << " " << 0 << " "<< 0<< " "<<  999999 << " " << 0 << endl;
      //ds.PutECalMJDB(channel.at(i), detname.at(i), 0,0,999999,0, cov);

    }else{

      cout <<  pos << " "<< channel.at(i)<< " " << detname.at(i).c_str() << " " << startrun << " " << endrun << " " << EnergyName << " " << 0 << " "<< 0<< " "<<  999999 << " " << 0 << endl;
      //ds.PutECalMJDB(channel.at(i), detname.at(i), 0,0,999999,0, cov);
    }
  }

  cov.clear();
  cov.push_back(0);
  cov.push_back(0);
  cov.push_back(0);
  cov.push_back(0);

  
  for(size_t i=0;i<pulserchannel.size();i++){
    cout << "000 " <<  pulserchannel.at(i)<< " " << "Pulser" << " " << startrun << " " << endrun << " " << EnergyName << " " << 0 << " "<< 0<< " "<<  999999 << " " << 0 << endl;
    
    //ds.PutECalMJDB(pulserchannel.at(i), "Pulser", 0,0,999999,0, cov);
  }


  vector<Double_t> UnknownChannel;
  

  //Unknown channel for DS1.
  UnknownChannel.push_back(674);
  UnknownChannel.push_back(675);
  UnknownChannel.push_back(676);
  UnknownChannel.push_back(677);
  //Pulser channel for DS5;
  //UnknownChannel.push_back(644);

  for(size_t i = 0;i<UnknownChannel.size();i++){
    cout << "000 " <<  UnknownChannel.at(i) << " " << "Unknown" << " " << startrun << " " << endrun << " " << EnergyName << " " 
	 <<  0 << " "<< 0<< " "<<  999999 << " " << 0 << endl;
    //ds.PutECalMJDB(UnknownChannel.at(i), "Unknown", 0,0,999999,0, cov);
  }


}


