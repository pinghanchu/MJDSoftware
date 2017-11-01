// 2016.2.18
// written by Pinghan Chu
// Following the logic in GATAutoCal.cc
// 
#include "GATAutoCal.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <string>
#include <stdlib.h>

Double_t IsNan(string Input){
  Double_t output = 0;
  if(Input == "nan" || Input == "-nan" || Input == "inf" || Input == "-inf"){
    output = 0;
  }else{
    output = atof(Input.c_str());
  }
  return output;
}


void GetPar(Int_t fStartRun, Int_t fEndRun, string fEName, string fPos, string fChannel, string fRefPar, vector<Double_t>* Par, vector<Double_t>* ParErr){
  ifstream fin2(Form("./%d_%d/parameters_%d_%d_%s_%s_%s.txt",fStartRun,fEndRun,fStartRun,fEndRun,fEName.c_str(),fPos.c_str(),fChannel.c_str()));
  string x1,x2;
  Double_t r1,r2;
  Int_t ind;
  string par;
  if(fin2.is_open()){
    while(!fin2.eof()){
      fin2 >> ind >> par >> x1 >> x2;      
      r1 = IsNan(x1);
      r2 = IsNan(x2);      
      if(par == fRefPar){
	Par->push_back(r1);
	ParErr->push_back(r2);
      }
    }
  }
}

void GetCal(Int_t fStartRun, Int_t fEndRun, string fEName, string fPos, string fChannel, vector<Double_t>* Par, vector<Double_t>* ParErr){
  ifstream fin2(Form("./%d_%d/calibration_%d_%d_%s_%s.txt",fStartRun,fEndRun,fStartRun,fEndRun,fEName.c_str(),fPos.c_str()));
  string x1,x2,x3,x4,x5,x6;
  Double_t r1,r2,r3,r4,r5,r6;
  Int_t run1,run2;
  string pos,chan,ename;
  if(fin2.is_open()){
    while(!fin2.eof()){
      fin2 >> run1 >> run2 >> ename >> pos >> chan>> x1 >> x2 >> x3 >> x4 >> x5 >> x6;
      r1 = IsNan(x1);
      r2 = IsNan(x2);
      r3 = IsNan(x3);
      r4 = IsNan(x4);
      r5 = IsNan(x5);
      r6 = IsNan(x6);
      if(fChannel == chan){
        Par->push_back(r1);
        ParErr->push_back(r2);
        Par->push_back(r3);
        ParErr->push_back(r4);
      }
    }
  }
}

void GetCov(Int_t fStartRun, Int_t fEndRun, string fEName, string fPos, string fChannel, vector<Double_t>* Par){
  ifstream fin2(Form("./%d_%d/cov_%d_%d_%s_%s.txt",fStartRun,fEndRun,fStartRun,fEndRun,fEName.c_str(),fPos.c_str()));
  string x1,x2,x3,x4,x5,x6;
  Double_t r1,r2,r3,r4,r5,r6;
  Int_t run1,run2;
  string pos,chan,ename;
  if(fin2.is_open()){
    while(!fin2.eof()){
      fin2 >> run1 >> run2 >> ename >> pos >> chan>> x1 >> x2 >> x3 >> x4;
      r1 = IsNan(x1);
      r2 = IsNan(x2);
      r3 = IsNan(x3);
      r4 = IsNan(x4);
      //cout << run1 << " " << run2 << " " << r1 << " " << r2 << " " << r4 <<endl;
      if(fChannel == chan){
        Par->push_back(r1);
        Par->push_back(r2);
        Par->push_back(r3);
        Par->push_back(r4);
      }
    }
  }
}

TGraphErrors* GetGraph(vector<Double_t> Px,vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr){
  const Int_t nPy = Py.size();

  Double_t A[nPy];
  Double_t AErr[nPy];
  Double_t B[nPy];
  Double_t BErr[nPy];

  for(Int_t j=0;j<nPy;j++){
    A[j]=Px.at(j);
    AErr[j] = PxErr.at(j);
    B[j]=Py.at(j);
    BErr[j] =PyErr.at(j);
  }
  TGraphErrors* fGraph = new TGraphErrors(nPy,A,B,AErr,BErr);
  fGraph->SetFillStyle(0);
  fGraph->SetFillColor(0);


  return fGraph;
}


int main(int argc, char** argv)
{
  if(argc != 3) {
    cout << "Usage: " << argv[0] << "[energy name] [run list]" << endl;
    return 1;
  }
  //Int_t fDataSet = atoi(argv[1]);
  string fEName = argv[1];
  //string fRefPar =  argv[3];
  string fInputName = argv[2];

  gROOT->ProcessLine(".x $GATDIR/MJDCalibration/MJDTalkPlotStyle.C");
  //cout << fInputName.c_str() << endl;
  ifstream fin(Form("%s",fInputName.c_str()));
  vector<Int_t> startrun;
  vector<Int_t> endrun;
  Int_t r1,r2;
  string r3;
  vector<string> fDS;
  if(fin.is_open()){
    while(!fin.eof()){
      //fin >> r1 >> r2 >> r3 >> r4;
      fin >> r1 >> r2 >> r3;
      startrun.push_back(r1);
      endrun.push_back(r2);
      fDS.push_back(r3);
    }
  }
  startrun.pop_back();
  endrun.pop_back();
  fDS.pop_back();
  vector<string> SumName;
  SumName.push_back("00");
  SumName.push_back("10");
  SumName.push_back("HG");
  vector<Int_t> CalEnr;
  CalEnr.push_back(238);
  CalEnr.push_back(240);
  CalEnr.push_back(277);
  CalEnr.push_back(300);
  CalEnr.push_back(583);
  CalEnr.push_back(727);
  CalEnr.push_back(785);
  CalEnr.push_back(860);
  CalEnr.push_back(2614);

  vector<Double_t> CalEnr1;
  CalEnr1.push_back(238.632);
  CalEnr1.push_back(240.986);
  CalEnr1.push_back(277.371);
  CalEnr1.push_back(300.087);
  CalEnr1.push_back(583.191);
  CalEnr1.push_back(727.330);
  CalEnr1.push_back(785.37);
  CalEnr1.push_back(860.557);
  CalEnr1.push_back(2614.533);

  vector<string> CalName;
  for(size_t i = 0;i<CalEnr.size();i++){
    string name = Form("#mu_{%d}",CalEnr.at(i));
    CalName.push_back(name);
  }
  vector<string> TauName;
  TauName.push_back("#tau_{0}");
  TauName.push_back("#tau_{1}");
  vector<string> HSName;
  HSName.push_back("h_{s}_{0}");
  HSName.push_back("h_{s}_{1}");
  HSName.push_back("h_{s}_{2}");
  vector<string> Ft0Name;
  Ft0Name.push_back("f_{t}_{0}");
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> calPar;
  vector<Double_t> calParErr;
  vector<Double_t> Delta;
  vector<Double_t> tauPar;
  vector<Double_t> tauParErr;
  vector<Double_t> hsPar;
  vector<Double_t> hsParErr;
  vector<Double_t> ft0Par;
  vector<Double_t> ft0ParErr;
  vector<Double_t> Tau;
  vector<Double_t> HS;
  vector<Double_t> Ft0;
  for(size_t i=0;i<startrun.size();i++){
    Delta.clear();
    Tau.clear();
    HS.clear();
    Ft0.clear();
    for(size_t j = 0;j<SumName.size();j++){
      //cout<< SumName.at(j) << endl;
      Par.clear();
      ParErr.clear();
      GetCal(startrun.at(i),endrun.at(i),fEName,"0",SumName.at(j),&Par,&ParErr);
      Double_t offset = Par.at(0);
      Double_t offsetunc= ParErr.at(0);
      Double_t slope = Par.at(1);
      Double_t slopeunc = ParErr.at(1);
      Par.clear();
      ParErr.clear();
      GetCov(startrun.at(i),endrun.at(i),fEName,"0",SumName.at(j),&Par);
      Double_t cov11 = Par.at(0);
      Double_t cov22 = Par.at(3);
      Double_t cov12 = Par.at(1);
      Double_t delta_mu = TMath::Sqrt(cov11+2039*cov12+2039*2039*cov22);
      //Double_t delta_mu = TMath::Sqrt
      calPar.clear();
      calParErr.clear();
      tauPar.clear();
      tauParErr.clear();
      hsPar.clear();
      hsParErr.clear();
      ft0Par.clear();
      ft0ParErr.clear();
      for(size_t k = 0;k<CalName.size();k++){
	Par.clear();
	ParErr.clear();
	GetPar(startrun.at(i),endrun.at(i),fEName,"000",SumName.at(j),CalName.at(k),&Par,&ParErr);
	Double_t Delta_k = TMath::Sqrt( pow((Par.at(0)-CalEnr1.at(k)),2)+pow(ParErr.at(0),2));
	Double_t delta_k = TMath::Sqrt(cov11+CalEnr1.at(k)*cov12+CalEnr1.at(k)*CalEnr1.at(k)*cov22);
	//Double_t delta_k = TMath::Sqrt(offsetunc*offsetunc+CalEnr1.at(k)*CalEnr1.at(k)*slopeunc*slopeunc);
	//cout << Delta_k << " " << delta_k << endl;
	Double_t ratio = Delta_k/delta_k;	
	calPar.push_back(ratio);
	//calParErr.push_back(ParErr.at(0));
      }
      
      for(size_t k=0;k<TauName.size();k++){
	Par.clear();
	ParErr.clear();
	GetPar(startrun.at(i),endrun.at(i),fEName,"000",SumName.at(j),TauName.at(k),&Par,&ParErr);
	tauPar.push_back(Par.at(0));
	tauParErr.push_back(ParErr.at(0));
      } 

      for(size_t k=0;k<HSName.size();k++){
        Par.clear();
        ParErr.clear();
        GetPar(startrun.at(i),endrun.at(i),fEName,"000",SumName.at(j),HSName.at(k),&Par,&ParErr);
        hsPar.push_back(Par.at(0));
        hsParErr.push_back(ParErr.at(0));
      }

      for(size_t k=0;k<Ft0Name.size();k++){
        Par.clear();
        ParErr.clear();
        GetPar(startrun.at(i),endrun.at(i),fEName,"000",SumName.at(j),Ft0Name.at(k),&Par,&ParErr);
        ft0Par.push_back(Par.at(0));
        ft0ParErr.push_back(ParErr.at(0));
      }

      Double_t tau0 = tauPar.at(0);
      Double_t tau1 = tauPar.at(1);
      Double_t hs0 = hsPar.at(0);
      Double_t hs1 = hsPar.at(1);
      Double_t hs2 = hsPar.at(2);
      Double_t ft0 = ft0Par.at(0);
      Double_t tauroi = tau0 + tau1*2039;
      Double_t hsroi = hs0/(2039*2039)+hs1*pow(2039,hs2);
      //cout << SumName.at(j) << " " << tau0 << " " << tau1 << " " << hs0 << " " << hs1 << " " << hs2 << endl;
      Tau.push_back(tauroi);
      HS.push_back(hsroi);
      Ft0.push_back(ft0);
      Double_t sum = 0;
      for(size_t k=0;k<calPar.size();k++){
	sum = sum + pow(calPar.at(k),2);
      }
      Double_t CorRatio = TMath::Sqrt(sum/((int)calPar.size()-2));
      Double_t delta_mu_NL = delta_mu*CorRatio;
      Delta.push_back(delta_mu_NL);
      //cout << "DS"<<fDS.at(i).c_str() << " " << SumName.at(j) << " " << delta_mu << " " << delta_mu_NL << endl;
    }
    //cout << "DS"<< fDS.at(i).c_str() << " & " << Tau.at(0) << " & " << Tau.at(1) << " & " << Tau.at(2) << "\\" << "\\" << endl;
    //cout << "DS"<< fDS.at(i).c_str() << " & " << HS.at(0) << " & " << HS.at(1) << " & " << HS.at(2) << "\\"<<"\\"<<endl;
    //cout << "DS"<< fDS.at(i).c_str() << " & " << Delta.at(0) << " & " << Delta.at(1) << " & " << Delta.at(2) << "\\"<<"\\"<<endl;
    cout << "DS"<< fDS.at(i).c_str() << " & " << Ft0.at(0) << " & " << Ft0.at(1) << " & " << Ft0.at(2) << "\\"<<"\\"<<endl;

  }
}


