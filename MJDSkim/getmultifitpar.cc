// 2016.2.18
// written by Pinghan Chu
// Following the logic in GATAutoCal.cc
// 
#include "MJDSkim.hh"
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


void GetPar(Int_t fStartRun, Int_t fEndRun, string fEName, Int_t fPos, Int_t fChannel, string fRefPar, vector<Double_t>* Par, vector<Double_t>* ParErr){
  ifstream fin2(Form("./%d_%d/parameters_%d_%d_%s_%d_%d.txt",fStartRun,fEndRun,fStartRun,fEndRun,fEName.c_str(),fPos,fChannel));
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
  if(argc != 5) {
    cout << "Usage: " << argv[0] << " [dataset] [energy name] [parameter] [run list]" << endl;
    return 1;
  }
  Int_t fDataSet = atoi(argv[1]);
  string fEName = argv[2];
  string fRefPar =  argv[3];
  string fInputName = argv[4];

  gROOT->ProcessLine(".x $GATDIR/MJDCalibration/MJDTalkPlotStyle.C");
  //cout << fInputName.c_str() << endl;
  ifstream fin(Form("%s",fInputName.c_str()));
  vector<Int_t> startrun;
  vector<Int_t> endrun;
  Int_t r1,r2,r3,r4;
  if(fin.is_open()){
    while(!fin.eof()){
      //fin >> r1 >> r2 >> r3 >> r4;
      fin >> r1 >> r2;
      startrun.push_back(r1);
      endrun.push_back(r2);
    }
  }
  startrun.pop_back();
  endrun.pop_back();
  vector<Int_t> Color;
  Color.push_back(2);
  Color.push_back(3);
  Color.push_back(4);
  Color.push_back(6);
  Color.push_back(7);
  Color.push_back(8);
  Color.push_back(9);
  Color.push_back(12);
  Color.push_back(28);
  Color.push_back(32);
  Color.push_back(37);
  Color.push_back(38);
  Color.push_back(41);
  Color.push_back(44);
  Color.push_back(46);
  Color.push_back(49);

  vector<Int_t> Style;
  Style.push_back(20);
  Style.push_back(24);
  Style.push_back(21);
  Style.push_back(25);
  Style.push_back(22);
  Style.push_back(26);
  Style.push_back(23);
  Style.push_back(32);
  Style.push_back(33);
  Style.push_back(27);
  Style.push_back(34);
  Style.push_back(28);
  Style.push_back(29);
  Style.push_back(30);
  
  ofstream fout("gainjump.txt");
  ofstream fgain("gain.txt",ios::app);
  ///////////////////////////////////////////
  //GATAutoCal ds(startrun.at(0),startrun.at(0));
  MJDSkim ds(fDataSet,0,0,1);
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos= ds.GetDetPosition();
  //vector<Int_t> fGoodBad = ds.GetGoodBad();
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> peak;
  vector<Double_t> peakerr;
  vector<Double_t> deltapeak;
  vector<Double_t> deltapeakerr;
  vector<Double_t> run;
  vector<Double_t> runerr;
  vector<Double_t> runtemp;
  vector<Double_t> runerrtemp;
  vector<Double_t> scale;
  vector<Double_t> offset;
  vector<Double_t> scaleerr;
  vector<Double_t> offseterr;
  //const Int_t channels = fChannel.size();
  string Pos;
  Int_t fPos;
  string FileName;
  string TitleName;
  cout.precision(12);
  vector<Double_t> sigma;
  for(size_t i=0;i<fChannel.size();i++){
    fPos = fCryo.at(i)*100+fStr.at(i)*10+fDetpos.at(i);
    Pos = Form("%d%d%d",fCryo.at(i),fStr.at(i),fDetpos.at(i));
    FileName = Form("DS%d_%s_%s_%d",fDataSet,fEName.c_str(),Pos.c_str(),fChannel.at(i));
    TitleName = Form("Data Set %d, C%dP%dD%d, Channel %d, Energy %s", fDataSet, fCryo.at(i),fStr.at(i),fDetpos.at(i),fChannel.at(i),fEName.c_str());
    scale.clear();
    offset.clear();
    peak.clear();
    peakerr.clear();
    deltapeak.clear();
    deltapeakerr.clear();
    run.clear();
    runerr.clear();
    runtemp.clear();
    runerrtemp.clear();
    //if(fGoodBad.at(i)>0){
      for(size_t j=0;j<startrun.size();j++){	
	Par.clear();
	ParErr.clear();
        Int_t run1 = startrun.at(j);
        Int_t run2 = endrun.at(j);
	GetPar(run1,run2,fEName,fPos,fChannel.at(i),fRefPar,&Par,&ParErr);
       	Int_t pars = Par.size();
	if(pars>0){
	  runtemp.push_back(startrun.at(j));
	  runerrtemp.push_back(endrun.at(j));
	  peak.push_back(Par.at(pars-1));
	  peakerr.push_back(ParErr.at(pars-1));
	}
      }
      //cout << scale.size() << " " << peak.size() << endl;
      if(peak.size()>0){
	for(size_t j=0;j<peak.size()-1;j++){
	  Double_t delta = peak.at(j)-2614.533;
	  
	  if(fChannel.at(i)%2==0 && abs(delta)<5){
	    sigma.push_back(delta);	  
	  }
	  if(abs(delta)>1){
	    cout << runtemp.at(j) << " "<< delta << endl;
	  }
	}
      }
      //}
  }
  //cout << sigma.size() << endl;
  Double_t sum0 = 0;
  Double_t sum1 = 0;
  Double_t sum2 = 0;
  Double_t n = 0;

  for(size_t ip=0;ip<sigma.size();ip++){
    sum0 = sum0 + sigma.at(ip);
    //sum1 = sum1 + pow(sigma.at(ip),2);
    //sum2 = sum2 + pow(sigma.at(ip),4);
    n = n +1;
  }
  Double_t ave0 = sum0/n;
  
  for(size_t ip=0;ip<sigma.size();ip++){
    //sum0 = sum0 + sigma.at(ip);
    sum1 = sum1 + pow(sigma.at(ip)-ave0,2);
    sum2 = sum2 + pow(sigma.at(ip)-ave0,4);
    //n = n +1;
    //cout << ip << " " << sigma.at(ip) << endl;
  }

  Double_t sig1 = sum1/n;
  Double_t sig2 = sum2/n;
  Double_t s1 = sqrt((sig2-(n-3)/(n-1)*sig1*sig1)/n);
  Double_t s2 = sqrt(pow((0.5*ave0),2)+sig1/(4*n));
  cout << "DS" << fDataSet << " " << ave0 << " "<< sig1 <<  " " << sig2 << " "<< n << " " << s1 <<  " " << s2 << endl;
  fgain << "DS" << fDataSet << " " << ave0 << " " << sig1 <<  " " << sig2 << " "<< n << " " << s1 << " " << s2 << endl;


}


