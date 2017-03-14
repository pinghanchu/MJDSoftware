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
#include "TFile.h"
#include <string>
#include <stdlib.h>

Double_t Flat(Double_t *v, Double_t *par){
  Double_t fitval = par[0];
  return fitval;
}


void GetCal(Int_t fStartRun, Int_t fEndRun, Int_t fChannel, string fEName, vector<Double_t>* Par){
  ifstream fin(Form("./List/cal/calibration_%d_%d.txt",fStartRun,fEndRun));

  Int_t cha,pos;
  string R[4];
  Double_t r[4];
  Int_t run1,run2;
  string ename;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> pos >> cha >> run1 >> run2 >> ename >> R[0] >> R[1] >> R[2] >>R[3];
      for(Int_t  i = 0;i<4;i++){
        if(R[i]=="nan" || R[i] =="-nan" || R[i] == "inf" || R[i] == "-inf"){
          r[i] = 0;
        }else{
          r[i] = atof(R[i].c_str());
        }
      }
      if(cha == fChannel && ename == fEName.c_str()){
	Par->push_back( r[0] );  // offset
	Par->push_back( r[1] );
	Par->push_back( r[2] ); // slope
	Par->push_back( r[3] );
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
  if(argc != 4 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [energy name]" << endl;
    return 1;
  }
  int fStartRun = atoi(argv[1]);
  int fEndRun = atoi(argv[2]);
  const char* EnergyName = argv[3];
  string fEName(EnergyName);
  //gStyle->SetOptStat(1100);
  //gStyle->SetOptFit(1111);
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  //TCanvas *c1 = new TCanvas("c1");
  ofstream fmiss(Form("./List/miss.cal_%d_%d.txt",fStartRun,fEndRun));
  ofstream falign(Form("./List/align.cal_%d_%d.txt",fStartRun,fEndRun));
  ofstream fbad(Form("./List/bad.cal_%d_%d.txt",fStartRun,fEndRun));
  ofstream foffseterr(Form("./List/offseterr.cal_%d_%d.txt",fStartRun,fEndRun));
  ifstream fin1("./List/runlist/cal.all.txt");
  Double_t offsetthreshold = 0.4;
  vector<Int_t> startrun;
  vector<Int_t> endrun;
  Int_t r1,r2,r3,r4;
  if(fin1.is_open()){
    while(!fin1.eof()){
      fin1 >> r1 >> r2 >> r3 >> r4 ;
      if(r1>= fStartRun && r2<= fEndRun){
	cout << r1 << " " << r2 << endl;
	startrun.push_back(r1);
	endrun.push_back(r2);
      }
    }
  }


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


  GATAutoCal ds(fStartRun,fStartRun);
  ds.SetEnergyName(EnergyName);
  string fDataSet = ds.GetDataSet();
  //  ds.SetParameters();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos= ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  const Int_t channels = fChannel.size();
  string Pos;
  string TitleName;
  string PlotName;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> offset;
  vector<Double_t> offseterr;
  vector<Double_t> slope;
  vector<Double_t> slopeerr;
  vector<Double_t> run;
  vector<Double_t> runerr;
  vector<Int_t> GoodBadIndex;

  TGraphErrors *fOffset[channels];
  TGraphErrors *fSlope[channels];
  for(Int_t i=0;i<channels;i++){
    Pos = Form("%d%d%d",fCryo.at(i),fStr.at(i),fDetpos.at(i));
    cout << Pos.c_str() << " " << fChannel.at(i) << " " << fGoodBad.at(i) << endl;
    if(fGoodBad.at(i)>0){
      //      c1->Update();
      run.clear();
      runerr.clear();
      offset.clear();
      offseterr.clear();
      slope.clear();
      slopeerr.clear();
      for(size_t j = 0 ; j<startrun.size();j++){
	Par.clear();        
	GetCal(startrun.at(j),endrun.at(j),fChannel.at(i),fEName,&Par);
	Int_t pars = Par.size();
	Double_t r1=(Double_t) (startrun.at(j)+endrun.at(j))/2.;
	Double_t r2=(Double_t) (-startrun.at(j)+endrun.at(j))/2.;

	run.push_back(r1);
	runerr.push_back(r2);

	if(Par.size()>0){	  
	  offset.push_back(Par.at(pars-4));
	  offseterr.push_back(Par.at(pars-3));
	  slope.push_back(Par.at(pars-2));
	  slopeerr.push_back(Par.at(pars-1));
	}else{
	  offset.push_back(9999999);
          offseterr.push_back(0);
          slope.push_back(0);
          slopeerr.push_back(0);
	  fmiss << startrun.at(j)<< " "<< endrun.at(j) << " " << Pos.c_str() << " " << fChannel.at(i) << endl;	
	}
      }
      /*
      if(offseterr.at(0)>0.1 || slopeerr.at(0)>0.1){
	falign << startrun.at(0) << " "<< endrun.at(0) << " "<< Pos.c_str() << " " << fChannel.at(i) << " " << slope.at(0) << " " << slopeerr.at(0)<< " " << offset.at(0) << " " << offseterr.at(0) << endl;
      }
      if(fDataSet =="P3JDY" && fChannel.at(i) == 663){
	if(slope.at(0)>3.){
	  fbad << startrun.at(0) << " " << endrun.at(0) << " " << Pos.c_str() << " " << fChannel.at(i) << endl;
	  slope.at(0) = 0;
	  slopeerr.at(0) = 0;
	  offset.at(0) = 999999;
	  offseterr.at(0) = 0;
	}
      }else if((fChannel.at(i)%2==0 && slope.at(0)>0.5) || (fChannel.at(i)%2==1 && slope.at(0)>1.5)){
	fbad << startrun.at(0) << " " << endrun.at(0) << " " << Pos.c_str() << " " << fChannel.at(i) << endl;
	slope.at(0) = 0;
	slopeerr.at(0) = 0;
	offset.at(0) = 999999;
	offseterr.at(0) = 0;
      }else if(abs(offset.at(0))>offsetthreshold){
	foffseterr << startrun.at(0) << " " << endrun.at(0) << " " << Pos.c_str() << " " << fChannel.at(i) << " " << offset.at(0) << endl;
      }
      */
	
      for(size_t j = 0 ; j<slope.size();j++){
	Double_t s0,s1,o0,o1;
	if(j==0){
	  s0 = slope.at(j+1);
	  s1 = s0;
	  o0 = offset.at(j+1);
	  o1 = o1;
	}else if(j==slope.size()-1){
	  s0 = slope.at(j-1);
	  s1 = s0;
	  o0 = offset.at(j-1);
	  o1 = o0;
	}else{
	  s0 = slope.at(j-1);
	  s1 = slope.at(j+1);
	  o0 = offset.at(j-1);
	  o1 = offset.at(j+1);
	}
	Double_t se = slopeerr.at(j);
	Double_t oe = offseterr.at(j);
	//cout << startrun.at(j) << " " << endrun.at(j) << " " << Pos.c_str() << " " << fChannel.at(i) << " " <<slope.at(j) << endl;
	if(fDataSet =="P3JDY" && fChannel.at(i) == 663){
	  if(slope.at(j)>3.){
	    fbad << startrun.at(j) << " " << endrun.at(j) << " " << Pos.c_str() << " " << fChannel.at(i) << " " << slope.at(j) << endl;
	    slope.at(j) = 0;
	    slopeerr.at(j) = 0;
	    offset.at(j) = 999999;
	    offseterr.at(j) = 0;
	  }
	}else if((fChannel.at(i)%2==0 && slope.at(j)>0.5) || (fChannel.at(i)%2==1 && slope.at(j)>1.5)){
	  fbad << startrun.at(j) << " " << endrun.at(j) << " " << Pos.c_str() << " " << fChannel.at(i) << " " << slope.at(j) << endl;
	  slope.at(j) = 0;
	  slopeerr.at(j) = 0;
	  offset.at(j) = 999999;
	  offseterr.at(j) = 0;
	}else if( abs(offset.at(j))>offsetthreshold) {
	  foffseterr << startrun.at(j) << " "<< endrun.at(j) << " "<< Pos.c_str() << " " << fChannel.at(i) << " " << offset.at(j) << " "
		     << endl;
	}else if((abs(slope.at(j)-s0)>0.1 && abs(slope.at(j)-s1)>0.1) || (abs(offset.at(j)-o0)>0.5 && abs(offset.at(j)-o1)>0.5) || oe>0.1 || se>0.1){
          falign << startrun.at(j) << " "<< endrun.at(j) << " "<< Pos.c_str() << " " << fChannel.at(i) << " " 
		 << slope.at(j) << " " << slopeerr.at(j) << " "  << offset.at(j) << " " << offseterr.at(j) << " "
		 << s0 << " " << o0 << " " << s1 << " " << o1 <<  endl;
	}
      }
      vector<Double_t> run1;
      vector<Double_t> runerr1;
      vector<Double_t> offset1;
      vector<Double_t> offseterr1;
      vector<Double_t> slope1;
      vector<Double_t> slopeerr1;
      for(size_t j = 0;j<slope.size();j++){
	ofstream fnew(Form("./List/cal/cal1_%d_%d.txt",startrun.at(j),endrun.at(j)),ios::app);
	fnew << Pos.c_str() << " " << fChannel.at(i) << " " << startrun.at(j) << " " << endrun.at(j) << " " << fEName.c_str() << " " 
	     << offset.at(j) << " "<< offseterr.at(j) << " "<< slope.at(j) << " " << slopeerr.at(j) << endl;	
	fnew.close();
	if(slope.at(j)==0 && offset.at(j) == 999999){
	}else{
	  slope1.push_back(slope.at(j));
	  offset1.push_back(offset.at(j));
	  slopeerr1.push_back(slopeerr.at(j));
	  offseterr1.push_back(offseterr.at(j));
	  run1.push_back(run.at(j));
	  runerr1.push_back(runerr.at(j));
	}
      }
      if(run.size()>0){
	GoodBadIndex.push_back(i);	//TF1 *fun = new TF1("fun",Flat,fStartRun,fEndRun,1);	
	fOffset[i] = GetGraph(run1,runerr1,offset1,offseterr1);
	Int_t ic = fCryo.at(i);
	Int_t is = fStr.at(i);
	Int_t id = fDetpos.at(i);
	cout << ic << is <<id << " "<< fChannel.at(i) << endl;
	fOffset[i]->SetTitle(Form("C%dP%dD%d, Channel %d;Run;Offset",ic,is,id,fChannel.at(i)));
	fOffset[i]->GetYaxis()->SetTitleOffset(1.5);
	fOffset[i]->SetMarkerStyle(2);
	fOffset[i]->SetMarkerSize(1);
	PlotName = Form("./Plot/offset_%s_%s_%d_%d_%d.pdf", fEName.c_str(),Pos.c_str(),fChannel.at(i), fStartRun,fEndRun);
	//fPeak[i]->SetMaximum(imax*1.01);
	//fPeak[i]->SetMinimum(imin*0.99);
	ds.PlotGraph(fOffset[i], PlotName);

        fSlope[i] = GetGraph(run1,runerr1,slope1,slopeerr1);
        fSlope[i]->SetTitle(Form("C%dP%dD%d, Channel %d;Run;Slope",ic,is,id,fChannel.at(i)));
        fSlope[i]->GetYaxis()->SetTitleOffset(1.5);
        fSlope[i]->SetMarkerStyle(2);
        fSlope[i]->SetMarkerSize(1);
        PlotName = Form("./Plot/slope_%s_%s_%d_%d_%d.pdf", fEName.c_str(),Pos.c_str(),fChannel.at(i), fStartRun,fEndRun);
        //fPeak[i]->SetMaximum(imax*1.01);
        //fPeak[i]->SetMinimum(imin*0.99);
        ds.PlotGraph(fSlope[i], PlotName);
      }
    }
  }
  
  TMultiGraph *mg1 = new TMultiGraph();
  TMultiGraph *mg2 = new TMultiGraph();
  TMultiGraph *mg3 = new TMultiGraph();
  TMultiGraph *mg4 = new TMultiGraph();

  TLegend *leg1 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg2 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg3 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg4 = new TLegend(0.15,0.55,0.85,0.85);


  gStyle->SetLegendBorderSize(0);

  for(size_t i=0;i<GoodBadIndex.size();i++){
    Int_t j = GoodBadIndex.at(i);
    Int_t ic = fCryo.at(j);
    Int_t is = fStr.at(j);
    Int_t id = fDetpos.at(j);
    Int_t hg = id*2;
    Int_t lg = hg+1;
    Int_t chan = fChannel.at(j);
    fOffset[j]->SetMarkerSize(1);
    fOffset[j]->SetMarkerColor(Color.at(is));
    fOffset[j]->SetLineColor(Color.at(is));
    fSlope[j]->SetMarkerSize(1);
    fSlope[j]->SetMarkerColor(Color.at(is));
    fSlope[j]->SetLineColor(Color.at(is));

    if(chan%2==0){
      fOffset[j]->SetMarkerStyle(Style.at(hg));
      fOffset[j]->SetLineStyle(0);
      fSlope[j]->SetMarkerStyle(Style.at(hg));
      fSlope[j]->SetLineStyle(0);
    }else{
      fOffset[j]->SetMarkerStyle(Style.at(lg));
      fOffset[j]->SetLineStyle(2);
      fSlope[j]->SetMarkerStyle(Style.at(lg));
      fSlope[j]->SetLineStyle(2);
    }

    TitleName = Form("#splitline{C%dP%dD%d}{%d}",ic,is,id,chan);
    if(chan%2 == 0){
      mg1->Add(fOffset[j],"ep");
      leg1->AddEntry(fOffset[j],Form("%s",TitleName.c_str()));
      mg3->Add(fSlope[j],"ep");
      leg3->AddEntry(fSlope[j],Form("%s",TitleName.c_str()));
    }else{
      mg2->Add(fOffset[j],"ep");
      leg2->AddEntry(fOffset[j],Form("%s",TitleName.c_str()));
      mg4->Add(fSlope[j],"ep");
      leg4->AddEntry(fSlope[j],Form("%s",TitleName.c_str()));
    }
  }
  gStyle->SetLegendTextSize(0.02);
  leg1->SetFillStyle(0);
  leg1->SetNColumns(6);
  leg1->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetNColumns(6);
  leg2->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetNColumns(6);
  leg3->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetNColumns(6);
  leg4->SetBorderSize(0);

  string TitleName1 = ";Run;Offset";
  string TitleName2 = ";Run;Offset";
  string TitleName3 = ";Run;Slope";
  string TitleName4 = ";Run;Slope";
  
  mg1->SetMaximum(0.3);
  mg2->SetMaximum(0.3);
  mg3->SetMaximum(0.5);
  mg4->SetMaximum(1.6);
  mg1->SetMinimum(-0.1);
  mg2->SetMinimum(-0.1);
  mg3->SetMinimum(0.35);
  mg4->SetMinimum(1.3);
  PlotName = Form("./Plot/offsethg_%s_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun);
  ds.PlotMultiGraph(mg1, leg1, "", PlotName,TitleName1);

  PlotName = Form("./Plot/offsetlg_%s_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun);
  ds.PlotMultiGraph(mg2, leg2, "", PlotName,TitleName2);

  PlotName = Form("./Plot/slopehg_%s_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun);
  ds.PlotMultiGraph(mg3, leg3, "", PlotName,TitleName3);

  PlotName = Form("./Plot/slopelg_%s_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun);
  ds.PlotMultiGraph(mg4, leg4, "", PlotName,TitleName4);

}


