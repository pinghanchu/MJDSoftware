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

Double_t Flat(Double_t *v, Double_t *par){
  Double_t fitval = par[0];
  return fitval;
}

Double_t IsNan(string Input){
  Double_t output = 0;
  if(Input == "nan" || Input == "-nan" || Input == "inf" || Input == "-inf"){
    output = 0;
  }else{
    output = atof(Input.c_str());
  }
  return output;
}

void GetCal(Int_t fStartRun, Int_t fEndRun, Int_t fChannel, string fEName, Int_t fIndex,vector<Double_t>* Par, vector<Double_t>* ParErr){

  //cout << "Read: "<<fStartRun << " " << fEndRun << " "<< fChannel << " " << fEName.c_str() << " " << fIndex << endl;
  ifstream fin2(Form("./List/cal/ezfit_%d_%d.txt",fStartRun,fEndRun));
  string x0,x1,x2,x3,x4;
  Double_t r1,r2,r3,r4;
  Int_t run1,run2;
  string cha;
  Int_t pos,index;
  string ename;
  Int_t chan;
  if(fin2.is_open()){
    while(!fin2.eof()){
      fin2 >> pos>> cha >> run1>>run2>>ename >> index >> x0 >> x1 >> x2 >>x3 >>x4;      
      r1 = IsNan(x1);
      r2 = IsNan(x2);
      r3 = IsNan(x3);
      r4 = IsNan(x4);
      chan = atoi(cha.c_str());
      //cout << "Load: " << cha<< " " << ename << " " << index << " " << r1 << " " << r2 << " " << r3 <<  " " << r4 << endl;
      if(chan == fChannel && ename == fEName && index==fIndex){	
      //if(cha == fChannel && index == fIndex){
	//cout << "Match: " << cha<< " " << ename << " " << index << " " << r1 << " " << r2 << " " << r3 <<  " " << r4 << endl;  
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
  if(argc != 5 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [energy name] [peak index]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  const char* EnergyName = argv[3];
  Int_t fIndex = atoi(argv[4]);
  string fEName(EnergyName);
  //gStyle->SetOptStat(1100);
  //gStyle->SetOptFit(1111);
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  ofstream fmiss(Form("./List/miss.ezfit_%d_%d.txt",fStartRun,fEndRun));
  ofstream falign(Form("./List/align.ezfit_%d_%d.txt",fStartRun,fEndRun));
  ifstream fin1("./List/runlist/cal.all.txt");
  vector<Int_t> startrun;
  vector<Int_t> endrun;
  Int_t r1,r2,r3,r4;
  if(fin1.is_open()){
    while(!fin1.eof()){
      fin1 >> r1 >> r2 >>r3 >>r4;
      if(r1>=fStartRun && r2 <=fEndRun){
	//cout << r1 << " " << r2 << endl;
	startrun.push_back(r1);
	endrun.push_back(r2);
      }
    }
  }
  if(startrun.at(startrun.size()-1)==startrun.at(startrun.size()-2)){
    startrun.pop_back();
    endrun.pop_back();
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

  ///////////////////////////////////////////
  GATAutoCal ds(fStartRun,fStartRun+1);
  ds.SetEnergyName(EnergyName);
  ds.SetParameters();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos= ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  const Int_t channels = fChannel.size();
  ofstream fpeak(Form("./List/ezfit_%d_%d.txt",fStartRun,fEndRun),ios::app);
  string Pos;
  string TitleName;
  string PlotName;
  vector<Double_t> Par;
  vector<Double_t> ParErr;
  vector<Double_t> peak;
  vector<Double_t> peakerr;
  vector<Double_t> fwhm;
  vector<Double_t> fwhmerr;
  vector<Double_t> run;
  vector<Double_t> runerr;
  vector<Int_t> GoodBadIndex;
  TGraphErrors *fPeak[channels];
  TGraphErrors *fPeakDelta[channels];
  TGraphErrors *fFWHM[channels];
  for(Int_t i=0;i<channels;i++){
    Pos = Form("%d%d%d",fCryo.at(i),fStr.at(i),fDetpos.at(i));
    cout << Pos.c_str() << " " << fChannel.at(i) << " " << fStartRun << " " << fEndRun << endl;
    ofstream fchan(Form("./List/cal/ezfitchan_%d_%d_%s_%d.txt",fStartRun,fEndRun,Pos.c_str(),fChannel.at(i)),ios::app);
    if(fGoodBad.at(i)>0){
      run.clear();
      runerr.clear();
      peak.clear();
      peakerr.clear();
      fwhm.clear();
      fwhmerr.clear();
      for(Int_t j=0;j<(int)startrun.size();j++){
        Int_t run1 = startrun.at(j);
        Int_t run2 = endrun.at(j);
	Par.clear();
	ParErr.clear();
	GetCal(run1,run2,fChannel.at(i),fEName,fIndex, &Par,&ParErr);
	Int_t pars = Par.size();
	if(Par.size()>0){
	  //cout << run1 << " " << run2 <<  " " << fEName << " " << fChannel.at(i) << " " << Par.at(pars-2) << " " << Par.at(pars-1) << endl; 
	  run.push_back(Par.at(pars-3));
	  peak.push_back(Par.at(pars-2));
	  runerr.push_back(ParErr.at(pars-3));
	  peakerr.push_back(ParErr.at(pars-2));
	  fwhm.push_back(Par.at(pars-1));
	  fwhmerr.push_back(ParErr.at(pars-1));
	}else{
	  fmiss << run1 << " " << run2 << " " << fEName.c_str() << " " << Pos.c_str() << " " << fChannel.at(i) << endl;
	}
      }
      /*
      for(size_t irun=1;irun<run.size()-1;irun++){
	Int_t run1 = run.at(irun)-runerr.at(irun);
	Int_t run2 = run.at(irun)+runerr.at(irun);
	Double_t peak1;
	Double_t peak2;
	if(irun == 0){
	  peak1 = peak.at(irun+1);
	  peak2 = peak2;
	}else if(irun==run.size()-1){
	  peak1 = peak.at(irun-1);
	  peak2 = peak1;
	}else{
	  peak1 = peak.at(irun-1);
	  peak2 = peak.at(irun+1);
	}
	if( (abs(peak.at(irun)-peak1) >30 && abs(peak.at(irun)-peak2)>30) || peakerr.at(irun)>1 ){
	  falign << run1 << " " << run2 << " " 
		 << Pos.c_str() << " " << fChannel.at(i) << " "<< fEName.c_str() << " " 
		 << peak.at(irun) << " " << peakerr.at(irun) << " "
		 << peak1<< " " << peak2<< " "
		 << endl;
	}
	if( (fChannel.at(i)%2==0 && peak.at(irun)< 6000) || (fChannel.at(i)%2==0 && peak.at(irun)< 1700)){
	  peak.at(irun) =0;
	  peakerr.at(irun) = 0;
	  fwhm.at(irun) = 0;
	  fwhmerr.at(irun)= 0; 
	}	
      }
      */
      for(size_t j = 0;j<peak.size();j++){
        ofstream fnew(Form("./List/cal/ezfit1_%d_%d.txt",startrun.at(j),endrun.at(j)),ios::app);
        fnew << Pos.c_str() << " " << fChannel.at(i) << " " << startrun.at(j) << " " << endrun.at(j) << " " << fEName.c_str() << " "
	     << fIndex << " " 
	     << peak.at(j) << " " << peakerr.at(j) << " " << fwhm.at(j) << " "<< fwhmerr.at(j) << endl;
        fnew.close();
        fchan << Pos.c_str() << " " << fChannel.at(i) << " " << startrun.at(j) << " " << endrun.at(j) << " " 
	      << fEName.c_str() << " "  << fIndex << " "
             << peak.at(j) << " " << peakerr.at(j) << " " << fwhm.at(j) << " "<< fwhmerr.at(j) << endl;
        fnew.close();
      }

      /*
      Double_t max = 0;
      Double_t min = 100000;
      
      for(size_t ip = 0;ip<peak.size();ip++){
	if(peak.at(ip)>max){
	  max = peak.at(ip);
	}
	if(peak.at(ip)<min){
	  min = peak.at(ip);
	}
      }
      fpeak << fStartRun << " " << fEndRun << " " << Pos.c_str() << " "<< fChannel.at(i) << " " << (max+min)/2 << " " << max-min << endl;
      */
      
      vector<Double_t> run1;
      vector<Double_t> runerr1;
      vector<Double_t> peak1;
      vector<Double_t> peakerr1;
      vector<Double_t> fwhm1;
      vector<Double_t> fwhmerr1;
      Double_t peakthresholdHG = 6000;
      Double_t peakthresholdLG = 1700;
      run1.clear();
      runerr1.clear();
      peak1.clear();
      peakerr1.clear();
      fwhm1.clear();
      fwhmerr1.clear();
      string Unit;
      for(size_t irun=0;irun<run.size();irun++){
	Int_t r1 = run.at(irun)-runerr.at(irun);
        Int_t r2 = run.at(irun)+runerr.at(irun);
  
	if(fEName  == "trapENFCal" || fEName == "trapECal" || fEName == "trapENMCal"){
	  Unit = "KeV";
	  run1.push_back(run.at(irun));
          runerr1.push_back(runerr.at(irun));
          peak1.push_back(peak.at(irun));
          peakerr1.push_back(peakerr.at(irun));	    
	  fwhm1.push_back(fwhm.at(irun));
	  if(fwhmerr.at(irun)>1){
	    fwhmerr1.push_back(sqrt(fwhm.at(irun)));
	  }else{
	    fwhmerr1.push_back(fwhmerr.at(irun));
	  }
	  if(peak.at(irun)< 2610 || peak.at(irun)>2620){
	    falign << r1 << " " << r2 << " " << fEName.c_str() << " " << Pos.c_str() << " " << fChannel.at(i) << " " << peak.at(irun)<< endl;	    
	  }
	}else if((fChannel.at(i)%2==0 &&peak.at(irun)>peakthresholdHG) || (fChannel.at(i)%2==1 && peak.at(irun)>peakthresholdLG)){
	  Unit = "ADC";
	  run1.push_back(run.at(irun));
	  runerr1.push_back(runerr.at(irun));
	  peak1.push_back(peak.at(irun));
	  peakerr1.push_back(peakerr.at(irun));
          fwhm1.push_back(fwhm.at(irun));
          fwhmerr1.push_back(fwhmerr.at(irun));
	}
      }

      if(run1.size()>0){
	TF1 *f = new TF1("f",Flat,fStartRun,fEndRun,1);
	GoodBadIndex.push_back(i);
	fPeak[i] = GetGraph(run1,runerr1,peak1,peakerr1);	
	fPeak[i]->Fit(f);
	Double_t mean = f->GetParameter(0);
	Int_t ic = fCryo.at(i);
	Int_t is = fStr.at(i);
	Int_t id = fDetpos.at(i);
	fPeak[i]->SetTitle(Form("C%dP%dD%d, Channel %d;Run;Peak(%s)",ic,is,id,fChannel.at(i),Unit.c_str()));
	fPeak[i]->GetYaxis()->SetTitleOffset(1.5);
	fPeak[i]->SetMarkerStyle(2);
	fPeak[i]->SetMarkerSize(1);
	PlotName = Form("./Plot/peak_%s_%s_%d_%d_%d_%d.pdf", fEName.c_str(),Pos.c_str(),fChannel.at(i), fStartRun,fEndRun,fIndex);
	ds.PlotGraph(fPeak[i], PlotName);

	for(size_t ip = 0;ip<peak1.size();ip++){
	  peak1.at(ip) = peak1.at(ip)/mean;
	  peakerr1.at(ip)=peakerr1.at(ip)/mean;
	}
	fPeakDelta[i] = GetGraph(run1,runerr1,peak1,peakerr1);
	fPeakDelta[i]->SetTitle(Form("C%dP%dD%d, Channel %d;Run;Peak(%s)",ic,is,id,fChannel.at(i),Unit.c_str()));
        fPeakDelta[i]->GetYaxis()->SetTitleOffset(1.5);
        fPeakDelta[i]->SetMarkerStyle(2);
        fPeakDelta[i]->SetMarkerSize(1);
        PlotName = Form("./Plot/peakdelta_%s_%s_%d_%d_%d_%d.pdf", fEName.c_str(),Pos.c_str(),fChannel.at(i), fStartRun,fEndRun,fIndex);
        ds.PlotGraph(fPeakDelta[i], PlotName);

	fFWHM[i] = GetGraph(run1,runerr1,fwhm1,fwhmerr1);
        fFWHM[i]->SetTitle(Form("C%dP%dD%d, Channel %d;Run;FWHM(%s)",ic,is,id,fChannel.at(i),Unit.c_str()));
        fFWHM[i]->GetYaxis()->SetTitleOffset(1.5);
        fFWHM[i]->SetMarkerStyle(2);
        fFWHM[i]->SetMarkerSize(1);
        PlotName = Form("./Plot/fwhm_%s_%s_%d_%d_%d_%d.pdf", fEName.c_str(),Pos.c_str(),fChannel.at(i), fStartRun,fEndRun,fIndex);
        ds.PlotGraph(fFWHM[i], PlotName);
      }
    }
  }


  TMultiGraph *mg1 = new TMultiGraph();
  TMultiGraph *mg2 = new TMultiGraph();
  TMultiGraph *mg3 = new TMultiGraph();
  TMultiGraph *mg4 = new TMultiGraph();
  TMultiGraph *mg5 = new TMultiGraph();
  TMultiGraph *mg6 = new TMultiGraph();

  TLegend *leg1 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg2 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg3 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg4 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg5 = new TLegend(0.15,0.55,0.85,0.85);
  TLegend *leg6 = new TLegend(0.15,0.55,0.85,0.85);

  gStyle->SetLegendBorderSize(0);

  for(size_t i=0;i<GoodBadIndex.size();i++){
    Int_t j = GoodBadIndex.at(i);
    Int_t ic = fCryo.at(j);
    Int_t is = fStr.at(j);
    Int_t id = fDetpos.at(j);
    Int_t hg = id*2;
    Int_t lg = hg+1;
    Int_t chan = fChannel.at(j);
    fPeak[j]->SetMarkerSize(1);
    fPeak[j]->SetMarkerColor(Color.at(is));
    fPeak[j]->SetLineColor(Color.at(is));
    fPeakDelta[j]->SetMarkerSize(1);
    fPeakDelta[j]->SetMarkerColor(Color.at(is));
    fPeakDelta[j]->SetLineColor(Color.at(is));
    fFWHM[j]->SetMarkerSize(1);
    fFWHM[j]->SetMarkerColor(Color.at(is));
    fFWHM[j]->SetLineColor(Color.at(is));

    if(chan%2==0){
      fPeak[j]->SetMarkerStyle(Style.at(hg));
      fPeak[j]->SetLineStyle(0);
      fPeakDelta[j]->SetMarkerStyle(Style.at(hg));
      fPeakDelta[j]->SetLineStyle(0);
      fFWHM[j]->SetMarkerStyle(Style.at(hg));
      fFWHM[j]->SetLineStyle(0);
    }else{
      fPeak[j]->SetMarkerStyle(Style.at(lg));
      fPeak[j]->SetLineStyle(2);
      fPeakDelta[j]->SetMarkerStyle(Style.at(lg));
      fPeakDelta[j]->SetLineStyle(2);
      fFWHM[j]->SetMarkerStyle(Style.at(lg));
      fFWHM[j]->SetLineStyle(2);
    }

    TitleName = Form("#splitline{C%dP%dD%d}{%d}",ic,is,id,chan);
    if(chan%2 == 0){
      mg1->Add(fPeak[j],"ep");
      leg1->AddEntry(fPeakDelta[j],Form("%s",TitleName.c_str()));
      mg2->Add(fPeakDelta[j],"ep");
      leg2->AddEntry(fPeakDelta[j],Form("%s",TitleName.c_str()));
      mg5->Add(fFWHM[j],"ep");
      leg5->AddEntry(fFWHM[j],Form("%s",TitleName.c_str()));
    }else{
      mg3->Add(fPeak[j],"ep");
      leg3->AddEntry(fPeak[j],Form("%s",TitleName.c_str()));
      mg4->Add(fPeakDelta[j],"ep");
      leg4->AddEntry(fPeakDelta[j],Form("%s",TitleName.c_str()));
      mg6->Add(fFWHM[j],"ep");
      leg6->AddEntry(fFWHM[j],Form("%s",TitleName.c_str()));
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
  
  leg5->SetFillStyle(0);
  leg5->SetNColumns(6);
  leg5->SetBorderSize(0);

  leg6->SetFillStyle(0);
  leg6->SetNColumns(6);
  leg6->SetBorderSize(0);

  string Unit;
  if(fEName  == "trapENFCal" || fEName == "trapECal" || fEName == "trapENMCal"){
    Unit = "KeV";
  }else{
    Unit = "ADC";
  }

  string TitleName1 = Form(";Run;Peak(%s)",Unit.c_str());
  string TitleName2 =";Run;Peak/Ave(Peak)";
  string TitleName3 = Form(";Run;FWHM(%s)",Unit.c_str());

  mg1->Draw();
  //mg1->SetMaximum(7500);
  //mg1->SetMinimum(6150);
  PlotName = Form("./Plot/peakhg_%s_%d_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun,fIndex);
  ds.PlotMultiGraph(mg1, leg1, "", PlotName,TitleName1);

  mg2->Draw();
  mg2->SetMaximum(1.002);
  mg2->SetMinimum(0.999);
  PlotName = Form("./Plot/peakdeltahg_%s_%d_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun,fIndex);
  ds.PlotMultiGraph(mg2, leg2, "", PlotName,TitleName2);

  mg3->Draw();
  //mg3->SetMaximum(2400);
  //mg3->SetMinimum(1800);
  
  PlotName = Form("./Plot/peaklg_%s_%d_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun,fIndex);
  ds.PlotMultiGraph(mg3, leg3, "", PlotName,TitleName1);


  mg4->Draw();
  mg4->SetMaximum(1.002);
  mg4->SetMinimum(0.999);
  PlotName = Form("./Plot/peakdeltalg_%s_%d_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun,fIndex);
  ds.PlotMultiGraph(mg4, leg4, "", PlotName,TitleName2);

  mg5->Draw();
  mg5->SetMaximum(5);
  mg5->SetMinimum(2);
  PlotName = Form("./Plot/fwhmhg_%s_%d_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun,fIndex);
  ds.PlotMultiGraph(mg5, leg5, "", PlotName,TitleName3);

  mg6->Draw();
  mg6->SetMaximum(6);
  mg6->SetMinimum(2);
  PlotName = Form("./Plot/fwhmlg_%s_%d_%d_%d.pdf", fEName.c_str(), fStartRun,fEndRun,fIndex);
  ds.PlotMultiGraph(mg6, leg6, "", PlotName,TitleName3);


}


