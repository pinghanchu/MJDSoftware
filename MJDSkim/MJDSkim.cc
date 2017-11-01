#include "MJDSkim.hh"
#include "GATAutoCal.hh"
#include "MJTRun.hh"
#include "MJAnalysisDoc.hh"
#include "MJAnalysisParameters.hh"
#include "GATWaveformBrowser.hh"
#include "MGTWaveform.hh"
#include "GATPeakShape.hh"
#include "GATMultiPeakFitter.hh"
#include "GATHybridMonteCarlo.hh"
#include "MJProvenance.hh"
#include "MJDatabase.hh"
//#include "MJDBUtilitiesx2.hh"
#include "KTreeFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include "TMatrixDSym.h"
#include "TPaveStats.h"
#include "Fit/FitResult.h"
#include "Fit/FitConfig.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TVirtualFFT.h"
#include <iostream>
#include <bitset>
#include <iomanip> 
//ClassImp(MJDSkim)

using namespace std;
using namespace katrin;
using namespace MJDB;
//using namespace GATPeakShapeFunction;

////////////////////////////////////////////////////////

//Set ChannelMap, mjdTree (from GATDataSet), CalibrationPeak, Initial Parameters
MJDSkim::MJDSkim(Int_t DataSet, Int_t SubSet1, Int_t SubSet2, Int_t IsCal) 
{
  fDataSet = DataSet;
  fSubSet = SubSet1;
  fSubSet2 = SubSet2;
  fSkimTree = new TChain("skimTree");
  fIsCal = IsCal;
  string path = "GAT-v01-07";
  if(fSubSet == fSubSet2){
    if(IsCal == 0){
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%d/%s/skimDS%d_%d.root",fDataSet,path.c_str(),fDataSet,fSubSet));
    }else{
      fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%dcal/%s/skimDS%d_run%d_small.root",fDataSet,path.c_str(),fDataSet,fSubSet));
    }
  }else{
    for(Int_t is = SubSet1;is<=SubSet2;is++){
      if(IsCal == 0){
	fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%d/%s/skimDS%d_%d.root",fDataSet,path.c_str(),fDataSet,is));
      }else{
	fSkimTree->Add(Form("$MJDDATADIR/surfmjd/analysis/skim/DS%dcal/%s/skimDS%d_run%d_small.root",fDataSet,path.c_str(),fDataSet,is));
      }
    }
  }
  fSkimTree->SetBranchStatus("*",1);
  
  fChannel = NULL;
  fP = NULL;
  fD = NULL;
  fC = NULL;
  fGain = NULL;
  fIsEnr = NULL;
  fIsNat = NULL;
  fIsGood = NULL;
  fTrapENFCal = NULL;
  ftOffset = NULL;
  fDCR = NULL;
  fAvsE = NULL;
  fTrapETailMin = NULL;
  fnX = NULL;
  fdtmu_s = NULL;
  fglobalTime = NULL;
  fSkimTree->SetBranchAddress("run", &fRun);
  fSkimTree->SetBranchAddress("iEvent", &fEvent);
  fSkimTree->SetBranchAddress("channel",&fChannel);
  fSkimTree->SetBranchAddress("P",&fP);
  fSkimTree->SetBranchAddress("D",&fD);
  fSkimTree->SetBranchAddress("C",&fC);
  fSkimTree->SetBranchAddress("gain",&fGain);
  fSkimTree->SetBranchAddress("isEnr",&fIsEnr);
  fSkimTree->SetBranchAddress("isNat",&fIsNat);
  fSkimTree->SetBranchAddress("isGood",&fIsGood);
  fSkimTree->SetBranchAddress("trapENFCalC", &fTrapENFCal);
  fSkimTree->SetBranchAddress("globalTime",&fglobalTime);
  fSkimTree->SetBranchAddress("localTime_s",&flocalTime_s);
  fSkimTree->SetBranchAddress("clockTime_s",&fclockTime_s);
  fSkimTree->SetBranchAddress("tOffset",&ftOffset);
  fSkimTree->SetBranchAddress("dcr99",&fDCR);
  fSkimTree->SetBranchAddress("avse",&fAvsE);
  fSkimTree->SetBranchAddress("trapETailMin",&fTrapETailMin);
  fSkimTree->SetBranchAddress("mHL",&fmH);
  fSkimTree->SetBranchAddress("nX",&fnX);
  fSkimTree->SetBranchAddress("muVeto",&fmuVeto);
  fSkimTree->SetBranchAddress("dtmu_s",&fdtmu_s);
  fSkimTree->SetBranchAddress("muTUnc",&fmuTUnc);
  fSkimTree->GetEntry(0);
  fEntries = fSkimTree->GetEntries();
  fStartRun = fRun;
  GATDataSet *fDS = new GATDataSet(fStartRun);
  if((fStartRun>0 && fStartRun<30000000) || (fStartRun > 60000000 && fStartRun<65000000)){
    fMap = fDS->GetChannelMap();
  }
  SetChannel();
}


void MJDSkim::SetChannel()
{
  //cout << fStartRun << endl;
  fMTChannel.clear();
  fCryo.clear();
  fString.clear();
  fDetector.clear();
  fDetectorName.clear();
  
  //DS0-DS3 (Mod1)
  if(fStartRun>=0 && fStartRun <18623){
    for(Int_t i =1;i<2;i++){
      for(Int_t j = 1;j<=7;j++){
	if(j==4){
	  for(Int_t k = 1;k<=5;k++){
	    fCryo.push_back(i);
	    fCryo.push_back(i);
	    fString.push_back(j);
	    fString.push_back(j);
	    fDetector.push_back(k);
	    fDetector.push_back(k);
	    fMTChannel.push_back(fMap->GetIDHi(i,j,k));
	    fMTChannel.push_back(fMap->GetIDLo(i,j,k));
	    fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
            fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	  }
	}else{
	  for(Int_t k = 1;k<=4;k++){
	    fCryo.push_back(i);
	    fCryo.push_back(i);
	    fString.push_back(j);
	    fString.push_back(j);
	    fDetector.push_back(k);
	    fDetector.push_back(k);
	    fMTChannel.push_back(fMap->GetIDHi(i,j,k));
	    fMTChannel.push_back(fMap->GetIDLo(i,j,k));
            fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
            fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	  }
	}
      }      
    }
  }

  //DS4 (Mod2)
  if(fStartRun>=60000001 && fStartRun <70000000){
    for(Int_t i =2;i<3;i++){
      for(Int_t j = 1;j<=7;j++){
        if(j==2 || j==4){
          for(Int_t k = 1;k<=5;k++){
	    fCryo.push_back(i);
	    fCryo.push_back(i);
            fString.push_back(j);
            fString.push_back(j);
            fDetector.push_back(k);
            fDetector.push_back(k);
            fMTChannel.push_back(fMap->GetIDHi(i,j,k));
            fMTChannel.push_back(fMap->GetIDLo(i,j,k));
	    fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
            fDetectorName.push_back(fMap->GetDetectorName(i,j,k));

          }
        }else if(j==3){
          for(Int_t k = 1;k<=3;k++){
	    fCryo.push_back(i);
	    fCryo.push_back(i);
            fString.push_back(j);
            fString.push_back(j);
            fDetector.push_back(k);
            fDetector.push_back(k);
            fMTChannel.push_back(fMap->GetIDHi(i,j,k));
            fMTChannel.push_back(fMap->GetIDLo(i,j,k));
	    fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
            fDetectorName.push_back(fMap->GetDetectorName(i,j,k));

          }
	}else{
          for(Int_t k = 1;k<=4;k++){
	    fCryo.push_back(i);
	    fCryo.push_back(i);
            fString.push_back(j);
            fString.push_back(j);
            fDetector.push_back(k);
            fDetector.push_back(k);
            fMTChannel.push_back(fMap->GetIDHi(i,j,k));
            fMTChannel.push_back(fMap->GetIDLo(i,j,k));
	    fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
            fDetectorName.push_back(fMap->GetDetectorName(i,j,k));

          }
        }
      }
    }
  }

  //DS5- (Mod1+Mod2)
  if(fStartRun>=18623 && fStartRun <30000000){
    for(Int_t i =1;i<3;i++){
      //Mod1
      if(i==1){
	for(Int_t j = 1;j<=7;j++){
	  if(j==4){
	    for(Int_t k = 1;k<=5;k++){
	      fCryo.push_back(i);
	      fCryo.push_back(i);
	      fString.push_back(j);
	      fString.push_back(j);
	      fDetector.push_back(k);
	      fDetector.push_back(k);
	      fMTChannel.push_back(fMap->GetIDHi(i,j,k));
	      fMTChannel.push_back(fMap->GetIDLo(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	    }
	  }else{
	    for(Int_t k = 1;k<=4;k++){
	      fCryo.push_back(i);
	      fCryo.push_back(i);
	      fString.push_back(j);
	      fString.push_back(j);
	      fDetector.push_back(k);
	      fDetector.push_back(k);
	      fMTChannel.push_back(fMap->GetIDHi(i,j,k));
	      fMTChannel.push_back(fMap->GetIDLo(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	    }
	  }
	}      
      }
      //Mod2
      if(i==2){
	for(Int_t j = 1;j<=7;j++){
	  if(j==2 || j==4){
	    for(Int_t k = 1;k<=5;k++){
	      fCryo.push_back(i);
	      fCryo.push_back(i);
	      fString.push_back(j);
	      fString.push_back(j);
	      fDetector.push_back(k);
	      fDetector.push_back(k);
	      fMTChannel.push_back(fMap->GetIDHi(i,j,k));
	      fMTChannel.push_back(fMap->GetIDLo(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));	      
	    }
	  }else if(j==3){
	    for(Int_t k = 1;k<=3;k++){
	      fCryo.push_back(i);
	      fCryo.push_back(i);
	      fString.push_back(j);
	      fString.push_back(j);
	      fDetector.push_back(k);
	      fDetector.push_back(k);
	      fMTChannel.push_back(fMap->GetIDHi(i,j,k));
	      fMTChannel.push_back(fMap->GetIDLo(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));	      
	    }
	  }else{
	    for(Int_t k = 1;k<=4;k++){
	      fCryo.push_back(i);
	      fCryo.push_back(i);
	      fString.push_back(j);
	      fString.push_back(j);
	      fDetector.push_back(k);
	      fDetector.push_back(k);
	      fMTChannel.push_back(fMap->GetIDHi(i,j,k));
	      fMTChannel.push_back(fMap->GetIDLo(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));	      
	    }
	  }
	}
      }
    }
  }
}


void MJDSkim::FindEvent(Int_t List, vector<Int_t>* Channel, vector<Double_t>* Energy){
  //////////////////////////////////////


  fSkimTree->GetEntry(List);
  Int_t run1 = fRun;
  Int_t event1 = fEvent;

  for(size_t i=0;i<fChannel->size();i++){	      
    Int_t chan1 = fChannel->at(i);    
    Double_t enr1 = fTrapENFCal->at(i);
    Double_t dcr1 = fDCR->at(i);
    Double_t avse1 = fAvsE->at(i);
    Double_t trapetailmin1 = fTrapETailMin->at(i);
    Double_t time1 = fglobalTime->AsDouble()+ ftOffset->at(i)/1e9;
    Double_t dtmu_s1 = fdtmu_s->at(i);

    if((enr1>5 && enr1 < 12000 && chan1%2==0) && 
       !(0.015<dcr1 && dcr1<0.018 && 49<enr1 && enr1<51) &&
       !(0.004<dcr1 && dcr1<0.006 && 49<enr1 && enr1<52) &&
       !(-0.006<dcr1 && dcr1<-0.004 && 48<enr1 && enr1<58) &&
       !(avse1>500 && (chan1 == 691 || chan1 == 1124)) &&
       !(trapetailmin1 >2 && trapetailmin1 <3.5 && enr1>191 && enr1<192)
       ){
      Channel->push_back(chan1);
      Energy->push_back(enr1);
    }
  }
}




void MJDSkim::FindDelayedEvent(Int_t List, Int_t Chan, Double_t fTime, string fOutputFile){
  //////////////////////////////////////

  ofstream fout(Form("%s",fOutputFile.c_str()),ios::app);
  fout << fixed << setprecision(5);

  vector<Int_t> Run2;
  vector<Int_t> Run3;
  vector<Int_t> List2;
  vector<Int_t> List3;
  vector<Int_t> Chan2;
  vector<Int_t> Chan3;
  vector<Double_t> Enr2;
  vector<Double_t> Enr3;
  vector<Double_t> Time2;
  vector<Double_t> Time3;
  vector<Int_t> Event2;
  vector<Int_t> Event3;
  vector<Double_t> Mu_s2;
  vector<Double_t> Mu_s3;
  vector<Double_t> DCR2;
  vector<Double_t> DCR3;
  vector<Double_t> AvsE2;
  vector<Double_t> AvsE3;
  vector<Double_t> TrapETailMin2;
  vector<Double_t> TrapETailMin3;
  cout.precision(15);



  fSkimTree->GetEntry(List);
  Int_t run1 = fRun;
  Int_t event1 = fEvent;

  for(size_t i=0;i<fChannel->size();i++){	      
    Int_t chan1 = fChannel->at(i);
    if(chan1 == Chan){
      Double_t enr1 = fTrapENFCal->at(i);
      Double_t dcr1 = fDCR->at(i);
      Double_t avse1 = fAvsE->at(i);
      Double_t trapetailmin1 = fTrapETailMin->at(i);
      Double_t time1 = fglobalTime->AsDouble()+ ftOffset->at(i)/1e9;
      Double_t dtmu_s1 = fdtmu_s->at(i);
      cout << fDataSet << "  " << fSubSet << " "<< List << " " << fRun << " " << fEvent << " " << time1 << " " << chan1 << " " << enr1 << endl;
      Int_t ii = List-1;
      
      fSkimTree->GetEntry(ii);
      Int_t run2 = fRun;
      Double_t time2 = fglobalTime->AsDouble();
      Int_t event2 = fEvent;
      
      while( (time1-time2)<fTime*20 && (time1-time2) >=0 && ii>0){
	for(size_t j=0;j<fChannel->size();j++){
	  Int_t chan2 = fChannel->at(j);
	  Double_t enr2 = fTrapENFCal->at(j);	    
	  Double_t dtmu_s2 = fdtmu_s->at(j);
	  Double_t dcr2 = fDCR->at(j);
	  Double_t avse2 = fAvsE->at(j);
	  Double_t trapetailmin2 = fTrapETailMin->at(j);
	  time2 =  fglobalTime->AsDouble() + ftOffset->at(j)/1e9;
	  if((enr2>5 && enr2 < 12000 && chan2%2==0) && 
	     !(0.015<dcr2 && dcr2<0.018 && 49<enr2 && enr2<51) &&
	     !(0.004<dcr2 && dcr2<0.006 && 49<enr2 && enr2<52) &&
	     !(-0.006<dcr2 && dcr2<-0.004 && 48<enr2 && enr2<58) &&
	     !(avse2>500 && (chan2 == 691 || chan2 == 1124)) &&
	     !(trapetailmin2 >2 && trapetailmin2 <3.5 && enr2>191 && enr2<192)
	     ){
	    Run2.push_back(run1);
	    List2.push_back(List);
	    Chan2.push_back(chan1);
	    Enr2.push_back(enr1);
	    Time2.push_back(time1);
	    Event2.push_back(event1);
	    Mu_s2.push_back(dtmu_s1);
	    DCR2.push_back(dcr1);
	    AvsE2.push_back(avse1);
	    TrapETailMin2.push_back(trapetailmin1);
	    Run3.push_back(run2);
	    List3.push_back(ii);
	    Chan3.push_back(chan2);
	    Enr3.push_back(enr2);			     
	    Time3.push_back(time2);
	    Event3.push_back(event2);
	    Mu_s3.push_back(dtmu_s2);
	    DCR3.push_back(dcr2);
	    AvsE3.push_back(avse2);
	    TrapETailMin3.push_back(trapetailmin2);
	  }	     
	  ii--;
	  fSkimTree->GetEntry(ii);
	  run2 = fRun;
	  time2 = fglobalTime->AsDouble();
	  event2 = fEvent;	  
	}
      }
    }
  }

  for(size_t i=0;i<List2.size();i++){
    fout << fDataSet << " " << fSubSet << " "
	 << List2.at(i) << " " << Run2.at(i) << " " << Event2.at(i) << " " << Time2.at(i) << " " 
	 << Chan2.at(i) << " " << Enr2.at(i) << " "
	 << DCR2.at(i) << " "  << AvsE2.at(i) << " " << TrapETailMin2.at(i) << " " << Mu_s2.at(i) << " "      
	 << List3.at(i) << " " << Run3.at(i) << " " << Event3.at(i) << " " << Time3.at(i) << " " 
	 << Chan3.at(i) << " " << Enr3.at(i) << " " 
	 << DCR3.at(i) << " "  << AvsE3.at(i) << " " << TrapETailMin3.at(i) << " " << Mu_s3.at(i) << " "
	 << Time2.at(i) - Time3.at(i) << endl;
  }
}



void MJDSkim::SearchDelayedEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the second gamma energy
  //fEnr2 : the first Q value
  //////////////////////////////////////
  //cout << fTime << endl;
  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout << fixed << setprecision(5);

  vector<Int_t> Run1;
  vector<Int_t> List1;
  vector<Int_t> Chan1;
  vector<Double_t> Enr1;
  vector<Double_t> DCR1;
  vector<Double_t> AvsE1;
  vector<Double_t> TrapETailMin1;
  vector<Double_t> Time1;
  for(size_t i=0;i<fEntries;i++){
    fSkimTree->GetEntry(i);
    Int_t run1 = fRun;
    for(size_t j=0;j<fChannel->size();j++){	      
      Int_t chan1 = fChannel->at(j);
      Double_t enr1 = fTrapENFCal->at(j);
      Double_t dcr1 = fDCR->at(j);
      Double_t avse1 = fAvsE->at(j);
      Double_t trapetailmin = fTrapETailMin->at(j);
      Double_t time1 = fglobalTime->AsDouble()+ ftOffset->at(j)/1e9;
      if( (abs(enr1-fEnr1)<20 && chan1%2==0) 
	  && 
	  !(0.015<dcr1 && dcr1<0.018 && 49<enr1 && enr1<51) && 
	  !(0.004<dcr1 && dcr1<0.006 && 49<enr1 && enr1<52) &&
	  !(-0.006<dcr1 && dcr1<-0.004 && 48<enr1 && enr1<58) &&
	  !(avse1>500 && (chan1 == 691 || chan1 == 1124)) &&
	  !(trapetailmin >2 && trapetailmin <3.5 && enr1>191 && enr1<192)
	 ){ 
	Run1.push_back(run1);
	List1.push_back(i);
	Chan1.push_back(chan1);
	Enr1.push_back(enr1);
	DCR1.push_back(dcr1);
	AvsE1.push_back(avse1);
	TrapETailMin1.push_back(trapetailmin);
	Time1.push_back(time1);
      }
    }
  }

  vector<Int_t> Run2;
  vector<Int_t> Run3;
  vector<Int_t> List2;
  vector<Int_t> List3;
  vector<Int_t> Chan2;
  vector<Int_t> Chan3;
  vector<Double_t> Enr2;
  vector<Double_t> Enr3;
  vector<Double_t> Time2;
  vector<Double_t> Time3;
  vector<Int_t> Event2;
  vector<Int_t> Event3;
  vector<Double_t> Mu_s2;
  vector<Double_t> Mu_s3;
  vector<Double_t> DCR2;
  vector<Double_t> DCR3;
  vector<Double_t> AvsE2;
  vector<Double_t> AvsE3;
  vector<Double_t> TrapETailMin2;
  vector<Double_t> TrapETailMin3;
  cout.precision(15);
  if((Int_t)List1.size()>0){
    for(size_t i=0;i<List1.size();i++){
      fSkimTree->GetEntry(List1.at(i));
      Int_t run1 = Run1.at(i);
      //Double_t time1 = flocalTime_s;
      Double_t time1 = Time1.at(i);
      Double_t time2 = time1;
      Double_t dtmu_s1 = fdtmu_s->at(0);

      Int_t event1 = fEvent;
      Int_t ii = List1.at(i)-1;
      fSkimTree->GetEntry(ii);
      Int_t run2 = fRun;
      time2 = fglobalTime->AsDouble();
      Int_t event2 = fEvent;
      Double_t dcr2 = DCR1.at(i);
      Double_t avse2 = AvsE1.at(i);
      Double_t trapetailmin2 = TrapETailMin1.at(i);
      //cout<< run1 << " "<< event1 << " "<< Enr1.at(i) << endl;
      while( (time1-time2)<fTime*20 && (time1-time2) >=0 && ii>0){
	for(size_t j=0;j<fChannel->size();j++){
	  Int_t chan2 = fChannel->at(j);
	  Double_t enr2 = fTrapENFCal->at(j);	    
	  Double_t dtmu_s2 = fdtmu_s->at(j);
	  Double_t dt = ftOffset->at(j)/1e9;
	  Double_t dcr3 = fDCR->at(j);
	  Double_t avse3 = fAvsE->at(j);
	  Double_t trapetailmin3 = fTrapETailMin->at(j);
	  //cout << run1 << " "<< event1 << " " << time1 << " "<< event2 << " "<< time2 << " " << dt << " " << chan2 << endl;
	  if( (enr2 < fEnr2 && enr2>5 && chan2%2==0 && (time1 - (time2+dt))>0) 
	      &&
	      !(0.015<dcr3 && dcr3<0.018 && 49<enr2 && enr2<51) &&
	      !(0.004<dcr3 && dcr3<0.006 && 49<enr2 && enr2<52) &&
	      !(-0.006<dcr3 && dcr3<-0.004 && 48<enr2 && enr2<58) &&
	      !(avse3>500 && (chan2 == 691 || chan2 == 1124)) &&
	      !(trapetailmin3 >2 && trapetailmin3 <3.5 && enr2>191 && enr2<192)
	     ){
	    Run2.push_back(run1);
	    List2.push_back(List1.at(i));
	    Chan2.push_back(Chan1.at(i));
	    Enr2.push_back(Enr1.at(i));
	    Time2.push_back(time1);
	    Event2.push_back(event1);
	    Mu_s2.push_back(dtmu_s1);
	    DCR2.push_back(dcr2);
	    AvsE2.push_back(avse2);
	    TrapETailMin2.push_back(trapetailmin2);
	    Run3.push_back(run2);
	    List3.push_back(ii);
	    Chan3.push_back(chan2);
	    Enr3.push_back(enr2);			     
	    Time3.push_back(time2+dt);
            Event3.push_back(event2);
            Mu_s3.push_back(dtmu_s2);
	    DCR3.push_back(dcr3);
	    AvsE3.push_back(avse3);
	    TrapETailMin3.push_back(trapetailmin3);
	    //cout << time1 << " " << time2 << " "<< time1-time2 << endl;
	  }
	}     
	ii--;
	fSkimTree->GetEntry(ii);
        run2 = fRun;
        time2 = fglobalTime->AsDouble();
        event2 = fEvent;
      }
    }
  }


  for(size_t i=0;i<List2.size();i++){ 


    fout << fDataSet << " " << fSubSet << " "
	 << List2.at(i) << " " << Run2.at(i) << " " << Event2.at(i) << " " << Time2.at(i) << " " 
	 << Chan2.at(i) << " " << Enr2.at(i) << " "
	 << DCR2.at(i) << " "  << AvsE2.at(i) << " " << TrapETailMin2.at(i) << " " << Mu_s2.at(i) << " "

	 << List3.at(i) << " " << Run3.at(i) << " " << Event3.at(i) << " " << Time3.at(i) << " " 
	 << Chan3.at(i) << " " << Enr3.at(i) << " " 
	 << DCR3.at(i) << " "  << AvsE3.at(i) << " " << TrapETailMin3.at(i) << " " << Mu_s3.at(i) << " "
	 << Time2.at(i) - Time3.at(i) << endl;
  }
}



void MJDSkim::SearchDelayedHighEnergyEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the second gamma energy
  //fEnr2 : the first Q value
  //////////////////////////////////////
  //cout << fTime << endl;
  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout << fixed << setprecision(5);

  vector<Int_t> Run1;
  vector<Int_t> List1;
  vector<Int_t> Chan1;
  vector<Double_t> Enr1;
  vector<Double_t> DCR1;
  vector<Double_t> AvsE1;
  vector<Double_t> TrapETailMin1;
  vector<Double_t> Time1;
  for(size_t i=0;i<fEntries;i++){
    fSkimTree->GetEntry(i);
    Int_t run1 = fRun;
    for(size_t j=0;j<fChannel->size();j++){	      
      Int_t chan1 = fChannel->at(j);
      Double_t enr1 = fTrapENFCal->at(j);
      Double_t dcr1 = fDCR->at(j);
      Double_t avse1 = fAvsE->at(j);
      Double_t trapetailmin = fTrapETailMin->at(j);
      Double_t time1 = fglobalTime->AsDouble()+ ftOffset->at(j)/1e9;
      if( (enr1>fEnr1 && enr1<10000 && chan1%2==0) 
	  && 
	  !(0.015<dcr1 && dcr1<0.018 && 49<enr1 && enr1<51) && 
	  !(0.004<dcr1 && dcr1<0.006 && 49<enr1 && enr1<52) &&
	  !(-0.006<dcr1 && dcr1<-0.004 && 48<enr1 && enr1<58) &&
	  !(avse1>500 && (chan1 == 691 || chan1 == 1124)) &&
	  !(trapetailmin >2 && trapetailmin <3.5 && enr1>191 && enr1<192)
	 ){ 
	Run1.push_back(run1);
	List1.push_back(i);
	Chan1.push_back(chan1);
	Enr1.push_back(enr1);
	DCR1.push_back(dcr1);
	AvsE1.push_back(avse1);
	TrapETailMin1.push_back(trapetailmin);
	Time1.push_back(time1);
      }
    }
  }

  vector<Int_t> Run2;
  vector<Int_t> Run3;
  vector<Int_t> List2;
  vector<Int_t> List3;
  vector<Int_t> Chan2;
  vector<Int_t> Chan3;
  vector<Double_t> Enr2;
  vector<Double_t> Enr3;
  vector<Double_t> Time2;
  vector<Double_t> Time3;
  vector<Int_t> Event2;
  vector<Int_t> Event3;
  vector<Double_t> Mu_s2;
  vector<Double_t> Mu_s3;
  vector<Double_t> DCR2;
  vector<Double_t> DCR3;
  vector<Double_t> AvsE2;
  vector<Double_t> AvsE3;
  vector<Double_t> TrapETailMin2;
  vector<Double_t> TrapETailMin3;
  cout.precision(15);
  if((Int_t)List1.size()>0){
    for(size_t i=0;i<List1.size();i++){
      fSkimTree->GetEntry(List1.at(i));
      Int_t run1 = Run1.at(i);
      //Double_t time1 = flocalTime_s;
      Double_t time1 = Time1.at(i);
      Double_t time2 = time1;
      Double_t dtmu_s1 = fdtmu_s->at(0);

      Int_t event1 = fEvent;
      Int_t ii = List1.at(i)-1;
      fSkimTree->GetEntry(ii);
      Int_t run2 = fRun;
      time2 = fglobalTime->AsDouble();
      Int_t event2 = fEvent;
      Double_t dcr2 = DCR1.at(i);
      Double_t avse2 = AvsE1.at(i);
      Double_t trapetailmin2 = TrapETailMin1.at(i);
      //cout<< run1 << " "<< event1 << " "<< Enr1.at(i) << endl;
      while( (time1-time2)<fTime*20 && (time1-time2) >=0 && ii>0){
	for(size_t j=0;j<fChannel->size();j++){
	  Int_t chan2 = fChannel->at(j);
	  Double_t enr2 = fTrapENFCal->at(j);	    
	  Double_t dtmu_s2 = fdtmu_s->at(j);
	  Double_t dt = ftOffset->at(j)/1e9;
	  Double_t dcr3 = fDCR->at(j);
	  Double_t avse3 = fAvsE->at(j);
	  Double_t trapetailmin3 = fTrapETailMin->at(j);
	  //cout << run1 << " "<< event1 << " " << time1 << " "<< event2 << " "<< time2 << " " << dt << " " << chan2 << endl;
	  if( (enr2 > fEnr2 && enr2<10000 && chan2%2==0 && (time1 - (time2+dt))>0) 
	      &&
	      !(0.015<dcr3 && dcr3<0.018 && 49<enr2 && enr2<51) &&
	      !(0.004<dcr3 && dcr3<0.006 && 49<enr2 && enr2<52) &&
	      !(-0.006<dcr3 && dcr3<-0.004 && 48<enr2 && enr2<58) &&
	      !(avse3>500 && (chan2 == 691 || chan2 == 1124)) &&
	      !(trapetailmin3 >2 && trapetailmin3 <3.5 && enr2>191 && enr2<192)
	     ){
	    Run2.push_back(run1);
	    List2.push_back(List1.at(i));
	    Chan2.push_back(Chan1.at(i));
	    Enr2.push_back(Enr1.at(i));
	    Time2.push_back(time1);
	    Event2.push_back(event1);
	    Mu_s2.push_back(dtmu_s1);
	    DCR2.push_back(dcr2);
	    AvsE2.push_back(avse2);
	    TrapETailMin2.push_back(trapetailmin2);
	    Run3.push_back(run2);
	    List3.push_back(ii);
	    Chan3.push_back(chan2);
	    Enr3.push_back(enr2);			     
	    Time3.push_back(time2+dt);
            Event3.push_back(event2);
            Mu_s3.push_back(dtmu_s2);
	    DCR3.push_back(dcr3);
	    AvsE3.push_back(avse3);
	    TrapETailMin3.push_back(trapetailmin3);
	    //cout << time1 << " " << time2 << " "<< time1-time2 << endl;
	  }
	}     
	ii--;
	fSkimTree->GetEntry(ii);
        run2 = fRun;
        time2 = fglobalTime->AsDouble();
        event2 = fEvent;
      }
    }
  }


  for(size_t i=0;i<List2.size();i++){ 


    fout << fDataSet << " " << fSubSet << " "
	 << List2.at(i) << " " << Run2.at(i) << " " << Event2.at(i) << " " << Time2.at(i) << " " 
	 << Chan2.at(i) << " " << Enr2.at(i) << " "
	 << DCR2.at(i) << " "  << AvsE2.at(i) << " " << TrapETailMin2.at(i) << " " << Mu_s2.at(i) << " "

	 << List3.at(i) << " " << Run3.at(i) << " " << Event3.at(i) << " " << Time3.at(i) << " " 
	 << Chan3.at(i) << " " << Enr3.at(i) << " " 
	 << DCR3.at(i) << " "  << AvsE3.at(i) << " " << TrapETailMin3.at(i) << " " << Mu_s3.at(i) << " "
	 << Time2.at(i) - Time3.at(i) << endl;
  }
}



void MJDSkim::SearchDelayedMultiplicityEvent(Double_t fEnr1, Double_t fEnr2, Double_t fTime, string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the second gamma energy
  //fEnr2 : the first Q value
  //////////////////////////////////////
  //cout << fTime << endl;
  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout << fixed << setprecision(5);

  vector<Int_t> Run1;
  vector<Int_t> List1;
  vector<Int_t> Chan1;
  vector<Double_t> Enr1;
  vector<Double_t> DCR1;
  vector<Double_t> AvsE1;
  vector<Double_t> TrapETailMin1;
  vector<Double_t> Time1;
  for(size_t i=0;i<fEntries;i++){
    fSkimTree->GetEntry(i);
    Int_t run1 = fRun;
    Int_t mHL = fmH;
    for(size_t j=0;j<fChannel->size();j++){	      
      Int_t chan1 = fChannel->at(j);
      Double_t enr1 = fTrapENFCal->at(j);
      Double_t dcr1 = fDCR->at(j);
      Double_t avse1 = fAvsE->at(j);
      Double_t trapetailmin = fTrapETailMin->at(j);
      Double_t time1 = fglobalTime->AsDouble()+ ftOffset->at(j)/1e9;
      if( (abs(enr1-fEnr1)<10 && chan1%2==0 && mHL>1) 
	  && 
	  !(0.015<dcr1 && dcr1<0.018 && 49<enr1 && enr1<51) && 
	  !(0.004<dcr1 && dcr1<0.006 && 49<enr1 && enr1<52) &&
	  !(-0.006<dcr1 && dcr1<-0.004 && 48<enr1 && enr1<58) &&
	  !(avse1>500 && (chan1 == 691 || chan1 == 1124)) &&
	  !(trapetailmin >2 && trapetailmin <3.5 && enr1>191 && enr1<192)
	 ){ 
	Run1.push_back(run1);
	List1.push_back(i);
	Chan1.push_back(chan1);
	Enr1.push_back(enr1);
	DCR1.push_back(dcr1);
	AvsE1.push_back(avse1);
	TrapETailMin1.push_back(trapetailmin);
	Time1.push_back(time1);
      }
    }
  }

  vector<Int_t> Run2;
  vector<Int_t> Run3;
  vector<Int_t> List2;
  vector<Int_t> List3;
  vector<Int_t> Chan2;
  vector<Int_t> Chan3;
  vector<Double_t> Enr2;
  vector<Double_t> Enr3;
  vector<Double_t> Time2;
  vector<Double_t> Time3;
  vector<Int_t> Event2;
  vector<Int_t> Event3;
  vector<Double_t> Mu_s2;
  vector<Double_t> Mu_s3;
  vector<Double_t> DCR2;
  vector<Double_t> DCR3;
  vector<Double_t> AvsE2;
  vector<Double_t> AvsE3;
  vector<Double_t> TrapETailMin2;
  vector<Double_t> TrapETailMin3;
  cout.precision(15);
  if((Int_t)List1.size()>0){
    for(size_t i=0;i<List1.size();i++){
      fSkimTree->GetEntry(List1.at(i));
      Int_t run1 = Run1.at(i);
      //Double_t time1 = flocalTime_s;
      Double_t time1 = Time1.at(i);
      Double_t time2 = time1;
      Double_t dtmu_s1 = fdtmu_s->at(0);

      Int_t event1 = fEvent;
      Int_t ii = List1.at(i)-1;
      fSkimTree->GetEntry(ii);
      Int_t run2 = fRun;
      time2 = fglobalTime->AsDouble();
      Int_t event2 = fEvent;
      Double_t dcr2 = DCR1.at(i);
      Double_t avse2 = AvsE1.at(i);
      Double_t trapetailmin2 = TrapETailMin1.at(i);
      //cout<< run1 << " "<< event1 << " "<< Enr1.at(i) << endl;
      while( (time1-time2)<fTime*20 && (time1-time2) >=0 && ii>0){
	for(size_t j=0;j<fChannel->size();j++){
	  Int_t chan2 = fChannel->at(j);
	  Double_t enr2 = fTrapENFCal->at(j);	    
	  Double_t dtmu_s2 = fdtmu_s->at(j);
	  Double_t dt = ftOffset->at(j)/1e9;
	  Double_t dcr3 = fDCR->at(j);
	  Double_t avse3 = fAvsE->at(j);
	  Double_t trapetailmin3 = fTrapETailMin->at(j);
	  //cout << run1 << " "<< event1 << " " << time1 << " "<< event2 << " "<< time2 << " " << dt << " " << chan2 << endl;
	  if( (enr2 < fEnr2 && enr2>5 && chan2%2==0 && (time1 - (time2+dt))>0) 
	      &&
	      !(0.015<dcr3 && dcr3<0.018 && 49<enr2 && enr2<51) &&
	      !(0.004<dcr3 && dcr3<0.006 && 49<enr2 && enr2<52) &&
	      !(-0.006<dcr3 && dcr3<-0.004 && 48<enr2 && enr2<58) &&
	      !(avse3>500 && (chan2 == 691 || chan2 == 1124)) &&
	      !(trapetailmin3 >2 && trapetailmin3 <3.5 && enr2>191 && enr2<192)
	     ){
	    Run2.push_back(run1);
	    List2.push_back(List1.at(i));
	    Chan2.push_back(Chan1.at(i));
	    Enr2.push_back(Enr1.at(i));
	    Time2.push_back(time1);
	    Event2.push_back(event1);
	    Mu_s2.push_back(dtmu_s1);
	    DCR2.push_back(dcr2);
	    AvsE2.push_back(avse2);
	    TrapETailMin2.push_back(trapetailmin2);
	    Run3.push_back(run2);
	    List3.push_back(ii);
	    Chan3.push_back(chan2);
	    Enr3.push_back(enr2);			     
	    Time3.push_back(time2+dt);
            Event3.push_back(event2);
            Mu_s3.push_back(dtmu_s2);
	    DCR3.push_back(dcr3);
	    AvsE3.push_back(avse3);
	    TrapETailMin3.push_back(trapetailmin3);
	    //cout << time1 << " " << time2 << " "<< time1-time2 << endl;
	  }
	}     
	ii--;
	fSkimTree->GetEntry(ii);
        run2 = fRun;
        time2 = fglobalTime->AsDouble();
        event2 = fEvent;
      }
    }
  }


  for(size_t i=0;i<List2.size();i++){ 


    fout << fDataSet << " " << fSubSet << " "
	 << List2.at(i) << " " << Run2.at(i) << " " << Event2.at(i) << " " << Time2.at(i) << " " 
	 << Chan2.at(i) << " " << Enr2.at(i) << " "
	 << DCR2.at(i) << " "  << AvsE2.at(i) << " " << TrapETailMin2.at(i) << " " << Mu_s2.at(i) << " "

	 << List3.at(i) << " " << Run3.at(i) << " " << Event3.at(i) << " " << Time3.at(i) << " " 
	 << Chan3.at(i) << " " << Enr3.at(i) << " " 
	 << DCR3.at(i) << " "  << AvsE3.at(i) << " " << TrapETailMin3.at(i) << " " << Mu_s3.at(i) << " "
	 << Time2.at(i) - Time3.at(i) << endl;
  }
}


void MJDSkim::SearchMuonCoinEvent(Double_t fEnr, Double_t fTime, string fOutputFile){
  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout.precision(15);
  vector<Int_t> Run1;
  vector<Int_t> Event1;
  vector<Int_t> Chan1;
  vector<Double_t> Enr1;
  vector<Double_t> DtMu1;
  for(size_t i=0;i<fEntries;i++){
    fSkimTree->GetEntry(i);
    Int_t run1 = fRun;
    Int_t event1 = fEvent;
    for(size_t j=0;j<fChannel->size();j++){	      
      Int_t chan1 = fChannel->at(j);
      Double_t enr1 = fTrapENFCal->at(j);
      Double_t dtmu1 = fdtmu_s->at(j);
      if(abs(enr1-fEnr)<20 && chan1%2==0 && fIsGood->at(j) == 1 && dtmu1<fTime*20){
	Run1.push_back(run1);
	Event1.push_back(event1);
	Chan1.push_back(chan1);
	Enr1.push_back(enr1);
	DtMu1.push_back(dtmu1);
      }
    }
  }
  for(size_t i=0;i<Run1.size();i++){
    fout << Run1.at(i) << " " << Event1.at(i)<< " " << Chan1.at(i) << " " << Enr1.at(i) << " " << DtMu1.at(i) << endl;
  }
}


void MJDSkim::SearchEnergyEvent(Double_t fEnr, Double_t fEnrWindow, string fOutputFile){
  //////////////////////////////////////
  //fEnr1 : the gamma energy
  //////////////////////////////////////
  ofstream fout(Form("%s",fOutputFile.c_str()));
  fout.precision(15);
  for(size_t i=0;i<fEntries;i++){
    fSkimTree->GetEntry(i);
    Int_t irun = fRun;
    Int_t ievent = fEvent;
    for(size_t j=0;j<fChannel->size();j++){
      Int_t chan1 = fChannel->at(j);
      Double_t enr1 = fTrapENFCal->at(j);
      Double_t dcr1 = fDCR->at(j);
      Double_t avse1 = fAvsE->at(j);
      Double_t trapetailmin = fTrapETailMin->at(j);
      Double_t time1 = fglobalTime->AsDouble()+ ftOffset->at(j)/1e9;
      Int_t ic = fC->at(j);
      Int_t ip = fP-> at(j);
      Int_t id = fD->at(j);
      string pos = Form("%d%d%d",ic,ip,id);
      if((abs(enr1-fEnr)<fEnrWindow && chan1%2==0) 
	  && 
	  !(0.015<dcr1 && dcr1<0.018 && 49<enr1 && enr1<51) && 
	  !(0.004<dcr1 && dcr1<0.006 && 49<enr1 && enr1<52) &&
	  !(-0.006<dcr1 && dcr1<-0.004 && 48<enr1 && enr1<58) &&
	  !(avse1>500 && (chan1 == 691 || chan1 == 1124)) &&
	  !(trapetailmin >2 && trapetailmin <3.5 && enr1>191 && enr1<192)
	 ){ 
	fout << irun << " "  << ievent << " " << time1 << " " << pos.c_str() << " "
	     << chan1 << " " << enr1 << " "
	     << dcr1 << " "<< avse1 << " " << trapetailmin << endl;	    
      }      
    }
  }
  cout << "SearchEnergy file is generated...." <<endl;
}


TH1D* MJDSkim::GetWaveform(Int_t fR,Int_t fEntry, Int_t fChan,Double_t fEnr){
  GATDataSet ds(fR);
  Int_t Entry = fEntry+1;
  TCut cut1 = Form("channel == %d", fChan);
  GATWaveformBrowser wb;
  wb.LoadWaveforms(ds.GetBuiltChain(), cut1, "fWaveforms", Entry);
  size_t i=wb.GetNWaveforms();
  if(i>0){
    //cout << i << endl;
    MGTWaveform *w1 = wb.GetWaveform(i-1).get();
    TH1D* h= (TH1D*)w1->GimmeHist();
    h->SetTitle(Form("Run=%d,Entry=%d, Channel=%d, Energy=%f(keV);Time(ns);ADC",fR,fEntry,fChan,fEnr));
    h->SetName(Form("waveform_%d_%d_%d",fR,fEntry,fChan));
    return h;
  }else{
    TH1D* h = NULL;
    return h;
    cout << "No this waveform!" << endl;
  }
}

TH1D* MJDSkim::TrapezoidalFilter(TH1D* hist, double RampTime, double FlatTime , double DecayTime )
{

  TH1D *h = (TH1D*)hist->Clone();
  h->Reset(0);
  Int_t entries = hist->GetEntries();
  vector<Double_t> anInput;
  for(Int_t i1 = 0;i1<entries;i1++){
    //cout << i1 << " "<< hist->GetBinContent(i1) << endl;
    anInput.push_back(hist->GetBinContent(i1));
  }

  // Arbitrarily choosing 1000 as decay constant
  double decayConstant = 0.0;
  if(DecayTime != 0) decayConstant = 1./(exp(1./DecayTime) - 1);
  // double decayConstant = 0;
  double rampStep = RampTime;
  double flatStep = FlatTime;
  double baseline = 0; // No baseline for now
  double norm = rampStep;
  if(decayConstant != 0)norm *= decayConstant;

  std::vector<double> fVector;
  std::vector<double> anOutput;
  if(fVector.size() != anInput.size()) 
    {
      fVector.resize(anInput.size());
      anOutput.resize(anInput.size());
    }

  fVector[0] = anInput[0] - baseline;
  anOutput[0] = (decayConstant+1.)*(anInput[0] - baseline);
  double scratch = 0.0;
  for(size_t i = 1; i < anInput.size(); i++)
    {
      // This is a little tricky with all the ternary operators, but it's faster
      // this way.  We must check the bounds.
      scratch = anInput[i]  - ((i>=rampStep) ? anInput[i-rampStep] : baseline)
	- ((i>=flatStep+rampStep) ? anInput[i-flatStep-rampStep] : baseline)
	+ ((i>=flatStep+2*rampStep) ? anInput[i-flatStep-2*rampStep] : baseline);  
      if(decayConstant != 0.0) 
	{
	  fVector[i] = fVector[i-1] + scratch; 
	  anOutput[i] = (anOutput[i-1] + fVector[i] + decayConstant*scratch);
	  // anOutput[i] /= norm;
	} 
      else anOutput[i] = anOutput[i-1] + scratch;
    }
  for(size_t i = 2*rampStep+flatStep; i < anInput.size(); i++) anOutput[i-(2*rampStep+flatStep)] = anOutput[i];
  anOutput.resize(anOutput.size()-(2*rampStep+flatStep));
  // Rescale event by normalization factor
  for(size_t i = 1; i < anOutput.size(); i++)anOutput[i] = anOutput[i]/norm;
  //cout << anOutput.size() << " "<< entries << endl;
  for(size_t i2 = 0;i2<anOutput.size();i2++){
    h->SetBinContent(i2,anOutput.at(i2));
  }
  return h;
  //return anOutput;
}

TH1D* MJDSkim::GetHistoSmooth(TH1D* hist, Int_t DeltaBin){
  TH1D *h = (TH1D*)hist->Clone();
  h->Reset(0);
  Int_t entries = hist->GetEntries();
  Double_t biny=0;
  Double_t biny1=0;
  Double_t biny2=0;
  Double_t ave1=0;
  Double_t ave2=0;
  for(Int_t i1 = DeltaBin;i1<entries-DeltaBin;i1++){
    Double_t fDummySum=0;
    Double_t fDummyAve=0;
    Int_t count = 0;
    for(Int_t i2 = i1-DeltaBin;i2<i1+DeltaBin;i2++){
      biny = hist->GetBinContent(i2);
      biny1 = hist->GetBinContent(i2-1);
      biny2 = hist->GetBinContent(i2+1);
      
      if(abs(biny-biny1)>10 || abs(biny-biny2)>10){
	fDummySum = fDummySum;	
      }else{
	fDummySum = fDummySum+biny;
	count++;
      }
    }
    fDummyAve = fDummySum/(count);
    h->SetBinContent(i1,fDummyAve);
    if(i1==DeltaBin){
      ave1 = fDummyAve;
    }
    if(i1==entries-DeltaBin-1){
      ave2 = fDummyAve;
    }
  }
  for(Int_t i = 0;i<DeltaBin;i++){
    h->SetBinContent(i,ave1);
  }
  for(Int_t i = entries-DeltaBin;i<entries;i++){
    h->SetBinContent(i,ave2);
  }
  h->Draw("wf.pdf");
  return h;
}


TH1D* MJDSkim::GetHistoDerivative(TH1D* hist, Int_t DeltaBin){
  TH1D *h = (TH1D*)hist->Clone();
  h->Reset(0);
  Int_t entries = hist->GetEntries();
  Double_t binx= hist->GetBinWidth(1);
  Double_t biny1=0;
  Double_t biny2=0;
  Double_t slope=0;
  Double_t slope1=0;
  Double_t slope2=0;
  for(Int_t i = DeltaBin;i<entries-DeltaBin;i++){
    biny1 = hist->GetBinContent(i-DeltaBin);
    biny2 = hist->GetBinContent(i+DeltaBin);
    slope = (biny2-biny1)/(binx*2*DeltaBin);
    h->SetBinContent(i,slope);
    if(i==DeltaBin){
      slope1 = slope;
    }
    if(i==entries-DeltaBin-1){
      slope2 = slope;
    }
  }
  for(Int_t i = 0;i<DeltaBin;i++){
    h->SetBinContent(i,slope1);
  }
  for(Int_t i = entries-DeltaBin;i<entries;i++){
    h->SetBinContent(i,slope2);
  }
  return h;
}


TH1* MJDSkim::GetHistoFFT(TH1D* hist){
  TH1 *h = NULL;
  TVirtualFFT::SetTransform(0);
  h = hist->FFT(h,"MAG");
  return h;
}


vector<Int_t> MJDSkim::Sort(vector<Double_t> X){
  Int_t n = X.size();
  const Int_t n1 = n;
  Double_t x1[n1];
  Int_t index[n1];
  for(Int_t i =0;i<n;i++){
    x1[i] = X.at(i);
  }
  TMath::Sort(n,x1,index,0);
  vector<Int_t> Index;
  for(Int_t i=0;i<n;i++){
    Index.push_back(index[i]);
  }
  return Index;
}

vector<Int_t> MJDSkim::Clean(vector<Double_t> X){

  vector<Int_t> Index;
  Index.push_back(0);
  for(size_t i = 1 ; i< X.size();i++){
    Int_t count = 0;
    for(size_t j=0;j<i;j++){
      if(X.at(i) == X.at(j)){
	count++;
      }
    }      
    if(count==0){
      Index.push_back(i);
    }
  }
  return Index;
}

Double_t MJDSkim::GetYValue(TH1D* hist, Double_t X){
  Double_t binwidth= hist->GetBinWidth(1);
  Double_t binx=0;
  Double_t biny=0;
  Double_t Y=0;
  for(size_t i = 0;i<hist->GetEntries();i++){
    binx = hist->GetBinCenter(i);
    biny = hist->GetBinContent(i);
    if(abs(X-binx)<binwidth){
      Y = biny;
    }
  }
  return Y;
}
Int_t MJDSkim::FindPeaks(TH1D* hist, Double_t Low, Double_t Up, Double_t Resolution, Double_t Sigma, Double_t Threshold, vector<Double_t>* fPositionX, vector<Double_t>* fPositionY){
  Double_t baseline = 0;
  Double_t blrms = 0;
  Int_t count =0;
  for(size_t i = 0; i< hist->GetEntries(); i++){
    if(hist->GetBinCenter(i)>1000 && hist->GetBinCenter(i)<4000){
      baseline = baseline+hist->GetBinContent(i);
      count++;
    }
  }
  baseline = baseline/(Double_t)count;
  blrms = 0;
  count = 0;
  for(size_t i = 0; i< hist->GetEntries(); i++){
    if(hist->GetBinCenter(i)>1000 && hist->GetBinCenter(i)<4000){
      blrms = blrms + pow((hist->GetBinContent(i)-baseline),2);
      count++;
    }
  }
  blrms = sqrt(blrms/(Double_t)count);
  Double_t Ymin = baseline+2*blrms;
  Double_t fResolution = Resolution;
  Double_t fSigma = Sigma;
  Double_t fThreshold = Threshold;
  TSpectrum *s = new TSpectrum(100,fResolution);
  Int_t npeaks = s->Search(hist,fSigma,"new",fThreshold);

  Double_t *xpeaks = s->GetPositionX();
  Double_t *ypeaks = s->GetPositionY();

  for(Int_t i=0;i<npeaks;i++){
    if(xpeaks[i]>Low && xpeaks[i]<Up && ypeaks[i]>Ymin){
      hist->GetXaxis()->SetRangeUser(xpeaks[i]-20,xpeaks[i]+20);
      double ymax = hist->GetMaximum();
      int xbin = hist->GetMaximumBin();
      double xmax = hist->GetBinCenter(xbin);
      fPositionX->push_back(xmax);
      fPositionY->push_back(ymax);
    }
  }
  npeaks = fPositionX->size();
  delete s;
  return npeaks;
}

Double_t MJDSkim::GetMax(TH1* hist, Double_t Low, Double_t Up){
  hist->GetXaxis()->SetRangeUser(Low,Up);
  Double_t xmax = hist->GetMaximum();
  return xmax;
}


void MJDSkim::IsPileUpTag(Int_t fList, vector<Int_t>* IsPileUp, vector<Double_t>* Ratio, vector<Double_t>* DeltaT){

  Double_t Xmin = 4000;
  Double_t Xmax = 19000;
  if((fRun >=14503 && fRun<=15892) || fRun>=25675){
    Xmax = 38000;
  }
  Double_t fResolution = 1;
  Double_t fThreshold = 0.01;
  Double_t fSigma = 5;

  vector<Double_t> xp;
  vector<Double_t> yp;
  vector<Double_t> xp1;
  vector<Double_t> yp1;
  vector<Double_t> xp2;
  vector<Double_t> yp2;

  fSkimTree->GetEntry(fList);
  
  ofstream fout(Form("pileup_%d.txt",fRun),ios::app);
  for(size_t j=0;j<fChannel->size();j++){
    xp.clear();
    yp.clear();
    xp1.clear();
    yp1.clear();
    xp2.clear();
    yp2.clear();
    
    Int_t chan1 = fChannel->at(j);
    Double_t enr1 = fTrapENFCal->at(j);
    Double_t dcr1 = fDCR->at(j);
    Double_t avse1 = fAvsE->at(j);
    Double_t trapetailmin = fTrapETailMin->at(j);
    Double_t time1 = fglobalTime->AsDouble()+ ftOffset->at(j)/1e9;
    
    if( enr1>40 && enr1<12000 && 
       !(0.015<dcr1 && dcr1<0.018 && 49<enr1 && enr1<51) &&
	!(0.004<dcr1 && dcr1<0.006 && 49<enr1 && enr1<52) &&
	!(-0.006<dcr1 && dcr1<-0.004 && 48<enr1 && enr1<58) &&
	!(avse1>500 && (chan1 == 691 || chan1 == 1124)) &&
	!(trapetailmin >2 && trapetailmin <3.5 && enr1>191 && enr1<192)
	){
      TH1D* h = MJDSkim::GetWaveform(fRun,fEvent,chan1,enr1);      
      Double_t maxY = h->GetMaximum();
      if(maxY<10000){
	TH1D* h2 = MJDSkim::GetHistoDerivative(h,10);
	TH1D* h4 = MJDSkim::GetHistoDerivative(h2,10);
	
	Int_t maxBin0 = h->GetMaximumBin();
	Double_t maxX0 = h->GetBinCenter(maxBin0);    
	Int_t nPeak1 = MJDSkim::FindPeaks(h2, Xmin, Xmax,fResolution,fSigma, fThreshold, &xp, &yp);

	if(nPeak1>1){
	  // Check if the peak is zero (<1e-3) in the 2nd derviative waveform
	  // Scan the position +- 20 bin
	  for(Int_t ip = 0;ip<(Int_t)xp.size();ip++){
	    cout << "1st x : " << xp.at(ip) << " 1st y: " << yp.at(ip) << endl;
	    Double_t yabsmin =  10;
	    for(Int_t is = 0;is<40;is++){
	      Double_t y = abs(MJDSkim::GetYValue(h4, xp.at(ip)-20+is));
	      if(y<yabsmin){
		yabsmin = y;
	      }
	    }
	    cout << "2nd y: " << yabsmin << " " << maxX0 << endl;
	    if(abs(yabsmin)< 1e-2 && xp.at(ip)< maxX0){
	      xp1.push_back(xp.at(ip));
	      yp1.push_back(yp.at(ip));
	    }
	  }
	  
	  if(xp1.size()>1){
	    //Pick the maximum two peaks
	    vector<Int_t> Index2 = MJDSkim::Sort(yp1);
	    for(size_t ip1 = Index2.size()-2;ip1<Index2.size();ip1++){
	      Int_t ii = Index2.at(ip1);
	      xp2.push_back(xp1.at(ii));
	      yp2.push_back(yp1.at(ii));
	    }
	    //Calculate the ratio and delta T
	    Double_t ratio = yp2.at(1)/yp2.at(0);
	    Double_t deltaT = xp2.at(0)-xp2.at(1);
	    cout << " ratio : " << ratio << " deltaT : "<<deltaT << endl;
	    IsPileUp->push_back(1);
	    Ratio->push_back(ratio);
	    DeltaT->push_back(deltaT);

	    fout << fDataSet << " " << fSubSet << " " << fList << " " 
	         << fRun << " " << fEvent << " " << time1 <<  " " << chan1 << " " << enr1 << " " 
		 << ratio << " " << deltaT << " " << avse1 << " " << dcr1 << " " << trapetailmin << endl;
	  }else{
	    IsPileUp->push_back(0);
	    Ratio->push_back(0);
	    DeltaT->push_back(0);
	  }
	}else{
	  IsPileUp->push_back(0);
	  Ratio->push_back(0);
	  DeltaT->push_back(0);
	}
	delete h2;
	delete h4;
      }else{
	IsPileUp->push_back(0);
	Ratio->push_back(0);
	DeltaT->push_back(0);
      }
      delete h;
    }else{
      IsPileUp->push_back(0);
      Ratio->push_back(0);
      DeltaT->push_back(0);
    }
  }
}


void MJDSkim::PileUpTree(const char* pathName){
  TFile *newfile = new TFile(Form("%spileup_%d_%d.root", pathName, fDataSet,fSubSet), "recreate");
  TTree *newtree = new TTree("pileupTree", "pile-up tag");

  vector<Int_t>* IsPileUp = NULL;
  vector<Double_t>* PileUpRatio = NULL;
  vector<Double_t>* PileUpDeltaT =NULL;
  vector<Double_t>* AvsE = NULL;
  vector<Double_t>* DCR = NULL;
  vector<Double_t>* TrapETailMin = NULL;

  vector<Int_t> IsPileUp1;
  vector<Double_t> PileUpRatio1;
  vector<Double_t> PileUpDeltaT1;

  newtree->Branch("IsPileUp",&IsPileUp);
  newtree->Branch("PileUpRatio", &PileUpRatio);
  newtree->Branch("PileUpDeltaT", &PileUpDeltaT);
  newtree->Branch("AvsE", &AvsE);
  newtree->Branch("DCR", &DCR);
  newtree->Branch("TrapETailMin",&TrapETailMin);
  Int_t Run;
  Int_t Event;
  vector<Int_t>* Channel = NULL;
  vector<Int_t>* C = NULL;
  vector<Int_t>* P = NULL;
  vector<Int_t>* D = NULL;
  vector<bool>* IsEnr = NULL;
  vector<bool>* IsNat = NULL;
  vector<bool>* IsGood = NULL;
  vector<Double_t>* TOffset = NULL;
  vector<Double_t>* Dtmu_s = NULL;
  vector<Double_t>* Energy = NULL;

  Double_t LocalTime;
  Double_t GlobalTime;
  Int_t mHClean;
  bool MuVeto;
  Double_t MuTUnc;

  newtree->Branch("Run", &Run);
  newtree->Branch("Event", &Event);
  newtree->Branch("Channel",&Channel);
  newtree->Branch("Energy",&Energy);
  newtree->Branch("P",&P);
  newtree->Branch("D",&D);
  newtree->Branch("C",&C);
  newtree->Branch("IsEnr",&IsEnr);
  newtree->Branch("IsNat",&IsNat);
  newtree->Branch("IsGood",&IsGood);
  newtree->Branch("LocalTime",&LocalTime);
  newtree->Branch("GlobalTime",&GlobalTime);
  newtree->Branch("TOffset",&TOffset);
  newtree->Branch("mHClean",&mHClean);
  newtree->Branch("MuVeto",&MuVeto);
  newtree->Branch("Dtmu_s",&Dtmu_s);
  newtree->Branch("MuTUnc",&MuTUnc);


  Int_t entries = 0;
  if(fIsCal == 0){
    entries = fSkimTree->GetEntries();
  }else{
    entries = 1000;
  }
  for(Int_t i=0;i<entries;i++){
    fSkimTree->GetEntry(i);
    //    cout << fRun << " " << fEvent << endl;
    if(fRun==21892 && fEvent == 59434){
      cout << fRun << " " << fEvent << endl;

    Double_t count = 0;


    IsPileUp1.clear();
    PileUpRatio1.clear();
    PileUpDeltaT1.clear();

    IsPileUp->clear();
    PileUpRatio->clear();
    PileUpDeltaT->clear();

    MJDSkim::IsPileUpTag(i,&IsPileUp1,&PileUpRatio1,&PileUpDeltaT1);
    
    Channel->clear();
    Energy->clear();
    P->clear();
    D->clear();
    C->clear();
    IsEnr->clear();
    IsNat->clear();
    IsGood->clear();
    TOffset->clear();
    Dtmu_s->clear();
    AvsE->clear();
    DCR->clear();
    TrapETailMin->clear();
    
    Run = fRun;
    Event = fEvent;
    LocalTime = flocalTime_s;
    GlobalTime = fglobalTime->AsDouble();
    mHClean = fmH;
    MuVeto = fmuVeto;
    MuTUnc = fmuTUnc;
    
    //cout << IsPileUp1.size() << " "<< fChannel->size() << endl;
    for(size_t j = 0;j<fChannel->size();j++){
      if(IsPileUp1.at(j)>0){
	
	IsPileUp->push_back(IsPileUp1.at(j));
	PileUpRatio->push_back(PileUpRatio1.at(j));
	PileUpDeltaT->push_back(PileUpDeltaT1.at(j));
		
	Channel->push_back(fChannel->at(j));
	Energy->push_back(fTrapENFCal->at(j));
	P->push_back(fP->at(j));
	D->push_back(fD->at(j));
	C->push_back(fC->at(j));
	IsEnr->push_back(fIsEnr->at(j));
	IsNat->push_back(fIsNat->at(j));
	IsGood->push_back(fIsGood->at(j));
	TOffset->push_back(ftOffset->at(j)/1e9);
	Dtmu_s->push_back(fdtmu_s->at(j));	       
	AvsE->push_back(fAvsE->at(j));
	DCR->push_back(fDCR->at(j));
	TrapETailMin->push_back(fTrapETailMin->at(j));
      }
    }
    newtree->Fill();      
    }
  }
  newtree->Write();
  cout << "Pile-up file is generated...." <<endl;
  delete newfile;
}

void MJDSkim::TimeDiffTree(const char* pathName){

  TFile *newfile = new TFile(Form("%stimediff_%d_%d.root", pathName, fDataSet,fSubSet), "recreate");
  TTree *newtree = new TTree("timediffTree", "");
  vector<Double_t>* TimeDiff = NULL;
  newtree->Branch("TimeDiff",&TimeDiff);
  Int_t entries = 0;
  if(fIsCal == 0){
    entries = fSkimTree->GetEntries();
  }else{
    entries = 1000;
  }
  for(Int_t i=0;i<entries;i++){
    Double_t count = 0;
    fSkimTree->GetEntry(i);
    TimeDiff->clear();
    if(fmH>1){
      for(size_t j = 0;j<fChannel->size();j++){	
	vector<Double_t> tOff;
	tOff.clear();
	for(size_t j = 0;j<fChannel->size();j++){
	  if(fChannel->at(j)%2==0){
	    tOff.push_back(ftOffset->at(j)/1e9);
	  }
	}
	vector<Int_t> tInd = MJDSkim::Sort(tOff);	
	vector<Double_t> tDiff;
	tDiff.clear();   
	if(tOff.size()>0){
	  Double_t t0 = tOff.at(0);
	  for(size_t j=0;j<tOff.size();j++){
	    Double_t t1 = tOff.at(j);
	    TimeDiff->push_back(t1-t0);
	    t0 = t1;
	  }
	}
      }
    }
    newtree->Fill();      
  }
  newtree->Write();
  cout << "TimeDiff file is generated...." <<endl;
  delete newfile;
}

TH1D* MJDSkim::FillHisto(TChain* mTree, string InputParaName, string OutputParaName,string CutName, Int_t Bin, Double_t Low, Double_t Up){
  TH1D *h = new TH1D(Form("%s", OutputParaName.c_str()),"", Bin, Low, Up);
  TCut cut1 = Form("%s",CutName.c_str());
  //TCanvas *c1 = new TCanvas("c1");
  mTree->Draw(Form("%s>>%s",InputParaName.c_str(),OutputParaName.c_str()),cut1);
  return h;
  cout << "histogram is generated" << endl;
}
