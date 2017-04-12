#include "GATAutoCal.hh"
#include "MAPP3END.hh"
#include "MAPP3JDY.hh"
#include "MAPP3KJR.hh"
#include "MAPP3LQG.hh"
#include "MAPP3LQK.hh"
#include "MJTRun.hh"
#include "MJAnalysisDoc.hh"
#include "MJAnalysisParameters.hh"
#include "GATWaveformBrowser.hh"
#include "GATPeakShape.hh"
#include "GATMultiPeakFitter.hh"
#include "GATHybridMonteCarlo.hh"
#include "MJProvenance.hh"
#include "MJDatabase.hh"
#include "MJDBUtilities.hh"
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
#include <iostream>
#include <bitset>

//ClassImp(GATAutoCal)

using namespace std;
using namespace katrin;
using namespace MJDB;
//using namespace GATPeakShapeFunction;
///////////////////////////////////////////////////////////////
// Fitting Function:
Double_t Flat(Double_t *par){
  Double_t fitval = par[0];
  return fitval;
}
Double_t Slope(Double_t *v, Double_t *par){
  Double_t fitval = par[0]+par[1]*v[0];
  return fitval;
}
Double_t Gaus(Double_t *v, Double_t *par){
  Double_t arg = 0;
  arg = (v[0] - (par[1]))/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}

Double_t TwoGaus(Double_t *v, Double_t *par){
  return Gaus(v,par)+Gaus(v,&par[3])+Slope(v,&par[6]);
}

Double_t GausSlope(Double_t *v, Double_t *par){
  return Gaus(v,par)+Slope(v,&par[3]);
}

Double_t SkewGaus(Double_t *v, Double_t *par){
  Double_t arg = 0;
  if (v[0]>=(par[1])){
    arg = (v[0] - (par[1]))/par[2];
  }
  if (v[0]<(par[1])){
    arg = (v[0] - (par[1]))/par[3];
  }
  
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}

Double_t SkewGausSlope(Double_t *v, Double_t *par){
   return SkewGaus(v,par)+Slope(v,&par[4]);
}

Double_t Reso(Double_t *v, Double_t *par){
  Double_t fitval = TMath::Sqrt(par[0]*par[0]+par[1]*par[1]*v[0]+par[2]*par[2]*v[0]*v[0]);
  return fitval;
}

Double_t Quad(Double_t *v, Double_t *par){
  Double_t fitval = 0;
  fitval = par[0]+par[1]*v[0]+par[2]*pow(v[0],2);
  return fitval;
}

Double_t Ramp(Double_t *v, Double_t *par){
  Double_t fitval = 0;
  fitval = pow((par[0]/v[0]),2)+pow(par[1],2)+pow(par[2]*v[0],2);
  return fitval;
}

////////////////////////////////////////////////////////

//Set ChannelMap, mjdTree (from GATDataSet), CalibrationPeak, Initial Parameters
GATAutoCal::GATAutoCal(Int_t StartRun, Int_t EndRun) : fDS(StartRun,EndRun)
{
  fStartRun = StartRun;
  fEndRun = EndRun;
  if((fStartRun>0 && fEndRun<30000000) || (fStartRun > 60000000 && fEndRun<65000000)){
    fMap = fDS.GetChannelMap();
    fMjdTree = fDS.GetGatifiedChain(false);
  }

  SetCalibrationPeak();
  SetParameters();

}

Double_t GATAutoCal::IsNan(string Input){
  Double_t output = 0;
  if(Input == "nan" || Input == "-nan" || Input == "inf" || Input == "-inf"){
    output = 0;
  }else{
    output = atof(Input.c_str());
  }
  return output;
}

//Set external mjdTree
void GATAutoCal::SetMjdTree(TChain *mjdTree){
  fMjdTree = mjdTree;
}

//Set Initial Parameters: isGood, isEnr, Pulser, PulserCal 
//The values are stored in MAPP3*.hh
void GATAutoCal::SetParameters()
{
  //DS0
  if(fStartRun >= 2337 && fEndRun < 8183){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad0[i]);
      fEnriched.push_back(OrtecBeGe1[i]);
      if(fEndRun<3464) fPulserCal.push_back(PulserCal01[i]);
      else fPulserCal.push_back(PulserCal02[i]);
      if(fEndRun<3464) fPulser.push_back(Pulser01[i]);
      else fPulser.push_back(Pulser02[i]);
    }
  }

  //DS1
  if(fStartRun >= 8722 && fEndRun<14503){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad1[i]);
      fEnriched.push_back(OrtecBeGe1[i]);
      fPulserCal.push_back(PulserCal1[i]);
      fPulser.push_back(Pulser1[i]);
    }
  }

  //DS2
  if(fStartRun >= 14503 && fEndRun<16836){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad2[i]);
      fEnriched.push_back(OrtecBeGe1[i]);
      fPulserCal.push_back(PulserCal1[i]);
      fPulser.push_back(Pulser1[i]);
    }
  }
  //DS3
  if(fStartRun >= 16836 && fEndRun<18589){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad3[i]);
      fEnriched.push_back(OrtecBeGe1[i]);
      fPulserCal.push_back(PulserCal1[i]);
      fPulser.push_back(Pulser1[i]);
    }
  }

  //DS4
  if(fStartRun >= 60000001 && fEndRun<=60002394){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad4[i]);
      fEnriched.push_back(OrtecBeGe2[i]);
      fPulserCal.push_back(PulserCal4[i]);
      fPulser.push_back(Pulser4[i]);
    }
  }

  //DS5
  if(fStartRun >= 18590 && fEndRun<5000000){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad51[i]);
      fEnriched.push_back(OrtecBeGe1[i]);
      fPulserCal.push_back(PulserCal51[i]);
      fPulser.push_back(Pulser51[i]);
    }
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad52[i]);
      fEnriched.push_back(OrtecBeGe2[i]);
      fPulserCal.push_back(PulserCal52[i]);
      fPulser.push_back(Pulser52[i]);
    }
  }

  //STC
  if(fStartRun >= 30000000 && fEndRun<40000000){
    for(Int_t i = 0;i<8;i++){
      fGoodBad.push_back(1);
      fEnriched.push_back(0);
      fPulserCal.push_back(0);
      fPulser.push_back(0);
    }
  }

  //DSPM
  if(fStartRun >= 45000000 && fEndRun<50000000){
    for(Int_t i = 0;i<proChannels;i++){
      cout << i << " " << GoodBadpro[i] << endl;
      fGoodBad.push_back(GoodBadpro[i]);
      fEnriched.push_back(0);
      fPulserCal.push_back(PulserCalpro[i]);
      fPulser.push_back(Pulserpro[i]);
    }
  }
  //Set Channel Parameters
  SetChannel();
  //Set Data Set Parameters
  SetDataSet();
}

//Set Run Range
void GATAutoCal::SetRunRange(Int_t StartRun,Int_t EndRun){
  fStartRun = StartRun;
  fEndRun = EndRun;
}

//Set Calibration Peak : default 238,583,727,2614
void GATAutoCal::SetCalibrationPeak(){
  const Int_t nReference = 18;  
  Double_t CalibrationPeak[nReference] = {74.97,77.107,84.5,87.349,89.9,115.183,215.5,238.6,277.358,300.087, 510.77, 
  583.191, 727.330,860.56,1592.5, 1620, 2104,2614.5};  
  vector<Int_t> index;
  index.push_back(7); //238 
  index.push_back(11);//583
  index.push_back(12);//727
  index.push_back(17);//2614
  Int_t nIndex = index.size();
  for(Int_t i = 0;i<nIndex;i++){
    fCalibrationPeak.push_back(CalibrationPeak[index.at(i)]);
  }
  fCalibrationPeaks = fCalibrationPeak.size();
}

//Set external gamma peak value
void GATAutoCal::SetCalibrationPeak(vector<Double_t> gammakeV){
  fCalibrationPeak.erase(fCalibrationPeak.begin(),fCalibrationPeak.end());
  Int_t nIndex = gammakeV.size();
  for(Int_t i = 0;i<nIndex;i++){
    fCalibrationPeak.push_back(gammakeV.at(i));
  }
  fCalibrationPeaks = fCalibrationPeak.size();
}

//Set data set name
void GATAutoCal::SetDataSet()
{
  if(fStartRun>0 && fEndRun<=111) fDataSet = "P3GKF";
  else if(fStartRun>=202 && fEndRun<=1041) fDataSet = "P3HUA";
  else if(fStartRun>=1042 && fEndRun<=2334) fDataSet = "P3JCJ";
  else if(fStartRun>=2337 && fEndRun<=8183) fDataSet = "P3JDY";
  else if(fStartRun>=8184 && fEndRun<=8721) fDataSet = "P3K93";
  else if(fStartRun>=8722 && fEndRun<18589) fDataSet = "P3KJR";
  else if(fStartRun>=18623 && fEndRun<5000000) fDataSet = "P3LQK";
  else if(fStartRun>=30000392 && fEndRun<=30000663) fDataSet = "P3CLR";
  else if(fStartRun>=30001706 && fEndRun<=30001921) fDataSet = "P3DCR";
  else if(fStartRun>=30001966 && fEndRun<=30002469) fDataSet = "P3DNQ";
  else if(fStartRun>=30004106 && fEndRun<=30004516) fDataSet = "P3EVV";
  else if(fStartRun>=30004968 && fEndRun<=30005020) fDataSet = "P3FPW";
  else if(fStartRun>=30005405 && fEndRun<=30005413) fDataSet = "P3FPV";
  else if(fStartRun>=60000000 && fEndRun<60000549) fDataSet = "P3LMF";
  else if(fStartRun>=60000550 && fEndRun<70000000) fDataSet = "P3LQG";
  else if(fStartRun>= 45000000&& fEndRun<=45008360) fDataSet = "P3END";
  else fDataSet = "NA";
}

//Set channel parameters: cryo, string, detector position, channel, detector name
void GATAutoCal::SetChannel()
{
  fChannel.clear();
  fCryo.clear();
  fString.clear();
  fDetector.clear();
  fDetectorName.clear();
  fPulserTagChannel.clear();
  //DS0-DS3 (Mod1)
  if(fStartRun>=0 && fEndRun <18623){
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
	    fChannel.push_back(fMap->GetIDHi(i,j,k));
	    fChannel.push_back(fMap->GetIDLo(i,j,k));
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
	    fChannel.push_back(fMap->GetIDHi(i,j,k));
	    fChannel.push_back(fMap->GetIDLo(i,j,k));
            fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
            fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	  }
	}
      }      
    }
    fStrings = 7;
  }

  //DS4 (Mod2)
  if(fStartRun>=60000001 && fEndRun <70000000){
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
            fChannel.push_back(fMap->GetIDHi(i,j,k));
            fChannel.push_back(fMap->GetIDLo(i,j,k));
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
            fChannel.push_back(fMap->GetIDHi(i,j,k));
            fChannel.push_back(fMap->GetIDLo(i,j,k));
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
            fChannel.push_back(fMap->GetIDHi(i,j,k));
            fChannel.push_back(fMap->GetIDLo(i,j,k));
	    fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
            fDetectorName.push_back(fMap->GetDetectorName(i,j,k));

          }
        }
      }
    }
    fStrings = 7;
  }

  //DS5- (Mod1+Mod2)
  if(fStartRun>=18623 && fEndRun <5000000){
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
	      fChannel.push_back(fMap->GetIDHi(i,j,k));
	      fChannel.push_back(fMap->GetIDLo(i,j,k));
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
	      fChannel.push_back(fMap->GetIDHi(i,j,k));
	      fChannel.push_back(fMap->GetIDLo(i,j,k));
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
	      fChannel.push_back(fMap->GetIDHi(i,j,k));
	      fChannel.push_back(fMap->GetIDLo(i,j,k));
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
	      fChannel.push_back(fMap->GetIDHi(i,j,k));
	      fChannel.push_back(fMap->GetIDLo(i,j,k));
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
	      fChannel.push_back(fMap->GetIDHi(i,j,k));
	      fChannel.push_back(fMap->GetIDLo(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));
	      fDetectorName.push_back(fMap->GetDetectorName(i,j,k));	      
	    }
	  }
	}
      }
    }
    fStrings = 7;
  }

  //STC
  if(fStartRun >=30000392 && fEndRun<=30000492){
    for(Int_t i = 0;i<8;i++){
      Int_t j = i/2;
      fChannel.push_back(i+48);
      fCryo.push_back(1);
      fString.push_back(1);
      fDetector.push_back(j+1);
      fDetectorName.push_back("");
    }
    fStrings = 1;
  }
  if(fStartRun >=30001706 && fEndRun<=30001921){
    for(Int_t i = 0;i<8;i++){
      Int_t j = i/2;
      fChannel.push_back(i+68);
      fCryo.push_back(1);
      fString.push_back(1);
      fDetector.push_back(j+1);
      fDetectorName.push_back("");
    }
    fStrings = 1;
  }
  if(fStartRun >=30001966 && fEndRun<=30002469){
    for(Int_t i = 0;i<8;i++){
      Int_t j = i/2;
      fChannel.push_back(i+66);
      fCryo.push_back(1);
      fString.push_back(1);
      fDetector.push_back(j+1);
      fDetectorName.push_back("");
    }
    fStrings = 1;
  }
  if(fStartRun >=30004106 && fEndRun<=30004516){
    for(Int_t i = 0;i<8;i++){
      Int_t j = i/2;
      fChannel.push_back(i+128);
      fCryo.push_back(1);
      fString.push_back(1);
      fDetector.push_back(j+1);
      fDetectorName.push_back("");
    }
    fStrings = 1;
  }
  if(fStartRun >=30004968 && fEndRun<=30005020){
    for(Int_t i = 0;i<8;i++){
      Int_t j = i/2;
      fChannel.push_back(i+176);
      fCryo.push_back(1);
      fString.push_back(1);
      fDetector.push_back(j+1);
      fDetectorName.push_back("");
    }
    fStrings = 1;
  }
  if(fStartRun >=30005405 && fEndRun<=30005414){
    for(Int_t i = 0;i<8;i++){
      Int_t j = i/2;
      fChannel.push_back(i+114);
      fCryo.push_back(1);
      fString.push_back(1);
      fDetector.push_back(j+1);
      fDetectorName.push_back("");
    }
    fStrings = 1;
  }
  
  //DSPM (Prototype)
  if(fStartRun >= 45000000 && fEndRun<50000000){
    for(Int_t i = 0;i<proChannels;i++){
      fChannel.push_back(Channelpro[i]);
      fDetectorName.push_back(DetectorNamepro[i]);
    }
    for(Int_t j = 0;j<=2;j++){
      if(j == 0){
	for(Int_t k =1;k<=4;k++){
	  fCryo.push_back(1);
	  fCryo.push_back(1);
	  fString.push_back(j);
	  fString.push_back(j);
	  fDetector.push_back(k);
	  fDetector.push_back(k);
	}
      }
      if(j == 1){
        for(Int_t k =1;k<=1;k++){
	  fCryo.push_back(1);
	  fCryo.push_back(1);
          fString.push_back(j);
          fString.push_back(j);
          fDetector.push_back(k);
          fDetector.push_back(k);
        }
      }
      if(j == 2){
        for(Int_t k =1;k<=5;k++){
	  fCryo.push_back(1);
	  fCryo.push_back(1);
          fString.push_back(j);
          fString.push_back(j);
          fDetector.push_back(k);
          fDetector.push_back(k);
        }
      }
    }
    fStrings = 3;
  }

  //Except STC and PM, Set Up Pulser Channels
  if((fStartRun>0 && fEndRun<30000000) || (fStartRun > 60000000 && fEndRun<65000000)){
    vector<uint32_t> pulsertagchan = fMap->GetPulserChanList();
    for(Int_t i =0;i<(Int_t)pulsertagchan.size();i++){
      Int_t chan = pulsertagchan[i];
      fPulserTagChannel.push_back(chan);
    }
  }

  fChannels = fChannel.size();
  fPulserTagChannels = fPulserTagChannel.size();
  fTotalChannels = fChannels+fPulserTagChannels;
}

//Get GAT Revision
Int_t GATAutoCal::GetGATRev(){
  Int_t gatrev = 0;
  //Except STC and PM
  if((fStartRun>0 && fEndRun<30000000) || (fStartRun > 60000000 && fEndRun<65000000)){
    GATDataSet ds(fStartRun);
    TChain *mjdT = ds.GetGatifiedChain(false);
    mjdT->SetBranchStatus("*",1);
    
    fMjdTree->SetBranchAddress("gatrev", &fGATRev);
    fMjdTree->GetEntry(0);
    gatrev = (Int_t)fGATRev;
  }
  return gatrev;
}

//Get Time Information from GAT
Double_t GATAutoCal::GetStartTime(){
  GATDataSet ds(fStartRun);
  TChain *mjdT = ds.GetGatifiedChain(false);
  vector<Double_t> *timestamp = NULL;
  mjdT->SetBranchStatus("*",1);
  mjdT->SetBranchAddress("timestamp", &timestamp);
  mjdT->SetBranchAddress("startTime",&fMTStartTime);
  mjdT->GetEntry(0);
  //Double_t starttime = fMTStartTime;
  Double_t starttime = timestamp->at(0);
  return starttime;
}

Double_t GATAutoCal::GetStartTimeMT(){
  GATDataSet ds(fStartRun);
  TChain *mjdT = ds.GetGatifiedChain(false);
  mjdT->SetBranchStatus("*",1);
  fMTTime = NULL;
  mjdT->SetBranchAddress("timeMT", &fMTTime);
  mjdT->GetEntry(0);
  Double_t time = fMTTime->at(0);
  return time;
}

Int_t GATAutoCal::GetStartDateMT(){
  GATDataSet ds(fStartRun);
  TChain *mjdT = ds.GetGatifiedChain(false);
  mjdT->SetBranchStatus("*",1);
  fMTDate = NULL;
  mjdT->SetBranchAddress("dateMT",&fMTDate);
  mjdT->GetEntry(0);
  Int_t date = fMTDate->at(0);
  return date;
}

Double_t GATAutoCal::GetStopTime(){
  GATDataSet ds(fEndRun);
  TChain *mjdT = ds.GetGatifiedChain(false);
  Int_t entries = mjdT->GetEntries();
  vector<Double_t> *timestamp = NULL;
  mjdT->SetBranchStatus("*",1);
  mjdT->SetBranchAddress("timestamp", &timestamp);
  mjdT->SetBranchAddress("stopTime",&fMTStopTime);
  mjdT->GetEntry(entries-1);
  //Double_t stoptime = fMTStopTime;
  Double_t stoptime = timestamp->at(0);
  return stoptime;
}

Double_t GATAutoCal::GetStopTimeMT(){
  GATDataSet ds(fEndRun);
  TChain *mjdT = ds.GetGatifiedChain(false);
  Int_t entries = mjdT->GetEntries();
  mjdT->SetBranchStatus("*",1);
  fMTTime = NULL;
  mjdT->SetBranchAddress("timeMT", &fMTTime);
  mjdT->GetEntry(entries-1);
  Double_t time = fMTTime->at(0);
  return time;
}

Int_t GATAutoCal::GetStopDateMT(){
  GATDataSet ds(fEndRun);
  TChain *mjdT = ds.GetGatifiedChain(false);
  Int_t entries = mjdT->GetEntries();
  mjdT->SetBranchStatus("*",1);
  fMTDate = NULL;
  mjdT->SetBranchAddress("dateMT",&fMTDate);
  mjdT->GetEntry(entries-1);
  Int_t date = fMTDate->at(0);
  return date;
}

//Fill Histogram of each channel
TH1D* GATAutoCal::FillHisto(TChain* mTree, string EnergyName, Int_t Channel, Int_t Bin, Double_t Low, Double_t Up){
  TH1D *h = new TH1D(Form("%s%d", EnergyName.c_str(),Channel),Form("Channel %d",Channel), Bin, Low, Up);
  TCut cut1;
  if((fStartRun == 9913 && fEndRun == 9926) || (fStartRun == 13071 && fEndRun == 13074)){
    cut1 = Form("channel == %d && (rawWFMax > 200 || rawWFMin < -20)",Channel);
  }else if(fStartRun >= 60002368 && fEndRun <= 60002372 && Channel==1175){
    cut1 = Form("channel == %d && (rawWFMax > 50 || rawWFMin < -200)",Channel);
  }else if(fStartRun >= 60002368 && fEndRun <= 60002372 && Channel==1106){
    cut1 = Form("channel == %d && (rawWFMax > 200 || rawWFMin < -200)",Channel);
  }else if(fStartRun >= 60002368 && fEndRun <= 60002372 && Channel==1107){
    cut1 = Form("channel == %d && (rawWFMax > 80 || rawWFMin < -50)",Channel); 
  }else if(fStartRun >= 60002368 && fEndRun <= 60002372){
    cut1 = Form("channel == %d && (rawWFMax > 400 || rawWFMin < -400)",Channel);
  }else if(fStartRun == 4171 && fEndRun == 4201){
    cut1 = Form("channel == %d",Channel);
  }else if(EnergyName == "trapENFBL"){
    cut1 = Form("channel == %d && trapENF>200 && RawWFblSlope>-0.0002 && RawWFblSlope<0.0002 && EventDC1Bits == 0",Channel);
  }else{
    cut1 = Form("channel == %d && EventDC1Bits == 0",Channel);
  }
  mTree->Draw(Form("%s>>%s%d",EnergyName.c_str(),EnergyName.c_str(),Channel),cut1);
  return h;
}

//Save into a ROOT file
void GATAutoCal::SaveHisto(TH1D* Hist, string FileName, string Option){
  TFile fhist(Form("%s",FileName.c_str()),Option.c_str());
  Hist->Write();
  cout << "Save histograms at: " << FileName.c_str() << endl;
}

//Load from a ROOT file
TH1D* GATAutoCal::LoadHisto(string FileName, string HistoName){
  TFile *fhist = TFile::Open(Form("%s",FileName.c_str()),"read");
  TH1D *h = (TH1D*)fhist->Get(Form("%s",HistoName.c_str()));
  return h;
  cout << "Load histograms from: "<< FileName.c_str() << endl;
}


Int_t GATAutoCal::MultiPeakFit(TH1F *Hist, Double_t scaleenergy, Double_t scaleamps,string FitName, vector<string>* ParName, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov, vector<Double_t>* CalPeak){

  ParName->clear();
  Par->clear();
  ParErr->clear();
  Cov->clear();

  TH1F* hist = (TH1F*)Hist->Clone();
  //TH1D* hist1 = (TH1D*)Hist->Clone();
  if(hist==NULL) {
    cout <<  "Histogram not found!" << endl;
    gSystem->Abort();
  }

  hist = (TH1F*) hist->Clone();
  GATMultiPeakFitter fitter;
  //Add the regions. This is AddRegion(low bin, high bin, list of peak energies in keV)
  fitter.AddRegion(200, 340, {238.632, 240.986, 277.371, 300.087});//4 peaks
  fitter.AddRegion(540, 620, {583.191}); //1 peaks
  fitter.AddRegion(716, 900, {727.330, 785.37, 860.557});//3 peaks
  //fitter.AddRegion(716, 900, {727.330, 860.557});//2 peaks
  fitter.AddRegion(2550, 2800, {2614.533}); //1 peak
 
  // Set up the boundary of each region
  vector<Double_t> Low;
  vector<Double_t> Up;
  Low.push_back(200);
  Low.push_back(540);
  Low.push_back(716);
  Low.push_back(2550);
  Up.push_back(340);
  Up.push_back(620);
  Up.push_back(900);
  Up.push_back(2800);

  // initial BG parameters
  fitter.SetBGPars(0, 3000, -9, 0);
  fitter.FixBGPar(0, 2, 0);
  fitter.SetBGPars(1, 500, 0, 0);
  fitter.FixBGPar(1, 1, 0);
  fitter.FixBGPar(1, 2, 0);
  fitter.SetBGPars(2, 250, 0, 0);
  fitter.FixBGPar(2, 1, 0);
  fitter.FixBGPar(2, 2, 0);
  fitter.SetBGPars(3, 4, 0, 0);
  fitter.FixBGPar(3, 1, 0);
  fitter.FixBGPar(3, 2, 0);


  // Start setting initial parameters
  //SetParFunction(parameter, function type, list of parameters)
  //Initial amplitudes for each peak. Limit these to >0
  fitter.SetParFunction(GATMultiPeakFitter::kAmp, GATMultiPeakFitter::kFree, {/*region 1*/ 160000, 12000, 8000, 10000, /*region2*/ 80000, /*region3*/ 15000, 3000, 10000, /*region4*/ 35000}); //25 total pars

  for(size_t i=0; i<fitter.NPeaks(); i++) fitter.LimitPar(GATMultiPeakFitter::kAmp, i, 10., (double) NAN);

  fitter.SetParFunction(GATMultiPeakFitter::kMu, GATMultiPeakFitter::kLinear, {0, 1.});
  fitter.SetParFunction(GATMultiPeakFitter::kSig, GATMultiPeakFitter::kRootQuad, {0.2, 0.017, 0.0002});

  fitter.SetParFunction(GATMultiPeakFitter::kFt, GATMultiPeakFitter::kConst, {0.3});
  fitter.SetParFunction(GATMultiPeakFitter::kTau, GATMultiPeakFitter::kLinear, {0., 0.0006});
  fitter.LimitPar(GATMultiPeakFitter::kFt, 0, 0., 1.);
  fitter.FixPar(GATMultiPeakFitter::kTau, 0, 0.);
  fitter.LimitPar(GATMultiPeakFitter::kTau, 1, 0.00001, 0.01);

  fitter.SetParFunction(GATMultiPeakFitter::kFht, GATMultiPeakFitter::kConst, {0.1});
  fitter.SetParFunction(GATMultiPeakFitter::kTauHT, GATMultiPeakFitter::kLinear, {0., 0.0006});
  fitter.LimitPar(GATMultiPeakFitter::kFht, 0, 0., 1.);
  fitter.LimitPar(GATMultiPeakFitter::kTauHT, 1, 0., 0.1);

  //comment the next 3 lines to turn on the HE tail
  fitter.SetParFunction(GATMultiPeakFitter::kTauHT, GATMultiPeakFitter::kConst, {0.5});
  fitter.FixPar(GATMultiPeakFitter::kFht, 0, 0);
  fitter.FixPar(GATMultiPeakFitter::kTauHT, 0, 0.5);
  fitter.SetParFunction(GATMultiPeakFitter::kHs, GATMultiPeakFitter::kStepHeightFun, {650, 0.1, -0.5});
  fitter.LimitPar(GATMultiPeakFitter::kHs, 2, -1, 0);

  fitter.ScaleEnergy(scaleenergy);
  fitter.ScaleAmps(scaleamps);

  fitter.SetHists(hist);

  // Set up a hybrid montecarlo fit to improve the initial parameters

  GATHybridMonteCarlo hmc;
  hmc.SetNLLFunc(dynamic_pointer_cast<ROOT::Math::IGradientFunctionMultiDim>(fitter.SetPoissonLLFCN()));
  // uncomment the following lines to generate a TTree with each MCMC step
  hmc.SetOutputFile(Form("hmc_%s.root",FitName.c_str()));
  hmc.SetRecordPaths();
  //hmc.SetOutputFile(("results_"+chname+".root").c_str());
  //hmc.SetRecordPaths();
  hmc.SetParameters(std::vector<double>(fitter.Parameters(), fitter.Parameters()+fitter.NPar()));
  // Set the HMC parameters. Step size is leapfrog step size. Step length is number of leapfrog steps for each MCMC step. NSteps is number of MCMC steps. Adapt step size and parameter scales will automatically adjust step size and individual parameter scales for each step
  hmc.SetStepSize(0.02);
  hmc.SetStepLength(50);
  hmc.SetNSteps(200);
  hmc.SetAdaptStepSize();
  hmc.SetAdaptParScales();
  hmc.SetLimits(fitter.GetParLimits());
  // output the random seed to be used. This might be useful for troubleshooting
  cout << "Random seed: " << hmc.GetCurrentSeed() << endl;
  // true -> verbose output that llists the negative loglikelihood at each step and the step size
  hmc.DoMCMC(true);
  // fitter.SetParameters(hmc.GetLikeliestPars().data());
  //perform a minuit fit using the most likely parameters found during the HMC as the initial parameters. Use minos error estimation
  fitter.SetParameters(hmc.GetLikeliestPars().data());  
  for(size_t i=0; i<fitter.NPeaks(); i++) fitter.LimitPar(GATMultiPeakFitter::kAmp, i, 10., (double) NAN);

  fitter.FitToHists(0., true, false, true);

  //release mu; get a better fitting result of each peak
  fitter.SetParFunction(GATMultiPeakFitter::kMu, GATMultiPeakFitter::kFree);
  fitter.FitToHists(0., true, false, true);
  Int_t isValid = fitter.GetResult().IsValid();
  if(isValid == 1){
    cout << "fitting is valid" << endl;
  }else{
    cout << "fitting is not valid" << endl;
  }

  for(size_t  i = 0;i<fitter.NPar();i++){
    ParName->push_back(fitter.ParameterName(i));
    Par->push_back(fitter.GetResult().Parameter(i));
    ParErr->push_back(fitter.GetResult().Error(i));
  }

  ofstream fout(Form("parameters_%s.txt",FitName.c_str()),ios::app);
  fout.precision(15);
  for(size_t i = 0;i<fitter.NPar();i++){
    fout << i<< " " <<  ParName->at(i).c_str() << " " << Par->at(i) << " " << ParErr->at(i) << endl;
  }

  vector<Double_t> calpeak = fitter.GetPeaks();
  for(size_t  i = 0;i<calpeak.size();i++){
    CalPeak->push_back(calpeak.at(i));
  }

  TCanvas *c1 = new TCanvas("c1");
  c1->SetLogy();
  for(size_t i=0; i<fitter.NRegions(); i++){
    hist->GetXaxis()->SetRangeUser(Low.at(i)*scaleenergy, Up.at(i)*scaleenergy);
    fitter.DrawRegion(i);
    c1->Print(Form("spectrum_%s_%d.pdf",FitName.c_str(),(Int_t)i));
  }

  return isValid;
}

Int_t GATAutoCal::LinearFit(vector<Double_t> Px, vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, string FitName, string TitleName, Double_t EnergyROI, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov){

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  
  //check if there are enough inputs
  Int_t nCal=Px.size();
  if(nCal<2){
    cout<<"Cannot fit a calibration with less than 2 points!"<<endl;
  }

  TGraphErrors fCalGraph;
  fCalGraph.Set(0);
  fCalGraph.SetTitle(Form("%s;%s;Energy (keV)",TitleName.c_str(),fEnergyName.c_str()));
  fCalGraph.SetMarkerStyle(8);
  fCalGraph.SetMarkerSize(2); 
  fCalGraph.SetMarkerColor(4);
  for(Int_t i = 0;i<nCal;i++){
    fCalGraph.SetPoint(i,Px.at(i),Py.at(i));
    fCalGraph.SetPointError(i,PxErr.at(i),PyErr.at(i));
  }


  TGraphErrors fDeltaGraph;
  fDeltaGraph.Set(0);
  fDeltaGraph.SetTitle(Form("%s;Energy (keV);#Delta E (keV)",TitleName.c_str()));
  fDeltaGraph.SetMarkerStyle(8);
  fDeltaGraph.SetMarkerSize(2);
  fDeltaGraph.SetMarkerColor(4);

  TF1 *fLinearFunction;
  fLinearFunction = new TF1("fLinearFunction",Slope, 0.9*Px.at(0),1.1*Px.at(nCal-1),2);
  fLinearFunction->SetLineColor(2);
  fLinearFunction->SetLineWidth(2);
  fLinearFunction->SetLineStyle(2);
  
  Double_t fScale = (Py.at(nCal-1)-Py.at(0))/(Px.at(nCal-1)-Px.at(0));
  Double_t fOffset = Py.at(0)-fScale*Px.at(0);
  
  fLinearFunction->SetParameters(fOffset,fScale);
  fCalGraph.Fit(fLinearFunction,"qremf");

  Double_t chisquare = fLinearFunction->GetChisquare();
  TFitResultPtr r = fCalGraph.Fit(fLinearFunction,"S");
  TMatrixDSym cov = r->GetCovarianceMatrix();
  Int_t isValid = r->IsValid();

  fScale=fLinearFunction->GetParameter(1);
  Double_t fScaleErr=fLinearFunction->GetParError(1);
  fOffset=fLinearFunction->GetParameter(0);
  Double_t fOffsetErr=fLinearFunction->GetParError(0);

  Par->push_back(fOffset); //#0
  ParErr->push_back(fOffsetErr); 
  Par->push_back(fScale); // #1
  ParErr->push_back(fScaleErr);

  Double_t LinearROI = (EnergyROI-fOffset)/fScale; 
  Double_t p[2];
  p[0] = -1/fScale;
  p[1] = -(EnergyROI-fOffset)/(fScale*fScale);

  Double_t errsum = 0;
  for(Int_t i = 0;i<2;i++){
    for(Int_t j = 0;j<2;j++){
      if(i==j){
        errsum = errsum + p[i]*TMath::Abs(cov(i,j))*p[j];
      }else{
        errsum = errsum + p[i]*cov(i,j)*p[j];
      }
    }
  }
  Double_t LinearROIErr = TMath::Sqrt(errsum);
  Par->push_back(LinearROI); //#2
  ParErr->push_back(LinearROIErr);
  Par->push_back(chisquare); //#3
  ParErr->push_back(0);

  for(Int_t k = 0;k<2;k++){
    for(Int_t j = 0;j<2;j++){
      if(k==j){
	Cov->push_back(abs(cov(k,j)));
      }else{
	Cov->push_back(cov(k,j));
      }
    }
  }

  ofstream ferror(Form("error_%s.txt",FitName.c_str()));
  //////////////////////////////////////////////
  for(Int_t i = 0;i<nCal;i++){
    Double_t diff =  Py.at(i)-(fOffset+fScale*Px.at(i));
    Double_t differr = TMath::Sqrt(pow(fOffsetErr,2)+pow(fScaleErr*Px.at(i),2)+pow(PxErr.at(i)*fScale,2));
    fDeltaGraph.SetPoint(i,Py.at(i),diff);
    fDeltaGraph.SetPointError(i,0,differr);
    if(abs(diff)>0.5){
      ferror << FitName.c_str() << " " << fEnergyName.c_str() << " " << i << " " << Py.at(i) << " " << Px.at(i) << " " << PxErr.at(i)<<  " " <<  diff << " " << differr << endl;
    }
  }
  
  TCanvas *c1 = new TCanvas("c1");
  fCalGraph.Draw();
  c1->Print(Form("linear_%s.pdf",FitName.c_str()));
  fDeltaGraph.Draw("AP");
  c1->Print(Form("delta_%s.pdf",FitName.c_str()));

  return isValid;
}


void GATAutoCal::SetUpProvenance(std::string title, std::string yourname,
				 Int_t gatrev,
				 //		 std::string detectorid,
				 int startrun, int endrun,
				 int coverstartrun, int coverendrun,
				 int startdate, int enddate,
				 int coverstartdate, int coverenddate)
{
  cout << gatrev << endl;
  //Construct static provenance
  fACProvenance.SetTitle(title);
  fACProvenance.SetPublisher(yourname);
  fACProvenance.SetRecordType(krtAutomatic); //Automatic calibration record
  fACProvenance.SetValid(true);              //Invalidate later, if necessary
  fACProvenance.SetSystemIdentifier(fDataSet); //from GATAutoCal
  //  fACProvenance.SetDetectorIdentifier(detectorid);

  //Dates
  //Start and end time of calibration data
  //time_t ltmpDate = MJDB::ISOStringToDate(startdate); //<-insert start date here
  //time_t utmpDate = MJDB::ISOStringToDate(enddate);//<-insert end date here
  fACProvenance.SetRunDate(MJDB::NumberToString(startdate),
		          MJDB::NumberToString(enddate)); 

  //Runs from which the calibrations are derived
  fACProvenance.SetRunIdentifier(MJDB::NumberToString(startrun),
				MJDB::NumberToString(endrun));


  //Coverage Dates - where this calibration can be applied
  //Start and End Date-time
  //  ltmpDate = MJDB::ISOStringToDate(startdate); //<-insert start date here
  if(coverendrun == 4999999){   
    fACProvenance.SetDateCoverage(MJDB::NumberToString(coverstartdate)); //no end date
  }else{
    fACProvenance.SetDateCoverage(MJDB::NumberToString(coverstartdate),MJDB::NumberToString(coverenddate));
  }

  //Run Coverage - Runs to which the calibration applies
  fACProvenance.SetRunCoverage(MJDB::NumberToString(coverstartrun),MJDB::NumberToString(coverendrun)); //no ending run #

  //What software generated the calibrations (fictional names for testing)
  fACProvenance.SetSoftwareSource(Form("GATAutoCal/GATRev%d",gatrev)); //<- GAT Version?
  fACProvenance.SetSoftwareCoverage("GAT/MJDBParameters");

  //Parameter source
  //fACProvenance.SetParametersSource(fEnergyName); //from the class
}

//This creates a single new record of energy calibration in the database
    void GATAutoCal::PutECalMJDB(int channel, std::string detectorid,
			     double scale, double scalerr,
			     double offset, double offerr,
			     std::vector<double> covariance)
{
  
  //REMOVE this after testing
  //This line adds the non-standard DB to the AnalysisDoc readout
  //fACAnalysisDB.SetAlternateDB("", "mjdbsandbox","","","");
  fACAnalysisDB.SetAlternateDB("", "mjdbatest","","","");
  
  MJDB::EnergyCalibration myE_On;
  fACProvenance.SetParametersSource(fEnergyName); //from the class
  fACProvenance.SetParametersType(myE_On.GetPType());//set the Parameter type
  fACProvenance.SetChannelIdentifier(MJDB::NumberToString(channel));
  fACProvenance.SetDetectorIdentifier(detectorid);


  myE_On.Scale.SetValue(scale);
  myE_On.Scale.SetUncertainty(scalerr);
  
  myE_On.Offset.SetValue(offset);
  myE_On.Offset.SetUncertainty(offerr);
  myE_On.Covariance.SetValueList(covariance);
  
  //Build a database record
  fACProvenance.SetDBProvenance(fACAnalysisDB);
  if (myE_On.SetDBValue(fACAnalysisDB)) {      //Make sure we save the parameters    
    int status = fACAnalysisDB.New_MJAnalysis_Record();
    if (status) {
      std::string statusMessage;
      fACAnalysisDB.GetDBStatus(statusMessage);
      cout << "\"" << statusMessage << "\"" << endl;
      return; //Error
    } 
  } else {
    //Reach here only if the parameter record is incomplete. Bad error
    cout << "Energy Calibration Parameter record not completely initialized."
	 <<"  Try again."
	 << endl;
  }
}
