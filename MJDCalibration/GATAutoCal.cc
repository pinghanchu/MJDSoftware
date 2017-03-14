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


GATAutoCal::GATAutoCal(Int_t StartRun, Int_t EndRun) : fDS(StartRun,EndRun)
{
  fStartRun = StartRun;
  fEndRun = EndRun;
  if((fStartRun>0 && fEndRun<30000000) || (fStartRun > 60000000 && fEndRun<65000000)){
    fMap = fDS.GetChannelMap();
    fMjdTree = fDS.GetGatifiedChain(false);
  }

  SetParameters();
  SetCalibrationPeak();
  for(Int_t i = 0;i<116;i++){
    TrapE[i] = NULL;
  }
}

void GATAutoCal::SetMjdTree(TChain *mjdTree){
  fMjdTree = mjdTree;
}
void GATAutoCal::SetParameters()
{
  if(fStartRun >= 2337 && fEndRun < 8183){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad0[i]);
      fEnriched.push_back(OrtecBeGe0[i]);
      if(fEndRun<3464) fPulserCal.push_back(PulserCal0[i]);
      else fPulserCal.push_back(PulserCal1[i]);
      if(fEndRun<3464) fPulser.push_back(Pulser0[i]);
      else fPulser.push_back(Pulser1[i]);
    }
  }
  //D1
  if(fStartRun >= 8722 && fEndRun<14503){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad3[i]);
      fEnriched.push_back(OrtecBeGe0[i]);
      fPulserCal.push_back(PulserCal3[i]);
      fPulser.push_back(Pulser3[i]);
    }
  }
  //DS2
  if(fStartRun >= 14503 && fEndRun<16836){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad32[i]);
      fEnriched.push_back(OrtecBeGe0[i]);
      fPulserCal.push_back(PulserCal3[i]);
      fPulser.push_back(Pulser3[i]);
    }
  }

  //DS3
  if(fStartRun >= 16836 && fEndRun<18589){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad4[i]);
      fEnriched.push_back(OrtecBeGe0[i]);
      fPulserCal.push_back(PulserCal3[i]);
      fPulser.push_back(Pulser3[i]);
    }
  }
  if(fStartRun >= 18590 && fEndRun<5000000){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad6[i]);
      fEnriched.push_back(OrtecBeGe0[i]);
      fPulserCal.push_back(PulserCal6[i]);
      fPulser.push_back(Pulser6[i]);
    }
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad7[i]);
      fEnriched.push_back(OrtecBeGe4[i]);
      fPulserCal.push_back(PulserCal7[i]);
      fPulser.push_back(Pulser7[i]);
    }
  }
  if(fStartRun >= 30000000 && fEndRun<40000000){
    for(Int_t i = 0;i<8;i++){
      fGoodBad.push_back(1);
      fEnriched.push_back(0);
      fPulserCal.push_back(0);
      fPulser.push_back(0);
    }
  }

  if(fStartRun >= 60000001 && fEndRun<60000550){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad5[i]);
      fEnriched.push_back(OrtecBeGe4[i]);
      fPulserCal.push_back(PulserCal5[i]);
      fPulser.push_back(Pulser5[i]);
    }
  }

  if(fStartRun >= 60000550 && fEndRun<70000000){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad5[i]);
      fEnriched.push_back(OrtecBeGe4[i]);
      fPulserCal.push_back(PulserCal5[i]);
      fPulser.push_back(Pulser5[i]);
    }
  }

  if(fStartRun >= 45000000 && fEndRun<50000000){
    for(Int_t i = 0;i<proChannels;i++){
      cout << i << " " << GoodBadpro[i] << endl;
      fGoodBad.push_back(GoodBadpro[i]);
      fEnriched.push_back(0);
      fPulserCal.push_back(PulserCalpro[i]);
      fPulser.push_back(Pulserpro[i]);
    }
  }
  SetChannel();
  SetDataSet();
}
void GATAutoCal::SetRunRange(Int_t StartRun,Int_t EndRun){
  fStartRun = StartRun;
  fEndRun = EndRun;
}
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

void GATAutoCal::SetCalibrationPeak(vector<Double_t> gammakeV){
  fCalibrationPeak.erase(fCalibrationPeak.begin(),fCalibrationPeak.end());
  Int_t nIndex = gammakeV.size();
  for(Int_t i = 0;i<nIndex;i++){
    fCalibrationPeak.push_back(gammakeV.at(i));
  }
  fCalibrationPeaks = fCalibrationPeak.size();
}

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

void GATAutoCal::SetChannel()
{
  fChannel.erase(fChannel.begin(),fChannel.end());
  fCryo.clear();
  fString.clear();
  fDetector.clear();
  fDetectorName.clear();
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

  if(fStartRun>=18623 && fEndRun <5000000){
    for(Int_t i =1;i<3;i++){
      
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
  if(fStartRun >=30000392 && fEndRun<=30000492){
    for(Int_t i = 0;i<8;i++){
      Int_t j = i/2;
      fChannel.push_back(i+48);
      fCryo.push_back(0);
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
      fCryo.push_back(0);
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
      fCryo.push_back(0);
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
      fCryo.push_back(0);
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
      fCryo.push_back(0);
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
      fCryo.push_back(0);
      fString.push_back(1);
      fDetector.push_back(j+1);
      fDetectorName.push_back("");
    }
    fStrings = 1;
  }
  
  if(fStartRun >= 45000000 && fEndRun<50000000){
    for(Int_t i = 0;i<proChannels;i++){
      fChannel.push_back(Channelpro[i]);
      cout << i << " " << Channelpro[i] <<endl;
    }
    for(Int_t j = 0;j<=2;j++){
      /*
      if(j == 0){
	for(Int_t k =1;k<=4;k++){
	  fCryo.push_back(1);
	  fCryo.push_back(1);
	  fString.push_back(j);
	  fString.push_back(j);
	  fDetector.push_back(k);
	  fDetector.push_back(k);
	  fDetectorName.push_back("");
	  fDetectorName.push_back("");
 
	}
      }
      if(j == 1){
        for(Int_t k =1;k<=1;k++){
	  fCryo.push_back(1);
	  fCryo.push_back(1);
          fString.push_back(j+1);
          fString.push_back(j+1);
          fDetector.push_back(k);
          fDetector.push_back(k);
          fDetectorName.push_back("");
          fDetectorName.push_back("");

        }
      }
      if(j == 1){
        for(Int_t k =1;k<=4;k++){
	  fCryo.push_back(1);
	  fCryo.push_back(1);
          fString.push_back(j+1);
          fString.push_back(j+1);
          fDetector.push_back(k);
          fDetector.push_back(k);
          fDetectorName.push_back("");
          fDetectorName.push_back("");

        }
      }
      */
     if(j == 0){
	for(Int_t k =1;k<=4;k++){
	  fCryo.push_back(0);
	  fCryo.push_back(0);
	  fString.push_back(j);
	  fString.push_back(j);
	  fDetector.push_back(k);
	  fDetector.push_back(k);
	  fDetectorName.push_back("");
	  fDetectorName.push_back("");
 
	}
      }
      if(j == 1){
        for(Int_t k =1;k<=1;k++){
	  fCryo.push_back(0);
	  fCryo.push_back(0);
          fString.push_back(j);
          fString.push_back(j);
          fDetector.push_back(k);
          fDetector.push_back(k);
          fDetectorName.push_back("");
          fDetectorName.push_back("");

        }
      }
      if(j == 2){
        for(Int_t k =1;k<=5;k++){
	  fCryo.push_back(0);
	  fCryo.push_back(0);
          fString.push_back(j);
          fString.push_back(j);
          fDetector.push_back(k);
          fDetector.push_back(k);
          fDetectorName.push_back("");
          fDetectorName.push_back("");
        }
      }
    }
    fStrings = 3;
  }
  fChannels = fChannel.size();
  cout << fChannels << endl;
  if((fStartRun>0 && fEndRun<30000000) || (fStartRun > 60000000 && fEndRun<65000000)){
    fPulserTagChannel.erase(fPulserTagChannel.begin(),fPulserTagChannel.end());
    vector<uint32_t> pulsertagchan = fMap->GetPulserChanList();
    for(Int_t i =0;i<(Int_t)pulsertagchan.size();i++){
      Int_t chan = pulsertagchan[i];
      fPulserTagChannel.push_back(chan);
    }
  }
  fPulserTagChannels = fPulserTagChannel.size();
  fTotalChannels = fChannels+fPulserTagChannels;
}

Int_t GATAutoCal::GetGATRev(){
  Int_t gatrev = 0;
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

Double_t GATAutoCal::GetStartTime(){
  GATDataSet ds(fStartRun);
  TChain *mjdT = ds.GetGatifiedChain(false);
  vector<Double_t> *timestamp = NULL;
  mjdT->SetBranchStatus("*",1);
  mjdT->SetBranchAddress("timestamp", &timestamp);
  mjdT->SetBranchAddress("startTime",&fMTStartTime);
  mjdT->GetEntry(0);
  Double_t starttime = fMTStartTime;
  //Double_t starttime = timestamp->at(0);
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
  Int_t tempRun = fEndRun;
  Int_t ifExist = 0;
  
  if(fEndRun < 45000000 || fEndRun > 60000000 ){
    if(ifstream(Form("./built/%s/OR_run%d.root",fDataSet.c_str(),tempRun))){
      ifExist = 1;
    }else{
      ifExist = 0;
    }

    while(ifExist == 0 && tempRun>=fStartRun){
      cout << tempRun <<  endl;
      tempRun--;
      if(ifstream(Form("./built/%s/OR_run%d.root",fDataSet.c_str(),tempRun))){
	ifExist = 1;
      }else{
	ifExist = 0;
      }
    }
  }
  
  GATDataSet ds(tempRun);
  TChain *mjdT = ds.GetGatifiedChain(false);
  Int_t entries = mjdT->GetEntries();
  vector<Double_t> *timestamp = NULL;
  mjdT->SetBranchStatus("*",1);
  mjdT->SetBranchAddress("timestamp", &timestamp);
  mjdT->SetBranchAddress("stopTime",&fMTStopTime);
  /*
  Double_t stoptime = 0;
  Int_t index = entries-1;
  while(stoptime == 0){
    index--;
    mjdT->GetEntry(index);
    Double_t stoptime = fMTStopTime;
    Double_t time_s = timestamp->at(0);
    cout << index << " " << stoptime << " " << time_s << endl;
    }*/
  mjdT->GetEntry(entries-1);
  //Double_t stoptime = timestamp->at(0);
  Double_t stoptime = fMTStopTime;

  return stoptime;
}

Double_t GATAutoCal::GetStopTimeMT(){
  Int_t tempRun = fEndRun;
  Int_t ifExist = 0;
  if(fEndRun < 45000000 || fEndRun > 60000000){
    if(ifstream(Form("./built/%s/OR_run%d.root",fDataSet.c_str(),tempRun))){
      ifExist = 1;
    }else{
      ifExist = 0;
    }
    while(ifExist == 0){
      cout << tempRun <<  endl;
      tempRun--;
      if(ifstream(Form("./built/%s/OR_run%d.root",fDataSet.c_str(),tempRun))){
	ifExist = 1;
      }else{
	ifExist = 0;
      }
    }
  }
  GATDataSet ds(tempRun);

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
  Int_t tempRun = fEndRun;
  Int_t ifExist = 0;
  if(fEndRun < 45000000 || fEndRun > 60000000){
    if(ifstream(Form("./built/%s/OR_run%d.root",fDataSet.c_str(),tempRun))){
      ifExist = 1;
    }else{
      ifExist = 0;
    }
    while(ifExist == 0){
      cout << tempRun <<  endl;
      tempRun--;
      if(ifstream(Form("./built/%s/OR_run%d.root",fDataSet.c_str(),tempRun))){
	ifExist = 1;
      }else{
	ifExist = 0;
      }
    }
  }
  GATDataSet ds(tempRun);
  TChain *mjdT = ds.GetGatifiedChain(false);
  Int_t entries = mjdT->GetEntries();
  mjdT->SetBranchStatus("*",1);
  fMTDate = NULL;
  mjdT->SetBranchAddress("dateMT",&fMTDate);
  mjdT->GetEntry(entries-1);
  Int_t date = fMTDate->at(0);
  return date;
}

void GATAutoCal::FillHisto(TChain* mTree, Int_t Bin1, Double_t Low1, Double_t Up1, Int_t Bin2, Double_t Low2, Double_t Up2){
  cout << fEnergyName << endl;
  string pos;
  for(Int_t i=0;i<fChannels;i++){
   
    TH1D *h;
    pos = Form("%d%d%d",fCryo.at(i),fString.at(i),fDetector.at(i));
    cout << pos.c_str() << " " << fChannel.at(i) << endl;
    if(i%2 == 0){
      h= new TH1D(Form("%s%s%d",fEnergyName.c_str(),pos.c_str(),fChannel.at(i)),Form("Channel %d (C%dP%dD%d) in Run %d - %d",fChannel.at(i),fCryo.at(i),fString.at(i),fDetector.at(i),fStartRun,fEndRun), Bin1, Low1,Up1);
    }else{
      h= new TH1D(Form("%s%s%d",fEnergyName.c_str(),pos.c_str(),fChannel.at(i)),Form("Channel %d (C%dP%dD%d) in Run %d - %d",fChannel.at(i),fCryo.at(i),fString.at(i),fDetector.at(i),fStartRun,fEndRun), Bin2, Low2,Up2);
    }
    TCut cut1;

    size_t found = fEnergyName.find("trapENFBL");
    cout << found << endl;
    if((fStartRun == 9913 && fEndRun == 9926) || (fStartRun == 13071 && fEndRun == 13074)){
      cut1 = Form("channel == %d && (rawWFMax > 200 || rawWFMin < -20)",fChannel.at(i));
    }else if(fStartRun >= 60002368 && fEndRun <= 60002372 && fChannel.at(i)==1175){
      cut1 = Form("channel == %d && (rawWFMax > 50 || rawWFMin < -200)",fChannel.at(i));
    }else if(fStartRun >= 60002368 && fEndRun <= 60002372 && fChannel.at(i)==1106){
      cut1 = Form("channel == %d && (rawWFMax > 200 || rawWFMin < -200)",fChannel.at(i));
    }else if(fStartRun >= 60002368 && fEndRun <= 60002372 && fChannel.at(i)==1107){
      cut1 = Form("channel == %d && (rawWFMax > 80 || rawWFMin < -50)",fChannel.at(i)); 
   }else if(fStartRun >= 60002368 && fEndRun <= 60002372){
      cut1 = Form("channel == %d && (rawWFMax > 400 || rawWFMin < -400)",fChannel.at(i));
    }else if(fStartRun == 4171 && fEndRun == 4201){
      cut1 = Form("channel == %d",fChannel.at(i));
    }else if(fEnergyName == "trapENFBL" || (found>0 && found < 10)){
      cut1 = Form("channel == %d && trapENF>200 && RawWFblSlope>-0.0002 && RawWFblSlope<0.0002",fChannel.at(i));
    }else{
      cut1 = Form("channel == %d",fChannel.at(i));
    }
    mTree->Draw(Form("%s>>%s%s%d",fEnergyName.c_str(),fEnergyName.c_str(),pos.c_str(),fChannel.at(i)),cut1);
    TrapE[i] = (TH1D*)h->Clone();
    TrapE[i]->SetTitle(Form("Channel %d (C%dP%dD%d) in Run %d - %d;%s;counts",fChannel.at(i),fCryo.at(i),fString.at(i),fDetector.at(i),fStartRun,fEndRun,fEnergyName.c_str()));
    delete h;
  }
}

void GATAutoCal::FillHistoCal(TChain* mTree, vector<Double_t> Offset, vector<Double_t> Slope, Int_t Bin1, Double_t Low1, Double_t Up1, Int_t Bin2, Double_t Low2, Double_t Up2){
  
  string pos;
  Double_t slope=1;
  Double_t offset=0;
  for(Int_t i=0;i<fChannels;i++){
    TH1D *h;
    pos = Form("%d%d%d",fCryo.at(i),fString.at(i),fDetector.at(i));
    slope = Slope.at(i);
    offset = Offset.at(i);
    if(i%2 == 0){
      h= new TH1D(Form("%sCal%s%d",fEnergyName.c_str(),pos.c_str(),fChannel.at(i)),Form("Channel %d (C%dP%dD%d) in Run %d - %d",fChannel.at(i),fCryo.at(i),fString.at(i),fDetector.at(i),fStartRun,fEndRun), Bin1, Low1,Up1);
    }else{
      h= new TH1D(Form("%sCal%s%d",fEnergyName.c_str(),pos.c_str(),fChannel.at(i)),Form("Channel %d (C%dP%dD%d) in Run %d - %d",fChannel.at(i),fCryo.at(i),fString.at(i),fDetector.at(i),fStartRun,fEndRun), Bin2, Low2,Up2);
    }
    TCut cut1;
    if((fStartRun == 9913 && fEndRun == 9926) || (fStartRun == 13071 && fEndRun == 13074)){
      cut1 = Form("channel == %d && (rawWFMax > 200 || rawWFMin < -20)",fChannel.at(i));
    }else if(fStartRun == 60002368 && fEndRun == 60002372){
      cut1 = Form("channel == %d && (rawWFMax > 1000 || rawWFMin < -1000)",fChannel.at(i));
    }else if(fStartRun == 4171 && fEndRun == 4201){
      cut1 = Form("channel == %d",fChannel.at(i));
    }else{
      cut1 = Form("channel == %d",fChannel.at(i));
    }
    mTree->Draw(Form("%s*%f+%f>>%sCal%s%d",fEnergyName.c_str(),slope, offset,fEnergyName.c_str(),pos.c_str(),fChannel.at(i)),cut1);
    TrapE[i] = (TH1D*)h->Clone();
    TrapE[i]->SetTitle(Form("Channel %d (C%dP%dD%d) in Run %d - %d;%s;counts",fChannel.at(i),fCryo.at(i),fString.at(i),fDetector.at(i),fStartRun,fEndRun,fEnergyName.c_str()));
    delete h;
  }
}

void GATAutoCal::FillHistoDC(TChain* mTree, vector<Int_t> Channel, vector<string> EnergyName, Int_t Bin1, Double_t Low1, Double_t Up1, Int_t Bin2, Double_t Low2, Double_t Up2){

  if(Channel.size() != EnergyName.size()){
    cout << " channel size is not equal to energy name size!" << endl;
  }

  vector<string> Pos;
  string pos;
  for(size_t j=0;j<Channel.size();j++){
    for(Int_t i=0;i<fChannels;i++){
      if(Channel.at(j) == fChannel.at(i)){
	pos = Form("%d%d%d",fCryo.at(i),fString.at(i),fDetector.at(i));
	Pos.push_back(pos);
      }
    }
  }

  Int_t str;
  Int_t detpos;
  for(size_t i=0;i<Channel.size();i++){
    TH1D *h;
    str = fString.at(i);
    detpos = fDetector.at(i);

    if(Channel.at(i)%2 == 0){
      h= new TH1D(Form("%s%s%d","trapENF",Pos.at(i).c_str(),Channel.at(i)),Form("Channel %d (P%dD%d) in Run %d - %d",Channel.at(i),str,detpos,fStartRun,fEndRun), Bin1, Low1,Up1);
    }else{
      h= new TH1D(Form("%s%s%d","trapENF",Pos.at(i).c_str(),Channel.at(i)),Form("Channel %d (P%dD%d) in Run %d - %d",Channel.at(i),str,detpos,fStartRun,fEndRun), Bin2, Low2,Up2);
    }
    TCut cut1;
    if((fStartRun == 9913 && fEndRun == 9926) || (fStartRun == 13071 && fEndRun == 13074)){
      cut1 = Form("channel == %d && (rawWFMax > 200 || rawWFMin < -20)",Channel.at(i));
    }else{
      cut1 = Form("channel == %d",Channel.at(i));
    }
    mTree->Draw(Form("%s>>%s%s%d",EnergyName.at(i).c_str(),"trapENF",Pos.at(i).c_str(),Channel.at(i)),cut1);
    TrapE[i] = (TH1D*)h->Clone();

    TrapE[i]->SetTitle(Form("Channel %d (P%dD%d) in Run %d - %d;%s;counts",Channel.at(i),str,detpos,fStartRun,fEndRun,fEnergyName.c_str()));
    delete h;
  }
}


void GATAutoCal::FillHistoDCCal(TChain* mTree, vector<string> EnergyName,vector<Double_t> Offset, vector<Double_t> Slope, Int_t Bin1, Double_t Low1, Double_t Up1, Int_t Bin2, Double_t Low2, Double_t Up2){
  
  string pos;
  Double_t slope=1;
  Double_t offset=0;
  for(Int_t i=0;i<fChannels;i++){
    TH1D *h;
    pos = Form("%d%d%d",fCryo.at(i),fString.at(i),fDetector.at(i));
    slope = Slope.at(i);
    offset = Offset.at(i);
    if(i%2 == 0){
      h= new TH1D(Form("%sCal%s%d","trapENF",pos.c_str(),fChannel.at(i)),Form("Channel %d (C%dP%dD%d) in Run %d - %d",fChannel.at(i),fCryo.at(i),fString.at(i),fDetector.at(i),fStartRun,fEndRun), Bin1, Low1,Up1);
    }else{
      h= new TH1D(Form("%sCal%s%d","trapENF",pos.c_str(),fChannel.at(i)),Form("Channel %d (C%dP%dD%d) in Run %d - %d",fChannel.at(i),fCryo.at(i),fString.at(i),fDetector.at(i),fStartRun,fEndRun), Bin2, Low2,Up2);
    }
    TCut cut1;
    if((fStartRun == 9913 && fEndRun == 9926) || (fStartRun == 13071 && fEndRun == 13074)){
      cut1 = Form("channel == %d && (rawWFMax > 200 || rawWFMin < -20)",fChannel.at(i));
    }else if(fStartRun == 60002368 && fEndRun == 60002372){
      cut1 = Form("channel == %d && (rawWFMax > 1000 || rawWFMin < -1000)",fChannel.at(i));
    }else if(fStartRun == 4171 && fEndRun == 4201){
      cut1 = Form("channel == %d",fChannel.at(i));
    }else{
      cut1 = Form("channel == %d",fChannel.at(i));
    }
    mTree->Draw(Form("%s*%f+%f>>%sCal%s%d",EnergyName.at(i).c_str(),slope, offset,"trapENF",pos.c_str(),fChannel.at(i)),cut1);
    cout << pos.c_str() << " " << fChannel.at(i) << " " << EnergyName.at(i) << " " << slope << " " << offset << endl;
    TrapE[i] = (TH1D*)h->Clone();
    TrapE[i]->SetTitle(Form("Channel %d (C%dP%dD%d) in Run %d - %d;%s;counts",fChannel.at(i),fCryo.at(i),fString.at(i),fDetector.at(i),fStartRun,fEndRun,fEnergyName.c_str()));
    delete h;
  }
}

void GATAutoCal::SaveFile(string PathName,string FileName){
  TFile fhist(Form("%s%s",PathName.c_str(),FileName.c_str()),"update");
  for(Int_t i=0;i<fChannels;i++){
    TrapE[i]->Write();
  }
  cout << "Save histograms at: " << PathName.c_str() << FileName.c_str() << endl;
}


void GATAutoCal::LoadFile(string PathName,string FileName){
  cout << "Load histograms from: "<< FileName.c_str() << endl;
  TFile *fhist = TFile::Open(Form("%s%s",PathName.c_str(),FileName.c_str()),"read");
  string pos;
  for(int i=0;i<fChannels;i++){
    pos = Form("%d%d%d",fCryo.at(i),fString.at(i),fDetector.at(i));
    cout << i << " " << fEnergyName.c_str() << " " << pos.c_str() << " " << fChannel.at(i) << endl;
    TH1D *h = (TH1D*)fhist->Get(Form("%s%s%d",fEnergyName.c_str(),pos.c_str(),fChannel.at(i)));
    TrapE[i] = (TH1D*)h->Clone();
    delete h;
  }
}

Double_t GATAutoCal::GetMaximumPeak(TH1D *Hist,Double_t Low, Double_t Up){
  TH1D *hc = (TH1D*)Hist->Clone();
  Double_t maxPeak = 0;
  hc->GetXaxis()->SetRangeUser(Low,Up);
  Int_t maxbin = hc->GetMaximumBin();
  Double_t maxX = hc->GetBinCenter(maxbin);
  maxPeak = maxX;
  return maxPeak;
}

vector<Double_t> GATAutoCal::GetRefPeak(Double_t LastPeak, Double_t FirstPeak){
  vector<Double_t> refpeak;
  Double_t ref = 0;
  for(Int_t i = 0;i<fCalibrationPeaks;i++){
    ref = (fCalibrationPeak.at(i)-238.6)/(fCalibrationPeak.at(fCalibrationPeaks-1)-238.6)*(LastPeak-FirstPeak)+FirstPeak;
    refpeak.push_back(ref);
  }
  return refpeak;
}

void GATAutoCal::GaussFit(TH1D *Hist, Double_t Mean, Double_t Window,vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName,string FileName, string TitleName){

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  Par->clear();
  ParErr->clear();  

  string PlotName1 = Form("%sgaus_%s.pdf",PathName.c_str(),FileName.c_str()); 
  string PlotName2 = Form("%sgausresidual_%s.pdf",PathName.c_str(),FileName.c_str()); 

  TCanvas *c1 = new TCanvas("c1");

  TH1D *h1 = (TH1D*)Hist->Clone();
  h1->SetTitle(Form("%s;",TitleName.c_str()));
 
  
  if(h1->GetEntries()>0){
    //mean = h1->GetMean();
    //max = h1->GetMaximum();
    //width = h1->GetRMS();
    //maxbin = h1->GetMaximumBin();
    //maxX = h1->GetBinCenter(maxbin);

    TFitResultPtr fun = h1->Fit("gaus","S");

    for(Int_t i =0;i<3;i++){
      cout << fun->Parameter(i) << endl;
      Par->push_back(fun->Parameter(i));
      ParErr->push_back(fun->ParError(i));
    }
    
    
    Double_t sigma1 = abs(fun->Parameter(2));
    Double_t sigma2 = abs(fun->Parameter(2));
    Double_t fwhm = TMath::Sqrt(2*TMath::Log(2.))*(sigma1+sigma2);
    Double_t sigmaerr1 = fun->ParError(2);
    Double_t sigmaerr2 = fun->ParError(2);
    Double_t fwhmerr = TMath::Sqrt(pow(sigmaerr1,2)+pow(sigmaerr2,2));
    Double_t chi = fun->Chi2();
    Double_t ndf = fun->Ndf();
    Par->push_back(fwhm); //3
    ParErr->push_back(fwhmerr);//3
    Par->push_back(chi);//4
    Par->push_back(ndf);//5
    
    TF1 *fun1 = new TF1("fun1",Gaus,Par->at(1)-10*Par->at(2),Par->at(1)+10*Par->at(2),3);
    for(Int_t i=0;i<3;i++){
      fun1->SetParameter(i,Par->at(i));
    }
 

    if(FileName !="NA"){
      h1->GetXaxis()->SetRangeUser(Par->at(1)-5.*abs(Par->at(2)), Par->at(1)+5.*abs(Par->at(2)));
      h1->SetMinimum(0);
      h1->SetMaximum(abs(Par->at(0))*1.1);
      h1->GetYaxis()->SetRangeUser(0,abs(Par->at(0))*1.1);
      h1->Draw("E1");
      c1->Print(Form("%s",PlotName1.c_str()));

      TH1D *h2 = (TH1D*)h1->Clone();
      h2->SetTitle(Form("%s;Residual",TitleName.c_str()));
      h2->GetXaxis()->SetRangeUser(Mean-Window,Mean+Window);

      const Int_t h2entries = h2->GetSize()-2;
      h2->SetMarkerStyle(2);
      h2->SetMarkerColor(2);
      
      for(Int_t j =0; j<h2entries;j++){
	Double_t bincenter = h1->GetBinCenter(j);
	Double_t funvalue = fun1->Eval(bincenter);
	Double_t histovalue = h1->GetBinContent(j);
	Double_t residual = histovalue-funvalue;
	h2->SetBinContent(j,residual);
	h2->SetBinError(j,TMath::Sqrt(histovalue));
      }
      
      h2->SetMinimum(-20);
      h2->SetMaximum(20);
      h2->Draw("E1");
      c1->Print(Form("%s", PlotName2.c_str()));      
    }    
    /*
    for(Int_t i =0;i<6;i++){
      Par->push_back(0);
      ParErr->push_back(0);
    }
    */
  }else{
    for(Int_t i =0;i<6;i++){
      Par->push_back(0);
      ParErr->push_back(0);
    }
  }
}


void GATAutoCal::SkewGaussFit(TH1D *Hist, Double_t Mean, Double_t Window,vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName, string FileName, string TitleName, string opt){

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  Par->clear();
  ParErr->clear();  

  string PlotName1 = Form("%sgaus_%s.pdf",PathName.c_str(),FileName.c_str()); 
  string PlotName2 = Form("%sgausresidual_%s.pdf",PathName.c_str(),FileName.c_str()); 

  TCanvas *c1 = new TCanvas("c1");

  TH1D *h1 = (TH1D*)Hist->Clone();
  h1->SetTitle(Form("%s;",TitleName.c_str()));

  Double_t max, width, maxX;
  Int_t maxbin;

  if(h1->GetEntries()>0){
    max = h1->GetMaximum();
    width = h1->GetRMS();
    maxbin = h1->GetMaximumBin();
    maxX = h1->GetBinCenter(maxbin);
    //width = 0.2;
    TF1 *fun;
    fun = new TF1("fun",SkewGaus,Mean-Window,Mean+Window,4);
    fun->SetParameter(0,max);
    fun->SetParameter(1,maxX);
    fun->SetParameter(2,width);
    fun->SetParameter(3,width);
    //fun->SetParLimits(2,0,width*10);
    //fun->SetParLimits(3,0,width*10);
    //fun->SetParameter(4,0);
    //fun->SetParameter(5,0);
    //Int_t icolor = h1->GetLineColor();
    //fun->SetLineColor(icolor);
    h1->Fit("fun");
    
    for(Int_t i =0;i<4;i++){
      string temp = Form("%f",fun->GetParameter(i));
      string temperr = Form("%f",fun->GetParError(i));
      if(temp == "nan" || temp == "-nan" || temp == "inf" || temp == "-inf"){
	Par->push_back(0);
      }else{
	Par->push_back(fun->GetParameter(i));
      }  
      if(temperr == "nan" || temperr == "-nan" || temperr == "inf" || temperr == "-inf"){
        ParErr->push_back(0);
      }else{
	ParErr->push_back(fun->GetParError(i));
      }
    }
    Double_t sigma1 = abs(fun->GetParameter(2));
    Double_t sigma2 = abs(fun->GetParameter(3));
    Double_t fwhm = TMath::Sqrt(2*TMath::Log(2.))*(sigma1+sigma2);
    Double_t sigmaerr1 = fun->GetParError(2);
    Double_t sigmaerr2 = fun->GetParError(3);
    Double_t fwhmerr = TMath::Sqrt(pow(sigmaerr1,2)+pow(sigmaerr2,2));
    Double_t chi = fun->GetChisquare();
    Double_t ndf = fun->GetNDF();
    Par->push_back(fwhm); //4
    ParErr->push_back(fwhmerr);//4
    Par->push_back(chi);//5
    Par->push_back(ndf);//6

    if(FileName !="NA"){
      //h1->GetXaxis()->SetRangeUser(fun->GetParameter(1)-5.*abs(fun->GetParameter(2)), fun->GetParameter(1)+5.*abs(fun->GetParameter(3)));
      h1->SetMinimum(0);
      h1->SetMaximum(abs(fun->GetParameter(0))*1.1);
      //h1->GetYaxis()->SetRangeUser(0,abs(fun->GetParameter(0))*1.1);
      h1->GetXaxis()->SetRangeUser(Mean-Window,Mean+Window);
      h1->Draw(Form("E1 %s",opt.c_str()));
      c1->Print(Form("%s",PlotName1.c_str()));
      /*
      TH1D *h2 = (TH1D*)h1->Clone();
      h2->SetTitle(Form("%s;Residual",TitleName.c_str()));
      h2->GetXaxis()->SetRangeUser(Mean-Window,Mean+Window);
      const Int_t h2entries = h2->GetSize()-2;
      h2->SetMarkerStyle(2);
      h2->SetMarkerColor(2);
      
      for(Int_t j =0; j<h2entries;j++){
	Double_t bincenter = h1->GetBinCenter(j);
	Double_t funvalue = fun->Eval(bincenter);
	Double_t histovalue = h1->GetBinContent(j);
	Double_t residual = histovalue-funvalue;
	h2->SetBinContent(j,residual);
	h2->SetBinError(j,TMath::Sqrt(histovalue));
      }
      
      h2->SetMinimum(-20);
      h2->SetMaximum(20);
      h2->Draw("E1");
      c1->Print(Form("%s",PlotName2.c_str()));      
      */
    }    
  }else{
    for(Int_t i =0;i<7;i++){
      Par->push_back(0);
      ParErr->push_back(0);
    }
  }
  //  delete c1;
}

void GATAutoCal::GaussSlopeFit(TH1D *Hist, Double_t Mean, Double_t Window,vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName,string FileName, string TitleName){

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  Par->clear();
  ParErr->clear();  


  string PlotName1 = Form("%sgaus_%s",PathName.c_str(),FileName.c_str()); 
  string PlotName2 = Form("%sgausresidual_%s",PathName.c_str(),FileName.c_str()); 

  TCanvas *c1 = new TCanvas("c1");

  TH1D *h1 = (TH1D*)Hist->Clone();
  h1->SetTitle(Form("%s;",TitleName.c_str()));
  h1->GetXaxis()->SetRangeUser(Mean-Window,Mean+Window);

  Double_t max, width, maxX;
  Int_t maxbin;

  if(h1->GetEntries()>0){
    TH1D *htemp = (TH1D*)h1->Clone();    
    //mean = h1->GetMean();
    max = h1->GetMaximum();
    maxbin = h1->GetMaximumBin();
    maxX = h1->GetBinCenter(maxbin);
    htemp->GetXaxis()->SetRangeUser(maxX*0.995,maxX*1.005);
    width=htemp->GetRMS();
    delete htemp;

    TF1 *fun1 = new TF1("fun1","pol1",Mean-Window,Mean+Window);
    //TF1 *fun2 = new TF1("fun2","gaus",Mean-0.2*Window,Mean+0.2*Window);
    TF1 *fun3 = new TF1("fun3","gaus+pol1(3)",Mean-Window,Mean+Window);
    Double_t par[5];
    h1->Fit("fun1");
    par[3]=fun1->GetParameter(0);
    par[4]=fun1->GetParameter(1);
    /*
    h1->Fit("fun2");
    par[0]=fun2->GetParameter(0);
    par[1]=fun2->GetParameter(1);
    par[2]=fun2->GetParameter(2);
    */
    fun3->SetParameter(3,par[3]);
    fun3->SetParameter(4,par[4]);
    fun3->SetParameter(0,max);
    fun3->SetParameter(1,maxX);
    fun3->SetParameter(2,width);
    h1->Fit("fun3","V");

    for(Int_t i =0;i<5;i++){
      cout << fun3->GetParameter(i) << endl;
      Par->push_back(fun3->GetParameter(i));
      ParErr->push_back(fun3->GetParError(i));
    }    
    Double_t sigma1 = abs(fun3->GetParameter(2));
    Double_t sigma2 = abs(fun3->GetParameter(2));
    Double_t fwhm = TMath::Sqrt(2*TMath::Log(2.))*(sigma1+sigma2);
    Double_t sigmaerr1 = fun3->GetParError(2);
    Double_t sigmaerr2 = fun3->GetParError(2);
    Double_t fwhmerr = TMath::Sqrt(pow(sigmaerr1,2)+pow(sigmaerr2,2));
    Double_t chi = fun3->GetChisquare();
    Double_t ndf = fun3->GetNDF();
    Par->push_back(fwhm); //5
    ParErr->push_back(fwhmerr);//5
    Par->push_back(chi);//6
    Par->push_back(ndf);//7
    

    if(FileName !="NA"){
      h1->GetXaxis()->SetRangeUser(Par->at(1)-5.*abs(Par->at(2)), Par->at(1)+5.*abs(Par->at(2)));
      h1->Draw("E1");
      fun3->Draw("same");
      c1->Print(Form("%s",PlotName1.c_str()));

      TH1D *h2 = (TH1D*)h1->Clone();
      h2->SetTitle(Form("%s;Residual",TitleName.c_str()));
      h2->GetXaxis()->SetRangeUser(Mean-Window,Mean+Window);

      const Int_t h2entries = h2->GetSize()-2;
      h2->SetMarkerStyle(2);
      h2->SetMarkerColor(2);
      
      for(Int_t j =0; j<h2entries;j++){
	Double_t bincenter = h1->GetBinCenter(j);
	Double_t funvalue = fun3->Eval(bincenter);
	Double_t histovalue = h1->GetBinContent(j);
	Double_t residual = histovalue-funvalue;
	h2->SetBinContent(j,residual);
	h2->SetBinError(j,TMath::Sqrt(histovalue));
      }
      
      h2->SetMinimum(-20);
      h2->SetMaximum(20);
      h2->Draw("E1");
      c1->Print(Form("%s", PlotName2.c_str()));      
    }    
  }else{
    for(Int_t i =0;i<8;i++){
      Par->push_back(0);
      ParErr->push_back(0);
    }
  }
}


void GATAutoCal::EzFit(TH1D *Hist, Double_t Mean, Double_t Window,Int_t ReBin, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName,string FileName, string TitleName){
  gStyle->SetOptStat(0);
  //gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  Par->clear();
  ParErr->clear();
  
  string FitName = Form("%sfit_%s.pdf",PathName.c_str(),FileName.c_str());
  string PeakName = Form("%speak_%s.pdf",PathName.c_str(),FileName.c_str());
  string PeakName1 = Form("%speaklog_%s.pdf",PathName.c_str(),FileName.c_str());
  //string PeakName2 = Form("%speak_%s.root",PathName.c_str(),FileName.c_str());


  string ResidualName = Form("%sresidual_%s.pdf",PathName.c_str(),FileName.c_str());
  
  TH1D *h1 = (TH1D*)Hist->Clone();
  h1->SetTitle("");
  h1->GetXaxis()->SetTitle("Energy (keV)");
  h1->GetYaxis()->SetTitle("Counts");
  h1->GetYaxis()->SetTitleOffset(1.5);
  h1->GetXaxis()->SetTitleFont(132);
  h1->GetYaxis()->SetTitleFont(132);
  h1->GetXaxis()->SetLabelFont(132);
  h1->GetYaxis()->SetLabelFont(132);
  h1->Rebin(ReBin);
  //h1->SetTitle(Form("%s",TitleName.c_str()));

  h1->GetXaxis()->SetRangeUser(Mean-10*Window,Mean+10*Window);
  string parerr;
  GATPeakShape ps;
  ps.Reset();
  ps.EZFit(h1,Mean,0,true);
  for(Int_t i = 0;i<11;i++){
    string par = Form("%f",ps.GetPar(i));
    if(par == "nan" || par == "-nan"){
      Par->push_back(0);
    }else{
      Par->push_back(ps.GetPar(i));
    }

    parerr = Form("%f",ps.GetFitResult().Error(i));
    if(parerr == "nan" || parerr == "-nan"){
      ParErr->push_back(0);
    }else{
      ParErr->push_back(ps.GetFitResult().Error(i));
    }
  }
  Double_t max = h1->GetMaximum();
  Double_t funmax = 0;
  Double_t mu = ps.GetPar(1);
  Double_t funmu = mu;
  for(Int_t j = 0;j<100000;j++){
    Double_t bincenter = mu+0.001*(j-50000);
    Double_t funvalue = ps(&bincenter);
    if(funvalue > funmax){
      funmax = funvalue;
      funmu = bincenter;
    }
  }
  vector<Double_t> bintemp;
  vector<Double_t> difftemp;
  bintemp.clear();
  difftemp.clear();
  Double_t fwhm1 = mu;
  Double_t fwhm2 = mu;
  for(Int_t j = 0;j<100000;j++){
    Double_t bincenter = funmu-0.001*j;
    Double_t funvalue = ps(&bincenter);
    if(abs(funvalue-funmax/2)<5.){
      bintemp.push_back(bincenter);
      difftemp.push_back(abs(funvalue-funmax/2.));
    }
  }
  Double_t difftempmin = 1000;
  for(size_t i1 = 0;i1<bintemp.size();i1++){
    if(difftemp.at(i1)<difftempmin){
      fwhm1 = bintemp.at(i1);
      difftempmin = difftemp.at(i1);
    }
  }
  bintemp.clear();
  difftemp.clear();
  for(Int_t j = 0;j<100000;j++){
    Double_t bincenter = funmu+0.001*j;
    Double_t funvalue = ps(&bincenter);
    if(abs(funvalue-funmax/2)<5.){
      bintemp.push_back(bincenter);
      difftemp.push_back(abs(funvalue-funmax/2.));
    }
  }
  
  difftempmin = 1000;
  for(size_t i1 = 0;i1<bintemp.size();i1++){
    if(difftemp.at(i1)<difftempmin){
      fwhm2 = bintemp.at(i1);
      difftempmin = difftemp.at(i1);
    }
  }
  Double_t fwhm = fwhm2-fwhm1;

  //  max = h1->GetBinContent(funmaxbin-1);

  Par->push_back(ps.GetCentroid());//11
  Par->push_back(ps.GetFWHM()); //12
  Par->push_back(max);//13
  Par->push_back(ps.GetFitResult().MinFcnValue());//14
  Par->push_back(ps.GetFitResult().Ndf()); //15
  Par->push_back(ps.GetFitResult().Edm());//16
  Par->push_back(ps.GetStdev());//17
  Par->push_back(ps.GetEMax()); //18
  Par->push_back(funmax); // 19
  Par->push_back(ps.GetRangeMin()); // 20
  Par->push_back(ps.GetRangeMax()); // 21
  Par->push_back(ps.GetFitResult().IsValid());//22
  Par->push_back(fwhm); //23
  parerr = Form("%f",ps.GetCentroidUncertainty());
  if(parerr == "nan" || parerr == "-nan"){
    ParErr->push_back(0);
  }else{
    ParErr->push_back(ps.GetCentroidUncertainty());
  }

  parerr = Form("%f",ps.GetFWHMUncertainty());
  if(parerr == "nan" || parerr == "-nan"){
    ParErr->push_back(0);
  }else{
    ParErr->push_back(ps.GetFWHMUncertainty());
  }

  ParErr->push_back(0);//13
  ParErr->push_back(0);//14
  ParErr->push_back(0);//15
  ParErr->push_back(0);//16
  ParErr->push_back(0);//17
  ParErr->push_back(0);//18
  ParErr->push_back(0);//19
  ParErr->push_back(0);//20
  ParErr->push_back(0);//21
  ParErr->push_back(0);//22
  ParErr->push_back(0);//23

  if(FileName !="NA"){
    TCanvas *c1 = new TCanvas("c1");
    c1->SetLogy(0); 

    TH1D *h2 = (TH1D*)h1->Clone();
    h2->SetTitle(Form("%s;Residual",TitleName.c_str()));
    const Int_t h2entries = h2->GetSize()-2;
    h2->SetMarkerStyle(2);
    h2->SetMarkerColor(2);
    for(Int_t j =0; j<h2entries;j++){
      Double_t bincenter = h1->GetBinCenter(j);
      Double_t funvalue = ps(&bincenter);
      Double_t histovalue = h1->GetBinContent(j);
      Double_t residual = histovalue-funvalue;
      h2->SetBinContent(j,residual);
      h2->SetBinError(j,TMath::Sqrt(histovalue));
    }

    /*
    Int_t kAxisTitleFont = 133; // 13 (Times New Roman) + 3 (size in pixels)
    Float_t kAxisTitleSize = 30;
    Float_t kAxisTitleOffset = 1.0;
    h1->SetTitleSize(kAxisTitleSize, "XYZ");
    h1->SetTitleFont(kAxisTitleFont, "XYZ");
    h1->GetXaxis()->SetTitleOffset(kAxisTitleOffset);
    h1->GetYaxis()->SetTitleOffset(kAxisTitleOffset);
    h2->SetTitleSize(kAxisTitleSize, "XYZ");
    h2->SetTitleFont(kAxisTitleFont, "XYZ");
    h2->GetXaxis()->SetTitleOffset(kAxisTitleOffset);
    h2->GetYaxis()->SetTitleOffset(kAxisTitleOffset);
    */
    //h1->GetXaxis()->SetRange(ps.GetRangeMin(), ps.GetRangeMax());
    //h2->GetXaxis()->SetRange(ps.GetRangeMin(), ps.GetRangeMax());

    h1->GetXaxis()->SetRangeUser(Mean-Window,Mean+Window);
    ps.DrawComponents(h1);
    c1->Print(Form("%s",FitName.c_str()));
    ps.Draw(h1);
    c1->Print(Form("%s",PeakName.c_str()));
    c1->SetLogy();
    ps.Draw(h1);
    c1->Print(Form("%s",PeakName1.c_str()));

    /*
    c1->Print(Form("%s",PeakName1.c_str()));
    ps.Draw(h1);
    c1->Print(Form("%s",PeakName2.c_str()));
    */
    h2->Draw("E1");
    c1->Print(Form("%s",ResidualName.c_str()));
    delete c1;
    delete h2;
  }  
  delete h1;
}

void GATAutoCal::NoBGFit(TH1D *Hist, Double_t Mean, Double_t Window,Int_t ReBin, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName,string FileName, string TitleName){

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  Par->clear();
  ParErr->clear();
  
  string FitName = Form("%sfit_%s.pdf",PathName.c_str(),FileName.c_str());
  string PeakName = Form("%speak_%s.pdf",PathName.c_str(),FileName.c_str());
  string ResidualName = Form("%sresidual_%s.pdf",PathName.c_str(),FileName.c_str());
  TH1D *h1 = (TH1D*)Hist->Clone();
  h1->Rebin(ReBin);
  h1->SetTitle(Form("%s",TitleName.c_str()));  
  h1->GetXaxis()->SetRangeUser(Mean-3*Window,Mean+3*Window);
  Double_t max = h1->GetMaximum();

  string parerr;
  GATPeakShape ps;
  //ps.SetRange(Mean-3*Window,Mean+3*Window);

  //for(size_t i =0;i<IniPar->size();i++){
  //  ps.SetParameter(i,IniPar->at(i));
  // }

  ps.EZFit(h1,Mean,0,true);
  ps.UseBG(false,false,false);
  ps.FitTo(h1);

  for(Int_t i = 0;i<11;i++){
    string par = Form("%f",ps.GetPar(i));
    if(par == "nan" || par == "-nan"){
      Par->push_back(0);
    }else{
      Par->push_back(ps.GetPar(i));
    }

    parerr = Form("%f",ps.GetFitResult().Error(i));
    if(parerr == "nan" || parerr == "-nan"){
      ParErr->push_back(0);
    }else{
      ParErr->push_back(ps.GetFitResult().Error(i));
    }
  }
  
  TH1D *h2 = (TH1D*)h1->Clone();
  h2->SetTitle(Form("%s;Residual",TitleName.c_str()));
  const Int_t h2entries = h2->GetSize()-2;
  h2->SetMarkerStyle(2);
  h2->SetMarkerColor(2);
  for(Int_t j =0; j<h2entries;j++){
    Double_t bincenter = h1->GetBinCenter(j);
    Double_t funvalue = ps(&bincenter);
    Double_t histovalue = h1->GetBinContent(j);
    Double_t residual = histovalue-funvalue;
    h2->SetBinContent(j,residual);
    h2->SetBinError(j,TMath::Sqrt(histovalue));
  }
  /*
  Int_t centroidbin=0;
  for(Int_t j =0; j<h1->GetSize()-2;j++){
    Double_t bincenter = h1->GetBinCenter(j);
    Double_t histovalue = h1->GetBinContent(j);
    if(abs(bincenter-Par->at(1))<h1->GetBinWidth(j)){
      centroidbin = j;
      max = histovalue;
    }
  }
  */
  Double_t funmax = 0;
  Double_t mu = ps.GetPar(1);
  for(Int_t j = 0;j<10000;j++){
    Double_t bincenter = mu+0.01*(j-5000);
    Double_t funvalue = ps(&bincenter);
    if(funvalue > funmax){
      funmax = funvalue;
    }
  }

  vector<Double_t> bintemp;
  vector<Double_t> difftemp;
  bintemp.clear();
  difftemp.clear();
  Double_t fwhm1 = mu;
  Double_t fwhm2 = mu;
  for(Int_t j = 0;j<10000;j++){
    Double_t bincenter = mu-0.01*j;
    Double_t funvalue = ps(&bincenter);
    if(abs(funvalue-funmax/2)<1.){
      bintemp.push_back(bincenter);
      difftemp.push_back(abs(funvalue-funmax/2.));
    }
  }
  Double_t difftempmin = 1000;
  for(size_t i1 = 0;i1<bintemp.size();i1++){
    if(difftemp.at(i1)<difftempmin){
      fwhm1 = bintemp.at(i1);
      difftempmin = difftemp.at(i1);
    }
  }
  bintemp.clear();
  difftemp.clear();
  for(Int_t j = 0;j<10000;j++){
    Double_t bincenter = mu+0.01*j;
    Double_t funvalue = ps(&bincenter);
    if(abs(funvalue-funmax/2)<1.){
      bintemp.push_back(bincenter);
      difftemp.push_back(abs(funvalue-funmax/2.));
    }
  }
  difftempmin = 1000;
  for(size_t i1 = 0;i1<bintemp.size();i1++){
    if(difftemp.at(i1)<difftempmin){
      fwhm2 = bintemp.at(i1);
      difftempmin = difftemp.at(i1);
    }
  }
  Double_t fwhm = fwhm2-fwhm1;


  Par->push_back(ps.GetCentroid());//11
  Par->push_back(ps.GetFWHM()); //12
  Par->push_back(max);//13
  Par->push_back(ps.GetFitResult().MinFcnValue());//14
  Par->push_back(ps.GetFitResult().Ndf()); //15
  Par->push_back(ps.GetFitResult().Edm());//16
  Par->push_back(ps.GetStdev());//17
  Par->push_back(ps.GetEMax()); //18
  Par->push_back(funmax); // 19
  Par->push_back(ps.GetRangeMin()); // 20
  Par->push_back(ps.GetRangeMax()); // 21
  Par->push_back(ps.GetFitResult().IsValid());//22
  Par->push_back(fwhm);//23

  parerr = Form("%f",ps.GetCentroidUncertainty());
  if(parerr == "nan" || parerr == "-nan"){
    ParErr->push_back(0);
  }else{
    ParErr->push_back(ps.GetCentroidUncertainty());
  }

  parerr = Form("%f",ps.GetFWHMUncertainty());
  if(parerr == "nan" || parerr == "-nan"){
    ParErr->push_back(0);
  }else{
    ParErr->push_back(ps.GetFWHMUncertainty());
  }

  ParErr->push_back(0);//13
  ParErr->push_back(0);//14
  ParErr->push_back(0);//15
  ParErr->push_back(0);//16
  ParErr->push_back(0);//17
  ParErr->push_back(0);//18
  ParErr->push_back(0);//19
  ParErr->push_back(0);//20
  ParErr->push_back(0);//21
  ParErr->push_back(0);//22
  ParErr->push_back(0);//23
  if(FileName !="NA"){
    TCanvas *c1 = new TCanvas("c1");
    c1->SetLogy(0); 

    TH1D *h2 = (TH1D*)h1->Clone();
    h2->SetTitle(Form("%s;Residual",TitleName.c_str()));
    const Int_t h2entries = h2->GetSize()-2;
    h2->SetMarkerStyle(2);
    h2->SetMarkerColor(2);
    for(Int_t j =0; j<h2entries;j++){
      Double_t bincenter = h1->GetBinCenter(j);
      Double_t funvalue = ps(&bincenter);
      Double_t histovalue = h1->GetBinContent(j);
      Double_t residual = histovalue-funvalue;
      h2->SetBinContent(j,residual);
      h2->SetBinError(j,TMath::Sqrt(histovalue));
    }

    /*
    Int_t kAxisTitleFont = 133; // 13 (Times New Roman) + 3 (size in pixels)
    Float_t kAxisTitleSize = 30;
    Float_t kAxisTitleOffset = 1.0;
    h1->SetTitleSize(kAxisTitleSize, "XYZ");
    h1->SetTitleFont(kAxisTitleFont, "XYZ");
    h1->GetXaxis()->SetTitleOffset(kAxisTitleOffset);
    h1->GetYaxis()->SetTitleOffset(kAxisTitleOffset);
    h2->SetTitleSize(kAxisTitleSize, "XYZ");
    h2->SetTitleFont(kAxisTitleFont, "XYZ");
    h2->GetXaxis()->SetTitleOffset(kAxisTitleOffset);
    h2->GetYaxis()->SetTitleOffset(kAxisTitleOffset);
    */
    //h1->GetXaxis()->SetRange(ps.GetRangeMin(), ps.GetRangeMax());
    //h2->GetXaxis()->SetRange(ps.GetRangeMin(), ps.GetRangeMax());

    ps.DrawComponents(h1);
    c1->Print(Form("%s",FitName.c_str()));
    ps.Draw(h1);
    c1->Print(Form("%s",PeakName.c_str()));
    h2->Draw("E1");
    c1->Print(Form("%s",ResidualName.c_str()));
    delete c1;
    delete h2;
  }  
  delete h1;
}

void GATAutoCal::MultiPeakFitterFull(TH1F *Hist, Double_t scaleenergy, Double_t scaleamps, string FileName, vector<string>* ParName, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov){
  ParName->clear();
  Par->clear();
  ParErr->clear();
  Cov->clear();
  gStyle->SetOptStat(0);
  TH1F* hist = (TH1F*)Hist->Clone();

  if(hist==NULL) {
    cout <<  "Histogram not found!" << endl;
    gSystem->Abort();
  }

  hist = (TH1F*) hist->Clone();
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle("Energy (keV)");
  hist->GetYaxis()->SetTitle("Counts");
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetTitleFont(132);
  hist->GetYaxis()->SetTitleFont(132);
  hist->GetXaxis()->SetLabelFont(132);
  hist->GetYaxis()->SetLabelFont(132);
  
  //Start setting up the fitter

  GATMultiPeakFitter fitter;

  //Add the regions. This is AddRegion(low bin, high bin, list of peak energies in keV)

  fitter.AddRegion(110, 120, {115.183});//1 peak
  fitter.AddRegion(200, 340, {215.983, 233.36, 238.632, 240.986, 252.61, 277.371, 288.2, 300.087, 328.03});//9 peaks
  fitter.AddRegion(420, 490, {452.98}); //1 peak
  fitter.AddRegion(530, 620, {549.73, 583.191}); //2 peaks
  fitter.AddRegion(716, 815, {722.04, 727.330, 763.130, 785.370});//4 peaks
  fitter.AddRegion(850, 1100, {860.557, 893.408, 952.12, 1078.62, 1093.9}); //5 peaks
  fitter.AddRegion(1450, 1570, {1512.7}); //1 peak
  fitter.AddRegion(1600, 1675, {1620.738}); //1 peak
  fitter.AddRegion(2550, 2800, {2614.533}); //1 peak
  //fitter.AddRegion(-10,10,{0}); // 1 Peak
  vector<Double_t> Low;
  Low.push_back(110);
  Low.push_back(200);
  Low.push_back(420);
  Low.push_back(530);
  Low.push_back(716);
  Low.push_back(850);
  Low.push_back(1450);
  Low.push_back(1600);
  Low.push_back(2550);
  //Low.push_back(-10);
  vector<Double_t> Up;
  Up.push_back(120);
  Up.push_back(340);
  Up.push_back(490);
  Up.push_back(620);
  Up.push_back(815);
  Up.push_back(1100);
  Up.push_back(1570);
  Up.push_back(1675);
  Up.push_back(2800);

  //Initial amplitudes for each peak. Limit these to >0
  
  fitter.SetParFunction(GATMultiPeakFitter::kAmp, GATMultiPeakFitter::kFree, {/*region1*/ 500, /*region 2*/ 1000, 200, 160000, 12000, 10000, 8000, 1000, 10000, 500, /*region4*/ 1100, /*region5*/ 400, 80000, /*region6*/ 200, 15000, 1000, 2500, /*region7*/ 10000, 800, 250, 1000, 300, /*region8*/ 500, /*region9*/ 2500, /*region10*/ 35000}); //26 total pars

  for(size_t i=0; i<fitter.NPeaks(); i++) fitter.LimitPar(GATMultiPeakFitter::kAmp, i, 0., (double) NAN);
  fitter.SetParFunction(GATMultiPeakFitter::kMu, GATMultiPeakFitter::kLinear, {0, 1.});
  fitter.SetParFunction(GATMultiPeakFitter::kSig, GATMultiPeakFitter::kRootQuad, {0.2, 0.017, 0.0002});
  fitter.SetParFunction(GATMultiPeakFitter::kFt, GATMultiPeakFitter::kConst, {0.3});
  fitter.SetParFunction(GATMultiPeakFitter::kTau, GATMultiPeakFitter::kLinear, {0., 0.0006});
  fitter.LimitPar(GATMultiPeakFitter::kFt, 0, 0., 1.);
  fitter.LimitPar(GATMultiPeakFitter::kTau, 1, 0.00001, 0.01);

  //uncomment the following 3 lines to turn off the LE tail
  //fitter.SetParFunction(GATMultiPeakFitter::kTau, GATMultiPeakFitter::kConst, {0.5});
  //fitter.FixPar(GATMultiPeakFitter::kFt, 0, 0);
  //fitter.FixPar(GATMultiPeakFitter::kTau, 0, 0.5);
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

  // initial BG parameters
  fitter.SetBGPars(0, 5000, 2, 0);
  fitter.FixBGPar(0, 2, 0);
  fitter.SetBGPars(1, 3000, -18, 0);
  fitter.SetBGPars(2, 1300, -6, 0);
  fitter.SetBGPars(3, 500, 0, 0);
  fitter.SetBGPars(4, 250, -0.2, 0);
  fitter.SetBGPars(5, 200, -0.2, 0);
  fitter.SetBGPars(6, 140, 0, 0);
  fitter.SetBGPars(7, 170, 0, 0);
  fitter.SetBGPars(8, 4, 0, 0);
  fitter.FixBGPar(8, 1, 0);
  fitter.FixBGPar(8, 2, 0);
  fitter.ScaleEnergy(scaleenergy);
  fitter.ScaleAmps(scaleamps); 
  // Set the histogram to be used for each region. This will set all histograms to be the same one
  fitter.SetHists(hist);
  // Set the histogram to be used for the 0E peak. This will not change any of the others. The first parameter is the index of the region
  //fitter.SetHist(0, hist0E);

  // Set up a hybrid montecarlo fit to improve the initial parameters
  GATHybridMonteCarlo hmc;
  hmc.SetNLLFunc(dynamic_pointer_cast<ROOT::Math::IGradientFunctionMultiDim>(fitter.SetPoissonLLFCN()));
  // uncomment the following lines to generate a TTree with each MCMC step
  hmc.SetOutputFile(Form("hmc_%s.root",FileName.c_str()));
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

  ofstream fout(Form("./List/parameters_%s.txt",FileName.c_str()),ios::app);
  fout.precision(15);
  TCanvas *c1 = new TCanvas("c1");

  //perform a minuit fit using the most likely parameters found during the HMC as the initial parameters. Use minos error estimation
  fitter.SetParameters(hmc.GetLikeliestPars().data());
  for(size_t i=0; i<fitter.NPeaks(); i++) fitter.LimitPar(GATMultiPeakFitter::kAmp, i, 10., (double) NAN);  
  fitter.FitToHists(0., true, false, true);
  fout << FileName.c_str() << " isValid: " << fitter.GetResult().IsValid() << endl;
  cout << "fitting is valid :" << fitter.GetResult().IsValid() << endl;
 
  for(size_t  i = 0;i<fitter.NPar();i++){
    ParName->push_back(fitter.ParameterName(i));
    Par->push_back(fitter.GetResult().Parameter(i));
    ParErr->push_back(fitter.GetResult().Error(i));
  }

  fout << fitter.GetResult().Parameter(25) << " " << fitter.GetResult().Error(25)<< " " << fitter.GetResult().Parameter(26) << " " << fitter.GetResult().Error(26) << endl;
  fout << fitter.GetResult().Parameter(27) << " " << fitter.GetResult().Error(27) << " " 
       << fitter.GetResult().Parameter(28) << " " << fitter.GetResult().Error(28) << " "
       << fitter.GetResult().Parameter(29) << " " << fitter.GetResult().Error(29) << endl;
  
  Double_t cov = 0;
  for(unsigned int k = 25;k<27;k++){
    for(unsigned int j = 25;j<27;j++){
      cov = fitter.GetResult().CovMatrix(k,j);
      if(k==j){
        fout << k << " " << j << " " << abs(cov) << endl;
      }else{
        fout << k << " " << j << " " << cov << endl;
      }
    }
  }

  cov = 0;
  for(unsigned int k = 27;k<30;k++){
    for(unsigned int j = 27;j<30;j++){
      cov = fitter.GetResult().CovMatrix(k,j);
      if(k==j){
        fout << k << " " << j << " " << abs(cov) << endl;
      }else{
	fout << k << " " << j << " " << cov << endl;
      }
    }
  }
  cov = 0;

  for(unsigned int k = 0;k<(unsigned int)fitter.NPar();k++){
    for(unsigned int j = 0;j<(unsigned int)fitter.NPar();j++){
      cov = fitter.GetResult().CovMatrix(k,j);
      Cov->push_back(cov);
    }
  }
  
  c1->SetLogy();
  fitter.DrawRegions();
  c1->Print(Form("./spectrum_%s.pdf",FileName.c_str()));
  c1->Print(Form("./spectrum_%s.png",FileName.c_str()));
  for(size_t i=0; i<fitter.NRegions(); i++){
    hist->GetXaxis()->SetRangeUser(Low.at(i)*scaleenergy, Up.at(i)*scaleenergy);
    cout << i << " " << Low.at(i) << " " << Up.at(i) << endl;
    fitter.DrawRegion(i);
    c1->Print(Form("./spectrum_%s_%d.png",FileName.c_str(),(Int_t)i));
    c1->Print(Form("./spectrum_%s_%d.pdf",FileName.c_str(),(Int_t)i));
  }
  c1->SetLogy(0);
  hist->GetXaxis()->SetRangeUser(2600,2625);
  fitter.DrawRegion(8);
  //c1->Print(Form("./spectrum_%s_2614.root",FileName.c_str()));
  c1->Print(Form("./spectrum_%s_2614.pdf",FileName.c_str()));
  c1->SetLogy();
  fitter.DrawRegion(8);
  //c1->Print(Form("./spectrum_%s_2614.root",FileName.c_str()));
  c1->Print(Form("./spectrumlog_%s_2614.pdf",FileName.c_str()));


  Double_t enr1=2039;
  Double_t fwhm1=0;
  Double_t fwhmunc1=0;
  Double_t enr2=2614;
  Double_t fwhm2;
  Double_t fwhmunc2;
  Double_t enr3=60;
  Double_t fwhm3;
  Double_t fwhmunc3;

  fitter.GetFWHM(enr1, fwhm1,  fwhmunc1);
  fout << FileName.c_str() << " " << enr1  << " " <<  fwhm1 << " " << fwhmunc1 << endl;
  fitter.GetFWHM(enr2,fwhm2,fwhmunc2);
  fout << FileName.c_str() << " " << enr2  << " " <<  fwhm2 << " " << fwhmunc3 << endl;
  fitter.GetFWHM(enr3,fwhm3,fwhmunc3);
  fout << FileName.c_str() << " " << enr3  << " " <<  fwhm3 << " " << fwhmunc3 << endl;

  vector<Double_t> calpeak = fitter.GetPeaks();
  for(size_t i = 0;i<calpeak.size();i++){
    ParName->push_back(Form("#Peak_%d",(Int_t)i));
    Par->push_back(calpeak.at(i));
    ParErr->push_back(0);
  }

  for(size_t i = 0;i<calpeak.size();i++){
    double ip;
    double iperr;
    fitter.GetCentroidOfPeak(i,ip,iperr);
    ParName->push_back(Form("#Mu_%d",(Int_t)i));
    Par->push_back(ip);
    ParErr->push_back(iperr);
  }
  //fitter.SetParFunction(GATMultiPeakFitter::kMu, GATMultiPeakFitter::kLinear, {0, 1.});

  fitter.SetParFunction(GATMultiPeakFitter::kMu, GATMultiPeakFitter::kFree);
  fitter.FitToHists(0., true, false, true);
  for(size_t i = 0;i<calpeak.size();i++){
    double ip;
    double iperr;
    fitter.GetCentroidOfPeak(i,ip,iperr);
    ParName->push_back(Form("#Mu_%d",(Int_t)i));
    Par->push_back(ip);
    ParErr->push_back(iperr);
  }
  

  /*
  for(size_t i = 0;i<Par->size();i++){
    fout << i << " " << Par->at(i) << " " << ParErr->at(i) << endl;
  }
  */
}

Int_t GATAutoCal::MultiPeakFitter(TH1F *Hist, Double_t scaleenergy, Double_t scaleamps,string FileName, vector<string>* ParName, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov){
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
  //  fitter.AddRegion(-10,10, {0});
  fitter.AddRegion(200, 340, {238.632, 240.986, 277.371, 300.087});//4 peaks
  fitter.AddRegion(540, 620, {583.191}); //1 peaks
  fitter.AddRegion(716, 900, {727.330, 785.37, 860.557});//3 peaks
  //fitter.AddRegion(716, 900, {727.330, 860.557});//2 peaks
  fitter.AddRegion(2550, 2800, {2614.533}); //1 peak
 
  // Start setting initial parameters
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
  hmc.SetOutputFile(Form("hmc_%s.root",FileName.c_str()));
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
  ofstream fout(Form("./List/parameters_%s.txt",FileName.c_str()),ios::app);
  fout.precision(15);
  TCanvas *c1 = new TCanvas("c1");

  fitter.SetParameters(hmc.GetLikeliestPars().data());  
  for(size_t i=0; i<fitter.NPeaks(); i++) fitter.LimitPar(GATMultiPeakFitter::kAmp, i, 10., (double) NAN);

  fitter.FitToHists(0., true, false, true);
  fout << FileName.c_str() << " isValid: " << fitter.GetResult().IsValid() << endl;
  cout << "fitting is valid :" << fitter.GetResult().IsValid() << endl;
  Int_t isValid = fitter.GetResult().IsValid();
  for(size_t  i = 0;i<fitter.NPar();i++){
    ParName->push_back(fitter.ParameterName(i));
    Par->push_back(fitter.GetResult().Parameter(i));
    ParErr->push_back(fitter.GetResult().Error(i));
  }
  /*
  Double_t fOffset;
  Double_t fSlope;
  fitter.GetCalibrationPars(fOffset, fSlope);
  ParName->push_back("Offset");
  ParName->push_back("Slope");
  Par->push_back(fOffset);
  Par->push_back(fSlope);
  ParErr->push_back(0);
  ParErr->push_back(0);
  */
  fout << fitter.GetResult().Parameter(9) << " " << fitter.GetResult().Error(9)<< " " << fitter.GetResult().Parameter(10) << " " << fitter.GetResult().Error(10) << endl;
  fout << fitter.GetResult().Parameter(11) << " " << fitter.GetResult().Error(11) << " " 
       << fitter.GetResult().Parameter(12) << " " << fitter.GetResult().Error(12) << " "
       << fitter.GetResult().Parameter(13) << " " << fitter.GetResult().Error(13) << endl;

  //TMatrixDSym cov;
  //fitter.GetResult().GetCovarianceMatrix(&cov);
  Double_t cov = 0;
  for(unsigned int k = 9;k<11;k++){
    for(unsigned int j = 9;j<11;j++){
      cov = fitter.GetResult().CovMatrix(k,j);
      if(k==j){
        fout << k << " " << j << " " << abs(cov) << endl;
      }else{
        fout << k << " " << j << " " << cov << endl;
      }
    }
  }

  cov = 0;
  for(unsigned int k = 11;k<14;k++){
    for(unsigned int j = 11;j<14;j++){
      cov = fitter.GetResult().CovMatrix(k,j);
      if(k==j){
        fout << k << " " << j << " " << abs(cov) << endl;
      }else{
	fout << k << " " << j << " " << cov << endl;
      }
    }
  }
  cov = 0;

  for(unsigned int k = 0;k<(unsigned int)fitter.NPar();k++){
    for(unsigned int j = 0;j<(unsigned int)fitter.NPar();j++){
      cov = fitter.GetResult().CovMatrix(k,j);
      fout << k*fitter.NPar()+j <<  " " << k << " " << j  << " " << cov << endl;
      Cov->push_back(cov);
    }
  }
  
  for(size_t i=0;i<9;i++){
    if(Par->at(i)<0){
      isValid = 0;
    }
  }
  /*
  ParName->push_back("isValid");
  Par->push_back(isValid);
  ParErr->push_back(0);
  */
  c1->SetLogy();
  for(size_t i=0; i<fitter.NRegions(); i++){
    hist->GetXaxis()->SetRangeUser(Low.at(i)*scaleenergy, Up.at(i)*scaleenergy);
    fitter.DrawRegion(i);
    c1->Print(Form("spectrum_%s_%d.pdf",FileName.c_str(),(Int_t)i));
  }
  return isValid;
}

Int_t GATAutoCal::MultiPeakFitter(TH1F *Hist, Double_t scaleenergy, Double_t scaleamps,string FileName, vector<string>* ParName, vector<Double_t>* Par, vector<Double_t>* ParErr, vector<Double_t>* Cov, TH1F *HistoE0){
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
  // Start setting initial parameters
  vector<Double_t> Low;
  vector<Double_t> Up;
  Low.push_back(200);
  Low.push_back(540);
  Low.push_back(716);
  Low.push_back(2550);
  Low.push_back(-3);
  
  Up.push_back(340);
  Up.push_back(620);
  Up.push_back(900);
  Up.push_back(2800);
  Up.push_back(3);
  // initial BG parameters
  //  fitter.AddRegion(-10,10, {0});

  fitter.AddRegion(Low.at(0), Up.at(0), {238.632, 240.986, 277.371, 300.087});//4 peaks
  fitter.AddRegion(Low.at(1), Up.at(1), {583.191}); //1 peaks
  fitter.AddRegion(Low.at(2), Up.at(2), {727.330, 785.37, 860.557});//3 peaks
  //fitter.AddRegion(716, 900, {727.330, 860.557});//2 peaks
  fitter.AddRegion(Low.at(3), Up.at(3), {2614.533}); //1 peak
  fitter.AddRegion(Low.at(4), Up.at(4),{0.}); // 1 peak

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
  //fitter.SetBGPars(4,0,0,0);
  //fitter.FixBGPar(4,1,0);
  //fitter.FixBGPar(4,2,0);

  // Start setting initial parameters  
  //SetParFunction(parameter, function type, list of parameters)

  
  //Initial amplitudes for each peak. Limit these to >0
  fitter.SetParFunction(GATMultiPeakFitter::kAmp, GATMultiPeakFitter::kFree, {/*region 1*/ 160000, 12000, 8000, 10000, /*region2*/ 80000, /*region3*/ 15000, 3000, 10000, /*region4*/ 35000,/*region 5*/ 35000}); //25 total pars

  for(size_t i=0; i<fitter.NPeaks(); i++) fitter.LimitPar(GATMultiPeakFitter::kAmp, i, 1., (double) NAN);
  fitter.SetParFunction(GATMultiPeakFitter::kMu, GATMultiPeakFitter::kLinear, {0., 1.});
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

  Int_t maxbin = HistoE0->GetMaximumBin();
  Double_t E0Int = HistoE0->Integral()*HistoE0->GetBinWidth(maxbin);
  cout << "How many peaks:" << fitter.NPeaks() << " Maximum Bin:" << maxbin << " " << HistoE0->Integral() << " " << endl; 
  //  fitter.SetPar(GATMultiPeakFitter::kAmp, fitter.NPeaks()-1, HistoE0->Integral()*HistoE0->GetBinWidth(maxbin));
  fitter.SetHist(fitter.NRegions()-1, HistoE0);
  fitter.SetPar(GATMultiPeakFitter::kAmp, fitter.NPeaks()-1, E0Int);

  // Set up a hybrid montecarlo fit to improve the initial parameters

  GATHybridMonteCarlo hmc;
  hmc.SetNLLFunc(dynamic_pointer_cast<ROOT::Math::IGradientFunctionMultiDim>(fitter.SetPoissonLLFCN()));
  // uncomment the following lines to generate a TTree with each MCMC step
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
  ofstream fout(Form("./List/parameters_%s.txt",FileName.c_str()),ios::app);
  fout.precision(15);
  TCanvas *c1 = new TCanvas("c1");

  fitter.SetParameters(hmc.GetLikeliestPars().data());  
  fitter.FitToHists(0., true, false, true);
  fout << FileName.c_str() << " isValid: " << fitter.GetResult().IsValid() << endl;
  cout << "fitting is valid :" << fitter.GetResult().IsValid() << endl;
  Int_t isValid = fitter.GetResult().IsValid();
  for(size_t  i = 0;i<fitter.NPar();i++){
    ParName->push_back(fitter.ParameterName(i));
    Par->push_back(fitter.GetResult().Parameter(i));
    ParErr->push_back(fitter.GetResult().Error(i));
  }
  /*
  Double_t fOffset;
  Double_t fSlope;
  fitter.GetCalibrationPars(fOffset, fSlope);
  ParName->push_back("Offset");
  ParName->push_back("Slope");
  Par->push_back(fOffset);
  Par->push_back(fSlope);
  ParErr->push_back(0);
  ParErr->push_back(0);
  */
  fout << fitter.GetResult().Parameter(9) << " " << fitter.GetResult().Error(9)<< " " << fitter.GetResult().Parameter(10) << " " << fitter.GetResult().Error(10) << endl;
  fout << fitter.GetResult().Parameter(11) << " " << fitter.GetResult().Error(11) << " " 
       << fitter.GetResult().Parameter(12) << " " << fitter.GetResult().Error(12) << " "
       << fitter.GetResult().Parameter(13) << " " << fitter.GetResult().Error(13) << endl;

  //TMatrixDSym cov;
  //fitter.GetResult().GetCovarianceMatrix(&cov);
  Double_t cov = 0;
  for(unsigned int k = 9;k<11;k++){
    for(unsigned int j = 9;j<11;j++){
      cov = fitter.GetResult().CovMatrix(k,j);
      if(k==j){
        fout << k << " " << j << " " << abs(cov) << endl;
      }else{
        fout << k << " " << j << " " << cov << endl;
      }
    }
  }

  cov = 0;
  for(unsigned int k = 11;k<14;k++){
    for(unsigned int j = 11;j<14;j++){
      cov = fitter.GetResult().CovMatrix(k,j);
      if(k==j){
        fout << k << " " << j << " " << abs(cov) << endl;
      }else{
	fout << k << " " << j << " " << cov << endl;
      }
    }
  }
  cov = 0;

  for(unsigned int k = 0;k<(unsigned int)fitter.NPar();k++){
    for(unsigned int j = 0;j<(unsigned int)fitter.NPar();j++){
      cov = fitter.GetResult().CovMatrix(k,j);
      fout << k*fitter.NPar()+j <<  " " << k << " " << j  << " " << cov << endl;
      Cov->push_back(cov);
    }
  }
  for(size_t i=0;i<9;i++){
    if(Par->at(i)<0){
      isValid = 0;
    }
  }
  /*
  ParName->push_back("isValid");
  Par->push_back(isValid);
  ParErr->push_back(0);
  */
  c1->SetLogy();
  for(size_t i=0; i<fitter.NRegions(); i++){
    hist->GetXaxis()->SetRangeUser(Low.at(i)*scaleenergy, Up.at(i)*scaleenergy);
    fitter.DrawRegion(i);
    c1->Print(Form("spectrum_%s_%d.pdf",FileName.c_str(),(Int_t)i));
  }
  return isValid;
}

void GATAutoCal::LinearFit(vector<Double_t> Px, vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName,string FileName,string TitleName, Double_t Energy){

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");

  string CalibrationName=Form("%scalibration_%s.pdf",PathName.c_str(),FileName.c_str());
  string ResidualName=Form("%scalibrationresidual_%s.pdf",PathName.c_str(),FileName.c_str());
  string DeltaName=Form("%scalibrationdelta_%s.pdf",PathName.c_str(),FileName.c_str());
  ofstream ferror(Form("linearfiterr_%s.txt",FileName.c_str()),ios::app);
  //ofstream fdelta(Form("delta_%s.txt",FileName.c_str()),ios::app);

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

  TGraphErrors fDeltaGraph;
  fDeltaGraph.Set(0);
  fDeltaGraph.SetTitle(Form("%s;Energy (keV);#Delta E (keV)",TitleName.c_str()));
  fDeltaGraph.SetMarkerStyle(8);
  fDeltaGraph.SetMarkerSize(2);
  fDeltaGraph.SetMarkerColor(4);
  
  TGraphErrors fResidual;
  fResidual.Set(0);
  fResidual.SetTitle(Form("%s;%s;Residual",TitleName.c_str(),fEnergyName.c_str()));
  fResidual.SetMarkerStyle(8);
  fResidual.SetMarkerSize(2);
  fResidual.SetMarkerColor(4);
  
  for(Int_t i = 0;i<nCal;i++){
    fCalGraph.SetPoint(i,Px.at(i),Py.at(i));
    fCalGraph.SetPointError(i,PxErr.at(i),PyErr.at(i));
  }

  TF1 *fLinearFunction;
  fLinearFunction = new TF1("fLinearFunction",Slope, 0.9*Px.at(0),1.1*Px.at(nCal-1),2);
  fLinearFunction->SetLineColor(2);
  fLinearFunction->SetLineWidth(2);
  fLinearFunction->SetLineStyle(2);

  TF1 *fLine = new TF1("fLine","[0]",0,1000000);
  fLine->SetParameter(0,0);
  fLine->SetLineColor(2);
  fLine->SetLineWidth(2);
  fLine->SetLineStyle(2);

  Double_t fScale = (Py.at(nCal-1)-Py.at(0))/(Px.at(nCal-1)-Px.at(0));
  Double_t fOffset = Py.at(0)-fScale*Px.at(0);
  fLinearFunction->SetParameters(fOffset,fScale);

  fCalGraph.Fit(fLinearFunction,"qremf");
  Double_t chisquare = fLinearFunction->GetChisquare();

  TFitResultPtr r = fCalGraph.Fit(fLinearFunction,"S");
  TMatrixDSym cov = r->GetCovarianceMatrix();
  
  fScale=fLinearFunction->GetParameter(1);
  Double_t fScaleErr=fLinearFunction->GetParError(1);
  fOffset=fLinearFunction->GetParameter(0);
  Double_t fOffsetErr=fLinearFunction->GetParError(0);

  Par->push_back(fOffset); //#0
  ParErr->push_back(fOffsetErr); 
  Par->push_back(fScale); // #1
  ParErr->push_back(fScaleErr);

  Double_t energyROI = Energy;
  Double_t LinearROI = (energyROI-fOffset)/fScale; 

  Double_t p[2];
  p[0] = -1/fScale;
  p[1] = -(energyROI-fOffset)/(fScale*fScale);

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
  for(Int_t k = 0;k<2;k++){
    for(Int_t j = 0;j<2;j++){
      if(k==j){
	ParErr->push_back(abs(cov(k,j)));
      }else{
	ParErr->push_back(cov(k,j));
      }
    }
  }
  //////////////////////////////////////////////
  if(FileName!="NA"){
    for(Int_t i = 0;i<nCal;i++){
      Double_t diff =  Py.at(i)-(fOffset+fScale*Px.at(i));
      Double_t differr = TMath::Sqrt(pow(fOffsetErr,2)+pow(fScaleErr*Px.at(i),2)+pow(PxErr.at(i)*fScale,2));
      fDeltaGraph.SetPoint(i,Py.at(i),diff);
      fDeltaGraph.SetPointError(i,0,differr);
      //fdelta << fEnergyName.c_str() << " " << Py.at(i) << " " << 0 << " " << diff << " " << differr << endl;
      if(abs(diff)>0.5){
	ferror << FileName.c_str() << " " << fEnergyName.c_str() << " " << i << " " << Px.at(i) << " " << PxErr.at(i)<<  " " <<  diff << " " << differr << endl;
      }
    }
    vector<Double_t> residual;
    vector<Double_t> residualerr;
    for(Int_t i = 0;i<nCal;i++){
      Double_t diff =  Px.at(i)-(Py.at(i)-fOffset)/fScale;
      Double_t differr = PxErr.at(i);
      residual.push_back(diff);
      residualerr.push_back(differr);
    }
    for(Int_t i = 0;i<nCal;i++){
      fResidual.SetPoint(i,Px.at(i),residual[i]);
      fResidual.SetPointError(i,PxErr[i],residualerr[i]);
    }    
    
    TCanvas *c1 = new TCanvas("c1");
    fCalGraph.Draw();
    c1->Print(Form("%s",CalibrationName.c_str()));

    fDeltaGraph.Draw("AP");
    fLine->Draw("same");
    c1->Print(Form("%s",DeltaName.c_str()));
    
    fResidual.Draw("AP");
    fLine->Draw("same");
    c1->Print(Form("%s",ResidualName.c_str()));
  }
}



void GATAutoCal::ResolutionFit(vector<Double_t> Px, vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName, string FileName, string TitleName, Double_t ROIEnergy){

  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");

  string ResolutionName=Form("%sresolution_%s.pdf",PathName.c_str(),FileName.c_str());
  string ResidualName=Form("%sresoresidual_%s.pdf",PathName.c_str(),FileName.c_str());
  ofstream fdelta(Form("resolution_%s.txt",FileName.c_str()),ios::app);
 
  const Int_t nCal=Px.size();
  if(nCal<2){
    cout<<"Cannot fit a calibration with less than 2 points!"<<endl;
  }
  Double_t A[nCal];
  Double_t B[nCal];
  Double_t AErr[nCal];
  Double_t BErr[nCal];

  Double_t maxX = 0;
  for(Int_t i = 0;i<nCal;i++){
    A[i] = Px.at(i);
    B[i] = Py.at(i);
    AErr[i] = PxErr.at(i);
    BErr[i] = PyErr.at(i);
    if(Px.at(i) > maxX){
      maxX = Px.at(i);
    }
  }

  TGraphErrors *fResGraph = new TGraphErrors(nCal,A,B,AErr,BErr);
  fResGraph->SetTitle(Form("%s; Energy (keV);FWHM (keV)",TitleName.c_str()));
  fResGraph->GetXaxis()->SetLimits(0.,maxX*1.1);
  fResGraph->SetMinimum(0);
  fResGraph->SetMarkerStyle(8);
  fResGraph->SetMarkerSize(2);
  fResGraph->SetMarkerColor(4);

  TGraphErrors fResidual;
  fResidual.Set(0);
  fResidual.SetTitle(Form("%s;Energy (keV);Residual",TitleName.c_str()));
  fResidual.SetMarkerStyle(8);
  fResidual.SetMarkerSize(2);
  fResidual.SetMarkerColor(4);

  TF1 *fResolutionFunction;
  fResolutionFunction = new TF1("fResolutionFunction",Reso,0,maxX*1.1,3);
  fResolutionFunction->SetRange(0,maxX*1.1);
  fResolutionFunction->SetLineColor(2);
  fResolutionFunction->SetLineWidth(2);
  fResolutionFunction->SetLineStyle(2);

  TF1 *fLine = new TF1("fLine","[0]",0,1000000);
  fLine->SetParameter(0,0);
  fLine->SetLineColor(2);
  fLine->SetLineWidth(2);
  fLine->SetLineStyle(2);

  fResGraph->Fit(fResolutionFunction,"qremf");
  TFitResultPtr r = fResGraph->Fit(fResolutionFunction,"S");
  TMatrixDSym cov = r->GetCovarianceMatrix();

  for(Int_t i = 0;i<3;i++){
    Par->push_back(fResolutionFunction->GetParameter(i)); // #0,#1,#2
    ParErr->push_back(fResolutionFunction->GetParError(i));
  }
  Double_t energyROI = ROIEnergy;
  Double_t ResoROI = TMath::Sqrt(Par->at(0)*Par->at(0)+Par->at(1)*Par->at(1)*energyROI+Par->at(2)*Par->at(2)*energyROI*energyROI);

  Double_t p[3];
  p[0] = 0.5/ResoROI*2*Par->at(0);
  p[1] = 0.5/ResoROI*2*Par->at(1)*energyROI;
  p[2] = 0.5/ResoROI*2*Par->at(2)*energyROI*energyROI;
  Double_t errsum = 0;
  for(Int_t i = 0;i<3;i++){
    for(Int_t j = 0;j<3;j++){
      if(i==j){
        errsum = errsum + p[i]*TMath::Abs(cov(i,j))*p[j];
      }else{
        errsum = errsum + p[i]*cov(i,j)*p[j];
      }
    }
  }
  Double_t ResoROIErr = TMath::Sqrt(errsum);
  Par->push_back(ResoROI); // #3
  ParErr->push_back(ResoROIErr);

  for(Int_t k = 0;k<3;k++){
    for(Int_t j = 0;j<3;j++){
      if(k==j){
	ParErr->push_back(abs(cov(k,j)));
      }else{
	ParErr->push_back(cov(k,j));
      }
    }
  }//#4-#12
  ////////////////////////////////////////////

  vector<Double_t> residual;
  vector<Double_t> residualerr;
  
  for(Int_t i = 0;i<nCal;i++){
    Double_t diff =  Py.at(i) - TMath::Sqrt(Par->at(0)*Par->at(0)+Par->at(1)*Par->at(1)*Px.at(i)+Par->at(2)*Par->at(2)*Px.at(i)*Px.at(i));;
    Double_t differr = PyErr.at(i);
    residual.push_back(diff);
    residualerr.push_back(differr);
    //fdelta << fStartRun << " " << fEndRun << " " << fEnergyName.c_str() << " " << Py.at(i) << " "<< PyErr.at(i) << " " << diff<< " " << differr << endl;
  }
  for(Int_t i = 0;i<nCal;i++){
    fResidual.SetPoint(i,Px.at(i),residual[i]);
    fResidual.SetPointError(i,PxErr[i],residualerr[i]);
  }

  if(PathName!="NA"){
    TCanvas *c1 = new TCanvas("c1");
    fResGraph->Draw("AP");
    c1->Print(Form("%s",ResolutionName.c_str()));

    fResidual.Draw("AP");
    fLine->Draw("same");
    c1->Print(Form("%s",ResidualName.c_str()));
  }
}

TGraphErrors* GATAutoCal::QuadFit(vector<Double_t> Px,vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName, string FileName, string TitleName){
  
  //gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  //cout << "Title Name : " << TitleName.c_str() <<  endl;
  string QuadName = Form("%squad_%s.pdf",PathName.c_str(),FileName.c_str());
 
  const Int_t nCal=Px.size();
  if(nCal<2){
    cout<<"Cannot fit a calibration with less than 2 points!"<<endl;
  }
  Double_t A[nCal];
  Double_t B[nCal];
  Double_t AErr[nCal];
  Double_t BErr[nCal];
  Double_t imin = 10000;
  Double_t imax =0;
  Double_t jmin = 10000;
  Double_t jminx = 0;
  for(Int_t i = 0;i<nCal;i++){
    A[i] = Px.at(i);
    B[i] = Py.at(i);
    AErr[i] = PxErr.at(i);
    BErr[i] = PyErr.at(i);
    if(Px.at(i)<imin){
      imin = Px.at(i);
    }
    if(Px.at(i)>imax){
      imax = Px.at(i);
    }
    if(Py.at(i)<jmin){
      jmin = Py.at(i);
      jminx = Px.at(i);
    }
  }

  TGraphErrors *fQuadGraph = new TGraphErrors(nCal,A,B,AErr,BErr);
  fQuadGraph->GetXaxis()->SetLimits(imin,imax);
  fQuadGraph->SetTitle(Form("%s;#tau'^{-1} (#mus^{-1});FWHM/Mean",TitleName.c_str()));
  fQuadGraph->SetMarkerStyle(8);
  fQuadGraph->SetMarkerColor(4);
  fQuadGraph->SetMarkerSize(2);
  fQuadGraph->SetLineColor(4);

  TF1 *fQuadFunction;
  fQuadFunction = new TF1("fQuadFunction",Quad,imin,imax,3);
  fQuadFunction->SetRange(0,imax*1.5);
  fQuadFunction->SetLineColor(2);
  fQuadFunction->SetLineWidth(2);
  fQuadFunction->SetLineStyle(2);
  TFitResultPtr r = fQuadGraph->Fit(fQuadFunction,"S");
  TMatrixDSym cov = r->GetCovarianceMatrix();
  for(Int_t i = 0;i<3;i++){
    Par->push_back(fQuadFunction->GetParameter(i));
    ParErr->push_back(fQuadFunction->GetParError(i));
  }
  if(FileName != "NA"){
    TCanvas *c1=new TCanvas("c1");
    fQuadGraph->GetYaxis()->SetTitleOffset(1.5);
    fQuadGraph->Draw("AP");    
    c1->Print(Form("%s",QuadName.c_str()));
  }
  Double_t xmin = fQuadFunction->GetMinimumX(0,imax*1.5,0.01);

  Par->push_back(jminx);
  Par->push_back(xmin);
  return fQuadGraph;
}



Double_t GATAutoCal::RampTimeFit(vector<Double_t> Px,vector<Double_t> PxErr, vector<Double_t> Py, vector<Double_t> PyErr, vector<Double_t>* Par, vector<Double_t>* ParErr, string PathName, string FileName, string TitleName){
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");

  string RampName = Form("%sramp_%s",PathName.c_str(), FileName.c_str());
  const Int_t nCal=Px.size();
  if(nCal<2){
    cout<<"Cannot fit a calibration with less than 2 points!"<<endl;
  }
  Double_t A[nCal];
  Double_t B[nCal];
  Double_t AErr[nCal];
  Double_t BErr[nCal];

  Double_t imin=10000;
  Double_t imax=0;
  Double_t jmin=10000;
  Double_t jminx = 0;
  for(Int_t i = 0;i<nCal;i++){
    A[i] = Px.at(i);
    B[i] = Py.at(i);
    AErr[i] = PxErr.at(i);
    BErr[i] = PyErr.at(i);
    if(Px.at(i)<imin){
      imin = Px.at(i);
    }
    if(Px.at(i)>imax){
      imax = Px.at(i);
    }
    if(Py.at(i)<jmin){
      jmin = Py.at(i);
      jminx = Px.at(i);
    }
  }

  TGraphErrors *fRampGraph = new TGraphErrors(nCal,A,B,AErr,BErr);
  fRampGraph->GetXaxis()->SetLimits(imin,imax);
  fRampGraph->SetTitle(Form("%s;#tau (us);FWHM/Mean",TitleName.c_str()));
  fRampGraph->SetMarkerStyle(8);
  fRampGraph->SetMarkerColor(4);
  fRampGraph->SetMarkerSize(2);
  fRampGraph->SetLineColor(4);

  TF1 *fRampFunction;
  fRampFunction = new TF1("fRampFunction",Ramp,imin,imax,3);
  fRampFunction->SetRange(0,imax*1.5);
  fRampFunction->SetLineColor(2);
  fRampFunction->SetLineWidth(2);
  fRampFunction->SetLineStyle(2);

  TFitResultPtr r = fRampGraph->Fit(fRampFunction,"S");
  TMatrixDSym cov = r->GetCovarianceMatrix();
  for(Int_t i = 0;i<3;i++){
    Par->push_back(fRampFunction->GetParameter(i));
    ParErr->push_back(fRampFunction->GetParError(i));
  }

  if(FileName != "NA"){
    TCanvas *c1 = new TCanvas("c1");
    fRampGraph->GetYaxis()->SetTitleOffset(1.5);
    fRampGraph->Draw("AP");    
    c1->Print(Form("%s",RampName.c_str()));
  }
  Double_t xmin = fRampFunction->GetMinimumX(0,imax*1.5,0.01);

  Par->push_back(jminx);
  Par->push_back(xmin);
  return xmin;
}

void GATAutoCal::PlotMultiGraph(TMultiGraph *MG, TLegend *Leg, string PathName, string FileName, string TitleName){
  //gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  cout << FileName.c_str() << endl;
  string PlotName = Form("%s%s",PathName.c_str(),FileName.c_str());
  TCanvas *c1 = new TCanvas("c1");
  MG->SetTitle(Form("%s",TitleName.c_str()));
  MG->Draw("AP");
  MG->GetYaxis()->SetTitleOffset(1.5);
  Leg->Draw();
  c1->Print(Form("%s",PlotName.c_str()));
}

void GATAutoCal::PlotGrid(vector<Int_t> Px, vector<Double_t> Py, vector<Double_t> PyErr, Double_t Ave, string PathName,string FileName,string TitleName){
  gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  string PlotName = Form("%s%s",PathName.c_str(),FileName.c_str());

  TCanvas *c1 = new TCanvas("c1");
  c1->SetGrid();
  const Int_t nPx = Px.size();
  string PosName[nPx];
  Int_t ic,is,id;
  for(size_t i=0;i<Px.size();i++){
    //cout << Px.at(i) << " " << ic<< " " << is << " " << id <<endl;
    ic = fCryo.at(Px.at(i));
    is = fString.at(Px.at(i));
    id = fDetector.at(Px.at(i));
    cout << Px.at(i) << " " << ic<< " " << is << " " << id <<endl;
 
    string detpos = Form("C%dP%dD%d",ic,is,id);
    PosName[i] = detpos;
  }

  Double_t A[nPx];
  Double_t AErr[nPx];
  Double_t B[nPx];
  Double_t BErr[nPx];

  for(Int_t j=0;j<nPx;j++){
    A[j]=j;
    AErr[j] =0;
    B[j]=Py.at(j);
    BErr[j] =PyErr.at(j);
  }
  TGraphErrors *fPeak;
  fPeak = new TGraphErrors(nPx,A,B,AErr,BErr);
  fPeak->SetMarkerStyle(8);
  fPeak->SetMarkerSize(2);
  fPeak->SetMarkerColor(4);

  TF1 *fFlatFun =new TF1("fFlatFun","[0]",0.,100.);
  fFlatFun->SetLineColor(2);
  fFlatFun->SetLineWidth(2);
  fFlatFun->SetLineStyle(2);
  fFlatFun->SetParameter(0,Ave);

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(fPeak,"P");
  mg->SetTitle(Form(";;%s (%s)",TitleName.c_str(),fEnergyName.c_str()));
  mg->Draw("AP");
  fFlatFun->Draw("same");
  
  TAxis *ax1 = mg->GetHistogram()->GetXaxis();
  Double_t x11 = ax1->GetBinLowEdge(3);
  Double_t x21 = ax1->GetBinUpEdge(100-2);
  mg->GetHistogram()->GetXaxis()->Set(nPx,x11,x21);

  for(Int_t j=0;j<nPx;j++){
    mg->GetHistogram()->GetXaxis()->SetBinLabel(j+1,PosName[j].c_str());
  }
  Float_t kAxisTitleSize = 18;
  mg->GetXaxis()->SetLabelSize(kAxisTitleSize);
  
  c1->Print(Form("/%s",PlotName.c_str()));
}

void GATAutoCal::PlotGrid(TMultiGraph *MG, TLegend *Leg, vector<Int_t> Px, vector<Double_t> Ave, Double_t Low, Double_t Up, string PathName,string FileName){

  //gROOT->ProcessLine(".x MJDTalkPlotStyle.C");
  string PlotName = Form("%s%s",PathName.c_str(),FileName.c_str()); 
  TCanvas *c1 = new TCanvas("c1", "c1",243,130,1000,600);
  //TCanvas *c1 = new TCanvas("c1");
  c1->SetGrid();
  const Int_t nPx = Px.size();
  string PosName[nPx];
  Int_t px;
  Int_t ic,is,id;
  string det;

  for(size_t i=0;i<Px.size();i++){
    px = (Int_t)Px.at(i);
    ic = px/100;
    is = (px-ic*100)/10;
    id = (px-ic*100-is*10);
    //ic = fCryo.at(px);
    //is = fString.at(px);
    //id = fDetector.at(px);
    //string detpos = Form("#splitline{C%dP%dD%d}{%d}",ic,is,id,fChannel.at(px));
    string detpos = Form("C%dP%dD%d",ic,is,id);
    PosName[i] = detpos;
  }

  const Int_t nAve = Ave.size();
  TF1 *fFlatFun[nAve];
  for(size_t ia = 0;ia<Ave.size();ia++){
    fFlatFun[ia] =new TF1("fFlatFun","[0]",0.,100.);
    fFlatFun[ia]->SetLineColor(ia+2);
    fFlatFun[ia]->SetLineWidth(2);
    fFlatFun[ia]->SetLineStyle(2);
    fFlatFun[ia]->SetParameter(0,Ave[ia]);
  }
  //MG->SetTitle(Form("DS0:1000-1400 keV;;%s",TitleName.c_str()));
  MG->SetMaximum(Up);
  MG->SetMinimum(Low);

  MG->Draw("AP");
  for(size_t ia=0;ia<Ave.size();ia++){
    fFlatFun[ia]->Draw("same");
  }
  TAxis *ax1 = MG->GetHistogram()->GetXaxis();
  Double_t x11 = ax1->GetBinLowEdge(4);
  Double_t x21 = ax1->GetBinUpEdge(97);
  MG->GetHistogram()->GetXaxis()->Set(nPx,x11,x21);

  for(Int_t j=0;j<nPx;j++){
    MG->GetHistogram()->GetXaxis()->SetBinLabel(j+1,PosName[j].c_str());
    //MG->GetHistogram()->GetXaxis()->SetBinLabel(MG->GetHistogram()->GetBin(j+1),PosName[j].c_str());
  }
  Float_t kAxisTitleSize = 20;
  MG->GetXaxis()->SetLabelSize(kAxisTitleSize);
  MG->GetXaxis()->SetLabelFont(5);
  MG->GetXaxis()->SetBit(TAxis::kLabelsVert);
  Leg->Draw();
  c1->Print(Form("%s",PlotName.c_str()));
}

void GATAutoCal::PlotSpectrum(TH1D *Hist,string FileName){
  TCanvas *c1 = new TCanvas("c1");
  TH1D *hc = (TH1D*)Hist->Clone();
  c1->SetLogy();
  hc->Draw();
  c1->Print(Form("%s",FileName.c_str()));
  delete c1;
}


void GATAutoCal::PlotSpectrum2(TH1D *H1,TH1D *H2,TLegend *Leg,string FileName){

  TCanvas *c1 = new TCanvas("c1");
  Double_t m1= H1->GetMaximum();
  Double_t m2= H2->GetMaximum();
  if(m1>m2){
    H1->SetMaximum(m1);
  }else{
    H1->SetMaximum(m2);
  }
  c1->SetLogy();

  H1->Draw();
  H2->Draw("same");
  Leg->Draw();

  c1->Print(Form("%s",FileName.c_str()));
}

void GATAutoCal::PlotGraph(TGraphErrors *Graph,string FileName){
  TCanvas *c1 = new TCanvas("c1");
  Graph->Draw("AP");
  c1->Print(Form("%s",FileName.c_str()));
  delete c1;
}
/*
void GATAutoCal::PlotMultiGraph(TMultiGraph *Graph, TLegend *Leg, string FileName){
  TCanvas *c1 = new TCanvas("c1");
  Graph->Draw("AP");
  Leg->Draw();
  c1->Print(Form("%s",FileName.c_str()));
  delete c1;
}
*/


TProfile* GATAutoCal::PeakProfileTime(TChain *fChain,string fCut,Double_t PeakChannel,Double_t PeakChannelWindow, string PathName,string FileName){
  TCanvas *c1 = new TCanvas("c1");
  Int_t windowbin = (Int_t) PeakChannelWindow*2;
  //GATDataSet ds(fStartRun,fEndRun);
  //TChain* chain = ds.GetGatifiedChain(0);
  Double_t starttime = GATAutoCal::GetStartTime();
  Double_t stoptime = GATAutoCal::GetStopTime();
  Double_t totaltime = (stoptime-starttime)/3600;
  cout << totaltime << endl;
  
  Int_t hourbin = (Int_t) totaltime + 1;
  Double_t totaltimeplus = (Double_t) hourbin;
  TH2F *h2 = new TH2F("h2","",hourbin,0,totaltimeplus,windowbin,PeakChannel-PeakChannelWindow,PeakChannel+PeakChannelWindow);
  fChain->Draw(Form("trapENFCal:(startTime-%f+timestamp*1e-8)/3600.>>h2",starttime),Form("(%s) && TMath::Abs(%s-%f)<%f",fCut.c_str(),fEnergyName.c_str(),PeakChannel,PeakChannelWindow));
  h2->Draw();
  c1->Print(Form("%splot_%s.pdf",PathName.c_str(),FileName.c_str()));
  Double_t meanY = h2->GetMean(2);
  Double_t rmsY = h2->GetRMS(2);

  TProfile *px = h2->ProfileX("px");

  px->GetYaxis()->SetRangeUser(meanY-rmsY*2,meanY+rmsY*2);
  if(PathName != "NA"){
    px->Draw();
    c1->Print(Form("%sprofiletime_%s.pdf",PathName.c_str(),FileName.c_str()));
  }
  return px;
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

  //Dates -> ask Pinghan
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

//This creates a single new record of peak locations in the database
void GATAutoCal::PutPeakLocationMJDB(Int_t channel, string detectorid,
				     vector<Double_t> &calpeak,
				     vector<Double_t> &calpeakerr,
				     vector<Double_t> &calwidth,
				     vector<Double_t> &calwidtherr
				     )
{
  //REMOVE this after testing
  //This line adds the non-standard DB to the AnalysisDoc readout
  fACAnalysisDB.SetAlternateDB("", "mjdbatest","","","");
  
  MJDB::PeakLocationList myPeakList;
  fACProvenance.SetParametersType(myPeakList.GetPType());//set the Parameter type
  fACProvenance.SetChannelIdentifier(MJDB::NumberToString(channel));
  fACProvenance.SetDetectorIdentifier(detectorid);

  MJDB::PeakLocation myPLoc;
  for(size_t i=0; i<calpeak.size(); i++) {
    myPLoc.PeakID.SetValue(i+1);
    myPLoc.Location.SetValue(calpeak.at(i));
    myPLoc.Location.SetUncertainty(calpeakerr.at(i));
    myPLoc.FWHM.SetValue(calwidth.at(i));
    myPLoc.FWHM.SetUncertainty(calwidtherr.at(i));
    myPeakList.AddLocation(myPLoc);
  }
  
  //Build a database record
  fACProvenance.SetDBProvenance(fACAnalysisDB);
  if (myPeakList.SetDBValue(fACAnalysisDB)) { //Make sure we save the parameters
    
    int status = fACAnalysisDB.New_MJAnalysis_Record();
    if (status) {
      std::string statusMessage;
      fACAnalysisDB.GetDBStatus(statusMessage);
      cout << "\"" << statusMessage << "\"" << endl;
      return; //Error
    } 
  } else {
    //Reach here only if the parameter record is incomplete. Bad error
    cout << "Peak Location Parameter is record not completely initialized."
	 << " Try again."
	 << endl;
  }
}
//This reads a single record of peak locations from the database
void GATAutoCal::GetPeakLocationMJDB(Int_t channel,
				     vector<Double_t> &calpeak,
				     vector<Double_t> &calpeakerr,
				     vector<Double_t> &calwidth,
				     vector<Double_t> &calwidtherr,
				     Int_t &calibrationpeaks)
{
  //The search returns an AnalysisDoc
  MJAnalysisDoc findResult;

  //Extra keys is a KTree, need a NullTree if there are no extra Keys
  katrin::KTree::KEmptyArray NullTree;
  katrin::KTree ExtraKey;
  

  //REMOVE this for real database operation
  fACAnalysisDB.SetAlternateDB("", "mjdbatest","","","");

  ExtraKey = NullTree;
  MJDB::PeakLocationList myPeakList;
  size_t Length = 
    findResult.GetAnalysisParameter((fStartRun),
				    (channel), 
				    myPeakList.GetPType(),
				    ExtraKey);
  
  
  //Loop over all the records returned from the DB
  for (size_t i=0; i< Length; i++) {
    MJAnalysisDoc temp = findResult[i];
    
    myPeakList.Clear();  //Clear away the previous iteration
    
    myPeakList.GetDBValue(temp);
    
    //This loop shows how to access the parameters
    calibrationpeaks = myPeakList.NumberOfLocations();
    
    for (size_t i=0; i<myPeakList.NumberOfLocations(); i++) {
      
      PeakLocation myPeakLoc;
      myPeakLoc = myPeakList.GetLocation(i);
      
      //Locations
      calpeak.push_back(myPeakLoc.Location.Value());
      calpeakerr.push_back(myPeakLoc.Location.Uncertainty());
      //FWHM
      calwidth.push_back(myPeakLoc.FWHM.Value());
      calwidtherr.push_back(myPeakLoc.FWHM.Uncertainty());
    }
  }
}


