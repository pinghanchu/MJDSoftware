// 2016.3.23
// Pulser Tagger written by Pinghan Chu (pinghan@gmail.com)
// This code can be only applied to background runs, but calibration runs.
// Currently, the code requires pulser peak location as an input. 
// In the future, the code will get the pulser peak location from database.
//
// 2016.4.5 Ported to GAT by J. Detwiler
// 2016.11.14 
// Modify the code in order to work the combined data of Mod1 and Mod2
//
#include "GATPulserTag.hh"
#include "MAPP3END.hh"
#include "MAPP3JDY.hh"
#include "MAPP3KJR.hh"
#include "MAPP3LQG.hh"
#include "MAPP3LQK.hh"
#include "MJTRun.hh"
#include "MJAnalysisDoc.hh"
#include "MJAnalysisParameters.hh"
#include "TMath.h"
#include "TStyle.h"
#include "TFile.h"
#include <iostream>
#include <bitset>

//ClassImp(GATPulserTag)

using namespace std;
using namespace katrin;
using namespace MJDB;


GATPulserTag::GATPulserTag(Int_t run) : fDS(run)
{
  fRun = run;
  fMap = fDS.GetChannelMap();
  fMjdTree = fDS.GetGatifiedChain(false);
  fMjdTree->SetBranchStatus("*",0);
  fMjdTree->SetBranchStatus("run",1);
  fMjdTree->SetBranchStatus("C",1);
  fMjdTree->SetBranchStatus("P",1);
  fMjdTree->SetBranchStatus("D",1);
  fMjdTree->SetBranchStatus("channel", 1);
  fMjdTree->SetBranchStatus("trapE", 1);
  fMjdTree->SetBranchStatus("timeinfo", 1);
  fMjdTree->SetBranchStatus("clockTime", 1);
  fMjdTree->SetBranchStatus("tOffset", 1);
  fMjdTree->SetBranchStatus("mH",1);
  fMjdTree->SetBranchStatus("mL",1);

  fMTChannel = NULL;
  fMTTrapE = NULL;
  fTimeInfo = NULL;
  fMTC = NULL;
  fMTP = NULL;
  fMTD = NULL;
  fMjdTree->SetBranchAddress("run", &fMTRun);
  fMjdTree->SetBranchAddress("C", &fMTC);
  fMjdTree->SetBranchAddress("P", &fMTP);
  fMjdTree->SetBranchAddress("D", &fMTD);
  fMjdTree->SetBranchAddress("channel", &fMTChannel);
  fMjdTree->SetBranchAddress("trapE", &fMTTrapE);
  fMjdTree->SetBranchAddress("timeinfo", &fTimeInfo);
  fMjdTree->SetBranchAddress("mH", &fmH);
  fMjdTree->SetBranchAddress("mL", &fmL);

  fMjdTree->GetEntry(0);
  fEntries = fMjdTree->GetEntries();

  SetParameters();

  TChain* bc = fDS.GetBuiltChain(false);
  MJTRun* runInfo = NULL; //new MJTRun();
  bc->SetBranchAddress("run",&runInfo);
  bc->GetEntry(0);
  bitset<32> event_type = runInfo->GetRunBits();
  bc->SetBranchAddress("run",NULL);

  if(event_type.test(5) || event_type.test(6)) fIsRadio = 1;
  else fIsRadio = 0;
}


void GATPulserTag::SetParameters()
{
  if(fRun >= 2337 && fRun < 8183){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad0[i]);
      if(fRun<3464) fPulser.push_back(Pulser0[i]);
      else if(fRun>3464 && fRun<4573) fPulser.push_back(Pulser1a[i]);
      else fPulser.push_back(Pulser1[i]);
    }
  }
  if(fRun >= 8722 && fRun<18589){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad3[i]);
      if(fRun>=16797 && fRun<=17009){
	fPulser.push_back(Pulser3a[i]);
      }else{
	fPulser.push_back(Pulser3[i]);
      }
    }
  }

  if(fRun >= 18623 && fRun<25159){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad6[i]);
      fPulser.push_back(Pulser6[i]);
    }
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad7[i]);
      fPulser.push_back(Pulser7[i]);
    }
  }

  if(fRun >= 25159 && fRun<25675){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad6[i]);
      fPulser.push_back(Pulser8[i]);
    }
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad7[i]);
      fPulser.push_back(Pulser9[i]);
    }
  }
  if(fRun >= 25675 && fRun<50000000){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad6[i]);
      fPulser.push_back(Pulser6[i]);
    }
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad7[i]);
      fPulser.push_back(Pulser7[i]);
    }
  }


  if(fRun >= 60000001 && fRun<60000550){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad4[i]);
      fPulser.push_back(Pulser5[i]);
    }
  }

  if(fRun >= 60000550 && fRun<70000000){
    for(Int_t i = 0;i<kChannels;i++){
      fGoodBad.push_back(GoodBad5[i]);
      fPulser.push_back(Pulser5[i]);
    }
  }

  if(fRun >= 45000000 && fRun<50000000){
    for(Int_t i = 0;i<proChannels;i++){
      fGoodBad.push_back(GoodBadpro[i]);
      fPulser.push_back(Pulserpro[i]);
    }
  }

  SetPulserWindow(0.1); // 10% of the pulser energy
  SetPulserTimeWindow(10); // 10 seconds time window
  SetChannel();
  SetDataSet();
}


void GATPulserTag::SetDataSet()
{
  if(fRun>0 && fRun<=111) fDataSet = "P3GKF";
  else if(fRun>=202 && fRun<=1041) fDataSet = "P3HUA";
  else if(fRun>=1042 && fRun<=2334) fDataSet = "P3JCJ";
  else if(fRun>=2337 && fRun<=8183) fDataSet = "P3JDY";
  else if(fRun>=8184 && fRun<=8721) fDataSet = "P3K93";
  else if(fRun>=8722 && fRun<18589) fDataSet = "P3KJR";
  else if(fRun>=18623 && fRun <= 25671) fDataSet = "P3LQK";
  else if(fRun>=25672 && fRun < 5000000) fDataSet = "P3LTP";
  else if(fRun>=60000000 && fRun<60000549) fDataSet = "P3LMF";
  else if(fRun>=60000550 && fRun<70000000) fDataSet = "P3LQG";
  else if(fRun>= 45000001&& fRun<=45008360) fDataSet = "P3END";
  else fDataSet = "NA";
}


void GATPulserTag::SetChannel()
{

  fChannel.erase(fChannel.begin(),fChannel.end());
  if(fRun>=0 && fRun<18623){
    for(Int_t i =1;i<2;i++){
      for(Int_t j = 1;j<=7;j++){
	if(j==4){
	  for(Int_t k = 1;k<=5;k++){
	    fChannel.push_back(fMap->GetIDHi(i,j,k));
	    fChannel.push_back(fMap->GetIDLo(i,j,k));
	  }
	}else{
	  for(Int_t k = 1;k<=4;k++){
	    fChannel.push_back(fMap->GetIDHi(i,j,k));
	    fChannel.push_back(fMap->GetIDLo(i,j,k));
	  }
	}
      }
    }
  }
  if(fRun>=60000001 && fRun <70000000){
    for(Int_t i =2;i<3;i++){
      for(Int_t j = 1;j<=7;j++){
        if(j==2 || j==4){
          for(Int_t k = 1;k<=5;k++){
            fChannel.push_back(fMap->GetIDHi(i,j,k));
            fChannel.push_back(fMap->GetIDLo(i,j,k));
          }
        }else if(j==3){
          for(Int_t k = 1;k<=3;k++){
            fChannel.push_back(fMap->GetIDHi(i,j,k));
            fChannel.push_back(fMap->GetIDLo(i,j,k));
          }
	}else{
          for(Int_t k = 1;k<=4;k++){
            fChannel.push_back(fMap->GetIDHi(i,j,k));
            fChannel.push_back(fMap->GetIDLo(i,j,k));
          }
        }
      }
    }
    fStrings = 7;
  }

  if(fRun>=18623 && fRun <50000000){
    for(Int_t i =1;i<3;i++){
      if(i==1){
        for(Int_t j = 1;j<=7;j++){
          if(j==4){
            for(Int_t k = 1;k<=5;k++){
              fChannel.push_back(fMap->GetIDHi(i,j,k));
              fChannel.push_back(fMap->GetIDLo(i,j,k));
            }
          }else{
            for(Int_t k = 1;k<=4;k++){
              fChannel.push_back(fMap->GetIDHi(i,j,k));
	      fChannel.push_back(fMap->GetIDLo(i,j,k));
            }
          }
        }
      }
      if(i==2){
        for(Int_t j = 1;j<=7;j++){
          if(j==2 || j==4){
            for(Int_t k = 1;k<=5;k++){
              fChannel.push_back(fMap->GetIDHi(i,j,k));
              fChannel.push_back(fMap->GetIDLo(i,j,k));
            }
          }else if(j==3){
            for(Int_t k = 1;k<=3;k++){
              fChannel.push_back(fMap->GetIDHi(i,j,k));
              fChannel.push_back(fMap->GetIDLo(i,j,k));
            }
          }else{
            for(Int_t k = 1;k<=4;k++){
              fChannel.push_back(fMap->GetIDHi(i,j,k));
              fChannel.push_back(fMap->GetIDLo(i,j,k));
	    }
	  }
	}
      }
    }
    fStrings = 14;
  }
  

  if(fRun >= 45000000 && fRun<50000000){
    for(Int_t i = 0;i<proChannels;i++){
      fChannel.push_back(Channelpro[i]);
    }
    fStrings = 3;
  }

  fPulserTagChannel.erase(fPulserTagChannel.begin(),fPulserTagChannel.end());
  vector<uint32_t> pulsertagchan = fMap->GetPulserChanList();
  for(Int_t i =0;i<(Int_t)pulsertagchan.size();i++){
    Int_t chan = pulsertagchan[i];
    fPulserTagChannel.push_back(chan);
  }
  fChannels = fChannel.size();
}


void GATPulserTag::SetPulser()
{
  fPulser.erase(fPulser.begin(),fPulser.end());
  for(Int_t i = 0; i<fChannels; i++){
    Int_t chan = fChannel[i];
    Double_t par = 0;
    Double_t parerr = 0;
    GetPulserPeakFromDB(fRun,fRun,chan, &par, &parerr);
    fPulser.push_back(par);  
  }
}


void GATPulserTag::GetPulserPeakFromDB(Int_t LrunNumber, Int_t UrunNumber, Int_t chanNumber, Double_t* par, Double_t* parerr)
{
  MJAnalysisDoc findResult;
  findResult.SetAlternateDB("", "mjdbsandbox", "", "", "");

  PulserPeak myPulser_Peak;
  std::vector<string> ptype_list;
  ptype_list.push_back(myPulser_Peak.GetPType());

  std::string pType;
  std::string pSource="";
  for (Int_t i=LrunNumber; i <= UrunNumber;i++) {
    for (size_t j=0; j<ptype_list.size(); j++) {
      pType = ptype_list[j];
      size_t Length = findResult.GetAnalysisParameter(i, chanNumber, pSource, pType);	

      for (size_t k=0; k< Length; k++) {
	MJAnalysisDoc temp = findResult[k];
	if (pType == myPulser_Peak.GetPType()) {
	  myPulser_Peak.GetDBValue(temp);
	  if(k==Length-1){
	    *par = myPulser_Peak.Location.Value();
	    *parerr = myPulser_Peak.Location.Uncertainty();
	  }
	}
      }
    }
  }
}


void GATPulserTag::Tag1(const vector<Int_t>& inputEntry, vector<Int_t>& outputEntry)
{
  map<int,int> detid;
  vector<Double_t> RunEnergy;
  vector<Double_t> RunEnergy1;
  vector<Double_t> RunEnergy2;
  for(int i=0; i<fChannels; i++){
    detid[fChannel.at(i)] = i;
    RunEnergy.push_back(0);
    RunEnergy1.push_back(0);
    RunEnergy2.push_back(0);
  }

  vector<Int_t> TagEntry;
  vector<Int_t> M1Entry;
  vector<Int_t> M2Entry;
  
  for(size_t i=0; i<inputEntry.size(); i++){
    fMjdTree->GetEntry(inputEntry[i]);
    Int_t ncha = fMTChannel->size();
    Int_t npulser = 0;

    for(Int_t j=0; j<fChannels; j++){
      RunEnergy.at(j) = 0;
    }
    for(Int_t j=0; j<ncha; j++) {
      Int_t cha = (*fMTChannel)[j];
      double energy = (*fMTTrapE)[j];
      RunEnergy.at(detid[cha]) = energy;
    }
    for(Int_t j=0; j<fChannels/2; j++) {
      Int_t k1= j*2;
      Int_t k2= k1+1;
      if( (RunEnergy.at(k1)>10 && abs(RunEnergy.at(k1)-fPulser.at(k1))< fPulser.at(k1)*fPulserWindow)||
	  (RunEnergy.at(k2)>10 && abs(RunEnergy.at(k2)-fPulser.at(k2))< fPulser.at(k2)*fPulserWindow)){
	npulser++;
      }
    }

    if(npulser>0) M1Entry.push_back(inputEntry[i]); // pulser candidate; E within the pulser raw E 10%
    else M2Entry.push_back(inputEntry[i]);
  }

  Int_t M1Entries = M1Entry.size();
  Double_t timewindow = 0;
  Int_t two = 0;
  Int_t one = 0;
  cout.precision(12);
  if(M1Entries >1){
    for(Int_t i=0;i<M1Entries;i++){
      fMjdTree->GetEntry(M1Entry[i]);
      Double_t t1 = fTimeInfo->clockTime/1e9;
      cout << i << " "<< t1 << endl; 
      for(Int_t j=0; j<fChannels; j++) RunEnergy1[j] = 0;
      one = 0;
      for(size_t j=0; j<fMTChannel->size(); j++){
	Int_t cha = (*fMTChannel)[j];
	if(cha == 677 && fRun <8183) one =1;
	Double_t energy = (*fMTTrapE)[j];
	RunEnergy1.at(detid[cha]) = energy;
      }

      timewindow = 0;
      two = 0;

      if(i<M1Entries-1){
	Int_t k = i+1;
	while(timewindow <fPulserTimeWindow && two < 1 && k<M1Entries){
	  fMjdTree->GetEntry(M1Entry[k]);
	  Double_t t2 = fTimeInfo->clockTime/1.e9;
	  timewindow = t2-t1;  
	  for(Int_t j=0; j<fChannels; j++) RunEnergy2[j] = 0;
	  for(size_t j=0; j<fMTChannel->size(); j++){
	    Int_t cha = (*fMTChannel)[j];
	    Double_t energy = (*fMTTrapE)[j];
	    RunEnergy2.at(detid[cha]) = energy;
	  }
	  for(Int_t j = 0; j<fChannels/2;j++){
	    Int_t k1 = j*2;
	    Int_t k2 = k1+1;
	    if((RunEnergy1.at(k1)>10 && abs(RunEnergy1.at(k1)-RunEnergy2.at(k1))<5) || 
	       (RunEnergy1.at(k2)>10 && abs(RunEnergy1.at(k2)-RunEnergy2.at(k2))<5)){
	      two++;
	    }
	  }
	  k++;
	}
      }

      timewindow = 0;
      if(i>1 && two < 1){
	Int_t k = i-1;
	while( timewindow < fPulserTimeWindow && two < 1 && k>=0){
	  fMjdTree->GetEntry(M1Entry[k]);
	  Double_t t2 = fTimeInfo->clockTime/1.e9;
	  timewindow = t1-t2;
	    
	  for(Int_t j=0;j<fChannels;j++) RunEnergy2[j] = 0;
	  for(size_t j=0; j<fMTChannel->size(); j++) {
	    Int_t cha = (*fMTChannel)[j];
	    Double_t energy = (*fMTTrapE)[j];
	    RunEnergy2.at(detid[cha]) = energy;
	  }
	  for(Int_t j = 0; j<fChannels/2;j++){
	    Int_t k1 = j*2;
	    Int_t k2 = k1+1;
	    if( (RunEnergy1.at(k1)>10 && abs(RunEnergy1.at(k1)-RunEnergy2.at(k1))<5)
		||(RunEnergy1.at(k2)>10 && abs(RunEnergy1.at(k2)-RunEnergy2.at(k2))<5)){
	      two++;
	    }
	  }
	  k--;
	}
      }
  
      if(two>0 || one>0) TagEntry.push_back(M1Entry[i]);
      else M2Entry.push_back(M1Entry[i]);
    }
  }

  if(M1Entry.size()==1) M2Entry.push_back(M1Entry[0]);

  if(M2Entry.size()>0){
    vector<Double_t> tempList;
    const Int_t nentry = M2Entry.size();
    for(Int_t i = 0;i<nentry;i++){
      tempList.push_back((Double_t)M2Entry[i]);
    }
    sort(tempList.begin(),tempList.end());
    outputEntry.resize(0);
    for(Int_t i =0;i<(Int_t)tempList.size();i++){
      outputEntry.push_back((Int_t)tempList[i]);
    }
  }
}


void GATPulserTag::Tag2(const vector<Int_t>& inputEntry, vector<Int_t>& outputEntry)
{
  if(fPulserTagChannel.size() == 0){
    outputEntry = inputEntry;
    return;
  }

  outputEntry.resize(0);
  for(size_t i=0; i<inputEntry.size(); i++) {
    fMjdTree->GetEntry(inputEntry[i]);
    Int_t npulser = 0;
    for(size_t j=0; j<fMTChannel->size(); j++){
      Int_t cha = (*fMTChannel)[j];
      for(Int_t k=0;k<(Int_t)fPulserTagChannel.size();k++){
        if(cha == fPulserTagChannel[k]){
          npulser++;
        }
      }
    }
    if(npulser == 0){
      outputEntry.push_back(inputEntry[i]);
    }
  }
}


void GATPulserTag::Tag3(const vector<Int_t>& inputEntry, vector<Int_t>& outputEntry)
{
  outputEntry.resize(0);
  vector<Double_t> RunEnergy;
  map<Int_t,Int_t> detid;
  for(Int_t i=0;i<fChannels;i++){
    detid[fChannel.at(i)]= i;
    RunEnergy.push_back(0);
  }
  vector<Int_t> TagEntry;
  vector<Int_t> M1Entry;
  vector<Int_t> M2Entry;

  for(size_t i=0; i<inputEntry.size(); i++) {
    fMjdTree->GetEntry(inputEntry[i]);
    Int_t npulser = 0;
    for(Int_t j=0; j<fChannels; j++) RunEnergy.at(j) = 0;
    for(size_t j=0; j<fMTChannel->size(); j++){
      Int_t cha = (*fMTChannel)[j];
      double energy = (*fMTTrapE)[j];
      RunEnergy.at(detid[cha]) = energy;
    }
    
    for(Int_t j=0;j<fChannels/2;j++){
      Int_t k1= j*2;
      Int_t k2= k1+1;
      if( (RunEnergy.at(k1)>10 && abs(RunEnergy.at(k1)-fPulser.at(k1))< fPulser.at(k1)*fPulserWindow)||
	  (RunEnergy.at(k2)>10 && abs(RunEnergy.at(k2)-fPulser.at(k2))< fPulser.at(k2)*fPulserWindow)){
	npulser++;
      }
    }

    if(npulser>1) M1Entry.push_back(inputEntry[i]); // pulser candidate; the E within the pulser raw E 10%
    else M2Entry.push_back(inputEntry[i]);
  }
  
  if(M2Entry.size()>0){
    const Int_t nentry = M2Entry.size();
    Int_t N1[nentry];
    Int_t N2[nentry];
    for(Int_t i = 0;i<nentry;i++){
      N1[i] = M2Entry[i];
    }
    TMath::Sort(nentry,N1,N2,0);
    for(Int_t i=0;i<nentry;i++){
      Int_t j= N2[i];
      outputEntry.push_back(N1[j]);
    }
  }
}

void GATPulserTag::Tag(const vector<Int_t>& inputEntry, vector<Int_t>& outputEntry)
{
  //Tag = Tag1 + Tag2
  if(fRun >= 4549){
    vector<Int_t> m1;
    Tag1(inputEntry, m1);
    Tag2(m1, outputEntry);
  }
  else Tag1(inputEntry, outputEntry);
}


void GATPulserTag::TagCorrect(vector<Int_t>& inputList, const vector<Int_t>& correctList)
{
  //Add corrected events
  for(size_t i = 0; i<correctList.size(); i++){
    inputList.push_back(correctList[i]);
  }
}


void GATPulserTag::GetTagList(const vector<Int_t>& inputEntry, vector<Int_t>& outputEntry)
{
  // Input the physical event list and output the pulser tag through 
  // all entries (0=physical event, 1=pulser event)
  outputEntry.resize(0);
  for(size_t i = 0; i<fEntries; i++) outputEntry.push_back(1);
  for(size_t i = 0; i<inputEntry.size(); i++) outputEntry[inputEntry[i]] = 0;
}


void GATPulserTag::GetPulserList(const vector<Int_t>& inputEntry, vector<Int_t>& PList)
{
  //Input the physical event list and output the pulser event list
  vector<Int_t> outputEntry;
  for(size_t i = 0; i<fEntries; i++) outputEntry.push_back(1);

  for(Int_t i = 0;i<(Int_t)inputEntry.size();i++){
    outputEntry[inputEntry[i]] = 0;
  }

  PList.resize(0);
  for(size_t i=0; i<fEntries; i++) {
    if(outputEntry[i] == 1) PList.push_back(i);
  }
}

//input the physical event list and output the pulser event of the channel list
void GATPulserTag::GetPulserChannelList(const vector<Int_t>& inputEntry, Int_t index, vector<Int_t>& PList)
{
  PList.resize(0);
  for(size_t i=0; i<inputEntry.size(); i++){
    fMjdTree->GetEntry(inputEntry[i]);
    for(size_t j=0; j<fMTChannel->size(); j++){
      Int_t chan = (*fMTChannel)[j];
      Int_t tE = (*fMTTrapE)[j];
      if(chan == fChannel[index] && abs(tE-fPulser[index])<fPulser[index]*0.1){
	PList.push_back(inputEntry[i]);
      }
    }
  }
}

//Input the pulser event list and output the time of the pulser event of the channel
void GATPulserTag::GetPulserTime(const vector<Int_t>& inputEntry, Int_t index, vector<Double_t>& PTime)
{
  PTime.resize(0);
  for(Int_t i=0;i<(Int_t)inputEntry.size();i++){
    fMjdTree->GetEntry(inputEntry[i]);
    for(size_t j = 0 ;j<fMTChannel->size();j++){
      Int_t chan = (*fMTChannel)[j];
      Double_t energy = (*fMTTrapE)[j];
      if(chan == fChannel[index] && abs(energy-fPulser[index])<fPulser[index]*0.1){
	Double_t ts = (fTimeInfo->clockTime + fTimeInfo->tOffset[j])/1e9;
	PTime.push_back(ts);
      }
    }
  }
}


Double_t GATPulserTag::GetPulserDiffTime(vector<Double_t> PulserTime)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  Double_t min=0;
  if(PulserTime.size()>0){
    vector<Double_t> DiffTime;
    Double_t difftime;
    for(Int_t i =1;i<(Int_t)PulserTime.size();i++){
      difftime = PulserTime[i] - PulserTime[i-1];
      DiffTime.push_back(difftime);
    }
    sort(DiffTime.begin(),DiffTime.end());
    
    Double_t diffdifftime;
    Double_t xi,xj;
    Double_t diffpulser = 1e-3;
    
    vector<Double_t> DiffTimeList;
    DiffTimeList.push_back(DiffTime[0]);
    for(Int_t i=1;i<(Int_t)DiffTime.size();i++){
      xi = DiffTime[i-1];
      xj = DiffTime[i];
      diffdifftime = xj-xi;
      if(diffdifftime>diffpulser){
	DiffTimeList.push_back(xj);
      }
    }
    const Int_t DiffTimeLists = DiffTimeList.size();
    Int_t DiffTimeNum[DiffTimeLists];
    for(Int_t j=0;j<DiffTimeLists;j++){
      DiffTimeNum[j] = 0;
    }
    
    for(Int_t i=0;i<(Int_t)DiffTime.size();i++){
      xi = DiffTime[i];
      for(Int_t j=0;j<DiffTimeLists;j++){
	if(abs(xi-DiffTimeList[j])<diffpulser){
	  DiffTimeNum[j] = DiffTimeNum[j]+1;
	}
      }
    }
    
    vector<Double_t> DiffTimeList2;
    vector<Int_t> DiffTimeNum2;
    for(Int_t j=0;j<DiffTimeLists;j++){
      if(DiffTimeNum[j]>=2){
	DiffTimeList2.push_back(DiffTimeList[j]);
	DiffTimeNum2.push_back(DiffTimeNum[j]);
      }
    }
    
    Double_t x2;
    Int_t y2;
    Double_t z2;
    if(DiffTimeList2.size()==1){
      min = DiffTimeList2[0];
    }else if(DiffTimeList2.size()==2){
      x2 = DiffTimeList2[1]/DiffTimeList2[0];
      y2 = (Int_t)x2;
      z2 = (Double_t) y2;
      if ((abs(x2-z2)<diffpulser) || (abs(x2-(z2+1))<diffpulser)){
	min = DiffTimeList2[0];
      }else{
	if(DiffTimeNum2[1]>DiffTimeNum2[0]){
	  min = DiffTimeList2[1];
	}else{
	  min = DiffTimeList2[0];
	}
      }
    }else{
      for(Int_t i = 0;i<(Int_t)DiffTimeList.size()-1;i++)
	{
	  x2 = DiffTimeList[i+1]/DiffTimeList[i];
	  y2 = (Int_t)x2;
	  z2 = (Double_t) y2;
	  if ((abs(x2-z2)<diffpulser) || (abs(x2-(z2+1))<diffpulser)){
	    min = DiffTimeList[i];
	    break;
	  }else{
	    min = DiffTimeList[i+1];
	  }
	}
    }
  }
  return min;
}


void GATPulserTag::GetDeviateEvent(const vector<Int_t>& inputEntry, 
                                   const vector<Double_t>& PulserTime, 
                                   Double_t PulserDiffTime,
                                   vector<Int_t>& DeviateList)
{
  DeviateList.resize(0);

  if(PulserTime.size()>1){
    for(Int_t i = 0;i<(Int_t)PulserTime.size();i++){
      if(i == 0){
	Double_t ts = PulserTime[i];
	Double_t tj = PulserTime[i+1];
	
	Double_t ratio = TMath::Abs(tj-ts)/PulserDiffTime;
	Int_t factor = (Int_t) ratio;
	Double_t delta = TMath::Abs(TMath::Abs(tj-ts)-PulserDiffTime*factor);

        if(delta>PulserDiffTime*0.5){
          delta = TMath::Abs(TMath::Abs(tj-ts)-PulserDiffTime*(factor+1));
        }

	if( delta > 1e-2 && delta < PulserDiffTime*0.5 ){
	  DeviateList.push_back(inputEntry[i]);
	}
      }


      if(i>0 && i<(Int_t)PulserTime.size()-1){
	Double_t ti = PulserTime[i-1];
	Double_t ts = PulserTime[i];
	Double_t tj = PulserTime[i+1];
	Double_t ratio1 = TMath::Abs(tj-ts)/PulserDiffTime;
	Int_t factor1 = (Int_t) ratio1;
	Double_t delta1 = TMath::Abs(TMath::Abs(tj-ts)-PulserDiffTime*factor1);
	Double_t ratio2 = TMath::Abs(ts-ti)/PulserDiffTime;
	Int_t factor2 = (Int_t) ratio2;
	Double_t delta2 = TMath::Abs(TMath::Abs(ts-ti)-PulserDiffTime*factor2);

        if(delta1>PulserDiffTime*0.5){
          delta1 = TMath::Abs(TMath::Abs(tj-ts)-PulserDiffTime*(factor1+1));
        }
	if(delta2>PulserDiffTime*0.5){
          delta2 = TMath::Abs(TMath::Abs(ts-ti)-PulserDiffTime*(factor2+1));
        }

	if( delta1 > 1e-2 && delta2 > 1e-2){
	  DeviateList.push_back(inputEntry[i]);
	}
      }

      if(i == (Int_t)PulserTime.size()-1){
	Double_t ti = PulserTime[i-1];
	Double_t ts = PulserTime[i];
	Double_t ratio = TMath::Abs(ts-ti)/PulserDiffTime;
	Int_t factor = (Int_t) ratio;
	Double_t delta = TMath::Abs(TMath::Abs(ts-ti)-PulserDiffTime*factor);

        if(delta>PulserDiffTime*0.5){
          delta = TMath::Abs(TMath::Abs(ts-ti)-PulserDiffTime*(factor+1));
        }
	if( delta > 1e-2){
	  DeviateList.push_back(inputEntry[i]);
	}
      }
    }
  }
}

void GATPulserTag::PhysicalEventList(vector<Int_t>& peList)
{
  vector<Int_t> inputList;
  for(size_t i=0; i<fEntries; i++) inputList.push_back(i);
  vector<Int_t> taglist;
  Tag(inputList, taglist);
  vector<Int_t> pulserlist;
  GetPulserList(taglist, pulserlist); 
  
  vector<Double_t> ptime;
  vector<Int_t> pchannellist;
  vector<Int_t> devlist;
  Double_t PulserTimeDiff = 0;
  for(Int_t i =0;i<fChannels;i++){
    ptime.clear();
    pchannellist.clear();
    devlist.clear();
    if(fGoodBad[i]==1){
      GetPulserChannelList(pulserlist, i, pchannellist);
      if(pchannellist.size()>1){
        GetPulserTime(pchannellist, i, ptime);
        PulserTimeDiff = GetPulserDiffTime(ptime);
        GetDeviateEvent(pchannellist, ptime, PulserTimeDiff, devlist);
        vector<Int_t> newlist1;
        Tag3(devlist, newlist1);
        if(fRun >= 4549){
          vector<Int_t> newlist2;
          Tag2(newlist1, newlist2);
          TagCorrect(taglist, newlist2);
        }
        else TagCorrect(taglist, newlist1);
      }
    }
  }
    
  vector<Double_t> R1;
  for(size_t i = 0; i<taglist.size(); i++) {
    R1.push_back((Double_t)taglist[i]);
  }
  sort(R1.begin(), R1.end());
  Double_t templist = -1;
  for(Int_t i = 0;i<(Int_t)R1.size();i++){
    if (R1[i] != templist){
      peList.push_back((Int_t)R1[i]);
    }
    templist = R1[i];
  }
}

void GATPulserTag::PulserTree(const char* pathName)
{
  vector<Int_t> physicaleventlist;
  PhysicalEventList(physicaleventlist);

  vector<Int_t> taglist;
  GetTagList(physicaleventlist, taglist);

  TFile *newfile = new TFile(Form("%spulser_%d.root", pathName, fRun), "recreate");
  TTree *newtree = new TTree("pulsertree", "pulser tag");

  Int_t Pulser;
  vector<Int_t>* PulserTagChan = NULL;
  newtree->Branch("Pulser",&Pulser);
  //newtree->Branch("PulserTagChan",&PulserTagChan);
  for(size_t ie=0; ie<taglist.size(); ie++){
    PulserTagChan->clear();
    fMjdTree->GetEntry(ie);
    Pulser = taglist[ie];
    /*
    if(Pulser>0){
      for(size_t ichan = 0;ichan<fMTChannel->size();ichan++){
	Int_t tagchan = fMap->GetInt(fMTChannel->at(ichan), MJTChannelMap::kPulserTagChan);
	PulserTagChan->push_back(tagchan);
      }
    }else{
      for(size_t ichan = 0;ichan<fMTChannel->size();ichan++){
        PulserTagChan->push_back(0);
      }
    }
    */
    newtree->Fill();
  }
  
  //newtree->Print();
  newtree->Write();
  cout << "Pulser file is generated...." <<endl;
  delete newfile;
}



void GATPulserTag::SaveTree(TChain* fChain, Int_t IsPulser,string PathName, string FileName, string PulserPath){
  TTree *oldtree = (TTree*)fChain->GetTree();
  TFile *newfile = new TFile(Form("%s%s.root",PathName.c_str(),FileName.c_str()),"recreate");
  TTree *newtree = oldtree->CloneTree(0);

  TFile *pulserfile = new TFile(Form("%spulser_%d.root", PulserPath.c_str(), fRun), "read");
  TTree *pulsertree = (TTree*)pulserfile->Get("pulsertree");
  Int_t Pulser;
  pulsertree->SetBranchAddress("Pulser", &Pulser);

  for(size_t i=0;i<fEntries;i++){
    oldtree->GetEntry(i);
    pulsertree->GetEntry(i);
    if(Pulser == IsPulser){
      newtree->Fill();
    }
  }
  newtree->Print();
  newtree->AutoSave();
  delete newfile;
}

