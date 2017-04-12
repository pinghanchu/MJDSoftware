// 2016.11.21 Pinghan Chu
// Fill histograms from the files on PDSF using GATDataSet();

#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TFile.h"
int main(int argc, char** argv)
{
  if(argc != 6 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run] [energy name] [Position = Mod*100+String*10+Det Position] [Input File] [Output File]" << endl;
    return 1;
  }

  Int_t fRun = atoi(argv[1]);
  string fEName = argv[2];
  const char* fEnergyName = Form("%s",fEName.c_str());
  Int_t fPos = atoi(argv[3]);
  string fInputFile = argv[4];
  string fOutputFile = argv[5];  

  GATAutoCal ds(fRun,fRun);
  ds.SetEnergyName(fEnergyName);
  string fDataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();

  TChain *mjdTree = new TChain("mjdTree");
  if(ifstream(fInputFile)){
    cout << "File " << fInputFile.c_str() << " exists!" << endl;
    mjdTree->Add(fInputFile.c_str());
  }else{
    cout << "File doesn't exist!" << endl;
  }

  vector<Int_t> index;
  for(size_t i = 0;i<fChannel.size();i++){
    Int_t pos = fCryo.at(i)*100+fStr.at(i)*10+fDetpos.at(i);
    if(pos == fPos){
      index.push_back(i);
    }
  }
  //const Int_t channels = index.size();

  Int_t Bin1,Bin2,Low1,Low2,Up1,Up2;
  
  if(fEName == "trapE" || fEName == "trapENM" || fEName == "trapENF"){
    Bin1 = 800000;
    Bin2 = 300000;
    Low1 = 0;
    Low2 = 0;
    Up1 = 8000;
    Up2 = 3000;
  }else if(fEName == "trapECal" || fEName == "trapENMCal" || fEName == "trapENFCal"){
    Bin1 = 300000;
    Bin2 = 300000;
    Low1 = 0;
    Low2 = 0;
    Up1 = 3000;
    Up2 = 3000;
  }else if(fEName == "trapENFBL" || (fRun>= 4171 && fRun<=4201) || (fRun>=9913 && fRun<= 9926) || (fRun>= 13071 && fRun<=13074) || (fRun>= 60002368 && fRun<= 60002372) ){
    Bin1 = 1000;
    Bin2 = 1000;
    Low1 = -10;
    Low2 = -10;
    Up1 = 10;
    Up2 = 10;
  }else{
    Bin1 = 0;
    Bin2 = 0;
    Low1 = 0;
    Low2 = 0;
    Up1 = 0;
    Up2 = 0;
    cout << "No this energy parameter!" <<endl;
  }
    
  Int_t Bin;
  Double_t Low = 0;
  Double_t Up = 0;

  TFile fhist(Form("%s",fOutputFile.c_str()),"update");

  for(size_t i=0;i<index.size();i++){
    Int_t ii = index.at(i);
    if(ii%2 == 0){
      Bin = Bin1;
      Low = Low1;
      Up  = Up1;
    }else{
      Bin = Bin2;
      Low = Low2;
      Up  = Up2;
    }    
    TH1D *h = ds.FillHisto(mjdTree, fEName, fChannel.at(ii),Bin,Low,Up);
    h->Write();
    delete h;
  }


}


