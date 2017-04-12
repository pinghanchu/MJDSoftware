// 2016.11.21 Pinghan Chu
// Combine histgram of each channel between startrun and endrun
//
#include "GATAutoCal.hh"
#include "TFile.h"

int main(int argc, char** argv)
{
  if(argc != 6 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [energy name] [Input Path] [Output File]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  string fEName = argv[3];
  const char* fEnergyName = Form("%s",fEName.c_str());
  string fInputPath = argv[4];
  string fOutputFile = argv[5];

  GATAutoCal ds(fStartRun,fStartRun);
  ds.SetEnergyName(fEnergyName);
  string fDataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();

  vector<Int_t> index;
  for(size_t i = 0;i<fChannel.size();i++){
    if(fGoodBad.at(i) == 1){
      index.push_back(i);
    }
  }

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
  }else if(fEName == "trapENFBL" || (fStartRun>= 4171 && fEndRun<=4201) || (fStartRun>=9913 && fEndRun<= 9926) || (fStartRun>= 13071 && fEndRun<=13074) || (fStartRun>= 60002368 && fEndRun<= 60002372) ){
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

  const Int_t channels = index.size();
  TH1D* Enr[channels];
  Int_t Bin;
  Double_t Low,Up;
  for(size_t i=0;i<index.size();i++){
    Enr[i] = NULL;
    Int_t ii = index.at(i);
    Int_t pos = fCryo.at(ii)*100+fStr.at(ii)*10+fDetpos.at(ii);

    if(ii%2 == 0){
      Bin = Bin1;
      Low = Low1;
      Up  = Up1;
    }else{
      Bin = Bin2;
      Low = Low2;
      Up  = Up2;
    }
    string fHistoName = Form("%s%d",fEName.c_str(),fChannel.at(ii));
    TH1D *h = new TH1D(Form("%s",fHistoName.c_str()), Form("C%dP%dD%d, Channel = %d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii),fChannel.at(ii)),Bin,Low,Up);
    string fInputFile = Form("%shist_%d_%d_%s_%d.root",fInputPath.c_str(),fStartRun,fEndRun,fEName.c_str(),pos);
    TH1D *h1 = ds.LoadHisto(fInputFile,fHistoName);
    h->Add(h1);
    delete h1;

    Enr[i] = (TH1D*)h->Clone();
    delete h;
  }
  TFile newfile(Form("%s",fOutputFile.c_str()),"update");
  for(size_t i = 0;i<index.size();i++){
    Enr[i]->Write();
  }
}


