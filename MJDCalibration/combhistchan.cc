// 2016.11.21 Pinghan Chu
// Combine histgram of each channel between startrun and endrun
//
#include "GATAutoCal.hh"
#include "TFile.h"

int main(int argc, char** argv)
{
  if(argc != 11 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [energy name] [Bin HG] [Low HG] [Up HG] [Bin LG] [Low LG] [Up LG] [channel]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  string fEName = argv[3];
  const char* EnergyName = Form("%s",fEName.c_str());
  Int_t BinHG = atoi(argv[4]);
  Int_t LowHG = atoi(argv[5]);
  Int_t UpHG = atoi(argv[6]);
  Int_t BinLG = atoi(argv[7]);
  Int_t LowLG = atoi(argv[8]);
  Int_t UpLG = atoi(argv[9]);
  Int_t fChan = atoi(argv[10]);

  GATAutoCal ds(fStartRun,fStartRun);
  ds.SetEnergyName(EnergyName);
  string DataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();

  Int_t pos;
  vector<Int_t> index;
  for(size_t i = 0;i<fChannel.size();i++){
    pos = fCryo.at(i)*10+fStr.at(i);
    if(fChannel.at(i) == fChan){
      index.push_back(i);
    }
  }
  const Int_t channels = index.size();
  TH1D *Energy[channels];

  string Pos;
  for(size_t i=0;i<index.size();i++){
    Int_t ii = index.at(i);
    Pos = Form("%d%d%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii));
    cout << ii << " " << fChannel.at(ii) << endl;
    if(fChannel.at(ii)%2==0){
      Energy[i] = new TH1D(Form("%s%s%d",fEName.c_str(),Pos.c_str(),fChannel.at(ii)),"",BinHG,LowHG,UpHG);
    }else{
      Energy[i] = new TH1D(Form("%s%s%d",fEName.c_str(),Pos.c_str(),fChannel.at(ii)),"",BinLG,LowLG,UpLG);
    }
  }

  for(Int_t irun = fStartRun;irun<=fEndRun;irun++){
    cout << "Add : hist_" << irun << ".root" << endl;
    TFile fhist(Form("./Hist/hist_%d.root",irun),"read");
    for(size_t  i =0;i<index.size();i++){
      Int_t ii = index.at(i);
      Pos = Form("%d%d%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii));
      TH1D* h = (TH1D*)fhist.Get(Form("%s%s%d",fEName.c_str(),Pos.c_str(),fChannel.at(ii)));
      Energy[i]->Add(h);
      delete h;
    }
    fhist.Close();
  }
   
  TFile newfile(Form("./Hist/hist_%d_%d.root",fStartRun,fEndRun),"update");
  for(size_t  i =0;i<index.size();i++){
    Int_t ii = index.at(i);
    cout << "Save : " << i << " " << fChannel.at(ii) << endl;
    Energy[i]->Write();
  }
  
}


