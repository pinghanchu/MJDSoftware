// 2016.11.21 Pinghan Chu
// Fill histograms from the files on PDSF using GATDataSet();

#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TFile.h"
int main(int argc, char** argv)
{
  if(argc != 10 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run] [energy name] [Bin HG] [Low HG] [Up HG] [Bin LG] [Low LG] [Up LG] [Position = Mod*10+String]" << endl;
    return 1;
  }
  Int_t fRun = atoi(argv[1]);
  string fEName = argv[2];
  const char* EnergyName = Form("%s",fEName.c_str());
  Int_t BinHG = atoi(argv[3]);
  Int_t LowHG = atoi(argv[4]);
  Int_t UpHG = atoi(argv[5]);
  Int_t BinLG = atoi(argv[6]);
  Int_t LowLG = atoi(argv[7]);
  Int_t UpLG = atoi(argv[8]);
  Int_t fPos = atoi(argv[9]);

  string FileName = Form("./Hist/hist_%d_%d.root",fRun,fPos);  

  GATAutoCal ds(fRun,fRun);
  ds.SetEnergyName(EnergyName);
  string fDataSet = ds.GetDataSet();
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fCryo = ds.GetCryo();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();

  TChain *mjdTree = new TChain("mjdTree");
  //string MotherName = Form("./gatified/%s/mjd_run%d.root",fDataSet.c_str(),fRun);
  string MotherName = Form("./Hist/gat/mjd_run%d.root",fRun);
  //string MotherName = Form("./mjd_run%d.root",fRun);
  if(ifstream(MotherName)){
    cout << "File mjd_run" << fRun << ".root exists!" << endl;
    mjdTree->Add(MotherName.c_str());
  }else{
    cout << "File doesn't exist!" << endl;
  }
  Int_t pos;
  vector<Int_t> index;
  for(size_t i = 0;i<fChannel.size();i++){
      pos = fCryo.at(i)*10+fStr.at(i);
      if(pos == fPos){
	index.push_back(i);
      }
  }
  const Int_t channels = index.size();
  TH1D *Energy[channels];

  Int_t Bin1 = BinHG;
  Int_t Bin2 = BinLG;
  Double_t Low1 = LowHG;
  Double_t Low2 = LowLG;
  Double_t Up1 = UpHG;
  Double_t Up2 = UpLG;
  Int_t Bin;
  Double_t Low = 0;
  Double_t Up = 0;
  string HistoName;
  string TitleName;
  string Pos;


  for(size_t i=0;i<index.size();i++){
    Int_t ii = index.at(i);
    Pos = Form("%d%d%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii));
    HistoName = Form("%s%s%d",fEName.c_str(), Pos.c_str(), fChannel.at(ii));
    TitleName = Form("C%dP%dD%d, channel =%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii), fChannel.at(ii));
    if(ii%2 == 0){
      Bin = Bin1;
      Low = Low1;
      Up  = Up1;
    }else{
      Bin = Bin2;
      Low = Low2;
      Up  = Up2;
    }
    
    TH1D *h = new TH1D(Form("%s",HistoName.c_str()),Form("%s",TitleName.c_str()),Bin, Low, Up);      
    TCut cut1;
    size_t found = fEName.find("trapENFBL"); 
 
    if((fRun >= 9913 && fRun <= 9926) || (fRun >= 13071 && fRun <= 13074)){
      cut1 = Form("channel == %d && (rawWFMax > 200 || rawWFMin < -20)",fChannel.at(ii));
    }else if(fRun >= 60002368 && fRun <= 60002372 && fChannel.at(ii)==1175){
      cut1 = Form("channel == %d && (rawWFMax > 50 || rawWFMin < -200)",fChannel.at(ii));
    }else if(fRun >= 60002368 && fRun <= 60002372 && fChannel.at(ii)==1106){
      cut1 = Form("channel == %d && (rawWFMax > 200 || rawWFMin < -200)",fChannel.at(ii));
    }else if(fRun >= 60002368 && fRun <= 60002372 && fChannel.at(ii)==1107){
      cut1 = Form("channel == %d && (rawWFMax > 80 || rawWFMin < -50)",fChannel.at(ii));
    }else if(fRun >= 60002368 && fRun <= 60002372){
      cut1 = Form("channel == %d && (rawWFMax > 400 || rawWFMin < -400)",fChannel.at(ii));
    }else if(fRun == 4171 && fRun == 4201){
      cut1 = Form("channel == %d",fChannel.at(ii));
    }else{
      cut1 = Form("channel == %d",fChannel.at(ii));
    }
    mjdTree->Draw(Form("%s>>%s",fEName.c_str(),HistoName.c_str()),cut1);
    Energy[i] = (TH1D*)h->Clone();
    delete h;
  }

  TFile fhist(Form("%s",FileName.c_str()),"update");
  
  for(size_t i =0;i<index.size();i++){
    Int_t ii = index.at(i);
    Pos = Form("%d%d%d",fCryo.at(ii),fStr.at(ii),fDetpos.at(ii));
    cout << "Save : " << fEName.c_str() << " " << Pos.c_str() << " " << fChannel.at(ii) << endl;
    Energy[i]->Write();
  }

}


