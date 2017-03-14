// 2016.11.21 Pinghan Chu
// Combine histgram of all high gain, all natural detectors, etc.
//
#include "GATAutoCal.hh"
#include "TFile.h"

int main(int argc, char** argv)
{
  if(argc != 9 ) {
    cout << "Usage: " << argv[0] << " [energy name] [Bin HG] [Low HG] [Up HG] [Bin LG] [Low LG] [Up LG] [List Path]" << endl;
    return 1;
  }
  string EName = argv[1];
  Int_t BinHG = atoi(argv[2]);
  Int_t LowHG = atoi(argv[3]);
  Int_t UpHG = atoi(argv[4]);
  Int_t BinLG = atoi(argv[5]);
  Int_t LowLG = atoi(argv[6]);
  Int_t UpLG = atoi(argv[7]);
  string PathName = argv[8];
  ifstream fin(PathName);
  Int_t run1,run2,run3,run4;
  vector<Int_t> Run1;
  vector<Int_t> Run2;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run1 >> run2 >> run3>>run4;
      Run1.push_back(run1);
      Run2.push_back(run2);
    }
  }
  Run1.pop_back();
  Run2.pop_back();

  vector<string> DataName;
  DataName.push_back("00"); // Natural HG
  DataName.push_back("01"); // Natural LG
  DataName.push_back("10"); // Enriched HG
  DataName.push_back("11"); // Enriched LG
  DataName.push_back("HG");
  DataName.push_back("LG");
  vector<string> TitleName;
  TitleName.push_back("Natural HG");
  TitleName.push_back("Natural LG");
  TitleName.push_back("Enriched HG");
  TitleName.push_back("Enrichd LG");
  TitleName.push_back("High Gain");
  TitleName.push_back("Low Gain");
  const Int_t nDataName = DataName.size();
  TH1D* Energy[nDataName];
  for(size_t i=0;i<DataName.size();i++){
    if((Int_t)i%2==0){
      Energy[i] = new TH1D(Form("%s%s",EName.c_str(),DataName.at(i).c_str()),Form("%s",TitleName.at(i).c_str()),BinHG,LowHG,UpHG);
    }else{
      Energy[i] = new TH1D(Form("%s%s",EName.c_str(),DataName.at(i).c_str()),Form("%s",TitleName.at(i).c_str()),BinLG,LowLG,UpLG);
    }
  }

  for(size_t j = 0;j<Run1.size();j++){
    TFile fhist(Form("./Hist/hist_%d_%d.root",Run1.at(j),Run2.at(j)),"read");
    cout << Run1.at(j) << " " << Run2.at(j) << endl;
    for(size_t  i =0;i<DataName.size();i++){
      TH1D* h1 = (TH1D*)fhist.Get(Form("%s%s",EName.c_str(),DataName.at(i).c_str()));
      Energy[i]->Add(h1);
      delete h1;
    }
    fhist.Close();
  } 
  TFile newfile(Form("./Hist/hist_%d_%d.root",Run1.at(0),Run2.at(Run2.size()-1)),"update");
  for(size_t  i =0;i<DataName.size();i++){
    Energy[i]->Write();
  }
  
}


