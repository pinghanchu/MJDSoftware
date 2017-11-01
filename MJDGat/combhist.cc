// 2016.11.21 Pinghan Chu
// Combine histgram of each channel between startrun and endrun
//
#include "MJDGat.hh"
#include "TFile.h"
#include "TROOT.h"
int main(int argc, char** argv)
{
  if(argc != 3 ) {
    cout << "Usage: " << argv[0] << "[energy name] [Input File]" << endl;
    return 1;
  }
  string fEName = argv[1];
  //const char* fEnergyName = Form("%s",fEName.c_str());
  string fInputFile = argv[2];
  string run,ds;
  ifstream fin(fInputFile);
  vector<Int_t> Run;
  vector<Int_t> DS;
  while(getline(fin,ds,',')){
    getline(fin,run,'\n');
    int run1 = stoi(run);
    int ds1  = stoi(ds);
    Run.push_back(run1);
    DS.push_back(ds1);
  }



  Int_t Bin = 300000;
  Double_t Low = 0;
  Double_t Up = 3000;
  
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
  
  const Int_t sums = DataName.size();
  TH1D* Sum[sums];
  
  for(size_t i=0;i<DataName.size();i++){
    Sum[i] = new TH1D(Form("%s%s",fEName.c_str(),DataName.at(i).c_str()),Form("%s",TitleName.at(i).c_str()),Bin,Low,Up);
  }

  /*
  for(size_t i = 0;i<fChannel.size();i++){
    if(fGoodBad.at(i) == 1){
      if(fChannel.at(i)%2==0){
	Sum[4]->Add(Enr[i]);
	if(fIsEnriched.at(i) ==0){
	  Sum[0]->Add(Enr[i]);
	}else{
	  Sum[2]->Add(Enr[i]);
	}
      }else{
	Sum[5]->Add(Enr[i]);
	if(fIsEnriched.at(i) ==0){
	  Sum[1]->Add(Enr[i]);
	}else{
	  Sum[3]->Add(Enr[i]);
	}
      }
    }
  }
  */
  
  for(size_t i=0;i<DataName.size();i++){
    for(size_t j = 0;j<Run.size();j++){
      cout <<DataName.at(i) << " "<<  Run.at(j) << endl;
      TFile fhist(Form("./data/hist_%d_%d.root",DS.at(j),Run.at(j)),"read");
      
      TH1D* h1 = (TH1D*)fhist.Get(Form("%s%s",fEName.c_str(),DataName.at(i).c_str()));
      Sum[i]->Add(h1);
      delete h1;
      fhist.Close();
    }
  }
  
  string fOutputFile = "alldata.root";
  TFile newfile(Form("%s",fOutputFile.c_str()),"update");
  for(size_t i = 0;i<DataName.size();i++){
    Sum[i]->Write();
  }
}


