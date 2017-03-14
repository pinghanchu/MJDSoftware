// 2016.11.21 Pinghan Chu
// Combine histgram of all high gain, all natural detectors, etc.
//
#include "GATAutoCal.hh"
#include "TFile.h"

int main(int argc, char** argv)
{
  if(argc != 10 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun number] [endrun number] [energy name] [Bin HG] [Low HG] [Up HG] [Bin LG] [Low LG] [Up LG]" << endl;
    return 1;
  }
  int startrun = atoi(argv[1]);
  int endrun = atoi(argv[2]);
  string EName = argv[3];
  Int_t BinHG = atoi(argv[4]);
  Int_t LowHG = atoi(argv[5]);
  Int_t UpHG = atoi(argv[6]);
  Int_t BinLG = atoi(argv[7]);
  Int_t LowLG = atoi(argv[8]);
  Int_t UpLG = atoi(argv[9]);
  GATAutoCal ac(startrun,startrun);
  vector<Int_t> channel = ac.GetChannel();
  vector<Int_t> cryo = ac.GetCryo();
  vector<Int_t> str = ac.GetString();
  vector<Int_t> detpos = ac.GetDetPosition();
  vector<Int_t> goodbad = ac.GetGoodBad();
  vector<Int_t> isenriched = ac.GetEnriched();

  string Pos;
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
    Pos = Form("%d%d%d",cryo.at(i),str.at(i),detpos.at(i));
    if((Int_t)i%2==0){
      Energy[i] = new TH1D(Form("%s%s",EName.c_str(),DataName.at(i).c_str()),Form("%s",TitleName.at(i).c_str()),BinHG,LowHG,UpHG);
    }else{
      Energy[i] = new TH1D(Form("%s%s",EName.c_str(),DataName.at(i).c_str()),Form("%s",TitleName.at(i).c_str()),BinLG,LowLG,UpLG);
    }
  }


  Int_t fPos;
  for(size_t  i =0;i<channel.size();i++){
    Pos = Form("%d%d%d",cryo.at(i),str.at(i),detpos.at(i));
    fPos = cryo.at(i)*10+str.at(i);
    TFile fhist(Form("./Hist/hist_%d_%d_%d.root",startrun,endrun,fPos),"read");
    cout << Pos.c_str() << " " << channel.at(i) << endl;
    TH1D* h1 = (TH1D*)fhist.Get(Form("%s%s%d",EName.c_str(),Pos.c_str(),channel.at(i)));
    if(goodbad.at(i)==1){
      if(channel.at(i)%2==0){
	Energy[4]->Add(h1);      
	if(isenriched.at(i) ==0){
	  Energy[0]->Add(h1);
	}else{
	  Energy[2]->Add(h1);
	}
      }else{
	Energy[5]->Add(h1);
        if(isenriched.at(i) ==0){
          Energy[1]->Add(h1);
        }else{
          Energy[3]->Add(h1);
        }
      }
    }
    delete h1;
    fhist.Close();
  }
 
  
  TFile newfile(Form("./Hist/hist_%d_%d.root",startrun,endrun),"update");
  for(size_t  i =0;i<DataName.size();i++){
    Energy[i]->Write();
  }
  
}


