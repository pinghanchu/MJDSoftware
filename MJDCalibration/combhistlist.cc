// 2016.11.21 Pinghan Chu
// Combine histgram of each channel from the run list in ./List/cal.list.txt
//
#include "GATAutoCal.hh"
#include "TFile.h"

int main(int argc, char** argv)
{
  if(argc != 9 ) {
    cout << "Usage: " << argv[0] << " [ename] [Bin HG] [Low HG] [Up HG] [Bin LG] [Low LG] [Up LG] [List Path] " << endl;
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
  //ifstream fin("./List/runlist/cal.DS2.txt");
  ifstream fin(PathName);
  Int_t run1,run2,r3,r4;
  vector<Int_t> Run1;
  vector<Int_t> Run2;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> run1 >> run2 >> r3 >> r4;
      Run1.push_back(run1);
      Run2.push_back(run2);
    }
  }
  Run1.pop_back();
  Run2.pop_back();

  GATAutoCal ac(Run1.at(0),Run1.at(0)+1);
  vector<Int_t> channel = ac.GetChannel();
  vector<Int_t> cryo = ac.GetCryo();
  vector<Int_t> str = ac.GetString();
  vector<Int_t> detpos = ac.GetDetPosition();
  const Int_t channels = channel.size();
  TH1D* Energy[channels];
  string Pos;
  for(size_t i=0;i<channel.size();i++){
    Pos = Form("%d%d%d",cryo.at(i),str.at(i),detpos.at(i));
    if(channel.at(i)%2==0){
      Energy[i] = new TH1D(Form("%s%s%d",EName.c_str(),Pos.c_str(),channel.at(i)),"",BinHG,LowHG,UpHG);
    }else{
      Energy[i] = new TH1D(Form("%s%s%d",EName.c_str(),Pos.c_str(),channel.at(i)),"",BinLG,LowLG,UpLG);
    }
  }
  for(size_t j = 0;j<Run1.size();j++){
    TFile fhist(Form("./Hist/hist_%d_%d.root",Run1.at(j),Run2.at(j)),"read");
    for(size_t  i =0;i<channel.size();i++){
      Pos = Form("%d%d%d",cryo.at(i),str.at(i),detpos.at(i));
      TH1D* h1 = (TH1D*)fhist.Get(Form("%s%s%d",EName.c_str(),Pos.c_str(),channel.at(i)));
      Energy[i]->Add(h1);      
      delete h1;
    }
    fhist.Close();
  }
  
  
  TFile newfile(Form("./Hist/hist_%d_%d.root",Run1.at(0),Run2.at(Run2.size()-1)),"update");
  for(size_t  i =0;i<channel.size();i++){
    Energy[i]->Write();
  }
  
}


