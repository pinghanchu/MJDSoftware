#include "GATAutoCal.hh"

int main(int argc, char** argv)
{
  if(argc != 4) {
    cout << "Usage: " << argv[0] << "[dataset] [startrun number] [endrun number]" << endl;
    return 1;
  }
  Int_t fDataSet = atoi(argv[1]);
  int startrun = atoi(argv[2]);
  int endrun = atoi(argv[3]);
  GATAutoCal ac1(startrun,startrun);
  string dataset = ac1.GetDataSet();
  Int_t index = endrun;
  Int_t ifExist = 0;

  if(ifstream(Form("./built/%s/OR_run%d.root",dataset.c_str(),index))){
    ifExist = 1;
  }else{
    ifExist = 0;
  }

  while(ifExist == 0){
    index--;
    if(ifstream(Form("./built/%s/OR_run%d.root",dataset.c_str(),index))){
      ifExist = 1;
    }else{
      ifExist = 0;
    }
  }
  GATAutoCal ac2(index,index);
  //cout << ac.GetGATRev() << endl;;
  ofstream fout(Form("./List/gettime_DS%d.txt",fDataSet),ios::app);
  Double_t starttime = ac1.GetStartTime();
  Double_t stoptime = ac2.GetStopTime();
  //Int_t index = endrun;
  //Int_t ifExist = 0;
  while(stoptime == 0 && index>startrun){
    index--;

    if(ifstream(Form("./built/%s/OR_run%d.root",dataset.c_str(),index))){
      //      ifExist = 1;
      GATAutoCal aci(index,index);
      stoptime = aci.GetStopTime();
    }else{
      //ifExist = 0;
    }
  }
  fout << startrun << " " << endrun << " " << index << " "<< starttime << " " << stoptime <<  " " << stoptime-starttime << endl;
}


