#include "GATAutoCal.hh"
#include "GATDataSet.hh"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH2F.h"


int main(int argc, char** argv)
{
  if(argc != 5 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [startrun] [endrun] [peak value] [peak window]" << endl;
    return 1;
  }
  Int_t fStartRun = atoi(argv[1]);
  Int_t fEndRun = atoi(argv[2]);
  Double_t fMeanValue = (Double_t) atoi(argv[3]);
  Double_t fWindow = (Double_t) atoi(argv[4]);
  
  GATAutoCal ds(fStartRun,fEndRun);
  const char* EnergyName = "trapENFCal";

  TChain *mjdTree = new TChain("mjdTree");
  //string MotherName = Form("./gatified/%s/mjd_run%d.root",fDataSet.c_str(),fRun);
  
  string MotherName;
  for(Int_t irun = fStartRun;irun<=fEndRun;irun++){
    MotherName= Form("./Hist/gat/mjd_run%d.root",irun);
    //string MotherName = Form("./mjd_run%d.root",fRun);
    if(ifstream(MotherName)){
      cout << "File mjd_run" << irun << ".root exists!" << endl;
      mjdTree->Add(MotherName.c_str());
    }else{
      cout << "File doesn't exist!" << endl;
    }
  }

  string PathName = "./Plot/";
  string FileName = Form("%d_%d",fStartRun,fEndRun);
  ds.SetEnergyName(EnergyName);
  vector<Int_t> fChannel = ds.GetChannel();
  vector<Int_t> fStr = ds.GetString();
  vector<Int_t> fDetpos = ds.GetDetPosition();
  vector<Int_t> fGoodBad = ds.GetGoodBad();
  ofstream fout(Form("./List/profiletime_%d_%d_%d_%d.txt",fStartRun,fEndRun,atoi(argv[3]),atoi(argv[4])),ios::app);
  fout.precision(15);

  string fCut;
  string temp;
  for(Int_t i = 0;i<(int)fChannel.size();i++){
    if(fGoodBad.at(i)>0 && fChannel.at(i)%2==0){
      temp = Form("channel == %d ||",fChannel.at(i));
      fCut.append(temp);
    }
  }
  fCut.erase(fCut.size()-1);
  fCut.erase(fCut.size()-1);
  cout << fCut.c_str() << endl;
  TProfile *p1 = ds.PeakProfileTime(mjdTree,fCut,fMeanValue,fWindow,PathName,FileName);
  Double_t sum = 0;
  Double_t ave = 0;
  Double_t count = 0;
  Int_t entries = p1->GetEntries();
  cout << entries << endl;
  for(Int_t j=0;j< entries;j++){
    if(p1->GetBinContent(j+1)>0){
      sum = sum + p1->GetBinContent(j+1);
      count++;
    }
  }
  ave = sum/count;
  
  for(Int_t j=0;j< entries;j++){
    if(p1->GetBinContent(j+1)>0){
      fout << j << " " << p1->GetBinContent(j+1) << " " << p1->GetBinError(j+1)<< " " << ave << " " << p1->GetBinContent(j+1)/ave << endl;
    }
  }
  //delete p1;
 
}


