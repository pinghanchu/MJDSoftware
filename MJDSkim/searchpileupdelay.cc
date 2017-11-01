#include "MJDSkim.hh"
#include "GATAutoCal.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 3) {
    cout << "Usage: " << argv[0] << " [input file] [output]" << endl;
    return 1;
  }
  string fInputFile = argv[1];
  string fOutputFile = argv[2];

  ifstream fin(Form("%s",fInputFile.c_str()));
  Int_t dataset,subset,list,run,event,chan;
  Double_t time, enr,ratio,deltaT,avse,dcr,trapetailmin;
  vector<Int_t> DataSet;
  vector<Int_t> SubSet;
  vector<Int_t> List;
  vector<Int_t> Chan;
  if(fin.is_open()){
    while(!fin.eof()){
      fin >> dataset >> subset >> list >> run >> event >> time >> chan >> enr >> ratio >> deltaT >> avse >> dcr >> trapetailmin;
      DataSet.push_back(dataset);
      SubSet.push_back(subset);
      List.push_back(list);
      Chan.push_back(chan);
    }
  }
  for(size_t i =0;i<DataSet.size();i++){
    MJDSkim ds(DataSet.at(i),SubSet.at(i),SubSet.at(i),0);
    ds.FindDelayedEvent(List.at(i),Chan.at(i),0.5,fOutputFile);
  }
  return 1;
}
