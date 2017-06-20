#include "MJDGat.hh"
#include "MJTChannelMap.hh"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run number]" << endl;
    return 1;
  }
  Int_t run = atoi(argv[1]);
  MJDGat ds(run);
  
  MJTChannelMap* map = ds.GetMap();
  vector<Int_t> channel = ds.GetChannel();
  vector<Int_t> cryo = ds.GetCryo();
  vector<Int_t> str = ds.GetString();
  vector<Int_t> position = ds.GetDetPosition();
  vector<Int_t> goodbad = ds.GetGoodBad();
  vector<string> detname = ds.GetDetectorName();
  vector<Int_t> pulserchannel = ds.GetPulserChannel();
  map->Print();
  //cout << channel.size() << " " << cryo.size() << endl;

  for(size_t i=0;i<channel.size();i++){
    cout << cryo.at(i)<<str.at(i)<<position.at(i) << " " << detname.at(i).c_str() << " " << channel.at(i) << endl;
  }
  for(size_t i=0;i<pulserchannel.size();i++){
    cout << "Pulser " << pulserchannel.at(i) << endl;
  }

}


