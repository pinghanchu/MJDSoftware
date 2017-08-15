#include "GATPulserTag.hh"
#include "TCanvas.h"

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [runnumber]" << endl;
    return 1;
  }
  
  Int_t run = atoi(argv[1]);
  const char* energyName = "trapE";
  const char* fileName = "";
  GATPulserTag pt(run);

  //Int_t isradio = pt.GetIsRadio();
  //  if(isradio == 1){
  //  cout << "calibration source run... quitting." << endl;
  //  return 0;
  //}else{
  pt.SetEnergyName(energyName);
  pt.PulserTree(fileName);
    //}
  return 0;
}


