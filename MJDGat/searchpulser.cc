#include "MJDGat.hh"
#include "TStyle.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2 || atoi(argv[1]) == 0) {
    cout << "Usage: " << argv[0] << " [run]" << endl;
    return 1;
  }

  Int_t fRun = atoi(argv[1]);
  MJDGat ds(fRun);
  ds.SearchPulserEvents();

}
