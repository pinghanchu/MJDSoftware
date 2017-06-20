#include "GATNonLinearityCorrector.hh"
#include "MJTGretina4DigitizerData.hh"
#include "MJTypes.hh"

using namespace std;

GATNonLinearityCorrector::GATNonLinearityCorrector(
  MGWFNonLinearityCorrector& nonLinCorrTransformer,
  const char* nlcDBPath, const char* postedWFNameToProcess, 
  const char* wfNameForPosting, bool useAuxWFs, bool skipBadWFs
) : 
  GATWaveformTransformer(
    nonLinCorrTransformer, postedWFNameToProcess, wfNameForPosting, useAuxWFs, skipBadWFs
  ), 
  fNLCTransformer(&nonLinCorrTransformer),
  fNLCMapDir(nlcDBPath)
{ }

void GATNonLinearityCorrector::LoadParameters(int ddID, int run)
{
  string path = fNLCMapDir;
  
  MGWFNonLinearityCorrectionMap* nlcMap = fNLCMaps[ddID];
  //This is a hack. Eventually we need to have these maps in a DB and pull them out properly.//

  if(nlcMap == NULL) {
    nlcMap = new MGWFNonLinearityCorrectionMap;
    int crate = MJUtil::GetCrate(ddID);
    int card = MJUtil::GetCard(ddID);
    int channel = MJUtil::GetChannel(ddID);
    if(run < 11339) path += "/Run0";
    else if(run < 11397) path += "/Run11339";
    else if(run < 18622) path += "/Run11397";
    else if(run < 18643) path += "/Run18623";
    else if(run < 18990) path += "/Run18643";
    else if(run < 19018) path += "/Run18990";
    else if(run < 19502) path += "/Run19018";
    else if(run < 6000000) path += "/Run19502";
    else if(run < 60000000) {
      //cout << "No NLC files for run " << run << endl;      
      return;
    }
    else if(run < 60002395) path += "/Run60000000";
    else if(run < 60002397) path += "/Run60002395";
    else if(run < 60002419) path += "/Run60002397";
    else {
      //cout << "No NLC files for run " << run << endl;      
      return;
    }
    if(crate == 1) {
      char upFileName[500];
      sprintf(upFileName,"%s/Crate%d_GRET%d_Ch%d_up1.dat", path.c_str(), crate, card, channel);
      char downFileName[500];
      sprintf(downFileName,"%s/Crate%d_GRET%d_Ch%d_up0.dat", path.c_str(), crate, card, channel);
      nlcMap->LoadFromUpDownFiles(upFileName, downFileName);
    }
    else if(crate == 2) {
      char combFileName[500];
      sprintf(combFileName,"%s/Crate%d_GRET%d_Ch%d_comb.dat", path.c_str(), crate, card, channel);
      nlcMap->LoadFromCombinedFile(combFileName);
    }
    else {
      cout << "GATNonLinearityCorrector::LoadParameters(" 
           << ddID << ", " << run << "): Error: got invalide crate number " 
           << crate << endl;
    }
    fNLCMaps[ddID] = nlcMap;
  }
  fNLCTransformer->SetNLCMap(nlcMap);
}

