// 2016.11.18 Pinghan Chu
// 1. Modify from the original GAT processor.
// 2. Retire most energy parameters. Energy parameters should be obtained in a separate processor.
// 3. Fix the bug for the noise power spectrum(Micah Buuck) and the pile-up event finder(Julieta Gruszko).
//
// 2016.12.20 Pinghan Chu
// 1. Retire rawWFftSlope, Offset, etc.
// 2. Retire rawWFwholeWFSlope, Offset, etc.
// 3. Keep only 50% and 97% point of waveforms.
// 4. Change trapENM(4-2.5-4). Add trapENMSample, trapENMMin, trapENMCal. Retire trapEMin, trapESample.
// 5. Add trapETailMin
// 6. Rename trirt100nsft10nsMax to triTrapMax and trirt100nsft10nsIntegralW to triTrapIntegralW. Retire triFilMin and toe.
// 7. Add pile-up parameters: fastTrapNLCWFsnRisingX, pileUpTrapNLCWFsnRisingX.
// 2017.1.31 Pinghan Chu
// 1. Add "RawWFftSlope" back
// 2. Add 1%, 10%, 90% point of waveforms
// 3. Change trapENM(4-2-4) and trapE(4-2-4)
//
// 2017.02.28 Clint Wiseman
// Add GATThresholdReader (takes a ROOT input file)
// Add a check that pulser and threshold files exist
// 2017.04.17 Pinghan Chu
// Keep parameters using database

#include <ctime> //clock_t
#include "TFile.h"
#include "TTree.h"
#include "GATWaveformTransformer.hh"
#include "GATWFExtremumProcessor.hh"
#include "GATPostedDataTreeWriter.hh"
#include "GATWFTimePointProcessor.hh"
#include "MGWFBaselineRemover.hh"
#include "GATWFLinearFitProcessor.hh"
#include "MGWFBySampleDerivative.hh"
#include "MGWFTrapSlopeFilter.hh"
#include "MGWFSmoother.hh"
#include "MGPCIntegral.hh"
#include "MGPCNoiseBands.hh"
#include "GATWFTransientNoiseTagger.hh"
#include "MGPC2ndDiffPeakValue.hh"
#include "MGPCNorm2ndDiffPeakValue.hh"
#include "MGPCSingleSampleFilter.hh"
#include "MGPCSingleSampleFilterSigned.hh"
#include "MGWFRCDifferentiation.hh"
#include "MGWFTrapezoidalFilter.hh"
#include "MGWFResidualTailSubtractor.hh"
#include "MGPCCountThresholdCrossings.hh"
//#include "GATEnergyCalibrationAscii.hh"
#include "GATEnergyCalibration.hh"
#include "GATRunSummary.hh"
#include "GATMultiplicityProcessor.hh"
#include "MJTChannelMap.hh"
#include "MJTChannelSettings.hh"
#include "GATPostedDataExpressionProcessor.hh"
#include "GATWFCalcParamProcessor.hh"
#include "GATWFCountCrossingsParamProcessor.hh"
#include "GATWFDCRParamProcessor.hh"
#include "GATWFDataCleaningProcessor.hh"
#include "GATDCBoxCut.hh"
#include "GATEventDataCleaningProcessor.hh"
#include "GATDCPulserFileTagger.hh"
#include "GATDCAoverETagger.hh"
#include "GATDCSlopeDeltaY.hh"
#include "GATDetInfoProcessor.hh"
#include "GATInputDataPoster.hh"
#include "MGWFSlopeCalculator.hh"
#include "MGDCRSlopeCalculator.hh"
#include "GATEventTimeInfoProcessor.hh"
#include "MJAnalysisDoc.hh"
#include "MJAnalysisParameters.hh"
#include <map>
#include "GATNonLinearityCorrector.hh"
#include "MGWFSavitzkyGolaySmoother.hh"
#include "MGWFStartTimeFinder.hh"
#include "GATWFFixedTimePickoffProcessor.hh"
#include "MGWFWindower.hh"
#include "MGWFMovingAverage.hh"
#include "GATUtils.hh"
#include "GATTrapezoidalFilter.hh"
#include "MJTRun.hh"
#include "MJTMSWaveform.hh"
#include "GATMSWFPoster.hh"
#include "GATThresholdReader.hh"

using namespace std;
using namespace CLHEP;
using namespace MJDB;

int main(int argc, char** argv)
{
  // Give usage info for no input (user must supply an input file)
  if(argc < 2 || argc > 2) {
    cout << "Usage: " << argv[0] << "[gat file]" << endl;
    return 1;
  }

  // start timer
  clock_t timer;
  timer = clock();

  //string calibrationMapFile = "";
  //if(argc == 4) calibrationMapFile = argv[3];

  // *******************************
  // Pull stuff out from input files
  // *******************************

  // open input file
  TFile* gatfile = TFile::Open(argv[1]);
  if(gatfile == NULL) return 0;

  // find mjdTree
  TTree* mjdTree = (TTree*) gatfile->Get("mjdTree");
  if(mjdTree == NULL) {
    cout << "mjdTree not found. Continue" <<endl;
  }
  
  //int runNumber = RunSettings->GetRunNumber();

  // **********************************
  // Set up data processing environment
  // **********************************

  // Set up access to Analysis Parameters Database
  // Processors can use this pointer for common access.
  // No constructor args means we are communicating with the default APDB
  MJAnalysisDoc* analysisDoc = new MJAnalysisDoc();
  //analysisDoc->SetAlternateDB("mjdb.phy.ornl.gov", "mjdbsandbox", "6984", "https", "mjd"); //use the sandbox DB

  // prepare selector and writer so we can add to them as we go
  TAMSelector selector;
  GATPostedDataTreeWriter mjdDataWriter("enr"); //the name mjd is also used at the end

  // ****************************************************
  // Copy stuff from the input to the output and intercom
  // ****************************************************

  // Add event type and timestamp to the intercom
  
  GATInputDataPoster inputDataPoster;
  inputDataPoster.PostInputLeafVector("trapE","trapE");
  inputDataPoster.PostInputLeafVector("trapENF","trapENF");
  inputDataPoster.PostInputLeafVector("trapENM","trapENM");
  selector.AddInput(&inputDataPoster);
  mjdDataWriter.AddInputLeafVector("trapE", "trapE");
  mjdDataWriter.AddInputLeafVector("trapENF", "trapENF");
  mjdDataWriter.AddInputLeafVector("trapENM", "trapENM");

  // ******************
  // Energy calibration
  // ******************

  // Energy calibration, "trapECal", "trapEMinCal", "trapENFCal"

  GATEnergyCalibration trapECalibration(analysisDoc, "trapE");
  trapECalibration.UseCalibrationFor("trapE");
  selector.AddInput(&trapECalibration);
  mjdDataWriter.AddPostedVector(trapECalibration.GetNameOfPostedVector());

  GATEnergyCalibration trapENMCalibration(analysisDoc, "trapENM");
  trapENMCalibration.UseCalibrationFor("trapENM"); // <- Need to update once the parameters uploaded to the database
  selector.AddInput(&trapENMCalibration);
  mjdDataWriter.AddPostedVector(trapENMCalibration.GetNameOfPostedVector());

  GATEnergyCalibration trapENFCalibration(analysisDoc, "trapENF");
  trapENFCalibration.UseCalibrationFor("trapENF");
  selector.AddInput(&trapENFCalibration);
  mjdDataWriter.AddPostedVector(trapENFCalibration.GetNameOfPostedVector());

  // ************
  // Now process!
  // ************

  selector.AddInput(&mjdDataWriter);
  mjdTree->Process(&selector);


  // *************************************
  // Finally, clean up and write to output
  // *************************************

  //  if(MCATree != NULL || ChanMap != NULL || ChanSettings != NULL){
  if(mjdTree != NULL){
    string inputFileName = gatfile->GetName();
    if(inputFileName.find("/") != string::npos) {
      inputFileName = inputFileName.substr(inputFileName.find_last_of("/")+1);
    }
    if(inputFileName.find("mjd") != string::npos) {
      inputFileName = inputFileName.substr(inputFileName.find_last_of("mjd")+1);
    }
    TFile* outputFileReopen=TFile::Open(("enr"+inputFileName).c_str(), "UPDATE");
    outputFileReopen->cd();
        
    outputFileReopen->Write("",TObject::kOverwrite);
    outputFileReopen->Close();
  }

  // stop timer
  timer = clock() - timer;
  cout << "Processing took " << ((float)timer)/CLOCKS_PER_SEC/60.0 << " minutes." << endl;

  delete analysisDoc;
  return 0;
}

