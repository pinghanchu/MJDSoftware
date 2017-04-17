// 2016.11.17 Pinghan Chu
// 1.A simple version of GAT. Only save rawWFMax, rawWFMin, trapE and trapENFDB.
// 2.Can process multiple sampling waveforms.
//
// 2017.1.10 Pinghan Chu
// 1. The energy parameters are trapE(4-1.5-4) (will change to 4-2-4 after reprocessing), trapENM(4-2-4) and trapENF(4-2.5-4)
// 2. Add the E=0 parameters, trapENFBL, using the first sample of a trap(4-1-4)
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
#include "GATEnergyCalibrationAscii.hh"
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

using namespace std;
using namespace CLHEP;
using namespace MJDB;

int main(int argc, char** argv)
{
  // Give usage info for no input (user must supply an input file)
  if(argc < 2 || argc > 3) {
    cout << "Usage: " << argv[0] << " [file] (calibrationMapFile)" << endl;
    return 1;
  }

  // start timer
  clock_t timer;
  timer = clock();

  string calibrationMapFile = "";
  if(argc == 3) calibrationMapFile = argv[2];


  // *******************************
  // Pull stuff out from input files
  // *******************************

  // open input file
  TFile* file = TFile::Open(argv[1]);
  if(file == NULL) return 0;

  // find MGTEvent tree in input file 
  TTree* MGTree = (TTree*) file->Get("MGTree");
  if(MGTree == NULL) {
    cout << "MGTree not found" << endl;
    return 0;
  }
  // find MGTEvent tree in input file
  TTree* MCATree = (TTree*) file->Get("MCATree");
  if(MCATree == NULL) {
    cout << "MCATree not found. Continue" << endl;
  }
  // find MJTChannelMap in input file
  MJTChannelMap* ChanMap = (MJTChannelMap*) file->Get("ChannelMap");
  if(ChanMap == NULL) {
    cout << "ChannelMap not found. Continue" << endl;
  }
  // find MJTChannelSettings in input file
  MJTChannelSettings* ChanSettings = (MJTChannelSettings*) file->Get("ChannelSettings");
  if(ChanSettings == NULL) {
    cout << "ChannelSettings not found. Continue" << endl;
  }

  // find MJTRun in input file
  bool bUseMS = false;
  MJTRun* RunSettings = (MJTRun*) file->Get("run");
  if(RunSettings == NULL) {
    cout << "MJTRun not found. Continue" << endl;
  }
  else bUseMS = RunSettings->GetUseMultisampling();


  // **********************************
  // Set up data processing environment
  // **********************************

  // Set up access to Analysis Parameters Database
  // Processors can use this pointer for common access.
  // No constructor args means we are communicating with the default APDB
  MJAnalysisDoc* analysisDoc = new MJAnalysisDoc();
  //  analysisDoc->SetAlternateDB("", "mjdbatest","","","");

  // prepare selector and writer so we can add to them as we go
  TAMSelector selector;
  GATPostedDataTreeWriter mjdDataWriter("mjd"); //the name mjd is also used at the end

  // ****************************************************
  // Copy stuff from the input to the output and intercom
  // ****************************************************

  // Add event type and timestamp to the intercom
  GATInputDataPoster inputDataPoster;
  inputDataPoster.PostInputLeafTrivial("fEventType","eventType");
  inputDataPoster.PostInputLeafVector("fDigitizerData.fTimeStamp","timestamp");
  inputDataPoster.PostInputLeafVector("fDigitizerData.fID","channel");
  selector.AddInput(&inputDataPoster);

  // add items to copy over from original tree
  mjdDataWriter.AddInputLeafTrivial("fRunNumber", "run");
  mjdDataWriter.AddInputLeafTrivial("fStartTime", "startTime");
  mjdDataWriter.AddInputLeafTrivial("fStopTime", "stopTime");
  mjdDataWriter.AddInputLeafVector("fDigitizerData.fEnergy", "energy");
  mjdDataWriter.AddInputLeafVector("fDigitizerData.fID", "channel");
  mjdDataWriter.AddInputLeafVector("fDigitizerData.fTimeStamp", "timestamp");

  // Add GAT version info to output tree
  TParameter<uint32_t> gatRevision("gatrev", strtol(GATUtils::GetGATRevision(), NULL, 16));
  inputDataPoster.PostObject(&gatRevision);
  mjdDataWriter.AddPostedTParameter("gatrev");


  // ********************************************************************
  // Add basic detector, multiplicity, and event time info to output tree
  // ********************************************************************
 
  // event time info processor: timeMT, dateMT
  GATEventTimeInfoProcessor eventTimeInfoProcessor;
  selector.AddInput(&eventTimeInfoProcessor);
  for(size_t i=0; i<eventTimeInfoProcessor.GetNPostedVectors(); i++) {
    mjdDataWriter.AddPostedVector(eventTimeInfoProcessor.GetNameOfPostedVector(i));
  }
 
  // detector info processor: detector ID, position / row, enr/nat, etc
  GATDetInfoProcessor detInfoProcessor;
  selector.AddInput(&detInfoProcessor);
  for(size_t i=0; i<detInfoProcessor.GetNPostedVectors(); i++) {
    mjdDataWriter.AddPostedVector(detInfoProcessor.GetNameOfPostedVector(i));
  }

  // (m,i,j)
  GATMultiplicityProcessor multProc;
  selector.AddInput(&multProc);
  mjdDataWriter.AddPostedObject(multProc.GetNameOfPostedObject());


  // *********************
  // Raw waveform analysis
  // *********************
 
  // Look for raw wf min and max before doing anything else
  GATWFExtremumProcessor rawWFMaxProc(NULL, false, true, "rawWFMax");
  rawWFMaxProc.GetExtremumFinder().SetLocalMaximumTime(20160.*ns);
  selector.AddInput(&rawWFMaxProc);
  mjdDataWriter.AddPostedVector(rawWFMaxProc.GetNameOfPostedVector());
  GATWFExtremumProcessor rawWFMinProc(NULL, false, false, "rawWFMin");
  rawWFMinProc.GetExtremumFinder().SetLocalMaximumTime(20160.*ns);
  selector.AddInput(&rawWFMinProc);
  mjdDataWriter.AddPostedVector(rawWFMinProc.GetNameOfPostedVector());


  // **************************************************
  // Prepare modified waveforms for subsequent analyses
  // **************************************************
  //
  // These just get posted to the intercom, they don't get written to the output file.

  GATMSWFPoster *wfPoster = NULL;
  MGWFWindower *wfWindower = NULL;
  GATWaveformTransformer *chopperProc = NULL;
  MGWFBaselineRemover *baselineRemover = NULL;
  GATWaveformTransformer *baselineRemoverProc = NULL;
  MGWFNonLinearityCorrector *nlcTransform = NULL;
  GATNonLinearityCorrector *nlcCorrector = NULL;
  const char* nlcDBPath = "/project/projectdirs/majorana/data/production/NLCDB";

  if(bUseMS)
    {
      // multisamplewf: unpacked multi-sampled waveform
      wfPoster = new GATMSWFPoster("multisamplewf");
      selector.AddInput(wfPoster);

      // blrwf: baseline-removed (raw) waveforms -- subtract average of first 200 samples    
      baselineRemover = new MGWFBaselineRemover();
      baselineRemover->SetBaselineSamples(200);
      baselineRemoverProc = new GATWaveformTransformer(*baselineRemover, "multisamplewf", "blrwf");
      selector.AddInput(baselineRemoverProc);

      // nlcwf: provide a version of waveforms that are corrected for ADC non-linearities
      nlcTransform = new MGWFNonLinearityCorrector();
      nlcTransform->SetTimeConstant_samples(140); // 1.4 us time constant for Radford time-lagged method
      nlcCorrector = new GATNonLinearityCorrector(*nlcTransform, nlcDBPath, "multisamplewf", "nlcwf");
      selector.AddInput(nlcCorrector);

    }
  else
    {
      // choppedwf: chop off those last few samples of waveform
      wfWindower = new MGWFWindower(0, 2016);
      chopperProc = new GATWaveformTransformer(*wfWindower, NULL, "choppedwf");
      selector.AddInput(chopperProc);  
      //mjdDataWriter.AddPostedObject("choppedwf");

      //mjdDataWriter.AddPostedObject("choppedwf");
      // blrwf: baseline-removed (raw) waveforms -- subtract average of first 200 samples  
      baselineRemover = new MGWFBaselineRemover();
      baselineRemover->SetBaselineSamples(200);
      baselineRemoverProc = new GATWaveformTransformer(*baselineRemover, "choppedwf", "blrwf");
      selector.AddInput(baselineRemoverProc);  
      //mjdDataWriter.AddPostedObject("blrwf");

      // nlcwf: provide a version of waveforms that are corrected for ADC non-linearities
      nlcTransform = new MGWFNonLinearityCorrector();
      nlcTransform->SetTimeConstant_samples(140); // 1.4 us time constant for Radford time-lagged method
      nlcCorrector = new GATNonLinearityCorrector(*nlcTransform, nlcDBPath, "choppedwf", "nlcwf");
      selector.AddInput(nlcCorrector);
      //mjdDataWriter.AddPostedObject("nlcwf");

    }


  // nlcblrwf: baseline-removed ADC-NL-corrected waveforms -- subtract average of first 400 samples
  MGWFBaselineRemover baselineRemover400;
  baselineRemover400.SetBaselineSamples(400);
  GATWaveformTransformer nlcBaselineRemoverProc(baselineRemover400, "nlcwf", "nlcblrwf");
  selector.AddInput(&nlcBaselineRemoverProc);
  //mjdDataWriter.AddPostedObject("nlcblrwf");

  GATWFLinearFitProcessor baselineFitterProc("bl");
  baselineFitterProc.GetFitter().SetFitSamples(200);
  selector.AddInput(&baselineFitterProc);
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedSlopesVector());
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedOffsetsVector());
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedSlopesUncVector());
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedOffsetsUncVector());
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedChi2sVector());
  /*
  //"RawWFftSlope", "RawWFftOffset", "RawWFftSlopeUnc", "RawWFftOffsetUnc", RawWFftChi2"
  GATWFLinearFitProcessor flattopFitterProc("ft");
  flattopFitterProc.GetFitter().SetFitSamples(200, 1400);
  selector.AddInput(&flattopFitterProc);
  mjdDataWriter.AddPostedVector(flattopFitterProc.GetNameOfPostedSlopesVector());
  mjdDataWriter.AddPostedVector(flattopFitterProc.GetNameOfPostedSlopesUncVector());
  mjdDataWriter.AddPostedVector(flattopFitterProc.GetNameOfPostedOffsetsVector());
  mjdDataWriter.AddPostedVector(flattopFitterProc.GetNameOfPostedOffsetsUncVector());
  mjdDataWriter.AddPostedVector(flattopFitterProc.GetNameOfPostedChi2sVector());
  */

  // d1wf and d2wf: 1st and 2nd derivatives of non-linearity-corrected waveforms
  /*
  MGWFBySampleDerivative bySampleDerivative;
  GATWaveformTransformer d1Proc(bySampleDerivative, "nlcwf", "d1wf");
  selector.AddInput(&d1Proc);
  GATWaveformTransformer d2Proc(bySampleDerivative, "d1wf", "d2wf");
  selector.AddInput(&d2Proc);
  */
  
  // ********************
  // Trapezoidal filters
  // ********************

  // shorter trapezoid for offline re-triggering
  GATTrapezoidalFilter triggerTrapFilter("nlcblrwf", "triggerTrap");
  triggerTrapFilter.GetTransform().SetDoNormalize();
  triggerTrapFilter.GetTransform().SetFitAndSubtractFlatBaseline();
  triggerTrapFilter.SetManualRampTime(1.0*us);
  triggerTrapFilter.SetManualFlatTime(1.5*us);
  triggerTrapFilter.SetManualPZCorrection(72.*us); // correct for RC decay only
  triggerTrapFilter.GetTransform().SetDoNotResize(false);
  selector.AddInput(&triggerTrapFilter);
  //mjdDataWriter.AddPostedObject("triggerTrap");  

  // start time finder on the trigger trapezoids
  MGWFStartTimeFinder trapStartTimeFinder;
  trapStartTimeFinder.SetupBLR(false);
  // "good" waveforms should have a max within 4-14 us. Only search within that range
  trapStartTimeFinder.GetExtremumFinder().SetLocalMinimumTime(0.*us);
  trapStartTimeFinder.GetExtremumFinder().SetLocalMaximumTime(10.*us);
  trapStartTimeFinder.RestrictToRegion(0, 1000);
  trapStartTimeFinder.SetThreshold(2.); // 2 ADC ~ 1 keV on high gain channels. Later pull this from APDB.
  GATWFCalcParamProcessor trapStartTimeProc(trapStartTimeFinder, "triggerTrap"); 
  selector.AddInput(&trapStartTimeProc);
  //mjdDataWriter.AddPostedVector(trapStartTimeProc.GetNameOfPostedVector());
 
  // best trap filters: original version with no NLC or charge trapping correction, 
  // and use max instead of fixed-time pickoff, trapE
  MGWFTrapezoidalFilter optTrapFilter;
  optTrapFilter.SetDoNormalize();
  optTrapFilter.SetRampTime(4.*us);
  optTrapFilter.SetFlatTime(2.*us);
  optTrapFilter.SetFitAndSubtractFlatBaseline();
  optTrapFilter.SetDecayConstant(72.*us);
  GATWaveformTransformer optTrapProc(optTrapFilter, "blrwf", "optTrapWFs");
  selector.AddInput(&optTrapProc);
  //mjdDataWriter.AddPostedObject("optTrapWFs");

  GATWFExtremumProcessor optTrapMaxProc("optTrapWFs", false, true, "trapE");
  selector.AddInput(&optTrapMaxProc);
  mjdDataWriter.AddPostedVector(optTrapMaxProc.GetNameOfPostedVector());

  // with NLC and use max instead of fixed-time pickoff, "trapENM"
  MGWFTrapezoidalFilter optTrapFilterNLC;
  optTrapFilterNLC.SetDoNormalize();
  optTrapFilterNLC.SetRampTime(4.*us);
  optTrapFilterNLC.SetFlatTime(2.*us);
  optTrapFilterNLC.SetFitAndSubtractFlatBaseline();
  optTrapFilterNLC.SetDecayConstant(72.*us);
  GATWaveformTransformer optTrapProcNLC(optTrapFilterNLC, "nlcblrwf", "optTrapWFsNLC");
  selector.AddInput(&optTrapProcNLC);
  GATWFExtremumProcessor optTrapMaxProcNLC("optTrapWFsNLC", false, true, "trapENM");
  selector.AddInput(&optTrapMaxProcNLC);
  mjdDataWriter.AddPostedVector(optTrapMaxProcNLC.GetNameOfPostedVector());

  // DB-driven optimal trap filter 
  // Optimal trap filter, trapENF
  GATTrapezoidalFilter dbTrapFilter("nlcblrwf", "dbTrap", analysisDoc, "trapENF", GATTrapezoidalFilter::kEffQTrapping);
  dbTrapFilter.GetTransform().SetDoNormalize();
  dbTrapFilter.GetTransform().SetFitAndSubtractFlatBaseline();
  dbTrapFilter.GetTransform().SetDoNotResize(false);
  dbTrapFilter.SetManualRampTime(4.0*us);
  dbTrapFilter.SetManualFlatTime(2.5*us);
  selector.AddInput(&dbTrapFilter);
  //mjdDataWriter.AddPostedObject();
  GATWFFixedTimePickoffProcessor trapENFDBProc(-7.0*us+4.0*us+2.0*us, trapStartTimeProc.GetNameOfPostedVector(), "dbTrap", "trapENF");
  selector.AddInput(&trapENFDBProc);
  mjdDataWriter.AddPostedVector(trapENFDBProc.GetNameOfPostedVector());

  // E=0 peak, trapENFBL
  GATTrapezoidalFilter dbTrapFilterBL("nlcblrwf", "dbTrapBL", analysisDoc, "trapENF", GATTrapezoidalFilter::kEffQTrapping);
  dbTrapFilterBL.GetTransform().SetDoNormalize();
  dbTrapFilterBL.GetTransform().SetFitAndSubtractFlatBaseline();
  dbTrapFilterBL.GetTransform().SetDoNotResize(false);
  dbTrapFilterBL.SetManualRampTime(4.0*us);
  dbTrapFilterBL.SetManualFlatTime(1.0*us);
  selector.AddInput(&dbTrapFilterBL);
  GATWFFixedTimePickoffProcessor trapENFDBProcBL(0.0*us,"", "dbTrapBL", "trapENFBL");
  selector.AddInput(&trapENFDBProcBL);
  mjdDataWriter.AddPostedVector(trapENFDBProcBL.GetNameOfPostedVector());
  
  // ******************
  // Energy calibration
  // ******************

  // Energy calibration, "trapECal", "trapEMinCal", "trapENFCal"
  //GATCalibrationMapAscii calibrationMap;

  cout << "Fill trapE" << endl;
  GATEnergyCalibration trapECalibration(analysisDoc, "trapE");
  trapECalibration.UseCalibrationFor("trapE");
  selector.AddInput(&trapECalibration);
  mjdDataWriter.AddPostedVector(trapECalibration.GetNameOfPostedVector());

  cout << "Fill trapENM" << endl;
  GATEnergyCalibration trapENMCalibration(analysisDoc, "trapENM");
  trapENMCalibration.UseCalibrationFor("trapENM");
  selector.AddInput(&trapENMCalibration);
  mjdDataWriter.AddPostedVector(trapENMCalibration.GetNameOfPostedVector());

  cout << "Fill trapENF" << endl;
  GATEnergyCalibration trapENFCalibration(analysisDoc, "trapENF");
  trapENFCalibration.UseCalibrationFor("trapENF");
  selector.AddInput(&trapENFCalibration);
  mjdDataWriter.AddPostedVector(trapENFCalibration.GetNameOfPostedVector());
  /*
  GATEnergyCalibration trapENFBLCalibration(analysisDoc, "trapENFBL");
  trapENFBLCalibration.UseCalibrationFor("trapENF");
  selector.AddInput(&trapENFBLCalibration);
  mjdDataWriter.AddPostedVector(trapENFBLCalibration.GetNameOfPostedVector());
  */
  // cout << trapENFBLCalibration.GetNameOfPostedVector() << endl;

  // ************
  // Now process!
  // ************

  selector.AddInput(&mjdDataWriter);
  MGTree->Process(&selector);


  // *************************************
  // Finally, clean up and write to output
  // *************************************

  if(MCATree != NULL || ChanMap != NULL || ChanSettings != NULL){
    string inputFileName = file->GetName();
    if(inputFileName.find("/") != string::npos) {
      inputFileName = inputFileName.substr(inputFileName.find_last_of("/")+1);
    }
    if(inputFileName.find("OR") != string::npos) {
      inputFileName = inputFileName.substr(inputFileName.find_last_of("OR")+1);
    }
    TFile* outputFileReopen=TFile::Open(("mjd"+inputFileName).c_str(), "UPDATE");
    outputFileReopen->cd();
    if(MCATree != NULL) {
      cout<<"Add MCA Tree"<<endl;
      TTree *MCATree2=(TTree*)MCATree->CopyTree("1");
      MCATree2->Write();
    }
    if(ChanMap != NULL) {
      cout<<"Add ChannelMap"<<endl;
      MJTChannelMap ChanMap2(*ChanMap);
      ChanMap2.Write();
    }
    if(ChanSettings != NULL) {
      cout<<"Add ChannelSettings"<<endl;
      MJTChannelSettings ChanSettings2(*ChanSettings);
      ChanSettings2.Write();
    }
    outputFileReopen->Write("",TObject::kOverwrite);
    outputFileReopen->Close();
  }
  
  // stop timer
  timer = clock() - timer;
  cout << "Processing took " << ((float)timer)/CLOCKS_PER_SEC/60.0 << " minutes." << endl;

  delete analysisDoc;
  return 0;
}

