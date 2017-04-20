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
  GATPostedDataTreeWriter mjdDataWriter("mjd"); //the name mjd is also used at the end

  // ****************************************************
  // Copy stuff from the input to the output and intercom
  // ****************************************************

  // Add event type and timestamp to the intercom
  GATInputDataPoster inputDataPoster;
  inputDataPoster.PostInputLeafTrivial("fEventType","eventType");
  inputDataPoster.PostInputLeafVector("fDigitizerData.fTimeStamp","timestamp");
  inputDataPoster.PostInputLeafVector("fDigitizerData.fID","channel");
  inputDataPoster.PostInputLeafVector("fDigitizerData.fIndex","index");
  selector.AddInput(&inputDataPoster);

  // add items to copy over from original tree
  //mjdDataWriter.AddInputLeafTrivial("fRunNumber", "run");
  //mjdDataWriter.AddInputLeafTrivial("fStartTime", "startTime");
  //mjdDataWriter.AddInputLeafTrivial("fStopTime", "stopTime");
  //mjdDataWriter.AddInputLeafTrivial("fStartClockTime", "startClockTime");
  //mjdDataWriter.AddInputLeafVector("fDigitizerData.fEnergy", "energy");
  //mjdDataWriter.AddInputLeafVector("fDigitizerData.fID", "channel");
  //mjdDataWriter.AddInputLeafVector("fDigitizerData.fIndex", "index");
  //mjdDataWriter.AddInputLeafVector("fDigitizerData.fTimeStamp", "timestamp");

  // Add GAT version info to output tree
  //TParameter<uint32_t> gatRevision("gatrev", strtol(GATUtils::GetGATRevision(), NULL, 16));
  //inputDataPoster.PostObject(&gatRevision);
  //mjdDataWriter.AddPostedTParameter("gatrev");


  // ********************************************************************
  // Add basic detector, multiplicity, and event time info to output tree
  // ********************************************************************
  /*
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

  */

  // *********************
  // Raw waveform analysis
  // *********************
  
  // Look for raw wf min and max before doing anything else, "rawWFMax", "rawWFMin"
  GATWFExtremumProcessor rawWFMaxProc(NULL, false, true, "rawWFMax");
  rawWFMaxProc.GetExtremumFinder().SetLocalMaximumTime(20160.*ns);
  selector.AddInput(&rawWFMaxProc);
  //mjdDataWriter.AddPostedVector(rawWFMaxProc.GetNameOfPostedVector());
  GATWFExtremumProcessor rawWFMinProc(NULL, false, false, "rawWFMin");
  rawWFMinProc.GetExtremumFinder().SetLocalMaximumTime(20160.*ns);
  selector.AddInput(&rawWFMinProc);
  //mjdDataWriter.AddPostedVector(rawWFMinProc.GetNameOfPostedVector());
  

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

  //if there is an auxiliary waveform(multiple sampling waveform)

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

      // blrwf: baseline-removed (raw) waveforms -- subtract average of first 200 samples
      baselineRemover = new MGWFBaselineRemover();
      baselineRemover->SetBaselineSamples(200);
      baselineRemoverProc = new GATWaveformTransformer(*baselineRemover, "choppedwf", "blrwf");
      selector.AddInput(baselineRemoverProc);

      // nlcwf: provide a version of waveforms that are corrected for ADC non-linearities
      nlcTransform = new MGWFNonLinearityCorrector();
      nlcTransform->SetTimeConstant_samples(140); // 1.4 us time constant for Radford time-lagged method
      nlcCorrector = new GATNonLinearityCorrector(*nlcTransform, nlcDBPath, "choppedwf", "nlcwf");
      selector.AddInput(nlcCorrector);
    }

  // nlcblrwf: baseline-removed ADC-NL-corrected waveforms -- subtract average of first 400 samples
  
  MGWFBaselineRemover baselineRemover400;
  baselineRemover400.SetBaselineSamples(400);
  GATWaveformTransformer nlcBaselineRemoverProc(baselineRemover400, "nlcwf", "nlcblrwf");
  selector.AddInput(&nlcBaselineRemoverProc);

  // d1wf and d2wf: 1st and 2nd derivatives of non-linearity-corrected waveforms
  MGWFBySampleDerivative bySampleDerivative;
  GATWaveformTransformer d1Proc(bySampleDerivative, "nlcwf", "d1wf");
  selector.AddInput(&d1Proc);
  GATWaveformTransformer d2Proc(bySampleDerivative, "d1wf", "d2wf");
  selector.AddInput(&d2Proc);


  // ********************************
  // Linear fits in different regions
  // ********************************

  // linear fits of baseline and flattop of raw waveforms, "RawWFblSlope", "RawWFblOffset", "RawWFblSlopeUnc", "RawWFblOffsetUnc", RawWFblChi2"
  /*
  GATWFLinearFitProcessor baselineFitterProc("bl");
  baselineFitterProc.GetFitter().SetFitSamples(200);
  selector.AddInput(&baselineFitterProc);
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedSlopesVector());
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedOffsetsVector());
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedSlopesUncVector());
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedOffsetsUncVector());
  mjdDataWriter.AddPostedVector(baselineFitterProc.GetNameOfPostedChi2sVector());


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
  /*
  // Linear fit of the whole waveform, "RawWFwholeWFLinFitSlope", "RawWFwholeWFLinFitOffset", "RawWFwholeWFLinFitSlopeUnc", "RawWFwholeWFLinFitOffsetUnc", RawWFwholeWFLinFitChi2"
  GATWFLinearFitProcessor wholeWFFitterProc("wholeWFLinFit");
  wholeWFFitterProc.GetFitter().SetFitSamples(2016);
  selector.AddInput(&wholeWFFitterProc);
  mjdDataWriter.AddPostedVector(wholeWFFitterProc.GetNameOfPostedSlopesVector());
  mjdDataWriter.AddPostedVector(wholeWFFitterProc.GetNameOfPostedSlopesUncVector());
  mjdDataWriter.AddPostedVector(wholeWFFitterProc.GetNameOfPostedOffsetsVector());
  mjdDataWriter.AddPostedVector(wholeWFFitterProc.GetNameOfPostedOffsetsUncVector());
  mjdDataWriter.AddPostedVector(wholeWFFitterProc.GetNameOfPostedChi2sVector());
  */

  // ***********
  // Time points
  // ***********
  
  // all time point calculators will use these, "blrwfFMRxx"
  vector<double> tpPercents; vector<string> tpLabels;
  //tpPercents.push_back(0.001); tpLabels.push_back("0p1");
  //tpPercents.push_back(0.01);  tpLabels.push_back("1");
  //tpPercents.push_back(0.03);  tpLabels.push_back("3");
  //tpPercents.push_back(0.1);   tpLabels.push_back("10");
  //tpPercents.push_back(0.2);   tpLabels.push_back("20");
  tpPercents.push_back(0.5);   tpLabels.push_back("50");
  //tpPercents.push_back(0.8);   tpLabels.push_back("80");
  //tpPercents.push_back(0.9);   tpLabels.push_back("90");
  //tpPercents.push_back(0.97);  tpLabels.push_back("97");
  //tpPercents.push_back(0.99);  tpLabels.push_back("99");


  // time points of baseline-subtracted, "blrwfFMRxx"
  GATWFTimePointProcessor* timePointProc = new GATWFTimePointProcessor("blrwf");
  for(size_t iPt = 0; iPt < tpPercents.size(); iPt++) {
    timePointProc->AddPoint(tpPercents[iPt], tpLabels[iPt].c_str());
  }
  selector.AddInput(timePointProc);
  //for(size_t iOut=0; iOut<timePointProc->GetNPostedVectors(); iOut++) {
    //string tpName = timePointProc->GetNameOfPostedVector(iOut);
    //if(tpName.find("FMR") != string::npos) mjdDataWriter.AddPostedVector(tpName.c_str());
  //}

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

  // start time finder on the trigger trapezoids, "triggerTrapt0"
  MGWFStartTimeFinder trapStartTimeFinder;
  trapStartTimeFinder.SetupBLR(false);
  // "good" waveforms should have a max within 4-14 us. Only search within that range
  trapStartTimeFinder.GetExtremumFinder().SetLocalMinimumTime(0.*us);
  trapStartTimeFinder.GetExtremumFinder().SetLocalMaximumTime(10.*us);
  trapStartTimeFinder.RestrictToRegion(0, 1000);
  trapStartTimeFinder.SetThreshold(2.); // 2 ADC ~ 1 keV on high gain channels. Later pull this from APDB.
  GATWFCalcParamProcessor trapStartTimeProc(trapStartTimeFinder, "triggerTrap");
  selector.AddInput(&trapStartTimeProc);
  mjdDataWriter.AddPostedVector(trapStartTimeProc.GetNameOfPostedVector());

  // best trap filters: original version with no NLC or charge trapping correction,
  // and use max instead of fixed-time pickoff, "trapE"
  MGWFTrapezoidalFilter optTrapFilter;
  optTrapFilter.SetDoNormalize();
  optTrapFilter.SetRampTime(4.*us);
  optTrapFilter.SetFlatTime(2.*us); // <- Will use 2.0*us once the new parameters uploaded to the database
  optTrapFilter.SetFitAndSubtractFlatBaseline();
  optTrapFilter.SetDecayConstant(72.*us);
  GATWaveformTransformer optTrapProc(optTrapFilter, "blrwf", "optTrapWFs");
  selector.AddInput(&optTrapProc);
  GATWFExtremumProcessor optTrapMaxProc("optTrapWFs", false, true, "trapE");
  selector.AddInput(&optTrapMaxProc);
  mjdDataWriter.AddPostedVector(optTrapMaxProc.GetNameOfPostedVector());

  // with NLC and use max instead of fixed-time pickoff, "trapENM", "trapENMSample", "trapENMMin"
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
  mjdDataWriter.AddPostedVector(optTrapMaxProcNLC.GetNameOfPostedVector(1)); // trapENM max sample
  GATWFExtremumProcessor optTrapMinProcNLC("optTrapWFsNLC", false, false, "trapENMMin");
  selector.AddInput(&optTrapMinProcNLC);
  mjdDataWriter.AddPostedVector(optTrapMinProcNLC.GetNameOfPostedVector());


  // DB-driven optimal trap filter, "trapENF"
  GATTrapezoidalFilter dbTrapFilter("nlcblrwf", "dbTrap", analysisDoc, "trapENF", GATTrapezoidalFilter::kEffQTrapping);
  dbTrapFilter.GetTransform().SetDoNormalize();
  dbTrapFilter.GetTransform().SetFitAndSubtractFlatBaseline();
  dbTrapFilter.GetTransform().SetDoNotResize(false);
  // currently using constant ramp and flat time for all channels.
  dbTrapFilter.SetManualRampTime(4.0*us);
  dbTrapFilter.SetManualFlatTime(2.5*us);
  selector.AddInput(&dbTrapFilter);
  GATWFFixedTimePickoffProcessor trapENFDBProc(-7.0*us+4.0*us+2.0*us, trapStartTimeProc.GetNameOfPostedVector(), "dbTrap", "trapENF");
  selector.AddInput(&trapENFDBProc);
  mjdDataWriter.AddPostedVector(trapENFDBProc.GetNameOfPostedVector());

  // the first sample of the short trap filter for the E=0 peak, "trapENFBL"
  GATTrapezoidalFilter dbTrapFilterBL("nlcblrwf", "dbTrapBL", analysisDoc, "trapENF", GATTrapezoidalFilter::kEffQTrapping);
  dbTrapFilterBL.GetTransform().SetDoNormalize();
  dbTrapFilterBL.GetTransform().SetFitAndSubtractFlatBaseline();
  dbTrapFilterBL.GetTransform().SetDoNotResize(false);
  dbTrapFilterBL.SetManualRampTime(4.0*us);
  dbTrapFilterBL.SetManualFlatTime(1.0*us);
  selector.AddInput(&dbTrapFilterBL);
  GATWFFixedTimePickoffProcessor trapENFDBProcBL(0.0*us, "", "dbTrapBL", "trapENFBL");
  selector.AddInput(&trapENFDBProcBL);
  mjdDataWriter.AddPostedVector(trapENFDBProcBL.GetNameOfPostedVector());

  //For the tail minimum cut, DON'T WANT TO DO PZ-CORRECT
  MGWFTrapezoidalFilter tailTrapFilter;
  tailTrapFilter.SetRampTime(4.*us);
  tailTrapFilter.SetFlatTime(1.5*us);
  GATWaveformTransformer tailTrapProc(tailTrapFilter, "blrwf", "tailTrapWFs");
  selector.AddInput(&tailTrapProc);
  GATWFExtremumProcessor optTrapTailMinProc("tailTrapWFs", false, false, "trapETailMin");
  optTrapTailMinProc.GetExtremumFinder().SetLocalMinimumTime(10000.*ns);
  selector.AddInput(&optTrapTailMinProc);
  mjdDataWriter.AddPostedVector(optTrapTailMinProc.GetNameOfPostedVector());
  /*
  // trap slope filters, "TSCurrentxxMax"
  vector<double> intTimes; vector<string> itLabels; vector<bool> tsitDoTPs;
  intTimes.push_back(50.*ns);  itLabels.push_back("50ns");  tsitDoTPs.push_back(false);
  intTimes.push_back(100.*ns); itLabels.push_back("100ns"); tsitDoTPs.push_back(false);
  intTimes.push_back(200.*ns); itLabels.push_back("200ns"); tsitDoTPs.push_back(false);
  size_t nTSFilters = intTimes.size();
  vector<MGWFTrapSlopeFilter*> tsFilters(nTSFilters);
  vector<GATWaveformTransformer*> tsCurrentProcs(nTSFilters);
  vector<GATWFExtremumProcessor*> tsCurrentMaxProcs(nTSFilters);
  vector<GATWFTimePointProcessor*> tsCurrentTimePointProcs(nTSFilters);
  for(size_t i=0; i<nTSFilters; i++) {
    tsFilters[i] = new MGWFTrapSlopeFilter;
    tsFilters[i]->SetPeakingTime(10.*ns);
    tsFilters[i]->SetIntegrationTime(intTimes[i]);
    tsFilters[i]->SetEvaluateMode(7);
    tsFilters[i]->OutputInternalParameter("s2");
    string tsWFs = "TSCurrent";
    tsWFs += itLabels[i];
    tsCurrentProcs[i] = new GATWaveformTransformer(*(tsFilters[i]), "blrwf", tsWFs.c_str());
    selector.AddInput(tsCurrentProcs[i]);
    tsCurrentMaxProcs[i] = new GATWFExtremumProcessor(tsWFs.c_str());
    //don't include final sample in maximum search
    tsCurrentMaxProcs[i]->GetExtremumFinder().SetLocalMaximumTime(20160.*ns);
    selector.AddInput(tsCurrentMaxProcs[i]);
    mjdDataWriter.AddPostedVector(tsCurrentMaxProcs[i]->GetNameOfPostedVector());

    // add timing when desired
    if(tsitDoTPs[i]) {
      tsCurrentTimePointProcs[i] = new GATWFTimePointProcessor(tsWFs.c_str());
      selector.AddInput(tsCurrentTimePointProcs[i]);
      for(size_t iPt = 0; iPt < tpPercents.size(); iPt++) {
        tsCurrentTimePointProcs[i]->AddPoint(tpPercents[iPt], tpLabels[iPt].c_str());
      }
      for(size_t iOut=0; iOut<tsCurrentTimePointProcs[i]->GetNPostedVectors(); iOut++) {
        mjdDataWriter.AddPostedVector(tsCurrentTimePointProcs[i]->GetNameOfPostedVector(iOut));
      }
    }
  }
  */
  /*
  //Time-point-based DCR slope "nlcblrwfSlope"
  MGWFSlopeCalculator tpSlopeCalculator;
  tpSlopeCalculator.SetAverageTime(1.*us);
  //tpSlopeCalculator.SetSecondStartSample(1900);
  GATWFDCRParamProcessor dcrTPSlopeProc(tpSlopeCalculator, "nlcblrwf"); // use adc non-lin corr wf's
  dcrTPSlopeProc.SetFirstPoint("blrwfFMR97", 200); // blrwfFMR97 should be "close enough" for nlcwf
  dcrTPSlopeProc.UseEndOfWaveform();
  selector.AddInput(&dcrTPSlopeProc);
  mjdDataWriter.AddPostedVector(dcrTPSlopeProc.GetNameOfPostedVector());

  //triangula filter, "triTrapMax" "triTrapIntegralW"
  MGWFTrapezoidalFilter triTrapFilter;
  triTrapFilter.SetDoNormalize();
  triTrapFilter.SetRampTime(100.*ns);
  triTrapFilter.SetFlatTime(10.*ns);
  triTrapFilter.SetFitAndSubtractFlatBaseline();
  //triTrapFilter.SetDecayConstant(72.*us);
  GATWaveformTransformer triTrapProc(triTrapFilter,"blrwf","triTrap");
  selector.AddInput(&triTrapProc);
  GATWFExtremumProcessor triTrapMaxProc("triTrap", false, true, "triTrapMax");
  //triTrapMaxProc.GetExtremumFinder().SetLocalMaximumTime(20160.*ns);
  selector.AddInput(&triTrapMaxProc);
  mjdDataWriter.AddPostedVector(triTrapMaxProc.GetNameOfPostedVector());

  MGPCIntegral triTrapFilterIntW;
  triTrapFilterIntW.SetRange(0,2000);
  GATWFCalcParamProcessor triTrapIntWProc(triTrapFilterIntW,"triTrap",false);
  selector.AddInput(&triTrapIntWProc);
  mjdDataWriter.AddPostedVector(triTrapIntWProc.GetNameOfPostedVector());
  */

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


  // *************************
  // Event level data cleaning
  // *************************
  /*
  // Pinghan's pulser-file-based event tagging
  GATDCPulserFileTagger pulserFileTagger("pinghan_pulser_tag");
  pulserFileTagger.SetFilePrefix("./pulser_");
  pulserFileTagger.SetFileSuffix(".root");
  pulserFileTagger.SetTreeName("pulsertree");
  pulserFileTagger.SetLeafName("Pulser");
  selector.AddInput(&pulserFileTagger);

  GATEventDataCleaningProcessor eventDCProc1("EventDC1");
  mjdDataWriter.AddPostedTParameter(eventDCProc1.GetNameOfPostedObject());
  eventDCProc1.SetBitParameter(0, "eventType");
  eventDCProc1.SetBitHighCut(0,0.5);
  eventDCProc1.SetBitParameter(1, "pinghan_pulser_tag");
  eventDCProc1.SetBitHighCut(1,0.5);
  selector.AddInput(&eventDCProc1);

  // Clint and Brian's automatic threshold finder
  GATThresholdReader threshFinder("threshKeV","threshSigma",runNumber);
  string threshName = Form("./thresh_%i.root",runNumber);
  ifstream threshFile(threshName);
  if (threshFile.good()) {
    selector.AddInput(&threshFinder);
    mjdDataWriter.AddPostedVector(threshFinder.GetNameOfThresholdVector());
    mjdDataWriter.AddPostedVector(threshFinder.GetNameOfSigmaVector());
  }
  else cout << "GATThresholdReader : No file found.  Proceeding without ...\n";
  */
  // *************************
  // Noise Power Spectrum
  // *************************

  // Do DFT of second derivative waveforms

  MGPCNoiseBands nb("MGPCNoiseBands");
  vector<pair<double,double> > bandBoundaries;
  bandBoundaries.push_back(make_pair(0*MHz,2*MHz));
  bandBoundaries.push_back(make_pair(2*MHz,10*MHz));
  bandBoundaries.push_back(make_pair(10*MHz,20*MHz));
  bandBoundaries.push_back(make_pair(5*MHz,30*MHz));
  bandBoundaries.push_back(make_pair(30*MHz,35*MHz));
  bandBoundaries.push_back(make_pair(48*MHz,50*MHz));
  bandBoundaries.push_back(make_pair(0*MHz,50*MHz));
  nb.SetBandBoundaries(bandBoundaries);

  GATWFTransientNoiseTagger noiseBandsProc(nb, analysisDoc, 3, ChanMap, "d2wf");
  selector.AddInput(&noiseBandsProc);
  mjdDataWriter.AddPostedVector(noiseBandsProc.GetNameOfNormalizedValuesVector());
  for(size_t iParam=0; iParam<noiseBandsProc.GetNPostedVectors(); iParam++) {
    mjdDataWriter.AddPostedVector(noiseBandsProc.GetNameOfPostedVector(iParam));
  }

  // *************************
  // Waveform level data cleaning
  // *************************

  MGPCSingleSampleFilterSigned wfSSFSi;
  wfSSFSi.SetRange(0,2016);
  wfSSFSi.SetNeighborN(50);
  GATWFCalcParamProcessor wfSSFSiproc(wfSSFSi,"blrwf");
  selector.AddInput(&wfSSFSiproc);
  mjdDataWriter.AddPostedVector(wfSSFSiproc.GetNameOfPostedVector());

  MGPCNorm2ndDiffPeakValue wfNDPV;
  GATWFCalcParamProcessor NdiffblrWproc(wfNDPV,"blrwf");
  selector.AddInput(&NdiffblrWproc);
  mjdDataWriter.AddPostedVector(NdiffblrWproc.GetNameOfPostedVector());

  GATDCBoxCut dc_SSPosSpikeBL("SSPosSpikeBL"); // wfDCBits 0
  dc_SSPosSpikeBL.SetInIsBad(true);
  dc_SSPosSpikeBL.AddDimension("blrwfSSFSigned",1,230,0,0);
  dc_SSPosSpikeBL.AddDimension("blrwfNorm2ndDiffPeakValue",1,1e-5,0,0);
  selector.AddInput(&dc_SSPosSpikeBL);

  GATDCBoxCut dc_SSNegSpikeBL("SSNegSpikeBL"); // wfDCBits 1
  dc_SSNegSpikeBL.SetInIsBad(true);
  dc_SSNegSpikeBL.AddDimension("blrwfSSFSigned",0,0,1,-200);
  dc_SSNegSpikeBL.AddDimension("blrwfNorm2ndDiffPeakValue",1,1e-5,0,0);
  selector.AddInput(&dc_SSNegSpikeBL);

  GATDCBoxCut dc_SSPosSpikePhys("SSPosSpikePhys"); // wfDCBits 2
  dc_SSPosSpikePhys.SetInIsBad(true);
  dc_SSPosSpikePhys.AddDimension("blrwfSSFSigned",1,190,0,0);
  dc_SSPosSpikePhys.AddDimension("blrwfNorm2ndDiffPeakValue",1,1.07e-6,1,1e-5);
  dc_SSPosSpikePhys.AddDimension("channel",1,593.5,1,594.5);
  selector.AddInput(&dc_SSPosSpikePhys);

  GATDCBoxCut dc_SSNegSpikePhys("SSNegSpikePhys"); // wfDCBits 3
  dc_SSNegSpikePhys.SetInIsBad(true);
  dc_SSNegSpikePhys.AddDimension("blrwfSSFSigned",0,0,1,-190);
  dc_SSNegSpikePhys.AddDimension("blrwfNorm2ndDiffPeakValue",1,1.07e-6,1,1e-5);
  dc_SSNegSpikePhys.AddDimension("channel",1,593.5,1,594.5);
  selector.AddInput(&dc_SSNegSpikePhys);

  // Late and Early Triggered WFs //////////////////////////////

  GATDCBoxCut dc_EarlyTrigger("EarlyTrigger"); //wfDCBits 4
  dc_EarlyTrigger.SetInIsBad(true);
  dc_EarlyTrigger.AddDimension("blrwfFMR50",0,0.,1,8000);
  dc_EarlyTrigger.AddDimension("trapENFCal",1,30,0,0.);
  selector.AddInput(&dc_EarlyTrigger);

  GATDCBoxCut dc_LateTrigger("LateTrigger"); //wfDCBits 5
  dc_LateTrigger.SetInIsBad(true);
  dc_LateTrigger.AddDimension("blrwfFMR50",1,12000.,0,0);
  dc_LateTrigger.AddDimension("trapENFCal",1,30,0,0.);
  selector.AddInput(&dc_LateTrigger);

  // Saturated WFs /////////////////////////////////////////////
  GATDCBoxCut dc_PosSatWFs("PosSaturatedWFs"); // wfDCBits 6
  dc_PosSatWFs.SetInIsBad(true);
  dc_PosSatWFs.AddDimension("rawWFMax",1,8190,1,8192);
  selector.AddInput(&dc_PosSatWFs);

  GATDCBoxCut dc_NegSatWFs("NegSaturatedWFs"); // wfDCBits 7
  dc_NegSatWFs.SetInIsBad(true);
  dc_NegSatWFs.AddDimension("rawWFMin",1,-8192,1,-8190);
  selector.AddInput(&dc_NegSatWFs);

  /*
  //Fast trap for transition layer multi-site detection
  MGWFTrapezoidalFilter fastTrap;
  fastTrap.SetDoNormalize();
  fastTrap.SetRampTime(.1*us); // FIXME: pull optimum from DB
  fastTrap.SetFlatTime(.3*us);
  fastTrap.SetFitAndSubtractFlatBaseline(); // FIXME: use resting baseline from DB
  fastTrap.SetDecayConstant(72.*us); // FIXME: use PZ const from DB
  GATWaveformTransformer fastTrapNLCProc(fastTrap, "nlcblrwf", "fastTrapNLCWFs");
  selector.AddInput(&fastTrapNLCProc);
  
  //For transition layer multi-site rejection, count the threshold crossings of a .1us ramp-time trap filter
  MGPCCountThresholdCrossings TLMultiSiteFinder;
  TLMultiSiteFinder.SetThreshold(10);
  TLMultiSiteFinder.SetRequiredDifference(.001);
  GATWFCountCrossingsParamProcessor TLMultiSiteProc(TLMultiSiteFinder, "fastTrapNLCWFs");
  selector.AddInput(&TLMultiSiteProc);
  mjdDataWriter.AddPostedVector(TLMultiSiteProc.GetNameOfPostedVector());

  //For pile-up rejection, count the threshold crossings of a .3us ramp-time trap filter
  MGWFTrapezoidalFilter pileUpTrap;
  pileUpTrap.SetDoNormalize();
  pileUpTrap.SetRampTime(.3*us); // FIXME: pull optimum from DB
  pileUpTrap.SetFlatTime(.3*us);
  pileUpTrap.SetFitAndSubtractFlatBaseline(); // FIXME: use resting baseline from DB
  pileUpTrap.SetDecayConstant(72.*us); // FIXME: use PZ const from DB
  GATWaveformTransformer pileUpTrapNLCProc(pileUpTrap, "nlcblrwf", "pileUpTrapNLCWFs");
  selector.AddInput(&pileUpTrapNLCProc);

  MGPCCountThresholdCrossings pileUpFinder;
  pileUpFinder.SetThreshold(20);
  pileUpFinder.SetRequiredDifference(.01);
  GATWFCountCrossingsParamProcessor pileUpProc(pileUpFinder, "pileUpTrapNLCWFs");
  selector.AddInput(&pileUpProc);
  //mjdDataWriter.AddPostedVector(pileUpProc.GetNameOfPostedVector());
  //mjdDataWriter.AddPostedVector(pileUpProc.GetNameOfPostedRisingTimeVector());

  GATDCBoxCut dc_PileUpWFs("PileUpWFs"); // wfDCBits 8
  dc_PileUpWFs.SetInIsBad(false);
  dc_PileUpWFs.AddDimension("pileUpTrapNLCWFsnRisingX",0, 0, 1, 2);//number of rising threshold crossings of pile-up trap filter
  selector.AddInput(&dc_PileUpWFs);
  //mjdDataWriter.AddPostedVector(dc_PileUpWFs.GetNameOfPostedVector());
*/
  /*
  // Set the bits
  GATWFDataCleaningProcessor wfdcProc("wfDC");
  mjdDataWriter.AddPostedVector(wfdcProc.GetNameOfPostedObject());
  wfdcProc.SetBitParameter(0,"SSPosSpikeBL");
  wfdcProc.SetBitHighCut(0,0.5);
  wfdcProc.SetBitParameter(1,"SSNegSpikeBL");
  wfdcProc.SetBitHighCut(1,0.5);
  wfdcProc.SetBitParameter(2,"SSPosSpikePhys");
  wfdcProc.SetBitHighCut(2,0.5);
  wfdcProc.SetBitParameter(3,"SSNegSpikePhys");
  wfdcProc.SetBitHighCut(3,0.5);

  wfdcProc.SetBitParameter(4,"EarlyTrigger");
  wfdcProc.SetBitHighCut(4,0.5);
  wfdcProc.SetBitParameter(5,"LateTrigger");
  wfdcProc.SetBitHighCut(5,0.5);

  wfdcProc.SetBitParameter(6,"PosSaturatedWFs");
  wfdcProc.SetBitHighCut(6,0.5);
  wfdcProc.SetBitParameter(7,"NegSaturatedWFs");
  wfdcProc.SetBitHighCut(7,0.5);
*/
  //wfdcProc.SetBitParameter(8, "PileUpWFs"); //Pile-up detection: trap filter rising threshold crossings counter
  //wfdcProc.SetBitHighCut(8,0.5);

  /*
  // Now that the bits are set we can add the data cleaning parameter values to the output tree if we wish
  for(size_t i=0; i<wfdcProc.GetNBits(); i++) {
    string name = wfdcProc.GetNameOfBit(i);
    if(name != "") mjdDataWriter.AddPostedVector(name.c_str());
  }
  //The Data cleaning processor must be added to the selector at the end!
  selector.AddInput(&wfdcProc);
*/
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
    TFile* outputFileReopen=TFile::Open(("enr"+inputFileName).c_str(), "UPDATE");
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

