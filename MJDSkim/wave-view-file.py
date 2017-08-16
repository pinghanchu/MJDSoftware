#!/usr/common/usg/software/python/2.7.9/bin/python
import sys
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import numpy as np

def main(argv):
    """Interactive-draw or rapid-draw waveforms that pass a given TCut.

       Bug: Doesn't always work with a TChain.  Add input files together
       with hadd and use a single TFile.
    """

    scanSpeed = 0.5
    opt1, opt2 = "", ""
    intMode, printWF, warpMode = False, False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    import matplotlib.pyplot as plt

    # Set input file and cuts
    inputFile = TFile("/projecta/projectdirs/majorana/users/bxyzhu/waveskim2/waveSkimDS0_55.root")
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    theCut = inputFile.Get("theCut").GetTitle()
    # Add additional cuts
    theCut += "&& trapENFCal > 10"

    # Print cut and events passing cut
    print "Using cut:\n",theCut,"\n"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."


    # Make a figure (only setting data in the loop is faster)
    fig = plt.figure(figsize=(10,7), facecolor='w')
    a1 = plt.subplot(111)
    a1.set_xlabel("time (ns)")
    a1.set_ylabel("ADC")
    p1, = a1.plot(np.ones(1), np.ones(1), color='blue', alpha=0.8, label='data')
    plt.show(block=False)

    # Loop over events
    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList !=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        nWFs = waveTree.MGTWaveforms.size()

        if (nWFs==0):
            print "Error - nWFs:",nWFs,"nChans",nChans
            continue

        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))


        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            energy = waveTree.trapENFCal.at(iH)
            wf = waveTree.MGTWaveforms.at(iH)
            enm = waveTree.trapENM.at(iH)
            signal = processWaveform(wf)
            waveBLSub = signal.GetWaveBLSub()
            waveTS = signal.GetTS()
            print "%d / %d  Run %d  nCh %d  chan %d  trapENFCal %.1f trapENM %.1f" % (iList,nList,run,nChans,chan,energy, enm)

            # fill the figure
            p1.set_ydata(waveBLSub)
            p1.set_xdata(waveTS)

            xmin, xmax = np.amin(waveTS), np.amax(waveTS)
            ymin, ymax = np.amin(waveBLSub), np.amax(waveBLSub)
            a1.set_xlim([xmin,xmax])
            a1.set_ylim([ymin-abs(0.1*ymin),ymax+abs(0.1*ymax)])
            plt.title("Run %d  Channel %d  Energy %.2f keV" % (run,chan,energy))

            plt.tight_layout()
            plt.pause(scanSpeed)

def baselineParameters(signalRaw):
    """ Finds basic parameters of baselines using first 500 samples. """
    rms = 0
    slope = 0
    baselineMean = 0
    baselineAveSq = 0
    for x in xrange(0,500):
        wfX = signalRaw[x]
        baselineMean += wfX
        baselineAveSq += wfX*wfX
    baselineMean /= 500.
    baselineAveSq /= 500.
    rms = np.sqrt( baselineAveSq - baselineMean*baselineMean);
    return rms, slope, baselineMean

class processWaveform:
    """ Handy class for auto-processing waveforms into various numpy arrays. """
    def __init__(self, wave, remLo=0, remHi=2):

        self.waveMGT = wave                            # input an MGTWaveform object
        self.offset = self.waveMGT.GetTOffset()        # time offset [ns]
        vec = wave.GetVectorData()
        npArr = np.fromiter(vec, dtype=np.double)      # raw numpy array

        hist = wave.GimmeUniqueHist()                  # get timestamp limits and make an array
        self.start = hist.GetXaxis().GetXmin() + 5     # add 5 ns to make it start at 0 (default is -5.)
        self.stop = hist.GetXaxis().GetXmax() + 5.
        self.binsPerNS = (self.stop - self.start) / hist.GetNbinsX()
        ts = np.arange(self.start,self.stop,self.binsPerNS)

        removeSamples = [] # resize the waveform
        if remLo > 0: removeSamples += [i for i in range(0,remLo+1)]
        if remHi > 0: removeSamples += [npArr.size - i for i in range(1,remHi+1)]
        # print "removing:",removeSamples

        self.ts = np.delete(ts,removeSamples)
        self.waveRaw = np.delete(npArr,removeSamples)   # force the size of the arrays to match

        self.noiseAvg,_,self.baseAvg = baselineParameters(self.waveRaw)
        self.waveBLSub = np.copy(self.waveRaw)
        self.waveBLSub[:] = [x - self.baseAvg for x in self.waveRaw] # subtract the baseline value

    # constants
    def GetOffset(self): return self.offset
    def GetStartTime(self): return self.start
    def GetStopTime(self): return self.stop
    def GetBins(self): return self.binsPerNS
    def GetBaseNoise(self): return self.baseAvg, self.noiseAvg
    def GetWindowIndex(self): return self.loWin, self.hiWin

    # arrays
    def GetTS(self): return self.ts
    def GetWaveRaw(self): return self.waveRaw
    def GetWaveBLSub(self): return self.waveBLSub

if __name__ == "__main__":
    main(sys.argv[1:])
