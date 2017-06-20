#ifndef GATNonLinearityCorrector_hh
#define GATNonLinearityCorrector_hh

/** Class Description:
 * This class applies an ADC non-linearity correction to a
 * TClonesArray of waveforms and puts the output in its protected internal
 * TClonesArray of waveforms. The input array can be the raw (default) or aux
 * waveforms in MGTEvent, or a published TClonesArray of waveforms. The internal
 * TClonesArray can be published for use by other processors, and/or
 * post-processed by a derived class.
 */

#include "GATWaveformTransformer.hh"
#include "MGWFNonLinearityCorrector.hh"
#include <map>

class GATNonLinearityCorrector : public GATWaveformTransformer 
{
  public:
    /// Constructor. By default runs on raw WFs in MGTEvent and posts nothing
    /** transform must be non-NULL.
     * Set nlcDBPath to the path to list (database) of non-linearity
     *   correction up/down files.
     * Set postedWFNameToProcess to non-NULL to run on a posted TClonesArray of
     *   waveforms.
     * Set wfNameForPosting to non-NULL to post the internal TClonesArray of
     *   waveforms to be accessed by other TAModules.
     * Set useAuxWFs to "true" when postedWFNameToProcess = wfNameForPosting =
     *   NULL to run on the auxiliary wf array in MGTEvent rather than the (raw)
     *   waveforms in the main array. */
    GATNonLinearityCorrector(MGWFNonLinearityCorrector& nonLinCorrTransformer,
                             const char* nlcDBPath,
                             const char* postedWFNameToProcess = NULL, 
                             const char* wfNameForPosting = NULL, 
                             bool useAuxWFs = false,
                             bool skipBadWFs = true);
    virtual ~GATNonLinearityCorrector() {}
		
  protected:
    virtual void LoadParameters(int ddID, int run);

  protected:
    MGWFNonLinearityCorrector* fNLCTransformer;
    std::string fNLCMapDir;
    std::map<int, MGWFNonLinearityCorrectionMap*> fNLCMaps;
};

#endif
