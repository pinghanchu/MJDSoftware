include $(MGDODIR)/buildTools/config.mk

# Give the list of applications, which must be the stems of cc files with 'main'.  There
# can be more than one.  In our example, this means there is a test.cc and a test1.cc
APPS = process_mjd_cal fillhist fillhistchan combhist combhistlist combhistsum combhistchan combhistsumlist calhist calhistchan calhistsum resohist resohistsum multihist multihistsum pulserfinder combpulser plotcalibration plotpulser plotpulserhisto plotezfit plotezfitsum plotezfitrun plotrun profile_time plotprofile plotcovariance plotresocurve plotresocurvesum plotresogrid2 plotroi plotroisum uploadlinearfit checklinearfit checkezfit gettime
# The next two lines are important 
SHLIB = 
ARCHIVE = 
TAMDIR ?= $(ROOTSYS)
# Include the correct flags, 
INCLUDEFLAGS = $(CLHEP_INCLUDE_FLAGS) -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms 
INCLUDEFLAGS += -I$(MGDODIR)/Majorana -I$(MGDODIR)/MJDB $(ROOT_INCLUDE_FLAGS) -I$(TAMDIR)/inc -I$(TAMDIR)/include -I$(MGDODIR)/Tabree
INCLUDEFLAGS += -I$(GATDIR)/BaseClasses -I$(GATDIR)/MGTEventProcessing -I$(GATDIR)/MGOutputMCRunProcessing -I$(GATDIR)/Analysis -I$(GATDIR)/MJDAnalysis -I$(GATDIR)/DCProcs
LIBFLAGS = -L$(MGDODIR)/lib -lMGDORoot -lMGDOBase -lMGDOTransforms -lMGDOMajorana -lMGDOGerdaTransforms -lMGDOMJDB -lMGDOTabree
LIBFLAGS += -L$(GATDIR)/lib -lGATBaseClasses -lGATMGTEventProcessing -lGATMGOutputMCRunProcessing -lGATAnalysis -lGATMJDAnalysis -lGATDCProcs $(ROOT_LIB_FLAGS) -lSpectrum -L$(TAMDIR)/lib -lTAM

include $(MGDODIR)/buildTools/BasicMakefile
