//!----------------------------------
//!  This class provides functions for the fitting of geane tracks
//!----------------------------------

#ifndef GEANEFITTINGUTILS_HH
#define GEANEFITTINGUTILS_HH

#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"

#include "gm2tracker/utils/GeaneParamUtils.hh"
#include "gm2tracker/utils/GeaneTrackUtils.hh"
#include "gm2tracker/utils/GeaneLRUtils.hh"
#include "gm2tracker/utils/LeftRightUtils.hh"

#include "gm2tracker/fitter/GeaneFitter.hh"
#include "gm2tracker/fitter/Geane5x5Fitter.hh"

namespace gm2strawtracker {

  class GeaneFittingUtils {
     
    public :

     GeaneFittingUtils(fhicl::ParameterSet const & p);

     int truthLRFit(gm2strawtracker::TrackDetailArtRecord & trackFitDetails);
     int wireFit(gm2strawtracker::TrackDetailArtRecord & trackFitDetails);
     int mainFit(gm2strawtracker::TrackDetailArtRecord & trackFitDetails);
     int fullSeqFit(gm2strawtracker::TrackDetailArtRecord & trackFitDetails);

     int fittingLoop(gm2strawtracker::TrackDetailArtRecord & trackFitDetails, int maxNumPasses, bool updateLR);

     int checkExtraneousFailureModes(gm2strawtracker::TrackDetailArtRecord & trackFitDetails);

     gm2strawtracker::GeaneParamUtils geaneParamUtils_;

    private : 

     gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;
     gm2strawtracker::GeaneLRUtils geaneLRUtils_;
     gm2strawtracker::LeftRightUtils LRUtils_;

     gm2strawtracker::GeaneFitter matrixCalculations_;
     gm2strawtracker::Geane5x5Fitter reducedMatrixCalculations_;

     std::string name_;
     string fitMode_;
     string fitterType_;

     double lockLowDCAs_;
     bool useNodes_;

     int numPassesWireFit_; // Number of passes to run for the wire fit - 3 is probably enough.
     int numPassesSeqFit_; // Same but for the full seq fit.

     double convergenceCriteria_;
    
     bool useTangentLR_;

  };

}

#endif
