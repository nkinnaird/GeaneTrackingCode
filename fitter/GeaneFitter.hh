//!--------------------------------------------------------------------------------
//! This is a least-squares minimizer in 3-dimensional space which runs off Geane objects.
//! Nick Kinnaird - nickknn@bu.edu
//!---------------------------------------------------------------------------------

#ifndef GEANEFitter_HH
#define GEANEFitter_HH

#include "TrackFitter.hh"

#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"

#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"
#include "gm2tracker/utils/GeaneTrackUtils.hh"

#include <Eigen/Dense>

namespace gm2strawtracker {

   class GeaneFitter {

      public :

        GeaneFitter(fhicl::ParameterSet const & p);

        double TrackCorrelation(gm2strawtracker::TrackDetailArtRecord & track);

        void makeHybridErrMat(gm2strawtracker::TrackDetailArtRecord & trackDetails, bool Uhybrid);
        double sequenceChecking(gm2strawtracker::TrackDetailArtRecord & trackDetails);

        void convertToGeVcm(Eigen::VectorXd& myVector);
        void convertToMeVmm(Eigen::VectorXd& myVector);
      
      private :

        bool matrixDebug_; // For turning on copious matrix verbose outputs.
        
        std::string name_;

        gm2geom::StrawTrackerGeometry sgeom_;
        gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;

        int maxNumPlanes;
        double wireMeasurementError; // uniform distribution standard deviation - not gaussian but approximated as such when fitting hits
        int numTrackParams;

        std::vector<unsigned int> skipLayers_; // Layers that we'll ignore in the fit

   };

}

#endif // GEANEFitter
