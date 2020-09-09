//!--------------------------------------------------------------------------------
//! This is a least-squares minimizer in 3-dimensional space which runs off Geane objects.
//! Nick Kinnaird - nickknn@bu.edu
//!---------------------------------------------------------------------------------
//! Edited to fit using 5x5 matrix inversion at each plane without considering material error
//! This reduced fitter is run without loading the Straws or StrawTrackerCadMesh in the RunGeane.fcl file
//! In GeaneFittingUtils.cc, call Geane5x5Fitter::TrackCorrelation() for eventChiSquared.
//! Calling the 5x5 fitter follows the same structure as the NxN fitter in GeaneFittingUtils.cc and its header file.
//! Deepak Sathyan - dsathyan@bu.edu

#ifndef GEANE5x5Fitter_HH
#define GEANE5x5Fitter_HH

#include "TrackFitter.hh"

#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"
#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"
#include "gm2tracker/utils/GeaneTrackUtils.hh"

#include <Eigen/Dense>

namespace gm2strawtracker {

   class Geane5x5Fitter {

      public :

        //! constructor
        Geane5x5Fitter(fhicl::ParameterSet const & p);

        double TrackCorrelation(gm2strawtracker::TrackDetailArtRecord& trackDetails);

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

#endif // GEANE5x5Fitter
