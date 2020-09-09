//!-----------------------------------------------------------------------------------------------------------------------
//! This is a base class for the various derived classes that support fitting routines in the gm2tracker framework.
//! Tools that are designed for fitting should inhert from it and implement the virtual method(s) :
//!
//!  --makeFittedTrack(), which takes in a track candidate and the type of fitting routine. The return is a 
//!    boolean which tells if the fit was successfully ( converges ).
//!
//!    Consider two possible cases: 
//!    First case, a routine can create a track candidate based on track finding algorithms and use this function 
//!    to fit the track, and the function stores the fit parameters onto the TrackCandidate Art Record.
//!    Second case, a routine refits the track candidate with a different fitting routine
//!    and the the original fitting parameters are replaced.
//!
//!  Note that this base class forces users to fit only track candidates. The primary reason is
//!  to prevent developers from creating the reconstructed tracks in the fitting rountines. Fitters are
//!  designed to only fit and not build objects. In addition, this will prevent a circular dependency
//!  among the folders in the gm2tracker package.
//!---------------------------------------------------------------------------------------------------------------------

#ifndef TRACKFITTER_HH
#define TRACKFITTER_HH

// art includes
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// root includes
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"

// constants includes
#include "gm2geom/common/Gm2Constants_service.hh"

// field includes
#include "gm2geom/fields/gm2FieldManager_service.hh"

// geometry
#include "gm2geom/coordSystems/CoordSystem.hh"
#include "gm2geom/coordSystems/CoordSystemsStoreData.hh"
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh"  

#include "gm2geom/Core/GeometryBase.hh"
#include "gm2geom/Core/Geometry_service.hh"
#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"

// data product includes
#include "gm2dataproducts/strawtracker/Helix.hh"
#include "gm2dataproducts/strawtracker/StrawSeedArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackCandidateArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"

// local
#include "LinearTrack.hh"
#include "CircularTrack.hh"
#include "HelicalTrack.hh"

#include "gm2tracker/utils/StrawObjectSorter.hh"


namespace gm2strawtracker {

   template <typename StateVector>

   class TrackFitter {

      public :

        // standard constructor
        TrackFitter();

        // make fit the track candidate and checks if the fit is successfully
        virtual bool makeFittedTrack( std::vector<gm2geom::CoordSystem3Vector>& points,
			              std::vector<gm2geom::CoordSystem3Vector>& errors,
                                      StateVector& trackState ) = 0;  

        // set the fhicl parameters
        inline void setFhiclCuts( fhicl::ParameterSet const & p ) {
            fhicl_ = p;
        }

        // get the fhicl parameters
        inline fhicl::ParameterSet getFhiclCuts() {
            return fhicl_;
        }

        // set the wire ID
	inline void setWireIDVector( gm2strawtracker::WireID const & wire ) {
            wires_.push_back( wire );
	}

	// get the wire ID
	inline std::vector< gm2strawtracker::WireID > getWires() {
	    return wires_;
	}

        // set the coordinate system data
        inline void setCoordSysData( gm2geom::CoordSystemsStoreData const & c ) {
            cs_ = c;
        }  

	// get the coordinate system data
	inline gm2geom::CoordSystemsStoreData getCoordSysData() {
            return cs_;
	}

	// get the geometry tool
	inline gm2geom::StrawTrackerGeometry getGeometryTool() {
            return geom_;
        }

        // get the constant service handle
        inline art::ServiceHandle<gm2geom::Gm2Constants> getConstantService() {
            return gm2consts_;
        }            

        // get the field service handle
	inline art::ServiceHandle<gm2geom::gm2FieldManager> getFieldService() {
            return gm2fields_;
	}

        // set the predicted states for the track candidate
        inline void setCandidateStates( TrackDetailStates& istates ) {
            states_.clear(); 
            states_ = istates; 
        }

        // get the predicted states for the track candidate
        inline TrackDetailStates getCandidateStates() { 
            return states_; 
        }

        // set the iteration number
        inline void setIterations( int& i ) {
           iterations_ = i;
        }

        // get the iteration number
        inline int getIterations() {
           return iterations_;
        }

        inline void setFieldType( int& i ) {
           fieldType_ = i;
	}

        // get the field type
	inline int getFieldType() {
           return fieldType_;
	}

      private :

        // a container to store the predicted states resulting from a Kalman fit
        TrackDetailStates states_;

        // iterations
        int iterations_;

        // what field type
	int fieldType_;

        // fhicl parameters
        fhicl::ParameterSet fhicl_;
  
	// wires
	std::vector< gm2strawtracker::WireID > wires_;

        // services tools
	gm2geom::CoordSystemsStoreData cs_;
        gm2geom::StrawTrackerGeometry geom_;
        art::ServiceHandle<gm2geom::Gm2Constants> gm2consts_;
        art::ServiceHandle<gm2geom::gm2FieldManager> gm2fields_;

   };

}

#endif // TRACKFITTER
