#ifndef SIMPLETRACKFINDINGUTILS_HH
#define SIMPLETRACKFINDINGUTILS_HH

//---------------------------------------------------------
// This fitter takes all of the hits and create a track
// candidate
// --------------------------------------------------------

// art includes
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// geometry
#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"

// coordinate system includes
#include "gm2geom/coordSystems/CoordSystem.hh"
#include "gm2geom/coordSystems/CoordSystemsStoreData.hh"
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh"  

// data product includes
#include "gm2dataproducts/strawtracker/StrawTimeIslandArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackCandidateArtRecord.hh"

// local
#include "gm2tracker/fitter/SimpleCircleFitter.hh"
#include "gm2tracker/utils/StrawObjectSorter.hh"

#include "boost/format.hpp"

namespace gm2strawtracker {

   class SimpleTrackFindingUtils {

      public :

       // default constructor
       SimpleTrackFindingUtils();

       // merge all hits on an island to form a track candidate
       TrackCandidateArtRecord findTrackCandidates( const StrawTimeIslandArtRecord& island ); 


       // set the starting (pivot) point for the circle fitter
       inline void setStartingPosition( gm2geom::CoordSystem3Vector& position ) {
         position_ = position;
       }

       // set the coordinate system data
       inline void setCoordSysData( gm2geom::CoordSystemsStoreData const & c ) {
         cs_ = c;
       }

       // set fhicl parameters
       inline void setFhiclCuts( fhicl::ParameterSet const & p ) {
         fhicl_ = p;
       }
 
      private :

       // helper tools
       SimpleCircleFitter fitter_;

       // message facility name
       std::string name_;

       // coordinate system tools
       gm2geom::CoordSystemsStoreData cs_;
 
       // fhicl parameters
       fhicl::ParameterSet fhicl_;

       // the starting position
       gm2geom::CoordSystem3Vector position_;
   };

}

#endif
