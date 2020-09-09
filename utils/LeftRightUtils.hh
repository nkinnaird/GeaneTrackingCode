//!----------------------------------
//!  This class provides useful functions for translating between UV and XY coordinates
//!----------------------------------

#ifndef LEFTRIGHTUTILS_HH
#define LEFTRIGHTUTILS_HH

#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"
#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"
#include "gm2tracker/utils/GeaneTrackUtils.hh"
#include <Eigen/Dense>
#include "fhiclcpp/ParameterSet.h"

namespace gm2strawtracker {

  class LeftRightUtils {
     
    public :
    
    LeftRightUtils();

    // return LR sides based on fit - mainly from wire fit
    void fillLRFromGeomAndTangent(gm2strawtracker::TrackDetailArtRecord & track);
    
    private : 

    gm2geom::StrawTrackerGeometry sgeom_;    
    gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;

    // fhicl parameters
    fhicl::ParameterSet fhicl_;
    int maxNumPlanes_;
  };
}

#endif
