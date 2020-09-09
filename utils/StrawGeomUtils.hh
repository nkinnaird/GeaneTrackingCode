//!----------------------------------
//!  This class provides useful functions for translating between UV and XY coordinates
//!----------------------------------

#ifndef STRAWGEOMUTILS_HH
#define STRAWGEOMUTILS_HH

#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"
#include "gm2dataproducts/strawtracker/StraightLineTrackArtRecord.hh"

namespace gm2strawtracker {

  class StrawGeomUtils {
     
    public :

     //Functions to convert between UV & XY coordinates - pass x,y in station or module coordinates
     void getUVfromXY(const gm2geom::StrawTrackerGeometry& geom, const double& x, const double& y, double& u, double& v);
     void getXYfromUV(const gm2geom::StrawTrackerGeometry& geom, double& x, double& y, const double& u, const double& v);
     double getUfromXY(const gm2geom::StrawTrackerGeometry& geom,const double& x, const double& y);
     double getVfromXY(const gm2geom::StrawTrackerGeometry& geom, const double& x, const double& y);

     // Functions to XY position for a given z position when using StraightLineTrackArtRecords
     void getTrackXY(const StraightLineTrackArtRecord& track, double& x, double& y, const double z);
     void getTrackTrueXY(const StraightLineTrackArtRecord& track, double& x, double& y, const double z);

  };

}

#endif
