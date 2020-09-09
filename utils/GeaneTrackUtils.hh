//!----------------------------------
//!  This class provides useful functions for translating between UV and XY coordinates as well as getting track parameters/sides
//!
//!  Added overloaded functions to handle the multiple track art record
//!----------------------------------

#ifndef GEANETRACKUTILS_HH
#define GEANETRACKUTILS_HH

#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"

#include <Eigen/Dense>

#include "artg4/gm2Geane/gm2GeanePropagator.hh"
#include "artg4/gm2Geane/gm2GeanePropagatorManager.hh"

namespace gm2strawtracker {

  class GeaneTrackUtils {
     
    public :

     GeaneTrackUtils();

     bool isUPlane(int planeNum); // true if U plane, false if V plane, where planeNum runs from 1 to 32

     // get predicted parameters from GEANEHits
     double  getPredMom (const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredMom (const gm2strawtracker::GEANEHits& track, int planeNum);

     double  getPredXMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredXMom(const gm2strawtracker::GEANEHits& track, int planeNum);

     double  getPredYMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredYMom(const gm2strawtracker::GEANEHits& track, int planeNum);

     double  getPredZMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredZMom(const gm2strawtracker::GEANEHits& track, int planeNum);

     double  getPredUMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredUMom(const gm2strawtracker::GEANEHits& track, int planeNum);

     double  getPredVMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredVMom(const gm2strawtracker::GEANEHits& track, int planeNum);

     double  getPredXPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredXPos(const gm2strawtracker::GEANEHits& track, int planeNum);

     double  getPredYPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredYPos(const gm2strawtracker::GEANEHits& track, int planeNum);

     double  getPredZPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredZPos(const gm2strawtracker::GEANEHits& track, int planeNum);
     
     double  getPredUPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredUPos(const gm2strawtracker::GEANEHits& track, int planeNum);

     double  getPredVPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum);
     double  getPredVPos(const gm2strawtracker::GEANEHits& track, int planeNum);


     int maxNumPlanes;

     Eigen::Matrix2d XYtoUVcoordinateTransformationMatrix;
     Eigen::MatrixXd JacobianToUV;

     gm2GeanePropagatorManager* g4emgr;

    private : 

     gm2geom::StrawTrackerGeometry sgeom_;

  };

}

#endif
