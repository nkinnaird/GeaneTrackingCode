#include "StrawGeomUtils.hh"

namespace gm2strawtracker {

  // Functions to convert between XY & UV
  // In module coordinates, Z is along beam direction, Y is vertical and X is in radial direction (from flobber to end of module)
  // Define positive U in the direction away from the flobber (closer to +ve x) and same for V
  void StrawGeomUtils::getUVfromXY(const gm2geom::StrawTrackerGeometry& geom, const double& x, const double& y, double& u, double& v){
    u = x*cos(geom.layerAngle) - y*sin(geom.layerAngle);
    v = x*cos(geom.layerAngle) + y*sin(geom.layerAngle);
    return;
  }
  
  void StrawGeomUtils::getXYfromUV(const gm2geom::StrawTrackerGeometry& geom, double& x, double& y, const double& u, const double& v){
    x =  u/(2*cos(geom.layerAngle)) + v/(2*cos(geom.layerAngle));
    y = -u/(2*sin(geom.layerAngle)) + v/(2*sin(geom.layerAngle));
    return;
  }

  double StrawGeomUtils::getUfromXY(const gm2geom::StrawTrackerGeometry& geom, const double& x, const double& y) {
    double u, v;
    getUVfromXY(geom, x,y,u,v);
    return u;
  }
  
  double StrawGeomUtils::getVfromXY(const gm2geom::StrawTrackerGeometry& geom, const double& x, const double& y) {
      double u, v;
      getUVfromXY(geom, x,y,u,v);
      return v;
  }

  void StrawGeomUtils::getTrackXY(const StraightLineTrackArtRecord& track, double& x, double& y, const double z){
    double t = (z - track.pos.z())/track.mom.z();
    x = track.pos.x() + t*track.mom.x();
    y = track.pos.y() + t*track.mom.y();
    return;
  }


  void StrawGeomUtils::getTrackTrueXY(const StraightLineTrackArtRecord& track, double& x, double& y, const double z){
    double t = (z - track.truePos.z())/track.trueMom.z();
    x = track.truePos.x() + t*track.trueMom.x();
    y = track.truePos.y() + t*track.trueMom.y();
    return;
  }

  
}//namespace
