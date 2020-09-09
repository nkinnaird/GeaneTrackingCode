#include "GeaneTrackUtils.hh"

namespace gm2strawtracker {

GeaneTrackUtils::GeaneTrackUtils() 
  : sgeom_()
{
  maxNumPlanes = sgeom_.getNumLayersPerStation()+1; // 33

  XYtoUVcoordinateTransformationMatrix(0,0)=cos(sgeom_.layerAngle);
  XYtoUVcoordinateTransformationMatrix(0,1)=-sin(sgeom_.layerAngle);
  XYtoUVcoordinateTransformationMatrix(1,0)=cos(sgeom_.layerAngle);
  XYtoUVcoordinateTransformationMatrix(1,1)=sin(sgeom_.layerAngle);

  JacobianToUV = Eigen::MatrixXd::Zero(5,5);
  JacobianToUV(0,0) = 1.;
  JacobianToUV.block(1,1,2,2) = XYtoUVcoordinateTransformationMatrix;
  JacobianToUV.block(3,3,2,2) = XYtoUVcoordinateTransformationMatrix;
}

bool GeaneTrackUtils::isUPlane(int planeNum){
  if (((planeNum-1)/2)%2 == 0) return true;
  else return false;
}

// get the total momentum using GEANE struct from TrackDetailArtRecord 
double GeaneTrackUtils::getPredMom (const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    if (planeNum == 0) return sqrt(getPredXMom(track, planeNum)*getPredXMom(track, planeNum) + getPredYMom(track, planeNum)*getPredYMom(track, planeNum) + getPredZMom(track, planeNum)*getPredZMom(track, planeNum));
    else return 1./track.geanePredictedParameters[0].at(planeNum);
}

// get the total momentum using GEANE struct from TrackArtRecord
double GeaneTrackUtils::getPredMom (const gm2strawtracker::GEANEHits& track, int planeNum){
    if (planeNum == 0) return sqrt(getPredXMom(track, planeNum)*getPredXMom(track, planeNum) + getPredYMom(track, planeNum)*getPredYMom(track, planeNum) + getPredZMom(track, planeNum)*getPredZMom(track, planeNum));
    else return 1./track.geanePredictedParameters[0].at(planeNum);
}


// get the momentum in the X-direction from TrackDetailArtRecord
double GeaneTrackUtils::getPredXMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(3);
    else return track.geanePredictedParameters[1].at(planeNum) * getPredZMom(track, planeNum);
}

// get the momentum in the X-direction from TrackArtRecord
double GeaneTrackUtils::getPredXMom(const gm2strawtracker::GEANEHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(3);
    else return track.geanePredictedParameters[1].at(planeNum) * getPredZMom(track, planeNum);
}

// get the momentum in the Y-direction from TrackDetailArtRecord
double GeaneTrackUtils::getPredYMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(4);
    else return track.geanePredictedParameters[2].at(planeNum) * getPredZMom(track, planeNum);
}

// get the momentum in the Y-direction from TrackArtRecord
double GeaneTrackUtils::getPredYMom(const gm2strawtracker::GEANEHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(4);
    else return track.geanePredictedParameters[2].at(planeNum) * getPredZMom(track, planeNum);
}

// get the momentum in the Z-direction from TrackDetailArtRecord
double GeaneTrackUtils::getPredZMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(5);
    else return getPredMom(track, planeNum)/sqrt(1 + track.geanePredictedParameters[1].at(planeNum)*track.geanePredictedParameters[1].at(planeNum) + track.geanePredictedParameters[2].at(planeNum)*track.geanePredictedParameters[2].at(planeNum));
}

// get the momentum in the Z-direction from TrackArtRecord
double GeaneTrackUtils::getPredZMom(const gm2strawtracker::GEANEHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(5);
    else return getPredMom(track, planeNum)/sqrt(1 + track.geanePredictedParameters[1].at(planeNum)*track.geanePredictedParameters[1].at(planeNum) + track.geanePredictedParameters[2].at(planeNum)*track.geanePredictedParameters[2].at(planeNum));
}

// get the momenum in the U-direction fromm TrackDetailArtRecord
double GeaneTrackUtils::getPredUMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    return XYtoUVcoordinateTransformationMatrix(0,0)*getPredXMom(track, planeNum) + XYtoUVcoordinateTransformationMatrix(0,1)*getPredYMom(track, planeNum);
}

// get the momenum in the U-direction fromm TrackArtRecord
double GeaneTrackUtils::getPredUMom(const gm2strawtracker::GEANEHits& track, int planeNum){
    return XYtoUVcoordinateTransformationMatrix(0,0)*getPredXMom(track, planeNum) + XYtoUVcoordinateTransformationMatrix(0,1)*getPredYMom(track, planeNum);
}

// get the momentum in the V-direction from TrackDetailArtRecord
double GeaneTrackUtils::getPredVMom(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    return XYtoUVcoordinateTransformationMatrix(1,0)*getPredXMom(track, planeNum) + XYtoUVcoordinateTransformationMatrix(1,1)*getPredYMom(track, planeNum);
}

// get the momentum in the V-direction from TrackArtRecord
double GeaneTrackUtils::getPredVMom(const gm2strawtracker::GEANEHits& track, int planeNum){
    return XYtoUVcoordinateTransformationMatrix(1,0)*getPredXMom(track, planeNum) + XYtoUVcoordinateTransformationMatrix(1,1)*getPredYMom(track, planeNum);
}

// get the position in the X-direction from TrackDetailArtRecord
double GeaneTrackUtils::getPredXPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(0);
    else return track.geanePredictedParameters[3][planeNum];
}

// get the position in the X-direction from TrackArtRecord
double GeaneTrackUtils::getPredXPos(const gm2strawtracker::GEANEHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(0);
    else return track.geanePredictedParameters[3][planeNum];
}

// get the position in the Y-direction from TrackDetailArtRecord
double GeaneTrackUtils::getPredYPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(1);
    else return track.geanePredictedParameters[4][planeNum];
}

// get the position in the Y-direction from TrackArtRecord
double GeaneTrackUtils::getPredYPos(const gm2strawtracker::GEANEHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(1);
    else return track.geanePredictedParameters[4][planeNum];
}

// get the position in the Z-direction from TrackDetailArtRecord
double GeaneTrackUtils::getPredZPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(2);
    else return track.planeZPositions.at(planeNum);
}

// get the position in the Z-direction from TrackArtRecord
double GeaneTrackUtils::getPredZPos(const gm2strawtracker::GEANEHits& track, int planeNum){
    if (planeNum == 0) return track.startingGeaneParameters.at(2);
    else return track.planeZPositions.at(planeNum);
}

// get the position in the U-direction from TrackDetailArtRecord
double GeaneTrackUtils::getPredUPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    return XYtoUVcoordinateTransformationMatrix(0,0)*getPredXPos(track, planeNum) + XYtoUVcoordinateTransformationMatrix(0,1)*getPredYPos(track, planeNum);
}

// get the position in the U-direction from TrackArtRecord
double GeaneTrackUtils::getPredUPos(const gm2strawtracker::GEANEHits& track, int planeNum){
    return XYtoUVcoordinateTransformationMatrix(0,0)*getPredXPos(track, planeNum) + XYtoUVcoordinateTransformationMatrix(0,1)*getPredYPos(track, planeNum);
}

// get the position in the V-direction from TrackDetailArtRecord
double GeaneTrackUtils::getPredVPos(const gm2strawtracker::GEANEDetailHits& track, int planeNum){
    return XYtoUVcoordinateTransformationMatrix(1,0)*getPredXPos(track, planeNum) + XYtoUVcoordinateTransformationMatrix(1,1)*getPredYPos(track, planeNum);
}

// get the position in the V-direction from TrackArtRecord
double GeaneTrackUtils::getPredVPos(const gm2strawtracker::GEANEHits& track, int planeNum){
    return XYtoUVcoordinateTransformationMatrix(1,0)*getPredXPos(track, planeNum) + XYtoUVcoordinateTransformationMatrix(1,1)*getPredYPos(track, planeNum);
}


}//namespace
