#include "LeftRightUtils.hh"

// Drift circle struct definition - need these to draw tangents between hits
struct DriftCircle{
  
  double uv; // U or V position
  double z;  // Z position
  double r;  // Radius
  
  DriftCircle() : uv(), z(), r(){}
  
  DriftCircle(const double uv, const double z, const double r)
    : uv(uv)
    , z(z)
    , r(r)
  {}
  
  bool operator<(const DriftCircle& circle)    {return z < circle.z;}
  
};
  
  
// Tangent defintion - for tangents between circles
struct Tangent {
    
  double m; // Gradient
  double c; // Intercept
  bool hit0_left; // LR ambiguity for hit 0
  bool hit1_left;
  int LRIndex;
    
  Tangent() : m(), c(), hit0_left(), hit1_left(), LRIndex() {}
    
  Tangent(const double m, const double c, const bool hit0_left, const bool hit1_left)
     : m(m)
     , c(c)
     , hit0_left(hit0_left)
     , hit1_left(hit1_left)
  {
     //set LR index based on other members
     if (hit0_left && hit1_left) LRIndex = 1;
     if (hit0_left && !hit1_left) LRIndex = 2;
     if (!hit0_left && hit1_left) LRIndex = 3;
     if (!hit0_left && !hit1_left) LRIndex = 4;
  }
}; // Tangent

vector<Tangent> calculateTangents(vector<DriftCircle>& circles) {
    
  vector<Tangent> calcTangents;
    
  // Get two circles
  DriftCircle* circle0 = &circles.at(0);
  DriftCircle* circle1 = &circles.at(1);
    
  // Parameters describing positions of circles - put in terms of u here, but equally applicable for v
  double du = circle1->uv - circle0->uv;
  double dz = circle1->z - circle0->z;
  double r0 = circle0->r;
  double r1 = circle1->r;
  double du_dz = du/dz;
  double centreDist = sqrt( du*du + dz*dz );
  
  // Have four tangents to 2 circles
  int nTangents = 4;
    
  // Unless we've passed just the wire values when there's only one
  if(r0 == 0 and r1 == 0) nTangents = 1;
    
  // Calculate each one by one in this four loop (shameless stolen from Wikipedia article)
  for (int iTang = 0; iTang < nTangents; iTang++){
      
    // These two parameters define which tangent we're taking
    int sgn;
    double radDiff;
    switch (iTang) {
      case 0:
	sgn = 1;
	radDiff = fabs(r1) - fabs(r0);
	break;
      case 1:
	sgn = -1;
	radDiff = fabs(r1) - fabs(r0);
	break;
      case 2:
	sgn = -1;
	radDiff = -fabs(r1) - fabs(r0);
	break;
      case 3:
	sgn = 1;
	radDiff = -fabs(r1) - fabs(r0);
	break;
    }
      
    // This sqrt shows up a lot - just redefined it to neaten things up
    double sqrtRad = sqrt(pow(centreDist/radDiff,2)-1);
    
    // Calculate gradient and intercept of line (U = m*Z + c)
    double m = (1 + sgn*du_dz*sqrtRad) / (sgn*sqrtRad - du_dz);
    double c = circle0->uv - m*circle0->z - fabs(r0)/(radDiff/(centreDist*centreDist) * (du-sgn*dz*sqrtRad));
    
    // Fix special case where radii are equal (and looking at external tangents)
    if (radDiff == 0){
	m = du_dz;
	c = circle0->uv - m*circle0->z - fabs(r0)/(sgn*dz/centreDist);
    }
      
    // Calculate whether we passed on the left or right side of each wire
    bool hit0_left = false;
    bool hit1_left = false;
    if (m*circle0->z + c > circle0->uv) hit0_left = true;
    if (m*circle1->z + c > circle1->uv) hit1_left = true;
      
    // Push back this gradient to vector we're going to return
    calcTangents.push_back(Tangent(m, c, hit0_left, hit1_left));
  }
    
  return calcTangents;
}


namespace gm2strawtracker {

  LeftRightUtils::LeftRightUtils() 
    : sgeom_()
    , geaneTrackUtils_()
    , fhicl_()
  {
    maxNumPlanes_ = sgeom_.getNumLayersPerStation()+1; // 33
  }
  
  void LeftRightUtils::fillLRFromGeomAndTangent(gm2strawtracker::TrackDetailArtRecord & track)
  {
    
    // fill with what we had from previous wire fit and overwrite with geom
    // if track hasn't been filled with the geometry this will not fill any right left info
    std::vector<gm2strawtracker::GeaneHitSide> hitSides = track.geaneHits.geaneHitSides; 

    std::vector<gm2strawtracker::GeaneHitSide> origSides;
    for (auto s : hitSides) {
      origSides.push_back(s);
    } 
   
    map<int, vector<Tangent> > tangentsMap;
    map<int, double> lowestDCAOfDoublet;
    vector<int> planesThatAreDoublets;
    map<int, double> planeAngleUV;
    
    //loop over and find doublets
    for (auto i : {1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31}){
      //bool isDoublet = false;
      bool isDoublet = (track.geaneHits.geaneHitSides.at(i) < 3 && track.geaneHits.geaneHitSides.at(i+1) < 3);
      if (isDoublet) {

	// get the drift circles using the calculated dca
	double zPosMeas1 = track.geaneHits.planeZPositions.at(i);
	double zPosMeas2 = track.geaneHits.planeZPositions.at(i+1);
	
	double wirePosUV1 = (geaneTrackUtils_.isUPlane(i))   ? track.geaneHits.geaneWireUVPositions[3][i]   : track.geaneHits.geaneWireUVPositions[4][i];
	double wirePosUV2 = (geaneTrackUtils_.isUPlane(i+1)) ? track.geaneHits.geaneWireUVPositions[3][i+1] : track.geaneHits.geaneWireUVPositions[4][i+1];
	
	DriftCircle dc1 = DriftCircle( wirePosUV1, zPosMeas1, track.measuredDCAs.at(i));
	DriftCircle dc2 = DriftCircle( wirePosUV2, zPosMeas2, track.measuredDCAs.at(i+1));
	
	std::vector<DriftCircle> dcv = {dc1, dc2};

	// get all 4 possible tangents given the drift circles
	tangentsMap[i] = calculateTangents(dcv);
	tangentsMap[i+1] = tangentsMap[i];
	
	double lowestDCA = track.measuredDCAs.at(i);
	if (track.measuredDCAs.at(i+1) < lowestDCA) lowestDCA = track.measuredDCAs.at(i+1);
	lowestDCAOfDoublet[i] = lowestDCA;
	lowestDCAOfDoublet[i+1] = lowestDCA;
	
	planesThatAreDoublets.push_back(i);
	planesThatAreDoublets.push_back(i+1);
      }
      //else{
      //	// otherwise push back unknown to hit side?
      //	hitSides.at(i) = gm2strawtracker::GeaneHitSide::gUnknown;
      //	hitSides.at(i+1) = gm2strawtracker::GeaneHitSide::gUnknown;
      //}
    }	
    
    for (auto planeNum : planesThatAreDoublets){

      auto geaneHitsOnTrack = track.geaneHits;

      double predMomU = geaneTrackUtils_.getPredUMom(geaneHitsOnTrack, planeNum);
      double predMomV = geaneTrackUtils_.getPredVMom(geaneHitsOnTrack, planeNum);
      double predMomX = geaneTrackUtils_.getPredXMom(geaneHitsOnTrack, planeNum);	
      
      planeAngleUV[planeNum] = (geaneTrackUtils_.isUPlane(planeNum))? atan2( predMomU, predMomX ) : atan2( predMomV, predMomX );
      
      double angleDiff = 3.14;
      int indexOfClosestTangent = 0;
      
      for (auto t : tangentsMap[planeNum]){
	if (fabs(atan(t.m) - planeAngleUV[planeNum]) < fabs(angleDiff)) {
	  angleDiff = fabs(atan(t.m) - planeAngleUV[planeNum]);
	  indexOfClosestTangent = t.LRIndex;
	}
      }

      gm2strawtracker::GeaneHitSide side = gm2strawtracker::GeaneHitSide::gUnknown;
      if (indexOfClosestTangent == 1){ //LL
	side = gm2strawtracker::GeaneHitSide::gLeft;
      }
      if (indexOfClosestTangent == 2) { //LR
	if (planeNum % 2) side = gm2strawtracker::GeaneHitSide::gLeft;
	else side = gm2strawtracker::GeaneHitSide::gRight;
      }
      if (indexOfClosestTangent == 3) { //RL
	if (planeNum % 2) side = gm2strawtracker::GeaneHitSide::gRight;
	else side = gm2strawtracker::GeaneHitSide::gLeft;
      }
      if (indexOfClosestTangent == 4) { //RR
	side = gm2strawtracker::GeaneHitSide::gRight;
      }
      
      //double wireDCA = fhicl_.get<double>("wireDCA", 0.5); //default to 500 microns
      
      //if wire fit (side) and geom (set in hitSides) match don't match set it to unknown
      if (side != hitSides.at(planeNum)) side = gm2strawtracker::GeaneHitSide::gUnknown;
      
      // set to wire centre if DCA less than wireDCA value
      //if (track.measuredDCAs.at(planeNum) < wireDCA) side = gm2strawtracker::GeaneHitSide::gCenter;
      
      hitSides.at(planeNum) = side;
      
    }

    //int iS = 0;
    //for (auto s : hitSides) {
    //  if (origSides.at(iS) != s) cout << "SIDE CHECK: " << iS << " had side: " << origSides.at(iS) << ", now: " << s << "\n";
    //  iS++;
    //}
    
    track.geaneHits.geaneHitSides = hitSides;
    return;
  }
  

}//namespace
