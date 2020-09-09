#include "GeaneDummyUtils.hh"

namespace gm2strawtracker {

GeaneDummyUtils::GeaneDummyUtils(fhicl::ParameterSet const & p) 
  : rseed_(p.get<int>("rseed",0))
  , xPosChange_(p.get<double>("xPosChange",5.0))
  , yPosChange_(p.get<double>("yPosChange",5.0))
  , xMomChange_(p.get<double>("xMomChange",100.0))
  , yMomChange_(p.get<double>("yMomChange",100.0))
  , zMomChange_(p.get<double>("zMomChange",100.0))
  , geaneTrackUtils_()
{
  r3.SetSeed(rseed_);
}

// Take input hits and put in map by eventNumInFill or it's painfully slow to associate digits
void GeaneDummyUtils::fillDummyHits(art::Handle<gm2truth::GhostDetectorArtRecordCollection> DummyDataHandle){

  // Clear anything already in our dummyHits map
  dummyHits.clear();

  // Loop over input art::Handle and split hits into map to make searching easier
  for(unsigned int index = 0; index < DummyDataHandle->size(); index++){
    auto &thisHit = DummyDataHandle->at(index);
    auto mapIter = dummyHits.find(thisHit.eventNumInFill);
    if(mapIter == dummyHits.end()){
      dummyHits[thisHit.eventNumInFill] = { pair<unsigned int,gm2truth::GhostDetectorArtRecord>(index,thisHit) };
    } else {
      mapIter->second.push_back(pair<unsigned int,gm2truth::GhostDetectorArtRecord>(index,thisHit));
    }
  }

}

void GeaneDummyUtils::clearDummyHolders(){
  sortedDummyPlaneHits.clear();
  dummyPlanesHit.clear();
}

void GeaneDummyUtils::getTrackDummyHits(art::Handle<gm2truth::GhostDetectorArtRecordCollection> DummyDataHandle, gm2strawtracker::TrackDetailArtRecord & trackArtRecord){

    std::vector< art::Ptr< gm2truth::GhostDetectorArtRecord> > dummyHitPtrs;

    for (int i = 0; i < int(sortedDummyPlaneHits.size()); ++i)
    {
        art::Ptr< gm2truth::GhostDetectorArtRecord > dHitPtr(DummyDataHandle, sortedDummyPlaneHits.at(i).second);
        dummyHitPtrs.push_back(dHitPtr);
    }

    trackArtRecord.dummyPlaneHits = dummyHitPtrs;  
}             

void GeaneDummyUtils::fillDummyHitInfo(gm2truth::GhostDetectorArtRecord dummyHit){

  string pvName = dummyHit.volumeName;

  std::string delimiter1 = "[";
  size_t pos = 0;
  
  std::vector< string > token;
  while ((pos = pvName.find(delimiter1)) != std::string::npos) {
    token.push_back(pvName.substr(0,pos-1));
    pvName.erase(0, pos + delimiter1.length());
  }
  token.push_back(pvName.substr(0, pvName.size()-1));
  token.erase(token.begin());

  dummystationNumber = atoi(token.at(0).c_str());
  dummymoduleNum = atoi(token.at(1).c_str());
  dummysubPlaneNum = atoi(token.at(2).c_str());

  dummytotalplaneNum = 1+(dummymoduleNum*4)+dummysubPlaneNum;

}

void GeaneDummyUtils::fillDigitHitInfo(const art::Ptr<gm2strawtracker::StrawDigitArtRecord> *digitHit){

    int moduleNum = (*digitHit)->wireID.getModule();
    int subPlaneNum = (*digitHit)->wireID.getLayer() + 2 * (*digitHit)->wireID.getView(); // layer number 0 to 1 in each layer, view number = 0 for u and 1 for v
    digitplaneNum = 1+(moduleNum*4)+subPlaneNum;

    digitTrackID = (*digitHit)->strawMCDigit.strawMCHits.at(0)->trackID;
    digitEventNumInFill = (*digitHit)->strawMCDigit.strawMCHits.at(0)->eventNumInFill;

}

void GeaneDummyUtils::fillTruthParams(int firstPlaneHit, bool zeroPlane){

    std::pair<gm2truth::GhostDetectorArtRecord, int> dummyPlaneReturn; // put these as members

    // Get the right set of dummy hits for this eventNumInFill
    auto& eventNumDummyHits = dummyHits.at(digitEventNumInFill);

    for( int index = 0; index < int(eventNumDummyHits.size()); index++ ) { // brute force loop over dummy plane hits, then associate them to the geane art record

      gm2truth::GhostDetectorArtRecord& dummyHit = eventNumDummyHits.at(index).second;

      if (dummyHit.parentTrackID != 1) continue; // skip if non-primary positron
      if (dummyHit.trackID != digitTrackID) continue; // make sure I'm looking at the same particle

      fillDummyHitInfo(dummyHit);

      if (dummystationNumber != digitstationNumber) continue;

      /////////////////////////////////////////////////////////////////////////////////////
      // code to determine how to check for correct dummy hit depending on what we're filling

      bool fillInfo = false;

      if(!zeroPlane && dummysubPlaneNum == 5) continue; // if filling non-zero plane params, skip 0 planes

      if (!zeroPlane && dummytotalplaneNum == digitplaneNum)
      {
         hitDummyPlaneAlready = (dummyPlanesHit.find(dummytotalplaneNum) != dummyPlanesHit.end());
         dummyPlanesHit.insert(dummytotalplaneNum); // used to determine whether the same dummy plane was hit more than once by the same particle

         fillInfo = true;
      } 
      else if (zeroPlane && dummysubPlaneNum == 5 && dummymoduleNum == int((firstPlaneHit-1)/4.)) // if filling zero plane params make sure we're on the right plane
      {
         fillInfo = true;
      } 

      /////////////////////////////////////////////////////////////////////////////////////

      if (fillInfo)
      {
         dummyPlaneReturn.first = dummyHit;
         dummyPlaneReturn.second = eventNumDummyHits.at(index).first;

        break;
      }
    }

    sortedDummyPlaneHits.emplace_back(dummyPlaneReturn);
}

void GeaneDummyUtils::createStartGuess(std::vector<double> *startingParameters){

    gm2geom::CoordSystem3Vector dHitPos(sortedDummyPlaneHits.at(0).first.position.x(), sortedDummyPlaneHits.at(0).first.position.y(), sortedDummyPlaneHits.at(0).first.position.z() , "world");
    gm2geom::CoordSystem3Vector planeIntersection = dHitPos.transform(cs_, stationStr);

    gm2geom::CoordSystem3Vector dHitMom(sortedDummyPlaneHits.at(0).first.momentum.x(), sortedDummyPlaneHits.at(0).first.momentum.y(), sortedDummyPlaneHits.at(0).first.momentum.z() , "world");
    gm2geom::CoordSystem3Vector momGlob = dHitMom.transform(cs_, stationStr, true);

    double initialYPosChange = 2.*(r3.Rndm()-0.5)*yPosChange_; // mm
    double initialXPosChange = 2.*(r3.Rndm()-0.5)*xPosChange_;

    double initialXMomChange = xMomChange_*2*(r3.Rndm()-0.5); // MeV
    double initialYMomChange = yMomChange_*2*(r3.Rndm()-0.5); 
    double initialZMomChange = zMomChange_*2*(r3.Rndm()-0.5); 


    startingParameters->push_back(planeIntersection.x() + initialXPosChange);
    startingParameters->push_back(planeIntersection.y() + initialYPosChange);
    startingParameters->push_back(planeIntersection.z());

    startingParameters->push_back(1.*momGlob.x() + initialXMomChange);
    startingParameters->push_back(1.*momGlob.y() + initialYMomChange);
    startingParameters->push_back(1.*momGlob.z() + initialZMomChange);

}

void GeaneDummyUtils::fillLRFromTruth(gm2strawtracker::TrackDetailArtRecord & trackDetailArtRecord) {
    TrackArtRecord track(trackDetailArtRecord);
    fillLRFromTruth(track);
} 

void GeaneDummyUtils::fillLRFromTruth(gm2strawtracker::TrackArtRecord & trackArtRecord) {

   trueHitSides.clear();
   trueHitSides.resize(geaneTrackUtils_.maxNumPlanes, gNA_side); // reset the sides whenever this is called

   for (int i = 0; i < int(trackArtRecord.trackPlanesHitList.size()); ++i) // loop over hit planes
   {
    auto planeNum = trackArtRecord.trackPlanesHitList.at(i);
    auto dummyHit = trackArtRecord.dummyPlaneHits.at(i+1); // i+1 since dummy planes start with a 0 plane

    gm2geom::CoordSystem3Vector dHitPos(dummyHit->position.x(), dummyHit->position.y(), dummyHit->position.z() , "world");
    gm2geom::CoordSystem3Vector planeIntersection = dHitPos.transform(cs_, stationStr);

    double planeUposition = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*planeIntersection.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,1)*planeIntersection.y();  
    double planeVposition = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,0)*planeIntersection.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,1)*planeIntersection.y();        

    if(geaneTrackUtils_.isUPlane(planeNum)){
      if (planeUposition >= trackArtRecord.geaneHits.geaneWireUVPositions[3][planeNum]) trueHitSides.at(planeNum) = gLeft;
      else trueHitSides.at(planeNum) = gRight;
    }
    else if (!geaneTrackUtils_.isUPlane(planeNum)){
      if (planeVposition >= trackArtRecord.geaneHits.geaneWireUVPositions[4][planeNum]) trueHitSides.at(planeNum) = gLeft;
      else trueHitSides.at(planeNum) = gRight;
    }                
  }
}

std::vector<int> GeaneDummyUtils::checkLRAgainstTruth(gm2strawtracker::TrackDetailArtRecord & trackArtRecord){
     TrackArtRecord track(trackArtRecord);
     auto wrongHitSides = checkLRAgainstTruth(track);
     return wrongHitSides;
}


std::vector<int> GeaneDummyUtils::checkLRAgainstTruth(gm2strawtracker::TrackArtRecord & trackArtRecord){

    std::vector<int> wrongHitSides(geaneTrackUtils_.maxNumPlanes, 1); // default to correct choice and fix if wrong

    for (int i = 0; i < int(trackArtRecord.geaneHits.geaneHitSides.size()); ++i)
    {
      if(trueHitSides.at(i) != trackArtRecord.geaneHits.geaneHitSides.at(i)) wrongHitSides.at(i) = -1;
      if(trackArtRecord.geaneHits.geaneHitSides.at(i) == gNA_side) wrongHitSides.at(i) = 0; // if no hit set to 0
    }

    return wrongHitSides;
}

  
}//namespace
