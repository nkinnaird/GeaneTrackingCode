#include "GeaneLRUtils.hh"

namespace gm2strawtracker {

GeaneLRUtils::GeaneLRUtils(fhicl::ParameterSet const & p) 
  : geaneTrackUtils_()
  , matrixCalculations_(p)
  , lockLowDCAs_(p.get<double>("lockLowDCAs",0))
  , wireCenterError(2.*lockLowDCAs_/sqrt(12))
{}


void GeaneLRUtils::lockSmallDCAErrors(gm2strawtracker::TrackDetailArtRecord & geaneTrack)
{
    for (int i = 0; i < int(geaneTrack.trackPlanesHitList.size()); ++i)
    {
      int planeNum = geaneTrack.trackPlanesHitList.at(i);
      if(geaneTrack.measuredDCAs.at(planeNum) < lockLowDCAs_) geaneTrack.dcaErrors.at(planeNum) = wireCenterError;
    }
}

std::vector<gm2strawtracker::GeaneHitSide> GeaneLRUtils::lockSmallDCACenters(gm2strawtracker::TrackDetailArtRecord & geaneTrack)
{
    std::vector<gm2strawtracker::GeaneHitSide> hitSides = geaneTrack.geaneHits.inputHitSides;

    for (int i = 0; i < int(geaneTrack.trackPlanesHitList.size()); ++i)
    {
      int planeNum = geaneTrack.trackPlanesHitList.at(i);
      if(geaneTrack.measuredDCAs.at(planeNum) < lockLowDCAs_) hitSides.at(planeNum) = gCenter;
    }

    return hitSides;
}


void GeaneLRUtils::fillLRFromFit(gm2strawtracker::TrackDetailArtRecord & track)
{
    std::vector<gm2strawtracker::GeaneHitSide> hitSides(geaneTrackUtils_.maxNumPlanes, gNA_side);
    auto geaneHitsOnTrack = track.geaneHits;

    for (int i = 0; i < int(track.trackPlanesHitList.size()); ++i) // loop over hit planes
    {
       int planeNum = track.trackPlanesHitList.at(i);

       if (geaneTrackUtils_.isUPlane(planeNum))
       {
           if (geaneTrackUtils_.getPredUPos(geaneHitsOnTrack, planeNum) >= track.geaneHits.geaneWireUVPositions[3][planeNum]) hitSides.at(planeNum) = gLeft;
           else hitSides.at(planeNum) = gRight;
       }
       else if (!geaneTrackUtils_.isUPlane(planeNum)){
          if (geaneTrackUtils_.getPredVPos(geaneHitsOnTrack, planeNum) >= track.geaneHits.geaneWireUVPositions[4][planeNum]) hitSides.at(planeNum) = gLeft;
          else hitSides.at(planeNum) = gRight;
      }
    }

    track.geaneHits.inputHitSides = hitSides; // this is needed since lockSmallDCACenters reads in from inputHitSides
    if (lockLowDCAs_) 
    {
       track.geaneHits.inputHitSides = lockSmallDCACenters(track);
    }

    track.geaneHits.geaneHitSides = track.geaneHits.inputHitSides;
}


void GeaneLRUtils::fillLRFromGeom(gm2strawtracker::TrackDetailArtRecord & track)
{
    std::vector<gm2strawtracker::GeaneHitSide> hitSides = track.geaneHits.geaneHitSides; // fill with what we had from previous fit and overwrite with geom

    StrawDigitPtrCollection hitDigits;
    std::for_each(track.strawDCADigits.begin(),track.strawDCADigits.end(),[&](auto dcaDigit){ hitDigits.push_back( dcaDigit->digit ); });

    std::sort(hitDigits.begin(), hitDigits.end(), StrawTrackerSort::StrawDigitArtPtrRowSorter()); // sort the digits by planeNum

    auto geaneHitsOnTrack = track.geaneHits;

    for (auto it = hitDigits.begin(); it != hitDigits.end()-1;  ++it )
    {
      auto& firstHit  = *it;
      auto& secondHit = *(it+1);

      int firstModuleNum   = firstHit->wireID.getModule();
      int firstViewNum     = firstHit->wireID.getView();
      int firstStrawNum    = firstHit->wireID.getWire();
      int firstSubPlaneNum = firstHit->wireID.getLayer() + 2 * firstViewNum;
      int firstPlaneNum    = 1+(firstModuleNum*4)+firstSubPlaneNum;

      int secondModuleNum   = secondHit->wireID.getModule();
      int secondViewNum     = secondHit->wireID.getView();
      int secondStrawNum    = secondHit->wireID.getWire();
      int secondSubPlaneNum = secondHit->wireID.getLayer() + 2 * secondViewNum;
      int secondPlaneNum    = 1+(secondModuleNum*4)+secondSubPlaneNum;

      if(secondPlaneNum <= firstPlaneNum) throw cet::exception("GeaneLRUtils") << "Digits not ordered by planeNum in LR utils. \n";

      // check to make sure the module and view nums are the same
      if(firstModuleNum == secondModuleNum && firstViewNum == secondViewNum)
      {
        if(firstViewNum == 0){

          //check here if angle less than .1 - for which the geom LR is 100% correct
          if( abs(std::atan2(geaneTrackUtils_.getPredUMom(geaneHitsOnTrack,firstPlaneNum), geaneTrackUtils_.getPredZMom(geaneHitsOnTrack,firstPlaneNum))) > 0.1)
          { 
              continue;
          }

          if(secondStrawNum - firstStrawNum >= 0){
             hitSides.at(firstPlaneNum)  = gLeft;
             hitSides.at(secondPlaneNum) = gRight;
          }
          else{
             hitSides.at(firstPlaneNum)  = gRight;
             hitSides.at(secondPlaneNum) = gLeft;
          }
        }
        else {

          //check here if angle less than .1
          if( abs(std::atan2(geaneTrackUtils_.getPredVMom(geaneHitsOnTrack,firstPlaneNum), geaneTrackUtils_.getPredZMom(geaneHitsOnTrack,firstPlaneNum))) > 0.1)
          { 
              continue;
          }

          if(secondStrawNum - firstStrawNum <= 0){
             hitSides.at(firstPlaneNum)  = gRight;
             hitSides.at(secondPlaneNum) = gLeft;
          }
          else{
             hitSides.at(firstPlaneNum)  = gLeft;
             hitSides.at(secondPlaneNum) = gRight;
          }
        }
      }

    }

    // this doesn't deal with locked centers which I think is okay, since geom trumps that
    track.geaneHits.geaneHitSides = hitSides;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void GeaneLRUtils::setUnknownSides(gm2strawtracker::TrackDetailArtRecord & geaneTrack)
{
    for (int i = 0; i < int(geaneTrack.trackPlanesHitList.size()); ++i)
    {
      int planeNum = geaneTrack.trackPlanesHitList.at(i);
      geaneTrack.geaneHits.inputHitSides.at(planeNum) = gUnknown; // fill input hit sides with unknown to start (non hits stay gNA_side)
    }

    if(lockLowDCAs_) geaneTrack.geaneHits.inputHitSides = lockSmallDCACenters(geaneTrack);
}

void GeaneLRUtils::fillLRFromNodes(gm2strawtracker::TrackDetailArtRecord & track)
{
    auto seeds = track.candidate->strawSeeds;

    for (auto& seed : seeds)
    {
       for (auto& node : seed->nodes)
       {
          auto& hit = node.strawDigitToCluster[0].first; 

          int moduleNum   = hit->wireID.getModule();
          int subPlaneNum = hit->wireID.getLayer() + 2 * hit->wireID.getView(); // layer number 0 to 1 in each layer, view number = 0 for u and 1 for v
          int planeNum    = 1+(moduleNum*4)+subPlaneNum;

          if(node.side == StrawSide::left)         track.geaneHits.inputHitSides.at(planeNum) = gLeft;
          else if(node.side == StrawSide::right)   track.geaneHits.inputHitSides.at(planeNum) = gRight;
          else if(node.side == StrawSide::na_side) track.geaneHits.inputHitSides.at(planeNum) = gUnknown;
          // should put another thing here if at some point the nodes.side can be centered

       } // loop over nodes
    } // loop over seeds

    if (lockLowDCAs_) 
    {  
       track.geaneHits.inputHitSides = lockSmallDCACenters(track); // this will overwrite the LR given from the nodes, but that should be okay
    }
}


std::vector<int> GeaneLRUtils::getUVUnkowns(gm2strawtracker::TrackDetailArtRecord & geaneTrack, bool Uplanes)
{
    std::vector<int> unknowns;
    for (int i = 0; i < geaneTrack.trackNumPlanesHit; ++i)
    {
       int planeNum = geaneTrack.trackPlanesHitList.at(i);

       if (Uplanes && geaneTrackUtils_.isUPlane(planeNum)) // uhit
       {
         if(geaneTrack.geaneHits.inputHitSides.at(planeNum) == gUnknown) unknowns.push_back(planeNum);
       }
       else if(!Uplanes && !geaneTrackUtils_.isUPlane(planeNum)) // vhit
       { 
         if(geaneTrack.geaneHits.inputHitSides.at(planeNum) == gUnknown) unknowns.push_back(planeNum);
       }
    }

    return unknowns;
}


std::vector<std::pair<double, int> > GeaneLRUtils::getTopSequences(gm2strawtracker::TrackDetailArtRecord & geaneTrack, bool Uplanes, int unkownsSize)
{
    std::vector<std::pair<double, int> > topUVSequences;

    // Call this method to perform actions that only need to be done once for all combinations.
    matrixCalculations_.makeHybridErrMat(geaneTrack,Uplanes); 

    //need to generalize this for varying errors per plane
    unsigned long sequence = 0;
    for (int i = 0; i < pow(2, unkownsSize); ++i)
    {
      sequence = i; 
      boost::dynamic_bitset<> mySequence(unkownsSize, sequence);

      /////////////////////////////////////////////////////////////////////////////////////
      std::vector<gm2strawtracker::GeaneHitSide> UVsides(geaneTrackUtils_.maxNumPlanes, gNA_side); // will be filled into GEANEHits

      int nUVHit = 0;
      for (int j = 0; j < geaneTrack.trackNumPlanesHit; ++j)
      {
          int planeNum = geaneTrack.trackPlanesHitList.at(j);

          if(Uplanes && !geaneTrackUtils_.isUPlane(planeNum)) UVsides.at(planeNum)      = gCenter; // If performing U sequences, set Vs to the center and vice versa. (Could lock these as well possibly...)
          else if(!Uplanes && geaneTrackUtils_.isUPlane(planeNum)) UVsides.at(planeNum) = gCenter;
          else if(geaneTrack.geaneHits.inputHitSides.at(planeNum) == gUnknown) { // fills unknowns with sequence choices
            if (mySequence[nUVHit] == 0) UVsides.at(planeNum)      = gLeft;
            else if (mySequence[nUVHit] == 1) UVsides.at(planeNum) = gRight;
            nUVHit++;
          }
          else UVsides.at(planeNum) = geaneTrack.geaneHits.inputHitSides.at(planeNum); // keeps locked ones locked
      }
      geaneTrack.geaneHits.geaneHitSides = UVsides;

      modifyMeasuredParams(geaneTrack); // call separate method for adding/subtracting dca to wire postions

      // calcMeasuredParams(track); // this is slower (+10 secs over ~160 tracks) and something is wrong with it, though it should give somewhat better results than using modifyMeasuredParams when working

      double sequenceChiSquared = matrixCalculations_.sequenceChecking(geaneTrack);

      if (topUVSequences.size() < 10)
      {
         topUVSequences.emplace_back(sequenceChiSquared, sequence);
         std::sort(topUVSequences.begin(), topUVSequences.end());
      }
      else if (sequenceChiSquared < topUVSequences.back().first)
      {
         topUVSequences.pop_back();
         topUVSequences.emplace_back(sequenceChiSquared, sequence);
         std::sort(topUVSequences.begin(), topUVSequences.end());
      }
    }

    return topUVSequences;
}


void GeaneLRUtils::modifyMeasuredParams(gm2strawtracker::TrackDetailArtRecord & geaneTrack)
{
    for (int i = 0; i < geaneTrack.trackNumPlanesHit; ++i)
    {
      // mySequence[planeNum] == 0 or 1 then + or - , Left or Right

      int planeNum = geaneTrack.trackPlanesHitList.at(i);
      double sign  = 20; // some random number for initialization

      if(geaneTrack.geaneHits.geaneHitSides.at(planeNum)      == gNA_side) throw cet::exception("GeaneLRUtils") << "Track side is na which shouldn't happen in modifyMeasuredParams. planeNum: " << planeNum << " \n";
      else if(geaneTrack.geaneHits.geaneHitSides.at(planeNum) == gUnknown) throw cet::exception("GeaneLRUtils") << "Track side is unknown which shouldn't happen in modifyMeasuredParams. planeNum: " << planeNum << " \n";
      else if(geaneTrack.geaneHits.geaneHitSides.at(planeNum) == gCenter)  sign =  0;
      else if(geaneTrack.geaneHits.geaneHitSides.at(planeNum) == gLeft)    sign = +1;
      else if(geaneTrack.geaneHits.geaneHitSides.at(planeNum) == gRight)   sign = -1;

      geaneTrack.geaneHits.geaneMeasuredParameters[3][planeNum] = geaneTrack.geaneHits.geaneWireUVPositions[3][planeNum] + sign*geaneTrack.measuredDCAs.at(planeNum);
      geaneTrack.geaneHits.geaneMeasuredParameters[4][planeNum] = geaneTrack.geaneHits.geaneWireUVPositions[4][planeNum] + sign*geaneTrack.measuredDCAs.at(planeNum);
      // nonhit U or V parameter won't be used in fitting math
    }

    return;
}


void GeaneLRUtils::getSequenceSides(std::vector<int> UVunknownsInThisEvent, int sequenceInt, std::vector<gm2strawtracker::GeaneHitSide> *sides)
{
    boost::dynamic_bitset<> mySequence(int(UVunknownsInThisEvent.size()), sequenceInt);
    for (int i = 0; i < int(UVunknownsInThisEvent.size()); ++i)
    {
        int planeNum = UVunknownsInThisEvent.at(i);
        if (mySequence[i] == 0) sides->at(planeNum)      = gLeft;
        else if (mySequence[i] == 1) sides->at(planeNum) = gRight;
    }
}


}//namespace
