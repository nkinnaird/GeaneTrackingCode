// header file
#include "SimpleTrackFindingUtils.hh"

// namespace
using namespace gm2strawtracker;

// standard constructor
SimpleTrackFindingUtils::SimpleTrackFindingUtils() 
   : fitter_()
   , name_("SimpleTrackFindingUtils")
   , cs_()
   , fhicl_()
   , position_()
{}

TrackCandidateArtRecord SimpleTrackFindingUtils::findTrackCandidates( const StrawTimeIslandArtRecord& island )
{
   mf::LogVerbatim(name_) << "Enter SimpleTrackFindingUtils::findTrackCandidates using the model "
                          << fhicl_.get<std::string>("trackFinderName");

   TrackCandidateArtRecord candidate;

   if( fhicl_.get<std::string>("trackFinderName") != "SimpleTrackFinder" ) {
     mf::LogWarning(name_) << "The track candidate cannot be formed.! Please check the input model name!\n";
     return candidate;
   }

   // get the digits
   auto digits = island.strawDigits;

   if( digits.size() > 32) {
     mf::LogWarning(name_) << "Too many digits were included in the island to form a track candidate!\n";
     return candidate;
   }

   // sort the digits
   std::sort(digits.begin(),digits.end(),StrawTrackerSort::StrawDigitArtPtrRowSorter());

   // containers to hold positions and errors (U,V,Z) coordinates
   std::vector< gm2geom::CoordSystem3Vector> positions;
   std::vector< gm2geom::CoordSystem3Vector> errors;

   // initialize
   auto dummyValue = -9999;  //TODO use NAN

   // set the starting position, is overwritten with the initial guess position
   gm2strawtracker::GeaneState geane;

   int iScallop = digits.front()->wireID.getStation();
   int iMod = digits.front()->wireID.getModule();

   WireID wire0(iScallop, iMod, gm2strawtracker::u_view , 0, 16);
   gm2geom::CoordSystem3Vector frontposition = wire0.getCentreInWorld(cs_);
   auto geaneFrame = Form("TrackerStation[%d]", iScallop);
   frontposition = frontposition.transform(cs_,geaneFrame);

   geane.XYZPosition = gm2geom::CoordSystem3Vector(0, 0, frontposition.z()-10, geaneFrame); // 10 is a hardcoded number for now, matching the placement of dummy planes in TrackerDummyPlane_service.cc

   // get positions
   for(auto& digit : digits) {
     auto trackerPosition      = digit.get()->wireID.getCentreInStation(cs_);
     gm2geom::CoordSystem3Vector position(0.,0.,0.,geaneFrame);
    
     if( digit.get()->wireID.getView() == gm2strawtracker::u_view ) {
       position.setX( trackerPosition.z() );
       position.setY( dummyValue );
       position.setZ( trackerPosition.x() );

     } else if( digit.get()->wireID.getView() == gm2strawtracker::v_view ) {
       position.setX( dummyValue );
       position.setY( trackerPosition.z() );
       position.setZ( trackerPosition.x() );
     } 

     positions.push_back( position );
   } // end of getting the position


   // create a track candidate
   candidate.type = TrackPatRecognition::SimplePatRec; 
   std::for_each(digits.begin(),digits.end(),[&](auto& dig){ candidate.strawDigits.push_back(dig); });

   // sanity check
   if( geane.XYZPosition.coordSystemName.find("Tracker") != std::string::npos ) {
   
     // fit the track candidate
     fitter_.makeFittedTrack(positions,errors,geane);

     // store the initial parameters to the track candidate
     candidate.geane = geane;
 
   } // end of sanity condition

   mf::LogVerbatim(name_) << "Exit SimpleTrackFindingUtils::findTrackCandidates\n";
   return candidate;
}
