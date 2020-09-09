//!----------------------------------
//!  This class provides useful functions for dealing with tracker dummy plane hits along with the GEANEArtRecord for a track.
//!----------------------------------

#ifndef GEANEDUMMYUTILS_HH
#define GEANEDUMMYUTILS_HH

#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"
#include "gm2dataproducts/mc/ghostdetectors/GhostDetectorArtRecord.hh"

#include "gm2geom/coordSystems/CoordSystem.hh"
#include "gm2geom/coordSystems/CoordSystemsStoreData.hh"
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh"  

#include "TRandom3.h"

#include <map>

#include "gm2tracker/utils/GeaneTrackUtils.hh"


namespace gm2strawtracker {

  class GeaneDummyUtils {
     
   public :

     GeaneDummyUtils(fhicl::ParameterSet const & p);

     void fillDummyHits(art::Handle<gm2truth::GhostDetectorArtRecordCollection> DummyDataHandle);
     void clearDummyHolders();

     void getTrackDummyHits(art::Handle<gm2truth::GhostDetectorArtRecordCollection> DummyDataHandle, gm2strawtracker::TrackDetailArtRecord & trackArtRecord);

     void fillCS(gm2geom::CoordSystemsStoreData cs) { cs_ = cs; };
     void fillStationStr(std::string stationString, int digStationNumber) { stationStr = stationString; 
                                                                    digitstationNumber = digStationNumber; };

     void fillDummyHitInfo(gm2truth::GhostDetectorArtRecord dummyHit);
     void fillDigitHitInfo(const art::Ptr<gm2strawtracker::StrawDigitArtRecord> *digitHit);

     void fillTruthParams(int firstPlaneHit, bool zeroPlane);

     void createStartGuess(std::vector<double> *startingParameters);

     void fillLRFromTruth(gm2strawtracker::TrackArtRecord & trackArtRecord); // fill trueHitSides
     void fillLRFromTruth(gm2strawtracker::TrackDetailArtRecord & trackDetailArtRecord); // fill trueHitSides

     std::vector<int> checkLRAgainstTruth(gm2strawtracker::TrackArtRecord & trackArtRecord); // check sides against truth
     std::vector<int> checkLRAgainstTruth(gm2strawtracker::TrackDetailArtRecord & trackDetailArtRecord); // check sides against truth

     bool hitDummyPlaneAlready = false;
     std::vector<gm2strawtracker::GeaneHitSide> trueHitSides; // truth sides

   private :

    std::map<int,std::vector<std::pair<unsigned int, gm2truth::GhostDetectorArtRecord> > > dummyHits;  // eventNumInFill, pair of index in main art collection + GhostDetector hit for this fill builder input

    std::vector< std::pair<gm2truth::GhostDetectorArtRecord, int> > sortedDummyPlaneHits; // pairs to record dummy plane art records and their associated indices in the collection that's passed in

    gm2geom::CoordSystemsStoreData cs_;

    string stationStr;

    int digitstationNumber;
    int digitplaneNum;
    int digitTrackID;
    int digitEventNumInFill;

    int dummystationNumber;
    int dummymoduleNum;
    int dummysubPlaneNum;
    int dummytotalplaneNum;

    std::set<int> dummyPlanesHit; // For checking to make sure there is only 1 dummy plane hit per plane, otherwise cancel the fitting.

    TRandom3 r3;
    int rseed_;

    double xPosChange_;
    double yPosChange_;

    double xMomChange_; // Variables for simulating the starting momentum guess. The truth momentum is taken and then gaussian smeared with these variables.
    double yMomChange_;
    double zMomChange_;

    gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;

  };

}

#endif
