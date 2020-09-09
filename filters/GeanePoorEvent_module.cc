////////////////////////////////////////////////////////////////////////
//
// Filter out reconstructed tracks with pulls greater than some cutoff.
// 
// 
//
// Nick Kinnaird 
//
////////////////////////////////////////////////////////////////////////

#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"


// art includes
#include "art/Framework/Core/EDFilter.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4PhysicalConstants.hh"

// data product includes
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"
#include "gm2dataproducts/mc/ghostdetectors/GhostDetectorArtRecord.hh"

// common util
#include "gm2util/common/dataModuleDefs.hh"

#include <Eigen/Dense>
#include "gm2tracker/utils/GeaneEigenStorageUtils.hh"

#include "gm2tracker/utils/GeaneTrackUtils.hh"


namespace gm2strawtracker {
  class GeanePoorEvent;
}

class gm2strawtracker::GeanePoorEvent : public art::EDFilter {

   public:
     explicit GeanePoorEvent(fhicl::ParameterSet const & p);

     bool filter(art::Event & e) override;
     bool beginRun(art::Run & r) override;


   private:

     //Input and output data collection labels
     std::string TrackModuleLabel_;
     std::string TrackInstanceName_;

     std::string poorGEANECollectionName_;

     //Set the poor event metric at any pull being > or < some threshold
     double pullThreshold_;

     //Helper tools and data container
     gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;
     gm2geom::CoordSystemsStoreData cs_;

};


gm2strawtracker::GeanePoorEvent::GeanePoorEvent(fhicl::ParameterSet const & pset) 
    : TrackModuleLabel_( pset.get<std::string>("GEANEModuleLabel", "wronggeanemodulelabel") )
    , TrackInstanceName_( pset.get<std::string>("GEANEInstanceName", "wronggeaneinstancename") )
    , poorGEANECollectionName_( pset.get<std::string>("poorGEANECollectionName", "poorGeaneArtRecords") )
    , pullThreshold_( pset.get<double>("pullThreshold", 5.0))
    , geaneTrackUtils_()
    , cs_()
{
    produces<gm2strawtracker::TrackArtRecordCollection>( poorGEANECollectionName_ );
}


bool gm2strawtracker::GeanePoorEvent::beginRun(art::Run & r)
{
   //Get coord systems
   cs_ = artg4::dataFromRunOrService<gm2geom::CoordSystemsStoreData, gm2geom::CoordSystemsStore>
         ( r, dataModuleDefs::coordSysModuleLabel(),dataModuleDefs::coordSysInstanceLabel() );

   if( cs_.size() == 0 ) {
      mf::LogWarning("GeanePoorEvent") << "This run does not contain any data associated with the coordinate system\n";
      return false;
   }

   return true;
}


bool gm2strawtracker::GeanePoorEvent::filter(art::Event & e) 
{
   mf::LogInfo info("GeanePoorEvent");
   info << "Enter GeanePoorEvent\n";

   // Get the reconstructed tracks
   art::Handle<gm2strawtracker::TrackArtRecordCollection> TrackDataHandle;
   bool foundTrackCollection = e.getByLabel(TrackModuleLabel_,TrackInstanceName_,TrackDataHandle);
   if( !foundTrackCollection ) {
     throw cet::exception("GeanePoorEvent") << "No Track Collection in this event (\"" << TrackModuleLabel_ << "\":\"" << TrackInstanceName_ << "\")\n";
   }

   // make a collection of poor tracks and GEANE hits 
   std::unique_ptr<gm2strawtracker::TrackArtRecordCollection> poorGEANEHits(new gm2strawtracker::TrackArtRecordCollection); 

   // Initialize
   int trackNum = 0;

   auto tracks = *TrackDataHandle;

   // Loop over tracks
   for( auto & track : tracks ) {
      trackNum++;

      // Apply Cuts
      if (track.failureMode != 0 ) 
      {
        continue;
      }
 
      // Get eigen object conversion
      auto out_track         = GeaneEigenStorageUtils::ReadEigenFromStorage(track);
      auto geaneHitsOnTrack  = out_track.geaneHits; 
 
      auto dummyHit = track.dummyPlaneHits.at(0);

      int stationNumber = track.candidate->strawDigits.at(0)->wireID.getStation();

      stringstream stationStream;
      stationStream << "TrackerStation[" << stationNumber << "]";
      string stationStr = stationStream.str();

      gm2geom::CoordSystem3Vector dHitPos(dummyHit->position.x(), dummyHit->position.y(), dummyHit->position.z() , "world");
      gm2geom::CoordSystem3Vector plane0Postion = dHitPos.transform(cs_, stationStr);

      double plane0Uposition = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*plane0Postion.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,1)*plane0Postion.y();
      double plane0Vposition = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,0)*plane0Postion.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,1)*plane0Postion.y();

      gm2geom::CoordSystem3Vector dHitMom(dummyHit->momentum.x(), dummyHit->momentum.y(), dummyHit->momentum.z() , "world");
      gm2geom::CoordSystem3Vector plane0Momentum = dHitMom.transform(cs_, stationStr, true);

      double plane0Umomentum = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*plane0Momentum.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,1)*plane0Momentum.y();
      double plane0Vmomentum = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,0)*plane0Momentum.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,1)*plane0Momentum.y();
                        
      double plane0Zmomentum = plane0Momentum.z();

      double plane0Totalmomentum = plane0Momentum.mag();


      double predicted0Momentum = geaneTrackUtils_.getPredMom(geaneHitsOnTrack, 0);

      double predicted0XMomentum = geaneTrackUtils_.getPredXMom(geaneHitsOnTrack, 0);
      double predicted0YMomentum = geaneTrackUtils_.getPredYMom(geaneHitsOnTrack, 0);
      double predicted0ZMomentum = geaneTrackUtils_.getPredZMom(geaneHitsOnTrack, 0);

      double predicted0UMomentum = geaneTrackUtils_.getPredUMom(geaneHitsOnTrack, 0);
      double predicted0VMomentum = geaneTrackUtils_.getPredVMom(geaneHitsOnTrack, 0);

      double predicted0YPosition = geaneTrackUtils_.getPredYPos(geaneHitsOnTrack, 0);
      double predicted0ZPosition = geaneTrackUtils_.getPredZPos(geaneHitsOnTrack, 0);

      double predicted0UPosition = geaneTrackUtils_.getPredUPos(geaneHitsOnTrack, 0);
      double predicted0VPosition = geaneTrackUtils_.getPredVPos(geaneHitsOnTrack, 0);


      Eigen::MatrixXd true0Covariance = geaneHitsOnTrack.covarianceTotalInverse.inverse();

      double pull1oP  = (1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0)));
      double pullPuPx = (predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1));
      double pullPvPx = (predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2));

      double pullU = (predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3)));
      double pullV = (predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4)));

      //check if pull values are above or below pull threshold
      if( pull1oP  <= -pullThreshold_ || pull1oP  >= pullThreshold_  || 
          pullPuPx <= -pullThreshold_ || pullPuPx >= pullThreshold_  ||
          pullPvPx <= -pullThreshold_ || pullPvPx >= pullThreshold_  ||
          pullU    <= -pullThreshold_ || pullU    >= pullThreshold_  ||
          pullV    <= -pullThreshold_ || pullV    >= pullThreshold_  )
      {
          info << "\n BadPull event: " << e.event() << " run: " << e.run() << " subrun: " << e.subRun() << " trackNum: " << trackNum << "\n";
          poorGEANEHits->push_back(track);
      } 

   } // end loop over tracks


   if(poorGEANEHits->size() > 0)
   {
     GeaneEigenStorageUtils::PrepareEigenForStorage(poorGEANEHits); 
     e.put(std::move(poorGEANEHits), poorGEANECollectionName_);

     return true; // still filter on whether a track within the event was poor or not, and make a separate collection of those poor events
   } 
   else
   {
     e.put(std::move(poorGEANEHits), poorGEANECollectionName_);
     return false;
   } 

}

DEFINE_ART_MODULE(gm2strawtracker::GeanePoorEvent)
