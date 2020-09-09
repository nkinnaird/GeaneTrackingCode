////////////////////////////////////////////////////////////////////////
// Class:       GeaneReco
// Module Type: producer
// File:        GeaneReco_module.cc
//
//
// Module for tracking using the error_propagation (GEANE) routines in 
// the Geant4 library, with many pieces built and added in - accesses several utils files.
// 
// Created by Nick Kinnaird. nickkinn@bu.edu
//
// Update from T. Walton -- no changes in track fitting
// update track fitting to work with changes in the modeling of 
// the data products
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "gm2dataproducts/strawtracker/StrawDigitArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackCandidateArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"
#include "gm2dataproducts/mc/ghostdetectors/GhostDetectorArtRecord.hh"

#include "gm2geom/strawtracker/WireID.hh"

#include "artg4/util/DataFromRunOrService.hh"

#include "gm2geom/coordSystems/CoordSystem.hh"
#include "gm2geom/coordSystems/CoordSystemsStoreData.hh"
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh"  
#include "gm2util/coordSystems/CoordSystemUtils.hh"

#include "gm2util/common/dataModuleDefs.hh"
#include "artg4/util/util.hh"

#include <memory>

#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"

#include "gm2tracker/utils/GeaneFittingUtils.hh"
#include "gm2tracker/utils/GeaneEigenStorageUtils.hh"
#include "gm2tracker/quality/TrackQuality_service.hh"

/////////////////////////////////////////////////////////////////////////////////////

// namespace
using namespace gm2strawtracker;

namespace gm2strawtracker {
  class GeaneReco;
}


// class
class gm2strawtracker::GeaneReco : public art::EDProducer {

public:
  explicit GeaneReco(fhicl::ParameterSet const & p);

  GeaneReco(GeaneReco const &) = delete;
  GeaneReco(GeaneReco &&) = delete;

  GeaneReco & operator = (GeaneReco const &) = delete;
  GeaneReco & operator = (GeaneReco &&) = delete;

  void produce(art::Event & e) override;
  void beginRun(art::Run & r ) override;

private:

  // data products module labels
  std::string candidateModuleLabel_;
  std::string candidateInstanceName_;

  std::string dcaDigitModuleLabel_;
  std::string dcaDigitInstanceLabel_;

  std::string csModuleLabel_;
  std::string csInstanceName_;

  std::string instanceLabel_;

  bool keepTrackDetailArtRecord_;
  bool keepTrackArtRecord_;

  gm2geom::CoordSystemsStoreData cs_;

  std::string name_;

  gm2geom::StrawTrackerGeometry sgeom_;

  gm2strawtracker::GeaneFittingUtils  geaneFittingUtils_;

  int numPassesWireFit_; // Number of passes to run for the wire fit - my initial thought is that 2 is fine.
  int numPassesMainFit_;
  
  gm2util::CoordSystemUtils csUtils_;

  string fitMode_;

  std::map< std::string, gm2geom::CoordSystemsStoreData > detCoordMap_;

  art::ServiceHandle<gm2strawtracker::TrackQuality> trackQuality_;
  bool onlyTrackGoodCandidates_;

  // Helper function
  int trackFitting(art::Event & e, 
                   gm2strawtracker::TrackDetailArtRecord & trackFitDetails, 
                   bool firstTrackInEvent);
};

// constructor
gm2strawtracker::GeaneReco::GeaneReco(fhicl::ParameterSet const & p)
  : candidateModuleLabel_( p.get<std::string>("candidateModuleLabel", dataModuleDefs::t0FinderModuleLabel()) )
  , candidateInstanceName_( p.get<std::string>("candidateInstanceName",dataModuleDefs::recoInstanceLabel()) )
  , dcaDigitModuleLabel_( p.get<std::string>("dcaDigitModuleLabel",dataModuleDefs::dcaDigitModuleLabel()) )
  , dcaDigitInstanceLabel_( p.get<std::string>("dcaDigitInstanceLabel",dataModuleDefs::recoInstanceLabel()) )
  , csModuleLabel_( p.get<std::string>("csModuleLabel",dataModuleDefs::coordSysModuleLabel()) )
  , csInstanceName_( p.get<std::string>("csInstanceName",dataModuleDefs::coordSysInstanceLabel()) )
  , instanceLabel_( p.get<std::string>("instanceLabel", dataModuleDefs::recoInstanceLabel()) )
  , keepTrackDetailArtRecord_( p.get<bool>("keepTrackDetailArtRecord", false) )
  , keepTrackArtRecord_( p.get<bool>("keepTrackArtRecord", true) )
  , cs_()
  , name_( "GeaneReco" )
  , sgeom_()
  , geaneFittingUtils_(p)
  , numPassesWireFit_(p.get<int>("numPassesWireFit", 2))
  , numPassesMainFit_(p.get<int>("numPassesMainFit", 2))
  , csUtils_()
  , fitMode_(p.get<string>("fitMode","badFitMode"))
  , detCoordMap_()
  , trackQuality_()
  , onlyTrackGoodCandidates_(p.get<bool>("onlyTrackGoodCandidates",false))
{
  if (keepTrackArtRecord_) {
     produces<gm2strawtracker::TrackArtRecordCollection>( instanceLabel_ );
  }

  if (keepTrackDetailArtRecord_) {
     produces<gm2strawtracker::TrackDetailArtRecordCollection>( instanceLabel_ );
  }
}


// begin of run
void gm2strawtracker::GeaneReco::beginRun(art::Run & r)
{
  // get the coordinate system data
  cs_ = artg4::dataFromRunOrService<gm2geom::CoordSystemsStoreData, gm2geom::CoordSystemsStore>(r,csModuleLabel_,csInstanceName_);
  if( cs_.size() == 0 ) { 
    throw cet::exception(name_) << "This run does not contain any data associated with the coordinate system\n";
    return;
  }

  std::vector<std::string> detNames;
  detNames.push_back("TrackerStation");
  for(auto s : sgeom_.whichScallopLocations) {
    for(unsigned int m = 0; m < sgeom_.getNumModulesPerStation(); ++m) {
       detNames.push_back(Form("Module%d:%d", s, m));
    }
  }

  csUtils_.setDetectorNames(detNames);
  csUtils_.detNameToCoordMap(cs_,detCoordMap_);
  sgeom_.setDetNameToCoords(detCoordMap_);

  // storing vector of z positions in beginRun for speed
  std::vector<double> geaneTrackerZPositions_;

  geaneTrackerZPositions_.push_back(-1); // add entry for plane 0 to match calls to planeNum
  for(unsigned int mod = 0; mod < sgeom_.getNumModulesPerStation(); mod++) {
    for(unsigned int view = 0; view < sgeom_.getNumViewsPerModule(); view++) {
      for(unsigned int layer = 0; layer < sgeom_.getNumLayersPerView(); layer++) {
        WireID straw(18, mod, StrawView(view), layer, 0); // tracker station 18, wire 0
        double truthPlaneTargetZ = straw.getCentreInWorld(cs_).transform(cs_, "TrackerStation[18]").z(); // Z points should be the same for all 3 trackers until alignment is included
        geaneTrackerZPositions_.push_back(truthPlaneTargetZ);
      }
    }
  }
  geaneFittingUtils_.geaneParamUtils_.setZPositions(geaneTrackerZPositions_);

  // pass cs to dummy utils for coord sys transforms  
  geaneFittingUtils_.geaneParamUtils_.dummyUtils_.fillCS(detCoordMap_.find("TrackerStation")->second); 
  geaneFittingUtils_.geaneParamUtils_.setCoordMap(detCoordMap_);
}

// produce geane tracks
void gm2strawtracker::GeaneReco::produce(art::Event & e)
{
  mf::LogInfo info(name_);

  // get the data
  art::Handle<gm2strawtracker::TrackCandidateArtRecordCollection> candidateDataHandle;
  bool success = e.getByLabel(candidateModuleLabel_,candidateInstanceName_,candidateDataHandle);
  if( !success ) {
    throw cet::exception(name_) << "Event " << e.id() << " does not contain any trackCandidate data \"" << candidateModuleLabel_ << ":" << candidateInstanceName_  <<"\"\n";
  }
  auto trackCandidates = *candidateDataHandle;

  // create track-fitting based data products
  std::unique_ptr<gm2strawtracker::TrackArtRecordCollection> trackPtrs(new gm2strawtracker::TrackArtRecordCollection);
  std::unique_ptr<gm2strawtracker::TrackDetailArtRecordCollection> trackDetailPtrs(new gm2strawtracker::TrackDetailArtRecordCollection);

  // sanity check
  if (trackCandidates.size() == 0) {
    if( keepTrackDetailArtRecord_ ) {
      GeaneEigenStorageUtils::PrepareEigenForStorage(trackDetailPtrs);
      e.put(std::move(trackDetailPtrs), instanceLabel_);
    }
  
    if( keepTrackArtRecord_ ) {
      GeaneEigenStorageUtils::PrepareEigenForStorage(trackPtrs);
      e.put(std::move(trackPtrs), instanceLabel_);
    }
     
    return;
  }

  // get coordinate system
  auto css = detCoordMap_.find("TrackerStation")->second;

  // loop over track candidates
  for (int nCand = 0; nCand < int(trackCandidates.size()); ++nCand)
  {
     info << "nCand: " << nCand << "\n";

     gm2strawtracker::TrackDetailArtRecord trackFitDetails; // create object for fitted track
     art::Ptr< gm2strawtracker::TrackCandidateArtRecord > candidatePtr(candidateDataHandle, nCand);

     //check quality
     if (onlyTrackGoodCandidates_ && !trackQuality_->goodCandidate(candidatePtr)){
       info << "candidate failed candidate quality, and onlyTrackGoodCandidates is on, continue\n";
       continue;
     }
     
     /////////////////////////////////////////////////////////////////////////////////////
     // For first iteration of tracking, use data objects defined by the track candidate 
     /////////////////////////////////////////////////////////////////////////////////////

     trackFitDetails.candidate     = candidatePtr; // attach candidate to geane art record at the start
     trackFitDetails.strawClusters = trackCandidates.at(nCand).strawClusters;
     trackFitDetails.island        = trackCandidates.at(nCand).island;
     trackFitDetails.time          = trackCandidates.at(nCand).t0;


     // get the straw DCA digits that are associated with the straw digits on the track candidate
     // a strawDigit may have more than one association
     art::FindManyP<StrawDCADigitArtRecord, int> dcaDigitAssns(trackCandidates.at(nCand).strawDigits, e, art::InputTag(dcaDigitModuleLabel_,dcaDigitInstanceLabel_));

     // container to store straw dca digits on the track candidate
     StrawDCADigitPtrCollection trackCandidateDCADigits;
     
     // loop over straw digits and get the straw dca digit
     for(unsigned int i = 0; i < trackCandidates.at(nCand).strawDigits.size(); ++i) {
       StrawDCADigitPtrCollection dcaDigits;
       std::vector<int const*> trackCandIds;
       dcaDigitAssns.get(i,dcaDigits,trackCandIds);

       auto dcaDigitIdItr = std::find_if(trackCandIds.begin(),trackCandIds.end(),[&](auto id){ return *id == trackCandidates.at(nCand).id; });
       if( dcaDigitIdItr != trackCandIds.end() ) {
         trackCandidateDCADigits.push_back( dcaDigits.at( std::distance(trackCandIds.begin(),dcaDigitIdItr) ) );
       }
     } 

     trackFitDetails.strawDCADigits = trackCandidateDCADigits;

     bool firstTrackInEvent = (nCand == 0) ? true : false;
     int fittingReturn = trackFitting(e, trackFitDetails, firstTrackInEvent); // fit the track

     /////////////////////////////////////////////////////////////////////////////////////

     if (fittingReturn != 0)
     {
        info << "Fitting failed in some way: " << fittingReturn << "\n";

        // if the track failed just make a blank geane art record with the failure mode and 
        // candidate and pas that back - at some point might want to make separate collection, and save some track fitting details
        gm2strawtracker::TrackDetailArtRecord fGEANETrack; 

        fGEANETrack.failureMode    = fittingReturn;
        fGEANETrack.candidate      = candidatePtr;
        fGEANETrack.strawClusters  = trackCandidates.at(nCand).strawClusters;
        fGEANETrack.strawDCADigits = trackCandidateDCADigits;
        fGEANETrack.island         = trackCandidates.at(nCand).island;
        fGEANETrack.time           = trackCandidates.at(nCand).t0;

        if (keepTrackDetailArtRecord_) {
          trackDetailPtrs->push_back(fGEANETrack);
        }

        if (keepTrackArtRecord_) {
           TrackArtRecord fTrack(fGEANETrack);
           trackPtrs->push_back(fTrack);
        }

        continue;
     }

     info << "Chi^2 from fitting: "  << trackFitDetails.chi2 << "\n";
     info << "pValue from fitting: " << trackFitDetails.pValue << "\n";

     // get station number for track
     int stationNumber = trackCandidates[nCand].strawDigits.at(0)->wireID.getStation();
     stringstream stationStream;
     stationStream << "TrackerStation[" << stationNumber << "]";
     std::string stationStr = stationStream.str();

     // loop through trackFitDetails objects and fill a
     // container to hold the track states
     gm2strawtracker::TrackDetailStates trackStates;

     for (int i = 0; i < int(trackFitDetails.trackPlanesHitList.size()); ++i)
     {
       int planeNum = trackFitDetails.trackPlanesHitList.at(i);

       // this world here is for the geanetracker world, which is different from the ring world
       gm2geom::CoordSystem3Vector statePositionGEANEWorld = gm2geom::CoordSystem3Vector(trackFitDetails.geaneHits.geanePredictedParameters[3].at(planeNum), 
                                                                                         trackFitDetails.geaneHits.geanePredictedParameters[4].at(planeNum), 
                                                                                         trackFitDetails.geaneHits.planeZPositions.at(planeNum), stationStr); 

       auto totalP = 1./trackFitDetails.geaneHits.geanePredictedParameters[0].at(planeNum);
       auto statePz = totalP/sqrt(1 + trackFitDetails.geaneHits.geanePredictedParameters[1].at(planeNum)*trackFitDetails.geaneHits.geanePredictedParameters[1].at(planeNum) 
                                    + trackFitDetails.geaneHits.geanePredictedParameters[2].at(planeNum)*trackFitDetails.geaneHits.geanePredictedParameters[2].at(planeNum));
       auto statePy = trackFitDetails.geaneHits.geanePredictedParameters[2].at(planeNum) * statePz;
       auto statePx = trackFitDetails.geaneHits.geanePredictedParameters[1].at(planeNum) * statePz;

       gm2geom::CoordSystem3Vector stateMomentumGEANEWorld = gm2geom::CoordSystem3Vector(statePx, statePy, statePz, stationStr);
      
       trackStates.positionVect.push_back( statePositionGEANEWorld.transform(css, "world") );
       trackStates.momentumVect.push_back( stateMomentumGEANEWorld.transform(css, "world", true) );
     }


     // get the  state position and momentum
     gm2geom::CoordSystem3Vector startStatePositionGEANEWorld = gm2geom::CoordSystem3Vector(trackFitDetails.geaneHits.startingGeaneParameters.at(0), 
                                                                                            trackFitDetails.geaneHits.startingGeaneParameters.at(1), 
                                                                                            trackFitDetails.geaneHits.startingGeaneParameters.at(2), stationStr);
     auto startPosition = startStatePositionGEANEWorld.transform(css,"world");

     gm2geom::CoordSystem3Vector startStateMomentumGEANEWorld = gm2geom::CoordSystem3Vector(trackFitDetails.geaneHits.startingGeaneParameters.at(3), 
                                                                                            trackFitDetails.geaneHits.startingGeaneParameters.at(4), 
                                                                                            trackFitDetails.geaneHits.startingGeaneParameters.at(5), stationStr);
     auto startMomentum = startStateMomentumGEANEWorld.transform(css, "world", true);

     trackStates.positionVect.insert(trackStates.positionVect.begin(),startPosition);
     trackStates.momentumVect.insert(trackStates.momentumVect.begin(),startMomentum);

     trackFitDetails.states.nstates = int(trackFitDetails.trackPlanesHitList.size()) + 1;
     trackFitDetails.states         = trackStates;
     trackFitDetails.momentum       = trackStates.momentumVect.front();
     trackFitDetails.position       = trackStates.positionVect.front();
     trackFitDetails.fitType        = static_cast<gm2strawtracker::FitType>(3);

     // add to collection
     if (keepTrackDetailArtRecord_) {
        trackDetailPtrs->push_back(trackFitDetails); 
     }

     if (keepTrackArtRecord_) {
        TrackArtRecord trackFit(trackFitDetails);
        trackPtrs->push_back(trackFit);
     }


     info << "\n End of candidate - " << fitMode_ << " successful \n ";
   } // end loop over candidates

   // add geane art records to events
   if (keepTrackDetailArtRecord_) {
      GeaneEigenStorageUtils::PrepareEigenForStorage(trackDetailPtrs);
      e.put(std::move(trackDetailPtrs), instanceLabel_);
   } 

   if (keepTrackArtRecord_) {
      GeaneEigenStorageUtils::PrepareEigenForStorage(trackPtrs);
      e.put(std::move(trackPtrs), instanceLabel_);
   }

   return;
}

// track fitting
int gm2strawtracker::GeaneReco::trackFitting(art::Event & e, 
                                             gm2strawtracker::TrackDetailArtRecord & trackFitDetails, 
                                             bool firstTrackInEvent)
{
   // Setup initial parameters before fitting
   int setupReturn = geaneFittingUtils_.geaneParamUtils_.setupParams(e, trackFitDetails, firstTrackInEvent);
   if(setupReturn != 0) return setupReturn;

   /////////////////////////////////////////////////////////////////////////////////////
   // Fit based on specified fit mode.
   /////////////////////////////////////////////////////////////////////////////////////

   int fitReturnInt = -1;

   if(fitMode_ == "truthLRFit")      fitReturnInt = geaneFittingUtils_.truthLRFit(trackFitDetails);
   else if(fitMode_ == "wireFit")    fitReturnInt = geaneFittingUtils_.wireFit(trackFitDetails);
   else if(fitMode_ == "mainFit")    fitReturnInt = geaneFittingUtils_.mainFit(trackFitDetails);
   else if(fitMode_ == "fullSeqFit") fitReturnInt = geaneFittingUtils_.fullSeqFit(trackFitDetails);
   else throw cet::exception(name_) << "Bad fit mode: " << fitMode_ << "\n";

   if(fitReturnInt != 0) return fitReturnInt;
   return geaneFittingUtils_.checkExtraneousFailureModes(trackFitDetails);
}

   
DEFINE_ART_MODULE(gm2strawtracker::GeaneReco)
