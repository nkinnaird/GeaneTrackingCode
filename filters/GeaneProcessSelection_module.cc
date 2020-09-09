////////////////////////////////////////////////////////////////////////
//
// Filter to separate out GEANEArtRecords based on mc truth processes.
// 
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
#include "gm2dataproducts/mc/actions/track/TrajectoryArtRecord.hh"

#include "artg4/util/DataFromRunOrService.hh"
#include "artg4/pluginActions/physicalVolumeStore/PhysicalVolumeStoreData.hh"
#include "artg4/pluginActions/physicalVolumeStore/physicalVolumeStore_service.hh"

// common util
#include "gm2util/common/dataModuleDefs.hh"

#include <Eigen/Dense>

namespace gm2strawtracker {
  class GeaneProcessSelection;
}

class gm2strawtracker::GeaneProcessSelection : public art::EDFilter {

   public:
     explicit GeaneProcessSelection(fhicl::ParameterSet const & p);
     virtual ~GeaneProcessSelection();

     bool filter(art::Event & e) override;
     bool beginRun(art::Run & r) override;


   private:

     // fhicl parameters
     std::string TrackModuleLabel_;
     std::string TrackInstanceName_;

     std::string trajectoryModuleLabel_;
     std::string trajectoryInstanceName_;

     std::string pvsModuleLabel_;
     std::string pvsInstanceName_;

     bool filterOnIonizations_;
     bool filterOnBrems_;
     bool filterOnTransportation_;
     bool filterOnOther_;

     gm2geom::CoordSystemsStoreData cs_;

     artg4::PhysicalVolumeStoreData pvs_;
};


gm2strawtracker::GeaneProcessSelection::GeaneProcessSelection(fhicl::ParameterSet const & pset) 
    : TrackModuleLabel_( pset.get<std::string>("TrackModuleLabel", "wronggeanemodulelabel") )
    , TrackInstanceName_( pset.get<std::string>("TrackInstanceName", "wronggeaneinstancename") )
    , trajectoryModuleLabel_( pset.get<std::string>("trajectoryModuleLabel", "artg4") )
    , trajectoryInstanceName_( pset.get<std::string>("trajectoryInstanceName", "") )
    , pvsModuleLabel_ ( pset.get<std::string>("pvsModuleLabel", "artg4") )
    , pvsInstanceName_ ( pset.get<std::string>("pvsInstanceName", "physicalVolumeStore") )
    , filterOnIonizations_(pset.get<bool>("filterOnIonizations", false))
    , filterOnBrems_(pset.get<bool>("filterOnBrems", false))
    , filterOnTransportation_(pset.get<bool>("filterOnTransportation", false))
    , filterOnOther_(pset.get<bool>("filterOnOther", false))
    , cs_()
{}


gm2strawtracker::GeaneProcessSelection::~GeaneProcessSelection()
{
   // Clean up dynamic memory and other resources here.
}

bool gm2strawtracker::GeaneProcessSelection::beginRun(art::Run & r)
{
    //Get coord systems
    cs_ = artg4::dataFromRunOrService<gm2geom::CoordSystemsStoreData, gm2geom::CoordSystemsStore>
          ( r, dataModuleDefs::coordSysModuleLabel(),dataModuleDefs::coordSysInstanceLabel() );
    if( cs_.size() == 0 ) {
      mf::LogWarning("GeaneProcessSelection") << "This run does not contain any data associated with the coordinate system\n";
    }

    pvs_ = artg4::dataFromRunOrService<artg4::PhysicalVolumeStoreData, artg4::PhysicalVolumeStoreService>(r, pvsModuleLabel_, pvsInstanceName_);
    if( pvs_.size() == 0 ) {
      mf::LogWarning("GeaneProcessSelection") << "This run does not contain any data associated with the physical volums store\n";
    }

  return true;
}


bool gm2strawtracker::GeaneProcessSelection::filter(art::Event & e) 
{
 
   mf::LogInfo info("GeaneProcessSelection");

   // Get track data 
   art::Handle<gm2strawtracker::TrackArtRecordCollection> TrackDataHandle;
   bool foundTrackcollection = e.getByLabel(TrackModuleLabel_,TrackInstanceName_,TrackDataHandle);
   if( ! foundTrackcollection ) {
     throw cet::exception("GeaneProcessSelection") << "No Trackcollection in this event (\"" << TrackModuleLabel_ << "\":\"" << TrackInstanceName_ << "\")\n";
     // return;
   }

   gm2strawtracker::TrackArtRecordCollection const& tracks = *TrackDataHandle;

   art::Handle<gm2truth::TrajectoryArtRecordCollection> trajDataHandle;
   bool success = e.getByLabel(trajectoryModuleLabel_, trajectoryInstanceName_, trajDataHandle);
   if( ! success ) {
     throw cet::exception("TrajectoryAna") << "No trajectories in this event (\"" << trajectoryModuleLabel_ << "\":\"" << trajectoryInstanceName_ << "\")\n";
     // return;
   }

   gm2truth::TrajectoryArtRecordCollection const& trajectories = *trajDataHandle;

   std::pair<int, double> largestEdepProcess; // process, edep
   largestEdepProcess.first  = -10;
   largestEdepProcess.second = 0.0;

   std::vector< std::pair<int, double> > eDepProcesses;
   eDepProcesses.emplace_back(2, 0.0);
   eDepProcesses.emplace_back(3, 0.0);
   eDepProcesses.emplace_back(92491, 0.0);
   eDepProcesses.emplace_back(-10, 0.0);

   // Loop over tracks
   for( auto & track : tracks ) { 

     // Get GEANE collection
     auto geaneHitsOnTrack = track.geaneHits;

     // Selection Cut
     if( track.failureMode != 0 ) 
     {
       return false;
     }

     // Get the Planes
     int firstPlane = track.trackPlanesHitList.front();
     int lastPlane  = track.trackPlanesHitList.back();


     // Loop over Truth Information
     for( auto traj : trajectories)
     {
        if (traj.trackID != 4) continue; // FIXME will not work for beam gun # 4 is the positron trajectory for the gas gun

        double totalEdep    = 0;
        double eDepBetwHits = 0;

        for (int i = 0; i < int(traj.points.size()); ++i)
        {
          auto point = traj.points.at(i);

          int stationNumber = track.candidate->strawDigits.at(0)->wireID.getStation();

          stringstream stationStream;
          stationStream << "TrackerStation[" << stationNumber << "]";
          string stationStr = stationStream.str();

          gm2geom::CoordSystem3Vector pointPos(point.position, "world");
          gm2geom::CoordSystem3Vector pointPosGeane = pointPos.transform(cs_, stationStr);

          double pointEdep = point.edep;
          int subtype      = point.postProcessSubType;

          if (pointPosGeane.z() > geaneHitsOnTrack.startingGeaneParameters.at(2) && pointPosGeane.z() < geaneHitsOnTrack.planeZPositions.at(lastPlane))
          {
             eDepBetwHits += pointEdep;

             if(pointEdep > largestEdepProcess.second){
              largestEdepProcess.first = subtype;
              largestEdepProcess.second = pointEdep;
             }

             if (subtype==2)                       eDepProcesses.at(0).second+=pointEdep;
             else if (subtype==3)                  eDepProcesses.at(1).second+=pointEdep;
             else if (subtype==92 || subtype==491) eDepProcesses.at(2).second+=pointEdep;
             else                                  eDepProcesses.at(3).second+=pointEdep;
          }

          totalEdep += pointEdep;
        }
      } // end loop over trajectories (2, 1 of which is used)

      std::sort(eDepProcesses.begin(), eDepProcesses.end(), [](auto &left, auto &right) {
        return left.second > right.second;
      });

      // Need to filter in a smarter way...
      if (eDepProcesses.at(0).first==2 && filterOnIonizations_)             return true;
      else if (eDepProcesses.at(0).first==3 && filterOnBrems_)              return true;
      else if (eDepProcesses.at(0).first==92491 && filterOnTransportation_) return true;
      else if (eDepProcesses.at(0).first==-10 && filterOnOther_)            return true;

   } // end loop over tracks

   return false;
}

DEFINE_ART_MODULE(gm2strawtracker::GeaneProcessSelection)
