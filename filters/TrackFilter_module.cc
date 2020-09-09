//Kill events that don't contain hits in any straws

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

// data product includes
#include "gm2dataproducts/mc/strawtracker/StrawArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"

// common util
#include "gm2util/common/dataModuleDefs.hh"

namespace gm2strawtracker {
  class TrackFilter;
}

class gm2strawtracker::TrackFilter : public art::EDFilter {

   public:
     explicit TrackFilter(fhicl::ParameterSet const & p);
     virtual ~TrackFilter();

     bool filter(art::Event & e) override;

   private:

     // fhicl parameters 
    std::string strawBuilderModuleLabel;
    std::string strawBuilderInstanceName;

    std::string TrackModuleLabel_;
    std::string TrackInstanceName_;

    bool filterOnStraws_;
    bool filterOnTracks_;

};


gm2strawtracker::TrackFilter::TrackFilter(fhicl::ParameterSet const & p) 
    : strawBuilderModuleLabel( p.get<std::string>("strawModuleLabel", dataModuleDefs::strawModuleLabel()) )
    , strawBuilderInstanceName( p.get<std::string>("strawInstanceName", dataModuleDefs::strawInstanceLabel()) )
    , TrackModuleLabel_( p.get<std::string>("TrackModuleLabel", dataModuleDefs::trackModuleLabel()) )
    , TrackInstanceName_( p.get<std::string>("TrackInstanceName", dataModuleDefs::recoInstanceLabel()) )
    , filterOnStraws_( p.get<bool>("filterOnStraws",false) )
    , filterOnTracks_( p.get<bool>("filterOnTracks",true) )
{
}


gm2strawtracker::TrackFilter::~TrackFilter()
{
   // Clean up dynamic memory and other resources here.
}


bool gm2strawtracker::TrackFilter::filter(art::Event & e) 
{

  gm2truth::StrawArtRecordCollection straws;

  //Get straw hits
  if(filterOnStraws_){
    art::Handle<gm2truth::StrawArtRecordCollection> strawDataHandle;
    bool foundStrawHits = e.getByLabel(strawBuilderModuleLabel,strawBuilderInstanceName,strawDataHandle);
    if( ! foundStrawHits ) {
      throw cet::exception("TrackFilter") << "No straw hit data in this event (\"" << strawBuilderModuleLabel << "\":\"" << strawBuilderInstanceName << "\")\n";
    }
    straws = *strawDataHandle; //Resolve handle to get collection
  }

  gm2strawtracker::TrackArtRecordCollection tracks;

  if(filterOnTracks_){
      art::Handle<gm2strawtracker::TrackArtRecordCollection> TrackDataHandle;
      bool foundTrackcollection = e.getByLabel(TrackModuleLabel_,TrackInstanceName_,TrackDataHandle);
      if( ! foundTrackcollection ) {
        throw cet::exception("TrackFilter") << "No Trackcollection in this event (\"" << TrackModuleLabel_ << "\":\"" << TrackInstanceName_ << "\")\n";
    }
    tracks = *TrackDataHandle;
  }


  //Reject event if no hits/tracks

  if(filterOnStraws_){
    if( straws.size() == 0 ) {
      return false;
    } 
    else {
      return true;
    }
  }

  if(filterOnTracks_){
    if( tracks.size() == 0 ) {
      return false;
    } 
    else {
      return true;
    }
  }

  // for(auto& track : tracks){
  //   if(filterOnTracks_){
  //     if(track.failureMode == 0 && track.pValue <= 0.01) return true;
  //     else return false;
  //   }
  // }

  return true;

}

DEFINE_ART_MODULE(gm2strawtracker::TrackFilter)
