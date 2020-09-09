//Module that selects track candidates based on number of hits per layer etc
//James Mott (Jan 2018)
//

// Include needed ART headers
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

//art records
#include "gm2dataproducts/strawtracker/StrawTimeIslandArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackCandidateArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"

//Straw Geometry
#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"
#include "gm2tracker/utils/StrawObjectSorter.hh"

//Coord systems
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh"  

//Util
#include "gm2util/common/dataModuleDefs.hh"

//C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <math.h> 

using std::vector;
using std::set;
using std::string;

namespace gm2strawtracker {

  //
  // Class declaration
  //
  class GeaneTrackSelector : public art::EDFilter {

  public:

    explicit GeaneTrackSelector(fhicl::ParameterSet const& pset);

    //Override desired art::EDFilter functions
    bool filter(art::Event& event) override;
    void beginJob() override;

  private:

    //Filter labels
    string TrackModuleLabel_;
    string TrackInstanceName_;
    string outputInstanceName_;

    //Parameters for choosing tracks
    int minHitsPerLayer_;
    int maxHitsPerLayer_;
    int minUorVHits_;
    int maxUorVHits_;
    int minUniqueModules_;
    int maxUniqueModules_;
    int minTracksPerTimeIsland_;
    int maxTracksPerTimeIsland_;
    int minLayersWithMultiHits_;
    int maxLayersWithMultiHits_;

    //Vectors of wires that we'll ignore when choosing tracks - fcl params and holder for module
    vector<short> ignoreModules_;
    vector<short> ignoreViews_;
    vector<short> ignoreLayers_;
    vector<short> ignoreWires_;
    set<WireID> ignoreStraws_;

    //Straw geometry
    gm2geom::StrawTrackerGeometry geom_;

    string name_;

  }; //End of class GeaneTrackSelector

  //
  // Class implementation
  //
  GeaneTrackSelector::GeaneTrackSelector(fhicl::ParameterSet const& pset)
    : TrackModuleLabel_( pset.get<string>("TrackModuleLabel",dataModuleDefs::trackModuleLabel()) )
    , TrackInstanceName_( pset.get<string>("TrackInstanceName",dataModuleDefs::recoInstanceLabel()) )
    , outputInstanceName_( pset.get<string>("outputInstanceName",dataModuleDefs::recoInstanceLabel()) )
    , minHitsPerLayer_( pset.get<int>("minHitsPerLayer", -1) )
    , maxHitsPerLayer_( pset.get<int>("maxHitsPerLayer", -1) )
    , minUorVHits_( pset.get<int>("minUorVHits", -1) )
    , maxUorVHits_( pset.get<int>("maxUorVHits", -1) )
    , minUniqueModules_( pset.get<int>("minUniqueModules", -1) )
    , maxUniqueModules_( pset.get<int>("maxUniqueModules", -1) )
    , minTracksPerTimeIsland_( pset.get<int>("minTracksPerTimeIsland", -1) )
    , maxTracksPerTimeIsland_( pset.get<int>("maxTracksPerTimeIsland", -1) )
    , minLayersWithMultiHits_( pset.get<int>("minLayersWithMultiHits", 0) )
    , maxLayersWithMultiHits_( pset.get<int>("maxLayersWithMultiHits", 32) )
    , ignoreModules_( pset.get<vector<short> >("ignoreModules", {}) )
    , ignoreViews_( pset.get<vector<short> >("ignoreViews", {}) )
    , ignoreLayers_( pset.get<vector<short> >("ignoreLayers", {}) )
    , ignoreWires_( pset.get<vector<short> >("ignoreWires", {}) ) 
    , ignoreStraws_()
    , geom_()
    , name_("GeaneTrackSelector")
  { 
    produces<TrackArtRecordCollection>(outputInstanceName_);
  }

  bool GeaneTrackSelector::filter(art::Event& event) {

    mf::LogVerbatim(name_) << "Event : " << event.id() << "\n";

    //Create new collection that we'll fill with our chosen track candidates
    std::unique_ptr<TrackArtRecordCollection> selectedTracks(new TrackArtRecordCollection);

    //Get the tracks
    art::Handle<gm2strawtracker::TrackArtRecordCollection> trackDataHandle;
    bool foundTracks = event.getByLabel(TrackModuleLabel_,TrackInstanceName_,trackDataHandle);
    if( ! foundTracks ) {
      throw cet::exception("GeaneTrackSelector") << "No tracks in this event (\"" << TrackModuleLabel_ << "\":\"" << TrackInstanceName_ << "\")\n";
      return false;
    }

    auto tracks = *trackDataHandle;

    // Flag to know whether to keep whole event
    bool acceptedTrackInEvent = false;

    // island number and track map?
    std::map<int, gm2strawtracker::TrackArtRecordCollection > islandNumToTrackMap;

    // Loop over tracks and use station & island number as the key to a map
    for( auto track : tracks ) {
      islandNumToTrackMap[ 1000*track.island->station + track.island->islandNumber ].push_back(track); 
    }

    // index
    int imap = 0;

    // Loop over hits in the track
    for( auto &ipair : islandNumToTrackMap ) {
 
      mf::LogVerbatim(name_) << "Island " << ipair.first << "\n";
      int nTracksThisIsland = ipair.second.size();

      // check any island conditions...

      // Create a map to store number of hits per layer (<global layer, number of hits>)
      std::map<unsigned int,int> islandLayerHits;
	
      // Fill map with zeros (so we can cut on minimums)
      for (unsigned int layer = 0; layer < geom_.getTotalNumLayers(); layer++){
	islandLayerHits[layer] = 0;
      }

      // Loop over digits in track and fill map of hits vs layer number
      for( auto digit : tracks.at(imap).island->strawDigits ) {
	// Only increment map counters if we're not ignoring this straw
	if(std::find(ignoreStraws_.begin(), ignoreStraws_.end(), digit->wireID) == ignoreStraws_.end()){
	  std::map<unsigned int,int>::iterator iter(islandLayerHits.find(geom_.getGlobalLayer(digit->wireID)));
	  iter->second++;
	}
      }

      int layersWithMultiHits = 0;
      for( auto& layerHitsIter : islandLayerHits ) if ( layerHitsIter.second > 1) layersWithMultiHits++;
      mf::LogVerbatim(name_) << "Layers with multiple hits: " << layersWithMultiHits << "\n";

      if (layersWithMultiHits < minLayersWithMultiHits_) continue;
      if (layersWithMultiHits > maxLayersWithMultiHits_) continue;
      mf::LogVerbatim(name_) << "this track passes multi layer hit checker... \n";
      
      // Loop over tracks in from this island
      for (auto track : ipair.second){

	mf::LogVerbatim(name_) << "number of tracks from this island: " << nTracksThisIsland << ", " << tracks.at(imap).candidate->island << " with: " <<tracks.at(imap).candidate->strawDigits.size() <<"\n";

	// filter on tracks with multiple tracks per island
	if ( maxTracksPerTimeIsland_ > 0 && nTracksThisIsland > maxTracksPerTimeIsland_) continue;
	if ( nTracksThisIsland < minTracksPerTimeIsland_) continue;

	mf::LogVerbatim(name_) << "this track passed multi selector... \n";

	// Create a map to store number of hits per layer (<global layer, number of hits>)
	std::map<unsigned int,int> layerHits;
	
	// Fill map with zeros (so we can cut on minimums)
	for (unsigned int layer = 0; layer < geom_.getTotalNumLayers(); layer++){
	  layerHits[layer] = 0;
	}

	// Loop over digits in track and fill map of hits vs layer number
	for( auto digit : tracks.at(imap).candidate->strawDigits ) {
	  // Only increment map counters if we're not ignoring this straw
	  if(std::find(ignoreStraws_.begin(), ignoreStraws_.end(), digit->wireID) == ignoreStraws_.end()){
	    std::map<unsigned int,int>::iterator iter(layerHits.find(geom_.getGlobalLayer(digit->wireID)));
	    iter->second++;
	  }
	}

	// Check number of hits in each layer passes minimum/maximum set by user - use flag to continue from track loop (not layerHits loop)
	bool passNHitsCheck = true;
	for( auto& layerHitsIter : layerHits ) {
	  if ( minHitsPerLayer_ > 0 and layerHitsIter.second < minHitsPerLayer_) passNHitsCheck = false;
	  if ( maxHitsPerLayer_ > 0 and layerHitsIter.second > maxHitsPerLayer_) passNHitsCheck = false;
	}
	if (!passNHitsCheck) continue;
	mf::LogVerbatim(name_) << "this track passes n hits checker... \n";

	// Count number of U/V layer hits
	int nULayerHits = 0;
	int nVLayerHits = 0; 
	for( auto& layerHitsIter : layerHits ) {
	  if ( (layerHitsIter.first/2) % 2 == 0 ) nULayerHits += layerHitsIter.second;
	  else                                    nVLayerHits += layerHitsIter.second;
	}
	
	// Check that we've got required numbers of U & V hits
	if( minUorVHits_ > 0 and (nULayerHits < minUorVHits_ or nVLayerHits < minUorVHits_)) continue;
	if( maxUorVHits_ > 0 and (nULayerHits > maxUorVHits_ or nVLayerHits > maxUorVHits_)) continue;
	
	mf::LogVerbatim(name_) << "this track passes n U and V checks... \n";

	// make a set of the unique module numbers
	std::set<int> modulesHit;
	for( auto digit : tracks.at(imap).candidate->strawDigits ) modulesHit.insert(digit->wireID.getModule());
	
	int nUniqueModules = int(modulesHit.size());
	
	mf::LogVerbatim(name_) << "n unique modules: " << nUniqueModules << " \n";

	// Check that we've got required numbers of unique modules
	if( minUniqueModules_ > 0 and (nUniqueModules < minUniqueModules_ )) continue;
	if( maxUniqueModules_ > 0 and (nUniqueModules > maxUniqueModules_ )) continue;
	
	mf::LogVerbatim(name_) << "this track passes unique modules... \n";

	// If we've got this far, then we're going to keep track
	selectedTracks->push_back(track);
	
	mf::LogVerbatim(name_) << "added this track to event \n";
	
	// Record number for stats (and set event-level flag)
	acceptedTrackInEvent = true;

	mf::LogVerbatim(name_) << "ACCEPTED TRACK : " << event.id() << " : station " << tracks.at(imap).candidate->island->station << " : island " << tracks.at(imap).candidate->island->islandNumber << "\n";
      } // for island

      // increment
      imap += 1;
    }

    //Put the raw digits in the art event
    event.put(std::move(selectedTracks),outputInstanceName_);

    //Filter overall events based on whether there was or wasn't a time island we want
    if(acceptedTrackInEvent){
      return true;
    } else {
      return false;
    }

  }//filter


  void GeaneTrackSelector::beginJob() {
  
    // Check that all vector sizes are OK sizes
    if (ignoreViews_.size() != ignoreModules_.size() and ignoreViews_.size() > 0){
      throw cet::exception("GeaneTrackSelector") << "Size of fcl parameter vectors are not the same: ignoreViews(" << ignoreViews_.size() << "), ignoreModules(" << ignoreModules_.size() << ")\n";
    } else if (ignoreLayers_.size() != ignoreViews_.size() and ignoreLayers_.size() > 0){
      throw cet::exception("GeaneTrackSelector") << "Size of fcl parameter vectors are not the same: ignoreLayers(" << ignoreLayers_.size() << "), ignoreViews(" << ignoreViews_.size() << ")\n";
    } else if (ignoreWires_.size() != ignoreLayers_.size() and ignoreWires_.size() > 0){
      throw cet::exception("GeaneTrackSelector") << "Size of fcl parameter vectors are not the same: ignoreWires(" << ignoreWires_.size() << "), ignoreLayers(" << ignoreLayers_.size() << ")\n";
    }

    // Put all ignored wires in set
    // If anything is blank in view->layer->wire chain, then fill all (and all below this point in chain)
    if(ignoreWires_.size() > 0){
      for(unsigned int i_wire = 0; i_wire < ignoreWires_.size(); i_wire++){
	WireID wire(0, ignoreModules_.at(i_wire), StrawView(ignoreViews_.at(i_wire)), ignoreLayers_.at(i_wire), ignoreWires_.at(i_wire));
	ignoreStraws_.insert(wire);
      }
    } else if (ignoreLayers_.size() > 0){
      for(unsigned int i_layer = 0; i_layer < ignoreLayers_.size(); i_layer++){
	for(unsigned int i_wire = 0; i_wire < geom_.getNumWiresPerLayer(); i_wire++){
	  WireID wire(0, ignoreModules_.at(i_layer), StrawView(ignoreViews_.at(i_layer)), ignoreLayers_.at(i_layer), i_wire);
	  ignoreStraws_.insert(wire);
	}
      }
    } else if(ignoreViews_.size() > 0){
      for(unsigned int i_view = 0; i_view < ignoreViews_.size(); i_view++){
	for(unsigned int i_layer = 0; i_layer < geom_.getNumLayersPerView(); i_layer++){
	  for(unsigned int i_wire = 0; i_wire < geom_.getNumWiresPerLayer(); i_wire++){
	    WireID wire(0, ignoreModules_.at(i_view), StrawView(ignoreViews_.at(i_view)), i_layer, i_wire);
	    ignoreStraws_.insert(wire);
	  }
	}
      }
    } else if (ignoreModules_.size() > 0){
      for(unsigned int i_mod = 0; i_mod < ignoreModules_.size(); i_mod++){
	for(unsigned int i_view = 0; i_view < geom_.getNumViewsPerModule(); i_view++){
	  for(unsigned int i_layer = 0; i_layer < geom_.getNumLayersPerView(); i_layer++){
	    for(unsigned int i_wire = 0; i_wire < geom_.getNumWiresPerLayer(); i_wire++){
	      WireID wire(0, ignoreModules_.at(i_mod), StrawView(i_view), i_layer, i_wire);
	      ignoreStraws_.insert(wire);
	    }
	  }
	}
      }
    }
  }


} // End of namespace gm2strawtracker

//
// Extras
//

// These are some necessary boilerplate for the ROOT persistency system
using gm2strawtracker::GeaneTrackSelector;
DEFINE_ART_MODULE(GeaneTrackSelector)
