////////////////////////////////////////////////////////////////////////
// Class:       GeaneExtrapolation
// Module Type: producer
// File:        GeaneExtrapolation_module.cc
//
// Generated at Tue Jan 27 14:11:23 2015 by Tammy Walton using artmod
// from cetpkgsupport v1_07_01.
//
// Modified by Saskia Charity on Jun 23
////////////////////////////////////////////////////////////////////////

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// cern 
#include "CLHEP/Units/SystemOfUnits.h"

// artg4 includes
#include "artg4/util/DataFromRunOrService.hh"
#include "artg4/util/util.hh"

// data product includes
#include "gm2dataproducts/mc/actions/track/TrajectoryArtRecord.hh"
#include "gm2dataproducts/mc/strawtracker/StrawArtRecord.hh"
#include "gm2dataproducts/mc/calo/CaloArtRecord.hh"
#include "gm2dataproducts/mc/actions/track/TrackingActionArtRecord.hh"
#include "gm2dataproducts/strawtracker/DecayVertexArtRecord.hh"

// geometry includes
#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"
#include "gm2geom/coordSystems/CoordSystem.hh"
#include "gm2geom/coordSystems/CoordSystemsStoreData.hh"
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh"  
#include "gm2util/coordSystems/CoordSystemUtils.hh"

// utils
#include "gm2util/common/dataModuleDefs.hh"
#include "gm2util/coordSystems/CoordSystem3VectorSort.hh"
#include "gm2tracker/utils/StrawObjectSorter.hh"
#include "gm2tracker/utils/GeaneExtrapolationUtils.hh"
#include "gm2tracker/utils/GeaneEigenStorageUtils.hh"
// c++ includes
//#include <utility>
#include<string>
#include<iostream>
#include<sstream>
#include<boost/format.hpp>

//root include
#include "TTree.h"
#include "TFile.h"

// namespace
using namespace gm2strawtracker;
using gm2geom::CoordSystem3Vector;

namespace gm2strawtracker {
  class GeaneExtrapolation;
}

class gm2strawtracker::GeaneExtrapolation : public art::EDProducer {

  public:
    explicit GeaneExtrapolation(fhicl::ParameterSet const & p);

  //GeaneExtrapolation(GeaneExtrapolation const &) = delete;
  //GeaneExtrapolation(GeaneExtrapolation &&) = delete;

  //GeaneExtrapolation & operator = (GeaneExtrapolation const &) = delete;
  //GeaneExtrapolation & operator = (GeaneExtrapolation &&) = delete;

    void produce(art::Event & e) override;
    void beginJob() override;
    void beginRun(art::Run &) override;
    void endRun(art::Run & r) override;
    void endJob() override;

  private:
    std::string decayVertexInstanceLabel_;
    std::string forwardsExtrapInstanceLabel_;
    std::string failedBackwardsExtrapInstanceLabel_;
    std::string failedForwardsExtrapInstanceLabel_;
    std::string trackModuleLabel_;
    std::string trackInstanceLabel_;
    std::string trajModuleLabel_;
    std::string trajInstanceLabel_;
    std::string caloModuleLabel_;
    std::string caloInstanceLabel_;
    std::string csModuleLabel_;
    std::string csInstanceLabel_;

    // name of the track extrapolation routine
    std::string extrapolater_;

    bool doBackwardsExtrap_;
    bool doForwardsExtrap_;
  
    // whether to check hit volumes
    bool checkVolumesHit_;
  
    // whether to go to decay vertex or true azimuth
    bool extrapolateToTrueAzimuth_;

    // straw tracker geometry
    gm2geom::StrawTrackerGeometry sgeom_;

    // helper function
    GeaneExtrapolationUtils utils_;
  //GeaneExtrapolationUtils utilsWorld_;
    GeaneFhiclOptions fclOptions_;
    gm2geom::CoordSystemsStoreData cs_;
    gm2util::CoordSystemUtils csUtils_;
    std::map< std::string, gm2geom::CoordSystemsStoreData > detCoordMap_;
  //art::ServiceHandle<artg4::DetectorHolderService> detectorHolder_;

    // station info
    std::set<int> stations_;

    // counters
    int counter_backwards_;
    int counter_forwards_;
    int counter_failed_backwards_;
    int counter_failed_forwards_;
    int counter_tracks_;
   
    // message facilty name
    std::string name_;
};

// standard constructor
gm2strawtracker::GeaneExtrapolation::GeaneExtrapolation(fhicl::ParameterSet const & p) 
   : decayVertexInstanceLabel_(p.get<std::string>("decayVertexInstanceLabel",dataModuleDefs::backExtrapInstanceLabel()) )
   , forwardsExtrapInstanceLabel_(p.get<std::string>("forwardsExtrapInstanceLabel",dataModuleDefs::forwardExtrapInstanceLabel()) )
   , failedBackwardsExtrapInstanceLabel_(p.get<std::string>("failedBackwardsExtrapInstanceLabel","failedBackwardsExtrapolation") )
   , failedForwardsExtrapInstanceLabel_(p.get<std::string>("failedForwardsExtrapInstanceLabel","failedForwardsExtrapolation") )
   , trackModuleLabel_(  p.get<std::string>("trackModuleLabel",dataModuleDefs::trackModuleLabel()) )
   , trackInstanceLabel_( p.get<std::string>("trackInstanceLabel",dataModuleDefs::recoInstanceLabel()) )
   , trajModuleLabel_(  p.get<std::string>("trajModuleLabel",dataModuleDefs::trajModuleLabel()) )
   , trajInstanceLabel_( p.get<std::string>("trajInstanceLabel",dataModuleDefs::trajInstanceLabel()))
   , caloModuleLabel_(  p.get<std::string>("caloModuleLabel","artg4") )
   , caloInstanceLabel_( p.get<std::string>("caloInstanceLabel","calorimeter"))
   , csModuleLabel_( p.get<std::string>("csModuleLabel",dataModuleDefs::coordSysModuleLabel()) )
   , csInstanceLabel_( p.get<std::string>("csInstanceLabel",dataModuleDefs::coordSysInstanceLabel()) )
   , extrapolater_( p.get<std::string>("extrapolater", "GeaneExtrapolater") )
   , doBackwardsExtrap_(p.get<bool>("doBackwardsExtrap",true))
   , doForwardsExtrap_(p.get<bool>("doForwardsExtrap",true))
   , checkVolumesHit_( p.get<bool>("checkVolumes", true) )
   , extrapolateToTrueAzimuth_( p.get<bool>("extrapolateToTrueAzimuth",false) )
   , sgeom_()
   , utils_()
   , fclOptions_(p)
   , cs_()
   , csUtils_()
   , detCoordMap_()
   , stations_{0,12,18}
   , name_("GeaneExtrapolation")
   
{
  //std::cout<<"ryan module constructor"<<std::endl;
  produces<DecayVertexArtRecordCollection>(decayVertexInstanceLabel_);
  produces<DecayVertexArtRecordCollection>(forwardsExtrapInstanceLabel_);
  produces<DecayVertexArtRecordCollection>(failedBackwardsExtrapInstanceLabel_);
  produces<DecayVertexArtRecordCollection>(failedForwardsExtrapInstanceLabel_);
  produces<std::vector<ExtrapolationStep>>("volumeSteps");
  utils_.setFhiclCuts( p );
  counter_backwards_ = 0;
  counter_forwards_ = 0;
  counter_failed_backwards_ = 0;
  counter_failed_forwards_ = 0;
  counter_tracks_ = 0;
}

void gm2strawtracker::GeaneExtrapolation::beginJob() {
  //std::cout<<"ryan module beginJob"<<std::endl;
}

void gm2strawtracker::GeaneExtrapolation::beginRun(art::Run & r)
{
  //std::cout<<"ryan module beginRun"<<std::endl;
  mf::LogInfo(name_) << "Start of beginRun, r.id() = " << r.id() << "\n";

  // get the coord systems
  cs_ = artg4::dataFromRunOrService<gm2geom::CoordSystemsStoreData, gm2geom::CoordSystemsStore>(r,csModuleLabel_,csInstanceLabel_);
  
  if( cs_.size() == 0 ){
    throw cet::exception(name_) << "This run does not contain any data associated with the coordinate system\n";
    mf::LogInfo(name_) << "This run does not contain any data associated with the coordinate system\n\n";
  }

  //set up the tracker system
  std::vector<std::string> detNames;
  detNames.push_back("TrackerStation");
  int nprob = 0;
  for(auto s : sgeom_.whichScallopLocations) {
    for(unsigned int m = 0; m < sgeom_.getNumModulesPerStation(); ++m) {
       std::cout <<"s variable "<< s <<"\n";
       detNames.push_back(Form("Module%d:%d", (int)s, m));
       std::cout << "Detector names " << detNames[nprob] << "\n";
       nprob++;
    }
  }

  csUtils_.setDetectorNames(detNames);
  csUtils_.detNameToCoordMap(cs_,detCoordMap_);
  sgeom_.setDetNameToCoords(detCoordMap_);

  // set up the detector coord system map
  // detNames.push_back("TrackerCalos");
  // csUtils_.setDetectorNames(detNames);
  // csUtils_.detNameToCoordMap(cs_,detCoordMap_);
  // sgeom_.setDetNameToCoords(detCoordMap_);

  // cs_ = detCoordMap_.find("TrackerCalos")->second;
  utils_.setCoordSysData(cs_);

  // pass cs to dummy utils for coord sys transforms  
  //  geaneFittingUtils_.geaneParamUtils_.dummyUtils_.fillCS(detCoordMap_.find("TrackerStation")->second); 
  utils_.setCoordMap(detCoordMap_);



  mf::LogInfo(name_) << "End of beginRun " << "\n\n";
}

void gm2strawtracker::GeaneExtrapolation::endRun(art::Run & r)
{
  //std::cout<<"ryan module endRun"<<std::endl;
  mf::LogInfo(name_) << "For run= " << r.id() << ", TOTAL TRACKS  = " << counter_tracks_ << "\n\n";
  mf::LogInfo(name_) << "For run= " << r.id() << ", TOTAL DECAY VERTICES  = " << counter_backwards_ << "\n\n";
  mf::LogInfo(name_) << "For run= " << r.id() << ", TOTAL FORWARDS EXTRAP = " << counter_forwards_ << "\n\n";
  mf::LogInfo(name_) << "For run= " << r.id() << ", TOTAL FAILED (BACKWARDS) = " << counter_failed_backwards_ << "\n\n";
  mf::LogInfo(name_) << "For run= " << r.id() << ", TOTAL FAILED (FORWARDS)  = " << counter_failed_forwards_ << "\n\n";
}

void gm2strawtracker::GeaneExtrapolation::endJob(){}


// produces the decay vertices
void gm2strawtracker::GeaneExtrapolation::produce(art::Event & e)
{
  //std::cout<<"ryan module produce"<<std::endl;
  //std::cout<<"Start of produce"<<std::endl;
  mf::LogInfo(name_) << "Enter GeaneExtrapolation::produce, e.event() = " << e.event() << "\n";
  
  // field service
  art::ServiceHandle<gm2geom::gm2FieldManager> fieldManager;
  
  // create a new collection of DecayVertexArtRecords that will be added to the event
  std::unique_ptr<DecayVertexArtRecordCollection> decayVertices(new DecayVertexArtRecordCollection);

  std::unique_ptr<std::vector<ExtrapolationStep>> volumeSteps(new std::vector<ExtrapolationStep>());

  // create a new collection of DecayVertexArtRecords that will be added to the event
  std::unique_ptr<DecayVertexArtRecordCollection> forwardsExtrapPoints(new DecayVertexArtRecordCollection);
  
  // failure modes
  std::unique_ptr<DecayVertexArtRecordCollection> failedBackwardsExtrapolationEvents(new DecayVertexArtRecordCollection);
  std::unique_ptr<DecayVertexArtRecordCollection> failedForwardsExtrapolationEvents(new DecayVertexArtRecordCollection);
  
  // success condition for the vertex reconstruction
  bool success = false;
  bool forwardSuccess = false;
  bool backwardSuccess = false;
  
  // get the reconstructed track
  art::Handle<gm2strawtracker::TrackArtRecordCollection> trackDataHandle;
  gm2strawtracker::TrackArtRecordCollection trackData;
  //art::Handle<gm2strawtracker::GEANEArtRecordCollection> tracksHandle;
  
  //bool success2 = false;
  //success2 = e.getByLabel(trackModuleLabel_,trackInstanceLabel_,tracksHandle);

  //if(!success2) {
  //  std::cout<<"No GEANEArtRecordCollection"<<std::endl;
  //}

  //gm2strawtracker::GEANEArtRecordCollection tracks = GeaneEigenStorageUtils::ReadEigenFromStorage(*tracksHandle);

  success = e.getByLabel(trackModuleLabel_,trackInstanceLabel_,trackDataHandle);
  
  if( !success ) {
    throw cet::exception(name_) << "Event " << e.id() << " does not contain any track data \"" << trackModuleLabel_ << ":" << trackInstanceLabel_  <<"\"\n";
  }
  
  trackData = *trackDataHandle;
  
  if (trackData.size() == 0 ){
    mf::LogWarning(name_) << "No tracks in this event" << "\n";
  }
  
  mf::LogInfo(name_) << "NUMBER OF TRACKS AT START OF PRODUCE: " << trackData.size() << "\n";

  Eigen::MatrixXd zeroMatrix = Eigen::MatrixXd::Zero(5,5);

  for(unsigned int itrack = 0; itrack < trackData.size(); ++itrack) {
    //std::cout<<"trackData.size: "<<trackData.size()<<std::endl;
    //std::cout<<"tracks.size: "<<tracks.size()<<std::endl;
    //std::cout<<"geaneErrorMatrices.size(): "<<tracks[itrack].geaneErrorMatrices.size()<<std::endl;
#if 0
    bool errorLoop = true;
    while(errorLoop){
      if(tracks[itrack].geaneErrorMatrices.size()==0){
	std::cout<<"No error matrices"<<std::endl;
	tracks[itrack].geaneErrorMatrices.push_back(zeroMatrix);
	break;
      }
      Eigen::MatrixXd initialError = tracks[itrack].geaneErrorMatrices.back();
      if(initialError == zeroMatrix){
	tracks[itrack].geaneErrorMatrices.pop_back();
      }
      else{
	errorLoop = false;
      }
    }
#endif
    //Eigen::MatrixXd initialError = tracks[itrack].geaneErrorMatrices.back();
    std::vector<double> covTotalInverseData;
    covTotalInverseData = trackData[itrack].geaneHits.covarianceTotalInverseData;
    Eigen::MatrixXd covMatrixUVInverse(5,5);
    for (unsigned int i(0); i<covTotalInverseData.size(); i++) {
      covMatrixUVInverse(i) = covTotalInverseData[i];
    }
    std::cout<<"Ryan UV Inverse error matrix "<<covMatrixUVInverse<<"\n\n";
    Eigen::MatrixXd covMatrixUVEigen = covMatrixUVInverse.inverse();

    //    Eigen::MatrixXd covMatrixXYEigen = geaneTrackUtils_.JacobianToUV.inverse() * covMatrixUVEigen * geaneTrackUtils_.JacobianToUV.transpose().inverse(); 
   Eigen::MatrixXd initialError = covMatrixUVEigen;
   std::cout<<"Ryan UV error matrix "<<initialError<<"\n\n";
    //  Eigen::MatrixXd initialError = covMatrixXYEigen;
    //Eigen::MatrixXd initialError = zeroMatrix;
    //std::cout<<"initialError: "<<initialError<<std::endl;
    counter_tracks_++;
    // initialize decay vertex to be created
    DecayVertexArtRecord decayVertex;
    DecayVertexMC mcDecayVertex;
    DecayVertexArtRecord vertexSteps;
    
    // initialize forwards extrap point (at calo) to be created
    DecayVertexArtRecord forwardsExtrapPoint;
    DecayVertexMC mcForwardsExtrapPoint;
        
    if ( stations_.find(trackData[itrack].island->station) == stations_.end() ) {
      mf::LogInfo(name_) << "Track not from a station of interest, continue" << "\n";
      continue;
    }

    if ( trackData[itrack].failureMode != 0) {
      mf::LogInfo(name_) << "trackData[itrack].geanePlanes->failureMode = " << trackData[itrack].failureMode << " so continue" << "\n";
      continue;
    }
    
    // reconstruct the decay vertex
    if (!extrapolateToTrueAzimuth_) {
      //std::cout<<"initialError: "<<initialError<<std::endl;
      //std::cout<<"volumeSteps size: "<<volumeSteps->size()<<std::endl;
      if (doBackwardsExtrap_) backwardSuccess = utils_.reconstructDecayVertex(trackData[itrack], decayVertex, fclOptions_, initialError, vertexSteps);
      //std::cout<<"vertexSteps.steps size: "<<vertexSteps.steps.size()<<std::endl;
      std::cout<<"TrackCounter: "<<counter_tracks_<<std::endl;
      mf::LogInfo(name_) << "1 backwardSuccess = " << backwardSuccess << "\n";
      if (doForwardsExtrap_)  forwardSuccess  = utils_.extrapolateToCalorimeter(trackData[itrack], forwardsExtrapPoint, fclOptions_);
      mf::LogInfo(name_) << "1 forwardSuccess = " << forwardSuccess << "\n";
    }
    
    else if (!e.isRealData()) {
      ExtrapolationStep firstStep;
      firstStep.position = trackData[itrack].position;
      firstStep.momentum = trackData[itrack].momentum;
      success = utils_.extrapolateToTrueAzimuth(firstStep, decayVertex, fclOptions_);
    }
    
    mf::LogInfo(name_) << "success = " << success << "\n";
    mf::LogInfo(name_) << "doBackwardsExtrap_ = " << doBackwardsExtrap_ << ", doForwardsExtrap_ = " << doForwardsExtrap_ << "\n";

    // decay vertex reconstruction condition
    if(doBackwardsExtrap_ && !backwardSuccess ) {
      mf::LogInfo(name_) << "failed reconstructDecayVertex with failure mode " << decayVertex.failureMode << "\n";
      failedBackwardsExtrapolationEvents->push_back(decayVertex);
      counter_failed_backwards_++;
      //continue; 
    }
    
    // forward extrap reconstruction condition
    if(doForwardsExtrap_ && !forwardSuccess ) {
      mf::LogInfo(name_) << "failed reconstructDecayVertex with failure mode " << forwardsExtrapPoint.failureMode << "\n";
      failedForwardsExtrapolationEvents->push_back(forwardsExtrapPoint);
      counter_failed_forwards_++;
      //continue; 
    }
    
    if (doBackwardsExtrap_ && backwardSuccess && (decayVertex.position.coordSystemName == "NULL" || decayVertex.momentum.coordSystemName == "NULL")) {
      mf::LogInfo(name_) << "decay vertex coord sys name == NULL; vertex was not correctly filled, continue" << "\n";
      decayVertex.failureMode = "coordSystemIsNull";
      failedBackwardsExtrapolationEvents->push_back(decayVertex);
      counter_failed_backwards_++;
      //continue;
    }
    
    if (doForwardsExtrap_ && forwardSuccess && (forwardsExtrapPoint.position.coordSystemName == "NULL" || forwardsExtrapPoint.momentum.coordSystemName == "NULL")) {
      mf::LogInfo(name_) << "forward extrap coord sys name == NULL; forward extrap point was not correctly filled, continue" << "\n";
      forwardsExtrapPoint.failureMode = "coordSystemIsNull";
      failedForwardsExtrapolationEvents->push_back(forwardsExtrapPoint);
      counter_failed_forwards_++;
      //continue;
    }
    //std::cout<<"module decayVertex.position: "<<decayVertex.position<<std::endl;
    //std::cout<<"module decayVertex.momentum: "<<decayVertex.momentum<<std::endl;

    std::cout<<"Ryan volumesHit: "<<std::endl;
    for(auto volume : decayVertex.volumesHit){
      std::cout<<volume<<std::endl;
    }

    mf::LogInfo(name_) << "decayVertex.position    = " << decayVertex.position << "\n";
    mf::LogInfo(name_) << "decayVertex.momentum    = " << decayVertex.momentum << "\n";
    mf::LogInfo(name_) << "decayVertex.time    = " << decayVertex.time << ", t0: " << trackData[itrack].time << " forwardsExtrapPoint.time: " << forwardsExtrapPoint.time << "\n";
    //std::cout<<"module uncorrected steps size: "<<decayVertex.uncorrectedSteps.size()<<std::endl;
    
    decayVertex.t0 = trackData[itrack].time;
    //decayVertex.trackID = trackData[itrack].geanePlanes->candidate->strawDigits.at(0)->strawMCDigit.strawMCHits.at(0)->trackID;
    //std::cout<<"ryan t0: "<<decayVertex.t0<<std::endl;
    //std::cout<<"ryan trackID: "<<decayVertex.trackID<<std::endl;
    //std::cout<<"ryan island: "<<trackData[itrack].island<<std::endl;
    //mf::LogInfo(name_) << "forwardsExtrapPoint.position    = " << forwardsExtrapPoint.position.transform(cs_,"CalorimeterNumber[18]") << "\n";
    //mf::LogInfo(name_) << "forwardsExtrapPoint.momentum    = " << forwardsExtrapPoint.momentum.transform(cs_,"CalorimeterNumber[18]",true) << "\n";
    //std::cout<<"Before isRealData"<<std::endl;
    // if using sim, then get the true decay pos/mom and save it to the art record
    if ( !e.isRealData() ) {
      int vertexTrackID(0);
      vertexTrackID = trackData[itrack].island->strawDigits.at(0)->strawMCDigit.strawMCHits.at(0)->trackID;
      mf::LogInfo(name_) << "vertexTrackID = " << vertexTrackID << "\n";
      
      if (doBackwardsExtrap_) {
	art::Handle<gm2truth::TrajectoryArtRecordCollection> trajDataHandle;
	gm2truth::TrajectoryArtRecordCollection trajData;
	success = e.getByLabel(trajModuleLabel_,trajInstanceLabel_,trajDataHandle);
	//std::cout<<"test before"<<std::endl;
	//std::cout<<"trajDataHandle: "<<trajDataHandle<<std::endl;
	auto trajectories = *trajDataHandle;
	//std::cout<<"trajDataHandle: "<<trajDataHandle<<std::endl;
	//std::cout<<"test after"<<std::endl;
	G4ThreeVector trueDecayPos; 
	G4ThreeVector trueDecayMom;
	double trueTimeAtDecay = 0.0;
	//std::cout<<"Before traj for loop"<<std::endl;
	for (auto traj : trajectories) {
	  if (traj.trackID == vertexTrackID) {
	    trueDecayPos = traj.getVertex().position;
	    trueDecayMom = traj.getVertex().momentum;
	    trueTimeAtDecay = traj.getVertex().time;
	  }
	  else continue;
	}
	
	mf::LogInfo(name_) << "trueDecayPos = " << trueDecayPos << "\n";
	mf::LogInfo(name_) << "trueDecayMom = " << trueDecayMom << "\n";
	std::cout<<"trueDecayPos: "<<trueDecayPos<<std::endl;
	std::cout<<"trueDecayMom: "<<trueDecayMom<<std::endl;
	//std::cout<<"trueDecayRadius: "<<sqrt(pow(trueDecayPos[0],2)+pow(trueDecayPos[2],2))<<std::endl;
	mcDecayVertex.position = trueDecayPos;
	mcDecayVertex.momentum = trueDecayMom;
	mcDecayVertex.time = trueTimeAtDecay;
	//std::cout<<"true block decayVertex.position: "<<decayVertex.position<<std::endl;
	//std::cout<<"true block decayVertex.momentum: "<<decayVertex.momentum<<std::endl;
      }
      //std::cout<<"before doForwardsExtrap"<<std::endl;
      if (doForwardsExtrap_) { 
	art::Handle<gm2truth::CaloArtRecordCollection> caloDataHandle;
	gm2truth::CaloArtRecordCollection caloData;
	success = e.getByLabel(caloModuleLabel_,caloInstanceLabel_,caloDataHandle);
	auto caloHits = *caloDataHandle;
	for (auto caloHit : caloHits) {
	  if (caloHit.trackID == vertexTrackID) {
	    G4ThreeVector trueCaloPos = {caloHit.x,caloHit.y,caloHit.z};
	    G4ThreeVector trueCaloMom = {caloHit.px,caloHit.py,caloHit.pz};
	    mcForwardsExtrapPoint.position = trueCaloPos;
	    mcForwardsExtrapPoint.momentum = trueCaloMom;
	  }
	}
      }
    }
    // add pointer to the track object to the vertex
    art::Ptr<gm2strawtracker::TrackArtRecord> trackPtr(trackDataHandle,itrack);
    
    // push the decayVertex back to the collection
    if (doBackwardsExtrap_ && backwardSuccess) {
      decayVertex.track  = trackPtr;
      decayVertex.island = trackData[itrack].island;
      if (!e.isRealData()) {
	decayVertex.mcVertex = mcDecayVertex;
      }
      decayVertices->push_back(decayVertex);
      for(uint i=0; i<vertexSteps.steps.size(); i++){
	volumeSteps->push_back(vertexSteps.steps[i]);
	  }
      //std::cout<<"volumeSteps size: "<<volumeSteps->size()<<std::endl;
    }

    if (doForwardsExtrap_  && forwardSuccess) {
      forwardsExtrapPoint.track  = trackPtr;
      forwardsExtrapPoint.island = trackData[itrack].island;
      if (!e.isRealData()) forwardsExtrapPoint.mcVertex = mcForwardsExtrapPoint;
      forwardsExtrapPoints->push_back(forwardsExtrapPoint);
    }
    
  } // loop over tracks

  //std::cout<<"after loop over tracks"<<std::endl;
  mf::LogInfo(name_) << "e.id() = " << e.id() << ", decayVertices->size() = " << decayVertices->size() << ", forwardsExtrapPoints->size() = " << forwardsExtrapPoints->size() << "\n";
  
  // count the number of extrapolated hits
  counter_backwards_ += (int)decayVertices->size();
  counter_forwards_  += (int)forwardsExtrapPoints->size();
  
  // store the collection of extrapolated hits to the event
  e.put(std::move(decayVertices),decayVertexInstanceLabel_);
  e.put(std::move(forwardsExtrapPoints),forwardsExtrapInstanceLabel_);
  e.put(std::move(failedBackwardsExtrapolationEvents),failedBackwardsExtrapInstanceLabel_);
  e.put(std::move(failedForwardsExtrapolationEvents),failedForwardsExtrapInstanceLabel_);
  //std::cout<<"volumeSteps size: "<<volumeSteps->size()<<std::endl;
  e.put(std::move(volumeSteps), "volumeSteps");
  mf::LogInfo(name_) << "Exit GeaneExtrapolation::produce, e.id() = " << e.id() << ", counter_backwards_ = " << counter_backwards_ << ", counter_forwards_ = " << counter_forwards_ << "\n";
  //std::cout<<"ryan module end of produce"<<std::endl;
}


DEFINE_ART_MODULE(gm2strawtracker::GeaneExtrapolation)
