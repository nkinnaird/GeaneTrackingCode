/////////////////////////////////////////////////////////////////////////
// Class:       GeaneExtrapolationAna
// Module Type: analyzer
// File:        GeaneExtrapolationAna_module.cc
// Author:      Ryan McCarthy
////////////////////////////////////////////////////////////////////////

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// data products
#include "gm2dataproducts/strawtracker/DecayVertexArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"
#include "gm2dataproducts/mc/actions/track/TrajectoryArtRecord.hh"
#include "gm2dataproducts/mc/strawtracker/StrawArtRecord.hh"
#include "gm2dataproducts/strawtracker/StrawTimeIslandArtRecord.hh"
#include "gm2dataproducts/mc/calo/CaloArtRecord.hh"
#include "gm2dataproducts/mc/ring/RingArtRecord.hh"
#include "gm2dataproducts/mc/actions/track/TrackingActionArtRecord.hh"

// common util
#include "gm2util/common/dataModuleDefs.hh"
#include "gm2util/common/RootManager.hh"
#include "gm2geom/common/StraightLineTools.hh"
#include "gm2tracker/utils/StrawObjectSorter.hh"

// geometry includes
#include "gm2geom/coordSystems/CoordSystem.hh"
#include "gm2geom/coordSystems/CoordSystemsStoreData.hh"
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh" 
#include "gm2geom/common/Gm2Constants_service.hh" 
#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"

// root 
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TTree.h"
#include "TFile.h"

// c++ includes
#include <typeinfo>
#include <cmath>
#include <math.h>
#include "boost/format.hpp"
#include <iostream>
#include <fstream>


// namespace
using namespace gm2strawtracker;
using gm2geom::CoordSystem3Vector;
using std::vector;
using std::string;

namespace gm2strawtracker {
  class GeaneExtrapolationAna;
}

class gm2strawtracker::GeaneExtrapolationAna : public art::EDAnalyzer {
  
public:
  explicit GeaneExtrapolationAna(fhicl::ParameterSet const & pset);
  virtual ~GeaneExtrapolationAna();
  
  void analyze(art::Event const & e) override;
  void beginRun(art::Run const & r) override;
  void beginJob() override;
  void endJob() override;
  
private:
  // helper functions
  void bookDecayVertexPlots(TDirectory* topDir);
  void bookFailedEventPlots(TDirectory* topDir);
  void bookWhichVolumePlots(TDirectory* topDir);
  void bookRecoCaloHitPlots(TDirectory* topDir);
  void bookTrackPlots(TDirectory* topDir);
  void bookTruthPlots(TDirectory* topDir);
  void bookCaloTruthPlots(TDirectory* topDir);
  void bookGeomGraphs(TDirectory* topDir);
  void makeProfilePlots(TDirectory* topDir);
  
  void fillDecayVertexHistograms( gm2strawtracker::DecayVertexArtRecordCollection const decayVertices, int eventNum, int runNumber, std::string dirName_, std::vector<std::vector<double>>& dataVector, std::vector<ExtrapolationStep> volumeSteps, gm2strawtracker::TrackArtRecordCollection trackStates);
  void fillFailedEventHistograms( gm2strawtracker::DecayVertexArtRecordCollection const failedBackwardsExtrapolationEvents,
				  gm2strawtracker::DecayVertexArtRecordCollection const failedForwardsExtrapolationEvents);
  void fillCaloHitHistograms( gm2strawtracker::DecayVertexArtRecordCollection const extrapolatedCaloHits, int eventNum, int runNumber );
  void fillCaloTruthHistograms( gm2strawtracker::DecayVertexArtRecordCollection const extrapolatedCaloHits, int eventNum, int runNumber );
  void fillTrackHistograms( gm2strawtracker::DecayVertexArtRecordCollection const decayVertices, int eventNum, int runNumber );
  void fillTruthHistograms( gm2strawtracker::DecayVertexArtRecordCollection const decayVertices, int eventNum, int runNumber, std::string dirName_ );
  
  // declare the data product module and instance names
  std::string decayVertexModuleLabel_;
  std::string decayVertexInstanceName_;
  std::string forwardsExtrapModuleLabel_;
  std::string forwardsExtrapInstanceName_;
  std::string failedBackwardsEventModuleLabel_;
  std::string failedBackwardsEventInstanceName_;
  std::string failedForwardsEventModuleLabel_;
  std::string failedForwardsEventInstanceName_;
  std::string csModuleLabel_;
  std::string csInstanceName_;
  
  // helper functions
  gm2geom::CoordSystemsStoreData cs_;
  
  // make truth plots?
  bool makeTruthPlots_;

  // momentum cut
  double p_low_;
  double p_high_;

  //bool
  bool useGeane;
  
  // distance cut
  double distCut_;
  
  // energy loss cut
  bool applyEnergyLossCut_;
  double energyLossCut_;

  // p value cut
  bool cutPoorPValues_;
  double pValueCut_;
  
  // steps?
  bool keepSteps_;
  bool plotSteps_;
  
  // root plotting members
  std::unique_ptr<TFile> outputRootFile_;
  std::unique_ptr<RootManager> rootManager_;
  std::string dirName_;
  
  // name for the messaging service
  std::string name_;
  
  // get the required services
  art::ServiceHandle<gm2geom::Gm2Constants> gm2consts_;

  // which stations to use
  std::set<int> stations_;
  gm2geom::StrawTrackerGeometry geom_;

  // counters
  int nEvt = 0;
  int radialPosOutliers = 0;
  int geaneRadialPosOutliers = 0;
  int radialPosWorldOutliers = 0;
  int geaneRadialPosWorldOutliers = 0;
  int verticalPosOutliers = 0;
  int geaneVerticalPosOutliers = 0;

  // failure modes
  std::map<std::string,int> failureModesBackwards_;
  std::map<std::string,int> failureModesForwards_;

  // volumes hit
  std::map<std::string,int> volumesHitByTracks;
  std::map<std::string,int> volumesHitByNonVacTracks;
  std::map<std::string,int> volumesHitByOutsideWorldTracks;

  std::map<std::string,int> volumesHitByTracksS0;
  std::map<std::string,int> volumesHitByNonVacTracksS0;
  std::map<std::string,int> volumesHitByOutsideWorldTracksS0;

  std::map<std::string,int> volumesHitByTracksS12;
  std::map<std::string,int> volumesHitByNonVacTracksS12;
  std::map<std::string,int> volumesHitByOutsideWorldTracksS12;

  std::map<std::string,int> volumesHitByTracksS18;
  std::map<std::string,int> volumesHitByNonVacTracksS18;
  std::map<std::string,int> volumesHitByOutsideWorldTracksS18;

  std::map<std::string,int> geaneVolumesHitByTracks;
  std::map<std::string,int> geaneVolumesHitByNonVacTracks;
  std::map<std::string,int> geaneVolumesHitByOutsideWorldTracks;

  std::map<std::string,int> geaneVolumesHitByTracksS0;
  std::map<std::string,int> geaneVolumesHitByNonVacTracksS0;
  std::map<std::string,int> geaneVolumesHitByOutsideWorldTracksS0;

  std::map<std::string,int> geaneVolumesHitByTracksS12;
  std::map<std::string,int> geaneVolumesHitByNonVacTracksS12;
  std::map<std::string,int> geaneVolumesHitByOutsideWorldTracksS12;

  std::map<std::string,int> geaneVolumesHitByTracksS18;
  std::map<std::string,int> geaneVolumesHitByNonVacTracksS18;
  std::map<std::string,int> geaneVolumesHitByOutsideWorldTracksS18;

  // constants
  double mm2m  = 0.001;
  double m2mm  = 100.0;
  double cm2m  = 0.1;
  double m2cm  = 10.0;
  double mm2cm = 0.1;
  double cm2mm = 10.0;

  // angle around ring of face of calos 12 and 18
  double caloTheta_0  = 0.20348; // TODO - check this number, currently just caloTheta_12 - (12 * ( (caloTheta_18 - caloTheta_12) / 6) )
  double caloTheta_12 = 3.34508; // rad
  double caloTheta_18 = 4.91588; // rad

  struct vertexData_t{
    double radialPosition;
    double verticalPosition;
    double radialMomentum;
    double verticalMomentum;
    double theta;
    double decayArcLength;
    double time;
    double totalMomentum;
    double worldRadialPos;
  };

  struct truthData_t{
    double trueRadialPos;
    double trueVerticalPos;
    double trueRadialMom;
    double trueVerticalMom;
    double trueTheta;
    double trueArcLength;
    double trueTotalMomentum;
    double trueWorldRadialPos;
  };
    
  struct diffData_t{
    double diffRadialPosition;
    double diffVerticalPosition;
    double diffRadialMomentum;
    double diffVerticalMomentum;
    double diffTheta;
    double diffDecayArcLength;
    double diffTime;
    double diffTotalMomentum;
    double diffWorldRadialPos;
  };
  
  vertexData_t vertexData;
  truthData_t truthData;
  diffData_t diffData;
  TTree* saskiaTree_;
  TTree* geaneTree_;
  TTree* truthTree_;
  TTree* diffTree_;


};


// standard constructor
gm2strawtracker::GeaneExtrapolationAna::GeaneExtrapolationAna(fhicl::ParameterSet const & pset) 
  : EDAnalyzer(pset)
  , decayVertexModuleLabel_( pset.get<std::string>("decayVertexModuleLabel",dataModuleDefs::vertexModuleLabel()) )
  , decayVertexInstanceName_( pset.get<std::string>("decayVertexInstanceName",dataModuleDefs::backExtrapInstanceLabel()) )
  , forwardsExtrapModuleLabel_( pset.get<std::string>( "forwardsExtrapModuleLabel",dataModuleDefs::vertexModuleLabel()) )
  , forwardsExtrapInstanceName_( pset.get<std::string>("forwardsExtrapInstanceName",dataModuleDefs::forwardExtrapInstanceLabel()) )
  , failedBackwardsEventModuleLabel_( pset.get<std::string>("failedBackwardsEventModuleLabel",dataModuleDefs::vertexModuleLabel()) )
  , failedBackwardsEventInstanceName_( pset.get<std::string>("failedBackwardsEventInstanceName","failedBackwardsExtrapolation") )
  , failedForwardsEventModuleLabel_( pset.get<std::string>("failedEventModuleLabel",dataModuleDefs::vertexModuleLabel()) )
  , failedForwardsEventInstanceName_( pset.get<std::string>("failedEventInstanceName","failedForwardsExtrapolation") )
  , csModuleLabel_( pset.get<std::string>("csModuleLabel",dataModuleDefs::coordSysModuleLabel()) )
  , csInstanceName_( pset.get<std::string>("csInstanceName",dataModuleDefs::coordSysInstanceLabel()) )
  , cs_()
  , makeTruthPlots_(pset.get<bool>("makeTruthPlots",true))
  , p_low_(pset.get<float>("lowMomentum",0.0) )
  , p_high_(pset.get<float>("highMomentum",2500.0) ) // MeV
  , distCut_(pset.get<float>("maxDistance",10000.0) ) // mm
  , applyEnergyLossCut_(pset.get<bool>("applyEnergyLossCut",false) )
  , energyLossCut_(pset.get<double>("energyLossCut",0.0) ) //MeV
  , cutPoorPValues_(pset.get<bool>("cutPoorPValues",false))
  , pValueCut_(pset.get<double>("pValueCut",0.005))
  , keepSteps_(pset.get<bool>("keepSteps",false))
  , plotSteps_(pset.get<bool>("plotSteps",false))
  , dirName_(pset.get<std::string>("dirName","Extrapolation"))
  , name_( "GeaneExtrapolationAna" )
  , geom_()
  , saskiaTree_()
  , geaneTree_()
  , truthTree_()
  , diffTree_()

{

  //better to load the station up from the geometry
  for(auto s : geom_.whichScallopLocations) {
    stations_.insert(s);
  }

}

//! standard destructor
gm2strawtracker::GeaneExtrapolationAna::~GeaneExtrapolationAna()
{
  //! Clean up dynamic memory and other resources here.
}

void gm2strawtracker::GeaneExtrapolationAna::beginRun(art::Run const & r)
{
  // get the coordinate system tools
  cs_ = artg4::dataFromRunOrService<gm2geom::CoordSystemsStoreData, gm2geom::CoordSystemsStore>(r,csModuleLabel_,csInstanceName_);
  if( cs_.size() == 0 )
    throw cet::exception("GeaneExtrapolationAna") << "This run does not contain any data associated with the coordinate system\n";
}


void gm2strawtracker::GeaneExtrapolationAna::analyze(art::Event const & e)
{
  mf::LogInfo info(name_);
  info << "Event= " << e.id() << " START \n";
  info << "e.run() = " << e.run() << "\n";
  
  int runNum = e.run();


  std::vector<double> saskiaRadialPositionVector = {};
  std::vector<double> saskiaVerticalPositionVector = {};
  std::vector<double> saskiaRadialMomentumVector = {};
  std::vector<double> saskiaVerticalMomentumVector = {};
  std::vector<double> saskiaThetaVector = {};
  std::vector<double> saskiaDecayArcLengthVector = {};
  std::vector<double> saskiaTimeVector = {};
  std::vector<double> saskiaTotalMomentumVector = {};
  std::vector<double> saskiaWorldRadialPosVector = {};
  std::vector<std::vector<double>> saskiaDataVector = {};
  saskiaDataVector.push_back(saskiaRadialPositionVector);
  saskiaDataVector.push_back(saskiaVerticalPositionVector);
  saskiaDataVector.push_back(saskiaRadialMomentumVector);
  saskiaDataVector.push_back(saskiaVerticalMomentumVector);
  saskiaDataVector.push_back(saskiaThetaVector);
  saskiaDataVector.push_back(saskiaDecayArcLengthVector);
  saskiaDataVector.push_back(saskiaTimeVector);
  saskiaDataVector.push_back(saskiaTotalMomentumVector);
  saskiaDataVector.push_back(saskiaWorldRadialPosVector);
 

  std::vector<double> ryanRadialPositionVector = {};
  std::vector<double> ryanVerticalPositionVector = {};
  std::vector<double> ryanRadialMomentumVector = {};
  std::vector<double> ryanVerticalMomentumVector = {};
  std::vector<double> ryanThetaVector = {};
  std::vector<double> ryanDecayArcLengthVector = {};
  std::vector<double> ryanTimeVector = {};
  std::vector<double> ryanTotalMomentumVector = {};
  std::vector<double> ryanWorldRadialPosVector = {};
  std::vector<std::vector<double>> ryanDataVector = {};
  ryanDataVector.push_back(ryanRadialPositionVector);
  ryanDataVector.push_back(ryanVerticalPositionVector);
  ryanDataVector.push_back(ryanRadialMomentumVector);
  ryanDataVector.push_back(ryanVerticalMomentumVector);
  ryanDataVector.push_back(ryanThetaVector);
  ryanDataVector.push_back(ryanDecayArcLengthVector);
  ryanDataVector.push_back(ryanTimeVector);
  ryanDataVector.push_back(ryanTotalMomentumVector);
  ryanDataVector.push_back(ryanWorldRadialPosVector);

  // success condition
  bool success1 = false;
  bool success2 = false;
  bool success3 = false;
  bool success4 = false;

  //! get decay vertices
  art::Handle<gm2strawtracker::DecayVertexArtRecordCollection> decayVertexDataHandle;
  success1 = e.getByLabel(decayVertexModuleLabel_,decayVertexInstanceName_,decayVertexDataHandle);

  std::cout<<"success1: "<<success1<<std::endl;
  
  art::Handle<gm2strawtracker::DecayVertexArtRecordCollection> geaneDecayVertexDataHandle;
  bool geaneSuccess1 = false;
  std::string const geaneModuleLabel = "geaneExtrapolater";
  std::string const geaneInstanceLabel = "backwards";
  geaneSuccess1 = e.getByLabel(geaneModuleLabel,geaneInstanceLabel,geaneDecayVertexDataHandle);

  std::cout<<"geaneSuccess1: "<<geaneSuccess1<<std::endl;

  art::Handle<gm2strawtracker::TrackArtRecordCollection> nickTrackStateDataHandle;
  bool nickSuccess1 = false;
  std::string const nickModuleLabel = "tracks";
  std::string const nickInstanceLabel = "straws";
  nickSuccess1 = e.getByLabel(nickModuleLabel,nickInstanceLabel,nickTrackStateDataHandle);
  std::cout<<"nickSuccess1: "<<nickSuccess1<<std::endl;
  
  art::Handle<std::vector<ExtrapolationStep>> volumeStepsHandle;
  bool geaneSuccess2 = false;
  std::string const stepsInstanceLabel = "volumeSteps";
  geaneSuccess2 = e.getByLabel(geaneModuleLabel,stepsInstanceLabel,volumeStepsHandle);
  std::cout<<"geaneSuccess2: "<<geaneSuccess2<<std::endl;

  //! get forwards extrap record
  art::Handle<gm2strawtracker::DecayVertexArtRecordCollection> forwardsExtrapDataHandle;
  success2 = e.getByLabel(forwardsExtrapModuleLabel_,forwardsExtrapInstanceName_,forwardsExtrapDataHandle);

  //! get events that fail the extrapolation
  art::Handle<gm2strawtracker::DecayVertexArtRecordCollection> failedBackwardsEventDataHandle;
  success3 = e.getByLabel(failedBackwardsEventModuleLabel_,failedBackwardsEventInstanceName_,failedBackwardsEventDataHandle);

  art::Handle<gm2strawtracker::DecayVertexArtRecordCollection> failedForwardsEventDataHandle;
  success4 = e.getByLabel(failedForwardsEventModuleLabel_,failedForwardsEventInstanceName_,failedForwardsEventDataHandle);

  if (!success1 && !success2 && !success3 && !success4) {
    mf::LogWarning(name_) << "Couldn't get any of the required art records" << "\n";
    return;
  }

  //! get the data
  gm2strawtracker::DecayVertexArtRecordCollection decayVertexData;
  gm2strawtracker::DecayVertexArtRecordCollection geaneDecayVertexData;
  gm2strawtracker::DecayVertexArtRecordCollection forwardsExtrapData;
  gm2strawtracker::DecayVertexArtRecordCollection failedBackwardsEvents;
  gm2strawtracker::DecayVertexArtRecordCollection failedForwardsEvents;
  std::vector<ExtrapolationStep> volumeSteps;
  gm2strawtracker::TrackArtRecordCollection nickTrackStateData;
  std::vector<gm2strawtracker::TrackStates> nickTrackStateCollection;

  if (success1) decayVertexData = *decayVertexDataHandle;
  if (geaneSuccess1) geaneDecayVertexData = *geaneDecayVertexDataHandle;
  if (nickSuccess1){
    nickTrackStateData = *nickTrackStateDataHandle;
    for( uint i=0; i<nickTrackStateData.size(); i++){
      nickTrackStateCollection.push_back(nickTrackStateData[i].states);
    }
  }					 
  if (geaneDecayVertexData.size() != 0){
    //std::cout<<"Geane decayVertex.position: "<<geaneDecayVertexData[0].position<<std::endl;
    //std::cout<<"Geane decayVertex.momentum: "<<geaneDecayVertexData[0].momentum<<std::endl;
  }

  std::cout<<"decayVertexData size: "<<decayVertexData.size()<<std::endl;
  std::cout<<"geaneDecayVertexData size: "<<geaneDecayVertexData.size()<<std::endl;
  //std::cout<<"saskia first vertex params: "<<decayVertexData[0].startMomentum<<" "<<decayVertexData[0].startPosition<<" "<<decayVertexData[0].t0<<" "<<decayVertexData[0].position<<std::endl;
  //std::cout<<"ryan first vertex params: "<<geaneDecayVertexData[0].startMomentum<<" "<<geaneDecayVertexData[0].startPosition<<" "<<geaneDecayVertexData[0].t0<<" "<<geaneDecayVertexData[0].position<<std::endl;
  
  if (geaneSuccess2) volumeSteps = *volumeStepsHandle;
  if (volumeSteps.size() != 0){
    //std::cout<<"volumeSteps size: "<<volumeSteps.size()<<std::endl;
  }

  if (success2) forwardsExtrapData = *forwardsExtrapDataHandle;
  if (success3) failedBackwardsEvents = *failedBackwardsEventDataHandle;
  if (success4) failedForwardsEvents  = *failedForwardsEventDataHandle;

  if (success1 && decayVertexData.size() == 0) return;
  if (success2 && forwardsExtrapData.size() == 0) return;

  if (!decayVertexData.empty()) {
    std::vector<ExtrapolationStep> dummyVolumeSteps;
    fillDecayVertexHistograms(decayVertexData, nEvt, runNum, "Extrapolation", saskiaDataVector, dummyVolumeSteps, nickTrackStateData);
    fillDecayVertexHistograms(geaneDecayVertexData, nEvt, runNum, "GeaneExtrapolation", ryanDataVector, volumeSteps, nickTrackStateData);
    fillTrackHistograms(decayVertexData, nEvt, runNum);
    if (!e.isRealData()){
      fillTruthHistograms(decayVertexData,nEvt,runNum, "Extrapolation");
      fillTruthHistograms(geaneDecayVertexData,nEvt,runNum, "GeaneExtrapolation");
    }
  }
  std::cout<<"After decayVertexData.empty"<<std::endl;
  if (!forwardsExtrapData.empty()) {
    fillCaloHitHistograms(forwardsExtrapData, nEvt, runNum);
    if (!e.isRealData()) fillCaloTruthHistograms(forwardsExtrapData, nEvt,runNum);
  }
  
  if ( ! (failedBackwardsEvents.empty() && failedForwardsEvents.empty()) ) fillFailedEventHistograms(failedBackwardsEvents,failedForwardsEvents);
  
  nEvt++;

  //std::cout<<"saskiaDataVector[0]: "<<saskiaDataVector[0].size()<<std::endl;
  //std::cout<<"saskiaDataVector[4]: "<<saskiaDataVector[4].size()<<std::endl;
  //std::cout<<"ryanDataVector[0]: "<<ryanDataVector[0].size()<<std::endl;
  //std::cout<<"ryanDataVector[4]: "<<ryanDataVector[4].size()<<std::endl;

  if(saskiaDataVector[0].size() != 0){
    if(ryanDataVector[0].size() != 0){
      std::cout<<"saskiaDataVector sizes: "<<saskiaDataVector[0].size()<< " " << saskiaDataVector[1].size() << " " << saskiaDataVector[2].size() << " " <<saskiaDataVector[3].size() << " " << saskiaDataVector[4].size() << " " << saskiaDataVector[5].size() << " " << saskiaDataVector[6].size() << " " << saskiaDataVector[7].size() << " " << saskiaDataVector[8].size();
      std::cout<<"ryanDataVector sizes: "<<ryanDataVector[0].size()<< " " << ryanDataVector[1].size() << " " << ryanDataVector[2].size() << " " <<ryanDataVector[3].size() << " " << ryanDataVector[4].size() << " " << ryanDataVector[5].size() << " " << ryanDataVector[6].size() << " " << ryanDataVector[7].size() << " " << ryanDataVector[8].size();

      for(unsigned int i=0; i < saskiaDataVector[0].size(); i++){
	std::cout<<"saskia radialPos: "<<saskiaDataVector[0][i]<<std::endl;
	std::cout<<"ryan radialPos: "<<ryanDataVector[0][i]<<std::endl;
	diffData.diffRadialPosition = saskiaDataVector[0][i] - ryanDataVector[0][i];
	diffData.diffVerticalPosition = saskiaDataVector[1][i] - ryanDataVector[1][i];
	diffData.diffRadialMomentum = saskiaDataVector[2][i] - ryanDataVector[2][i];
	diffData.diffVerticalMomentum = saskiaDataVector[3][i] - ryanDataVector[3][i];
	std::cout<<"verticalMomDiff: "<<diffData.diffVerticalMomentum<<std::endl;
	diffData.diffTheta = saskiaDataVector[4][i] - ryanDataVector[4][i];
	diffData.diffDecayArcLength = saskiaDataVector[5][i] - ryanDataVector[5][i];
	diffData.diffTime = saskiaDataVector[6][i] - ryanDataVector[6][i];
	diffData.diffTotalMomentum = saskiaDataVector[7][i] - ryanDataVector[7][i];
	diffData.diffWorldRadialPos = saskiaDataVector[8][i] - ryanDataVector[8][i];
	
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialPosDiff")->Fill(diffData.diffRadialPosition);
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialMomDiff")->Fill(diffData.diffRadialMomentum);
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialMomDiffZoom")->Fill(diffData.diffRadialMomentum);
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialMomDiffOriginZoom")->Fill(diffData.diffRadialMomentum);
	rootManager_->Get<TH2*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialPosDiff_vs_radialMomDiff")->Fill(diffData.diffRadialMomentum, diffData.diffRadialPosition);
	rootManager_->Get<TH2*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialPosDiff_vs_radialMomDiff")->SetOption("COLZ");
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_verticalPosDiff")->Fill(diffData.diffVerticalPosition);
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_verticalMomDiff")->Fill(diffData.diffVerticalMomentum);

	diffTree_->Fill();
      }
    }
  }

#if 0
  bool saskiaIsIncreasing = true;
  bool ryanIsIncreasing = true;
  uint j = 0;
  DecayVertexArtRecord saskiaCurrentVertex;
  DecayVertexArtRecord ryanCurrentVertex;

  for (unsigned int i=0; i < decayVertexData.size(); i++){
    if((i>0) && (j>0)){
      int saskiaT0Diff = decayVertexData[i].t0 - decayVertexData[i-1].t0;
      if(saskiaT0Diff>0){
	saskiaIsIncreasing = true;
      }
      else{
	saskiaIsIncreasing = false;
      }
      int ryanT0Diff = geaneDecayVertexData[j].t0 - geaneDecayVertexData[j-1].t0;
      if(ryanT0Diff>0){
	ryanIsIncreasing = true;
      }
      else{
	ryanIsIncreasing = false;
      }
      if((saskiaIsIncreasing && ryanIsIncreasing) || (!saskiaIsIncreasing && !ryanIsIncreasing)){
	if(decayVertexData[i].t0 == geaneDecayVertexData[j].t0){
	  saskiaCurrentVertex = decayVertexData[i];
	  ryanCurrentVertex = geaneDecayVertexData[j];
	  if((saskiaCurrentVertex.startMomentum==ryanCurrentVertex.startMomentum)&&(saskiaCurrentVertex.startPosition==ryanCurrentVertex.startPosition)){
	    std::cout<<"matches: true"<<std::endl;
	  }
	  else{
	    std::cout<<"matches: false"<<std::endl;
	  }
	}
	else if(decayVertexData[i].t0 > geaneDecayVertexData[j].t0){
	  i--;
	}
	else{
	  j--;
	}	  
      }
      else if(saskiaIsIncreasing & !ryanIsIncreasing){
	if(decayVertexData[i].t0 == geaneDecayVertexData[j].t0){
	  saskiaCurrentVertex = decayVertexData[i];
	  ryanCurrentVertex = geaneDecayVertexData[j];
	  if((saskiaCurrentVertex.startMomentum==ryanCurrentVertex.startMomentum)&&(saskiaCurrentVertex.startPosition==ryanCurrentVertex.startPosition)){
	    std::cout<<"matches: true"<<std::endl;
	  }
	  else{
	    std::cout<<"matches: false"<<std::endl;
	  }
	}
	else if(decayVertexData[i].t0 > geaneDecayVertexData[j].t0){
	  j--;
	}
	else{
	  j--;
	}
      }
      else{
	if(decayVertexData[i].t0 == geaneDecayVertexData[j].t0){
	  saskiaCurrentVertex = decayVertexData[i];
	  ryanCurrentVertex = geaneDecayVertexData[j];
	  if((saskiaCurrentVertex.startMomentum==ryanCurrentVertex.startMomentum)&&(saskiaCurrentVertex.startPosition==ryanCurrentVertex.startPosition)){
	    std::cout<<"matches: true"<<std::endl;
	  }
	  else{
	    std::cout<<"matches: false"<<std::endl;
	  }
	}
	else if(decayVertexData[i].t0 > geaneDecayVertexData[j].t0){
	  i--;
	}
	else{
	  i--;
	}
      }

    }//if i and j above 0
    else{
      if(decayVertexData[i].t0 == geaneDecayVertexData[j].t0){
	saskiaCurrentVertex = decayVertexData[i];
	ryanCurrentVertex = geaneDecayVertexData[j];
	if((saskiaCurrentVertex.startMomentum==ryanCurrentVertex.startMomentum)&&(saskiaCurrentVertex.startPosition==ryanCurrentVertex.startPosition)){
	    std::cout<<"matches: true"<<std::endl;
	  }
	  else{
	    std::cout<<"matches: false"<<std::endl;
	  }
      }
      else if(decayVertexData[i].t0 > geaneDecayVertexData[j].t0){
	i--;
      }
      else{
	j--;
      }
    }

    
    j++;
  }//loop over decayVertices
#endif
#if 0 
  for( auto saskiaDecayVertex : decayVertexData ){
    saskiaDecayVertex.t0 = saskiaDecayVertex.track->geanePlanes->candidate->t0;
    saskiaDecayVertex.startPosition = saskiaDecayVertex.track->states.at(0)->position;
    saskiaDecayVertex.startMomentum = saskiaDecayVertex.track->momentum;
  }

  for( auto ryanDecayVertex : geaneDecayVertexData ){
    ryanDecayVertex.t0 = ryanDecayVertex.track->geanePlanes->candidate->t0;
    ryanDecayVertex.startPosition = ryanDecayVertex.track->states.at(0)->position;
    ryanDecayVertex.startMomentum = ryanDecayVertex.track->momentum;
  }
#endif
  std::cout<<"Before matches"<<std::endl;
  int matches = 0;
  std::cout<<"decayVertexData size: "<<decayVertexData.size()<<std::endl;
  std::cout<<"geaneDecayVertexData size: "<<geaneDecayVertexData.size()<<std::endl;
  for( auto saskiaDecayVertex : decayVertexData ){
    for( auto ryanDecayVertex : geaneDecayVertexData ){
      std::cout<<"decayVertexData cs: "<<saskiaDecayVertex.position.coordSystemName<<std::endl;
      std::cout<<"geaneDecayVertexData cs: "<<ryanDecayVertex.position.coordSystemName<<std::endl;
      std::cout<<"decayVertexData start cs: "<<saskiaDecayVertex.startPosition.coordSystemName<<std::endl;
      std::cout<<"geaneDecayVertexData start cs: "<<ryanDecayVertex.startPosition.coordSystemName<<std::endl;
      if((saskiaDecayVertex.startMomentum==ryanDecayVertex.startMomentum)&&(saskiaDecayVertex.startPosition==ryanDecayVertex.startPosition)&&(saskiaDecayVertex.t0==ryanDecayVertex.t0)){
	matches++;
	std::cout<<"matches: "<<matches<<std::endl;
	std::cout<<"saskiaDecayVertex startMomentum: "<<saskiaDecayVertex.startMomentum<<std::endl;
	std::cout<<"saskiaDecayVertex startPosition: "<<saskiaDecayVertex.startPosition<<std::endl;
	std::cout<<"saskiaDecayVertex t0: "<<saskiaDecayVertex.t0<<std::endl;
	std::cout<<"ryanDecayVertex startMomentum: "<<ryanDecayVertex.startMomentum<<std::endl;
	std::cout<<"ryanDecayVertex startPosition: "<<ryanDecayVertex.startPosition<<std::endl;
	std::cout<<"ryanDecayVertex t0: "<<ryanDecayVertex.t0<<std::endl;
	
	G4ThreeVector ryanTangentPointPos{ryanDecayVertex.position.x()/mm, ryanDecayVertex.position.y()/mm, ryanDecayVertex.position.z()/mm};
	G4ThreeVector ryanTangentPointMom{-ryanDecayVertex.momentum.x()/MeV, -ryanDecayVertex.momentum.y()/MeV, -ryanDecayVertex.momentum.z()/MeV};
	G4ThreeVector saskiaTangentPointPos{saskiaDecayVertex.position.x()/mm, saskiaDecayVertex.position.y()/mm, saskiaDecayVertex.position.z()/mm};
	G4ThreeVector saskiaTangentPointMom{saskiaDecayVertex.momentum.x()/MeV, saskiaDecayVertex.momentum.y()/MeV, saskiaDecayVertex.momentum.z()/MeV};
	G4ThreeVector ryanUncorrectedTangentPointPos{ryanDecayVertex.uncorrectedPosition.x(), ryanDecayVertex.uncorrectedPosition.y(), ryanDecayVertex.uncorrectedPosition.z()};
	G4ThreeVector ryanUncorrectedTangentPointMom{ryanDecayVertex.uncorrectedPosition.x(), ryanDecayVertex.uncorrectedPosition.y(), ryanDecayVertex.uncorrectedPosition.z()};
	std::cout<<"ryanTangentPointMom: "<<ryanTangentPointMom<<std::endl;
	std::cout<<"saskiaTangentPointMom: "<<saskiaTangentPointMom<<std::endl;

	double ryanTangentPointRadialPos = gm2consts_->ComputeRhat( &ryanTangentPointPos );
	double ryanUncorrectedTangentPointRadialPos = gm2consts_->ComputeRhat( &ryanUncorrectedTangentPointPos );
	double saskiaTangentPointRadialPos = gm2consts_->ComputeRhat( &saskiaTangentPointPos );
	double ryanWorldRadialPos = sqrt(pow(ryanTangentPointPos[0],2)+pow(ryanTangentPointPos[2],2));
	double saskiaWorldRadialPos = sqrt(pow(saskiaTangentPointPos[0],2)+pow(saskiaTangentPointPos[2],2));
	double ryanTangentPointVerticalPos = gm2consts_->ComputeVhat( &ryanTangentPointPos );
	double ryanUncorrectedTangentPointVerticalPos = gm2consts_->ComputeVhat( &ryanUncorrectedTangentPointPos );
	double saskiaTangentPointVerticalPos = gm2consts_->ComputeVhat( &saskiaTangentPointPos );
	double ryanTangentPointRadialMom = gm2consts_->ComputePrhat( &ryanTangentPointPos, &ryanTangentPointMom );
	double ryanUncorrectedTangentPointRadialMom = gm2consts_->ComputePrhat( &ryanUncorrectedTangentPointPos, &ryanUncorrectedTangentPointMom );
	double saskiaTangentPointRadialMom = gm2consts_->ComputePrhat( &saskiaTangentPointPos, &saskiaTangentPointMom );
	double ryanTangentPointVerticalMom = gm2consts_->ComputePvhat( &ryanTangentPointPos, &ryanTangentPointMom) * ryanTangentPointMom.mag();
	double ryanUncorrectedTangentPointVerticalMom = gm2consts_->ComputePvhat( &ryanUncorrectedTangentPointPos, &ryanUncorrectedTangentPointMom) * ryanUncorrectedTangentPointMom.mag();
	double saskiaTangentPointVerticalMom = gm2consts_->ComputePvhat( &saskiaTangentPointPos, &saskiaTangentPointMom) * saskiaTangentPointMom.mag();
	double ryanTangentPointTheta = gm2consts_->ComputeTheta( &ryanTangentPointPos );
	double saskiaTangentPointTheta = gm2consts_->ComputeTheta( &saskiaTangentPointPos );
	
	int worseRadialPos = 0;
	int worseRadialMom = 0;
	int worseVerticalPos = 0;
	int worseVerticalMom = 0;

	std::cout<<"ryanDecayVertex.uncorrectedSteps size: "<<ryanDecayVertex.uncorrectedSteps.size()<<std::endl;

	if( abs(saskiaTangentPointRadialPos - ryanUncorrectedTangentPointRadialPos) < abs(saskiaTangentPointRadialPos - ryanTangentPointRadialPos) ){
	  worseRadialPos++;
	}
	
	if( abs(saskiaTangentPointRadialMom - ryanUncorrectedTangentPointRadialMom) < abs(saskiaTangentPointRadialMom - ryanTangentPointRadialMom) ){
	  worseRadialMom++;
	}

	if( abs(saskiaTangentPointVerticalPos - ryanUncorrectedTangentPointVerticalPos) < abs(saskiaTangentPointVerticalPos - ryanTangentPointVerticalPos) ){
	  worseVerticalPos++;
	}

	if( abs(saskiaTangentPointVerticalMom - ryanUncorrectedTangentPointVerticalMom) < abs(saskiaTangentPointVerticalMom - ryanTangentPointVerticalMom) ){
	  worseVerticalMom++;
	}

	std::cout<<"worseRadialPos: "<<worseRadialPos<<std::endl;
	std::cout<<"worseRadialMom: "<<worseRadialMom<<std::endl;
	std::cout<<"worseVerticalPos: "<<worseVerticalPos<<std::endl;
	std::cout<<"worseVerticalMom: "<<worseVerticalMom<<std::endl;

	for( auto step : ryanDecayVertex.steps ){
	  std::cout<<"step position: "<<step.position<<std::endl;
	  std::cout<<"step momentum: "<<step.momentum<<std::endl<<std::endl;
	}
	for( auto step : ryanDecayVertex.uncorrectedSteps ){
	  std::cout<<"uncorrectedStep position: "<<step.position<<std::endl;
	  std::cout<<"uncorrectedStep momentum: "<<step.momentum<<std::endl<<std::endl;
	}
	for( auto step : saskiaDecayVertex.steps ){
	  std::cout<<"saskia step position: "<<step.position<<std::endl;
	  std::cout<<"saskia step momentum: "<<step.momentum<<std::endl<<std::endl;
	}

	int ryanStationNum = ryanDecayVertex.island->station;
	int saskiaStationNum = saskiaDecayVertex.island->station;
	double ryanThetaDiff_0 = caloTheta_0 - ryanTangentPointTheta;
	double saskiaThetaDiff_0 = caloTheta_0 - saskiaTangentPointTheta;
	double ryanThetaDiff_12 = caloTheta_12 - ryanTangentPointTheta;
	double saskiaThetaDiff_12 = caloTheta_12 - saskiaTangentPointTheta;
	double ryanThetaDiff_18 = caloTheta_18 - ryanTangentPointTheta;
	double saskiaThetaDiff_18 = caloTheta_18 - saskiaTangentPointTheta;
	double ryanThetaDiff(0.0);
	if ( ryanStationNum == 0 ) ryanThetaDiff = ryanThetaDiff_0; 
	if ( ryanStationNum == 12 ) ryanThetaDiff = ryanThetaDiff_12; 
	if ( ryanStationNum == 18 ) ryanThetaDiff = ryanThetaDiff_18; 
	double saskiaThetaDiff(0.0);
	if (saskiaStationNum == 0 ) saskiaThetaDiff = saskiaThetaDiff_0;
	if (saskiaStationNum == 12 ) saskiaThetaDiff = saskiaThetaDiff_12;
	if (saskiaStationNum == 18 ) saskiaThetaDiff = saskiaThetaDiff_18;
	
	double ryanDecayArcLength = ryanThetaDiff * (gm2consts_->R_magic() + ryanTangentPointRadialPos) * mm2m;
	double saskiaDecayArcLength = saskiaThetaDiff * (gm2consts_->R_magic() + saskiaTangentPointRadialPos) * mm2m;
	auto ryanTime = ryanDecayVertex.time * 1e-3;
	auto saskiaTime = saskiaDecayVertex.time * 1e-3;
	double ryanTotalMomentum = ryanDecayVertex.momentum.mag();
	double saskiaTotalMomentum = saskiaDecayVertex.momentum.mag();

	diffData.diffRadialPosition = saskiaTangentPointRadialPos - ryanTangentPointRadialPos;
	diffData.diffVerticalPosition = saskiaTangentPointVerticalPos - ryanTangentPointVerticalPos;
	diffData.diffRadialMomentum = saskiaTangentPointRadialMom - ryanTangentPointRadialMom;
	diffData.diffVerticalMomentum = saskiaTangentPointVerticalMom + ryanTangentPointVerticalMom;
	diffData.diffTheta = saskiaTangentPointTheta - ryanTangentPointTheta;
	diffData.diffDecayArcLength = saskiaDecayArcLength - ryanDecayArcLength;
	diffData.diffTime = saskiaTime - ryanTime;
	diffData.diffTotalMomentum = saskiaTotalMomentum - ryanTotalMomentum;
	diffData.diffWorldRadialPos = saskiaWorldRadialPos - ryanWorldRadialPos;
	
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialPosDiff")->Fill(diffData.diffRadialPosition);
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialMomDiff")->Fill(diffData.diffRadialMomentum);
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialMomDiffZoom")->Fill(diffData.diffRadialMomentum);
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialMomDiffOriginZoom")->Fill(diffData.diffRadialMomentum);
	rootManager_->Get<TH2*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialPosDiff_vs_radialMomDiff")->Fill(diffData.diffRadialMomentum, diffData.diffRadialPosition);
	rootManager_->Get<TH2*>("GeaneExtrapolation/vertices/allStations/allEvents","h_radialPosDiff_vs_radialMomDiff")->SetOption("COLZ");
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_verticalPosDiff")->Fill(diffData.diffVerticalPosition);
	rootManager_->Get<TH1*>("GeaneExtrapolation/vertices/allStations/allEvents","h_verticalMomDiff")->Fill(diffData.diffVerticalMomentum);
	rootManager_->Get<TH2*>("GeaneExtrapolation/vertices/allStations/allEvents","h_verticalPosDiff_vs_verticalMomDiff")->Fill(diffData.diffVerticalMomentum, diffData.diffVerticalPosition);
	rootManager_->Get<TH2*>("GeaneExtrapolation/vertices/allStations/allEvents","h_verticalPosDiff_vs_radialPosDiff")->Fill(diffData.diffRadialPosition, diffData.diffVerticalPosition);

	diffTree_->Fill();
      }
    }
  }
  
     
  


} // analyze

void gm2strawtracker::GeaneExtrapolationAna::fillFailedEventHistograms(  gm2strawtracker::DecayVertexArtRecordCollection const failedBackwardsEvents
									,gm2strawtracker::DecayVertexArtRecordCollection const failedForwardsEvents)
{
  if (!failedBackwardsEvents.empty()) {
    for (auto failedEvent : failedBackwardsEvents) {
      failureModesBackwards_[failedEvent.failureMode]++;
    }
  }

  if (!failedForwardsEvents.empty()) {
    for (auto failedEvent : failedForwardsEvents) {
      failureModesForwards_[failedEvent.failureMode]++;
    }
  }
}


void gm2strawtracker::GeaneExtrapolationAna::fillDecayVertexHistograms( gm2strawtracker::DecayVertexArtRecordCollection const decayVertices,
									int eventNum,
									int runNumber, std::string dirName_, std::vector<std::vector<double>>& dataVector, std::vector<ExtrapolationStep> volumeSteps, gm2strawtracker::TrackArtRecordCollection trackStates)
{

  mf::LogInfo info(name_);
  mf::LogDebug debug(name_);
  
  info << "Event number " << eventNum << " has " << decayVertices.size() << " vertices\n";
  
  // check there are vertices and trajectories
  if ( decayVertices.empty() ) return;
  
  if (dirName_ == "GeaneExtrapolation"){
    TDirectory* topDir = rootManager_->GetDir(dirName_,true);
    TDirectory* whichVolumes = rootManager_->GetDir(topDir,"whichVolumes",true);

    for ( auto decayVertex : decayVertices ){
      std::vector<gm2geom::CoordSystem3Vector> impactVectors = decayVertex.impactVectors;
      for ( auto impactVector : impactVectors ){
	rootManager_->Get<TH2*>(whichVolumes,"h_impactVectors_2D") ->Fill(impactVector[0], impactVector[2]);
      }
    }

    int count = 0;
    for(auto state: trackStates){
      count++;
      if(count<200){
	rootManager_->Get<TH3*>(whichVolumes,"h_world_3D") ->Fill(state.position[0],state.position[2],state.position[1]);
      }
      if(state.position[0] < -4000){
	rootManager_->Get<TH3*>(whichVolumes,"h_vacuumChamberCadMesh_3D") ->Fill(state.position[0],state.position[2],state.position[1]);
      }
      else{
	rootManager_->Get<TH3*>(whichVolumes,"h_vacuumChamberCadMesh_3D_2") ->Fill(state.position[0],state.position[2],state.position[1]);
      }
    }
    for (auto volumeStep : volumeSteps) {
#if 0
      double hitTheta = atan(volumeStep.position[2]/volumeStep.position[0]);
      double hitRadius = sqrt(pow(volumeStep.position[0],2)+pow(volumeStep.position[2],2));
      G4ThreeVector hitPosition(volumeStep.position[0], volumeStep.position[1], volumeStep.position[2]);
      double hitBeamRadius = gm2consts_->ComputeRhat( &hitPosition );
      rootManager_->Get<TH3*>(whichVolumes,"h_world_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);

      //std::cout<<"volumeStep.volume: "<<volumeStep.volume<<std::endl;
      if (volumeStep.volume == "supportPostLVShell"){
	rootManager_->Get<TH3*>(whichVolumes,"h_supportPostLVShell_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	if(volumeStep.position[0] < 400){
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVShell_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]);
	}
	else if(volumeStep.position[0] < 500){
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVShell_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-7);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-7);
	}
	else if(volumeStep.position[0] < 700){
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVShell_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-16);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-16);
	}
	else if(volumeStep.position[0] < 800){
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVShell_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-28);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-28);
	}
	else{
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVShell_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-43);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-43);
	}
	//rootManager_->Get<TH3*>(whichVolumes,"h_supportPostLVBoth_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
      }
      if (volumeStep.volume == "supportPostLV"){
	rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLV") ->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVUnshifted_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]);
	if(volumeStep.position[0] < 400){
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLV_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]);
	}
	else if(volumeStep.position[0] < 500){
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLV_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-7);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-7);
	}
	else if(volumeStep.position[0] < 700){
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLV_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-16);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-16);
	}
	else if(volumeStep.position[0] < 800){
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLV_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-28);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-28);
	}
	else{
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLV_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-43);
	  rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLVBoth_2D") ->Fill(volumeStep.position[0],volumeStep.position[2]-43);
	}
	rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLV_y") ->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_supportPostLV_beam") ->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_supportPostLV_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	  //rootManager_->Get<TH3*>(whichVolumes,"h_supportPostLVBoth_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_supportPostLV_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	  //rootManager_->Get<TH3*>(whichVolumes,"h_supportPostLVBoth_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "bellowsRail"){
	rootManager_->Get<TH2*>(whichVolumes,"h_bellowsRail") ->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_bellowsRail_y") ->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_bellowsRail_beam") ->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_bellowsRail_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_bellowsRail_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "VacuumChamberCadMesh"){
	rootManager_->Get<TH2*>(whichVolumes,"h_vacuumChamberCadMesh") ->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_vacuumChamberCadMesh_y") ->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_vacuumChamberCadMesh_x") ->Fill(hitRadius,volumeStep.position[0]);
	rootManager_->Get<TH2*>(whichVolumes,"h_vacuumChamberCadMesh_beam") ->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_vacuumChamberCadMesh_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_vacuumChamberCadMesh_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "trolleyRail"){
	rootManager_->Get<TH2*>(whichVolumes,"h_trolleyRail") ->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_trolleyRail_y") ->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_trolleyRail_beam") ->Fill(hitTheta,hitBeamRadius);
	std::cout<<"volumeStep trolleyRail: "<<volumeStep.position[0]<<" "<<volumeStep.position[2]<<" "<<volumeStep.position[1]<<std::endl;
	rootManager_->Get<TH3*>(whichVolumes,"h_trolleyRail_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
      }
      if (volumeStep.volume == "xtal"){
	rootManager_->Get<TH2*>(whichVolumes,"h_xtal")->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_xtal_y")->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_xtal_beam")->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_xtal_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_xtal_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "StationNumber"){
	rootManager_->Get<TH2*>(whichVolumes,"h_StationNumber")->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_StationNumber_y")->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_StationNumber_beam")->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_StationNumber_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_StationNumber_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "PbF2Bounding"){
	rootManager_->Get<TH2*>(whichVolumes,"h_PbF2Bounding")->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_PbF2Bounding_y")->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_PbF2Bounding_beam")->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_PbF2Bounding_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_PbF2Bounding_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "Calorimeter"){
	rootManager_->Get<TH2*>(whichVolumes,"h_Calorimeter")->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_Calorimeter_y")->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_Calorimeter_beam")->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_Calorimeter_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_Calorimeter_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "insideCalo"){
	rootManager_->Get<TH2*>(whichVolumes,"h_insideCalo")->Fill(hitTheta,hitRadius); 
	rootManager_->Get<TH2*>(whichVolumes,"h_insideCalo_y")->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_insideCalo_beam")->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_insideCalo_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_insideCalo_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "photodetector"){
	rootManager_->Get<TH2*>(whichVolumes,"h_photodetector")->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_photodetector_y")->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_photodetector_beam")->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_photodetector_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_photodetector_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "outsideStorageRing"){
	std::cout<<"outsideStorageRing position: "<<volumeStep.position<<std::endl;
	rootManager_->Get<TH2*>(whichVolumes,"h_outsideStorageRing")->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_outsideStorageRing_y")->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_outsideStorageRing_beam")->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_outsideStorageRing_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_outsideStorageRing_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "frontWrapping"){
	rootManager_->Get<TH2*>(whichVolumes,"h_frontWrapping")->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_frontWrapping_y")->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_frontWrapping_beam")->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_frontWrapping_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_frontWrapping_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
      if (volumeStep.volume == "backWrapping_P"){
	rootManager_->Get<TH2*>(whichVolumes,"h_backWrapping")->Fill(hitTheta,hitRadius);
	rootManager_->Get<TH2*>(whichVolumes,"h_backWrapping_y")->Fill(hitRadius,volumeStep.position[1]);
	rootManager_->Get<TH2*>(whichVolumes,"h_backWrapping_beam")->Fill(hitTheta,hitBeamRadius);
	if(volumeStep.position[0] < -4000){
	  rootManager_->Get<TH3*>(whichVolumes,"h_backWrapping_3D") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
	else{
	  rootManager_->Get<TH3*>(whichVolumes,"h_backWrapping_3D_2") ->Fill(volumeStep.position[0],volumeStep.position[2],volumeStep.position[1]);
	}
      }
#endif
    }
  }




  // initialize histogram name
  std::string name = "";

  for (auto decayVertex : decayVertices) { 
    
    info << "decayVertex.position = " << decayVertex.position << "\n";
    info << "decayVertex.mcVertex.position = " << decayVertex.mcVertex.position << "\n";
    //std::cout<<"decayVertex.position: "<<decayVertex.position<<std::endl;
    //std::cout<<"decayVertex.momentum: "<<decayVertex.momentum<<std::endl;


    if (decayVertex.position.coordSystemName == "NULL") continue;
    if (decayVertex.position.mag() == 0) continue;
    if (decayVertex.momentum.mag() == 0) continue;
    if (keepSteps_ && decayVertex.steps.empty()) {
      continue;
    }

    int stationNum = decayVertex.island->station;
    if (stations_.find(stationNum) == stations_.end()) continue;

    // tangent point vertical and radial pos/mom
    G4ThreeVector tangentPointPos{decayVertex.position.x()/mm, decayVertex.position.y()/mm, decayVertex.position.z()/mm};
    G4ThreeVector tangentPointMom{-decayVertex.momentum.x()/MeV, -decayVertex.momentum.y()/MeV, -decayVertex.momentum.z()/MeV};
    vertexData.totalMomentum = decayVertex.momentum.mag();
    double tangentPointRadialPos   = gm2consts_->ComputeRhat( &tangentPointPos );
    vertexData.radialPosition = tangentPointRadialPos;
    double worldRadialPos = sqrt(pow(tangentPointPos[0],2)+pow(tangentPointPos[2],2));
    vertexData.worldRadialPos = worldRadialPos;
    double tangentPointVerticalPos = gm2consts_->ComputeVhat( &tangentPointPos );
    vertexData.verticalPosition = tangentPointVerticalPos;
    double tangentPointRadialMom   = gm2consts_->ComputePrhat( &tangentPointPos, &tangentPointMom );
    vertexData.radialMomentum = tangentPointRadialMom;
    double tangentPointVerticalMom = gm2consts_->ComputePvhat( &tangentPointPos, &tangentPointMom ) * decayVertex.momentum.mag();
    vertexData.verticalMomentum = tangentPointVerticalMom;
    double tangentPointTheta       = gm2consts_->ComputeTheta( &tangentPointPos ); 
    vertexData.theta = tangentPointTheta;
    //std::cout<<"Theta: "<<tangentPointTheta<<std::endl;
    double thetaDiff_0  = caloTheta_0  - tangentPointTheta;
    double thetaDiff_12 = caloTheta_12 - tangentPointTheta;
    double thetaDiff_18 = caloTheta_18 - tangentPointTheta;
    
    double thetaDiff(0.0);
    if (stationNum == 0 ) thetaDiff = thetaDiff_0;
    if (stationNum == 12) thetaDiff = thetaDiff_12;
    if (stationNum == 18) thetaDiff = thetaDiff_18;
    
    double decayArcLength = thetaDiff * (gm2consts_->R_magic() + tangentPointRadialPos) * mm2m;
    vertexData.decayArcLength = decayArcLength;

    //
    // Fill histograms
    //

    // access the time of the event
    
    //auto track  = decayVertex.track;
    //auto island = track->island;
    //auto time   = island->meanTime * 1e-3; // us
    auto time   = decayVertex.time * 1e-3; // us
    vertexData.time = time;

    if(dirName_ == "GeaneExtrapolation"){
      geaneTree_->Fill();
      //std::cout<<"GeaneExtrapolation block vector: "<<dataVector[0].size()<<std::endl;
      //std::cout<<"tangentPointRadialPos: "<<tangentPointRadialPos<<std::endl;
      
      if(tangentPointRadialPos<-45 || tangentPointRadialPos>45){
	geaneRadialPosOutliers++;
	//std::cout<<"geaneRadialPosOutliers: "<<geaneRadialPosOutliers<<std::endl;
      }
      if(worldRadialPos<7067 || worldRadialPos>7157){
	geaneRadialPosWorldOutliers++;
	//std::cout<<"geaneRadialPosWorldOutliers: "<<geaneRadialPosWorldOutliers<<std::endl;
      }
      if(tangentPointVerticalPos<-45 || tangentPointVerticalPos>45){
	geaneVerticalPosOutliers++;
	//std::cout<<"geaneVerticalPosOutliers: "<<geaneVerticalPosOutliers<<std::endl;
      }

      dataVector[0].push_back(tangentPointRadialPos);
      dataVector[1].push_back(tangentPointVerticalPos);
      dataVector[2].push_back(tangentPointRadialMom);
      dataVector[3].push_back(tangentPointVerticalMom);
      dataVector[4].push_back(tangentPointTheta);
      dataVector[5].push_back(decayArcLength);
      dataVector[6].push_back(time);
      dataVector[7].push_back(decayVertex.momentum.mag());
      dataVector[8].push_back(worldRadialPos);
      //std::cout<<"GeaneExtrapolation block vector: "<<dataVector[0].size()<<std::endl;
    }
    if(dirName_ == "Extrapolation"){
      saskiaTree_->Fill();
      //std::cout<<"Extrapolation block vector: "<<dataVector[0].size()<<std::endl;
      //std::cout<<"tangentPointRadialPos: "<<tangentPointRadialPos<<std::endl;

      if(tangentPointRadialPos<-45 || tangentPointRadialPos>45){
	radialPosOutliers++;
	//std::cout<<"radialPosOutliers: "<<radialPosOutliers<<std::endl;
      }
      if(worldRadialPos<7067 || worldRadialPos>7157){
	radialPosWorldOutliers++;
	//std::cout<<"radialPosWorldOutliers: "<<radialPosWorldOutliers<<std::endl;
      }
      if(tangentPointVerticalPos<-45 || tangentPointVerticalPos>45){
	verticalPosOutliers++;
	//std::cout<<"verticalPosOutliers: "<<verticalPosOutliers<<std::endl;
      }

      dataVector[0].push_back(tangentPointRadialPos);
      dataVector[1].push_back(tangentPointVerticalPos);
      dataVector[2].push_back(tangentPointRadialMom);
      dataVector[3].push_back(tangentPointVerticalMom);
      dataVector[4].push_back(tangentPointTheta);
      dataVector[5].push_back(decayArcLength);
      dataVector[6].push_back(time);
      dataVector[7].push_back(decayVertex.momentum.mag());
      dataVector[8].push_back(worldRadialPos);
      //std::cout<<"Extrapolation block vector: "<<dataVector[0].size()<<std::endl;
    }

    //
    // make a directory containing plots with cuts applied
    //

    bool passVolCut = false;
    bool passPvalueCut = false;

    std::vector<std::string> subDirs;
    subDirs.push_back("/vertices/allStations/allEvents");
    subDirs.push_back(Form("/vertices/station%02d/allEvents",stationNum));

    if (decayVertex.track.get()->pValue > pValueCut_) {
      passPvalueCut = true;
    }

    // pass volume cut only
    if (!decayVertex.hitVolume) {    
      passVolCut = true;
    }

    if (passPvalueCut && passVolCut) {
      subDirs.push_back(Form("/vertices/allStations/pValue>%.3f_and_noVolumesHit",pValueCut_));
      subDirs.push_back(Form("/vertices/station%02d/pValue>%.3f_and_noVolumesHit",stationNum,pValueCut_));
    }
   
    // check which volumes are hit
    if (decayVertex.hitVolume) {
      if(dirName_ == "GeaneExtrapolation"){
	for (auto hitVol : decayVertex.volumesHit) {
	  geaneVolumesHitByTracks[hitVol]++;
	  if (stationNum == 0 ) geaneVolumesHitByTracksS0 [hitVol]++;
	  if (stationNum == 12) geaneVolumesHitByTracksS12[hitVol]++;
	  if (stationNum == 18) geaneVolumesHitByTracksS18[hitVol]++;
	}
      }
      else{
	for (auto hitVol : decayVertex.volumesHit) {
	  volumesHitByTracks[hitVol]++;
	  if (stationNum == 0 ) volumesHitByTracksS0 [hitVol]++;
	  if (stationNum == 12) volumesHitByTracksS12[hitVol]++;
	  if (stationNum == 18) volumesHitByTracksS18[hitVol]++;
	}
      }
    }

    if (decayVertex.hitVolume && std::find(decayVertex.volumesHit.begin(),decayVertex.volumesHit.end(),"VacuumChamberCadMesh") == decayVertex.volumesHit.end()) {
      if(dirName_ == "GeaneExtrapolation"){
	for (auto hitVol : decayVertex.volumesHit) {
	  geaneVolumesHitByNonVacTracks[hitVol]++;
	  if (stationNum == 0 ) geaneVolumesHitByNonVacTracksS0 [hitVol]++;
	  if (stationNum == 12) geaneVolumesHitByNonVacTracksS12[hitVol]++;
	  if (stationNum == 18) geaneVolumesHitByNonVacTracksS18[hitVol]++;
	}
      }
      else{
	for (auto hitVol : decayVertex.volumesHit) {
	  volumesHitByNonVacTracks[hitVol]++;
	  if (stationNum == 0 ) volumesHitByNonVacTracksS0 [hitVol]++;
	  if (stationNum == 12) volumesHitByNonVacTracksS12[hitVol]++;
	  if (stationNum == 18) volumesHitByNonVacTracksS18[hitVol]++;
	}
      }
    }
    
    if (decayVertex.hitVolume && std::find(decayVertex.volumesHit.begin(),decayVertex.volumesHit.end(),"outsideStorageRing") != decayVertex.volumesHit.end()) {
      if(dirName_ == "GeaneExtrapolation"){
	for (auto hitVol : decayVertex.volumesHit) {
	  volumesHitByOutsideWorldTracks[hitVol]++;
	  if (stationNum == 0 ) geaneVolumesHitByOutsideWorldTracksS0 [hitVol]++;
	  if (stationNum == 12) geaneVolumesHitByOutsideWorldTracksS12[hitVol]++;
	  if (stationNum == 18) geaneVolumesHitByOutsideWorldTracksS18[hitVol]++;
	}
      }
      else{
	for (auto hitVol : decayVertex.volumesHit) {
	  volumesHitByOutsideWorldTracks[hitVol]++;
	  if (stationNum == 0 ) volumesHitByOutsideWorldTracksS0 [hitVol]++;
	  if (stationNum == 12) volumesHitByOutsideWorldTracksS12[hitVol]++;
	  if (stationNum == 18) volumesHitByOutsideWorldTracksS18[hitVol]++;
	}
      }
    }

    for (auto subDirName : subDirs) {
      subDirName = dirName_ + subDirName;
      //if(useGeane) subDirName = "GeaneExtrapolation" + subDirName;
      info << "filling histograms for dir " << subDirName << "\n";
      
      //auto subDir = rootManager_->GetDir(dir,cutDirName.c_str());
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_numSteps")->Fill(decayVertex.steps.size());
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_radialPos")  ->Fill(tangentPointRadialPos);
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_worldRadialPos")  ->Fill(worldRadialPos);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_worldRadialPos_vs_theta") ->Fill(tangentPointTheta, worldRadialPos);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_worldRadialPos_vs_theta") ->SetOption("COLZ");
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_theta_vs_p") ->Fill(decayVertex.momentum.mag(), tangentPointTheta);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_theta_vs_p") ->SetOption("COLZ");
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_verticalPos")->Fill(tangentPointVerticalPos);
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_azimuth")  ->Fill(tangentPointTheta);
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_decayArcLength")  ->Fill(decayArcLength);
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_mom")->Fill(tangentPointMom.mag());
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_radialMom")  ->Fill(tangentPointRadialMom);
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_verticalMom")->Fill(tangentPointVerticalMom);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_vertexPosSpread")->Fill(tangentPointRadialPos,tangentPointVerticalPos);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_vertexMomSpread")->Fill(tangentPointRadialMom,tangentPointVerticalMom);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_vertexPosSpread")->SetOption("COLZ");
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_vertexMomSpread")->SetOption("COLZ");
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_zPos_vs_xPos")->Fill(tangentPointPos.z(),tangentPointPos.x()); // world coords
      rootManager_->Get<TH1*>(subDirName.c_str(),"h_meanTime")->Fill(time);
	
      double g2period = (TMath::TwoPi() / gm2consts_->omegaAMagic()) * 1e-3; 
      double fracTime = time / g2period;
      int fracTimeInt = int(fracTime);
      double moduloTime = (fracTime - fracTimeInt) * g2period; // get rid of weird start bit in data
      
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_verticalPos_vs_time_fine")->Fill(time,tangentPointVerticalPos);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_radialPos_vs_time_fine")->Fill(time,tangentPointRadialPos);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_verticalPos_vs_time")->Fill(time,tangentPointVerticalPos);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_radialPos_vs_time")->Fill(time,tangentPointRadialPos);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_verticalMom_vs_time")->Fill(time,tangentPointVerticalMom);
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_radialMom_vs_time")->Fill(time,tangentPointRadialMom);
      
      // cut out first 50 us for modulo plot
      if (time > 50) {
	rootManager_->Get<TH2*>(subDirName.c_str(),"h_verticalPos_vs_time_mod")->Fill(moduloTime,tangentPointVerticalPos);
	rootManager_->Get<TH2*>(subDirName.c_str(),"h_radialPos_vs_time_mod")->Fill(moduloTime,tangentPointRadialPos);
      }
      
      rootManager_->Get<TProfile*>(subDirName.c_str(),"h_avgVerticalPos_vs_runNum")->Fill(runNumber,tangentPointVerticalPos);
      rootManager_->Get<TProfile*>(subDirName.c_str(),"h_avgRadialPos_vs_runNum")->Fill(runNumber,tangentPointRadialPos);
      
      rootManager_->Get<TH2*>(subDirName.c_str(),"h_decayArcLength_vs_mom")->Fill(tangentPointMom.mag(),decayArcLength);
	
    }
    info << "End of loop over decay vertices" << "\n";
   
  } // loop over vertices
  

  info << "End of fill histograms" << "\n\n";
  return;
}


void gm2strawtracker::GeaneExtrapolationAna::fillTruthHistograms( gm2strawtracker::DecayVertexArtRecordCollection const vertices,
								  int eventNum,
								  int runNumber, std::string dirName_)
{

  std::string subDirName = "vertices/truthComparison";
  auto truthDir = rootManager_->GetDir((TDirectory*)rootManager_->GetDir(dirName_),subDirName);
  //if(useGeane) auto truthDir = rootManager_->GetDir((TDirectory*)rootManager_->GetDir("GeaneExtrapolation"),subDirName);

  for (auto vertex : vertices) {
    
    int stationNum = vertex.island->station;
    if (stations_.find(stationNum) == stations_.end()) {
      mf::LogInfo(name_) << "continuing as station " << stationNum << " is not in stations_" << "\n";
      continue;
    }

    CoordSystem3Vector trueDecayPos = {vertex.mcVertex.position,"world"};
    CoordSystem3Vector trueDecayMom = {vertex.mcVertex.momentum,"world"};
    
    mf::LogInfo(name_) << "trueDecayPos = " << trueDecayPos << "\n";
    mf::LogInfo(name_) << "trueDecayMom = " << trueDecayMom << "\n";

    mf::LogInfo(name_) << "extrapDecayPos = " << vertex.position << "\n";
    mf::LogInfo(name_) << "extrapDecayMom = " << vertex.momentum << "\n";
    
    if (vertex.position.coordSystemName == "NULL" || vertex.momentum.coordSystemName == "NULL") {
      continue;
    }
    
    if (trueDecayPos.coordSystemName == "NULL" || trueDecayPos.mag() == 0) {
      continue;
    }
    
    if (trueDecayMom.coordSystemName == "NULL" || trueDecayMom.mag() == 0) {
      continue;
    }
    
    double track_momentum_true = trueDecayMom.mag();
    
    G4ThreeVector trueDecayPosVect = trueDecayPos.getVector();
    G4ThreeVector trueDecayMomVect = trueDecayMom.getVector();
    truthData.trueTotalMomentum = trueDecayMom.mag();
    double trueDecayPos_radial   = gm2consts_->ComputeRhat( &trueDecayPosVect );
    truthData.trueRadialPos = trueDecayPos_radial;
    double trueWorldRadialPos = sqrt(pow(trueDecayPos[0],2)+pow(trueDecayPos[2],2));
    truthData.trueWorldRadialPos = trueWorldRadialPos;
    double trueDecayPos_vertical = gm2consts_->ComputeVhat( &trueDecayPosVect );
    truthData.trueVerticalPos = trueDecayPos_vertical;
    double trueDecayPos_theta    = gm2consts_->ComputeTheta( &trueDecayPosVect );
    truthData.trueTheta = trueDecayPos_theta;
    double trueDecayMom_radial   = gm2consts_->ComputePrhat( &trueDecayPosVect, &trueDecayMomVect );
    truthData.trueRadialMom = trueDecayMom_radial;
    truthData.trueVerticalMom = gm2consts_->ComputePvhat( &trueDecayPosVect, &trueDecayMomVect) * trueDecayMomVect.mag(); 
    // tangent point variables
    G4ThreeVector recoDecayPos = vertex.position.getVector();
    G4ThreeVector recoDecayMom = vertex.momentum.getVector();
    
    double recoDecayPos_radial = gm2consts_->ComputeRhat( &recoDecayPos );
    double recoDecayPos_vertical = gm2consts_->ComputeVhat( &recoDecayPos );
    double recoDecayPos_theta    = gm2consts_->ComputeTheta( &recoDecayPos );

    // decay arc length, true
    double trueThetaDiff_0  = caloTheta_0  - trueDecayPos_theta;
    double trueThetaDiff_12 = caloTheta_12 - trueDecayPos_theta;
    double trueThetaDiff_18 = caloTheta_18 - trueDecayPos_theta;
    
    double trueThetaDiff(0.0);
    if (stationNum == 0 ) trueThetaDiff = trueThetaDiff_0;
    if (stationNum == 12) trueThetaDiff = trueThetaDiff_12;
    if (stationNum == 18) trueThetaDiff = trueThetaDiff_18;
    
    double trueDecayArcLength = trueThetaDiff * (gm2consts_->R_magic() + trueDecayPos_radial) * mm2m;
    truthData.trueArcLength = trueDecayArcLength;
    truthTree_->Fill();

    // decay arc length, reco
    double recoThetaDiff_0  = caloTheta_0  - recoDecayPos_theta;
    double recoThetaDiff_12 = caloTheta_12 - recoDecayPos_theta;
    double recoThetaDiff_18 = caloTheta_18 - recoDecayPos_theta;
    
    double recoThetaDiff(0.0);
    if (stationNum == 0 ) recoThetaDiff = recoThetaDiff_0;
    if (stationNum == 12) recoThetaDiff = recoThetaDiff_12;
    if (stationNum == 18) recoThetaDiff = recoThetaDiff_18;
    
    double recoDecayArcLength = recoThetaDiff * (gm2consts_->R_magic() + recoDecayPos_radial) * mm2m;
    
    std::vector<std::string> subDirList;
    
    subDirList.push_back("allEvents");
    if (vertex.hitVolume) subDirList.push_back("hitVolume");
    if (!vertex.hitVolume) subDirList.push_back("noVolumes");
    
    for (auto subDirName : subDirList) {
      
      TDirectory* subDir = (TDirectory*)rootManager_->GetDir(truthDir,subDirName.c_str(),true);
     
      rootManager_->Get<TH2F*>(subDir,"h_zx")->Fill(trueDecayPos.z(),trueDecayPos.x());
      rootManager_->Get<TH2F*>(subDir,"h_zx")->SetOption("COLZ");
      
      rootManager_->Get<TH2F*>(subDir,"h_truePrhat_vs_mom")->Fill(trueDecayMom_radial,track_momentum_true);
      rootManager_->Get<TProfile*>(subDir,"h_truePrhat_vs_mom_prof")->Fill(trueDecayMom_radial,track_momentum_true);

      rootManager_->Get<TH1*>(subDir,"h_verticalPos_diff")->Fill(trueDecayPos_vertical - recoDecayPos_vertical);
      rootManager_->Get<TH1*>(subDir,"h_radialPos_diff")->Fill(trueDecayPos_radial - recoDecayPos_radial);
      rootManager_->Get<TH1*>(subDir,"h_azimuth_diff")->Fill(trueDecayPos_theta - recoDecayPos_theta);
      rootManager_->Get<TH1*>(subDir,"h_decayArcLength_diff")->Fill(trueDecayArcLength - recoDecayArcLength);

      rootManager_->Get<TH1*>(subDir,"h_verticalPos")->Fill(trueDecayPos_vertical);
      rootManager_->Get<TH1*>(subDir,"h_radialPos")->Fill(trueDecayPos_radial);
      rootManager_->Get<TH2*>(subDir,"h_worldRadialPos_vs_theta")->Fill(trueDecayPos_theta, trueWorldRadialPos);
      rootManager_->Get<TH2*>(subDir,"h_theta_vs_p")->Fill(vertex.momentum.mag(), trueDecayPos_theta);
      rootManager_->Get<TH2*>(subDir,"h_theta_vs_p")->SetOption("COLZ");
      rootManager_->Get<TH1*>(subDir,"h_azimuth")->Fill(trueDecayPos_theta);
      rootManager_->Get<TH1*>(subDir,"h_decayArcLength")->Fill(trueDecayArcLength);
      
      rootManager_->Get<TH2*>(subDir,"h_true_vs_tangent_radial")->Fill(recoDecayPos_radial,trueDecayPos_radial);
      rootManager_->Get<TH2*>(subDir,"h_true_vs_tangent_vertical")->Fill(recoDecayPos_vertical,trueDecayPos_vertical);
      
      // as function of momentum
      rootManager_->Get<TH2*>(subDir,"h_verticalPos_diff_vs_mom")->Fill(vertex.momentum.mag(),trueDecayPos_vertical - recoDecayPos_vertical);
      rootManager_->Get<TH2*>(subDir,"h_radialPos_diff_vs_mom")->Fill(vertex.momentum.mag(),trueDecayPos_radial - recoDecayPos_radial);
      rootManager_->Get<TH2*>(subDir,"h_decayArcLength_diff_vs_mom")->Fill(vertex.momentum.mag(),trueDecayArcLength - recoDecayArcLength);
      
      rootManager_->Get<TProfile*>(subDir,"h_verticalPos_diff_vs_mom_prof")->Fill(vertex.momentum.mag(),trueDecayPos_vertical - recoDecayPos_vertical);
      rootManager_->Get<TProfile*>(subDir,"h_radialPos_diff_vs_mom_prof")->Fill(vertex.momentum.mag(),trueDecayPos_radial - recoDecayPos_radial);
      rootManager_->Get<TProfile*>(subDir,"h_decayArcLength_diff_vs_mom_prof")->Fill(vertex.momentum.mag(),trueDecayArcLength - recoDecayArcLength);
      
      rootManager_->Get<TH1*>(subDir,"h_mom_true")->Fill(trueDecayMom.mag());
      rootManager_->Get<TH1*>(subDir,"h_mom_diff")->Fill(trueDecayMom.mag() - vertex.momentum.mag());
      if(dirName_ == "GeaneExtrapolation"){
	rootManager_->Get<TH1*>(subDir,"h_py_diff")->Fill(trueDecayMom.y() - vertex.momentum.y());
      }
      else{
	rootManager_->Get<TH1*>(subDir,"h_py_diff")->Fill(trueDecayMom.y() + vertex.momentum.y()); // Momenta are in different directions
      }
      rootManager_->Get<TH2*>(subDir,"h_py_diff_vs_mom")->Fill(vertex.momentum.mag(),trueDecayMom.y() + vertex.momentum.y());
      rootManager_->Get<TH2*>(subDir,"h_decayArcLength_vs_mom")->Fill(trueDecayMom.mag(),trueDecayArcLength);
    }

  }
}

void gm2strawtracker::GeaneExtrapolationAna::fillCaloTruthHistograms( gm2strawtracker::DecayVertexArtRecordCollection const extrapolatedCaloHits,
								      int eventNum,
								      int runNumber )
{

  TDirectory* caloDir = rootManager_->GetDir(dirName_.c_str(),"forwards",true);
  TDirectory* caloTruthDir = rootManager_->GetDir(caloDir,"caloComparison",true);

  int itrack(0);
  for (auto extrapolatedCaloHit : extrapolatedCaloHits) {

    if (extrapolatedCaloHit.position.coordSystemName == "NULL") {
      continue;
    }

    if (extrapolatedCaloHit.mcVertex.position.mag() == 0) {
      continue;
    }

    int stationNum = extrapolatedCaloHit.track->island->station;
    std::string stationName = Form("CalorimeterNumber[%2d]",stationNum);

    gm2geom::CoordSystem3Vector forwardsExtrapPos   = extrapolatedCaloHit.position.transform(cs_,stationName);
    gm2geom::CoordSystem3Vector forwardsExtrapMom   = extrapolatedCaloHit.momentum.transform(cs_,stationName,true);

    gm2geom::CoordSystem3Vector trueCaloHitPos_world = {extrapolatedCaloHit.mcVertex.position,"world"};
    gm2geom::CoordSystem3Vector trueCaloHitMom_world = {extrapolatedCaloHit.mcVertex.momentum,"world"};

    gm2geom::CoordSystem3Vector trueCaloHitPos = trueCaloHitPos_world.transform(cs_,stationName);
    gm2geom::CoordSystem3Vector trueCaloHitMom = trueCaloHitMom_world.transform(cs_,stationName,true);
    
    mf::LogInfo(name_) << "trueCaloHitPos = " << trueCaloHitPos << "\n";
    mf::LogInfo(name_) << "forwardsExtrapPos = " << forwardsExtrapPos << "\n";

    rootManager_->Get<TH1*>(caloTruthDir,"h_trueCaloHitPos_x")->Fill(trueCaloHitPos.x());
    rootManager_->Get<TH1*>(caloTruthDir,"h_trueCaloHitPos_y")->Fill(trueCaloHitPos.y());
    rootManager_->Get<TH1*>(caloTruthDir,"h_trueCaloHitPos_z")->Fill(trueCaloHitPos.z());
    rootManager_->Get<TH1*>(caloTruthDir,"h_trueCaloHitMom_x")->Fill(trueCaloHitMom.x());
    rootManager_->Get<TH1*>(caloTruthDir,"h_trueCaloHitMom_y")->Fill(trueCaloHitMom.y());
    rootManager_->Get<TH1*>(caloTruthDir,"h_trueCaloHitMom_z")->Fill(trueCaloHitMom.z());
    rootManager_->Get<TH1*>(caloTruthDir,"h_trueCaloHitMom")  ->Fill(trueCaloHitMom.mag());
    rootManager_->Get<TH2*>(caloTruthDir,"h_TrueCaloHitMom_vs_xpos")->Fill(trueCaloHitPos.x(),trueCaloHitMom.mag());
    rootManager_->Get<TH2*>(caloTruthDir,"h_TrueCaloHitMom_vs_ypos")->Fill(trueCaloHitPos.y(),trueCaloHitMom.mag());
    rootManager_->Get<TH2*>(caloTruthDir,"h_trueCaloHitPos_xy")->Fill(trueCaloHitPos.x(),trueCaloHitPos.y());
    
    // true and reco comparison plots
    rootManager_->Get<TH1*>(caloTruthDir,"h_caloXpos_diff")->Fill(trueCaloHitPos.x() - forwardsExtrapPos.x());
    rootManager_->Get<TH1*>(caloTruthDir,"h_caloYpos_diff")->Fill(trueCaloHitPos.y() - forwardsExtrapPos.y());
    rootManager_->Get<TH1*>(caloTruthDir,"h_caloZpos_diff")->Fill(trueCaloHitPos.z() - forwardsExtrapPos.z());

    itrack++;
  }
}

void gm2strawtracker::GeaneExtrapolationAna::fillTrackHistograms( gm2strawtracker::DecayVertexArtRecordCollection const decayVertices,
								  int eventNum,
								  int runNumber )
{
  
  mf::LogInfo(name_) << "in fillTrackHistograms" << "\n";
  mf::LogInfo(name_) << "decayVertices.size() = " << decayVertices.size() << "\n";

  std::string dirName = dirName_ + "/fittedTracks";
  
  int stationNumber = 0;

  for (auto vertex : decayVertices) {

    //auto lastTrackState = vertex.track->states.back();
    //auto firsTrackState = vertex.track->states.at(0);
    auto lastTrackPosition = vertex.track->states.positionVect.back();
    auto lastTrackMomentum = vertex.track->states.momentumVect.back();
    auto firstTrackPosition = vertex.track->states.positionVect.front();
    auto firstTrackMomentum = vertex.track->states.momentumVect.front();
    
    stationNumber = vertex.island->station;
    if (stations_.find(stationNumber) == stations_.end()) {
      continue;
    }
    
    //gm2geom::CoordSystem3Vector posAtLastState = lastTrackState->position;
    //gm2geom::CoordSystem3Vector momAtLastState = lastTrackState->momentum;
    gm2geom::CoordSystem3Vector posAtLastState = lastTrackPosition;
    gm2geom::CoordSystem3Vector momAtLastState = lastTrackMomentum;
    //gm2geom::CoordSystem3Vector posAtFirsState = firsTrackState->position;
    //gm2geom::CoordSystem3Vector momAtFirsState = firsTrackState->momentum;
    gm2geom::CoordSystem3Vector posAtFirsState = firstTrackPosition;
    gm2geom::CoordSystem3Vector momAtFirsState = firstTrackMomentum;

    gm2geom::CoordSystem3Vector posAtLastState_module = posAtLastState.transform(cs_,Form("strawModuleNumber[%02d][7]",stationNumber));
    gm2geom::CoordSystem3Vector momAtLastState_module = momAtLastState.transform(cs_,Form("strawModuleNumber[%02d][7]",stationNumber),true);
    
    gm2geom::CoordSystem3Vector posAtFirsState_module = posAtFirsState.transform(cs_,Form("strawModuleNumber[%02d][0]",stationNumber));
    gm2geom::CoordSystem3Vector momAtFirsState_module = momAtFirsState.transform(cs_,Form("strawModuleNumber[%02d][0]",stationNumber),true);
    
    rootManager_->Get<TH1*>(dirName.c_str(),"h_xMom_atModule7")->Fill(momAtLastState_module.x());
    rootManager_->Get<TH1*>(dirName.c_str(),"h_xMom_atModule0")->Fill(momAtFirsState_module.x());
    
    gm2geom::CoordSystem3Vector posAtFirsState_station = posAtFirsState.transform(cs_,Form("TrackerStation[%02d]",stationNumber));
    gm2geom::CoordSystem3Vector posAtLastState_station = posAtLastState.transform(cs_,Form("TrackerStation[%02d]",stationNumber));
    
    // radial pos
    
    double radialPos_first = sqrt( pow(posAtFirsState.x(),2) + pow(posAtFirsState.z(),2) );
    double radialPos_last  = sqrt( pow(posAtLastState.x(),2) + pow(posAtLastState.z(),2) );
    
    //double xdiff = posAtFirsState_station.x() - posAtLastState_station.x();
    double xdiff = radialPos_first - radialPos_last;
    
    // << "xdiff = " << xdiff << "\n";
    
    rootManager_->Get<TH1*>(dirName.c_str(),"h_xdiff")->Fill(xdiff);
    
  }

}

void gm2strawtracker::GeaneExtrapolationAna::fillCaloHitHistograms( gm2strawtracker::DecayVertexArtRecordCollection const extrapolatedCaloHits,
								    int eventNum,
								    int runNumber )
{

  mf::LogInfo  info(name_);
  mf::LogDebug debug(name_);

  info << "extrapolatedCaloHits.size() = " << extrapolatedCaloHits.size() << "\n";

  std::string statName = "";
  
  TDirectory* dir = rootManager_->GetDir((TDirectory*)rootManager_->GetDir(dirName_), "forwards");

  TDirectory* passCut = rootManager_->GetDir(dir,"passCut",true);
  TDirectory* failCut = rootManager_->GetDir(dir,"failCut",true);

  gm2geom::CoordSystem3Vector recoCaloHitPos_world;
  gm2geom::CoordSystem3Vector recoCaloHitMom_world;

  gm2geom::CoordSystem3Vector recoCaloHitPos_calo;
  gm2geom::CoordSystem3Vector recoCaloHitMom_calo;

  for (auto recoCaloHit : extrapolatedCaloHits) {

    int stationNum = recoCaloHit.track->island->station;
    statName = Form("CalorimeterNumber[%2d]",stationNum);

    recoCaloHitPos_world = recoCaloHit.position;
    recoCaloHitMom_world = recoCaloHit.momentum;
    
    if (recoCaloHitPos_world.coordSystemName == "NULL") {
      // << "calo hit reco failed, continue" << "\n";
      continue;
    }
    
    recoCaloHitPos_calo  = recoCaloHitPos_world.transform(cs_,statName);
    recoCaloHitMom_calo  = recoCaloHitMom_world.transform(cs_,statName,true); // true = mom-like
    
    std::string dirName = "forwards";
    
    rootManager_->Get<TH1*>(dir,"h_forwardsExtrapTrackPosX")->Fill(recoCaloHitPos_calo.x());
    rootManager_->Get<TH1*>(dir,"h_forwardsExtrapTrackPosY")->Fill(recoCaloHitPos_calo.y());
    rootManager_->Get<TH1*>(dir,"h_forwardsExtrapTrackPosZ")->Fill(recoCaloHitPos_calo.z());

    rootManager_->Get<TH1*>(dir,"h_forwardsExtrapTrackMomX")->Fill(recoCaloHitMom_calo.x());
    rootManager_->Get<TH1*>(dir,"h_forwardsExtrapTrackMomY")->Fill(recoCaloHitMom_calo.y());
    rootManager_->Get<TH1*>(dir,"h_forwardsExtrapTrackMomZ")->Fill(recoCaloHitMom_calo.z());

    rootManager_->Get<TH1*>(dir,"h_forwardsExtrapTrackMom" )->Fill(recoCaloHitMom_calo.mag());

    G4ThreeVector mom = recoCaloHitMom_world.getVector();
    G4ThreeVector pos = recoCaloHitPos_world.getVector();

    double radialMom = gm2consts_->ComputePrhat(&pos,&mom) * recoCaloHitMom_world.mag();
    rootManager_->Get<TH1*>(dir,"h_forwardsExtrapTrackMom_radial")->Fill(radialMom);

    // angle of incidence to calo face
    double incidentAngle = atan2(-recoCaloHitMom_calo.x(),recoCaloHitMom_calo.mag()) * (180 / M_PI);
    double incidentVerticalAngle = atan2(-recoCaloHitMom_calo.y(),recoCaloHitMom_calo.mag()) * (180 / M_PI);
 
    rootManager_->Get<TH1*>(dir,"h_incidentAngle")->Fill(incidentAngle);
    rootManager_->Get<TH2*>(dir,"h_incidentAngle_vs_mom")->Fill(recoCaloHitMom_calo.mag(),incidentAngle);

    rootManager_->Get<TH1*>(dir,"h_incidentVerticalAngle")->Fill(incidentVerticalAngle);
    rootManager_->Get<TH2*>(dir,"h_incidentVerticalAngle_vs_mom")->Fill(recoCaloHitMom_calo.mag(),incidentVerticalAngle);

    if (recoCaloHitMom_calo.mag() < 2500) {
      rootManager_->Get<TH1*>(passCut,"h_forwardsExtrapTrackMom_radial")->Fill(radialMom * recoCaloHitMom_world.mag()); // +ve = radially out of ring
      rootManager_->Get<TH1*>(passCut,"h_forwardsExtrapTrackMomX")->Fill(recoCaloHitMom_calo.x()); // +ve = radially out of ring
      rootManager_->Get<TH1*>(passCut,"h_forwardsExtrapTrackMomY")->Fill(recoCaloHitMom_calo.y());
      rootManager_->Get<TH1*>(passCut,"h_forwardsExtrapTrackMomZ")->Fill(recoCaloHitMom_calo.z());
      rootManager_->Get<TH1*>(passCut,"h_incidentAngle")->Fill(incidentAngle);
      rootManager_->Get<TH2*>(passCut,"h_incidentAngle_vs_mom")->Fill(recoCaloHitMom_calo.mag(),incidentAngle);
      rootManager_->Get<TH1*>(passCut,"h_incidentVerticalAngle")->Fill(incidentVerticalAngle);
      rootManager_->Get<TH2*>(passCut,"h_incidentVerticalAngle_vs_mom")->Fill(recoCaloHitMom_calo.mag(),incidentVerticalAngle);
    }

    else {
      rootManager_->Get<TH1*>(failCut,"h_forwardsExtrapTrackMom_radial")->Fill(radialMom * recoCaloHitMom_world.mag()); // +ve = radially out of ring
      rootManager_->Get<TH1*>(failCut,"h_forwardsExtrapTrackMomX")->Fill(recoCaloHitMom_calo.x()); // +ve = radially out of ring
      rootManager_->Get<TH1*>(failCut,"h_forwardsExtrapTrackMomY")->Fill(recoCaloHitMom_calo.y());
      rootManager_->Get<TH1*>(failCut,"h_forwardsExtrapTrackMomZ")->Fill(recoCaloHitMom_calo.z());

      rootManager_->Get<TH1*>(failCut,"h_incidentAngle")->Fill(incidentAngle);
      rootManager_->Get<TH2*>(failCut,"h_incidentAngle_vs_mom")->Fill(recoCaloHitMom_calo.mag(),incidentAngle);

      rootManager_->Get<TH1*>(failCut,"h_incidentVerticalAngle")->Fill(incidentVerticalAngle);
      rootManager_->Get<TH2*>(failCut,"h_incidentVerticalAngle_vs_mom")->Fill(recoCaloHitMom_calo.mag(),incidentVerticalAngle);
    }

    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackPos_xy")->Fill(recoCaloHitPos_calo.x(),recoCaloHitPos_calo.y());
    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackPos_xy")->SetOption("COLZ");
    
    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackMom_xy")->Fill(recoCaloHitMom_calo.x(),recoCaloHitMom_calo.y());
    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackMom_xy")->SetOption("COLZ");
    
    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackMom_zx")->Fill(recoCaloHitMom_calo.z(),recoCaloHitMom_calo.x());
    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackMom_zx")->SetOption("COLZ");
    
    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackMom_zy")->Fill(recoCaloHitMom_calo.z(),recoCaloHitMom_calo.y());
    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackMom_zy")->SetOption("COLZ");

    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackMom_vs_xpos")->Fill(recoCaloHitPos_calo.x(),recoCaloHitMom_calo.mag());
    rootManager_->Get<TH2*>(dir,"h_forwardsExtrapTrackMom_vs_ypos")->Fill(recoCaloHitPos_calo.y(),recoCaloHitMom_calo.mag());
    
    double meanTime = recoCaloHit.track->island->meanTime;
    
    rootManager_->Get<TH1*>(dir,"h_meanIslandTime")->Fill(meanTime*1e-3);

    //
    // only plot events that will hit calo
    //
  
    double xmin_station = -120.0;
    double xmax_station = 120.0;
    double ymin_station = -80.0;
    double ymax_station = 80.0;

    if ( recoCaloHitPos_calo.x() >= xmin_station && recoCaloHitPos_calo.x() <= xmax_station ) {
      if (recoCaloHitPos_calo.y() >= ymin_station && recoCaloHitPos_calo.y() <= ymax_station) {

	dirName = dirName_ + "/forwards/inCaloRegion";
	
	rootManager_->Get<TH1*>(dirName.c_str(),"h_forwardsExtrapTrackPosX")->Fill(recoCaloHitPos_calo.x());
	rootManager_->Get<TH1*>(dirName.c_str(),"h_forwardsExtrapTrackPosY")->Fill(recoCaloHitPos_calo.y());
	rootManager_->Get<TH1*>(dirName.c_str(),"h_forwardsExtrapTrackPosZ")->Fill(recoCaloHitPos_calo.z());
	rootManager_->Get<TH1*>(dirName.c_str(),"h_forwardsExtrapTrackMomX")->Fill(recoCaloHitMom_calo.x());
	rootManager_->Get<TH1*>(dirName.c_str(),"h_forwardsExtrapTrackMomY")->Fill(recoCaloHitMom_calo.y());
	rootManager_->Get<TH1*>(dirName.c_str(),"h_forwardsExtrapTrackMomZ")->Fill(recoCaloHitMom_calo.z());
	rootManager_->Get<TH1*>(dirName.c_str(),"h_forwardsExtrapTrackMom" )->Fill(recoCaloHitMom_calo.mag());

	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackPos_xy")->Fill(recoCaloHitPos_calo.x(),recoCaloHitPos_calo.y());
	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackPos_xy")->SetOption("COLZ");
	
	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackMom_xy")->Fill(recoCaloHitMom_calo.x(),recoCaloHitMom_calo.y());
	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackMom_xy")->SetOption("COLZ");
	
	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackMom_zx")->Fill(recoCaloHitMom_calo.z(),recoCaloHitMom_calo.x());
	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackMom_zx")->SetOption("COLZ");
	
	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackMom_zy")->Fill(recoCaloHitMom_calo.z(),recoCaloHitMom_calo.y());
	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackMom_zy")->SetOption("COLZ");
	
	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackMom_vs_xpos")->Fill(recoCaloHitPos_calo.x(),recoCaloHitMom_calo.mag());
	rootManager_->Get<TH2*>(dirName.c_str(),"h_forwardsExtrapTrackMom_vs_ypos")->Fill(recoCaloHitPos_calo.y(),recoCaloHitMom_calo.mag());

	rootManager_->Get<TH1*>(dirName.c_str(),"h_incidentAngle")->Fill(incidentAngle);
	rootManager_->Get<TH2*>(dirName.c_str(),"h_incidentAngle_vs_mom")->Fill(recoCaloHitMom_calo.mag(),incidentAngle);

	rootManager_->Get<TH1*>(dirName.c_str(),"h_meanIslandTime")->Fill(meanTime*1e-3);
	
      }
    }
  }
}

void gm2strawtracker::GeaneExtrapolationAna::beginJob()
{
  mf::LogInfo info(name_);
  info << "Enter GeaneExtrapolationAna::beginJob\n";

  // create a root file and manager
  art::ServiceHandle<art::TFileService> tfs;
  auto& outputRootFile_ = tfs->file();
  rootManager_.reset( new RootManager("GeaneExtrapolationAnaPlots",&outputRootFile_) );  

  
  //tree
  saskiaTree_ = tfs->make<TTree>("saskiaData","Data from extrapolation");
  geaneTree_ = tfs->make<TTree>("geaneData","Data from extrapolation");
  truthTree_ = tfs->make<TTree>("truthData","Data from extrapolation");
  diffTree_ = tfs->make<TTree>("diffData", "Data from extrapolation");
  saskiaTree_->Branch("saskiaData", &vertexData, "radialPosition/D:verticalPosition:radialMomentum:verticalMomentum:theta:decayArcLength:time:totalMomentum:worldRadialPos"); 
  geaneTree_->Branch("geaneData", &vertexData, "radialPosition/D:verticalPosition:radialMomentum:verticalMomentum:theta:decayArcLength:time:totalMomentum:worldRadialPos");
  truthTree_->Branch("truthData", &truthData, "trueRadialPos/D:trueVerticalPos:trueRadialMom:trueVerticalMom:trueTheta:trueArcLength:trueTotalMomentum:trueWorldRadialPos");
  diffTree_->Branch("diffData", &diffData, "diffRadialPosition/D:diffVerticalPosition:diffRadialMomentum:diffVerticalMomentum:diffTheta:diffDecayArcLength:diffTime:diffTotalMomentum:diffWorldRadialPos");


  // create directories
  std::vector<std::string> topDirList = {"Extrapolation", "GeaneExtrapolation"};
  for (auto dirName_ : topDirList){
    auto topDir    = rootManager_->GetDir(dirName_,true); // true -> create if doesn't exist
    auto vertexDir = rootManager_->GetDir(topDir,"vertices",true);
    auto failureModeDir = rootManager_->GetDir(topDir,"failureModes",true);
    auto whichVolumes = rootManager_->GetDir(topDir,"whichVolumes",true);
    auto caloDir   = rootManager_->GetDir(topDir,"forwards",true);
    auto trackDir  = rootManager_->GetDir(topDir,"fittedTracks",true);
    auto truthDir  = rootManager_->GetDir(vertexDir,"truthComparison",true);
    auto caloTruthDir = rootManager_->GetDir(caloDir,"caloComparison",true);
    auto geomDir   = rootManager_->GetDir(topDir,"strawGeometry",true);

    bookDecayVertexPlots(vertexDir);
    bookFailedEventPlots(failureModeDir);
    bookWhichVolumePlots(whichVolumes);
    bookRecoCaloHitPlots(caloDir);
    bookTrackPlots(trackDir);
    bookTruthPlots(truthDir);
    bookCaloTruthPlots(caloTruthDir);
    bookGeomGraphs(geomDir);

  }

  info << "Exit GeaneExtrapolationAna::beginJob\n";
 
}

void gm2strawtracker::GeaneExtrapolationAna::bookTruthPlots( TDirectory* dir )
{

  std::vector<std::string> subDirList = {"allEvents","hitVolume","noVolumes"};

  for (auto subDirName : subDirList) {
   
    TDirectory* subDir = rootManager_->GetDir(dir,subDirName.c_str(),true);
    
    // book decay vertex and truth comparison plots
    rootManager_->Add(subDir, new TH1F("h_verticalPos_diff","True - tangent point;#Deltay[mm];",200,-100,100));
    rootManager_->Add(subDir, new TH1F("h_radialPos_diff" ,"True - tangent point;#Deltax[mm];",200,-100,100));
    rootManager_->Add(subDir, new TH1F("h_azimuth_diff" ,"True - tangent point;#Delta#theta[rad];",200,-1,1));
    rootManager_->Add(subDir, new TH1F("h_decayArcLength_diff" ,"True - tangent point;True - tangent decay arc length [m];",100,-2,2));

    rootManager_->Add(subDir, new TH1F("h_verticalPos","True vertical pos;y [mm];",100,-50,50));
    rootManager_->Add(subDir, new TH1F("h_radialPos" ,"True radial pos;x [mm];",100,-50,50));
    rootManager_->Add(subDir, new TH2F("h_worldRadialPos_vs_theta", "Truth;worldRadialPos [mm]; Theta",30,3,6,400,7000,7200));
    rootManager_->Add(subDir, new TH2F("h_theta_vs_p","Truth;|p|[MeV/c];Theta",1000,0,3000,12,0,6));
    rootManager_->Add(subDir, new TH1F("h_azimuth" ,"True theta;#theta [rad];",100*M_PI,0,2*M_PI));
    rootManager_->Add(subDir, new TH1F("h_decayArcLength" ,"True azimuth;True decay arc length [m];",250,0,10));
    
    rootManager_->Add(subDir, new TH2F("h_verticalPos_diff_vs_mom","True - tangent point;Momentum [MeV];#Deltay[mm]",35,0,3500,50,-50,50));
    rootManager_->Add(subDir, new TH2F("h_radialPos_diff_vs_mom"  ,"True - tangent point;Momentum [MeV];#Deltax[mm]",35,0,3500,50,-50,50));
    rootManager_->Add(subDir, new TH2F("h_azimuth_diff_vs_mom"  ,"True - tangent point;Momentum [MeV];#Delta#theta[rad]",35,0,3500,200,-10,10));
    rootManager_->Add(subDir, new TH2F("h_decayArcLength_diff_vs_mom"  ,"True - tangent point;Momentum [MeV];True - tangent decay arc length [m]",35,0,3500,100,-2,2));
    
    rootManager_->Add(subDir, new TProfile("h_verticalPos_diff_vs_mom_prof","True - tangent point;Momentum [MeV];#Deltay[mm]",35,0,3500));
    rootManager_->Add(subDir, new TProfile("h_radialPos_diff_vs_mom_prof"  ,"True - tangent point;Momentum [MeV];#Deltax[mm]",35,0,3500));
    rootManager_->Add(subDir, new TProfile("h_decayArcLength_diff_vs_mom_prof"  ,"True - tangent point;Momentum [MeV];True - tangent decay arc length [m]",35,0,3500));
    
    rootManager_->Add(subDir, new TH2F("h_zx","Truth;z [mm];x[mm]",750,-7500,7500,750,-7500,7500));
    rootManager_->Add(subDir, new TH2F("h_truePrhat_vs_mom","Truth;Momentum [MeV];prhat",35,0,3500,1000,-0.1,0.1));
    rootManager_->Add(subDir, new TProfile("h_truePrhat_vs_mom_prof","Truth;Momentum [MeV];Average prhat",35,0,3500,-0.1,0.1));
    
    // for unfolding
    rootManager_->Add(subDir, new TH2F("h_true_vs_tangent_radial"  ,";Tangency radial pos [mm];true radial pos [mm]",200,-100,100,200,-100,100));
    rootManager_->Add(subDir, new TH2F("h_true_vs_tangent_vertical",";Tangency vertical pos [mm];true vertical pos [mm]",200,-100,100,200,-100,100));

    // momentum
    rootManager_->Add(subDir, new TH1F("h_mom_diff","True - tangent point;#Delta |p| [MeV];",1000,-500,500));
    rootManager_->Add(subDir, new TH1F("h_mom_true","True momentum;|p| [MeV];",350,0,3500));
    
    // vertical momentum
    rootManager_->Add(subDir, new TH1F("h_py_diff","True - tangent point;#Delta py [MeV];",1000,-500,500));
    rootManager_->Add(subDir, new TH2F("h_py_diff_vs_mom","True - tangent point;Momentum [MeV];#Delta py [MeV];",350,0,3500,1000,-500,500));

    // arc length vs momentum
    rootManager_->Add(subDir, new TH2F("h_decayArcLength_vs_mom",";Momentum [MeV];Decay arc length [m]",128,0,3200,250,0,10));
  }  
}   
  
void gm2strawtracker::GeaneExtrapolationAna::bookCaloTruthPlots( TDirectory* dir )
{

  rootManager_->Add(dir, new TH1F("h_trueCaloHitPos_x",";x pos [mm];",120,-600,600));
  rootManager_->Add(dir, new TH1F("h_trueCaloHitPos_y",";y pos [mm];",120,-600,600));
  rootManager_->Add(dir, new TH1F("h_trueCaloHitPos_z",";z pos [mm];",100,-75,-65));
  rootManager_->Add(dir, new TH1F("h_trueCaloHitMom_x",";x mom [MeV];",200,-1000,1000));
  rootManager_->Add(dir, new TH1F("h_trueCaloHitMom_y",";y mom [MeV];",200,-1000,1000));
  rootManager_->Add(dir, new TH1F("h_trueCaloHitMom_z",";z mom [MeV];",350,0,3500));
  rootManager_->Add(dir, new TH1F("h_trueCaloHitMom",";|p| [MeV];",350,0,3500));

  rootManager_->Add(dir, new TH2F("h_TrueCaloHitMom_vs_xpos",";x pos [mm];Momentum [MeV]",60,-300,300,40,0,4000));
  rootManager_->Add(dir, new TH2F("h_TrueCaloHitMom_vs_ypos",";y pos [mm];Momentum [MeV]",40,-200,200,40,0,4000));

  rootManager_->Add(dir, new TH2F("h_trueCaloHitPos_xy",";x pos [mm];y pos [mm]",60,-200,400,30,-150,150));

  // true and reco comparison plots
  rootManager_->Add(dir, new TH1F("h_caloXpos_diff","True - extrapolated calo pos;#Deltax[mm];",2000,-5,5));
  rootManager_->Add(dir, new TH1F("h_caloYpos_diff","True - extrapolated calo pos;#Deltay[mm];",2000,-5,5));
  rootManager_->Add(dir, new TH1F("h_caloZpos_diff","True - extrapolated calo pos;#Deltaz[mm];",2000,-5,5));

}   

void gm2strawtracker::GeaneExtrapolationAna::bookTrackPlots( TDirectory* dir )
{
  rootManager_->Add(dir, new TH1F("h_xMom_atModule7","Module 7;p_{x} [MeV];N tracks",1000,-500,500));
  rootManager_->Add(dir, new TH1F("h_xMom_atModule0","Module 0;p_{x} [MeV];N tracks",1000,-500,500));
  rootManager_->Add(dir, new TH1F("h_xdiff",";Module 0 x pos - Module 7 x pos [mm];",100,-1000,1000));
}

void gm2strawtracker::GeaneExtrapolationAna::bookRecoCaloHitPlots(TDirectory* dir)
{
  
  TDirectory* passCut = rootManager_->GetDir(dir,"passCut",true);
  TDirectory* failCut = rootManager_->GetDir(dir,"failCut",true);

  rootManager_->Add(dir, new TH1F("h_forwardsExtrapTrackPosX","Forwards extrapolation to calo;extrapolated track pos_{x} [mm];",120,-600,600));
  rootManager_->Add(dir, new TH1F("h_forwardsExtrapTrackPosY","Forwards extrapolation to calo;extrapolated track pos_{y} [mm];",120,-600,600));
  rootManager_->Add(dir, new TH1F("h_forwardsExtrapTrackPosZ","Forwards extrapolation to calo;extrapolated track pos_{z} [mm];",10,-70,-71));//1100,1200));
  rootManager_->Add(dir, new TH1F("h_forwardsExtrapTrackMomX","Forwards extrapolation to calo;extrapolated track p_{x} [MeV];" ,200,-1000,1000));
  rootManager_->Add(dir, new TH1F("h_forwardsExtrapTrackMomY","Forwards extrapolation to calo;extrapolated track p_{y} [MeV];" ,200,-1000,1000));
  rootManager_->Add(dir, new TH1F("h_forwardsExtrapTrackMomZ","Forwards extrapolation to calo;extrapolated track p_{z} [MeV];" ,35,0,3500));
  rootManager_->Add(dir, new TH1F("h_forwardsExtrapTrackMom" ,"Forwards extrapolation to calo;extrapolated track |p| [MeV];"   ,35,0,3500));

  rootManager_->Add(dir, new TH1F("h_forwardsExtrapTrackMom_radial",";extrapolated track radial momentum [MeV];" ,400,-500,3500));

  rootManager_->Add(failCut, new TH1F("h_forwardsExtrapTrackMom_radial",";extrapolated track radial momentum [MeV];" ,400,-500,3500));
  rootManager_->Add(failCut, new TH1F("h_forwardsExtrapTrackMomX","|p| > 2500 MeV;extrapolated track p_{x} [MeV];" ,200,-1000,1000));
  rootManager_->Add(failCut, new TH1F("h_forwardsExtrapTrackMomY","|p| > 2500 MeV;extrapolated track p_{y} [MeV];" ,200,-1000,1000));
  rootManager_->Add(failCut, new TH1F("h_forwardsExtrapTrackMomZ","|p| > 2500 MeV;extrapolated track p_{z} [MeV];" ,350,0,3500));

  rootManager_->Add(passCut, new TH1F("h_forwardsExtrapTrackMom_radial",";extrapolated track radial momentum [MeV];" ,400,-500,3500));
  rootManager_->Add(passCut, new TH1F("h_forwardsExtrapTrackMomX","|p| < 2500 MeV;extrapolated track p_{x} [MeV];" ,200,-1000,1000));
  rootManager_->Add(passCut, new TH1F("h_forwardsExtrapTrackMomY","|p| < 2500 MeV;extrapolated track p_{y} [MeV];" ,200,-1000,1000));
  rootManager_->Add(passCut, new TH1F("h_forwardsExtrapTrackMomZ","|p| < 2500 MeV;extrapolated track p_{z} [MeV];" ,350,0,3500));

  rootManager_->Add(dir, new TH2F("h_forwardsExtrapTrackPos_xy","Forwards extrapolation to calo;pos_{x} [mm];pos_{y} [mm]",60,-600,600,30,-150,150));
  rootManager_->Add(dir, new TH2F("h_forwardsExtrapTrackMom_xy","Forwards extrapolation to calo;p_{x} [MeV];p_{y} [MeV]",60,-400,800,30,-150,150));
  rootManager_->Add(dir, new TH2F("h_forwardsExtrapTrackMom_zx","Forwards extrapolation to calo;p_{z} [MeV];p_{x} [MeV]",35,0,3500,60,-400,800));
  rootManager_->Add(dir, new TH2F("h_forwardsExtrapTrackMom_zy","Forwards extrapolation to calo;p_{z} [MeV];p_{y} [MeV]",35,0,3500,30,-150,150));

  rootManager_->Add(dir, new TH2F("h_forwardsExtrapTrackMom_vs_xpos","Forwards extrapolation to calo;Extrapolated x pos [mm];Track momentum [MeV]",60,-300,300,40,0,4000));
  rootManager_->Add(dir, new TH2F("h_forwardsExtrapTrackMom_vs_ypos","Forwards extrapolation to calo;Extrapolated y pos [mm];Track momentum [MeV]",40,-200,200,40,0,4000));

  rootManager_->Add(dir, new TH1F("h_meanIslandTime","Forwards extraplation to calo;Mean island time [#mus];N tracks",100,0,1000));

  rootManager_->Add(dir, new TH1F("h_incidentAngle","Angle of incidence;Incident angle [degrees];",180,-45,45));
  rootManager_->Add(dir, new TH2F("h_incidentAngle_vs_mom","Angle of incidence vs momentum;Track momentum [MeV];Incident angle [degrees]",400,0,4000,50,-10,40));

  rootManager_->Add(passCut, new TH1F("h_incidentAngle","Angle of incidence;Incident angle [degrees];",180,-45,45));
  rootManager_->Add(passCut, new TH2F("h_incidentAngle_vs_mom","Angle of incidence vs momentum;Track momentum [MeV];Incident angle [degrees]",300,0,3000,90,-45,45));

  rootManager_->Add(failCut, new TH1F("h_incidentAngle","Angle of incidence;Incident angle [degrees];",180,-45,45));
  rootManager_->Add(failCut, new TH2F("h_incidentAngle_vs_mom","Angle of incidence vs momentum;Track momentum [MeV];Incident angle [degrees]",200,2000,4000,50,-10,40));

  rootManager_->Add(dir, new TH1F("h_incidentVerticalAngle","Angle of incidence (vertical);Incident angle [degrees];",180,-45,45));
  rootManager_->Add(dir, new TH2F("h_incidentVerticalAngle_vs_mom","Angle of incidence vs momentum;Track momentum [MeV];Incident angle [degrees]",400,0,4000,900,-45,45));

  rootManager_->Add(passCut, new TH1F("h_incidentVerticalAngle","Angle of incidence (vertical);Incident angle [degrees];",180,-45,45));
  rootManager_->Add(passCut, new TH2F("h_incidentVerticalAngle_vs_mom","Angle of incidence vs momentum;Track momentum [MeV];Incident angle [degrees]",300,0,3000,900,-45,45));

  rootManager_->Add(failCut, new TH1F("h_incidentVerticalAngle","Angle of incidence (vertical);Incident angle [degrees];",180,-45,45));
  rootManager_->Add(failCut, new TH2F("h_incidentVerticalAngle_vs_mom","Angle of incidence vs momentum;Track momentum [MeV];Incident angle [degrees]",200,2000,4000,900,-45,45));
  
  TDirectory* cut = rootManager_->GetDir(dir,"inCaloRegion",true);
  rootManager_->Add(cut, new TH1F("h_forwardsExtrapTrackPosX","Forwards extrapolation to calo;extrapolated track pos_{x} [mm];",120,-600,600));
  rootManager_->Add(cut, new TH1F("h_forwardsExtrapTrackPosY","Forwards extrapolation to calo;extrapolated track pos_{y} [mm];",120,-600,600));
  rootManager_->Add(cut, new TH1F("h_forwardsExtrapTrackPosZ","Forwards extrapolation to calo;extrapolated track pos_{z} [mm];",100,-75,-65));//1100,1200));

  rootManager_->Add(cut, new TH1F("h_forwardsExtrapTrackMomX","Forwards extrapolation to calo;extrapolated track p_{x} [MeV];" ,200,-1000,1000));
  rootManager_->Add(cut, new TH1F("h_forwardsExtrapTrackMomY","Forwards extrapolation to calo;extrapolated track p_{y} [MeV];" ,200,-1000,1000));
  rootManager_->Add(cut, new TH1F("h_forwardsExtrapTrackMomZ","Forwards extrapolation to calo;extrapolated track p_{z} [MeV];" ,350,0,3500));
  rootManager_->Add(cut, new TH1F("h_forwardsExtrapTrackMom" ,"Forwards extrapolation to calo;extrapolated track |p| [MeV];"   ,350,0,3500));

  rootManager_->Add(cut, new TH1F("h_forwardsExtrapTrackMom_radial",";extrapolated track radial momentum [MeV];",400,-1000,1000));

  rootManager_->Add(cut, new TH2F("h_forwardsExtrapTrackPos_xy","Forwards extrapolation to calo;pos_{x} [mm];pos_{y} [mm]",120,-200,400,60,-150,150));
  rootManager_->Add(cut, new TH2F("h_forwardsExtrapTrackMom_xy","Forwards extrapolation to calo;p_{x} [MeV];p_{y} [MeV]",60,-800,400,30,-150,150));
  rootManager_->Add(cut, new TH2F("h_forwardsExtrapTrackMom_zx","Forwards extrapolation to calo;p_{z} [MeV];p_{x} [MeV]",35,0,3500,60,-400,800));
  rootManager_->Add(cut, new TH2F("h_forwardsExtrapTrackMom_zy","Forwards extrapolation to calo;p_{z} [MeV];p_{y} [MeV]",35,0,3500,30,-150,150));

  rootManager_->Add(cut, new TH2F("h_forwardsExtrapTrackMom_vs_xpos","Forwards extrapolation to calo;Extrapolated x pos [mm];Track momentum [MeV]",60,-300,300,40,0,4000));
  rootManager_->Add(cut, new TH2F("h_forwardsExtrapTrackMom_vs_ypos","Forwards extrapolation to calo;Extrapolated y pos [mm];Track momentum [MeV]",40,-200,200,40,0,4000));

  rootManager_->Add(cut, new TH1F("h_incidentAngle","Angle of incidence;Incident angle [degrees];",180,-45,45));
  rootManager_->Add(cut, new TH2F("h_incidentAngle_vs_mom","Angle of incidence vs momentum;Track momentum [MeV];Incident angle [degrees]",400,0,4000,90,-45,45));

  rootManager_->Add(cut, new TH1F("h_incidentVerticalAngle","Angle of incidence (vertical);Incident angle [degrees];",180,-45,45));
  rootManager_->Add(cut, new TH2F("h_incidentVerticalAngle_vs_mom","Angle of incidence vs momentum;Track momentum [MeV];Incident angle [degrees]",400,0,4000,900,-45,45));

  rootManager_->Add(cut, new TH1F("h_meanIslandTime","Forwards extraplation to calo;Mean island time [#mus];N tracks",100,0,1000));

}

void gm2strawtracker::GeaneExtrapolationAna::bookFailedEventPlots(TDirectory* dir)
{
  rootManager_->Add(dir, new TH1F("h_failureModesBackwards","Failure modes for backwards extrap",15,0,15));
  rootManager_->Add(dir, new TH1F("h_failureModesForwards","Failure modes for forwards extrap",15,0,15));
}

void gm2strawtracker::GeaneExtrapolationAna::bookWhichVolumePlots(TDirectory* dir)
{
  rootManager_->Add(dir, new TH1F("h_volumesHit","Volumes hit by all events",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHit_noVacuumChamber","Volumes hit by all events that do not hit vac chamber",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHit_outsideWorld","Volumes hit by all events that go on to leave the SR",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHitS0","Volumes hit by all events, station 0",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHit_noVacuumChamberS0","Volumes hit by all events that do not hit vac chamber, station 0",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHit_outsideWorldS0","Volumes hit by all events that go on to leave the SR, station 0",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHitS12","Volumes hit by all events, station 12",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHit_noVacuumChamberS12","Volumes hit by all events that do not hit vac chamber, station 12",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHit_outsideWorldS12","Volumes hit by all events that go on to leave the SR, station 12",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHitS18","Volumes hit by all events, station 18",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHit_noVacuumChamberS18","Volumes hit by all events that do not hit vac chamber, station 18",30,0,30));
  rootManager_->Add(dir, new TH1F("h_volumesHit_outsideWorldS18","Volumes hit by all events that go on to leave the SR, station 18",30,0,30));

  rootManager_->Add(dir, new TH3F("h_world_3D","Position of steps that any volume",1600,-7400,7400,1600,-7400,7400,50,-100,100));

  rootManager_->Add(dir, new TH2F("h_impactVectors_2D","Impact vectors for tracks near supportPost",100,-10,10,100,-10,10));

  rootManager_->Add(dir, new TH2F("h_supportPostLVShell_2D","Global z vs x position for supportPostLVShell hits",1000,320,890,1000,-7080,-6960));
  rootManager_->Add(dir, new TH2F("h_supportPostLVBoth_2D","Global z vs x position for supportPostLV and shell hits",1000,320,890,1000,-7080,-6960));
  rootManager_->Add(dir, new TH3F("h_supportPostLVShell_3D","Position of steps that hit supportPostLVShell",500,200,1000,500,-7400,-6800,200,-100,100));
  //rootManager_->Add(dir, new TH3F("h_supportPostLVBoth_3D","Position of steps that hit supportPostLV and shell",500,200,1000,500,-7400,-6800,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_supportPostLV","Global radius vs theta for steps that hit supportPostLV",28,-7,7,400,7000,7200));
  rootManager_->Add(dir, new TH2F("h_supportPostLV_2D","Global z vs x position for supportPostLV hits",1000,320,890,1000,-7080,-6960));
  rootManager_->Add(dir, new TH2F("h_supportPostLVUnshifted_2D","Global z vs x position for supportPostLV hits",1000,300,900,1000,-7080,-6960));
  rootManager_->Add(dir, new TH2F("h_supportPostLV_y","Global y vs radius for steps that hit supportPostLV",400,7000,7200,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_supportPostLV_beam","Beam radius vs theta for steps that hit supportPostLV",28,-7,7,400,-100,100));
  rootManager_->Add(dir, new TH3F("h_supportPostLV_3D","Position of steps that hit supportPostLV",500,-7400,-6900,500,-1400,0,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_supportPostLV_3D_2","Position of steps that hit supportPostLV",500,200,1000,500,-7400,-6800,200,-100,100));
  

  rootManager_->Add(dir, new TH2F("h_bellowsRail","Global radius vs theta for steps that hit bellowsRail",28,-7,7,400,7000,7200));
  rootManager_->Add(dir, new TH2F("h_bellowsRail_y","Global y vs radius for steps that hit bellowsRail",400,7000,7200,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_bellowsRail_beam","Beam radius vs theta for steps that hit bellowsRail",28,-7,7,400,-100,100));
  rootManager_->Add(dir, new TH3F("h_bellowsRail_3D","Position of steps that hit bellowsRail",500,-7400,-6800,500,-400,500,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_bellowsRail_3D_2","Position of steps that hit bellowsRail",500,-500,500,500,-7400,-6800,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_vacuumChamberCadMesh","Global radius vs theta for steps that hit vacuumChamberCadMesh",28,-7,7,740,6800,7170));
  rootManager_->Add(dir, new TH2F("h_vacuumChamberCadMesh_y","Global y vs radius for steps that hit vacuumChamberCadMesh",800,6800,7200,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_vacuumChamberCadMesh_x","Global x vs radius for steps that hit vacuumChamberCadMesh",800,6800,7200,800,-7200,7200));
  rootManager_->Add(dir, new TH2F("h_vacuumChamberCadMesh_beam","Beam radius vs theta for steps that hit vacuumChamberCadMesh",28,-7,7,300,-200,-50));
  rootManager_->Add(dir, new TH3F("h_vacuumChamberCadMesh_3D","Position of steps that hit vacuumChamberCadMesh",500,-7400,-6300,500,-1000,1000,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_vacuumChamberCadMesh_3D_2","Position of steps that hit vacuumChamberCadMesh",500,-3000,1000,500,-7400,-6400,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_trolleyRail","Global radius vs theta for steps that hit trolleyRail",28,-7,7,400,7000,7200));
  rootManager_->Add(dir, new TH2F("h_trolleyRail_y","Global y vs radius for steps that hit trolleyRail",400,7000,7200,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_trolleyRail_beam","Beam radius vs theta for steps that hit trolleyRail",28,-7,7,400,-100,100));
  rootManager_->Add(dir, new TH3F("h_trolleyRail_3D","Position of steps that hit trolleyRail",500,-7400,7400,500,-7400,7400,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_xtal","Global radius vs theta for steps that hit xtal",28,-7,7,460,6800,7030));
  rootManager_->Add(dir, new TH2F("h_xtal_y","Global y vs radius for steps that hit xtal",460,6800,7030,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_xtal_beam","Beam radius vs theta for steps that hit xtal",28,-7,7,120,-140,-80));
  rootManager_->Add(dir, new TH3F("h_xtal_3D","Position of steps that hit xtal",500,-7200,-6500,500,-200,600,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_xtal_3D_2","Position of steps that hit xtal",500,-2200,0,500,-7400,-6200,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_StationNumber","Global radius vs theta for steps that hit StationNumber",28,-7,7,440,6800,7020));
  rootManager_->Add(dir, new TH2F("h_StationNumber_y","Global y vs radius for steps that hit StationNumber",440,6800,7020,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_StationNumber_beam","Beam radius vs theta for steps that hit StationNumber",28,-7,7,60,-120,-90));
  rootManager_->Add(dir, new TH3F("h_StationNumber_3D","Position of steps that hit StationNumber",500,-7200,-6600,500,0,2250,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_StationNumber_3D_2","Position of steps that hit StationNumber",500,-2200,0,500,-7400,-6400,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_PbF2Bounding","Global radius vs theta for steps that hit PbF2Bounding",28,-7,7,460,6800,7030));
  rootManager_->Add(dir, new TH2F("h_PbF2Bounding_y","Global y vs radius for steps that hit PbF2Bounding",460,6800,7030,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_PbF2Bounding_beam","Beam radius vs theta for steps that hit PbF2Bounding",28,-7,7,60,-120,-90));
  rootManager_->Add(dir, new TH3F("h_PbF2Bounding_3D","Position of steps that hit PbF2Bounding",500,-7200,-6500,500,0,2500,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_PbF2Bounding_3D_2","Position of steps that hit PbF2Bounding",500,-2400,0,500,-7400,-6200,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_Calorimeter","Global radius vs theta for steps that hit Calorimeter",28,-7,7,460,6800,7030));
  rootManager_->Add(dir, new TH2F("h_Calorimeter_y","Global y vs radius for steps that hit Calorimeter",460,6800,7030,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_Calorimeter_beam","Beam radius vs theta for steps that hit Calorimeter",28,-7,7,80,-120,-80));
  rootManager_->Add(dir, new TH3F("h_Calorimeter_3D","Position of steps that hit Calorimeter",500,-7100,-6500,500,0,2500,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_Calorimeter_3D_2","Position of steps that hit Calorimeter",500,-2200,0,500,-7400,-6400,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_insideCalo","Global radius vs theta for steps that hit insideCalo",28,-7,7,460,6800,7030));
  rootManager_->Add(dir, new TH2F("h_insideCalo_y","Global y vs radius for steps that hit insideCalo",460,6800,7030,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_insideCalo_beam","Beam radius vs theta for steps that hit insideCalo",28,-7,7,80,-120,-80));
  rootManager_->Add(dir, new TH3F("h_insideCalo_3D","Position of steps that hit insideCalo",500,-7100,-6600,500,-400,400,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_insideCalo_3D_2","Position of steps that hit insideCalo",500,-2600,200,500,-7400,-6400,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_photodetector","Global radius vs theta for steps that hit photodetector",28,-7,7,420,6800,7010));
  rootManager_->Add(dir, new TH2F("h_photodetector_y","Global y vs radius for steps that hit photodetector",420,6800,7010,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_photodetector_beam","Beam radius vs theta for steps that hit photodetector",28,-7,7,800,-200,200));
  rootManager_->Add(dir, new TH3F("h_photodetector_3D","Position of steps that hit photodetector",500,-7200,-6700,500,-400,500,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_photodetector_3D_2","Position of steps that hit photodetector",500,-400,0,500,-7400,-6400,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_outsideStorageRing","Global radius vs theta for steps that hit outsideStorageRing",28,-7,7,1200,6800,7400));
  rootManager_->Add(dir, new TH2F("h_outsideStorageRing_y","Global y vs radius for steps that hit outsideStorageRing",400,7000,7200,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_outsideStorageRing_beam","Beam radius vs theta for steps that hit outsideStorageRing",28,-7,7,800,-200,200));
  rootManager_->Add(dir, new TH3F("h_outsideStorageRing_3D","Position of steps that hit outsideStorageRing",500,-7400,-4000,500,-7400,7400,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_outsideStorageRing_3D_2","Position of steps that hit outsideStorageRing",500,-4000,7400,500,-7400,7400,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_frontWrapping","Global radius vs theta for steps that hit frontWrapping",28,-7,7,460,6800,7030));
  rootManager_->Add(dir, new TH2F("h_frontWrapping_y","Global y vs radius for steps that hit frontWrapping",460,6800,7030,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_frontWrapping_beam","Beam radius vs theta for steps that hit frontWrapping",28,-7,7,80,-120,-80));
  rootManager_->Add(dir, new TH3F("h_frontWrapping_3D","Position of steps that hit frontWrapping",500,-7100,-6500,500,200,600,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_frontWrapping_3D_2","Position of steps that hit frontWrapping",500,-600,-200,500,-7400,-6400,200,-100,100));

  rootManager_->Add(dir, new TH2F("h_backWrapping","Global radius vs theta for steps that hit backWrapping",28,-7,7,440,6800,7020));
  rootManager_->Add(dir, new TH2F("h_backWrapping_y","Global y vs radius for steps that hit backWrapping",440,6800,7020,400,-100,100));
  rootManager_->Add(dir, new TH2F("h_backWrapping_beam","Beam radius vs theta for steps that hit backWrapping",28,-7,7,60,-120,-90));
  rootManager_->Add(dir, new TH3F("h_backWrapping_3D","Position of steps that hit backWrapping",500,-7100,-6500,500,-200,600,200,-100,100));
  rootManager_->Add(dir, new TH3F("h_backWrapping_3D_2","Position of steps that hit backWrapping",500,-600,-200,500,-7400,-6400,200,-100,100));


}

void gm2strawtracker::GeaneExtrapolationAna::bookDecayVertexPlots(TDirectory* dir) 
{

  // book histograms
  std::string name = "";
  std::string title = "";

  //int runNumMin = 1000;
  //int runNumMax = 2000;
  //int nRuns = runNumMax - runNumMin;

  // book the same plots for all the various quality cuts: Momentum pass/fail, Hit volume (no vols hit) etc
  
  std::vector<std::string> subDirList = { "allStations/allEvents"
					  ,Form("allStations/pValue>%.3f_and_noVolumesHit",pValueCut_)
                                        };
  
  for (auto station : stations_) {
    subDirList.push_back(Form("station%02d/allEvents",station));
    subDirList.push_back(Form("station%02d/pValue>%.3f_and_noVolumesHit",station,pValueCut_));
  }
  
  for (auto subDirName : subDirList) {
    
    TDirectory* subDir = rootManager_->GetDir(dir,subDirName.c_str(),true);
    
    rootManager_->Add(subDir, new TH1F("h_numSteps"   ,"Number of Steps;Steps Taken;", 8000, 0, 8000));

    // radial, vertical and azimuthal pos
    rootManager_->Add(subDir, new TH1F("h_radialPos"  ,"Tangent point;Radial pos [mm];"  ,400,-200,200));
    rootManager_->Add(subDir, new TH1F("h_worldRadialPos"  ,"Tangent point;World Radial pos [mm];"  ,400,7000,7200));
    rootManager_->Add(subDir, new TH1F("h_verticalPos","Tangent point;Vertical pos [mm];",400,-200,200));
    rootManager_->Add(subDir, new TH1F("h_azimuth"    ,"Tangent point;theta [rad];"    ,1280,-6.4,6.4));
    rootManager_->Add(subDir, new TH1F("h_decayArcLength"    ,"Tangent point;Decay arc length [m];",250,0,10));
    rootManager_->Add(subDir, new TH2F("h_worldRadialPos_vs_theta", "worldRadialPos [mm]; Theta",30,3,6,400,7000,7200));
    rootManager_->Add(subDir, new TH2F("h_theta_vs_p","Tangent point;|p| [MeV/c];Theta",1000,0,3000,12,0,6));
    rootManager_->Add(subDir, new TH1F("h_radialPosDiff", "Tangent point;Radial pos[mm];", 100, -0.3, 0.3));
    rootManager_->Add(subDir, new TH1F("h_radialMomDiff", "Tangent point;Radial mom[MeV/c];", 500, -.00005, .00005));
    rootManager_->Add(subDir, new TH1F("h_radialMomDiffZoom", "Tangent point;Radial mom[MeV/c];", 500, -.00003, .000005));
    rootManager_->Add(subDir, new TH1F("h_radialMomDiffOriginZoom", "Tangent point;Radial mom[MeV/c];", 500, -.0000001, .0000001));
    rootManager_->Add(subDir, new TH2F("h_radialPosDiff_vs_radialMomDiff", "Tangent point;radialMomDiff[MeV/c];radialPosDiff[mm]",800,-.001,.0002,800,-.3,.1));
    rootManager_->Add(subDir, new TH1F("h_verticalPosDiff", "Tangent point;Vertical pos[mm];", 200, -.3, .3));
    rootManager_->Add(subDir, new TH1F("h_verticalMomDiff", "Tangent point;Vertical mom[MeV/c];", 280, -.002, .002));
    rootManager_->Add(subDir, new TH2F("h_verticalPosDiff_vs_verticalMomDiff", "Tangent point;Vertical mom[MeV/c];Vertical pos[mm];",400,-.002,.002,400,-.5,.5));
    rootManager_->Add(subDir, new TH2F("h_verticalPosDiff_vs_radialPosDiff", "Tangent point;Radial pos[mm];Vertical pos[mm];",800,-0.3,0.3,800,-0.3,0.3));

    // beam spot position profile
    rootManager_->Add(subDir, new TH2F("h_vertexPosSpread","Tangent point;Vertex pos_{radial} [mm] ;Vertex pos_{vertical} [mm]" ,200,-100,100,200,-75,75));
    rootManager_->Add(subDir, new TH1F("h_mom","Tangent point;|p| [MeV];",128,0,3200));
    rootManager_->Add(subDir, new TH1F("h_radialMom"  ,"Tangent point;p_{rad} [MeV];" ,200,-0.0000001,0.0000001));
    rootManager_->Add(subDir, new TH1F("h_verticalMom","Tangent point;p_{vert} [MeV];",100,-500,500));
    
    // momentum profile
    rootManager_->Add(subDir, new TH2F("h_vertexMomSpread","Tangent point;Vertex mom_{radial} [MeV];Vertex mom_{vertical} [MeV]",100,-0.1,0.1,100,-0.1,0.1));
    
    // y vs y' (e+)
    rootManager_->Add(subDir, new TH2F("h_y_yprime","y vs y';y [mm];y' [mm]",1000,-100,100,1000,-100,100));
    
    // world zx plot
    rootManager_->Add(subDir, new TH2F("h_zPos_vs_xPos","Tangent point;Vertex z pos [mm];Vertex x pos [mm]",750,-7500,7500,750,-7500,7500));
    
    // n events per run
    //rootManager_->Add(subDir, new TH1F("h_numTracksPerRun","N tracks per run;Run number;N tracks",nRuns,runNumMin,runNumMax));
    
    // mean time
    rootManager_->Add(subDir, new TH1F("h_meanTime",";Mean time [us];",6000,0,6000*0.148936));

    rootManager_->Add(subDir, new TH2F("h_radialPos_vs_time_fine"  ,"Radial pos vs mean time;Track time [us];Radial pos [mm]",20000,0,200,140,-70,70));
    rootManager_->Add(subDir, new TH2F("h_verticalPos_vs_time_fine"  ,"Vertical pos vs mean time;Track time [us];Vertical pos [mm]",20000,0,200,140,-70,70));
    rootManager_->Add(subDir, new TH2F("h_radialPos_vs_time"  ,"Radial pos vs mean time;Mean island time [us];Radial pos [mm]",6000,0,6000*0.148936,500,-250,250));
    rootManager_->Add(subDir, new TH2F("h_verticalPos_vs_time","Vertical pos vs mean time;Mean island time [us];Vertical pos [mm]",6000,0,6000*0.148936,500,-250,250));

    rootManager_->Add(subDir, new TH2F("h_radialMom_vs_time"  ,"Radial mom vs mean time;Mean island time [us];Radial mom [MeV]",6000,0,6000*0.148936,100,-0.1,0.1));
    rootManager_->Add(subDir, new TH2F("h_verticalMom_vs_time","Vertical mom vs mean time;Mean island time [us];Vertical mom [MeV]",6000,0,6000*0.148936,500,-100,100));

    // modulo
    double g2period = (TMath::TwoPi() / gm2consts_->omegaAMagic()) * 1e-3; // us 
    rootManager_->Add(subDir, new TH2F("h_radialPos_vs_time_mod"  ,"Radial pos vs mean time % g2period;Mean island time % g2period [us];Radial pos [mm]",300,0,g2period,500,-250,250));
    rootManager_->Add(subDir, new TH2F("h_verticalPos_vs_time_mod","Vertical pos vs mean time % g2period;Mean island time % g2period [us];Vertical pos [mm]",300,0,g2period,500,-250,250));
    
    rootManager_->Add(subDir, new TProfile("h_avgVerticalPos_vs_runNum","Average vertical pos vs run number;Run number;Avg vertical pos [mm]",3000,7000,10000,-200,200));
    rootManager_->Add(subDir, new TProfile("h_avgRadialPos_vs_runNum","Average radial pos vs run number;Run number;Avg radial pos [mm]",3000,7000,10000,-200,200));
    
    // arc length vs momentum
    rootManager_->Add(subDir, new TH2F("h_decayArcLength_vs_mom",";Momentum [MeV];Decay arc length [m]",128,0,3200,250,0,10));

  }
}

void gm2strawtracker::GeaneExtrapolationAna::bookGeomGraphs(TDirectory* dir) 
{

  TGraph* g_tracker0 = new TGraph;
  g_tracker0->SetName("strawGeometry_trackerCoords_station0");
  rootManager_->Add(dir, g_tracker0);
  
  TGraph* g_world0 = new TGraph;
  g_world0->SetName("strawGeometry_worldCoords_station0");
  rootManager_->Add(dir, g_world0);

  TGraph* g_tracker12 = new TGraph;
  g_tracker12->SetName("strawGeometry_trackerCoords_station12");
  rootManager_->Add(dir, g_tracker12);
  
  TGraph* g_world12 = new TGraph;
  g_world12->SetName("strawGeometry_worldCoords_station12");
  rootManager_->Add(dir, g_world12);

  TGraph* g_tracker18 = new TGraph;
  g_tracker18->SetName("strawGeometry_trackerCoords_station18");
  rootManager_->Add(dir, g_tracker18);
  
  TGraph* g_world18 = new TGraph;
  g_world18->SetName("strawGeometry_worldCoords_station18");
  rootManager_->Add(dir, g_world18);


}

void gm2strawtracker::GeaneExtrapolationAna::makeProfilePlots(TDirectory* dir) {
  
  std::vector<std::string> subDirList = { "allStations/allEvents"
					  ,Form("allStations/pValue>%.3f_and_noVolumesHit",pValueCut_)
                                        };
  
  for (auto station : stations_) {
    subDirList.push_back(Form("station%02d/allEvents",station));
    subDirList.push_back(Form("station%02d/pValue>%.3f_and_noVolumesHit",station,pValueCut_));
  }
  
  for (auto subDirName : subDirList) {
    
    TDirectory* subDir = rootManager_->GetDir(dir,subDirName.c_str(),true);
    TH2F* h_verticalPos_vs_time = rootManager_->Get<TH2F*>(subDir,"h_verticalPos_vs_time");
    TH2F* h_radialPos_vs_time   = rootManager_->Get<TH2F*>(subDir,"h_radialPos_vs_time");
    TH2F* h_verticalMom_vs_time = rootManager_->Get<TH2F*>(subDir,"h_verticalMom_vs_time");
    TH2F* h_radialMom_vs_time   = rootManager_->Get<TH2F*>(subDir,"h_radialMom_vs_time");
    TH2F* h_verticalPos_vs_time_mod = rootManager_->Get<TH2F*>(subDir,"h_verticalPos_vs_time_mod");
    TH2F* h_radialPos_vs_time_mod   = rootManager_->Get<TH2F*>(subDir,"h_radialPos_vs_time_mod");

    TProfile* tp_verticalPos_vs_time = (TProfile*)h_verticalPos_vs_time->ProfileX(); 
    TProfile* tp_radialPos_vs_time = (TProfile*)h_radialPos_vs_time->ProfileX(); 
    TProfile* tp_verticalMom_vs_time = (TProfile*)h_verticalMom_vs_time->ProfileX(); 
    TProfile* tp_radialMom_vs_time = (TProfile*)h_radialMom_vs_time->ProfileX(); 

    TProfile* tp_verticalPos_vs_time_mod = (TProfile*)h_verticalPos_vs_time_mod->ProfileX(); 
    TProfile* tp_radialPos_vs_time_mod = (TProfile*)h_radialPos_vs_time_mod->ProfileX(); 

    tp_verticalPos_vs_time->SetTitle("Tangent point"); 
    tp_verticalPos_vs_time->GetXaxis()->SetTitle("Time [us]"); 
    tp_verticalPos_vs_time->GetYaxis()->SetTitle("Average vertical pos [mm]"); 

    tp_radialPos_vs_time->SetTitle("Tangent point"); 
    tp_radialPos_vs_time->GetXaxis()->SetTitle("Time [us]"); 
    tp_radialPos_vs_time->GetYaxis()->SetTitle("Average radial pos [mm]"); 

    tp_verticalMom_vs_time->SetTitle("Tangent point"); 
    tp_verticalMom_vs_time->GetXaxis()->SetTitle("Time [us]");
    tp_verticalMom_vs_time->GetYaxis()->SetTitle("Average vertical mom [MeV]"); 

    tp_radialMom_vs_time->SetTitle("Tangent point"); 
    tp_radialMom_vs_time->GetXaxis()->SetTitle("Time [us]");
    tp_radialMom_vs_time->GetYaxis()->SetTitle("Average radial mom [MeV]");

    tp_verticalPos_vs_time_mod->SetTitle("Tangent point"); 
    tp_verticalPos_vs_time_mod->GetXaxis()->SetTitle("Time % g2period [us]"); 
    tp_verticalPos_vs_time_mod->GetYaxis()->SetTitle("Average vertical pos [mm]"); 

    tp_radialPos_vs_time_mod->SetTitle("Tangent point"); 
    tp_radialPos_vs_time_mod->GetXaxis()->SetTitle("Time % g2period [us]"); 
    tp_radialPos_vs_time_mod->GetYaxis()->SetTitle("Average radial pos [mm]"); 
    
    rootManager_->Add(subDir,tp_verticalPos_vs_time);
    rootManager_->Add(subDir,tp_radialPos_vs_time);
    rootManager_->Add(subDir,tp_verticalMom_vs_time);
    rootManager_->Add(subDir,tp_radialMom_vs_time);
    rootManager_->Add(subDir,tp_verticalPos_vs_time_mod);
    rootManager_->Add(subDir,tp_radialPos_vs_time_mod);
    
  }

}

void gm2strawtracker::GeaneExtrapolationAna::endJob()
{
  mf::LogInfo info(name_);
  info << "Enter GeaneExtrapolationAna::endJob\n";
  std::vector<std::string> topDirList = {"Extrapolation", "GeaneExtrapolation"};
  for (auto dirName_ : topDirList){
    TDirectory* topDir = rootManager_->GetDir(dirName_,true);
    TDirectory* vertexDir = rootManager_->GetDir(topDir,"vertices",true);
    makeProfilePlots(vertexDir);

    TDirectory* whichVolumes = rootManager_->GetDir(topDir,"whichVolumes",true);

    TH1F* h_volumesHit = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHit");
    TH1F* h_volumesHit_noVac = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHit_noVacuumChamber");
    TH1F* h_volumesHit_outsideWorld = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHit_outsideWorld");
    TH1F* h_volumesHitS0 = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHitS0");
    TH1F* h_volumesHit_noVacS0 = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHit_noVacuumChamberS0");
    TH1F* h_volumesHit_outsideWorldS0 = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHit_outsideWorldS0");
    TH1F* h_volumesHitS12 = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHitS12");
    TH1F* h_volumesHit_noVacS12 = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHit_noVacuumChamberS12");
    TH1F* h_volumesHit_outsideWorldS12 = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHit_outsideWorldS12");
    TH1F* h_volumesHitS18 = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHitS18");
    TH1F* h_volumesHit_noVacS18 = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHit_noVacuumChamberS18");
    TH1F* h_volumesHit_outsideWorldS18 = (TH1F*)rootManager_->Get(whichVolumes,"h_volumesHit_outsideWorldS18");
  
    int i(0);
    for (std::map<std::string,int>::iterator it=volumesHitByTracks.begin(); it!=volumesHitByTracks.end(); ++it) {
      h_volumesHit->SetBinContent(i+1,it->second);
      h_volumesHit->GetXaxis()->SetBinLabel(i+1,it->first.c_str());
      i++;
    }
    
    int j(0);
    for (std::map<std::string,int>::iterator it=volumesHitByNonVacTracks.begin(); it!=volumesHitByNonVacTracks.end(); ++it) {
      h_volumesHit_noVac->SetBinContent(j+1,it->second);
      h_volumesHit_noVac->GetXaxis()->SetBinLabel(j+1,it->first.c_str());
      j++;
    }
    
    j = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByOutsideWorldTracks.begin(); it!=volumesHitByOutsideWorldTracks.end(); ++it) {
      h_volumesHit_outsideWorld->SetBinContent(j+1,it->second);
      h_volumesHit_outsideWorld->GetXaxis()->SetBinLabel(j+1,it->first.c_str());
      j++;
    }

    //S0
    i = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByTracksS0.begin(); it!=volumesHitByTracksS0.end(); ++it) {
      h_volumesHitS0->SetBinContent(i+1,it->second);
      h_volumesHitS0->GetXaxis()->SetBinLabel(i+1,it->first.c_str());
      i++;
    }
  
    j = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByNonVacTracksS0.begin(); it!=volumesHitByNonVacTracksS0.end(); ++it) {
      h_volumesHit_noVacS0->SetBinContent(j+1,it->second);
      h_volumesHit_noVacS0->GetXaxis()->SetBinLabel(j+1,it->first.c_str());
      j++;
    }
  
    j = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByOutsideWorldTracksS0.begin(); it!=volumesHitByOutsideWorldTracksS0.end(); ++it) {
      h_volumesHit_outsideWorldS0->SetBinContent(j+1,it->second);
      h_volumesHit_outsideWorldS0->GetXaxis()->SetBinLabel(j+1,it->first.c_str());
      j++;
    }
  
    //S12
    i = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByTracksS12.begin(); it!=volumesHitByTracksS12.end(); ++it) {
      h_volumesHitS12->SetBinContent(i+1,it->second);
      h_volumesHitS12->GetXaxis()->SetBinLabel(i+1,it->first.c_str());
      i++;
    }
    
    j = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByNonVacTracksS12.begin(); it!=volumesHitByNonVacTracksS12.end(); ++it) {
      h_volumesHit_noVacS12->SetBinContent(j+1,it->second);
      h_volumesHit_noVacS12->GetXaxis()->SetBinLabel(j+1,it->first.c_str());
      j++;
    }
    
    j = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByOutsideWorldTracksS12.begin(); it!=volumesHitByOutsideWorldTracksS12.end(); ++it) {
      h_volumesHit_outsideWorldS12->SetBinContent(j+1,it->second);
      h_volumesHit_outsideWorldS12->GetXaxis()->SetBinLabel(j+1,it->first.c_str());
      j++;
    }
    
    //S18
    i = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByTracksS18.begin(); it!=volumesHitByTracksS18.end(); ++it) {
      h_volumesHitS18->SetBinContent(i+1,it->second);
      h_volumesHitS18->GetXaxis()->SetBinLabel(i+1,it->first.c_str());
      i++;
    }
  
    j = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByNonVacTracksS18.begin(); it!=volumesHitByNonVacTracksS18.end(); ++it) {
      h_volumesHit_noVacS18->SetBinContent(j+1,it->second);
      h_volumesHit_noVacS18->GetXaxis()->SetBinLabel(j+1,it->first.c_str());
      j++;
    }

    j = 0;
    for (std::map<std::string,int>::iterator it=volumesHitByOutsideWorldTracksS18.begin(); it!=volumesHitByOutsideWorldTracksS18.end(); ++it) {
      h_volumesHit_outsideWorldS18->SetBinContent(j+1,it->second);
      h_volumesHit_outsideWorldS18->GetXaxis()->SetBinLabel(j+1,it->first.c_str());
      j++;
    }
  }
  // fill failure mode histo
  auto dir = rootManager_->GetDir(dirName_,"failureModes",true);
  TH1F* h_failureModesBackwards = rootManager_->Get<TH1F*>(dir,"h_failureModesBackwards"); 
  TH1F* h_failureModesForwards  = rootManager_->Get<TH1F*>(dir,"h_failureModesForwards"); 
  
  int ibin(0);
  for (auto failureMode : failureModesBackwards_) {
    h_failureModesBackwards->SetBinContent(ibin+1,failureMode.second);
    h_failureModesBackwards->GetXaxis()->SetBinLabel(ibin+1,failureMode.first.c_str());
    ibin++;
  }

  int jbin(0);
  for (auto failureMode : failureModesForwards_) {
    h_failureModesForwards->SetBinContent(ibin+1,failureMode.second);
    h_failureModesForwards->GetXaxis()->SetBinLabel(jbin+1,failureMode.first.c_str());
    jbin++;
  }
  TDirectory* topDir = rootManager_->GetDir(dirName_,true);
  TDirectory* geomDir = rootManager_->GetDir(topDir,"strawGeometry",true);

  gm2strawtracker::StrawView view = StrawView::u_view;

  std::string graphName = "";
  
  for (auto station : stations_) {

    if (cs_.find(Form("TrackerStation[%d]",station)) == cs_.end()) continue;

    for (auto i : {0,1,2,3,4,5,6,7}){
      WireID wire0(station,i,view,0,0 );
      WireID wire31(station,i,view,0,31 );
      
      gm2geom::CoordSystem3Vector v0_w  = wire0.getCentreInWorld(cs_);
      gm2geom::CoordSystem3Vector v31_w = wire31.getCentreInWorld(cs_);
      
      gm2geom::CoordSystem3Vector v0_s  = v0_w.transform(cs_,  Form("TrackerStation[%d]",station));
      gm2geom::CoordSystem3Vector v31_s = v31_w.transform(cs_, Form("TrackerStation[%d]",station));
      
      graphName = Form("strawGeometry_trackerCoords_station%d",station);

      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->SetPoint(rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->GetN(), v0_s.x(), v0_s.z());
      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->SetPoint(rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->GetN(), v31_s.x(), v31_s.z());
      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->SetTitle(Form("Straw module locations, station %d",station));
      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->GetXaxis()->SetTitle("Station x pos [mm]");
      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->GetYaxis()->SetTitle("Station z pos [mm]");

      graphName = Form("strawGeometry_worldCoords_station%d",station);
      
      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->SetPoint(rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->GetN(), v0_w.z(), v0_w.x());
      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->SetPoint(rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->GetN(), v31_w.z(), v31_w.x());
      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->SetTitle(Form("Straw module locations, station %d",station));
      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->GetXaxis()->SetTitle("World z pos [mm]");
      rootManager_->Get<TGraph*>(geomDir,graphName.c_str())->GetYaxis()->SetTitle("World x pos [mm]");
    }
  }

  rootManager_->ClearEmpty(rootManager_->GetDir(dirName_));

  info  << "Exit GeaneExtrapolationAna::endJob\n";
}

DEFINE_ART_MODULE(gm2strawtracker::GeaneExtrapolationAna)
