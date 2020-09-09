
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

// Include needed ART headers
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

//Geant4
#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4PhysicalConstants.hh"

//art records
#include "gm2dataproducts/strawtracker/StrawDigitArtRecord.hh"
#include "gm2dataproducts/strawtracker/StrawTimeIslandArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackCandidateArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"
#include "gm2dataproducts/mc/ghostdetectors/GhostDetectorArtRecord.hh"

//Utils
#include "gm2geom/common/Gm2Constants_service.hh"
#include "gm2util/common/dataModuleDefs.hh"
#include "gm2util/common/RootManager.hh"

//C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <math.h> 

#include <Eigen/Dense>
#include "gm2tracker/utils/GeaneEigenStorageUtils.hh"

#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"
#include "gm2tracker/utils/GeaneTrackUtils.hh"
#include "gm2tracker/utils/GeaneDummyUtils.hh"

#include "TMath.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TLegend.h"


namespace gm2strawtracker {

//
// Class declaration
//
class GeaneSingleEventViewer : public art::EDAnalyzer {

  public:

    explicit GeaneSingleEventViewer(fhicl::ParameterSet const& pset);

    //Override desired art::EDAnalyzer functions
    void analyze(const art::Event& event ) override;
    void beginJob() override;
    void beginRun(art::Run const & r) override;
    void endJob() override;

  private:

    std::string name_;

    //Producer labels
    std::string TrackModuleLabel_;
    std::string TrackInstanceName_;

    std::string DummyModuleLabel_;
    std::string DummyInstanceName_;

    //ROOT plotting members
    std::unique_ptr<RootManager> rootManager_;
    std::string dirName_;

    bool specifyEvents_;
    std::vector<unsigned int> eventNum_;

    double energyLossCut_;

    //Helper tools
    gm2geom::StrawTrackerGeometry sgeom_;
    gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;
    gm2strawtracker::GeaneDummyUtils geaneDummyUtils_;
    gm2geom::CoordSystemsStoreData cs_;


    bool foundDummyCollection_;

    void BookTotalRootObjects(TDirectory* dir);
    void BookEventRootObjects(TDirectory* dir);
    void AnalyzeEvent(const art::Event& event, const gm2strawtracker::TrackDetailArtRecord& geaneTrack, int& trackNum, int& station);


}; //End of class GeaneSingleEventViewer


// Class implementation
GeaneSingleEventViewer::GeaneSingleEventViewer(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset)
    , name_( "GeaneSingleEventViewer" )
    , TrackModuleLabel_( pset.get<std::string>("TrackModuleLabel", "mainGEANEFit") )
    , TrackInstanceName_( pset.get<std::string>("TrackInstanceName", "mainGEANEtracks") )
    , DummyModuleLabel_( pset.get<std::string>("DummyModuleLabel","artg4") )
    , DummyInstanceName_( pset.get<std::string>("DummyInstanceName","trackerdummyplane") )
    , rootManager_()
    , dirName_( pset.get<std::string>("dirName","Top") )
    , specifyEvents_(pset.get<bool>("specifyEvents", false)) 
    , eventNum_(pset.get<std::vector<unsigned int>>("eventNum"))
    , energyLossCut_( pset.get<double>("energyLossCut", 1000000.))
    , sgeom_()
    , geaneTrackUtils_()
    , geaneDummyUtils_(pset)
    , cs_()
{ 
     // auto eventsDirectory = rootManager_->GetDir(dirName_,"IndividualEvents",true);
}

void GeaneSingleEventViewer::analyze(const art::Event& event) 
{
    mf::LogInfo info(name_);
    info << "Enter GeaneSingleEventViewer\n";

    //
    // Get data from art record
    //

    // Get the TrackDetailArtRecord
    art::Handle<gm2strawtracker::TrackDetailArtRecordCollection> TrackDataHandle;
    bool foundTrackCollection = event.getByLabel(TrackModuleLabel_,TrackInstanceName_,TrackDataHandle);
    if( !foundTrackCollection ) {
      throw cet::exception(name_) << "No Track Collection in this event (\"" << TrackModuleLabel_ << "\":\"" << TrackInstanceName_ << "\")\n";
      return;
    }

    // Get the track dummy planes art record
    art::Handle<gm2truth::GhostDetectorArtRecordCollection> DummyDataHandle;
    foundDummyCollection_ = event.getByLabel(DummyModuleLabel_,DummyInstanceName_,DummyDataHandle);
    if( !foundDummyCollection_ ) {
      mf::LogWarning(name_) << "No Dummy collection in this event (\"" << DummyModuleLabel_ << "\":\"" << DummyInstanceName_ << "\")\n";
    }

    //
    // Record stats
    //

    // Move to next event if we've not got any tracks in this event
    if( (*TrackDataHandle).size() == 0 ) return;

    // Apply event selection
    // FIXME -- doesn't currently work with more than 1 track per event
    if( specifyEvents_ ) {
      if (std::find(eventNum_.begin(), eventNum_.end(),  event.event()) == eventNum_.end() ){
        info << "Skipping event = " << event.event() << ", waiting for events: "; 
        for (auto i : eventNum_) info << i << ", ";
        info << "\n";
        return;
      }
    }

    info << "\tAnalyzing event " << event.event() << "\n";


    //
    // Make plots
    //

    int trackNum = 0;

    for(auto& track : *TrackDataHandle) {
      trackNum++;
      if (track.failureMode != 0 ) continue;

      int station = track.candidate->strawDigits.at(0)->wireID.getStation();

      info << "Event number: " << event.event() << " track number: " << trackNum << " with chi2 = " << track.chi2 << " pvalue = " << track.pValue << std::endl;

      AnalyzeEvent(event,track,trackNum,station);
    }

    info << "Exit GeaneSingleEventViewer\n";

}//analyze

 
// begin job 
void GeaneSingleEventViewer::beginJob() {

    //
    // Set up plotting helpers
    //
    
    // Create a ROOT file and manager 
    art::ServiceHandle<art::TFileService> tfs;
    auto& outputRootFile = tfs->file();
    rootManager_.reset( new RootManager(name_, &outputRootFile) ); 

    // Create directories
    auto topDir = rootManager_->GetDir(&outputRootFile,dirName_,true); //true -> create if doesn't exist

    // Book histograms
    auto totalDir = rootManager_->GetDir(topDir,"Total",true);
    BookTotalRootObjects(totalDir);

    auto eventsDirectory = rootManager_->GetDir(topDir,"IndividualEvents",true);

    if (specifyEvents_)
    {
      for (auto evt : eventNum_){
         stringstream ss;
         ss << evt;
         std::string fName = ss.str();

         auto eventDir = rootManager_->GetDir(eventsDirectory,fName,true);
         BookEventRootObjects(eventDir);
      }
    }

}//beginJob



// begin run
void GeaneSingleEventViewer::beginRun(art::Run const & r) {

  // Get coord systems
  cs_ = artg4::dataFromRunOrService<gm2geom::CoordSystemsStoreData, gm2geom::CoordSystemsStore>
        ( r, dataModuleDefs::coordSysModuleLabel(),dataModuleDefs::coordSysInstanceLabel() );

  if( cs_.size() == 0 ) {
    mf::LogWarning(name_) << "This run does not contain any data associated with the coordinate system\n";
  }

  geaneDummyUtils_.fillCS(cs_); // pass cs to dummy utils for coord sys transforms
}

// end job
void GeaneSingleEventViewer::endJob() {}

// book root objects
void GeaneSingleEventViewer::BookTotalRootObjects(TDirectory* totalDir) {

   auto allpullDir = rootManager_->GetDir(totalDir,"allpulls",true);
   rootManager_->Add( allpullDir, new TH1F( "1/P Truth Pull", "; delta(1/P)/sigma; Events", 100, -10, 10) );
   rootManager_->Add( allpullDir, new TH1F( "PuPz Truth Pull", "; delta(PuPz)/sigma; Events", 100, -10, 10) );
   rootManager_->Add( allpullDir, new TH1F( "PvPz Truth Pull", "; delta(PvPz)/sigma; Events", 100, -10, 10) );
   rootManager_->Add( allpullDir, new TH1F( "U Truth Pull", "; delta(U)/sigma; Events", 100, -10, 10) );
   rootManager_->Add( allpullDir, new TH1F( "V Truth Pull", "; delta(V)/sigma; Events", 100, -10, 10) );
}

// book root objects
void GeaneSingleEventViewer::BookEventRootObjects(TDirectory* eventDir) {

   auto resDir         = rootManager_->GetDir(eventDir,"residuals",true);
   auto momDir         = rootManager_->GetDir(eventDir,"momentums",true);
   auto posDir         = rootManager_->GetDir(eventDir,"positions",true);
   auto posDirUV       = rootManager_->GetDir(eventDir,"positionsUV",true);    
   auto posDirUVplanes = rootManager_->GetDir(eventDir,"positionsUVplanes",true);    
   auto otherDir       = rootManager_->GetDir(eventDir,"other",true);
   auto pullDir        = rootManager_->GetDir(eventDir,"pulls",true);

   TGraph* g       = new TGraph();
   TMultiGraph* mg = new TMultiGraph();


   /////////////////////////////////////////////////////////////////////////////////////
   //residuals
   /////////////////////////////////////////////////////////////////////////////////////

   g->SetLineColor(1);
   g->SetLineWidth(2);
   g->SetMarkerStyle(20);
   g->SetName("chi2PerPlane");
   g->SetTitle("chi2 vs Plane; planeNum; chi2");
   rootManager_->Add( eventDir, g->Clone());

   g->SetName("predtruthUVresidualPerPlane");
   g->SetTitle("UV predtruth residual vs Plane; planeNum; residual (UV mm)");
   rootManager_->Add( resDir, g->Clone());

   g->SetName("predmeasUVresidualPerPlane");
   g->SetTitle("UV predmeas residual vs Plane; planeNum; residual (UV mm)");
   rootManager_->Add( resDir, g->Clone());   

   g->SetName("predwireUVresidualPerPlane");
   g->SetTitle("UV predwire residual vs Plane; planeNum; residual (UV mm)");
   rootManager_->Add( resDir, g->Clone());    

   g->SetName("meastruthUVresidualPerPlane");
   g->SetTitle("UV meastruth residual vs Plane; planeNum; residual (UV mm)");
   rootManager_->Add( resDir, g->Clone());          

   g->SetName("predtruthZresidualPerPlane");
   g->SetTitle("X predtruth residual vs Plane; planeNum; residual (X mm)");
   rootManager_->Add( resDir, g->Clone());


   /////////////////////////////////////////////////////////////////////////////////////
   // momentum
   /////////////////////////////////////////////////////////////////////////////////////   

   g->SetLineColor(1);
   g->SetName("predMomPerPosition");
   g->SetTitle("Predicted Momentum vs Position; planeX (mm); momentum (MeV)");
   rootManager_->Add( momDir, g->Clone());

   g->SetName("predMomPerPlane");
   g->SetTitle("Predicted Momentum vs Plane; planeNum; momentum (MeV)");
   rootManager_->Add( momDir, g->Clone());

   g->SetLineColor(2);
   g->SetName("trueMomPerPosition");
   g->SetTitle("Truth Momentum vs Position; planeX (mm); momentum (MeV)");
   rootManager_->Add( momDir, g->Clone());

   g->SetName("trueMomPerPlane");
   g->SetTitle("Truth Momentum vs Plane; planeNum; momentum (MeV)");
   rootManager_->Add( momDir, g->Clone());


   mg->SetName("momPerPosition");
   mg->SetTitle("Pred and Truth Mom per positions; planeX (mm); momentum (MeV)");
   rootManager_->Add( momDir, mg->Clone());

   mg->SetName("momPerPlane");
   mg->SetTitle("Pred and Truth Mom per planes; planeNum; momentum (MeV)");
   rootManager_->Add( momDir, mg->Clone());

   /////////////////////////////////////////////////////////////////////////////////////
   // position
   /////////////////////////////////////////////////////////////////////////////////////

   g->SetLineColor(1);
   g->SetName("xyPosPred");
   g->SetTitle("Predicted xy positions; x (mm); y (mm)");
   rootManager_->Add( posDir, g->Clone());

   g->SetLineColor(2);
   g->SetName("xyPosTruth");
   g->SetTitle("Truth xy positions; x (mm); y (mm)");
   rootManager_->Add( posDir, g->Clone());

   mg->SetName("xyPos");
   mg->SetTitle("xy positions; x (mm); y (mm)");
   rootManager_->Add( posDir, mg->Clone());

   g->SetLineColor(1);
   g->SetName("xzPosPred");
   g->SetTitle("Predicted xz positionss; x (mm); z (mm)");
   rootManager_->Add( posDir, g->Clone());

   g->SetLineColor(2);
   g->SetName("xzPosTruth");
   g->SetTitle("Truth xz positions; x (mm); z (mm)");
   rootManager_->Add( posDir, g->Clone());

   mg->SetName("xzPos");
   mg->SetTitle("xz positions; x (mm); z (mm)");
   rootManager_->Add( posDir, mg->Clone());

   g->SetLineColor(1);
   g->SetName("zyPosPred");
   g->SetTitle("Predicted zy positions; z (mm); y (mm)");
   rootManager_->Add( posDir, g->Clone());

   g->SetLineColor(2);
   g->SetName("zyPosTruth");
   g->SetTitle("Truth zy positions; z (mm); y (mm)");
   rootManager_->Add( posDir, g->Clone());

   mg->SetName("zyPos");
   mg->SetTitle("zy positions; z (mm); y (mm)");
   rootManager_->Add( posDir, mg->Clone());

   /////////////////////////////////////////////////////////////////////////////////////
   // position UV
   /////////////////////////////////////////////////////////////////////////////////////

   g->SetLineColor(1);
   g->SetName("zuPosPred");
   g->SetTitle("Predicted zu positions; z (mm); u (mm)");
   rootManager_->Add( posDirUV, g->Clone());

   g->SetLineColor(2);
   g->SetName("zuPosTruth");
   g->SetTitle("Truth zu positions; z (mm); u (mm)");
   rootManager_->Add( posDirUV, g->Clone());

   g->SetLineColor(4);
   g->SetName("zuPosMeas");
   g->SetTitle("Meas zu positions; z (mm); u (mm)");
   rootManager_->Add( posDirUV, g->Clone());

   g->SetLineColor(6);
   g->SetName("zuPosWire");
   g->SetTitle("Wire zu positions; z (mm); u (mm)");
   rootManager_->Add( posDirUV, g->Clone());

   mg->SetName("zuPos");
   mg->SetTitle("zu positions; z (mm); u (mm)");
   rootManager_->Add( posDirUV, mg->Clone());


   g->SetLineColor(1);
   g->SetName("zvPosPred");
   g->SetTitle("Predicted zv positionss; z (mm); v (mm)");
   rootManager_->Add( posDirUV, g->Clone());

   g->SetLineColor(2);
   g->SetName("zvPosTruth");
   g->SetTitle("Truth zv positions; z (mm); v (mm)");
   rootManager_->Add( posDirUV, g->Clone());

   g->SetLineColor(4);
   g->SetName("zvPosMeas");
   g->SetTitle("Meas zv positions; z (mm); v (mm)");
   rootManager_->Add( posDirUV, g->Clone());

   g->SetLineColor(6);
   g->SetName("zvPosWire");
   g->SetTitle("Wire zv positions; z (mm); v (mm)");
   rootManager_->Add( posDirUV, g->Clone());

   mg->SetName("zvPos");
   mg->SetTitle("zv positions; z (mm); v (mm)");
   rootManager_->Add( posDirUV, mg->Clone());


   g->SetLineColor(1);
   g->SetName("zuPosPredplanes");
   g->SetTitle("Predicted zu positions; planeNum; u (mm)");
   rootManager_->Add( posDirUVplanes, g->Clone());

   g->SetLineColor(2);
   g->SetName("zuPosTruthplanes");
   g->SetTitle("Truth zu positions; planeNum; u (mm)");
   rootManager_->Add( posDirUVplanes, g->Clone());

   g->SetLineColor(4);
   g->SetName("zuPosMeasplanes");
   g->SetTitle("Meas zu positions; planeNum; u (mm)");
   rootManager_->Add( posDirUVplanes, g->Clone());

   g->SetLineColor(6);
   g->SetName("zuPosWireplanes");
   g->SetTitle("Wire zu positions; planeNum; u (mm)");
   rootManager_->Add( posDirUVplanes, g->Clone());

   mg->SetName("zuPosplanes");
   mg->SetTitle("zu positions; planeNum; u (mm)");
   rootManager_->Add( posDirUVplanes, mg->Clone());


   g->SetLineColor(1);
   g->SetName("zvPosPredplanes");
   g->SetTitle("Predicted zv positionss; planeNum; v (mm)");
   rootManager_->Add( posDirUVplanes, g->Clone());

   g->SetLineColor(2);
   g->SetName("zvPosTruthplanes");
   g->SetTitle("Truth zv positions; planeNum; v (mm)");
   rootManager_->Add( posDirUVplanes, g->Clone());

   g->SetLineColor(4);
   g->SetName("zvPosMeasplanes");
   g->SetTitle("Meas zv positions; planeNum; v (mm)");
   rootManager_->Add( posDirUVplanes, g->Clone());

   g->SetLineColor(6);
   g->SetName("zvPosWireplanes");
   g->SetTitle("Wire zv positions; planeNum; v (mm)");
   rootManager_->Add( posDirUVplanes, g->Clone());

   mg->SetName("zvPosplanes");
   mg->SetTitle("zv positions; planeNum; v (mm)");
   rootManager_->Add( posDirUVplanes, mg->Clone());

   /////////////////////////////////////////////////////////////////////////////////////
   // other
   /////////////////////////////////////////////////////////////////////////////////////

   g->SetName("dcaMeas");
   g->SetTitle("Measured dca; planeNum; dca (mm)");
   rootManager_->Add( otherDir, g->Clone());

   g->SetName("UVerrors");
   g->SetTitle("UV errors; planeNum; error (mm)");
   rootManager_->Add( otherDir, g->Clone());

   g->SetName("correctLRChoicePerPlane");
   g->SetTitle("correctLRChoicePerPlane; planeNum; correct or not");
   rootManager_->Add( otherDir, g->Clone());

   g->SetName("LRchoiceTimesDCA");
   g->SetTitle("LRchoiceTimesDCA; planeNum; correct or not * dca");
   rootManager_->Add( otherDir, g->Clone());     

   /////////////////////////////////////////////////////////////////////////////////////
   // pulls
   /////////////////////////////////////////////////////////////////////////////////////

   rootManager_->Add( pullDir, new TH1F( "1/P Truth Pull", "; delta(1/P)/sigma; Events", 100, -10, 10) );
   rootManager_->Add( pullDir, new TH1F( "PuPz Truth Pull", "; delta(PuPz)/sigma; Events", 100, -10, 10) );
   rootManager_->Add( pullDir, new TH1F( "PvPz Truth Pull", "; delta(PvPz)/sigma; Events", 100, -10, 10) );
   rootManager_->Add( pullDir, new TH1F( "U Truth Pull", "; delta(U)/sigma; Events", 100, -10, 10) );
   rootManager_->Add( pullDir, new TH1F( "V Truth Pull", "; delta(V)/sigma; Events", 100, -10, 10) );
}


// analyze event
void GeaneSingleEventViewer::AnalyzeEvent(const art::Event& event, const gm2strawtracker::TrackDetailArtRecord& geaneTrack, int& trackNum, int& station) {

   mf::LogInfo info(name_);

   //Resolve handle to get collection - now with eigen object conversion
   auto track = GeaneEigenStorageUtils::ReadEigenFromStorage(geaneTrack); 
   auto geaneHitsOnTrack = track.geaneHits;

   // Track dummy Collection
   if(foundDummyCollection_) {
     double energyDiff = track.dummyPlaneHits.at(0)->momentum.mag() - track.dummyPlaneHits.back()->momentum.mag();
     if (energyDiff > energyLossCut_){
        info << "LOSSCATCH  event " << event.event() << " with greater energy loss than momentum tolerance: " << energyDiff << "\n";
        return;
     } 
   }

   auto eventsDirectory = rootManager_->GetDir(dirName_,"IndividualEvents",true);

   stringstream ss;
   ss << event.event();
   if(!specifyEvents_) ss << "-" << event.subRun() << "-" << event.run() << "-" << trackNum;
   std::string fName = ss.str();

   auto eventDir = rootManager_->GetDir(eventsDirectory,fName,true);

   auto resDir         = rootManager_->GetDir(eventDir,"residuals",true);
   auto momDir         = rootManager_->GetDir(eventDir,"momentums",true);
   auto posDir         = rootManager_->GetDir(eventDir,"positions",true);
   auto posDirUV       = rootManager_->GetDir(eventDir,"positionsUV",true);    
   auto posDirUVplanes = rootManager_->GetDir(eventDir,"positionsUVplanes",true);    
   auto otherDir       = rootManager_->GetDir(eventDir,"other",true);
   auto pullDir        = rootManager_->GetDir(eventDir,"pulls",true);


   if(!specifyEvents_) BookEventRootObjects(eventDir);

   stringstream stationStream;
   stationStream << "TrackerStation[" << station << "]";
   string stationStr = stationStream.str();

   if(foundDummyCollection_){
     geaneDummyUtils_.fillStationStr(stationStr, station); // pass station number to dummy utils
   }

   for (int i = 0; i < int(track.trackPlanesHitList.size()); ++i)
   {
      int planeNum = track.trackPlanesHitList.at(i);

      double truthMom; // have to declare these outside the dummy check statement
      double xPosTruth; 
      double yPosTruth; 
      double zPosTruth; 
      double uPosTruth; 
      double vPosTruth; 

      if(foundDummyCollection_){
        auto dummyHit = track.dummyPlaneHits.at(i+1);

        gm2geom::CoordSystem3Vector pHitPos(dummyHit->position.x(), dummyHit->position.y(), dummyHit->position.z() , "world");
        gm2geom::CoordSystem3Vector planeIntersection = pHitPos.transform(cs_, stationStr);

        gm2geom::CoordSystem3Vector pHitMom(dummyHit->momentum.x(), dummyHit->momentum.y(), dummyHit->momentum.z() , "world");
        gm2geom::CoordSystem3Vector momGlob = pHitMom.transform(cs_, stationStr, true);

        truthMom  = momGlob.mag();
        xPosTruth = planeIntersection.x();
        yPosTruth = planeIntersection.y();
        zPosTruth = planeIntersection.z();

        uPosTruth = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*xPosTruth + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,1)*yPosTruth;
        vPosTruth = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,0)*xPosTruth + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,1)*yPosTruth;
      }

      double zPosMeas = geaneHitsOnTrack.planeZPositions.at(planeNum);
      double uPosMeas = geaneHitsOnTrack.geaneMeasuredParameters[3][planeNum];
      double vPosMeas = geaneHitsOnTrack.geaneMeasuredParameters[4][planeNum];

      double uPosWire = geaneHitsOnTrack.geaneWireUVPositions[3][planeNum];
      double vPosWire = geaneHitsOnTrack.geaneWireUVPositions[4][planeNum];

      double predMom  = geaneTrackUtils_.getPredMom(geaneHitsOnTrack, planeNum);
      double xPosPred = geaneTrackUtils_.getPredXPos(geaneHitsOnTrack, planeNum);
      double yPosPred = geaneTrackUtils_.getPredYPos(geaneHitsOnTrack, planeNum);
      double zPosPred = geaneTrackUtils_.getPredZPos(geaneHitsOnTrack, planeNum);
      double uPosPred = geaneTrackUtils_.getPredUPos(geaneHitsOnTrack, planeNum);
      double vPosPred = geaneTrackUtils_.getPredVPos(geaneHitsOnTrack, planeNum);

      rootManager_->Get<TGraph*>( eventDir, "chi2PerPlane")->SetPoint(rootManager_->Get<TGraph*>( eventDir, "chi2PerPlane")->GetN(), planeNum , track.chi2Planes.at(i));
      if(track.chi2Planes.at(i) < 0) info << "Negative plane chi2 for event: " << event.event() << " with chi2: " << track.chi2Planes.at(i) << "\n";

          
      double predtruthUVresidual;
      double predmeasUVresidual;
      double predwireUVresidual;
      double meastruthUVresidual;

      if (geaneTrackUtils_.isUPlane(planeNum))
      {
        predmeasUVresidual = uPosPred - uPosMeas;
        predwireUVresidual = uPosPred - geaneHitsOnTrack.geaneWireUVPositions[3][planeNum];

        if(foundDummyCollection_){
          predtruthUVresidual = uPosPred - uPosTruth;
          meastruthUVresidual = uPosMeas - uPosTruth;               
        }

      } 
      else
      {
        predmeasUVresidual = vPosPred - vPosMeas;
        predwireUVresidual = vPosPred - geaneHitsOnTrack.geaneWireUVPositions[4][planeNum];

        if(foundDummyCollection_){
          predtruthUVresidual = vPosPred - vPosTruth;
          meastruthUVresidual = vPosMeas - vPosTruth;
        }
      } 

          
      rootManager_->Get<TGraph*>( resDir, "predmeasUVresidualPerPlane")->SetPoint(rootManager_->Get<TGraph*>( resDir, "predmeasUVresidualPerPlane")->GetN(), planeNum , predmeasUVresidual);
      rootManager_->Get<TGraph*>( resDir, "predwireUVresidualPerPlane")->SetPoint(rootManager_->Get<TGraph*>( resDir, "predwireUVresidualPerPlane")->GetN(), planeNum , predwireUVresidual);
                    
      if(foundDummyCollection_){
        rootManager_->Get<TGraph*>( resDir, "predtruthUVresidualPerPlane")->SetPoint(rootManager_->Get<TGraph*>( resDir, "predtruthUVresidualPerPlane")->GetN(), planeNum , predtruthUVresidual);
        rootManager_->Get<TGraph*>( resDir, "meastruthUVresidualPerPlane")->SetPoint(rootManager_->Get<TGraph*>( resDir, "meastruthUVresidualPerPlane")->GetN(), planeNum , meastruthUVresidual);

        rootManager_->Get<TGraph*>( resDir, "predtruthZresidualPerPlane")->SetPoint(rootManager_->Get<TGraph*>( resDir, "predtruthZresidualPerPlane")->GetN(), planeNum , zPosPred - zPosTruth);
        if(abs(zPosPred - zPosTruth) > .00001) info << "Z residual very large for event: " << event.event() << " trackNum: " << trackNum << "\n";
      }


      rootManager_->Get<TGraph*>( momDir, "predMomPerPosition")->SetPoint(rootManager_->Get<TGraph*>( momDir, "predMomPerPosition")->GetN(), zPosPred , predMom);
      rootManager_->Get<TGraph*>( momDir, "predMomPerPlane")->SetPoint(rootManager_->Get<TGraph*>( momDir, "predMomPerPlane")->GetN(), planeNum , predMom);

      rootManager_->Get<TGraph*>( posDir, "xyPosPred")->SetPoint(rootManager_->Get<TGraph*>( posDir, "xyPosPred")->GetN(), xPosPred , yPosPred);
      rootManager_->Get<TGraph*>( posDir, "xzPosPred")->SetPoint(rootManager_->Get<TGraph*>( posDir, "xzPosPred")->GetN(), xPosPred , zPosPred);
      rootManager_->Get<TGraph*>( posDir, "zyPosPred")->SetPoint(rootManager_->Get<TGraph*>( posDir, "zyPosPred")->GetN(), zPosPred , yPosPred);


      if(foundDummyCollection_){          
        rootManager_->Get<TGraph*>( momDir, "trueMomPerPosition")->SetPoint(rootManager_->Get<TGraph*>( momDir, "trueMomPerPosition")->GetN(), zPosTruth , truthMom);
        rootManager_->Get<TGraph*>( momDir, "trueMomPerPlane")->SetPoint(rootManager_->Get<TGraph*>( momDir, "trueMomPerPlane")->GetN(), planeNum , truthMom);

        rootManager_->Get<TGraph*>( posDir, "xyPosTruth")->SetPoint(rootManager_->Get<TGraph*>( posDir, "xyPosTruth")->GetN(), xPosTruth , yPosTruth);
        rootManager_->Get<TGraph*>( posDir, "xzPosTruth")->SetPoint(rootManager_->Get<TGraph*>( posDir, "xzPosTruth")->GetN(), xPosTruth , zPosTruth);
        rootManager_->Get<TGraph*>( posDir, "zyPosTruth")->SetPoint(rootManager_->Get<TGraph*>( posDir, "zyPosTruth")->GetN(), zPosTruth , yPosTruth);
      }


      if (geaneTrackUtils_.isUPlane(planeNum))          
      {
        rootManager_->Get<TGraph*>( posDirUV, "zuPosPred")->SetPoint(rootManager_->Get<TGraph*>( posDirUV, "zuPosPred")->GetN(), zPosPred , uPosPred);
        rootManager_->Get<TGraph*>( posDirUV, "zuPosMeas")->SetPoint(rootManager_->Get<TGraph*>( posDirUV, "zuPosMeas")->GetN(), zPosMeas , uPosMeas);
        rootManager_->Get<TGraph*>( posDirUV, "zuPosWire")->SetPoint(rootManager_->Get<TGraph*>( posDirUV, "zuPosWire")->GetN(), zPosMeas , uPosWire);
        if(foundDummyCollection_) rootManager_->Get<TGraph*>( posDirUV, "zuPosTruth")->SetPoint(rootManager_->Get<TGraph*>( posDirUV, "zuPosTruth")->GetN(), zPosTruth , uPosTruth);
      }
      else
      {
        rootManager_->Get<TGraph*>( posDirUV, "zvPosPred")->SetPoint(rootManager_->Get<TGraph*>( posDirUV, "zvPosPred")->GetN(), zPosPred , vPosPred);
        rootManager_->Get<TGraph*>( posDirUV, "zvPosMeas")->SetPoint(rootManager_->Get<TGraph*>( posDirUV, "zvPosMeas")->GetN(), zPosMeas , vPosMeas);
        rootManager_->Get<TGraph*>( posDirUV, "zvPosWire")->SetPoint(rootManager_->Get<TGraph*>( posDirUV, "zvPosWire")->GetN(), zPosMeas , vPosWire);
        if(foundDummyCollection_) rootManager_->Get<TGraph*>( posDirUV, "zvPosTruth")->SetPoint(rootManager_->Get<TGraph*>( posDirUV, "zvPosTruth")->GetN(), zPosTruth , vPosTruth);
      }


      if (geaneTrackUtils_.isUPlane(planeNum))          
      {
        rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosPredplanes")->SetPoint(rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosPredplanes")->GetN(), planeNum , uPosPred);
        rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosMeasplanes")->SetPoint(rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosMeasplanes")->GetN(), planeNum , uPosMeas);
        rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosWireplanes")->SetPoint(rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosWireplanes")->GetN(), planeNum , uPosWire);
        if(foundDummyCollection_) rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosTruthplanes")->SetPoint(rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosTruthplanes")->GetN(), planeNum , uPosTruth);
      }
      else
      {
        rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosPredplanes")->SetPoint(rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosPredplanes")->GetN(), planeNum , vPosPred);
        rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosMeasplanes")->SetPoint(rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosMeasplanes")->GetN(), planeNum , vPosMeas);
        rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosWireplanes")->SetPoint(rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosWireplanes")->GetN(), planeNum , vPosWire);
        if(foundDummyCollection_) rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosTruthplanes")->SetPoint(rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosTruthplanes")->GetN(), planeNum , vPosTruth);
      }


      rootManager_->Get<TGraph*>( otherDir, "dcaMeas")->SetPoint(rootManager_->Get<TGraph*>( otherDir, "dcaMeas")->GetN(), planeNum , track.measuredDCAs.at(planeNum));
      rootManager_->Get<TGraph*>( otherDir, "UVerrors")->SetPoint(rootManager_->Get<TGraph*>( otherDir, "UVerrors")->GetN(), planeNum , track.UVerrors.at(planeNum));

   } // for loop over hit planes
        
   rootManager_->Get<TMultiGraph*>(posDir, "xyPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDir, "xyPosPred")->Clone());
   rootManager_->Get<TMultiGraph*>(posDir, "xzPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDir, "xzPosPred")->Clone());
   rootManager_->Get<TMultiGraph*>(posDir, "zyPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDir, "zyPosPred")->Clone());

   if(foundDummyCollection_){
     rootManager_->Get<TMultiGraph*>(posDir, "xyPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDir, "xyPosTruth")->Clone());
     rootManager_->Get<TMultiGraph*>(posDir, "xzPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDir, "xzPosTruth")->Clone());
     rootManager_->Get<TMultiGraph*>(posDir, "zyPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDir, "zyPosTruth")->Clone());
   }


   rootManager_->Get<TMultiGraph*>(posDirUV, "zuPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUV, "zuPosPred")->Clone());
   rootManager_->Get<TMultiGraph*>(posDirUV, "zuPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUV, "zuPosMeas")->Clone());
   rootManager_->Get<TMultiGraph*>(posDirUV, "zuPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUV, "zuPosWire")->Clone());
   if(foundDummyCollection_) rootManager_->Get<TMultiGraph*>(posDirUV, "zuPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUV, "zuPosTruth")->Clone());


   rootManager_->Get<TMultiGraph*>(posDirUV, "zvPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUV, "zvPosPred")->Clone());
   rootManager_->Get<TMultiGraph*>(posDirUV, "zvPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUV, "zvPosMeas")->Clone());
   rootManager_->Get<TMultiGraph*>(posDirUV, "zvPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUV, "zvPosWire")->Clone());
   if(foundDummyCollection_) rootManager_->Get<TMultiGraph*>(posDirUV, "zvPos")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUV, "zvPosTruth")->Clone());


   rootManager_->Get<TMultiGraph*>(posDirUVplanes, "zuPosplanes")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosPredplanes")->Clone());
   rootManager_->Get<TMultiGraph*>(posDirUVplanes, "zuPosplanes")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosMeasplanes")->Clone());
   rootManager_->Get<TMultiGraph*>(posDirUVplanes, "zuPosplanes")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosWireplanes")->Clone());
   if(foundDummyCollection_) rootManager_->Get<TMultiGraph*>(posDirUVplanes, "zuPosplanes")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUVplanes, "zuPosTruthplanes")->Clone());


   rootManager_->Get<TMultiGraph*>(posDirUVplanes, "zvPosplanes")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosPredplanes")->Clone());
   rootManager_->Get<TMultiGraph*>(posDirUVplanes, "zvPosplanes")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosMeasplanes")->Clone());
   rootManager_->Get<TMultiGraph*>(posDirUVplanes, "zvPosplanes")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosWireplanes")->Clone());
   if(foundDummyCollection_) rootManager_->Get<TMultiGraph*>(posDirUVplanes, "zvPosplanes")->Add((TGraph*)rootManager_->Get<TGraph*>( posDirUVplanes, "zvPosTruthplanes")->Clone());



   rootManager_->Get<TMultiGraph*>(momDir, "momPerPosition")->Add((TGraph*)rootManager_->Get<TGraph*>( momDir, "predMomPerPosition")->Clone());
   if(foundDummyCollection_) rootManager_->Get<TMultiGraph*>(momDir, "momPerPosition")->Add((TGraph*)rootManager_->Get<TGraph*>( momDir, "trueMomPerPosition")->Clone());

   rootManager_->Get<TMultiGraph*>(momDir, "momPerPlane")->Add((TGraph*)rootManager_->Get<TGraph*>( momDir, "predMomPerPlane")->Clone());
   if(foundDummyCollection_) rootManager_->Get<TMultiGraph*>(momDir, "momPerPlane")->Add((TGraph*)rootManager_->Get<TGraph*>( momDir, "trueMomPerPlane")->Clone());


   if(foundDummyCollection_){
     auto dummyHit = track.dummyPlaneHits.at(0);

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


     double predicted0Momentum  = geaneTrackUtils_.getPredMom(geaneHitsOnTrack, 0);
     double predicted0XMomentum = geaneTrackUtils_.getPredXMom(geaneHitsOnTrack, 0);
     double predicted0YMomentum = geaneTrackUtils_.getPredYMom(geaneHitsOnTrack, 0);
     double predicted0ZMomentum = geaneTrackUtils_.getPredZMom(geaneHitsOnTrack, 0);
     double predicted0UMomentum = geaneTrackUtils_.getPredUMom(geaneHitsOnTrack, 0);
     double predicted0VMomentum = geaneTrackUtils_.getPredVMom(geaneHitsOnTrack, 0);
     double predicted0YPosition = geaneTrackUtils_.getPredYPos(geaneHitsOnTrack, 0);
     double predicted0ZPosition = geaneTrackUtils_.getPredZPos(geaneHitsOnTrack, 0);
     double predicted0UPosition = geaneTrackUtils_.getPredUPos(geaneHitsOnTrack, 0);
     double predicted0VPosition = geaneTrackUtils_.getPredVPos(geaneHitsOnTrack, 0);

     Eigen::MatrixXd true0Covariance = geaneHitsOnTrack.covarianceTotalInverse.inverse(); // covarianceTotal is inverted from what it should be in the main code for speed improvements

     // errors are from the sigma0 covariance matrix 
     rootManager_->Get<TH1F*>( pullDir, "1/P Truth Pull" )->Fill((1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0)))); // using diagonal elements as the errors - should I be using a sum over columns? (for each row)
     rootManager_->Get<TH1F*>( pullDir, "PuPz Truth Pull" )->Fill((predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1))); 
     rootManager_->Get<TH1F*>( pullDir, "PvPz Truth Pull" )->Fill((predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2))); 
     rootManager_->Get<TH1F*>( pullDir, "U Truth Pull" )->Fill((predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3)))); 
     rootManager_->Get<TH1F*>( pullDir, "V Truth Pull" )->Fill((predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4)))); 

     auto topDir     = rootManager_->GetDir(dirName_,true);
     auto totalDir   = rootManager_->GetDir(topDir,"Total",true);
     auto allpullDir = rootManager_->GetDir(totalDir,"allpulls",true);

     rootManager_->Get<TH1F*>( allpullDir, "1/P Truth Pull" )->Fill((1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0)))); // using diagonal elements as the errors - should I be using a sum over columns? (for each row)
     rootManager_->Get<TH1F*>( allpullDir, "PuPz Truth Pull" )->Fill((predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1))); 
     rootManager_->Get<TH1F*>( allpullDir, "PvPz Truth Pull" )->Fill((predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2))); 
     rootManager_->Get<TH1F*>( allpullDir, "U Truth Pull" )->Fill((predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3)))); 
     rootManager_->Get<TH1F*>( allpullDir, "V Truth Pull" )->Fill((predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4)))); 
   } // foundDummyCollection 0 plane


   // check LR sides against truth
   if(foundDummyCollection_){
     geaneDummyUtils_.fillLRFromTruth(track);
     std::vector<int> wrongHitSides = geaneDummyUtils_.checkLRAgainstTruth(track);

     for (int i = 0; i < int(wrongHitSides.size()); ++i)
     {
        rootManager_->Get<TGraph*>( otherDir, "correctLRChoicePerPlane")->SetPoint(rootManager_->Get<TGraph*>( otherDir, "correctLRChoicePerPlane")->GetN(), i , wrongHitSides.at(i));
     }

     for (int i = 0; i < int(wrongHitSides.size()); ++i)
     {
        rootManager_->Get<TGraph*>( otherDir, "LRchoiceTimesDCA")->SetPoint(rootManager_->Get<TGraph*>( otherDir, "LRchoiceTimesDCA")->GetN(), i , wrongHitSides.at(i)*track.measuredDCAs.at(i));
     }
  }

  return;

} // end AnalyzeEvent method


} // End of namespace gm2strawtracker



//
// Extras
//

// These are some necessary boilerplate for the ROOT persistency system
using gm2strawtracker::GeaneSingleEventViewer;
DEFINE_ART_MODULE(GeaneSingleEventViewer)
