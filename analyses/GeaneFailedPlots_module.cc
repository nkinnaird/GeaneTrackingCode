//Plots for events that fail Geane fitting somehow.

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
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"
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

namespace gm2strawtracker {

  //
  // Class declaration
  //
  class GeaneFailedPlots : public art::EDAnalyzer {

  public:

    explicit GeaneFailedPlots(fhicl::ParameterSet const& pset);

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

    //Define what data products to summarise


    //ROOT plotting members
    std::unique_ptr<RootManager> rootManager_;
    std::string dirName_;

     //Helper tools
    gm2geom::CoordSystemsStoreData cs_;

    //Stats

    void BookHistograms(TDirectory* dir);
    void FillHistograms(const art::Event& event, art::Handle<gm2strawtracker::TrackArtRecordCollection> tracksHandle, art::Handle<gm2truth::GhostDetectorArtRecordCollection> dummyHandle, TDirectory* trackDir);

  }; //End of class GeaneFailedPlots


  //
  // Class implementation
  //

  GeaneFailedPlots::GeaneFailedPlots(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset)
    , name_( "GeaneFailedPlots" )
    , TrackModuleLabel_( pset.get<std::string>("TrackModuleLabel","") )  
    , TrackInstanceName_( pset.get<std::string>("TrackInstanceName","") )
    , DummyModuleLabel_( pset.get<std::string>("DummyModuleLabel","artg4") )
    , DummyInstanceName_( pset.get<std::string>("DummyInstanceName","trackerdummyplane") )
    , rootManager_()
    , dirName_( pset.get<std::string>("dirName","Failed") )
    , cs_()
  {}
  

  void GeaneFailedPlots::analyze(const art::Event& event) {
    
    //Get pointers to TDirectories for use later
    auto topDir = rootManager_->GetDir(dirName_);

    //
    // Get data from art record
    //

    //Get TrackArtRecordCollection
      art::Handle<gm2strawtracker::TrackArtRecordCollection> TrackDataHandle;
      bool foundTrackcollection = event.getByLabel(TrackModuleLabel_,TrackInstanceName_,TrackDataHandle);
      if( ! foundTrackcollection ) {
        throw cet::exception("GeaneFailedPlots") << "No Trackcollection in this event (\"" << TrackModuleLabel_ << "\":\"" << TrackInstanceName_ << "\")\n";
        return;
    }

      art::Handle<gm2truth::GhostDetectorArtRecordCollection> DummyDataHandle;
      bool foundDummycollection = event.getByLabel(DummyModuleLabel_,DummyInstanceName_,DummyDataHandle);
      if( ! foundDummycollection ) {
        mf::LogWarning(name_) << "No Dummy collection in this event (\"" << DummyModuleLabel_ << "\":\"" << DummyInstanceName_ << "\")\n";
        // return;
    }

    //
    // Record stats
    //

    //Move to next event if we've not got any tracks in this event
    if( (*TrackDataHandle).size() == 0 ) return;


    //
    // Make plots
    //
    FillHistograms(event, TrackDataHandle, DummyDataHandle, topDir);

  }//analyze

  
  void GeaneFailedPlots::beginJob() {

    //
    // Set up plotting helpers
    //
    
    //Create a ROOT file and manager 
    art::ServiceHandle<art::TFileService> tfs;
    auto& outputRootFile = tfs->file();
    rootManager_.reset( new RootManager(name_, &outputRootFile) ); 

    //Create directories
    auto topDir = rootManager_->GetDir(&outputRootFile,dirName_,true); //true -> create if doesn't exist

    //
    // Book histograms
    //

    BookHistograms(topDir);


  }//beginJob



  void GeaneFailedPlots::beginRun(art::Run const & r) {

    //Get coord systems
    cs_ = artg4::dataFromRunOrService<gm2geom::CoordSystemsStoreData, gm2geom::CoordSystemsStore>
          ( r, dataModuleDefs::coordSysModuleLabel(),dataModuleDefs::coordSysInstanceLabel() );
    if( cs_.size() == 0 ) {
      mf::LogWarning(name_) << "This run does not contain any data associated with the coordinate system\n";
    }

  }//beginRun


  void GeaneFailedPlots::endJob() {
    // can clear empty histograms here if I want
  }//endJob




  void GeaneFailedPlots::BookHistograms(TDirectory* dir){

   // std::string fmodes[14] = {"worked", "neg or nan Chi2", "Chi2 diverging", "dummy plane hit > 1", "digit hit > 1", "start pos > target", "geant4e source error", "step length 0", ">7 iterations", "circle fit nan", "", "non primary mchit", "neg chi2 on plane", "x residual large"};
   std::string fmodes[14] = {"worked", "neg or nan Chi2", "Chi2 diverging", "dummy plane hit > 1", "digit hit > 1", "start pos > target", "geant4e source error", "step length 0", "step length 1000", "n steps > 100", "fullSeqFit none converged", "non primary mchit", "tracking in air", "x residual non-zero"};

    //Create directories
   rootManager_->Add( dir, new TH1F( "Worked", "; Worked; Events", 1, 0, 1) );

   rootManager_->Add( dir, new TH1F( "FailureMode", "; Failure Mode; Events", 14, 0, 14) );
   
   for (int bini = 1; bini < 15; ++bini)
   {
        rootManager_->Get<TH1F*>( dir, "FailureMode" )->GetXaxis()->SetBinLabel(bini, fmodes[bini-1].c_str());
   }

  }


  void GeaneFailedPlots::FillHistograms(const art::Event& event, art::Handle<gm2strawtracker::TrackArtRecordCollection> tracksHandle, art::Handle<gm2truth::GhostDetectorArtRecordCollection> dummyHandle, TDirectory* trackDir){

    mf::LogTrace(name_) << "Begin filling histograms. \n";

    // Histogram directories

    // Loop over each track
    gm2strawtracker::TrackArtRecordCollection const& tracks = *tracksHandle; //Resolve handle to get collection

    for( auto & track : tracks ) {

      if (track.failureMode == 0 ) rootManager_->Get<TH1F*>( trackDir, "Worked" )->Fill(track.failureMode);
      else rootManager_->Get<TH1F*>( trackDir, "FailureMode" )->Fill(track.failureMode);
      // else rootManager_->Get<TH1F*>( trackDir, "FailureMode" )->Fill(fmodes[track.failureMode].c_str(), 1);

    } // end auto loop

    mf::LogTrace(name_) << "End filling histograms. \n";

} // end FillHistograms

} // End of namespace gm2strawtracker

// These are some necessary boilerplate for the ROOT persistency system
using gm2strawtracker::GeaneFailedPlots;
DEFINE_ART_MODULE(GeaneFailedPlots)
