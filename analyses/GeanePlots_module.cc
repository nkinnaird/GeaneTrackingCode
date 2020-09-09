//Geane fitting results plots

// Include needed ART headers
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"
#include "gm2tracker/utils/GeaneTrackUtils.hh"

//art records
#include "gm2dataproducts/strawtracker/StrawDigitArtRecord.hh"
#include "gm2dataproducts/strawtracker/StrawDCADigitArtRecord.hh"
#include "gm2dataproducts/strawtracker/StrawTimeIslandArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackCandidateArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"
#include "gm2dataproducts/mc/ghostdetectors/GhostDetectorArtRecord.hh"

//Utils
#include "gm2geom/common/Gm2Constants_service.hh"
#include "gm2util/common/dataModuleDefs.hh"
#include "gm2util/common/RootManager.hh"
#include "gm2util/coordSystems/CoordSystemUtils.hh"
#include "gm2tracker/quality/TrackQuality_service.hh"

//C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <math.h> 

#include <Eigen/Dense>
#include "gm2tracker/utils/GeaneEigenStorageUtils.hh"
#include "gm2tracker/utils/TrackDetailArtRecordUtils.hh"

namespace gm2strawtracker {

  //
  // Class declaration
  //
  class GeanePlots : public art::EDAnalyzer {

  public:

    explicit GeanePlots(fhicl::ParameterSet const& pset);

    //Override desired art::EDAnalyzer functions
    void analyze(const art::Event& event ) override;
    void beginJob() override;
    void beginRun(art::Run const & r) override;
    void endRun(art::Run const & r) override;

  private:

    std::string name_;

    //Producer labels
    std::string TrackModuleLabel_;
    std::string TrackInstanceName_;

    std::string dcaDigitModuleLabel_;
    std::string dcaDigitInstanceLabel_;

    std::string DummyModuleLabel_;
    std::string DummyInstanceName_;

    std::string trajectoryModuleLabel_;
    std::string trajectoryInstanceName_;

    //ROOT plotting members
    std::unique_ptr<RootManager> rootManager_;
    std::string dirName_;

    //Helper tools
    gm2geom::CoordSystemsStoreData cs_;
    gm2geom::StrawTrackerGeometry sgeom_;
    gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;
    art::ServiceHandle<gm2strawtracker::TrackQuality> trackQuality_;

    //Variables that we cut on  
    bool applyTrackQuality_;
    double pValueCut_;
    int numPlanesHitCut_;
    double energyLossCut_;
    vector<double> timeWindow_;
    vector<double> momWindow_;
    vector<int> stations_;

    // Keep track of number of tracks and number of events to work out tracks per event for this run
    int tracksInRun_;
    int eventsInRun_;

    // For setting max/min values of plots (larger for wire fit)
    string fitMode_;

    // Flag for whether we want to make truth plots (we don't know whether it's data or MC at the beginRun stage).
    bool makeTruthPlots_;

    // Flag for grabbing the detail track art record
    bool useTrackDetailArtRecord_;

    gm2util::CoordSystemUtils csUtils_;
    std::map< std::string, gm2geom::CoordSystemsStoreData > detCoordMap_;

    // Data and MC plots
    void BookHistograms(TDirectory* dir);
    void FillHistograms(const art::Event& event);
    void FillPlaneHistos(int numPlanesHit, int stationNum, const gm2strawtracker::TrackDetailArtRecord& track, double& t0);

    // MC-only plots
    void BookTruthHistograms(TDirectory* dir);
    void FillTruthHistograms(const art::Event& event);
    void FillPlaneTruthHistos(int numPlanesHit, int stationNum, const gm2strawtracker::TrackDetailArtRecord& track, art::Ptr<gm2truth::GhostDetectorArtRecord>& dummyHit);
    void FillGlobalPullHistos(int stationNum, const gm2strawtracker::TrackDetailArtRecord& track, art::Ptr<gm2truth::GhostDetectorArtRecord>& dummyHit);

  }; //End of class GeanePlots


  //
  // Class implementation
  //

  GeanePlots::GeanePlots(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset)
    , name_( "GeanePlots" )
    , TrackModuleLabel_( pset.get<std::string>("TrackModuleLabel", dataModuleDefs::trackModuleLabel()) )
    , TrackInstanceName_( pset.get<std::string>("TrackInstanceName", dataModuleDefs::recoInstanceLabel()) )
    , dcaDigitModuleLabel_( pset.get<std::string>("dcaDigitModuleLabel",dataModuleDefs::dcaDigitModuleLabel()) )
    , dcaDigitInstanceLabel_( pset.get<std::string>("dcaDigitInstanceLabel",dataModuleDefs::digitInstanceLabel()) )
    , DummyModuleLabel_( pset.get<std::string>("DummyModuleLabel",dataModuleDefs::strawBuilderModuleLabel()) )
    , DummyInstanceName_( pset.get<std::string>("DummyInstanceName","trackerdummyplane") )
    , trajectoryModuleLabel_( pset.get<std::string>("trajectoryModuleLabel", dataModuleDefs::trajBuilderModuleLabel()) )
    , trajectoryInstanceName_( pset.get<std::string>("trajectoryInstanceName", dataModuleDefs::trajBuilderInstanceLabel()) )
    , rootManager_()
    , dirName_( pset.get<std::string>("dirName","TrackSummary") )
    , cs_()
    , sgeom_()
    , geaneTrackUtils_()
    , applyTrackQuality_( pset.get<bool>("applyTrackQuality", true) )
    , pValueCut_( pset.get<double>("pValueCut", 0.0) )
    , numPlanesHitCut_( pset.get<int>("numPlanesHitCut", 0))
    , energyLossCut_( pset.get<double>("energyLossCut", 1000000.)) // default energy loss cut large
    , timeWindow_(pset.get<vector<double> >("timeWindow",{}))
    , momWindow_(pset.get<vector<double> >("momWindow",{}))
    , stations_(pset.get<vector<int> >("stations",{}))
    , tracksInRun_(0)
    , eventsInRun_(0)
    , fitMode_(pset.get<string>("fitMode","badFitMode"))
    , makeTruthPlots_(pset.get<bool>("makeTruthPlots",true))
    , useTrackDetailArtRecord_(pset.get<bool>("useTrackDetailArtRecord",false))
    , csUtils_()
    , detCoordMap_()
  {}
  

  void GeanePlots::beginJob() {

    //Create a ROOT file and manager 
    art::ServiceHandle<art::TFileService> tfs;
    auto& outputRootFile = tfs->file();
    rootManager_.reset( new RootManager(name_, &outputRootFile) ); 

    //Create directory structure (do it here so that they're in reasonable order)
    auto topDir = rootManager_->GetDir(&outputRootFile,dirName_,true); //true -> create if doesn't exist
    rootManager_->GetDir(topDir,"RunInfo",true);
    rootManager_->GetDir(topDir,"FitResults",true);
    if(makeTruthPlots_){
      rootManager_->GetDir(topDir,"TruthParameters",true);
      rootManager_->GetDir(topDir,"TruthFitComparison",true);
    }
    rootManager_->GetDir(topDir,"PerPlane",true);
    rootManager_->GetDir(topDir,"nPlanesHit",true);
    rootManager_->GetDir(topDir,"nPlanesHit/chiSquaredsPerPlanesHit",true);
    rootManager_->GetDir(topDir,"nPlanesHit/pValuePerPlanesHit",true);
    rootManager_->GetDir(topDir,"Iterations",true);
    rootManager_->GetDir(topDir,"TrackLength",true);
    rootManager_->GetDir(topDir,"energyLoss",true);
    rootManager_->GetDir(topDir,"Digits",true);
    // Book histograms
    BookHistograms(topDir);
    if(makeTruthPlots_) BookTruthHistograms(topDir);


  }//beginJob


  void GeanePlots::beginRun(art::Run const & r) {

    //Get coord systems
    cs_ = artg4::dataFromRunOrService<gm2geom::CoordSystemsStoreData, gm2geom::CoordSystemsStore>
      ( r, dataModuleDefs::coordSysModuleLabel(),dataModuleDefs::coordSysInstanceLabel() );
    if( cs_.size() == 0 ) {
      mf::LogWarning(name_) << "This run does not contain any data associated with the coordinate system\n";
    }

    // Add extra coordinate system maps to speed up transforms
    std::vector<std::string> detNames;
    detNames.push_back("TrackerStation");
    detNames.push_back("TrackerModule");
    for(auto s : sgeom_.whichScallopLocations) {
      for(unsigned int m = 0; m < sgeom_.getNumModulesPerStation(); ++m) {
        detNames.push_back(Form("Module%d:%d", s, m));
      }
    }
    csUtils_.setDetectorNames(detNames);
    csUtils_.detNameToCoordMap(cs_,detCoordMap_);
    sgeom_.setDetNameToCoords(detCoordMap_);

    // Reset tracks and events
    eventsInRun_ = 0;
    tracksInRun_ = 0;

  }//beginRun


  void GeanePlots::analyze(const art::Event& event) {
    
    // Set flag for whether to fill truth plots
    if(event.isRealData()) makeTruthPlots_ = false;

    //Fill plots
    FillHistograms(event);

    // Increment event counter
    eventsInRun_++;

  }//analyze


  void GeanePlots::endRun(art::Run const & r) {

    auto runInfoDirectory = rootManager_->GetDir(rootManager_->GetDir(dirName_),"RunInfo");
    auto tgTracks = rootManager_->Get<TGraph*>(runInfoDirectory, "TracksPerEvent" );
    if(eventsInRun_ > 0) tgTracks->SetPoint(tgTracks->GetN(), r.run(), double(tracksInRun_)/eventsInRun_);
    
  }

  void GeanePlots::BookHistograms(TDirectory* dir){

    // Set scale for some plots
    double boundUV = (fitMode_ == "wireFit") ? 2.535 : 1.5;

    // Run Info (run number and number of tracks per run)
    auto runInfoDirectory = rootManager_->GetDir(dir,"RunInfo");
    rootManager_->Add( runInfoDirectory, new TH1F( "Run", ";Run Number; Tracks", 3000, 6999.5, 9999.5) );
    TGraph* tgTracks = new TGraph();
    tgTracks->SetName( "TracksPerEvent" );
    tgTracks->SetTitle(";Run Number;No. Tracks per Fill");
    tgTracks->SetMarkerStyle(20);
    tgTracks->SetLineWidth(0);
    rootManager_->Add( runInfoDirectory, tgTracks);

    // Fit results (at plane 0)
    auto fitResultsDir = rootManager_->GetDir(dir,"FitResults");
    rootManager_->Add( fitResultsDir, new TH1F( "P", "; P [MeV]; Tracks", 360, 0, 3600) );
    rootManager_->Add( fitResultsDir, new TH1F( "Pu", "; P_{U} [MeV]; Tracks", 500, -1000, 1000) );
    rootManager_->Add( fitResultsDir, new TH1F( "Pv", "; P_{V} [MeV]; Tracks", 500, -1000, 1000) );
    rootManager_->Add( fitResultsDir, new TH1F( "U", "; U [mm]; Tracks", 500, -1000, 1000) );
    rootManager_->Add( fitResultsDir, new TH1F( "V", "; V [mm]; Tracks", 500, -1000, 1000) );
    rootManager_->Add( fitResultsDir, new TH1F( "Px", "; P_{x} [MeV]; Tracks", 500, -1000, 1000) );
    rootManager_->Add( fitResultsDir, new TH1F( "Py", "; P_{y} [MeV]; Tracks", 500, -250, 250) );
    rootManager_->Add( fitResultsDir, new TH1F( "Pz", "; P_{z} [MeV]; Tracks", 360, 0, 3600) );
    rootManager_->Add( fitResultsDir, new TH1F( "Pr", "; P_{r} [MeV]; Tracks", 500, -1000, 1000) );
    rootManager_->Add( fitResultsDir, new TH1F( "X", "; X [mm]; Tracks", 400, -400, 400) );
    rootManager_->Add( fitResultsDir, new TH1F( "Y", "; Y [mm]; Tracks", 400, -100, 100) );
    rootManager_->Add( fitResultsDir, new TH1F( "Z", "; Z [mm]; Tracks", 500, -200, 800) );
    rootManager_->Add( fitResultsDir, new TH1F( "R", "; R [mm]; Tracks", 600, 6600, 7200) );
    rootManager_->Add( fitResultsDir, new TH2F( "EntrancePoint", "; Module X [mm]; Module Y", 500, -75, 175, 400, -100, 100) );
    rootManager_->Add( fitResultsDir, new TH1F( "Chi2", ";#chi^{2}; Events", 200, 0, 100) );
    rootManager_->Add( fitResultsDir, new TH1F( "pValues", "; pValue; Events", 200, 0, 1) );
    rootManager_->Add( fitResultsDir, new TH1F( "Times", ";Track t_{0} (us); Events", 6000, 0, 6000*0.148936) );
    rootManager_->Add( fitResultsDir, new TH1F( "Times_gt_1800MeV", ";Track t_{0} (us); Events", 6000, 0, 6000*0.148936) );
    rootManager_->Add( fitResultsDir, new TH1F( "Times_fastRotation", ";Track t_{0} (us); Events", 36000, 0, 6000*0.148936) );
    rootManager_->Add( fitResultsDir, new TH2F( "Time_vs_Pr", ";Track t_{0} (us) Events; P_{r} [MeV]", 6000, 0, 6000*0.148936, 300, -800, 400) );
    rootManager_->Add( fitResultsDir, new TH2F( "Time_vs_P" , ";Track t_{0} (us) Events; P [MeV]",     6000, 0, 6000*0.148936, 360, 0, 3600) );
    
    // Number of iterations
    auto iterationsDir = rootManager_->GetDir(dir,"Iterations");
    rootManager_->Add( iterationsDir, new TH1F( "numIterations", "; numIterations; Events", 11, 0, 11) );
    rootManager_->Add( iterationsDir, new TH2F( "numIterations vs Chi2", "; numIterations; #chi^{2}", 8, 0, 8, 100, 0, 100) );

    // Number of planes hit
    auto planesHitDir = rootManager_->GetDir(dir, "nPlanesHit");
    rootManager_->Add( planesHitDir, new TH1F( "nPlanesHit", ";Planes Hit; Events", 33, 0, 33) );
    rootManager_->Add( planesHitDir, new TH1F( "nUplanesHit", ";U Planes Hit; Events", 17, 0, 17) );
    rootManager_->Add( planesHitDir, new TH1F( "nVplanesHit", ";V Planes Hit; Events", 17, 0, 17) );
    rootManager_->Add( planesHitDir, new TH2F( "nU vs nV", ";U Planes Hit;V Planes Hit", 17, 0, 17, 17, 0, 17) );
    rootManager_->Add( planesHitDir, new TH2F( "nPlanesHit vs P", ";P;Planes Hit", 110, 0, 3300, 33, 0, 33) );
    rootManager_->Add( planesHitDir, new TH2F( "nUplanesHit vs P", ";P;U Planes Hit", 110, 0, 3300, 17, 0, 17) );
    rootManager_->Add( planesHitDir, new TH2F( "nVplanesHit vs P", ";P;V Planes Hit", 110, 0, 3300, 17, 0, 17) );

    // Energy loss
    auto energyLossDir = rootManager_->GetDir(dir, "energyLoss");
    rootManager_->Add( energyLossDir, new TH1F( "predEnergyLoss", "; MeV lost; Events", 10000, 0, 1000));
    rootManager_->Add( energyLossDir, new TH1F( "predEnergyLoss2MeVBound", "; MeV lost; Events", 2000, 0, 2));
    rootManager_->Add( energyLossDir, new TH1F( "predEnergyLossRev", "; MeV lost; Events", 10000, -1000, 0));
    rootManager_->Add( energyLossDir, new TH1F( "predEnergyLoss2MeVBoundRev", "; MeV lost; Events", 2000, -2, 0));

    /////////////////////////////////////////////////////////////////////////////////////    
    // Per-plane plots (loop over number of planes)
    for (int iPlane = 1; iPlane < geaneTrackUtils_.maxNumPlanes; iPlane++) {

      auto planeNDirectory = rootManager_->GetDir(dir,Form("PerPlane/Plane%d",iPlane),true);

      // Fit/predicted parameters
      auto fitParametersDirectory = rootManager_->GetDir(planeNDirectory, "Fit Parameters", true);
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("P Fit Plane %d", iPlane), "; P; Events", 110, 0, 3300) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("Pu Fit Plane %d", iPlane), "; Pu; Events", 100, -800, 800) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("Pv Fit Plane %d", iPlane), "; Pv; Events", 100, -800, 800) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("Px Fit Plane %d", iPlane), "; Px; Events", 100, -800, 800) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("Py Fit Plane %d", iPlane), "; Py; Events", 100, -200, 200) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("Pz Fit Plane %d", iPlane), "; Pz; Events", 110, 0, 3300) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("Pr Fit Plane %d", iPlane), "; Pr; Events", 100, -800, 800) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("X Fit Plane %d", iPlane), "; X; Events", 100, -400, 400) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("Y Fit Plane %d", iPlane), "; Y; Events", 100, -100, 100) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("Z Fit Plane %d", iPlane), "; Z; Events", 100, -200, 800) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("R Fit Plane %d", iPlane), "; R; Events", 600, 6600, 7200) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("1P Fit Plane %d", iPlane), "; 1/P; Events", 110, 0, 0.0033) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("PuPz Fit Plane %d", iPlane), "; Pu/Px; Events", 100, -.4, .4) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("PvPz Fit Plane %d", iPlane), "; Pv/Px; Events", 100, -.4, .4) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("U Fit Plane %d", iPlane), "; U; Events", 100, -1000, 1000) );
      rootManager_->Add( fitParametersDirectory, new TH1F( Form("V Fit Plane %d", iPlane), "; V; Events", 100, -1000, 1000) );

      // Measured plots
      auto measuredResidualDirectory = rootManager_->GetDir(planeNDirectory, "Measure Residuals", true);
      rootManager_->Add( measuredResidualDirectory, new TH1F( Form("UVresidualsMeasPred Plane %d", iPlane), "; UV residual p-m (mm); Events", 100, -boundUV, boundUV) );
      rootManager_->Add( measuredResidualDirectory, new TH2F( Form("UVresidualsMeasPredvsDCA Plane %d", iPlane), "; DCA (um); UV residual p-m (mm)", 450, -1000, 3500, 300, -boundUV, boundUV) );
      
      // Pull plots
      auto measurePullDirectory = rootManager_->GetDir(planeNDirectory, "Measure Pulls", true);
      rootManager_->Add( measurePullDirectory, new TH1F( Form("UV Measure Pull Plane %d", iPlane), "; #Delta(UV)/#sigma; Events", 100, -10, 10) );
      rootManager_->Add( measurePullDirectory, new TH2F( Form("UV Measure Pull vs DCA Plane %d", iPlane), ";DCA (um);#Delta(UV)/#sigma", 450, -1000, 3500, 100, -10, 10) );

      auto mPullNumDenomDir = rootManager_->GetDir(measurePullDirectory, "Num Denom", true);
      rootManager_->Add( mPullNumDenomDir, new TH1F( Form("Meas Pull Numerator Plane %d", iPlane), "; #Delta(UV); Events", 100, -boundUV, boundUV) );
      rootManager_->Add( mPullNumDenomDir, new TH1F( Form("Meas Pull Denominator Plane %d", iPlane), "; #sigma; Events", 100, 0, boundUV) );
      rootManager_->Add( mPullNumDenomDir, new TH1F( Form("Meas Pull sigmaMeas Plane %d", iPlane), "; #sigma_{Meas}; Events", 100, 0, boundUV) );
      rootManager_->Add( mPullNumDenomDir, new TH1F( Form("Meas Pull sigmaFit Plane %d", iPlane), "; #sigma_{Fit}; Events", 100, 0, boundUV) );
      
      // Plots as a function of distance along track
      auto distMeasurePullDir = rootManager_->GetDir(dir,Form("TrackLength/Dist%d/Measure Pulls",iPlane),true);
      rootManager_->Add( distMeasurePullDir, new TH1F( Form("UV Measure Pull Dist %d", iPlane), "; #Delta(UV)/#sigma; Events", 100, -10, 10) );
      
      auto distMPullNumDenomDir = rootManager_->GetDir(distMeasurePullDir, "Num Denom", true);
      rootManager_->Add( distMPullNumDenomDir, new TH1F( Form("Meas Pull Numerator Dist %d", iPlane), "; #Delta(UV); Events", 100, -boundUV, boundUV) );
      rootManager_->Add( distMPullNumDenomDir, new TH1F( Form("Meas Pull Denominator Dist %d", iPlane), "; #sigma; Events", 100, 0, boundUV) );
      rootManager_->Add( distMPullNumDenomDir, new TH1F( Form("Meas Pull sigmaMeas Dist %d", iPlane), "; #sigma_{Meas}; Events", 100, 0, boundUV) );
      rootManager_->Add( distMPullNumDenomDir, new TH1F( Form("Meas Pull sigmaFit Dist %d", iPlane), "; #sigma_{Fit}; Events", 100, 0, boundUV) );

      // Chi-squared and p-value based on number of planes hit
      auto chi2PlanesHitDirectory = rootManager_->GetDir(planesHitDir,"chiSquaredsPerPlanesHit");
      rootManager_->Add( chi2PlanesHitDirectory, new TH1F( Form("ChiSquareds Planes Hit %d", iPlane), ";#chi^{2}; Events", 100, 0, 100) );
      auto pValPlanesHitDirectory = rootManager_->GetDir(planesHitDir,"pValuePerPlanesHit");
      rootManager_->Add( pValPlanesHitDirectory, new TH1F( Form("pValue Planes Hit %d", iPlane), "; p Value; Events", 100, 0, 1) );

    } // end loop over 33 planes

    // Overall residual and pull plots
    auto planePlotsDirectory = rootManager_->GetDir(dir,"PerPlane");
    rootManager_->Add( planePlotsDirectory, new TH1F( "UVresidualsMeasPred", "; UV residual p-m (mm); Events", 100, -boundUV, boundUV) );
    rootManager_->Add( planePlotsDirectory, new TH2F( "UVresidualsMeasPredvsDCA", "; DCA (um); UV residual p-m (mm)", 450, -1000, 3500, 300, -boundUV, boundUV) );
    rootManager_->Add( planePlotsDirectory, new TH1F( "UV Measure Pull", "; #Delta(UV)/#sigma; Events", 100, -10, 10) );
    rootManager_->Add( planePlotsDirectory, new TH2F( "UV Measure Pull vs DCA", ";DCA (um);#Delta(UV)/#sigma", 450, -1000, 3500, 100, -10, 10) );

    // Overall residual and pull plots - vs distance along track
    auto distPlotsDirectory = rootManager_->GetDir(dir,"TrackLength");
    rootManager_->Add( distPlotsDirectory, new TH2F( "UV Meas Pull vs Dist", "; dist(mm); UV meas pull", 200, 0, 1000, 100, -10, 10) );
    rootManager_->Add( distPlotsDirectory, new TH2F( "Num Meas Pull vs Dist", "; dist(mm); #Delta(UV)", 200, 0, 1000, 100, -boundUV, boundUV) );
    rootManager_->Add( distPlotsDirectory, new TH2F( "Denom Meas Pull vs Dist", "; dist(mm); #sigma", 200, 0, 1000, 100, 0, boundUV) );
    rootManager_->Add( distPlotsDirectory, new TH2F( "sigmaMeas Meas Pull vs Dist", "; dist(mm); #sigma_{Meas}", 200, 0, 1000, 100, 0, boundUV) );
    rootManager_->Add( distPlotsDirectory, new TH2F( "sigmaFit Meas Pull vs Dist", "; dist(mm); #sigma_{Fit}", 200, 0, 1000, 100, 0, boundUV) );

    /////////////////////////////////////////////////////////////////////////////////////

    // Individual straw digit hit info
    auto digitsDir = rootManager_->GetDir(dir,"Digits");
    rootManager_->Add(digitsDir, new TH1F("DriftTime", ";Reco drift time [ns]", 400, -100, 150.) );
    rootManager_->Add(digitsDir, new TH1F("DCA", ";Reco DCA [um]", 400, -1.e3, 3.5e3) );
    rootManager_->Add(digitsDir, new TH1F("DCAError", ";Reco DCA Error [um]", 400, 0, 400) );
    rootManager_->Add(digitsDir, new TH1F("Width", ";Straw digit width [ns]", 96, -10., 110.) );

  }


  void GeanePlots::BookTruthHistograms(TDirectory* dir){

    auto planePlotsDirectory = rootManager_->GetDir(dir,"PerPlane");

    auto energyLossDir = rootManager_->GetDir(dir, "energyLoss");

    rootManager_->Add( energyLossDir, new TH1F( "energyLoss", "; MeV lost; Events", 10000, 0, 1000));
    rootManager_->Add( energyLossDir, new TH1F( "energyLoss2MeVBound", "; MeV lost; Events", 2000, 0, 2));
    rootManager_->Add( energyLossDir, new TH1F( "energyLossRev", "; MeV lost; Events", 10000, -1000, 0));
    rootManager_->Add( energyLossDir, new TH1F( "energyLoss2MeVBoundRev", "; MeV lost; Events", 2000, -2, 0));
    rootManager_->Add( energyLossDir, new TH1F( "energyLossFraction", "; true Eloss / pred Eloss; Events", 1000, 0, 1));

    double boundUV = (fitMode_ == "wireFit") ? 2.535 : 1.5;
    
    // Overall truth parameters
    auto truthDirectory = rootManager_->GetDir(dir,"TruthParameters");
    rootManager_->Add( truthDirectory, new TH1F( "P", "; P; Events", 360, 0, 3600) );
    rootManager_->Add( truthDirectory, new TH1F( "Pu", "; Pu; Events", 500, -800, 800) );
    rootManager_->Add( truthDirectory, new TH1F( "Pv", "; Pv; Events", 500, -800, 800) );
    rootManager_->Add( truthDirectory, new TH1F( "U", "; U; Events", 500, -1000, 1000) );
    rootManager_->Add( truthDirectory, new TH1F( "V", "; V; Events", 500, -1000, 1000) );
    
    rootManager_->Add( truthDirectory, new TH1F( "Px", "; Px; Events", 500, -1000, 1000) );
    rootManager_->Add( truthDirectory, new TH1F( "Py", "; Py; Events", 500, -250, 250) );
    rootManager_->Add( truthDirectory, new TH1F( "Pz", "; Pz; Events", 360, 0, 3600) );
    rootManager_->Add( truthDirectory, new TH1F( "Pr", "; Pr; Events", 500, -1000, 1000) );
    
    rootManager_->Add( truthDirectory, new TH1F( "X", "; X; Events", 400, -400, 400) );
    rootManager_->Add( truthDirectory, new TH1F( "Y", "; Y; Events", 400, -100, 100) );
    rootManager_->Add( truthDirectory, new TH1F( "Z", "; Z; Events", 500, -200, 800) );
    rootManager_->Add( truthDirectory, new TH1F( "R", "; R; Events", 600, 6600, 7200) );
    rootManager_->Add( truthDirectory, new TH2F( "EntrancePoint", "; Module X [mm]; Module Y", 500, -75, 175, 400, -100, 100) );

    // Per plane truth parameters
    for (int iPlane = 1; iPlane < geaneTrackUtils_.maxNumPlanes; iPlane++) {

      auto planeNDirectory = rootManager_->GetDir(planePlotsDirectory,Form("Plane%d",iPlane));
      auto measuredResidualDirectory = rootManager_->GetDir(planeNDirectory, "Measure Residuals");

      auto truthResidualDirectory = rootManager_->GetDir(planeNDirectory, "Truth Residuals", true);
      auto absoluteDir = rootManager_->GetDir(truthResidualDirectory, "Absolute", true);
      auto relativeDir = rootManager_->GetDir(truthResidualDirectory, "Relative", true);

      auto truthParametersDirectory = rootManager_->GetDir(planeNDirectory, "Truth Parameters", true);

      /////////////////////////////////////////////////////////////////////////////////////
      // True to measured residuals
      rootManager_->Add( measuredResidualDirectory, new TH1F( Form("UVresidualsMeasTruth Plane %d", iPlane), "; UV residual m-t (mm); Events", 100, -boundUV, boundUV) );
      rootManager_->Add( measuredResidualDirectory, new TH2F( Form("UVresidualsMeasTruthvsDCA Plane %d", iPlane), "; DCA (um); UV residual m-t (mm)", 450, -1000, 3500, 300, -boundUV, boundUV) );

      /////////////////////////////////////////////////////////////////////////////////////
      // Truth parameters
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("P Truth Plane %d", iPlane), "; P; Events", 110, 0, 3300) );

      rootManager_->Add( truthParametersDirectory, new TH1F( Form("Pu Truth Plane %d", iPlane), "; Pu; Events", 100, -800, 800) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("Pv Truth Plane %d", iPlane), "; Pv; Events", 100, -800, 800) );

      rootManager_->Add( truthParametersDirectory, new TH1F( Form("Px Truth Plane %d", iPlane), "; Px; Events", 100, -800, 800) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("Py Truth Plane %d", iPlane), "; Py; Events", 100, -200, 200) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("Pz Truth Plane %d", iPlane), "; Pz; Events", 110, 0, 3300) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("Pr Truth Plane %d", iPlane), "; Pr; Events", 100, -800, 800) );

      rootManager_->Add( truthParametersDirectory, new TH1F( Form("X Truth Plane %d", iPlane), "; X; Events", 100, -400, 400) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("Y Truth Plane %d", iPlane), "; Y; Events", 100, -100, 100) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("Z Truth Plane %d", iPlane), "; Z; Events", 100, -200, 800) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("R Truth Plane %d", iPlane), "; R; Events", 600, 6600, 7200) );

      rootManager_->Add( truthParametersDirectory, new TH1F( Form("1P Truth Plane %d", iPlane), "; 1/P; Events", 110, 0, 0.0033) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("PuPz Truth Plane %d", iPlane), "; Pu/Px; Events", 100, -.4, .4) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("PvPz Truth Plane %d", iPlane), "; Pv/Px; Events", 100, -.4, .4) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("U Truth Plane %d", iPlane), "; U; Events", 100, -1000, 1000) );
      rootManager_->Add( truthParametersDirectory, new TH1F( Form("V Truth Plane %d", iPlane), "; V; Events", 100, -1000, 1000) );

      /////////////////////////////////////////////////////////////////////////////////////
      // residuals

      // absolute
      rootManager_->Add( absoluteDir, new TH1F( Form("P Truth Residual Abs Plane %d", iPlane), "; #Delta(P); Events", 100, -300, 300) );
      rootManager_->Add( absoluteDir, new TH2F( Form("P Truth Residual Abs vs P Plane %d", iPlane), ";  P; #Delta(P)", 110, 0, 3300, 100, -150, 150));

      rootManager_->Add( absoluteDir, new TH1F( Form("Pu Truth Residual Abs Plane %d", iPlane), "; #Delta(Pu); Events", 100, -100, 100) );
      rootManager_->Add( absoluteDir, new TH1F( Form("Pv Truth Residual Abs Plane %d", iPlane), "; #Delta(Pv); Events", 100, -100, 100) );

      rootManager_->Add( absoluteDir, new TH1F( Form("Px Truth Residual Abs Plane %d", iPlane), "; #Delta(Px); Events", 100, -150, 150) );
      rootManager_->Add( absoluteDir, new TH1F( Form("Pz Truth Residual Abs Plane %d", iPlane), "; #Delta(Pz); Events", 100, -100, 100) );

      rootManager_->Add( absoluteDir, new TH1F( Form("Py Truth Residual Abs Plane %d", iPlane), "; #Delta(Py); Events", 100, -100, 100) );
      rootManager_->Add( absoluteDir, new TH2F( Form("Py Truth Residual Abs vs Py Plane %d", iPlane), "; Py; #Delta(Py)", 100, -100, 100, 100, -100, 100) );

      rootManager_->Add( absoluteDir, new TH1F( Form("X Truth Residual Abs Plane %d", iPlane), "; #Delta(X); Events", 1001, -5, 5) );
      rootManager_->Add( absoluteDir, new TH1F( Form("Y Truth Residual Abs Plane %d", iPlane), "; #Delta(Y); Events", 100, -5, 5) );
      rootManager_->Add( absoluteDir, new TH1F( Form("Z Truth Residual Abs Plane %d", iPlane), "; #Delta(Z); Events", 100, -.001001, 0.001001) );

      rootManager_->Add( absoluteDir, new TH1F( Form("1P Truth Residual Abs Plane %d", iPlane), "; #Delta(1/P); Events", 100, -0.0001, 0.0001) );
      rootManager_->Add( absoluteDir, new TH2F( Form("1P Truth Residual Abs vs P Plane %d", iPlane), "; P; #Delta(1/P)", 110, 0, 3300, 100, -0.0001, 0.0001) );

      rootManager_->Add( absoluteDir, new TH1F( Form("PuPz Truth Residual Abs Plane %d", iPlane), "; #Delta(Pu/Px); Events", 100, -.01, .01) );
      rootManager_->Add( absoluteDir, new TH1F( Form("PvPz Truth Residual Abs Plane %d", iPlane), "; #Delta(Pv/Px); Events", 100, -.01, .01) );
      rootManager_->Add( absoluteDir, new TH1F( Form("U Truth Residual Abs Plane %d", iPlane), "; #Delta(U); Events", 100, -boundUV, boundUV) );
      rootManager_->Add( absoluteDir, new TH1F( Form("V Truth Residual Abs Plane %d", iPlane), "; #Delta(V); Events", 100, -boundUV, boundUV) );

      /////////////////////////////////////////////////////////////////////////////////////
      // relative
      rootManager_->Add( relativeDir, new TH1F( Form("P Truth Residual Rel Plane %d", iPlane), "; #Delta(P)/P; Events", 100, -.2, .2) );
      rootManager_->Add( relativeDir, new TH2F( Form("P Truth Residual Rel vs P Plane %d", iPlane), "; P; #Delta(P)/P;",110, 0, 3300,  100, -.2, .2));

      rootManager_->Add( relativeDir, new TH1F( Form("Px Truth Residual Rel Plane %d", iPlane), "; #Delta(Px)/Px; Events", 100, -.2, .2) );
      rootManager_->Add( relativeDir, new TH1F( Form("Pz Truth Residual Rel Plane %d", iPlane), "; #Delta(Pz)/Pz; Events", 100, -.2, .2) );

      rootManager_->Add( relativeDir, new TH1F( Form("Py Truth Residual Rel Plane %d", iPlane), "; #Delta(Py)/Py; Events", 100, -.2, .2) );
      rootManager_->Add( relativeDir, new TH2F( Form("Py Truth Residual Rel vs Py Plane %d", iPlane), "; Py; #Delta(Py)/Py", 100, -100, 100,  100, -.2, .2) );

      /////////////////////////////////////////////////////////////////////////////////////
      // pull plots on planes
      auto truthPullDirectory = rootManager_->GetDir(planeNDirectory, "Truth Pulls", true);

      rootManager_->Add( truthPullDirectory, new TH1F( Form("1P Truth Pull Plane %d", iPlane), "; #Delta(1/P)/#sigma; Events", 100, -10, 10) );
      rootManager_->Add( truthPullDirectory, new TH1F( Form("PuPz Truth Pull Plane %d", iPlane), "; #Delta(Pu/Px)/#sigma; Events", 100, -10, 10) );
      rootManager_->Add( truthPullDirectory, new TH1F( Form("PvPz Truth Pull Plane %d", iPlane), "; #Delta(Pv/Px)/#sigma; Events", 100, -10, 10) );
      rootManager_->Add( truthPullDirectory, new TH1F( Form("U Truth Pull Plane %d", iPlane), "; #Delta(U)/#sigma; Events", 100, -10, 10) );
      rootManager_->Add( truthPullDirectory, new TH1F( Form("V Truth Pull Plane %d", iPlane), "; #Delta(V)/#sigma; Events", 100, -10, 10) );

      /////////////////////////////////////////////////////////////////////////////////////
      // plots as a function of distance
      
      auto distPlotsDirectory = rootManager_->GetDir(dir,"TrackLength");
      
      stringstream sss;
      sss << "Dist" << iPlane;
      
      auto distNDirectory = rootManager_->GetDir(distPlotsDirectory,sss.str(),true);
      
      auto distTruthPullDir = rootManager_->GetDir(distNDirectory, "Truth Pulls", true);
      
      rootManager_->Add( distTruthPullDir, new TH1F( Form("1P Truth Pull Dist %d", iPlane), "; #Delta(1/P)/#sigma; Events", 100, -10, 10) );
      rootManager_->Add( distTruthPullDir, new TH1F( Form("PuPz Truth Pull Dist %d", iPlane), "; #Delta(Pu/Px)/#sigma; Events", 100, -10, 10) );
      rootManager_->Add( distTruthPullDir, new TH1F( Form("PvPz Truth Pull Dist %d", iPlane), "; #Delta(Pv/Px)/#sigma; Events", 100, -10, 10) );
      rootManager_->Add( distTruthPullDir, new TH1F( Form("U Truth Pull Dist %d", iPlane), "; #Delta(U)/#sigma; Events", 100, -10, 10) );
      rootManager_->Add( distTruthPullDir, new TH1F( Form("V Truth Pull Dist %d", iPlane), "; #Delta(V)/#sigma; Events", 100, -10, 10) );
      
    } // loop over 32 planes

    /////////////////////////////////////////////////////////////////////////////////////
    // only fill 0 planes hit pulls, which means fill for all tracks
    auto truthPullDirectory = rootManager_->GetDir(dir, "TruthFitComparison");

    double bound1overP = 2.5e-4;
    double boundPuvoverPx = 0.022;

    //pulls
    rootManager_->Add( truthPullDirectory, new TH1F( "1P Truth Pull", "; #Delta(1/P)/#sigma; Events", 100, -10, 10) );
    rootManager_->Add( truthPullDirectory, new TH1F( "PuPz Truth Pull", "; #Delta(Pu/Px)/#sigma; Events", 100, -10, 10) );
    rootManager_->Add( truthPullDirectory, new TH1F( "PvPz Truth Pull", "; #Delta(Pv/Px)/#sigma; Events", 100, -10, 10) );
    rootManager_->Add( truthPullDirectory, new TH1F( "U Truth Pull", "; #Delta(U)/#sigma; Events", 100, -10, 10) );
    rootManager_->Add( truthPullDirectory, new TH1F( "V Truth Pull", "; #Delta(V)/#sigma; Events", 100, -10, 10) );

    auto pullPullDir = rootManager_->GetDir(truthPullDirectory, "Pull Pull", true);

    rootManager_->Add( pullPullDir, new TH2F( "1P vs PuPz Pull", "; pull(1/P); pull(Pu/Px)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "1P vs PvPz Pull", "; pull(1/P); pull(Pv/Px)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "1P vs U Pull", "; pull(1/P); pull(U)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "1P vs V Pull", "; pull(1/P); pull(V)", 100, -10, 10, 100, -10, 10) );

    rootManager_->Add( pullPullDir, new TH2F( "PuPz vs 1P Pull", "; pull(Pu/Px); pull(1/P)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "PuPz vs PvPz Pull", "; pull(Pu/Px); pull(Pv/Px)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "PuPz vs U Pull", "; pull(Pu/Px); pull(U)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "PuPz vs V Pull", "; pull(Pu/Px); pull(V)", 100, -10, 10, 100, -10, 10) );

    rootManager_->Add( pullPullDir, new TH2F( "PvPz vs 1P Pull", "; pull(Pv/Px); pull(1/P)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "PvPz vs PuPz Pull", "; pull(Pv/Px); pull(Pu/Px)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "PvPz vs U Pull", "; pull(Pv/Px); pull(U)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "PvPz vs V Pull", "; pull(Pv/Px); pull(V)", 100, -10, 10, 100, -10, 10) );

    rootManager_->Add( pullPullDir, new TH2F( "U vs 1P Pull", "; pull(U); pull(1/P)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "U vs PuPz Pull", "; pull(U); pull(Pu/Px)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "U vs PvPz Pull", "; pull(U); pull(Pv/Px)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "U vs V Pull", "; pull(U); pull(V)", 100, -10, 10, 100, -10, 10) );

    rootManager_->Add( pullPullDir, new TH2F( "V vs 1P Pull", "; pull(V); pull(1/P)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "V vs PuPz Pull", "; pull(V); pull(Pu/Px)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "V vs PvPz Pull", "; pull(V); pull(Pv/Px)", 100, -10, 10, 100, -10, 10) );
    rootManager_->Add( pullPullDir, new TH2F( "V vs U Pull", "; pull(V); pull(U)", 100, -10, 10, 100, -10, 10) );

    auto NumDenomDir = rootManager_->GetDir(truthPullDirectory, "Num Denom", true);
      
    rootManager_->Add( NumDenomDir, new TH1F( "1P Truth Pull Numerator", "; #Delta(1/P); Events", 100, -bound1overP, bound1overP) );
    rootManager_->Add( NumDenomDir, new TH1F( "PuPz Truth Pull Numerator", "; #Delta(Pu/Px); Events", 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumDenomDir, new TH1F( "PvPz Truth Pull Numerator", "; #Delta(Pv/Px); Events", 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumDenomDir, new TH1F( "U Truth Pull Numerator", "; #Delta(U); Events", 100, -boundUV, boundUV) );
    rootManager_->Add( NumDenomDir, new TH1F( "V Truth Pull Numerator", "; #Delta(V); Events", 100, -boundUV, boundUV) );

    rootManager_->Add( NumDenomDir, new TH1F( "1P Truth Pull Denominator", "; #sigma; Events", 100, 0, bound1overP) );
    rootManager_->Add( NumDenomDir, new TH1F( "PuPz Truth Pull Denominator", "; #sigma; Events", 100, 0, boundPuvoverPx) );
    rootManager_->Add( NumDenomDir, new TH1F( "PvPz Truth Pull Denominator", "; #sigma; Events", 100, 0, boundPuvoverPx) );
    rootManager_->Add( NumDenomDir, new TH1F( "U Truth Pull Denominator", "; #sigma; Events", 100, 0, boundUV) );
    rootManager_->Add( NumDenomDir, new TH1F( "V Truth Pull Denominator", "; #sigma; Events", 100, 0, boundUV) );

    rootManager_->Add( NumDenomDir, new TH2F( "Denom vs Num 1P", "; #Delta(1/P); #sigma 1/P", 100, -bound1overP, bound1overP, 100, 0, bound1overP) );
    rootManager_->Add( NumDenomDir, new TH2F( "Denom vs Num PuPz", "; #Delta(Pu/Px); #sigma Pu/Px", 100, -boundPuvoverPx, boundPuvoverPx, 100, 0, boundPuvoverPx) );
    rootManager_->Add( NumDenomDir, new TH2F( "Denom vs Num PvPz", "; #Delta(Pv/Px); #sigma Pv/Px", 100, -boundPuvoverPx, boundPuvoverPx, 100, 0, boundPuvoverPx) );
    rootManager_->Add( NumDenomDir, new TH2F( "Denom vs Num U", "; #Delta(U); #sigma U", 100, -boundUV, boundUV, 100, 0, boundUV) );
    rootManager_->Add( NumDenomDir, new TH2F( "Denom vs Num V", "; #Delta(V); #sigma V", 100, -boundUV, boundUV, 100, 0, boundUV) );

    rootManager_->Add( NumDenomDir, new TH2F( "Denom 1P vs P", "; P; #sigma 1/P", 110, 0, 3300, 100, 0, bound1overP) );
    rootManager_->Add( NumDenomDir, new TH2F( "Denom PuPz vs P", "; P; #sigma Pu/Px", 110, 0, 3300, 100, 0, boundPuvoverPx) );
    rootManager_->Add( NumDenomDir, new TH2F( "Denom PvPz vs P", "; P; #sigma Pv/Px", 110, 0, 3300, 100, 0, boundPuvoverPx) );
    rootManager_->Add( NumDenomDir, new TH2F( "Denom U vs P", "; P; #sigma U", 110, 0, 3300, 100, 0, boundUV) );
    rootManager_->Add( NumDenomDir, new TH2F( "Denom V vs P", "; P; #sigma V", 110, 0, 3300, 100, 0, boundUV) );

    auto fitErrorVsGuessDir = rootManager_->GetDir(truthPullDirectory, "Fit Error Vs Guess", true);

    rootManager_->Add( fitErrorVsGuessDir, new TH2F( "Sigma 1P vs Sdiff 1P", "; diff 1/P; #sigma 1/P", 100, -0.001, 0.001,  100, 0, bound1overP) );
    rootManager_->Add( fitErrorVsGuessDir, new TH2F( "Sigma PUPz vs Sdiff PUPz", "; diff PU/PX; #sigma PU/PX", 100, -.40, .40,  100, 0, boundPuvoverPx) );
    rootManager_->Add( fitErrorVsGuessDir, new TH2F( "Sigma PVPz vs Sdiff PVPz", "; diff PV/PX; #sigma PV/PX", 100, -.40, .40,  100, 0, boundPuvoverPx) );
    rootManager_->Add( fitErrorVsGuessDir, new TH2F( "Sigma U vs Sdiff U", "; diff U; #sigma U", 100, -8, 8,  100, 0, boundUV) );
    rootManager_->Add( fitErrorVsGuessDir, new TH2F( "Sigma V vs Sdiff V", "; diff V; #sigma V", 100, -8, 8,  100, 0, boundUV) );

    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma 1P over Guess Diff", "; #sigma/guessDiff; Events", 100, -40, 40) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma PUPz over Guess Diff", "; #sigma/guessDiff; Events", 100, -40, 40) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma PVPz over Guess Diff", "; #sigma/guessDiff; Events", 100, -40, 40) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma U over Guess Diff", "; #sigma/guessDiff; Events", 100, -40, 40) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma V over Guess Diff", "; #sigma/guessDiff; Events", 100, -40, 40) );

    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Guess Diff 1P over Sigma", "; guessDiff/#sigma; Events", 200, -100, 100) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Guess Diff PUPz over Sigma", "; guessDiff/#sigma; Events", 200, -100, 100) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Guess Diff PVPz over Sigma", "; guessDiff/#sigma; Events", 200, -100, 100) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Guess Diff U over Sigma", "; guessDiff/#sigma; Events", 200, -100, 100) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Guess Diff V over Sigma", "; guessDiff/#sigma; Events", 200, -100, 100) );

    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma 1P over Guess Diff Max", "; #sigma/guessDiffMax; Events", 100, -10, 10) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma PUPz over Guess Diff Max", "; #sigma/guessDiffMax; Events", 100, -10, 10) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma PVPz over Guess Diff Max", "; #sigma/guessDiffMax; Events", 100, -10, 10) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma U over Guess Diff Max", "; #sigma/guessDiffMax; Events", 100, -10, 10) );
    rootManager_->Add( fitErrorVsGuessDir, new TH1F( "Sigma V over Guess Diff Max", "; #sigma/guessDiffMax; Events", 100, -10, 10) );

    auto DenomVsDenomDir = rootManager_->GetDir(truthPullDirectory, "Fit Error Vs Fit Error", true);

    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom 1P vs 1P", "; #sigma 1/P; #sigma 1/P", 100, 0, bound1overP, 100, 0, bound1overP) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom 1P vs PuPz", "; #sigma 1/P; #sigma Pu/Px", 100, 0, bound1overP, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom 1P vs PvPz", "; #sigma 1/P; #sigma Pv/Px", 100, 0, bound1overP, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom 1P vs U", "; #sigma 1/P; #sigma U", 100, 0, bound1overP, 100, 0, boundUV) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom 1P vs V", "; #sigma 1/P; #sigma V", 100, 0, bound1overP, 100, 0, boundUV) );

    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PuPz vs 1P", "; #sigma Pu/Px; #sigma 1/P", 100, 0, boundPuvoverPx, 100, 0, bound1overP) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PuPz vs PuPz", "; #sigma Pu/Px; #sigma Pu/Px", 100, 0, boundPuvoverPx, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PuPz vs PvPz", "; #sigma Pu/Px; #sigma Pv/Px", 100, 0, boundPuvoverPx, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PuPz vs U", "; #sigma Pu/Px; #sigma U", 100, 0, boundPuvoverPx, 100, 0, boundUV) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PuPz vs V", "; #sigma Pu/Px; #sigma V", 100, 0, boundPuvoverPx, 100, 0, boundUV) );

    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PvPz vs 1P", "; #sigma Pv/Px; #sigma 1/P", 100, 0, boundPuvoverPx, 100, 0, bound1overP) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PvPz vs PuPz", "; #sigma Pv/Px; #sigma Pu/Px", 100, 0, boundPuvoverPx, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PvPz vs PvPz", "; #sigma Pv/Px; #sigma Pv/Px", 100, 0, boundPuvoverPx, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PvPz vs U", "; #sigma Pv/Px; #sigma U", 100, 0, boundPuvoverPx, 100, 0, boundUV) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom PvPz vs V", "; #sigma Pv/Px; #sigma V", 100, 0, boundPuvoverPx, 100, 0, boundUV) );

    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom U vs 1P", "; #sigma U; #sigma 1/P", 100, 0, boundUV, 100, 0, bound1overP) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom U vs PuPz", "; #sigma U; #sigma Pu/Px", 100, 0, boundUV, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom U vs PvPz", "; #sigma U; #sigma Pv/Px", 100, 0, boundUV, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom U vs U", "; #sigma U; #sigma U", 100, 0, boundUV, 100, 0, boundUV) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom U vs V", "; #sigma U; #sigma V", 100, 0, boundUV, 100, 0, boundUV) );

    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom V vs 1P", "; #sigma V; #sigma 1/P", 100, 0, boundUV, 100, 0, bound1overP) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom V vs PuPz", "; #sigma V; #sigma Pu/Px", 100, 0, boundUV, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom V vs PvPz", "; #sigma V; #sigma Pv/Px", 100, 0, boundUV, 100, 0, boundPuvoverPx) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom V vs U", "; #sigma V; #sigma U", 100, 0, boundUV, 100, 0, boundUV) );
    rootManager_->Add( DenomVsDenomDir, new TH2F( "Denom V vs V", "; #sigma V; #sigma V", 100, 0, boundUV, 100, 0, boundUV) );

    auto NumVsNumDir = rootManager_->GetDir(truthPullDirectory, "Delta Vs Delta", true);

    rootManager_->Add( NumVsNumDir, new TH2F( "Num 1P vs 1P", "; delta 1/P; delta 1/P", 100, -bound1overP, bound1overP, 100, -bound1overP, bound1overP) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num 1P vs PuPz", "; delta 1/P; delta Pu/Px", 100, -bound1overP, bound1overP, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num 1P vs PvPz", "; delta 1/P; delta Pv/Px", 100, -bound1overP, bound1overP, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num 1P vs U", "; delta 1/P; delta U", 100, -bound1overP, bound1overP, 100, -boundUV, boundUV) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num 1P vs V", "; delta 1/P; delta V", 100, -bound1overP, bound1overP, 100, -boundUV, boundUV) );

    rootManager_->Add( NumVsNumDir, new TH2F( "Num PuPz vs 1P", "; delta Pu/Px; delta 1/P", 100, -boundPuvoverPx, boundPuvoverPx, 100, -bound1overP, bound1overP) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num PuPz vs PuPz", "; delta Pu/Px; delta Pu/Px", 100, -boundPuvoverPx, boundPuvoverPx, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num PuPz vs PvPz", "; delta Pu/Px; delta Pv/Px", 100, -boundPuvoverPx, boundPuvoverPx, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num PuPz vs U", "; delta Pu/Px; delta U", 100, -boundPuvoverPx, boundPuvoverPx, 100, -boundUV, boundUV) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num PuPz vs V", "; delta Pu/Px; delta V", 100, -boundPuvoverPx, boundPuvoverPx, 100, -boundUV, boundUV) );

    rootManager_->Add( NumVsNumDir, new TH2F( "Num PvPz vs 1P", "; delta Pv/Px; delta 1/P", 100, -boundPuvoverPx, boundPuvoverPx, 100, -bound1overP, bound1overP) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num PvPz vs PuPz", "; delta Pv/Px; delta Pu/Px", 100, -boundPuvoverPx, boundPuvoverPx, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num PvPz vs PvPz", "; delta Pv/Px; delta Pv/Px", 100, -boundPuvoverPx, boundPuvoverPx, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num PvPz vs U", "; delta Pv/Px; delta U", 100, -boundPuvoverPx, boundPuvoverPx, 100, -boundUV, boundUV) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num PvPz vs V", "; delta Pv/Px; delta V", 100, -boundPuvoverPx, boundPuvoverPx, 100, -boundUV, boundUV) );

    rootManager_->Add( NumVsNumDir, new TH2F( "Num U vs 1P", "; delta U; delta 1/P", 100, -boundUV, boundUV, 100, -bound1overP, bound1overP) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num U vs PuPz", "; delta U; delta Pu/Px", 100, -boundUV, boundUV, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num U vs PvPz", "; delta U; delta Pv/Px", 100, -boundUV, boundUV, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num U vs U", "; delta U; delta U", 100, -boundUV, boundUV, 100, -boundUV, boundUV) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num U vs V", "; delta U; delta V", 100, -boundUV, boundUV, 100, -boundUV, boundUV) );

    rootManager_->Add( NumVsNumDir, new TH2F( "Num V vs 1P", "; delta V; delta 1/P", 100, -boundUV, boundUV, 100, -bound1overP, bound1overP) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num V vs PuPz", "; delta V; delta Pu/Px", 100, -boundUV, boundUV, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num V vs PvPz", "; delta V; delta Pv/Px", 100, -boundUV, boundUV, 100, -boundPuvoverPx, boundPuvoverPx) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num V vs U", "; delta V; delta U", 100, -boundUV, boundUV, 100, -boundUV, boundUV) );
    rootManager_->Add( NumVsNumDir, new TH2F( "Num V vs V", "; delta V; delta V", 100, -boundUV, boundUV, 100, -boundUV, boundUV) );
      
    // Meas - True plots
    rootManager_->Add( planePlotsDirectory, new TH1F( "UVresidualsMeasTruth", "; UV residual m-t (mm); Events", 100, -boundUV, boundUV) );
    rootManager_->Add( planePlotsDirectory, new TH2F( "UVresidualsMeasTruthvsDCA", "; DCA (um); UV residual m-t (mm)", 450, -1000, 3500, 300, -boundUV, boundUV) );
    
    /////////////////////////////////////////////////////////////////////////////////////
    auto trajPlotsDir = rootManager_->GetDir(dir,"TrajPlots",true);

    rootManager_->Add( trajPlotsDir, new TH1F( "cosAngle", "; cosAngle; Events", 1010, 0, 1.01));
    rootManager_->Add( trajPlotsDir, new TH2F( "cosAngleMomDiff", "; cosAngle; Mom", 1010, 0, 1.01, 101, -.1, 10));
    rootManager_->Add( trajPlotsDir, new TH2F( "cosAngleEdep", "; cosAngle; Edep", 1010, 0, 1.01, 101, -.1, 10));
    
    rootManager_->Add( trajPlotsDir, new TH1F( "maxCosAngle", "; maxCosAngle; Events", 1010, 0, 1.01));
    
    /////////////////////////////////////////////////////////////////////////////////////
  }


  void GeanePlots::FillHistograms(const art::Event& event){

    mf::LogInfo info(name_);
    mf::LogTrace(name_) << "Begin filling histograms. \n";

    // Get histogram directories
    auto dir = rootManager_->GetDir(dirName_);
    auto runInfoDirectory = rootManager_->GetDir(dir,"RunInfo");
    auto fitResultsDirectory = rootManager_->GetDir(dir,"FitResults");
    auto iterationsDir = rootManager_->GetDir(dir,"Iterations");
    auto planesHitDir = rootManager_->GetDir(dir, "nPlanesHit");
    auto chi2PlanesHitDirectory = rootManager_->GetDir(planesHitDir,"chiSquaredsPerPlanesHit");
    auto pValPlanesHitDirectory = rootManager_->GetDir(planesHitDir,"pValuePerPlanesHit");
    auto energyLossDir = rootManager_->GetDir(dir, "energyLoss");

    //Get Track Fitting Results
    art::Handle<gm2strawtracker::TrackArtRecordCollection> tracksHandle;
    bool foundTrackCollection = event.getByLabel(TrackModuleLabel_,TrackInstanceName_,tracksHandle);
    if( ! foundTrackCollection ) {
      throw cet::exception("GeanePlots") << "No Trackcollection in this event (\"" << TrackModuleLabel_ << "\":\"" << TrackInstanceName_ << "\")\n";
      return;
    }
   
    art::Handle<gm2strawtracker::TrackDetailArtRecordCollection> tracksDetailHandle;
    if( useTrackDetailArtRecord_ ) {
      foundTrackCollection = event.getByLabel(TrackModuleLabel_,TrackInstanceName_,tracksDetailHandle);
      if( ! foundTrackCollection ) {
        throw cet::exception("GeanePlots") << "No Trackcollection in this event (\"" << TrackModuleLabel_ << "\":\"" << TrackInstanceName_ << "\")\n";
        return;
      }
    } 
 
    //Don't fill anything if we've not got any tracks in this event
    if( !useTrackDetailArtRecord_ && (*tracksHandle).size() == 0 )       return;
    if( useTrackDetailArtRecord_  && (*tracksDetailHandle).size() == 0 ) return; 

    //Get true Dummy Plane Records
    art::Handle<gm2truth::GhostDetectorArtRecordCollection> dummyHandle;
    bool foundDummyCollection = event.getByLabel(DummyModuleLabel_,DummyInstanceName_,dummyHandle);
    if( ! foundDummyCollection ) {
      mf::LogWarning(name_) << "No Dummy collection in this event (\"" << DummyModuleLabel_ << "\":\"" << DummyInstanceName_ << "\")\n";
    }
    
    //Get true trajectory Records
    art::Handle<gm2truth::TrajectoryArtRecordCollection> trajHandle;
    bool foundTrajCollection = event.getByLabel(trajectoryModuleLabel_, trajectoryInstanceName_, trajHandle);
    if( ! foundTrajCollection ) {
       mf::LogWarning(name_) << "No trajectories in this event (\"" << trajectoryModuleLabel_ << "\":\"" << trajectoryInstanceName_ << "\")\n";
    }


    const gm2strawtracker::TrackArtRecordCollection & trackArtRecordData = tracksHandle->size() == 0  ? gm2strawtracker::TrackArtRecordCollection() : *tracksHandle;
    auto tracksData = useTrackDetailArtRecord_  ? *tracksDetailHandle : TrackDetailArtRecordUtils::FillTrackDetailArtRecord(trackArtRecordData);
  
    if( tracksData.size() == 0 ) {
      throw cet::exception("GeanePlots") << "Something went wrong when using the TrackDetailArtRecordUtis!\n";
      return;
    }

    // Loop over each track and decide whether to cut
    for( auto & trackData : tracksData ) {

      // Resolve handle to get collection - now with eigen object conversion
      auto track            = GeaneEigenStorageUtils::ReadEigenFromStorage(trackData);
      auto geaneHitsOnTrack = track.geaneHits;

      // Apply cuts
      if (applyTrackQuality_ && !trackQuality_->goodTrack(track)) continue;
      if (track.failureMode != 0 ) continue;
      if (track.pValue < pValueCut_) continue;
      if (track.trackNumPlanesHit < numPlanesHitCut_) continue;
      if (timeWindow_.size() == 2 && (track.candidate->t0 < timeWindow_[0] or track.candidate->t0 > timeWindow_[1])) continue;

      auto startingTrackPx = track.geaneHits.startingGeaneParameters.at(3);
      auto startingTrackPy = track.geaneHits.startingGeaneParameters.at(4);
      auto startingTrackPz = track.geaneHits.startingGeaneParameters.at(5);
 
      double startingTrackMomentum = sqrt(startingTrackPx*startingTrackPx + startingTrackPy*startingTrackPy + startingTrackPz*startingTrackPz );
      if (momWindow_.size() == 2 && (startingTrackMomentum < momWindow_[0] or startingTrackMomentum > momWindow_[1])) continue;

      int thisStation = track.candidate->strawDigits.at(0)->wireID.getStation();
      if (stations_.size() > 0 && (std::find(stations_.begin(),stations_.end(), thisStation) == stations_.end())) continue;

      double energyDiff = 0;

      if(foundDummyCollection){
	energyDiff = track.dummyPlaneHits.at(0)->momentum.mag() - track.dummyPlaneHits.back()->momentum.mag();
	if (energyDiff > energyLossCut_){
	  info << "LOSSCATCH: event " << event.event() << " with greater energy loss than momentum tolerance: " << energyDiff << "\n";
	  continue;
	}
      }

      rootManager_->Get<TH1F*>( runInfoDirectory, "Run" )->Fill(event.run());
      rootManager_->Get<TH1F*>( fitResultsDirectory, "Chi2" )->Fill(track.chi2);
      rootManager_->Get<TH1F*>( fitResultsDirectory, "pValues" )->Fill(track.pValue);
      rootManager_->Get<TH1F*>( fitResultsDirectory, "Times" )->Fill( (track.candidate->t0) * 1.e-3 ); // ns to us
      if(startingTrackMomentum > 1800) rootManager_->Get<TH1F*>( fitResultsDirectory, "Times_gt_1800MeV" )->Fill( (track.candidate->t0) * 1.e-3 ); // ns to us
      rootManager_->Get<TH1F*>( fitResultsDirectory, "Times_fastRotation" )->Fill( (track.candidate->t0) * 1.e-3 ); // ns to us

      rootManager_->Get<TH1F*>( chi2PlanesHitDirectory, Form("ChiSquareds Planes Hit %d", track.trackNumPlanesHit) )->Fill(track.chi2); 
      rootManager_->Get<TH1F*>( pValPlanesHitDirectory, Form("pValue Planes Hit %d", track.trackNumPlanesHit) )->Fill(track.pValue); 

      /////////////////////////////////////////////////////////////////////////////////////

      double endTrackMomentum = geaneTrackUtils_.getPredMom(track.geaneHits, track.trackLastPlaneHit);
      double predEnergyLoss   = startingTrackMomentum - endTrackMomentum;

      rootManager_->Get<TH1F*>( energyLossDir, "predEnergyLoss" )->Fill(predEnergyLoss);
      rootManager_->Get<TH1F*>( energyLossDir, "predEnergyLoss2MeVBound" )->Fill(predEnergyLoss);
      rootManager_->Get<TH1F*>( energyLossDir, "predEnergyLossRev" )->Fill(-predEnergyLoss);
      rootManager_->Get<TH1F*>( energyLossDir, "predEnergyLoss2MeVBoundRev" )->Fill(-predEnergyLoss);

      if(foundDummyCollection and makeTruthPlots_){
	rootManager_->Get<TH1F*>( energyLossDir, "energyLoss" )->Fill(energyDiff);
	rootManager_->Get<TH1F*>( energyLossDir, "energyLoss2MeVBound" )->Fill(energyDiff);
	rootManager_->Get<TH1F*>( energyLossDir, "energyLossRev" )->Fill(-energyDiff);
	rootManager_->Get<TH1F*>( energyLossDir, "energyLoss2MeVBoundRev" )->Fill(-energyDiff);

	rootManager_->Get<TH1F*>( energyLossDir, "energyLossFraction" )->Fill(energyDiff/predEnergyLoss);
      }


      /////////////////////////////////////////////////////////////////////////////////////
      rootManager_->Get<TH1F*>(iterationsDir, "numIterations")->Fill(track.numIterations);
      rootManager_->Get<TH2F*>(iterationsDir, "numIterations vs Chi2")->Fill(track.numIterations, track.chi2);

      /////////////////////////////////////////////////////////////////////////////////////

      int numUHit = 0;
      int numVHit = 0;
 
      auto trackCandidateT0 = track.candidate->t0;

      for (int i = 0; i < int(track.trackPlanesHitList.size()); ++i) {// loop over hit planes
        
	int planeNum = track.trackPlanesHitList.at(i);
	if (geaneTrackUtils_.isUPlane(planeNum)) numUHit++;
	else numVHit++;

	FillPlaneHistos(planeNum, thisStation, track, trackCandidateT0);

	/////////////////////////////////////////////////////////////////////////////////////                
	if (foundDummyCollection and makeTruthPlots_){
	  auto dummyHit = track.dummyPlaneHits.at(i+1); // i+1 since dummy plane hits is size N+1
	  FillPlaneTruthHistos(planeNum, thisStation, track, dummyHit);
	} // dummy hit if statement
      } // end loop over track planes hit

      rootManager_->Get<TH1F*>( planesHitDir, "nPlanesHit" )->Fill(track.trackNumPlanesHit);

      rootManager_->Get<TH1F*>( planesHitDir, "nUplanesHit" )->Fill(numUHit);
      rootManager_->Get<TH1F*>( planesHitDir, "nVplanesHit" )->Fill(numVHit);
      rootManager_->Get<TH2F*>( planesHitDir, "nU vs nV" )->Fill(numUHit, numVHit);

      rootManager_->Get<TH2F*>(planesHitDir, "nPlanesHit vs P")->Fill(startingTrackMomentum, track.trackNumPlanesHit);
      rootManager_->Get<TH2F*>(planesHitDir, "nUplanesHit vs P")->Fill(startingTrackMomentum, numUHit);
      rootManager_->Get<TH2F*>(planesHitDir, "nVplanesHit vs P")->Fill(startingTrackMomentum, numVHit);

      // Fill overall (plane 0) plots
      FillPlaneHistos(0, thisStation, track, trackCandidateT0);
      if (foundDummyCollection and makeTruthPlots_){
	auto dummyHit = track.dummyPlaneHits.at(0);
	FillPlaneTruthHistos(0, thisStation, track, dummyHit);
	FillGlobalPullHistos(thisStation, track, dummyHit);
      }

      /////////////////////////////////////////////////////////////////////////////////////

      // Straw digit in track information
      auto digitsDir = rootManager_->GetDir(dir,"Digits");

      // Get dcaDigits associated with hitDigits (through art:Assns in event)
      // Use FindMany as there may be multiple dcaDigits associated to digit if it's used in multiple track candidates
      art::FindManyP<StrawDCADigitArtRecord, int> digitAssns(track.candidate->strawDigits, event, art::InputTag(dcaDigitModuleLabel_,dcaDigitInstanceLabel_));

      // Loop over associated digits
      for(unsigned int i_digit = 0; i_digit < track.candidate->strawDigits.size(); i_digit++){

	// Fill vector of dcaDigits from Assn
	std::vector<art::Ptr<StrawDCADigitArtRecord> > dcaDigits;
	std::vector<int const*> trackCandIds;
	digitAssns.get(i_digit,dcaDigits,trackCandIds);

	// Get digit ptr from this dca digit - find the one with the right track candidate ID (in case there's more than one)
	unsigned int dcaDigitIndex = -1;
	for(unsigned int iDCADig = 0; iDCADig < trackCandIds.size(); iDCADig++){
	  if(*(trackCandIds.at(iDCADig)) == track.candidate->id){
	    dcaDigitIndex = iDCADig;
	    break;
	  }
	}
	if(dcaDigitIndex < 0) throw cet::exception(name_) << "Event " << event.id() << ": trackCandidate with id " << track.candidate->id << " was not found in list of ids returned from Assn\n";
	auto& dcaDigit = dcaDigits.at(dcaDigitIndex);

	// Fill information
	rootManager_->Get<TH1F*>(digitsDir, "DriftTime")->Fill(dcaDigit->driftTime);
	rootManager_->Get<TH1F*>(digitsDir, "DCA")->Fill(1e3*dcaDigit->dca);
	rootManager_->Get<TH1F*>(digitsDir, "DCAError")->Fill(1e3*dcaDigit->dcaError);
	rootManager_->Get<TH1F*>(digitsDir, "Width")->Fill(dcaDigit->digit->hitWidth);
      }

      /////////////////////////////////////////////////////////////////////////////////////
      // Trajectory Plots
/*
      if (foundTrajCollection and makeTruthPlots_){

	auto trajPlotsDir = rootManager_->GetDir(dir,"TrajPlots");

	gm2truth::TrajectoryArtRecordCollection const& trajectories = *trajHandle;
	
	for( auto traj : trajectories) {
	  if (traj.trackID != 4 || traj.pdgID != -11) continue;

	  double maxCosAngle = 1.0;

	  for (int i = 0; i < int(traj.points.size()); ++i){
	    auto point = traj.points.at(i);

	    int stationNumber = track.candidate->strawDigits.at(0)->wireID.getStation();

	    stringstream stationStream;
	    stationStream << "TrackerStation[" << stationNumber << "]";
	    string stationStr = stationStream.str();

	    gm2geom::CoordSystem3Vector pointPos(point.position, "world");

	    auto css = detCoordMap_.find("TrackerStation")->second;
	    gm2geom::CoordSystem3Vector pointPosGeane = pointPos.transform(css, stationStr);

	    if (pointPosGeane.x() > track.startingGeaneParameters.at(0) && pointPosGeane.x() < track.planeXPositions.at(track.trackPlanesHitList.back())){
	      // G4cout << "Step momentum: " << point.momentum << G4endl;

	      if (i != 0){
		G4ThreeVector prevMomentum = traj.points.at(i-1).momentum;
		G4ThreeVector currMomentum = point.momentum;

		double cosAngle = (prevMomentum.x()*currMomentum.x() + prevMomentum.y()*currMomentum.y() + prevMomentum.z()*currMomentum.z()) / (prevMomentum.mag() * currMomentum.mag());

		std::cout.precision(17);
		if (cosAngle >= 1.00001) {
		  std::cout << "prevMomentum x: " << prevMomentum.x() << " " << "currMomentum x: " << currMomentum.x() << std::endl;
		  std::cout << "prevMomentum y: " << prevMomentum.y() << " " << "currMomentum y: " << currMomentum.y() << std::endl;
		  std::cout << "prevMomentum z: " << prevMomentum.z() << " " << "currMomentum z: " << currMomentum.z() << std::endl;
		  std::cout << "prevMomentum mag: " << prevMomentum.mag() << " " << "currMomentum mag: " << currMomentum.mag() << std::endl;
		  std::cout << "Cos(angle) = " << cosAngle << std::endl;
		  std::cout << "Event: " << event.event() << " subRun: " << event.subRun() << std::endl;
		  info << "cosAngle > 1.00001 " << " Event: " << event.event() << " subRun: " << event.subRun() << std::endl;
		}


		if (cosAngle < maxCosAngle) maxCosAngle = cosAngle; // larger angle for smaller cosAngle

		rootManager_->Get<TH1F*>( trajPlotsDir, "cosAngle")->Fill(cosAngle);
		rootManager_->Get<TH2F*>( trajPlotsDir, "cosAngleMomDiff")->Fill(cosAngle, prevMomentum.mag()-currMomentum.mag());
		rootManager_->Get<TH2F*>( trajPlotsDir, "cosAngleEdep")->Fill(cosAngle, point.edep);

	      }
	    }
	  }
	  rootManager_->Get<TH1F*>( trajPlotsDir, "maxCosAngle")->Fill(maxCosAngle);

	}
      } // end traj plots
*/


      tracksInRun_++;

    } // end auto loop over tracks

    mf::LogTrace(name_) << "End filling histograms. \n";

  } // end FillHistograms main method

  void GeanePlots::FillPlaneHistos(int planeNum, int stationNum, const gm2strawtracker::TrackDetailArtRecord & track, double & trackCandidateT0){

    mf::LogTrace(name_) << "Begin filling plane histograms. \n";

    auto & geaneHitsOnTrack = track.geaneHits;

    stringstream stationStream;
    stationStream << "TrackerStation[" << stationNum << "]";
    string stationStr = stationStream.str();

    // Upstream module coordsystem
    int firstModule = track.candidate->upstreamDigit->wireID.getModule();

    auto planePlotsDirectory = rootManager_->GetDir(rootManager_->GetDir(dirName_), "PerPlane");

    /////////////////////////////////////////////////////////////////////////////////////
    // predicted parameters

    double predictedMomentum = geaneTrackUtils_.getPredMom(geaneHitsOnTrack, planeNum);
    double predictedXMomentum = geaneTrackUtils_.getPredXMom(geaneHitsOnTrack, planeNum);
    double predictedYMomentum = geaneTrackUtils_.getPredYMom(geaneHitsOnTrack, planeNum);
    double predictedZMomentum = geaneTrackUtils_.getPredZMom(geaneHitsOnTrack, planeNum);
    double predictedUMomentum = geaneTrackUtils_.getPredUMom(geaneHitsOnTrack, planeNum);
    double predictedVMomentum = geaneTrackUtils_.getPredVMom(geaneHitsOnTrack, planeNum);
    double predictedXPosition = geaneTrackUtils_.getPredXPos(geaneHitsOnTrack, planeNum);
    double predictedYPosition = geaneTrackUtils_.getPredYPos(geaneHitsOnTrack, planeNum);
    double predictedZPosition = geaneTrackUtils_.getPredZPos(geaneHitsOnTrack, planeNum);
    double predictedUPosition = geaneTrackUtils_.getPredUPos(geaneHitsOnTrack, planeNum);
    double predictedVPosition = geaneTrackUtils_.getPredVPos(geaneHitsOnTrack, planeNum);

    /////////////////////////////////////////////////////////////////////////////////////
    // fill plots
    if (planeNum > 0) {

      auto measuredResidualDirectory = rootManager_->GetDir(planePlotsDirectory, Form("Plane%d/Measure Residuals",planeNum));

      if(geaneTrackUtils_.isUPlane(planeNum)){
	rootManager_->Get<TH1F*>( measuredResidualDirectory, Form("UVresidualsMeasPred Plane %d", planeNum) )->Fill(predictedUPosition - geaneHitsOnTrack.geaneMeasuredParameters[3].at(planeNum));
	rootManager_->Get<TH2F*>( measuredResidualDirectory, Form("UVresidualsMeasPredvsDCA Plane %d", planeNum) )->Fill(1e3*track.measuredDCAs.at(planeNum), predictedUPosition - geaneHitsOnTrack.geaneMeasuredParameters[3].at(planeNum));
	rootManager_->Get<TH1F*>( planePlotsDirectory, "UVresidualsMeasPred" )->Fill(predictedUPosition - geaneHitsOnTrack.geaneMeasuredParameters[3].at(planeNum));
	rootManager_->Get<TH2F*>( planePlotsDirectory, "UVresidualsMeasPredvsDCA" )->Fill(1e3*track.measuredDCAs.at(planeNum), predictedUPosition - geaneHitsOnTrack.geaneMeasuredParameters[3].at(planeNum));
      } else {
	rootManager_->Get<TH1F*>( measuredResidualDirectory, Form("UVresidualsMeasPred Plane %d", planeNum) )->Fill(predictedVPosition - geaneHitsOnTrack.geaneMeasuredParameters[4].at(planeNum));
	rootManager_->Get<TH2F*>( measuredResidualDirectory, Form("UVresidualsMeasPredvsDCA Plane %d", planeNum) )->Fill(1e3*track.measuredDCAs.at(planeNum), predictedVPosition - geaneHitsOnTrack.geaneMeasuredParameters[4].at(planeNum));
	rootManager_->Get<TH1F*>( planePlotsDirectory, "UVresidualsMeasPred" )->Fill(predictedVPosition - geaneHitsOnTrack.geaneMeasuredParameters[4].at(planeNum));
	rootManager_->Get<TH2F*>( planePlotsDirectory, "UVresidualsMeasPredvsDCA" )->Fill(1e3*track.measuredDCAs.at(planeNum), predictedVPosition - geaneHitsOnTrack.geaneMeasuredParameters[4].at(planeNum));
      }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // parameters

    // Calculate radial fit results
    gm2geom::CoordSystem3Vector predPosStation(predictedXPosition, predictedYPosition, predictedZPosition, stationStr);
    gm2geom::CoordSystem3Vector predMomStation(predictedXMomentum, predictedYMomentum, predictedZMomentum, stationStr);
     
    auto css = detCoordMap_.find("TrackerStation")->second;
    gm2geom::CoordSystem3Vector predPosWorld = predPosStation.transform(css,"world");
    gm2geom::CoordSystem3Vector predMomWorld = predMomStation.transform(css,"world",true); //true -> momentum type transform

    css = detCoordMap_.find("TrackerModule")->second;
    gm2geom::CoordSystem3Vector predPosModule = predPosWorld.transform(css,Form("TrackerModule[%d][%d]",stationNum,firstModule));

    // Get radial position and momentum
    predPosWorld.setY(0); // In world coordinates, XZ is radial plane so we get rid of Y component
    double radialPosition = predPosWorld.mag();
    double radialMomentum = predMomWorld.dot(predPosWorld.unit());

    // fit
    if(planeNum > 0){
      auto fitParametersDirectory = rootManager_->GetDir(planePlotsDirectory, Form("Plane%d/Fit Parameters",planeNum));

      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("P Fit Plane %d", planeNum) )->Fill( predictedMomentum );

      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("Pu Fit Plane %d", planeNum) )->Fill( predictedUMomentum );
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("Pv Fit Plane %d", planeNum) )->Fill( predictedVMomentum );

      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("Px Fit Plane %d", planeNum) )->Fill( predictedXMomentum );
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("Py Fit Plane %d", planeNum) )->Fill( predictedYMomentum );
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("Pz Fit Plane %d", planeNum) )->Fill( predictedZMomentum );
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("Pr Fit Plane %d", planeNum) )->Fill( radialMomentum );
      
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("X Fit Plane %d", planeNum) )->Fill( predictedXPosition );
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("Y Fit Plane %d", planeNum) )->Fill( predictedYPosition );
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("Z Fit Plane %d", planeNum) )->Fill( predictedZPosition );
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("R Fit Plane %d", planeNum) )->Fill( radialPosition );

      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("1P Fit Plane %d", planeNum) )->Fill(1./predictedMomentum);
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("PuPz Fit Plane %d", planeNum) )->Fill(predictedUMomentum/predictedZMomentum ); 
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("PvPz Fit Plane %d", planeNum) )->Fill(predictedVMomentum/predictedZMomentum );                    
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("U Fit Plane %d", planeNum) )->Fill(predictedUPosition ); 
      rootManager_->Get<TH1F*>( fitParametersDirectory, Form("V Fit Plane %d", planeNum) )->Fill(predictedVPosition );  

    } else {
      auto fitResultsDir = rootManager_->GetDir(rootManager_->GetDir(dirName_), "FitResults");
      rootManager_->Get<TH1F*>( fitResultsDir, "P")->Fill( predictedMomentum );
      rootManager_->Get<TH1F*>( fitResultsDir, "Pu")->Fill( predictedUMomentum );
      rootManager_->Get<TH1F*>( fitResultsDir, "Pv")->Fill( predictedVMomentum );
      rootManager_->Get<TH1F*>( fitResultsDir, "U")->Fill(predictedUPosition ); 
      rootManager_->Get<TH1F*>( fitResultsDir, "V")->Fill(predictedVPosition );  
      rootManager_->Get<TH1F*>( fitResultsDir, "Px")->Fill( predictedXMomentum );
      rootManager_->Get<TH1F*>( fitResultsDir, "Py")->Fill( predictedYMomentum );
      rootManager_->Get<TH1F*>( fitResultsDir, "Pz")->Fill( predictedZMomentum );
      rootManager_->Get<TH1F*>( fitResultsDir, "Pr")->Fill( radialMomentum );
      rootManager_->Get<TH1F*>( fitResultsDir, "X")->Fill( predictedXPosition );
      rootManager_->Get<TH1F*>( fitResultsDir, "Y")->Fill( predictedYPosition );
      rootManager_->Get<TH1F*>( fitResultsDir, "Z")->Fill( predictedZPosition );
      rootManager_->Get<TH1F*>( fitResultsDir, "R")->Fill( radialPosition );
      rootManager_->Get<TH2F*>( fitResultsDir, "Time_vs_Pr")->Fill( trackCandidateT0 * 1.e-3, radialMomentum );// ns to us
      rootManager_->Get<TH2F*>( fitResultsDir, "Time_vs_P")->Fill( trackCandidateT0 * 1.e-3, predictedMomentum );// ns to us
      rootManager_->Get<TH2F*>( fitResultsDir, "EntrancePoint")->Fill( predPosModule.getX(), predPosModule.getY() );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // this code copied below, if making changes change both          
    Eigen::MatrixXd true0Covariance = geaneHitsOnTrack.covarianceTotalInverse.inverse();
    Eigen::MatrixXd planeCovarianceProp;


    if (planeNum != 0){
      int distInt = planeNum - 4*int((track.trackFirstPlaneHit-1)/4);

      double distZ = predictedZPosition - geaneHitsOnTrack.startingGeaneParameters.at(2);
      // Eigen::MatrixXd planeMaterialPropError = (geaneTrackUtils_.JacobianToUV * geaneHitsOnTrack.geaneErrorMatrices[planeNum] * geaneTrackUtils_.JacobianToUV.transpose()); // was .inverse() and should be .transpose(), but it still shouldn't be needed at all
      planeCovarianceProp = geaneHitsOnTrack.extendedTransportMatrixBegToEnd[planeNum]*true0Covariance*geaneHitsOnTrack.extendedTransportMatrixBegToEnd[planeNum].transpose();

      auto topDir = rootManager_->GetDir(dirName_);
      auto distPlotsDirectory = rootManager_->GetDir(topDir,"TrackLength");

      stringstream sss;
      sss << "Dist" << distInt;

      auto distNDirectory = rootManager_->GetDir(distPlotsDirectory,sss.str());
      auto distMeasurePullDir = rootManager_->GetDir(distNDirectory, "Measure Pulls");
      auto distMPullNumDenomDir = rootManager_->GetDir(distMeasurePullDir, "Num Denom");

      double predPos = 0.0;
      int j = 50;
      if(geaneTrackUtils_.isUPlane(planeNum)){
	j=3;
	predPos = predictedUPosition;
      } else {
	j=4;
	predPos = predictedVPosition;
      } 

      double measError = track.UVerrors.at(planeNum)*track.UVerrors.at(planeNum); // + 100*planeMaterialPropError(j,j); // this last part shouldn't be needed
      double fitError = 100.*planeCovarianceProp(j,j);
      // there is some reason to believe that the material error should instead be added to the fit error, instead of the measured error, but then the resulting pull distributions are off

      double pullNum = predPos - geaneHitsOnTrack.geaneMeasuredParameters[j].at(planeNum);
      double pullDenom = sqrt(measError - fitError);

      rootManager_->Get<TH1F*>( distMeasurePullDir, Form("UV Measure Pull Dist %d", distInt))->Fill( pullNum/pullDenom );

      rootManager_->Get<TH1F*>( distMPullNumDenomDir,  Form("Meas Pull Numerator Dist %d", distInt))->Fill(pullNum);
      rootManager_->Get<TH1F*>( distMPullNumDenomDir,  Form("Meas Pull Denominator Dist %d", distInt))->Fill(pullDenom);
      rootManager_->Get<TH1F*>( distMPullNumDenomDir,  Form("Meas Pull sigmaMeas Dist %d", distInt))->Fill(sqrt(measError));
      rootManager_->Get<TH1F*>( distMPullNumDenomDir,  Form("Meas Pull sigmaFit Dist %d", distInt))->Fill(sqrt(fitError));

      auto measurePullDirectory = rootManager_->GetDir(planePlotsDirectory, Form("Plane%d/Measure Pulls",planeNum));
      auto mPullNumDenomDir = rootManager_->GetDir(measurePullDirectory, "Num Denom");

      rootManager_->Get<TH1F*>( measurePullDirectory, Form("UV Measure Pull Plane %d", planeNum))->Fill( pullNum/pullDenom );
      rootManager_->Get<TH2F*>( measurePullDirectory, Form("UV Measure Pull vs DCA Plane %d", planeNum))->Fill(1e3*track.measuredDCAs.at(planeNum), pullNum/pullDenom);

      rootManager_->Get<TH1F*>( planePlotsDirectory, "UV Measure Pull" )->Fill( pullNum/pullDenom );
      rootManager_->Get<TH2F*>( planePlotsDirectory, "UV Measure Pull vs DCA" )->Fill(1e3*track.measuredDCAs.at(planeNum), pullNum/pullDenom);

      rootManager_->Get<TH1F*>( mPullNumDenomDir,  Form("Meas Pull Numerator Plane %d", planeNum))->Fill(pullNum);
      rootManager_->Get<TH1F*>( mPullNumDenomDir,  Form("Meas Pull Denominator Plane %d", planeNum))->Fill(pullDenom);
      rootManager_->Get<TH1F*>( mPullNumDenomDir,  Form("Meas Pull sigmaMeas Plane %d", planeNum))->Fill(sqrt(measError));
      rootManager_->Get<TH1F*>( mPullNumDenomDir,  Form("Meas Pull sigmaFit Plane %d", planeNum))->Fill(sqrt(fitError));

      rootManager_->Get<TH2F*>( distPlotsDirectory, "UV Meas Pull vs Dist")->Fill(distZ, pullNum/pullDenom);
      rootManager_->Get<TH2F*>( distPlotsDirectory, "Num Meas Pull vs Dist")->Fill(distZ, pullNum);
      rootManager_->Get<TH2F*>( distPlotsDirectory, "Denom Meas Pull vs Dist")->Fill(distZ, pullDenom);
      rootManager_->Get<TH2F*>( distPlotsDirectory, "sigmaMeas Meas Pull vs Dist")->Fill(distZ, sqrt(measError));
      rootManager_->Get<TH2F*>( distPlotsDirectory, "sigmaFit Meas Pull vs Dist")->Fill(distZ, sqrt(fitError));

    } // planeNum != 0

    mf::LogTrace(name_) << "End of filling plane histograms. \n";

  }

  void GeanePlots::FillPlaneTruthHistos(int planeNum, int stationNum, const gm2strawtracker::TrackDetailArtRecord & track, art::Ptr<gm2truth::GhostDetectorArtRecord> & dummyHit){

    mf::LogTrace(name_) << "Begin filling plane truth histograms. \n";

    auto geaneHitsOnTrack = track.geaneHits;

    stringstream stationStream;
    stationStream << "TrackerStation[" << stationNum << "]";
    string stationStr = stationStream.str();

    // Upstream module coordsystem
    int firstModule = track.candidate->upstreamDigit->wireID.getModule();

    auto planePlotsDirectory = rootManager_->GetDir(rootManager_->GetDir(dirName_), "PerPlane");

    gm2geom::CoordSystem3Vector pHitPos(dummyHit->position.x(), dummyHit->position.y(), dummyHit->position.z() , "world");
    auto css = detCoordMap_.find("TrackerStation")->second;
    gm2geom::CoordSystem3Vector planeIntersection = pHitPos.transform(css, stationStr);

    double planeUposition = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*planeIntersection.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,1)*planeIntersection.y();
    double planeVposition = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,0)*planeIntersection.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,1)*planeIntersection.y();

    double planeXposition = planeIntersection.x();
    double planeYposition = planeIntersection.y();
    double planeZposition = planeIntersection.z();


    gm2geom::CoordSystem3Vector pHitMom(dummyHit->momentum.x(), dummyHit->momentum.y(), dummyHit->momentum.z() , "world");
    gm2geom::CoordSystem3Vector momGlob = pHitMom.transform(css, stationStr, true);

    double planeUmomentum = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*momGlob.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,1)*momGlob.y();
    double planeVmomentum = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,0)*momGlob.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,1)*momGlob.y();
                        
    double planeXmomentum = momGlob.x();
    double planeYmomentum = momGlob.y();
    double planeZmomentum = momGlob.z();

    double planeTotalmomentum = momGlob.mag();

    /////////////////////////////////////////////////////////////////////////////////////
    // predicted parameters
    
    double predictedMomentum = geaneTrackUtils_.getPredMom(geaneHitsOnTrack, planeNum);
    double predictedXMomentum = geaneTrackUtils_.getPredXMom(geaneHitsOnTrack, planeNum);
    double predictedYMomentum = geaneTrackUtils_.getPredYMom(geaneHitsOnTrack, planeNum);
    double predictedZMomentum = geaneTrackUtils_.getPredZMom(geaneHitsOnTrack, planeNum);
    double predictedUMomentum = geaneTrackUtils_.getPredUMom(geaneHitsOnTrack, planeNum);
    double predictedVMomentum = geaneTrackUtils_.getPredVMom(geaneHitsOnTrack, planeNum);
    double predictedXPosition = geaneTrackUtils_.getPredXPos(geaneHitsOnTrack, planeNum);
    double predictedYPosition = geaneTrackUtils_.getPredYPos(geaneHitsOnTrack, planeNum);
    double predictedZPosition = geaneTrackUtils_.getPredZPos(geaneHitsOnTrack, planeNum);
    double predictedUPosition = geaneTrackUtils_.getPredUPos(geaneHitsOnTrack, planeNum);
    double predictedVPosition = geaneTrackUtils_.getPredVPos(geaneHitsOnTrack, planeNum);


    /////////////////////////////////////////////////////////////////////////////////////
    // fill plots

    if (planeNum != 0) {

      auto measuredResidualDirectory = rootManager_->GetDir(planePlotsDirectory, Form("Plane%d/Measure Residuals",planeNum));

      if(geaneTrackUtils_.isUPlane(planeNum)){
	rootManager_->Get<TH1F*>( measuredResidualDirectory, Form("UVresidualsMeasTruth Plane %d", planeNum) )->Fill(geaneHitsOnTrack.geaneMeasuredParameters[3].at(planeNum) - planeUposition);
	rootManager_->Get<TH1F*>( planePlotsDirectory, "UVresidualsMeasTruth" )->Fill(geaneHitsOnTrack.geaneMeasuredParameters[3].at(planeNum) - planeUposition);
	rootManager_->Get<TH2F*>( measuredResidualDirectory, Form("UVresidualsMeasTruthvsDCA Plane %d", planeNum) )->Fill(1e3*track.measuredDCAs.at(planeNum), geaneHitsOnTrack.geaneMeasuredParameters[3].at(planeNum) - planeUposition);
	rootManager_->Get<TH2F*>( planePlotsDirectory, "UVresidualsMeasTruthvsDCA" )->Fill(1e3*track.measuredDCAs.at(planeNum), geaneHitsOnTrack.geaneMeasuredParameters[3].at(planeNum) - planeUposition);
      } else {
	rootManager_->Get<TH1F*>( measuredResidualDirectory, Form("UVresidualsMeasTruth Plane %d", planeNum) )->Fill(geaneHitsOnTrack.geaneMeasuredParameters[4].at(planeNum) - planeVposition);
	rootManager_->Get<TH1F*>( planePlotsDirectory, "UVresidualsMeasTruth" )->Fill(geaneHitsOnTrack.geaneMeasuredParameters[4].at(planeNum) - planeVposition);
	rootManager_->Get<TH2F*>( measuredResidualDirectory, Form("UVresidualsMeasTruthvsDCA Plane %d", planeNum) )->Fill(1e3*track.measuredDCAs.at(planeNum), geaneHitsOnTrack.geaneMeasuredParameters[4].at(planeNum) - planeVposition);
	rootManager_->Get<TH2F*>( planePlotsDirectory, "UVresidualsMeasTruthvsDCA" )->Fill(1e3*track.measuredDCAs.at(planeNum), geaneHitsOnTrack.geaneMeasuredParameters[4].at(planeNum) - planeVposition);
      }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // parameters

    // Calculate radial truth data
    gm2geom::CoordSystem3Vector truePosStation(planeXposition, planeYposition, planeZposition, stationStr);
    gm2geom::CoordSystem3Vector trueMomStation(planeXmomentum, planeYmomentum, planeZmomentum, stationStr);
      
    gm2geom::CoordSystem3Vector truePosWorld = truePosStation.transform(css,"world");
    gm2geom::CoordSystem3Vector trueMomWorld = trueMomStation.transform(css,"world",true); //true -> momentum type transform

    css = detCoordMap_.find("TrackerModule")->second;
    gm2geom::CoordSystem3Vector truePosModule = truePosWorld.transform(css,Form("TrackerModule[%d][%d]",stationNum,firstModule));

    // Get radial position and momentum
    truePosWorld.setY(0); // In world coordinates, XZ is radial plane so we get rid of Y component
    double radialPosition = truePosWorld.mag();
    double radialMomentum = trueMomWorld.dot(truePosWorld.unit());

    // truth
    if(planeNum > 0){
      auto truthParametersDirectory = rootManager_->GetDir(planePlotsDirectory, Form("Plane%d/Truth Parameters",planeNum));

      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("P Truth Plane %d", planeNum) )->Fill( planeTotalmomentum );
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("Pu Truth Plane %d", planeNum) )->Fill( planeUmomentum );
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("Pv Truth Plane %d", planeNum) )->Fill( planeVmomentum );

      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("Px Truth Plane %d", planeNum) )->Fill( planeXmomentum );
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("Py Truth Plane %d", planeNum) )->Fill( planeYmomentum );
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("Pz Truth Plane %d", planeNum) )->Fill( planeZmomentum );
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("Pr Truth Plane %d", planeNum) )->Fill( radialMomentum );
                   
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("X Truth Plane %d", planeNum) )->Fill( planeXposition );
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("Y Truth Plane %d", planeNum) )->Fill( planeYposition );
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("Z Truth Plane %d", planeNum) )->Fill( planeZposition );
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("R Truth Plane %d", planeNum) )->Fill( radialPosition );

      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("1P Truth Plane %d", planeNum) )->Fill(1./planeTotalmomentum);
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("PuPz Truth Plane %d", planeNum) )->Fill(planeUmomentum/planeZmomentum ); 
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("PvPz Truth Plane %d", planeNum) )->Fill(planeVmomentum/planeZmomentum ); 
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("U Truth Plane %d", planeNum) )->Fill(planeUposition ); 
      rootManager_->Get<TH1F*>( truthParametersDirectory, Form("V Truth Plane %d", planeNum) )->Fill(planeVposition );                        


      /////////////////////////////////////////////////////////////////////////////////////
      // residuals fit - truth                       

      auto truthResidualDirectory = rootManager_->GetDir(planePlotsDirectory, Form("Plane%d/Truth Residuals",planeNum));

      // absolute
      auto absoluteDir = rootManager_->GetDir(truthResidualDirectory, "Absolute");
      rootManager_->Get<TH1F*>( absoluteDir, Form("P Truth Residual Abs Plane %d", planeNum) )->Fill( predictedMomentum - planeTotalmomentum );
      rootManager_->Get<TH2F*>( absoluteDir, Form("P Truth Residual Abs vs P Plane %d", planeNum) )->Fill( predictedMomentum , predictedMomentum - planeTotalmomentum);

      rootManager_->Get<TH1F*>( absoluteDir, Form("Pu Truth Residual Abs Plane %d", planeNum) )->Fill( predictedUMomentum - planeUmomentum );
      rootManager_->Get<TH1F*>( absoluteDir, Form("Pv Truth Residual Abs Plane %d", planeNum) )->Fill( predictedVMomentum - planeVmomentum );

      rootManager_->Get<TH1F*>( absoluteDir, Form("Px Truth Residual Abs Plane %d", planeNum) )->Fill( predictedXMomentum - planeXmomentum );
      rootManager_->Get<TH1F*>( absoluteDir, Form("Pz Truth Residual Abs Plane %d", planeNum) )->Fill( predictedZMomentum - planeZmomentum );

      rootManager_->Get<TH1F*>( absoluteDir, Form("Py Truth Residual Abs Plane %d", planeNum) )->Fill( predictedYMomentum - planeYmomentum );
      rootManager_->Get<TH2F*>( absoluteDir, Form("Py Truth Residual Abs vs Py Plane %d", planeNum) )->Fill( predictedYMomentum, predictedYMomentum - planeYmomentum );
                   
      // G4cout << G4endl << "subPlaneNum number: " << subPlaneNum << G4endl;
      // G4cout << G4endl << "X residual: " << (predictedXPosition - planeXposition) << G4endl;

      rootManager_->Get<TH1F*>( absoluteDir, Form("X Truth Residual Abs Plane %d", planeNum) )->Fill( predictedXPosition - planeXposition );
      rootManager_->Get<TH1F*>( absoluteDir, Form("Y Truth Residual Abs Plane %d", planeNum) )->Fill( predictedYPosition - planeYposition );
      rootManager_->Get<TH1F*>( absoluteDir, Form("Z Truth Residual Abs Plane %d", planeNum) )->Fill( predictedZPosition - planeZposition );

      rootManager_->Get<TH1F*>( absoluteDir, Form("1P Truth Residual Abs Plane %d", planeNum) )->Fill(1./predictedMomentum - 1./planeTotalmomentum );
      rootManager_->Get<TH2F*>( absoluteDir, Form("1P Truth Residual Abs vs P Plane %d", planeNum) )->Fill( predictedMomentum , 1./predictedMomentum - 1./planeTotalmomentum );

      rootManager_->Get<TH1F*>( absoluteDir, Form("PuPz Truth Residual Abs Plane %d", planeNum) )->Fill(predictedUMomentum/predictedZMomentum - planeUmomentum/planeZmomentum ); 
      rootManager_->Get<TH1F*>( absoluteDir, Form("PvPz Truth Residual Abs Plane %d", planeNum) )->Fill(predictedVMomentum/predictedZMomentum - planeVmomentum/planeZmomentum ); 
      rootManager_->Get<TH1F*>( absoluteDir, Form("U Truth Residual Abs Plane %d", planeNum) )->Fill(predictedUPosition - planeUposition ); 
      rootManager_->Get<TH1F*>( absoluteDir, Form("V Truth Residual Abs Plane %d", planeNum) )->Fill(predictedVPosition - planeVposition ); 
    
      /////////////////////////////////////////////////////////////////////////////////////
      //relative
      auto relativeDir = rootManager_->GetDir(truthResidualDirectory, "Relative");

      rootManager_->Get<TH1F*>( relativeDir, Form("P Truth Residual Rel Plane %d", planeNum) )->Fill( (predictedMomentum - planeTotalmomentum)/predictedMomentum );
      rootManager_->Get<TH2F*>( relativeDir, Form("P Truth Residual Rel vs P Plane %d", planeNum) )->Fill( predictedMomentum , (predictedMomentum - planeTotalmomentum)/predictedMomentum );

      rootManager_->Get<TH1F*>( relativeDir, Form("Px Truth Residual Rel Plane %d", planeNum) )->Fill( (predictedXMomentum - planeXmomentum)/predictedXMomentum );
      rootManager_->Get<TH1F*>( relativeDir, Form("Pz Truth Residual Rel Plane %d", planeNum) )->Fill( (predictedZMomentum - planeZmomentum)/predictedZMomentum );

      rootManager_->Get<TH1F*>( relativeDir, Form("Py Truth Residual Rel Plane %d", planeNum) )->Fill( (predictedYMomentum - planeYmomentum)/predictedYMomentum );
      rootManager_->Get<TH2F*>( relativeDir, Form("Py Truth Residual Rel vs Py Plane %d", planeNum) )->Fill( predictedYMomentum , (predictedYMomentum - planeYmomentum)/predictedYMomentum );
    } else {

      auto truthDir = rootManager_->GetDir(rootManager_->GetDir(dirName_),"TruthParameters");
      rootManager_->Get<TH1F*>( truthDir, "P" )->Fill( planeTotalmomentum );
      rootManager_->Get<TH1F*>( truthDir, "Pu" )->Fill( planeUmomentum );
      rootManager_->Get<TH1F*>( truthDir, "Pv" )->Fill( planeVmomentum );
      rootManager_->Get<TH1F*>( truthDir, "U" )->Fill(planeUposition ); 
      rootManager_->Get<TH1F*>( truthDir, "V" )->Fill(planeVposition );                        

      rootManager_->Get<TH1F*>( truthDir, "Px" )->Fill( planeXmomentum );
      rootManager_->Get<TH1F*>( truthDir, "Py" )->Fill( planeYmomentum );
      rootManager_->Get<TH1F*>( truthDir, "Pz" )->Fill( planeZmomentum );
      rootManager_->Get<TH1F*>( truthDir, "Pr" )->Fill( radialMomentum );
                   
      rootManager_->Get<TH1F*>( truthDir, "X" )->Fill( planeXposition );
      rootManager_->Get<TH1F*>( truthDir, "Y" )->Fill( planeYposition );
      rootManager_->Get<TH1F*>( truthDir, "Z" )->Fill( planeZposition );
      rootManager_->Get<TH1F*>( truthDir, "R" )->Fill( radialPosition );
      rootManager_->Get<TH2F*>( truthDir, "EntrancePoint")->Fill( truePosModule.getX(), truePosModule.getY() );

    }

    /////////////////////////////////////////////////////////////////////////////////////
    // this code copied above, if making changes change both
    Eigen::MatrixXd true0Covariance = geaneHitsOnTrack.covarianceTotalInverse.inverse();
    Eigen::MatrixXd planeCovarianceProp;

    if (planeNum != 0){
      int distInt = planeNum - 4*int((track.trackFirstPlaneHit-1)/4);

      // Eigen::MatrixXd planeMaterialPropError = (geaneTrackUtils_.JacobianToUV * geaneHitsOnTrack.geaneErrorMatrices[planeNum] * geaneTrackUtils_.JacobianToUV.transpose()); // was .inverse() and should be .transpose(), but it still shouldn't be needed at all
      planeCovarianceProp = geaneHitsOnTrack.extendedTransportMatrixBegToEnd[planeNum]*true0Covariance*geaneHitsOnTrack.extendedTransportMatrixBegToEnd[planeNum].transpose();

      double tPull1PDenom = (1./1000)*sqrt(planeCovarianceProp(0,0)); // - planeMaterialPropError(0,0)); // this last part shouldn't be needed
      double tPullPuPzDenom = sqrt(planeCovarianceProp(1,1)); // - planeMaterialPropError(1,1));
      double tPullPvPzDenom = sqrt(planeCovarianceProp(2,2)); // - planeMaterialPropError(2,2));
      double tPullUDenom = 10.*sqrt(planeCovarianceProp(3,3)); // - planeMaterialPropError(3,3));
      double tPullVDenom = 10.*sqrt(planeCovarianceProp(4,4)); // - planeMaterialPropError(4,4));

      auto truthPullDirectory = rootManager_->GetDir(planePlotsDirectory, Form("Plane%d/Truth Pulls",planeNum));
      rootManager_->Get<TH1F*>( truthPullDirectory, Form("1P Truth Pull Plane %d", planeNum))->Fill((1./predictedMomentum - 1./planeTotalmomentum )/tPull1PDenom);
      rootManager_->Get<TH1F*>( truthPullDirectory, Form("PuPz Truth Pull Plane %d", planeNum))->Fill((predictedUMomentum/predictedZMomentum - planeUmomentum/planeZmomentum )/tPullPuPzDenom);
      rootManager_->Get<TH1F*>( truthPullDirectory, Form("PvPz Truth Pull Plane %d", planeNum))->Fill((predictedVMomentum/predictedZMomentum - planeVmomentum/planeZmomentum )/tPullPvPzDenom);
      rootManager_->Get<TH1F*>( truthPullDirectory, Form("U Truth Pull Plane %d", planeNum))->Fill((predictedUPosition - planeUposition )/tPullUDenom);
      rootManager_->Get<TH1F*>( truthPullDirectory, Form("V Truth Pull Plane %d", planeNum))->Fill((predictedVPosition - planeVposition )/tPullVDenom);

      auto distPlotsDirectory = rootManager_->GetDir(rootManager_->GetDir(dirName_),"TrackLength");

      stringstream sss;
      sss << "Dist" << distInt;

      auto distNDirectory = rootManager_->GetDir(distPlotsDirectory,sss.str());
      auto distTruthPullDir = rootManager_->GetDir(distNDirectory, "Truth Pulls");

      rootManager_->Get<TH1F*>( distTruthPullDir, Form("1P Truth Pull Dist %d", distInt))->Fill((1./predictedMomentum - 1./planeTotalmomentum )/tPull1PDenom);
      rootManager_->Get<TH1F*>( distTruthPullDir, Form("PuPz Truth Pull Dist %d", distInt))->Fill((predictedUMomentum/predictedZMomentum - planeUmomentum/planeZmomentum )/tPullPuPzDenom);
      rootManager_->Get<TH1F*>( distTruthPullDir, Form("PvPz Truth Pull Dist %d", distInt))->Fill((predictedVMomentum/predictedZMomentum - planeVmomentum/planeZmomentum )/tPullPvPzDenom);
      rootManager_->Get<TH1F*>( distTruthPullDir, Form("U Truth Pull Dist %d", distInt))->Fill((predictedUPosition - planeUposition )/tPullUDenom);
      rootManager_->Get<TH1F*>( distTruthPullDir, Form("V Truth Pull Dist %d", distInt))->Fill((predictedVPosition - planeVposition )/tPullVDenom);

    }
  }


  void GeanePlots::FillGlobalPullHistos(int stationNum, const gm2strawtracker::TrackDetailArtRecord& trackDetail, art::Ptr<gm2truth::GhostDetectorArtRecord>& dummyHit){

    mf::LogTrace(name_) << "Begin filling pull histograms. \n";

    stringstream stationStream;
    stationStream << "TrackerStation[" << stationNum << "]";
    string stationStr = stationStream.str();

    // get the GEANE hits on the track
    auto track = trackDetail.geaneHits;

    Eigen::MatrixXd true0Covariance = track.covarianceTotalInverse.inverse(); // covarianceTotal is inverted from what it should be in the main code for speed improvements
    auto truthPullDirectory = rootManager_->GetDir(rootManager_->GetDir(dirName_), "TruthFitComparison");

    // truth parameters
    gm2geom::CoordSystem3Vector dHitPos(dummyHit->position.x(), dummyHit->position.y(), dummyHit->position.z() , "world");
    auto css = detCoordMap_.find("TrackerStation")->second;
    gm2geom::CoordSystem3Vector plane0Postion = dHitPos.transform(css, stationStr);

    double plane0Uposition = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*plane0Postion.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,1)*plane0Postion.y();
    double plane0Vposition = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,0)*plane0Postion.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,1)*plane0Postion.y();

    gm2geom::CoordSystem3Vector dHitMom(dummyHit->momentum.x(), dummyHit->momentum.y(), dummyHit->momentum.z() , "world");
    gm2geom::CoordSystem3Vector plane0Momentum = dHitMom.transform(css, stationStr, true);

    double plane0Umomentum = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*plane0Momentum.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,1)*plane0Momentum.y();
    double plane0Vmomentum = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,0)*plane0Momentum.x() + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,1)*plane0Momentum.y();
                       
    double plane0Zmomentum = plane0Momentum.z();
    double plane0Totalmomentum = plane0Momentum.mag();

    // predicted parameters
    double predicted0Momentum = geaneTrackUtils_.getPredMom(track, 0);
    double predicted0ZMomentum = geaneTrackUtils_.getPredZMom(track, 0);
    double predicted0UMomentum = geaneTrackUtils_.getPredUMom(track, 0);
    double predicted0VMomentum = geaneTrackUtils_.getPredVMom(track, 0);
    double predicted0UPosition = geaneTrackUtils_.getPredUPos(track, 0);
    double predicted0VPosition = geaneTrackUtils_.getPredVPos(track, 0);

    /////////////////////////////////////////////////////////////////////////////////////
    // fill plots
    /////////////////////////////////////////////////////////////////////////////////////

    // errors are from the sigma0 covariance matrix 
    rootManager_->Get<TH1F*>( truthPullDirectory, "1P Truth Pull" )->Fill((1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0)))); // using diagonal elements as the errors - should I be using a sum over columns? (for each row)
    rootManager_->Get<TH1F*>( truthPullDirectory, "PuPz Truth Pull" )->Fill((predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1))); 
    rootManager_->Get<TH1F*>( truthPullDirectory, "PvPz Truth Pull" )->Fill((predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2))); 
    rootManager_->Get<TH1F*>( truthPullDirectory, "U Truth Pull" )->Fill((predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3)))); 
    rootManager_->Get<TH1F*>( truthPullDirectory, "V Truth Pull" )->Fill((predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4)))); 

    auto pullPullDir = rootManager_->GetDir(truthPullDirectory, "Pull Pull");
    rootManager_->Get<TH2F*>( pullPullDir, "1P vs PuPz Pull" )->Fill((1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0))), (predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>( pullPullDir, "1P vs PvPz Pull" )->Fill((1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0))), (predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>( pullPullDir, "1P vs U Pull" )->Fill((1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0))), (predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3))));
    rootManager_->Get<TH2F*>( pullPullDir, "1P vs V Pull" )->Fill((1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0))), (predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4))));

    rootManager_->Get<TH2F*>( pullPullDir, "PuPz vs 1P Pull" )->Fill((predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1)), (1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0))));
    rootManager_->Get<TH2F*>( pullPullDir, "PuPz vs PvPz Pull" )->Fill((predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1)), (predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>( pullPullDir, "PuPz vs U Pull" )->Fill((predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1)), (predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3))));
    rootManager_->Get<TH2F*>( pullPullDir, "PuPz vs V Pull" )->Fill((predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1)), (predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4))));

    rootManager_->Get<TH2F*>( pullPullDir, "PvPz vs 1P Pull" )->Fill((predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2)), (1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0))));
    rootManager_->Get<TH2F*>( pullPullDir, "PvPz vs PuPz Pull" )->Fill((predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2)), (predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>( pullPullDir, "PvPz vs U Pull" )->Fill((predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2)), (predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3))));
    rootManager_->Get<TH2F*>( pullPullDir, "PvPz vs V Pull" )->Fill((predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2)), (predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4))));

    rootManager_->Get<TH2F*>( pullPullDir, "U vs 1P Pull" )->Fill((predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3))), (1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0))));
    rootManager_->Get<TH2F*>( pullPullDir, "U vs PuPz Pull" )->Fill((predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3))), (predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>( pullPullDir, "U vs PvPz Pull" )->Fill((predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3))), (predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>( pullPullDir, "U vs V Pull" )->Fill((predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3))), (predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4))));

    rootManager_->Get<TH2F*>( pullPullDir, "V vs 1P Pull" )->Fill((predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4))), (1./predicted0Momentum - 1./plane0Totalmomentum )/(1./1000.*sqrt(true0Covariance(0,0))));
    rootManager_->Get<TH2F*>( pullPullDir, "V vs PuPz Pull" )->Fill((predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4))), (predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum )/sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>( pullPullDir, "V vs PvPz Pull" )->Fill((predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4))), (predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum )/sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>( pullPullDir, "V vs U Pull" )->Fill((predicted0VPosition - plane0Vposition )/(10.*sqrt(true0Covariance(4,4))), (predicted0UPosition - plane0Uposition )/(10.*sqrt(true0Covariance(3,3))));

    auto NumDenomDir = rootManager_->GetDir(truthPullDirectory, "Num Denom");
    rootManager_->Get<TH1F*>( NumDenomDir, "1P Truth Pull Numerator" )->Fill(1./predicted0Momentum - 1./plane0Totalmomentum ); // using diagonal elements as the errors - should I be using a sum over columns? (for each row)
    rootManager_->Get<TH1F*>( NumDenomDir, "PuPz Truth Pull Numerator" )->Fill(predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum ); 
    rootManager_->Get<TH1F*>( NumDenomDir, "PvPz Truth Pull Numerator" )->Fill(predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum ); 
    rootManager_->Get<TH1F*>( NumDenomDir, "U Truth Pull Numerator" )->Fill(predicted0UPosition - plane0Uposition ); 
    rootManager_->Get<TH1F*>( NumDenomDir, "V Truth Pull Numerator" )->Fill(predicted0VPosition - plane0Vposition ); 

    rootManager_->Get<TH1F*>( NumDenomDir, "1P Truth Pull Denominator" )->Fill(1./1000.*sqrt(true0Covariance(0,0))); // using diagonal elements as the errors - should I be using a sum over columns? (for each row)
    rootManager_->Get<TH1F*>( NumDenomDir, "PuPz Truth Pull Denominator" )->Fill(sqrt(true0Covariance(1,1))); 
    rootManager_->Get<TH1F*>( NumDenomDir, "PvPz Truth Pull Denominator" )->Fill(sqrt(true0Covariance(2,2))); 
    rootManager_->Get<TH1F*>( NumDenomDir, "U Truth Pull Denominator" )->Fill(10.*sqrt(true0Covariance(3,3))); 
    rootManager_->Get<TH1F*>( NumDenomDir, "V Truth Pull Denominator" )->Fill(10.*sqrt(true0Covariance(4,4))); 

    rootManager_->Get<TH2F*>( NumDenomDir, "Denom vs Num 1P" )->Fill(1./predicted0Momentum - 1./plane0Totalmomentum , 1./1000.*sqrt(true0Covariance(0,0)));
    rootManager_->Get<TH2F*>( NumDenomDir, "Denom vs Num PuPz" )->Fill(predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum , sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>( NumDenomDir, "Denom vs Num PvPz" )->Fill(predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum , sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>( NumDenomDir, "Denom vs Num U" )->Fill(predicted0UPosition - plane0Uposition , 10.*sqrt(true0Covariance(3,3)));
    rootManager_->Get<TH2F*>( NumDenomDir, "Denom vs Num V" )->Fill(predicted0VPosition - plane0Vposition , 10.*sqrt(true0Covariance(4,4)));

    rootManager_->Get<TH2F*>( NumDenomDir, "Denom 1P vs P" )->Fill(predicted0Momentum, 1./1000.*sqrt(true0Covariance(0,0)));
    rootManager_->Get<TH2F*>( NumDenomDir, "Denom PuPz vs P" )->Fill(predicted0Momentum, sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>( NumDenomDir, "Denom PvPz vs P" )->Fill(predicted0Momentum, sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>( NumDenomDir, "Denom U vs P" )->Fill(predicted0Momentum, 10.*sqrt(true0Covariance(3,3)));
    rootManager_->Get<TH2F*>( NumDenomDir, "Denom V vs P" )->Fill(predicted0Momentum, 10.*sqrt(true0Covariance(4,4)));


    auto DenomVsDenomDir = rootManager_->GetDir(truthPullDirectory, "Fit Error Vs Fit Error");

    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom 1P vs 1P" )->Fill((1./1000.*sqrt(true0Covariance(0,0))), (1./1000.*sqrt(true0Covariance(0,0))));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom 1P vs PuPz" )->Fill((1./1000.*sqrt(true0Covariance(0,0))), sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom 1P vs PvPz" )->Fill((1./1000.*sqrt(true0Covariance(0,0))), sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom 1P vs U" )->Fill((1./1000.*sqrt(true0Covariance(0,0))), 10.*sqrt(true0Covariance(3,3)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom 1P vs V" )->Fill((1./1000.*sqrt(true0Covariance(0,0))), 10.*sqrt(true0Covariance(4,4)));

    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PuPz vs 1P" )->Fill(sqrt(true0Covariance(1,1)), (1./1000.*sqrt(true0Covariance(0,0))));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PuPz vs PuPz" )->Fill(sqrt(true0Covariance(1,1)), sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PuPz vs PvPz" )->Fill(sqrt(true0Covariance(1,1)), sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PuPz vs U" )->Fill(sqrt(true0Covariance(1,1)), 10.*sqrt(true0Covariance(3,3)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PuPz vs V" )->Fill(sqrt(true0Covariance(1,1)), 10.*sqrt(true0Covariance(4,4)));

    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PvPz vs 1P" )->Fill(sqrt(true0Covariance(2,2)), (1./1000.*sqrt(true0Covariance(0,0))));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PvPz vs PuPz" )->Fill(sqrt(true0Covariance(2,2)), sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PvPz vs PvPz" )->Fill(sqrt(true0Covariance(2,2)), sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PvPz vs U" )->Fill(sqrt(true0Covariance(2,2)), 10.*sqrt(true0Covariance(3,3)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom PvPz vs V" )->Fill(sqrt(true0Covariance(2,2)), 10.*sqrt(true0Covariance(4,4)));

    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom U vs 1P" )->Fill(10.*sqrt(true0Covariance(3,3)), (1./1000.*sqrt(true0Covariance(0,0))));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom U vs PuPz" )->Fill(10.*sqrt(true0Covariance(3,3)), sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom U vs PvPz" )->Fill(10.*sqrt(true0Covariance(3,3)), sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom U vs U" )->Fill(10.*sqrt(true0Covariance(3,3)), 10.*sqrt(true0Covariance(3,3)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom U vs V" )->Fill(10.*sqrt(true0Covariance(3,3)), 10.*sqrt(true0Covariance(4,4)));

    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom V vs 1P" )->Fill(10.*sqrt(true0Covariance(4,4)), (1./1000.*sqrt(true0Covariance(0,0))));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom V vs PuPz" )->Fill(10.*sqrt(true0Covariance(4,4)), sqrt(true0Covariance(1,1)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom V vs PvPz" )->Fill(10.*sqrt(true0Covariance(4,4)), sqrt(true0Covariance(2,2)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom V vs U" )->Fill(10.*sqrt(true0Covariance(4,4)), 10.*sqrt(true0Covariance(3,3)));
    rootManager_->Get<TH2F*>(DenomVsDenomDir, "Denom V vs V" )->Fill(10.*sqrt(true0Covariance(4,4)), 10.*sqrt(true0Covariance(4,4)));

    auto NumVsNumDir = rootManager_->GetDir(truthPullDirectory, "Delta Vs Delta");

    rootManager_->Get<TH2F*>(NumVsNumDir, "Num 1P vs 1P" )->Fill(1./predicted0Momentum - 1./plane0Totalmomentum, 1./predicted0Momentum - 1./plane0Totalmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num 1P vs PuPz" )->Fill(1./predicted0Momentum - 1./plane0Totalmomentum, predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num 1P vs PvPz" )->Fill(1./predicted0Momentum - 1./plane0Totalmomentum, predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num 1P vs U" )->Fill(1./predicted0Momentum - 1./plane0Totalmomentum, predicted0UPosition - plane0Uposition);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num 1P vs V" )->Fill(1./predicted0Momentum - 1./plane0Totalmomentum, predicted0VPosition - plane0Vposition);

    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PuPz vs 1P" )->Fill(predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum, 1./predicted0Momentum - 1./plane0Totalmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PuPz vs PuPz" )->Fill(predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum, predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PuPz vs PvPz" )->Fill(predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum, predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PuPz vs U" )->Fill(predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum, predicted0UPosition - plane0Uposition);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PuPz vs V" )->Fill(predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum, predicted0VPosition - plane0Vposition);

    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PvPz vs 1P" )->Fill(predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum, 1./predicted0Momentum - 1./plane0Totalmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PvPz vs PuPz" )->Fill(predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum, predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PvPz vs PvPz" )->Fill(predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum, predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PvPz vs U" )->Fill(predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum, predicted0UPosition - plane0Uposition);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num PvPz vs V" )->Fill(predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum, predicted0VPosition - plane0Vposition);

    rootManager_->Get<TH2F*>(NumVsNumDir, "Num U vs 1P" )->Fill(predicted0UPosition - plane0Uposition, 1./predicted0Momentum - 1./plane0Totalmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num U vs PuPz" )->Fill(predicted0UPosition - plane0Uposition, predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num U vs PvPz" )->Fill(predicted0UPosition - plane0Uposition, predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num U vs U" )->Fill(predicted0UPosition - plane0Uposition, predicted0UPosition - plane0Uposition);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num U vs V" )->Fill(predicted0UPosition - plane0Uposition, predicted0VPosition - plane0Vposition);

    rootManager_->Get<TH2F*>(NumVsNumDir, "Num V vs 1P" )->Fill(predicted0VPosition - plane0Vposition, 1./predicted0Momentum - 1./plane0Totalmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num V vs PuPz" )->Fill(predicted0VPosition - plane0Vposition, predicted0UMomentum/predicted0ZMomentum - plane0Umomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num V vs PvPz" )->Fill(predicted0VPosition - plane0Vposition, predicted0VMomentum/predicted0ZMomentum - plane0Vmomentum/plane0Zmomentum);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num V vs U" )->Fill(predicted0VPosition - plane0Vposition, predicted0UPosition - plane0Uposition);
    rootManager_->Get<TH2F*>(NumVsNumDir, "Num V vs V" )->Fill(predicted0VPosition - plane0Vposition, predicted0VPosition - plane0Vposition);

    mf::LogTrace(name_) << "End filling pull histograms. \n";

  }

} // namespace gm2strawtracker

DEFINE_ART_MODULE(gm2strawtracker::GeanePlots)
