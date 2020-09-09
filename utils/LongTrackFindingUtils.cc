// header file
#include "LongTrackFindingUtils.hh"

// namespace
using namespace gm2strawtracker;

// standard constructor
LongTrackFindingUtils::LongTrackFindingUtils()
   : ClusterSelectingUtils_() 
   , LSQFitter_()
   , CircleFitter_()
   , geom_()
   , name_("LongTrackFindingUtils")
   , cs_()
   , csMap_()
   , direction_("forward")
   , nmodules_(8)
   , multipleTrackCandidates_(false) 
   , removeClusters_(true)
   , makePlots_(false)
   , wireCentreInWorld_()
   , wireCentreInModule_()
   , wireCentreInStation_()
{
}

LongTrackFindingUtils::LongTrackFindingUtils(bool makePlots)
   : ClusterSelectingUtils_() 
   , LSQFitter_()
   , CircleFitter_()
   , geom_()
   , name_("LongTrackFindingUtils")
   , cs_()
   , csMap_()
   , direction_("forward")
   , nmodules_(8)
   , multipleTrackCandidates_(false) 
   , removeClusters_(true)
   , wirePositions_()
   , counter_(0)
   , makePlots_(makePlots)
   , wireCentreInWorld_()
   , wireCentreInModule_()
   , wireCentreInStation_()
{
  //book debug file
  if (makePlots_){
    std::string outputRootFileName = "LongTrackDebugPlots.root";
    outputRootFile_.reset( new TFile(outputRootFileName.c_str(),"recreate") );
    rm_.reset( new RootManager("LongTrackDebugPlots",outputRootFile_.get()) );
  }
}

void LongTrackFindingUtils::setFhiclCuts (fhicl::ParameterSet const & p ){
  trackFinderName_ = p.get<std::string>("trackFinderName");
  longTrackFinderParams_ = p.get<fhicl::ParameterSet>("longTrackFinderParams");
  minSeedsPerFit_ = longTrackFinderParams_.get<int>("minSeedsPerFit");
  maxSeedsToBeMultiTrackCand_ = longTrackFinderParams_.get<int>("maxSeedsToBeMultiTrackCand");
  shareSeeds_ = longTrackFinderParams_.get<bool>("shareSeeds");
  maxOverlapSeeds_ = longTrackFinderParams_.get<int>("maxOverlapSeeds");
  minPlanesPerFit_ = longTrackFinderParams_.get<int>("minPlanesPerFit");
  applyClusterSelection_ = longTrackFinderParams_.get<bool>("applyClusterSelection");
  candidateSeedsMaxTime_ = longTrackFinderParams_.get<double>("candidateSeedsMaxTime"); 
  candidateSeedsStrictTime_ = longTrackFinderParams_.get<double>("candidateSeedsStrictTime");
  candidateSeedsMaxWire_ = longTrackFinderParams_.get<int>("candidateSeedsMaxWire");
  maxWireSeparationForNeighbours_ = longTrackFinderParams_.get<int>("maxWireSeparationForNeighbours", 32 );
  maxWireSeparation_ = longTrackFinderParams_.get<int>("maxWireSeparation"); 
  candidateTracksMaxTime_ =   longTrackFinderParams_.get<double>("candidateTracksMaxTime");
  maxTimeSeparation_ = longTrackFinderParams_.get<double>("maxTimeSeparation");
  avgSeedRMSTime_ = longTrackFinderParams_.get<double>("avgSeedRMSTime");
  requiredFitConverges_ = longTrackFinderParams_.get<bool>("requiredFitConverges");

  fastFitterParams_      = p.get<fhicl::ParameterSet>("fastFitterParams");
  fitterName_ = fastFitterParams_.get<std::string>("fitterName");
  ClusterSelectingUtils_.setFhiclCuts( p );
  
  
} 

LongTrackFindingUtils::~LongTrackFindingUtils(){
  // close the root file
  if (makePlots_){
    rm_->WriteToFile();
    outputRootFile_->Close();
  }
}

TrackCandidateArtRecordCollection LongTrackFindingUtils::findTrackCandidates( StrawSeedPtrCollection& seeds )
{
  
  mf::LogDebug(name_) << "Enter LongTrackFindingUtils::findTrackCandidates using the model "
			 << trackFinderName_;
  
   TrackCandidateArtRecordCollection candidates;

   if( trackFinderName_ == "LongTrackFinder" ) {
     candidates = mergeSeedsToFormTrackCandidates(seeds);
   } else {
     mf::LogWarning(name_) << "This track finder name is not defined cannot form track candidates\n";
   } 

   mf::LogDebug(name_) << "Exit LongTrackFindingUtils::findTrackCandidates\n";
   return candidates;
}

//!======================================================================================================
//! merge seeds to form track candidates 
//! -- this function uses the LSQ minimizer to quantify the merging
//! -- for now an uniform magnetic field is assumed
//!======================================================================================================
TrackCandidateArtRecordCollection LongTrackFindingUtils::mergeSeedsToFormTrackCandidates( StrawSeedPtrCollection& seeds )
{

   counter_++;
   if (makePlots_ && counter_ > 10) makePlots_= false;
   mf::LogDebug(name_) << "Enter LongTrackFindingUtils::mergeSeedsToFormTrackCandidates with nseeds= " << seeds.size() << "\n"; 

   if( fitterName_ == "LSQMinimizer" ) { 
     // pass fhicl parameters to the fitter
     LSQFitter_.setFhiclCuts( fastFitterParams_ );

     // pass coordinate system tools to fitter
     LSQFitter_.setCoordSysData( cs_ );
   }

   //for straw plot
   if(wirePositions_.size() == 0){
     for(auto & wireID : geom_.getWireIDs()){
       if (geom_.getGlobalStation(wireID) != 2) continue;
       //std::string coordName = Form("Module%d:%d", wireID.getStation(), wireID.getModule());
       //gm2geom::CoordSystemsStoreData* newcs_ = &csMap_.find(coordName)->second;
       gm2geom::CoordSystem3Vector centreInWorld = wireCentreInWorld_.at(geom_.getGlobalWire(wireID));
       TEllipse* wire = new TEllipse(centreInWorld.z(), centreInWorld.x(),0.5,0.5);
       wire->SetLineColor(1);
       wirePositions_.push_back(wire);
     }
   }

   // container to store processed seeds and created tracks
   TrackCandidateArtRecordCollection preTrackCandidates, processTrackCandidates, trackCandidates;

   // check the size of the seed container 
   if( (int)seeds.size() < minSeedsPerFit_ ) {
     mf::LogDebug(name_) << "\tCannot create a track from nseeds= " << seeds.size() << "\n"; 
     return trackCandidates;
   }

   // more seeds in the event, more likely multi-particles are on the island 
   multipleTrackCandidates_ = (int)seeds.size() >= maxSeedsToBeMultiTrackCand_ ? true : false;

   // container to store seeds for re-organization
   std::vector< StrawSeedPtrCollection > seedsByModule;
   std::vector< int > indices; 
   organizeSeedsByModule(seeds,seedsByModule,indices);

   //flag up modules with high activity and remove them from the map
   vector<int> modulesToRemove;
   // loop over all modules
   for (auto i : indices){
     if (seedsByModule.at(i).size() > 2){
       //if neighbouring hits on same layer remove from selection for now
       bool remove = false;

       // for each module, make a map of wire numbers hit in each layer 0,1,2,3 = U0,U1,V0,V1
       std::map<int, vector<int> > layerAndWireNumbersMap;
       for (auto &seed : seedsByModule.at(i) ) {
	 for (auto &dig : seed->strawDigits) {
	   int viewNum = ( dig.get()->wireID.getView() == gm2strawtracker::u_view )? 0 : 1; 
	   int layerNum = 2 * viewNum + dig.get()->wireID.getLayer();
	   layerAndWireNumbersMap[ layerNum ].push_back( dig.get()->wireID.getWire());
	 }
       }

       //sorts wire numbers in ascending order
       for (auto iLay : {0,1,2,3}){
	 std::sort(layerAndWireNumbersMap[iLay].begin(), layerAndWireNumbersMap[iLay].end());

	 int prevWire = -1;
	 for (auto w: layerAndWireNumbersMap[iLay]) {
	   if (w - prevWire < 1) remove = true;
	   prevWire = w;
	 }
       }
       if (remove) modulesToRemove.push_back(i);
     }
   }
   
   for (auto i : modulesToRemove) mf::LogDebug(name_) << "module to remove: " << i << "\n";

   // Make Histogram that we'll add TMarkers to for wires and this island, only for debugging
   TH1F *strawMapPlot = NULL; 
   //TH1F* 
   if (makePlots_){
     strawMapPlot = new TH1F("plot1",";Ring z [mm];Ring x [mm]", 100, -7025, -6660);
     strawMapPlot->GetYaxis()->SetRangeUser(200,1400);
     strawMapPlot->SetStats(0);
     // Add ellipses for wire positions
     for(auto& ellipse : wirePositions_ ) strawMapPlot->GetListOfFunctions()->Add(ellipse);
   }
  
   // determine the direction
   int beginIndex = indices.front();
   int endIndex   = indices.back();

   direction_  = seedsByModule[beginIndex].size() <= seedsByModule[endIndex].size() ? "forward" : "backward";
   auto mindex = direction_ == "forward" ? beginIndex : endIndex;

   // get the seeds in the first registered module
   StrawSeedPtrCollection initialSeeds = seedsByModule[mindex];

   // container to hold used and unused seeds
   StrawSeedPtrCollection usedSeeds, unusedSeeds;
 
   //DebugPlots
   TH1F *hinitial = NULL; 
   auto eventFolderName = Form("event_%d", counter_);
   if (makePlots_){
     hinitial = (TH1F*) strawMapPlot->Clone();
     TH1F* hignored = (TH1F*) strawMapPlot->Clone();
     for (auto& seedCol : seedsByModule) {
       for (auto& seed : seedCol )  {
	 hinitial = strawPlot(seed->strawDigits, hinitial);
	 if (std::find( modulesToRemove.begin(), modulesToRemove.end(), seed.get()->frontModule ) != modulesToRemove.end()) {	 
	   hignored = strawPlot(seed->strawDigits, hignored, 3);
	 }
       }
     }
     hinitial->SetName("initialHits");
     hignored->SetName("ignoredHits");
     rm_->Add(eventFolderName, hinitial->Clone());
     rm_->Add(eventFolderName, hignored->Clone());
   }

   // processed all seeds into track candidates
   mf::LogDebug(name_) << "\tProcessed seeds in track candidates starting at module= " 
                          << initialSeeds[0].get()->frontModule << " and seeds size= " 
                          << initialSeeds.size() 
                          << ", searching in the direction= " << direction_
                          << "\n";

   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! BUILD TRACK CANDIDATES
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   int iSeed(0);
   while( !initialSeeds.empty() ) {
     iSeed++;
     // create a collection of track candidates
     TrackCandidateArtRecordCollection buildCandidates;
      
     // get reference seed
     auto refSeed   = *initialSeeds.begin();
     auto refModule = refSeed.get()->frontModule;

     //DebugPlots
     std::string seedFolderName( Form("/seed%d", iSeed) );
     if (makePlots_){
       TH1F* initSeed = (TH1F*) strawMapPlot->Clone();
       initSeed = strawPlot(refSeed->strawDigits, hinitial, 3);
       initSeed->SetName("initialSeed");
       rm_->Add(eventFolderName + seedFolderName, initSeed->Clone());
     }
      
     // create at track candidate
     TrackCandidateArtRecord candidate;
     candidate.strawSeeds.push_back( refSeed );
     mf::LogDebug(name_) << "\tFirst seed= [ " << refSeed.get()->frontModule << ", " << refSeed.get()->backModule << " ]\n";      

     // add to collection
     buildCandidates.push_back( candidate );

     // remove seed from container
     initialSeeds.erase( initialSeeds.begin() );

     // build candidates
     int keepBuildingIndex(0);
     bool keepBuilding = true;
     while( keepBuilding ) {
       keepBuildingIndex++;
       // set reference seed
       auto refSeed = buildCandidates.back().strawSeeds.back();

       // safety if the seed is composed of clusters on different modules
       int moduleIndex = 1;
       if( refSeed.get()->frontModule != refSeed.get()->backModule )
         moduleIndex = 2; 

       // get the number of the neighboring modules 
       int neighModule = direction_ == "forward" ? refModule+moduleIndex : refModule-moduleIndex;
   
       //if ignore layer/view/module
       while (std::find( modulesToRemove.begin(), modulesToRemove.end(), neighModule ) != modulesToRemove.end()) {
	 mf::LogDebug(name_) << "removal test, skipping mod " << refSeed.get()->frontModule + moduleIndex << "\n";
	 moduleIndex++;
	 neighModule = direction_ == "forward" ? refModule+moduleIndex : refModule-moduleIndex;
       }

       if( neighModule == nmodules_ || neighModule < 0 ) {
         keepBuilding = false;
         break;
       }

       // get neighboring seeds
       StrawSeedPtrCollection neighSeeds;
       auto seedVect = seedsByModule.at( neighModule );
       mf::LogDebug(name_) << "\t\tnumber of neighboring seeds for reference module (" << refModule << "): " << seedVect.size() << "\n";

       // determine whether to find the best neighboring seeds in a high occupancy region of seeds
       bool isManyCloseSeedsInSameModule = false;
       if( seedVect.size() > 3 ) {
	 //additional check to see if these seeds are anywhere near the reference seed
	 int nNearSeed = 0;
	 
	 //TODO work out the average correctly, for now use front or back wire number
	 auto refWireAvg = (direction_=="forward")? refSeed.get()->strawDigits.back().get()->wireID.getWire() : refSeed.get()->strawDigits.front().get()->wireID.getWire();
	 for (auto seed : seedVect){
	   auto candWireAvg = direction_=="forward" ? seed.get()->strawDigits.back().get()->wireID.getWire() : seed.get()->strawDigits.front().get()->wireID.getWire();
	   //TODO make steerable

	   //if (fabs(refWireAvg - candWireAvg) < maxWireNumbersBetweenSeeds_) nNearSeed++;
	   if (fabs(refWireAvg - candWireAvg) < 4) nNearSeed++;
	 }
	 
	 isManyCloseSeedsInSameModule = (nNearSeed > 3) && overlappingSeeds(seedVect);
	 if( isManyCloseSeedsInSameModule ) {
	   keepBuilding = false;
	   break;
	 }
       }


       // determine whether to allow track candidates to share seeds
       if( shareSeeds_ && seedVect.size() > 1 && !preTrackCandidates.empty() ) {
         mf::LogDebug(name_) << "\t\tdetermine whether to allow track candidates to share seeds\n";
         checkSeedIsOnTrackCandidate(seedVect,preTrackCandidates);
       }
 
       // get the avg time of the track candidate
       double timeDiff = 9999.;
       if( buildCandidates.back().strawSeeds.size() >= 2 ) {
         auto trackSeeds = buildCandidates.back().strawSeeds;
         timeDiff = getTimeDifferenceBetweenSeeds(trackSeeds);
       } 

       // get the best neighboring seeds
       // seedVect contains all seeds in neighboring module, refseed the seed in module before and neighSeeds gets filled in the function
       getNeighboringSeeds(seedVect,refSeed,neighSeeds,timeDiff,isManyCloseSeedsInSameModule);
       if( neighSeeds.empty() ) {
         keepBuilding = false;
         break;
       }

       //DebugPlots, all possible neighbouring seeds
       if (makePlots_){
	 int neighIndex(0);
	 for (auto &nSeed: seedVect){
	   neighIndex++;
	   TH1F* neighSeed = (TH1F*) strawMapPlot->Clone();
	   neighSeed = strawPlot(nSeed->strawDigits, hinitial, 3);
	   neighSeed->SetName(Form("neighSeed_%d_%d", keepBuildingIndex, neighIndex));
	   rm_->Add(eventFolderName + seedFolderName, neighSeed->Clone());
	 }
       }
       
       // create new container for bookkeeping and organization
       unsigned int index = 0;
       unsigned int nsize = buildCandidates.size()*neighSeeds.size();
       TrackCandidateArtRecordCollection candidateVect(nsize);

       // make copy of each candidate
       for(unsigned int i = 0; i < buildCandidates.size(); ++i) {
         TrackCandidateArtRecord copyCandidate = buildCandidates.at(i);

         for(unsigned int j = 0; j < neighSeeds.size(); ++j,++index) {
           candidateVect.at(index) = copyCandidate;
         }
       }

       // copy assignment
       buildCandidates.clear();
       buildCandidates = candidateVect;

       // merge seeds
       unsigned int icol = 0;
       for(unsigned int row = 0; row < buildCandidates.size(); ++row) {
         mf::LogDebug(name_) << "\t\tadding seed: [ " << neighSeeds[icol].get()->frontModule 
                                << ", " << neighSeeds[icol].get()->backModule << " ]\n";

         buildCandidates[row].strawSeeds.push_back( neighSeeds[icol] );
         icol += 1;
         if( icol == neighSeeds.size() )
            icol = 0;
       } 

       // increment
       refModule = neighModule;
     } // end loop over building track candidate using this initial seed, e.g. keep building loop
     
     // added candidates to collection
     for(auto cand : buildCandidates) {
       preTrackCandidates.push_back( cand );

       //DebugPlots
       if (makePlots_){
	 std::string candFolderName( Form("/candidate_%d", int(preTrackCandidates.size()) ) );
	 TH1F* hcand = (TH1F*) strawMapPlot->Clone();
	 for(auto& iseed : cand.strawSeeds) hcand = strawPlot(iseed->strawDigits, hcand, 3);
	 hcand->SetName("preCandidate");
	 rm_->Add(eventFolderName + candFolderName, hcand->Clone());
       }
       for(auto iseed : cand.strawSeeds) {
	 if( std::find(usedSeeds.begin(),usedSeeds.end(),iseed) == usedSeeds.end() )
           usedSeeds.push_back( iseed );
       } //end of seed bookkeeping
       
     } // end of finding candidates  
     
     // check if more track candidates can be made
     if( initialSeeds.empty() ) {
       for(auto& iseed : seeds) {
         if( std::find(usedSeeds.begin(),usedSeeds.end(),iseed) == usedSeeds.end() ){
	   //don't push back unusable seeds, do that later...
	   if (std::find( modulesToRemove.begin(), modulesToRemove.end(), iseed.get()->frontModule ) == modulesToRemove.end())
	     unusedSeeds.push_back( iseed );
	 }
       }
       
       mf::LogDebug(name_) << "\t  build track candidates= " << preTrackCandidates.size() 
                              << ", with unused seeds= " << unusedSeeds.size() << "\n";
       
       if( (int)unusedSeeds.size() >= minSeedsPerFit_  ) {
         std::sort(unusedSeeds.begin(),unusedSeeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter());
	 
         organizeSeedsByModule(unusedSeeds,seedsByModule,indices);

         beginIndex = indices.front();
         endIndex   = indices.back();

         direction_  = seedsByModule[beginIndex].size() <= seedsByModule[endIndex].size() ? "forward" : "backward";
         mindex = direction_ == "forward" ? beginIndex : endIndex;
	 
         initialSeeds = seedsByModule[ mindex];
         unusedSeeds.clear();
	 
         mf::LogDebug(name_) << "\tProcessed seeds in track candidates starting at module= " 
                                << initialSeeds[0].get()->frontModule << " and seeds size= " 
                                << initialSeeds.size() 
                                << ", searching in the direction= " << direction_
                                << "\n";
       } // end of unused seed check
     } // end of checking for more track candidates
     
   } // end loop over creating track candidates from all seeds, e.g. initial seed loop
   
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! REMOVE preTRACK CANDIDATES WITH ONE SEED
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   auto ipreCand = preTrackCandidates.begin();
   while( ipreCand != preTrackCandidates.end() ) {
      if( (*ipreCand).strawSeeds.size() == 1 ) {
        unusedSeeds.push_back( (*ipreCand).strawSeeds.front() );
        ipreCand = preTrackCandidates.erase(ipreCand);
      } else {
        ++ipreCand;
      }
   }

   mf::LogDebug(name_) << "\tAfter building and removal, found pre- track candidates= " << preTrackCandidates.size() 
                          << ", with unused seeds= " << unusedSeeds.size() << "\n";
   
   //DebugPlots
   if (makePlots_){
     std::string fName( "/AfterSingleSeedRemoval" );
     int itmp(0);
     for (auto &cand: preTrackCandidates){
       itmp++;
       TH1F* hcand = (TH1F*) strawMapPlot->Clone();
       for(auto& iseed : cand.strawSeeds) hcand = strawPlot(iseed->strawDigits, hcand, 3);
       hcand->SetName(Form("preCandidate_%d", itmp));
       rm_->Add(eventFolderName + fName, hcand->Clone());
     }
   }
   
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! SORT THE SEEDS OF THE preTRACK CANDIDATES
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   for(auto& candidate : preTrackCandidates) {
     std::sort(candidate.strawSeeds.begin(),candidate.strawSeeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter());
     mf::LogDebug(name_) << "\t\tpre- track candidate's number of seeds= " << candidate.strawSeeds.size() << "\n";
   }

   std::sort(unusedSeeds.begin(),unusedSeeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter());

   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! MERGE NEIGHBORING UN-USED SEEDS TO preTRACK CANDIDATES
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   // check if unused seeds can be added onto the track candidates
   mf::LogDebug(name_) << "\tMerge neighboring unused seeds onto the pre-track candidates\n";

   // determine if the seeds are good candidates for merging -- if the track candidate runs into a high occupancy module
   // do not attempt to add those seeds on the track candidate -- this should be done at a later stage
   if( !unusedSeeds.empty() ) {
     organizeSeedsByModule(unusedSeeds,seedsByModule,indices);
   }
 
   auto iunusedSeed = unusedSeeds.begin();
   while( iunusedSeed != unusedSeeds.end() ) {

      auto mergeSeed   = false;
      auto seed        = *iunusedSeed;
      auto frontModule = seed.get()->frontModule; 

      if( seedsByModule.at(frontModule).size() == 1 ) { 
        bool mergeBrokenTracks = true;

        for(unsigned int i = 0; i < preTrackCandidates.size(); ++i) {
          if( trackCandidateSeedNeighbors(preTrackCandidates[i],seed,mergeBrokenTracks) ) {
            mergeSeed = true;
            preTrackCandidates[i].strawSeeds.push_back( seed );
            std::sort(preTrackCandidates[i].strawSeeds.begin(),preTrackCandidates[i].strawSeeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter() );
            mf::LogDebug(name_) << "\t\ti= " << i << ", added unused seed to track candidate= " << preTrackCandidates[i].strawSeeds.size() << "\n";
          }
        }
      } // end of deciding to enter a high or low occupancy region of the tracker

      if( mergeSeed ) 
        iunusedSeed = unusedSeeds.erase(iunusedSeed);
      else 
        iunusedSeed++;
   } // end of merging un-used seeds to track candidates
   
   //DebugPlots
   if (makePlots_){
     std::string fName = "/AfterMergeOfNeighbouringUnusedSeeds";
     int itmp = 0;
     for (auto &cand: preTrackCandidates){
       itmp++;
       TH1F* hcand = (TH1F*) strawMapPlot->Clone();
       for(auto& iseed : cand.strawSeeds) hcand = strawPlot(iseed->strawDigits, hcand, 3);
       hcand->SetName(Form("preCandidate_%d", itmp));
       rm_->Add(eventFolderName + fName, hcand->Clone());
     }
   }

   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! MERGE OVERLAPPING UN-USED SEEDS TO preTRACK CANDIDATES
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   // check if unused seeds can be added onto the track candidates
   mf::LogDebug(name_) << "\tMerge overlapping unused seeds ( " << unusedSeeds.size() << " ) onto the pre-track candidates\n";

   // re-organize the input list of seeds, this is needed to determine if the unused seed
   // is in a high occupancy region, we only want to merge seeds in isolated regions
   if( !unusedSeeds.empty() ) {
     organizeSeedsByModule(seeds,seedsByModule,indices);
   }

   auto junusedSeed = unusedSeeds.begin();
   while( junusedSeed != unusedSeeds.end() ) {

      bool mergeSeed = false;
      auto seed      = *junusedSeed;

      auto moduleIdx     = seed.get()->frontModule;
      auto seedsInRegion = seedsByModule.at(moduleIdx);
      
      bool isSeedInHighOccupancyRegion = seedsInRegion.size() > 2 ? true : false;

      if( !isSeedInHighOccupancyRegion ) {
        for(unsigned int i = 0; i < preTrackCandidates.size(); ++i) {
          if( trackCandidateSeedOverlaps(preTrackCandidates[i],seed) ) {
            mergeSeed = true;
            preTrackCandidates[i].strawSeeds.push_back( seed );
            std::sort(preTrackCandidates[i].strawSeeds.begin(),preTrackCandidates[i].strawSeeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter() );
            mf::LogDebug(name_) << "\t\ti= " << i << ", added unused seed to track candidate= " << preTrackCandidates[i].strawSeeds.size() << "\n";
          }
        }
      } // end of checking region occupancy
 
      if( mergeSeed ) {
        junusedSeed = unusedSeeds.erase(junusedSeed);
      } else {
        junusedSeed++;
      }

   } // end of merging un-used seeds to track candidates

   //DebugPlots
   if (makePlots_){
     std::string fName = "/AfterMergeOfOverlappingUnusedSeeds";
     int itmp = 0;
     for (auto &cand: preTrackCandidates){
       itmp++;
       TH1F* hcand = (TH1F*) strawMapPlot->Clone();
       for(auto& iseed : cand.strawSeeds) hcand = strawPlot(iseed->strawDigits, hcand, 3);
       hcand->SetName(Form("preCandidate_%d", itmp));
       rm_->Add(eventFolderName + fName, hcand->Clone());
     }
   }
   
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! MERGE BROKEN pre-TRACK CANDIDATES
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   // check if the unused seeds or track candidates are part of a broken track candidate
   mf::LogDebug(name_) << "\tMerge broken pre-track candidates= " << preTrackCandidates.size() << ", with unused seeds= " << unusedSeeds.size() << "\n";

   for(unsigned int i = 0; i < preTrackCandidates.size(); ++i) {
     auto refCandidate = preTrackCandidates[i];

     for(unsigned int j = i; j < preTrackCandidates.size(); ++j) {
       auto candidate = preTrackCandidates[j];

       if( candidate == refCandidate )
         continue;

       auto overlaps  = trackSeedsOverlaps(refCandidate.strawSeeds,candidate.strawSeeds);
       auto neighbors = trackCandidatesNeighbors(refCandidate,candidate);
       auto isGap     = isGapBetweenTrackCandidates(refCandidate,candidate,unusedSeeds);

       if( !overlaps && neighbors && isGap ) {
         refCandidate.strawSeeds.insert(refCandidate.strawSeeds.end(),candidate.strawSeeds.begin(),candidate.strawSeeds.end());
         std::sort(refCandidate.strawSeeds.begin(),refCandidate.strawSeeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter());
         mf::LogDebug(name_) << "\t\tmerged track onto track candidate with seeds= " << refCandidate.strawSeeds.size() << "\n";
       } 
     }

     processTrackCandidates.push_back( refCandidate );
   }  // end of merging broken tracks

   //DebugPlots - here we have moved to process track candidates
   if (makePlots_){
     std::string fName = "/AfterMergeOfBrokenTrack";
     int itmp = 0;
     for (auto &cand: processTrackCandidates){
       itmp++;
       TH1F* hcand = (TH1F*) strawMapPlot->Clone();
       for(auto& iseed : cand.strawSeeds) hcand = strawPlot(iseed->strawDigits, hcand, 3);
       hcand->SetName(Form("preCandidate_%d", itmp));
       rm_->Add(eventFolderName + fName, hcand->Clone());
     }
   }


   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! REMOVE FAKE TRACK CANDIDATES
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   // remove unneccessary tracks
   mf::LogDebug(name_) << "\tRemove duplicated track candidates ( total= " << processTrackCandidates.size() << " )\n";

   while( !processTrackCandidates.empty() ) {

     // get reference track candidate
     auto refCandidate = *processTrackCandidates.begin();

     // added processed track candidates
     trackCandidates.push_back( refCandidate );

     // remove track candidate from container
     auto icand = processTrackCandidates.erase( processTrackCandidates.begin() );

     // remove overlapping tracks
     while( icand != processTrackCandidates.end() ) {

       // check if candidate can be removed
       StrawSeedPtrCollection seedNotOnRefCandidate;     
       for(auto iseed : icand->strawSeeds) {
         auto it = std::find(trackCandidates.back().strawSeeds.begin(),trackCandidates.back().strawSeeds.end(),iseed);
         if( it == trackCandidates.back().strawSeeds.end() )
           seedNotOnRefCandidate.push_back( iseed );
       }

       // seed condition
       if( (int)seedNotOnRefCandidate.size() > maxOverlapSeeds_ ) {
         ++icand;
         continue;
       }

       // determine to merge track candidates
       bool mergeCandidates = false;
       for(auto iseed : seedNotOnRefCandidate) {
         if( trackCandidateSeedShareCluster(*icand,iseed) )
           mergeCandidates = true;
       }

       if( mergeCandidates ) {
         std::for_each(seedNotOnRefCandidate.begin(),seedNotOnRefCandidate.end(),[&](auto iseed){ trackCandidates.back().strawSeeds.push_back( iseed ); });
         icand = processTrackCandidates.erase( icand );
       } else {
         ++icand;
       }
     }
   } // end remove overlapping track candidates

   //DebugPlots - here we have moved to track candidates
   if (makePlots_){
     std::string fName = "/AfterRemoveOverlappingTracks";
     int itmp = 0;
     for (auto &cand: trackCandidates){
       itmp++;
       TH1F* hcand = (TH1F*) strawMapPlot->Clone();
       for(auto& iseed : cand.strawSeeds) hcand = strawPlot(iseed->strawDigits, hcand, 3);
       hcand->SetName(Form("preCandidate_%d", itmp));
       rm_->Add(eventFolderName + fName, hcand->Clone());
     }
   }
 
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! CHECK PLANE CONDITION
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   // determine if candidate meets the plane condition
   mf::LogDebug(name_) << "\tDetermine if the track candidates ( " << trackCandidates.size() << " ) meet the plane condition\n";

   auto icand = trackCandidates.begin();
   while( icand != trackCandidates.end() ) {

     auto candidate = *icand;
     std::vector< std::string > planeNames;

     for(auto seed : candidate.strawSeeds) {
       for(auto clus : seed.get()->strawClusters) {
         std::string name( Form("Station%dModule%dPlane%d", clus.get()->station, clus.get()->module, clus.get()->view) );
         if( std::find(planeNames.begin(),planeNames.end(),name) == planeNames.end() ) 
           planeNames.push_back( name );
       }
     }

     if( int(planeNames.size()) < minPlanesPerFit_ ) {
       auto it = std::find(trackCandidates.begin(),trackCandidates.end(),candidate);
       if( it != trackCandidates.end() ) {
         icand = trackCandidates.erase(it);
         mf::LogDebug(name_) << "\t\tremoved candidate with planes= " << planeNames.size() << "\n"; 
       }
       
     } else {
       ++icand;
     }
   }

   //DebugPlots
   if (makePlots_){
     std::string fName = "/AfterPlaneCondition";
     int itmp = 0;
     for (auto &cand: trackCandidates){
       itmp++;
       TH1F* hcand = (TH1F*) strawMapPlot->Clone();
       for(auto& iseed : cand.strawSeeds) hcand = strawPlot(iseed->strawDigits, hcand, 3);
       hcand->SetName(Form("preCandidate_%d", itmp));
       rm_->Add(eventFolderName + fName, hcand->Clone());
     }
   }

   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! REMOVE CLUSTERS FROM THE TRACK CANDIDATE THAT DONT BELONG ON
   //! THE TRACK
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   //perform cluster selection here
   if (applyClusterSelection_) {
     mf::LogDebug(name_) << "\tPerform the clustering selection on the track candidates= " << trackCandidates.size() << "\n";
     ClusterSelectingUtils_.trackGrouping( trackCandidates, removeClusters_ );
   }

   //DebugPlots
   if (makePlots_){
     std::string fName = "/AfterRemovalOfClusters";
     int itmp = 0;
     for (auto &cand: trackCandidates){
       itmp++;
       TH1F* hcand = (TH1F*) strawMapPlot->Clone();
       for(auto& iseed : cand.strawSeeds) hcand = strawPlot(iseed->strawDigits, hcand, 3);
       hcand->SetName(Form("preCandidate_%d", itmp));
       rm_->Add(eventFolderName + fName, hcand->Clone());
     }
   }

   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   //!
   //! GET THE INITIAL FIT PARAMETERS
   //!
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   // get fit parameters
   mf::LogDebug(name_) << "\tGet the initial fit parameters of the track candidates\n";

   for(auto& candidate : trackCandidates) {

     if( fitterName_ == "LSQMinimizer" )
       callLSQFitter(candidate);
     else if( fitterName_ == "SimpleCircleFitter" )
       callSimpleCircleFitter(candidate);
 
   } // end loop over track candidates

   //DebugPlots
   if (makePlots_){
     std::string fName = "/FinalCandidates";
     int itmp = 0;
     for (auto &cand: trackCandidates){
       itmp++;
       TH1F* hcand = (TH1F*) strawMapPlot->Clone();
       for(auto& iseed : cand.strawSeeds) hcand = strawPlot(iseed->strawDigits, hcand, 3);
       hcand->SetName(Form("preCandidate_%d", itmp));
       rm_->Add(eventFolderName + fName, hcand->Clone());
     }
   }
   mf::LogDebug(name_) << "Exit LongTrackFindingUtils::mergeSeedsToFormTrackCandidates with n track candidates= " << trackCandidates.size() << "\n";
   return trackCandidates;
}


//--------------------------------------------------------------------------
// reconstructed the track's momentum based on the curvature
//--------------------------------------------------------------------------
void LongTrackFindingUtils::getTrackMomentum( TrackCandidateArtRecord& track )
{
   mf::LogDebug(name_) << "Enter LongTrackFindingUtils::getTrackMomentum\n";

   if( !track.helixFit ) 
     return;

   if( track.helix.kappa.value == 0 ) {
     mf::LogWarning(name_) << "The curvature has a zero value, cannot reconstruct the momentum\n";
     return;
   }

   /*
    The helix parameters were computed using the ring coordinate system. 
    */    

   const double pt  = (1.0 / track.helix.kappa.value)*1000.; 
   const double phi = track.helix.phi;

   mf::LogDebug(name_) << "\tpt= " << pt << ", turning angle= " << phi << "\n";

   // set momentum in units of MeV
   track.momentum.setX( pt*(-sin(track.helix.phi0.value+phi)) );
   track.momentum.setZ( pt*cos(track.helix.phi0.value+phi) );
   track.momentum.setY( pt*track.helix.tlambda.value );
   track.momentum.coordSystemName = "world";

   mf::LogDebug(name_) << "\thelix= " << track.helix << "\n";
   mf::LogDebug(name_) << "\tmomentum= " << track.momentum << "\n";
   mf::LogDebug(name_) << "Exit LongTrackFindingUtils::getTrackMomentum\n";
}


//-----------------------------------------------------------
// organize the seeds
//-----------------------------------------------------------
void LongTrackFindingUtils::organizeSeedsByModule( StrawSeedPtrCollection& seeds, std::vector< StrawSeedPtrCollection>& seedsByModule, std::vector< int >& indices )
{
   seedsByModule.clear();
   indices.clear();

   std::vector< StrawSeedPtrCollection> tmp(nmodules_);
   for(auto& seed : seeds) {
     int idx = seed.get()->frontModule;

     tmp.at(idx).push_back( seed );

     if( std::find(indices.begin(),indices.end(),idx) == indices.end() )
       indices.push_back( idx );
   }

   std::sort(indices.begin(),indices.end(),[](auto i, auto j){ return i < j; });

   seedsByModule = tmp;

   tmp.clear();
   return;
}


//-----------------------------------------------------------
// select the most probable neighboring seed
//-----------------------------------------------------------
void LongTrackFindingUtils::getNeighboringSeeds( StrawSeedPtrCollection& seeds, art::Ptr<StrawSeedArtRecord>& iseed, StrawSeedPtrCollection& oseeds,
                                                 double avgTimeTrackCandidate, bool isManyCloseSeedsInSameModule )
{
   auto wireSeparation = 9999;
   auto timeSeparation = avgTimeTrackCandidate + 25;

   //------ get the reference digit information
   double refTimeAvg = 0;
   std::for_each(iseed.get()->strawDigits.begin(),
                 iseed.get()->strawDigits.end(),
                 [&](auto& idig){ refTimeAvg += idig.get()->calTime; });

   refTimeAvg = refTimeAvg/double(iseed.get()->strawDigits.size());

   auto idigits  = iseed.get()->strawDigits;
   std::sort(idigits.begin(),idigits.end(),StrawTrackerSort::StrawDigitArtPtrRowSorter());

   //TODO U and V can be quite different wireIDs, get separate averages
   float refWireAvg = 0.;
   for (auto &dig : idigits) refWireAvg += float(dig.get()->wireID.getWire());
   refWireAvg /= float(idigits.size());

   //------- get the neighboring digit information
   std::vector< DataObject > neighDataVect;
   for(auto& seed : seeds) {
     DataObject data;
     data.seed = seed;

     double neighTimeAvg = 0;
     std::for_each(seed.get()->strawDigits.begin(),
                   seed.get()->strawDigits.end(),
                   [&](auto& idig){ neighTimeAvg += idig.get()->calTime; });
     data.meanTime = neighTimeAvg/double(seed.get()->strawDigits.size());

     auto idigits  = seed.get()->strawDigits;
     std::sort(idigits.begin(),idigits.end(),StrawTrackerSort::StrawDigitArtPtrRowSorter());

     //TODO meanWire is an int
     float neighWireAvg = 0.;
     for (auto &dig : idigits) neighWireAvg += float(dig.get()->wireID.getWire());
     neighWireAvg /= idigits.size();
     data.meanWire = neighWireAvg;

     neighDataVect.push_back( data );
   }

   //----- get geometric information
   bool candidateSeedsVeryCloseInTime = false;
   bool candidateSeedsCloseInTime    = false;
   bool candidateSeedsCloseInSpace   = false;

   for(unsigned int i = 0; i < neighDataVect.size(); ++i) {
     for(unsigned int j = 0; j < neighDataVect.size(); ++j) {
       if( i == j ) continue;
       if( fabs(neighDataVect[j].meanTime-neighDataVect[i].meanTime) <= candidateSeedsMaxTime_ ) 
         candidateSeedsCloseInTime    = true;
       if( fabs(neighDataVect[j].meanTime-neighDataVect[i].meanTime) <= candidateSeedsStrictTime_ )
         candidateSeedsVeryCloseInTime = true;
       if( fabs(neighDataVect[j].meanWire-neighDataVect[i].meanWire) <= candidateSeedsMaxWire_ )
         candidateSeedsCloseInSpace   = true;
     }
   } 

   //--------------------------------------------------------------
   //
   //  FIND THE BEST NEIGHBORING SEED
   //
   //--------------------------------------------------------------  

   // if we think there is only one candidate then add this neighboring seed to the output candidate
   if( seeds.size() == 1 && !multipleTrackCandidates_ ) {
     auto neighWireAvg = neighDataVect[0].meanWire;
     if( direction_ == "forward" && refWireAvg-neighWireAvg >= 0 ) 
       oseeds.push_back( seeds.front() );
     else if( direction_ == "backward" && refWireAvg-neighWireAvg <= 0 )
       oseeds.push_back( seeds.front() );

   } else {
     //multipleTrackCandidates_ = true;
 
     oseeds.clear();
     art::Ptr<StrawSeedArtRecord> oseed;
     bool addedSeed = false;
     
     for(auto& n : neighDataVect) {
       auto neighWireAvg = n.meanWire;
       auto neighTimeAvg = n.meanTime;
       auto seed         = n.seed;
       auto digits       = seed.get()->strawDigits;
       
       // IF THERE ARE 2 NEIGHBORING SEEDS (ON SAME MODULE) THAT ARE CLOSE DO SPECIAL CASE
       if( candidateSeedsCloseInTime && candidateSeedsCloseInSpace ) { //TODO USE FITTER
         bool seedsOverlap = overlappingSeeds(iseed,seed);
         mf::LogDebug(name_) << "\t\tThe seeds are close in time and space ( overlap= " << seedsOverlap << " )\n";
	 
	 /*
	   if( seedsOverlap ) //Turned this off the fitter should removed hits that don't belong
           continue;
	 */
       }
       
       // TODO: make the diffrerence in incorrect direction steerable
       bool wireMoveDirection = false;
       if( direction_ == "forward" && refWireAvg-neighWireAvg >= -2 ) 
         wireMoveDirection = true;
       else if( direction_ == "backward" && refWireAvg-neighWireAvg <= 2 ) 
         wireMoveDirection = true;
       
       if( !wireMoveDirection )
         continue;
       
       //IF NEIGHBORING SEEDS ARE NOT CLOSE IN TIME, AND THE AVG TIME IS CLOSE TO THE REFERENCE SEED AVG TIME
       if( !candidateSeedsCloseInTime ) {
	 if( fabs(refTimeAvg-neighTimeAvg) < timeSeparation ) {
	   //ADDED CHECK IF THEY ARE CLOSE IN WIRE NUMBER
	   if( fabs(neighWireAvg - refWireAvg) < maxWireSeparationForNeighbours_ ) {
	     timeSeparation = fabs(refTimeAvg-neighTimeAvg);
	     oseed          = seed;
	     addedSeed      = true;
	   }
	 }
       } 
       else { //THE POTENTIAL NEIGHBORING SEEDS ARE CLOSE IN TIME, SO NEED TO DISTINGUISH
	 auto diff = direction_ == "forward" ? refWireAvg-neighWireAvg : neighWireAvg-refWireAvg;
	 //if( candidateSeedsVeryCloseInTime )
	 //	 wireSeparation = longTrackFinderParams_.get<int>("candidateSeedsMaxWire");
	 if( diff <= wireSeparation ) {
	   
	   if( fabs(neighWireAvg - refWireAvg) < maxWireSeparationForNeighbours_) {
	     mf::LogDebug(name_) << "\t\t adding seed" << "\n";
	     wireSeparation = diff;
	     oseed          = seed;
	     addedSeed      = true;
	     if( candidateSeedsVeryCloseInTime ){
	       //need to clear first so that the seed size is only 1
	       oseeds.clear();
	       oseeds.push_back( oseed );
	     }
	   }
	 }
       }
       
     } // end loop over seeds
     
     if( addedSeed ) {
       if( oseeds.size() == 0 )
         oseeds.push_back( oseed );
       else if( oseeds.size() > 1 ) 
         oseeds.clear();
     } // end of added seeds

   } // end of selecting seeds  
   
   return;
}

 

//-----------------------------------------------------------
// check if the two seeds are neighbors
//-----------------------------------------------------------
bool LongTrackFindingUtils::neighboringSeeds( art::Ptr<StrawSeedArtRecord>& seed1, art::Ptr<StrawSeedArtRecord>& seed2 )
{
   int wire1 = seed1.get()->strawDigits.front().get()->wireID.getWire();
   int wire2 = seed2.get()->strawDigits.front().get()->wireID.getWire();

   if( fabs(wire1-wire2) < maxWireSeparation_ ) 
     return true;

   return false;
}



//------------------------------------------------------------------------------
// check if the seeds associated with track1 and track2 overlap
//------------------------------------------------------------------------------
bool LongTrackFindingUtils::trackSeedsOverlaps( StrawSeedPtrCollection track1, StrawSeedPtrCollection track2 )
{
   for(auto iseed : track1) {
     for(auto jseed : track2) {
       if( jseed == iseed ) return true;
     }
   }

   return false;
}


//------------------------------------------------------------------------------
// check if the track candidates are neighbored
//------------------------------------------------------------------------------
bool LongTrackFindingUtils::trackCandidatesNeighbors( TrackCandidateArtRecord& refCandidate, TrackCandidateArtRecord& candidate )
{
   auto seeds = candidate.strawSeeds;
   std::sort(seeds.begin(),seeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter());
   auto upstreamSeed   = seeds.front();
   auto downstreamSeed = seeds.back();

   bool mergeBrokenTracks = true;

   bool isUpstreamNeighbor   = trackCandidateSeedNeighbors(refCandidate,upstreamSeed,mergeBrokenTracks);
   bool isDownstreamNeighbor = trackCandidateSeedNeighbors(refCandidate,downstreamSeed,mergeBrokenTracks);

   auto refSeeds = refCandidate.strawSeeds;
   std::sort(refSeeds.begin(),refSeeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter());

   auto refFrontWire = refSeeds.front().get()->strawDigits.front().get()->wireID.getWire();
   auto refBackWire  = refSeeds.back().get()->strawDigits.back().get()->wireID.getWire();

   auto refFrontMod = refSeeds.front().get()->strawDigits.front().get()->wireID.getModule();
   auto refBackMod  = refSeeds.back().get()->strawDigits.back().get()->wireID.getModule();

   auto canFrontMod = upstreamSeed.get()->strawDigits.front().get()->wireID.getModule();
   auto canBackMod = downstreamSeed.get()->strawDigits.front().get()->wireID.getModule();

   //ADD IN SANITY CHECK TO MAKE SURE WE DON'T ADD 2 PERFECTLY GOOD TRACKS.
   // IN ORDER TO BE A BROKEN TRACK THERE MUST NOT BE MORE MODULES IN TOTAL THAN 8
   int nModulesRef = abs(refBackMod - refFrontMod);
   int nModulesCan = abs(canBackMod - canFrontMod);
   //mf::LogDebug(name_) << "\t\t NEIGH CHECK, n ref mods: " << nModulesRef << " n can mods: " << nModulesCan << "\n";
   if ( (nModulesRef + nModulesCan) > 8) return false;

   if( isUpstreamNeighbor ) {
     auto updown1 = fabs(refBackWire-upstreamSeed.get()->strawDigits.front().get()->wireID.getWire());
     auto updown2 = fabs(refBackWire-downstreamSeed.get()->strawDigits.back().get()->wireID.getWire());

     if( updown1 <= updown2 )
       isUpstreamNeighbor = true;
     else 
       isUpstreamNeighbor = false;

   } else if( isDownstreamNeighbor ) {
     auto updown1 = fabs(refFrontWire-upstreamSeed.get()->strawDigits.front().get()->wireID.getWire());
     auto updown2 = fabs(refFrontWire-downstreamSeed.get()->strawDigits.back().get()->wireID.getWire());

     if( updown2 <= updown1 )
       isDownstreamNeighbor = true;
     else 
       isDownstreamNeighbor = false;
   }

   auto sizeRefCandidate    = 0.;
   auto avgTimeRefCandidate = 0.;
   for(auto seed : refSeeds) {
     std::for_each(seed.get()->strawDigits.begin(),
                   seed.get()->strawDigits.end(),
                   [&](auto dig){ avgTimeRefCandidate += dig.get()->calTime; sizeRefCandidate += 1; });
   }
   avgTimeRefCandidate = avgTimeRefCandidate/double(sizeRefCandidate);

   auto sizeNeighCandidate    = 0.;
   auto avgTimeNeighCandidate = 0.;
   for(auto seed : seeds) {
     std::for_each(seed.get()->strawDigits.begin(),
                   seed.get()->strawDigits.end(),
                   [&](auto dig){ avgTimeNeighCandidate += dig.get()->calTime; sizeNeighCandidate += 1; });
   }
   avgTimeNeighCandidate = avgTimeNeighCandidate/double(sizeNeighCandidate);

   if( fabs(avgTimeRefCandidate-avgTimeNeighCandidate) <= candidateTracksMaxTime_) { 
     if( isUpstreamNeighbor || isDownstreamNeighbor )
       return true;
   }

   return false;
}


//------------------------------------------------------------------------------
// check if there is a gap between the neighboring track candidates
//------------------------------------------------------------------------------
bool LongTrackFindingUtils::isGapBetweenTrackCandidates( TrackCandidateArtRecord& refCandidate, TrackCandidateArtRecord& candidate, StrawSeedPtrCollection& unusedSeeds )
{
   if( unusedSeeds.empty() )
     return true;

   int refUpstreamModule   = refCandidate.strawSeeds.front().get()->frontModule;
   int refDownstreamModule = refCandidate.strawSeeds.back().get()->backModule;

   int neighUpstreamModule   = candidate.strawSeeds.front().get()->frontModule;
   int neighDownstreamModule = candidate.strawSeeds.back().get()->backModule;

   for(auto& seed : unusedSeeds) {
     int seedModule = seed.get()->frontModule;
     if( fabs(refDownstreamModule-seedModule) == 1 && fabs(neighUpstreamModule-seedModule) == 1 )
       return false;
     else if( fabs(refUpstreamModule-seedModule) == 1 && fabs(neighDownstreamModule-seedModule) == 1 )
       return false;
   }

   return true;
}

 

//------------------------------------------------------------------------------------
// check if a seed neighbors the track candidate at the front or back
//------------------------------------------------------------------------------------
bool LongTrackFindingUtils::trackCandidateSeedNeighbors( TrackCandidateArtRecord& track, art::Ptr<StrawSeedArtRecord>& seed, bool mergeBrokenTracks )
{
   auto trackSeeds = track.strawSeeds;
   std::sort(trackSeeds.begin(),trackSeeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter());

   auto seedfront  = seed.get()->frontModule;
   auto seedback   = seed.get()->backModule;

   auto trackfront = trackSeeds.front().get()->frontModule;
   auto trackback  = trackSeeds.back().get()->backModule;

   auto seedfrontWire  = seed.get()->strawDigits.front().get()->wireID.getWire();
   auto seedbackWire   = seed.get()->strawDigits.back().get()->wireID.getWire();

   auto trackbackWire  = trackSeeds.back().get()->strawDigits.back().get()->wireID.getWire();
   auto trackfrontWire = trackSeeds.front().get()->strawDigits.front().get()->wireID.getWire(); 
 
   auto seedfrontTime  = seed.get()->strawDigits.front().get()->calTime;
   auto seedbackTime   = seed.get()->strawDigits.back().get()->calTime;

   auto trackbackTime  = trackSeeds.back().get()->strawDigits.back().get()->calTime;
   auto trackfrontTime = trackSeeds.front().get()->strawDigits.front().get()->calTime; 

   //ADDED CHECK IF module of seed is exactly 1 after end of track NEIGHBOURING
   if ( (seedfront - trackback) != 1 && (seedback - trackback ) != 1 ) return false;

   if( seedfront-trackback == 1 && fabs(seedfrontTime-trackbackTime) <= maxTimeSeparation_ ) {
     if( trackbackWire-seedfrontWire >= 0 && trackbackWire-seedfrontWire <= maxWireSeparation_ ) 
       return true;
   }
   else if( trackfront-seedback == 1 && fabs(seedbackTime-trackfrontTime) <= maxTimeSeparation_ ) {
     if( trackfrontWire-seedbackWire <= 0 && fabs(trackfrontWire-seedbackWire) <= maxWireSeparation_ )
       return true;
   }

   bool overlaps = trackCandidateSeedShareCluster(track,seed);
   if( mergeBrokenTracks && !overlaps ) {

     //TODO in future use the fitter to determine if the seed should be added to the track candidate
     if( fabs(trackbackWire-seedfrontWire) > 2*maxWireSeparation_ )
       return false; 
     if( fabs(trackfrontWire-seedbackWire) > 2*maxWireSeparation_ )
       return false;

     if( seedfront-trackback == 2 )
       return true;
     else if( trackfront-seedback == 2 )
       return true;
   } 

   return false;
}


//------------------------------------------------------------------------------------
// check if a seed overlaps a track candidate 
//------------------------------------------------------------------------------------
bool LongTrackFindingUtils::trackCandidateSeedOverlaps( TrackCandidateArtRecord& track, art::Ptr<StrawSeedArtRecord>& seed )
{
   auto trackSeeds = track.strawSeeds;
   std::sort(trackSeeds.begin(),trackSeeds.end(),StrawTrackerSort::StrawSeedPtrModuleSorter());

   art::Ptr<StrawSeedArtRecord> tseed;
   if( trackSeeds.front().get()->frontModule == seed.get()->frontModule )
     tseed = trackSeeds.front();
   else if( trackSeeds.back().get()->frontModule == seed.get()->frontModule )
     tseed = trackSeeds.back();
   else return false;

   auto count = std::count_if(trackSeeds.begin(),trackSeeds.end(),[&](auto& iseed){ return iseed.get()->frontModule == tseed.get()->frontModule; });
   if( count != 1 )
     return false;

   auto seedWithMinDigits = tseed.get()->strawDigits.size() <= seed.get()->strawDigits.size() ? tseed : seed;
   auto seedWithMaxDigits = tseed.get()->strawDigits.size() >  seed.get()->strawDigits.size() ? tseed : seed;

   unsigned int matched = 0;
   for(auto& idigit : seedWithMinDigits.get()->strawDigits) {
     auto imod  = idigit.get()->wireID.getModule();
     auto iview = idigit.get()->wireID.getView();
     auto ilay  = idigit.get()->wireID.getLayer();
     auto iwire = idigit.get()->wireID.getWire();

     for(auto& jdigit : seedWithMaxDigits.get()->strawDigits) {
       auto jmod  = jdigit.get()->wireID.getModule();
       auto jview = jdigit.get()->wireID.getView();
       auto jlay  = jdigit.get()->wireID.getLayer();
       auto jwire = jdigit.get()->wireID.getWire();

       if( imod == jmod && iview == jview && ilay == jlay && fabs(iwire-jwire) <= 1 ) {
         matched += 1;
         break;
        }
     } // end loop over seed with maximum number of digits
   } // end loop over seed with minimum number of digits


   if( matched == 0 )
     return false;

   auto avgTimeDiff = getTimeDifferenceBetweenSeeds(trackSeeds);
 
   auto avgTrackTime    = 0.;
   auto avgTrackTimeDen = 0.;
   for(auto& trackSeed : trackSeeds) {
     std::for_each(trackSeed.get()->strawDigits.begin(),
                   trackSeed.get()->strawDigits.end(),
                   [&](auto& dig){ avgTrackTime += dig.get()->calTime; });
     avgTrackTimeDen += double(trackSeed.get()->strawDigits.size());
   }
   avgTrackTime = avgTrackTime/avgTrackTimeDen;

   auto avgSeedTime = 0.;
   std::for_each(seed.get()->strawDigits.begin(),
                 seed.get()->strawDigits.end(),
                 [&](auto& dig){ avgSeedTime += dig.get()->calTime; });
   avgSeedTime = avgSeedTime/double(seed.get()->strawDigits.size());
   
   if( fabs(avgTrackTime-avgSeedTime) < (avgTimeDiff+ avgSeedRMSTime_) )
     return true;

   return false;
}
 

//------------------------------------------------------------------------------------------ 
// check if a seed and the track candidate share a cluster
//------------------------------------------------------------------------------------------ 
bool LongTrackFindingUtils::trackCandidateSeedShareCluster( TrackCandidateArtRecord& track, art::Ptr<StrawSeedArtRecord>& seed )
{
  bool shareCluster = false;
  for(auto icluster : seed.get()->strawClusters ) {
    for(auto trackSeed : track.strawSeeds) {
      for(auto jcluster : trackSeed.get()->strawClusters) {
        if( *jcluster.get() == *icluster.get() )
          shareCluster = true;
      }
    }
  }
  return shareCluster;
}

//------------------------------------------------------------------------------------------ 
// determine if the digits on the track candidate are reasonable
//------------------------------------------------------------------------------------------ 
bool LongTrackFindingUtils::checkDigitsOnTrackCandidate( TrackCandidateArtRecord& track )
{
   bool pass = true;

   StrawDigitPtrCollection digits;
   for(auto& seed : track.strawSeeds) {
     for(auto& dig : seed.get()->strawDigits) {
       if( std::find(digits.begin(),digits.end(),dig) == digits.end() )
         digits.push_back( dig );
     }
   }

   auto refDigit = *digits.begin();
   auto refWire  = refDigit.get()->wireID.getWire();
   auto idigit   = digits.erase( digits.begin() );

   while( idigit != digits.end() ) {
     
      auto curDigit = *idigit;
      auto curWire  = curDigit.get()->wireID.getWire();

      if( fabs(curWire-refWire) < maxWireSeparation_ ) {
        refDigit = curDigit;
        refWire  = curWire;
        idigit   = digits.erase( idigit );

      } else {
        pass = false;
        break;
      }
   }

   return pass;
}

//----------------------------------------------------------------------------------------
// determine if the seeds are overlapping the same module
//----------------------------------------------------------------------------------------
bool LongTrackFindingUtils::overlappingSeeds( StrawSeedPtrCollection& seeds )
{
   bool overlap = false;
   for(auto& iseed : seeds) {
     auto iclusters = iseed.get()->strawClusters;

     for(auto& jseed : seeds) {
       if( jseed == iseed ) continue;
       auto jclusters = jseed.get()->strawClusters;

       for(auto jcluster : jclusters) {
         if( std::find(iclusters.begin(),iclusters.end(),jcluster) != iclusters.end() ) 
           overlap = true;
       }
     }
   }

   mf::LogDebug(name_) << "\t\tAre the seeds overlapping? " << overlap << "\n";
   return overlap;
}


//----------------------------------------------------------------------------------------
// determine if the reference and neighboring seeds overlap in wire regions
//----------------------------------------------------------------------------------------
bool LongTrackFindingUtils::overlappingSeeds( art::Ptr<StrawSeedArtRecord>& refSeed, art::Ptr<StrawSeedArtRecord>& neighSeed )
{
   bool overlap = false;
   
   if( fabs(refSeed.get()->frontModule-neighSeed.get()->frontModule) != 1 )
     return false;
   if( fabs(refSeed.get()->backModule-neighSeed.get()->backModule)   != 1 )
     return false;

   for(auto& refDigit : refSeed.get()->strawDigits) {
     auto found = std::find_if(neighSeed.get()->strawDigits.begin(),
                               neighSeed.get()->strawDigits.end(),
                               [&](auto& neighDigit)
                               { return fabs(neighDigit.get()->wireID.getWire()-refDigit.get()->wireID.getWire()) == 0; });

     if( found != neighSeed.get()->strawDigits.end() )
       overlap = true;
   }

   return overlap;
}


//----------------------------------------------------------------------------------------
// remove seeds that are already living on a track candidate
//----------------------------------------------------------------------------------------
void LongTrackFindingUtils::checkSeedIsOnTrackCandidate( StrawSeedPtrCollection& seeds, TrackCandidateArtRecordCollection& candidates )
{
   auto iseed = seeds.begin();
   while( iseed != seeds.end() ) {
     bool removeSeed = false;

     for(auto candidate : candidates) {
       auto found = std::find(candidate.strawSeeds.begin(),candidate.strawSeeds.end(),*iseed);
       if( found != candidate.strawSeeds.end() )
         removeSeed = true;
     }

     if( removeSeed )
       iseed = seeds.erase(iseed);
     else 
       ++iseed;
   }

   return;
}
 

//----------------------------------------------------------------------------------------
// get the averge time difference between the seeds on the track candidate
//----------------------------------------------------------------------------------------
double LongTrackFindingUtils::getTimeDifferenceBetweenSeeds( StrawSeedPtrCollection& seeds )
{
   double timeDiff    = 0.; 
   double timeDiffDen = 0.;

   if( seeds.size() < 2 )
     return timeDiff;

   auto itime = 0.;
   auto iseed = seeds.at(0);
   std::for_each(iseed.get()->strawDigits.begin(),
                 iseed.get()->strawDigits.end(),
                 [&](auto& idig){ itime += idig.get()->calTime; });
   itime /= double(iseed.get()->strawDigits.size());

   for(unsigned int j = 1; j < seeds.size(); ++j) {
     auto jtime = 0.;
     auto jseed = seeds[j];
     std::for_each(jseed.get()->strawDigits.begin(),
                   jseed.get()->strawDigits.end(),
                   [&](auto& idig){ jtime += idig.get()->calTime; });
     jtime /= double(jseed.get()->strawDigits.size());


     timeDiff    += fabs(itime-jtime);
     timeDiffDen += 1.;
     itime        = jtime;
   }  

   timeDiff = timeDiff/timeDiffDen;
   return timeDiff;
}


//----------------------------------------------------------------------------------------
// call the least squares minimizer fitter
//----------------------------------------------------------------------------------------
void LongTrackFindingUtils::callLSQFitter( TrackCandidateArtRecord& candidate )
{
   std::vector<gm2geom::CoordSystem3Vector> measurements;
   std::vector<gm2geom::CoordSystem3Vector> errors;
   gm2strawtracker::Helix helix;

   for(auto seed : candidate.strawSeeds) {
     for(auto node : seed.get()->nodes) {
       measurements.push_back( node.position.transform(cs_,"world") );
       errors.push_back( node.posError );
     }
   }

   bool fitConverges  = LSQFitter_.makeFittedTrack(measurements,errors,helix); 
   bool keepCandidate = requiredFitConverges_ && !fitConverges ? false : true;

   mf::LogDebug(name_) << "\tdoes the fit converges= " << fitConverges 
                          << ", chi2= " << candidate.chi2 
                          << ", keep candidate= " << keepCandidate << "\n";

   candidate.helix    = helix;
   candidate.helixFit = true;
   candidate.type     = TrackPatRecognition::LongPatRec;
}

//TWAT
TH1F* LongTrackFindingUtils::strawPlot(StrawDigitPtrCollection digits, TH1F* bgPlot, int fillColor ){
  
  //bgPlot gets drawn over

  //auto eventFolderName = Form("event_%d", counter_);
  TH1F* hdigits = (TH1F*) bgPlot->Clone();
  hdigits->SetName("twatname");

  // get digits from candidate
  //for (auto& digit : cand.get()->strawDigits ){
  for (auto& digit : digits ){
    mf::LogDebug(name_) << "\t\t" << digit->wireID << "\n";
    
    //std::string coordName = Form("Module%d:%d", digit->wireID.getStation(), digit->wireID.getModule());
    //gm2geom::CoordSystemsStoreData newcs_ = csMap_.find(coordName)->second;

    //gm2geom::CoordSystem3Vector centreInWorld = digit->wireID.getCentreInWorld(newcs_);    
    gm2geom::CoordSystem3Vector centreInWorld = wireCentreInWorld_.at(geom_.getGlobalWire(digit->wireID));
    TEllipse* hitEllipse = new TEllipse(centreInWorld.z(), centreInWorld.x(),2.5,2.5);
    hitEllipse->SetFillColor(fillColor);
    hitEllipse->SetLineColor(fillColor);
    hdigits->GetListOfFunctions()->Add(hitEllipse->Clone());
    delete hitEllipse;
  }
  
  return (TH1F*) hdigits->Clone();
}



//----------------------------------------------------------------------------------------
// call the simple circle fast fitter
//----------------------------------------------------------------------------------------
void LongTrackFindingUtils::callSimpleCircleFitter( TrackCandidateArtRecord& candidate )
{
   StrawDigitPtrCollection digits;

   if( candidate.strawDigits.empty() ) {
     std::for_each(candidate.strawSeeds.begin(),
                   candidate.strawSeeds.end(),
                   [&](auto& seed)
                   { digits.insert(digits.end(),seed.get()->strawDigits.begin(),seed.get()->strawDigits.end()); });
   } else {
     digits = candidate.strawDigits;
   }
              
   std::sort(digits.begin(),digits.end(),StrawTrackerSort::StrawDigitArtPtrRowSorter());

   std::vector< gm2geom::CoordSystem3Vector> positions;
   std::vector< gm2geom::CoordSystem3Vector> errors;

   auto dummyValue = -9999;  //TODO use NAN
   gm2strawtracker::GeaneState geane;

   WireID wire0(digits.front()->wireID.getStation(),
                digits.front()->wireID.getModule(),
                gm2strawtracker::u_view , 0, 16);

   

   //old std::string coordName = Form("Module%d:%d", wire0.getStation(), wire0.getModule());
   //old gm2geom::CoordSystemsStoreData* cssStraw             = &csMap_.find(coordName)->second;
   //old auto frontPosition        = wire0.getCentreInStation(*cssStraw);
   auto frontPosition = wireCentreInStation_.at(geom_.getGlobalWire(wire0));
   auto geaneFrame           = frontPosition.getCoordSystemName(); //same as: Form("TrackerStation[%d]", digits.front()->wireID.getStation());
   
   geane.XYZPosition = gm2geom::CoordSystem3Vector(0, 0, frontPosition.z()-10, geaneFrame); // 10 is a hardcoded number for now, matching the placement of dummy planes in TrackerDummyPlane_service.cc

   //int currentMod = wire0.getModule();
   for(auto& digit : digits) {
     
     ////speed increase: all digits are in the same station, so the cssMap is only refreshed if the module changes
     //if (digit->wireID.getModule() != currentMod){
     //  coordName = Form("Module%d:%d", digit->wireID.getStation(), digit->wireID.getModule());
     //  cssStraw = &csMap_.find(coordName)->second;
     //  currentMod = digit->wireID.getModule();
     //}
     
     //old auto trackerPosition      = digit.get()->wireID.getCentreInStation(*cssStraw);
     auto trackerPosition = wireCentreInStation_.at(geom_.getGlobalWire(digit.get()->wireID));

     gm2geom::CoordSystem3Vector position(0.,0.,0.,geaneFrame);

    if( digit.get()->wireID.getView() == gm2strawtracker::u_view ) {
       position.setX( trackerPosition.z() );
       position.setY( dummyValue );
       position.setZ( trackerPosition.x() );

     } else if( digit.get()->wireID.getView() == gm2strawtracker::v_view ) {
      position.setX( dummyValue );
      position.setY( trackerPosition.z() );
       position.setZ( trackerPosition.x() );
    } 
    
    positions.push_back( position );
   } // end of getting the position
   
   
   CircleFitter_.makeFittedTrack(positions,errors,geane);
   candidate.geane = geane;
   candidate.type  = TrackPatRecognition::LongPatRec; 

   mf::LogDebug(name_) << "\t\tGeane: " << geane << "\n";

   return; 
}



