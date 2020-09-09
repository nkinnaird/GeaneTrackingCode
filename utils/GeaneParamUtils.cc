#include "GeaneParamUtils.hh"

namespace gm2strawtracker {

GeaneParamUtils::GeaneParamUtils(fhicl::ParameterSet const & p) 
  : dummyUtils_(p)
  , name_( "GeaneParamUtils" )
  , fitMode_(p.get<string>("fitMode","badFitMode"))
  , DummyModuleLabel_( p.get<std::string>("DummyModuleLabel",dataModuleDefs::strawBuilderModuleLabel()) )
  , DummyInstanceName_( p.get<std::string>("DummyInstanceName","trackerdummyplane") )
  , theTarget(0)
  , sgeom_()
  , noHit(900000.)
  , G4EVERBOSE_(p.get<int>("G4EVERBOSE",0))
  , trackingVerbose_(p.get<int>("trackingVerbose",5))
  , onlyPrimaryPositrons_(p.get<bool>("onlyPrimaryPositrons"))
  , geaneTrackUtils_()
  , useCircleGuess_(p.get<bool>("useCircleGuess",false))
  , skipLayers_(p.get<vector<unsigned int> >("skipLayers",{}))
{
   mf::LogDebug(name_) << "Begin GeaneParamUtils Constructor" << "\n";

   // Make sure to perform all of this only once before iterating through all tracks.
   art::ServiceHandle<artg4::DetectorHolderService> dh;
   mf::LogDebug(name_) << std::endl << "Created Detector Holder Service" << std::endl;

   if(!dh->isInitialized()){
      dh->initialize();
      mf::LogDebug(name_) << std::endl << "Initialized Detector Holder Service" << std::endl;

      dh->constructAllLVs();
      mf::LogDebug(name_) << std::endl << "Constructed LVs" << std::endl;
   }

   G4VPhysicalVolume* world = dh->worldPhysicalVolume();
   mf::LogDebug(name_) << std::endl << "Grabbed World Physical Volume" << std::endl;

   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////

   // Set global transportation field manager directly from new gm2FieldManager_service (before was taking fields directly from arc volume)
   art::ServiceHandle<gm2geom::gm2FieldManager> fMgr;
   G4TransportationManager::GetTransportationManager()->SetFieldManager(fMgr->GetUnifiedFieldManager()); 
   G4TransportationManager::GetTransportationManager()->GetFieldManager()->SetChordFinder(fMgr->GetUnifiedFieldManager()->GetChordFinder());

   // Initialize the GEANT4e manager 
   geaneTrackUtils_.g4emgr  = gm2GeanePropagatorManager::GetErrorPropagatorManager();
   g4edata = G4ErrorPropagatorData::GetErrorPropagatorData();

   // set error_propagation verbose level
   g4edata->SetVerbose(G4EVERBOSE_);
   mf::LogDebug(name_) << " G4 verbose level is : " <<  g4edata->verbose() << std::endl;
   mf::LogDebug(name_) << std::endl << "Initializing GEANE and setting WORLD" << world << std::endl;

   geaneTrackUtils_.g4emgr->SetUserInitialization(world);
   geaneTrackUtils_.g4emgr->InitGeant4e();

   mf::LogDebug(name_) << std::endl << "GEANE Initialized" << std::endl;

   // set the tracking verbose level
   std::ostringstream oss;
   oss << "/tracking/verbose " << trackingVerbose_;
   std::string var = oss.str();

   G4UImanager::GetUIpointer()->ApplyCommand(var);

   /*
    * change step length with this - can make steerable - maybe don't want to, and just let geant take care of it
    *
   G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength .01 mm"); 
    */

   wireDiamError = 2.*sgeom_.outerRadiusOfTheGas/sqrt(12);

   mf::LogDebug(name_) << "End GeaneParamUtils Constructor" << "\n";
}

// set up the GEANE parameters
int GeaneParamUtils::setupParams(art::Event & e, gm2strawtracker::TrackDetailArtRecord& trackFitDetails, bool firstTrackInEvent)
{
  mf::LogInfo info(name_); 

  art::Handle<gm2truth::GhostDetectorArtRecordCollection> DummyDataHandle;
  bool foundDummycollection = e.getByLabel(DummyModuleLabel_,DummyInstanceName_,DummyDataHandle);
  if( !foundDummycollection ) {
    mf::LogWarning(name_) << "No Dummy collection in this event (\"" << DummyModuleLabel_ << "\":\"" << DummyInstanceName_ << "\")\n";
  } else if(firstTrackInEvent){
    dummyUtils_.fillDummyHits(DummyDataHandle); // pass dummy hits to utils file
  }

  if(!foundDummycollection && fitMode_ == "truthLRFit"){
    throw cet::exception(name_) << "No dummy collection but trying to do a fit involving truth!\n";
    return -2;
  }

  // Initialization of vectors, arrays, matrices, etc.
  std::vector <std::vector<double> > wireParamMeasured;
  wireParamMeasured.resize(5, std::vector<double>(geaneTrackUtils_.maxNumPlanes, noHit)); // Initialize 2D vector with 5 rows and maxNumPlanes columns filled with noHit values.

  std::vector<int> planesHitInThisEvent;
  std::set<int> planesHitSet; // For checking to make sure there is only 1 digit per plane, otherwise cancel the fitting.

  std::vector<double> dcaMeasured(geaneTrackUtils_.maxNumPlanes, noHit);

  std::vector<double> strawDiamErrors(geaneTrackUtils_.maxNumPlanes, wireDiamError);
  std::vector<double> dcaErrors(geaneTrackUtils_.maxNumPlanes, noHit);

  std::vector<gm2strawtracker::GeaneHitSide> inputHitSides(geaneTrackUtils_.maxNumPlanes, gNA_side); // hit sides which are filled in different ways

  // get station number for track
  stationNumber = trackFitDetails.strawDCADigits.at(0)->digit->wireID.getStation();
  stringstream stationStream;
  stationStream << "TrackerStation[" << stationNumber << "]";
  stationStr = stationStream.str();

  dummyUtils_.fillStationStr(stationStr, stationNumber); // pass station number to dummy utils

  int firstPlaneHit = 99; // Which is not the 0 plane
  int numPlanesHit  = 0;
  int lastPlaneHit  = 0;

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  gm2strawtracker::StrawDCADigitPtrCollection hitDCADigits; // temp container

  for (auto& dcahit : trackFitDetails.strawDCADigits) // loop through digits before sorting them
  {
    auto & hit = dcahit->digit;

    if(!e.isRealData() && onlyPrimaryPositrons_){
      for (int i = 0; i < int(hit->strawMCDigit.strawMCHits.size()); ++i)
      {
        if(hit->strawMCDigit.strawMCHits[i]->parent_ID != 1 || hit->strawMCDigit.strawMCHits[i]->pdg_ID != -11) 
        {
          mf::LogWarning(name_) << "\n" << "MCHits include non-primary positron - drop event.\n";
          return 11; // if hits include non-primary positron then drop it
        }
      }
    }

    // Skip requested layers (don't include in fit)
    if(std::find(skipLayers_.begin(), skipLayers_.end(), sgeom_.getGlobalLayer(hit->wireID)) != skipLayers_.end()) continue;
    hitDCADigits.push_back(dcahit);
  }


  // the way the dummy hits are recorded, these digits need to be ordered
  // perhaps there is some other order dependence as well, but I removed what I could find
  // if I sort my own hitDCADigits container, that ends up sorting the digits in the TrackCandidateArtRecord through the pointers as well

  std::sort(trackFitDetails.strawDCADigits.begin(), trackFitDetails.strawDCADigits.end(), StrawTrackerSort::StrawDCADigitArtPtrRowSorter());
  std::sort(hitDCADigits.begin(), hitDCADigits.end(), StrawTrackerSort::StrawDCADigitArtPtrRowSorter()); // sort the digits by planeNum

  if (hitDCADigits.size() < 5)  {
    mf::LogWarning(name_) << "Track candidates digits size is too small to fit\n";
    return 20;
  }

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////


  // Loop over associated straw dca digits
  for(unsigned int i_digit = 0; i_digit < hitDCADigits.size(); i_digit++){

    // copy
    auto & dcaDigit = hitDCADigits.at(i_digit);

    // Get digit ptr from this dca digit
    auto& hit = dcaDigit->digit;
    numPlanesHit++;
    int moduleNum = hit->wireID.getModule();
    int subPlaneNum = hit->wireID.getLayer() + 2 * hit->wireID.getView(); // layer number 0 to 1 in each layer, view number = 0 for u and 1 for v
    int planeNum = 1+(moduleNum*4)+subPlaneNum;
   
    planesHitInThisEvent.push_back(planeNum);
    planesHitSet.insert(planeNum); // for checking same plane isn't hit twice
   
    if (i_digit == 0)
    { 
      firstPlaneHit = planeNum; // update just once at the start
      if (foundDummycollection){
         dummyUtils_.clearDummyHolders(); // Clear dummy plane holders before re-filling for this track
         dummyUtils_.fillDigitHitInfo(&hit);
         dummyUtils_.fillTruthParams(firstPlaneHit, true); // call once for 0th plane params
      } 
    }
    lastPlaneHit = planeNum; // update each time
   
    mf::LogTrace(name_) << "Particle hit module: " << moduleNum << " subplane number: " << subPlaneNum << " total plane number: " << planeNum << "\n";
    if(!e.isRealData()) mf::LogTrace(name_) << " trackID " << hit->strawMCDigit.strawMCHits[0]->trackID << " parentID: " << hit->strawMCDigit.strawMCHits[0]->parent_ID << "\n";

    // fill measured vectors with hit info
    auto cssTrackerWorld = detCoordMap_.find("TrackerStation")->second;
    auto cssStraw        = detCoordMap_.find(Form("Module%d:%d", hit->wireID.getStation(), hit->wireID.getModule()))->second;
    auto wireGeaneTrackerPosition = hit->wireID.getCentreInWorld(cssStraw).transform(cssTrackerWorld,stationStr);
    auto wireGeanePosition        = geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*wireGeaneTrackerPosition.x(); // y coordinate is 0 so omit

    wireParamMeasured[3][planeNum] = wireGeanePosition; // only one of the U or V values will be used, U and V are the same for the wires with y=0
    wireParamMeasured[4][planeNum] = wireGeanePosition;

    dcaMeasured[planeNum] = dcaDigit->dca; // Fill dca values here for later use
    dcaErrors[planeNum]   = dcaDigit->dcaError;

    inputHitSides.at(planeNum) = gCenter; // default sides to the centers at the start (non hit layers default to na side)

    if(foundDummycollection){
      dummyUtils_.fillDigitHitInfo(&hit);
      dummyUtils_.fillTruthParams(firstPlaneHit, false);
    }

  } // hit digits iterator loop 

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  if (planesHitSet.size() != planesHitInThisEvent.size())
  {
    info << "Individual plane(s) hit more than once in digits - skipping event." << "\n";
    return 4;
  }

  if(foundDummycollection){
    if(dummyUtils_.hitDummyPlaneAlready)
    {
      info << "Single dummy plane was hit more than once for some reason - break out and fail." << "\n";
      return 3;
    }
  }

  mf::LogTrace(name_) << "Number of planes hit in this event: " << numPlanesHit << "\n"; 
  mf::LogTrace(name_) << " first plane hit " << firstPlaneHit << "\n";

  mf::LogTrace(name_) << " Filled into parameter measured and extra arrays: " << "\n";
  for (int planeNum = 0; planeNum < geaneTrackUtils_.maxNumPlanes; ++planeNum)
  {
    mf::LogTrace(name_) << "planeNum: " << planeNum << "\n"
                        << "wire u: "   << wireParamMeasured[3][planeNum] << "\n"
                        << "wire v: "   << wireParamMeasured[4][planeNum] << "\n" 
                        << "dca: "      << dcaMeasured[planeNum] << "\n" 
                        << "dca error: " << dcaErrors[planeNum] << "\n"
                        << "\n";
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // Set the initial starting trajectory.
  // These units are in MeV mm.
  // Can take the starting plane values from truth (if it exists), and add changes to similulate a non-perfect starting trajectory from some initial weak trackfitting algorithm - or can take guess from circle fitter.
  // 
  // TODO : add another condition for track refinement

  std::vector<double> startingParameters;
  if( !useCircleGuess_ && foundDummycollection){
    dummyUtils_.createStartGuess(&startingParameters);
  }
  else if ( useCircleGuess_ )
  {
    startingParameters.push_back(trackFitDetails.candidate->geane.XYZPosition.x());
    startingParameters.push_back(trackFitDetails.candidate->geane.XYZPosition.y());
    startingParameters.push_back(trackFitDetails.candidate->geane.XYZPosition.z());
    startingParameters.push_back(trackFitDetails.candidate->geane.momentum.x());
    startingParameters.push_back(trackFitDetails.candidate->geane.momentum.y()); 
    startingParameters.push_back(trackFitDetails.candidate->geane.momentum.z());
  }
  else{
    throw cet::exception(name_) << "Something wrong with starting guess fcl configuration for initial guess. \n";    
  }


  for (auto i : startingParameters)
  {
    if(std::isnan(i)){
      info << "\n" << "Circle fit returned nan value for a starting parameter. \n";
      return -1;
    } 
  }

  // Fill track starting parameters 
  trackFitDetails.geaneHits.startingGeaneParameters = startingParameters; 

  /////////////////////////////////////////////////////////////////////////////////////
  // Fill track variables that haven't been filled yet.

  trackFitDetails.geaneHits.geaneMeasuredParameters = wireParamMeasured; // geaneMeasuredParameters is filled in angularCorrection or modifyMeasuredParams
  trackFitDetails.trackNumPlanesHit       = numPlanesHit;
  trackFitDetails.trackFirstPlaneHit      = firstPlaneHit;
  trackFitDetails.trackLastPlaneHit       = lastPlaneHit;
  trackFitDetails.trackPlanesHitList      = planesHitInThisEvent;

  trackFitDetails.measuredDCAs    = dcaMeasured;
  trackFitDetails.geaneHits.geaneWireUVPositions = wireParamMeasured;

  trackFitDetails.UVerrors        = strawDiamErrors; // set UV errors to wire errors for default just for filling - they will be updated later
  trackFitDetails.uncorrErrors    = strawDiamErrors;
  trackFitDetails.strawDiamErrors = strawDiamErrors;
  trackFitDetails.dcaErrors       = dcaErrors;

  trackFitDetails.geaneHits.geaneHitSides = inputHitSides;
  trackFitDetails.geaneHits.inputHitSides = inputHitSides;
  /////////////////////////////////////////////////////////////////////////////////////

  if(foundDummycollection){
    dummyUtils_.getTrackDummyHits(DummyDataHandle,trackFitDetails); // dummy hits at wire planes
    dummyUtils_.fillLRFromTruth(trackFitDetails);
    if(fitMode_ == "truthLRFit"){
       trackFitDetails.geaneHits.geaneHitSides = dummyUtils_.trueHitSides;
    }
  } // if tracker dummy planes exist


  /////////////////////////////////////////////////////////////////////////////////////
  // set tracing vectors for Geane error propagation (parallel to wire planes)
  // need to define these in world coordinates

  // set tracker Y vector as vertical
  tracingWY.setX(0);
  tracingWY.setY(1);
  tracingWY.setZ(0);

  // tracker X vector as radially out along straws
  WireID straw1(stationNumber, 0, gm2strawtracker::u_view, 0, 0);
  WireID straw2(stationNumber, 0, gm2strawtracker::u_view, 0, 31);

  auto cssStraw = detCoordMap_.find(Form("Module%d:%d", stationNumber, 0))->second;
  auto firstV   = straw1.getCentreInWorld(cssStraw);
  auto secondV  = straw2.getCentreInWorld(cssStraw);
  auto thirdV   = secondV - firstV;

  tracingVX.setX(thirdV.x()/thirdV.mag());
  tracingVX.setY(0);
  tracingVX.setZ(thirdV.z()/thirdV.mag());

  auto crossV = tracingVX.cross(tracingWY);

  // set normal vector (equal to Z in tracker frame) for target planes
  surfNorm.setX(crossV.x());
  surfNorm.setY(crossV.y());
  surfNorm.setZ(crossV.z());

  return 0;
}


// get and fill the error propagation matrix
int GeaneParamUtils::errorProp(gm2strawtracker::TrackDetailArtRecord & track){ 

  mf::LogInfo info(name_);

  int ierr = 0; // geant4 error propagation failure int

  std::vector<double> wirePlaneZPositions;
  wirePlaneZPositions.resize(geaneTrackUtils_.maxNumPlanes, noHit);

  // First vector is param number (rows), second is plane number (columns)
  std::vector <std::vector<double> > trackParamPredicted; // What GEANT4E calculates as the average predicted params from propagation from tracing vectors defined as above - 1/p, px/pz, py/pz, x, y for now.
  trackParamPredicted.resize(5, std::vector<double>(geaneTrackUtils_.maxNumPlanes, noHit)); // Initialize things to noHit to start - noHit isn't actually used for checking things anymore.

  // starting parameters is defined in the geane frame (now equal to the tracker frames) and needs to be put into the world frame for propagation
  gm2geom::CoordSystem3Vector startingPosCoord(track.geaneHits.startingGeaneParameters[0], track.geaneHits.startingGeaneParameters[1], track.geaneHits.startingGeaneParameters[2], stationStr);
  gm2geom::CoordSystem3Vector startingMomCoord(track.geaneHits.startingGeaneParameters[3], track.geaneHits.startingGeaneParameters[4], track.geaneHits.startingGeaneParameters[5], stationStr);

  auto worldStartingPosCoord = startingPosCoord.transform(detCoordMap_.find("TrackerStation")->second, "world");
  auto worldStartingMomCoord = startingMomCoord.transform(detCoordMap_.find("TrackerStation")->second, "world", true);

  double startingXPos = worldStartingPosCoord.x(); // Grab starting params from trackDetails to propagate the particle.
  double startingYPos = worldStartingPosCoord.y();
  double startingZPos = worldStartingPosCoord.z();
  double startingXMom = worldStartingMomCoord.x();
  double startingYMom = worldStartingMomCoord.y();
  double startingZMom = worldStartingMomCoord.z();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<gm2GeaneMatrix> transferMatrices; // set of transfer/transport matrices, one for each plane, describes transport from one plane to the next plane. Surface system
  std::vector<gm2GeaneMatrix> transferMatricesFreeSystem; // free coordinate system
  std::vector<gm2GeaneMatrix> errorMatrices; // set of error matrices, one for each plane - surface system

  transferMatrices.clear(); // Clear vectors at the beginning of each pass so things can be push_backed correctly into the vectors - can probably be removed now that the pass loop is external, but leaving for now.
  transferMatricesFreeSystem.clear();
  errorMatrices.clear();


  G4ThreeVector xv3( startingXPos, startingYPos, startingZPos ); // Define starting params for geant4e
  G4ThreeVector pv3( startingXMom, startingYMom, startingZMom );

  wirePlaneZPositions.at(0) = track.geaneHits.startingGeaneParameters[2];

  theG4ErrorMode = G4ErrorMode_PropForwards; // Propagate forward from plane to plane.

  // Need to provide error matrix at the start of the trajectory state. Units for the error (and transport) matrices are in GeV cm.
  gm2GeaneTrajErr error( 5, 0 ); // starts empty

  /////////////////////////////////////////////////////////////////////////////////////
  // starting error notes
  //
  // When this is not included I get perfect results (for no material) but concepetually it seems like it would be better to include an error/spread in intial parameters before fitting.
  // When this is included the fit overestimates the errors and the p values rise towards one. (when error is set to covarianceTotal from the previous iteration)
  //
  // Might also want to include spread in initial parameters just before the 1st iteration, from the initial quick fitter. (Circle fitter)
  //
  /////////////////////////////////////////////////////////////////////////////////////

  // This for a free trajectory state with 1/p, lambda, phi, y perp, and z perp.
  gm2GeaneFreeTrajState* myFreeTrajState 
    = new gm2GeaneFreeTrajState("e+", xv3, pv3, error );

  std::vector<gm2GeaneMatrix> SC2SDTransformationMatrix; // Matrices for transformation to and from free system and plane/surface system. These will be filled with per-step Jacobian matrices, not per-plane matrices.
  std::vector<gm2GeaneMatrix> SC2SDTransformationMatrixInverse;

  gm2GeaneMatrix fillTransformationMatrix0(5,0); // 0th step transformation matrix.

  // This for a surface trajectory state with 1/p, v', w', v, and w. This needs to be orthogonal as currently written in the GEANT4 source code.
  // This creates the surface trajectory state from the free state.
  gm2GeaneSurfaceTrajState* mySurfaceTrajState
    = new gm2GeaneSurfaceTrajState(*myFreeTrajState, tracingVX, tracingWY, fillTransformationMatrix0); // Passing it the empty transformation matrix here fills it with that matrix for the starting plane point.

   SC2SDTransformationMatrix.push_back(fillTransformationMatrix0);

  int failure; // Int to pass and see if the inversion failed.

  // In this case, the 0 here does not refer to the 0 plane, but the starting point plane. This is necessary for the transport matrix coordinate transformation later. 
  SC2SDTransformationMatrixInverse.push_back(SC2SDTransformationMatrix[0].inverse(failure)); 

  geaneTrackUtils_.g4emgr = gm2GeanePropagatorManager::GetErrorPropagatorManager();

  ////////////////////////////////////////////////////////////////////////////////////

  gm2GeaneMatrix zeroMatrix = myFreeTrajState->GetTransfMat(); // Transport matrix for the free state when no steps have been made is a zero matrix.
  mf::LogTrace(name_) << "Check that this is an initial zero matrix: " << zeroMatrix << "\n";

  // Free state system has this particular method while the surface system does not. On plane 0 it is a 0 matrix which will be excluded in the GEANE formulation later on.
  // There is no transfer matrix to plane 0 (or the starting plane) because there is no plane before it.

  transferMatrices.push_back(zeroMatrix); // Have to fill with an object to start otherwise it seg faults. Surface transfer matrices can also just be filled with the 0 matrix.
  transferMatricesFreeSystem.push_back(zeroMatrix);

  errorMatrices.push_back(mySurfaceTrajState->GetError()); // Starts off with 0 matrix as well, but writing more naturally.

  /////////////////////////////////////////////////////////////////////////////////////

  // adjusts first and last modules included in the tracking when layers are skipped (fitting will still only consider hit layers when doing the minimization)

  int firstModuleHit = int((track.trackFirstPlaneHit-1)/4.);
  int lastModuleHit  = int((track.trackLastPlaneHit-1)/4.);

  int firstPlaneOfFirstModuleHit = firstModuleHit * 4 + 1;
  int lastPlaneOfLastModuleHit   = (lastModuleHit+1) * 4;

  if(skipLayers_.size() > 0) // loop through skipped layers, and adjust beginning module and ending module of tracking accordingly (eg. if a single hit in the first module was dropped, resulting in the tracking starting a module later than it should)
  {
    for (uint i = 0; i < skipLayers_.size(); ++i)
    {
       int stationN = int(skipLayers_.at(i)/32); // 0, 1, 2
       int skipLayerStation = -1; // 0, 12, 18

       if(stationN == 0) skipLayerStation = 0;
       else if (stationN == 1) skipLayerStation = 12;
       else if (stationN == 2) skipLayerStation = 18;
       else throw cet::exception(name_) << "Something weird with comparing skipped layer station number. \n";    

       if(skipLayerStation != stationNumber) continue; // the skipLayers vector holds layer numbers for mulitple tracking stations (0-95)

       int skippedLayer = (skipLayers_.at(i)+1)%32;
       int skippedLayerModule = int((skippedLayer-1)/4.);

       if (skippedLayerModule < firstModuleHit)
       {
         firstPlaneOfFirstModuleHit = skippedLayerModule * 4 + 1;
         firstModuleHit = skippedLayerModule;
       } 
       else if (skippedLayerModule > lastModuleHit)
       {
         lastPlaneOfLastModuleHit = (skippedLayerModule+1) * 4;
         lastModuleHit = skippedLayerModule;
       } 
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////////////
  // trace from plane to plane - traces from first module hit to last module hit, including intermediate planes where no measurements were made, and skips the outside ones

  for (int planeNum = 1; planeNum < geaneTrackUtils_.maxNumPlanes; ++planeNum)
  {
     mf::LogTrace(name_) << "\n" << "Tracking to plane: " << planeNum << "\n";

     if (planeNum < firstPlaneOfFirstModuleHit) // Skip over non-hit planes in front of the first plane of the first module hit. Starting track parameters are always defined in front of this plane.
     {
        mf::LogTrace(name_) << "The first plane hit was: " << track.trackFirstPlaneHit << " in module: " << int((track.trackFirstPlaneHit-1)/4.) << " Skipping planenum " << planeNum << " in non-hit module. \n";
        transferMatrices.push_back(zeroMatrix); // Fill matrix arrays with zero arrays for skipped planes.
        transferMatricesFreeSystem.push_back(zeroMatrix);
        errorMatrices.push_back(zeroMatrix); // Same for error matrices.

        SC2SDTransformationMatrix.push_back(zeroMatrix);
        SC2SDTransformationMatrixInverse.push_back(zeroMatrix);

        continue;
     }

     if (planeNum > lastPlaneOfLastModuleHit) // Skip over planes after the last plane of the last module hit for speed reasons.
     {
        mf::LogTrace(name_) << "The last plane hit was: " << track.trackLastPlaneHit << " in module: " << int((track.trackLastPlaneHit-1)/4.)  << " Skipping planenum " << planeNum << " in non-hit module.\n";
        transferMatrices.push_back(zeroMatrix); // Fill matrix arrays with zero arrays for skipped planes.
        transferMatricesFreeSystem.push_back(zeroMatrix);
        errorMatrices.push_back(zeroMatrix); // Same for error matrices.

        SC2SDTransformationMatrix.push_back(zeroMatrix);
        SC2SDTransformationMatrixInverse.push_back(zeroMatrix);

        continue;
     }

     // declare this variable separately for the cases where the first plane hit != 1, so then the previous plane is plane 0 - manually set below. 
     // Since missed intermediate planes are still traced to, this doesn't have any use beyond that first hit and is defined obviously.
     int previousPlaneHit = planeNum-1; 

     // for the tracing part, all intermediate planes are tracked to, even if the plane wasn't hit in the input data
     // Want the previous step matrices to be defined on the 0 'starting' plane. Eg. If plane 7 was the first plane hit, then plane 0, the starting plane, is defined at a small distance in front of the 2nd module, etc.
     if (planeNum == firstPlaneOfFirstModuleHit) previousPlaneHit = 0;

     // Keep track of how many steps geant4 takes between planes - usually it's a small number dependent on the intermediate material/volume regions.
     int stepNumBetweenPlanes = 0; 

     // Re-Build and set my own target here without the build target method for plane to plane propagation.
     double truthPlaneTargetX=0.;
     double truthPlaneTargetY=0.;
     double truthPlaneTargetZ=0.;

     wirePlaneZPositions.at(planeNum) = geaneTrackerZPositions_.at(planeNum); // Filled every passNum loop which isn't necessary since the x positions don't change.

     // get point on target plane from a wire center
     short moduleNum = int((planeNum-1)/4.);
     short wireLayerNum = 0;

     gm2strawtracker::StrawView wireView;
     if (((planeNum-1)/2)%2 == 0) wireView = gm2strawtracker::u_view;
     else  wireView = gm2strawtracker::v_view;

     if (planeNum % 2 == 0) wireLayerNum = 1;

     WireID targetStraw(stationNumber, moduleNum, wireView, wireLayerNum, 16);

     auto cssStraw = detCoordMap_.find(Form("Module%d:%d", stationNumber, moduleNum))->second;
     auto tSInWorld = targetStraw.getCentreInWorld(cssStraw);

     truthPlaneTargetX = tSInWorld.x();
     truthPlaneTargetY = tSInWorld.y();
     truthPlaneTargetZ = tSInWorld.z();

     mf::LogTrace(name_) << "Starting position of free state (world): " << myFreeTrajState->GetPosition() << "\n";
     mf::LogTrace(name_) << "Starting momentum of free state (world): " << myFreeTrajState->GetMomentum() << "\n";

     /////////////////////////////////////////////////////////////////////////////////////
     // If the Geant SD recorded a too large of a starting point (323.1 mm instead of 323 mm), skip the event. Or skip if the propagation somehow jumps past the target without triggering it.
     //  Check this by converting coordinates to the tracker frame and checking Z values
     gm2geom::CoordSystem3Vector currentCoord(myFreeTrajState->GetPosition().x(), myFreeTrajState->GetPosition().y(), myFreeTrajState->GetPosition().z(), "world");
     gm2geom::CoordSystem3Vector targetCoord(truthPlaneTargetX, truthPlaneTargetY, truthPlaneTargetZ, "world");

     auto trackerCurrentCoord = currentCoord.transform(detCoordMap_.find("TrackerStation")->second, stationStr);
     auto trackerTargetCoord = targetCoord.transform(detCoordMap_.find("TrackerStation")->second, stationStr);

     if (trackerCurrentCoord.z() >= trackerTargetCoord.z()){ 
        mf::LogTrace(name_) << "Returned by Bad starting position. " << "\n";
        return 5; 
     }

     /////////////////////////////////////////////////////////////////////////////////////
     // free memory before reassigning
     if (theTarget){
      	delete theTarget;
      	theTarget = nullptr;
     }

     G4Point3D surfPos(truthPlaneTargetX,truthPlaneTargetY,truthPlaneTargetZ);
     theTarget = new gm2GeanePlaneSurfaceTarget(surfNorm, surfPos ); //surfNorm set in setupParams method
     g4edata->SetTarget( theTarget ); // Set the target to track to.

     ////////////////////////////////////////////////////////////////////////////////////
     // Propagate until G4ErrorTarget is reached step by step.

     geaneTrackUtils_.g4emgr->InitTrackPropagation();

     bool moreEvt = TRUE; // boolean for more steps to be taken before reaching target
     while( moreEvt )
     {
         mf::LogTrace(name_) << "//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << "\n";
         mf::LogTrace(name_) << "Step Propagating to Plane number: " << planeNum << " stepNumBetweenPlanes: " << stepNumBetweenPlanes << "\n";        

         mf::LogTrace(name_) << "\n" << "Propagate One Step with Free State" << "\n";

         ierr = geaneTrackUtils_.g4emgr->PropagateOneStep( myFreeTrajState, theG4ErrorMode );

         mf::LogTrace(name_) << " Step length from PropagateOneStep: " << myFreeTrajState->GetG4Track()->GetStepLength() << "\n";

         /////////////////////////////////////////////////////////////////////////////////////
         // Some failure modes - Can I move these out somewhere? Or fix them?
         // ierr is the failure integer which is non-zero when step length is super small, 
         // propagation is exactly along z, or a few other things. In this case transport matrices 
         // and things will get set to 0 and later on the program will seg fault, so break out 
         // and skip this event if it's non-zero here.
         if (ierr != 0){
            info << "ierr = " << ierr << " is non-zero, skipping event. \n";
            return 6;
         }

         // See some issues where geant returns a step length of 0 - sometimes it reaches the target and doesn't register, sometimes it keeps returning steps of 0 length.      
         if (myFreeTrajState->GetG4Track()->GetStepLength() == 0.0){
            moreEvt = 0;
            info << "Step length was seen to be 0 - Geant is having step length issues - return out and fail. \n";
            // break; // boundary is only seen to be reached a step later, which then changes the last transport matrix incorrectly, but can break in this case and salvage tracking since it does reach the target, however root cause is still unknown - return fail for now
            return 7;
         }
         // Some steps return the max step length and then the tracking fails - unknown root cause
         else if(myFreeTrajState->GetG4Track()->GetStepLength() == 1000){
              info << "Step length was seen to be 1000, and maxed out - Geant is having step length issues - return out and fail. \n";
              return 8;
         }

         // Some events end up trying to track in air which causes failures - possibly because of a really bad initial guess
         // It's possible some tracks might start in air but be fixed after an iteration, so if this causes too many failures, then remove it
         if (myFreeTrajState->GetG4Track()->GetMaterial()->GetName() == "G4_AIR"){
            info << "Particle tracking has ended up in air somehow, which shoudn't happen - return out and fail. \n";
            return 12;
         }

	 /////////////////////////////////////////////////////////////////////////////////////
         // Use step by step propagation to properly save the transfer matrices.

         gm2GeaneMatrix fillTransformationMatrix(5,0); // Remake this object everytime for each step before adding to the vector.

	 delete mySurfaceTrajState; mySurfaceTrajState = nullptr; // Free up memory before reassiging

         // Create a new surface trajectory state each step of the way from the propagating free trajectory state. The passed empty matrix gets filled with the Jacobian transformation matrix.
         mySurfaceTrajState = new gm2GeaneSurfaceTrajState(*myFreeTrajState, tracingVX, tracingWY, fillTransformationMatrix); 

         if (stepNumBetweenPlanes!=0) // Only want the coordinate transformation matrices on the planes, not the individual steps - so remove the last then add the most recent.
         {
            SC2SDTransformationMatrix.pop_back();
            SC2SDTransformationMatrixInverse.pop_back();
         }

         SC2SDTransformationMatrix.push_back(fillTransformationMatrix);

         int failureInt;
         SC2SDTransformationMatrixInverse.push_back((SC2SDTransformationMatrix.back()).inverse(failureInt));

         mf::LogTrace(name_) << "Plane number is: " << planeNum << " step num is: " << stepNumBetweenPlanes << " SC2SDTransformationMatrix sizes are: " 
                             << "\n" << SC2SDTransformationMatrix.size() << " " << SC2SDTransformationMatrixInverse.size() << "\n";

        // Fill transport matrices.
        // Transfer matrices uses push_back since they are only filled for hit planes.
        if (stepNumBetweenPlanes==0) 
        { 
           // For matrix multiplication, it should be: R10 = A1 * T10 * A0^-1, where R and T are the transport matrices in separate coordinate systems from plane 0 to 1, and A is the transformation between the two. 
           transferMatricesFreeSystem.push_back(myFreeTrajState->GetTransfMat()); // Gets the transport or transfer matrix for the last step. If there is only a single step between planes this is fine.
           transferMatrices.push_back(SC2SDTransformationMatrix.at(planeNum)*(myFreeTrajState->GetTransfMat())*SC2SDTransformationMatrixInverse.at(previousPlaneHit));              
        }
        else 
        {
           transferMatricesFreeSystem.at(planeNum) = myFreeTrajState->GetTransfMat()*transferMatricesFreeSystem.at(planeNum);

           // Here I gather up the transport matrices for the free system for all steps between planes, then I multiply on the right by the 
           // previous plane coordinate transformation matrix, and on the left by the next step coordinate transformation matrix. It keeps 
           // looping until the next step transformation matrix is the same as the next plane transformation matrix.
           transferMatrices.at(planeNum) = (SC2SDTransformationMatrix.at(planeNum)*(transferMatricesFreeSystem.at(planeNum))*SC2SDTransformationMatrixInverse.at(previousPlaneHit)); 
           // Have to be careful with plane to plane transport vs step by step transport. This is the right way to do it I believe.
        }

        stepNumBetweenPlanes++; // Increment the step number between planes. (The 1st step is the same as stepNumBetweenPlanes==0.)

        //---- Check if target is reached
        if( geaneTrackUtils_.g4emgr->GetPropagator()->CheckIfLastStep( myFreeTrajState->GetG4Track() )) 
        {
           geaneTrackUtils_.g4emgr->GetPropagator()->InvokePostUserTrackingAction( myFreeTrajState->GetG4Track() );  
           moreEvt = 0;
           mf::LogTrace(name_) << "STEP_BY_STEP propagation: Last Step " << "\n";
        }

	if(stepNumBetweenPlanes > 100){
	  info << "stepNumBetweenPlanes (" << stepNumBetweenPlanes << ") > 100.  Particle is likely spiralling - return out and fail. \n";
	  return 9;
	}

     } // end of while


     /////////////////////////////////////////////////////////////////////////////////////
     // After getting to the target, once again create the surface trajectory state on the target.

     mf::LogTrace(name_) << "\n" << "Create end Surface State for plane parameter gathering" << "\n";

     gm2GeaneMatrix unneededMatrix(5,0);
	    
     delete mySurfaceTrajState; mySurfaceTrajState = nullptr; // Free up memory before reassiging

     // The idea is to propagate in steps with the free state, then at the end on the plane turn that into a surface state and grab associated parameters and the error to fill predicted objects.
     mySurfaceTrajState = new gm2GeaneSurfaceTrajState(*myFreeTrajState, tracingVX, tracingWY, unneededMatrix); 

     trackParamPredicted[0].at(planeNum) = mySurfaceTrajState->GetParameters().GetInvP();

     // need to get position values from world frame and convert to tracker frame, just grabbing VW geane parameters may be off in magnitude
     G4Point3D predPosFromSurface = mySurfaceTrajState->GetPosition();

     gm2geom::CoordSystem3Vector worldPredPos(predPosFromSurface.x(), predPosFromSurface.y(), predPosFromSurface.z(), "world");
     auto trackerPredPos = worldPredPos.transform(detCoordMap_.find("TrackerStation")->second, stationStr);

     double predForwardMom = sqrt( 1./trackParamPredicted[0].at(planeNum)*1./trackParamPredicted[0].at(planeNum) - mySurfaceTrajState->GetParameters().GetPV()*mySurfaceTrajState->GetParameters().GetPV() - mySurfaceTrajState->GetParameters().GetPW()*mySurfaceTrajState->GetParameters().GetPW() );
            
     trackParamPredicted[1].at(planeNum) = (mySurfaceTrajState->GetParameters().GetPV()) / predForwardMom; // I think these are good
     trackParamPredicted[2].at(planeNum) = (mySurfaceTrajState->GetParameters().GetPW()) / predForwardMom;

     trackParamPredicted[3].at(planeNum) = trackerPredPos.x();
     trackParamPredicted[4].at(planeNum) = trackerPredPos.y();


     mf::LogTrace(name_) << " Predicted Track Parameters Plane: " << planeNum
            << " param 0: " << trackParamPredicted[0].at(planeNum)
            << " param 1: " << trackParamPredicted[1].at(planeNum)
            << " param 2: " << trackParamPredicted[2].at(planeNum)
            << " param 3: " << trackParamPredicted[3].at(planeNum)
            << " param 4: " << trackParamPredicted[4].at(planeNum)
            << "\n";

     gm2GeaneTrajErr errorEnd(5,0);
     errorEnd = mySurfaceTrajState->GetError();
     errorMatrices.push_back(errorEnd);

     mf::LogTrace(name_) << "\n" << "Error matrix plane: " << planeNum << "\n" << errorMatrices[planeNum] << "\n";
     mf::LogTrace(name_) << "\n" << "Transport matrix free plane: " << planeNum << "\n" << transferMatricesFreeSystem[planeNum] << "\n";
     mf::LogTrace(name_) << "\n" << "Transport matrix surface plane: " << planeNum << "\n" << transferMatrices[planeNum] << "\n";


  }  // End planeNum loop
      
  /////////////////////////////////////////////////////////////////////////////////////
  // Convert gm2GeaneMatrix objects to Eigen here now for the gm2dataproduct where there is no geant

  std::vector<Eigen::MatrixXd> myTransferMatricesGeVcm; // same names as objects in method below for same info
  std::vector<Eigen::MatrixXd> errorMatricesGeVcm;

   for(int ipl=0;ipl<geaneTrackUtils_.maxNumPlanes;ipl++) {
      myTransferMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));
      errorMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));

     for(int i=0;i<5;i++) {
       for(int j=0; j<5;j++) {
          myTransferMatricesGeVcm[ipl](i,j) = transferMatrices[ipl][i][j]*1.; // GEVCM (Don't scale them to MeV mm.)
          errorMatricesGeVcm[ipl](i,j) = errorMatrices[ipl][i][j]*1.; // GEVCM
       }
     } 
  } 

  track.geaneHits.geaneTransportMatrices   = myTransferMatricesGeVcm; // Fill trackDetails with predicted/propagated objects. 
  track.geaneHits.geaneErrorMatrices       = errorMatricesGeVcm;
  track.geaneHits.geanePredictedParameters = trackParamPredicted;
  track.geaneHits.planeZPositions          = wirePlaneZPositions; // Filled every passNum loop with the same values - should reduce.

  // Tidy up memory
  delete myFreeTrajState;
  delete mySurfaceTrajState;
  delete theTarget; theTarget = nullptr;

  return 0;

} // end errorPropagation method


void GeaneParamUtils::calcMeasuredParams(gm2strawtracker::TrackDetailArtRecord & trackFitDetails)
{

  // this method calculates measured paramaters for the geane track, where the measured parameters are the the wire postions +- the angle corrected dcas
  // also angle corrects the dca error
  double* emField = new double[6];

  for (int i = 0; i < int(trackFitDetails.trackPlanesHitList.size()); ++i) // loop over hit planes
  {

    int planeNum = trackFitDetails.trackPlanesHitList.at(i);

    double predictedZPosition = trackFitDetails.geaneHits.planeZPositions.at(planeNum);
    double predictedYPosition = trackFitDetails.geaneHits.geanePredictedParameters[4].at(planeNum);
    double predictedXPosition = trackFitDetails.geaneHits.geanePredictedParameters[3].at(planeNum);

    // double predPoint[4] = {predictedXPosition, predictedYPosition, predictedZPosition, 1}; // 1 for the time - might need to change with time varying fields 

    gm2geom::CoordSystem3Vector pointForFieldTracker(predictedXPosition, predictedYPosition, predictedZPosition, stationStr);
    auto pointForFieldWorld = pointForFieldTracker.transform(detCoordMap_.find("TrackerStation")->second, "world");
    double predPoint[4] = {pointForFieldWorld.x(), pointForFieldWorld.y(), pointForFieldWorld.z(), trackFitDetails.time};

    G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField()->GetFieldValue(predPoint, emField);

    double xField = emField[0]/tesla; 
    double yField = emField[1]/tesla;
    double zField = emField[2]/tesla; 

    double predictedMomentum = 1./trackFitDetails.geaneHits.geanePredictedParameters[0].at(planeNum);

    double predictedZMomentum = predictedMomentum/sqrt(1 + trackFitDetails.geaneHits.geanePredictedParameters[1].at(planeNum)*trackFitDetails.geaneHits.geanePredictedParameters[1].at(planeNum) + 
                                                           trackFitDetails.geaneHits.geanePredictedParameters[2].at(planeNum)*trackFitDetails.geaneHits.geanePredictedParameters[2].at(planeNum));
    double predictedYMomentum = trackFitDetails.geaneHits.geanePredictedParameters[2].at(planeNum) * predictedZMomentum;
    double predictedXMomentum = trackFitDetails.geaneHits.geanePredictedParameters[1].at(planeNum) * predictedZMomentum;

    double parMomentum; // get momentum and field components parallel and perpendicular to U and V measurement axes respectively, take dot products of U and V like vectors with field and momentum vectors
    double perpField;

    gm2geom::CoordSystem3Vector fieldVector(xField, yField, zField, "world");
    auto trackerFieldVector = fieldVector.transform(detCoordMap_.find("TrackerStation")->second, stationStr, true);

    if (geaneTrackUtils_.isUPlane(planeNum))
    {   
       // these are the correct dot products explicitly written out, there's a better way to do this
       // axis up and to the right, negative U, negative shouldn't actually matter since this is always squared
       parMomentum = -(geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,0)*predictedXMomentum + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(0,1)*predictedYMomentum); 
       perpField = std::cos(7.5*pi/180)*trackerFieldVector.y() + std::sin(7.5*pi/180)*trackerFieldVector.x();
    }
    else
    {
       parMomentum = -(geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,0)*predictedXMomentum + geaneTrackUtils_.XYtoUVcoordinateTransformationMatrix(1,1)*predictedYMomentum); // axis down and to the right, negative V
       perpField = std::cos(7.5*pi/180)*trackerFieldVector.y() - std::sin(7.5*pi/180)*trackerFieldVector.x();
    }

    double circleMomentum = sqrt(parMomentum*parMomentum + predictedZMomentum*predictedZMomentum);

    double approxplaneFactorCorrection = sqrt(1 - (parMomentum/predictedMomentum)*(parMomentum/predictedMomentum) );

    double radius = 1000*(circleMomentum*1e6*1./299792458)/abs(perpField); // convert from MeV to eV, divide by c, then divide by the e*field, then multiple by 1000 to get mm

    double dca = trackFitDetails.measuredDCAs[planeNum];
    // double dcaCorrected = dca/approxplaneFactorCorrection; // straight line correction to dca

    double rightside = -radius * sqrt(1 - (parMomentum/predictedMomentum)*(parMomentum/predictedMomentum)) + sqrt( dca*dca + 2 * dca * radius + (1 - (parMomentum/predictedMomentum)*(parMomentum/predictedMomentum)) * radius * radius);
    double leftside  = +radius * sqrt(1 - (parMomentum/predictedMomentum)*(parMomentum/predictedMomentum)) - sqrt( dca*dca - 2 * dca * radius + (1 - (parMomentum/predictedMomentum)*(parMomentum/predictedMomentum)) * radius * radius);

    // sometimes tracks with poor initial guesses (usually with low numbers of hits) can result in predicted parameters in regions of zero field
    // just default positions to uncorrected dcas and fit on that before coming back around again
    if (yField == 0 && xField == 0 && zField == 0)
    {
       rightside = dca;
       leftside  = dca;
    }

    /////////////////////////////////////////////////////////////////////////////////////

    if (geaneTrackUtils_.isUPlane(planeNum))
    {
       if (trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gLeft) trackFitDetails.geaneHits.geaneMeasuredParameters[3][planeNum] = trackFitDetails.geaneHits.geaneWireUVPositions[3][planeNum] + leftside;
       else if (trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gRight) trackFitDetails.geaneHits.geaneMeasuredParameters[3][planeNum] = trackFitDetails.geaneHits.geaneWireUVPositions[3][planeNum] - rightside; // add or subtract correction based on truth
       else if (trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gCenter) trackFitDetails.geaneHits.geaneMeasuredParameters[3][planeNum] = trackFitDetails.geaneHits.geaneWireUVPositions[3][planeNum];
       else if (trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gNA_side || trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gUnknown){
              throw cet::exception(name_) << "Track sides are na or unknown which shouldn't happen in angularCorrection. planeNum: " << planeNum << " Uside: " << trackFitDetails.geaneHits.geaneHitSides.at(planeNum) << " \n";
       }
    }
    else if(!geaneTrackUtils_.isUPlane(planeNum))
    {
       if (trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gLeft) trackFitDetails.geaneHits.geaneMeasuredParameters[4][planeNum] = trackFitDetails.geaneHits.geaneWireUVPositions[4][planeNum] + leftside; 
       else if (trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gRight) trackFitDetails.geaneHits.geaneMeasuredParameters[4][planeNum] = trackFitDetails.geaneHits.geaneWireUVPositions[4][planeNum] - rightside; // add or subtract correction based on truth
       else if (trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gCenter) trackFitDetails.geaneHits.geaneMeasuredParameters[4][planeNum] = trackFitDetails.geaneHits.geaneWireUVPositions[4][planeNum];
       else if (trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gNA_side || trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gUnknown){
            throw cet::exception(name_) << "Track sides are na or unknown which shouldn't happen in angularCorrection. planeNum: " << planeNum << " Vside: " << trackFitDetails.geaneHits.geaneHitSides.at(planeNum) << " \n";
       }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // for errors, it seems fine to use the approximate straight line correction - higher order becomes much more difficult and is probably unnessecary 

    if(!(trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gCenter))   trackFitDetails.UVerrors[planeNum] = trackFitDetails.uncorrErrors[planeNum]/approxplaneFactorCorrection;
    else if(trackFitDetails.geaneHits.geaneHitSides.at(planeNum) == gCenter) trackFitDetails.UVerrors[planeNum] = trackFitDetails.uncorrErrors[planeNum]; // if hit is set at the center of the wire, don't correct error

  }

  delete[] emField;

  return;
}
  

}//namespace
