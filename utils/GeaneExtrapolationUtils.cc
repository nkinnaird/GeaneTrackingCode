//! header file
#include "GeaneExtrapolationUtils.hh"
//! namespace
using namespace gm2strawtracker;
using gm2geom::CoordSystem3Vector;

//! standard constructor
GeaneExtrapolationUtils::GeaneExtrapolationUtils()
   : name_("GeaneExtrapolationUtils")
   , cs_()
   , fhicl_()
   , gm2consts_()
   , fieldManager_()
   , nav_(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking())
   , geomUtils_()
   , sgeom_()
   , theTarget(0)
   , theTarget2(0)
   , theTarget3(0)
   , theTarget4(0)  
   , geaneTrackUtils_()
   , world()
{
  art::ServiceHandle<artg4::DetectorHolderService> dh;

  if(!dh->isInitialized()){
    dh->initialize();
    dh->constructAllLVs();
  }

  G4VPhysicalVolume* world = dh->worldPhysicalVolume();

  art::ServiceHandle<gm2geom::gm2FieldManager> fMgr;
  G4TransportationManager::GetTransportationManager()->SetFieldManager(fMgr->GetUnifiedFieldManager());
  G4TransportationManager::GetTransportationManager()->GetFieldManager()->SetChordFinder(fMgr->GetUnifiedFieldManager()->GetChordFinder());

  geaneTrackUtils_.g4emgr = gm2GeanePropagatorManager::GetErrorPropagatorManager();
  g4edata = G4ErrorPropagatorData::GetErrorPropagatorData();

  g4edata->SetVerbose(0);

  geaneTrackUtils_.g4emgr->SetUserInitialization(world);
  geaneTrackUtils_.g4emgr->InitGeant4e();

  //G4Navigator* nav_ = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  nav_->SetWorldVolume(world);

  //std::ostringstream oss;
  //oss << "/tracking/verbose " << trackingVerbose_;
  //std::string var = oss.str();
  //G4UImanager::GetUIpointer()->ApplyCommand(var);  
} 


// reconstruct the track vertex (decay) position
bool GeaneExtrapolationUtils::reconstructDecayVertex( const TrackArtRecord &track, DecayVertexArtRecord& decayVertex, const GeaneFhiclOptions& fclOptions, Eigen::MatrixXd initialError, DecayVertexArtRecord& vertexSteps) 
{
   std::cout <<"Ryan UV Error XY "<<initialError<<"\n\n";

  //qqq
   Eigen::MatrixXd covMatrixXYEigen = geaneTrackUtils_.JacobianToUV.inverse() * initialError * geaneTrackUtils_.JacobianToUV.transpose().inverse(); 
   initialError = covMatrixXYEigen;
   std::cout <<"Ryan Initial Error XY "<<covMatrixXYEigen<<"\n\n";
    //  Eigen::MatrixXd initialError = covMatrixXYEigen;

  mf::LogInfo(name_) << "in GeaneExtrapolationUtils::reconstructDecayVertex for trackArtRecord" << "\n";

  
  printf("In Extrapolation utils code \n");
  if (decayVertex.volumesHit.size() !=0) decayVertex.volumesHit.clear();
  
  // initialize
  bool success = true;

  // get the required fhicl parameters
  bool useTruth    = fclOptions.useTruth;
  bool useOnlyPDPs = fclOptions.useOnlyPDPs;
  int minNumPlanes = fclOptions.minNumPlanes;

  std::string extrapolater = fclOptions.extrapolater;

  // momentum cut
  double pmin = fclOptions.pmin;
  double pmax = fclOptions.pmax;
  
  // p value cut
  bool cutPoorPValues = fclOptions.cutPoorPValues;
  double pValueCut = fclOptions.pValueCut;
  

  if(extrapolater == "GeaneExtrapolater") {
    printf("We have chosen the Geane extrapolator \n");
    // initialize information to get from TrackArtRecord
    art::Ptr<gm2strawtracker::StrawTimeIslandArtRecord> timeIsl;
    gm2strawtracker::StrawDigitPtrCollection digits;
    art::Ptr<gm2strawtracker::StrawDigitArtRecord> startDigit;
    gm2strawtracker::StrawMCDigit startMCDigit;
    art::Ptr<gm2truth::StrawArtRecord> startHit;
    
    timeIsl = track.island;
    digits  = timeIsl->strawDigits;
    
    std::vector<double> covTotalInverseData;
    CoordSystem3Vector startMomentum;
    CoordSystem3Vector startPosition;
    double startTime(0);
    double station(0);

    std::vector< gm2strawtracker::ExtrapolationStep > extrapolationSteps;
    gm2strawtracker::ExtrapolationStep newStep;
    
    // if extrapolating backwards, flip the hit's momentum
    if (useTruth) {
      int previousModuleNumber = 9; // higher than can be
      for (auto dig : digits) {
	if(useOnlyPDPs) {
	  if (dig->strawMCDigit.strawMCHits.at(0)->particle_name != "e+" || dig->strawMCDigit.strawMCHits.at(0)->parent_ID != 1) {
	    decayVertex.failureMode = "simParticleNotPDP";
	    success = false;
	    return success;
	  }
	}
	for (auto hit : dig->strawMCDigit.strawMCHits) {
	  mf::LogDebug(name_) << "hit->moduleNumber = " << hit->moduleNumber << ", previousModuleNumber = " << previousModuleNumber << "\n";
	  if (hit->moduleNumber < previousModuleNumber) {
	    startHit = hit;
	    previousModuleNumber = hit->moduleNumber;
	  }
	}
      }
      
      // set the initial position and momentum for the extrapolation
      //startPosition.set(startHit->x_global , startHit->y_global , startHit->z_global , "world");
      //startMomentum.set(startHit->px_global, startHit->py_global, startHit->pz_global, "world");
      startTime = startHit->time;
      station = startHit->stationNumber;
      
      startPosition = startHit->position.transform(cs_, "world");
      startMomentum = startHit->momentum.transform(cs_, "world", true);
      
    }
    
    else {
      // if not using truth, use the reconstructed track momentum/position etc
      startMomentum = track.momentum;
      startPosition = track.states.positionVect.front();
      startTime = track.time;
      station = track.island->station;

      if (startMomentum.mag() < pmin || startMomentum.mag() > pmax) {
	decayVertex.failureMode = "failedMomCut";
	success = false;
	return success;
      }
      
      mf::LogInfo(name_) << "startMomentum = " << startMomentum << "\n";
      mf::LogInfo(name_) << "startPosition = " << startPosition << "\n";
      decayVertex.startMomentum = startMomentum;
      decayVertex.startPosition = startPosition;
      std::cout<<"decayVertex.startPosition: "<<decayVertex.startPosition<<std::endl;
      std::cout<<"decayVertex.startMomentum: "<<decayVertex.startMomentum<<std::endl;
      // get the covariance matrix of the track
      if( track.trackPlanesHitList.size() ) {
	covTotalInverseData = track.geaneHits.covarianceTotalInverseData;
      
	
	// check the p-value of the track
	if (cutPoorPValues) {
	  if (track.pValue < pValueCut) {
	    mf::LogInfo(name_) << "pValue = " << track.pValue << ", which is < " << pValueCut << " so return" << "\n";
	    decayVertex.failureMode = "poorPValue";
	    success = false;
	    return success;
	  }
	}
      }
      
      //TrackMatrix5x5 covMatrix = track.covMatrix;
      //if (covMatrix != 0) {
      //  std::cout << covMatrix << "\n";
      //}
      
      //Eigen::MatrixXd true0Covariance = track.covarianceTotalInverse.inverse(); 
      //Eigen::MatrixXd covMeVmm (5,5);
      //Eigen::MatrixXd variances(5,5);
      //variances.setZero();
      //variances(0,0) = sqrt(true0Covariance(0,0)); 
      //variances(1,1) = sqrt(true0Covariance(1,1)); 
      //variances(2,2) = sqrt(true0Covariance(2,2)); 
      //variances(3,3) = sqrt(true0Covariance(3,3)); 
      //variances(4,4) = sqrt(true0Covariance(4,4)); 
      //Eigen::MatrixXd corrMatrix(5,5);
      //corrMatrix = variances.inverse() * true0Covariance * variances.inverse();
      //covMeVmm(0,0) = true0Covariance(0,0) * (1/1000)*(1/1000);
      //covMeVmm(0,1) = true0Covariance(0,1) * (1/1000);
      //covMeVmm(0,2) = true0Covariance(0,2) * (1/1000);
      //covMeVmm(0,3) = true0Covariance(0,3) * (1/1000) * 10;
      //covMeVmm(0,4) = true0Covariance(0,4) * (1/1000) * 10;
      //covMeVmm(1,0) = true0Covariance(1,0) * (1/1000);
      //covMeVmm(2,0) = true0Covariance(2,0) * (1/1000);
      //covMeVmm(3,0) = true0Covariance(3,0) * (1/1000) * 10;
      //covMeVmm(4,0) = true0Covariance(4,0) * (1/1000) * 10;
      //covMeVmm(1,3) = true0Covariance(1,3) * 10;
      //covMeVmm(1,4) = true0Covariance(1,4) * 10;
      //covMeVmm(2,3) = true0Covariance(2,3) * 10;
      //covMeVmm(2,4) = true0Covariance(2,4) * 10;
      //covMeVmm(3,1) = true0Covariance(3,1) * 10;
      //covMeVmm(4,1) = true0Covariance(4,1) * 10;
      //covMeVmm(3,2) = true0Covariance(3,2) * 10;
      //covMeVmm(4,2) = true0Covariance(4,2) * 10;
      //covMeVmm(3,3) = true0Covariance(3,3) * 10 * 10;
      //covMeVmm(3,4) = true0Covariance(3,4) * 10 * 10;
      //covMeVmm(4,3) = true0Covariance(4,3) * 10 * 10;
      //covMeVmm(4,4) = true0Covariance(4,4) * 10 * 10;
      
      //	corrMatrix.print();
      
      // check that there are enough track states
      if (track.trackNumPlanesHit < minNumPlanes) {
	mf::LogInfo(name_) << "track.states.size() = " << track.trackNumPlanesHit << " which is < minNumPlanes = " << minNumPlanes << " so return false." << "\n";
	decayVertex.failureMode = "poorPValue";
	success = false;
	return success;
      }
    }
    
    // set the information for the extrapolation step
    newStep.momentum = startMomentum;
    newStep.position = startPosition;
    newStep.time = startTime;

    G4double initialPosition[4] = {startPosition[0], startPosition[1], startPosition[2], 0.0};
    G4double Bfield[3] = {0.0, 0.0, 0.0};
    fieldManager_->GetFieldValue(initialPosition, Bfield);
    printf("Finished printing start position and momentum\n");
    printf("Extrapolation set up complete. Now do extrapolation \n");

    // reconstruct vertex given the pos/mom/time of start point
    success = reconstructDecayVertex(newStep,decayVertex,station,fclOptions, initialError, vertexSteps);
  }
  
  else {
    decayVertex.failureMode = "extrapolaterTypeNotSet";
    success = false;
  }

  return success;

}

bool GeaneExtrapolationUtils::reconstructDecayVertex( const ExtrapolationStep &startPoint, DecayVertexArtRecord& decayVertex, int stationNumber, const GeaneFhiclOptions& fclOptions, Eigen::MatrixXd initialError, DecayVertexArtRecord& vertexSteps) 
{
  
  mf::LogInfo(name_) << "in GeaneExtrapolationUtils::reconstructDecayVertex for ExtrapolationStep" << "\n";
  
  printf("We reconstruct! Decay vertex \n");

  double extrapolationTime = startPoint.time;


  std::string extrapolationDirection = "backwards";
  double totalDistance = fclOptions.totalDistance;
  unsigned int nSteps = fclOptions.nSteps;
  std::string extrapolater = fclOptions.extrapolater;
  double stepSize = fclOptions.stepSize;
  double smallStepSize = fclOptions.smallStepSize;

  ExtrapolationStep prevStep;
  ExtrapolationStep afterStep;

  //bool oldCode = true;
  bool success = true;

  std::string direction = "backwards";

  // initialize empty vector of extrapolation steps to be saved to the decay vertex
  std::vector< gm2strawtracker::ExtrapolationStep > extrapolationSteps;
  std::vector< gm2strawtracker::ExtrapolationStep > uncorrectedExtrapolationSteps;
  std::vector< gm2strawtracker::ExtrapolationStep > extrapolationStepsInVolume;
  std::vector< gm2strawtracker::ExtrapolationStep > stepsInVolume;
  // set the step size
  if (smallStepSize == 0.0) stepSize = totalDistance / nSteps;
  else nSteps = totalDistance / smallStepSize;
  mf::LogInfo(name_) << "stepSize = " << stepSize << ", nSteps = " << nSteps << "\n";  
  
  if(extrapolater == "GeaneExtrapolater"){

    printf("Choosing GeaneExtrapolater\n");

    //Track parameters
    gm2geom::CoordSystem3Vector startPosition = startPoint.position;
    gm2geom::CoordSystem3Vector startMomentum = startPoint.momentum;
    gm2strawtracker::ExtrapolationStep timeStep;

    //Extrapolater parameters
    gm2strawtracker::ExtrapolationStep firstStep;
    int numsteps = 0;
    stringstream stationStream;
    std::string stationStr;
    stationStr = stationStream.str();

    //Surface vectors - first the vertical
    tracingWY.setX(0);  //tracing vectors defined in .hh
    tracingWY.setY(1);
    tracingWY.setZ(0);
    std::cout<< "tracingWY "<<tracingWY << "\n";

    //u_view defined where
    printf("station number %d \n",stationNumber);
    WireID straw1(stationNumber, 0,  gm2strawtracker::u_view, 0,0);
    printf("After straw1 \n");
    WireID straw2(stationNumber, 0,  gm2strawtracker::u_view, 0,31);
    printf("After straw2 \n");


    // for (auto iter = detCoordMap_.begin()+1; iter != detCoordMap_.end(); iter++)
    // {
    //   cout << "Key: " << iter->first;
    //   cout << " v1: " <<  straw1.getCentreInWorld(iter->second) << endl;
    //   // cout << "Key: " << iter->first << " instancename: " <<  iter->second.coordInstanceName_ << endl;
    // }


std::cout << std::endl << "Debug Got Here 1" << std::endl;



    //detCoordMap_ defined in .hh
    auto cssStraw = detCoordMap_.find(Form("Module%d:%d", stationNumber, 0))->second;
    printf("After cssstraw \n");
    auto firstV = straw1.getCentreInWorld(cssStraw);
    std::cout<< "firstV "<< firstV << "\n";

std::cout << std::endl << "Debug Got Here 2" << std::endl;


    auto secondV = straw2.getCentreInWorld(cssStraw);
    std::cout<< "secondV "<< secondV << "\n";
    auto thirdV = secondV - firstV;
    tracingVX.setX(thirdV.x()/thirdV.mag());
    tracingVX.setY(0);
    tracingVX.setZ(thirdV.z()/thirdV.mag());

    std::cout<< "tracingVX "<<tracingVX << "\n";

    auto crossV = tracingVX.cross(tracingWY);

    std::cout<< "crossV "<<crossV << "\n";

    //    set the normal vectors (Z in tracker frame) for target planes (surfnorm defined in hh
    surfNorm.setX(crossV.x());
    surfNorm.setY(crossV.y());
    surfNorm.setZ(crossV.z());


    //Break track momentum and position vectors into components so we can convert from CoordSystem3Vector to G4ThreeVector
    double startingXPos = startPosition[0];
    double startingYPos = startPosition[1];
    double startingZPos = startPosition[2];
    double startingXMom = -1*startMomentum[0];
    double startingYMom = -1*startMomentum[1];
    double startingZMom = -1*startMomentum[2];
    G4ThreeVector xv3( startingXPos, startingYPos, startingZPos );
    G4ThreeVector pv3( startingXMom, startingYMom, startingZMom );


    cout << endl << endl << "/////////////////////////////////////////////////////////////////////////////////////" << endl << endl;
    cout << endl << endl << "/////////////////////////////////////////////////////////////////////////////////////" << endl << endl;

    cout << "Error input: " << endl;

    //Create a surface state to go with this position, momentum
    gm2GeaneTrajErr errorss( 5, 0);
    //Error matrix carries over from the track fit
    for(int i=0;i<5;i++){
      for(int j=0;j<5;j++){
	      errorss[i][j] = initialError(i,j);
        cout << initialError(i,j) << " ";
      }
      cout << endl;
    }



    gm2GeaneSurfaceTrajState *yourSurfaceTrajState = new gm2GeaneSurfaceTrajState( "e+", xv3, pv3, tracingVX, tracingWY, errorss);
    // gm2GeaneSurfaceTrajState *yourSurfaceTrajState = new gm2GeaneSurfaceTrajState( "gamma", xv3, pv3, tracingVX, tracingWY, errorss);
    gm2GeaneFreeTrajState *myFreeTrajState = new gm2GeaneFreeTrajState(*yourSurfaceTrajState);

      gm2GeaneMatrix fillTransformationMatrix0temp(5,0); // 0th step transformation matrix.
    gm2GeaneSurfaceTrajState *testSurfaceState = new gm2GeaneSurfaceTrajState(*myFreeTrajState, tracingVX, tracingWY, fillTransformationMatrix0temp);

    gm2GeaneFreeTrajState *myFreeTrajState2 = new gm2GeaneFreeTrajState(*testSurfaceState);

      gm2GeaneMatrix fillTransformationMatrix3temp(5,0); // 0th step transformation matrix.
    gm2GeaneSurfaceTrajState *testSurfaceState3 = new gm2GeaneSurfaceTrajState(*myFreeTrajState2, tracingVX, tracingWY, fillTransformationMatrix3temp);

    cout << "surface state 1: " << endl;
    cout << *yourSurfaceTrajState;

    cout << "surface state 2: " << endl;
    cout << *testSurfaceState;

    cout << "surface state 3: " << endl;
    cout << *testSurfaceState3;

    // cout << "surface state 1 error: " << endl;
    // cout << yourSurfaceTrajState->GetError() << endl;

    // cout << "surface state 2 error: " << endl;
    // cout << testSurfaceState->GetError() << endl;

    // cout << "surface state 3 error: " << endl;
    // cout << testSurfaceState3->GetError() << endl;

    // cout << "surface state 1 parameters: " << yourSurfaceTrajState->GetParameters() << endl;
    // cout << "surface state 2 parameters: " << testSurfaceState->GetParameters() << endl;

    // cout << "surface state 1 vector V: " << yourSurfaceTrajState->GetVectorV() << endl;
    // cout << "surface state 2 vector V: " << testSurfaceState->GetVectorV() << endl;

    // cout << "surface state 1 vector W: " << yourSurfaceTrajState->GetVectorW() << endl;
    // cout << "surface state 2 vector W: " << testSurfaceState->GetVectorW() << endl;


    cout << endl << endl << "/////////////////////////////////////////////////////////////////////////////////////" << endl << endl;
    cout << endl << endl << "/////////////////////////////////////////////////////////////////////////////////////" << endl << endl;



    //    gm2GeaneTrajErr error( 5, 0 );
    //   gm2GeaneFreeTrajState* myFreeTrajState = new gm2GeaneFreeTrajState("e+", xv3, pv3, error );
    theG4ErrorMode = G4ErrorMode_PropBackwards;
    std::cout<<"initial Error surface: "<<yourSurfaceTrajState->GetError() <<std::endl;
    std::cout<<"initial Error free: "<<myFreeTrajState->GetError() <<std::endl;

    //    gm2GeaneTrajErr startingError(5);
    //std::cout<<"startingError: "<<startingError<<std::endl;
    //    for(int i=0;i<5;i++){
    //  for(int j=0;j<5;j++){
    //	startingError[i][j] = initialError(i,j);
    //  }
    // }
    std::cout<<"startingError: "<<myFreeTrajState->GetError()<<std::endl;
    //myFreeTrajState->SetError(startingError);

    std::vector<gm2GeaneMatrix> transferMatrices;
    std::vector<gm2GeaneMatrix> targetTransferMatrices;
    std::vector<gm2GeaneMatrix> errorMatrices;
    std::vector<gm2GeaneMatrix> targetErrorMatrices;

    errorMatrices.push_back(myFreeTrajState->GetError());
    transferMatrices.push_back(myFreeTrajState->GetTransfMat());
    printf("Pushing initial Transfer Matrix \n");
    std::cout<<"Transfer Matrix: "<<myFreeTrajState->GetTransfMat()<<std::endl;

    int stepNumBetweenPlanes = 0;
    
    //Failure integer. Event will be skipped if this is non-zero
    int ierr = 0;

    //Boolean used to signal beginning and end of first while loop, the extrapolation to a cylinder target
    bool stepping1 = true;

    std::string volumeHit = "";
    std::string volumeHitBefore = "";
    bool hitVol = false;
    decayVertex.volumesHit.clear();


    //If a target is already set, delete it
    //if(theTarget){
    //delete theTarget;
    //theTarget = nullptr;
    //}

    //Cylinder parameters (values are in mm)
    G4double cylradius = 7062.0*mm;
    G4ThreeVector cyltrans = G4ThreeVector(0,0,0);
    G4RotationMatrix cylrotm = G4RotationMatrix();
    cylrotm.rotateX(90.*deg);
    cylrotm.rotateY(0.*deg);
    cylrotm.rotateZ(0.*rad);    
    //End of parameter initialization

    //Create cylinder target
    theTarget = new gm2GeaneCylSurfaceTarget(cylradius, cyltrans, cylrotm);

    //Set cylinder as target
    g4edata->SetTarget( theTarget );

    geaneTrackUtils_.g4emgr->InitTrackPropagation();

    measuredTrackLength = 0.0;

    //First while loop that extrapolates backwards to our cylinder target which represents the uniform magnetic field inside the storage ring
    while( stepping1 ) {
      //Break out of loop after max number of steps
      auto currentPosition = myFreeTrajState->GetPosition();
      auto currentMomentum = myFreeTrajState->GetMomentum();

      //std::cout<<"currentPosition: "<<currentPosition<<std::endl;
      //std::cout<<"currentMomentum: "<<currentMomentum<<std::endl<<std::endl;

      if(numsteps == 1000){
	stepping1 = false;
      }

      //Increment runner variable for number of steps
      numsteps++;

      G4ThreeVector beforePos = currentPosition;

      //The heart of the extrapolater, extrapolates one step so we loop this until the target is reached
      ierr = geaneTrackUtils_.g4emgr->PropagateOneStep( myFreeTrajState, theG4ErrorMode );

      errorMatrices.push_back(myFreeTrajState->GetError());
      transferMatrices.push_back(myFreeTrajState->GetTransfMat());
      printf("Pushing Main Stepping 1 Transfer Matrix \n");
    std::cout<<"Transfer Matrix: "<<myFreeTrajState->GetTransfMat()<<std::endl;

      //Check for errors
      if (ierr != 0){
	mf::LogTrace(name_) << "ierr = " << ierr << " is non-zero, skipping event and returning False. \n";
	std::cout << "ierr = " << ierr << " is non-zero, skipping event and returning False." << std::endl;
	success = false;
	return success;
      }
#if 0
      if (myFreeTrajState->GetG4Track()->GetStepLength() == 0.0){
	std::cout<<"stepLength==0"<<std::endl;
	success = false;
	return success;
      }

      if (myFreeTrajState->GetG4Track()->GetStepLength() == 1000){
	std::cout<<"stepLength==1000"<<std::endl;
	success = false;
	return success;
      }

      if(myFreeTrajState->GetG4Track()->GetMaterial()->GetName() == "G4_AIR"){
	std::cout<<"G4_air"<<std::endl;
	success = false;
	return success;
      }
#endif
      currentPosition = myFreeTrajState->GetPosition();
      currentMomentum = myFreeTrajState->GetMomentum();
      G4ThreeVector g4StepPos(currentPosition[0], currentPosition[1], currentPosition[2]);
      CoordSystem3Vector stepPos;
      stepPos.set(currentPosition[0], currentPosition[1], currentPosition[2]);
      CoordSystem3Vector stepMom;
      stepMom.set(currentMomentum[0], currentMomentum[1], currentMomentum[2]);

      ExtrapolationStep currentStep;
      currentStep.position = stepPos;
      currentStep.momentum = stepMom;
      extrapolationSteps.push_back(currentStep);
      uncorrectedExtrapolationSteps.push_back(currentStep);
      

      if(sqrt(pow(currentPosition[0],2)+pow(currentPosition[2],2)) > 7320){
	std::cout<<"Outside ring so track killed"<<std::endl;
	success = true;
	decayVertex.position = stepPos;
	decayVertex.momentum = stepMom;
	decayVertex.position.coordSystemName = "world";
	decayVertex.momentum.coordSystemName = "world";
	decayVertex.time = extrapolationTime;
	vertexSteps.steps = stepsInVolume;
	return success;
      }

      G4TouchableHandle theTouchableBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetTouchableHandle();
      G4TouchableHandle theTouchable = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetTouchableHandle();
      const G4NavigationHistory *history = theTouchable->GetHistory();
      const G4NavigationHistory *historyBefore = theTouchableBefore->GetHistory();
      G4int depth = history->GetDepth();
      G4int depthBefore = historyBefore->GetDepth();

      volumeHit = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetPhysicalVolume()->GetName();
      volumeHitBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetName();
      G4VPhysicalVolume* physicalVolumeHit = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetPhysicalVolume();
      G4VPhysicalVolume* physicalVolumeHitBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetPhysicalVolume();

      if(volumeHit=="supportPostLVShell"){
	G4AffineTransform transform = history->GetTransform(depth-1).Inverse();
	G4ThreeVector postLocation = transform.TransformPoint(physicalVolumeHit->GetTranslation());
	std::cout<<"postLocation: "<<postLocation<<std::endl;
	CoordSystem3Vector postLocationVector;
	postLocationVector.set(postLocation[0], postLocation[1], postLocation[2]);
	CoordSystem3Vector impactVector = stepPos - postLocationVector;
	std::cout<<"impactVector: "<<impactVector<<std::endl;
	decayVertex.impactVectors.push_back(impactVector);
      }

      if(volumeHitBefore=="supportPostLVShell"){
	G4AffineTransform transform = historyBefore->GetTransform(depthBefore-1).Inverse();
	G4ThreeVector postLocation = transform.TransformPoint(physicalVolumeHitBefore->GetTranslation());
	std::cout<<"postLocation: "<<postLocation<<std::endl;
	CoordSystem3Vector postLocationVector;
	postLocationVector.set(postLocation[0], postLocation[1], postLocation[2]);
	CoordSystem3Vector impactVector = stepPos - postLocationVector;
	std::cout<<"impactVector: "<<impactVector<<std::endl;
	decayVertex.impactVectors.push_back(impactVector);
      }
 
      std::cout<<"physicalVolumeHit:"<<volumeHit<<std::endl;
      std::cout<<"physicalVolumeHit translation: "<<physicalVolumeHit->GetTranslation()<<std::endl;
      hitVol = didStepHitVolume(volumeHit);
      if(hitVol){
	decayVertex.hitVolume = true;
	ExtrapolationStep volumeStep;
	volumeStep.position = stepPos;
	volumeStep.momentum = stepMom;
	volumeStep.volume = volumeHit;
	if(myFreeTrajState->GetG4Track()->GetStepLength() != 0.0){
	  stepsInVolume.push_back(volumeStep);
	}
	if(std::find(decayVertex.volumesHit.begin(),decayVertex.volumesHit.end(),volumeHit) == decayVertex.volumesHit.end()) {
	  decayVertex.volumesHit.push_back(volumeHit);
	  decayVertex.uncorrectedVolumesHit.push_back(volumeHit);
	}
      }
      else if (!decayVertex.hitVolume) {
	decayVertex.hitVolume = false;
      }

      G4ThreeVector afterPos = currentPosition;
      G4ThreeVector stepDistance = afterPos - beforePos;
      extrapolationTime -= 1e9 * (0.001*stepDistance.mag()) / gm2consts_->cLight();
      
      //Increment second runner variable
      stepNumBetweenPlanes++;

      //Check to see if we hit target and end stepping1
      if( geaneTrackUtils_.g4emgr->GetPropagator()-> CheckIfLastStep( myFreeTrajState->GetG4Track() ))
	{
	  targetTransferMatrices.push_back(myFreeTrajState->GetTransfMat());
      printf("Hit target 1 Push Transfer Matrix \n");
    std::cout<<"Target Transfer Matrix: "<<myFreeTrajState->GetTransfMat()<<std::endl;

	  targetErrorMatrices.push_back(myFreeTrajState->GetError());
	  geaneTrackUtils_.g4emgr->GetPropagator()->InvokePostUserTrackingAction( myFreeTrajState->GetG4Track() );
	  stepping1 = false;
	}

    }//end of stepping1

    //Extrapolation Math Parameters
    auto finalPosition = myFreeTrajState->GetPosition();
    auto finalMomentum = myFreeTrajState->GetMomentum();

    std::cout<<"RM: Final Position on inner cylinder: "<<finalPosition<<"\n\n";
    finalPosition[1] = 0.0;   // mm
    finalMomentum[1] = 0.0;   // MeV/c
    
    long double eCharge;
    eCharge = 1.60217662e-19; // C

    long double uniformBField;
    long double radiusOfCurviture;

    uniformBField = 1.451; // T
    radiusOfCurviture = finalMomentum.mag() / (eCharge * uniformBField); // (MeV/c) / (C*T)
    radiusOfCurviture = (5.344286e-19) * radiusOfCurviture; //mm
    
    G4ThreeVector ri{(-1*finalMomentum[2]), 0.0, finalMomentum[0]}; //(MeV/c)
    ri = (ri * radiusOfCurviture) / finalMomentum.mag(); //mm
    
    G4ThreeVector Rs(0.0, 0.0, 0.0);
    Rs = finalPosition - ri; //mm
    
    G4ThreeVector Rf(0.0, 0.0, 0.0);
    Rf = (Rs.mag() + radiusOfCurviture) * (Rs / Rs.mag()); //mm

    long double trackLength;
    G4ThreeVector displacementVector(0.0, 0.0, 0.0);
    displacementVector = Rf - finalPosition; //mm
    trackLength = displacementVector.mag(); //mm
    trackLength = trackLength / radiusOfCurviture; //unitless
    trackLength = pow(trackLength, 2); // unitless
    trackLength = acos(1-(.5*trackLength)); //unitless
    trackLength = radiusOfCurviture * trackLength; //mm

    //if(theTarget2){
    //delete theTarget2;
    //theTarget2 = nullptr;
    //}
    
    std::cout<<"stepping2 trackLength: "<<trackLength<<std::endl;

    //Create track length target
    theTarget2 = new gm2GeaneTrackLengthTarget(trackLength+(myFreeTrajState->GetG4Track()->GetTrackLength()));

    //Set track length target
    g4edata->SetTarget( theTarget2 );

    //Boolean for second while loop, stepping2
    bool stepping2;
    stepping2 = true;
    bool foundMin = false;
    bool pastTangent = false;
    //Reset runner variables
    numsteps = 0;
    stepNumBetweenPlanes = 0;
    
    while(stepping2){
      auto currentPosition = myFreeTrajState->GetPosition();
      auto currentMomentum = myFreeTrajState->GetMomentum();

      //std::cout<<"currentPosition: "<<currentPosition<<std::endl;
      //std::cout<<"currentMomentum: "<<currentMomentum<<std::endl<<std::endl;


      G4ThreeVector g4Position{currentPosition[0],currentPosition[1],currentPosition[2]};
      G4ThreeVector g4Momentum{currentMomentum[0],currentMomentum[1],currentMomentum[2]};
      double stepRadialMomBefore = gm2consts_->ComputePrhat(&g4Position,&g4Momentum);
      G4ThreeVector beforePos = currentPosition;
      if(numsteps == 1000){
	stepping2 = false;
	break;
      }

#if 0
      if((((cos(theta)*currentMomentum[0])+(sin(theta)*currentMomentum[2]))/sqrt(pow(currentMomentum[0],2)+pow(currentMomentum[2],2))) < .001){
	std::cout<<"At tangent point"<<std::endl;
	stepping2 = false;
	break;
      }
#endif
      //Increment runner variable
      numsteps++;

      //Propagate particle one step
      ierr = geaneTrackUtils_.g4emgr->PropagateOneStep( myFreeTrajState, theG4ErrorMode );

      errorMatrices.push_back(myFreeTrajState->GetError());
      transferMatrices.push_back(myFreeTrajState->GetTransfMat());
       printf("Pushing Main Stepping 2 Transfer Matrix \n");
    std::cout<<"Transfer Matrix: "<<myFreeTrajState->GetTransfMat()<<std::endl;
     if(ierr != 0){
	success = false;
	return success;
      }
#if 0
      if (myFreeTrajState->GetG4Track()->GetStepLength() == 0.0){
	success = false;
	return success;
      }

      if (myFreeTrajState->GetG4Track()->GetStepLength() == 1000){
	success = false;
	return success;
      }
      
      if(myFreeTrajState->GetG4Track()->GetMaterial()->GetName() == "G4_AIR"){
	success = false;
	return success;
      }
#endif
      currentPosition = myFreeTrajState->GetPosition();
      currentMomentum = myFreeTrajState->GetMomentum(); 
      G4ThreeVector g4PositionAfter{currentPosition[0],currentPosition[1],currentPosition[2]};
      G4ThreeVector g4MomentumAfter{currentMomentum[0],currentMomentum[1],currentMomentum[2]};
      double stepRadialMomAfter = gm2consts_->ComputePrhat(&g4PositionAfter,&g4MomentumAfter);

      if(fabs(stepRadialMomAfter) < fabs(stepRadialMomBefore)){
	foundMin = true;
      }
      
      if((foundMin) && (fabs(stepRadialMomAfter) > fabs(stepRadialMomBefore))){
	pastTangent = true;
      }

      G4ThreeVector g4StepPos(currentPosition[0], currentPosition[1], currentPosition[2]);
      CoordSystem3Vector stepPos;
      stepPos.set(currentPosition[0], currentPosition[1], currentPosition[2]);
      CoordSystem3Vector stepMom;
      stepMom.set(currentMomentum[0], currentMomentum[1], currentMomentum[2]);
 
      ExtrapolationStep currentStep;
      currentStep.position = stepPos;
      currentStep.momentum = stepMom;
      extrapolationSteps.push_back(currentStep);
      uncorrectedExtrapolationSteps.push_back(currentStep);
      


      if(sqrt(pow(currentPosition[0],2)+pow(currentPosition[2],2)) > 7320){
	std::cout<<"Outside ring so track killed"<<std::endl;
	success = true;
	decayVertex.position = stepPos;
	decayVertex.momentum = stepMom;
	decayVertex.position.coordSystemName = "world";
	decayVertex.momentum.coordSystemName = "world";
	decayVertex.time = extrapolationTime;
	vertexSteps.steps = stepsInVolume;
	return success;
      }

      G4TouchableHandle theTouchable = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetTouchableHandle();
      G4TouchableHandle theTouchableBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetTouchableHandle();
      const G4NavigationHistory *history = theTouchable->GetHistory();
      const G4NavigationHistory *historyBefore = theTouchableBefore->GetHistory();
      G4int depth = history->GetDepth();
      G4int depthBefore = historyBefore->GetDepth();

      volumeHit = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetPhysicalVolume()->GetName();
      volumeHitBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetName();

      G4VPhysicalVolume* physicalVolumeHit = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetPhysicalVolume();
      G4VPhysicalVolume* physicalVolumeHitBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetPhysicalVolume();

      if(volumeHit=="supportPostLVShell"){
	G4AffineTransform transform = history->GetTransform(depth-1).Inverse();
	G4ThreeVector postLocation = transform.TransformPoint(physicalVolumeHit->GetTranslation());
	std::cout<<"postLocation: "<<postLocation<<std::endl;
	CoordSystem3Vector postLocationVector;
	postLocationVector.set(postLocation[0], postLocation[1], postLocation[2]);
	CoordSystem3Vector impactVector = stepPos - postLocationVector;
	std::cout<<"impactVector: "<<impactVector<<std::endl;
	decayVertex.impactVectors.push_back(impactVector);
      }

      if(volumeHitBefore=="supportPostLVShell"){
	G4AffineTransform transform = historyBefore->GetTransform(depthBefore-1).Inverse();
	G4ThreeVector postLocation = transform.TransformPoint(physicalVolumeHitBefore->GetTranslation());
	std::cout<<"postLocation: "<<postLocation<<std::endl;
	CoordSystem3Vector postLocationVector;
	postLocationVector.set(postLocation[0], postLocation[1], postLocation[2]);
	CoordSystem3Vector impactVector = stepPos - postLocationVector;
	std::cout<<"impactVector: "<<impactVector<<std::endl;
	decayVertex.impactVectors.push_back(impactVector);
      }

      hitVol = didStepHitVolume(volumeHit);
      if(hitVol){
	decayVertex.hitVolume = true;
	ExtrapolationStep volumeStep;
	volumeStep.position = stepPos;
	volumeStep.momentum = stepMom;
	volumeStep.volume = volumeHit;
	if(myFreeTrajState->GetG4Track()->GetStepLength() != 0.0){
	  stepsInVolume.push_back(volumeStep);
	}
	if(std::find(decayVertex.volumesHit.begin(),decayVertex.volumesHit.end(),volumeHit) == decayVertex.volumesHit.end()) {
	  decayVertex.volumesHit.push_back(volumeHit);
	  decayVertex.uncorrectedVolumesHit.push_back(volumeHit);
	}
      }
      else if (!decayVertex.hitVolume) {
	decayVertex.hitVolume = false;
      }


      G4ThreeVector afterPos = currentPosition;
      G4ThreeVector stepDistance = afterPos - beforePos;
      extrapolationTime -= 1e9 * (0.001*stepDistance.mag()) / gm2consts_->cLight();

      //Increment second runner variable
      stepNumBetweenPlanes++;

      //Check to see if we hit target and end stepping2
      if( geaneTrackUtils_.g4emgr->GetPropagator()->CheckIfLastStep( myFreeTrajState->GetG4Track() ))
	{
	  targetTransferMatrices.push_back(myFreeTrajState->GetTransfMat());
          printf("Hit target stepping 2. Pushing Stepping 2 Transfer Matrix \n");
         std::cout<<"Target Transfer Matrix: "<<myFreeTrajState->GetTransfMat()<<std::endl;

	  targetErrorMatrices.push_back(myFreeTrajState->GetError());
	  geaneTrackUtils_.g4emgr->GetPropagator()->InvokePostUserTrackingAction( myFreeTrajState->GetG4Track() );
	  stepping2 = false;
	}

    }//end of stepping2

    auto correctionPosition = myFreeTrajState->GetPosition();
    auto correctionMomentum = myFreeTrajState->GetMomentum();

    G4Vector3D stepping3Momentum = -1*correctionMomentum;
    G4Vector3D stepping3Position = correctionPosition;
    
    correctionPosition[1] = 0.0;
    correctionMomentum[1] = 0.0;

    radiusOfCurviture = correctionMomentum.mag() / (eCharge * uniformBField);
    radiusOfCurviture = (5.344286e-19) * radiusOfCurviture;

    ri[0] = -1*correctionMomentum[2];
    ri[1] = 0.0;
    ri[2] = correctionMomentum[0];
    ri = (ri * radiusOfCurviture) / correctionMomentum.mag();

    Rs = correctionPosition - ri;

    Rf = (Rs.mag() + radiusOfCurviture) * (Rs / Rs.mag());

    displacementVector = Rf - correctionPosition;
    trackLength = displacementVector.mag();
    trackLength = trackLength / radiusOfCurviture;
    trackLength = pow(trackLength, 2);
    trackLength = acos(1-(.5*trackLength));
    trackLength = radiusOfCurviture * trackLength;
    
    std::cout<<"stepping3 trackLength: "<<trackLength<<std::endl;

    //if(theTarget3){
    //delete theTarget3;
    //theTarget3 = nullptr;
    //}

    theTarget3 = new gm2GeaneTrackLengthTarget(trackLength+(myFreeTrajState->GetG4Track()->GetTrackLength()));

    g4edata->SetTarget( theTarget3 );
 
    bool stepping3;
    stepping3 = true;

    std::cout<<"Stepping2 past tangent: "<<pastTangent<<std::endl;
    if(pastTangent){
      myFreeTrajState->SetMomentum(stepping3Momentum);
      myFreeTrajState->GetG4Track()->SetMomentumDirection(-1*(myFreeTrajState->GetG4Track()->GetMomentumDirection()));

      theG4ErrorMode = G4ErrorMode_PropForwards;
    }

    numsteps = 0;
    stepNumBetweenPlanes = 0;

    foundMin = false;
    pastTangent = false;

    while(stepping3){
      numsteps++;
      auto currentPosition = myFreeTrajState->GetPosition();
      auto currentMomentum = myFreeTrajState->GetMomentum();

      //std::cout<<"currentPosition: "<<currentPosition<<std::endl;
      //std::cout<<"currentMomentum: "<<currentMomentum<<std::endl<<std::endl;

      G4ThreeVector g4Position{currentPosition[0],currentPosition[1],currentPosition[2]};
      G4ThreeVector g4Momentum{currentMomentum[0],currentMomentum[1],currentMomentum[2]};
      double stepRadialMomBefore = gm2consts_->ComputePrhat(&g4Position,&g4Momentum);

      G4ThreeVector beforePos = currentPosition;

      ierr = geaneTrackUtils_.g4emgr->PropagateOneStep( myFreeTrajState, theG4ErrorMode );

      errorMatrices.push_back(myFreeTrajState->GetError());
      transferMatrices.push_back(myFreeTrajState->GetTransfMat());
      printf("Pushing Main Stepping 3 Transfer Matrix \n");
      std::cout<<"Transfer Matrix: "<<myFreeTrajState->GetTransfMat()<<std::endl;

      if(ierr != 0){
	std::cout<<"ierr != 0"<<std::endl;
	success = false;
	return success;
      }
#if 0
      if (myFreeTrajState->GetG4Track()->GetStepLength() == 0.0){
	std::cout<<"StepLength == 0.0"<<std::endl;
	success = false;
	return success;
      }

      if (myFreeTrajState->GetG4Track()->GetStepLength() == 1000){
	std::cout<<"StepLength == 1000"<<std::endl;
	success = false;
	return success;
      }
      
      if(myFreeTrajState->GetG4Track()->GetMaterial()->GetName() == "G4_AIR"){
	std::cout<<"In G4_AIR"<<std::endl;
	success = false;
	return success;
      }
#endif
      stepNumBetweenPlanes++;
      currentPosition = myFreeTrajState->GetPosition();
      currentMomentum = myFreeTrajState->GetMomentum();
      G4ThreeVector g4PositionAfter{currentPosition[0],currentPosition[1],currentPosition[2]};
      G4ThreeVector g4MomentumAfter{currentMomentum[0],currentMomentum[1],currentMomentum[2]};
      double stepRadialMomAfter = gm2consts_->ComputePrhat(&g4PositionAfter,&g4MomentumAfter);

      if(fabs(stepRadialMomAfter) < fabs(stepRadialMomBefore)){
	foundMin = true;
      }
      
      if((foundMin) && (fabs(stepRadialMomAfter) > fabs(stepRadialMomBefore))){
	pastTangent = true;
      }

      G4ThreeVector g4StepPos(currentPosition[0], currentPosition[1], currentPosition[2]);
      CoordSystem3Vector stepPos;
      stepPos.set(currentPosition[0], currentPosition[1], currentPosition[2]);
      CoordSystem3Vector stepMom;
      stepMom.set(currentMomentum[0], currentMomentum[1], currentMomentum[2]);

      ExtrapolationStep currentStep;
      currentStep.position = stepPos;
      currentStep.momentum = stepMom;
      extrapolationSteps.push_back(currentStep);
      uncorrectedExtrapolationSteps.push_back(currentStep);


      if(sqrt(pow(currentPosition[0],2)+pow(currentPosition[2],2)) > 7320){
	std::cout<<"Outside ring so track killed"<<std::endl;
	success = true;
	decayVertex.position = stepPos;
	decayVertex.momentum = stepMom;
	decayVertex.position.coordSystemName = "world";
	decayVertex.momentum.coordSystemName = "world";
	decayVertex.time = extrapolationTime;
	vertexSteps.steps = stepsInVolume;
	return success;
      }

      G4TouchableHandle theTouchable = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetTouchableHandle();
      G4TouchableHandle theTouchableBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetTouchableHandle();
      const G4NavigationHistory *history = theTouchable->GetHistory();
      const G4NavigationHistory *historyBefore = theTouchableBefore->GetHistory();
      G4int depth = history->GetDepth();
      G4int depthBefore = historyBefore->GetDepth();

      volumeHit = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetPhysicalVolume()->GetName();
      volumeHitBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetName();

      G4VPhysicalVolume* physicalVolumeHit = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetPhysicalVolume();
      G4VPhysicalVolume* physicalVolumeHitBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetPhysicalVolume();

      if(volumeHit=="supportPostLVShell"){
	G4AffineTransform transform = history->GetTransform(depth-1).Inverse();
	G4ThreeVector postLocation = transform.TransformPoint(physicalVolumeHit->GetTranslation());
	std::cout<<"postLocation: "<<postLocation<<std::endl;
	CoordSystem3Vector postLocationVector;
	postLocationVector.set(postLocation[0], postLocation[1], postLocation[2]);
	CoordSystem3Vector impactVector = stepPos - postLocationVector;
	std::cout<<"impactVector: "<<impactVector<<std::endl;
	decayVertex.impactVectors.push_back(impactVector);
      }

      if(volumeHitBefore=="supportPostLVShell"){
	G4AffineTransform transform = historyBefore->GetTransform(depthBefore-1).Inverse();
	G4ThreeVector postLocation = transform.TransformPoint(physicalVolumeHitBefore->GetTranslation());
	std::cout<<"postLocation: "<<postLocation<<std::endl;
	CoordSystem3Vector postLocationVector;
	postLocationVector.set(postLocation[0], postLocation[1], postLocation[2]);
	CoordSystem3Vector impactVector = stepPos - postLocationVector;
	std::cout<<"impactVector: "<<impactVector<<std::endl;
	decayVertex.impactVectors.push_back(impactVector);
      }

      hitVol = didStepHitVolume(volumeHit);
      if(hitVol){
	decayVertex.hitVolume = true;
	ExtrapolationStep volumeStep;
	volumeStep.position = stepPos;
	volumeStep.momentum = stepMom;
	volumeStep.volume = volumeHit;
	if(myFreeTrajState->GetG4Track()->GetStepLength() != 0.0){
	  stepsInVolume.push_back(volumeStep);
	}
	if(std::find(decayVertex.volumesHit.begin(),decayVertex.volumesHit.end(),volumeHit) == decayVertex.volumesHit.end()) {
	  decayVertex.volumesHit.push_back(volumeHit);
	  decayVertex.uncorrectedVolumesHit.push_back(volumeHit);
	}
      }
      else if (!decayVertex.hitVolume) {
	decayVertex.hitVolume = false;
      }

      G4ThreeVector afterPos = currentPosition;
      G4ThreeVector stepDistance = afterPos - beforePos;
      extrapolationTime += 1e9 * (0.001*stepDistance.mag()) / gm2consts_->cLight();

      if( geaneTrackUtils_.g4emgr->GetPropagator()->CheckIfLastStep( myFreeTrajState->GetG4Track() ))
	{
	  targetTransferMatrices.push_back(myFreeTrajState->GetTransfMat());
          printf("Pushing Stepping 3 Transfer Matrix \n");
          std::cout<<"Target Transfer Matrix: "<<myFreeTrajState->GetTransfMat()<<std::endl;

	  targetErrorMatrices.push_back(myFreeTrajState->GetError());
	  geaneTrackUtils_.g4emgr->GetPropagator()->InvokePostUserTrackingAction( myFreeTrajState->GetG4Track() );
	  stepping3 = false;
	}
      
    }//end of stepping3

    correctionPosition = myFreeTrajState->GetPosition();
    correctionMomentum = myFreeTrajState->GetMomentum();

    CoordSystem3Vector uncorrectedPosition;
    uncorrectedPosition.set(correctionPosition[0],correctionPosition[1],correctionPosition[2]);
    CoordSystem3Vector uncorrectedMomentum;
    uncorrectedMomentum.set(correctionMomentum[0],correctionMomentum[1],correctionMomentum[2]);

    decayVertex.uncorrectedPosition = uncorrectedPosition;
    decayVertex.uncorrectedMomentum = uncorrectedMomentum;
    std::cout<<"uncorrectedExtrapolationSteps size: "<<uncorrectedExtrapolationSteps.size()<<std::endl;
    decayVertex.uncorrectedSteps = uncorrectedExtrapolationSteps;
    std::cout<<"uncorrectedSteps size: "<<decayVertex.uncorrectedSteps.size()<<std::endl;

    G4Vector3D stepping4Momentum = -1*correctionMomentum;
    G4Vector3D stepping4Position = correctionPosition;

    correctionPosition[1] = 0.0;
    correctionMomentum[1] = 0.0;

    radiusOfCurviture = correctionMomentum.mag() / (eCharge * uniformBField);
    radiusOfCurviture = (5.344286e-19) * radiusOfCurviture;

    ri[0] = -1*correctionMomentum[2];
    ri[1] = 0.0;
    ri[2] = correctionMomentum[0];
    ri = (ri * radiusOfCurviture) / correctionMomentum.mag();

    Rs = correctionPosition - ri;

    Rf = (Rs.mag() + radiusOfCurviture) * (Rs / Rs.mag());

    displacementVector = Rf - correctionPosition;
    trackLength = displacementVector.mag();
    trackLength = trackLength / radiusOfCurviture;
    trackLength = pow(trackLength, 2);
    trackLength = acos(1-(.5*trackLength));
    trackLength = radiusOfCurviture * trackLength;
    
    std::cout<<"stepping4 trackLength: "<<trackLength<<std::endl;

    theTarget4 = new gm2GeaneTrackLengthTarget(trackLength+(myFreeTrajState->GetG4Track()->GetTrackLength()));

    g4edata->SetTarget( theTarget4 );
 
    bool stepping4;
    stepping4 = true;

    std::cout<<"Stepping3 past tangent: "<<pastTangent<<std::endl;
    if(pastTangent){
      myFreeTrajState->SetMomentum(stepping4Momentum);
      myFreeTrajState->GetG4Track()->SetMomentumDirection(-1*(myFreeTrajState->GetG4Track()->GetMomentumDirection()));
      
      if(theG4ErrorMode == G4ErrorMode_PropForwards){
	theG4ErrorMode = G4ErrorMode_PropBackwards;
      }
      else{
	theG4ErrorMode = G4ErrorMode_PropForwards;
      }
    }
    
    numsteps = 0;
    stepNumBetweenPlanes = 0;

    while(stepping4){
      numsteps++;
      auto currentPosition = myFreeTrajState->GetPosition();
      auto currentMomentum = myFreeTrajState->GetMomentum();

      //std::cout<<"currentPosition: "<<currentPosition<<std::endl;
      //std::cout<<"currentMomentum: "<<currentMomentum<<std::endl<<std::endl;


      G4ThreeVector beforePos = currentPosition;

      ierr = geaneTrackUtils_.g4emgr->PropagateOneStep( myFreeTrajState, theG4ErrorMode );

      errorMatrices.push_back(myFreeTrajState->GetError());
      transferMatrices.push_back(myFreeTrajState->GetTransfMat());
      printf("Pushing Main Stepping 4 Transfer Matrix \n");
      std::cout<<"Transfer Matrix: "<<myFreeTrajState->GetTransfMat()<<std::endl;

      if(ierr != 0){
	std::cout<<"ierr != 0"<<std::endl;
	success = false;
	return success;
      }
#if 0
      if (myFreeTrajState->GetG4Track()->GetStepLength() == 0.0){
	std::cout<<"StepLength == 0.0"<<std::endl;
	success = false;
	return success;
      }

      if (myFreeTrajState->GetG4Track()->GetStepLength() == 1000){
	std::cout<<"StepLength == 1000"<<std::endl;
	success = false;
	return success;
      }
      
      if(myFreeTrajState->GetG4Track()->GetMaterial()->GetName() == "G4_AIR"){
	std::cout<<"In G4_AIR"<<std::endl;
	success = false;
	return success;
      }
#endif
      stepNumBetweenPlanes++;
      currentPosition = myFreeTrajState->GetPosition();
      currentMomentum = myFreeTrajState->GetMomentum();
      
      G4ThreeVector g4StepPos(currentPosition[0], currentPosition[1], currentPosition[2]);
      CoordSystem3Vector stepPos;
      stepPos.set(currentPosition[0], currentPosition[1], currentPosition[2]);
      CoordSystem3Vector stepMom;
      stepMom.set(currentMomentum[0], currentMomentum[1], currentMomentum[2]);

      ExtrapolationStep currentStep;
      currentStep.position = stepPos;
      currentStep.momentum = stepMom;
      extrapolationSteps.push_back(currentStep);
      


      if(sqrt(pow(currentPosition[0],2)+pow(currentPosition[2],2)) > 7320){
	std::cout<<"Outside ring so track killed"<<std::endl;
	success = true;
	decayVertex.position = stepPos;
	decayVertex.momentum = stepMom;
	decayVertex.position.coordSystemName = "world";
	decayVertex.momentum.coordSystemName = "world";
	decayVertex.time = extrapolationTime;
	vertexSteps.steps = stepsInVolume;
	return success;
      }

      G4TouchableHandle theTouchable = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetTouchableHandle();
      G4TouchableHandle theTouchableBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetTouchableHandle();
      const G4NavigationHistory *history = theTouchable->GetHistory();
      const G4NavigationHistory *historyBefore = theTouchableBefore->GetHistory();
      G4int depth = history->GetDepth();
      G4int depthBefore = historyBefore->GetDepth();

      volumeHit = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetPhysicalVolume()->GetName();
      volumeHitBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetName();

      G4VPhysicalVolume* physicalVolumeHit = myFreeTrajState->GetG4Track()->GetStep()->GetPostStepPoint()->GetPhysicalVolume();
      G4VPhysicalVolume* physicalVolumeHitBefore = myFreeTrajState->GetG4Track()->GetStep()->GetPreStepPoint()->GetPhysicalVolume();

      if(volumeHit=="supportPostLVShell"){
	G4AffineTransform transform = history->GetTransform(depth-1).Inverse();
	G4ThreeVector postLocation = transform.TransformPoint(physicalVolumeHit->GetTranslation());
	std::cout<<"postLocation: "<<postLocation<<std::endl;
	CoordSystem3Vector postLocationVector;
	postLocationVector.set(postLocation[0], postLocation[1], postLocation[2]);
	CoordSystem3Vector impactVector = stepPos - postLocationVector;
	std::cout<<"impactVector: "<<impactVector<<std::endl;
	decayVertex.impactVectors.push_back(impactVector);
      }

      if(volumeHitBefore=="supportPostLVShell"){
	G4AffineTransform transform = historyBefore->GetTransform(depthBefore-1).Inverse();
	G4ThreeVector postLocation = transform.TransformPoint(physicalVolumeHitBefore->GetTranslation());
	std::cout<<"postLocation: "<<postLocation<<std::endl;
	CoordSystem3Vector postLocationVector;
	postLocationVector.set(postLocation[0], postLocation[1], postLocation[2]);
	CoordSystem3Vector impactVector = stepPos - postLocationVector;
	std::cout<<"impactVector: "<<impactVector<<std::endl;
	decayVertex.impactVectors.push_back(impactVector);
      }

      hitVol = didStepHitVolume(volumeHit);
      if(hitVol){
	decayVertex.hitVolume = true;
	ExtrapolationStep volumeStep;
	volumeStep.position = stepPos;
	volumeStep.momentum = stepMom;
	volumeStep.volume = volumeHit;
	if(myFreeTrajState->GetG4Track()->GetStepLength() != 0.0){
	  stepsInVolume.push_back(volumeStep);
	}
	if(std::find(decayVertex.volumesHit.begin(),decayVertex.volumesHit.end(),volumeHit) == decayVertex.volumesHit.end()) {
	  decayVertex.volumesHit.push_back(volumeHit);
	}
      }
      else if (!decayVertex.hitVolume) {
	decayVertex.hitVolume = false;
      }

      G4ThreeVector afterPos = currentPosition;
      G4ThreeVector stepDistance = afterPos - beforePos;
      extrapolationTime += 1e9 * (0.001*stepDistance.mag()) / gm2consts_->cLight();

      if( geaneTrackUtils_.g4emgr->GetPropagator()->CheckIfLastStep( myFreeTrajState->GetG4Track() ))
	{
	  targetTransferMatrices.push_back(myFreeTrajState->GetTransfMat());
	  printf("Last step stepping 4. Push transfer matrix \n");
          std::cout<<"Target Transfer Matrix: "<<myFreeTrajState->GetTransfMat()<<std::endl;
	  targetErrorMatrices.push_back(myFreeTrajState->GetError());
	  geaneTrackUtils_.g4emgr->GetPropagator()->InvokePostUserTrackingAction( myFreeTrajState->GetG4Track() );
	  stepping4 = false;
	}
      
    }//end of stepping4

    auto currentPosition = myFreeTrajState->GetPosition();
    auto currentMomentum = myFreeTrajState->GetMomentum();

    std::cout<<"myFreeTrajState->GetMomentum() "<< myFreeTrajState->GetMomentum()<<std::endl;

    G4Vector3D normalToFinalPlane;
    G4Vector3D parallelToFinalPlane;
    G4Vector3D radialVect;

    normalToFinalPlane[0] = currentMomentum[0];
    normalToFinalPlane[1] = currentMomentum[1];
    normalToFinalPlane[2] = currentMomentum[2];
    normalToFinalPlane = normalToFinalPlane / normalToFinalPlane.mag();

    radialVect[0] = currentPosition[0];
    radialVect[1] = currentPosition[1];
    radialVect[2] = currentPosition[2];

    radialVect = radialVect / radialVect.mag();

    parallelToFinalPlane = radialVect.cross(normalToFinalPlane);

    std::cout<<"NormalToFinalPlane"<< normalToFinalPlane << std::endl;
    std::cout<<"ParallelToFinalPlane"<< parallelToFinalPlane << std::endl;
    std::cout<<"radialVect"<< radialVect << std::endl;



    gm2GeaneMatrix fillTransformationMatrix(5,0);
    std::cout << "fillTransformationMatrix (B)"<< fillTransformationMatrix << std::endl;
    gm2GeaneSurfaceTrajState* mySurfaceTrajState = new gm2GeaneSurfaceTrajState(*myFreeTrajState, radialVect, parallelToFinalPlane, fillTransformationMatrix);
    std::cout << "fillTransformationMatrix (A)"<< fillTransformationMatrix << std::endl;



    targetErrorMatrices.push_back(mySurfaceTrajState->GetError());
    errorMatrices.push_back(mySurfaceTrajState->GetError());
    targetTransferMatrices.push_back(fillTransformationMatrix);
    transferMatrices.push_back(fillTransformationMatrix);

    gm2GeaneMatrix FinalProductMatrix = transferMatrices[1];
    int ntm = transferMatrices.size();
    printf("Number of transfer matrices of all kinds %d \n",ntm);
    for(uint i=2; i<transferMatrices.size();i++) {
      FinalProductMatrix = transferMatrices[i]*FinalProductMatrix;
    }

    std::cout<<"Final Transfer Matrix"<<std::endl;
    std::cout<<FinalProductMatrix<<std::endl;



    std::cout<<"targetErrorMatrices: "<<std::endl;
    for(uint i=0; i<targetErrorMatrices.size();i++){
      std::cout<<targetErrorMatrices[i]<<std::endl;
    }

    std::cout<<"errorMatrices: "<<std::endl;
    for(uint i=0; i<errorMatrices.size();i++){
      std::cout<<errorMatrices[i]<<std::endl;
    }

    std::cout<<"targetTransferMatrices: "<<std::endl;
    for(uint i=0; i<targetTransferMatrices.size();i++){
      printf(" i %d \n",i);
      std::cout<<targetTransferMatrices[i]<<std::endl;
    }


    std::cout<<"transferMatrices: "<<std::endl;
    for(uint i=0; i<transferMatrices.size();i++){
      printf(" i %d \n",i);
      std::cout<<transferMatrices[i]<<std::endl;
    }

    gm2geom::CoordSystem3Vector vectorPosition;
    gm2geom::CoordSystem3Vector vectorMomentum;

    vectorPosition[0] = currentPosition[0];
    vectorPosition[1] = currentPosition[1];
    vectorPosition[2] = currentPosition[2];

    vectorMomentum[0] = currentMomentum[0];
    vectorMomentum[1] = currentMomentum[1];
    vectorMomentum[2] = currentMomentum[2];

    decayVertex.position = vectorPosition;
    decayVertex.momentum = vectorMomentum;
    decayVertex.position.coordSystemName = "world";
    decayVertex.momentum.coordSystemName = "world";
    decayVertex.time = extrapolationTime;
    decayVertex.steps = extrapolationSteps;

    vertexSteps.steps = stepsInVolume;

    std::cout<<"Ryan decayVertex.position: "<<decayVertex.position<<std::endl;
    std::cout<<"Ryan decayVertex.momentum: "<<decayVertex.momentum<<std::endl;

    gm2GeaneTrajErr errorEnd(5,0);
    errorEnd = myFreeTrajState->GetError();
#if 0 
    Eigen::MatrixXd vertexError;
    for(int i=0;i<5;i++){
      for(int j=0;j<5;j++){
	vertexError(i,j) = errorEnd(i,j);
      }
    }
    decayVertex.error = vertexError;
#endif
    //std::cout<<"errorEnd: "<<errorEnd<<std::endl;
    for(std::vector<std::string>::const_iterator i = decayVertex.volumesHit.begin(); i != decayVertex.volumesHit.end(); ++i){
      std::cout<< "Volume Hit: "<< *i << std::endl;
    }
    std::cout<<"ryan volumesHit: "<<std::endl;
    for(auto volume : decayVertex.volumesHit){
      std::cout<<volume<<std::endl;
    }
    
#if 0
    G4double vertexInvP = 1 / currentMomentum.mag();
    G4double vertexDipAngle = atan( currentMomentum[1] / sqrt(pow(currentMomentum[0],2)+pow(currentMomentum[2],2)) );
    G4double vertexAzimuthalAngle = atan ( currentMomentum[2] / currentMomentum[0] );
    G4double vertexYPerp = -1*sin(vertexAzimuthalAngle)*currentPosition[0] + cos(vertexAzimuthalAngle)*currentPosition[2];
    G4double vertexZPerp = -1*sin(vertexDipAngle)*cos(vertexAzimuthalAngle)*currentPosition[0]+(-1)*sin(vertexDipAngle)*sin(vertexAzimuthalAngle)*currentPosition[2]+cos(vertexDipAngle)*currentPosition[1];
#endif

    //delete myFreeTrajState;
    //delete theTarget;
    //delete theTarget2;
    //delete theTarget3;
    //theTarget = nullptr;
    //theTarget2 = nullptr;
    //theTarget3 = nullptr;
    
    return success;

    } //end of GeaneExtrapolater block
  
  else {
    decayVertex.failureMode = "extrapolaterTypeNotSet";
    success = false;
  }
  
  if (decayVertex.position.coordSystemName == "NULL" || decayVertex.momentum.coordSystemName == "NULL") {
    mf::LogInfo(name_) << "decayVertex.position = " << decayVertex.position << "\n";
    mf::LogInfo(name_) << "decayVertex.momentum = " << decayVertex.momentum << "\n";
    //mf::LogInfo(name_) << "decayVertex.time = " << decayVertex.time << "\n";
    mf::LogInfo(name_) << "returning false as decay vertex was not correctly filled" << "\n";
    decayVertex.failureMode = "vertexVectorsNotFilled";
    success = false;
  }
  
  mf::LogInfo(name_) << "Returning from GeaneExtrapolationUtils::ReconstructDecayVertex with success = " << success << "\n";
  return success;
}

bool GeaneExtrapolationUtils::extrapolateToCalorimeter( const TrackArtRecord &track, DecayVertexArtRecord& extrapolatedCaloHit , const GeaneFhiclOptions& fclOptions) 
{
  
  mf::LogInfo(name_) << "in utils::extrapolateToCalorimeter(track,hit)" << "\n";

  bool success = false;
  bool useTruth    = fclOptions.useTruth;
  bool useOnlyPDPs = fclOptions.useOnlyPDPs;
  int minNumPlanes = fclOptions.minNumPlanes;

  std::string extrapolater = fclOptions.extrapolater;

  // momentum cut
  double pmin = fclOptions.pmin;
  double pmax = fclOptions.pmax;
  
  // p value cut
  bool cutPoorPValues = fclOptions.cutPoorPValues;
  double pValueCut = fclOptions.pValueCut;

  // initialize information to get from TrackArtRecord
  art::Ptr<gm2strawtracker::StrawTimeIslandArtRecord> timeIsl;
  gm2strawtracker::StrawDigitPtrCollection digits;
  art::Ptr<gm2strawtracker::StrawDigitArtRecord> startDigit;
  gm2strawtracker::StrawMCDigit startMCDigit;
  art::Ptr<gm2truth::StrawArtRecord> startHit;
  
  timeIsl = track.island;
  digits  = timeIsl->strawDigits;
  
  CoordSystem3Vector startMomentum;
  CoordSystem3Vector startPosition;
  double startTime(0);
  double station(0);
    
  if (useTruth) {
    
    mf::LogInfo(name_) << "Using truth" << "\n";
    
    if (digits.size() == 0) {
      mf::LogInfo(name_) << "no digits, returning false" << "\n";
      extrapolatedCaloHit.failureMode = "noDigits";
      success = false;
      return success;
    }
    
    int previousModuleNumber = -1; // smaller than can be
    for (auto dig : digits) {
      
      // if this digit has no mc information then return
      if (!dig.get()->hasMCHits()) {
	extrapolatedCaloHit.failureMode = "noMCHits";
	success = false;
	return success;
      }
      
      if(useOnlyPDPs) {
	if (dig->strawMCDigit.strawMCHits.at(0)->particle_name != "e+" || dig->strawMCDigit.strawMCHits.at(0)->parent_ID != 1) {
	  extrapolatedCaloHit.failureMode = "simParticleNotPDP";
	  success = false;
	  return success;
	}
      }
      for (auto hit : dig->strawMCDigit.strawMCHits) {
	if (hit->moduleNumber > previousModuleNumber) {
	  startHit = hit;
	  previousModuleNumber = hit->moduleNumber;
	}
      }
    }
    
    mf::LogDebug(name_) << "previousModuleNumber = " << previousModuleNumber << "\n";
    mf::LogDebug(name_) << "startHit->moduleNumber = " << startHit->moduleNumber << "\n";
    
    int minModuleNum = 0;
    if (startHit->moduleNumber < minModuleNum) {
      mf::LogInfo(name_) << "Last dig in track was not in module " << minModuleNum << "6 or 7, return" << "\n";
      extrapolatedCaloHit.failureMode = "lastDigitInWrongModule";
      success = false;
      return success;
    }
    
    //startPosition.set(startHit.get()->x_global, startHit.get()->y_global, startHit.get()->z_global,"world");
    //startMomentum.set(startHit.get()->px_global, startHit.get()->py_global, startHit.get()->pz_global,"world");
    startTime = startHit->time;
    station = startHit->stationNumber;

    startPosition = startHit->position.transform(cs_, "world");
    startMomentum = startHit->momentum.transform(cs_, "world", true);
    
  } // useTruth
  
  else {
    
    mf::LogInfo(name_) << "forwards, useTruth = " << useTruth << "\n";
    
    // check hit enough planes
    if (track.trackNumPlanesHit < minNumPlanes) {
      mf::LogInfo(name_) << "track.states.size() = " << track.trackNumPlanesHit << " which is < minNumPlanes = " << minNumPlanes << " so return false." << "\n";
      extrapolatedCaloHit.failureMode = "notEnoughTrackingPlanes";
      success = false;
      return success;
    }
    
    // check the p-value of the track
    if (cutPoorPValues) {
      if (track.pValue < pValueCut) {
	mf::LogInfo(name_) << "pValue = " << track.pValue << ", which is < " << pValueCut << " so return" << "\n";
	extrapolatedCaloHit.failureMode = "poorPValue";
	success = false;
	return success;
      }
    }
    
    startMomentum = track.states.momentumVect.back();
    startPosition = track.states.positionVect.back();
    startTime = track.island->meanTime;
    
    if (startMomentum.mag() < pmin || startMomentum.mag() > pmax) {
      extrapolatedCaloHit.failureMode = "failedMomCut";
      success = false;
      return success;
    }
    station = track.island->station;
  }

  ExtrapolationStep startPoint;
  startPoint.position = startPosition;
  startPoint.momentum = startMomentum;
  startPoint.time = startTime;

  mf::LogInfo(name_) << "startPoint.position = " << startPosition << "\n";
  mf::LogInfo(name_) << "startPoint.position = " << startMomentum << "\n";
  mf::LogInfo(name_) << "startPoint.time = " << startTime << "\n";
  mf::LogInfo(name_) << "station = " << station << "\n";  

  success = extrapolateToCalorimeter(startPoint, extrapolatedCaloHit, station, fclOptions);
  if (extrapolatedCaloHit.position.coordSystemName == "NULL" || extrapolatedCaloHit.momentum.coordSystemName == "NULL") success = false;
  return success;
}

bool GeaneExtrapolationUtils::extrapolateToCalorimeter( const ExtrapolationStep &startPoint, DecayVertexArtRecord& extrapolatedCaloHit , int stationNum, const GeaneFhiclOptions& fclOptions ) 
{

  mf::LogInfo(name_) << "in utils::extrapolateToCalorimeter(startPoint,hit) " << "\n";

  std::string extrapolationDirection = "forwards";

  double totalDistance = fclOptions.totalDistance;
  unsigned int nSteps = fclOptions.nSteps;
  std::string extrapolater = fclOptions.extrapolater;
  bool keepSteps = fclOptions.keepSteps;
  double targetZpos = fclOptions.targetZpos;
  double stepSize = fclOptions.stepSize;
  double smallStepSize = fclOptions.smallStepSize;

  double position_err[3] = {0.0,0.0,0.0};
  int station = stationNum;

  ExtrapolationStep prevStep;
  ExtrapolationStep afterStep;

  bool oldCode = true;
  bool success = true;
  bool leftFieldMap = false;
  // initialize empty vector of extrapolation steps to be saved to the decay vertex
  std::vector< gm2strawtracker::ExtrapolationStep > extrapolationSteps;
  std::vector< gm2strawtracker::ExtrapolationStep > extrapolationStepsInVolume;

  // set the step size
  if (smallStepSize == 0.0) stepSize = totalDistance / nSteps;
  else nSteps = totalDistance / smallStepSize;
  stepSize = smallStepSize;

  mf::LogDebug(name_) << "stepSize = " << stepSize << ", nSteps = " << nSteps << "\n";
  
  if (extrapolater == "GeaneExtrapolater") {
    
    gm2strawtracker::ExtrapolationStep firstStep;
    firstStep.position = startPoint.position;
    firstStep.momentum = startPoint.momentum;
    firstStep.time = startPoint.time;

    mf::LogDebug(name_) << "firstStep.position = " << firstStep.position << "\n";
    mf::LogDebug(name_) << "firstStep.momentum = " << firstStep.momentum << "\n";
    mf::LogDebug(name_) << "firstStep.position.transform(cs_,CalorimeterNumber[18]) = " << firstStep.position.transform(cs_,"CalorimeterNumber[18]") << "\n";
    mf::LogDebug(name_) << "firstStep.momentum.transform(cs_,CalorimeterNumber[18]) = " << firstStep.momentum.transform(cs_,"CalorimeterNumber[18]",true) << "\n";
    mf::LogDebug(name_) << "station = " << station << "\n";  

    // push step back to the collection to be saved to the vertex
    extrapolationSteps.push_back(firstStep);

    mf::LogDebug(name_) << "about to loop over steps (forwards)" << "\n";
    for (unsigned int i(1); i<nSteps+1; i++) {
      
      mf::LogDebug(name_) << "forwards step i = " << i << "\n";

      double local_pos_err[3] = {0.0,0.0,0.0};

      mf::LogDebug(name_) << "test0 extrapolationSteps[i-1].position = " << extrapolationSteps[i-1].position << "\n";
      mf::LogDebug(name_) << "test0 extrapolationSteps[i-1].momentum = " << extrapolationSteps[i-1].momentum << "\n";
      mf::LogDebug(name_) << "test0 extrapolationSteps[i-1].position.transform(cs_,CalorimeterNumber[18]) = " << extrapolationSteps[i-1].position.transform(cs_,"CalorimeterNumber[18]") << "\n";
      mf::LogDebug(name_) << "test0 extrapolationSteps[i-1].momentum.transform(cs_,CalorimeterNumber[18],true) = " << extrapolationSteps[i-1].momentum.transform(cs_,"CalorimeterNumber[18]",true) << "\n";
  
      mf::LogDebug(name_) << "extrapolationDirection = " << extrapolationDirection << "\n";

      gm2strawtracker::ExtrapolationStep step = oldCode ? oldGeaneExtrapolation(extrapolationSteps[i-1], stepSize, extrapolationDirection,leftFieldMap) : GeaneExtrapolation(extrapolationSteps[i-1], stepSize, extrapolationDirection, local_pos_err,leftFieldMap);
      
      mf::LogDebug(name_) << "test1 step.position = " << step.position << "\n";
      mf::LogDebug(name_) << "test1 step.momentum = " << step.momentum << "\n";

      if(leftFieldMap) {
	mf::LogInfo(name_) << "track left field map, returning false" << "\n";
	extrapolatedCaloHit.failureMode = "leftFieldMap";
	success = false;
	return success;
      }

      if(step.position.mag() == 0 && step.momentum.mag() == 0) {
	mf::LogDebug(name_) << "step.position.mag() = " << step.position.mag() << " and step.momentum.mag() = " << step.momentum.mag() << " so return false." << "\n";
	extrapolatedCaloHit.failureMode = "returnedFromRKWithEmptyVectors";
	success = false;
	return success;
      }
      
      //step.time = updateTime(extrapolationSteps[i-1],stepSize,extrapolationDirection);

      // errors for step position
      position_err[0] = std::sqrt(position_err[0]*position_err[0] + local_pos_err[0]*local_pos_err[0]);
      position_err[1] = std::sqrt(position_err[0]*position_err[0] + local_pos_err[0]*local_pos_err[1]);
      position_err[2] = std::sqrt(position_err[0]*position_err[0] + local_pos_err[0]*local_pos_err[2]);
      
      // initialize momentum and position for step - will compare momentum of next step to this one when calculating tangent point
      gm2geom::CoordSystem3Vector eStepMom_before;
      gm2geom::CoordSystem3Vector eStepPos_before;
      //      double eStepTime_before;
      
      // extrapolate one step
      extrapolationSteps.push_back(step);
      
      bool useCaloCoords = true; // TODO in future only use the calo coords here, but until we can drop a new coord system store in, use the tracker station coords instead if running over files that do not have calo coord systems loaded.
      
      std::string coordSysName = "";
      double endZpos = 0.0;

      mf::LogDebug(name_) << "step.position = " << step.position << "\n";
      mf::LogDebug(name_) << "step.momentum = " << step.momentum << "\n";
      mf::LogDebug(name_) << "step.position.transform(cs_,CalorimeterNumber[18]) = " << step.position.transform(cs_,"CalorimeterNumber[18]") << "\n";
      mf::LogDebug(name_) << "step.momentum.transform(cs_,CalorimeterNumber[18],true) = " << step.momentum.transform(cs_,"CalorimeterNumber[18]",true) << "\n";
  
      if (useCaloCoords) {
	coordSysName = Form("CalorimeterNumber[%02d]", station);
	endZpos = targetZpos;
      }
      
      else {
	coordSysName = Form("TrackerStation[%02d]",station);
	endZpos = -1094.6;
      }
      
      auto stepPosCalo = step.position.transform(cs_,coordSysName);
      auto stepMomCalo = step.momentum.transform(cs_,coordSysName,true);

      mf::LogDebug(name_) << "endZpos = " << endZpos << "\n";
      mf::LogDebug(name_) << "stepPosCalo = " << stepPosCalo << "\n";
      
      if (stepPosCalo.z() == endZpos) { 
	if (keepSteps) extrapolatedCaloHit.steps = extrapolationSteps;
	extrapolatedCaloHit.position = step.position;
	extrapolatedCaloHit.momentum = step.momentum;
	mf::LogDebug(name_) << "extrapolatedCaloHit.position = " << extrapolatedCaloHit.position << "\n";
	mf::LogDebug(name_) << "extrapolatedCaloHit.position.transform(cs_,coordSysName) = " << extrapolatedCaloHit.position.transform(cs_,coordSysName) << "\n";
	return success;
      }
      
      else if ( (useCaloCoords && (stepPosCalo.z() > endZpos)) || (!useCaloCoords && stepPosCalo.z() > endZpos )) {

	extrapolatedCaloHit.position = step.position;
	extrapolatedCaloHit.momentum = step.momentum;

	mf::LogDebug(name_) << "extrapolatedCaloHit.position.transform(cs_,coordSysName) = " << extrapolatedCaloHit.position.transform(cs_,coordSysName) << "\n";
	
	if (keepSteps) extrapolatedCaloHit.steps = extrapolationSteps;
	
	// code below is for extrapolating a partial step to get the position and momentum at the target z plane
	// rather than close to it
	
	ExtrapolationStep prevStep = extrapolationSteps[i-1];
	CoordSystem3Vector prevStepPos = prevStep.position;
	CoordSystem3Vector prevStepPosCalo = prevStepPos.transform(cs_,coordSysName);
	
	double dist1 = fabs(prevStepPosCalo.z() - endZpos);
	double dist2 = fabs(endZpos - stepPosCalo.z());
	double frac = dist1 / (dist1 + dist2);
	double stepSizeFrac = stepSize * frac;
	
	ExtrapolationStep step_corrected;
	step_corrected = oldCode ? oldGeaneExtrapolation(prevStep,stepSizeFrac,extrapolationDirection,leftFieldMap) : GeaneExtrapolation(prevStep,stepSizeFrac,extrapolationDirection,local_pos_err,leftFieldMap);
	//step_corrected.time = updateTime(prevStep,stepSizeFrac,extrapolationDirection);

	G4ThreeVector pos_corr = step_corrected.position.getVector();
	G4ThreeVector mom_corr = step_corrected.momentum.getVector();
	
	extrapolatedCaloHit.position = step_corrected.position;
	extrapolatedCaloHit.momentum = step_corrected.momentum;
	//extrapolatedCaloHit.time = step_corrected.time;

	if(leftFieldMap) {
	  mf::LogInfo(name_) << "track left field map, returning false" << "\n";
	  extrapolatedCaloHit.failureMode = "leftFieldMap";
	  success = false;
	  return success;
	}
		
	if (extrapolatedCaloHit.position.coordSystemName == "NULL" || extrapolatedCaloHit.momentum.coordSystemName == "NULL") {
	  mf::LogDebug(name_) << "extrapolatedCaloHit.position = " << extrapolatedCaloHit.position << "\n";
	  mf::LogDebug(name_) << "extrapolatedCaloHit.momentum = " << extrapolatedCaloHit.momentum << "\n";
	  //mf::LogDebug(name_) << "extrapolatedCaloHit.time = " << extrapolatedCaloHit.time << "\n";
	  mf::LogDebug(name_) << "returning false as decay vertex was not correctly filled" << "\n";
	  extrapolatedCaloHit.failureMode = "vertexVectorsNotFilled";
	  success = false;
	}
	
	else {
	  mf::LogInfo(name_) << "Returning from utils::extrapolateToCalorimeter" << "\n";
	  success = true;
	  return success;
	}
	
      }
      else {
	continue;
      }
    }
  }

  if (extrapolatedCaloHit.position.coordSystemName == "NULL" || extrapolatedCaloHit.momentum.coordSystemName == "NULL") {
    mf::LogInfo(name_) << "extrapolatedCaloHit.position = " << extrapolatedCaloHit.position << "\n";
    mf::LogInfo(name_) << "extrapolatedCaloHit.momentum = " << extrapolatedCaloHit.momentum << "\n";
    //mf::LogInfo(name_) << "extrapolatedCaloHit.time = " << extrapolatedCaloHit.time << "\n";
    mf::LogInfo(name_) << "returning false as decay vertex was not correctly filled" << "\n";
    extrapolatedCaloHit.failureMode = "vertexVectorsNotFilled";
    success = false;
  }
  
  mf::LogInfo(name_) << "Returning from GeaneExtrapolationUtils::ExtrapolateToCalorimeter with success = " << success << "\n";
  return success;
}

bool GeaneExtrapolationUtils::extrapolateToTrueAzimuth( const ExtrapolationStep &startPoint, DecayVertexArtRecord& decayVertex, const GeaneFhiclOptions& fclOptions)
{
  // function for extrapolating to the true decay azimuth rather than the tangent point. only in simulation.
  bool success = true;
  bool leftFieldMap = false;

  std::string direction = "backwards"; // must be backwards for true azimuth
  
  double stepSize = fclOptions.stepSize;
  
  // vector of extrapolation steps
  std::vector<ExtrapolationStep> extrapolationSteps;
  
  gm2geom::CoordSystem3Vector startPosition = startPoint.position;
  gm2geom::CoordSystem3Vector startMomentum = startPoint.momentum;
  
  // do one step forwards to get the correct momentum for the first backwards step
  ExtrapolationStep firstStep;
  ExtrapolationStep startStep; 

  double local_pos_err[3] = {0.0,0.0,0.0};
 
  // extrapolate forwards one step
  firstStep = GeaneExtrapolation(startPoint,stepSize,"forwards",local_pos_err,leftFieldMap);
  
  if(firstStep.position.mag() == 0 && firstStep.momentum.mag() == 0) {
    success = false;
    return success;
  }
  
  // extrapolate backwards one step to get the correct momentum for the start point
  startStep.momentum = -firstStep.momentum;
  startStep.position = startPosition;

  // need to transform here as plane could be defined in any coords - just need to be in world for individual steps
  double extrapolationDistance = fclOptions.totalDistance;
  int nsteps = extrapolationDistance / stepSize;

  // extrapolate one step
  extrapolationSteps.push_back(startStep);
  
  for (int i(1); i<nsteps+1; i++) {
    
    // extrapolate one step
    gm2strawtracker::ExtrapolationStep step = GeaneExtrapolation(extrapolationSteps[i-1], stepSize, direction, local_pos_err,leftFieldMap);
    if(step.position.mag() == 0 && step.momentum.mag() == 0) {
      success = false;
      return success;
    }
    
    G4ThreeVector truePosVect = decayVertex.mcVertex.position;
    G4ThreeVector stepPosVect = {step.position.x(),step.position.y(),step.position.z()};

    double trueAzimuth = gm2consts_->ComputeTheta(&truePosVect);
    double stepAzimuth = gm2consts_->ComputeTheta(&stepPosVect);

    if ( stepAzimuth == trueAzimuth || stepAzimuth < trueAzimuth) {
      decayVertex.position = step.position;
      decayVertex.momentum = step.momentum;
      if (decayVertex.position.coordSystemName == "NULL" || decayVertex.momentum.coordSystemName == "NULL") {
	success = false;
      }
      return success;
    }

    extrapolationSteps.push_back(step);
  }// end loop over steps

  decayVertex.position = extrapolationSteps.back().position;
  decayVertex.momentum = extrapolationSteps.back().momentum;

  if (decayVertex.position.coordSystemName == "NULL" || decayVertex.momentum.coordSystemName == "NULL") {
    success = false;
  }

  return success;
}

bool GeaneExtrapolationUtils::extrapolateToPlane( const ExtrapolationStep &startPoint, const double targetZ, std::string coordSysName, ExtrapolationStep &endPoint, std::string direction, double stepSize)
{
  bool success = true;
  bool leftFieldMap = false;
  
  // vector of extrapolation steps
  std::vector<ExtrapolationStep> extrapolationSteps;
  
  // don't need to transform here - RKExtrapolation sets pos/mom to world coords
  gm2geom::CoordSystem3Vector startPosition = startPoint.position;//.transform(cs_,"world");
  gm2geom::CoordSystem3Vector startMomentum = startPoint.momentum;//.transform(cs_,"world",true);
  
  // do one step forwards to get the correct momentum for the first backwards step
  ExtrapolationStep firstStep;
  ExtrapolationStep startStep;
  
  double local_pos_err[3] = {0.0,0.0,0.0};
  
  // extrapolate forwards one step
  if (direction == "backwards") {
    firstStep = GeaneExtrapolation(startPoint,stepSize,"forwards",local_pos_err,leftFieldMap);
  }

  else if (direction == "forwards") {
    gm2strawtracker::ExtrapolationStep backStep;
    backStep.position = startPoint.position;
    backStep.momentum = -startPoint.momentum;
    firstStep = GeaneExtrapolation(backStep,stepSize,"backwards",local_pos_err,leftFieldMap);
  }
  
  if(firstStep.position.mag() == 0 && firstStep.momentum.mag() == 0) {
    success = false;
    return success;
  }
  
  // extrapolate backwards one step to get the correct momentum for the start point
  startStep.momentum = -firstStep.momentum;
  startStep.position = startPosition;

  // need to transform here as plane could be defined in any coords - just need to be in world for individual steps
  double extrapolationDistance = fabs(targetZ - startPoint.position.transform(cs_,coordSysName).z());
  int nsteps = extrapolationDistance / stepSize;

  // extrapolate one step
  extrapolationSteps.push_back(startStep);
  
  for (int i(1); i<nsteps+1; i++) {
    
    // extrapolate one step
    gm2strawtracker::ExtrapolationStep step = GeaneExtrapolation(extrapolationSteps[i-1], stepSize, direction, local_pos_err,leftFieldMap);
    if(step.position.mag() == 0 && step.momentum.mag() == 0) {
      success = false;
      return success;
    }

    // RKExtrapolation function returns pos/mom in world coord system.
    // Transform into specified coord system for checking hit plane.
    CoordSystem3Vector stepPos = step.position.transform(cs_,coordSysName);
    CoordSystem3Vector stepMom = step.momentum.transform(cs_,coordSysName,true);
    
    if (stepPos.z() == targetZ) { // TODO get this from geometry not hard-coded
      endPoint.position = stepPos;
      endPoint.momentum = stepMom;
      //      debug << "in utils : endPoint.position = " << endPoint.position << "\n";
      if (endPoint.position.coordSystemName == "NULL" || endPoint.momentum.coordSystemName == "NULL") {
	success = false;
      }
      return success;
    }
    
    else if (stepPos.z() < targetZ) {
      // need to transform here as only transformed stepPos before
      endPoint.position = extrapolationSteps[i-1].position.transform(cs_,coordSysName);
      endPoint.momentum = extrapolationSteps[i-1].momentum.transform(cs_,coordSysName,true);
      if (endPoint.position.coordSystemName == "NULL" || endPoint.momentum.coordSystemName == "NULL") {
	success = false;
      }
      return success;
    }

    else {
      extrapolationSteps.push_back(step);
      continue;
    }

  }// end loop over steps

  mf::LogInfo(name_) << "At end of loop - never reached the target plane, returning false" << "\n";
  success = false;
  return success;
}

gm2strawtracker::ExtrapolationStep GeaneExtrapolationUtils::oldGeaneExtrapolation(gm2strawtracker::ExtrapolationStep step, double stepSize, std::string extrapolationDirection, bool leftFieldMap, bool forceUniformField) {

  //------------------------------------------------
  // note that Barry's algorithm sets the
  // constants using units:
  // cm, kGauss, Gev/c
  //------------------------------------------------
      
  //initialize
  gm2strawtracker::ExtrapolationStep eStep;

  // get the charge of the particle that made the track
  int Charge = 1;
  
  // magnetic field: 1.5 T = 15 kG
  double Bfield[3] = {0.0, 0.0, 0.0};
  double Bfield_new[3] = {0.0, 0.0, 0.0};
  
  // need to pass in initial position, direction cosines of track and magnitude of momentum
  // -- if extrapolationDirection is backwards, these should be negative
  gm2geom::CoordSystem3Vector stepPos_world = step.position.transform(cs_, "world");
  gm2geom::CoordSystem3Vector stepMom_world = step.momentum.transform(cs_, "world");
  
  double dirCosine_x = stepMom_world.x()  / stepMom_world.mag();
  double dirCosine_y = stepMom_world.y()  / stepMom_world.mag();
  double dirCosine_z = stepMom_world.z()  / stepMom_world.mag();
  
  double pTot = stepMom_world.mag(); 
  
  double mm2cm = 0.1;
  double cm2mm = 10.0;
  double mev2gev = 0.001; 
  double gev2mev = 1000.0; 
  double convert2kgauss = 10000.0;// going from Geant B field to kGauss, appears that geant is in mT??

  pTot *= mev2gev;
  
  double Vin[7] = {double(stepPos_world.x() * mm2cm),
		   double(stepPos_world.y() * mm2cm),
		   double(stepPos_world.z() * mm2cm), 
		   double(dirCosine_x), 
		   double(dirCosine_y), 
		   double(dirCosine_z), 
		   double(pTot)};
  
  // this will be the 'track vector' output from the algorithm
  double Vout[7];
  
  // B field vector
  double BFieldVect[4];

  // temporary momentum and position arrays
  double xyzt[3];
  double xyz[3];
  
  // arrays that store the sec of the angle of the track in each direction
  double Secxs[4];
  double Secys[4];
  double Seczs[4];
  
  // constants in units cm, GeV/c and kGauss
  int Maxit = 1992;
  int Maxcut = 11;
  double Ec = 2.9979251e-4; // if this is in cm/s it should be e-6
  double Dlt = 1.0e-4;
  double Dlt32 = Dlt/32.0;
  
  double Pisqua = 0.986960440109e+01;
  
  double zzz = 0.0;
  double xxx = 0.0;
  double yyy = 0.0;

  // initialize distance elements
  double Dxt = 0.0;
  double Dyt = 0.0;
  double Dzt = 0.0;
  
  int Iter = 0;
  int Ncut = 0;
  
  // initially, set the final position/momentum equal to the initial
  for(int j=1; j<=7; j++) {
      Vout[j-1] = Vin[j-1];
  } 
  
  double Pinv   = Ec * Charge / Vin[6];

  // AUG 25
  if (extrapolationDirection == "backwards") Pinv *= -1; 
  
  // tolerance constant
  double Tl = 0.0;
  double H  = stepSize;
  
  int Iloop = 0;
  do {
    Iloop = 1;
    int Icon1 = 1;
    int Icon2 = 1;
    int Icon3 = 1;
    
    // H can only be so different from the input step size
    // check that it is within tolerance Tl of step size
    // if not, set equal to the difference 
    double Rest  = stepSize - Tl;
    if(fabs(H) > fabs(Rest)) {
      H = Rest;
    }

    xxx = Vout[0];
    yyy = Vout[1];
    zzz = Vout[2];    
    
    const G4double trackPos2[4] = {double(xxx*cm2mm), double(yyy*cm2mm), double(zzz*cm2mm), 0.0};
    fieldManager_->GetFieldValue(trackPos2, Bfield); 
    
    if (std::isnan(Bfield[0]) || std::isnan(Bfield[1]) || std::isnan(Bfield[2])) {
      mf::LogInfo(name_) << "Bfield is nan --> have left field map region. Returning false from RKExtrapolation" << "\n";
      leftFieldMap = true;
      eStep.position.set(Vout[0] , Vout[1] , Vout[2] , "world");
      eStep.momentum.set(Vout[3]*Vout[6], Vout[4]*Vout[6], Vout[5]*Vout[6], "world");
      return eStep;
    }

    if (forceUniformField) {
      Bfield[0] = 0.0;
      Bfield[1] = 0.00145127;
      Bfield[2] = 0.0;
    }

    BFieldVect[0] = Bfield[0]*convert2kgauss;
    BFieldVect[1] = Bfield[1]*convert2kgauss;
    BFieldVect[2] = Bfield[2]*convert2kgauss;

    // --------------------------------
    //    Start of integration
    // --------------------------------
    
    // give intermediate position xyz the values of the output vector
    xyz[0]  = Vout[0];
    xyz[1]  = Vout[1];
    xyz[2]  = Vout[2];
    
    // A, B, C are the x y z components of the momentum
    double A = Vout[3];
    double B = Vout[4];
    double C = Vout[5];
    
    // make the 'step size' H smaller, and decrease the size of Pinv accordingly
    double H2 = H/2.;
    double H4 = H/4.;
    double Ph = Pinv * H;
    double Ph2= Ph/2;
    
    // calculate the first guess of the angle of the track
    // this is effectively:
    //    momentum.cross(BFieldVect) * Pinv * H/2
    // so the cross product of the momentum and b field, divided by the magnitude of the momentum, for that step
    
    Secxs[0] = (B * BFieldVect[2] - C * BFieldVect[1]) * Ph2;
    Secys[0] = (C * BFieldVect[0] - A * BFieldVect[2]) * Ph2;
    Seczs[0] = (A * BFieldVect[1] - B * BFieldVect[0]) * Ph2;
    
    // angle squared
    double Ang2 = (Secxs[0]*Secxs[0] + Secys[0]*Secys[0] + Seczs[0]*Seczs[0]);
    
    if(Ang2 > Pisqua) {
      eStep.position.set(0.0,0.0,0.0,"");
      eStep.momentum.set(0.0,0.0,0.0,"");
      return eStep;
    }
    
    // small distance elements to move the track along a little bit
    // H2 * A etc are the momentum components along half the step size
    // then add a smaller bit of the step (H/4) in direction defined by the angles
    // Secxs[0] etc
    Dxt    = H2 * A + H4 * Secxs[0];
    Dyt    = H2 * B + H4 * Secys[0];
    Dzt    = H2 * C + H4 * Seczs[0];
    
    // update the positions by adding this small disance element to the current
    // intermediate position
    xyzt[0]= xyz[0] + Dxt;
    xyzt[1]= xyz[1] + Dyt;
    xyzt[2]= xyz[2] + Dzt;
    
    // ------------------------------------
    // Second intermediate point
    // ------------------------------------
    
    // make an estimate of how far the track has moved in the previous step
    double Est = fabs(Dxt)+ fabs(Dyt)+ fabs(Dzt);
    
    // if this is bigger than the step size, either use a helix if we have performed too many cuts
    // already, or just set H = H/2 
    if(Est > H) {
      Ncut = Ncut + 1;
      if(Ncut > Maxcut) {
	eStep.position.set(0.0,0.0,0.0,"");
	eStep.momentum.set(0.0,0.0,0.0,"");
	return eStep;
      }
      H = H/2;
      Icon1 = 0;
    }
    // this will be true if Est was not > H (gets set to 0 inside that if statement)
    if(Icon1 == 1) {
      
      xxx = xyzt[0];
      yyy = xyzt[1];
      zzz = xyzt[2];    

      const G4double trackPos3[4] = {double(xxx*cm2mm), double(yyy*cm2mm), double(zzz*cm2mm), 0.0};
      fieldManager_->GetFieldValue(trackPos3, Bfield); 

      if (std::isnan(Bfield[0]) || std::isnan(Bfield[1]) || std::isnan(Bfield[2])) {
	mf::LogInfo(name_) << "Bfield is nan --> have left field map region. Returning false from RKExtrapolation" << "\n";
	leftFieldMap = true;
	eStep.position.set(Vout[0] , Vout[1] , Vout[2] , "world");
	eStep.momentum.set(Vout[3]*Vout[6], Vout[4]*Vout[6], Vout[5]*Vout[6], "world");
	return eStep;
      }

      if (forceUniformField) {
	Bfield[0] = 0.0;
	Bfield[1] = 0.00145127;
	Bfield[2] = 0.0;
      }
      
      BFieldVect[0] = Bfield[0]*convert2kgauss;
      BFieldVect[1] = Bfield[1]*convert2kgauss;
      BFieldVect[2] = Bfield[2]*convert2kgauss;
      
      //--std::cout << "\tBFieldVect[0] = " << BFieldVect[0] << "\tBFieldVect[1] = " << BFieldVect[1] << "\tBFieldVect[2] = " << BFieldVect[2] << "\n\n";
      
      // update the momentum components
      double At     = A + Secxs[0];
      double Bt     = B + Secys[0];
      double Ct     = C + Seczs[0];
      
      //--std::cout << "\tAt = " << At << "\tBt = " << Bt << "\tCt = " << Ct << "\n\n";
      
      // now update the next angles
      // this is 
      //        (momentum + secxyz).cross(BFieldVect) * (Pinv * H / 2)
      Secxs[1] = (Bt * BFieldVect[2] - Ct * BFieldVect[1]) * Ph2;
      Secys[1] = (Ct * BFieldVect[0] - At * BFieldVect[2]) * Ph2;
      Seczs[1] = (At * BFieldVect[1] - Bt * BFieldVect[0]) * Ph2;
      
      //--std::cout << "\tSecxs[1] = " << Secxs[1] << "\tSecys[1] = " << Secys[1] << "\tSeczs[1] = " << Seczs[1] << "\n\n";
      
      // add the new angle to the momentum
      At     = A + Secxs[1];
      Bt     = B + Secys[1];
      Ct     = C + Seczs[1];
      
      //--std::cout << "\tAt = " << At << "\tBt = " << Bt << "\tCt = " << Ct << "\n\n";
      
      // update the angles again
      // using the new momentum components
      Secxs[2] = (Bt * BFieldVect[2] - Ct * BFieldVect[1]) * Ph2;
      Secys[2] = (Ct * BFieldVect[0] - At * BFieldVect[2]) * Ph2;
      Seczs[2] = (At * BFieldVect[1] - Bt * BFieldVect[0]) * Ph2;
      
      //--std::cout << "\tSecxs[2] = " << Secxs[2] << "\tSecys[2] = " << Secys[2] << "\tSeczs[2] = " << Seczs[2] << "\n\n";
      
      // set the distance to extrapolate along
      Dxt    = H * (A + Secxs[2]);
      Dyt    = H * (B + Secys[2]);
      Dzt    = H * (C + Seczs[2]);
      
      //--std::cout << "\tDxt = " << Dxt << "\tDyt = " << Dyt << "\tDzt = " << Dzt << "\n\n";
      
      // update the intermediate position by adding the distance element
      xyzt[0]  = xyz[0] + Dxt;
      xyzt[1]  = xyz[1] + Dyt;
      xyzt[2]  = xyz[2] + Dzt;
      
      //--std::cout << "\txyzt[0] = " << xyzt[0] << "\txyzt[1] = " << xyzt[1] << "\txyzt[2] = " << xyzt[2] << "\n";
      
      // update the momentum
      At     = A + 2*Secxs[2];
      Bt     = B + 2*Secys[2];
      Ct     = C + 2*Seczs[2];
      
      //--std::cout << "\tAt = " << At << "\tBt = " << Bt << "\tCt = " << Ct << "\n\n";
      
      // make a new estimate of how far we have got
      // using the new distance elements
      Est = fabs(Dxt)+ fabs(Dyt) + fabs(Dzt);
      
      //--std::cout << "\tEst = " << Est << "\n\n";
      
      // check if this is more than twice the step size, if so either use helix
      // or set the step size H = H/2
      // basically convergence check
      if(Est > (2.0*fabs(H))) {
	//--std::cout << "\tEst = " << Est << " which is > " << 2.0*fabs(H) << "\n\n";
	Ncut = Ncut + 1;
	if(Ncut > Maxcut) {

	  eStep.position.set(0.0,0.0,0.0,"world");
	  eStep.momentum.set(0.0,0.0,0.0,"world");
	  return eStep;
	}
	H = H/2;
	Icon2 = 0;
      }
      
      if(Icon2 == 1) {

	xxx = xyzt[0];
	yyy = xyzt[1];
	zzz = xyzt[2];    

	const G4double trackPos4[4] = {double(xxx*cm2mm), double(yyy*cm2mm), double(zzz*cm2mm), 0.0};
	fieldManager_->GetFieldValue(trackPos4, Bfield); 

	if (std::isnan(Bfield[0]) || std::isnan(Bfield[1]) || std::isnan(Bfield[2])) {
	  mf::LogInfo(name_) << "Bfield is nan --> have left field map region. Returning false from RKExtrapolation" << "\n";
	  leftFieldMap = true;
	  eStep.position.set(Vout[0] , Vout[1] , Vout[2] , "world");
	  eStep.momentum.set(Vout[3]*Vout[6], Vout[4]*Vout[6], Vout[5]*Vout[6], "world");
	  return eStep;
	}
	
	if (forceUniformField) {
	  Bfield[0] = 0.0;
	  Bfield[1] = 0.00145127;
	  Bfield[2] = 0.0;
	}
      
	BFieldVect[0] = Bfield[0]*convert2kgauss;
	BFieldVect[1] = Bfield[1]*convert2kgauss;
	BFieldVect[2] = Bfield[2]*convert2kgauss;
	
	//--std::cout << "\tBFieldVect[0] = " << BFieldVect[0] << "\tBFieldVect[1] = " << BFieldVect[1] << "\tBFieldVect[2] = " << BFieldVect[2] << "\n\n";
	
	if (sqrt(Bfield[0]*Bfield[0] + Bfield[1]*Bfield[1] + Bfield[2]*Bfield[2]) == 0){
	  fieldManager_->GetFieldValue(trackPos2, Bfield_new); 
	  
	  if (std::isnan(Bfield[0]) || std::isnan(Bfield[1]) || std::isnan(Bfield[2])) {
	    mf::LogInfo(name_) << "Bfield is nan --> have left field map region. Returning false from RKExtrapolation" << "\n";
	    leftFieldMap = true;
	    eStep.position.set(Vout[0] , Vout[1] , Vout[2] , "world");
	    eStep.momentum.set(Vout[3]*Vout[6], Vout[4]*Vout[6], Vout[5]*Vout[6], "world");
	    return eStep;
	  }
	  
	  if (forceUniformField) {
	    Bfield[0] = 0.0;
	    Bfield[1] = 0.00145127;
	    Bfield[2] = 0.0;
	  }

	  BFieldVect[0] = Bfield_new[0]*convert2kgauss;
	  BFieldVect[1] = Bfield_new[1]*convert2kgauss;
	  BFieldVect[2] = Bfield_new[2]*convert2kgauss;
	  
	  if (sqrt(Bfield_new[0]*Bfield_new[0] + Bfield_new[1]*Bfield_new[1] + Bfield_new[2]*Bfield_new[2]) == 0){
	    //--std::cout << "\tGone beyond mapped storage ring field region, setting to (0, 0, 15) kGauss." << "\n";
	    BFieldVect[0] = 0.0;
	    BFieldVect[1] = 14.5127;
	    BFieldVect[2] = 0.0;
	  } 
	}
	
	xyz[2] = xyz[2] + (C + (Seczs[0] + Seczs[1] + Seczs[2]) * (1/3)) * H;
	xyz[1] = xyz[1] + (B + (Secys[0] + Secys[1] + Secys[2]) * (1/3)) * H;
	xyz[0] = xyz[0] + (A + (Secxs[0] + Secxs[1] + Secxs[2]) * (1/3)) * H;
	
	//--std::cout << "\txyz[2] = " << xyz[2] << "\txyz[1] = " << xyz[1] << "\txyz[0] = " << xyz[0] << "\n\n";
	
	Secxs[3] = (Bt*BFieldVect[2] - Ct*BFieldVect[1])* Ph2;
	Secys[3] = (Ct*BFieldVect[0] - At*BFieldVect[2])* Ph2;
	Seczs[3] = (At*BFieldVect[1] - Bt*BFieldVect[0])* Ph2;
	
	A      = A + (Secxs[0]+Secxs[3]+ 2*(Secxs[1]+Secxs[2]))*(1.0/3.0);
	B      = B + (Secys[0]+Secys[3]+ 2*(Secys[1]+Secys[2]))*(1.0/3.0);
	C      = C + (Seczs[0]+Seczs[3]+ 2*(Seczs[1]+Seczs[2]))*(1.0/3.0);
	
	Est = fabs(Secxs[0]+Secxs[3] - (Secxs[1]+Secxs[2]))
	  + fabs(Secys[0]+Secys[3] - (Secys[1]+Secys[2]))
	  + fabs(Seczs[0]+Seczs[3] - (Seczs[1]+Seczs[2]));
	
	if(Est > Dlt && fabs(H) > 1.E-4) {
	  Ncut = Ncut + 1;
	  if(Ncut > Maxcut) {

	    eStep.position.set(0.0,0.0,0.0,"");
	    eStep.momentum.set(0.0,0.0,0.0,"");
	    return eStep;
	    
	  }
	  H = H/2;
	  Icon3 = 0; 
	}
	
	if(Icon3 == 1) {
	  Iter = Iter + 1;
	  Ncut = 0;
	  if(Iter > Maxit) {
	    eStep.position.set(0.0,0.0,0.0,"");
	    eStep.momentum.set(0.0,0.0,0.0,"");
	    return eStep;
	  }
	  Tl = Tl + H;
	  if(Est < Dlt32) {
	    H = 2*H;
	  }
	  
	  double Cba = 1/sqrt(A*A + B*B + C*C);
	  Vout[0] = xyz[0];
	  Vout[1] = xyz[1];
	  Vout[2] = xyz[2];
	  Vout[3] = Cba*A;
	  Vout[4] = Cba*B;
	  Vout[5] = Cba*C;
	  	  
	  eStep.position.set(Vout[0] , Vout[1] , Vout[2] , "world");
	  eStep.momentum.set(Vout[3]*Vout[6], Vout[4]*Vout[6], Vout[5]*Vout[6], "world");
	  
	  Rest = stepSize - Tl;
	  if(stepSize < 0.0) {
	    Rest = -Rest;
	  }
	  if(Rest < 1.0e-5*fabs(stepSize)) {
	    Iloop = 0;
	  }
	}                // Icon3
      }                  // Icon2
    }                    // Icon1 
  } while (Iloop == 1);
  
  //convert back to Mev
  eStep.momentum.set( eStep.momentum.x()*gev2mev, eStep.momentum.y()*gev2mev, eStep.momentum.z()*gev2mev );
  
  //convert back to mm
  eStep.position.set( eStep.position.x()*cm2mm, eStep.position.y()*cm2mm, eStep.position.z()*cm2mm );

  return eStep;
}


gm2strawtracker::ExtrapolationStep GeaneExtrapolationUtils::GeaneExtrapolation(gm2strawtracker::ExtrapolationStep step, double stepSize, std::string extrapolationDirection, double pos_err[3], bool leftFieldMap, bool forceUniformField)
{

  //------------------------------------------------
  // this algorithm sets the
  // constants using units:
  // cm, kGauss, Gev/c
  //------------------------------------------------
  
  //  std::cout << "Begin RKExtrapolation" << "\n";

  //initialize
  gm2strawtracker::ExtrapolationStep eStep;
  
  // get the charge of the particle that made the track
  int Charge = 1; // e+
  
  // magnetic field: 1.5 T = 15 kG
  double Bfield[3] = {0.0, 0.0, 0.0};
  
  // initial pos and mom of step - must be world coords
  gm2geom::CoordSystem3Vector stepPos_world = step.position.transform(cs_,"world");
  gm2geom::CoordSystem3Vector stepMom_world = step.momentum.transform(cs_,"world",true);
  
  // step direction cosines
  double dirCosine_x = stepMom_world.x()  / stepMom_world.mag();
  double dirCosine_y = stepMom_world.y()  / stepMom_world.mag();
  double dirCosine_z = stepMom_world.z()  / stepMom_world.mag();
  
  // magnitude of initial momentum
  double pTot = stepMom_world.mag(); 
  
  // conversions
  double mm2cm = 0.1;
  double cm2mm = 10.0;
  double mev2gev = 0.001; 
  double gev2mev = 1000.0; 
  double convert2kgauss = 10000.0;// going from Geant B field to kGauss, geant is in mT

  // convert mom magnitude to GeV for algorithm
  pTot *= mev2gev;

  // initial step vector - pos and mom components and mom magnitude
  double Vin[7] = {double(stepPos_world.x() * mm2cm),
		   double(stepPos_world.y() * mm2cm),
		   double(stepPos_world.z() * mm2cm), 
		   double(dirCosine_x), 
		   double(dirCosine_y), 
		   double(dirCosine_z), 
		   double(pTot)};
  
  // output track vector
  double Vout[7];
  
  // B Field vector
  double BFieldVect[4];
  
  // forces Bfield to be {0, 1.45127, 0} - TODO make fhicl parameter
  //  bool forceUniformField = false;
  
  // temporary momentum and position arrays
  double intermediate_pos[3]; 
  double intermediate_mom[3];
  double track_pos[3];
  double track_mom[3];
  
  // arrays that store the updated angle in each direction
  double A[4];
  double B[4];
  double C[4];
  
  // constants in units cm, GeV/c and kGauss
  int Maxit    = 1992;
  int Maxcut   = 11;
  double Ec    = 2.9979251e-4; // if this is in cm/s it should be e-6, speed of light
  double Dlt   = 1.0e-4;
  double Dlt32 = Dlt/32.0;
  
  // initialize distance elements
  double Dxt = 0.0;
  double Dyt = 0.0;
  double Dzt = 0.0;
  
  // conversion checks
  int Iter = 0;
  int Ncut = 0;
  
  // initially, set the final position/momentum equal to the initial
  for(int j=1; j<=7; j++) Vout[j-1] = Vin[j-1];
  
  // inverse momentum
  double Pinv = Ec * Charge / Vin[6];
  
  // flip charge so that track bends correct way around B field if direction is backwards (towards decay vertex)
  if (extrapolationDirection == "backwards") Pinv *= -1; 
  
  // tolerance constant
  double Tl = 0.0;//fhicl_.get<double>("tolerance",0.0);
  
  bool toleranceTesting = false;
  if (toleranceTesting) Tl = 0.1 * stepSize;
  
  // step size that gets overwritten - if step does not converge within tolerance constant
  // make step size smaller
  double H = stepSize;

  int Iloop = 0;
  do {
    Iloop = 1;
    
    // convergence checks
    int Icon1 = 1; 
    int Icon2 = 1; 
    int Icon3 = 1; 
    
    // check if within tolerance of end of step
    double Rest  = H - Tl;
    if(fabs(H) > fabs(Rest)) {
      H = Rest;
    }
    
     // current pos = output pos
    track_pos[0]  = Vout[0];
    track_pos[1]  = Vout[1];
    track_pos[2]  = Vout[2];
    
    // set the momentum components
    track_mom[0] = Vout[3]; 
    track_mom[1] = Vout[4];
    track_mom[2] = Vout[5];
    
//    std::cout << "START" << "\n";
//    std::cout << "track_pos = [" << track_pos[0] << ", " << track_pos[1] << ", " << track_pos[2] << "]\n";
//    std::cout << "track_mom = [" << track_mom[0] << ", " << track_mom[1] << ", " << track_mom[2] << "]\n";
    
    // position for field lookup
    const G4double initialTrackPos[4] = {track_pos[0]*cm2mm,track_pos[1]*cm2mm,track_pos[2]*cm2mm,0.0};

    // initial field lookup
    fieldManager_->GetFieldValue(initialTrackPos, Bfield); 
    
    if (std::isnan(Bfield[0]) || std::isnan(Bfield[1]) || std::isnan(Bfield[2])) {
      mf::LogInfo(name_) << "Bfield is nan --> have left field map region. Returning false from RKExtrapolation" << "\n";
      leftFieldMap = true;
      eStep.position.set(Vout[0] , Vout[1] , Vout[2] , "world");
      eStep.momentum.set(Vout[3]*Vout[6], Vout[4]*Vout[6], Vout[5]*Vout[6], "world");
      return eStep;
    }
    
    if (forceUniformField) {
      Bfield[0] = 0.0;
      Bfield[1] = 0.00145127;
      Bfield[2] = 0.0;
    }
    
    // convert B field components to kgauss for algorithm - fieldManager returns in mT
    BFieldVect[0] = Bfield[0]*convert2kgauss;
    BFieldVect[1] = Bfield[1]*convert2kgauss;
    BFieldVect[2] = Bfield[2]*convert2kgauss;
    
    // --------------------------------
    //    Start of integration
    // --------------------------------
    
    // calculate the first guess of the angle of the track
    // this is effectively:
    //    momentum.cross(BFieldVect) * Pinv * H/2
    // so the cross product of the momentum and b field, divided by the magnitude of the momentum, for that step
    
    // these are the initial direction components - 'K' in RK Nystrom algorithm literature
    A[0] = (track_mom[1] * BFieldVect[2] - track_mom[2] * BFieldVect[1]) * ((Pinv*H)/2);
    B[0] = (track_mom[2] * BFieldVect[0] - track_mom[0] * BFieldVect[2]) * ((Pinv*H)/2);
    C[0] = (track_mom[0] * BFieldVect[1] - track_mom[1] * BFieldVect[0]) * ((Pinv*H)/2);

    double Ang2 = (A[0]*A[0] + B[0]*B[0] + C[0]*C[0]);
    
    // -------------------------------------
    // If angle is too big just use a helix
    // (TODO implement this, for now return)
    // -------------------------------------
    if(Ang2 > pow(TMath::Pi(),2)) {
       eStep.momentum.set(0.0,0.0,0.0,"world");
       eStep.position.set(0.0,0.0,0.0,"world");
       return eStep;
     }

    // update the track position using the direction components of the initial momentum
    // and the direction of the angle between the initial momentum and the B field at
    // the start position
    Dxt    = (H/2) * track_mom[0] + (H/4) * A[0];
    Dyt    = (H/2) * track_mom[1] + (H/4) * B[0];
    Dzt    = (H/2) * track_mom[2] + (H/4) * C[0];

    // update the positions by adding this small disance element to the current
    // intermediate position
    intermediate_pos[0]= track_pos[0] + Dxt;
    intermediate_pos[1]= track_pos[1] + Dyt;
    intermediate_pos[2]= track_pos[2] + Dzt;

    // ------------------------------------
    // Second intermediate point
    // ------------------------------------
    
    // make an estimate of how far the track has moved in the previous step
    double Est = fabs(Dxt)+ fabs(Dyt)+ fabs(Dzt);
    
    // if this is bigger than the step size, either use a helix if we have performed too many cuts
    // already, or just set H = H/2 
    if(Est > H) {
      Ncut = Ncut + 1;
      if(Ncut > Maxcut) {
	 eStep.position.set(0.0,0.0,0.0,"world");
	 eStep.momentum.set(0.0,0.0,0.0,"world");
	 return eStep;
      }
      H = H/2;
      // set to 0 to go back round loop
      Icon1 = 0;
    }

    if(Icon1 == 1) {

      // field lookup
      const G4double trackPos3[4] = {intermediate_pos[0]*cm2mm, intermediate_pos[1]*cm2mm, intermediate_pos[2]*cm2mm,0.0};
      fieldManager_->GetFieldValue(trackPos3, Bfield);

      if (std::isnan(Bfield[0]) || std::isnan(Bfield[1]) || std::isnan(Bfield[2])) {
	mf::LogInfo(name_) << "Bfield is nan --> have left field map region. Returning false from RKExtrapolation" << "\n";
	leftFieldMap = true;
	eStep.position.set(Vout[0] , Vout[1] , Vout[2] , "world");
	eStep.momentum.set(Vout[3]*Vout[6], Vout[4]*Vout[6], Vout[5]*Vout[6], "world");
	return eStep;
      }
      
      if (forceUniformField) {
	Bfield[0] = 0.0;
	Bfield[1] = 0.00145127;
	Bfield[2] = 0.0;
      }
      
      // convert to kgauss
      BFieldVect[0] = Bfield[0]*convert2kgauss;
      BFieldVect[1] = Bfield[1]*convert2kgauss;
      BFieldVect[2] = Bfield[2]*convert2kgauss;
      
      // update the momentum components
      intermediate_mom[0] = track_mom[0] + A[0]; // this is the momentum direction as calculated at START POS
      intermediate_mom[1] = track_mom[1] + B[0]; 
      intermediate_mom[2] = track_mom[2] + C[0];
      
      // now update the next angles
      // this is 
      //        (momentum + A).cross(BFieldVect) * (Pinv * H / 2)
      // update with angle from FIRST INTERMEDIATE STEP momentum (extra angle)
      A[1] = (intermediate_mom[1] * BFieldVect[2] - intermediate_mom[2] * BFieldVect[1]) * ((Pinv*H)/2); 
      B[1] = (intermediate_mom[2] * BFieldVect[0] - intermediate_mom[0] * BFieldVect[2]) * ((Pinv*H)/2);
      C[1] = (intermediate_mom[0] * BFieldVect[1] - intermediate_mom[1] * BFieldVect[0]) * ((Pinv*H)/2);
      
      // add the new angle to the momentum
      intermediate_mom[0] = track_mom[0] + A[1]; // THIS IS THE TOTAL MOMENTUM DIRECTION COMPONENT TAKING INTO ACCOUNT START POINT MOM AND INTERMEDIATE POINT MOM
      intermediate_mom[1] = track_mom[1] + B[1];
      intermediate_mom[2] = track_mom[2] + C[1];

      // update the angles again using the new momentum components
      A[2] = (intermediate_mom[1] * BFieldVect[2] - intermediate_mom[2] * BFieldVect[1]) * ((Pinv*H)/2);
      B[2] = (intermediate_mom[2] * BFieldVect[0] - intermediate_mom[0] * BFieldVect[2]) * ((Pinv*H)/2);
      C[2] = (intermediate_mom[0] * BFieldVect[1] - intermediate_mom[1] * BFieldVect[0]) * ((Pinv*H)/2);
      
      // set the distance to extrapolate along
      Dxt = H * (track_mom[0] + A[2]);
      Dyt = H * (track_mom[1] + B[2]);
      Dzt = H * (track_mom[2] + C[2]);
      
      // *** NOTE this point is still evaluated at H/2 *** 

      // update the intermediate position by adding the distance element to the current track position
      intermediate_pos[0]  = track_pos[0] + Dxt;
      intermediate_pos[1]  = track_pos[1] + Dyt;
      intermediate_pos[2]  = track_pos[2] + Dzt;
      
      // update the momentum again using the latest angles
      intermediate_mom[0] = track_mom[0] + 2*A[2];
      intermediate_mom[1] = track_mom[1] + 2*B[2];
      intermediate_mom[2] = track_mom[2] + 2*C[2];

      // make a new estimate of how far we have got using the new distance elements
      Est = fabs(Dxt)+ fabs(Dyt) + fabs(Dzt);

      // check if this is more than twice the step size, if so either use helix
      // or set the step size H = H/2
      // basically convergence check
      if(Est > (2.0*fabs(H))) {
	
	Ncut = Ncut + 1;
	
	if(Ncut > Maxcut) {
	  eStep.position.set(0.0,0.0,0.0,"world");
	  eStep.momentum.set(0.0,0.0,0.0,"world");
	  return eStep;
	}
	
	// decrease step size again
	H = H/2;
	
	// set to 0 to go back round loop
	Icon2 = 0;
	Iloop = 0;
      }
      
      // if returned successfully from previous loop
      if(Icon2 == 1) {

	// field lookup
	const G4double trackPos4[4] = {intermediate_pos[0]*cm2mm, intermediate_pos[1]*cm2mm, intermediate_pos[2]*cm2mm, 0.0};
	fieldManager_->GetFieldValue(trackPos4, Bfield);
	
	if (std::isnan(Bfield[0]) || std::isnan(Bfield[1]) || std::isnan(Bfield[2])) {
	  mf::LogInfo(name_) << "Bfield is nan --> have left field map region. Returning false from RKExtrapolation" << "\n";
	  leftFieldMap = true;
	  eStep.position.set(Vout[0] , Vout[1] , Vout[2] , "world");
	  eStep.momentum.set(Vout[3]*Vout[6], Vout[4]*Vout[6], Vout[5]*Vout[6], "world");
	  return eStep;
	}
	
	if (forceUniformField) {
	  Bfield[0] = 0.0;
	  Bfield[1] = 0.00145127;
	  Bfield[2] = 0.0;
	}
	
	// convert to kgauss
	BFieldVect[0] = Bfield[0]*convert2kgauss;
	BFieldVect[1] = Bfield[1]*convert2kgauss;
	BFieldVect[2] = Bfield[2]*convert2kgauss;
	
	// update the track position with the latest angles -- Fn in literature
	track_pos[0] = track_pos[0] + (track_mom[0] + (A[0] + A[1] + A[2]) * (1/3)) * H;
	track_pos[1] = track_pos[1] + (track_mom[1] + (B[0] + B[1] + B[2]) * (1/3)) * H;
	track_pos[2] = track_pos[2] + (track_mom[2] + (C[0] + C[1] + C[2]) * (1/3)) * H;

	// update the angles again 
	A[3] = (intermediate_mom[1]*BFieldVect[2] - intermediate_mom[2]*BFieldVect[1])* ((Pinv*H)/2);
	B[3] = (intermediate_mom[2]*BFieldVect[0] - intermediate_mom[0]*BFieldVect[2])* ((Pinv*H)/2);
	C[3] = (intermediate_mom[0]*BFieldVect[1] - intermediate_mom[1]*BFieldVect[0])* ((Pinv*H)/2);

	// update the momentum with the latest angles
	// THIS IS Gn IN LITERATURE
	track_mom[0] = track_mom[0] + (A[0]+A[3]+ 2*(A[1]+A[2]))*(1.0/3.0);
	track_mom[1] = track_mom[1] + (B[0]+B[3]+ 2*(B[1]+B[2]))*(1.0/3.0);
	track_mom[2] = track_mom[2] + (C[0]+C[3]+ 2*(C[1]+C[2]))*(1.0/3.0);

	// distance travelled along track step
	Est = (fabs(A[0]+A[3] - (A[1]+A[2])) +
	       fabs(B[0]+B[3] - (B[1]+B[2])) +
	       fabs(C[0]+C[3] - (C[1]+C[2])) );
	
	// convergence check
	if(Est > Dlt && fabs(H) > 1.e-4) {
	  Ncut = Ncut + 1;
	  
	  // use helix
	  if(Ncut > Maxcut) {
	    eStep.position.set(0.0,0.0,0.0,"world");
	    eStep.momentum.set(0.0,0.0,0.0,"world");
	    return eStep;
	  }
	  
	  // decrease step size and go around again
	  H = H/2;
	  Icon3 = 0; 
	}
	
	if(Icon3 == 1) {
	  // count number of iteration
	  Iter = Iter + 1;
	  Ncut = 0;
	  // if too many iterations, use helix
	  if(Iter > Maxit) {
	    eStep.position.set(0.0,0.0,0.0,"world");
	    eStep.momentum.set(0.0,0.0,0.0,"world");
	    return eStep;
	  }
	  
	  // adjust tolerance constant using current step size
	  Tl = Tl + H;
	  // check if distance along step is less than Dlt32, if it is increase step size
	  if(Est < Dlt32) {
	    H = 2*H;
	  }
	  
	  // 1/P
	  // double track_mom_mag_inv = 1/stepMom_world.mag();
	  
	  // fill the output track vector
	  Vout[0] = track_pos[0];
	  Vout[1] = track_pos[1];
	  Vout[2] = track_pos[2];
	  Vout[3] = track_mom[0];
	  Vout[4] = track_mom[1];
	  Vout[5] = track_mom[2];
	  
	  eStep.position.set(Vout[0],Vout[1],Vout[2],"world");
	  eStep.momentum.set(Vout[3]*Vout[6],Vout[4]*Vout[6],Vout[5]*Vout[6],"world"); 

	  ///
	  /// now update the Jacobian with the derivatives for the step
	  ///
	  /*
	  // unit
	  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(7,7);
	  J(0,0) = 1;
	  J(1,1) = 1;
	  J(2,2) = 1;
	  J(6,6) = Pinv;
	  
	  // TODO get this per track from GeaneArtRecord
	  Eigen::MatrixXd covMatrix = Eigen::MatrixXd::Zero(5,5);
	  covMatrix(0,0) = 1.68124e-05;
	  covMatrix.row(1) << 0, 1.93178e-06, 1.59488e-06, -0.000134473, -9.90522e-05;
	  covMatrix.row(2) << 0, 1.59488e-06, 2.17225e-06, -8.39095e-05, -0.000165608;
	  covMatrix.row(3) << -0, -0.000134473, -8.39095e-05,    0.0199097,   0.00602013;
	  covMatrix.row(4) << -0, -9.90522e-05, -0.000165608,   0.00602013,    0.0233725;
	  
	  //double lambda = Pinv * H/2;
	  //double lambda_prime = (lambda**3)

	  //	  propagateJacobian(Vout,A,B,C,H,covMatrix,J);
	  */
	  // set the local errors for the step
	  double err_x = H * H * std::fabs(A[0] - A[1] - A[2] + A[3]);
	  double err_y = H * H * std::fabs(B[0] - B[1] - B[2] + B[3]);
	  double err_z = H * H * std::fabs(C[0] - C[1] - C[2] + C[3]);
	  
	  //	   double err_px = (1/track_mom_mag_inv) * H * std::fabs(A[0] - A[1] - A[2] + A[3]);
	  //	   double err_py = (1/track_mom_mag_inv) * H * std::fabs(B[0] - B[1] - B[2] + B[3]);
	  //	   double err_pz = (1/track_mom_mag_inv) * H * std::fabs(C[0] - C[1] - C[2] + C[3]);
	  
	  pos_err[0] = err_x;
	  pos_err[1] = err_y;
	  pos_err[2] = err_z;
	  
	  // final convergence check
	  Rest = H - Tl;
	  
	  // can't have negative step size
	  if(H < 0.0) {
	    Rest = -Rest;
	  }
	  // check if converged or need to go round again
	  if(Rest < 1.0e-5*fabs(H)) {
	    Iloop = 0;
	  }
	}                // Icon3
      }                  // Icon2
    }                    // Icon1 
  } while (Iloop == 1);

//  std::cout << "END OF RKEXTRAPOLATION LOOP" << "\n";
//  std::cout << "---------------------------" << "\n";
    
  //convert back to Mev
  eStep.momentum.set( eStep.momentum.x()*gev2mev, eStep.momentum.y()*gev2mev, eStep.momentum.z()*gev2mev );
  
  //convert back to mm
  eStep.position.set( eStep.position.x()*cm2mm, eStep.position.y()*cm2mm, eStep.position.z()*cm2mm );
  return eStep;

}

void gm2strawtracker::GeaneExtrapolationUtils::propagateJacobian(double trackParameters[7], double A[4], double B[4], double C[4], double H, Eigen::MatrixXd covMatrix, Eigen::MatrixXd J)
{
  /*
  std::cout << "check J(2,2) = " << J(2,2) << "\n";
  std::cout << "check covMatrix(1,1) = " << covMatrix(1,1) << "\n";
  std::cout << "check covMatrix(2,1) = " << covMatrix(2,1) << "\n";

  std::cout << "H = " << H << "\n";
  
  for (int i(0); i<4; i++) {
    std::cout << "A[" << i << " ] = " << A[i] << "\n";
    std::cout << "B[" << i << " ] = " << B[i] << "\n";
    std::cout << "C[" << i << " ] = " << C[i] << "\n";
  }

  //
  // BEGIN DIFFERENTIATION  
  //


  double u[4]          = {trackParameters[0],trackParameters[1],trackParameters[2]}; // position
  double u_prime[4]    = {trackParameters[3],trackParameters[4],trackParameters[5]}; // direction

  double u_prime2_x[4] = {A[0],A[1],A[2],A[3]};
  double u_prime2_y[4] = {B[0],B[1],B[2],B[3]};
  double u_prime2_z[4] = {C[0],C[1],C[2],C[3]};

  // dF_n/du_n = 1 + H/3(dA0/du_n + dA1/du_n + dA2/du_n) etc
  
  // Ak matrices = differentiation of u_prime wrt u_prime
  // = dA0_dpx --> dA0_dP, dA1_dpx --> dA1_dP for A, B and C
  double dA0[4];
  double dB0[4];
  double dC0[4];

  double dA1[4];
  double dB1[4];
  double dC1[4];

  double dA2[4];
  double dB2[4];
  double dC2[4];

  double dA3[4];
  double dB3[4];
  double dC3[4];
  
  //dA0[0] = 0;
  //dA0[1] = lambda * BField[2];
  //dA0[2] = -lambda * Bfield[1];
  //dA0[3] = A[0] / lambda;
  //
  //dB0[0] = -lambda * Bfield[2];
  //dB0[1] = 0;
  //dB0[2] = 
  
  // Ck matrices = differentiation of u_prime wrt u
  // = dA0_dx --> dA0_dLambda, dA1_dx --> dA1_dLambda for A, B and C

  

  // transform covMatrix elements into global coords (7x7 matrix)
  
  double d_pu_px = covMatrix(1,1);
  double d_pv_px = covMatrix(2,2);
  double d_px_pz(0);
  double d_py_pz(0);
  
  // transform the errors on pu and pv into px and py errors
  geomUtils_.getXYfromUV(sgeom_,d_px_pz,d_py_pz,d_pu_px,d_pv_px);
  std::cout << "d_px_pz = " << d_px_pz << "\n";
  std::cout << "d_py_pz = " << d_py_pz << "\n";
  */

}


void GeaneExtrapolationUtils::extrapolateForwardsAndBackwards(gm2strawtracker::ExtrapolationStep& startingPoint, gm2strawtracker::ExtrapolationStep& finalPoint, std::vector<ExtrapolationStep>& forwardsSteps, std::vector<ExtrapolationStep>& backwardsSteps, double stepSize, bool useOldCode, bool forceUniformField, double distance){

  mf::LogInfo info(name_);

  bool leftFieldMap = false;

  // transform everything to world coordinates
  startingPoint.position = startingPoint.position.transform(cs_,"world");
  startingPoint.momentum = startingPoint.momentum.transform(cs_,"world",true);

  // check which field starting in
  const G4double initialPos[4] = {startingPoint.position.x(),startingPoint.position.y(),startingPoint.position.z(),0.0};
  double Bfield[3] = {0.0,0.0,0.0};
  fieldManager_->GetFieldValue(initialPos, Bfield);
  
  mf::LogDebug(name_) << "Field at start pos for forward/backwards plots = {" << Bfield[0] << ", " << Bfield[1] << ", " << Bfield[2] << "}\n";

  // initialize local errors
  double errs[3] = {0.0,0.0,0.0};
  
  gm2strawtracker::ExtrapolationStep initialForwardsStep;
  initialForwardsStep.position = startingPoint.position;
  initialForwardsStep.momentum = startingPoint.momentum;

  // initialize vector of forwards steps
  gm2strawtracker::ExtrapolationStep firstStepForwards = useOldCode ? oldGeaneExtrapolation(initialForwardsStep,stepSize,"forwards",leftFieldMap,forceUniformField) : GeaneExtrapolation(initialForwardsStep,stepSize,"forwards",errs,leftFieldMap,forceUniformField);
  
  // sanity check
  gm2strawtracker::ExtrapolationStep sanityStep;
  sanityStep.position = initialForwardsStep.position;
  sanityStep.momentum = -firstStepForwards.momentum;
  gm2strawtracker::ExtrapolationStep sanityStepBackwards = useOldCode ? oldGeaneExtrapolation(sanityStep,stepSize,"backwards",leftFieldMap,forceUniformField) : GeaneExtrapolation(sanityStep,stepSize,"backwards",errs,leftFieldMap,forceUniformField);
  
  firstStepForwards.position = startingPoint.position;
  firstStepForwards.momentum = startingPoint.momentum;

  forwardsSteps.push_back(firstStepForwards);
  
  double radius = (startingPoint.momentum.mag() * 1e6 * 1e3)  / (1.45127 * 3e8);

  // extrapolate backwards for n steps
  if (distance == -1.0) {
    distance = M_PI * 2 * radius * 0.1;
  }

  int nsteps = distance / stepSize;
    
  for (int istep(1); istep<nsteps+1; istep++) {

    gm2strawtracker::ExtrapolationStep nextStepForwards = useOldCode ? oldGeaneExtrapolation(forwardsSteps.at(istep-1),stepSize,"forwards",leftFieldMap,forceUniformField) : GeaneExtrapolation(forwardsSteps.at(istep-1),stepSize,"forwards",errs,leftFieldMap,forceUniformField);
    if (nextStepForwards.position.mag() == 0 && nextStepForwards.momentum.mag() == 0) {
      mf::LogDebug(name_) << "failed extrapolation" << "\n";
      return;
    }
    forwardsSteps.push_back(nextStepForwards);
  }

  // now go backwards the same distance
  gm2strawtracker::ExtrapolationStep firstStepBackwards;

  firstStepBackwards.position = forwardsSteps.at(forwardsSteps.size()-1).position;
  firstStepBackwards.momentum = -forwardsSteps.at(forwardsSteps.size()-2).momentum;
  backwardsSteps.push_back(firstStepBackwards);
  
  for (int istep(1); istep<nsteps+1; istep++) {
    gm2strawtracker::ExtrapolationStep nextStepBackwards = useOldCode ? oldGeaneExtrapolation(backwardsSteps.at(istep-1),stepSize,"backwards",forceUniformField) : GeaneExtrapolation(backwardsSteps.at(istep-1),stepSize,"backwards",errs,forceUniformField);
    if (nextStepBackwards.position.mag() == 0 && nextStepBackwards.momentum.mag() == 0) {
      return;
    }
    backwardsSteps.push_back(nextStepBackwards);
  }
  
  finalPoint = backwardsSteps.back();
  
  // make sure in world coords
  finalPoint.position.transform(cs_,"world");
  finalPoint.momentum.transform(cs_,"world",true);
  
  //std::cout << "finalPoint.position = " << finalPoint.position << "\n";
  //std::cout << "finalPoint.momentum = " << finalPoint.momentum << "\n\n";
}

void GeaneExtrapolationUtils::extrapolateBackwardsAndForwards(gm2strawtracker::ExtrapolationStep& startingPoint, gm2strawtracker::ExtrapolationStep& finalPoint, double stepSize, double distance){

  // transform everything to world coordinates
  startingPoint.position = startingPoint.position.transform(cs_,"world");
  startingPoint.momentum = startingPoint.momentum.transform(cs_,"world",true);

  // check which field starting in
  const G4double initialPos[4] = {startingPoint.position.x(),startingPoint.position.y(),startingPoint.position.z(),0.0};
  double Bfield[3] = {0.0,0.0,0.0};
  fieldManager_->GetFieldValue(initialPos, Bfield);
  
  mf::LogDebug(name_) << "Field at start pos for forward/backwards plots = {" << Bfield[0] << ", " << Bfield[1] << ", " << Bfield[2] << "}\n";
  
  // extrapolate in same direction as the startingPoint momentum then flip around and extrapolate the other way for same distance
  // - default is backwards then forwards, default coordSysName = "world"
  double errs[3] = {0.0,0.0,0.0};
  bool leftFieldMap = false;

  // do one forwards step to get the correct initial backwards momentum
  gm2strawtracker::ExtrapolationStep initialForwardsStep = GeaneExtrapolation(startingPoint,stepSize,"forwards",errs,leftFieldMap);
  if (initialForwardsStep.position.mag() == 0 && initialForwardsStep.momentum.mag() == 0) {
    mf::LogInfo(name_) << "went into helix condition and returned false" << "\n";
    return;
  }

  mf::LogDebug(name_) << "initialForwardsStep.position = " << initialForwardsStep.position << "\n";
  mf::LogDebug(name_) << "initialForwardsStep.momentum = " << initialForwardsStep.momentum << "\n";

  // initialize the first backwards step
  gm2strawtracker::ExtrapolationStep initialBackwardsStep;
  initialBackwardsStep.position = startingPoint.position;//initialForwardsStep.position;
  initialBackwardsStep.momentum = -initialForwardsStep.momentum;
  
  mf::LogDebug(name_) << "initialBackwardsStep.position = " << initialBackwardsStep.position << "\n";
  mf::LogDebug(name_) << "initialBackwardsStep.momentum = " << initialBackwardsStep.momentum << "\n";

  // initialize vector of backwards steps
  std::vector<gm2strawtracker::ExtrapolationStep> backwardsSteps;

  backwardsSteps.push_back(initialBackwardsStep);

  // extrapolate backwards for n steps
  int nsteps = distance / stepSize;
  
  for (int istep(1); istep<nsteps+1; istep++) {
    gm2strawtracker::ExtrapolationStep nextStepBackwards = GeaneExtrapolation(backwardsSteps.at(istep-1),stepSize,"backwards",errs,leftFieldMap);
    if (nextStepBackwards.position.mag() == 0 && nextStepBackwards.momentum.mag() == 0) {
      mf::LogInfo(name_) << "went into helix condition and returned false" << "\n";
      return;
    }
    
    mf::LogDebug(name_) << "nextStepBackwards.position = " << nextStepBackwards.position << "\n";
    mf::LogDebug(name_) << "nextStepBackwards.momentum = " << nextStepBackwards.momentum << "\n";
    
    const G4double nextBackwardsPos[4] = {nextStepBackwards.position.x(),nextStepBackwards.position.y(),nextStepBackwards.position.z(),0.0};
    fieldManager_->GetFieldValue(nextBackwardsPos, Bfield);
    
    backwardsSteps.push_back(nextStepBackwards);
  }
  
  // do one extra backwards step
  gm2strawtracker::ExtrapolationStep extraBackwardsStep = GeaneExtrapolation(backwardsSteps.back(),stepSize,"backwards",errs,leftFieldMap);
  if (extraBackwardsStep.position.mag() == 0 && extraBackwardsStep.momentum.mag() == 0) {
    mf::LogInfo(name_) << "went into helix condition and returned false" << "\n";
    return;
  }
  
  mf::LogDebug(name_) << "extraBackwardsStep.position = " << extraBackwardsStep.position << "\n";
  mf::LogDebug(name_) << "extraBackwardsStep.momentum = " << extraBackwardsStep.momentum << "\n";
  
  // now go forwards
  std::vector<gm2strawtracker::ExtrapolationStep> forwardsSteps; 
  gm2strawtracker::ExtrapolationStep firstStepForwards;
  firstStepForwards.position = backwardsSteps.back().position;
  firstStepForwards.momentum = -extraBackwardsStep.momentum;
  forwardsSteps.push_back(firstStepForwards);
  
  mf::LogDebug(name_) << "firstStepForwards.position = " << firstStepForwards.position << "\n";
  mf::LogDebug(name_) << "firstStepForwards.momentum = " << firstStepForwards.momentum << "\n";
  
  //  for (int istep(1); istep<nsteps+1; istep++) {
  for (int istep(1); istep<nsteps+1; istep++) {
    gm2strawtracker::ExtrapolationStep nextStepForwards = GeaneExtrapolation(forwardsSteps.at(istep-1),stepSize,"forwards",errs,leftFieldMap);
    if (nextStepForwards.position.mag() == 0 && nextStepForwards.momentum.mag() == 0) {
      mf::LogInfo(name_) << "went into helix condition and returned false" << "\n";
      return;
    }

    mf::LogDebug(name_) << "nextStepForwards.position = " << nextStepForwards.position << "\n";
    mf::LogDebug(name_) << "nextStepForwards.momentum = " << nextStepForwards.momentum << "\n";
    forwardsSteps.push_back(nextStepForwards);

    const G4double nextForwardsPos[4] = {nextStepForwards.position.x(),nextStepForwards.position.y(),nextStepForwards.position.z(),0.0};
    //    double Bfield[3] = {0.0,0.0,0.0};
    fieldManager_->GetFieldValue(nextForwardsPos, Bfield);
  }
  
  finalPoint = forwardsSteps.back();
  
  // make sure in world coords
  finalPoint.position.transform(cs_,"world");
  finalPoint.momentum.transform(cs_,"world",true);
    
}


bool GeaneExtrapolationUtils::GeaneExtrapolationUtils::didStepHitVolume(std::string givenVolume) {
  std::string volume = givenVolume;
  bool checkOnlyOneVolume = false;
  
  std::vector<std::string> expectedVolumes = {"SingleStraw","strawModule","ArcSection","world"}; // expect that all tracks hit these volumes
  std::vector<std::string> expectedVolumesForwards = {"SingleStraw","strawModule","ArcSection","world","xtal","Calorimeter","insideCalo","VacuumChamberCadMesh"}; // expect that all tracks hit these volumes
  std::vector<std::string> duplicateVolumeNames = {"mylarWindowPV","SingleStraw","strawModule","VacuumChamberCadMesh","xtal","quadField","Calorimeter","PbF2Bounding","insideCalo","trolleyRail","quadPlate","Collimator","StationNumber","frontWrapping","ArcSection","bellowsRail","quadStandoff","opticalCoupling","photodetector"}; // many different volumes with these names but don't care specifically which one track hits
#if 0
  // check if the step position is still in the world region. LocateGlobalPointAndSetup will seg fault if no volume found.
  if (fabs(stepPos.x()) > 10000 || fabs(stepPos.y()) > 3500  || fabs(stepPos.z()) > 10000 ) {
    mf::LogInfo(name_) << "Step pos = " << stepPos << " which is outside world volume; set volume = outsideWorld and return true for didStepHitVolume." << "\n";
    volume = "outsideWorld";
    return true;
  }
  if (nav_->GetWorldVolume()) volume = nav_->LocateGlobalPointAndSetup(stepPos)->GetName();
#endif
  for (auto vol : duplicateVolumeNames) {
    if (volume.find(vol) != std::string::npos) volume = vol;
  }
  std::string direction = "backwards";
  //std::cout<<"After for loop"<<std::endl;
  if (!checkOnlyOneVolume) {
    if (direction != "forwards" && (std::find(expectedVolumes.begin(),expectedVolumes.end(),volume) == expectedVolumes.end())) {
      mf::LogDebug(name_) << "from map of bad volumes, step hit volume " << volume << ", return true\n";
      //std::cout<<"returning true"<<std::endl;
      return true;
    }
    else if (direction == "forwards" && (std::find(expectedVolumesForwards.begin(),expectedVolumesForwards.end(),volume) == expectedVolumesForwards.end())) {
      //std::cout<<"returning true"<<std::endl;
      return true;
    }
    else {
      mf::LogDebug(name_) << "volume " << volume << " is not in the list of bad volumes so return false " << "\n";
      //std::cout<<"returning false"<<std::endl;
      return false;
    }
  } 
  
  else {
    std::string volumeToCheck = "mylarWindowPV";
    //std::cout<<"mylarWindowPV block"<<std::endl;
    if (volume.find(volumeToCheck) != std::string::npos) {
      mf::LogDebug(name_) << "only checking if hit " << volumeToCheck << ", and have hit " << volume << "\n";
      return true;
    }
    else return false;
  }
}
