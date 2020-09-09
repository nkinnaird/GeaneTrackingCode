// Module that takes tracks and rotates them around the y-axis (which should have a rotational symmetry).
// Before our bug fix, this had big peaks in the chi2 of fit track at +/- z axes.

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "artg4/material/Materials.hh"

// Geane & Geant4 includes
#include "artg4/services/DetectorHolder_service.hh"
#include "artg4/gm2Geane/gm2GeanePropagator.hh"
#include "Geant4/G4ErrorPropagatorData.hh"
#include "artg4/gm2Geane/gm2GeanePropagatorManager.hh"
#include "artg4/gm2Geane/gm2GeanePlaneSurfaceTarget.hh"
#include "artg4/gm2Geane/gm2GeaneFreeTrajState.hh"
#include "artg4/gm2Geane/gm2GeaneSurfaceTrajState.hh"
#include "artg4/gm2Geane/gm2GeaneMatrix.hh"
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4FieldManager.hh"
#include "Geant4/G4UniformMagField.hh"
#include "Geant4/G4TransportationManager.hh"
#include "gm2geom/fields/gm2FieldManager_service.hh"

// Output TFile includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "gm2util/common/RootManager.hh"
#include "TRandom3.h"
#include "TGraph.h"

class GeaneRotationTester;

class GeaneRotationTester : public art::EDAnalyzer {

public:

  explicit GeaneRotationTester(fhicl::ParameterSet const & p);

  void analyze(const art::Event & e) override;
  void beginJob() override;

  int    PropagateTrack(double truePhi, G4ThreeVector startPos, G4ThreeVector startMom, std::vector<Eigen::VectorXd>& predictedParams, std::vector<Eigen::MatrixXd>& transferMatrices, std::vector<Eigen::MatrixXd>& errorMatrices);
  double TrackIteration(double truePhi, G4ThreeVector& startPos, G4ThreeVector& startMom, std::vector<Eigen::VectorXd> predictedParams, std::vector<Eigen::MatrixXd> transferMatrices, std::vector<Eigen::MatrixXd> errorMatrices, std::vector<double> measHitsH, std::vector<double> measHitsV);

private :

  gm2GeanePropagatorManager* g4emgr_;
  G4ErrorPropagatorData* g4edata_; 

  //ROOT plotting members
  std::string name_;
  std::unique_ptr<RootManager> rootManager_;
  std::string dirName_;

  // Some setup params
  double By_;
  int    nPlanes_;
  double planeSep_;
  double measRes_;
};

GeaneRotationTester::GeaneRotationTester(fhicl::ParameterSet const & p)
  : art::EDAnalyzer(p)
  , name_("GeaneRotationTester")
  , rootManager_()
  , dirName_( p.get<std::string>("dirName","GeaneRotationTester") )
{

  // Instantiate everything for GEANE
  art::ServiceHandle<artg4::DetectorHolderService> dh;
  dh->initialize();
  dh->constructAllLVs();
  G4VPhysicalVolume* world = dh->worldPhysicalVolume();
  world->GetLogicalVolume()->UpdateMaterial(artg4Materials::Vacuum());

  // Set magnetic field in y direction
  By_ = 1.45*tesla;
  G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0.,By_,0.));
  G4TransportationManager::GetTransportationManager()->GetFieldManager()->SetDetectorField(magField);
  G4TransportationManager::GetTransportationManager()->GetFieldManager()->CreateChordFinder(magField);

  // Initialize the GEANT4e manager 
  g4emgr_ = gm2GeanePropagatorManager::GetErrorPropagatorManager();
  g4emgr_->SetUserInitialization(world);
  g4emgr_->InitGeant4e();

  // Initialize the GEANT4e Propagator Data 
  g4edata_ = G4ErrorPropagatorData::GetErrorPropagatorData();

  // "Detector" config
  nPlanes_ = 8;
  planeSep_ = 10*cm;
  measRes_ = 0.15*mm; 
}


void GeaneRotationTester::beginJob(){

  // Grab output file
  art::ServiceHandle<art::TFileService> tfs;
  auto& outputRootFile = tfs->file();
  rootManager_.reset( new RootManager(name_, &outputRootFile) ); 
  
  //Create directory 
  rootManager_->GetDir(&outputRootFile,dirName_,true); //true -> create if doesn't exist

}


void GeaneRotationTester::analyze(const art::Event & e){

  // Graphs of position and momentum - make this for first track (which has same measured offsets)
  TGraph* tgX = new TGraph();
  tgX->SetMarkerStyle(20);
  tgX->SetMarkerSize(0.4);
  tgX->SetName("X");
  tgX->SetTitle(";#phi_{z};x [mm]");

  TGraph* tgY = new TGraph();
  tgY->SetMarkerStyle(20);
  tgY->SetMarkerSize(0.4);
  tgY->SetName("Y");
  tgY->SetTitle(";#phi_{z};y [mm]");

  TGraph* tgZ = new TGraph();
  tgZ->SetMarkerStyle(20);
  tgZ->SetMarkerSize(0.4);
  tgZ->SetName("Z");
  tgZ->SetTitle(";#phi_{z};z [mm]");

  TGraph* tgP = new TGraph();
  tgP->SetMarkerStyle(20);
  tgP->SetMarkerSize(0.4);
  tgP->SetName("P");
  tgP->SetTitle(";#phi_{z};p [MeV]");

  TGraph* tgPX = new TGraph();
  tgPX->SetMarkerStyle(20);
  tgPX->SetMarkerSize(0.4);
  tgPX->SetName("PX");
  tgPX->SetTitle(";#phi_{z};p_{x} [MeV]");

  TGraph* tgPY = new TGraph();
  tgPY->SetMarkerStyle(20);
  tgPY->SetMarkerSize(0.4);
  tgPY->SetName("PY");
  tgPY->SetTitle(";#phi_{z};p_{y} [MeV]");

  TGraph* tgPZ = new TGraph();
  tgPZ->SetMarkerStyle(20);
  tgPZ->SetMarkerSize(0.4);
  tgPZ->SetName("PZ");
  tgPZ->SetTitle(";#phi_{z};p_{z} [MeV]");

  // Histogram of track fit details (Chi2 vs angle)
  TH1D* chi2 = new TH1D("chi2",";Track #chi^{2}",500,0,50);
  TH1D* pVal = new TH1D("pVal",";Track p-value",500,0,1);
  TH2D* chi2VsLambda = new TH2D("chi2VsLambda",";#phi_{z};Track #chi^{2}",381,-100.5,280.5,100,0,50);

  // Some basic running config
  int nAnglePts = 360*1;
  int tracksPerPt = 10000;

  // Track config parameters
  double trueMom = 1500*MeV;

  // Calculate true hit positions w.r.t particle direction (these are same in our rotating frame)
  vector<double> trueHitH, trueHitV;
  for(int iPlane = 0; iPlane < nPlanes_; iPlane++){
    double planeX = (iPlane+1)*planeSep_;
    double r_c = (trueMom/By_) / 299.792458;
    double hitH = By_ > 0 ? r_c - sqrt(r_c*r_c - planeX*planeX) : 0;
    trueHitV.push_back(0);
    trueHitH.push_back(hitH);
  }

  // Loop over angle w.r.t. x-axis (denoted phi here)
  for(int i = 0; i < nAnglePts; i++){

    // Set track angle w.r.t x-axis
    double truePhi = i*TMath::TwoPi()/nAnglePts;
    if(i % 10 == 0) std::cout << "Phi = " << truePhi << "..." << std::endl;

    // Random number generator - recreate with same seed each time so we get same tracks at each angle
    TRandom3 rng(12345);

    // Fire many tracks at this angle
    for(int iTrack = 0; iTrack < tracksPerPt; iTrack++){
      
      // Comment measHits out here if you want to fire same track and hits
      // Holders for hit information (horizontal and vertical in plane perpendicular to truePhi).
      vector<double> measHitH, measHitV;

      // Smear truth to get measured hits
      for(int iPlane = 0; iPlane < nPlanes_; iPlane++){
       	measHitH.push_back(trueHitH.at(iPlane) + rng.Gaus(0,measRes_));
       	measHitV.push_back(trueHitV.at(iPlane) + rng.Gaus(0,measRes_));
      }

      // Guess a starting track - fire from origin at angle phi
      G4ThreeVector startPos( 0, 0, 0 );
      G4ThreeVector startMom( trueMom*cos(truePhi)*MeV, 0, trueMom*sin(truePhi)*MeV);  // In xz plane at angle truePhi

      // Iterate twice (one correction)
      double chiSquared;
      bool skipTrack = false;
      for(int iteration = 0; iteration < 2; iteration++){

	// Propagate track to get predicted parameters, transfer and error matrices
	std::vector<Eigen::VectorXd> predictedParams;
	std::vector<Eigen::MatrixXd> transferMatrices, errorMatrices;
	int error = PropagateTrack(truePhi, startPos, startMom, predictedParams, transferMatrices, errorMatrices);
	if(error > 0) {
	  std::cout << "error = " << error << std::endl;
	  skipTrack = true;
	  break;
	}

	// Iterate track (linearly correct starting params to minimize chi2)
	chiSquared = TrackIteration(truePhi, startPos, startMom, predictedParams, transferMatrices, errorMatrices, measHitH, measHitV);

      }     
      
      if(skipTrack) continue;

      // Chi-squared after fit plot
      chi2->Fill(chiSquared);
      pVal->Fill(TMath::Prob(chiSquared,nPlanes_*2-5)); // 5 pars
      chi2VsLambda->Fill(180*(truePhi-TMath::Pi()/2.)/TMath::Pi(),chiSquared);

      // Put parameters in graphs for first track
      if(iTrack == 0){
	//	std::cout << "phi = " << truePhi << "\tFilling with : " << startPos << " : " << startMom << std::endl;
	tgX->SetPoint(tgX->GetN(), 180*(truePhi-TMath::Pi()/2.)/TMath::Pi(), startPos.x());
	tgY->SetPoint(tgY->GetN(), 180*(truePhi-TMath::Pi()/2.)/TMath::Pi(), startPos.y());
	tgZ->SetPoint(tgZ->GetN(), 180*(truePhi-TMath::Pi()/2.)/TMath::Pi(), startPos.z());
	tgP->SetPoint(tgP->GetN(), 180*(truePhi-TMath::Pi()/2.)/TMath::Pi(), startMom.mag());
	tgPX->SetPoint(tgPX->GetN(), 180*(truePhi-TMath::Pi()/2.)/TMath::Pi(), startMom.x());
	tgPY->SetPoint(tgPY->GetN(), 180*(truePhi-TMath::Pi()/2.)/TMath::Pi(), startMom.y());
	tgPZ->SetPoint(tgPZ->GetN(), 180*(truePhi-TMath::Pi()/2.)/TMath::Pi(), startMom.z());
      }
 
    } // Loop over tracks

  } // Loop over angles

  // Add histograms to TFile
  auto topDir = rootManager_->GetDir(dirName_);
  rootManager_->Add(topDir,chi2);
  rootManager_->Add(topDir,pVal);
  rootManager_->Add(topDir,chi2VsLambda);
  rootManager_->Add(topDir,tgX);
  rootManager_->Add(topDir,tgY);
  rootManager_->Add(topDir,tgZ);
  rootManager_->Add(topDir,tgP);
  rootManager_->Add(topDir,tgPX);
  rootManager_->Add(topDir,tgPY);
  rootManager_->Add(topDir,tgPZ);

  return;
}


int GeaneRotationTester::PropagateTrack(double truePhi, G4ThreeVector startPos, G4ThreeVector startMom, std::vector<Eigen::VectorXd>& predictedParams, std::vector<Eigen::MatrixXd>& transferMatrices, std::vector<Eigen::MatrixXd>& errorMatrices){

  std::vector<gm2GeaneMatrix> freeTransferMatrices; // free coordinate system
  std::vector<gm2GeaneMatrix> surfTransferMatrices; // set of transfer/transport matrices, one for each plane, describes transport from one plane to the next plane. Surface system
  std::vector<gm2GeaneMatrix> surfErrorMatrices; // set of error matrices, one for each plane - surface system
  std::vector<gm2GeaneMatrix> SC2SDTransformationMatrix; // Matrices for transformation to and from free system and plane/surface system. These will be filled with per-step Jacobian matrices, not per-plane matrices.
  std::vector<gm2GeaneMatrix> SC2SDTransformationMatrixInverse;

  // Free state (that we will propagate)
  gm2GeaneFreeTrajState* myFreeTrajState = new gm2GeaneFreeTrajState("e+", startPos, startMom );

  // A surface state. The passed empty matrix gets filled with the Jacobian transformation matrix.
  G4Normal3D surfNorm(-cos(truePhi),0,-sin(truePhi));
  G4ThreeVector surfVert(0,1,0);
  G4ThreeVector surfHoriz = surfVert.cross(surfNorm);
  gm2GeaneMatrix jacobMatrix(5,0);
  gm2GeaneSurfaceTrajState* startingSurfState = new gm2GeaneSurfaceTrajState(*myFreeTrajState, surfVert, surfHoriz, jacobMatrix);

  // Push back the Jacobian and the inverse
  SC2SDTransformationMatrix.push_back(jacobMatrix);
  int tmp; // Int to pass and see if the inversion failed.
  SC2SDTransformationMatrixInverse.push_back(SC2SDTransformationMatrix.at(0).inverse(tmp));
  
  // Push back other matrices
  freeTransferMatrices.push_back(myFreeTrajState->GetTransfMat());   // Transport matrix for the free state when no steps have been made is a zero matrix.
  surfTransferMatrices.push_back(myFreeTrajState->GetTransfMat()); // Surface transfer matrices also just filled with the 0 matrix.
  surfErrorMatrices.push_back(startingSurfState->GetError());        // Starts off with 0 matrix as well, but writing more naturally.
 
  // Loop over planes that we'll propagate between - set at plane distance from origin normal to angle truePhi
  for(int iPlane = 1; iPlane <= nPlanes_; iPlane++){
	
    // Track bends in field so won't be perpendicular to this plane
    G4Point3D surfPos = surfNorm*iPlane*planeSep_*-1;
    G4ErrorTarget* target = new gm2GeanePlaneSurfaceTarget(surfNorm, surfPos );
    g4edata_->SetTarget( target );

    // Propagate track to surface target
    g4emgr_->InitTrackPropagation();

    bool moreSteps = true;
    int stepNum = 0;
    while(moreSteps){

      int ierr = g4emgr_->PropagateOneStep( myFreeTrajState, G4ErrorMode::G4ErrorMode_PropForwards );
	  
      // ierr is the failure integer which is non-zero when step length is super small, propagation is exactly along z, or a few other things.
      // In this case transport matrices and things will get set to 0 and later on the program will seg fault, so break out and skip this event if it's non-zero here.
      if(ierr != 0){
	std::cout << "ierr = " << ierr << " is non-zero, skipping event." << std::endl;
	return ierr;
      }
	  
      // See some issues where geant returns a step length of 0
      if (myFreeTrajState->GetG4Track()->GetStepLength() == 0.0){
	std::cout << "Step length was seen to be 0 - Geant is having step length issues - return out and fail. \n";
	return 10;
      }
	  
      // Some steps return the max step length and then the tracking fails - unknown root cause
      else if(myFreeTrajState->GetG4Track()->GetStepLength() == 1000){
	std::cout << "Step length was seen to be 1000, and maxed out - Geant is having step length issues - return out and fail. \n";
	return 11;
      }
	  
      // --- Check if target is reached
      if( g4emgr_->GetPropagator()->CheckIfLastStep( myFreeTrajState->GetG4Track() )) {
	g4emgr_->GetPropagator()->InvokePostUserTrackingAction( myFreeTrajState->GetG4Track() );  
	moreSteps = false;
      }

      // Use step by step propagation to properly save the transfer matrices.
      gm2GeaneMatrix stepSurfTransfMat(5,0);
      gm2GeaneSurfaceTrajState* stepSurfState = new gm2GeaneSurfaceTrajState(*myFreeTrajState, surfVert, surfHoriz, stepSurfTransfMat);
      // Only want the transformation matrices on the planes, not the individual steps - so remove the last then add the most recent.
      if (stepNum > 0){
	SC2SDTransformationMatrix.pop_back();
	SC2SDTransformationMatrixInverse.pop_back();
      }
      SC2SDTransformationMatrix.push_back(stepSurfTransfMat);
      SC2SDTransformationMatrixInverse.push_back((SC2SDTransformationMatrix.back()).inverse(tmp));

      // Fill transport matrices.
      if(stepNum == 0) { 
	// For matrix multiplication, it should be: R10 = A1 * T10 * A0^-1, where R and T are the transport matrices in separate coordinate systems from plane 0 to 1, and A is the transformation between the two. 
	freeTransferMatrices.push_back(myFreeTrajState->GetTransfMat()); // Gets the transport or transfer matrix for the last step. If there is only a single step between planes this is fine.
	surfTransferMatrices.push_back(SC2SDTransformationMatrix.at(iPlane)*(myFreeTrajState->GetTransfMat())*SC2SDTransformationMatrixInverse.at(iPlane-1));              
      }
      else {
	freeTransferMatrices.at(iPlane) = myFreeTrajState->GetTransfMat()*freeTransferMatrices.at(iPlane);
	// Here I gather up the transport matrices for the free system for all steps between planes, then I multiply on the right by the previous plane coordinate transformation matrix, and on the left by the next step coordinate transformation matrix. It keeps looping until the next step transformation matrix is the same as the next plane transformation matrix.
	surfTransferMatrices.at(iPlane) = (SC2SDTransformationMatrix.at(iPlane)*(freeTransferMatrices.at(iPlane))*SC2SDTransformationMatrixInverse.at(iPlane-1));
      }
      
      stepNum++; // Increment the step number
      delete stepSurfState;

    } // end of while

    // After getting to the target, create the surface trajectory state on the target and get error matrix
    gm2GeaneSurfaceTrajState* mySurfTrajState = new gm2GeaneSurfaceTrajState(*myFreeTrajState, surfVert, surfHoriz, jacobMatrix);

    // Get predicted parameters
    double invP = mySurfTrajState->GetParameters().GetInvP();
    double pV = mySurfTrajState->GetParameters().GetPV();
    double pW = mySurfTrajState->GetParameters().GetPW();
    double pX = sqrt( 1./invP * 1./invP - pV*pV - pW*pW);
    double V = mySurfTrajState->GetParameters().GetV();
    double W = mySurfTrajState->GetParameters().GetW();
    predictedParams.push_back(Eigen::VectorXd::Zero(5));
    predictedParams.at(iPlane-1)[0] = invP;
    predictedParams.at(iPlane-1)[1] = pV/pX;
    predictedParams.at(iPlane-1)[2] = pW/pX;
    predictedParams.at(iPlane-1)[3] = V;
    predictedParams.at(iPlane-1)[4] = W;

    // Get error matrix
    surfErrorMatrices.push_back(mySurfTrajState->GetError());

    // Delete the target and surface trajectory states
    delete target;
    delete mySurfTrajState;

  } // Loop over planes

  // Loop from 0 to nPlanes (inclusive) to fill eigen objects with these matrices
  for(int iPlane = 0; iPlane <= nPlanes_; iPlane++){
    // Convert these matrices to Eigen form
    transferMatrices.push_back(Eigen::MatrixXd::Zero(5,5));
    errorMatrices.push_back(Eigen::MatrixXd::Zero(5,5));
    for(int i = 0; i < 5; i++) {
      for(int j = 0; j < 5; j++) {
	transferMatrices.at(iPlane)(i,j) = surfTransferMatrices.at(iPlane)[i][j]*1.; // GEVCM (Don't scale them to MeV mm.)
	errorMatrices.at(iPlane)(i,j) = surfErrorMatrices.at(iPlane)[i][j]*1.; // GEVCM
      }
    } 
  }

  // Delete free state
  delete myFreeTrajState;

  return 0;
}


double GeaneRotationTester::TrackIteration(double truePhi, G4ThreeVector& startPos, G4ThreeVector& startMom, std::vector<Eigen::VectorXd> predictedParams, std::vector<Eigen::MatrixXd> transferMatrices, std::vector<Eigen::MatrixXd> errorMatrices, std::vector<double> measHitsH, std::vector<double> measHitsV ){

  ///////////////////////////////////////////////////////////////////////////////////// 
  // PAY ATTENTION TO THIS !!!!!!
  // NOTE NOTE NOTE: mainDiagonalErrorMatrices is shifted down 1 from errorMatrices
  // the latter has a 0 matrix for the 0 plane, while the former doesn't consider the fictional 0 plane at all
  // this now is sort of run around by looping through the trackPlanesHitList array 

  std::vector<Eigen::MatrixXd> mainDiagonalErrorMatrices;
  std::vector<Eigen::MatrixXd> mainDiagonalErrorMatricesInverse;

  for (int iPlane = 0; iPlane < nPlanes_; iPlane++) { 
    mainDiagonalErrorMatrices.push_back(std::numeric_limits<double>::infinity() * Eigen::MatrixXd::Identity(5,5)); // Fill infinities into non-measured elements.   
    mainDiagonalErrorMatrices.at(iPlane)(3,3) = 0.1*measRes_*0.1*measRes_;
    mainDiagonalErrorMatrices.at(iPlane)(4,4) = 0.1*measRes_*0.1*measRes_;
    //mainDiagonalErrorMatrices.at(iPlane)(4,4) = errorMatrices.at(iPlane+1)(4,4) +0.1*measRes_*0.1*measRes_;
    //mainDiagonalErrorMatrices.at(iPlane)(3,3) = errorMatrices.at(iPlane+1)(3,3) +0.1*measRes_*0.1*measRes_;

    mainDiagonalErrorMatricesInverse.push_back(Eigen::MatrixXd::Zero(5,5));
    mainDiagonalErrorMatricesInverse[iPlane](3,3) = 1./(0.1*measRes_*0.1*measRes_);
    mainDiagonalErrorMatricesInverse[iPlane](4,4) = 1./(0.1*measRes_*0.1*measRes_);
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // Multiply transport matrices together. Be careful with Transport matrix = 0 matrix for the starting plane.
  std::vector<Eigen::MatrixXd> transportMatrixBegToEnd;
  for (int iPlane = 0; iPlane <= nPlanes_; ++iPlane){
    transportMatrixBegToEnd.push_back(Eigen::MatrixXd::Zero(5,5)); // Initialize matrices.
  }
  for (int iPlane = 1; iPlane <= nPlanes_; ++iPlane){
    if(iPlane == 1){
      transportMatrixBegToEnd.at(1) = transferMatrices.at(1);
    } else {
      // Multiply previous full transport matrices on the left by the next single gap transport matrix.
      transportMatrixBegToEnd.at(iPlane) = transferMatrices.at(iPlane)*transportMatrixBegToEnd.at(iPlane-1);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////      
  // Form the overall 5x5 covariance matrix, equation 33 in the geane manual. (uninverted)
  // Change how covariance total inverse is calculated to equation 27 in geane manual
  //covarianceTotalInv = (combinedTransportMatrices.transpose())*reducedMatrixInverse*combinedTransportMatrices;
  Eigen::MatrixXd covarianceTotalInv(5,5);
  covarianceTotalInv.setZero();

  // At each hit plane planenum, add on transport matrix transpose(0, iPlane) times reducedErrorAtPlaneInverse[iPlane] times transport matrix(0, iPlane)
  for (int iPlane = 0; iPlane < nPlanes_; ++iPlane){
    int planeNum = iPlane+1;
    covarianceTotalInv += transportMatrixBegToEnd.at(planeNum).transpose()*mainDiagonalErrorMatricesInverse.at(iPlane)*transportMatrixBegToEnd.at(planeNum);
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // Convert measured parameter arrays into Eigen objects and YZ measured to UV measured:
  std::vector<Eigen::VectorXd> measParams;
  for(int iPlane = 0; iPlane < nPlanes_; iPlane++) {
    measParams.push_back(Eigen::VectorXd::Zero(5));
    measParams.at(iPlane)(3) = measHitsV.at(iPlane);
    measParams.at(iPlane)(4) = measHitsH.at(iPlane);
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // Calculate the residuals at every step: - Far right side of equation 23, 24, or 26 in the geane manual paper.
  std::vector<Eigen::VectorXd> residuals;
  for(int iPlane = 0; iPlane < nPlanes_; iPlane++) {
    residuals.push_back(measParams.at(iPlane)- predictedParams.at(iPlane));
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // Chi-squared calculation
  double chi2= 0;
  for(int iPlane = 0; iPlane < nPlanes_; ++iPlane) {
    chi2 += pow(residuals.at(iPlane)(3)/measRes_, 2);
    chi2 += pow(residuals.at(iPlane)(4)/measRes_, 2);
  }
  
  // GEVCM Convert residuals from MeV mm to GeV cm to mulitply against GeV cm transport and error matrices.
  for (int iPlane = 0; iPlane < nPlanes_; ++iPlane){
    residuals.at(iPlane)(0) = residuals.at(iPlane)(0) * 1.0e3;
    residuals.at(iPlane)(1) = residuals.at(iPlane)(1) * 1.0e0;
    residuals.at(iPlane)(2) = residuals.at(iPlane)(2) * 1.0e0;
    residuals.at(iPlane)(3) = residuals.at(iPlane)(3) * 1.0e-1;
    residuals.at(iPlane)(4) = residuals.at(iPlane)(4) * 1.0e-1;
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // Calculate correction to starting track

  // Sum over all planes: --- This is the right side of equation 26 or 32 from the geane manual paper. 
  Eigen::VectorXd trialTrajectoryAllPlanes(5);
  //trialTrajectoryAllPlanes = combinedTransportMatrices.transpose()*reducedMatrixInverse*combinedResiduals;
  trialTrajectoryAllPlanes.setZero();
  for(int iPlane = 0; iPlane < nPlanes_; ++iPlane){
    int planeNum = iPlane+1;
    trialTrajectoryAllPlanes += transportMatrixBegToEnd.at(planeNum).transpose()*mainDiagonalErrorMatricesInverse.at(iPlane)*residuals.at(iPlane);
  }  

  // Calculate the solution:   this is deltaPsi from equation 26 or 32 of the geane manual paper
  // Solves the equation Ax = b for x, where A = covarianceTotalInv, b = trialTrajectoryAllPlanes, and x = deltaPsi
  Eigen::VectorXd deltaPsi;
  deltaPsi = covarianceTotalInv.colPivHouseholderQr().solve(trialTrajectoryAllPlanes); 

  /////////////////////////////////////////////////////////////////////////////////////
  // Apply correction to starting parameters
  
  // Convert back to MeV mm units for the main code to use.
  deltaPsi(0) = deltaPsi(0) * 1.0e-3; 
  deltaPsi(1) = deltaPsi(1) * 1.0e0; 
  deltaPsi(2) = deltaPsi(2) * 1.0e0;
  deltaPsi(3) = deltaPsi(3) * 1.0e1;
  deltaPsi(4) = deltaPsi(4) * 1.0e1;

  // startPos & startMom are in world coordinates.  Change them to detector planes here (by rotating by truePhi) so we can add deltaPsi
  G4ThreeVector localPos( startPos.x()*cos(truePhi)+startPos.z()*sin(truePhi), startPos.y(), -startPos.x()*sin(truePhi)+startPos.z()*cos(truePhi) );
  G4ThreeVector localMom( startMom.x()*cos(truePhi)+startMom.z()*sin(truePhi), startMom.y(), -startMom.x()*sin(truePhi)+startMom.z()*cos(truePhi) );

  // Calculate new momentum and position (using deltaPsi) and update localMom
  double recoMomentum = 1./(deltaPsi[0]+(1./localMom.mag()));
  double pvpOVERpup = deltaPsi[1] + localMom.y()/localMom.x();
  double pwpOVERpup = deltaPsi[2] + localMom.z()/localMom.x();
  double updatedPX = sqrt( (recoMomentum*recoMomentum) / (1 + pvpOVERpup*pvpOVERpup + pwpOVERpup*pwpOVERpup) ); // This is only correct for an orthogonal system. 
  localMom.setX(updatedPX);
  localMom.setY(pvpOVERpup*updatedPX);
  localMom.setZ(pwpOVERpup*updatedPX);
  localPos.setY(localPos.y()+deltaPsi[3]);
  localPos.setZ(localPos.z()+deltaPsi[4]);

  // Rotate back to global system and pass back corrected values
  startMom.setX(localMom.x()*cos(truePhi) - localMom.z()*sin(truePhi));
  startMom.setY(localMom.y());
  startMom.setZ(localMom.x()*sin(truePhi) + localMom.z()*cos(truePhi));
  startPos.setX(localPos.x()*cos(truePhi) - localPos.z()*sin(truePhi));
  startPos.setY(localPos.y());
  startPos.setZ(localPos.x()*sin(truePhi) + localPos.z()*cos(truePhi));

  return chi2;
}


   
DEFINE_ART_MODULE(GeaneRotationTester)
