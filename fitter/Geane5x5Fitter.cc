//! header
#include "Geane5x5Fitter.hh"

using namespace gm2strawtracker;
using namespace std; 

//! constructor
Geane5x5Fitter::Geane5x5Fitter(fhicl::ParameterSet const & p) 
    : matrixDebug_(p.get<bool>("matrixDebug",false))
    , name_("Geane5x5Fitter")
    , sgeom_()
    , geaneTrackUtils_()
    , numTrackParams(5)
    , skipLayers_(p.get<vector<unsigned int> >("skipLayers",{}))
{
  maxNumPlanes = geaneTrackUtils_.maxNumPlanes;
  wireMeasurementError = 2*sgeom_.outerRadiusOfTheGas/sqrt(12); // this is used for the sequence checking stuff - specifically the hybrid error matrix
}

double Geane5x5Fitter::TrackCorrelation(gm2strawtracker::TrackDetailArtRecord& track)
{
  auto geaneHitsOnTrackDetails = track.geaneHits;

  double startingXPos = geaneHitsOnTrackDetails.startingGeaneParameters[0];
  double startingYPos = geaneHitsOnTrackDetails.startingGeaneParameters[1];
  double startingZPos = geaneHitsOnTrackDetails.startingGeaneParameters[2];
  double startingXMom = geaneHitsOnTrackDetails.startingGeaneParameters[3];
  double startingYMom = geaneHitsOnTrackDetails.startingGeaneParameters[4];
  double startingZMom = geaneHitsOnTrackDetails.startingGeaneParameters[5];

  double totalStartMomentum = sqrt ((startingXMom)*(startingXMom) 
                              + (startingYMom)*(startingYMom)
                              + (startingZMom)*(startingZMom));


  std::vector<Eigen::VectorXd> paramMeasuredEigen;
  std::vector<Eigen::VectorXd> paramPredictedEigen; // Already in XY
  std::vector<Eigen::VectorXd> paramPredictedInUVEigen;
  std::vector<Eigen::VectorXd> residuals;

  Eigen::MatrixXd covarianceTotalInv(5,5);

  Eigen::VectorXd deltaPsiNought;
  Eigen::VectorXd deltaPsiNoughtXY;

  Eigen::VectorXd trialTrajectoryAllPlanes(5);

  std::vector<Eigen::MatrixXd> transportMatrixBegToEnd;

  std::vector<Eigen::MatrixXd> mainDiagonalErrorMatrices;
  std::vector<Eigen::MatrixXd> mainDiagonalErrorMatricesInverse;

  double unmeasured = 1.e300; // large value for measurement uncertainty of V in U straws, etc. Is used purely as a placeholder for which rows and columns should be removed before multiplications.

  std::vector<Eigen::MatrixXd> myTransferMatricesGeVcm;
  std::vector<Eigen::MatrixXd> errorMatricesGeVcm;

  for (int ipl = 0; ipl < maxNumPlanes; ++ipl) // Use this loop to convert XY transport and error matrices to VU.
  {
    myTransferMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));
    errorMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));  

    myTransferMatricesGeVcm[ipl] = geaneTrackUtils_.JacobianToUV * geaneHitsOnTrackDetails.geaneTransportMatrices[ipl] * geaneTrackUtils_.JacobianToUV.inverse();
    errorMatricesGeVcm[ipl]      = geaneTrackUtils_.JacobianToUV * geaneHitsOnTrackDetails.geaneErrorMatrices[ipl] * geaneTrackUtils_.JacobianToUV.transpose();
  }

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
     mf::LogTrace(name_) << "\n";
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       mf::LogTrace(name_) << "Transport matrices in GeV cm XY space " << ipl << "\n";
       mf::LogTrace(name_) << (geaneHitsOnTrackDetails.geaneTransportMatrices)[ipl];
       mf::LogTrace(name_) << "\n \n";
     }

     mf::LogTrace(name_) << "\n";
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       mf::LogTrace(name_) << "Transport matrices in GeV cm U V space " << ipl << "\n";
       mf::LogTrace(name_) << myTransferMatricesGeVcm[ipl];
       mf::LogTrace(name_) << "\n \n";
     }

     mf::LogTrace(name_) << "\n";
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       mf::LogTrace(name_) << "Propagated error matrices in GeV cm U V space" << ipl << "\n";
       mf::LogTrace(name_) << errorMatricesGeVcm[ipl];
       mf::LogTrace(name_) << "\n \n";
     }
  }
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////// 

  // PAY ATTENTION TO THIS !!!!!!
  // NOTE NOTE NOTE: mainDiagonalErrorMatrices is shifted down 1 from errorMatricesGeVcm
  // the latter has a 0 matrix for the 0 plane, while the former doesn't consider the fictional 0 plane at all
  // this now is sort of run around by looping through the trackPlanesHitList array 
  /////////////////////////////////////////////////////////////////////////////////////

  for (int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl) // Fill diagonal errors as well as check which planes were hit.
  { 
     mainDiagonalErrorMatrices.push_back(unmeasured * Eigen::MatrixXd::Identity(5,5)); // Fill infinities into non-measured elements.   
     mainDiagonalErrorMatricesInverse.push_back(Eigen::MatrixXd::Zero(5,5));

     int planeNum = track.trackPlanesHitList.at(ipl);

     if (matrixDebug_)
     { 
        mf::LogTrace(name_) << "Plane: " << planeNum << " Param measured 3 and 4: " << (geaneHitsOnTrackDetails.geaneMeasuredParameters)[3][planeNum] 
                            << " " << (geaneHitsOnTrackDetails.geaneMeasuredParameters)[4][planeNum] << "\n";           
     }

     // removed W diagonal terms. Only using measurement errors
     if (geaneTrackUtils_.isUPlane(planeNum)) // uhit
     {
        if (matrixDebug_) mf::LogTrace(name_) << " Hit U. " << "\n";
        mainDiagonalErrorMatrices[ipl](3,3) = .1*.1*track.UVerrors.at(planeNum)*track.UVerrors.at(planeNum); // grabbing different per plane errors
     }
     else // vhit
     {
        if (matrixDebug_) mf::LogTrace(name_) << " Hit V. " << "\n";
        mainDiagonalErrorMatrices[ipl](4,4) = .1*.1*track.UVerrors.at(planeNum)*track.UVerrors.at(planeNum); // grabbing different per plane errors
     }
  }

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
     mf::LogTrace(name_) << "\n";
     for(int ipl=0;ipl<track.trackNumPlanesHit;ipl++) {
       mf::LogTrace(name_) << "Pre inversion diagonal error matrices " << track.trackPlanesHitList.at(ipl) << "\n";
       mf::LogTrace(name_) << mainDiagonalErrorMatrices[ipl];
       mf::LogTrace(name_) << "\n \n";
     }
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // Multiply transport matrices together. Be careful with Transport matrix = 0 matrix for the starting plane.
  /////////////////////////////////////////////////////////////////////////////////////

  // adjusts first layer/module hit so that transport matrices are properly multiplied against each other
  int firstModuleHit             = int((track.trackFirstPlaneHit-1)/4.);
  int firstPlaneOfFirstModuleHit = firstModuleHit * 4 + 1;
  int stationNumber              = track.island->station;

  // loop through skipped layers, and adjust beginning module of tracking accordingly (eg. if a single hit in the first module was dropped, resulting in the tracking starting a module later than it should)
  if (skipLayers_.size() > 0) 
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

        if(skippedLayerModule < firstModuleHit)
        {
           firstPlaneOfFirstModuleHit = skippedLayerModule * 4 + 1;
           firstModuleHit = skippedLayerModule;
        } 
     }
  }

  for (int ipl = 0; ipl < maxNumPlanes; ++ipl)
  { 
     transportMatrixBegToEnd.push_back(Eigen::MatrixXd::Zero(5,5)); // Initialize matrices.
  }

  // Fill the first plane of first module hit matrix, ignore the previous plane matrices since they are just 0's and won't be included in the sum later on.
  transportMatrixBegToEnd[firstPlaneOfFirstModuleHit] = myTransferMatricesGeVcm[firstPlaneOfFirstModuleHit]; 

  for (int ipl = (firstPlaneOfFirstModuleHit+1); ipl < maxNumPlanes; ++ipl)
  {
     transportMatrixBegToEnd[ipl] = myTransferMatricesGeVcm[ipl]*transportMatrixBegToEnd[ipl-1]; // Multiply previous full transport matrices on the left by the next single gap transport matrix.
  }

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
    mf::LogTrace(name_) << "\n";
    for(int ipl=0;ipl<maxNumPlanes;ipl++) {
        mf::LogTrace(name_) << "transport matrix beg to end Plane " << ipl << "\n";
        mf::LogTrace(name_) << transportMatrixBegToEnd[ipl];
        mf::LogTrace(name_) << "\n \n";
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////

  // Change the way reduction of error matrices is done.
  // Grab the non-infinite error corresponding to the hit and use it as the only non-zero value in the inverse error matrix
  // Inverse error matrices stored in vector mainDiagonalErrorMatricesInverse<Eigen::MatrixXd>
  std::vector<double> reducedErrorAtPlane;
  std::vector<double> reducedErrorAtPlaneInverse;
  for (int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl)
  {
     int planeNum = track.trackPlanesHitList.at(ipl);

     if (geaneTrackUtils_.isUPlane(planeNum)) // uhit
     {
        double errorUHit = mainDiagonalErrorMatrices[ipl](3,3);
        mainDiagonalErrorMatricesInverse[ipl](3,3) = 1./errorUHit;
        reducedErrorAtPlane.push_back(errorUHit);
        reducedErrorAtPlaneInverse.push_back(1./errorUHit);
     }
     else // vhit
     {
        double errorVHit = mainDiagonalErrorMatrices[ipl](4,4); 
        mainDiagonalErrorMatricesInverse[ipl](4,4) = 1./errorVHit;
        reducedErrorAtPlane.push_back(errorVHit);
        reducedErrorAtPlaneInverse.push_back(1./errorVHit);
      }
  }
      
  /////////////////////////////////////////////////////////////////////////////////////      
  // Form the overall 5x5 covariance matrix, equation 33 in the geane manual. (uninverted)
  // Change how covariance total inverse is calculated to equation 27 in geane manual
  //covarianceTotalInv = (combinedTransportMatrices.transpose())*reducedMatrixInverse*combinedTransportMatrices;
  covarianceTotalInv.setZero();
  for (int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl)
  {
     // At each hit plane planenum, add on transport matrix transpose(0, ipl) times reducedErrorAtPlaneInverse[ipl] times trasnport matrix(0, ipl)
     int planeNum        = track.trackPlanesHitList.at(ipl);
     covarianceTotalInv += transportMatrixBegToEnd[planeNum].transpose()*mainDiagonalErrorMatricesInverse[ipl]*transportMatrixBegToEnd[planeNum];
  }

  if (matrixDebug_)
  {
    mf::LogTrace(name_) << "\n" << "covariance matrix out front, should be 5x5: " << "\n" << covarianceTotalInv << "\n";
  } 

  /////////////////////////////////////////////////////////////////////////////////////
  // Convert parameters arrays into Eigen objects and XY measured to UV measured:   
  //====================
  for (int ipl=0;ipl<maxNumPlanes;ipl++) 
  {
     paramMeasuredEigen.push_back(Eigen::VectorXd::Zero(5));
     paramPredictedEigen.push_back(Eigen::VectorXd::Zero(5));
     paramPredictedInUVEigen.push_back(Eigen::VectorXd::Zero(5));

     for(int i=0;i<5;i++) {
        paramMeasuredEigen[ipl](i)  = (geaneHitsOnTrackDetails.geaneMeasuredParameters)[i][ipl];
        paramPredictedEigen[ipl](i) = (geaneHitsOnTrackDetails.geanePredictedParameters)[i][ipl];
     } 
  }

  for (int ipl=0;ipl<maxNumPlanes;ipl++) 
  {
     paramPredictedInUVEigen[ipl] = geaneTrackUtils_.JacobianToUV * paramPredictedEigen[ipl]; 
  } 

  /////////////////////////////////////////////////////////////////////////////////////
  // Now we can calculate the next guess:
  //====================================
  //
  // Calculate the residuals at every step: - Far right side of equation 23, 24, or 26 in the geane manual paper.
  //======================================

  for (int ipl=0;ipl<maxNumPlanes;ipl++) 
  {
     residuals.push_back(Eigen::VectorXd::Zero(5));
     residuals[ipl] = paramMeasuredEigen[ipl]-paramPredictedInUVEigen[ipl];
  }

  // GEVCM Convert residuals from MeV mm to GeV cm to multiply against GeV cm transport and error matrices.
  for (int ipl = 1; ipl < maxNumPlanes; ++ipl)
  {
     convertToGeVcm(residuals[ipl]);
  } 

  if (matrixDebug_)
  {
    mf::LogTrace(name_) << "\n";
    for (int ipl = 1; ipl < maxNumPlanes; ++ipl){
       mf::LogTrace(name_) << "params measured Plane " << ipl << "\n";
       mf::LogTrace(name_) << "Eigen objects: " << "\n" << paramMeasuredEigen[ipl] << std::endl;
       mf::LogTrace(name_) << "\n \n";
    } 

    mf::LogTrace(name_) << "\n";
    for (int ipl = 1; ipl < maxNumPlanes; ++ipl){
       mf::LogTrace(name_) << "params predicted Plane " << ipl << "\n";
       mf::LogTrace(name_) << "Eigen objects: " << "\n" << paramPredictedEigen[ipl] << "\n";
       mf::LogTrace(name_) << " And from XY to UV: " << "\n" << paramPredictedInUVEigen[ipl]; // Also Eigen
       mf::LogTrace(name_) << "\n \n";
    } 

    mf::LogTrace(name_) << "\n";
    for (int ipl = 1; ipl < maxNumPlanes; ++ipl){
       mf::LogTrace(name_) << "residuals Plane " << ipl << "\n";
       mf::LogTrace(name_) << residuals[ipl];
       mf::LogTrace(name_) << "\n \n";
    } 
  }

  // Now multiply out the right side of equation 26 (within the sum).
  //======================================================/
  // Sum over all planes: --- This is the right side of equation 26 or 32 from the geane manual paper. 
  //====================

  //trialTrajectoryAllPlanes = combinedTransportMatrices.transpose()*reducedMatrixInverse*combinedResiduals;
  trialTrajectoryAllPlanes.setZero();
  for (int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl)
  {
     int planeNum = track.trackPlanesHitList.at(ipl);
     trialTrajectoryAllPlanes += transportMatrixBegToEnd[planeNum].transpose()*mainDiagonalErrorMatricesInverse[ipl]*residuals[planeNum];
  }  

  // Calculate the solution:   this is deltaPsiNought from equation 26 or 32 of the geane manual paper
  //=======================
  // Solves the equation Ax = b for x, where A = covarianceTotalInv, b = trialTrajectoryAllPlanes, and x = deltaPsiNought.

  deltaPsiNought = covarianceTotalInv.colPivHouseholderQr().solve(trialTrajectoryAllPlanes); 

  if (matrixDebug_)
  {
     mf::LogTrace(name_) << "\n starting parameter deltas \n";
     for(int i=0;i<5;i++) {
       mf::LogTrace(name_) << "parameter " << i << "\n";
       mf::LogTrace(name_) << deltaPsiNought[i] << "\n"; 
     } 
  }

  track.chi2Planes.clear(); // clear first on each iteration
  double eventChiSquared = 0;
  
  for(int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl)
  {
    int planeNum = track.trackPlanesHitList.at(ipl);
    double planeChiSquared = 0;

    if(geaneTrackUtils_.isUPlane(planeNum)) // uhit 
    {
      planeChiSquared  = residuals[planeNum](3)*reducedErrorAtPlaneInverse[ipl]*residuals[planeNum](3);
      eventChiSquared += planeChiSquared;
    }
    else 
    {
      planeChiSquared  = residuals[planeNum](4)*reducedErrorAtPlaneInverse[ipl]*residuals[planeNum](4);
      eventChiSquared += planeChiSquared;  
    }

    track.chi2Planes.push_back(planeChiSquared);
  }
  
  mf::LogTrace(name_) << "Chi^2 for the track is: " << eventChiSquared << "\n";

  // I can multiply the Jacobian against deltaPsiNought here to immediately convert back to XY for easier use in the main code. Careful with inverses.
  deltaPsiNoughtXY = geaneTrackUtils_.JacobianToUV.inverse() * deltaPsiNought;

  Eigen::VectorXd deltaStartingTrack = Eigen::VectorXd::Zero(5);

  deltaStartingTrack(0) = deltaPsiNoughtXY(0);
  deltaStartingTrack(1) = deltaPsiNoughtXY(1);
  deltaStartingTrack(2) = deltaPsiNoughtXY(2);
  deltaStartingTrack(3) = deltaPsiNoughtXY(3);
  deltaStartingTrack(4) = deltaPsiNoughtXY(4);

  convertToMeVmm(deltaStartingTrack); // Convert back to MeV mm units for the main code to use.

  mf::LogTrace(name_) << "\n starting parameter deltas in XY MeV mm \n";
  for (int i=0;i<5;i++) {
      mf::LogTrace(name_) << "parameter " << i << "\n";
      mf::LogTrace(name_) << deltaStartingTrack[i] << "\n"; 
  } 

  double reconstructedTotalMomentum = 1./(deltaStartingTrack[0]+(1./totalStartMomentum));

  double pvpOVERpup = deltaStartingTrack[1] + startingXMom/startingZMom;
  double pwpOVERpup = deltaStartingTrack[2] + startingYMom/startingZMom;
  
  startingZMom = sqrt( reconstructedTotalMomentum*reconstructedTotalMomentum / (1 + pvpOVERpup*pvpOVERpup + pwpOVERpup*pwpOVERpup) ); // This is only correct for an orthogonal system.
  startingXMom = pvpOVERpup * startingZMom; // This is the new startingZMom.
  startingYMom = pwpOVERpup * startingZMom;

  startingXPos = startingXPos + deltaStartingTrack[3];
  startingYPos = startingYPos + deltaStartingTrack[4];
  startingZPos = startingZPos + 0.; // There are no returned variables that allow for changing the starting z position. 

  // Change the track details to include the newest starting parameters.
  geaneHitsOnTrackDetails.startingGeaneParameters[0] = startingXPos;
  geaneHitsOnTrackDetails.startingGeaneParameters[1] = startingYPos;
  geaneHitsOnTrackDetails.startingGeaneParameters[2] = startingZPos;
  geaneHitsOnTrackDetails.startingGeaneParameters[3] = startingXMom;
  geaneHitsOnTrackDetails.startingGeaneParameters[4] = startingYMom;
  geaneHitsOnTrackDetails.startingGeaneParameters[5] = startingZMom;

  geaneHitsOnTrackDetails.extendedTransportMatrixBegToEnd = transportMatrixBegToEnd;
  geaneHitsOnTrackDetails.covarianceTotalInverse          = covarianceTotalInv; // inverted
  geaneHitsOnTrackDetails.paramPredictedInUVEigen         = paramPredictedInUVEigen;

  track.geaneHits = geaneHitsOnTrackDetails;
  return eventChiSquared;
}

void Geane5x5Fitter::convertToGeVcm(Eigen::VectorXd& myVector){
     myVector(0) = myVector(0) * 1.0e3;
     myVector(1) = myVector(1) * 1.0e0;
     myVector(2) = myVector(2) * 1.0e0;
     myVector(3) = myVector(3) * 1.0e-1;
     myVector(4) = myVector(4) * 1.0e-1;
}

void Geane5x5Fitter::convertToMeVmm(Eigen::VectorXd& myVector){
     myVector(0) = myVector(0) * 1.0e-3; 
     myVector(1) = myVector(1) * 1.0e0; 
     myVector(2) = myVector(2) * 1.0e0;
     myVector(3) = myVector(3) * 1.0e1;
     myVector(4) = myVector(4) * 1.0e1;
}


