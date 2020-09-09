//! header
#include "GeaneFitter.hh"

using namespace gm2strawtracker;
using namespace std;

//! constructor
GeaneFitter::GeaneFitter(fhicl::ParameterSet const & p) 
    : matrixDebug_(p.get<bool>("matrixDebug",false))
    , name_("GeaneFitter")
    , sgeom_()
    , geaneTrackUtils_()
    , numTrackParams(5)
    , skipLayers_(p.get<vector<unsigned int> >("skipLayers",{}))
{
  maxNumPlanes = geaneTrackUtils_.maxNumPlanes;
  wireMeasurementError = 2*sgeom_.outerRadiusOfTheGas/sqrt(12); // this is used for the sequence checking stuff - specifically the hybrid error matrix
}

// the track correlation function
double GeaneFitter::TrackCorrelation(gm2strawtracker::TrackDetailArtRecord & track)
{
  auto geaneHitsOnTrackDetails = track.geaneHits;

  double startingXPos = (geaneHitsOnTrackDetails).startingGeaneParameters[0];
  double startingYPos = (geaneHitsOnTrackDetails).startingGeaneParameters[1];
  double startingZPos = (geaneHitsOnTrackDetails).startingGeaneParameters[2];
  double startingXMom = (geaneHitsOnTrackDetails).startingGeaneParameters[3];
  double startingYMom = (geaneHitsOnTrackDetails).startingGeaneParameters[4];
  double startingZMom = (geaneHitsOnTrackDetails).startingGeaneParameters[5];

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

  Eigen::VectorXd trialTrajectoryAllPlanes;

  std::vector<Eigen::MatrixXd> transportMatrixBegToEnd;
  std::vector<Eigen::MatrixXd> mainDiagonalErrorMatrices;

  // large value for measurement uncertainty of V in U straws, etc. Is used purely as a placeholder for which rows and columns should be removed before multiplications.
  double unmeasured = 1.e300; 

  std::vector<Eigen::MatrixXd> myTransferMatricesGeVcm;
  std::vector<Eigen::MatrixXd> errorMatricesGeVcm;

  // std::vector<Eigen::MatrixXd> errorMatricesBetweenPlanes; // For use with banded inversion.

  /////////////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXd combinedResiduals(track.trackNumPlanesHit); // Reduced to size N x 1
  Eigen::MatrixXd combinedTransportMatrices(track.trackNumPlanesHit,numTrackParams); // Reduced to size N x 5

  // Eigen::MatrixXd combinedTotalErrorCorrelationMatrix(numTrackParams*(maxNumPlanes-1),numTrackParams*(maxNumPlanes-1)); // Full size, 160 x 160
  Eigen::MatrixXd combinedTotalErrorCorrelationMatrix(numTrackParams*track.trackNumPlanesHit,numTrackParams*track.trackNumPlanesHit); // 5N x 5N before being reduced to NxN

  combinedTransportMatrices.setZero();
  combinedTotalErrorCorrelationMatrix.setZero();

  /////////////////////////////////////////////////////////////////////////////////////
  // GeV cm are the better units for the matrix multiplication to work out, otherwise numerical errors can start ocurring with large transport matrices.
  // Commented out because this part was moved to end of errorPropagation method
  // for(int ipl=0;ipl<maxNumPlanes;ipl++) {
  //    myTransferMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));
  //    errorMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));
  //    for(int i=0;i<5;i++) {
  //      for(int j=0; j<5;j++) {
  //        myTransferMatricesGeVcm[ipl](i,j) = (geaneHitsOnTrackDetails.geaneTransportMatrices)[ipl][i][j]*1.; // GEVCM (Don't scale them to MeV mm.)
  //        errorMatricesGeVcm[ipl](i,j)      = (geaneHitsOnTrackDetails.geaneErrorMatrices)[ipl][i][j]*1.; // GEVCM
  //      }
  //    } 
  //  } 

  for (int ipl = 0; ipl < maxNumPlanes; ++ipl) // Use this loop to convert XY transport and error matrices to VU.
  {
     myTransferMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));
     errorMatricesGeVcm.push_back(Eigen::MatrixXd::Zero(5,5));  

     myTransferMatricesGeVcm[ipl] = geaneTrackUtils_.JacobianToUV * geaneHitsOnTrackDetails.geaneTransportMatrices[ipl] * geaneTrackUtils_.JacobianToUV.inverse();
     errorMatricesGeVcm[ipl]      = geaneTrackUtils_.JacobianToUV * geaneHitsOnTrackDetails.geaneErrorMatrices[ipl] * geaneTrackUtils_.JacobianToUV.transpose();
  }

  // for (int ipl = 0; ipl < maxNumPlanes; ++ipl) // Use this loop to convert XY transport and error matrices to VU.
  // {
  //    myTransferMatricesGeVcm[ipl] =  geaneTrackUtils_.JacobianToUV * myTransferMatricesGeVcm[ipl] * geaneTrackUtils_.JacobianToUV.inverse();
  //    errorMatricesGeVcm[ipl]      =  geaneTrackUtils_.JacobianToUV * errorMatricesGeVcm[ipl] * geaneTrackUtils_.JacobianToUV.transpose();
  // }

  // For use with banded inversion if it's ever implemented.
  // errorMatricesBetweenPlanes.push_back(Eigen::MatrixXd::Zero(5,5)); // Put a zero matrix into the zeroth position.
  // for (int ipl = 1; ipl < maxNumPlanes; ++ipl)
  // {
  //    errorMatricesBetweenPlanes.push_back(errorMatricesGeVcm[ipl]-errorMatricesGeVcm[ipl-1]); // Construct error matrices between planes, I think subtraction is okay. 
  // }

  if (matrixDebug_)
  {
     mf::LogTrace(name_) << "\n";
     for(int ipl=0;ipl<maxNumPlanes;ipl++) {
       mf::LogTrace(name_) << "Transport matrices in GeV cm X Y space " << ipl << "\n";
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
  // PAY ATTENTION TO THIS !!!!!!
  // NOTE NOTE NOTE: mainDiagonalErrorMatrices is shifted down 1 from errorMatricesGeVcm
  // the latter has a 0 matrix for the 0 plane, while the former doesn't consider the fictional 0 plane at all
  // this now is sort of run around by looping through the trackPlanesHitList array 
  /////////////////////////////////////////////////////////////////////////////////////

  for (int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl) // Fill diagonal errors as well as check which planes were hit.
  { 
     mainDiagonalErrorMatrices.push_back(unmeasured * Eigen::MatrixXd::Identity(5,5)); // Fill infinities into non-measured elements.     

     int planeNum = track.trackPlanesHitList.at(ipl);

     if (matrixDebug_)
     {
       mf::LogTrace(name_) << "Plane: " << planeNum << " Param measured 3 and 4: " << (geaneHitsOnTrackDetails.geaneMeasuredParameters)[3][planeNum] 
                           << " " << (geaneHitsOnTrackDetails.geaneMeasuredParameters)[4][planeNum] << "\n";
     }

     if (geaneTrackUtils_.isUPlane(planeNum)) // uhit
     {
        if (matrixDebug_) mf::LogTrace(name_) << " Hit U. " << "\n";
        mainDiagonalErrorMatrices[ipl](3,3) = errorMatricesGeVcm[planeNum](3,3) + .1*.1*track.UVerrors.at(planeNum)*track.UVerrors.at(planeNum); // grabbing different per plane errors
     }
     else // vhit
     {
        if (matrixDebug_) mf::LogTrace(name_) << " Hit V. " << "\n";
        mainDiagonalErrorMatrices[ipl](4,4) = errorMatricesGeVcm[planeNum](4,4) + .1*.1*track.UVerrors.at(planeNum)*track.UVerrors.at(planeNum); // grabbing different per plane errors
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

  for (int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl) 
  {
     combinedTotalErrorCorrelationMatrix.block(numTrackParams*(ipl),numTrackParams*(ipl),numTrackParams,numTrackParams) = mainDiagonalErrorMatrices[ipl];
  }


  if (matrixDebug_)
  {
    mf::LogTrace(name_) << "\n" << "Combined 5N x 5N error correlation matrix just diagonals top left block: " << "\n" << combinedTotalErrorCorrelationMatrix.topLeftCorner(15,15) << "\n";
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // Multiply transport matrices together. Be careful with Transport matrix = 0 matrix for the starting plane.
  /////////////////////////////////////////////////////////////////////////////////////

  // adjusts first layer/module hit so that transport matrices are properly multiplied against each other
  int firstModuleHit             = int((track.trackFirstPlaneHit-1)/4.);
  int firstPlaneOfFirstModuleHit = firstModuleHit * 4 + 1;
  int stationNumber              = track.island->station;

  // loop through skipped layers, and adjust beginning module of tracking accordingly (eg. if a single hit in the first module was dropped, resulting in the tracking starting a module later than it should)
  if(skipLayers_.size() > 0) 
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
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////
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

  for (int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl)
  {
     int planeNum = track.trackPlanesHitList.at(ipl);

     // Only grab the row in the transport matrix corresponding to what will multiply against the u or v hit.
     if (geaneTrackUtils_.isUPlane(planeNum)) // uhit
     {
        combinedTransportMatrices.block(ipl,0,1,numTrackParams) = transportMatrixBegToEnd[planeNum].block(3,0,1,numTrackParams); 
     }
     else // vhit
     {
        combinedTransportMatrices.block(ipl,0,1,numTrackParams) = transportMatrixBegToEnd[planeNum].block(4,0,1,numTrackParams);
     }
  }

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
     mf::LogTrace(name_) << "Combined N x 5 transport matrix: " << "\n" << combinedTransportMatrices << "\n";
  }
  /////////////////////////////////////////////////////////////////////////////////////

  for (int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl)
  {
     for (int index = ipl+1; index < int(track.trackPlanesHitList.size()); ++index)
     {
        int planeNum = track.trackPlanesHitList.at(ipl);
        Eigen::MatrixXd tempMatrix = (transportMatrixBegToEnd[track.trackPlanesHitList.at(index)]*(transportMatrixBegToEnd[planeNum].inverse())*errorMatricesGeVcm[planeNum]);
        Eigen::MatrixXd topDiagonalMatrixFill = Eigen::MatrixXd::Zero(5,5);

        topDiagonalMatrixFill.bottomRightCorner<2,2>() = tempMatrix.bottomRightCorner<2,2>();
        combinedTotalErrorCorrelationMatrix.block(numTrackParams*(ipl),numTrackParams*(index),numTrackParams,numTrackParams) = topDiagonalMatrixFill;
        combinedTotalErrorCorrelationMatrix.block(numTrackParams*(index),numTrackParams*(ipl),numTrackParams,numTrackParams) = topDiagonalMatrixFill.transpose();
      }
  }

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
     mf::LogTrace(name_) << "\n" << "Combined 5N x 5N error correlation matrix top left block: " << "\n" << combinedTotalErrorCorrelationMatrix.topLeftCorner(15,15) << "\n";
     mf::LogTrace(name_) << "\n" << "Combined 5N x 5N error correlation matrix bottom right block: " << "\n" << combinedTotalErrorCorrelationMatrix.bottomRightCorner(15,15) << "\n";
     mf::LogTrace(name_) << "\n" << "Number of rows and cols: " << combinedTotalErrorCorrelationMatrix.rows() << " " << combinedTotalErrorCorrelationMatrix.cols() << "\n";
  }     
  /////////////////////////////////////////////////////////////////////////////////////

  // Reduce matrix to make inversions faster.
  // will have to include a -1 if numPlanesHit ever includes some 0 plane - which it never should at this point
  Eigen::MatrixXd reducedMatrix((track.trackNumPlanesHit),(track.trackNumPlanesHit)); 

  /////////////////////////////////////////////////////////////////////////////////////
  int i = 0;
  while(i < int(track.trackPlanesHitList.size())) // remove rows and cols from big error correlation matrix
  {
      if (combinedTotalErrorCorrelationMatrix(i,i) == unmeasured){
          
        combinedTotalErrorCorrelationMatrix.block(i, 0, combinedTotalErrorCorrelationMatrix.rows()-(i+1), combinedTotalErrorCorrelationMatrix.cols()) = 
        combinedTotalErrorCorrelationMatrix.block(i+1, 0, combinedTotalErrorCorrelationMatrix.rows()-(i+1), combinedTotalErrorCorrelationMatrix.cols()); // remove row
  
        combinedTotalErrorCorrelationMatrix.block(0, i, combinedTotalErrorCorrelationMatrix.rows(), combinedTotalErrorCorrelationMatrix.cols()-(i+1)) = 
        combinedTotalErrorCorrelationMatrix.block(0, i+1, combinedTotalErrorCorrelationMatrix.rows(), combinedTotalErrorCorrelationMatrix.cols()-(i+1)); // remove col

        combinedTotalErrorCorrelationMatrix.conservativeResize(combinedTotalErrorCorrelationMatrix.rows()-1,combinedTotalErrorCorrelationMatrix.cols()-1);
      }
      else i++;
  }

  // There is some error that crops up occasionally where there is one too many rows/cols at the very end that doesn't get properly deleted. 
  // Why it doesn't appear for more last U hits I'm not sure. I can resize here just in case to get rid of those occurences.
  combinedTotalErrorCorrelationMatrix.conservativeResize(int(track.trackPlanesHitList.size()),int(track.trackPlanesHitList.size())); 

  // Just feed it into a new variable since that's what I already use right now
  reducedMatrix = combinedTotalErrorCorrelationMatrix; 

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
    mf::LogTrace(name_) << "\n" << "Combined error matrix size before placing in reduced matrix: " << "\n" << combinedTotalErrorCorrelationMatrix.rows() << " " << combinedTotalErrorCorrelationMatrix.cols() << "\n";
    mf::LogTrace(name_) << "\n" << "Reduced matrix: " << "\n" << reducedMatrix << "\n";
  } 
  /////////////////////////////////////////////////////////////////////////////////////
 
  Eigen::MatrixXd reducedMatrixInverse = reducedMatrix.inverse();

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
    mf::LogTrace(name_) << "\n" << "Reduced matrix inverse: " << "\n" << reducedMatrixInverse << "\n";
  } 
  /////////////////////////////////////////////////////////////////////////////////////
 
  // Form the overall 5x5 covariance matrix, equation 33 in the geane manual. (uninverted)
  covarianceTotalInv = (combinedTransportMatrices.transpose())*reducedMatrixInverse*combinedTransportMatrices;

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
    mf::LogTrace(name_) << "\n" << "covariance matrix out front, should be 5x5: " << "\n" << covarianceTotalInv << "\n";
  } 
  /////////////////////////////////////////////////////////////////////////////////////


  // Convert parameters arrays into Eigen objects and XY measured to UV measured:   
  for (int ipl=0;ipl<maxNumPlanes;ipl++) {
     paramMeasuredEigen.push_back(Eigen::VectorXd::Zero(5));
     paramPredictedEigen.push_back(Eigen::VectorXd::Zero(5));
     paramPredictedInUVEigen.push_back(Eigen::VectorXd::Zero(5));

     for (int i=0;i<5;i++) {
        paramMeasuredEigen[ipl](i)  = (geaneHitsOnTrackDetails.geaneMeasuredParameters)[i][ipl];
        paramPredictedEigen[ipl](i) = (geaneHitsOnTrackDetails.geanePredictedParameters)[i][ipl];
     } 
  }

  for (int ipl=0;ipl<maxNumPlanes;ipl++) {
     paramPredictedInUVEigen[ipl] = geaneTrackUtils_.JacobianToUV * paramPredictedEigen[ipl]; 
  } 

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Now we can calculate the next guess:
  // ====================================
  //
  // Calculate the residuals at every step: - Far right side of equation 23, 24, or 26 in the geane manual paper.
  // ======================================
  //           
  /////////////////////////////////////////////////////////////////////////////////////

  for (int ipl=0;ipl<maxNumPlanes;ipl++) {
     residuals.push_back(Eigen::VectorXd::Zero(5));
     residuals[ipl] = paramMeasuredEigen[ipl]-paramPredictedInUVEigen[ipl];
  }

  // GEVCM Convert residuals from MeV mm to GeV cm to multiply against GeV cm transport and error matrices.
  for (int ipl = 1; ipl < maxNumPlanes; ++ipl)
  {
     convertToGeVcm(residuals[ipl]);
  } 

  /////////////////////////////////////////////////////////////////////////////////////
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
  /////////////////////////////////////////////////////////////////////////////////////
 
  // Now multiply out the right side of equation 26 (within the sum).
  for (int ipl = 0; ipl < int(track.trackPlanesHitList.size()); ++ipl)
  {
     int planeNum = track.trackPlanesHitList.at(ipl);
     if (geaneTrackUtils_.isUPlane(planeNum)) // uhit
     {
       combinedResiduals.block(ipl,0,1,1) = residuals[planeNum].block(3,0,1,1);
     }
     else
     {
       combinedResiduals.block(ipl,0,1,1) = residuals[planeNum].block(4,0,1,1);
     }
  }

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
    mf::LogTrace(name_) << "Combined residuals: " << "\n" << combinedResiduals << "\n";
  } 
  /////////////////////////////////////////////////////////////////////////////////////


  // Sum over all planes: --- This is the right side of equation 26 or 32 from the geane manual paper.
  trialTrajectoryAllPlanes = combinedTransportMatrices.transpose()*reducedMatrixInverse*combinedResiduals;  

  // Calculate the solution:   this is deltaPsiNought from equation 26 or 32 of the geane manual paper
  // Solves the equation Ax = b for x, where A = covarianceTotalInv, b = trialTrajectoryAllPlanes, and x = deltaPsiNought.
  deltaPsiNought = covarianceTotalInv.colPivHouseholderQr().solve(trialTrajectoryAllPlanes); 

  /////////////////////////////////////////////////////////////////////////////////////
  if (matrixDebug_)
  {
     mf::LogTrace(name_) << "\n starting parameter deltas \n";
     for (int i=0;i<5;i++) {
        mf::LogTrace(name_) << "parameter " << i << "\n";
        mf::LogTrace(name_) << deltaPsiNought[i] << "\n"; 
     } 
  }
  /////////////////////////////////////////////////////////////////////////////////////
 
  double eventChiSquared = combinedResiduals.transpose()*reducedMatrixInverse*combinedResiduals;

  // calculate individual plane chi2s - does not include correlations - and so sum of per plane chi2s does not equal total chi2
  track.chi2Planes.clear(); // clear first on each iteration

  for (int i = 0; i < track.trackNumPlanesHit; ++i)
  {
     track.chi2Planes.push_back( combinedResiduals.transpose().segment(i,1) * reducedMatrixInverse(i,i) * combinedResiduals.segment(i,1) );
  }

  mf::LogTrace(name_) << "Chi^2 for the track is: " << eventChiSquared << "\n";

  // I can multiply the Jacobian against deltaPsiNought here to immediately convert back to XY for easire use in the main code. Careful with inverses.
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
  startingZPos = startingZPos + 0.; // There are no returned variables that allow for changing the starting Z position. 

  // Change the track details to include the newest starting parameters.
  geaneHitsOnTrackDetails.startingGeaneParameters[0] = startingXPos;
  geaneHitsOnTrackDetails.startingGeaneParameters[1] = startingYPos;
  geaneHitsOnTrackDetails.startingGeaneParameters[2] = startingZPos;
  geaneHitsOnTrackDetails.startingGeaneParameters[3] = startingXMom;
  geaneHitsOnTrackDetails.startingGeaneParameters[4] = startingYMom;
  geaneHitsOnTrackDetails.startingGeaneParameters[5] = startingZMom;

  geaneHitsOnTrackDetails.extendedTransportMatrixBegToEnd            = transportMatrixBegToEnd;
  geaneHitsOnTrackDetails.extendedCombinedTransportMatricesTranspose = combinedTransportMatrices.transpose(); // Have to later transpose this since this is saved untransposed - might want to fix.

  geaneHitsOnTrackDetails.extendedReducedMatrix        = reducedMatrix;
  geaneHitsOnTrackDetails.extendedReducedMatrixInverse = reducedMatrixInverse;
  geaneHitsOnTrackDetails.covarianceTotalInverse       = covarianceTotalInv; // inverted
  geaneHitsOnTrackDetails.paramPredictedInUVEigen      = paramPredictedInUVEigen;

  track.geaneHits = geaneHitsOnTrackDetails;
  return eventChiSquared;
}

void GeaneFitter::makeHybridErrMat(gm2strawtracker::TrackDetailArtRecord & trackDetails, bool Uhybrid)
{
  // this matrix has the wire measurement errors on all diagonals before this method is called (both times)
  trackDetails.geaneHits.extendedModifiedReducedMatrix = trackDetails.geaneHits.extendedReducedMatrix; 

  for (int i = 0; i < trackDetails.geaneHits.extendedModifiedReducedMatrix.rows(); ++i) // matrix rows and trackPlanesHitList should have the same size
  {
     int planeNum = trackDetails.trackPlanesHitList.at(i);
     if ((geaneTrackUtils_.isUPlane(planeNum) && Uhybrid) // uplane and uhybrid
         || (!geaneTrackUtils_.isUPlane(planeNum) && !Uhybrid)) // vplane and vhybrid 
     {
        // replace wire/straw errors with hit errors
        trackDetails.geaneHits.extendedModifiedReducedMatrix(i,i) = 
                trackDetails.geaneHits.extendedModifiedReducedMatrix(i,i) - .1*wireMeasurementError*.1*wireMeasurementError + .1*trackDetails.UVerrors.at(planeNum)*.1*trackDetails.UVerrors.at(planeNum); 
     }
  }
  trackDetails.geaneHits.extendedModifiedReducedMatrix = trackDetails.geaneHits.extendedModifiedReducedMatrix.inverse();
}

double GeaneFitter::sequenceChecking(gm2strawtracker::TrackDetailArtRecord & trackDetails)
{
  std::vector<Eigen::VectorXd> newParamMeasuredEigen;
  std::vector<Eigen::VectorXd> newParamPredictedEigen;

  std::vector<Eigen::VectorXd> newResiduals;
  Eigen::VectorXd newCombinedResiduals(trackDetails.trackNumPlanesHit);

  Eigen::VectorXd newTrialTrajectoryAllPlanes;

  for (int ipl=0;ipl<maxNumPlanes;ipl++) {
     newParamMeasuredEigen.push_back(Eigen::VectorXd::Zero(5));

     for (int i=0;i<5;i++) {
        newParamMeasuredEigen[ipl](i) = (trackDetails.geaneHits.geaneMeasuredParameters)[i][ipl];
     } 
  }

  for (int ipl=0;ipl<maxNumPlanes;ipl++) {
     newResiduals.push_back(Eigen::VectorXd::Zero(5));
     newResiduals[ipl] = newParamMeasuredEigen[ipl]-trackDetails.geaneHits.paramPredictedInUVEigen[ipl];
  }

  for (int ipl = 1; ipl < maxNumPlanes; ++ipl)
  {
     convertToGeVcm(newResiduals[ipl]);
  } 

  for (int ipl = 0; ipl < int(trackDetails.trackPlanesHitList.size()); ++ipl)
  {
     int planeNum = trackDetails.trackPlanesHitList.at(ipl);

     if (geaneTrackUtils_.isUPlane(planeNum)) // uhit
     {
        newCombinedResiduals.block(ipl,0,1,1) = newResiduals[planeNum].block(3,0,1,1);
     }
     else // vhit
     {
       newCombinedResiduals.block(ipl,0,1,1) = newResiduals[planeNum].block(4,0,1,1);
     }
  }

  /////////////////////////////////////////////////////////////////////////////////////
  // seems that I could stop here for a faster sequence check - seems to get really close for some events and not for others - study more
  // double sequenceChiSquaredtest = newCombinedResiduals.transpose()*(*geaneHitsOnTrackDetails).extendedModifiedReducedMatrix*newCombinedResiduals;
  // return sequenceChiSquaredtest;
  /////////////////////////////////////////////////////////////////////////////////////

  newTrialTrajectoryAllPlanes        = trackDetails.geaneHits.extendedCombinedTransportMatricesTranspose*(trackDetails.geaneHits.extendedReducedMatrixInverse)*newCombinedResiduals;  
  Eigen::VectorXd deltaStartingTrack = (trackDetails.geaneHits.covarianceTotalInverse).colPivHouseholderQr().solve(newTrialTrajectoryAllPlanes);

  std::vector<Eigen::VectorXd> paramPredictedChanges;
  for (int ipl=0;ipl<maxNumPlanes;ipl++) {
     paramPredictedChanges.push_back(Eigen::VectorXd::Zero(5));
     paramPredictedChanges[ipl] = (trackDetails.geaneHits.extendedTransportMatrixBegToEnd)[ipl]*(deltaStartingTrack);
  }

  std::vector<Eigen::VectorXd> paramPredictedGeVcm;
  for (int ipl = 0; ipl < maxNumPlanes; ++ipl)
  {
     paramPredictedGeVcm.push_back(Eigen::VectorXd::Zero(5));
     paramPredictedGeVcm[ipl] = (trackDetails.geaneHits.paramPredictedInUVEigen)[ipl];
     convertToGeVcm(paramPredictedGeVcm[ipl]); // need these in GeV cm
  }

  for (int ipl = 0; ipl < maxNumPlanes; ++ipl)
  {
     newParamPredictedEigen.push_back(Eigen::VectorXd::Zero(5));
     newParamPredictedEigen[ipl] = paramPredictedGeVcm[ipl]+paramPredictedChanges[ipl];
  }

  for (int ipl = 1; ipl < maxNumPlanes; ++ipl)
  {
     convertToGeVcm(newParamMeasuredEigen[ipl]);
  } 

  newResiduals.clear(); // remake new residuals
  for (int ipl=0;ipl<maxNumPlanes;ipl++) {
     newResiduals.push_back(Eigen::VectorXd::Zero(5));
     newResiduals[ipl] = newParamMeasuredEigen[ipl]-newParamPredictedEigen[ipl];
  }

  for (int ipl = 0; ipl < int(trackDetails.trackPlanesHitList.size()); ++ipl)
  {
     int planeNum = trackDetails.trackPlanesHitList.at(ipl);

     if (geaneTrackUtils_.isUPlane(planeNum)) // uhit
     {
        newCombinedResiduals.block(ipl,0,1,1) = newResiduals[planeNum].block(3,0,1,1); // make residuals N x 1 sizes
     }
     else // vhit
     {
        newCombinedResiduals.block(ipl,0,1,1) = newResiduals[planeNum].block(4,0,1,1);
     }
  }

  double sequenceChiSquared = newCombinedResiduals.transpose()*trackDetails.geaneHits.extendedModifiedReducedMatrix*newCombinedResiduals; // The extended modified reduced matrix here is already inverted.
  return sequenceChiSquared;
}

void GeaneFitter::convertToGeVcm(Eigen::VectorXd& myVector){
     myVector(0) = myVector(0) * 1.0e3;
     myVector(1) = myVector(1) * 1.0e0;
     myVector(2) = myVector(2) * 1.0e0;
     myVector(3) = myVector(3) * 1.0e-1;
     myVector(4) = myVector(4) * 1.0e-1;
}

void GeaneFitter::convertToMeVmm(Eigen::VectorXd& myVector){
     myVector(0) = myVector(0) * 1.0e-3; 
     myVector(1) = myVector(1) * 1.0e0; 
     myVector(2) = myVector(2) * 1.0e0;
     myVector(3) = myVector(3) * 1.0e1;
     myVector(4) = myVector(4) * 1.0e1;
}


