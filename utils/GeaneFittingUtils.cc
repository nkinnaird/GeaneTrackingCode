#include "GeaneFittingUtils.hh"

namespace gm2strawtracker {

GeaneFittingUtils::GeaneFittingUtils(fhicl::ParameterSet const & p) 
  : geaneParamUtils_(p)
  , geaneTrackUtils_()
  , geaneLRUtils_(p)
  , matrixCalculations_(p)
  , reducedMatrixCalculations_(p)
  , name_( "GeaneFittingUtils" )
  , fitMode_(p.get<string>("fitMode","badFitMode"))
  , fitterType_(p.get<string>("fitterType"))
  , lockLowDCAs_(p.get<double>("lockLowDCAs",0))
  , useNodes_(p.get<bool>("useNodes"))
  , numPassesWireFit_(p.get<int>("numPassesWireFit"))
  , numPassesSeqFit_(p.get<int>("numPassesSeqFit"))
  , convergenceCriteria_(p.get<double>("convergenceCriteria"))
  , useTangentLR_(p.get<bool>("useTangentLR"))
{}

//----------------------
// perform true leff/right fit
//----------------------
int GeaneFittingUtils::truthLRFit(gm2strawtracker::TrackDetailArtRecord & trackArtRecord)
{
  int passesToPerform = 10;
  trackArtRecord.uncorrErrors = trackArtRecord.dcaErrors; // if fitting truth then set straw errors to digit dca errors
  bool updateLR = false;

  return fittingLoop(trackArtRecord, passesToPerform, updateLR);
}


//----------------------
// perform wire fit
//----------------------
int GeaneFittingUtils::wireFit(gm2strawtracker::TrackDetailArtRecord & trackArtRecord)
{
  int passesToPerform = 10;
  bool updateLR = false;

  return fittingLoop(trackArtRecord, passesToPerform, updateLR);
}


//------------------------
// peform main fit
//------------------------
int GeaneFittingUtils::mainFit(gm2strawtracker::TrackDetailArtRecord & trackArtRecord)
{
  int passesToPerform = numPassesWireFit_;
  bool updateLR = false;

  // first wire fit
  int firstWireFit = fittingLoop(trackArtRecord, passesToPerform, updateLR);
  if(firstWireFit != 0) return firstWireFit;

  // updateLRFromFit will lock the centers if dcas are small
  if (lockLowDCAs_) {
     geaneLRUtils_.lockSmallDCAErrors(trackArtRecord); // updates dca errors
  }

  // main fit
  passesToPerform = 10;
  updateLR = true;

  trackArtRecord.uncorrErrors = trackArtRecord.dcaErrors;
  return fittingLoop(trackArtRecord, passesToPerform, updateLR);
}


//------------------------
// peform full sequence fit
//------------------------
int GeaneFittingUtils::fullSeqFit(gm2strawtracker::TrackDetailArtRecord & trackArtRecord)
{
  int passesToPerform = numPassesWireFit_;
  bool updateLR = false;

  //first wire fit
  int firstWireFit = fittingLoop(trackArtRecord, passesToPerform, updateLR);
  if(firstWireFit != 0) return firstWireFit;

  passesToPerform = numPassesSeqFit_;
  trackArtRecord.uncorrErrors = trackArtRecord.dcaErrors; // set errors to digit errors

  if(lockLowDCAs_) geaneLRUtils_.lockSmallDCAErrors(trackArtRecord);

  if(useNodes_) geaneLRUtils_.fillLRFromNodes(trackArtRecord);
  else geaneLRUtils_.setUnknownSides(trackArtRecord);

  /////////////////////////////////////////////////////////////////////////////////////
  // fill with geom if desired -- this will overwrite the unknown sides with geom LR if angle condition is met
  // maybe put an if statement with a fcl parameter for whether to turn this on or not? for now default to on

  trackArtRecord.geaneHits.geaneHitSides = trackArtRecord.geaneHits.inputHitSides; // set hit sides first since method below uses that as input
  geaneLRUtils_.fillLRFromGeom(trackArtRecord);
  if (useTangentLR_) LRUtils_.fillLRFromGeomAndTangent(trackArtRecord); // additional check that tangent matches geom
  trackArtRecord.geaneHits.inputHitSides = trackArtRecord.geaneHits.geaneHitSides; // rewrite input sides since that's what's used below


  /////////////////////////////////////////////////////////////////////////////////////

  std::vector<int> UunknownsInThisEvent = geaneLRUtils_.getUVUnkowns(trackArtRecord, true);
  std::vector<int> VunknownsInThisEvent = geaneLRUtils_.getUVUnkowns(trackArtRecord, false);

  std::vector<std::pair<double, int> > topUSequences = geaneLRUtils_.getTopSequences(trackArtRecord, true,  int(UunknownsInThisEvent.size()));
  std::vector<std::pair<double, int> > topVSequences = geaneLRUtils_.getTopSequences(trackArtRecord, false, int(VunknownsInThisEvent.size()));

  std::vector<std::pair<double, gm2strawtracker::TrackDetailArtRecord> > topTotalSequences; // best sequences with corresponding chi2, U sequence, V sequence and TrackDetailArtRecord

  /////////////////////////////////////////////////////////////////////////////////////

  // loop here over top best 4 or so chi squared tracks in U and V and their combinations with their associated sequences, want to then save the best chi squared at the end
  for (int k = 0; k <= 4; ++k)
  {
    int i = 0;
    int j = k;

    while(i+j <= k && j>=0)
    {
        if(i >= int(topUSequences.size()) || j >= int(topVSequences.size())){ 
          j--;
          i++;
          continue;
        }

        gm2strawtracker::TrackDetailArtRecord mainFitGeaneRecord = trackArtRecord; // this geaneArtRecord should be good to go as a base for the main fitting (starting params, wire positions, etc)

        std::vector<gm2strawtracker::GeaneHitSide> fullSides(geaneTrackUtils_.maxNumPlanes, gNA_side); // will be filled into TrackDetailArtRecord
        fullSides = trackArtRecord.geaneHits.inputHitSides;

        geaneLRUtils_.getSequenceSides(UunknownsInThisEvent, topUSequences.at(i).second, &fullSides); // updates fullSides based on sequences
        geaneLRUtils_.getSequenceSides(VunknownsInThisEvent, topVSequences.at(j).second, &fullSides);

        mainFitGeaneRecord.geaneHits.geaneHitSides = fullSides;

        /////////////////////////////////////////////////////////////////////////////////////
        // perform fitting on sequence

        int fitReturn = fittingLoop(mainFitGeaneRecord, passesToPerform, updateLR); // loop through different methods in fitting
        if(fitReturn != 0){
            j--;
            i++;        
            continue; // continue to next sequence but don't drop the event
        } 

        // fitReturn = geaneFittingUtils_.checkExtraneousFailureModes(&mainFitGeaneRecord); // might be necessary at some point to check the failure modes per sequence instead of just once for the minimum chi2 track
        // if(fitReturn != 0 ) return fitReturn;

        mf::LogTrace(name_) << "\n End of seq i j: " << i << " " << j << "\n\n";

        topTotalSequences.emplace_back(mainFitGeaneRecord.chi2, mainFitGeaneRecord);

        j--;
        i++;

        std::sort(topTotalSequences.begin(), topTotalSequences.end(), [](auto &left, auto &right) {
                 return left.first <= right.first;
                 });

        if (topTotalSequences.size() >= 3 && (topTotalSequences.at(topTotalSequences.size()-1).first - topTotalSequences.at(topTotalSequences.size()-2).first) > 10){
           mf::LogInfo(name_) << "Large jump in chi2s - break \n";
           goto chi2jump;
        }

    } // end while
  } // end loop over best chi squareds

  /////////////////////////////////////////////////////////////////////////////////////
  chi2jump:
  /////////////////////////////////////////////////////////////////////////////////////

  if(topTotalSequences.size() == 0 ) return 10; // if no sequences work return here, geaneArtRecord stays as the wire fit but won't get saved anywhere

  std::sort(topTotalSequences.begin(), topTotalSequences.end(), [](auto &left, auto &right) {
        return left.first <= right.first;
  });

  (trackArtRecord) = topTotalSequences.at(0).second; // fill back into the primary TrackDetailArtRecord that will be passed out and filled into the art file (through address)

  // for (int i = 0; i < int(topTotalSequences.size()); ++i)
  // {
  //   G4cout << "Top sequences i:  " << i << " chi2: " << std::get<0>(topTotalSequences.at(i)) << G4endl;
  // }

  return 0; // return a fitted track - might be a poor fit but didn't fail

}

//----------------------
// fitting loop
//----------------------
int GeaneFittingUtils::fittingLoop(gm2strawtracker::TrackDetailArtRecord & trackArtRecord, int maxNumPasses, bool updateLR)
{
  mf::LogInfo info(name_); 

  double eventChiSquared = 0;
  double prevChiSquared  = 0;
    
  for (int passNum = 1; passNum <= maxNumPasses; ++passNum)
  {

     mf::LogTrace(name_) << "\n" << "" << fitMode_ << " pass number: " << passNum << "\n";
     int errProp = geaneParamUtils_.errorProp(trackArtRecord);

     if (errProp != 0){
       info << "\n" << "Error propagation geant routines failed somehow. failureMode: " << errProp << "\n";
       return errProp;
     }

     if(updateLR == true) geaneLRUtils_.fillLRFromFit(trackArtRecord); // this helps the wire->mainFit converge for more tracks

     geaneParamUtils_.calcMeasuredParams(trackArtRecord); // this sets geaneMeasuredParameters and updates errors based on particle momentum at each plane

     if (fitterType_ == "fullFitter") 
     {
        eventChiSquared = matrixCalculations_.TrackCorrelation(trackArtRecord); // Straws and StrawTrackerCadMesh must be loaded in fcl file
     }
     else if(fitterType_ == "fastFitter") {
        eventChiSquared = reducedMatrixCalculations_.TrackCorrelation(trackArtRecord); // No Straws or StrawTrackerCadMesh loaded in fcl file
     }

     if (!(eventChiSquared > 0)){
       info << "\n" << fitMode_ << " returned negative or nan Chi^2: " << eventChiSquared << "\n";
       return 1; 
     }

     if (passNum > 1 && (eventChiSquared > prevChiSquared+10) ){
       info << "\n" << fitMode_ << " fit is not converging. Current chi2: " << eventChiSquared << " prev chi2: " << prevChiSquared << "\n";
       return 2;
     }

     // fill track details - will be overwritten later with further fitting which is good
     trackArtRecord.numIterations = passNum;
     trackArtRecord.chi2          = eventChiSquared;
     trackArtRecord.dof           = trackArtRecord.trackNumPlanesHit - 5; // dof = N-4 with no field and N-5 with a field
     trackArtRecord.chi2DoF       = trackArtRecord.chi2 /(trackArtRecord.trackNumPlanesHit-5);
     trackArtRecord.pValue        = TMath::Prob(eventChiSquared,trackArtRecord.trackNumPlanesHit-5);

     // the rest of the track details should be filled properly already
     mf::LogTrace(name_) << "\n" << fitMode_ << " end of pass: " << passNum  << "\n" << "\n";

     if (passNum > 1 && abs(eventChiSquared-prevChiSquared) < convergenceCriteria_ ){
        info << "\n" << fitMode_ << " converged " << " on pass: " << passNum << "\n";
        return 0;
     }

     prevChiSquared = eventChiSquared; // fill for next pass

  } // end of for loop

  info << "\n" << fitMode_ << " fittingLoop reached end without converging to .1 in chi2 or failing. " << "\n";
  return 0;
}

//---------------------------------------------------------------------------------------------------------------------
// method for checking weird, different failure modes, can expand, reduce, or remove in the future if desired
//---------------------------------------------------------------------------------------------------------------------
int GeaneFittingUtils::checkExtraneousFailureModes(gm2strawtracker::TrackDetailArtRecord & trackFitDetails)
{

  /////////////////////////////////////////////////////////////////////////////////////
  // Before saw a very small number of events whose x state positions and x targets did not match up to exactly or near exactly zero
  // This doesn't make much since because the target and state x positions are set from the same source
  // Supposedly geant was messing up somewhere, or transforming the wire centers was off
  // Deleting the failure mode but leaving the method for now
  /////////////////////////////////////////////////////////////////////////////////////
 
  return 0;
}


}//namespace
