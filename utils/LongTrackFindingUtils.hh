//!----------------------------------
//!  This class provides useful functions for merging seeds to form track candidates
//!----------------------------------

#ifndef LONGTRACKFINDINGUTILS_HH
#define LONGTRACKFINDINGUTILS_HH

// art includes
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// geometry
#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"

// coordinate system includes
#include "gm2geom/coordSystems/CoordSystem.hh"
#include "gm2geom/coordSystems/CoordSystemsStoreData.hh"
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh"  

// data product includes
#include "gm2dataproducts/strawtracker/StrawSeedArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackCandidateArtRecord.hh"

//For DebugPlotting
#include "gm2util/common/RootManager.hh"
#include "TEllipse.h"

// local
#include "gm2tracker/fitter/LeastSquaresFitter.hh"
#include "gm2tracker/fitter/SimpleCircleFitter.hh"  
#include "gm2tracker/utils/StrawObjectSorter.hh"
#include "gm2tracker/utils/TrackFindingUtils.hh"

#include "boost/format.hpp"

namespace gm2strawtracker {

   class LongTrackFindingUtils {

      public :

       // default constructor
       LongTrackFindingUtils();

       // constructor to make plots
       LongTrackFindingUtils(bool makePlots);

       // DebugPlottingdefault destructor
       ~LongTrackFindingUtils();

       // merge seeds to form track candidates -- this function uses the LSQ minimizer to determine if the merging should occur
       TrackCandidateArtRecordCollection findTrackCandidates( StrawSeedPtrCollection& seeds ); 

       // reconstructed the track's momentum based on the curvature
       void getTrackMomentum( TrackCandidateArtRecord& track );

       inline void setWireVectors( std::vector<gm2geom::CoordSystem3Vector>& wireCentreInWorld, std::vector<gm2geom::CoordSystem3Vector>& wireCentreInModule, std::vector<gm2geom::CoordSystem3Vector>& wireCentreInStation ){
	 wireCentreInWorld_ = wireCentreInWorld;
	 wireCentreInModule_ = wireCentreInModule;
	 wireCentreInStation_ = wireCentreInStation;
       }

       void setFhiclCuts (fhicl::ParameterSet const & p ); 

      private :

       // merge seeds to form track candidates -- this function uses the LSQ minimizer to determine the quality of the track candidate
       TrackCandidateArtRecordCollection mergeSeedsToFormTrackCandidates( StrawSeedPtrCollection& seeds ); 

       // organize the seeds
       void organizeSeedsByModule( StrawSeedPtrCollection& seeds, std::vector< StrawSeedPtrCollection>& seedsByModule, std::vector< int >& indices );

       // select the most probable neighboring seed
       void getNeighboringSeeds( StrawSeedPtrCollection& seeds, art::Ptr<StrawSeedArtRecord>& iseed, StrawSeedPtrCollection& oseeds,
                                 double avgTimeTrackCandidate, bool isManyCloseSeedsInSameModule );

       // checks if two seeds are neighbors
       bool neighboringSeeds( art::Ptr<StrawSeedArtRecord>& seed1, art::Ptr<StrawSeedArtRecord>& seed2 );

       // checks if two seeds share the same seed(s)
       bool trackSeedsOverlaps( StrawSeedPtrCollection track1, StrawSeedPtrCollection track2 );

       // check if the track candidates are neighbor
       bool trackCandidatesNeighbors( TrackCandidateArtRecord& refCandidate, TrackCandidateArtRecord& candidate );

       // check if there is a gap between the neighboring track candidates
       bool isGapBetweenTrackCandidates( TrackCandidateArtRecord& refCandidate, TrackCandidateArtRecord& candidate, StrawSeedPtrCollection& unusedSeeds );

       // check if a seed neighbors the track candidate at the front or back
       bool trackCandidateSeedNeighbors( TrackCandidateArtRecord& track, art::Ptr<StrawSeedArtRecord>& seed, bool mergeBrokenTracks = false );

       // check if a seed overlaps the track candidate 
       bool trackCandidateSeedOverlaps( TrackCandidateArtRecord& track, art::Ptr<StrawSeedArtRecord>& seed );

       // check if a seed and the track candidate share a cluster
       bool trackCandidateSeedShareCluster( TrackCandidateArtRecord& track, art::Ptr<StrawSeedArtRecord>& seed );

       // determine if the digits on the track candidate are reasonable
       bool checkDigitsOnTrackCandidate( TrackCandidateArtRecord& track );

       // determine if the seeds are overlapping the same module
       bool overlappingSeeds( StrawSeedPtrCollection& seeds );

       // determine if the seed over in module and wire regions
       bool overlappingSeeds( art::Ptr<StrawSeedArtRecord>& refSeed, art::Ptr<StrawSeedArtRecord>& neighSeed );

       // remove seeds that are already living on a track candidate
       void checkSeedIsOnTrackCandidate( StrawSeedPtrCollection& seeds, TrackCandidateArtRecordCollection& candidates );

       // get the averge time difference between the seeds on the track candidate
       double getTimeDifferenceBetweenSeeds( StrawSeedPtrCollection& seeds );

       // call the least squares minimizer fitter
       void callLSQFitter( TrackCandidateArtRecord& candidate );

       // call the simple circle fast fitter
       void callSimpleCircleFitter( TrackCandidateArtRecord& candidate );
     
       // for DebugPlots
       TH1F* strawPlot(StrawDigitPtrCollection digits, TH1F* bgPlot, int fillColor = 2 );

       // helper tools
       TrackFindingUtils  ClusterSelectingUtils_;
       LeastSquaresFitter LSQFitter_;
       SimpleCircleFitter CircleFitter_;

       // geometry helper tool
       gm2geom::StrawTrackerGeometry geom_;

       // message facility name
       std::string name_;

       // leaving full coord system for LSQ Fitter
       gm2geom::CoordSystemsStoreData cs_;

       // coordinate system helper tool
       std::map< std::string, gm2geom::CoordSystemsStoreData > csMap_;
 
       // fhicl parameters
       fhicl::ParameterSet fhicl_;
       fhicl::ParameterSet longTrackFinderParams_;
       fhicl::ParameterSet fastFitterParams_; 

       // variables
       std::string direction_;
       int  nmodules_;
       bool multipleTrackCandidates_;
       bool removeClusters_;

       // data helper organization
       struct DataObject {
          int    meanWire;
          double meanTime;
          art::Ptr<StrawSeedArtRecord> seed;
       };

       // DebugPlotting root plotting members
       std::unique_ptr<TFile> outputRootFile_;
       std::unique_ptr<RootManager> rm_;

       //DebugPlotting TMarkers to hold wire positions
       vector<TEllipse*> wirePositions_;
       int counter_;

       int maxWireSeparationForNeighbours_;
       bool makePlots_;

       string trackFinderName_;
       int minSeedsPerFit_;
       int maxSeedsToBeMultiTrackCand_;
       bool shareSeeds_;
       int maxOverlapSeeds_;
       int minPlanesPerFit_;
       bool applyClusterSelection_;
       double candidateSeedsMaxTime_;
       double candidateSeedsStrictTime_;
       int candidateSeedsMaxWire_;
       int maxWireSeparation_;
       double candidateTracksMaxTime_;
       double maxTimeSeparation_;
       double avgSeedRMSTime_;
       bool requiredFitConverges_;

       string fitterName_;

       std::vector<gm2geom::CoordSystem3Vector> wireCentreInWorld_;    
       std::vector<gm2geom::CoordSystem3Vector> wireCentreInModule_;    
       std::vector<gm2geom::CoordSystem3Vector> wireCentreInStation_;    
   };

}

#endif
