//!-----------------------------------------------------------------------------
//!  This class provides functions for the track extrapolation 
//!-----------------------------------------------------------------------------

#ifndef GEANEEXTRAPOLATIONUTILS_HH
#define GEANEEXTRAPOLATIONUTILS_HH

// art includes
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"

// constants includes
#include "gm2geom/common/Gm2Constants_service.hh"

// Eigen
#include <Eigen/Dense>

// utils
//#include "gm2util/common/RootManager.hh"
#include "gm2geom/common/StraightLineTools.hh"

// root includes
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TMath.h"

// c++ includes
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>

// field stuff
#include "gm2geom/fields/gm2FieldManager_service.hh"
#include "gm2geom/calo/CalorimeterGeometry.hh" // need the calo geometry for forwards extrapolation
#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"

// geometry
#include "gm2geom/coordSystems/CoordSystem.hh"
#include "gm2geom/coordSystems/CoordSystemsStoreData.hh"
#include "gm2geom/coordSystems/CoordSystemsStore_service.hh"
#include "gm2geom/coordSystems/CoordSystem3Vector.hh"  

// data product includes
#include "gm2dataproducts/mc/strawtracker/StrawArtRecord.hh"
#include "gm2dataproducts/mc/actions/track/TrackingActionArtRecord.hh"
#include "gm2dataproducts/strawtracker/DecayVertexArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"
#include "gm2dataproducts/strawtracker/StrawTimeIslandArtRecord.hh"
#include "gm2dataproducts/strawtracker/StrawDigitArtRecord.hh"

// geant4 stuff for volume checking
#include "Geant4/G4Navigator.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4TransportationManager.hh"
#include "artg4/services/DetectorHolder_service.hh"
#include "Geant4/G4ErrorPropagatorData.hh"
#include "Geant4/G4ErrorSurfaceTarget.hh"
#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4SteppingManager.hh"
#include "Geant4/G4EventManager.hh"
#include "Geant4/G4TrackingManager.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4StateManager.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4PhysicalVolumeStore.hh"
#include "Geant4/G4UnitsTable.hh"
#include "Geant4/G4RegionStore.hh"
#include "Geant4/G4UserLimits.hh"
#include "Geant4/G4Region.hh"
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4UIterminal.hh"

//gm2 Geane includes
#include "artg4/gm2Geane/gm2GeanePropagator.hh"
#include "artg4/gm2Geane/gm2GeanePropagatorManager.hh"
#include "artg4/gm2Geane/gm2GeaneFreeTrajState.hh"
#include "artg4/gm2Geane/gm2GeaneSurfaceTrajState.hh"
#include "artg4/gm2Geane/gm2GeaneCylSurfaceTarget.hh"
#include "artg4/gm2Geane/gm2GeaneTrackLengthTarget.hh"
#include "gm2tracker/utils/GeaneTrackUtils.hh"
#include "GeaneParamUtils.hh"
#include "artg4/gm2Geane/gm2GeanePropagator.hh"
#include "artg4/gm2Geane/gm2GeaneFreeTrajState.hh"
#include "artg4/gm2Geane/gm2GeaneSurfaceTrajState.hh"
#include "artg4/gm2Geane/gm2GeaneGeomVolumeTarget.hh"


// local
#include "gm2tracker/utils/StrawObjectSorter.hh"
#include "gm2tracker/utils/StrawGeomUtils.hh"
#include "gm2tracker/fitter/SimpleCircleFitter.hh"

// c++ stuff
#include "boost/format.hpp"
#include "TRandom.h"
#include "TMatrix.h"

namespace gm2strawtracker {

  struct GeaneFhiclOptions { 
    
  public : 
    
    // default constructor
    GeaneFhiclOptions() 
      : totalDistance()
      , nSteps()
      , stepsBeyondMin()
      , doBackwardsExtrap()
      , doForwardsExtrap()
      , extrapolater()
      , keepSteps()
      , onlyKeepStepsInVolume()
      , checkVolumes()
      , keepHitVolumes()
      , returnAtHitVolume()
      , returnAtTangentPoint()
      , useLargeStepsInStorageField()
      , useSmallStepsAtLowRadius()
      , useSmallStepsNearTangentPoint()
      , targetZpos()
      , stepSize()
      , smallStepSize()
      , largeStepSize()
      , initialRadius()
      , distanceBeyondTangentPoint()
      , useTruth()
      , useOnlyPDPs()
      , minNumPlanes()
      , cutPoorPValues()
      , pValueCut()
      , pmin()
      , pmax()
    {} 

    double totalDistance;
    unsigned int nSteps;
    int stepsBeyondMin;
    bool doBackwardsExtrap;
    bool doForwardsExtrap;
    std::string extrapolater;
    bool keepSteps;
    bool onlyKeepStepsInVolume;
    bool checkVolumes;
    bool keepHitVolumes;
    bool returnAtHitVolume;
    bool returnAtTangentPoint;
    bool useLargeStepsInStorageField;
    bool useSmallStepsAtLowRadius;
    bool useSmallStepsNearTangentPoint;
    double targetZpos;
    double stepSize;
    double smallStepSize;
    double largeStepSize;
    double initialRadius;
    double distanceBeyondTangentPoint;
    bool useTruth;
    bool useOnlyPDPs;
    int minNumPlanes;
    bool cutPoorPValues;
    double pValueCut;
    double pmin;
    double pmax;

    // copy constructor

#ifndef __ROOTCLING__
    
    GeaneFhiclOptions(fhicl::ParameterSet const & p) 
      : totalDistance(p.get<double>("extrapolationDistance"))
      , nSteps(p.get<int>("numberOfSteps"))
      , stepsBeyondMin(p.get<int>("stepsBeyondMin"))
      , doBackwardsExtrap(p.get<bool>("doBackwardsExtrap",true))
      , doForwardsExtrap(p.get<bool>("doForwardsExtrap",true))
      , extrapolater(p.get<std::string>("extrapolater","GeaneExtrapolater"))
      , keepSteps(p.get<bool>("keepSteps"))
      , onlyKeepStepsInVolume(p.get<bool>("onlyKeepStepsInVolume"))
      , checkVolumes(p.get<bool>("checkVolumes"))
      , keepHitVolumes(p.get<bool>("keepHitVolumes"))
      , returnAtHitVolume(p.get<bool>("returnAtHitVolume"))
      , returnAtTangentPoint(p.get<bool>("returnAtTangentPoint"))
      , useLargeStepsInStorageField(p.get<bool>("useLargeStepsInStorageField"))
      , useSmallStepsAtLowRadius(p.get<bool>("useSmallStepsAtLowRadius"))
      , useSmallStepsNearTangentPoint(p.get<bool>("useSmallStepsNearTangentPoint"))
      , targetZpos(p.get<double>("targetZpos"))
      , stepSize(p.get<double>("stepSize"))
      , smallStepSize(p.get<double>("smallStepSize"))
      , largeStepSize(p.get<double>("largeStepSize"))
      , initialRadius(p.get<double>("initialRadius"))
      , distanceBeyondTangentPoint(p.get<double>("distanceBeyondTangentPoint"))
      , useTruth(p.get<bool>("useTruth"))
      , useOnlyPDPs(p.get<bool>("useOnlyPDPs"))
      , minNumPlanes(p.get<int>("minNumPlanes"))
      , cutPoorPValues(p.get<bool>("cutPoorPValues",false))
      , pValueCut(p.get<double>("pValueCut",0.005))
      , pmin(p.get<double>("pmin"))
      , pmax(p.get<double>("pmax"))
    {}

#endif
    

  };

   class GeaneExtrapolationUtils {

      public :

       // default constructor
       GeaneExtrapolationUtils();

       // decay vertex reconstruction function
     bool reconstructDecayVertex( const TrackArtRecord& track, DecayVertexArtRecord& decayVertex, const GeaneFhiclOptions& fclOptions, Eigen::MatrixXd initialError, DecayVertexArtRecord& vertexSteps); 
     bool reconstructDecayVertex( const gm2strawtracker::ExtrapolationStep& startPoint, DecayVertexArtRecord& decayVertex, int stationNum, const GeaneFhiclOptions& fclOptions, Eigen::MatrixXd initialError, DecayVertexArtRecord& vertexSteps); 
     
       // forwards extrapolation to calo function - these call above functions but with more useful names
       bool extrapolateToCalorimeter( const TrackArtRecord& track, DecayVertexArtRecord& decayVertex, const GeaneFhiclOptions& fclOptions ); 
       bool extrapolateToCalorimeter( const gm2strawtracker::ExtrapolationStep& startPoint, DecayVertexArtRecord& decayVertex, int stationNum, const GeaneFhiclOptions& fclOptions ); 
     
       // for extrapolating to the true decay position
       bool extrapolateToTrueAzimuth( const ExtrapolationStep &startPoint, DecayVertexArtRecord& decayVertex, const GeaneFhiclOptions& fclOptions);

       // for extrapolating to a plane intersecting the given 'target position'
       bool extrapolateToPlane( const gm2strawtracker::ExtrapolationStep& startPoint, const double targetZ, std::string coordSysName, gm2strawtracker::ExtrapolationStep& endPoint, std::string direction = "backwards", double stepSize = 0.01);

     void extrapolateForwardsAndBackwards(gm2strawtracker::ExtrapolationStep& startingPoint, gm2strawtracker::ExtrapolationStep& finalPoint, std::vector<ExtrapolationStep>& forwardsSteps, std::vector<ExtrapolationStep>& backwardsSteps, double stepSize, bool useOldCode = false, bool forceUniformField = false, double distance = -1.0);

       void extrapolateBackwardsAndForwards(gm2strawtracker::ExtrapolationStep& startingPoint, gm2strawtracker::ExtrapolationStep& finalPoint, double stepSize, double distance);
    
       // find out if the step hit a solid volume
       bool didStepHitVolume(std::string givenVolume);


       void setCoordMap(std::map< std::string, gm2geom::CoordSystemsStoreData >& coordMap) { detCoordMap_ = coordMap; };

 
       // set the coordinate system data
       void setCoordSysData( gm2geom::CoordSystemsStoreData const & c ) {
         cs_ = c;
       } 
     
       // set fhicl parameters
       void setFhiclCuts( fhicl::ParameterSet const & p ) {
         fhicl_ = p;
       }
     
       void setDetectorHolder( art::ServiceHandle<artg4::DetectorHolderService> const & dh ) {
         detectorHolder_ = dh;
       }
     
       void setupNavigator( G4Navigator* nav) {
         nav_ = nav;
       }
 
      private :

       // for updating the time of each step
       double updateTime(gm2strawtracker::ExtrapolationStep step, double stepSize, bool backwardsFlag);
     
     gm2strawtracker::ExtrapolationStep GeaneExtrapolation( gm2strawtracker::ExtrapolationStep step, double stepSize, std::string extrapolationDirection, double pos_err[3], bool leftFieldMap, bool forceUniformField = false);
     
     gm2strawtracker::ExtrapolationStep oldGeaneExtrapolation( gm2strawtracker::ExtrapolationStep step, double stepSize, std::string extrapolationDirection, bool leftFieldMap, bool forceUniformField = false);

       // Jacobian propagation
       void propagateJacobian(double trackParameters[7], double A[4], double B[4], double C[4], double H, Eigen::MatrixXd covMatrix, Eigen::MatrixXd J);

       // message facility name
       std::string name_;
       
       // coordinate system helper tool
       gm2geom::CoordSystemsStoreData cs_;

       // fhicl parameters
       fhicl::ParameterSet fhicl_;

       // circle fitter
       SimpleCircleFitter circleFitter_;

       // services tools
       art::ServiceHandle<gm2geom::Gm2Constants> gm2consts_;
       art::ServiceHandle<gm2geom::gm2FieldManager> fieldManager_;    
       art::ServiceHandle<artg4::DetectorHolderService> detectorHolder_;
       G4Navigator* nav_;

       // utils 
       gm2strawtracker::StrawGeomUtils geomUtils_;
       gm2geom::StrawTrackerGeometry sgeom_;

       // Geane Related objects
       G4ErrorTarget* theTarget;
       G4ErrorTarget* theTarget2;
       G4ErrorTarget* theTarget3;
       G4ErrorTarget* theTarget4;
       G4ErrorMode theG4ErrorMode;
       G4ErrorPropagatorData* g4edata;
       gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;
       G4VPhysicalVolume* world;
       std::map< std::string, gm2geom::CoordSystemsStoreData > detCoordMap_;
       int stationNumber;

       G4ThreeVector tracingVX; // Set tracing vectors to be orthogonal in wire planes, Y vertical and 'X' radially out, and transform to UV later. Necessary because of what's in the GEANT4 source code.
       G4ThreeVector tracingWY;
       G4Normal3D surfNorm;




       G4ThreeVector initialRadius;
       G4ThreeVector initialMomentum;
       long double eCharge;
       G4double p;
       G4double B;
       G4double rc;
       G4ThreeVector rs;
       G4ThreeVector finalRadius;
       G4ThreeVector displacement;
       G4double trackLength;
       G4double rxp;
       G4double radialMomentum;
       G4double scaledRadialMomentum;
       G4double theta;
       G4double verticalMomentum;
       G4double endRadius;
       G4double px;
       G4double pz;
       G4double verticalAngle;
       G4double stepLength;
       G4double verticalPosition;
       G4double measuredTrackLength;
       G4ThreeVector beforeStepPosition;
       G4ThreeVector afterStepPosition;
       G4ThreeVector deltaStepPosition;
   };
  
}

#endif
