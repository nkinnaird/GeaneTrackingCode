//!----------------------------------
//!  This class provides the methods for geant error propagation and correcting measured parameters
//!----------------------------------

#ifndef GEANEPARAMUTILS_HH
#define GEANEPARAMUTILS_HH

/////////////////////////////////////////////////////////////////////////////////////

// art includes

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "artg4/services/DetectorHolder_service.hh"

/////////////////////////////////////////////////////////////////////////////////////

// gm2Geane includes

#include "artg4/gm2Geane/gm2GeanePropagator.hh"
#include "Geant4/G4ErrorPropagatorData.hh"
#include "artg4/gm2Geane/gm2GeanePropagatorManager.hh"
#include "Geant4/G4ErrorPlaneSurfaceTarget.hh"
#include "artg4/gm2Geane/gm2GeaneFreeTrajState.hh"
#include "artg4/gm2Geane/gm2GeaneSurfaceTrajState.hh"
#include "artg4/gm2Geane/gm2GeaneMatrix.hh"
#include "artg4/gm2Geane/gm2GeanePlaneSurfaceTarget.hh"

/////////////////////////////////////////////////////////////////////////////////////

// Geant4 and field includes

#include "Geant4/G4UImanager.hh"

#include "Geant4/G4FieldManager.hh"
#include "Geant4/G4TransportationManager.hh"

#include "gm2geom/fields/gm2FieldManager_service.hh"

/////////////////////////////////////////////////////////////////////////////////////

#include "gm2geom/strawtracker/StrawTrackerGeometry.hh"

#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"
#include "gm2dataproducts/strawtracker/GeaneState.hh"

#include "gm2dataproducts/strawtracker/StrawDigitArtRecord.hh"
#include "gm2dataproducts/strawtracker/StrawDCADigitArtRecord.hh"

#include "gm2tracker/utils/GeaneTrackUtils.hh"
#include "gm2tracker/utils/GeaneDummyUtils.hh"

#include "gm2tracker/utils/StrawObjectSorter.hh"

#include "gm2util/common/dataModuleDefs.hh"

#include <Eigen/Dense>

namespace gm2strawtracker {

  class GeaneParamUtils {
     
    public :

     GeaneParamUtils(fhicl::ParameterSet const & p);

     void setZPositions(std::vector<double> positions) { geaneTrackerZPositions_ = positions; };
     void setCoordMap(std::map< std::string, gm2geom::CoordSystemsStoreData >& coordMap) { detCoordMap_ = coordMap; };

     int setupParams(art::Event & e, gm2strawtracker::TrackDetailArtRecord & trackFitDetails, bool firstTrackInEvent);

     int errorProp(gm2strawtracker::TrackDetailArtRecord & trackDetails);
     void calcMeasuredParams(gm2strawtracker::TrackDetailArtRecord & trackDetails);

     gm2strawtracker::GeaneDummyUtils dummyUtils_;
     std::vector<double> geaneTrackerZPositions_;

    private : 

     std::string name_;
     string fitMode_;

     std::string DummyModuleLabel_;
     std::string DummyInstanceName_;

     // GEANE related objects
     G4ErrorTarget* theTarget;
     G4ErrorMode theG4ErrorMode;
     G4ErrorPropagatorData* g4edata; 

     gm2geom::StrawTrackerGeometry sgeom_;

     double noHit; // Unphysical double value to signify the lack of a hit within a tracker plane, now only used to initialize values
     double wireDiamError; // uniform distribution standard deviation - not gaussian but approximated as such when fitting hits

     int G4EVERBOSE_; // Geant4 error_propagation verbose statements - the majority of these only work with the proper debug build of the geant4 error_propagation package.
     int trackingVerbose_; // Geant4 tracking level verbose statements.

     bool onlyPrimaryPositrons_;

     gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;

     std::map< std::string, gm2geom::CoordSystemsStoreData > detCoordMap_;

     bool useCircleGuess_;

     std::string stationStr;
     int stationNumber;

     G4ThreeVector tracingVX; // Set tracing vectors to be orthogonal in wire planes, Y vertical and 'X' radially out, and transform to UV later. Necessary because of what's in the GEANT4 source code.
     G4ThreeVector tracingWY;
     G4Normal3D surfNorm;

     std::vector<unsigned int> skipLayers_; // Layers that we'll ignore in the fit

  };

}

#endif
