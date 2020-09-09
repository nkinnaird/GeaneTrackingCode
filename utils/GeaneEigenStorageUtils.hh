//!----------------------------------
//!  This class provides ability to store/read data from Eigen matrices and vectors for the GEANEArtRecord
//!----------------------------------

#ifndef GEANEEIGENSTORAGEUTILS_HH
#define GEANEEIGENSTORAGEUTILS_HH

#include "art/Framework/Principal/Event.h"
#include "gm2dataproducts/strawtracker/TrackArtRecord.hh"
#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"
#include "gm2dataproducts/strawtracker/DecayVertexArtRecord.hh"

namespace gm2strawtracker {
  
  class GeaneEigenStorageUtils {
    
    public :

    static void PrepareEigenForStorage(std::unique_ptr<DecayVertexArtRecordCollection>& recordsPtr);

    // Function to take eigen objects from art records and convert to flat array of doubles for storage
    static void PrepareEigenForStorage(std::unique_ptr<TrackDetailArtRecordCollection>& recordsPtr);
    static void PrepareEigenForStorage(std::unique_ptr<TrackArtRecordCollection>& recordsPtr);

    // Function to take read flat array of double from art records and populate Eigen objects
    static TrackArtRecordCollection ReadEigenFromStorage(const TrackArtRecordCollection& records);
    static TrackDetailArtRecordCollection ReadEigenFromStorage(const TrackDetailArtRecordCollection& records);
    static DecayVertexArtRecordCollection ReadEigenFromStorage(const DecayVertexArtRecordCollection& records);

    static TrackArtRecord ReadEigenFromStorage(const TrackArtRecord& record);
    static TrackDetailArtRecord ReadEigenFromStorage(const TrackDetailArtRecord& record);
    static DecayVertexArtRecord ReadEigenFromStorage(const DecayVertexArtRecord& record);


  };

}

#endif
