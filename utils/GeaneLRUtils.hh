//!----------------------------------
//!  This class provides useful functions for LR fitting of Geane tracks
//!----------------------------------

#ifndef GEANELRUTILS_HH
#define GEANELRUTILS_HH

#include "gm2dataproducts/strawtracker/TrackDetailArtRecord.hh"
#include "gm2tracker/utils/GeaneTrackUtils.hh"
#include "gm2tracker/fitter/GeaneFitter.hh"

#include "gm2tracker/utils/StrawObjectSorter.hh"

#include "boost/dynamic_bitset.hpp"

namespace gm2strawtracker {

  class GeaneLRUtils {
     
    public :

     GeaneLRUtils(fhicl::ParameterSet const & p);

     // updates errors if dca is less than some value
     void lockSmallDCAErrors(gm2strawtracker::TrackDetailArtRecord & geaneTrack);

     // locks sides to the center of the wire if dca is less than some value
     std::vector<gm2strawtracker::GeaneHitSide> lockSmallDCACenters(gm2strawtracker::TrackDetailArtRecord & geaneTrack);

     // fill LR sides based on fit - used after wire fitting
     void fillLRFromFit(gm2strawtracker::TrackDetailArtRecord & geaneTrack);

     // fill LR sides based on geom and angle
     void fillLRFromGeom(gm2strawtracker::TrackDetailArtRecord & geaneTrack);

     /////////////////////////////////////////////////////////////////////////////////////
     // methods below here are for the fullSeqFit mode
     /////////////////////////////////////////////////////////////////////////////////////
 
     // set all UV hits to unkown to start before locking (also lock small dcas in this method)
     void setUnknownSides(gm2strawtracker::TrackDetailArtRecord & geaneTrack);

     // set UV sides from the node sides
     void fillLRFromNodes(gm2strawtracker::TrackDetailArtRecord & geaneTrack);

     //get U or V planes which have sides that are unknown (not locked in some way)
     std::vector<int> getUVUnkowns(gm2strawtracker::TrackDetailArtRecord & geaneTrack, bool Uplanes);

     //get the best U or V sequences using the fast/approximate matrix checker
     std::vector<std::pair<double, int> > getTopSequences(gm2strawtracker::TrackDetailArtRecord & geaneTrack, bool Uplanes, int unkownsSize);

     // adds or subtracts dca to wire position based on particular sequence
     void modifyMeasuredParams(gm2strawtracker::TrackDetailArtRecord & track); 

     // sets sides for full fit track based on U and V sequences
     void getSequenceSides(std::vector<int> UVunknownsInThisEvent, int sequenceInt, std::vector<gm2strawtracker::GeaneHitSide> *sides);


    private : 

     gm2strawtracker::GeaneTrackUtils geaneTrackUtils_;
     gm2strawtracker::GeaneFitter matrixCalculations_;

     double lockLowDCAs_;
     double wireCenterError; // uniform distribution standard deviation - for when the hits are close to the wire


  };

}

#endif
