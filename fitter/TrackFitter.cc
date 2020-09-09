//! include header
#include "TrackFitter.hh"

using namespace gm2strawtracker;

//!-----------------------------------------------------------
//! standard constructor 
//!-----------------------------------------------------------
template< typename StateVector >
TrackFitter<StateVector>::TrackFitter() 
   : states_()
   , iterations_(0)
   , fieldType_(0)
   , fhicl_()
   , wires_()
   , cs_()
   , geom_()
   , gm2consts_()
   , gm2fields_()	
 {}

// put explicit template instantiation at here, so the compilier knows all 
// about the various implementation details TODO figure out why this is needed

template class TrackFitter<gm2strawtracker::Helix>;
template class TrackFitter<gm2strawtracker::GeaneState>;
