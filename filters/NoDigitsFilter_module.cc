// Kill events that don't contain digits
// This is an easy way to remove laser events

// art includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// data product includes
#include "gm2dataproducts/strawtracker/StrawDigitArtRecord.hh"

// common util
#include "gm2util/common/dataModuleDefs.hh"

namespace gm2strawtracker {
  class NoDigitsFilter;
}

class gm2strawtracker::NoDigitsFilter : public art::EDFilter {

   public:
     explicit NoDigitsFilter(fhicl::ParameterSet const & p);

     bool filter(art::Event & e) override;

   private:

    std::string digitModuleLabel_;
    std::string digitInstanceLabel_;

};


gm2strawtracker::NoDigitsFilter::NoDigitsFilter(fhicl::ParameterSet const & p) 
  : digitModuleLabel_( p.get<std::string>("digitModuleLabel",dataModuleDefs::deadTimeModuleLabel()) )
  , digitInstanceLabel_( p.get<std::string>("digitInstanceLabel",dataModuleDefs::digitInstanceLabel()) )
{}

bool gm2strawtracker::NoDigitsFilter::filter(art::Event & e) {

  //Get digits hits
  art::Handle<gm2strawtracker::StrawDigitArtRecordCollection> digitDataHandle;
  bool foundStrawHits = e.getByLabel(digitModuleLabel_,digitInstanceLabel_,digitDataHandle);
  if( ! foundStrawHits ) {
    return false;
  } else {
    return true;
  }

}

DEFINE_ART_MODULE(gm2strawtracker::NoDigitsFilter)
