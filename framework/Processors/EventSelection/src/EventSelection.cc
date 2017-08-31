
#include "Processors/EventSelection/interface/EventSelection.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FillerConstants.h"
#include "DataFormats/interface/Lepton.h"

namespace TAna {
namespace EventSelection {


bool passEventFilters(const EventReader& reader_event) {
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_goodVertices) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_globalTightHalo2016Filter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_HBHENoiseFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_HBHENoiseIsoFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_EcalDeadCellTriggerPrimitiveFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_eeBadScFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::AnaTM_badMuons) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::AnaTM_badChargedHadrons) ) return false;
    if(reader_event.goodVtx == 0) return false;
    return true;
}




}

}



