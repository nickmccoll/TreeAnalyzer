
#include "DataFormats/interface/Muon.h"
#include "TreeReaders/interface/FillerConstants.h"

namespace TAna {
//--------------------------------------------------------------------------------------------------
void Muon::addMuonInfo(float dbRelISO, ASTypes::size16 id) {
    _dBRelISO=  dbRelISO;
    _id      =  id      ;
}
//--------------------------------------------------------------------------------------------------


float Muon::dbRelISO   () const {return _dBRelISO;}
bool  Muon::passSoftID () const {return FillerConstants::doesPass(_id,FillerConstants::MUID_SOFT);}
bool  Muon::passLooseID() const {return FillerConstants::doesPass(_id,FillerConstants::MUID_LOOSE);}
bool  Muon::passMedID  () const {return FillerConstants::doesPass(_id,FillerConstants::MUID_MED);}
bool  Muon::passMed16ID() const {return FillerConstants::doesPass(_id,FillerConstants::MUID_MED16);}
bool  Muon::passTightID() const {return FillerConstants::doesPass(_id,FillerConstants::MUID_TIGHT);}
bool  Muon::passHighPT () const {return FillerConstants::doesPass(_id,FillerConstants::MUID_HIGHPT);}
}
