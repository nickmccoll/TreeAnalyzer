
#include "DataFormats/interface/Electron.h"
#include "TreeReaders/interface/FillerConstants.h"

namespace TAna {
//--------------------------------------------------------------------------------------------------
void Electron::addElectronInfo(float scEta, float mvaID, float eaRelISO, ASTypes::size16 id) {
    _scEta   =  scEta   ;
    _mvaID   =  mvaID   ;
    _eaRelISO=  eaRelISO;
    _id      =  id      ;
}
//--------------------------------------------------------------------------------------------------
float Electron::scEta      () const {return _scEta;}
float Electron::mvaID      () const {return _mvaID;}
float Electron::eaRelISO   () const {return _eaRelISO;}
bool  Electron::passVetoID () const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_VETO);}
bool  Electron::passLooseID() const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_LOOSE);}
bool  Electron::passMedID  () const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_MED);}
bool  Electron::passTightID() const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_TIGHT);}
bool  Electron::passHEEPID () const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_HEEP);}
//--------------------------------------------------------------------------------------------------
bool  Electron::passVetoID_noISO () const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_NOISO_VETO);}
bool  Electron::passLooseID_noISO() const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_NOISO_LOOSE);}
bool  Electron::passMedID_noISO  () const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_NOISO_MED);}
bool  Electron::passTightID_noISO() const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_NOISO_TIGHT);}
bool  Electron::passHEEPID_noISO () const {return FillerConstants::doesPass(_id,FillerConstants::ELID_CUT_NOISO_HEEP);}
}
