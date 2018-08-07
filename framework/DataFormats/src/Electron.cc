
#include "DataFormats/interface/Electron.h"
#include "TreeReaders/interface/FillerConstants.h"

namespace TAna {
//--------------------------------------------------------------------------------------------------
void Electron::addElectronInfo(float scEta, float mvaID, ASTypes::size8 mvaIDCat, float eaRelISO, ASTypes::size16 id, float sc_act_o_pt, float sc_dr_act) {
    _scEta   =  scEta   ;
    _mvaID   =  mvaID   ;
    _mvaIDCat=  mvaIDCat   ;
    _eaRelISO=  eaRelISO;
    _id      =  id      ;
    _sc_act_o_pt = sc_act_o_pt;
    _sc_dr_act = sc_dr_act;
}
//--------------------------------------------------------------------------------------------------
float Electron::scEta      () const {return _scEta;}
float Electron::mvaID      () const {return _mvaID;}
size  Electron::mvaIDCat   () const {return _mvaIDCat;}
float Electron::eaRelISO   () const {return _eaRelISO;}
float Electron::sc_act_o_pt() const {return _sc_act_o_pt;}
float Electron::sc_dr_act  () const {return _sc_dr_act;}
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
//dont like this....but I don't think it will get anywhere
bool  Electron::passMVA80ID      () const { return mvaID() >= (mvaIDCat() == 0 ?  0.941 : (mvaIDCat() == 1  ?  0.899 : 0.758));}
bool  Electron::passMVA90ID      () const { return mvaID() >= (mvaIDCat() == 0 ?  0.837 : (mvaIDCat() == 1  ?  0.715 : 0.357));}
}
