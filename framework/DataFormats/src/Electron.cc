
#include "DataFormats/interface/Electron.h"
#include "TreeReaders/interface/FillerConstants.h"

using namespace FillerConstants;
namespace TAna {
//--------------------------------------------------------------------------------------------------
void Electron::addElectronInfo(float scEta,float scE, float mvaID,
        float miniIsoFP, ASTypes::size16 id, float sc_act_o_pt, float sc_dr_act){
    _scEta   =  scEta   ;
    _scE   =  scE   ;
    _mvaID   =  mvaID   ;
    _miniIsoFP=  miniIsoFP   ;
    _id      =  id      ;
    _sc_act_o_pt = sc_act_o_pt;
    _sc_dr_act = sc_dr_act;
}
//--------------------------------------------------------------------------------------------------
float Electron::scEta      () const {return _scEta;}
float Electron::mvaID      () const {return _mvaID;}
float Electron::miniIsoFP   () const {return _miniIsoFP;}
float Electron::sc_act_o_pt() const {return _sc_act_o_pt;}
float Electron::sc_dr_act  () const {return _sc_dr_act;}
//--------------------------------------------------------------------------------------------------
bool  Electron::passLooseID()    const{return doesPass(_id,ELID_cut_loose );}
bool  Electron::passMedID  ()    const{return doesPass(_id,ELID_cut_medium);}
bool  Electron::passTightID()    const{return doesPass(_id,ELID_cut_tight );}
bool  Electron::passHEEPID ()    const{return doesPass(_id,ELID_heep);}
bool  Electron::passMVAHZZ()     const{return doesPass(_id,ELID_mva_wpHZZ);}
bool  Electron::passMVALooseID() const{return doesPass(_id,ELID_mva_wpLoose);}
bool  Electron::passMVA80ID()    const{return doesPass(_id,ELID_mva_wp80);}
bool  Electron::passMVA90ID()    const{return doesPass(_id,ELID_mva_wp90);}
//--------------------------------------------------------------------------------------------------
bool  Electron::passLooseID_noIso()    const{return doesPass(_id,ELID_cut_loose_noIso );}
bool  Electron::passMedID_noIso  ()    const{return doesPass(_id,ELID_cut_medium_noIso);}
bool  Electron::passTightID_noIso()    const{return doesPass(_id,ELID_cut_tight_noIso );}
bool  Electron::passHEEPID_noIso ()    const{return doesPass(_id,ELID_heep_noIso);}
bool  Electron::passMVALooseID_noIso() const{return doesPass(_id,ELID_mva_wpLoose_noIso);}
bool  Electron::passMVA80ID_noIso()    const{return doesPass(_id,ELID_mva_wp80_noIso);}
bool  Electron::passMVA90ID_noIso()    const{return doesPass(_id,ELID_mva_wp90_noIso);}
}
