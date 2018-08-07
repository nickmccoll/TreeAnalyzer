#ifndef TREEANALYZER_DATAFORMATS_ELECTRON_H
#define TREEANALYZER_DATAFORMATS_ELECTRON_H

#include "DataFormats/interface/Lepton.h"
namespace TAna {


//--------------------------------------------------------------------------------------------------
// Electron: Class for electrons
//--------------------------------------------------------------------------------------------------
class Electron : public Lepton
{
public :
    Electron() : Lepton(false) {}

    template <class InputCoordSystem>
    Electron(const ROOT::Math::LorentzVector<InputCoordSystem> &mom,
            const int idx,
            const ASTypes::int8 q, const  float d0,
            const  float dz,const  float sip3D,const  float miniIso, const float dRnorm, const float lepAct_o_pt)
        : Lepton(mom, idx,q,d0,dz,sip3D,miniIso,false, dRnorm, lepAct_o_pt) {}
    ~Electron() {}

    void addElectronInfo(float scEta, float mvaID, ASTypes::size8 mvaIDCat, float eaRelISO, ASTypes::size16 id, float sc_act_o_pt, float sc_dr_act);

    float scEta      () const;
    float mvaID      () const;
    size  mvaIDCat   () const;
    float eaRelISO   () const;
    float eaRelISO   () const;
    float sc_act_o_pt() const;
    float sc_dr_act  () const;

    bool  passVetoID () const;
    bool  passLooseID() const;
    bool  passMedID  () const;
    bool  passTightID() const;
    bool  passHEEPID () const;
    bool  passMVA80ID() const;
    bool  passMVA90ID() const;

    bool  passVetoID_noISO () const;
    bool  passLooseID_noISO() const;
    bool  passMedID_noISO  () const;
    bool  passTightID_noISO() const;
    bool  passHEEPID_noISO () const;

protected :
    float           _scEta      =0;
    float           _mvaID      =0;
    float           _eaRelISO   =0;
    ASTypes::size16 _id         =0;
    ASTypes::size8  _mvaIDCat   =0;
    float           _sc_act_o_pt=0;
    float           _sc_dr_act  =0;

};

typedef std::vector<Electron> ElectronCollection;
}
#endif
