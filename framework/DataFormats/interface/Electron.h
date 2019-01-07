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
            const  float dz,const  float sip3D)
        : Lepton(mom, idx,q,d0,dz,sip3D,false) {}
    ~Electron() {}

    void addElectronInfo(float scEta,float scE, float mvaID,
            float miniIsoFP, ASTypes::size16 id, float sc_act_o_pt, float sc_dr_act);

    float scEta      () const;
    float scE        () const;
    float mvaID      () const;
    float miniIsoFP  () const;
    float sc_act_o_pt() const;
    float sc_dr_act  () const;

    bool  passLooseID()    const;
    bool  passMedID  ()    const;
    bool  passTightID()    const;
    bool  passHEEPID ()    const;
    bool  passMVAHZZ()     const;
    bool  passMVALooseID() const;
    bool  passMVA80ID()    const;
    bool  passMVA90ID()    const;

    bool  passLooseID_noIso()    const;
    bool  passMedID_noIso  ()    const;
    bool  passTightID_noIso()    const;
    bool  passHEEPID_noIso ()    const;
    bool  passMVALooseID_noIso() const;
    bool  passMVA80ID_noIso()    const;
    bool  passMVA90ID_noIso()    const;

protected :
    float           _scEta      =0;
    float           _scE        =0;
    float           _mvaID      =0;
    float           _miniIsoFP  =0;
    ASTypes::size16 _id         =0;
    float           _sc_act_o_pt=0;
    float           _sc_dr_act  =0;

};

typedef std::vector<Electron> ElectronCollection;
}
#endif
