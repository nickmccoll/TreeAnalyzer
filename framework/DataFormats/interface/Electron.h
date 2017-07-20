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
            const  float dz,const  float sip3D,const  float miniIso)
        : Lepton(mom, idx,q,d0,dz,sip3D,miniIso,false) {}
    ~Electron() {}

    void addElectronInfo(float scEta, float mvaID, float eaRelISO, ASTypes::size16 id);

    float scEta      () const;
    float mvaID      () const;
    float eaRelISO   () const;

    bool  passVetoID () const;
    bool  passLooseID() const;
    bool  passMedID  () const;
    bool  passTightID() const;
    bool  passHEEPID () const;

    bool  passVetoID_noISO () const;
    bool  passLooseID_noISO() const;
    bool  passMedID_noISO  () const;
    bool  passTightID_noISO() const;
    bool  passHEEPID_noISO () const;

protected :
    float           _scEta     =0;
    float           _mvaID     =0;
    float           _eaRelISO  =0;
    ASTypes::size16 _id        =0;

};

typedef std::vector<Electron> ElectronCollection;
}
#endif
