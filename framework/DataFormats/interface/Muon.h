#ifndef TREEANALYZER_DATAFORMATS_MUON_H
#define TREEANALYZER_DATAFORMATS_MUON_H

#include "DataFormats/interface/Lepton.h"
namespace TAna {
//--------------------------------------------------------------------------------------------------
// Muon: Class for muons
//--------------------------------------------------------------------------------------------------
class Muon : public Lepton
{
public :
    Muon() : Lepton(true) {}

    template <class InputCoordSystem>
    Muon(const ROOT::Math::LorentzVector<InputCoordSystem> &mom,
            const int idx,
            const ASTypes::int8 q, const  float d0,
            const  float dz,const  float sip3D)
        : Lepton(mom, idx,q,d0,dz,sip3D,true) {}
    ~Muon() {}

    void setMuonInfo(ASTypes::size id);

    ASTypes::size  id () const {return _id;}
    bool  passSoftID () const;
    bool  passLooseID() const;
    bool  passMedID  () const;
    bool  passTightID() const;
    bool  passHighPT () const;

protected :
    ASTypes::size _id        =0;
};
typedef std::vector<Muon> MuonCollection;
}
#endif
