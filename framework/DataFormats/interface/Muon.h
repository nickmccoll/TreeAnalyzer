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
    Muon() {}

    template <class InputCoordSystem>
    Muon(const ROOT::Math::LorentzVector<InputCoordSystem> &mom,
            const int idx,
            const ASTypes::int8 q, const  float d0,
            const  float dz,const  float sip3D,const  float miniIso)
        : Lepton(mom, idx,q,d0,dz,sip3D,miniIso) {}
    ~Muon() {}

    void addMuonInfo(float dbRelISO, ASTypes::size16 id);

    float dbRelISO   () const;
    bool  passSoftID () const;
    bool  passLooseID() const;
    bool  passMedID  () const;
    bool  passMed16ID() const;
    bool  passTightID() const;
    bool  passHighPT () const;

protected :
    float           _dBRelISO  =0;
    ASTypes::size16 _id        =0;
};
typedef std::vector<Muon> MuonCollection;
}
#endif
