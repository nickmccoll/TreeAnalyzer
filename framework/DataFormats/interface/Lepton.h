#ifndef TREEANALYZER_DATAFORMATS_LEPTON_H
#define TREEANALYZER_DATAFORMATS_LEPTON_H

#include "DataFormats/interface/Momentum.h"
namespace TAna {


//--------------------------------------------------------------------------------------------------
// Lepton: Base class for leptons
//--------------------------------------------------------------------------------------------------
class Lepton : public IndexedMomentumF
{
public :
    Lepton() {}

    template <class InputCoordSystem>
    Lepton(const ROOT::Math::LorentzVector<InputCoordSystem> &mom,
            const int idx,
            const ASTypes::int8 q, const  float d0,
            const  float dz,const  float sip3D,const  float miniIso)
        : IndexedMomentumF(mom, idx), _q(q), _d0(d0),_dz(dz),_sip3D(sip3D),_miniIso(miniIso) {}
    ~Lepton() {}

    int   q       () const {return _q      ;}
    float d0      () const {return _d0     ;}
    float dz      () const {return _dz     ;}
    float sip3D   () const {return _sip3D  ;}
    float miniIso () const {return _miniIso;}


protected :
    ASTypes::int8  _q         = 0;
    float          _d0        = 0;
    float          _dz        = 0;
    float          _sip3D     = 0;
    float          _miniIso   = 0;

};
}
#endif
