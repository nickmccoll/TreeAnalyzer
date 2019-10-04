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
    Lepton(const bool muon) : _muon(muon) {}

    template <class InputCoordSystem>
    Lepton(const ROOT::Math::LorentzVector<InputCoordSystem> &mom,
            const int idx,
            const ASTypes::int8 q, const  float d0,
            const  float dz,const  float sip3D, const bool muon)
        : IndexedMomentumF(mom, idx), _q(q), _d0(d0),_dz(dz),_sip3D(sip3D), _muon(muon) {}

    ~Lepton() {}

    void setIsos(  float miniIso,  float pfIso, float ptRel, float ptRatio, float ttHMVA)
    { _miniIso=miniIso;_pfIso=pfIso;_ptRel=ptRel;_ptRatio=ptRatio;_ttHMVA=ttHMVA;};
    void setSysts (float dRnorm, const float lepAct_o_pt){_dRnorm =dRnorm;_lepAct_o_pt=lepAct_o_pt;}

    int   q           () const {return _q      ;}
    float d0          () const {return _d0     ;}
    float dz          () const {return _dz     ;}
    float sip3D       () const {return _sip3D  ;}

    float miniIso     () const {return _miniIso;}
    float pfIso       () const {return _pfIso;}
    float ptRel       () const {return _ptRel;}
    float ptRatio     () const {return _ptRatio;}
    float ttHMVA      () const {return _ttHMVA;}

    bool isMuon       () const {return _muon;}
    bool isElectron   () const {return !_muon;}

    float dRnorm      () const {return _dRnorm ;}
    float lepAct_o_pt () const {return _lepAct_o_pt;}

    bool  passInclID  () const {return true;}
    float inclIso     () const {return 0.0;}


protected :
    ASTypes::int8  _q           = 0;
    float          _d0          = 0;
    float          _dz          = 0;
    float          _sip3D       = 0;
    float          _miniIso     = 0;
    float          _pfIso       = 0;
    float          _ptRel       = 0;
    float          _ptRatio     = 0;
    float          _ttHMVA     = 0;

    bool           _muon        = false;

    float 		   _dRnorm      = 0;
    float          _lepAct_o_pt = 0;


};
}
#endif
