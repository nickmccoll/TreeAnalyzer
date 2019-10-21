#ifndef TREEANALYZER_DATAFORMATS_JET_H
#define TREEANALYZER_DATAFORMATS_JET_H

#include "DataFormats/interface/Momentum.h"

namespace TAna {


//--------------------------------------------------------------------------------------------------
// GenJet...no difference wrt IndexedMomentum
//--------------------------------------------------------------------------------------------------
typedef IndexedMomentumF GenJet;
typedef std::vector<GenJet> GenJetCollection;

//--------------------------------------------------------------------------------------------------
// BaseRecoJet: Base class for jets, subjets, and fatjets
//--------------------------------------------------------------------------------------------------
class BaseRecoJet : public IndexedMomentumF
{
public :
    BaseRecoJet() {}

    template <class InputCoordSystem>
    BaseRecoJet(const ROOT::Math::LorentzVector<InputCoordSystem> &mom,
            const int idx,const float toRaw)
        : IndexedMomentumF(mom, idx),_toRaw(toRaw) {}
    virtual ~BaseRecoJet() {}

    void addMCInfo(const  ASTypes::int8 hadronFlv,
            const  ASTypes::int8 partonFlv,const float JECUnc);
    void addBTagging( const float deep_csv);

    void setRawFactor(float nrw) {_toRaw=nrw;}
    float toRawFactor() const {return _toRaw;}
    float deep_csv() const {return _deep_csv;}
    int  hadronFlv() const {return _hadronFlv;}
    int  partonFlv() const {return _partonFlv;}
    float  jecUnc() const {return _JECUnc;}

   const MomentumF rawMom() const;


protected :
    float          _toRaw     = 0;
    float          _deep_csv  = 0;
    ASTypes::int8  _hadronFlv = 0;
    ASTypes::int8  _partonFlv = 0;
    float          _JECUnc    = 0;

};
//--------------------------------------------------------------------------------------------------
// Jet: standard RecoJet
//--------------------------------------------------------------------------------------------------
class Jet : public BaseRecoJet
{
public :
    Jet() {}

    template <class InputCoordSystem>
    Jet(const ROOT::Math::LorentzVector<InputCoordSystem> &mom,
            const int idx,
            const float toRaw,
            const ASTypes::size8 jetID)
        : BaseRecoJet(mom, idx,toRaw), _jetID(jetID){}
    virtual ~Jet() {}

    void addMCInfo(const  ASTypes::int8 hadronFlv,
            const  ASTypes::int8 partonFlv,const float JECUnc, GenJet *gj);
    void addBTagging( const float deep_csv, const float csv, const float deep_flavor);

    bool passPUID() const;
    bool passTightID() const;
    bool passTightNoLepID()  const;
    float csv() const {return _csv;}
    float deep_flavor() const {return _deep_flavor;}

    ASTypes::size8 jetID() const {return _jetID;}

    const GenJet  *genJet()        const { return _gj;  }
    GenJet        *genJet()        { return _gj;  }



protected :
    ASTypes::size8 _jetID     = 0;
    float          _csv       = 0;
    float          _deep_flavor       = 0;
    GenJet  *      _gj        = 0;

};

typedef std::vector<Jet>   JetCollection;
}
#endif
