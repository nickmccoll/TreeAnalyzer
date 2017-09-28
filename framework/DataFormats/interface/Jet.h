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
            const int idx,
            const float csv, const  ASTypes::int8 hadronFlv = 0,
            const  ASTypes::int8 partonFlv=0)
        : IndexedMomentumF(mom, idx), _csv(csv), _hadronFlv(hadronFlv)
          ,_partonFlv(partonFlv) {}
    ~BaseRecoJet() {}

    float csv() const {return _csv;}
    int  hadronFlv() const {return _hadronFlv;}
    int  partonFlv() const {return _partonFlv;}


protected :
    float          _csv       = 0;
    ASTypes::int8  _hadronFlv = 0;
    ASTypes::int8  _partonFlv = 0;

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
            const float csv, const  ASTypes::size8 jetID, const  ASTypes::int8 hadronFlv = 0,
            const  ASTypes::int8 partonFlv=0, GenJet *gj = 0)
        : BaseRecoJet(mom, idx,csv,hadronFlv,partonFlv),_jetID(jetID),_gj(gj) {}
    ~Jet() {}


    bool passPUID() const;
    bool passLooseID() const;
    bool passTightID() const;

    const GenJet  *genJet()        const { return _gj;  }
    GenJet        *genJet()        { return _gj;  }



protected :
    ASTypes::size8 _jetID     = 0;
    GenJet  *      _gj        = 0;

};

typedef std::vector<Jet>   JetCollection;
}
#endif
