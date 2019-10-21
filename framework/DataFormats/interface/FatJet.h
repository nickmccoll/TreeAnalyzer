#ifndef TREEANALYZER_DATAFORMATS_FATJET_H
#define TREEANALYZER_DATAFORMATS_FATJET_H

#include "DataFormats/interface/Jet.h"

namespace TAna {


//--------------------------------------------------------------------------------------------------
// GenFatJet...no difference wrt GenJet
//--------------------------------------------------------------------------------------------------
typedef GenJet GenFatJet;
typedef std::vector<GenFatJet> GenFatJetCollection;

//--------------------------------------------------------------------------------------------------
// SubJet...we add raw momentum on top of BaseRecoJet
//--------------------------------------------------------------------------------------------------
typedef BaseRecoJet SubJet;
typedef std::vector<SubJet> SubJetCollection;

//--------------------------------------------------------------------------------------------------
// FatJet: standard RecoFatJet
//--------------------------------------------------------------------------------------------------
class FatJet : public Jet
{
public :
    FatJet() {}

    template <class InputCoordSystem>
    FatJet(const ROOT::Math::LorentzVector<InputCoordSystem> &mom,
            const int idx,
            const float toRaw,
            const ASTypes::size8 jetID)
            : Jet(mom, idx,toRaw,jetID){}
    virtual ~FatJet() {}

    void addFJBtagging( const float bbt, const float deep_MDZHbb, const float deep_MDHbb
            , const float deep_Hbb);
    void addWTaging(const float deep_W);
    void addSubStructure(const float tau1, const float tau2, const float sdUp, const float sdDown);

    void addSubJet(const SubJet& sj);

    float     bbt()       const;
    float     deep_MDZHbb()       const;
    float     deep_MDHbb ()       const;
    float     deep_Hbb   ()       const;
    float     deep_W     ()       const;
    float     tau1()      const;
    float     tau2()      const;
    float     tau2otau1() const;
    float     rawSDMassUp()   const;
    float     rawSDMassDown() const;

    ASTypes::size      nSubJets()  const;
    MomentumF sdMom()     const;
    MomentumF rawSdMom()  const;


    const SubJet&     subJet(const ASTypes::size idx)  const;
    SubJet&           subJet(const ASTypes::size idx);
    const std::vector<SubJet>& subJets() const;

protected :
    float _tau1 =0;
    float _tau2 =0;
    float _sdUp =0;
    float _sdDown =0;
    float _bbt          =0;
    float _deep_MDZHbb  =0;
    float _deep_MDHbb   =0;
    float _deep_Hbb     =0;
    float _deep_W       =0;

    std::vector<SubJet> _sjs;


};
typedef std::vector<FatJet>   FatJetCollection;
}
#endif
