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
            const float toRaw=0, const float csv=0, const float deep_csv=0)
            : Jet(mom, idx,toRaw,csv,deep_csv){}
    virtual ~FatJet() {}

    void addExtraInfo(const float jetID,
            const float bbt, const float tau1, const float tau2);

    void addSubJet(const SubJet& sj);

    float     bbt()       const;
    float     tau1()      const;
    float     tau2()      const;
    float     tau3()      const;
    float     tau2otau1() const;
    float     tau3otau1() const;
    float     tau3otau2() const;
    float     ecfb1()     const;
    float     ecfb2()     const;

    ASTypes::size      nSubJets()  const;
    MomentumF sdMom()     const;
    MomentumF rawSdMom()  const;


    const SubJet&     subJet(const ASTypes::size idx)  const;
    SubJet&           subJet(const ASTypes::size idx);
    const std::vector<SubJet>& subJets() const;

protected :

    float _bbt  =0;
    float _tau1 =0;
    float _tau2 =0;
    std::vector<SubJet> _sjs;


};
typedef std::vector<FatJet>   FatJetCollection;
}
#endif
