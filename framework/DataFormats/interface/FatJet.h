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
            const float toRaw, const float csv, const  ASTypes::size8 jetID, const  ASTypes::int8 hadronFlv = 0,
            const  ASTypes::int8 partonFlv=0,const float JECUnc =0, GenJet *gj = 0)
            : Jet(mom, idx,toRaw,csv,jetID,hadronFlv,partonFlv,JECUnc,gj){}
    ~FatJet() {}

    void addFatJetInfo(const float bbt, const float tau1, const float tau2, const float tau3);
    void addSubJet(const SubJet& sj);

    float     bbt()       const;
    float     tau1()      const;
    float     tau2()      const;
    float     tau3()      const;
    float     tau2otau1() const;
    float     tau3otau1() const;
    float     tau3otau2() const;

    size      nSubJets()  const;
    MomentumF sdMom()     const;
    MomentumF rawSdMom()  const;


    const SubJet&     subJet(const size idx)  const;
    SubJet&           subJet(const size idx);
    const std::vector<SubJet>& subJets() const;

protected :

    float _bbt  =0;
    float _tau1 =0;
    float _tau2 =0;
    float _tau3 =0;
    std::vector<SubJet> _sjs;


};
typedef std::vector<FatJet>   FatJetCollection;
}
#endif
