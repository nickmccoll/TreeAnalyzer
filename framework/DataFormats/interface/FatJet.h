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
class SubJet : public BaseRecoJet
{
public :
    SubJet() {}
    template <class InputCoordSystem>
    SubJet(const ROOT::Math::LorentzVector<InputCoordSystem> &mom,
            const int idx,
            const float csv, const  ASTypes::int8 hadronFlv = 0,
            const  ASTypes::int8 partonFlv=0) : BaseRecoJet(mom, idx,csv,hadronFlv,partonFlv) {}
    ~SubJet() {}

    template <class InputCoordSystem>
    void setRawMomentum(const ROOT::Math::LorentzVector<InputCoordSystem> &mom) {_rawMom=mom;}
    const MomentumF& rawMom() const {return _rawMom;}

protected:
    MomentumF      _rawMom       ;
};
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
            const float csv, const  ASTypes::size8 jetID, const  ASTypes::int8 hadronFlv = 0,
            const  ASTypes::int8 partonFlv=0, GenJet *gj = 0)
            : Jet(mom, idx,csv,jetID,hadronFlv,partonFlv,gj){}
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
    float     minSJCSV()  const;
    float     maxSJCSV()  const;
    MomentumF sdMom()     const;
    MomentumF rawSdMom()  const;


    const SubJet&     subJet(const size idx)  const;
    SubJet&           subJet(const size idx);

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
