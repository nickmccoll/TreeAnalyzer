
#include "DataFormats/interface/FatJet.h"

namespace TAna {
//--------------------------------------------------------------------------------------------------
void FatJet::addExtraInfo(const float jetID,
        const float bbt, const float tau1, const float tau2){
    Jet::addExtraInfo(jetID);
    _bbt  = bbt  ;
    _tau1 = tau1 ;
    _tau2 = tau2 ;
    }
//--------------------------------------------------------------------------------------------------
void FatJet::addSubJet(const SubJet& sj) {_sjs.push_back(sj);}
//--------------------------------------------------------------------------------------------------
float     FatJet::bbt()       const{return _bbt;}
float     FatJet::tau1()      const{return _tau1;}
float     FatJet::tau2()      const{return _tau2;}
float     FatJet::tau2otau1() const{return _tau1 == 0 ? 99 : _tau2/_tau1;}
//--------------------------------------------------------------------------------------------------
size     FatJet::nSubJets()  const{return _sjs.size();}
MomentumF FatJet::sdMom()     const{
    MomentumF sd;
    for(const auto& sj : _sjs) { sd.p4() += sj.p4(); }
    return sd;
}
MomentumF FatJet::rawSdMom()  const{
    MomentumF sd;
    for(const auto& sj : _sjs) { sd.p4() += sj.rawMom().p4(); }
    return sd;
}
//--------------------------------------------------------------------------------------------------
const SubJet&     FatJet::subJet(const size idx)  const {
    if(idx >= nSubJets())
        throw std::out_of_range("FatJet::subJet() -> Not a valid SubJet idx!)");
    return _sjs[idx];
}
SubJet&           FatJet::subJet(const size idx) {
    if(idx >= nSubJets())
        throw std::out_of_range("FatJet::subJet() -> Not a valid SubJet idx!)");
    return _sjs[idx];
}
const std::vector<SubJet>& FatJet::subJets() const {return _sjs;}
}
