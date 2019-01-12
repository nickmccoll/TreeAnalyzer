
#include "DataFormats/interface/FatJet.h"

namespace TAna {
//--------------------------------------------------------------------------------------------------
void FatJet::addExtraInfo(const float jetID,
        const float bbt, const float tau1, const float tau2, const float tau3,
        const float ecfb1,const float ecfb2){
    Jet::addExtraInfo(jetID);
    _bbt  = bbt  ;
    _tau1 = tau1 ;
    _tau2 = tau2 ;
    _tau3 = tau3 ;
    _ecfb1= ecfb1;
    _ecfb2= ecfb2;
}
//--------------------------------------------------------------------------------------------------
void FatJet::addSubJet(const SubJet& sj) {_sjs.push_back(sj);}
//--------------------------------------------------------------------------------------------------
float     FatJet::bbt()       const{return _bbt;}
float     FatJet::tau1()      const{return _tau1;}
float     FatJet::tau2()      const{return _tau2;}
float     FatJet::tau3()      const{return _tau3;}
float     FatJet::tau2otau1() const{return _tau1 == 0 ? 99 : _tau2/_tau1;}
float     FatJet::tau3otau1() const{return _tau1 == 0 ? 99 : _tau3/_tau1;}
float     FatJet::tau3otau2() const{return _tau2 == 0 ? 99 : _tau3/_tau2;}
float     FatJet::ecfb1()     const{return _ecfb1;}
float     FatJet::ecfb2()     const{return _ecfb2;}
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
