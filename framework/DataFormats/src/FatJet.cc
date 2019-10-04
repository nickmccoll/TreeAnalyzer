
#include "DataFormats/interface/FatJet.h"

namespace TAna {
//--------------------------------------------------------------------------------------------------
void FatJet::addFJBtagging(const float bbt, const float deep_MDZHbb, const float deep_MDHbb
        , const float deep_Hbb){
    _bbt         = bbt          ;
    _deep_MDZHbb = deep_MDZHbb  ;
    _deep_MDHbb  = deep_MDHbb   ;
    _deep_Hbb    = deep_Hbb     ;

}
//--------------------------------------------------------------------------------------------------
void FatJet::addWTaging(const float deep_W){
    _deep_W=deep_W;
}
//--------------------------------------------------------------------------------------------------
void FatJet::addSubStructure(const float tau1, const float tau2,
        const float sdUp, const float sdDown){
    _tau1   =tau1  ;
    _tau2   =tau2  ;
    _sdUp   =sdUp  ;
    _sdDown =sdDown;
}
//--------------------------------------------------------------------------------------------------
void FatJet::addSubJet(const SubJet& sj) {_sjs.push_back(sj);}
//--------------------------------------------------------------------------------------------------
float     FatJet::bbt        ()       const{return _bbt        ;}
float     FatJet::deep_MDZHbb()       const{return _deep_MDZHbb;}
float     FatJet::deep_MDHbb ()       const{return _deep_MDHbb ;}
float     FatJet::deep_Hbb   ()       const{return _deep_Hbb   ;}
float     FatJet::deep_W     ()       const{return _deep_W   ;}
float     FatJet::tau1()      const{return _tau1;}
float     FatJet::tau2()      const{return _tau2;}
float     FatJet::tau2otau1() const{return _tau1 == 0 ? 99 : _tau2/_tau1;}
float     FatJet::rawSDMassUp()   const{return _sdUp;}
float     FatJet::rawSDMassDown() const{return _sdDown;}
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
