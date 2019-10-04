
#include "DataFormats/interface/Jet.h"
#include "Configuration/interface/FillerConstants.h"

namespace TAna {

//--------------------------------------------------------------------------------------------------
// BaseRecoJet
//--------------------------------------------------------------------------------------------------
const MomentumF BaseRecoJet::rawMom() const { return p4()*toRawFactor();}
//--------------------------------------------------------------------------------------------------
bool Jet::passPUID()    const {
    return pt() > 50 || FillerConstants::doesPass(_jetID,FillerConstants::JETID_PU_T);}
//--------------------------------------------------------------------------------------------------
bool Jet::passTightID() const {
    return FillerConstants::doesPass(_jetID,FillerConstants::JETID_TIGHT);}
//--------------------------------------------------------------------------------------------------
bool Jet::passTightNoLepID() const {
    return FillerConstants::doesPass(_jetID,FillerConstants::JETID_TIGHTNOLEP);}
//--------------------------------------------------------------------------------------------------
void BaseRecoJet::addMCInfo(const  ASTypes::int8 hadronFlv,
        const  ASTypes::int8 partonFlv,const float JECUnc){
    _hadronFlv =hadronFlv;
    _partonFlv =partonFlv;
    _JECUnc    =JECUnc   ;
}
//--------------------------------------------------------------------------------------------------
void BaseRecoJet::addBTagging( const float deep_csv) {
    _deep_csv = deep_csv;
}

//--------------------------------------------------------------------------------------------------
// Jet
//--------------------------------------------------------------------------------------------------
void Jet::addMCInfo(const  ASTypes::int8 hadronFlv,
        const  ASTypes::int8 partonFlv,const float JECUnc, GenJet *gj){
    BaseRecoJet::addMCInfo(hadronFlv,partonFlv,JECUnc);
    _gj = gj;
}
void Jet::addBTagging( const float deep_csv, const float csv, const float deep_flavor){
    BaseRecoJet::addBTagging(deep_csv);
    _csv =csv;
    _deep_flavor =deep_flavor;
}

}
