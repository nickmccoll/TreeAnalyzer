
#include "DataFormats/interface/Jet.h"
#include "Configuration/interface/FillerConstants.h"

namespace TAna {

//--------------------------------------------------------------------------------------------------
const MomentumF BaseRecoJet::rawMom() const { return p4()*toRawFactor();}
bool Jet::passPUID()    const {return FillerConstants::doesPass(_jetID,FillerConstants::JETID_PU_T);}
bool Jet::passLooseID() const {return FillerConstants::doesPass(_jetID,FillerConstants::JETID_LOOSE);}
bool Jet::passTightID() const {return FillerConstants::doesPass(_jetID,FillerConstants::JETID_TIGHT);}
bool Jet::passTightNoLepID() const {return FillerConstants::doesPass(_jetID,FillerConstants::JETID_TIGHTNOLEP);}
void BaseRecoJet::addMCInfo(const  ASTypes::int8 hadronFlv,
        const  ASTypes::int8 partonFlv,const float JECUnc){
    _hadronFlv =hadronFlv;
    _partonFlv =partonFlv;
    _JECUnc    =JECUnc   ;
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
void Jet::addMCInfo(const  ASTypes::int8 hadronFlv,
        const  ASTypes::int8 partonFlv,const float JECUnc, GenJet *gj){
    BaseRecoJet::addMCInfo(hadronFlv,partonFlv,JECUnc);
    _gj = gj;
}
void Jet::addExtraInfo(const float jetID){
    _jetID = jetID;
}

}
