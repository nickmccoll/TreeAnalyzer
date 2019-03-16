
#include "DataFormats/interface/Muon.h"
#include "Configuration/interface/FillerConstants.h"
using namespace FillerConstants;
namespace TAna {
//--------------------------------------------------------------------------------------------------
void Muon::setMuonInfo( ASTypes::size id) {
    _id      =  id      ;
}
//--------------------------------------------------------------------------------------------------

bool  Muon::passSoftID () const {return doesPass(_id,MUID_SoftCutBasedId);}
bool  Muon::passLooseID() const {return doesPass(_id,MUID_CutBasedIdLoose);}
bool  Muon::passMedID  () const {return doesPass(_id,MUID_CutBasedIdMedium);}
bool  Muon::passTightID() const {return doesPass(_id,MUID_CutBasedIdTight);}
bool  Muon::passHighPT () const {return doesPass(_id,MUID_CutBasedIdGlobalHighPt);}
}
