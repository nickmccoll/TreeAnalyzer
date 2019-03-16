#ifndef CONFIGURATION_INTERFACE_READERCONSTANTS_H
#define CONFIGURATION_INTERFACE_READERCONSTANTS_H

#include "AnalysisSupport/Utilities/interface/Types.h"
#include <vector>
#include <string>

// Forward declarations
namespace TAna {
class FatJet;
namespace FatJetSelHelpers {
typedef  bool (FatJet::*fjFunBool)() const;
}
}


namespace TAna {

struct FatJetParameters{
    //parameters to select jets on the loadFatJets func
    float cand_minPT     =-1;
    float cand_maxETA    =-1;
    FatJetSelHelpers::fjFunBool fjJetID =0 ;

    //parameters applied to both wjj and hbb AND dilep hbb
    float sj_minPT       =-1; //cand sel
    float sj_maxETA      =-1; //cand sel
    float sj_minBTagPT       =-1;  //For counting btags in pass sel
    float sj_maxBTagETA      =-1;  //For counting btags in pass sel

    //parameters applied to wjj candidate sel
    float wjj_maxLepDR   =-1;
    float wjj_minPT      =-1;
    float wjj_minSJs     =-1;

    //parameters applied to hbb candidate sel
    float hbb_minLepDPhi =-1;
    float hbb_minPT      =-1;
    float hbb_minSJs     =-1;

    //parameters applied to dilep Hbb candidate sel
    float hbbLL_minDphiBBLL = -1;
    float hbbLL_minPT       = -1;
    float hbbLL_minSJs      = -1;
    float hbbLL_minDRbbLL   = -1;
};


struct ParameterSet {
    FatJetParameters fatJets;
};

namespace ReaderConstants{
ParameterSet set2017Parameters();
}
}


#endif
