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
class Muon;
class Electron;
namespace LepSelHelpers {
    typedef  bool (Muon::*muFunBool)() const;
    typedef  bool (Electron::*elFunBool)() const;
    typedef  float (Muon::*muFunFloat)() const;
    typedef  float (Electron::*elFunFloat)() const;
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

struct LeptonParameters {
    float el_minPT  ;
    float el_maxETA ;
    float el_maxDZ  ;
    float el_maxD0  ;
    float el_maxSip3D  ;
    float el_maxISO ;
    LepSelHelpers::elFunBool el_getID ;
    LepSelHelpers::elFunFloat el_getISO;

    float mu_minPT  ;
    float mu_maxETA ;
    float mu_maxDZ  ;
    float mu_maxD0  ;
    float mu_maxSip3D  ;
    float mu_maxISO ;
    LepSelHelpers::muFunBool mu_getID ;
    LepSelHelpers::muFunFloat mu_getISO;
};


struct ParameterSet {
    FatJetParameters fatJets;
    LeptonParameters leptons;
};

namespace ReaderConstants{
ParameterSet set2017Parameters();
}
}


#endif
