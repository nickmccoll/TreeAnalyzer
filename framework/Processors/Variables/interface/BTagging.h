
#ifndef PROCESSORS_VARIABLES_BTAGGING_H
#define PROCESSORS_VARIABLES_BTAGGING_H

#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Configuration/interface/ReaderConstants.h"
namespace TAna {

namespace BTagging{
enum  BTAGWP { BTAG_INCL, BTAG_L, BTAG_M, BTAG_T};

template<typename Jet>
bool passJetBTagWP(const JetParameters& params, const Jet& jet){
    return (jet.*params.getJetBTagVal)() >= params.jetBTagWP ;
}

template<typename Jet>
bool passSubjetBTagLWP(const JetParameters& params, const Jet& subjet){
    return (subjet.*params.getSJBTagVal)() >= params.sjBTagLWP ;
}

template<typename Jet>
bool passSubjetBTagMWP(const JetParameters& params, const Jet& subjet){
    return (subjet.*params.getSJBTagVal)() >= params.sjBTagMWP ;
}


enum CSVSJ_CAT {CSVSJ_INCL, //inclusive
                CSVSJ_FF,  // no Loose b-tags
                CSVSJ_LF,  // One loose, not med
                CSVSJ_LL,  // Two loose, not med
                CSVSJ_MF,  // One med, not loose
                CSVSJ_ML,  // One med, one loose not med
                CSVSJ_MM,  // Two med
};

template<typename Jet>
CSVSJ_CAT getCSVSJCat(const JetParameters& params, const std::vector<Jet>& subjets){
    int nMB = 0;
    int nLB = 0;
    for(const auto& j : subjets) {
        if(j.pt() < params.minBtagJetPT) continue;
        if(params.maxBTagJetETA > 0 && j.absEta() > params.maxBTagJetETA) continue;
        if(passSubjetBTagLWP(params,j)) nLB++ ;
        if(passSubjetBTagMWP(params,j)) nMB++ ;
    }

    if(nLB == 0) return CSVSJ_FF;
    if(nMB == 0){
        if(nLB == 1) return CSVSJ_LF;
        return CSVSJ_LL;
    }
    if(nMB == 1) {
        if(nLB == 1) return CSVSJ_MF;
        return CSVSJ_ML;
    }
    return CSVSJ_MM;
};

enum  FLAVOR { FLV_B, FLV_C, FLV_L};

template<typename Jet>
FLAVOR jetFlavor(const Jet& jet) {
    switch(jet.hadronFlv()){
    case 5:
        return FLV_B;
    case 4:
        return FLV_C;
    default:
        return FLV_L;
    }
}
}


}

#endif

