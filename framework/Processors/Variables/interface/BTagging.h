
#ifndef PROCESSORS_VARIABLES_BTAGGING_H
#define PROCESSORS_VARIABLES_BTAGGING_H

#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
namespace TAna {

namespace BTagging{
enum  CSVWP { CSV_INCL, CSV_L, CSV_M, CSV_T};
const float CSVWP_VALS[] = {-100,0.5426,0.8484,0.9535};

template<typename Jet>
bool isLooseCSVTagged(const Jet& jet) { return jet.csv() >=  CSVWP_VALS[CSV_L];}

template<typename Jet>
bool isMediumCSVTagged(const Jet& jet) { return jet.csv() >=  CSVWP_VALS[CSV_M];}

template<typename Jet>
bool isTightCSVTagged(const Jet& jet) { return jet.csv() >=  CSVWP_VALS[CSV_T];}

template<typename Jet>
bool csvTagged(const Jet& jet, CSVWP wp) { return jet.csv() >=  CSVWP_VALS[wp];}

template<typename Jet>
bool csvNotTagged(const Jet& jet, CSVWP wp) { return jet.csv() <  CSVWP_VALS[wp];}


enum CSVSJ_CAT {CSVSJ_INCL, //inclusive
                CSVSJ_FF,  // no Loose b-tags
                CSVSJ_LF,  // One loose, not med
                CSVSJ_LL,  // Two loose, not med
                CSVSJ_MF,  // One med, not loose
                CSVSJ_ML,  // One med, one loose not med
                CSVSJ_MM,  // Two med
};

template<typename Jet>
CSVSJ_CAT getCSVSJCat(const std::vector<Jet>& subjets, const float minPT, const float maxAETA){
    int nMB = 0;
    int nLB = 0;
    for(const auto& j : subjets) {
        if(j.pt() < minPT) continue;
        if(maxAETA > 0 && j.absEta() > maxAETA) continue;
        if(isLooseCSVTagged(j)) nLB++ ;
        if(isMediumCSVTagged(j)) nMB++ ;
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

enum  BBTWP { BBT_INCL, BBT_L, BBT_M1, BBT_M2, BBT_T};
const float BBT_VALS[] = {-100,0.3,0.6,0.8,0.9};
}


}

#endif

