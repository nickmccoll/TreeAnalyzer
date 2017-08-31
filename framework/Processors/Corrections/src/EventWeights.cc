
#include "Processors/Corrections/interface/EventWeights.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FillerConstants.h"
namespace TAna {
namespace EventWeights {
float calcNormalizedEventWeight(const EventReader& reader_event, const float cross, const float numE, const float lumi) {
    return (reader_event.weight > 0 ? 1.0 : -1.0) *lumi * cross *1000 /numE;
}
float getNormalizedEventWeight(const EventReader& reader_event, const float cross, const float numE, const float lumi) {
    if(reader_event.normWeightLoaded) return lumi*reader_event.normWeight;
    if(cross < 0 ||numE < 0) return reader_event.weight;
    return calcNormalizedEventWeight(reader_event,cross,numE,lumi);
}
float get4bXSecLimit(size mass) {
    static const std::vector<std::pair<size,float>> sigMassesAnd4bXSecLimit = {
            std::pair<size,float>(600, 120),
            std::pair<size,float>(800, 50),
            std::pair<size,float>(1000,20),
            std::pair<size,float>(1200,10),
            std::pair<size,float>(1400,6),
            std::pair<size,float>(1600,5),
            std::pair<size,float>(1800,4),
            std::pair<size,float>(2000,3),
            std::pair<size,float>(2500,2),
            std::pair<size,float>(3000,2),
            std::pair<size,float>(3500,2),
            std::pair<size,float>(4000,2),
            std::pair<size,float>(4500,2)
    };
    static const size nSigMassesAnd4bXSecLimit= sigMassesAnd4bXSecLimit.size();
    for(unsigned int iB = 0; iB < nSigMassesAnd4bXSecLimit; ++iB){
        if(mass == sigMassesAnd4bXSecLimit[iB].first) return sigMassesAnd4bXSecLimit[iB].second;
        if(mass < sigMassesAnd4bXSecLimit[iB].first  ){
            if(iB == 0) return  sigMassesAnd4bXSecLimit[iB].second;
            return sigMassesAnd4bXSecLimit[iB-1].second +
                    float(mass)*(sigMassesAnd4bXSecLimit[iB].second - sigMassesAnd4bXSecLimit[iB-1].second)
                    /float(sigMassesAnd4bXSecLimit[iB].first - sigMassesAnd4bXSecLimit[iB-1].first);
        }
    }
    return sigMassesAnd4bXSecLimit[nSigMassesAnd4bXSecLimit-1].second;
}


bool passEventFilters(const EventReader& reader_event) {
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_goodVertices) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_globalTightHalo2016Filter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_HBHENoiseFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_HBHENoiseIsoFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_EcalDeadCellTriggerPrimitiveFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_eeBadScFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::AnaTM_badMuons) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::AnaTM_badChargedHadrons) ) return false;
    if(reader_event.goodVtx == 0) return false;
    return true;
}
}

}



