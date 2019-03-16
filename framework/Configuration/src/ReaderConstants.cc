

#include "Configuration/interface/ReaderConstants.h"
#include "DataFormats/interface/FatJet.h"
namespace TAna {

namespace ReaderConstants{
ParameterSet set2017Parameters() {
    ParameterSet paramSet;

    paramSet.fatJets.cand_minPT     = 50                   ;
    paramSet.fatJets.cand_maxETA    = 2.4                  ;
    paramSet.fatJets.fjJetID        = &FatJet::passTightID ;
    paramSet.fatJets.sj_minPT       = 20 ;
    paramSet.fatJets.sj_maxETA      = 2.4;
    paramSet.fatJets.sj_minBTagPT    = 30 ;
    paramSet.fatJets.sj_maxBTagETA  = 2.4;
    paramSet.fatJets.wjj_maxLepDR   = 1.2    ;
    paramSet.fatJets.wjj_minPT      = 50     ;
    paramSet.fatJets.wjj_minSJs     = 2      ;
    paramSet.fatJets.hbb_minLepDPhi = 2.0    ;
    paramSet.fatJets.hbb_minPT      = 200    ;
    paramSet.fatJets.hbb_minSJs     = 2      ;
    paramSet.fatJets.hbbLL_minDphiBBLL = 2.0;
    paramSet.fatJets.hbbLL_minPT       = 200;
    paramSet.fatJets.hbbLL_minSJs      = 2;
    paramSet.fatJets.hbbLL_minDRbbLL   = 2.0;


    return paramSet;
}
}
}
