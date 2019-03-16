

#include "Configuration/interface/ReaderConstants.h"
#include "DataFormats/interface/FatJet.h"
#include "DataFormats/interface/Electron.h"
#include "DataFormats/interface/Muon.h"
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

    paramSet.leptons.el_minPT   = 20  ;
    paramSet.leptons.el_maxETA  = 2.5 ;
    paramSet.leptons.el_maxDZ   = 0.1 ;
    paramSet.leptons.el_maxD0   = 0.05;
    paramSet.leptons.el_maxSip3D   = 4   ;
    paramSet.leptons.el_maxISO  = 0.1 ;
    paramSet.leptons.el_getID   = &Electron::passTightID_noIso;
    paramSet.leptons.el_getISO  = &Electron::miniIso;

    paramSet.leptons.mu_minPT   = 20  ;
    paramSet.leptons.mu_maxETA  = 2.4 ;
    paramSet.leptons.mu_maxDZ   = 0.1 ;
    paramSet.leptons.mu_maxD0   = 0.05;
    paramSet.leptons.mu_maxSip3D   = 4   ;
    paramSet.leptons.mu_maxISO  = 0.2 ;
    paramSet.leptons.mu_getID   = &Muon::passMedID;
    paramSet.leptons.mu_getISO  = &Muon::miniIso;


    return paramSet;
}
}
}
