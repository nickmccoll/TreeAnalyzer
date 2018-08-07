
#include "Processors/EventSelection/interface/EventSelection.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FillerConstants.h"
#include "DataFormats/interface/Lepton.h"

namespace TAna {
namespace EventSelection {


bool passEventFilters(const EventReader& reader_event) {
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_goodVertices) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_globalTightHalo2016Filter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_HBHENoiseFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_HBHENoiseIsoFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_EcalDeadCellTriggerPrimitiveFilter) ) return false;
    if(reader_event.realData && !FillerConstants::doesPass(reader_event.metFilters,FillerConstants::Flag_eeBadScFilter) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::AnaTM_badMuons) ) return false;
    if(!FillerConstants::doesPass(reader_event.metFilters,FillerConstants::AnaTM_badChargedHadrons) ) return false;
    if(reader_event.goodVtx == 0) return false;
    return true;
}

bool passTriggerSuite(const EventReader& reader_event   ) {
    auto passTrig = [&](FillerConstants::Triggers trig) -> bool {
        return FillerConstants::doesPass(reader_event.triggerAccepts,trig);
    };

    const bool passMuon = passTrig(FillerConstants::HLT_IsoMu24) || passTrig(FillerConstants::HLT_IsoTkMu24)
            || passTrig(FillerConstants::HLT_TkMu50)|| passTrig(FillerConstants::HLT_Mu50);
    const bool passEle = passTrig(FillerConstants::HLT_Ele27_WPTight_Gsf)
            || passTrig(FillerConstants::HLT_Ele45_WPLoose_Gsf)|| passTrig(FillerConstants::HLT_Ele115_CaloIdVT_GsfTrkIdT);
    const bool passJetHT = passTrig(FillerConstants::HLT_PFHT800) || passTrig(FillerConstants::HLT_PFHT900) ||
            passTrig(FillerConstants::HLT_AK8PFJet450) || passTrig(FillerConstants::HLT_AK8PFJet360_TrimMass30);
    const bool passMET = passTrig(FillerConstants::HLT_PFMETNoMu110_PFMHTNoMu110_IDTight) || passTrig(FillerConstants::HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);

    if(reader_event.realData){
        const bool passMCross = (reader_event.run < 274954) ?  passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
        const bool passECross = (reader_event.run < 274954) ?  passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
        switch(reader_event.dataset){
        case FillerConstants::SINGLEMU:
            return (passMuon || passMCross);
        case FillerConstants::SINGLEE:
            return (passEle || passECross) && !(passMuon || passMCross);
        case FillerConstants::JETHT:
            return passJetHT && !(passEle || passECross) && !(passMuon || passMCross);
        case FillerConstants::MET:
            return passMET && !passJetHT && !(passEle || passECross) && !(passMuon || passMCross);
        default:
            return false;
        }
    } else {
        const bool passMCross = passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
        const bool passECross = passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
        return (passMuon || passMCross) || (passEle || passECross)  || passJetHT || passMET;
    }
}

bool passMuonTriggerSuite(const EventReader& reader_event   ) {
    auto passTrig = [&](FillerConstants::Triggers trig) -> bool {
        return FillerConstants::doesPass(reader_event.triggerAccepts,trig);
    };

    const bool passMuon = passTrig(FillerConstants::HLT_IsoMu24) || passTrig(FillerConstants::HLT_IsoTkMu24)
            || passTrig(FillerConstants::HLT_TkMu50)|| passTrig(FillerConstants::HLT_Mu50);
    const bool passJetHT = passTrig(FillerConstants::HLT_PFHT800) || passTrig(FillerConstants::HLT_PFHT900) ||
            passTrig(FillerConstants::HLT_AK8PFJet450) || passTrig(FillerConstants::HLT_AK8PFJet360_TrimMass30);
    const bool passMET = passTrig(FillerConstants::HLT_PFMETNoMu110_PFMHTNoMu110_IDTight) || passTrig(FillerConstants::HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);

    if(reader_event.realData){
        const bool passMCross = (reader_event.run < 274954) ?  passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
        switch(reader_event.dataset){
        case FillerConstants::SINGLEMU:
            return (passMuon || passMCross);
        case FillerConstants::JETHT:
            return passJetHT  && !(passMuon || passMCross);
        case FillerConstants::MET:
            return passMET && !passJetHT  && !(passMuon || passMCross);
        default:
            return false;
        }
    } else {
        const bool passMCross = passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
        return (passMuon || passMCross)  || passJetHT || passMET;
    }
}

bool passElectronTriggerSuite(const EventReader& reader_event   ) {
    auto passTrig = [&](FillerConstants::Triggers trig) -> bool {
        return FillerConstants::doesPass(reader_event.triggerAccepts,trig);
    };

    const bool passEle = passTrig(FillerConstants::HLT_Ele27_WPTight_Gsf)
            || passTrig(FillerConstants::HLT_Ele45_WPLoose_Gsf)|| passTrig(FillerConstants::HLT_Ele115_CaloIdVT_GsfTrkIdT);
    const bool passJetHT = passTrig(FillerConstants::HLT_PFHT800) || passTrig(FillerConstants::HLT_PFHT900) ||
            passTrig(FillerConstants::HLT_AK8PFJet450) || passTrig(FillerConstants::HLT_AK8PFJet360_TrimMass30);
    const bool passMET = passTrig(FillerConstants::HLT_PFMETNoMu110_PFMHTNoMu110_IDTight) || passTrig(FillerConstants::HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);

    if(reader_event.realData){
        const bool passECross = (reader_event.run < 274954) ?  passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
        switch(reader_event.dataset){
        case FillerConstants::SINGLEE:
            return (passEle || passECross);
        case FillerConstants::JETHT:
            return passJetHT && !(passEle || passECross);
        case FillerConstants::MET:
            return passMET && !passJetHT && !(passEle || passECross);
        default:
            return false;
        }
    } else {
        const bool passECross = passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
        return  (passEle || passECross)  || passJetHT || passMET;
    }
}

bool passTriggerPreselection(const EventReader& reader_event,const float ht, const std::vector<const Lepton    *>& selectedLeptons ){
    if(!passTriggerSuite(reader_event)) return false;

    const float minHT = 400;
    const float minElePT = 30;
    const float minMuPT = 26;

    if(ht < minHT) return false;

    float maxElePT = 0;
    float maxMuPT = 0;
    for(const auto * l : selectedLeptons ) {
        if(l->isMuon()) maxMuPT = std::max(maxMuPT, l->pt());
        else maxElePT = std::max(maxElePT, l->pt());
    }

    if(maxElePT < minElePT && maxMuPT < minMuPT ) return false;
    return true;
};



}

}



