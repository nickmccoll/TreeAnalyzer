
#include "Processors/EventSelection/interface/EventSelection.h"
#include "TreeReaders/interface/EventReader.h"
#include "Configuration/interface/FillerConstants.h"
#include "DataFormats/interface/Lepton.h"

using namespace FillerConstants;

namespace TAna {
namespace EventSelection {
//--------------------------------------------------------------------------------------------------
bool passTriggerSuite2017(const EventReader& reader_event   ) {
    auto passTrig = [&](FillerConstants::Triggers_2017 trig) -> bool {
        return FillerConstants::doesPass(*reader_event.triggerAccepts,trig);
    };

    const bool passMuon = passTrig(HLT17_IsoMu27) || passTrig(HLT17_Mu50)
            || passTrig(HLT17_Mu15_IsoVVVL_PFHT450);
    const bool passEle =
              passTrig(HLT17_Ele32_WPTight_Gsf)
            || passTrig(HLT17_Ele32_WPTight_Gsf_L1DoubleEG)
            || passTrig(HLT17_Ele35_WPTight_Gsf)
            || passTrig(HLT17_Ele115_CaloIdVT_GsfTrkIdT)
            || passTrig(HLT17_Ele15_IsoVVVL_PFHT450)
            || passTrig(HLT17_Ele28_eta2p1_WPTight_Gsf_HT150)
            || passTrig(HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned)
            || passTrig(HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);

    const bool passMET = passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight)
            || passTrig(HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)
            || passTrig(HLT17_PFMETNoMu140_PFMHTNoMu140_IDTight)
			|| passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight)
			|| passTrig(HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60)
			|| passTrig(HLT17_PFMETTypeOne140_PFMHT140_IDTight)
			|| passTrig(HLT17_PFMET120_PFMHT120_IDTight)
			|| passTrig(HLT17_PFMET120_PFMHT120_IDTight_PFHT60)
			|| passTrig(HLT17_PFMET140_PFMHT140_IDTight);

    const bool passJetHT = passTrig(HLT17_PFHT1050) || passTrig(HLT17_AK8PFJet500) ||
            passTrig(HLT17_AK8PFJet400_TrimMass30) || passTrig(HLT17_AK8PFHT850_TrimMass50);
    const bool passPhoton = passTrig(HLT17_Photon200);

    if(reader_event.realData){
        switch(*reader_event.dataset){
        case PD_SingleMuon:
            return passMuon;
        case PD_SingleElectron:
            return passEle && !(passMuon);
        case PD_JetHT:
            return passJetHT && !(passEle || passMuon);
        case PD_MET:
            return passMET && !(passJetHT ||passEle || passMuon);
        case PD_SinglePhoton:
            return passPhoton && !(passMET || passJetHT ||passEle || passMuon);
        default:
            throw std::invalid_argument(
                    "EventSelection::passTriggerSuite2017 -> "
                    "This is not a valid primary dataset");
        }
    } else {
        return (passMuon || passEle  || passJetHT || passMET || passPhoton);
    }
    return true;
}
//--------------------------------------------------------------------------------------------------
bool alwaysTrue(const EventReader& reader_event   ) {
    if(*reader_event.dataEra){}
    return true;
}
//--------------------------------------------------------------------------------------------------
bool passEventFilters(const EventParameters& params,const EventReader& reader_event) {
    for(const auto& filter : reader_event.realData  ? params.dataFilters : params.mcFilters){
        if(!FillerConstants::doesPass(*reader_event.metFilters,filter) ) return false;
    }
    if(*(reader_event.goodVtx) == 0) return false;
    return true;
}
//--------------------------------------------------------------------------------------------------
bool passTriggerPreselection(const EventParameters& params,
        const EventReader& reader_event,const float ht,
        const std::vector<const Lepton    *>& selectedLeptons ){
    if(!(*params.passTrigger)(reader_event)) return false;

    if(ht < params.minHT) return false;

    float maxElePT = 0;
    float maxMuPT = 0;
    for(const auto * l : selectedLeptons ) {
        if(l->isMuon()) maxMuPT = std::max(maxMuPT, l->pt());
        else maxElePT = std::max(maxElePT, l->pt());
    }

    if(maxElePT < params.minTriggerEl && maxMuPT < params.minTriggerMu ) return false;
    return true;
};



}

}


//2016 ->
//bool passTriggerSuite(const EventReader& reader_event   ) {
//    auto passTrig = [&](FillerConstants::Triggers trig) -> bool {
//        return FillerConstants::doesPass(*reader_event.triggerAccepts,trig);
//    };
//
//    const bool passMuon = passTrig(FillerConstants::HLT_IsoMu24) || passTrig(FillerConstants::HLT_IsoTkMu24)
//            || passTrig(FillerConstants::HLT_TkMu50)|| passTrig(FillerConstants::HLT_Mu50);
//    const bool passEle = passTrig(FillerConstants::HLT_Ele27_WPTight_Gsf)
//            || passTrig(FillerConstants::HLT_Ele45_WPLoose_Gsf)|| passTrig(FillerConstants::HLT_Ele115_CaloIdVT_GsfTrkIdT);
//    const bool passJetHT = passTrig(FillerConstants::HLT_PFHT800) || passTrig(FillerConstants::HLT_PFHT900) ||
//            passTrig(FillerConstants::HLT_AK8PFJet450) || passTrig(FillerConstants::HLT_AK8PFJet360_TrimMass30);
//    const bool passMET = passTrig(FillerConstants::HLT_PFMETNoMu110_PFMHTNoMu110_IDTight) || passTrig(FillerConstants::HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
//
//    if(reader_event.realData){
//        const bool passMCross = (*reader_event.run < 274954) ?  passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
//        const bool passECross = (*reader_event.run < 274954) ?  passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
//        switch(*reader_event.dataset){
//        case FillerConstants::SINGLEMU:
//            return (passMuon || passMCross);
//        case FillerConstants::SINGLEE:
//            return (passEle || passECross) && !(passMuon || passMCross);
//        case FillerConstants::JETHT:
//            return passJetHT && !(passEle || passECross) && !(passMuon || passMCross);
//        case FillerConstants::MET:
//            return passMET && !passJetHT && !(passEle || passECross) && !(passMuon || passMCross);
//        default:
//            return false;
//        }
//    } else {
//        const bool passMCross = passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
//        const bool passECross = passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
//        return (passMuon || passMCross) || (passEle || passECross)  || passJetHT || passMET;
//    }
//    return true;
//}

//bool passMuonTriggerSuite(const EventReader& reader_event   ) {
//    auto passTrig = [&](FillerConstants::Triggers trig) -> bool {
//        return FillerConstants::doesPass(*reader_event.triggerAccepts,trig);
//    };
//
//    const bool passMuon = passTrig(FillerConstants::HLT_IsoMu24) || passTrig(FillerConstants::HLT_IsoTkMu24)
//            || passTrig(FillerConstants::HLT_TkMu50)|| passTrig(FillerConstants::HLT_Mu50);
//    const bool passJetHT = passTrig(FillerConstants::HLT_PFHT800) || passTrig(FillerConstants::HLT_PFHT900) ||
//            passTrig(FillerConstants::HLT_AK8PFJet450) || passTrig(FillerConstants::HLT_AK8PFJet360_TrimMass30);
//    const bool passMET = passTrig(FillerConstants::HLT_PFMETNoMu110_PFMHTNoMu110_IDTight) || passTrig(FillerConstants::HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
//
//    if(reader_event.realData){
//        const bool passMCross = (*reader_event.run < 274954) ?  passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
//        switch(*reader_event.dataset){
//        case FillerConstants::SINGLEMU:
//            return (passMuon || passMCross);
//        case FillerConstants::JETHT:
//            return passJetHT  && !(passMuon || passMCross);
//        case FillerConstants::MET:
//            return passMET && !passJetHT  && !(passMuon || passMCross);
//        default:
//            return false;
//        }
//    } else {
//        const bool passMCross = passTrig(FillerConstants::HLT_Mu15_IsoVVVL_PFHT400);
//        return (passMuon || passMCross)  || passJetHT || passMET;
//    }
//    return true;
//}

//bool passElectronTriggerSuite(const EventReader& reader_event   ) {
//    auto passTrig = [&](FillerConstants::Triggers trig) -> bool {
//        return FillerConstants::doesPass(*reader_event.triggerAccepts,trig);
//    };
//
//    const bool passEle = passTrig(FillerConstants::HLT_Ele27_WPTight_Gsf)
//            || passTrig(FillerConstants::HLT_Ele45_WPLoose_Gsf)|| passTrig(FillerConstants::HLT_Ele115_CaloIdVT_GsfTrkIdT);
//    const bool passJetHT = passTrig(FillerConstants::HLT_PFHT800) || passTrig(FillerConstants::HLT_PFHT900) ||
//            passTrig(FillerConstants::HLT_AK8PFJet450) || passTrig(FillerConstants::HLT_AK8PFJet360_TrimMass30);
//    const bool passMET = passTrig(FillerConstants::HLT_PFMETNoMu110_PFMHTNoMu110_IDTight) || passTrig(FillerConstants::HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
//
//    if(reader_event.realData){
//        const bool passECross = (*reader_event.run < 274954) ?  passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT350):   passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
//        switch(*reader_event.dataset){
//        case FillerConstants::SINGLEE:
//            return (passEle || passECross);
//        case FillerConstants::JETHT:
//            return passJetHT && !(passEle || passECross);
//        case FillerConstants::MET:
//            return passMET && !passJetHT && !(passEle || passECross);
//        default:
//            return false;
//        }
//    } else {
//        const bool passECross = passTrig(FillerConstants::HLT_Ele15_IsoVVVL_PFHT400);
//        return  (passEle || passECross)  || passJetHT || passMET;
//    }
//    return true;
//}

