#ifndef ANALYSISTREEMAKER_TREEFILLERS_FILLERCONSTANTS_H
#define ANALYSISTREEMAKER_TREEFILLERS_FILLERCONSTANTS_H

#include "AnalysisSupport/Utilities/interface/Types.h"
#include <vector>
#include <string>
namespace FillerConstants{

template <class storage>
void addPass(storage& passList, const ASTypes::size passed ) { passList |= (storage(1) << passed);}
template <class storage>
void removePass(storage& passList, const ASTypes::size passed ) { passList &= ~(storage(1) << passed);}
template <class storage>
bool doesPass(const storage passList, const ASTypes::size checkPassed ) {return  passList & (storage(1) << checkPassed);};

enum DataEra   {NOERA, ERA_2016,ERA_2017,ERA_2018};
const std::vector<std::string> DataEraNames = { "none","2016","2017","2018"};
enum DataRun   {NODATARUN, RUN2016A,RUN2016B,RUN2016C,RUN2016D,RUN2016E,RUN2016F,RUN2016G,RUN2016H,
                           RUN2017A,RUN2017B,RUN2017C,RUN2017D,RUN2017E,RUN2017F,RUN2017G,RUN2017H,
                           RUN2018A,RUN2018B,RUN2018C,RUN2018D,RUN2018E};
const std::vector<std::string> DataRunNames =
    {"NONE"," Run2016A","Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H",
        "Run2017A","Run2017B","Run2017C","Run2017D","Run2017E","Run2017F","Run2017G","Run2017H",
        "Run2018A","Run2018B","Run2018C","Run2018D","Run2018E"};
enum Dataset   {
                NODATASET,
                PD_BTagCSV        ,
                PD_BTagMu         ,
                PD_Charmonium     ,
                PD_DoubleEG       ,
                PD_DoubleMuon     ,
                PD_HTMHT          ,
                PD_JetHT          ,
                PD_MET            ,
                PD_MuOnia         ,
                PD_MuonEG         ,
                PD_SingleElectron ,
                PD_SingleMuon     ,
                PD_SinglePhoton   ,
                PD_Tau            ,
                PD_EGamma         ,
};
const std::vector<std::string> DatasetNames = { "NONE","BTagCSV","BTagMu","Charmonium","DoubleEG","DoubleMuon","HTMHT","JetHT","MET","MuOnia","MuonEG","SingleElectron","SingleMuon","SinglePhoton","Tau","EGamma"};
enum MCProcess {NOPROCESS, SIGNAL,TTBAR,WJETS,ZJETS,SINGLET,DIBOSON,TTX,QCD,HX};
const std::vector<std::string> MCProcessNames = { "NONE","signal","ttbar","wjets","zjets","singlet","diboson","ttX","qcd","hx"};



enum JetIDStatus { JETID_PU, JETID_LOOSE, JETID_TIGHT};

enum ElectronID {ELID_CUT_VETO
                ,ELID_CUT_LOOSE
                ,ELID_CUT_MED
                ,ELID_CUT_TIGHT
                ,ELID_CUT_HEEP
                ,ELID_CUT_NOISO_VETO
                ,ELID_CUT_NOISO_LOOSE
                ,ELID_CUT_NOISO_MED
                ,ELID_CUT_NOISO_TIGHT
                ,ELID_CUT_NOISO_HEEP
};

enum ElectronRECOStatus {ELRECO_TrckDrv
                        ,ELRECO_ECALDrv
};

enum MuonID   { MUID_SOFT
               ,MUID_LOOSE
               ,MUID_MED
               ,MUID_TIGHT
               ,MUID_HIGHPT
               ,MUID_MED16
};

//From
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
//https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
enum METFilters{
     Flag_HBHENoiseFilter
    ,Flag_HBHENoiseIsoFilter
    ,Flag_CSCTightHaloFilter
    ,Flag_CSCTightHaloTrkMuUnvetoFilter
    ,Flag_CSCTightHalo2015Filter
    ,Flag_globalTightHalo2016Filter
    ,Flag_globalSuperTightHalo2016Filter
    ,Flag_HcalStripHaloFilter
    ,Flag_hcalLaserEventFilter
    ,Flag_EcalDeadCellTriggerPrimitiveFilter
    ,Flag_EcalDeadCellBoundaryEnergyFilter
    ,Flag_ecalBadCalibFilter
    ,Flag_goodVertices
    ,Flag_eeBadScFilter
    ,Flag_ecalLaserCorrFilter
    ,Flag_trkPOGFilters
    ,Flag_chargedHadronTrackResolutionFilter
    ,Flag_muonBadTrackFilter
    ,Flag_BadChargedCandidateFilter
    ,Flag_BadPFMuonFilter
    ,Flag_BadChargedCandidateSummer16Filter
    ,Flag_BadPFMuonSummer16Filter
    ,Flag_trkPOG_manystripclus53X
    ,Flag_trkPOG_toomanystripclus53X
    ,Flag_trkPOG_logErrorTooManyClusters
    ,Flag_METFilters
    ,Flag_NFilters
};




enum Triggers_2017 {
 HLT17_PFHT500_PFMET100_PFMHT100_IDTight
,HLT17_PFHT700_PFMET85_PFMHT85_IDTight
,HLT17_PFHT800_PFMET75_PFMHT75_IDTight
,HLT17_AK8PFHT850_TrimMass50
,HLT17_AK8PFJet400_TrimMass30
,HLT17_AK8PFJet500
,HLT17_PFHT1050
,HLT17_PFMET120_PFMHT120_IDTight
,HLT17_PFMET120_PFMHT120_IDTight_PFHT60
,HLT17_PFMET140_PFMHT140_IDTight
,HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight
,HLT17_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60
,HLT17_PFMETNoMu140_PFMHTNoMu140_IDTight
,HLT17_PFMETTypeOne120_PFMHT120_IDTight
,HLT17_PFMETTypeOne120_PFMHT120_IDTight_PFHT60
,HLT17_PFMETTypeOne140_PFMHT140_IDTight
,HLT17_Ele115_CaloIdVT_GsfTrkIdT
,HLT17_Ele15_IsoVVVL_PFHT450
,HLT17_Ele28_eta2p1_WPTight_Gsf_HT150
,HLT17_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned
,HLT17_Ele32_WPTight_Gsf
,HLT17_Ele32_WPTight_Gsf_L1DoubleEG
,HLT17_Ele35_WPTight_Gsf
,HLT17_Ele50_CaloIdVT_GsfTrkIdT_PFJet165
,HLT17_IsoMu27
,HLT17_Mu15_IsoVVVL_PFHT450
,HLT17_Mu50
,HLT17_Photon200
,HLT17_NTrig
};




}


#endif
