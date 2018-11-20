#ifndef ANALYSISTREEMAKER_TREEFILLERS_FILLERCONSTANTS_H
#define ANALYSISTREEMAKER_TREEFILLERS_FILLERCONSTANTS_H

#include "AnalysisSupport/Utilities/interface/Types.h"
#include <vector>
#include <string>
namespace FillerConstants{

template <class storage, class type>
void addPass(storage& passList, const type passed ) { passList |= passed;}
template <class storage, class type>
void removePass(storage& passList, const type passed ) { passList &= ~passed;}
template <class storage, class type>
bool doesPass(const storage passList, const type checkPassed ) {return  checkPassed & passList;};

enum DataEra   {NOERA, ERA_2016,ERA_2017,ERA_2018};
const std::vector<std::string> DataEraNames = { "none","2016","2017","2018"};
enum DataRun   {NODATARUN, RUN2016A,RUN2016B,RUN2016C,RUN2016D,RUN2016E,RUN2016F,RUN2016G,RUN2016H,
                           RUN2017A,RUN2017B,RUN2017C,RUN2017D,RUN2017E,RUN2017F,RUN2017G,RUN2017H,
                           RUN2018A,RUN2018B,RUN2018C,RUN2018D,RUN2018E};
const std::vector<std::string> DataRunNames =
    {"none"," Run2016A","Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H",
        "Run2017A","Run2017B","Run2017C","Run2017D","Run2017E","Run2017F","Run2017G","Run2017H",
        "Run2018A","Run2018B","Run2018C","Run2018D","Run2018E"};
enum Dataset   {NODATASET, SINGLEE,SINGLEMU, JETHT,MET};
const std::vector<std::string> DatasetNames = { "none","data_e","data_mu","data_jetht","data_met"};
enum MCProcess {NOPROCESS, SIGNAL,TTBAR,WJETS,ZJETS,SINGLET,DIBOSON,TTX,QCD,HX};
const std::vector<std::string> MCProcessNames = { "none","signal","ttbar","wjets","zjets","singlet","diboson","ttX","qcd","hx"};



enum JetIDStatus { JETID_PU = (1 << 0), JETID_LOOSE = (1 << 1), JETID_TIGHT  = (1 << 2)};

enum ElectronID {ELID_CUT_VETO         = (1 << 0)
                ,ELID_CUT_LOOSE        = (1 << 1)
                ,ELID_CUT_MED          = (1 << 2)
                ,ELID_CUT_TIGHT        = (1 << 3)
                ,ELID_CUT_HEEP         = (1 << 4)
                ,ELID_CUT_NOISO_VETO   = (1 << 5)
                ,ELID_CUT_NOISO_LOOSE  = (1 << 6)
                ,ELID_CUT_NOISO_MED    = (1 << 7)
                ,ELID_CUT_NOISO_TIGHT  = (1 << 8)
                ,ELID_CUT_NOISO_HEEP   = (1 << 9)
};

enum ElectronRECOStatus {ELRECO_TrckDrv         = (1 << 0)
                        ,ELRECO_ECALDrv         = (1 << 1)
};

enum MuonID   { MUID_SOFT   = (1 << 0)
               ,MUID_LOOSE  = (1 << 1)
               ,MUID_MED    = (1 << 2)
               ,MUID_TIGHT  = (1 << 3)
               ,MUID_HIGHPT = (1 << 4)
               ,MUID_MED16  = (1 << 5)
};

//From
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
//https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
enum METFilters{
    Flag_HBHENoiseFilter                     = (1 <<0)
    ,Flag_HBHENoiseIsoFilter                 = (1 <<1)
    ,Flag_CSCTightHaloFilter                 = (1 <<2)
    ,Flag_CSCTightHaloTrkMuUnvetoFilter      = (1 <<3)
    ,Flag_CSCTightHalo2015Filter             = (1 <<4)
    ,Flag_globalTightHalo2016Filter          = (1 <<5)
    ,Flag_globalSuperTightHalo2016Filter     = (1 <<6)
    ,Flag_HcalStripHaloFilter                = (1 <<7)
    ,Flag_hcalLaserEventFilter               = (1 <<8)
    ,Flag_EcalDeadCellTriggerPrimitiveFilter = (1 <<9)
    ,Flag_EcalDeadCellBoundaryEnergyFilter   = (1 <<10)
    ,Flag_ecalBadCalibFilter                 = (1 <<11)
    ,Flag_goodVertices                       = (1 <<12)
    ,Flag_eeBadScFilter                      = (1 <<13)
    ,Flag_ecalLaserCorrFilter                = (1 <<14)
    ,Flag_trkPOGFilters                      = (1 <<15)
    ,Flag_chargedHadronTrackResolutionFilter = (1 <<15)
    ,Flag_muonBadTrackFilter                 = (1 <<17)
    ,Flag_BadChargedCandidateFilter          = (1 <<18)
    ,Flag_BadPFMuonFilter                    = (1 <<19)
    ,Flag_BadChargedCandidateSummer16Filter  = (1 <<20)
    ,Flag_BadPFMuonSummer16Filter            = (1 <<21)
    ,Flag_trkPOG_manystripclus53X            = (1 <<22)
    ,Flag_trkPOG_toomanystripclus53X         = (1 <<23)
    ,Flag_trkPOG_logErrorTooManyClusters     = (1 <<24)
    ,Flag_METFilters                         = (1 <<25)
};
const std::vector<std::string> metFilterStrings = {
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_CSCTightHaloFilter",
        "Flag_CSCTightHaloTrkMuUnvetoFilter",
        "Flag_CSCTightHalo2015Filter",
        "Flag_globalTightHalo2016Filter",
        "Flag_globalSuperTightHalo2016Filter",
        "Flag_HcalStripHaloFilter",
        "Flag_hcalLaserEventFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_EcalDeadCellBoundaryEnergyFilter",
        "Flag_ecalBadCalibFilter",
        "Flag_goodVertices",
        "Flag_eeBadScFilter",
        "Flag_ecalLaserCorrFilter",
        "Flag_trkPOGFilters",
        "Flag_chargedHadronTrackResolutionFilter",
        "Flag_muonBadTrackFilter",
        "Flag_BadChargedCandidateFilter",
        "Flag_BadPFMuonFilter",
        "Flag_BadChargedCandidateSummer16Filter",
        "Flag_BadPFMuonSummer16Filter",
        "Flag_trkPOG_manystripclus53X",
        "Flag_trkPOG_toomanystripclus53X",
        "Flag_trkPOG_logErrorTooManyClusters",
        "Flag_METFilters"
};






enum Triggers{
    //single mu
     HLT_TkMu50                            = (ASTypes::size64(1) << 0)
    ,HLT_Mu50                              = (ASTypes::size64(1) << 1)
    ,HLT_Mu45_eta2p1                       = (ASTypes::size64(1) << 2)
    ,HLT_IsoMu24                           = (ASTypes::size64(1) << 3)
    ,HLT_IsoTkMu24                         = (ASTypes::size64(1) << 4)
    ,HLT_Mu15_IsoVVVL_PFHT350              = (ASTypes::size64(1) << 5)
    ,HLT_Mu15_IsoVVVL_PFHT400              = (ASTypes::size64(1) << 6)
    ,HLT_Mu15_IsoVVVL_PFHT600              = (ASTypes::size64(1) << 7)
    ,HLT_IsoMu16_eta2p1_MET30              = (ASTypes::size64(1) << 8)
    //single e
    ,HLT_Ele27_WPTight_Gsf                 = (ASTypes::size64(1) << 9)
    ,HLT_Ele45_WPLoose_Gsf                 = (ASTypes::size64(1) << 10)
    ,HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = (ASTypes::size64(1) << 11)
    ,HLT_Ele105_CaloIdVT_GsfTrkIdT         = (ASTypes::size64(1) << 12)
    ,HLT_Ele115_CaloIdVT_GsfTrkIdT         = (ASTypes::size64(1) << 13)
    ,HLT_Ele15_IsoVVVL_PFHT350             = (ASTypes::size64(1) << 14)
    ,HLT_Ele15_IsoVVVL_PFHT400             = (ASTypes::size64(1) << 15)
    ,HLT_Ele15_IsoVVVL_PFHT600             = (ASTypes::size64(1) << 16)
    //jetHT
    ,HLT_AK8PFHT650_TrimR0p1PT0p03Mass50   = (ASTypes::size64(1) << 17)
    ,HLT_AK8PFHT700_TrimR0p1PT0p03Mass50   = (ASTypes::size64(1) << 18)
    ,HLT_AK8PFJet320                       = (ASTypes::size64(1) << 19)
    ,HLT_AK8PFJet400                       = (ASTypes::size64(1) << 20)
    ,HLT_AK8PFJet450                       = (ASTypes::size64(1) << 21)
    ,HLT_AK8PFJet500                       = (ASTypes::size64(1) << 22)
    ,HLT_AK8PFJet360_TrimMass30            = (ASTypes::size64(1) << 23)
    ,HLT_AK8PFJet400_TrimMass30            = (ASTypes::size64(1) << 24)
    ,HLT_PFHT125                           = (ASTypes::size64(1) << 25)
    ,HLT_PFHT200                           = (ASTypes::size64(1) << 26)
    ,HLT_PFHT250                           = (ASTypes::size64(1) << 27)
    ,HLT_PFHT300                           = (ASTypes::size64(1) << 28)
    ,HLT_PFHT350                           = (ASTypes::size64(1) << 29)
    ,HLT_PFHT400                           = (ASTypes::size64(1) << 30)
    ,HLT_PFHT475                           = (ASTypes::size64(1) << 31)
    ,HLT_PFHT600                           = (ASTypes::size64(1) << 32)
    ,HLT_PFHT650_WideJetMJJ900DEtaJJ1p5    = (ASTypes::size64(1) << 33)
    ,HLT_PFHT650_WideJetMJJ950DEtaJJ1p5    = (ASTypes::size64(1) << 34)
    ,HLT_PFHT650                           = (ASTypes::size64(1) << 35)
    ,HLT_PFHT800                           = (ASTypes::size64(1) << 36)
    ,HLT_PFHT900                           = (ASTypes::size64(1) << 37)
    //met
    ,HLT_PFMETNoMu110_PFMHTNoMu110_IDTight = (ASTypes::size64(1) << 38)
    ,HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = (ASTypes::size64(1) << 39)
    ,HLT_PFMET110_PFMHT110_IDTight         = (ASTypes::size64(1) << 40)
    ,HLT_PFMET120_PFMHT120_IDTight         = (ASTypes::size64(1) << 41)
};

const std::vector<std::string> triggerStrings = {
        "HLT_TkMu50"                                ,
        "HLT_Mu50"                                  ,
        "HLT_Mu45_eta2p1"                           ,
        "HLT_IsoMu24"                               ,
        "HLT_IsoTkMu24"                             ,
        "HLT_Mu15_IsoVVVL_PFHT350"                  ,
        "HLT_Mu15_IsoVVVL_PFHT400"                  ,
        "HLT_Mu15_IsoVVVL_PFHT600"                  ,
        "HLT_IsoMu16_eta2p1_MET30"                  ,
        "HLT_Ele27_WPTight_Gsf"                     ,
        "HLT_Ele45_WPLoose_Gsf"                     ,
        "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"     ,
        "HLT_Ele105_CaloIdVT_GsfTrkIdT"             ,
        "HLT_Ele115_CaloIdVT_GsfTrkIdT"             ,
        "HLT_Ele15_IsoVVVL_PFHT350"                 ,
        "HLT_Ele15_IsoVVVL_PFHT400"                 ,
        "HLT_Ele15_IsoVVVL_PFHT600"                 ,
        "HLT_AK8PFHT650_TrimR0p1PT0p03Mass50"       ,
        "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50"       ,
        "HLT_AK8PFJet320"                           ,
        "HLT_AK8PFJet400"                           ,
        "HLT_AK8PFJet450"                           ,
        "HLT_AK8PFJet500"                           ,
        "HLT_AK8PFJet360_TrimMass30"                ,
        "HLT_AK8PFJet400_TrimMass30"                ,
        "HLT_PFHT125"                               ,
        "HLT_PFHT200"                               ,
        "HLT_PFHT250"                               ,
        "HLT_PFHT300"                               ,
        "HLT_PFHT350"                               ,
        "HLT_PFHT400"                               ,
        "HLT_PFHT475"                               ,
        "HLT_PFHT600"                               ,
        "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5"        ,
        "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5"        ,
        "HLT_PFHT650"                               ,
        "HLT_PFHT800"                               ,
        "HLT_PFHT900"                               ,
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight"     ,
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"     ,
        "HLT_PFMET110_PFMHT110_IDTight"             ,
        "HLT_PFMET120_PFMHT120_IDTight"

};




}


#endif
