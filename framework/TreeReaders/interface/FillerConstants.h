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


enum DataRun   {NODATARUN, RUN2016A,RUN2016B,RUN2016C,RUN2016D,RUN2016E,RUN2016F,RUN2016G,RUN2016H};
enum Dataset   {NODATASET, SINGLEE,SINGLEMU, JETHT,MET};
const std::vector<std::string> DatasetNames = { "none","singlee","singlemu","jetht","met"};
enum MCProcess {NOPROCESS, SIGNAL,TTBAR,WJETS,ZJETS,SINGLET,DIBOSON,TTX,QCD};
const std::vector<std::string> MCProcessNames = { "none","signal","ttbar","wjets","zjets","singlet","diboson","ttX","qcd"};



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
enum MuonID   { MUID_SOFT   = (1 << 0)
               ,MUID_LOOSE  = (1 << 1)
               ,MUID_MED    = (1 << 2)
               ,MUID_TIGHT  = (1 << 3)
               ,MUID_HIGHPT = (1 << 4)
               ,MUID_MED16  = (1 << 5)
};

enum METFilters{
     Flag_goodVertices                      = (1 << 0)
    ,Flag_globalTightHalo2016Filter         = (1 << 1)
    ,Flag_HBHENoiseFilter                   = (1 << 2)
    ,Flag_HBHENoiseIsoFilter                = (1 << 3)
    ,Flag_EcalDeadCellTriggerPrimitiveFilter= (1 << 4)
    ,Flag_eeBadScFilter                     = (1 << 5)
    ,Flag_CSCTightHaloFilter                = (1 << 6)
    ,Flag_CSCTightHalo2015Filter            = (1 << 7)
    ,Flag_trackingFailureFilter             = (1 << 8)
    ,Flag_trkPOGFilters                     = (1 << 9)
    ,Flag_ecalLaserCorrFilter               = (1 << 10)
    ,Flag_hcalLaserEventFilter              = (1 << 11)
    ,Flag_badMuons                          = (1 << 12)
    ,Flag_duplicateMuons                    = (1 << 13)
     //Ones we add ourselves
     //bad muon filters
    ,AnaTM_badMuons                         = (1 << 14)
    ,AnaTM_badChargedHadrons                = (1 << 15)
    //ECAL slew rate
    ,AnaTM_dupECALClusters                  = (1 << 16)  //true if duplicates are present..bad
    ,AnaTM_hitsNotReplaced                  = (1 << 17)  //true of not empty...bad
};
const std::vector<std::string> metFilterStrings = {
        "Flag_goodVertices",
        "Flag_globalTightHalo2016Filter",
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_eeBadScFilter",
        "Flag_CSCTightHaloFilter",
        "Flag_CSCTightHalo2015Filter",
        "Flag_trackingFailureFilter",
        "Flag_trkPOGFilters",
        "Flag_ecalLaserCorrFilter",
        "Flag_hcalLaserEventFilter",
        "Flag_badMuons",
        "Flag_duplicateMuons",
        "AnaTM_badMuons",
        "AnaTM_badChargedHadrons",
        "AnaTM_dupECALClusters",
        "AnaTM_hitsNotReplaced"
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
