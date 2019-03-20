

#include "Configuration/interface/ReaderConstants.h"
#include "DataFormats/interface/FatJet.h"
#include "DataFormats/interface/Electron.h"
#include "DataFormats/interface/Muon.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Corrections/interface/BTagCalibrationStandalone.h"
namespace TAna {

namespace ReaderConstants{
ParameterSet set2017Parameters() {
    ParameterSet paramSet;
    paramSet.event.lumi = 41.53; //https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM
    paramSet.event.mcFilters = {
            FillerConstants::Flag_goodVertices,
            FillerConstants::Flag_globalSuperTightHalo2016Filter,
            FillerConstants::Flag_HBHENoiseFilter,
            FillerConstants::Flag_HBHENoiseIsoFilter,
            FillerConstants::Flag_EcalDeadCellTriggerPrimitiveFilter,
            FillerConstants::Flag_BadPFMuonFilter,
            FillerConstants::Flag_BadChargedCandidateFilter,
            FillerConstants::Flag_ecalBadCalibFilter
    };
    paramSet.event.dataFilters = paramSet.event.mcFilters;
    paramSet.event.dataFilters.push_back(FillerConstants::Flag_eeBadScFilter);
    paramSet.event.passTrigger = &EventSelection::passTriggerSuite2017;
    paramSet.event.minHT         =400;
    paramSet.event.minTriggerEl  =30;
    paramSet.event.minTriggerMu  =26;


    paramSet.fatJets.cand_minPT     = 50                   ;
    paramSet.fatJets.cand_maxETA    = 2.4                  ;
    paramSet.fatJets.fjJetID        = &FatJet::passTightID ;
    paramSet.fatJets.sj_minPT       = 20 ;
    paramSet.fatJets.sj_maxETA      = 2.4;
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



    paramSet.jets.minJetPT      = 30;
    paramSet.jets.maxJetETA     = 2.4;
    paramSet.jets.passJetID     = &Jet::passTightID;

    paramSet.jets.CSV_WP        = {-100,0.5803,0.8838,0.9693};
    paramSet.jets.DeepCSV_WP    = {-100,0.1522,0.4941,0.8001};
    paramSet.jets.minBtagJetPT  = 30;
    paramSet.jets.maxBTagJetETA = 2.4;

    paramSet.jets.getJetBTagVal = &BaseRecoJet::csv;
    paramSet.jets.jetBTagWP     = paramSet.jets.CSV_WP [BTagging::BTAG_M];
    paramSet.jets.getSJBTagVal =&BaseRecoJet::csv;
    paramSet.jets.sjBTagLWP = paramSet.jets.CSV_WP [BTagging::BTAG_L];
    paramSet.jets.sjBTagMWP = paramSet.jets.CSV_WP [BTagging::BTAG_M];

    paramSet.jets.jetBtagCorrSFFile ="corrections/CSVv2_Moriond17_B_H.csv";
    paramSet.jets.jetBtagCorrEffFile ="corrections/ak4_csvEff.root";
    paramSet.jets.jetBtagCorrWP       = paramSet.jets.CSV_WP      ;
    paramSet.jets.jetBtagCorrGetBTagVal = paramSet.jets.getJetBTagVal;
    paramSet.jets.sjBtagCorrWP         = paramSet.jets.CSV_WP      ;
    paramSet.jets.sjBtagCorrGetBTagVal = paramSet.jets.getSJBTagVal;
    paramSet.jets.sjBtagCorrSFFile     = "corrections/subjet_CSVv2_Moriond17_B_H.csv";
    paramSet.jets.sjBtagCorrEffFile    = "corrections/sj_csvEff.root";


    return paramSet;
}
}
}
