

#include "Configuration/interface/ReaderConstants.h"
#include "DataFormats/interface/FatJet.h"
#include "DataFormats/interface/Electron.h"
#include "DataFormats/interface/Muon.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Corrections/interface/BTagCalibrationStandalone.h"
namespace TAna {

namespace ReaderConstants{
ParameterSet setCommonParameters() {
    ParameterSet paramSet;
    paramSet.event.lumi = 1; //https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM
    paramSet.event.mcFilters = {
            FillerConstants::Flag_goodVertices,
            FillerConstants::Flag_globalSuperTightHalo2016Filter,
            FillerConstants::Flag_HBHENoiseFilter,
            FillerConstants::Flag_HBHENoiseIsoFilter,
            FillerConstants::Flag_EcalDeadCellTriggerPrimitiveFilter,
            FillerConstants::Flag_BadPFMuonFilter
    };
    paramSet.event.dataFilters = paramSet.event.mcFilters;
    paramSet.event.dataFilters.push_back(FillerConstants::Flag_eeBadScFilter);
    paramSet.event.passTrigger = &EventSelection::passTriggerSuite2017;
    paramSet.event.minHT         =400;
    paramSet.event.minTriggerEl  =30;
    paramSet.event.minTriggerMu  =27;
    paramSet.event.leptonCorrSFFile = "corrections/triggerSF_17.root";
    paramSet.event.puCorrSFFile = "corrections/puSF_17.root";
    paramSet.event.ttbarXSecSF_1000toInf_nLep0 = 0.000718573;
    paramSet.event.ttbarXSecSF_1000toInf_nLep1 = 0.000739261;
    paramSet.event.ttbarXSecSF_1000toInf_nLep2 = 0.000523889;
    paramSet.event.ttbarXSecSF_700to1000_nLep0 = 0.00123404;
    paramSet.event.ttbarXSecSF_700to1000_nLep1 = 0.00129566;
    paramSet.event.ttbarXSecSF_700to1000_nLep2 = 0.000767859;
    paramSet.event.doTTBarStitching = true;

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

    paramSet.leptons.el_minPT   = 30  ;
    paramSet.leptons.el_maxETA  = 1.479 ;
    paramSet.leptons.el_maxDZ   = 0.1 ;
    paramSet.leptons.el_maxD0   = 0.05;
    paramSet.leptons.el_maxSip3D   = 4   ;
    paramSet.leptons.el_maxISO  = 0.2 ;
    paramSet.leptons.el_getID   = &Electron::passMVA90ID_noIso;
    paramSet.leptons.el_getISO  = &Electron::miniIso;

    paramSet.leptons.mu_minPT   = 27  ;
    paramSet.leptons.mu_maxETA  = 2.4 ;
    paramSet.leptons.mu_maxDZ   = 0.1 ;
    paramSet.leptons.mu_maxD0   = 0.05;
    paramSet.leptons.mu_maxSip3D   = 4   ;
    paramSet.leptons.mu_maxISO  = 0.2 ;
    paramSet.leptons.mu_getID   = &Muon::passMedID;
    paramSet.leptons.mu_getISO  = &Muon::miniIso;

    paramSet.dileptons.el_minPT   = 10  ;
    paramSet.dileptons.el_maxETA  = 2.5 ;
    paramSet.dileptons.el_maxDZ   = 0.1 ;
    paramSet.dileptons.el_maxD0   = 0.05;
    paramSet.dileptons.el_maxSip3D   = 9999   ;
    paramSet.dileptons.el_maxISO  = 0.15 ;
    paramSet.dileptons.el_getID1   = &Electron::passMedID_noIso;
    paramSet.dileptons.el_getID2   = &Electron::passMedID_noIso;
    paramSet.dileptons.el_getISO  = &Electron::miniIso;

    paramSet.dileptons.mu_minPT   = 10  ;
    paramSet.dileptons.mu_maxETA  = 2.4 ;
    paramSet.dileptons.mu_maxDZ   = 0.1 ;
    paramSet.dileptons.mu_maxD0   = 0.05;
    paramSet.dileptons.mu_maxSip3D   = 9999   ;
    paramSet.dileptons.mu_maxISO  = 0.15 ;
    paramSet.dileptons.mu_getID1   = &Muon::passLooseID;
    paramSet.dileptons.mu_getID2   = &Muon::passLooseID;
    paramSet.dileptons.mu_getISO  = &Muon::miniIso;

    paramSet.jets.minJetPT      = 20;
    paramSet.jets.maxJetETA     = 2.4;
    paramSet.jets.passJetID     = &Jet::passTightID;

    paramSet.jets.CSV_WP        = {-100,0.5803,0.8838,0.9693};
    paramSet.jets.DeepCSV_WP    = {-100,0.1522,0.4941,0.8001};
    paramSet.jets.minBtagJetPT  = 20;
    paramSet.jets.maxBTagJetETA = 2.4;

    paramSet.jets.getJetBTagVal = &BaseRecoJet::deep_csv;
    paramSet.jets.jetBTagWP     = paramSet.jets.DeepCSV_WP [BTagging::BTAG_M];
    paramSet.jets.getSJBTagVal =&BaseRecoJet::deep_csv;
    paramSet.jets.sjBTagLWP = paramSet.jets.DeepCSV_WP [BTagging::BTAG_L];
    paramSet.jets.sjBTagMWP = paramSet.jets.DeepCSV_WP [BTagging::BTAG_M];

    paramSet.jets.jetBtagCorrSFFile ="corrections/CSVv2_Moriond17_B_H.csv";
    paramSet.jets.jetBtagCorrEffFile ="corrections/ak4_csvEff.root";
    paramSet.jets.jetBtagCorrWP       = paramSet.jets.DeepCSV_WP      ;
    paramSet.jets.jetBtagCorrGetBTagVal = paramSet.jets.getJetBTagVal;
    paramSet.jets.sjBtagCorrWP         = paramSet.jets.DeepCSV_WP      ;
    paramSet.jets.sjBtagCorrGetBTagVal = paramSet.jets.getSJBTagVal;
    paramSet.jets.sjBtagCorrSFFile     = "corrections/subjet_CSVv2_Moriond17_B_H.csv";
    paramSet.jets.sjBtagCorrEffFile    = "corrections/sj_csvEff.root";

    paramSet.hww.posMETParErr     = 0.067;
    paramSet.hww.negMETParErr     = 0.16;
    paramSet.hww.metPerpErr       = 31;
    paramSet.hww.jetErr           = 0.11;
    paramSet.hww.onWlnuMeanJet    = 30;
    paramSet.hww.offWlnuMeanJet   = 80;
    paramSet.hww.onWlnuMeanWlnu   = 80;
    paramSet.hww.offWlnuMeanWlnu  = 41;
    paramSet.hww.offWlnuPosWlnuErr= 5;
    paramSet.hww.offWnluNegWlnuErr= 16;
    paramSet.hww.onWlnuWlnuErr    = 2.3;
    paramSet.hww.onWlnuHWWErr     = 8.3;
    paramSet.hww.offWlnuHWWErr    = 3;

    paramSet.hww.ptCorB        = 0.9663;
    paramSet.hww.ptCorM        = -0.00001013;
    paramSet.hww.liFileName    ="variables/hhSol_templates.root";
    paramSet.hww.bkgLiFileName ="variables/hhSol_bkgTemplates.root";

    paramSet.hww.dilepInvMassGuess = 55;

    return paramSet;
}
ParameterSet set2016Parameters() {
    ParameterSet paramSet = setCommonParameters();
    paramSet.event.lumi = 35.9;
    paramSet.jets.jer_AK4CHS_resFile    ="corrections/JER/Summer16_25nsV1_MC/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
    paramSet.jets.jer_AK4CHS_sfFile     ="corrections/JER/Summer16_25nsV1_MC/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
    paramSet.jets.jer_AK8Puppi_resFile  ="corrections/JER/Summer16_25nsV1_MC/Summer16_25nsV1_MC_PtResolution_AK8PFPuppi.txt";
    paramSet.jets.jer_AK8Puppi_sfFile   ="corrections/JER/Summer16_25nsV1_MC/Summer16_25nsV1_MC_SF_AK8PFPuppi.txt";
    return paramSet;
}
ParameterSet set2017Parameters() {
    ParameterSet paramSet = setCommonParameters();
    paramSet.event.lumi = 41.5;
    paramSet.jets.jer_AK4CHS_resFile    ="corrections/JER/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK4PFchs.txt";
    paramSet.jets.jer_AK4CHS_sfFile     ="corrections/JER/Fall17_V3_MC/Fall17_V3_MC_SF_AK4PFchs.txt";
    paramSet.jets.jer_AK8Puppi_resFile  ="corrections/JER/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK8PFPuppi.txt";
    paramSet.jets.jer_AK8Puppi_sfFile   ="corrections/JER/Fall17_V3_MC/Fall17_V3_MC_SF_AK8PFPuppi.txt";
    paramSet.event.mcFilters.push_back(FillerConstants::FLAG_ecalBadCalibFilterUpdate);
    paramSet.event.dataFilters.push_back(FillerConstants::FLAG_ecalBadCalibFilterUpdate);
    return paramSet;
}
ParameterSet set2018Parameters() {
    ParameterSet paramSet = setCommonParameters();
    paramSet.event.lumi = 59.7;
    paramSet.jets.jer_AK4CHS_resFile    ="corrections/JER/Autumn18_V7_MC/Autumn18_V7_MC_PtResolution_AK4PFchs.txt";
    paramSet.jets.jer_AK4CHS_sfFile     ="corrections/JER/Autumn18_V7_MC/Autumn18_V7_MC_SF_AK4PFchs.txt";
    paramSet.jets.jer_AK8Puppi_resFile  ="corrections/JER/Autumn18_V7_MC/Autumn18_V7_MC_PtResolution_AK8PFPuppi.txt";
    paramSet.jets.jer_AK8Puppi_sfFile   ="corrections/JER/Autumn18_V7_MC/Autumn18_V7_MC_SF_AK8PFPuppi.txt";
    paramSet.event.mcFilters.push_back(FillerConstants::FLAG_ecalBadCalibFilterUpdate);
    paramSet.event.dataFilters.push_back(FillerConstants::FLAG_ecalBadCalibFilterUpdate);
    return paramSet;
}
}
}
