

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
    paramSet.event.minHT         =400;
    paramSet.event.minTriggerEl  =30;
    paramSet.event.minTriggerMu  =27;
    paramSet.event.leptonCorrSFFile = "";
    paramSet.event.puCorrSFFile = "";
    paramSet.event.ttbarXSecSF_1000toInf_nLep0 = -1;
    paramSet.event.ttbarXSecSF_1000toInf_nLep1 = -1;
    paramSet.event.ttbarXSecSF_1000toInf_nLep2 = -1;
    paramSet.event.ttbarXSecSF_700to1000_nLep0 = -1;
    paramSet.event.ttbarXSecSF_700to1000_nLep1 = -1;
    paramSet.event.ttbarXSecSF_700to1000_nLep2 = -1;
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
    paramSet.dileptons.el_maxSip3D   = 99999   ;
    paramSet.dileptons.el_maxISO  = 0.15 ;
    paramSet.dileptons.el_getID1   = &Electron::passMedID_noIso;
    paramSet.dileptons.el_getID2   = &Electron::passMedID_noIso;
    paramSet.dileptons.el_getISO  = &Electron::miniIso;

    paramSet.dileptons.mu_minPT   = 10  ;
    paramSet.dileptons.mu_maxETA  = 2.4 ;
    paramSet.dileptons.mu_maxDZ   = 0.1 ;
    paramSet.dileptons.mu_maxD0   = 0.05;
    paramSet.dileptons.mu_maxSip3D   = 99999   ;
    paramSet.dileptons.mu_maxISO  = 0.15 ;
    paramSet.dileptons.mu_getID1   = &Muon::passLooseID;
    paramSet.dileptons.mu_getID2   = &Muon::passLooseID;
    paramSet.dileptons.mu_getISO  = &Muon::miniIso;

    paramSet.jets.minJetPT      = 20;
    paramSet.jets.maxJetETA     = 2.4;
    paramSet.jets.passJetID     = &Jet::passTightID;

    paramSet.jets.CSV_WP        = {-100,100,100,100};
    paramSet.jets.DeepCSV_WP    = {-100,100,100,100};
    paramSet.jets.DeepFlavor_WP = {-100,100,100,100};
    paramSet.jets.minBtagJetPT  = 30;
    paramSet.jets.maxBTagJetETA = 2.4;

    paramSet.jets.getJetBTagVal = &BaseRecoJet::deep_flavor;
    paramSet.jets.jetBTagWP     = paramSet.jets.DeepFlavor_WP [BTagging::BTAG_M];
    paramSet.jets.getSJBTagVal = &BaseRecoJet::deep_csv;
    paramSet.jets.sjBTagLWP = paramSet.jets.DeepCSV_WP [BTagging::BTAG_L];
    paramSet.jets.sjBTagMWP = paramSet.jets.DeepCSV_WP [BTagging::BTAG_M];

    paramSet.jets.jetBtagCorrSFFile ="corrections/CSVv2_Moriond17_B_H.csv";
    paramSet.jets.jetBtagCorrEffFile ="corrections/ak4_csvEff.root";
    paramSet.jets.jetBtagCorrWP       = paramSet.jets.DeepFlavor_WP      ;
    paramSet.jets.jetBtagCorrGetBTagVal = paramSet.jets.getJetBTagVal;
    paramSet.jets.sjBtagCorrWP         = paramSet.jets.DeepCSV_WP      ;
    paramSet.jets.sjBtagCorrGetBTagVal = paramSet.jets.getSJBTagVal;
    paramSet.jets.sjBtagCorrSFFile     = "corrections/subjet_CSVv2_Moriond17_B_H.csv";
    paramSet.jets.sjBtagCorrEffFile    = "corrections/sj_csvEff.root";

    paramSet.hww.dilepInvMassGuess = 55;

    return paramSet;
}

void setJetWorkingPoints(JetParameters& jetParams, std::vector<float> wps) {
    jetParams.jetBTagWP = wps [BTagging::BTAG_M];
    jetParams.jetBtagCorrWP = wps;
}

void setSubjetWorkingPoints(JetParameters& jetParams, std::vector<float> wps) {
    jetParams.sjBTagLWP = wps [BTagging::BTAG_L];
    jetParams.sjBTagMWP = wps [BTagging::BTAG_M];
    jetParams.sjBtagCorrWP  = wps;
}

ParameterSet set2016Parameters() {
    ParameterSet paramSet = setCommonParameters();
    paramSet.event.lumi = 35.9;
    paramSet.event.puCorrSFFile = "corrections/pileup/puSF_2016.root";
    paramSet.event.leptonCorrSFFile = "corrections/trigger/triggerSF_2016.root";
    paramSet.event.passTrigger = &EventSelection::passTriggerSuite2016;
    paramSet.event.ttbarXSecSF_700to1000_nLep0 = 0.00159334;
    paramSet.event.ttbarXSecSF_700to1000_nLep1 = 0.00135322;
    paramSet.event.ttbarXSecSF_700to1000_nLep2 = 0.00078532;
    paramSet.event.ttbarXSecSF_1000toInf_nLep0 = 0.000707199;
    paramSet.event.ttbarXSecSF_1000toInf_nLep1 = 0.000655208;
    paramSet.event.ttbarXSecSF_1000toInf_nLep2 = 0.000478046;

    paramSet.leptons.el_SFFile          ="corrections/lepton/electron_1l_2016_SF.root";
    paramSet.leptons.mu_SFFile          ="corrections/lepton/muon_1l_2016_SF.root";
    paramSet.dileptons.el_SFFile        ="corrections/lepton/electron_2l_2016_SF.root";
    paramSet.dileptons.mu_SFFile        ="corrections/lepton/muon_2l_2016_SF.root";

    paramSet.jets.jer_AK4CHS_resFile    ="corrections/JER/Summer16_25nsV1_MC/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
    paramSet.jets.jer_AK4CHS_sfFile     ="corrections/JER/Summer16_25nsV1_MC/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
    paramSet.jets.jer_AK8Puppi_resFile  ="corrections/JER/Summer16_25nsV1_MC/Summer16_25nsV1_MC_PtResolution_AK8PFPuppi.txt";
    paramSet.jets.jer_AK8Puppi_sfFile   ="corrections/JER/Summer16_25nsV1_MC/Summer16_25nsV1_MC_SF_AK8PFPuppi.txt";

    paramSet.jets.DeepCSV_WP    = {-100,0.2217,0.6321,0.8953}; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
    paramSet.jets.DeepFlavor_WP = {-100,0.0614,0.3093,0.7221}; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
    setJetWorkingPoints(paramSet.jets,paramSet.jets.DeepFlavor_WP);
    setSubjetWorkingPoints(paramSet.jets,paramSet.jets.DeepCSV_WP);

    paramSet.hww.ptCorB        = 0.989509;
    paramSet.hww.ptCorM        = -1.50788e-05;
    paramSet.hww.liFileName    ="variables/hhSol_templates_2016.root";

    return paramSet;
}
ParameterSet set2017Parameters() {
    ParameterSet paramSet = setCommonParameters();
    paramSet.event.lumi = 41.5;
    paramSet.event.puCorrSFFile = "corrections/pileup/puSF_2017.root";
    paramSet.event.leptonCorrSFFile = "corrections/trigger/triggerSF_2017.root";
    paramSet.event.passTrigger = &EventSelection::passTriggerSuite2017;
    paramSet.event.ttbarXSecSF_700to1000_nLep0 = 0.00122835;
    paramSet.event.ttbarXSecSF_700to1000_nLep1 = 0.00129256;
    paramSet.event.ttbarXSecSF_700to1000_nLep2 = 0.00076952;
    paramSet.event.ttbarXSecSF_1000toInf_nLep0 = 0.000719423;
    paramSet.event.ttbarXSecSF_1000toInf_nLep1 = 0.000741457;
    paramSet.event.ttbarXSecSF_1000toInf_nLep2 = 0.00052646;

    paramSet.leptons.el_SFFile          ="corrections/lepton/electron_1l_2017_SF.root";
    paramSet.leptons.mu_SFFile          ="corrections/lepton/muon_1l_2017_SF.root";
    paramSet.dileptons.el_SFFile        ="corrections/lepton/electron_2l_2017_SF.root";
    paramSet.dileptons.mu_SFFile        ="corrections/lepton/muon_2l_2017_SF.root";

    paramSet.jets.jer_AK4CHS_resFile    ="corrections/JER/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK4PFchs.txt";
    paramSet.jets.jer_AK4CHS_sfFile     ="corrections/JER/Fall17_V3_MC/Fall17_V3_MC_SF_AK4PFchs.txt";
    paramSet.jets.jer_AK8Puppi_resFile  ="corrections/JER/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK8PFPuppi.txt";
    paramSet.jets.jer_AK8Puppi_sfFile   ="corrections/JER/Fall17_V3_MC/Fall17_V3_MC_SF_AK8PFPuppi.txt";
    paramSet.event.mcFilters.push_back(FillerConstants::FLAG_ecalBadCalibFilterUpdate);
    paramSet.event.dataFilters.push_back(FillerConstants::FLAG_ecalBadCalibFilterUpdate);

    paramSet.jets.CSV_WP        = {-100,0.5803,0.8838,0.9693}; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    paramSet.jets.DeepCSV_WP    = {-100,0.1522,0.4941,0.8001}; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    paramSet.jets.DeepFlavor_WP = {-100,0.0521,0.3033,0.7489}; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    setJetWorkingPoints(paramSet.jets,paramSet.jets.DeepFlavor_WP);
    setSubjetWorkingPoints(paramSet.jets,paramSet.jets.DeepCSV_WP);

    paramSet.hww.ptCorB        = 0.966024;
    paramSet.hww.ptCorM        = -8.90763e-06;
    paramSet.hww.liFileName    ="variables/hhSol_templates_2017.root";

    return paramSet;
}
ParameterSet set2018Parameters() {
    ParameterSet paramSet = setCommonParameters();
    paramSet.event.lumi = 59.7;
    paramSet.event.puCorrSFFile = "corrections/pileup/puSF_2018.root";
    paramSet.event.leptonCorrSFFile = "corrections/trigger/triggerSF_2018.root";
    paramSet.event.passTrigger = &EventSelection::passTriggerSuite2018;
    paramSet.event.ttbarXSecSF_700to1000_nLep0 = 0.000748296;
    paramSet.event.ttbarXSecSF_700to1000_nLep1 = 0.000776814;
    paramSet.event.ttbarXSecSF_700to1000_nLep2 = 0.000800626;
    paramSet.event.ttbarXSecSF_1000toInf_nLep0 = 0.00050157;
    paramSet.event.ttbarXSecSF_1000toInf_nLep1 = 0.000514053;
    paramSet.event.ttbarXSecSF_1000toInf_nLep2 = 0.000517094;

    paramSet.leptons.el_SFFile          ="corrections/lepton/electron_1l_2018_SF.root";
    paramSet.leptons.mu_SFFile          ="corrections/lepton/muon_1l_2018_SF.root";
    paramSet.dileptons.el_SFFile        ="corrections/lepton/electron_2l_2018_SF.root";
    paramSet.dileptons.mu_SFFile        ="corrections/lepton/muon_2l_2018_SF.root";

    paramSet.jets.jer_AK4CHS_resFile    ="corrections/JER/Autumn18_V7_MC/Autumn18_V7_MC_PtResolution_AK4PFchs.txt";
    paramSet.jets.jer_AK4CHS_sfFile     ="corrections/JER/Autumn18_V7_MC/Autumn18_V7_MC_SF_AK4PFchs.txt";
    paramSet.jets.jer_AK8Puppi_resFile  ="corrections/JER/Autumn18_V7_MC/Autumn18_V7_MC_PtResolution_AK8PFPuppi.txt";
    paramSet.jets.jer_AK8Puppi_sfFile   ="corrections/JER/Autumn18_V7_MC/Autumn18_V7_MC_SF_AK8PFPuppi.txt";
    paramSet.event.mcFilters.push_back(FillerConstants::FLAG_ecalBadCalibFilterUpdate);
    paramSet.event.dataFilters.push_back(FillerConstants::FLAG_ecalBadCalibFilterUpdate);

    paramSet.jets.DeepCSV_WP    = {-100,0.1241,0.4184,0.7527}; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    paramSet.jets.DeepFlavor_WP = {-100,0.0494,0.2770,0.7264}; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    setJetWorkingPoints(paramSet.jets,paramSet.jets.DeepFlavor_WP);
    setSubjetWorkingPoints(paramSet.jets,paramSet.jets.DeepCSV_WP);

    paramSet.hww.ptCorB        = 0.965834;
    paramSet.hww.ptCorM        = -8.75461e-06;
    paramSet.hww.liFileName    ="variables/hhSol_templates_2018.root";

    return paramSet;
}
}
}
