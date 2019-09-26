#ifndef CONFIGURATION_INTERFACE_READERCONSTANTS_H
#define CONFIGURATION_INTERFACE_READERCONSTANTS_H

#include "AnalysisSupport/Utilities/interface/Types.h"
#include "Configuration/interface/FillerConstants.h"
#include <vector>
#include <string>

namespace TAna {

//--------------------------------------------------------------------------------------------------
class EventReader;
namespace EventSelection {
typedef bool (*TriggerSel)(const EventReader& reader_event);
}
struct EventParameters{
    float lumi = -1;
    bool  doTTBarStitching = true;

    std::vector<FillerConstants::METFilters> dataFilters;
    std::vector<FillerConstants::METFilters> mcFilters;
    EventSelection::TriggerSel  passTrigger =0 ;

    float minHT;
    float minTriggerEl;
    float minTriggerMu;
    float ttbarXSecSF_1000toInf_nLep0;
    float ttbarXSecSF_1000toInf_nLep1;
    float ttbarXSecSF_1000toInf_nLep2;
    float ttbarXSecSF_700to1000_nLep0;
    float ttbarXSecSF_700to1000_nLep1;
    float ttbarXSecSF_700to1000_nLep2;

    std::string leptonCorrSFFile;
    std::string puCorrSFFile;
};

//--------------------------------------------------------------------------------------------------
class FatJet;
namespace FatJetSelHelpers {
typedef  bool (FatJet::*fjFunBool)() const;
}
struct FatJetParameters{
    //parameters to select jets on the loadFatJets func
    float cand_minPT     =-1;
    float cand_maxETA    =-1;
    FatJetSelHelpers::fjFunBool fjJetID =0 ;

    //parameters applied to both wjj and hbb AND dilep hbb
    float sj_minPT       =-1; //cand sel
    float sj_maxETA      =-1; //cand sel

    //parameters applied to wjj candidate sel
    float wjj_maxLepDR   =-1;
    float wjj_minPT      =-1;
    float wjj_minSJs     =-1;

    //parameters applied to hbb candidate sel
    float hbb_minLepDPhi =-1;
    float hbb_minPT      =-1;
    float hbb_minSJs     =-1;

    //parameters applied to dilep Hbb candidate sel
    float hbbLL_minDphiBBLL = -1;
    float hbbLL_minPT       = -1;
    float hbbLL_minSJs      = -1;
    float hbbLL_minDRbbLL   = -1;
};

//--------------------------------------------------------------------------------------------------
class Muon;
class Electron;
namespace LeptonProcessor {
    typedef  bool (Muon::*muFunBool)() const;
    typedef  bool (Electron::*elFunBool)() const;
    typedef  float (Muon::*muFunFloat)() const;
    typedef  float (Electron::*elFunFloat)() const;
}
struct LeptonParameters {
    float el_minPT  ;
    float el_maxETA ;
    float el_maxDZ  ;
    float el_maxD0  ;
    float el_maxSip3D  ;
    float el_maxISO ;
    LeptonProcessor::elFunBool el_getID ;
    LeptonProcessor::elFunFloat el_getISO;

    float mu_minPT  ;
    float mu_maxETA ;
    float mu_maxDZ  ;
    float mu_maxD0  ;
    float mu_maxSip3D  ;
    float mu_maxISO ;
    LeptonProcessor::muFunBool mu_getID ;
    LeptonProcessor::muFunFloat mu_getISO;
};

//--------------------------------------------------------------------------------------------------
namespace DileptonProcessor {
    typedef  bool (Muon::*muFunBool)() const;
    typedef  bool (Electron::*elFunBool)() const;
    typedef  float (Muon::*muFunFloat)() const;
    typedef  float (Electron::*elFunFloat)() const;
}
struct DileptonParameters {
    float el_minPT  ;
    float el_maxETA ;
    float el_maxDZ  ;
    float el_maxD0  ;
    float el_maxSip3D  ;
    float el_maxISO ;
    DileptonProcessor::elFunBool el_getID1 ;
    DileptonProcessor::elFunBool el_getID2 ;
    DileptonProcessor::elFunFloat el_getISO;

    float mu_minPT  ;
    float mu_maxETA ;
    float mu_maxDZ  ;
    float mu_maxD0  ;
    float mu_maxSip3D  ;
    float mu_maxISO ;
    DileptonProcessor::muFunBool mu_getID1 ;
    DileptonProcessor::muFunBool mu_getID2 ;
    DileptonProcessor::muFunFloat mu_getISO;

};
//--------------------------------------------------------------------------------------------------

class BaseRecoJet;
class Jet;
namespace BTagging {
    typedef  float (BaseRecoJet::*getBTagVal)() const;
}
namespace JetProcessor{
typedef  bool (Jet::*jetFun)() const;
}
struct JetParameters {


    float minJetPT      ;
    float maxJetETA     ;
    JetProcessor::jetFun    passJetID;

    std::vector<float> CSV_WP;
    std::vector<float> DeepCSV_WP;
    float minBtagJetPT  ;
    float maxBTagJetETA ;

    BTagging::getBTagVal getJetBTagVal;
    float jetBTagWP;
    BTagging::getBTagVal getSJBTagVal;
    float sjBTagLWP;
    float sjBTagMWP;

    std::vector<float>   jetBtagCorrWP;
    BTagging::getBTagVal jetBtagCorrGetBTagVal;
    std::string          jetBtagCorrSFFile;
    std::string          jetBtagCorrEffFile;
    std::vector<float>   sjBtagCorrWP;
    BTagging::getBTagVal sjBtagCorrGetBTagVal;
    std::string          sjBtagCorrSFFile;
    std::string          sjBtagCorrEffFile;


};

struct HWWParameters {
    // CHI parameters
    double posMETParErr     ;
    double negMETParErr     ;
    double metPerpErr       ;
    double jetErr           ;
    double onWlnuMeanJet    ;
    double offWlnuMeanJet   ;
    double onWlnuMeanWlnu   ;
    double offWlnuMeanWlnu  ;
    double offWlnuPosWlnuErr;
    double offWnluNegWlnuErr;
    double onWlnuWlnuErr    ;
    double onWlnuHWWErr     ;
    double offWlnuHWWErr    ;

    //LI parameters
    double ptCorB           ;
    double ptCorM           ;
    std::string   liFileName;
    std::string   bkgLiFileName;

    double dilepInvMassGuess;

};

//--------------------------------------------------------------------------------------------------
struct ParameterSet {
    EventParameters  event;
    FatJetParameters fatJets;
    LeptonParameters leptons;
    DileptonParameters dileptons;
    JetParameters    jets;
    HWWParameters    hww;
};
//--------------------------------------------------------------------------------------------------
namespace ReaderConstants{
ParameterSet setCommonParameters();
ParameterSet set2016Parameters();
ParameterSet set2017Parameters();
ParameterSet set2018Parameters();
}
}


#endif
