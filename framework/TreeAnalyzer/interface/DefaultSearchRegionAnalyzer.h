#ifndef TREEANALYZER_TREEANALYZER_DEFAULTSEARCHREGIONANALYZER_H
#define TREEANALYZER_TREEANALYZER_DEFAULTSEARCHREGIONANALYZER_H

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "DataFormats/interface/Momentum.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Configuration/interface/FillerConstants.h"
#include "Configuration/interface/ReaderConstants.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Variables/interface/DileptonSelection.h"

namespace TAna{
class EventReader       ;
class GenParticleReader ;
class ElectronReader    ;
class MuonReader        ;
class JetReader         ;
class FatJetReader      ;

class FatJetProcessor   ;
class TriggerScaleFactors;
class PUScaleFactors;
class LeptonScaleFactors;
class JetBTagScaleFactors;
class SubJetBTagScaleFactors;
class HbbFatJetScaleFactors;
class TopPTWeighting;
class JERCorrector;
class JESUncShifter;
class METUncShifter;
class HSolverChi;
class HSolverLi;

class Jet               ;
class FatJet            ;
class Lepton            ;

//--------------------------------------------------------------------------------------------------
// DefaultSearchRegionAnalyzer
//--------------------------------------------------------------------------------------------------
class DefaultSearchRegionAnalyzer : public BaseTreeAnalyzer {
public:
    //corrections that can be applied
    enum Corrections {CORR_XSEC
                     ,CORR_TRIG
                     ,CORR_PU
                     ,CORR_LEP
                     ,CORR_SJBTAG
                     ,CORR_AK4BTAG
                     ,CORR_SDMASS
                     ,CORR_TOPPT
                     ,CORR_JER
                     ,CORR_JES
                     ,CORR_MET
    };

    enum LepChannels {NOCHANNEL, SINGLELEP, DILEP};
    LepChannels lepChan = NOCHANNEL;

    DefaultSearchRegionAnalyzer(std::string fileName, std::string treeName,
            int treeInt, ASTypes::size randomSeed =0);

    virtual ~DefaultSearchRegionAnalyzer();

    //turn on, off correction
    void resetCorr() ;
    bool isCorrOn(Corrections corr) const;
    void turnOnCorr(Corrections corr);
    void turnOffCorr(Corrections corr);

    //default variable loading
    virtual void loadVariables() override;

    //check configuration
    virtual void checkConfig();

    //Set all the reader constants
    virtual void setupParameters();

    //fills class member variables....can be run before child runEvent
    virtual bool runEvent() override;

    bool isSignal() const {return mcProc == FillerConstants::SIGNAL;}

    std::shared_ptr<EventReader      > reader_event    ;
    std::shared_ptr<GenParticleReader> reader_genpart  ;
    std::shared_ptr<ElectronReader   > reader_electron ;
    std::shared_ptr<MuonReader       > reader_muon     ;
    std::shared_ptr<FatJetReader     > reader_fatjet   ;
    std::shared_ptr<FatJetReader     > reader_fatjet_noLep;
    std::shared_ptr<JetReader        > reader_jet        ;

    FillerConstants::MCProcess mcProc = FillerConstants::NOPROCESS;
    int             signal_mass=0;
    TString         smpName  = "";

    std::vector<const Jet*> jets_ht;
    float ht;

    std::vector<const Jet*> jets;
    std::vector<const Jet*> jets_HbbV;
    int nMedBTags = 0;
    int nMedBTags_HbbV = 0;

    ASTypes::size   corrections=0; //list of corrections to carry out
    float           weight     =0; //std weight after all corrections

    bool            passEventFilters = false;
    bool            passTriggerPreselection = false;
    bool            passTriggerPreselection2l = false;

    DiHiggsEvent    diHiggsEvt;
    SMDecayEvent    smDecayEvt;

    std::vector<const Lepton    *> selectedLeptons;
    const Lepton *                 selectedLepton=0;

    std::vector<const Lepton    *> selectedDileptons;
    const Lepton *                 dilep1=0;
    const Lepton *                 dilep2=0;

    const FatJet*              wjjCand     =0;
    const FatJet*              hbbCand     =0;
    BTagging::CSVSJ_CAT        hbbCSVCat   = BTagging::CSVSJ_INCL;
    BTagging::CSVSJ_CAT        wjjCSVCat   = BTagging::CSVSJ_INCL;

    float                      wwDM        = 0;
    float                      hwwChi      = 0;
    float                      hwwLi       = 0;

    MomentumF                  neutrino           ;
    MomentumF                  wlnu               ;
    MomentumF                  wqq                ;
    MomentumF                  hWW                ;
    MomentumF                  hh                 ;
    MomentumF                  hh_chi             ;
    MomentumF                  hh_basic           ;
    float                      hbbMass     =0     ;


    float llMass = 0;
    float llDR = 0;
    float llMetDphi = 0;

    ParameterSet parameters;
    FillerConstants::DataEra lastEra = FillerConstants::NOERA; //See if we need to load new params.


    std::unique_ptr<FatJetProcessor>        fjProc     ;
    std::unique_ptr<TriggerScaleFactors>    trigSFProc ;
    std::unique_ptr<PUScaleFactors>         puSFProc ;
    std::unique_ptr<LeptonScaleFactors>     leptonSFProc ;
    std::unique_ptr<JetBTagScaleFactors>    ak4btagSFProc ;
    std::unique_ptr<SubJetBTagScaleFactors> sjbtagSFProc ;
    std::unique_ptr<HbbFatJetScaleFactors>  hbbFJSFProc ;
    std::unique_ptr<TopPTWeighting>         topPTProc ;
    std::unique_ptr<JERCorrector>         JERProc ;
    std::unique_ptr<JESUncShifter>        JESUncProc ;
    std::unique_ptr<METUncShifter>          METUncProc;
    std::unique_ptr<HSolverChi>             hSolverChi;
    std::unique_ptr<HSolverLi>             hSolverLi;
};
}
#endif
