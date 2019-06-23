#ifndef TREEANALYZER_TREEANALYZER_DEFAULTSEARCHREGIONANALYZER_H
#define TREEANALYZER_TREEANALYZER_DEFAULTSEARCHREGIONANALYZER_H

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "DataFormats/interface/Momentum.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Configuration/interface/FillerConstants.h"
#include "Configuration/interface/ReaderConstants.h"
#include "Processors/Variables/interface/BTagging.h"


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
class HiggsSolver;

class Jet               ;
class FatJet            ;
class Lepton            ;

//--------------------------------------------------------------------------------------------------
// DefaultSearchRegionAnalyzer
//--------------------------------------------------------------------------------------------------
class DefaultSearchRegionAnalyzer : public BaseTreeAnalyzer {
public:
    //corrections that can be applied
    enum Corrections {CORR_XSEC    =(1<<0)
                     ,CORR_TRIG    =(1<<1)
                     ,CORR_PU      =(1<<2)
                     ,CORR_LEP     =(1<<3)
                     ,CORR_SJBTAG  =(1<<4)
                     ,CORR_AK4BTAG =(1<<5)
                     ,CORR_SDMASS  =(1<<6)
                     ,CORR_TOPPT   =(1<<7)
                     ,CORR_JER     =(1<<8)
                     ,CORR_JES     =(1<<9)
                     ,CORR_MET     =(1<<10)
    };

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
    std::shared_ptr<JetReader        > reader_jet_chs    ;

    FillerConstants::MCProcess mcProc = FillerConstants::NOPROCESS;
    int             signal_mass=0;
    TString         smpName  = "";

    std::vector<const Jet*> jets;
    std::vector<const Jet*> jets_HbbV;
    std::vector<const Jet*> jets_HbbV_2l;
    int nMedBTags = 0;
    int nMedBTags_HbbV = 0;
    int nMedBTags_HbbV_2l = 0;
    float ht = 0;

    std::vector<const Jet*> jets_chs_30;
    float                   ht_chs_30    =0;
    std::vector<const Jet*> jets_chs_20;
    float                   ht_chs_20    =0;

    std::vector<const Jet*> jets_puppi_30;
    float                   ht_puppi_30    =0;
    std::vector<const Jet*> jets_puppi_20;
    float                   ht_puppi_20    =0;

    std::vector<const Jet*> jets_NoDilepOverlap;
    float                   ht_NoDilepOverlap = 0;

    ASTypes::size   corrections=0; //list of corrections to carry out
    float           weight     =0; //std weight after all corrections

    bool            passEventFilters = false;
    bool            passTriggerPreselection = false;

    DiHiggsEvent    diHiggsEvt;
    SMDecayEvent    smDecayEvt;
    float           genMtt = 0;
	int nLepsTT = -1;

    std::vector<const Lepton    *> selectedLeptons;
    const Lepton *                 selectedLepton=0;

    std::vector<const Lepton    *> selectedDileptons;
    enum DilepChan {ee, mumu, emu};
    DilepChan dilepChan;
    std::map<DilepChan,TString> dilepMap = { {ee, "_ee_"}, {mumu,"_mumu_"}, {emu,"_emu_"}};

    const FatJet*              wjjCand     =0;
    const FatJet*              hbbCand     =0;
    BTagging::CSVSJ_CAT        hbbCSVCat   = BTagging::CSVSJ_INCL;
    BTagging::CSVSJ_CAT        wjjCSVCat   = BTagging::CSVSJ_INCL;

    const FatJet*       hbbCand_2l       =0;
    BTagging::CSVSJ_CAT hbbCSVCat_2l     = BTagging::CSVSJ_INCL;

    float                      wwDM        = 0;
    float                      hwwChi      = 0;
    bool                       passWWDM    = false;


    MomentumF                  neutrino           ;
    MomentumF                  wlnu               ;
    MomentumF                  wqq                ;
    MomentumF                  hWW                ;
    MomentumF                  hh                 ;
    MomentumF                  hh_old             ;
    float                      hbbMass     =0     ;

    MomentumF                  hWW_2l             ;
    MomentumF                  hh_2l              ;
    float                      hbbMass_2l  =0     ;

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
    std::unique_ptr<JERCorrector>         JERAK4PuppiProc ;
    std::unique_ptr<JERCorrector>         JERAK4CHSProc ;
    std::unique_ptr<JERCorrector>         JERAK8PuppiProc ;
    std::unique_ptr<JESUncShifter>        JESUncProc ;
    std::unique_ptr<METUncShifter>          METUncProc;
    std::unique_ptr<HiggsSolver>            higgsSolver;
};
}
#endif
