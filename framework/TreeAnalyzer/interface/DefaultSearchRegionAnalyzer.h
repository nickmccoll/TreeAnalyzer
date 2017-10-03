#ifndef TREEANALYZER_TREEANALYZER_DEFAULTSEARCHREGIONANALYZER_H
#define TREEANALYZER_TREEANALYZER_DEFAULTSEARCHREGIONANALYZER_H

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "DataFormats/interface/Momentum.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"



namespace TAna{
class EventReader       ;
class GenParticleReader ;
class ElectronReader    ;
class MuonReader        ;
class JetReader         ;
class FatJetReader      ;

class FatJetProcessor   ;
class LeptonProcessor   ;
class TriggerScaleFactors;
class PUScaleFactors;
class LeptonScaleFactors;
class JetBTagScaleFactors;
class SubJetBTagScaleFactors;

class Jet               ;
class FatJet            ;
class Lepton            ;

class DefaultSearchRegionAnalyzer : public BaseTreeAnalyzer {
public:
    //corrections that can be applied
    enum Corrections {CORR_XSEC =(1<<0)
                     ,CORR_TRIG =(1<<1)
                     ,CORR_PU   =(1<<2)
                     ,CORR_LEP  =(1<<3)
                     ,CORR_SJBTAG  =(1<<4)
                     ,CORR_AK4BTAG  =(1<<5)
    };

    DefaultSearchRegionAnalyzer(std::string fileName, std::string treeName, int treeInt);

    virtual ~DefaultSearchRegionAnalyzer();

    //all of the processor/variable constructing
    virtual void setupProcessors(std::string fileName);

    //turn on, off correction
    void resetCorr() ;
    bool isCorrOn(Corrections corr) const;
    void turnOnCorr(Corrections corr);
    void turnOffCorr(Corrections corr);

    //default variable loading
    virtual void loadVariables() override;

    //check configuration
    virtual void checkConfig();

    //fills class member variables....can be run before child runEvent
    virtual bool runEvent() override;

    std::shared_ptr<EventReader      > reader_event    ;
    std::shared_ptr<GenParticleReader> reader_genpart  ;
    std::shared_ptr<ElectronReader   > reader_electron ;
    std::shared_ptr<MuonReader       > reader_muon     ;
    std::shared_ptr<FatJetReader     > reader_fatjet   ;
    std::shared_ptr<JetReader        > reader_jet      ;
    std::shared_ptr<JetReader        > reader_jetwlep  ;


    int             signal_mass=0;
    TString         smpName  = "";

    std::vector<const Jet*> jets;

    std::vector<const Jet*> jets_wlep;
    float                   ht_wlep    =0;

    ASTypes::size   corrections=0; //list of corrections to carry out
    float           weight     =0; //std weight after all corrections

    bool            passEventFilters = false;
    bool            passTriggerPreselection = false;

    DiHiggsEvent    diHiggsEvt;
    SMDecayEvent    smDecayEvt;

    std::vector<const Lepton    *> selectedLeptons;
    const Lepton *                 selectedLepton=0;


    std::vector<const FatJet*> fatjetCands;
    const FatJet*              wjjCand     =0;
    const FatJet*              hbbCand     =0;
    bool                       passWjjSel  = false;
    bool                       passHbbSel  = false;
    bool                       passHbbTSel = false;
    MomentumF                  neutrino           ;
    MomentumF                  hh                 ;

    std::unique_ptr<FatJetProcessor>     fjProc     ;
    std::unique_ptr<LeptonProcessor>     leptonProc ;
    std::unique_ptr<TriggerScaleFactors> trigSFProc ;
    std::unique_ptr<PUScaleFactors>      puSFProc ;
    std::unique_ptr<LeptonScaleFactors>  leptonSFProc ;
    std::unique_ptr<JetBTagScaleFactors>    ak4btagSFProc ;
    std::unique_ptr<SubJetBTagScaleFactors>    sjbtagSFProc ;
};
}
#endif
