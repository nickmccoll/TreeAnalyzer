#ifndef TREEANALYZER_TREEANALYZER_DEFAULTSEARCHREGIONANALYZER_H
#define TREEANALYZER_TREEANALYZER_DEFAULTSEARCHREGIONANALYZER_H

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "DataFormats/interface/Momentum.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"



namespace TAna{
class EventReader       ;
class GenParticleReader ;
class ElectronReader    ;
class MuonReader        ;
class JetReader         ;
class JetReader         ;
class FatJetReader      ;

class FatJetProcessor   ;
class LeptonProcessor   ;

class FatJet            ;
class Lepton            ;

class DefaultSearchRegionAnalyzer : public BaseTreeAnalyzer {
public:

    DefaultSearchRegionAnalyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
        setupProcessors(fileName);
    }

    //all of the processor/variable constructing
    virtual void setupProcessors(std::string fileName);

    //default variable loading
    virtual void loadVariables() override;

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

    float           ht_wlep    =0;
    float           weight     =0;
    bool            passEventFilters = 0;
    DiHiggsEvent    diHiggsEvt;

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

    std::unique_ptr<FatJetProcessor> fjProc     ;
    std::unique_ptr<LeptonProcessor> leptonProc ;
};
}
#endif
