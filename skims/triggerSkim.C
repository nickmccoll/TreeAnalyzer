
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "TreeReaders/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/EventSelection/interface/EventSelection.h"

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
        leptonProc .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonProc);
    }

    virtual BaseEventAnalyzer * setupEventAnalyzer() override {return new CopierEventAnalyzer();}

    virtual void bookOutputVariables() override {}

    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );
        reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData(),false);      load(reader_jetwlep     );
    }

    bool runEvent() override {
        auto jets = JetKinematics::selectObjects(reader_jetwlep->jets,30);
        const float ht_wlep = JetKinematics::ht(jets);
        std::vector<const Lepton    *> selectedLeptons;

        if(reader_electron && reader_muon){
            selectedLeptons = leptonProc->getLeptons(*reader_event,*reader_muon,*reader_electron);
        }
        const bool passTriggerPreselection= EventSelection::passTriggerPreselection(*reader_event,ht_wlep,selectedLeptons);
        if(!passTriggerPreselection) return false;
        return true;

    }

    std::shared_ptr<EventReader      > reader_event    ;
    std::shared_ptr<ElectronReader   > reader_electron ;
    std::shared_ptr<MuonReader       > reader_muon     ;
    std::shared_ptr<JetReader        > reader_jetwlep  ;
    std::unique_ptr<LeptonProcessor> leptonProc ;
};

#endif


void triggerSkim(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_ALL);
    a.analyze();
}
void triggerSkim(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_ALL);
    a.analyze();
}
