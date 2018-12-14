#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "TreeReaders/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        if(diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        TString sn;
        const GenParticle* lep1 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1;
        const GenParticle* lep2 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1;
        int lep1id = lep1->pdgId();
        int lep2id = lep2->pdgId();

        // throw away dilepton events with taus and same-sign leptons, then record the dilep channel
        if (lep1id<0 == lep2id<0) return false;

        if (abs(lep1id) == 15 || abs(lep2id) == 15) return false;
        else if (abs(lep1id) == 11 && abs(lep2id) == 11) sn += "_ee";
        else if (abs(lep1id) == 13 && abs(lep2id) == 13) sn += "_mumu";
        else if ((abs(lep1id)==13 && abs(lep2id)==11) || ((abs(lep1id)==11 && abs(lep2id)==13))) sn += "_emu";

        // lepton pt cuts
        bool passPt1 = lep1->absPdgId() == 13 ? lep1->pt() > 26 : lep1->pt() > 30;
        bool passPt2 = lep2->pt() > 10;
        if (!(passPt1 && passPt2)) return false;



        return true;
    }
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void compareMonoDiLepStats(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void compareMonoDiLepStats(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
