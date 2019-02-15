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

#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


#include "TSystem.h"
using namespace TAna;

// Macro Description:

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
    }



    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(diHiggsEvt.type <  DiHiggsEvent::MU) return false;

        auto wqq =  diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4();

        plotter.getOrMake1DPre(smpName,"deltaRlW",";#Delta#it{R}(lepton, W#rightarrowq#bar{q})",160,0,1.6)->
                Fill(std::abs(PhysicsUtilities::deltaR(wqq,diHiggsEvt.w1_d1->p4())),weight);

        plotter.getOrMake1DPre(smpName,"deltaPhiHH",";|#Delta#phi(H#rightarrowWW*, H#rightarrowb#bar{b})|",140,1.74,3.14)->
                Fill(std::fabs(PhysicsUtilities::deltaPhi(diHiggsEvt.hbb->p4(),diHiggsEvt.hww->p4())),weight);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void plotPaperGenPlots(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void plotPaperGenPlots(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
