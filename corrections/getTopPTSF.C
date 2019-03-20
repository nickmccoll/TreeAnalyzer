
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Corrections/interface/EventWeights.h"


using namespace TAna;
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
        resetCorr();
        turnOnCorr(CORR_XSEC);
        turnOnCorr(CORR_PU);

    }


    void loadVariables() override  {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        checkConfig();
    }



    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(mcProc != FillerConstants::TTBAR) return false;

        plotter.getOrMake1D("nTops",";# of tops",5,-0.5,4.5)->Fill(smDecayEvt.topDecays.size(), weight);
        float truePT(0);
        if(smDecayEvt.topDecays.size() == 2){
            truePT = std::sqrt( smDecayEvt.topDecays[0].top->pt() * smDecayEvt.topDecays[1].top->pt());
        }
        plotter.getOrMake1D("true_pt",";<top pt>",100,0,2000)->Fill(truePT, weight);

        const float avgPT= topPTProc->getAvgPT(smDecayEvt);
        const float corrNoNorm = topPTProc->getCorrectionNoNorm(mcProc,smDecayEvt);
        const float corr = topPTProc->getCorrection(mcProc,smDecayEvt);

        plotter.getOrMake1D("van_pt",";<top pt>",100,0,2000)->Fill(avgPT, weight);
        plotter.getOrMake1D("cor_pt",";<top pt>",100,0,2000)->Fill(avgPT, corrNoNorm*weight);
        plotter.getOrMake1D("corNorm_pt",";<top pt>",100,0,2000)->Fill(avgPT, corr*weight);

        plotter.getOrMake1D("counts","; van / cor / corNorm",3,-0.5,2.5)->Fill(0., weight);
        plotter.getOrMake1D("counts","; van / cor / corNorm",3,-0.5,2.5)->Fill(1., corrNoNorm*weight);
        plotter.getOrMake1D("counts","; van / cor / corNorm",3,-0.5,2.5)->Fill(2., corr*weight);

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getTopPTSF(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getTopPTSF(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
