#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DileptonSearchRegionAnalyzer.h"
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
#include "Processors/Variables/interface/JetKinematics.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/LeptonSelection.h"

#include "TSystem.h"
using namespace TAna;

class Skimmer : public DileptonSearchRegionAnalyzer {
public:

    Skimmer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DileptonSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){}

    virtual BaseEventAnalyzer * setupEventAnalyzer() override {return new CopierEventAnalyzer();}
    virtual void bookOutputVariables() override {}

    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event);
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiJet",isRealData());       load(reader_fatjet);
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiJet",isRealData(),false); load(reader_jet);
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon);

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");          load(reader_genpart   );
        }
        checkConfig();
    }

    bool runEvent() override {
        if(!DileptonSearchRegionAnalyzer::runEvent()) return false;

        // require Hbb cand, which itself requires two selected Dileptons to be fille in DileptonSearchRegionAnalyzer
        if(!hbbCand) return false;
        reader_event->genWeights->clear();

        return true;
    }
};

#endif

void skimDileptonSel(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Skimmer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_LOADED);
    a.analyze();
}
void skimDileptonSel(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Skimmer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_LOADED);
    a.analyze();
}
