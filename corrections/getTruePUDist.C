
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
#include "TreeReaders/interface/EventReader.h"
using namespace TAna;
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){


        ASTypes::size  newCorr=0;
        FillerConstants::addPass(newCorr,CORR_XSEC);
        corrections = newCorr;
    }


    void loadVariables() override  {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        plotter.getOrMake1DPre(smpName,"numTruePUInt",";numTruePUInt",75,0,75)->Fill(reader_event->nTruePUInts,weight);
        if(reader_event->process == FillerConstants::SIGNAL ) plotter.getOrMake1DPre("signal","numTruePUInt",";numTruePUInt",75,0,75)->Fill(reader_event->nTruePUInts,weight);
        plotter.getOrMake1DPre("all","numTruePUInt",";numTruePUInt",75,0,75)->Fill(reader_event->nTruePUInts,weight);
        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getTruePUDist(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getTruePUDist(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
