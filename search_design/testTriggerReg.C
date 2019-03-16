
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
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(selectedLeptons.size() != 1) return false;
        TString prefix = smpName;
        if(passHbbTSel){
            prefix += "_passTight";
        } else if(passHbbSel){
            prefix += "_passLoose";
        } else {
            prefix += "_fail";
        }

        plotter.getOrMake1DPre(prefix, "lepPT"  ,";lepton p_{T} [GeV]; arbitrary units",500,0,500)->Fill( selectedLepton->pt(),weight);
        plotter.getOrMake1DPre(prefix, "ht"  ,";H_{T} [GeV]; arbitrary units",2000,0,2000)->Fill( ht_wlep,weight);

        double rebinHT[] = {0,200,300,400,450,500,550,600,800,1000,1200,1400,1800,2000,4000,8000,10000};
        int nRHT = 16;

        double rebinLeppt[] = {20,25,30,35,40,50,75,100,125,200,300,500,1000};
        int nRL = 12;

        plotter.getOrMake2DPre(prefix, "htvlepPT"  ,";lepton p_{T} [GeV]; H_{T} [GeV]; arbitrary units",nRL,rebinLeppt,nRHT,rebinHT)->Fill( selectedLepton->pt(),ht_wlep,weight);
        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void testTriggerReg(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void testTriggerReg(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
