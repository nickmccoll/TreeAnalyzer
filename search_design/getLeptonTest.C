
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FatJetReader.h"

using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {
    }

//    bool runEvent() override {
//        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
//        if(selectedLeptons.size() != 1) return false;
//        if(diHiggsEvt.type < DiHiggsEvent::MU) return false;
//
//        const float minHT = 400;
//        const float minElePT = 30;
//        const float minMuPT = 26;
//
//        if(ht_chs < minHT) return false;
//
//        float maxElePT = 0;
//        float maxMuPT = 0;
//        for(const auto * l : selectedLeptons ) {
//            if(l->isMuon()) maxMuPT = std::max(maxMuPT, l->pt());
//            else maxElePT = std::max(maxElePT, l->pt());
//        }
//        if(maxElePT < minElePT && maxMuPT < minMuPT ) return false;
//
//
//        std::string prefix = std::string(diHiggsEvt.type == DiHiggsEvent::MU ? "mu" : "e") + "_"+ smpName.Data();
//
//        plotter.getOrMake1DPre(prefix.c_str(),"incl",";HT",36,400,4000)->Fill(ht_chs,weight);
//
//        if(!EventSelection::passTriggerSuite(*reader_event)) return false;
//
//        plotter.getOrMake1DPre(prefix.c_str(),"pass",";HT",36,400,4000)->Fill(ht_chs,weight);
//
//        return true;
//    }

    bool runEvent() override {
        bool passPre = true;
        if(!DefaultSearchRegionAnalyzer::runEvent()) passPre = false;
        if(diHiggsEvt.type < DiHiggsEvent::MU) return false;
        if(!passTriggerPreselection) passPre = false;
//        if(!passEventFilters) passPre = false;
        if(selectedLeptons.size() != 1) passPre = false;
        if(!hbbCand)  passPre = false;
        if(!wjjCand)  passPre = false;
        if(!passPre) return false;
        if(wwDM>=125)  passPre = false;
        if(hWW.pt()/hh.mass() <= 0.3)  passPre = false;
        if(nMedBTags_HbbV >0) passPre=false;
        if(wjjCand->tau2otau1() >=0.75) passPre = false;
        if(hbbCSVCat < 4) passPre = false;
        if(!passPre) return false;

        double drHbb = PhysicsUtilities::deltaR(*diHiggsEvt.hbb,*hbbCand);
        double drbb = std::max(PhysicsUtilities::deltaR(*diHiggsEvt.b1,*hbbCand),PhysicsUtilities::deltaR(*diHiggsEvt.b2,*hbbCand));

        plotter.getOrMake1DPre("","drhbb",";drhbb",500,0,5)->Fill(drHbb,weight);
        plotter.getOrMake1DPre("","drbb",";drbb",500,0,5)->Fill(drbb,weight);


        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getLeptonTest(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void getLeptonTest(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
