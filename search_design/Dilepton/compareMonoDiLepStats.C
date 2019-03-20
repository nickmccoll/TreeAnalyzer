#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "Configuration/interface/FillerConstants.h"

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

        // take only dilepton and single electron or muon events
        if(!(diHiggsEvt.type == DiHiggsEvent::DILEP || diHiggsEvt.type == DiHiggsEvent::E || diHiggsEvt.type == DiHiggsEvent::MU)) return false;

        // event filter and trigger preselection for all dilepton events
        if (!passEventFilters) return false;
        if (!passTriggerPreselection) return false;

        // throw away dilepton events with taus and same-sign leptons, then record the dilep channel
        bool isdilepEE = false;
        bool isdilepEMU = false;
        bool isdilepMUMU = false;
        if(diHiggsEvt.type == DiHiggsEvent::DILEP) {
        	int lep1id = diHiggsEvt.w1_d1->pdgId();
        	int lep2id = diHiggsEvt.w2_d1->pdgId();

        	if (lep1id<0 == lep2id<0) return false;

            if (abs(lep1id) == 15 || abs(lep2id) == 15) return false;
            else if (abs(lep1id) == 11 && abs(lep2id) == 11) isdilepEE = true;
            else if (abs(lep1id) == 13 && abs(lep2id) == 13) isdilepMUMU = true;
            else if ((abs(lep1id)==13 && abs(lep2id)==11) || ((abs(lep1id)==11 && abs(lep2id)==13))) isdilepEMU = true;
            else {
            	std::cout<<"Error: d1 not a charged lepton"<<std::endl;
            	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
            }
        }
        // require leading lepton to pass our pt WP
        if (selectedLepton->pt() < (selectedLepton->isMuon() ? 26:30)) return false;

        // separate into mono and dilepton cases and fill statistics
        if (diHiggsEvt.type > DiHiggsEvent::DILEP) {
            plotter.getOrMake1DPre(smpName, "evts",";Cat",5,0,5)->Fill(0.5, weight);
        } else if (diHiggsEvt.type == DiHiggsEvent::DILEP) {
        	plotter.getOrMake1DPre(smpName, "evts",";Cat",5,0,5)->Fill(1.5, weight);
        	if (isdilepEE) plotter.getOrMake1DPre(smpName, "evts",";Cat",5,0,5)->Fill(2.5, weight);
        	else if (isdilepMUMU) plotter.getOrMake1DPre(smpName, "evts",";Cat",5,0,5)->Fill(3.5, weight);
        	else if (isdilepEMU) plotter.getOrMake1DPre(smpName, "evts",";Cat",5,0,5)->Fill(4.5, weight);
        	else printf("something went wrong in dilep cat\n");
        } else std::cout<<"More than 2 sel Leps"<<std::endl;

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
