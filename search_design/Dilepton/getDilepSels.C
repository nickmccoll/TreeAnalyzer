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
#include "Processors/Corrections/interface/FatJetScaleFactors.h"

#include "TSystem.h"
using namespace TAna;


class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        // take only dilepton events
        if(diHiggsEvt.type != DiHiggsEvent::DILEP) return false;

        TString sN = smpName;
        const GenParticle *lep1 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1);
        const GenParticle *lep2 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1);

        // take only e or mu events (get rid of taus)
        if (lep1->absPdgId()==15 || lep2->absPdgId()==15) return false;

        bool validHbb = false;
        if (reader_fatjet) {
        	for (const auto& fj : reader_fatjet->jets) {
        		// require at least one medium subjet btag (LMT Btag Cat)
        		bool hasLMTbtag = false;
        		for (const auto& sj : fj.subJets()) {
        			if (sj.csv() >= 0.8484) hasLMTbtag = true;
        		}
        		if (!hasLMTbtag) continue;

        		// require the correct SD mass within a window
        		float hbbmass = hbbFJSFProc->getCorrSDMass(&fj);
        		if (hbbmass < 100 || hbbmass > 150) continue;

        		validHbb = true;
        		sN += "_validfj";
        		break;
        	}
        }
        // fill for the inclusive, 1lep_pt>26, 2lep_pt>15&5, and the former two with ht>400 additionally (5 cases)
        plotter.getOrMake1DPre(sN, "evts", ";cats",5,0,5)->Fill(0.5, weight);
        if (lep1->pt() > 26) {
            plotter.getOrMake1DPre(sN, "evts", ";cats",5,0,5)->Fill(1.5, weight);
            if (ht_chs > 400) plotter.getOrMake1DPre(sN, "evts", ";cats",5,0,5)->Fill(3.5, weight);
        }
        if (lep1->pt() > 15 && lep2->pt() > 5) {
            plotter.getOrMake1DPre(sN, "evts", ";cats",5,0,5)->Fill(2.5, weight);
            if (ht_chs > 400) plotter.getOrMake1DPre(sN, "evts", ";cats",5,0,5)->Fill(4.5, weight);
        }
        return true;
    }
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void getDilepSels(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getDilepSels(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
