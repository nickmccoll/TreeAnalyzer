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

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if (ht_chs < 400) return false;

        // take only dilepton events
        if(diHiggsEvt.type != DiHiggsEvent::DILEP) return false;

        TString sN = smpName;
        const GenParticle *lep1 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1);
        const GenParticle *lep2 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1);

        // take only e or mu events (get rid of taus)
        if (lep1->absPdgId()==15 || lep2->absPdgId()==15) return false;

        double minDR = 99;
        int idx = -1;
        if (reader_fatjet) {
			for (const auto& fj : reader_fatjet->jets) {
				double dr = PhysicsUtilities::deltaR(*diHiggsEvt.hbb, fj);
				if (dr < minDR) {
	        		bool hasLMTbtag = false;
	        		for (const auto& sj : fj.subJets()) {
	        			if (sj.csv() >= 0.8484) hasLMTbtag = true;
	        		}
	        		if (hasLMTbtag) {
						minDR = dr;
						idx = fj.index();
	        		}
				}
			}
        }
        if (idx == -1) return false;
        if (minDR > 0.8) return false;
        if (reader_fatjet->jets[idx].pt() < 200) return false;

        if (lep1->absPdgId() == 13) {
        	if (lep1->pt() < 26) return false;
        } else {
        	if (lep1->pt() < 30) return false;
        }
        plotter.getOrMake1DPre(sN,"gen2_pt",";p_{T}",100,0,500)->Fill(lep2->pt());
        plotter.getOrMake1DPre(sN,"gen1_pt",";p_{T}",100,0,500)->Fill(lep1->pt());

        return true;
    }
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void getDilepSels(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void getDilepSels(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
