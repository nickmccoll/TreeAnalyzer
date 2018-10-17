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
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"

#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(ht_chs < 400) return false;

        if(diHiggsEvt.type != DiHiggsEvent::E) return false;
        const GenParticle* genlep1 = diHiggsEvt.w1_d1;
        if (genlep1->absPdgId() != 11) {
        	printf("genlep1 not an electron\n");
        	return false;
        }
        if (genlep1->pt() < 30) return false;
    	plotter.getOrMake1DPre(smpName+"e_gen","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

    	const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,2.5);
    	double maxDR = 0.1;
    	int e_idx = -1;
    	for (const auto& e : electrons) {
    		double dr = PhysicsUtilities::deltaR(*genlep1,*e);
    		if (dr < maxDR) {
    			maxDR = dr;
    			e_idx = e->index();
    		}
    	}
    	if (maxDR == 0.1) return false;
    	if (electrons[e_idx]->pt() < 30) {
    		printf("reco electron below pt threshold\n");
    		return false;
    	}
    	plotter.getOrMake1DPre(smpName+"e_reco","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        return true;
    }

    void write(TString fileName){
    	plotter.write(fileName);
    	double den = plotter.getOrMake1DPre(smpName+"e_gen","evts",";M_{X}",50,600,4600)->GetEntries();
    	double num = plotter.getOrMake1DPre(smpName+"e_reco","evts",";M_{X}",50,600,4600)->GetEntries();
    	printf("eff = %f\n",num/den);
    }
    HistGetter plotter;
};

#endif

void checkRecoE(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void checkRecoE(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
