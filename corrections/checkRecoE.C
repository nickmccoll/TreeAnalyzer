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
using namespace std;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    }
    void plotRecoVars(TString sn, const Electron* el, int idx) {
    	if ((*reader_electron->hcalOverEcal)[idx] < 0.15) return;
    	if ((*reader_electron->ecalDrivenSeed)[idx] == 0) return;

    	float HoEtower = (*reader_electron->hcalOverEcalBc)[idx];
    	float HoEcone = (*reader_electron->hcalOverEcal)[idx];
    	float dPhi_sc = (*reader_electron->dPhi_sc)[idx];
    	float dEta_sc = (*reader_electron->dEta_sc)[idx];
    	float dEta_seed = (*reader_electron->dEta_seed)[idx];
    	float sigmaIetaIeta = (*reader_electron->sigmaIetaIeta)[idx];
    	float full5x5_sigmaIetaIeta = (*reader_electron->full5x5_sigmaIetaIeta)[idx];
    	float e1x5 = (*reader_electron->e1x5)[idx];
    	float e5x5 = (*reader_electron->e5x5)[idx];

    	plotter.getOrMake1DPre(sn,"HoE",";H/E (cone)",200,0,1)->Fill(HoEcone,weight);
    	plotter.getOrMake1DPre(sn,"dPhi_sc",";|#Delta#phi_{in}|",100,-0.05,0.05)->Fill(dPhi_sc,weight);
    	plotter.getOrMake1DPre(sn,"dEta_sc",";|#Delta#eta_{in}|",100,-0.03,0.03)->Fill(dEta_sc,weight);
    	plotter.getOrMake1DPre(sn,"dEta_seed",";|#Delta#eta_{in}^{seed}|",100,-0.02,0.02)->Fill(dEta_seed,weight);

    	plotter.getOrMake1DPre(sn,"sigmaIetaIeta",";#sigma_{i#etai#eta}",100,0,0.02)->Fill(sigmaIetaIeta,weight);
    	plotter.getOrMake1DPre(sn,"full5x5_sigmaIetaIeta",";full 5x5 #sigma_{i#etai#eta}",100,0,0.02)->Fill(full5x5_sigmaIetaIeta,weight);
    	plotter.getOrMake1DPre(sn,"e1x5",";E",100,0,2000)->Fill(e1x5,weight);
    	plotter.getOrMake1DPre(sn,"e5x5",";E",100,0,2000)->Fill(e5x5,weight);

    	plotter.getOrMake1DPre(sn,"d0",";",100,0,0.05)->Fill(el->d0(),weight);
    	plotter.getOrMake1DPre(sn,"e1x5_o_e5x5",";E",100,0,1)->Fill(e1x5/e5x5,weight);
    	plotter.getOrMake1DPre(sn,"slurm",";",100,0,0.0001)->Fill(abs(1/el->E() - 1/el->p()),weight);
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
    	if (maxDR >= 0.1) return false;
    	if (reader_electron->electrons[e_idx].pt() < 30) {
    		printf("reco electron below pt threshold\n");
    		return false;
    	}
    	// electron is matched and passes reconstruction
    	const Electron* el = &reader_electron->electrons[e_idx];

    	plotter.getOrMake1DPre(smpName+"e_reco","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
    	if (el->absEta() <= 1.479) plotRecoVars(smpName+"_barrel_",el,e_idx);
    	else plotRecoVars(smpName+"_endcap_",el,e_idx);
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
