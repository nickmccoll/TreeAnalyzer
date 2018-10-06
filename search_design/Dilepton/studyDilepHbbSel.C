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
#include "TFile.h"
#include "TH1.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    }

    void plotSpectra(TString sn) {
    	plotter.getOrMake1DPre(sn,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process >= FillerConstants::ZJETS && reader_event->process <= FillerConstants::TTX )
            smpName = "other";
        TString sn = "";

        // denominator selection: take Dilepton GEN events and require two RECO leptons passing our current criteria
        if(diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        if (selectedLeptons.size() < 2) return false;
        if (selectedLepton->pt() < (selectedLepton->isMuon() ? 26:30)) return false;

        // order the genleptons by pt, then keep only e and mu events with opposite-sign leptons
        const GenParticle *genlep1 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1);
        const GenParticle *genlep2 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1);
        int lep1id = diHiggsEvt.w1_d1->pdgId();
        int lep2id = diHiggsEvt.w2_d1->pdgId();

        // throw away dilepton events with taus and same-sign leptons, then record the dilep channel
        if (lep1id<0 == lep2id<0) return false;
        if (abs(lep1id) == 15 || abs(lep2id) == 15) return false;
        else if (abs(lep1id) == 11 && abs(lep2id) == 11) sn += "ee";
        else if (abs(lep1id) == 13 && abs(lep2id) == 13) sn += "mumu";
        else if ((abs(lep1id)==13 && abs(lep2id)==11) || ((abs(lep1id)==11 && abs(lep2id)==13))) sn += "emu";
        else {
            std::cout<<"Error: d1 not a charged lepton"<<std::endl;
            ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
        }
        plotSpectra(sn+"_baseline");

        // GEN ACCEPTANCE: impose some Hbb selection on the gen object and the gen b quarks
        if (diHiggsEvt.hbb->pt() < 200 || diHiggsEvt.hbb->absEta() > 2.4) return false;
        if (diHiggsEvt.b1->pt() < 20 || diHiggsEvt.b2->pt() < 20 || diHiggsEvt.b1->absEta() > 2.4 || diHiggsEvt.b2->absEta() > 2.4) return false;

        plotSpectra(sn+"_genacc");

        // find a match in the RECO FJ collection to the GEN Hbb -> save index of the matched FJ
        double minDR = 99;
        int idx = -1;
        if (reader_fatjet) {
			for (const auto& fj : reader_fatjet->jets) {
				double dr = PhysicsUtilities::deltaR(*diHiggsEvt.hbb, fj);
				if (dr < minDR) {
					minDR = dr;
					idx = fj.index();
				}
			}
        }
        if (minDR > 0.8) return false;
        if (reader_fatjet->jets[idx].pt() < 200) return false;
        std::cout<<"idx = "<<idx<<std::endl;
        plotSpectra(sn+"_FJmatch");

        // make Hbb selection based on distance from dilepton system -> save index of selected FJ
        const MomentumF recodilepton = selectedLepton->p4() + selectedLeptons[1]->p4();
        const MomentumF recohww = HiggsSolver::getInvisible(reader_event->met,recodilepton);
        std::vector<int> indices;
        printf("indices = ");
        for (const auto& fj : reader_fatjet->jets) {
            if (abs(PhysicsUtilities::deltaPhi(fj,recodilepton)) > 2.0 && PhysicsUtilities::deltaR(fj,*selectedLepton) > 0.8
                && PhysicsUtilities::deltaR(fj,*selectedLeptons[1]) > 0.8) {
            	indices.push_back(fj.index());
            	printf("%d, ",fj.index());
            }
        }
        printf("\n");
        int idx2 = -1;
        double maxpt = 200;
        if (indices.size() == 0) return false;
        else if (indices.size() > 1) {
        	for (int k : indices) {
        		if (reader_fatjet->jets[k].pt() > maxpt) {
//        			printf("fj %d pt = %f; max is %f\n",k,reader_fatjet->jets[k].pt(),maxpt);
        			maxpt = reader_fatjet->jets[k].pt();
        			idx2 = k;
//        			printf("idx2 is %d\n",idx2);
        		}
        	}
        }
        else idx2 = indices[0];
        // if the two indices match, then it is a good selection
        printf("idx2 = %d\n\n",idx2);
        if (idx == idx2) plotSpectra(sn+"_HbbSel");
        else {
        	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
        	for (const auto& fj : reader_fatjet->jets) {
        		printf("fatjet %d: (E= %f pT=   %f eta= %f phi=  %f)",fj.index(),fj.E(),fj.pt(),fj.eta(),fj.phi());
        		printf(" ----> dR = %f    dPhi = %f\n\n",PhysicsUtilities::deltaR(fj,recodilepton), PhysicsUtilities::deltaPhi(fj,recodilepton));
        	}
    		printf("recodilepton: (E= %f pT=   %f eta= %f phi=  %f)\n\n",recodilepton.E(),recodilepton.pt(),recodilepton.eta(),recodilepton.phi());
    		printf("recohww: (E= %f pT=   %f eta= %f phi=  %f)\n\n",recohww.E(),recohww.pt(),recohww.eta(),recohww.phi());
        }
        return true;
    }

    void write(TString fileName){
    	plotter.write(fileName);
    	float ee_num = plotter.getOrMake1DPre("ee_HbbSel","evts",";M_{X}",50,600,4600)->GetEntries();
    	float ee_den = plotter.getOrMake1DPre("ee_FJmatch","evts",";M_{X}",50,600,4600)->GetEntries();
    	float emu_num = plotter.getOrMake1DPre("emu_HbbSel","evts",";M_{X}",50,600,4600)->GetEntries();
    	float emu_den = plotter.getOrMake1DPre("emu_FJmatch","evts",";M_{X}",50,600,4600)->GetEntries();
    	float mumu_num = plotter.getOrMake1DPre("mumu_HbbSel","evts",";M_{X}",50,600,4600)->GetEntries();
    	float mumu_den = plotter.getOrMake1DPre("mumu_FJmatch","evts",";M_{X}",50,600,4600)->GetEntries();

    	printf("\nMass = %d GeV\n",signal_mass);
    	printf("emu eff = %f\n",emu_num/emu_den);
    	printf("ee eff = %f\n",ee_num/ee_den);
    	printf("mumu eff = %f\n",mumu_num/mumu_den);
    }
    HistGetter plotter;
};

#endif

void studyDilepHbbSel(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void studyDilepHbbSel(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
