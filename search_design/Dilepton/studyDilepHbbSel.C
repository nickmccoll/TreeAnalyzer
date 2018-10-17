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

    bool findHbbCand(const Lepton* lep1, const Lepton* lep2, int genidx) {
    	// lambda function to determine if a FJ is LMT b-tagged (has at least one subjet that passes medium CSV WP)
    	auto isBtag = [&] (const FatJet* fj) {
    		bool hasBtag = false;
    		for (const auto& sj : fj->subJets()) {
    			if (sj.csv() > 0.8484) hasBtag = true;
    		}
    		return hasBtag;
    	};
        // only consider the top two fatjets in pt, provided they are above 200 GeV, then take furthest
    	double minDPhi = 2.0;
    	double minPt = 200;
    	int idx = -1;
    	//
        const MomentumF recodilepton = lep1->p4() + lep2->p4();
        std::vector<const FatJet*> fatjets;
        for (const auto& fj : reader_fatjet->jets) {
        	if (fj.pt() > minPt) fatjets.push_back(&fj);
        }
        std::sort(fatjets.begin(),fatjets.end(), PhysicsUtilities::greaterPTDeref<FatJet>());

        if (fatjets.size() == 0) return false;
        else if (fatjets.size() == 1) {
        	bool separatedFJ = abs(PhysicsUtilities::deltaPhi(*fatjets[0],recodilepton)) > 2.0 && PhysicsUtilities::deltaR(*fatjets[0],*lep1) > 0.8
        			&& PhysicsUtilities::deltaR(*fatjets[0],*lep2) > 0.8;
        	if (separatedFJ && isBtag(fatjets[0])) idx = fatjets[0]->index();

        } else {
            double fj_dr = 0;
            for (int k=0; k<2; k++) {
            	bool separatedFJ = abs(PhysicsUtilities::deltaPhi(*fatjets[k],recodilepton)) > 2.0 && PhysicsUtilities::deltaR(*fatjets[k],*lep1) > 0.8
                	    && PhysicsUtilities::deltaR(*fatjets[k],*lep2) > 0.8;

                if (separatedFJ && isBtag(fatjets[k])) {
                    if (PhysicsUtilities::deltaR(*fatjets[k],recodilepton) > fj_dr) {
                	    fj_dr = PhysicsUtilities::deltaR(*fatjets[k],recodilepton);
                	    idx = fatjets[k]->index();
                	}
                }
            }
        }
        if (idx == -1) return false;

        printf("idx = %d\n\n",idx);
        if (idx != genidx) {
        	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
        	for (const auto& fj : reader_fatjet->jets) {
        		printf("fatjet %d: (E= %f pT=   %f eta= %f phi=  %f)",fj.index(),fj.E(),fj.pt(),fj.eta(),fj.phi());
        		printf(" ----> dR = %f    dPhi = %f: subJet CSV: ",PhysicsUtilities::deltaR(fj,recodilepton), PhysicsUtilities::deltaPhi(fj,recodilepton));
        		for (const auto& sj : fj.subJets()) {
        			printf("%f, ",sj.csv());
        		}
        		printf("\n\n");
        	}
    		printf("recodilepton: (E= %f pT=   %f eta= %f phi=  %f)\n\n",recodilepton.E(),recodilepton.pt(),recodilepton.eta(),recodilepton.phi());
//    		printf("recohww: (E= %f pT=   %f eta= %f phi=  %f)\n\n",recohww.E(),recohww.pt(),recohww.eta(),recohww.phi());
    		return false;
        }
        return true;
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
//        if (genlep1->pt() < genlep1->absPdgId() == 13 ? 26 : 30) return false;

        plotSpectra(sn+"_genacc");
        if (!reader_fatjet) return false;

        // find a match in the RECO FJ collection to the GEN Hbb -> save index of the matched FJ
        double minDR = 99;
        int genidx = -1;
        if (reader_fatjet) {
			for (const auto& fj : reader_fatjet->jets) {
				double dr = PhysicsUtilities::deltaR(*diHiggsEvt.hbb, fj);
				if (dr < minDR) {
					minDR = dr;
					genidx = fj.index();
				}
			}
        }
        if (minDR > 0.8) return false;
        if (reader_fatjet->jets[genidx].pt() < 200) return false;
        std::cout<<"genidx = "<<genidx<<std::endl;
        plotSpectra(sn+"_FJmatch");

        // get the index of the selected FatJet given the two selectedLeptons
        bool goodHbbCand = findHbbCand(selectedLepton,selectedLeptons[1],genidx);
        if (goodHbbCand) {
        	plotSpectra(sn+"_HbbSel");
        	double genlep2_pt = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1->pt() : diHiggsEvt.w1_d1->pt();
        	plotter.getOrMake1DPre(sn+"_HbbSel","gen2_pt",";p_{T}",40,0,400)->Fill(genlep2_pt,weight);
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
