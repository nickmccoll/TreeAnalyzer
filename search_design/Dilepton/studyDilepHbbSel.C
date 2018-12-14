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

    void plotSpectra(TString sn, const FatJet hbbjet, const Lepton* lep1, const Lepton* lep2) {
    	const MomentumF dilepmom = lep1->p4() + lep2->p4();
    	double dR_bbll = PhysicsUtilities::deltaR(hbbjet,dilepmom);
    	double dphi_bbll = PhysicsUtilities::deltaPhi(hbbjet,dilepmom);

    	plotter.getOrMake1DPre(sn,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
    	plotter.getOrMake1DPre(sn,"ptbb",";p_{T}",400,0,4000)->Fill(hbbjet.pt(),weight);
    	plotter.getOrMake1DPre(sn,"dR_bbll",";#DeltaR_{bb,ll}",100,0,7)->Fill(dR_bbll,weight);
    	plotter.getOrMake1DPre(sn,"dphi_bbll",";#Delta#Phi_{bb,ll}",100,0,3.14)->Fill(abs(dphi_bbll),weight);
    }

    const FatJet* findHbbCand(const Lepton* lep1, const Lepton* lep2, double minDPhi, double minDR, bool takeFurthestFJ) {
    	// lambda function to determine if a FJ is LMT b-tagged (has at least one subjet that passes medium CSV WP)
    	auto isBtag = [&] (const FatJet* fj) {
    		bool hasBtag = false;
    		for (const auto& sj : fj->subJets()) {
    			if (sj.csv() > 0.8484) hasBtag = true;
    		}
    		return hasBtag;
    	};
    	// lambda function to determine if a FJ has two SJs each with pt > 20 and eta < 2.4
    	auto hasGoodSJs = [&] (const FatJet* fj) {
    		bool goodSJs = false;
    		int nGoodSJ = 0;
    		for (const auto& sj : fj->subJets()) {
    			if (sj.pt() > 20 && sj.absEta() < 2.4) nGoodSJ++;
    		}
    		if (nGoodSJ > 1) goodSJs = true;
    		return goodSJs;
    	};
    	double minPt = 200;
    	int idx = -1;

    	// only consider fatjets with pt > 200 and sort them by pt
        const MomentumF recodilepton = lep1->p4() + lep2->p4();
        std::vector<const FatJet*> fatjets;
        for (const auto& fj : reader_fatjet->jets) {
        	if (fj.pt() > minPt) fatjets.push_back(&fj);
        }
        std::sort(fatjets.begin(),fatjets.end(), PhysicsUtilities::greaterPTDeref<FatJet>());

        // make sure any considered fatjet is B-tagged, has 2 good subJets, is sufficiently far from dileptons. If several, take furthest FJ
        if (fatjets.size() == 0) return 0;
        else if (fatjets.size() == 1) {
        	bool separatedFJ = abs(PhysicsUtilities::deltaPhi(*fatjets[0],recodilepton)) > minDPhi && abs(PhysicsUtilities::deltaR(*fatjets[0],recodilepton)) > minDR
        	                   && PhysicsUtilities::deltaR(*fatjets[0],*lep1) > 0.8 && PhysicsUtilities::deltaR(*fatjets[0],*lep2) > 0.8;
        	if (separatedFJ && isBtag(fatjets[0]) && hasGoodSJs(fatjets[0])) idx = fatjets[0]->index();
        } else {
            double dr = 0;
            double pt = 0;
            for (int k=0; k<2; k++) {
            	bool separatedFJ = abs(PhysicsUtilities::deltaPhi(*fatjets[k],recodilepton)) > minDPhi && abs(PhysicsUtilities::deltaR(*fatjets[k],recodilepton)) > minDR
            	        	                   && PhysicsUtilities::deltaR(*fatjets[k],*lep1) > 0.8 && PhysicsUtilities::deltaR(*fatjets[k],*lep2) > 0.8;

                if (separatedFJ && isBtag(fatjets[k]) && hasGoodSJs(fatjets[k])) {
                	if (takeFurthestFJ) {
                        if (PhysicsUtilities::deltaR(*fatjets[k],recodilepton) > dr) {
                	        dr = PhysicsUtilities::deltaR(*fatjets[k],recodilepton);
                	        idx = fatjets[k]->index();
                	    }
                	} else {
                		if (fatjets[k]->pt() > pt) {
                			pt = fatjets[k]->pt();
                			idx = fatjets[k]->index();
                		}
                	}
                }
            }
        }
        if (idx < 0) return 0;

        // debug statements
/*        printf("idx = %d\n\n",idx);
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
    		return 0;
        } */
        return &reader_fatjet->jets[idx];
    }
    bool passLepSel(const Lepton* lep) {
    	double dz = 0.1;
    	double d0 = 0.05;
    	double sip = 4.0;

    	bool passIP = lep->dz() < dz && lep->d0() < d0 && lep->sip3D() < sip;
    	if (!passIP) return false;

    	bool passID = false;
    	bool passISO = false;
    	if (lep->isMuon()) {
    		passISO = lep->miniIso() < 0.2;
    		passID = ((const Muon*)lep)->passMed16ID();
    	} else {
    		passISO = lep->miniIso() < 0.2;
    		passID = ((const Electron*)lep)->passMedID_noISO();
    	}
    	if (passID && passISO) return true;
    	else return false;
    }
    void testHbbSel(TString sn, const FatJet* hbbjet, const Lepton* lep1, const Lepton* lep2) {
    	static const std::vector<double> dphis = {1.6,1.8,2.0,2.2,2.4};
    	static const std::vector<double> drs = {1.6,1.8,2.0,2.2,2.4};
    	static const std::vector<bool> bools = {true,false};

    	plotSpectra(sn+"_beforeSel",*hbbjet,lep1,lep2);
    	double DPHI = PhysicsUtilities::deltaPhi(*hbbjet,lep1->p4()+lep2->p4());
    	double DR = PhysicsUtilities::deltaR(*hbbjet,lep1->p4()+lep2->p4());
    	for (const bool takeFurthest : bools) {
    		for (const auto& dphi : dphis) {
    			for (const auto& dr : drs) {
    				if (dr < dphi) continue;
    				bool matchPass = abs(DPHI) > dphi && DR > dr;
    				if (!matchPass) continue;

    				TString histname = TString::Format("_dphi_%.1f_dr_%.1f",dphi,dr);
    				histname.ReplaceAll(".","p");
    				if (takeFurthest) histname += "_selDR";
    				else histname += "_selPT";

    				plotSpectra(sn+histname+"_matchPass",*hbbjet,lep1,lep2);
    				const FatJet* recoHbb = findHbbCand(lep1,lep2,dphi,dr,takeFurthest);
    				if (recoHbb && recoHbb->index() == hbbjet->index()) plotSpectra(sn+histname+"_goodSel",*recoHbb,lep1,lep2);
    			}
    		}
    	}
    }
    const Lepton* getMatchedLepton(const GenParticle* genLep,const std::vector<const Muon*> muons, const std::vector<const Electron*> electrons, double maxDR, bool chargeMatch) {
       if(genLep->absPdgId() == ParticleInfo::p_muminus){
           double nearestDR =10;
           int idx = -1;
       	   for (const auto& mu : muons) {
       	       if (chargeMatch) {
       		       if ((mu->q() > 0) == (genLep->pdgId() > 0)) continue;
       	       }
       	       double dr = PhysicsUtilities::deltaR(*genLep,*mu);
       	       if (dr < nearestDR) {
       	    	   nearestDR = dr;
           	       idx = mu->index();
       	       }
       	   }
       	   if (nearestDR > maxDR) return 0;
       	   if (idx < 0) return 0;
       	   return &reader_muon->muons[idx];
       } else {
           double nearestDR =10;
           int idx = -1;
       	   for (const auto& el : electrons) {
       	       if (chargeMatch) {
       		       if ((el->q() > 0) == (genLep->pdgId() > 0)) continue;
       	       }
       	       double dr = PhysicsUtilities::deltaR(*genLep,*el);
       	       if (dr < nearestDR) {
       	    	   nearestDR = dr;
           	       idx = el->index();
       	       }
       	   }
       	   if (nearestDR > maxDR) return 0;
       	   if (idx < 0) return 0;
       	   return &reader_electron->electrons[idx];
       }
    }
    const FatJet* getMatchedFJ(){
        double nearestDR = 99;
        int idx = -1;
        if (reader_fatjet) {
            for (const auto& fj : reader_fatjet->jets) {
        		double dr = PhysicsUtilities::deltaR(*diHiggsEvt.hbb, fj);
        		if (dr < nearestDR) {
        			nearestDR = dr;
        			idx = fj.index();
        		}
        	}
        }
        if (nearestDR > 0.4) return 0;
        if (idx < 0) return 0;
        else return &reader_fatjet->jets[idx];
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        TString sn = smpName;

        // denominator selection: take Dilepton GEN events
        if(diHiggsEvt.type != DiHiggsEvent::DILEP) return false;

        // order the genleptons by pt, then keep only e and mu events with opposite-sign leptons
        const GenParticle *genlep1 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1);
        const GenParticle *genlep2 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1);
        int lep1id = diHiggsEvt.w1_d1->pdgId();
        int lep2id = diHiggsEvt.w2_d1->pdgId();

        // throw away dilepton events with taus and same-sign leptons, then record the dilep channel
        if (lep1id<0 == lep2id<0) return false;
        if (abs(lep1id) == 15 || abs(lep2id) == 15) return false;
        else if (abs(lep1id) == 11 && abs(lep2id) == 11) sn += "_ee";
        else if (abs(lep1id) == 13 && abs(lep2id) == 13) sn += "_mumu";
        else if ((abs(lep1id)==13 && abs(lep2id)==11) || ((abs(lep1id)==11 && abs(lep2id)==13))) sn += "_emu";
        else {
            std::cout<<"Error: d1 not a charged lepton"<<std::endl;
            ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
        }
    	plotter.getOrMake1DPre(sn+"_gendilep","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

        // GEN ACCEPTANCE: impose some Hbb selection on the gen hbb object, the gen b quarks, and the gen leptons
        if (diHiggsEvt.hbb->pt() < 200 || diHiggsEvt.hbb->absEta() > 2.4) return false;
        if (diHiggsEvt.b1->pt() < 20 || diHiggsEvt.b2->pt() < 20 || diHiggsEvt.b1->absEta() > 2.4 || diHiggsEvt.b2->absEta() > 2.4) return false;
        bool passGenLepPT = (genlep1->absPdgId() == 13 ? genlep1->pt() > 26 : genlep1->pt() > 30) && (genlep2->pt() > 10);
        if (!passGenLepPT) return false;

    	plotter.getOrMake1DPre(sn+"_genacc","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        if (!reader_fatjet) return false;

        // find a match in the RECO FJ collection to the GEN Hbb -> save index of the matched FJ
        const FatJet* matchedHbb = getMatchedFJ();
        if (!matchedHbb) return false;
    	plotter.getOrMake1DPre(sn+"_FJmatch","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

        const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
        const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,2.5);
        const Lepton* matchLep1 = getMatchedLepton(genlep1,muons,electrons,0.1,true);
        const Lepton* matchLep2 = getMatchedLepton(genlep2,muons,electrons,0.1,true);
        if (!(matchLep1 && matchLep2)) return false;
    	plotter.getOrMake1DPre(sn+"_slurm","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        bool passRecoPt = matchLep1->isMuon() ? matchLep1->pt() > 26 : matchLep1->pt() > 30;
        if (!passRecoPt) return false;
    	plotter.getOrMake1DPre(sn+"_churm","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

        testHbbSel(sn,matchedHbb,matchLep1,matchLep2);
        return true;
    }

    void write(TString fileName){
    	plotter.write(fileName);
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
