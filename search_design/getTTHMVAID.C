
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "Processors/EventSelection/interface/EventSelection.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/DileptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/Hww2lSolver.h"

#include "Processors/Variables/interface/SignalHelper.h"

using namespace TAna;
using namespace std;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

        turnOffCorr(CORR_TRIG);
        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
//        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);

        fjProc2      .reset(new FatJetProcessor ());
    }

    const FatJet* getMatchedFJ(const MomentumF& genJet, bool doHbb, const std::vector<FatJet>& fatjets, const std::vector<FatJet>& fatjets_nolep) {
    	double nearestDR = 10;
    	if (genJet.pt() < 200) return 0;
    	if (doHbb){
    		if (!fatjets.size()) return 0;
        	int idx = PhysicsUtilities::findNearestDR(genJet,fatjets,nearestDR,0.2);
        	if (idx < 0) return 0;
        	else return &fatjets[idx];
    	} else {
    		if (!fatjets_nolep.size()) return 0;
        	int idx = PhysicsUtilities::findNearestDR(genJet,fatjets_nolep,nearestDR,0.2);
        	if (idx < 0) return 0;
        	else return &fatjets_nolep[idx];
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

    void plot(TString sn, const Lepton *lep1, const Lepton *lep2, bool is1l) {

    	TString chS = is1l ? "1l" : "2l";
		plotter.getOrMake1DPre(sn+"_"+chS,"ht",";",100,0,4000)->Fill(ht,weight);

    	TString s1 = lep1->isMuon() ? "mu1" :"e1";
		plotter.getOrMake1DPre(sn+"_"+chS+"_"+s1,"tthMVA",";",100,-1,1)->Fill(lep1->ttHMVA(),weight);

    	bool passStd1 = lep1->isMuon() ? ((Muon*)lep1)->passLooseID() : ((Electron*)lep1)->passMedID_noIso();
    	bool passtth1 = lep1->isMuon() ? lep1->ttHMVA() > 0.85 : lep1->ttHMVA() > 0.8;

		if(passStd1) {
			plotter.getOrMake1DPre(sn+"_"+chS+"_passSTD_"+s1,"pt",";",100,0,1000)->Fill(lep1->pt(),weight);
			if (is1l) plotter.getOrMake1DPre(sn+"_"+chS+"_passSTD","pt",";",100,0,1000)->Fill(lep1->pt(),weight);
		}
		if(passtth1) {
			plotter.getOrMake1DPre(sn+"_"+chS+"_passTTH_"+s1,"pt",";",100,0,1000)->Fill(lep1->pt(),weight);
			if (is1l) plotter.getOrMake1DPre(sn+"_"+chS+"_passTTH","pt",";",100,0,1000)->Fill(lep1->pt(),weight);
		}

		if(!is1l) {
	    	TString s2 = lep2->isMuon() ? "mu2" :"e2";
			plotter.getOrMake1DPre(sn+"_"+chS+"_"+s2,"tthMVA",";",100,-1,1)->Fill(lep2->ttHMVA(),weight);

	    	bool passStd2 = lep2->isMuon() ? ((Muon*)lep2)->passLooseID() : ((Electron*)lep2)->passMedID_noIso();
	    	bool passtth2 = lep2->isMuon() ? lep2->ttHMVA() > 0.85 : lep2->ttHMVA() > 0.8;

			if(passStd2) {
				plotter.getOrMake1DPre(sn+"_"+chS+"_passSTD_"+s2,"pt",";",100,0,1000)->Fill(lep2->pt(),weight);
				if(passStd1) plotter.getOrMake1DPre(sn+"_"+chS+"_passSTD","pt",";",100,0,1000)->Fill(lep2->pt(),weight);
			}
			if(passtth2) {
				plotter.getOrMake1DPre(sn+"_"+chS+"_passTTH_"+s2,"pt",";",100,0,1000)->Fill(lep2->pt(),weight);
				if(passtth1) plotter.getOrMake1DPre(sn+"_"+chS+"_passTTH","pt",";",100,0,1000)->Fill(lep2->pt(),weight);
			}
		}

    }

    const GenParticle *getGenLepFromTau(const GenParticle* tau) {
    	if (tau->absPdgId() != ParticleInfo::p_tauminus) return 0;
    	if (!ParticleInfo::isLastInChain(tau)) {
    		cout<<"tau not last in chain"<<endl;
    		return 0;
    	}

    	const GenParticle* gp=0;
    	int nlepsfromtau = 0;
    	for (unsigned int k = 0; k < tau->numberOfDaughters(); k++) {
    		const auto* dau = tau->daughter(k);
    		if (dau->absPdgId() == ParticleInfo::p_muminus || dau->absPdgId() == ParticleInfo::p_eminus) {

    			if (nlepsfromtau >= 1) {
    				cout<<"already a lep in the tau decay!!"<<endl;
    				continue;
    			}
    			gp = dau;
    			nlepsfromtau++;
    		}
    	}
    	if (!gp) return 0;
    	if (nlepsfromtau > 1) {
    		cout<<"Found "<<nlepsfromtau<<" leps in tau decay!!!!"<<endl;
    		return 0;
    	}

    	return gp;
    }

    bool passLepSelNoID(const Lepton* lep1, const Lepton* lep2) {
    	bool passIP1 = fabs(lep1->dz()) < 0.1 && fabs(lep1->d0()) < 0.05;
    	bool passIP2 = fabs(lep2->dz()) < 0.1 && fabs(lep2->d0()) < 0.05;

    	bool passIso = lep1->miniIso() < 0.15 && lep2->miniIso() < 0.15;

    	bool pass = passIP1 && passIP2 && passIso;

    	if (pass) return true;
    	else      return false;
    }

    void doSignalDilep(TString sn) {
    	const GenParticle* glep1_0 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1;
    	const GenParticle* glep2_0 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1;

    	if (glep1_0->pdgId()<0 == glep2_0->pdgId()<0) return;
    	const GenParticle *glep1=0;
    	const GenParticle *glep2=0;

    	if (glep1_0->absPdgId() == 15) glep1 = getGenLepFromTau(glep1_0);
    	else glep1 = glep1_0;

    	if (glep2_0->absPdgId() == 15) glep2 = getGenLepFromTau(glep2_0);
    	else glep2 = glep2_0;

    	if (!glep1 || !glep2) return;

//    	if (glep1->pt() < 10 || glep2->pt() < 10) return;
//    	bool goodPt1 = glep1->absPdgId() == 13 ? (glep1->pt() > 27) : (glep1->pt() > 30);
//    	bool goodPt2 = glep2->absPdgId() == 13 ? (glep2->pt() > 27) : (glep2->pt() > 30);
//    	if (!goodPt1 && !goodPt2) return;

    	const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
    	const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,2.5);
    	const Lepton* matchLep1 = getMatchedLepton(glep1,muons,electrons,0.1,true);
    	const Lepton* matchLep2 = getMatchedLepton(glep2,muons,electrons,0.1,true);

    	if (!matchLep1 || !matchLep2) return;
    	if ((matchLep1->isMuon() == matchLep2->isMuon()) && (matchLep1->index() == matchLep2->index())) return; // discard if these are the same RECO lep

    	// RECO lepton pt cuts
//    	if (matchLep1->pt() < 10 || matchLep2->pt() < 10) return;
//    	goodPt1 = matchLep1->isMuon() ? (matchLep1->pt() > 27) : (matchLep1->pt() > 30);
//    	goodPt2 = matchLep2->isMuon() ? (matchLep2->pt() > 27) : (matchLep2->pt() > 30);
//    	if (!goodPt1 && !goodPt2) return;

//    	if (!passLepSelNoID(matchLep1,matchLep2)) return;

    	plot(sn,matchLep1,matchLep2,false);
    }

    void doSignalSingleLep(TString sn) {
    	const GenParticle* glep0 = 0;
    	if (ParticleInfo::isLepton(diHiggsEvt.w1_d1->absPdgId())) glep0 = diHiggsEvt.w1_d1;
    	else if (ParticleInfo::isLepton(diHiggsEvt.w2_d1->absPdgId())) glep0 = diHiggsEvt.w2_d1;

    	if (!glep0) return;
    	const GenParticle *glep=0;

    	if (glep0->absPdgId() == 15) glep = getGenLepFromTau(glep0);
    	else glep = glep0;

    	if (!glep) return;

//    	bool goodPt = glep->absPdgId() == 13 ? (glep->pt() > 27) : (glep->pt() > 30);
//    	if (!goodPt) return;

    	const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
    	const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,1.479);
    	const Lepton* matchLep1 = getMatchedLepton(glep,muons,electrons,0.1,true);
    	if (!matchLep1) return;

    	// RECO lepton pt cuts
//    	goodPt = matchLep1->isMuon() ? (matchLep1->pt() > 27) : (matchLep1->pt() > 30);
//    	if (!goodPt) return;

//    	if (matchLep1->miniIso() > 0.2) return;
//    	if (matchLep1->sip3D() > 4 || std::fabs(matchLep1->dz()) > 0.1 || std::fabs(matchLep1->d0()) > 0.05) return;

    	plot(sn,matchLep1,0,true);
    }

    void testSigHelper(TString sn, SignalHelper& helper) {
    	bool passGenLepPt = false;
    	bool passRecoLepPt = false;
    	bool passIsoIP = false;

//    	cout<<"begin"<<endl;
    	if (helper.type == DiHiggsEvent::DILEP) {
    		if (!helper.recolep1 || !helper.recolep2 || !helper.genlep1 || !helper.genlep2) return;

    		bool hasHighPt = (helper.genlep1->absPdgId() == 13 ? helper.genlep1->pt() > 27 : helper.genlep1->pt() > 30)
    				|| (helper.genlep2->absPdgId() == 13 ? helper.genlep2->pt() > 27 : helper.genlep2->pt() > 30);
    		passGenLepPt = hasHighPt && (helper.genlep1->pt() > 10 && helper.genlep2->pt() > 10);

    		hasHighPt = (helper.recolep1->isMuon() ? helper.recolep1->pt() > 27 : helper.recolep1->pt() > 30)
    				|| (helper.recolep2->isMuon() ? helper.recolep2->pt() > 27 : helper.recolep2->pt() > 30);
    		passRecoLepPt = hasHighPt && (helper.recolep1->pt() > 10 && helper.recolep2->pt() > 10);

    		passIsoIP = helper.recolep1->miniIso() < 0.15 && helper.recolep2->miniIso() < 0.15
    				&& helper.recolep1->sip3D() < 4.0 && fabs(helper.recolep1->dz()) < 0.1 && fabs(helper.recolep1->d0()) < 0.05
					&& helper.recolep2->sip3D() < 4.0 && fabs(helper.recolep2->dz()) < 0.1 && fabs(helper.recolep2->d0()) < 0.05;
    	} else {
    		if (!helper.recolep1 || !helper.genlep1) return;

    		passGenLepPt = (helper.genlep1->absPdgId() == 13 ? helper.genlep1->pt() > 27 : helper.genlep1->pt() > 30);
    		passRecoLepPt = (helper.recolep1->isMuon() ? helper.recolep1->pt() > 27 : helper.recolep1->pt() > 30);

    		passIsoIP = helper.recolep1->miniIso() < 0.2 && helper.recolep1->sip3D() < 4.0
    				&& fabs(helper.recolep1->dz()) < 0.1 && fabs(helper.recolep1->d0()) < 0.05;
    	}
//    	cout<<"dbg0"<<endl;

    	if (!passGenLepPt || !passRecoLepPt || !passIsoIP) return;
//    	cout<<"dbg1"<<endl;

    	plot(sn,helper.recolep1,helper.recolep2,helper.type == DiHiggsEvent::DILEP ? false : true);
//    	cout<<"dbg2"<<endl;

    }

    void testBkg(TString sn, bool doStd) {

    	auto passElSelNoID1 = [&](const Electron* el) {
    		if (el->pt() < 30) return false;
    		if (el->absEta() > 1.479) return false;
    		if (el->sip3D() > 4.0) return false;
    		if (fabs(el->dz()) > 0.1) return false;
    		if (fabs(el->d0()) > 0.05) return false;
    		if (el->miniIso() > 0.2) return false;
    	};

    	auto passElSelNoID2 = [&](const Electron* el) {
    		if (el->pt() < 10) return false;
    		if (el->absEta() > 2.5) return false;
    		if (fabs(el->dz()) > 0.1) return false;
    		if (fabs(el->d0()) > 0.05) return false;
    		if (el->miniIso() > 0.15) return false;
    	};

    	auto passMuSelNoID1 = [&](const Electron* el) {
    		if (el->pt() < 27) return false;
    		if (el->absEta() > 2.4) return false;
    		if (el->sip3D() > 4.0) return false;
    		if (fabs(el->dz()) > 0.1) return false;
    		if (fabs(el->d0()) > 0.05) return false;
    		if (el->miniIso() > 0.2) return false;
    	};

    	auto passMuSelNoID2 = [&](const Electron* el) {
    		if (el->pt() < 10) return false;
    		if (el->absEta() > 2.4) return false;
    		if (fabs(el->dz()) > 0.1) return false;
    		if (fabs(el->d0()) > 0.05) return false;
    		if (el->miniIso() > 0.15) return false;
    	};

    	auto getLeptons = [&](bool is1l) {
        	std::vector<const Electron*> els;
        	std::vector<const Muon*> mus;

        	for (const auto& el : reader_electron->electrons) {
        		if (is1l) {
        			if (!passElSelNoID1(&el)) continue;
        		} else {
        			if (!passElSelNoID2(&el)) continue;
        		}
    			if (doStd) {
    				if (!el.passMVA90ID_noIso()) continue;
    			} else {
    				if (el.ttHMVA() < 0.8) continue;
    			}
        		els.push_back(&el);
        	}
        	for (const auto& mu : reader_muon->muons) {
        		if (is1l) {
        			if (!passMuSelNoID1(&mu)) continue;
        		} else {
        			if (!passMuSelNoID2(&mu)) continue;
        		}
    			if (doStd) {
    				if (!mu.passMedID()) continue;
    			} else {
    				if (mu.ttHMVA() < 0.85) continue;
    			}
        		mus.push_back(&mu);
        	}

        	std::vector<const Lepton*> leps;
            leps.reserve(els.size() + mus.size());
            for(const auto* l : els) leps.push_back(l);
            for(const auto* l : mus) leps.push_back(l);
            std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Lepton>());
            return leps;
    	};

    	vector<const Lepton*> singleleps = getLeptons(true);
    	vector<const Lepton*> dileps = getLeptons(false);

    	bool goodDilepPt = false;
    	if (dileps.size() >= 2) {
    		goodDilepPt = (dileps[0]->pt() > (dileps[0]->isMuon() ? 27 : 30))
    				|| (dileps[1]->pt() > (dileps[1]->isMuon() ? 27 : 30));
    	}

    	TString chS;
    	if (dileps.size() == 2 && goodDilepPt) chS = "2l";
    	else if (dileps.size() < 2) chS = "1l";
    	else return;

    	plotter.getOrMake1DPre(sn+"_"+chS,"ht",";ht",1000,0,4000)->Fill(ht,weight);
//    	if (!doStd) {
//        	plotter.getOrMake1DPre(sn+"_"+chS,"ht",";ht",1000,0,4000)->Fill(ht,weight);
//    	}

    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht < 400) return false;
	    if(!EventSelection::passTriggerSuite2017(*reader_event)) return false;

	    TString sn = smpName;
        if (mcProc == FillerConstants::TTBAR && smDecayEvt.nLepsTT >= 0 && smDecayEvt.nLepsTT <= 2) {
        	sn += TString::Format("%d",smDecayEvt.nLepsTT);
        }

        if (isSignal()) {
            SignalHelper sigInfo(diHiggsEvt,reader_muon,reader_electron);
        	sigInfo.minElRecoPt = 10;
        	sigInfo.minMuRecoPt = 10;
        	sigInfo.maxMuRecoEta = 2.4;

        	if (diHiggsEvt.type == DiHiggsEvent::DILEP) {
//        		doSignalDilep(sn);
        		sigInfo.maxElRecoEta = 2.5;
        		sigInfo.setRecoLeptons(0.1);
        		testSigHelper(sn,sigInfo);
        	} else if (diHiggsEvt.type >= DiHiggsEvent::TAU_MU) {
//        		doSignalSingleLep(sn);
        		sigInfo.maxElRecoEta = 1.479;
        		sigInfo.setRecoLeptons(0.1);
        		testSigHelper(sn,sigInfo);
        	}

        } else {
        	testBkg(sn+"_passSTD",true);
        	testBkg(sn+"_passTTH",false);
        }

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<FatJetProcessor>        fjProc2     ;

};

#endif

void getTTHMVAID(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
