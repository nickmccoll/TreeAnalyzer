
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
    				&& fabs(helper.recolep1->dz()) < 0.1 && fabs(helper.recolep1->d0()) < 0.05
					&& fabs(helper.recolep2->dz()) < 0.1 && fabs(helper.recolep2->d0()) < 0.05;
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
    		return true;
    	};

    	auto passElSelNoID2 = [&](const Electron* el) {
    		if (el->pt() < 10) return false;
    		if (el->absEta() > 2.5) return false;
    		if (fabs(el->dz()) > 0.1) return false;
    		if (fabs(el->d0()) > 0.05) return false;
    		if (el->miniIso() > 0.15) return false;
    		return true;
    	};

    	auto passMuSelNoID1 = [&](const Muon* mu) {
    		if (mu->pt() < 27) return false;
    		if (mu->absEta() > 2.4) return false;
    		if (mu->sip3D() > 4.0) return false;
    		if (fabs(mu->dz()) > 0.1) return false;
    		if (fabs(mu->d0()) > 0.05) return false;
    		if (mu->miniIso() > 0.2) return false;
    		return true;
    	};

    	auto passMuSelNoID2 = [&](const Muon* mu) {
    		if (mu->pt() < 10) return false;
    		if (mu->absEta() > 2.4) return false;
    		if (fabs(mu->dz()) > 0.1) return false;
    		if (fabs(mu->d0()) > 0.05) return false;
    		if (mu->miniIso() > 0.15) return false;
    		return true;
    	};

    	auto getLeptons = [&](bool is1l, bool doID) {
        	std::vector<const Electron*> els;
        	std::vector<const Muon*> mus;

        	for (const auto& el : reader_electron->electrons) {
        		if (is1l) {
        			if (!passElSelNoID1(&el)) continue;
        		} else {
        			if (!passElSelNoID2(&el)) continue;
        		}
        		if(doID) {
        			if (doStd) {
        				if (is1l) {
        					if (!el.passMVA90ID_noIso()) continue;
        				} else {
        					if (!el.passMedID_noIso()) continue;
        				}
        			} else {
        				if (el.ttHMVA() < 0.8) continue;
        			}
        		}
        		els.push_back(&el);
        	}
        	for (const auto& mu : reader_muon->muons) {
        		if (is1l) {
        			if (!passMuSelNoID1(&mu)) continue;
        		} else {
        			if (!passMuSelNoID2(&mu)) continue;
        		}
        		if(doID) {
        			if (doStd) {
        				if (is1l) {
        					if (!mu.passMedID()) continue;
        				} else {
        					if (!mu.passLooseID()) continue;
        				}
        			} else {
        				if (mu.ttHMVA() < 0.85) continue;
        			}
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

    	vector<const Lepton*> singleleps = getLeptons(true,true);
    	vector<const Lepton*> dileps = getLeptons(false,true);

    	bool goodDilepPt = false;
    	if (dileps.size() >= 2) {
    		goodDilepPt = (dileps[0]->pt() > (dileps[0]->isMuon() ? 27 : 30))
    				|| (dileps[1]->pt() > (dileps[1]->isMuon() ? 27 : 30));
    	}

    	TString chS;
    	if (dileps.size() == 2 && goodDilepPt) chS = "2l";
    	else if (dileps.size() < 2 && singleleps.size()) chS = "1l";
    	else return;

    	plotter.getOrMake1DPre(sn+"_"+chS,"ht",";ht",1000,0,4000)->Fill(ht,weight);
    	if (!doStd) {
    		if (chS == "1l") {
        		TString lS = (singleleps[0]->isMuon() ? "m1" : "e1");
            	plotter.getOrMake1DPre(sn+"_"+chS+"_"+lS,"tthMVA",";",100,-1,1)->Fill(singleleps[0]->ttHMVA(),weight);
    		} else {
        		TString lS1 = (dileps[0]->isMuon() ? "m1" : "e1");
        		TString lS2 = (dileps[1]->isMuon() ? "m2" : "e2");
            	plotter.getOrMake1DPre(sn+"_"+chS+"_"+lS1,"tthMVA",";",100,-1,1)->Fill(dileps[0]->ttHMVA(),weight);
            	plotter.getOrMake1DPre(sn+"_"+chS+"_"+lS2,"tthMVA",";",100,-1,1)->Fill(dileps[1]->ttHMVA(),weight);
    		}
    	}
    }

    void validateBkg(TString sn) {
    	if (lepChan == NOCHANNEL) return;
    	TString chS;
    	if (lepChan == SINGLELEP) {
    		chS = "1l";
    		TString lS = (selectedLepton->isMuon() ? "m1" : "e1");
        	plotter.getOrMake1DPre(sn+"_"+chS+"_"+lS,"tthMVA",";",100,-1,1)->Fill(selectedLepton->ttHMVA(),weight);
    	} else if (lepChan == DILEP) {
    		chS = "2l";
    		TString lS1 = (dilep1->isMuon() ? "m1" : "e1");
    		TString lS2 = (dilep2->isMuon() ? "m2" : "e2");
        	plotter.getOrMake1DPre(sn+"_"+chS+"_"+lS1,"tthMVA",";",100,-1,1)->Fill(dilep1->ttHMVA(),weight);
        	plotter.getOrMake1DPre(sn+"_"+chS+"_"+lS2,"tthMVA",";",100,-1,1)->Fill(dilep2->ttHMVA(),weight);
    	}
    	plotter.getOrMake1DPre(sn+"_"+chS,"ht",";ht",1000,0,4000)->Fill(ht,weight);
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
        	validateBkg(sn+"_true");
        }

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getTTHMVAID(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
