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
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Variables/interface/Hww2lSolver.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/DileptonSelection.h"
#include "Processors/Variables/interface/FatJetSelection.h"

#include "TSystem.h"
using namespace TAna;
using namespace std;

typedef float (Lepton::*lepFunFloat)() const;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    	//        turnOffCorr(CORR_TRIG);
    	//        turnOffCorr(CORR_PU  );
    	        turnOffCorr(CORR_LEP );
    	        turnOffCorr(CORR_SJBTAG);
    	        turnOffCorr(CORR_AK4BTAG);
    	//        turnOffCorr(CORR_SDMASS);
    	        turnOffCorr(CORR_TOPPT);
    	        turnOffCorr(CORR_JER);
    }

    void plotSpectra(TString sn, const Lepton* recolep1, const Lepton* recolep2, const FatJet* hbb) {
    	if (!hbb) return;
		if (BTagging::getCSVSJCat(parameters.jets,hbb->subJets()) < BTagging::CSVSJ_MF) return;
		auto bbmass = isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(hbb) : hbb->sdMom().mass();
		if (bbmass < 30 || bbmass > 210) return;

    	MomentumF dilepMOM = recolep1->p4() + recolep2->p4();
    	double dR_ll = PhysicsUtilities::deltaR(*recolep1,*recolep2);
    	double dPhi_ll = PhysicsUtilities::deltaPhi(*recolep1,*recolep2);

    	MomentumF hww = Hww2lSolver::getSimpleHiggsMom(dilepMOM,reader_event->met,40);
    	double mhh = (hbb->p4() + hww.p4()).mass();
    	if (mhh < 700) return;

    	auto plt = [&](TString pre) {
        	plotter.getOrMake1DPre(pre,"Mll",";m_{ll}",100,0,200)->Fill(dilepMOM.mass(),weight);
        	plotter.getOrMake1DPre(pre,"dR_ll",";#DeltaR_{ll}",50,0,5)->Fill(dR_ll,weight);
        	plotter.getOrMake1DPre(pre,"dPhi_ll",";#Delta#Phi_{ll}",50,-3.14,3.14)->Fill(dPhi_ll,weight);
        	plotter.getOrMake1DPre(pre,"pt1",";p_{T} lep1",100,0,1000)->Fill(recolep1->pt(),weight);
        	plotter.getOrMake1DPre(pre,"pt2",";p_{T} lep2",100,0,1000)->Fill(recolep2->pt(),weight);
        	plotter.getOrMake1DPre(pre,"ht",";H_{T}",5000,0,5000)->Fill(ht_puppi,weight);
        	plotter.getOrMake1DPre(pre,"mhh",";M_{HH}",3000,0,3000)->Fill(mhh,weight);
        	plotter.getOrMake1DPre(pre,"met",";MET",100,0,1000)->Fill(reader_event->met.pt(),weight);

        	plotter.getOrMake1DPre(pre,"sip1",";sip",500,0,100)->Fill(recolep1->sip3D(),weight);
        	plotter.getOrMake1DPre(pre,"sip2",";sip",500,0,100)->Fill(recolep2->sip3D(),weight);
    	};

    	plt(sn);
//    	plt("bkg");

    	if (dR_ll > 1.6) return;
    	if (reader_event->met.pt() < 40) return;
    	if (dilepMOM.mass() < 12 || dilepMOM.mass() > 75) return;

        const auto jetsVeto = PhysicsUtilities::selObjsD(jets,
                [&](const Jet* j){return  PhysicsUtilities::deltaR2(*j,*hbb ) >= 1.2*1.2;});
        const auto nBtags = PhysicsUtilities::selObjsMomD(jetsVeto,
                parameters.jets.minBtagJetPT,parameters.jets.maxBTagJetETA,
                [&](const Jet* j){return BTagging::passJetBTagWP(parameters.jets,*j);} ).size();

    	if (nBtags > 0) return;
    	if (PhysicsUtilities::absDeltaPhi(dilepMOM,reader_event->met.p4()) > TMath::PiOver2()) return;

    	plt(sn+"_pass2lSel");
    }
    void printDebugInfo(TString sn, const GenParticle* genlep1, const GenParticle* genlep2, int idx1, int idx2) {
    	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
		std::cout << sn << std::endl;
		printf("gen1 = %i --> idx1 = %i; gen2 = %i --> idx2 = %i\n",genlep1->pdgId(),idx1,genlep2->pdgId(),idx2);
    	for (const auto& mu : reader_muon->muons) {
    		printf("muon %d (%i): (E= %f pT= %f eta= %f phi=  %f)\n",mu.index(),mu.q()*(-13),mu.E(),mu.pt(),mu.eta(),mu.phi());
    	}
    	for (const auto& el : reader_electron->electrons) {
    		printf("electron %d (%i): (E= %f pT= %f eta= %f phi=  %f)\n",el.index(),el.q()*(-11),el.E(),el.pt(),el.eta(),el.phi());
    	}
		printf("\n");
    }

	TString getDilepChan(const Lepton* lep1, const Lepton* lep2) {
		if (lep1->isMuon() && lep2->isMuon()) return "_mumu_";
		else if (lep1->isElectron() && lep2->isElectron()) return "_ee_";
		else return "_emu_";
	}

	bool passSel(int relax, const Lepton* lep1, const Lepton* lep2) {
		bool pass = false;
		bool passIP1 = fabs(lep1->d0()) < 0.05 && fabs(lep1->dz()) < 0.1 && fabs(lep1->sip3D()) < 9999.0;
		bool passIP2 = fabs(lep2->d0()) < 0.05 && fabs(lep2->dz()) < 0.1 && fabs(lep2->sip3D()) < 9999.0;
		bool passIso1 = lep1->miniIso() < 0.15;
		bool passIso2 = lep2->miniIso() < 0.15;

		bool passID1 = lep1->isElectron() ? ((Electron*)lep1)->passMedID_noIso() : ((Muon*)lep1)->passLooseID();
		bool passID2 = lep2->isElectron() ? ((Electron*)lep2)->passMedID_noIso() : ((Muon*)lep2)->passLooseID();
		bool passID = passID1 && passID2;

//		if (lep1->isMuon() && lep2->isMuon() && passID == true) {
//			bool pass1 = (((Muon*)lep1)->passMedID() && lep1->pt() > 27);
//			bool pass2 = (((Muon*)lep2)->passMedID() && lep2->pt() > 27);
//
//			if (!pass1 && !pass2) passID = false;
//		}

		// 0 = relax IP; 1 = relax ID; 2 = relax ISO ////////
		if (relax==0)      pass = passID && passIso1 && passIso2;
		else if (relax==1) pass = passIP1 && passIP2 && passIso1 && passIso2;
		else if (relax==2) pass = passIP1 && passIP2 && passID;
		else {cout<<"passSel arg needs to be 0, 1, or 2"<<endl;}

		return pass;
	}

	void testISO(TString sn, bool isSignal, const Lepton *sigLep1=0, const Lepton *sigLep2=0, const FatJet* hbbJet=0) {

    	vector<double> isoWPs = {0.1, 0.15, 0.2, 0.25, 0.3, 0.4};
    	static const vector<LeptonProcessor::elFunFloat> elIsos = {&Electron::miniIso, &Electron::trackerIso, &Electron::pfIso};
    	static const vector<LeptonProcessor::muFunFloat> muIsos = {&Muon::miniIso, &Muon::trackerIso, &Muon::pfIso};

    	vector<TString> isos = {"miniIso","trkIso","pfIso"};
    	DileptonParameters param = parameters.dileptons;

		TString eVal, mVal;
		TString id = "passIP_ID_eMVA_mL_inclIso";

		if (isSignal) {
	    	if (!passSel(2,sigLep1,sigLep2)) return;
			plotSpectra(sn+"_passIP_ID_eMVA_mL_inclIso",sigLep1,sigLep2,hbbJet);
			plotSpectra(sn+getDilepChan(sigLep1,sigLep2)+"passIP_ID_eMVA_mL_inclIso",sigLep1,sigLep2,hbbJet);
		} else {
			param.mu_getISO = &Muon::inclIso;
			param.el_getISO = &Electron::inclIso;
    		const auto leps = DileptonProcessor::getLeptons(param,*reader_muon,*reader_electron);
    		if (leps.size() == 2) {
    			const Lepton *lep1 = leps.front();
    			const Lepton *lep2 = leps.back();
	    		float pt = lep1->isMuon() ? 27 : 30;
	    		if (lep1->pt() < pt) return;

    			fjProc.reset(new FatJetProcessor ());
    			fjProc->loadDilepFatJet(parameters.fatJets,*reader_fatjet,lep1,lep2);
    			const FatJet *bbjet = fjProc->getDilepHbbCand();

    			if (bbjet) {
    				plotSpectra(sn+"_"+id,lep1,lep2,bbjet);
    				plotSpectra(sn+getDilepChan(lep1,lep2)+id,lep1,lep2,bbjet);
    			}
    		}
		}

    	for(unsigned int iS = 0; iS < isos.size(); iS++) {
    		for(unsigned int im=0; im<isoWPs.size(); im++) for(unsigned int ie=0; ie<isoWPs.size(); ie++) {

	    		eVal = TString::Format("%.2f",isoWPs[ie]); eVal.ReplaceAll(".","p");
	    		mVal = TString::Format("%.2f",isoWPs[im]); mVal.ReplaceAll(".","p");
	    		id = "passIP_ID_eMVA_mL_"+isos[iS]+"_e"+eVal+"_m"+mVal;

    			if (isSignal) {

    				bool passIso1 = false;
    				bool passIso2 = false;

					passIso1 = sigLep1->isMuon() ? (((Muon*)sigLep1)->*muIsos[iS])() < isoWPs[im] : (((Electron*)sigLep1)->*elIsos[iS])() < isoWPs[ie];
					passIso2 = sigLep2->isMuon() ? (((Muon*)sigLep2)->*muIsos[iS])() < isoWPs[im] : (((Electron*)sigLep2)->*elIsos[iS])() < isoWPs[ie];

    				if (passIso1) plotter.getOrMake1DPre(sn+getDilepChan(sigLep1,sigLep2)+id,"ht_lep1",";H_{T}",5000,0,5000)->Fill(ht_puppi,weight);
    				if (passIso2) plotter.getOrMake1DPre(sn+getDilepChan(sigLep1,sigLep2)+id,"ht_lep2",";H_{T}",5000,0,5000)->Fill(ht_puppi,weight);
    				if (passIso1 && passIso2) {
    					plotSpectra(sn+"_"+id,sigLep1,sigLep2,hbbJet);
    					plotSpectra(sn+getDilepChan(sigLep1,sigLep2)+id,sigLep1,sigLep2,hbbJet);
    				}

    			} else {
    				param.el_maxISO = isoWPs[ie];
    				param.mu_maxISO = isoWPs[im];
    				param.mu_getISO = muIsos[iS];
    				param.el_getISO = elIsos[iS];

    	    		const auto leps = DileptonProcessor::getLeptons(param,*reader_muon,*reader_electron);
        	    	if (leps.size() >  2) plotter.getOrMake1DPre(sn+"_"+id+"_LepsGt2","ht",";H_{T}",5000,0,5000)->Fill(ht_puppi,weight);
    	    		if (leps.size() != 2) continue;
    	    		const Lepton *lep1 = leps.front();
    	    		const Lepton *lep2 = leps.back();
    	    		float pt = lep1->isMuon() ? 27 : 30;
    	    		if (lep1->pt() < pt) continue;

    	    	    fjProc.reset(new FatJetProcessor ());
    	            fjProc->loadDilepFatJet(parameters.fatJets,*reader_fatjet,lep1,lep2);
    	            const FatJet *bbjet = fjProc->getDilepHbbCand();

        			if (bbjet) {
        				plotSpectra(sn+"_"+id,lep1,lep2,bbjet);
        				plotSpectra(sn+getDilepChan(lep1,lep2)+id,lep1,lep2,bbjet);
        			}
    			}
    		}
    	}
	}

	void testID(TString sn, bool isSignal, const Lepton *sigLep1=0, const Lepton *sigLep2=0, const FatJet* hbbJet=0) {

    	static const vector<LeptonProcessor::elFunBool> elIds = {&Electron::passLooseID_noIso, &Electron::passMedID_noIso, &Electron::passTightID_noIso, &Electron::passMVA90ID_noIso};
    	static const vector<LeptonProcessor::muFunBool> muIds = {&Muon::passLooseID, &Muon::passMedID, &Muon::passTightID};
    	vector<TString> elNames = {"L","M","T","MVA"};
    	vector<TString> muNames = {"L","M","T"};

    	DileptonParameters param = parameters.dileptons;
    	TString eid1, eid2, mid1, mid2;
    	TString id = "passIP_ID_incl_miniIso0p2";

		if (isSignal) {
	    	if (!passSel(1,sigLep1,sigLep2)) return;
			plotSpectra(sn+"_"+id,sigLep1,sigLep2,hbbJet);
			plotSpectra(sn+getDilepChan(sigLep1,sigLep2)+id,sigLep1,sigLep2,hbbJet);
		} else {
			param.el_getID1 = &Electron::passInclID;
			param.el_getID2 = &Electron::passInclID;
			param.mu_getID1 = &Muon::passInclID;
			param.mu_getID2 = &Muon::passInclID;

	    	const auto leps = DileptonProcessor::getLeptons(param,*reader_muon,*reader_electron);
	    	if (leps.size() == 2) {
	    		const Lepton *lep1 = leps.front();
	    		const Lepton *lep2 = leps.back();
	    		float pt = lep1->isMuon() ? 27 : 30;
	    		if (lep1->pt() < pt) return;

	    		fjProc.reset(new FatJetProcessor ());
	    		fjProc->loadDilepFatJet(parameters.fatJets,*reader_fatjet,lep1,lep2);
	    		const FatJet *bbjet = fjProc->getDilepHbbCand();

    			if (bbjet) {
    				plotSpectra(sn+"_"+id,lep1,lep2,bbjet);
    				plotSpectra(sn+getDilepChan(lep1,lep2)+id,lep1,lep2,bbjet);
    			}
	    	}
		}

    	for(unsigned int im1=0; im1<muIds.size(); im1++) for(unsigned int ie1=0; ie1<elIds.size(); ie1++)
    		for(unsigned int im2=0; im2<muIds.size(); im2++) for(unsigned int ie2=0; ie2<elIds.size(); ie2++) {

	    	eid1 = elNames[ie1];
	    	mid1 = muNames[im1];
	    	eid2 = elNames[ie2];
	    	mid2 = muNames[im2];
	    	id = "passIP_ID_e1_"+eid1+"_e2_"+eid2+"_m1_"+mid1+"_m2_"+mid2+"_miniIso0p2";

    		if (isSignal) {

    			bool passId1 = false;
    			bool passId2 = false;
    			passId1 = sigLep1->isMuon() ? (((Muon*)sigLep1)->*muIds[im1])() : (((Electron*)sigLep1)->*elIds[ie1])();
    			passId2 = sigLep2->isMuon() ? (((Muon*)sigLep2)->*muIds[im2])() : (((Electron*)sigLep2)->*elIds[ie2])();

    			if (passId1) plotter.getOrMake1DPre(sn+getDilepChan(sigLep1,sigLep2)+id,"ht_lep1",";H_{T}",5000,0,5000)->Fill(ht_puppi,weight);
    			if (passId2) plotter.getOrMake1DPre(sn+getDilepChan(sigLep1,sigLep2)+id,"ht_lep2",";H_{T}",5000,0,5000)->Fill(ht_puppi,weight);
    			if (passId1 && passId2) {
    				plotSpectra(sn+"_"+id,sigLep1,sigLep2,hbbJet);
    				plotSpectra(sn+getDilepChan(sigLep1,sigLep2)+id,sigLep1,sigLep2,hbbJet);
    			}

    		} else {
    			param.el_getID1 = elIds[ie1];
    			param.el_getID2 = elIds[ie2];
    			param.mu_getID1 = muIds[im1];
    			param.mu_getID2 = muIds[im2];

    	    	const auto leps = DileptonProcessor::getLeptons(param,*reader_muon,*reader_electron);
    	    	if (leps.size() >  2) plotter.getOrMake1DPre(sn+"_"+id+"_LepsGt2","ht",";H_{T}",5000,0,5000)->Fill(ht_puppi,weight);
    	    	if (leps.size() != 2) continue;
    	    	const Lepton *lep1 = leps.front();
    	    	const Lepton *lep2 = leps.back();
	    		float pt = lep1->isMuon() ? 27 : 30;
	    		if (lep1->pt() < pt) continue;

    	    	bool goodSel = true;
    	    	if (lep1->isMuon()) {
    	    		if (!(((Muon*)lep1)->*muIds[im1])()) goodSel = false;
    	    	} else {
    	    		if (!(((Electron*)lep1)->*elIds[ie1])()) goodSel = false;
    	    	}
    	    	if (lep2->isMuon()) {
    	    		if (!(((Muon*)lep2)->*muIds[im2])()) goodSel = false;
    	    	} else {
    	    		if (!(((Electron*)lep2)->*elIds[ie2])()) goodSel = false;
    	    	}
    	    	if (!goodSel) continue;

    	    	fjProc.reset(new FatJetProcessor ());
    	        fjProc->loadDilepFatJet(parameters.fatJets,*reader_fatjet,lep1,lep2);
    	        const FatJet *bbjet = fjProc->getDilepHbbCand();

    			if (bbjet) {
    				plotSpectra(sn+"_"+id,lep1,lep2,bbjet);
    				plotSpectra(sn+getDilepChan(lep1,lep2)+id,lep1,lep2,bbjet);
    			}
    		}
    	}
	}

	void testSpecialID(TString sn, bool isSignal, const Lepton *sigLep1=0, const Lepton *sigLep2=0, const FatJet* hbbJet=0) {

		bool passId = false;
		TString id = "passIP_ID_e1_MVA_e2_MVA_m1_S_m2_S_miniIso0p2";

		if (isSignal) {
			if (sigLep1->isElectron() || sigLep2->isElectron()) return;
	    	if (!passSel(1,sigLep1,sigLep2)) return;

			if ( ((Muon*)sigLep1)->passLooseID() && ((Muon*)sigLep2)->passLooseID() ) {
				bool pass21 = ((Muon*)sigLep1)->passMedID() && sigLep1->pt() > 27;
				bool pass22 = ((Muon*)sigLep2)->passMedID() && sigLep2->pt() > 27;

				if (pass21 || pass22) plotSpectra(sn+getDilepChan(sigLep1,sigLep2)+id,sigLep1,sigLep2,hbbJet);
			}
		} else {
	    	DileptonParameters param = parameters.dileptons;
			param.mu_getID1 = &Muon::passLooseID;
			param.mu_getID2 = &Muon::passLooseID;

	    	const auto leps = DileptonProcessor::getLeptons(param,*reader_muon,*reader_electron);
	    	if (leps.size() != 2) return;
	    	const Lepton *lep1 = leps.front();
	    	const Lepton *lep2 = leps.back();
    		float pt = lep1->isMuon() ? 27 : 30;
    		if (lep1->pt() < pt) return;

	    	if (lep1->isElectron() || lep2->isElectron()) return;

			bool pass21 = ((Muon*)lep1)->passMedID() && lep1->pt() > 27;
			bool pass22 = ((Muon*)lep2)->passMedID() && lep2->pt() > 27;

			if (!pass21 && !pass22) return;

	    	fjProc.reset(new FatJetProcessor ());
	        fjProc->loadDilepFatJet(parameters.fatJets,*reader_fatjet,lep1,lep2);
	        const FatJet *bbjet = fjProc->getDilepHbbCand();

			if (bbjet) {
				plotSpectra(sn+getDilepChan(lep1,lep2)+id,lep1,lep2,bbjet);
			}

		}
	}

	void testIP(TString sn, bool isSignal, const Lepton *sigLep1=0, const Lepton *sigLep2=0, const FatJet* hbbJet=0) {

		vector<double> dzWPs = {0.05,0.1,0.15,0.2,0.25,0.3,0.4,9999};
		vector<double> d0WPs = {0.05,0.1,0.15,0.2,0.25,0.3,0.4,9999};
		vector<double> sipWPs = {1,2,3,4,5,6,9999};
		static const vector<lepFunFloat> ipVars = {&Lepton::dz, &Lepton::d0, &Lepton::sip3D};
		vector<vector<double>> wps = {dzWPs, d0WPs, sipWPs};
		vector<TString> ipstrs = {"dz","d0","sip"};

    	DileptonParameters param = parameters.dileptons;
    	TString constIdIso = "ID_eMVA_mM_miniIso0p2";
    	TString ipid = "inclIP_";

		if (isSignal) {
	    	if (!passSel(0,sigLep1,sigLep2)) return;
			plotSpectra(sn+"_"+ipid+constIdIso,sigLep1,sigLep2,hbbJet);
			plotSpectra(sn+getDilepChan(sigLep1,sigLep2)+ipid+constIdIso,sigLep1,sigLep2,hbbJet);
		} else {
			param.mu_maxDZ = 9999;
			param.mu_maxD0 = 9999;
			param.mu_maxSip3D = 9999;
			param.el_maxDZ = 9999;
			param.el_maxD0 = 9999;
			param.el_maxSip3D = 9999;

	    	const auto leps = DileptonProcessor::getLeptons(param,*reader_muon,*reader_electron);
	    	if (leps.size() == 2) {
	    		const Lepton *lep1 = leps.front();
	    		const Lepton *lep2 = leps.back();

	    		fjProc.reset(new FatJetProcessor ());
	    		fjProc->loadDilepFatJet(parameters.fatJets,*reader_fatjet,lep1,lep2);
	    		const FatJet *bbjet = fjProc->getDilepHbbCand();

    			if (bbjet) {
    				auto bbmass = isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(bbjet) : bbjet->sdMom().mass();
    				if (bbmass > 30 && bbmass < 210) {
    					plotSpectra(sn+"_"+ipid+constIdIso,lep1,lep2,bbjet);
    					plotSpectra(sn+getDilepChan(lep1,lep2)+ipid+constIdIso,lep1,lep2,bbjet);
    				}
    			}
	    	}
		}
    	for (unsigned int var = 0; var < ipVars.size(); var++) {

			if (isSignal) {
	    		if (!passIPsel(var,sigLep1,sigLep2)) continue;
				plotSpectra(sn+"_incl_"+ipstrs[var]+"_"+constIdIso,
						sigLep1,sigLep2,hbbJet);
				plotSpectra(sn+getDilepChan(sigLep1,sigLep2)+"incl_"+ipstrs[var]+"_"+constIdIso,
						sigLep1,sigLep2,hbbJet);
			}

    		vector<double>& wps_ = wps[var];
    		for (unsigned int wp = 0; wp < wps_.size(); wp++) {

    			ipid = ipstrs[var]+TString::Format("_%.2f_",wps_[wp]); ipid.ReplaceAll(".","p");
    			if (isSignal) {

        			if(std::fabs( ((sigLep1)->*ipVars[var])() ) > wps_[wp]) continue;
        			if(std::fabs( ((sigLep2)->*ipVars[var])() ) > wps_[wp]) continue;

        			plotSpectra(sn+"_"+ipid+constIdIso,sigLep1,sigLep2,hbbJet);
        			plotSpectra(sn+getDilepChan(sigLep1,sigLep2)+ipid+constIdIso,sigLep1,sigLep2,hbbJet);

    			} else {
    				if (var==0) {
    					param.el_maxDZ = wps_[wp];
    					param.mu_maxDZ = wps_[wp];
    				} else if (var == 1) {
    					param.el_maxD0 = wps_[wp];
    					param.mu_maxD0 = wps_[wp];
    				} else if (var == 2) {
    					param.el_maxSip3D = wps_[wp];
    					param.mu_maxSip3D = wps_[wp];
    				} else {cout<<"var should not be anything other than 0, 1, or 2!"<<endl;}

        	    	const auto leps = DileptonProcessor::getLeptons(param,*reader_muon,*reader_electron);
        	    	if (leps.size() != 2) continue;
        	    	const Lepton *lep1 = leps.front();
        	    	const Lepton *lep2 = leps.back();
    	    		float pt = lep1->isMuon() ? 27 : 30;
    	    		if (lep1->pt() < pt) continue;

        	    	fjProc.reset(new FatJetProcessor ());
        	        fjProc->loadDilepFatJet(parameters.fatJets,*reader_fatjet,lep1,lep2);
        	        const FatJet *bbjet = fjProc->getDilepHbbCand();

        			if (bbjet) {
        				auto bbmass = isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(bbjet) : bbjet->sdMom().mass();
        				if (bbmass > 30 && bbmass < 210) {
        					plotSpectra(sn+"_"+ipid+constIdIso,lep1,lep2,bbjet);
        					plotSpectra(sn+getDilepChan(lep1,lep2)+ipid+constIdIso,lep1,lep2,bbjet);
        				}
        			}
    			}
    		}
    	}
	}

	bool passIPsel(int relax, const Lepton* lep1, const Lepton* lep2) {
		bool pass = false;

		bool passDZ1 = fabs(lep1->dz()) < 0.1;
		bool passDZ2 = fabs(lep2->dz()) < 0.1;
		bool passD01 = fabs(lep1->d0()) < 0.05;
		bool passD02 = fabs(lep2->d0()) < 0.05;
		bool passSip1 = fabs(lep1->sip3D()) < 9999.0;
		bool passSip2 = fabs(lep2->sip3D()) < 9999.0;

		if (relax==0)      pass = passD01 && passD02 && passSip1 && passSip2;
		else if (relax==1) pass = passDZ1 && passDZ2 && passSip1 && passSip2;
		else if (relax==2) pass = passDZ1 && passDZ2 && passD01 && passD02;

		return pass;
	}

    const FatJet *getSigHbb(const Lepton *lep1, const Lepton *lep2) {

    	auto isBtag = [&] (const FatJet* fj) {
    		bool hasBtag = false;
    		for (const auto& sj : fj->subJets()) {
    			if (sj.deep_csv() >= parameters.jets.DeepCSV_WP[2]) hasBtag = true;
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

    	auto sepLep = [&] (const FatJet* fj)->bool {
    		bool sep = true;
    		if (PhysicsUtilities::deltaR2(*fj,*lep1) < 0.8*0.8) sep = false;
    		if (PhysicsUtilities::deltaR2(*fj,*lep2) < 0.8*0.8) sep = false;
    		if (PhysicsUtilities::absDeltaPhi(*fj,(lep1->p4()+lep2->p4())) < 2.0) sep = false;
    		return sep;
    	};

    	const FatJet *hbbjet = getMatchedFJ(diHiggsEvt.hbb->p4(),true,reader_fatjet->jets,reader_fatjet_noLep->jets);

    	if(!hbbjet) return 0;
    	if(!hasGoodSJs(hbbjet) || !isBtag(hbbjet)) return 0;
    	if(!sepLep(hbbjet)) return 0;

    	return hbbjet;
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

    	const GenParticle *ptcl = gp;
    	return ptcl;
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht_puppi < 400) return false;
        TString sn = smpName;
        if (mcProc == FillerConstants::TTBAR && smDecayEvt.nLepsTT >= 0 && smDecayEvt.nLepsTT <= 2) {
        	sn += TString::Format("%d",smDecayEvt.nLepsTT);
        }


        // SIGNAL
        if(isSignal() && diHiggsEvt.type == DiHiggsEvent::DILEP) {

			plotter.getOrMake1DPre(sn+"_full_signal","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			const GenParticle* glep1_0 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1;
	        const GenParticle* glep2_0 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1;


			if (glep1_0->pdgId()<0 == glep2_0->pdgId()<0) return false;

			const GenParticle *glep1=0;
			const GenParticle *glep2=0;

			if (glep1_0->absPdgId() == 15) glep1 = getGenLepFromTau(glep1_0);
			else glep1 = glep1_0;

			if (glep2_0->absPdgId() == 15) glep2 = getGenLepFromTau(glep2_0);
			else glep2 = glep2_0;

			if (!glep1 || !glep2) return false;

	        TString llChan = "_emu";
	        if (glep1->absPdgId()==13 && glep2->absPdgId()==13) llChan = "_mumu";
	        else if (glep1->absPdgId()==11 && glep2->absPdgId()==11) llChan = "_ee";

			plotter.getOrMake1DPre(sn+"_gen_dilep","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			plotter.getOrMake1DPre(sn+llChan+"_gen_dilep","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			// GEN lepton pt cuts
            bool passPt1 = glep1->absPdgId() == 13 ? (glep1->pt() > 27) : (glep1->pt() > 30);
            bool passPt2 = (glep2->pt() > 10);
            if (!(passPt1 && passPt2)) return false;

			plotter.getOrMake1DPre(sn+"_gen_dilep_passPt","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			plotter.getOrMake1DPre(sn+llChan+"_gen_dilep_passPt","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

            const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,1.479);
            const Lepton* matchLep1 = getMatchedLepton(glep1,muons,electrons,0.1,true);
            const Lepton* matchLep2 = getMatchedLepton(glep2,muons,electrons,0.1,true);

            if (!matchLep1 || !matchLep2) return false;
			plotter.getOrMake1DPre(sn+"_foundLepMatch","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			if ((matchLep1->isMuon() == matchLep2->isMuon()) && (matchLep1->index() == matchLep2->index())) return false; // discard if these are the same RECO lep
			plotter.getOrMake1DPre(sn+"_goodLepMatch","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			plotter.getOrMake1DPre(sn+llChan+"_goodLepMatch","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			// RECO lepton pt cuts
			bool passRecoPt1 = matchLep1->isMuon() ? (matchLep1->pt() > 27) : (matchLep1->pt() > 30);
			bool passRecoPt2 = (matchLep2->pt() > 10);
			if (!passRecoPt1 || !passRecoPt2) return false;
			plotter.getOrMake1DPre(sn+"_passRecoPt","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
			plotter.getOrMake1DPre(sn+llChan+"_passRecoPt","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			const FatJet* matchHbb = getSigHbb(matchLep1,matchLep2);
			if (!matchHbb) return false;
			plotter.getOrMake1DPre(sn+"_validHbb","evts",";M_{X}",5000,0,5000)->Fill(signal_mass,weight);
			plotter.getOrMake1DPre(sn+getDilepChan(matchLep1,matchLep2)+"validHbb","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			plotSpectra(sn+"_base",matchLep1,matchLep2,matchHbb);
			plotSpectra(sn+getDilepChan(matchLep1,matchLep2)+"base",matchLep1,matchLep2,matchHbb);

//			testIP(sn,isSignal(),matchLep1,matchLep2,matchHbb);
			testID(sn,isSignal(),matchLep1,matchLep2,matchHbb);
			testISO(sn,isSignal(),matchLep1,matchLep2,matchHbb);

			testSpecialID(sn,isSignal(),matchLep1,matchLep2,matchHbb);


        }
        // BKG
        if (!isSignal()) {
        	plotter.getOrMake1DPre(sn+"_baseline","evts",";H_{T}",3000,0,3000)->Fill(ht_puppi,weight);

//			testIP (sn,isSignal(),0,0,0);
			testID (sn,isSignal(),0,0,0);
			testISO(sn,isSignal(),0,0,0);
			testSpecialID(sn,isSignal(),0,0,0);
        }
        return true;
    }

    std::unique_ptr<FatJetProcessor> fjProc;
    void write(TString fileName){
    	plotter.write(fileName);
    }
    HistGetter plotter;
};

#endif

void getDileptonSelection(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
