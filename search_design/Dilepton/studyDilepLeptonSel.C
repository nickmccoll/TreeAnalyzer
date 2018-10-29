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
    struct RecoLepSelection {
    	std::vector<const Lepton*> selectedLeps;
    	std::vector<const Lepton*> filteredLeps;
    };

    void plotSpectra(TString sn, const Lepton* recolep1, const Lepton* recolep2, const FatJet* hbbjet) {
    	MomentumF dilepMOM = recolep1->p4() + recolep2->p4();
    	MomentumF bbllMOM = dilepMOM.p4() + hbbjet->p4();
    	double dR_ll = PhysicsUtilities::deltaR(*recolep1,*recolep2);
    	double dPhi_ll = PhysicsUtilities::deltaPhi(*recolep1,*recolep2);
    	double dR_bbll = PhysicsUtilities::deltaR(dilepMOM,*hbbjet);

    	plotter.getOrMake1DPre(sn,"Mll",";m_{ll}",100,0,200)->Fill(dilepMOM.mass(),weight);
    	plotter.getOrMake1DPre(sn,"Mbbll",";m_{bbll}",100,0,200)->Fill(bbllMOM.mass(),weight);
    	plotter.getOrMake1DPre(sn,"dR_ll",";#DeltaR_{ll}",50,0,5)->Fill(dR_ll,weight);
    	plotter.getOrMake1DPre(sn,"dPhi_ll",";#Delta#Phi_{ll}",50,-3.14,3.14)->Fill(dPhi_ll,weight);
    	plotter.getOrMake1DPre(sn,"pt1",";p_{T} lep1",100,0,1000)->Fill(recolep1->pt(),weight);
    	plotter.getOrMake1DPre(sn,"pt2",";p_{T} lep2",100,0,1000)->Fill(recolep2->pt(),weight);
    	plotter.getOrMake1DPre(sn,"ptbb",";p_{T} bb",100,0,1000)->Fill(hbbjet->pt(),weight);
    	plotter.getOrMake1DPre(sn,"dR_bbll",";#DeltaR_{bb,ll}",100,0,5)->Fill(dR_bbll,weight);
    }
    void printDebugInfo(TString sn, const GenParticle* genlep1, const GenParticle* genlep2, int idx1, int idx2, std::vector<const Lepton*> leps) {
    	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
		std::cout << sn << std::endl;
		printf("gen1 = %i --> idx1 = %i; gen2 = %i --> idx2 = %i\n",genlep1->pdgId(),idx1,genlep2->pdgId(),idx2);
    	for (const auto& lep : leps) {
    		printf("lepton %d (%i): (E= %f pT=   %f eta= %f phi=  %f)\n",lep->index(),lep->isMuon() ? (-1)*lep->q()*13:(-1)*lep->q()*11,lep->E(),lep->pt(),lep->eta(),lep->phi());
    	}
		printf("\n");
    }
    const FatJet* findHbbCand(const Lepton* lep1, const Lepton* lep2) {
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
    			if (sj.pt() > 20 && sj.eta() < 2.4) nGoodSJ++;
    		}
    		if (nGoodSJ > 2) goodSJs = true;
    		return goodSJs;
    	};
        // only consider the top two fatjets in pt, provided they are above 200 GeV, then take furthest
    	double minDPhi = 2.0;
    	double minPt = 200;
    	int idx = -1;

        const MomentumF recodilepton = lep1->p4() + lep2->p4();
        std::vector<const FatJet*> fatjets;
        for (const auto& fj : reader_fatjet->jets) {
        	if (fj.pt() > minPt) fatjets.push_back(&fj);
        }
        std::sort(fatjets.begin(),fatjets.end(), PhysicsUtilities::greaterPTDeref<FatJet>());

        if (fatjets.size() == 1) {
        	bool separatedFJ = abs(PhysicsUtilities::deltaPhi(*fatjets[0],recodilepton)) > 2.0 && PhysicsUtilities::deltaR(*fatjets[0],*lep1) > 0.8
        			&& PhysicsUtilities::deltaR(*fatjets[0],*lep2) > 0.8;
        	if (separatedFJ && isBtag(fatjets[0]) && hasGoodSJs(fatjets[0])) idx = fatjets[0]->index();

        } else if (fatjets.size() > 1) {
            double fj_dr = 0;
            for (int k=0; k<2; k++) {
            	bool separatedFJ = abs(PhysicsUtilities::deltaPhi(*fatjets[k],recodilepton)) > 2.0 && PhysicsUtilities::deltaR(*fatjets[k],*lep1) > 0.8
                	    && PhysicsUtilities::deltaR(*fatjets[k],*lep2) > 0.8;

                if (separatedFJ && isBtag(fatjets[k]) && hasGoodSJs(fatjets[k])) {
                    if (PhysicsUtilities::deltaR(*fatjets[k],recodilepton) > fj_dr) {
                	    fj_dr = PhysicsUtilities::deltaR(*fatjets[k],recodilepton);
                	    idx = fatjets[k]->index();
                	}
                }
            }
        }
        if (idx < 0) return 0;
        return &reader_fatjet->jets[idx];
    }
    bool matchDileptons(const Lepton* matchReco1, const Lepton* matchReco2, const Lepton* recoCand1, const Lepton* recoCand2) {

    	// check if the two selected reco leptons are matches to the DR-matched leptons (degeneracy = 2)
    	if (matchReco1 && matchReco2 && recoCand1 && recoCand2) {
			if (matchReco1->isMuon() == recoCand1->isMuon() && matchReco2->isMuon() == recoCand2->isMuon()) {
				if (matchReco1->index() == recoCand1->index() && matchReco2->index() == recoCand2->index()) return true;
			} else if (matchReco1->isMuon() == recoCand2->isMuon() && matchReco2->isMuon() == recoCand1->isMuon()) {
				if (matchReco1->index() == recoCand2->index() && matchReco2->index() == recoCand1->index()) return true;
			}
    	}
    	return false;
    }
    const Lepton * getMatchedLepton(const GenParticle& genLepton,const std::vector<const Muon *> muons, const std::vector<const Electron*> electrons){
       if(genLepton.absPdgId() == ParticleInfo::p_muminus){
           double nearestDR =10;
           int idx = PhysicsUtilities::findNearestDRDeref(genLepton,muons,nearestDR,0.2);
           if(idx < 0) return 0;
           else return muons[idx];
       } else {
           double nearestDR =10;
           int idx = PhysicsUtilities::findNearestDRDeref(genLepton,electrons,nearestDR,0.2);
           if(idx < 0) return 0;
           else return electrons[idx];
       }
    }
    Analyzer::RecoLepSelection selectRecoLeptons(const std::vector<const Muon*> muons, const std::vector<const Electron*> electrons) {
    	RecoLepSelection LepInfo;

    	// selection on lepton1
    	float minPt1_mu = 26;
    	float minPt1_e = 30;

    	// collect the muons and electrons together and then sort the leptons by pt
    	std::vector<const Lepton*> leps;
    	for (const auto* mu : muons) leps.push_back(mu);
    	for (const auto* e : electrons) leps.push_back(e);
        std::sort(leps.begin(),leps.end(), PhysicsUtilities::greaterPTDeref<Lepton>());

        LepInfo.filteredLeps = leps;
        std::vector<const Lepton*> selected_leps;
    	if (leps.size() < 2) return LepInfo;

        // take highest pt lepton as the primary lepton, provided it passes pt cut
        if (leps[0]->isMuon()) {
        	if (leps[0]->pt() < minPt1_mu) return LepInfo;
        	else selected_leps.push_back(leps[0]);
        } else if (leps[0]->isElectron()) {
        	if (leps[0]->pt() < minPt1_e) return LepInfo;
        	else selected_leps.push_back(leps[0]);
        } else printf("Error: lep not identified as e or mu\n");

        // look for next highest pt lepton with opposite charge
        int lep2_index = -1;
        for (unsigned long k=1; k<leps.size(); k++){
        	if (leps[k]->q() == leps[0]->q()) continue; // disregard same-sign leptons
            lep2_index = k;
            break;
        }
        if (lep2_index < 1) return LepInfo;
        else {
        	selected_leps.push_back(leps[lep2_index]);
        	LepInfo.selectedLeps = selected_leps;
        	return LepInfo;
        }
    }
    bool passIP(const Lepton* lep1, const Lepton* lep2) {
    	double maxSIP = 4.0;
    	double maxD0 = 0.05;
    	double maxDZ = 0.1;
    	if (lep1->sip3D() > maxSIP || lep2->sip3D() > maxSIP) return false;
    	if (lep1->d0() > maxD0 || lep2->d0() > maxD0) return false;
    	if (lep1->dz() > maxDZ || lep2->dz() > maxDZ) return false;
    	return true;
    }
    void testISO(TString sn, const Lepton* lep1, const Lepton* lep2, const FatJet* hbbjet) {
    	bool passID_1 = lep1->isMuon() ? ((const Muon*)lep1)->passMed16ID() : ((const Electron*)lep1)->passMedID_noISO();
    	bool passID_2 = lep2->isMuon() ? ((const Muon*)lep2)->passMed16ID() : ((const Electron*)lep2)->passMedID_noISO();

    	if (passID_1 && passID_2) {
    		std::vector<double> isoVals = {0.1,0.2,0.3};
    		std::vector<TString> isoTypes = {"miniIso_","relIso_"};
    		sn += "ID_eT_muM_";

    		for (unsigned long type=0; type<isoTypes.size(); type++) {
    			if (type == 0) {
    				for (const auto& e_iso : isoVals) {
    					for (const auto& mu_iso : isoVals) {
    						bool passIso1 = lep1->isMuon() ? lep1->miniIso() < mu_iso : lep1->miniIso() < e_iso;
    						bool passIso2 = lep2->isMuon() ? lep2->miniIso() < mu_iso : lep2->miniIso() < e_iso;

							TString name = isoTypes[type]+TString::Format("e%.1f_mu%.1f_",e_iso,mu_iso);
							name.ReplaceAll(".","");
							if (passIso1) plotter.getOrMake1DPre(sn+name,"evts_lep1",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
							if (passIso2) plotter.getOrMake1DPre(sn+name,"evts_lep2",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
    						if (passIso1 && passIso2) {
    							plotSpectra(sn+name,lep1,lep2,hbbjet);
    							plotter.getOrMake1DPre(sn+name,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

    						}
    					}
    				}
				} else if (type == 1) {
    				for (const auto& e_iso : isoVals) {
    					for (const auto& mu_iso : isoVals) {
    						bool passIso1 = lep1->isMuon() ? ((const Muon*)lep1)->dbRelISO() < mu_iso : ((const Electron*)lep1)->eaRelISO() < e_iso;
    						bool passIso2 = lep2->isMuon() ? ((const Muon*)lep2)->dbRelISO() < mu_iso : ((const Electron*)lep2)->eaRelISO() < e_iso;

							TString name = isoTypes[type]+TString::Format("e%.1f_mu%.1f_",e_iso,mu_iso);
							name.ReplaceAll(".","");
							if (passIso1) plotter.getOrMake1DPre(sn+name,"evts_lep1",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
							if (passIso2) plotter.getOrMake1DPre(sn+name,"evts_lep2",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
    						if (passIso1 && passIso2) {
    							plotSpectra(sn+name,lep1,lep2,hbbjet);
    							plotter.getOrMake1DPre(sn+name,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
    						}
    					}
    				}
				}
    		}
    	}
    }
    void testID(TString sn, const Lepton* lep1, const Lepton* lep2, const FatJet* hbbjet) {
    	bool passISO_1 = lep1->isMuon() ? lep1->miniIso() < 0.2 : lep1->miniIso() < 0.1;
    	bool passISO_2 = lep2->isMuon() ? lep2->miniIso() < 0.2 : lep2->miniIso() < 0.1;
    	TString suf = "miniIso_e01_mu02_";

    	if (passISO_1 && passISO_2) {
    		std::vector<TString> ids = {"L","M","T","H"};
    		for (unsigned long iMu=0; iMu<ids.size(); iMu++) {
    			for (unsigned long iE=0; iE<ids.size(); iE++) {
					bool passId1 = false;
					bool passId2 = false;

					if (lep1->isMuon()) {
						if (iMu==0) passId1 = ((const Muon*)lep1)->passLooseID();
						else if (iMu==1) passId1 = ((const Muon*)lep1)->passMed16ID();
						else if (iMu==2) passId1 = ((const Muon*)lep1)->passTightID();
						else if (iMu==3) passId1 = ((const Muon*)lep1)->passHighPT();
					} else {
						if (iE==0) passId1 = ((const Electron*)lep1)->passLooseID_noISO();
						else if (iE==1) passId1 = ((const Electron*)lep1)->passMedID_noISO();
						else if (iE==2) passId1 = ((const Electron*)lep1)->passTightID_noISO();
						else if (iE==3) passId1 = ((const Electron*)lep1)->passHEEPID_noISO();
					}
					if (lep2->isMuon()) {
						if (iMu==0) passId2 = ((const Muon*)lep2)->passLooseID();
						else if (iMu==1) passId2 = ((const Muon*)lep2)->passMed16ID();
						else if (iMu==2) passId2 = ((const Muon*)lep2)->passTightID();
						else if (iMu==3) passId2 = ((const Muon*)lep2)->passHighPT();
					} else {
						if (iE==0) passId2 = ((const Electron*)lep2)->passLooseID_noISO();
						else if (iE==1) passId2 = ((const Electron*)lep2)->passMedID_noISO();
						else if (iE==2) passId2 = ((const Electron*)lep2)->passTightID_noISO();
						else if (iE==3) passId2 = ((const Electron*)lep2)->passHEEPID_noISO();
					}

					TString name = "ID_e_"+ids[iE]+"_mu_"+ids[iMu]+"_";
					if (passId1) plotter.getOrMake1DPre(sn+name+suf,"evts_lep1",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
					if (passId2) plotter.getOrMake1DPre(sn+name+suf,"evts_lep2",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
					if (passId1 && passId2) {
						plotSpectra(sn+name+suf,lep1,lep2,hbbjet);
						plotter.getOrMake1DPre(sn+name+suf,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
					}
    			}
    		}
    	}
    }
    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(ht_chs < 400) return false;
        TString sn = smpName;

        // SIGNAL
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type == DiHiggsEvent::DILEP) {

			// order the genleptons by pt, then keep only e and mu events with opposite-sign leptons
			const GenParticle *genlep1 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1);
			const GenParticle *genlep2 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1);
			int lep1id = genlep1->pdgId();
			int lep2id = genlep2->pdgId();

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
			plotter.getOrMake1DPre(sn+"_baseline_","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

			// get reco leptons that pass Preselection maxEta cut and minPt cut
	    	float minPt2_mu = 10;
	    	float minPt2_e = 10;
	    	float maxEta_e = 2.5;
	    	float maxEta_mu = 2.4;
	    	const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,minPt2_mu,maxEta_mu);
	    	const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,minPt2_e,maxEta_e);

			// Select the reco dileptons
        	RecoLepSelection LepInfo = selectRecoLeptons(muons,electrons);
        	std::vector<const Lepton*> selectedLeps = LepInfo.selectedLeps;

        	if (selectedLeps.size() < 2) return false;
        	if (selectedLeps.size() > 2) printf("WARNING: Selected more than 2 leps\n");

        	// Check if the selected RecoLeps are the same as the ones matched to GenParticles
        	const Lepton* match1 = getMatchedLepton(*diHiggsEvt.w1_d1,muons,electrons);
        	const Lepton* match2 = getMatchedLepton(*diHiggsEvt.w2_d1,muons,electrons);
        	if (!matchDileptons(match1, match2, selectedLeps[0], selectedLeps[1])) return false;

        	// Check if a valid Hbb candidate is found in the event
        	const FatJet* hbbjet = findHbbCand(selectedLeps[0],selectedLeps[1]);
        	if (!hbbjet) return false;

        	sn += "_ept10_mupt10_";
			plotter.getOrMake1DPre(sn+"lepmatch_","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        	// filter selectedLeps thru IP cuts
        	if (!passIP(selectedLeps[0],selectedLeps[1])) return false;
        	testISO(sn,selectedLeps[0],selectedLeps[1],hbbjet);
        	testID(sn,selectedLeps[0],selectedLeps[1],hbbjet);
        }
        // BKG
        if (reader_event->process != FillerConstants::SIGNAL) {
			plotter.getOrMake1DPre(sn+"_baseline_","evts",";M_{X}",50,600,4600)->Fill(ht_chs,weight);

			// get reco leptons that pass Preselection maxEta cut and minPt cut
	    	float minPt2_mu = 10;
	    	float minPt2_e = 10;
	    	float maxEta_e = 2.5;
	    	float maxEta_mu = 2.4;
	    	const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,minPt2_mu,maxEta_mu);
	    	const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,minPt2_e,maxEta_e);

			// Select the reco dileptons
        	RecoLepSelection LepInfo = selectRecoLeptons(muons,electrons);
        	std::vector<const Lepton*> selectedLeps = LepInfo.selectedLeps;

        	if (selectedLeps.size() < 2) return false;
        	if (selectedLeps.size() > 2) printf("WARNING: Selected more than 2 leps\n");

        	// Check if a valid Hbb candidate is found in the event
        	const FatJet* hbbjet = findHbbCand(selectedLeps[0],selectedLeps[1]);
        	if (!hbbjet) return false;

			plotter.getOrMake1DPre(sn+"_passTrigPreSel","evts",";M_{X}",50,600,4600)->Fill(ht_chs,weight);
            if (selectedLeps[0]->isMuon() && selectedLeps[1]->isMuon()) sn += "_mumu";
            else if (selectedLeps[0]->isElectron() && selectedLeps[1]->isElectron()) sn += "_ee";
            else sn += "_emu";
			plotter.getOrMake1DPre(sn+"_passTrigPreSel","evts",";M_{X}",50,600,4600)->Fill(ht_chs,weight);

        	if (!passIP(selectedLeps[0],selectedLeps[1])) return false;
        	testISO(sn,selectedLeps[0],selectedLeps[1],hbbjet);
        	testID(sn,selectedLeps[0],selectedLeps[1],hbbjet);
        }
        return true;
    }

    void write(TString fileName){
    	plotter.write(fileName);
/*		if (reader_event->process == FillerConstants::SIGNAL) {
		    float ee_num = plotter.getOrMake1DPre(TString::Format("m%i_ee_lepmatch_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float ee_den = plotter.getOrMake1DPre(TString::Format("m%i_ee_baseline_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float emu_num = plotter.getOrMake1DPre(TString::Format("m%i_emu_lepmatch_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float emu_den = plotter.getOrMake1DPre(TString::Format("m%i_emu_baseline_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float mumu_num = plotter.getOrMake1DPre(TString::Format("m%i_mumu_lepmatch_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float mumu_den = plotter.getOrMake1DPre(TString::Format("m%i_mumu_baseline_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();

			printf("\nMass = %d GeV\n",signal_mass);
			printf("emu eff = %f\n",emu_num/emu_den);
			printf("ee eff = %f\n",ee_num/ee_den);
			printf("mumu eff = %f\n",mumu_num/mumu_den);
    	} else {
		    float ee_num = plotter.getOrMake1DPre(smpName+"_ee_passLepSel_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02","evts",";M_{X}",50,600,4600)->GetEntries();
			float emu_num = plotter.getOrMake1DPre(smpName+"_emu_passLepSel_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02","evts",";M_{X}",50,600,4600)->GetEntries();
			float mumu_num = plotter.getOrMake1DPre(smpName+"_mumu_passLepSel_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02","evts",";M_{X}",50,600,4600)->GetEntries();
			float num = plotter.getOrMake1DPre(smpName+"_passLepSel_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02","evts",";M_{X}",50,600,4600)->GetEntries();
			float den = plotter.getOrMake1DPre(smpName+"_baseline_ept10_mupt10_eID_L_muID_L_eISO_01_muISO_02","evts",";M_{X}",50,600,4600)->GetEntries();

			std::cout<< smpName <<std::endl;
			printf("eff = %f\n",num/den);
			printf("emu eff = %f\n",emu_num/den);
			printf("ee eff = %f\n",ee_num/den);
			printf("mumu eff = %f\n",mumu_num/den);
    	}
    	*/
    }
    HistGetter plotter;
};

#endif

void studyDilepLeptonSel(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void studyDilepLeptonSel(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
