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

    void plotAcc(TString sn) {
		plotter.getOrMake1DPre(sn,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
    }

    void plotSpectra(TString sn, const Lepton* recolep1, const Lepton* recolep2) {
    	MomentumF dilepMOM = recolep1->p4() + recolep2->p4();
    	MomentumF recoHww = dilepMOM.p4() + reader_event->met.p4();
    	double dR_ll = PhysicsUtilities::deltaR(*recolep1,*recolep2);
    	double dPhi_ll = PhysicsUtilities::deltaPhi(*recolep1,*recolep2);

    	plotter.getOrMake1DPre(sn,"Mll",";m_{ll}",100,0,200)->Fill(dilepMOM.mass(),weight);
    	plotter.getOrMake1DPre(sn,"Mww",";m_{WW}",100,0,200)->Fill(recoHww.mass(),weight);
    	plotter.getOrMake1DPre(sn,"dR_ll",";#DeltaR_{ll}",50,0,5)->Fill(dR_ll,weight);
    	plotter.getOrMake1DPre(sn,"dPhi_ll",";#Delta#Phi_{ll}",50,-3.14,3.14)->Fill(dPhi_ll,weight);
    	plotter.getOrMake1DPre(sn,"pt1",";p_{T} lep1",100,0,1000)->Fill(recolep1->pt(),weight);
    	plotter.getOrMake1DPre(sn,"pt2",";p_{T} lep2",100,0,1000)->Fill(recolep2->pt(),weight);

/*    	if (reader_event->process == FillerConstants::SIGNAL) {
        	double genlep2_pt = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1->pt() : diHiggsEvt.w1_d1->pt();
        	plotter.getOrMake1DPre(sn,"genpt2",";p_{T}",40,0,200)->Fill(genlep2_pt,weight);
    	}
    	*/
    }

    void addPlots(TString sn, const Lepton* recolep1, const Lepton* recolep2) {
    	sn += "_ept10_mupt10";
    	plotAcc(sn);
    	plotSpectra(sn,recolep1,recolep2);

    	TString id_abbv, iso_abbv;
    	std::vector<double> iso_vals = {0.1, 0.2, 0.3};
    	std::vector<TString> ids = {"loose","medium","tight","high"};

    	for (const auto& id : ids) {
    		if (id=="loose") id_abbv = "L";
    		else if (id=="medium") id_abbv = "M";
    		else if (id=="tight") id_abbv = "T";
    		else id_abbv = "H";

    		if (passID(recolep1,id,id) && passID(recolep2,id,id)) {
    			plotAcc(sn+"_eID_"+id_abbv+"_muID_"+id_abbv);
    			plotSpectra(sn+"_eID_"+id_abbv+"_muID_"+id_abbv,recolep1,recolep2);

        		for (const auto& iso : iso_vals) {
/*        			std::cout<<iso<<std::endl;
        			std::cout<< (iso == 0.1) << std::endl;
        			std::cout<<(iso==0.2)<<std::endl;
        			std::cout<<(iso==0.3)<<std::endl;
*/
        			if (iso == 0.1) iso_abbv = "01";
        			else if (iso == 0.2) iso_abbv = "02";
        			else iso_abbv = "03";

        			if (passISO(recolep1,iso,iso,0) && passISO(recolep2,iso,iso,0)) {
            			plotAcc(sn+"_eID_"+id_abbv+"_muID_"+id_abbv+"_eMISO_"+iso_abbv+"_muMISO_"+iso_abbv);
            			plotSpectra(sn+"_eID_"+id_abbv+"_muID_"+id_abbv+"_eMISO_"+iso_abbv+"_muMISO_"+iso_abbv,recolep1,recolep2);
        			}
        			if (passISO(recolep1,iso,iso,1) && passISO(recolep2,iso,iso,1)) {
            			plotAcc(sn+"_eID_"+id_abbv+"_muID_"+id_abbv+"_eRISO_"+iso_abbv+"_muRISO_"+iso_abbv);
            			plotSpectra(sn+"_eID_"+id_abbv+"_muID_"+id_abbv+"_eRISO_"+iso_abbv+"_muRISO_"+iso_abbv,recolep1,recolep2);
        			}
        		}
    		}
    	}
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

    bool findHbbCand(const Lepton* lep1, const Lepton* lep2) {
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
        return true;
    }

    bool matchDileptons(TString sn, const GenParticle* gen1, const GenParticle* gen2, const Lepton* reco1, const Lepton* reco2, std::vector<const Lepton*> leps) {
    	const auto allmuons = PhysicsUtilities::selObjsMom(reader_muon->muons,0,99);
    	const auto allelectrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,0,99);

    	// find a match to each GenLepton in the reco muon and electron collections by looking for closest in DR
    	double nearestDR1 = 10;
    	double nearestDR2 = 10;
    	int idx1 = 99;
    	int idx2 = 99;
    	if (gen1->absPdgId() == 13) {
    		idx1 = PhysicsUtilities::findNearestDRDeref(*gen1,allmuons,nearestDR1);
    	} else {
    		idx1 = PhysicsUtilities::findNearestDRDeref(*gen1,allelectrons,nearestDR1);
    	}
    	if (gen2->absPdgId() == 13) {
    		idx2 = PhysicsUtilities::findNearestDRDeref(*gen2,allmuons,nearestDR2);
    	} else {
    		idx2 = PhysicsUtilities::findNearestDRDeref(*gen2,allelectrons,nearestDR2);
    	}

    	int reco1_id = 0;
    	int reco2_id = 0;

    	if (reco1->isMuon()) reco1_id = (-1)*reco1->q()*13;
    	else reco1_id = (-1)*reco1->q()*11;

    	if (reco2->isMuon()) reco2_id = (-1)*reco2->q()*13;
    	else reco2_id = (-1)*reco2->q()*11;

    	// check if the two selected reco leptons are matches to the genleptons (degeneracy = 2)
    	if (gen1->pdgId() == reco1_id && gen2->pdgId() == reco2_id) {
    		if (idx1 == reco1->index() && idx2 == reco2->index()) return true;
    	    else {
//    		    printDebugInfo(sn,gen1,gen2,idx1,idx2,leps);
    		    return false;
    	    }
    	} else if (gen1->pdgId() == reco2_id && gen2->pdgId() == reco1_id) {
    		if (idx1 == reco2->index() && idx2 == reco1->index()) return true;
    		else {
//    		    printDebugInfo(sn,gen1,gen2,idx1,idx2,leps);
    		    return false;
    		}
    	} else {
//		    printDebugInfo(sn,gen1,gen2,idx1,idx2,leps);
		    return false;
    	}
    }

    Analyzer::RecoLepSelection selectRecoLeptons() {
    	RecoLepSelection LepInfo;

    	// selection on leptons
    	float minPt1_mu = 26;
    	float minPt2_mu = 10;
    	float minPt1_e = 30;
    	float minPt2_e = 10;

    	float maxEta_e = 2.5;
    	float maxEta_mu = 2.4;

    	float minDileptonMass = 0;
    	float maxDileptonMass = 9999;
    	float maxDileptonDR = 99;

    	// get reco leptons that pass the maxEta cut and minPt cut
    	const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,minPt2_mu,maxEta_mu);
    	const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,minPt2_e,maxEta_e);

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
        	if (PhysicsUtilities::deltaR2(*leps[0],*leps[k]) > maxDileptonDR*maxDileptonDR) continue; // cut on dR_ll

        	// require that the dilepton mass be within some window
        	const MomentumF dileptonMOM = leps[0]->p4() + leps[k]->p4();
        	if (dileptonMOM.mass() < minDileptonMass || dileptonMOM.mass() > maxDileptonMass) continue;

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

    bool passID(const Lepton* lep, TString cat_mu, TString cat_e) {
    	if (lep->isMuon()) {
    		if (cat_mu == "loose") {
    			if (!((const Muon*)lep)->passLooseID()) return false;
    		} else if (cat_mu == "medium") {
    			if (!((const Muon*)lep)->passMed16ID()) return false;
    		} else if (cat_mu == "tight") {
    			if (!((const Muon*)lep)->passTightID()) return false;
    		} else if (cat_mu == "high") {
    			if (!((const Muon*)lep)->passHighPT()) return false;
    		} else {
    			printf("Error with ID cat provided");
    		}
    	} else if (lep->isElectron()) {
    		if (cat_e == "loose") {
    			if (!((const Electron*)lep)->passLooseID_noISO()) return false;
    		} else if (cat_e == "medium") {
    			if (!((const Electron*)lep)->passMedID_noISO()) return false;
    		} else if (cat_e == "tight") {
    			if (!((const Electron*)lep)->passTightID_noISO()) return false;
    		} else if (cat_e == "high") {
    			if (!((const Electron*)lep)->passHEEPID_noISO()) return false;
    		} else {
    			printf("Error with ID cat provided");
    		}
    	} else {
    		printf("lep in id not muon or electron\n");
    	}
    	return true;
    }

    bool passISO(const Lepton* lep, float maxIso_mu, float maxIso_e, int iso_type) {
    	if (iso_type == 0) {
    	    if (lep->isMuon()) {
    		    if (lep->miniIso() > maxIso_mu) return false;
    	    } else if (lep->isElectron()) {
    	    	if (lep->miniIso() > maxIso_e) return false;
    	    } else {
    		    printf("lep in iso not muon or electron\n");
    	    }
    	} else if (iso_type == 1) {
    	    if (lep->isMuon()) {
    		    if (((const Muon*)lep)->dbRelISO() > maxIso_mu) return false;
    	    } else if (lep->isElectron()) {
    	    	if (((const Electron*)lep)->eaRelISO() > maxIso_e) return false;
    	    } else {
    		    printf("lep in iso not muon or electron\n");
    	    }
    	} else std::cout << "Error choosing isolation type" << std::endl;
    	return true;
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
			plotAcc(sn+"_baseline");

			// Select the reco dileptons
        	RecoLepSelection LepInfo = selectRecoLeptons();
        	std::vector<const Lepton*> selectedLeps = LepInfo.selectedLeps;

        	if (selectedLeps.size() < 2) return false;
        	if (selectedLeps.size() > 2) printf("WARNING: Selected more than 2 leps\n");

        	// Check if the selected RecoLeps match with the ones associated to the GenParticles
        	if (!matchDileptons(sn,genlep1, genlep2, selectedLeps[0], selectedLeps[1], LepInfo.filteredLeps)) return false;

        	// Check if a valid Hbb candidate is found in the event
        	if (!findHbbCand(selectedLeps[0], selectedLeps[1])) return false;

        	plotAcc(sn+"_lepmatch");
        	addPlots(sn,selectedLeps[0],selectedLeps[1]);
        }
        // BKG
        if (reader_event->process != FillerConstants::SIGNAL) {
        	plotAcc(sn+"_baseline");

			// Select the reco dileptons
        	RecoLepSelection LepInfo = selectRecoLeptons();
        	std::vector<const Lepton*> selectedLeps = LepInfo.selectedLeps;

        	if (selectedLeps.size() < 2) return false;
        	if (selectedLeps.size() > 2) printf("WARNING: Selected more than 2 leps\n");

        	// Check if a valid Hbb candidate is found in the event
        	if (!findHbbCand(selectedLeps[0], selectedLeps[1])) return false;

            plotAcc(sn+"_passLepSel");
            if (selectedLeps[0]->isMuon() && selectedLeps[1]->isMuon()) sn += "_mumu";
            else if (selectedLeps[0]->isElectron() && selectedLeps[1]->isElectron()) sn += "_ee";
            else sn += "_emu";
            plotAcc(sn+"_passLepSel");
        	addPlots(sn,selectedLeps[0],selectedLeps[1]);
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
