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
		plotter.getOrMake1DPre(sn+"_ept15_mupt15","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
    }

    void plotSpectra(TString sn, const Lepton* recolep1, const Lepton* recolep2) {
    	const MomentumF dilepMOM = recolep1->p4() + recolep2->p4();
    	double dR_ll = PhysicsUtilities::deltaR(*recolep1,*recolep2);

    	plotter.getOrMake1DPre(sn+"_ept15_mupt15","Mll",";m_{ll}",100,0,150)->Fill(dilepMOM.mass(),weight);
    	plotter.getOrMake1DPre(sn+"_ept15_mupt15","dR_ll",";#DeltaR_{ll}",50,0,5)->Fill(dR_ll,weight);
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
    	float minPt2_mu = 15;
    	float minPt1_e = 30;
    	float minPt2_e = 15;
    	float maxEta_e = 2.5;
    	float maxEta_mu = 2.4;
    	float minDileptonMass = 0;
    	float maxDileptonMass = 999;
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

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(ht_chs < 400) return false;
        if(reader_event->process >= FillerConstants::ZJETS && reader_event->process <= FillerConstants::TTX )
            smpName = "other";
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
        	plotAcc(sn+"_lepmatch");
        	plotSpectra(sn,selectedLeps[0],selectedLeps[1]);
        }
        // BKG
        if (reader_event->process != FillerConstants::SIGNAL) {
        	plotAcc(sn+"_baseline");

			// Select the reco dileptons
        	RecoLepSelection LepInfo = selectRecoLeptons();
        	std::vector<const Lepton*> selectedLeps = LepInfo.selectedLeps;

        	if (selectedLeps.size() < 2) return false;
        	if (selectedLeps.size() > 2) printf("WARNING: Selected more than 2 leps\n");

            plotAcc(sn+"_passLepSel");
            if (selectedLeps[0]->isMuon() && selectedLeps[1]->isMuon()) sn += "_mumu";
            else if (selectedLeps[0]->isElectron() && selectedLeps[1]->isElectron()) sn += "_ee";
            else sn += "_emu";
            plotAcc(sn+"_passLepSel");
        	plotSpectra(sn,selectedLeps[0],selectedLeps[1]);
        }
        return true;
    }

    void write(TString fileName){
    	plotter.write(fileName);
		if (reader_event->process == FillerConstants::SIGNAL) {
		    float ee_num = plotter.getOrMake1DPre(TString::Format("m%i_ee_lepmatch_ept15_mupt15",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float ee_den = plotter.getOrMake1DPre(TString::Format("m%i_ee_baseline_ept15_mupt15",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float emu_num = plotter.getOrMake1DPre(TString::Format("m%i_emu_lepmatch_ept15_mupt15",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float emu_den = plotter.getOrMake1DPre(TString::Format("m%i_emu_baseline_ept15_mupt15",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float mumu_num = plotter.getOrMake1DPre(TString::Format("m%i_mumu_lepmatch_ept15_mupt15",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();
			float mumu_den = plotter.getOrMake1DPre(TString::Format("m%i_mumu_baseline_ept15_mupt15",signal_mass),"evts",";M_{X}",50,600,4600)->GetEntries();

			printf("\nMass = %d GeV\n",signal_mass);
			printf("emu eff = %f\n",emu_num/emu_den);
			printf("ee eff = %f\n",ee_num/ee_den);
			printf("mumu eff = %f\n",mumu_num/mumu_den);
    	} else {
		    float ee_num = plotter.getOrMake1DPre(smpName+"_ee_passLepSel_ept15_mupt15","evts",";M_{X}",50,600,4600)->GetEntries();
			float emu_num = plotter.getOrMake1DPre(smpName+"_emu_passLepSel_ept15_mupt15","evts",";M_{X}",50,600,4600)->GetEntries();
			float mumu_num = plotter.getOrMake1DPre(smpName+"_mumu_passLepSel_ept15_mupt15","evts",";M_{X}",50,600,4600)->GetEntries();
			float num = plotter.getOrMake1DPre(smpName+"_passLepSel_ept15_mupt15","evts",";M_{X}",50,600,4600)->GetEntries();
			float den = plotter.getOrMake1DPre(smpName+"_baseline_ept15_mupt15","evts",";M_{X}",50,600,4600)->GetEntries();

			std::cout<< smpName <<std::endl;
			printf("eff = %f\n",num/den);
			printf("emu eff = %f\n",emu_num/den);
			printf("ee eff = %f\n",ee_num/den);
			printf("mumu eff = %f\n",mumu_num/den);
    	}
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
