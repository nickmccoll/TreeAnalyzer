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
#include "Processors/Variables/interface/LeptonSelection.h"

#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){}

    void plotSpectra(TString sn, const Lepton* recolep1, const Lepton* recolep2, const FatJet* hbb) {
    	MomentumF dilepMOM = recolep1->p4() + recolep2->p4();
    	MomentumF bbllMOM = dilepMOM.p4() + hbb->p4();
    	double dR_ll = PhysicsUtilities::deltaR(*recolep1,*recolep2);
    	double dPhi_ll = PhysicsUtilities::deltaPhi(*recolep1,*recolep2);
    	double dR_bbll = PhysicsUtilities::deltaR(dilepMOM,*hbb);

    	plotter.getOrMake1DPre(sn,"Mll",";m_{ll}",100,0,200)->Fill(dilepMOM.mass(),weight);
    	plotter.getOrMake1DPre(sn,"Mbbll",";m_{bbll}",200,0,3000)->Fill(bbllMOM.mass(),weight);
    	plotter.getOrMake1DPre(sn,"dR_ll",";#DeltaR_{ll}",50,0,5)->Fill(dR_ll,weight);
    	plotter.getOrMake1DPre(sn,"dPhi_ll",";#Delta#Phi_{ll}",50,-3.14,3.14)->Fill(dPhi_ll,weight);
    	plotter.getOrMake1DPre(sn,"pt1",";p_{T} lep1",100,0,1000)->Fill(recolep1->pt(),weight);
    	plotter.getOrMake1DPre(sn,"pt2",";p_{T} lep2",100,0,1000)->Fill(recolep2->pt(),weight);
    	plotter.getOrMake1DPre(sn,"ptbb",";p_{T} bb",150,0,1500)->Fill(hbb->pt(),weight);
    	plotter.getOrMake1DPre(sn,"dR_bbll",";#DeltaR_{bb,ll}",100,0,6)->Fill(dR_bbll,weight);
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
    		if (nGoodSJ > 1) goodSJs = true;
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
    void testISO(TString sn, const Lepton* matchLep1, const Lepton* matchLep2, bool isSignal) {
    	LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);
        proc.lepSelParams.mu_minPT = 10;
        proc.lepSelParams.el_minPT = 10;
        TString constID = "_ept10_mupt10_ID_eT_muM_";

        static const std::vector<float>                     isoWPs = {0.1,0.2,0.3};
        static const std::vector<LepSelHelpers::muFunFloat> muIsoTypes = {&Muon::miniIso,&Muon::dbRelISO};
        static const std::vector<LepSelHelpers::elFunFloat> elIsoTypes = {&Electron::miniIso,&Electron::eaRelISO};
        static const std::vector<TString>                   isoTypeName = {"miniIso","relIso"};

        for (unsigned long j=0; j<isoTypeName.size(); j++) {
        	proc.lepSelParams.mu_getISO = muIsoTypes[j];
        	proc.lepSelParams.el_getISO = elIsoTypes[j];

        	for (const auto& wp_mu : isoWPs) {
        		for (const auto& wp_e : isoWPs) {
                	proc.lepSelParams.mu_maxISO = wp_mu;
                	proc.lepSelParams.el_maxISO = wp_e;

                	TString name = TString::Format(isoTypeName[j]+"_e%.1f_mu%.1f_",wp_e,wp_mu);
                	name.ReplaceAll(".","");
                	const auto leps = proc.getLeptons(*reader_event,*reader_muon,*reader_electron);
                	if (leps.size() != 2) continue;
                	if (leps[0]->q() == leps[1]->q()) continue;

                	const Lepton* lep1 = leps.front();
                	const Lepton* lep2 = leps[1];
                	if (isSignal) {
                		if (!matchDileptons(matchLep1,matchLep2,lep1,lep2)) {
                			std::cout<<"Leptons do not match"<<std::endl;
                			continue;
                		}
                	}
                	const FatJet* hbbjet = findHbbCand(lep1,lep2);
                	if (!hbbjet) continue;

            		TString chan;
            		if (!isSignal) {
            			if (lep1->isMuon() && lep2->isMuon()) chan = "_mumu";
            			else if (lep1->isElectron() && lep2->isElectron()) chan = "_ee";
            			else chan = "_emu";
            		}
                	plotter.getOrMake1DPre(sn+chan+constID+name,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
                    plotSpectra(sn+chan+constID+name,lep1,lep2,hbbjet);
        		}
        	}
        }
    }
    void testID(TString sn, const Lepton* matchLep1, const Lepton* matchLep2, bool isSignal) {
    	LeptonProcessor proc;
        DefaultLeptonSelections::setDefaultLeptonProcessor(proc);
        proc.lepSelParams.mu_minPT = 10;
        proc.lepSelParams.el_minPT = 10;
        TString kinsel = "_ept10_mupt10_";
        TString suf = "_miniIso_e01_mu02_";

        static const std::vector<LepSelHelpers::muFunBool> muIdTypes = {&Muon::passLooseID,&Muon::passMed16ID,&Muon::passTightID,&Muon::passHighPT};
        static const std::vector<LepSelHelpers::elFunBool> elIdTypes = {&Electron::passLooseID_noISO,&Electron::passMedID_noISO,
        																&Electron::passTightID_noISO,&Electron::passHEEPID_noISO};
        static const std::vector<TString> ids = {"L","M","T","H"};

        for (unsigned long j=0; j<muIdTypes.size(); j++) {
        	for (unsigned long k=0; k<elIdTypes.size(); k++) {
        		proc.lepSelParams.mu_getID = muIdTypes[j];
        		proc.lepSelParams.el_getID = elIdTypes[k];

        		TString name = "ID_e_"+ids[k]+"_mu_"+ids[j];
        		const auto leps = proc.getLeptons(*reader_event,*reader_muon,*reader_electron);
        		if (leps.size() != 2) continue;
        		if (leps[0]->q() == leps[1]->q()) continue;

        		const Lepton* lep1 = leps.front();
        		const Lepton* lep2 = leps[1];
        		if (isSignal) {
        		    if (!matchDileptons(matchLep1,matchLep2,lep1,lep2)) {
        			    std::cout<<"Leptons do not match"<<std::endl;
        			    continue;
        		    }
        		}
        		const FatJet* hbbjet = findHbbCand(lep1,lep2);
        		if (!hbbjet) continue;

        		TString chan;
        		if (!isSignal) {
        			if (lep1->isMuon() && lep2->isMuon()) chan = "_mumu";
        			else if (lep1->isElectron() && lep2->isElectron()) chan = "_ee";
        			else chan = "_emu";
        		}
        		plotter.getOrMake1DPre(sn+chan+kinsel+name+suf,"evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        		plotSpectra(sn+chan+kinsel+name+suf,lep1,lep2,hbbjet);
        	}
        }
    }
    void testEachLep(TString sn, const Lepton* matchLep1, const Lepton* matchLep2, bool isSignal) {
    	LeptonProcessor proc;
    	DefaultLeptonSelections::setDefaultLeptonProcessor(proc);
        proc.lepSelParams.mu_minPT = 10;
        proc.lepSelParams.el_minPT = 10;
        proc.lepSelParams.el_maxISO = 9999;
        proc.lepSelParams.mu_maxISO = 9999;
        proc.lepSelParams.el_getID = &Electron::passVetoID_noISO;
        proc.lepSelParams.mu_getID = &Muon::passSoftID;

        sn += "_ept10_mupt10_";
        TString stdiso = "_miniIso_e01_mu02_";

        static const std::vector<TString> ids = {"L","M","T","H"};
        static const std::vector<double> isos = {0.1,0.2,0.3};

        const auto leps = proc.getLeptons(*reader_event,*reader_muon,*reader_electron);
        const Lepton* lep1=nullptr;
        const Lepton* lep2=nullptr;
        if (leps.size() > 0) lep1 = leps.front();
        if (leps.size() > 1) {
        	for (unsigned long k=1; k<leps.size(); k++) {
        		if (lep1->q() == leps[k]->q()) continue;
        		lep2 = leps[k];
        	}
        }
        if (lep1 && lep2) {
            if (isSignal) {
            	if (!matchDileptons(matchLep1,matchLep2,lep1,lep2)) return;
            }
    		bool passStdId1 = lep1->isMuon() ? ((const Muon*)lep1)->passMed16ID() : ((const Electron*)lep1)->passTightID_noISO();
    		bool passStdId2 = lep2->isMuon() ? ((const Muon*)lep2)->passMed16ID() : ((const Electron*)lep2)->passTightID_noISO();
    		bool passStdIso1 = lep1->isMuon() ? lep1->miniIso() < 0.2 : lep1->miniIso() < 0.1;
    		bool passStdIso2 = lep2->isMuon() ? lep2->miniIso() < 0.2 : lep2->miniIso() < 0.1;
    		// ID
        	for (unsigned long id=0; id<ids.size(); id++) {
        		if (passStdId2 && passStdIso1 && passStdIso2) {
        			bool goodId1 = false;
        			if (id==0) goodId1 = lep1->isMuon() ? ((const Muon*)lep1)->passLooseID() : ((const Electron*)lep1)->passLooseID_noISO();
        			else if (id==1) goodId1 = lep1->isMuon() ? ((const Muon*)lep1)->passMed16ID() : ((const Electron*)lep1)->passMedID_noISO();
        			else if (id==2) goodId1 = lep1->isMuon() ? ((const Muon*)lep1)->passTightID() : ((const Electron*)lep1)->passTightID_noISO();
        			else if (id==3) goodId1 = lep1->isMuon() ? ((const Muon*)lep1)->passHighPT() : ((const Electron*)lep1)->passHEEPID_noISO();

        			if (goodId1) plotter.getOrMake1DPre(sn+"ID_eT_muM"+stdiso,"_evts_lep1_"+ids[id],";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        		}
        		if (passStdId1 && passStdIso1 && passStdIso2) {
        			bool goodId2 = false;
        			if (id==0) goodId2 = lep2->isMuon() ? ((const Muon*)lep2)->passLooseID() : ((const Electron*)lep2)->passLooseID_noISO();
        			else if (id==1) goodId2 = lep2->isMuon() ? ((const Muon*)lep2)->passMed16ID() : ((const Electron*)lep2)->passMedID_noISO();
        			else if (id==2) goodId2 = lep2->isMuon() ? ((const Muon*)lep2)->passTightID() : ((const Electron*)lep2)->passTightID_noISO();
        			else if (id==3) goodId2 = lep2->isMuon() ? ((const Muon*)lep2)->passHighPT() : ((const Electron*)lep2)->passHEEPID_noISO();

        			if (goodId2) plotter.getOrMake1DPre(sn+"ID_eT_muM"+stdiso,"_evts_lep2_"+ids[id],";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        		}
        	}
        	// ISO
        	for (const auto& iso : isos) {
    			TString isoname = TString::Format("%.1f",iso);
    			isoname.ReplaceAll(".","");
        		if (passStdIso2 && passStdId1 && passStdId2) {
        			if (lep1->miniIso() < iso) plotter.getOrMake1DPre(sn+"ID_eT_muM"+stdiso,"_evts_lep1_"+isoname,";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        		}
        		if (passStdIso1 && passStdId1 && passStdId2) {
        			if (lep2->miniIso() < iso) plotter.getOrMake1DPre(sn+"ID_eT_muM"+stdiso,"_evts_lep2_"+isoname,";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        		}
        	}
        }
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
    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(ht_chs < 400) return false;
        TString sn = smpName;

        // SIGNAL
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type == DiHiggsEvent::DILEP) {

            const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,2.5);
            const auto* matchLep1 = getMatchedLepton(*diHiggsEvt.w1_d1,muons,electrons);
            const auto* matchLep2 = getMatchedLepton(*diHiggsEvt.w2_d1,muons,electrons);

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
			plotter.getOrMake1DPre(sn+"_baseline_","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

        	testISO(sn,matchLep1,matchLep2,true);
        	testID(sn,matchLep1,matchLep2,true);
        	testEachLep(sn,matchLep1,matchLep2,true);
        }
        // BKG
        if (reader_event->process != FillerConstants::SIGNAL) {
        	plotter.getOrMake1DPre(sn+"_baseline_","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);
        	testISO(sn,0,0,false);
        	testID(sn,0,0,false);
        	testEachLep(sn,0,0,false);
        }
        return true;
    }

    void write(TString fileName){
    	plotter.write(fileName);
    }
    HistGetter plotter;
};

#endif

void DileptonSelection(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void DileptonSelection(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
