
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


using namespace TAna;
using namespace std;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

//        turnOffCorr(CORR_TRIG);
//        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
//        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);
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

    void plot(TString sn, const Lepton *lep1, const Lepton *lep2, const FatJet *hbb) {

    	MomentumF dilepMOM = lep1->p4() + lep2->p4();
    	double dR_ll = PhysicsUtilities::deltaR(*lep1,*lep2);
    	double dPhi_ll = PhysicsUtilities::absDeltaPhi(*lep1,*lep2);
    	double dEta_ll = PhysicsUtilities::absDeltaEta(*lep1,*lep2);
    	double dPhi_metll = PhysicsUtilities::absDeltaPhi(reader_event->met,dilepMOM);
    	double mll = dilepMOM.mass();

    	MomentumF ww = Hww2lSolver::getSimpleHiggsMom(lep1->p4()+lep2->p4(),reader_event->met,40);
    	MomentumF xhh = ww.p4() + hbb->p4();
		auto bbmass = isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(hbb) : hbb->sdMom().mass();

    	auto plt = [&](TString pre) {
    		plotter.getOrMake1DPre(pre,"pt1",";p_{T}",100,0,2000)->Fill(lep1->pt(),weight);
    		plotter.getOrMake1DPre(pre,"pt2",";p_{T}",100,0,2000)->Fill(lep2->pt(),weight);
    		plotter.getOrMake1DPre(pre,"eta1",";#eta",100,-2.4,2.4)->Fill(lep1->eta(),weight);
    		plotter.getOrMake1DPre(pre,"eta2",";#eta",100,-2.4,2.4)->Fill(lep2->eta(),weight);
    		plotter.getOrMake1DPre(pre,"ptll",";p_{T}",100,0,2000)->Fill(dilepMOM.pt(),weight);
    		plotter.getOrMake1DPre(pre,"ptbb",";p_{T}",100,0,2000)->Fill(hbb->pt(),weight);
        	plotter.getOrMake1DPre(pre,"Mww",";M_{WW}",100,0,400)->Fill(ww.mass(),weight);
        	plotter.getOrMake1DPre(pre,"ptww",";p_{T,WW}",100,0,2000)->Fill(ww.pt(),weight);
        	plotter.getOrMake1DPre(pre,"ptwwomhh",";p_{T,WW}/M_{HH}",100,0,1)->Fill(ww.pt()/xhh.mass(),weight);
        	plotter.getOrMake1DPre(pre,"Mbb",";M_{bb}",100,30,210)->Fill(bbmass,weight);
        	plotter.getOrMake1DPre(pre,"Mhh",";M_{HH}",100,700,4500)->Fill(xhh.mass(),weight);
        	plotter.getOrMake1DPre(pre,"dRll",";#DeltaR_{ll}",100,0,4.0)->Fill(dR_ll,weight);
        	plotter.getOrMake1DPre(pre,"dPhill",";#Delta#phi_{ll}",100,0,TMath::Pi())->Fill(dPhi_ll,weight);
        	plotter.getOrMake1DPre(pre,"dEtall",";#Delta#eta_{ll}",100,0,4.0)->Fill(dEta_ll,weight);
        	plotter.getOrMake1DPre(pre,"dPhi_metll",";#Delta#Phi_{met,ll}",100,0,TMath::Pi())->Fill(dPhi_metll,weight);
        	plotter.getOrMake1DPre(pre,"Mll",";M_{ll}",100,0,200)->Fill(mll,weight);
        	plotter.getOrMake1DPre(pre,"met",";MET",100,0,2000)->Fill(reader_event->met.pt(),weight);
        	plotter.getOrMake1DPre(pre,"metOht",";MET",100,0,1)->Fill(reader_event->met.pt()/ht_puppi,weight);
        	plotter.getOrMake1DPre(pre,"metOrht",";MET",100,0,20)->Fill(reader_event->met.pt()/sqrt(ht_puppi),weight);
        	plotter.getOrMake1DPre(pre,"metOmhh",";MET/M_{HH}",100,0,1)->Fill(reader_event->met.pt()/xhh.mass(),weight);
        	plotter.getOrMake1DPre(pre,"mtww",";M_{T}",100,0,2000)->Fill(ww.mt(),weight);
        	plotter.getOrMake1DPre(pre,"yww",";y",100,-5,5)->Fill(ww.rap(),weight);

        	plotter.getOrMake2DPre(pre,"ptwwXmww",";p_{T,WW};M_{WW}",100,0,2000,100,0,400)->Fill(ww.pt(),ww.mass(),weight);
    	};

    	plt(sn);

    	bool passMll = mll > 6 && mll < 75;
    	bool passDR = dR_ll < 1.0;
    	bool passDPhi = dPhi_metll < TMath::PiOver2();
    	bool passMET = reader_event->met.pt()/xhh.mass() > 0.05;

    	if (passMll) plt(sn+"_Mll5to75");
    	if (dR_ll < 1.6)  plt(sn+"_dRll1p6");

    	if (!passMll) return;

    	if (passDR && passDPhi)  plt(sn+"_relaxMET");
    	if (passDR && passMET)   plt(sn+"_relaxDPhi");
    	if (passMET && passDPhi) plt(sn+"_relaxDR");

    	if (passDR && passMET && passDPhi) plt(sn+"_fullSel");
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

    bool passLepSel(const Lepton* lep1, const Lepton* lep2) {
    	bool passIP1 = fabs(lep1->dz()) < 0.1 && fabs(lep1->d0()) < 0.05;
    	bool passIP2 = fabs(lep2->dz()) < 0.1 && fabs(lep2->d0()) < 0.05;

    	bool passIso = lep1->miniIso() < 0.15 && lep2->miniIso() < 0.15;

    	bool passId1 = lep1->isMuon() ? ((Muon*)lep1)->passLooseID() : ((Electron*)lep1)->passMedID_noIso();
    	bool passId2 = lep2->isMuon() ? ((Muon*)lep2)->passLooseID() : ((Electron*)lep2)->passMedID_noIso();

    	bool pass = passIP1 && passIP2 && passIso && passId1 && passId2;

    	if (pass) return true;
    	else      return false;
    }

    void makePlots(TString sn, TString regS, bool doBtagRegs, const Lepton* lep1, const Lepton* lep2, const FatJet* hbb) {

    	TString llS = DileptonProcessor::getDilepStr(lep1,lep2);
    	plot(sn+"_"+regS,lep1,lep2,hbb);
    	plot(sn+llS+regS,lep1,lep2,hbb);

    	if (doBtagRegs) {
    		TString bS = "";
    		if      (hbbCSVCat_2l == BTagging::CSVSJ_MF) bS = "_bL";
    		else if (hbbCSVCat_2l == BTagging::CSVSJ_ML) bS = "_bM";
    		else if (hbbCSVCat_2l == BTagging::CSVSJ_MM) bS = "_bT";
    		else return;

        	plot(sn+"_"+regS+bS,lep1,lep2,hbb);
        	plot(sn+llS+regS+bS,lep1,lep2,hbb);
    	}
    }

    void doNormal(TString sn) {
        if (selectedDileptons.size() != 2) return;
        if (!hbbCand_2l) return;
        if (hbbMass_2l < 30 || hbbMass_2l > 210) return;
        if (hh_2l.mass() < 700) return;

        bool passPt1 = selectedDileptons[0]->isMuon() ? selectedDileptons[0]->pt() > 27 : selectedDileptons[0]->pt() > 30;
        if (!passPt1) return;

		// plot dif regions
		if (nMedBTags_HbbV_2l == 0 && hbbCSVCat_2l >= BTagging::CSVSJ_MF) {
			makePlots(sn,"SR",true,dilep1,dilep2,hbbCand_2l);
		}
		if (nMedBTags_HbbV_2l  > 0 && hbbCSVCat_2l >= BTagging::CSVSJ_MF) {
			makePlots(sn,"TopCR",true,dilep1,dilep2,hbbCand_2l);
		}
		if (nMedBTags_HbbV_2l == 0 && hbbCSVCat_2l == BTagging::CSVSJ_FF) {
			makePlots(sn,"QgCR",false,dilep1,dilep2,hbbCand_2l);
		}
    }

    void doGenMatching(TString sn) {
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
    	plotter.getOrMake1DPre(sn+"_gen_dilep","evts",";M_{X}",50,600,4600)->Fill(signal_mass,weight);

    	bool passPt1 = glep1->absPdgId() == 13 ? (glep1->pt() > 27) : (glep1->pt() > 30);
    	bool passPt2 = (glep2->pt() > 10);
    	if (!(passPt1 && passPt2)) return;

    	const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,10,2.4);
    	const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,10,2.5);
    	const Lepton* matchLep1 = getMatchedLepton(glep1,muons,electrons,0.1,true);
    	const Lepton* matchLep2 = getMatchedLepton(glep2,muons,electrons,0.1,true);

    	if (!matchLep1 || !matchLep2) return;
    	if ((matchLep1->isMuon() == matchLep2->isMuon()) && (matchLep1->index() == matchLep2->index())) return; // discard if these are the same RECO lep

    	// RECO lepton pt cuts
    	bool passRecoPt1 = matchLep1->isMuon() ? (matchLep1->pt() > 27) : (matchLep1->pt() > 30);
    	bool passRecoPt2 = (matchLep2->pt() > 10);
    	if (!passRecoPt1 || !passRecoPt2) return;
    	if (!passLepSel(matchLep1,matchLep2)) return;

    	const FatJet* matchHbb = getSigHbb(matchLep1,matchLep2);
    	if (!matchHbb) return;
    	plotter.getOrMake1DPre(sn+"_validHbb","evts",";M_{X}",5000,0,5000)->Fill(signal_mass,weight);

    	auto bbmass = isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(matchHbb) : matchHbb->sdMom().mass();
    	if (bbmass < 30 || bbmass > 210) return;
    	if (!passLepSel(matchLep1,matchLep2)) return;

    	const auto bbCat = BTagging::getCSVSJCat(parameters.jets,matchHbb->subJets());

    	const auto jetsVeto = PhysicsUtilities::selObjsD(jets,
    		        [&](const Jet* j){return  PhysicsUtilities::deltaR2(*j,*matchHbb ) >= 1.2*1.2;});
    	const auto nBtags = PhysicsUtilities::selObjsMomD(jetsVeto,
    	            parameters.jets.minBtagJetPT,parameters.jets.maxBTagJetETA,
    		        [&](const Jet* j){return BTagging::passJetBTagWP(parameters.jets,*j);} ).size();

    	if (nBtags == 0 && bbCat >= BTagging::CSVSJ_MF) {
    		makePlots(sn+"_genMatch","SR",true,matchLep1,matchLep2,matchHbb);
    	}
    	if (nBtags > 0  && bbCat >= BTagging::CSVSJ_MF) {
    		makePlots(sn+"_genMatch","TopCR",true,matchLep1,matchLep2,matchHbb);
    	}
    	if (nBtags == 0 && bbCat == BTagging::CSVSJ_FF) {
    		makePlots(sn+"_genMatch","QgCR",false,matchLep1,matchLep2,matchHbb);
    	}

    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht_puppi < 400) return false;
	    if(!EventSelection::passTriggerSuite2017(*reader_event)) return false;

	    TString sn = smpName;
        if (mcProc == FillerConstants::TTBAR && smDecayEvt.nLepsTT >= 0 && smDecayEvt.nLepsTT <= 2) {
        	sn += TString::Format("%d",smDecayEvt.nLepsTT);
        }

        if (isSignal()) {
        	if (diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        	doGenMatching(sn);
        }

        doNormal(sn);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getSelVars(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
