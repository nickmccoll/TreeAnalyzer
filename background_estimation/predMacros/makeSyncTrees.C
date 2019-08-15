
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Corrections/interface/TriggerScaleFactors.h"
#include "Processors/Corrections/interface/BTagScaleFactors.h"
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "Processors/Corrections/interface/JetAndMETCorrections.h"

#include "TPRegexp.h"
using namespace std;
using namespace TAna;
using namespace CorrHelp;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {}

    virtual BaseEventAnalyzer * setupEventAnalyzer() override {return new CopierEventAnalyzer();}

    void loadVariables() override {
        reader_event       =loadReader<EventReader>   ("event",isRealData());
        reader_fatjet      =loadReader<FatJetReader>  ("ak8PuppiJet",isRealData());
        reader_fatjet_noLep=loadReader<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false,false);
        reader_jet_chs     =loadReader<JetReader>     ("ak4Jet",isRealData());
        reader_jet         =loadReader<JetReader>     ("ak4PuppiJet",isRealData(),false);
        reader_electron    =loadReader<ElectronReader>("electron",true);
        reader_muon        =loadReader<MuonReader>    ("muon");

        if(!isRealData()){
            reader_genpart =loadReader<GenParticleReader>   ("genParticle");
        }

        checkConfig();
    }

    virtual void bookOutputVariables() override {

        outTree->addSingle(event,  "",  "event");
        outTree->addSingle(lumiSec,  "",  "ls");
        outTree->addSingle(run,  "",  "run");
        outTree->addSingle(n_mu,  "",  "n_presel_mu");
        outTree->addSingle(n_el,  "",  "n_presel_ele");
        outTree->addSingle(n_ak4,  "",  "n_presel_ak4Jet");
        outTree->addSingle(n_ak8,  "",  "n_presel_ak8Jet");
        outTree->addSingle(n_ak8LS,  "",  "n_presel_ak8lsJet");
        outTree->addSingle(pu_weight,  "",  "PU_weight");
        outTree->addSingle(mc_weight,  "",  "MC_weight");
        outTree->addSingle(met_pt,  "",  "PFMET");
        outTree->addSingle(met_phi,  "",  "PFMETphi");

        outTree->addSingle(mu1_pt,  "",  "mu1_pt");
        outTree->addSingle(mu1_eta,  "",  "mu1_eta");
        outTree->addSingle(mu1_phi,  "",  "mu1_phi");
        outTree->addSingle(mu1_E,  "",  "mu1_E");
        outTree->addSingle(mu1_charge,  "",  "mu1_charge");
        outTree->addSingle(mu1_miniIso,  "",  "mu1_miniRelIso");
        outTree->addSingle(mu1_sip,  "",  "mu1_sip3D");
        outTree->addSingle(mu1_dz,  "",  "mu1_dz");
        outTree->addSingle(mu1_dxy,  "",  "mu1_dxy");
        outTree->addSingle(mu1_dxyAbs,  "",  "mu1_dxyAbs");

        outTree->addSingle(mu2_pt,  "",  "mu2_pt");
        outTree->addSingle(mu2_eta,  "",  "mu2_eta");
        outTree->addSingle(mu2_phi,  "",  "mu2_phi");
        outTree->addSingle(mu2_E,  "",  "mu2_E");
        outTree->addSingle(mu2_charge,  "",  "mu2_charge");
        outTree->addSingle(mu2_miniIso,  "",  "mu2_miniRelIso");
        outTree->addSingle(mu2_sip,  "",  "mu2_sip3D");
        outTree->addSingle(mu2_dz,  "",  "mu2_dz");
        outTree->addSingle(mu2_dxy,  "",  "mu2_dxy");
        outTree->addSingle(mu2_dxyAbs,  "",  "mu2_dxyAbs");

        outTree->addSingle(el1_pt,  "",  "ele1_pt");
        outTree->addSingle(el1_eta,  "",  "ele1_eta");
        outTree->addSingle(el1_phi,  "",  "ele1_phi");
        outTree->addSingle(el1_E,  "",  "ele1_E");
        outTree->addSingle(el1_charge,  "",  "ele1_charge");
        outTree->addSingle(el1_miniIso,  "",  "ele1_miniRelIso");
        outTree->addSingle(el1_sip,  "",  "ele1_sip3D");
        outTree->addSingle(el1_dz,  "",  "ele1_dz");
        outTree->addSingle(el1_dxy,  "",  "ele1_dxy");
        outTree->addSingle(el1_dxyAbs,  "",  "ele1_dxyAbs");
        outTree->addSingle(el1_MVAID,  "",  "ele1_mvaIdFall17noIso");

        outTree->addSingle(el2_pt,  "",  "ele2_pt");
        outTree->addSingle(el2_eta,  "",  "ele2_eta");
        outTree->addSingle(el2_phi,  "",  "ele2_phi");
        outTree->addSingle(el2_E,  "",  "ele2_E");
        outTree->addSingle(el2_charge,  "",  "ele2_charge");
        outTree->addSingle(el2_miniIso,  "",  "ele2_miniIso");
        outTree->addSingle(el2_sip,  "",  "ele2_sip3D");
        outTree->addSingle(el2_dz,  "",  "ele2_dz");
        outTree->addSingle(el2_dxy,  "",  "ele2_dxy");
        outTree->addSingle(el2_dxyAbs,  "",  "ele2_dxyAbs");
        outTree->addSingle(el2_MVAID,  "",  "ele2_mvaIdFall17noIso");

        outTree->addSingle(ak4jet1_pt,  "",  "ak4Jet1_pt");
        outTree->addSingle(ak4jet1_eta,  "",  "ak4Jet1_eta");
        outTree->addSingle(ak4jet1_phi,  "",  "ak4Jet1_phi");
        outTree->addSingle(ak4jet1_E,  "",  "ak4Jet1_E");
        outTree->addSingle(ak4jet1_csv,  "",  "ak4Jet1_CSV");
        outTree->addSingle(ak4jet2_pt,  "",  "ak4Jet2_pt");
        outTree->addSingle(ak4jet2_eta,  "",  "ak4Jet2_eta");
        outTree->addSingle(ak4jet2_phi,  "",  "ak4Jet2_phi");
        outTree->addSingle(ak4jet2_E,  "",  "ak4Jet2_E");
        outTree->addSingle(ak4jet2_csv,  "",  "ak4Jet2_CSV");
        outTree->addSingle(ak4jet3_pt,  "",  "ak4Jet3_pt");
        outTree->addSingle(ak4jet3_eta,  "",  "ak4Jet3_eta");
        outTree->addSingle(ak4jet3_phi,  "",  "ak4Jet3_phi");
        outTree->addSingle(ak4jet3_E,  "",  "ak4Jet3_E");
        outTree->addSingle(ak4jet3_csv,  "",  "ak4Jet3_CSV");
        outTree->addSingle(ak4jet4_pt,  "",  "ak4Jet4_pt");
        outTree->addSingle(ak4jet4_eta,  "",  "ak4Jet4_eta");
        outTree->addSingle(ak4jet4_phi,  "",  "ak4Jet4_phi");
        outTree->addSingle(ak4jet4_E,  "",  "ak4Jet4_E");
        outTree->addSingle(ak4jet4_csv,  "",  "ak4Jet4_CSV");

        outTree->addSingle(ak8jet1_pt,  "",  "ak8Jet1_pt");
        outTree->addSingle(ak8jet1_eta,  "",  "ak8Jet1_eta");
        outTree->addSingle(ak8jet1_phi,  "",  "ak8Jet1_phi");
        outTree->addSingle(ak8jet1_E,  "",  "ak8Jet1_E");
        outTree->addSingle(ak8jet1_mSD,  "",  "ak8Jet1_msoftdrop");
        outTree->addSingle(ak8jet1_tau1,  "",  "ak8Jet1_tau1");
        outTree->addSingle(ak8jet1_tau2,  "",  "ak8Jet1_tau2");
        outTree->addSingle(ak8jet1_sj0_pt,  "",  "ak8Jet1_subjet0_pt");
        outTree->addSingle(ak8jet1_sj0_eta,  "",  "ak8Jet1_subjet0_eta");
        outTree->addSingle(ak8jet1_sj0_phi,  "",  "ak8Jet1_subjet0_phi");
        outTree->addSingle(ak8jet1_sj0_csv,  "",  "ak8Jet1_subjet0_CSV");
        outTree->addSingle(ak8jet1_sj1_pt,  "",  "ak8Jet1_subjet1_pt");
        outTree->addSingle(ak8jet1_sj1_eta,  "",  "ak8Jet1_subjet1_eta");
        outTree->addSingle(ak8jet1_sj1_phi,  "",  "ak8Jet1_subjet1_phi");
        outTree->addSingle(ak8jet1_sj1_csv,  "",  "ak8Jet1_subjet1_CSV");

        outTree->addSingle(ak8jet2_pt,  "",  "ak8Jet2_pt");
        outTree->addSingle(ak8jet2_eta,  "",  "ak8Jet2_eta");
        outTree->addSingle(ak8jet2_phi,  "",  "ak8Jet2_phi");
        outTree->addSingle(ak8jet2_E,  "",  "ak8Jet2_E");
        outTree->addSingle(ak8jet2_mSD,  "",  "ak8Jet2_msoftdrop");
        outTree->addSingle(ak8jet2_tau1,  "",  "ak8Jet2_tau1");
        outTree->addSingle(ak8jet2_tau2,  "",  "ak8Jet2_tau2");
        outTree->addSingle(ak8jet2_sj0_pt,  "",  "ak8Jet2_subjet0_pt");
        outTree->addSingle(ak8jet2_sj0_eta,  "",  "ak8Jet2_subjet0_eta");
        outTree->addSingle(ak8jet2_sj0_phi,  "",  "ak8Jet2_subjet0_phi");
        outTree->addSingle(ak8jet2_sj0_csv,  "",  "ak8Jet2_subjet0_CSV");
        outTree->addSingle(ak8jet2_sj1_pt,  "",  "ak8Jet2_subjet1_pt");
        outTree->addSingle(ak8jet2_sj1_eta,  "",  "ak8Jet2_subjet1_eta");
        outTree->addSingle(ak8jet2_sj1_phi,  "",  "ak8Jet2_subjet1_phi");
        outTree->addSingle(ak8jet2_sj1_csv,  "",  "ak8Jet2_subjet1_CSV");

        outTree->addSingle(ak8LSjet1_pt,  "",  "ak8lsJet1_pt");
        outTree->addSingle(ak8LSjet1_eta,  "",  "ak8lsJet1_eta");
        outTree->addSingle(ak8LSjet1_phi,  "",  "ak8lsJet1_phi");
        outTree->addSingle(ak8LSjet1_E,  "",  "ak8lsJet1_E");
        outTree->addSingle(ak8LSjet1_mSD,  "",  "ak8lsJet1_msoftdrop");
        outTree->addSingle(ak8LSjet1_tau1,  "",  "ak8lsJet1_tau1");
        outTree->addSingle(ak8LSjet1_tau2,  "",  "ak8lsJet1_tau2");
        outTree->addSingle(ak8LSjet1_sj0_pt,  "",  "ak8lsJet1_subjet0_pt");
        outTree->addSingle(ak8LSjet1_sj0_eta,  "",  "ak8lsJet1_subjet0_eta");
        outTree->addSingle(ak8LSjet1_sj0_phi,  "",  "ak8lsJet1_subjet0_phi");
        outTree->addSingle(ak8LSjet1_sj0_csv,  "",  "ak8lsJet1_subjet0_CSV");
        outTree->addSingle(ak8LSjet1_sj1_pt,  "",  "ak8lsJet1_subjet1_pt");
        outTree->addSingle(ak8LSjet1_sj1_eta,  "",  "ak8lsJet1_subjet1_eta");
        outTree->addSingle(ak8LSjet1_sj1_phi,  "",  "ak8lsJet1_subjet1_phi");
        outTree->addSingle(ak8LSjet1_sj1_csv,  "",  "ak8lsJet1_subjet1_CSV");

        outTree->addSingle(ak8LSjet2_pt,  "",  "ak8lsJet2_pt");
        outTree->addSingle(ak8LSjet2_eta,  "",  "ak8lsJet2_eta");
        outTree->addSingle(ak8LSjet2_phi,  "",  "ak8lsJet2_phi");
        outTree->addSingle(ak8LSjet2_E,  "",  "ak8lsJet2_E");
        outTree->addSingle(ak8LSjet2_mSD,  "",  "ak8lsJet2_msoftdrop");
        outTree->addSingle(ak8LSjet2_tau1,  "",  "ak8lsJet2_tau1");
        outTree->addSingle(ak8LSjet2_tau2,  "",  "ak8lsJet2_tau2");
        outTree->addSingle(ak8LSjet2_sj0_pt,  "",  "ak8lsJet2_subjet0_pt");
        outTree->addSingle(ak8LSjet2_sj0_eta,  "",  "ak8lsJet2_subjet0_eta");
        outTree->addSingle(ak8LSjet2_sj0_phi,  "",  "ak8lsJet2_subjet0_phi");
        outTree->addSingle(ak8LSjet2_sj0_csv,  "",  "ak8lsJet2_subjet0_CSV");
        outTree->addSingle(ak8LSjet2_sj1_pt,  "",  "ak8lsJet2_subjet1_pt");
        outTree->addSingle(ak8LSjet2_sj1_eta,  "",  "ak8lsJet2_subjet1_eta");
        outTree->addSingle(ak8LSjet2_sj1_phi,  "",  "ak8lsJet2_subjet1_phi");
        outTree->addSingle(ak8LSjet2_sj1_csv,  "",  "ak8lsJet2_subjet1_CSV");
    }

    bool passMuPresel(const Muon* mu) {
    	if (mu->pt() < 5) return false;
    	if (mu->absEta() > 2.4) return false;
    	if (fabs(mu->dz()) > 0.1) return false;
    	if (fabs(mu->d0()) > 0.05) return false;
    	if (fabs(mu->sip3D()) > 8) return false;
    	if (!mu->passLooseID()) return false;
    	if (mu->miniIso() > 0.4) return false;
    	return true;
    }

    bool passElPresel(const Electron* el, vector<const Muon*> muons) {
    	if (el->pt() < 7) return false;
    	if (el->absEta() > 2.5) return false;
    	if (fabs(el->dz()) > 0.1) return false;
    	if (fabs(el->d0()) > 0.05) return false;
    	if (fabs(el->sip3D()) > 8) return false;
    	if (!el->passMVA90ID_noIso()) return false;
    	if (el->miniIso() > 0.4) return false;
    	if ((reader_electron->missInnerHits)[el->index()] > 1) return false;

    	bool doesOverlap = false;
    	if (muons.size()) {
        	for (const auto& mu : muons) {
        		if (PhysicsUtilities::deltaR2(*el,*mu) < 0.3*0.3) doesOverlap = true;
        	}
    	}
    	if (doesOverlap) return false;

    	return true;
    }

    bool passAk4Presel(const Jet* jet, vector<const Electron*> els, vector<const Muon*> mus) {
    	if (jet->pt() < 25) return false;
    	if (jet->absEta() > 2.4) return false;
    	if (!jet->passTightID()) return false;

    	bool doesOverlap = false;
    	if (mus.size()) {
        	for (const auto& mu : mus) {
        		if (PhysicsUtilities::deltaR2(*jet,*mu) < 0.4*0.4) doesOverlap = true;
        	}
    	}
    	if (els.size()) {
        	for (const auto& el : els) {
        		if (PhysicsUtilities::deltaR2(*jet,*el) < 0.4*0.4) doesOverlap = true;
        	}
    	}
    	if (doesOverlap) return false;

    	return true;
    }

    bool passAk8Presel(const FatJet* fj, vector<const Electron*> els, vector<const Muon*> mus) {
    	if (!fj->passTightID()) return false;
    	if (fj->pt() < 200) return false;
    	if (fj->absEta() > 2.4) return false;
    	if (fj->nSubJets() != 2) return false;

    	if (fj->subJets()[0].pt() < 20 || fj->subJets()[0].absEta() > 2.4) return false;
    	if (fj->subJets()[1].pt() < 20 || fj->subJets()[1].absEta() > 2.4) return false;
    	if (fj->tau2() / fj->tau1() > 0.75) return false;
        if (fj->sdMom().mass() < 30 || fj->sdMom().mass() > 210) return false;

    	bool doesOverlap = false;
    	if (mus.size()) {
        	for (const auto& mu : mus) {
        		if (PhysicsUtilities::deltaR2(*fj,*mu) < 0.8*0.8) doesOverlap = true;
        	}
    	}

    	if (els.size()) {
        	for (const auto& el : els) {
        		if (PhysicsUtilities::deltaR2(*fj,*el) < 0.8*0.8) doesOverlap = true;
        	}
    	}

    	if (doesOverlap) return false;

    	return true;
    }

    bool passAk8LSPresel(const FatJet* fj, vector<const Electron*> els, vector<const Muon*> mus) {
    	if (!fj->passLooseID()) return false;
    	if (fj->pt() < 100) return false;
    	if (fj->absEta() > 2.4) return false;
    	if (fj->nSubJets() != 2) return false;

    	if (fj->subJets()[0].pt() < 20 || fj->subJets()[0].absEta() > 2.4) return false;
    	if (fj->subJets()[1].pt() < 20 || fj->subJets()[1].absEta() > 2.4) return false;
    	if (fj->tau2() / fj->tau1() > 0.75) return false;
        if (fj->sdMom().mass() < 30 || fj->sdMom().mass() > 210) return false;

    	bool doesOverlap = false;
    	if (mus.size()) {
        	for (const auto& mu : mus) {
        		if (PhysicsUtilities::deltaR2(*fj,*mu) < 1.2*1.2) doesOverlap = true;
        	}
    	}
    	if (els.size()) {
        	for (const auto& el : els) {
        		if (PhysicsUtilities::deltaR2(*fj,*el) < 1.2*1.2) doesOverlap = true;
        	}
    	}

    	if (!doesOverlap) return false;

    	return true;
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
//        if(!passEventFilters) return false;

        vector<const Electron*> electrons;
        vector<const Muon*> muons;
        vector<const Jet*> ak4jets;
        vector<const FatJet*> ak8jets;
        vector<const FatJet*> ak8LSjets;

        if (reader_muon) {
            for (const auto& mu : reader_muon->muons) {
            	if (passMuPresel(&mu)) muons.push_back(&mu);
            }
        }

        if (reader_electron) {
            for (const auto& el : reader_electron->electrons) {
            	if (passElPresel(&el,muons)) electrons.push_back(&el);
            }
        }

        if (reader_jet_chs) {
            for (const auto& j : reader_jet_chs->jets) {
            	if (passAk4Presel(&j,electrons,muons)) ak4jets.push_back(&j);
            }
        }

        if (reader_fatjet) {
            for (const auto& j : reader_fatjet->jets) {
            	if (passAk8Presel(&j,electrons,muons)) ak8jets.push_back(&j);
            }
        }

        if (reader_fatjet_noLep) {
            for (const auto& j : reader_fatjet_noLep->jets) {
            	if (passAk8LSPresel(&j,electrons,muons)) ak8LSjets.push_back(&j);
            }
        }

        std::sort(muons.begin(),muons.end(), PhysicsUtilities::greaterPTDeref<Muon>());
        std::sort(electrons.begin(),electrons.end(), PhysicsUtilities::greaterPTDeref<Electron>());
        std::sort(ak4jets.begin(),ak4jets.end(), PhysicsUtilities::greaterPTDeref<Jet>());
        std::sort(ak8jets.begin(),ak8jets.end(), PhysicsUtilities::greaterPTDeref<FatJet>());
        std::sort(ak8LSjets.begin(),ak8LSjets.end(), PhysicsUtilities::greaterPTDeref<FatJet>());

        event = size64(reader_event->event.val());
        lumiSec = size(reader_event->lumi.val());

        n_mu = size(muons.size());
        n_el = size(electrons.size());
        n_ak4 = size(ak4jets.size());
        n_ak8 = size(ak8jets.size());
        n_ak8LS = size(ak8LSjets.size());

        pu_weight = float(puSFProc->getCorrection(*reader_event->nTruePUInts,CorrHelp::NOMINAL));
        mc_weight = float(*reader_event->genWeight);

        met_pt = reader_event->met.pt();
        met_phi = reader_event->met.phi();


        if (muons.size()) {
        	numE_preselMu++;

        	const auto mu1 = muons[0];
            mu1_pt = mu1->pt();
            mu1_eta = mu1->eta();
            mu1_phi = mu1->phi();
            mu1_E = mu1->E();
            mu1_charge = mu1->q();
            mu1_miniIso = mu1->miniIso();
            mu1_sip = mu1->sip3D();
            mu1_dz = mu1->dz();
            mu1_dxy = mu1->d0();
            mu1_dxyAbs = fabs(mu1->d0());

            if (muons.size() >= 2) {
            	const auto mu2 = muons[1];
                mu2_pt = mu2->pt();
                mu2_eta = mu2->eta();
                mu2_phi = mu2->phi();
                mu2_E = mu2->E();
                mu2_charge = mu2->q();
                mu2_miniIso = mu2->miniIso();
                mu2_sip = mu2->sip3D();
                mu2_dz = mu2->dz();
                mu2_dxy = mu2->d0();
                mu2_dxyAbs = fabs(mu2->d0());
            }
        }

        if (electrons.size()) {
        	numE_preselEl++;

        	const auto el1 = electrons[0];
            el1_pt = el1->pt();
            el1_eta = el1->eta();
            el1_phi = el1->phi();
            el1_E = el1->E();
            el1_charge = el1->q();
            el1_miniIso = el1->miniIso();
            el1_sip = el1->sip3D();
            el1_dz = el1->dz();
            el1_dxy = el1->d0();
            el1_dxyAbs = fabs(el1->d0());
            el1_MVAID = el1->passMVA90ID_noIso();

            if (electrons.size() >= 2) {
            	const auto el2 = electrons[1];
                el2_pt = el2->pt();
                el2_eta = el2->eta();
                el2_phi = el2->phi();
                el2_E = el2->E();
                el2_charge = el2->q();
                el2_miniIso = el2->miniIso();
                el2_sip = el2->sip3D();
                el2_dz = el2->dz();
                el2_dxy = el2->d0();
                el2_dxyAbs = fabs(el2->d0());
                el2_MVAID = el2->passMVA90ID_noIso();

            }
        }

        if (ak4jets.size()) {
        	numE_preselAk4++;

            ak4jet1_pt = ak4jets[0]->pt();
            ak4jet1_eta = ak4jets[0]->eta();
            ak4jet1_phi = ak4jets[0]->phi();
            ak4jet1_E = ak4jets[0]->E();
            ak4jet1_csv = ak4jets[0]->deep_csv();
            if (ak4jets.size() >= 2) {
                ak4jet2_pt = ak4jets[1]->pt();
                ak4jet2_eta = ak4jets[1]->eta();
                ak4jet2_phi = ak4jets[1]->phi();
                ak4jet2_E = ak4jets[1]->E();
                ak4jet2_csv = ak4jets[1]->deep_csv();
            }
            if (ak4jets.size() >= 3) {
                ak4jet3_pt = ak4jets[2]->pt();
                ak4jet3_eta = ak4jets[2]->eta();
                ak4jet3_phi = ak4jets[2]->phi();
                ak4jet3_E = ak4jets[2]->E();
                ak4jet3_csv = ak4jets[2]->deep_csv();
            }
            if (ak4jets.size() >= 4) {
                ak4jet4_pt = ak4jets[3]->pt();
                ak4jet4_eta = ak4jets[3]->eta();
                ak4jet4_phi = ak4jets[3]->phi();
                ak4jet4_E = ak4jets[3]->E();
                ak4jet4_csv = ak4jets[3]->deep_csv();
            }
        }

        if (ak8jets.size()) {
        	numE_preselAk8++;

            ak8jet1_pt = ak8jets[0]->pt();
            ak8jet1_eta = ak8jets[0]->eta();
            ak8jet1_phi = ak8jets[0]->phi();
            ak8jet1_E = ak8jets[0]->E();
            ak8jet1_mSD = ak8jets[0]->sdMom().mass();
            ak8jet1_tau1 = ak8jets[0]->tau1();
            ak8jet1_tau2 = ak8jets[0]->tau2();
            ak8jet1_sj0_pt = ak8jets[0]->subJets()[0].pt();
            ak8jet1_sj0_eta = ak8jets[0]->subJets()[0].eta();
            ak8jet1_sj0_phi = ak8jets[0]->subJets()[0].phi();
            ak8jet1_sj0_csv = ak8jets[0]->subJets()[0].deep_csv();
            ak8jet1_sj1_pt = ak8jets[0]->subJets()[1].pt();
            ak8jet1_sj1_eta = ak8jets[0]->subJets()[1].eta();
            ak8jet1_sj1_phi = ak8jets[0]->subJets()[1].phi();
            ak8jet1_sj1_csv = ak8jets[0]->subJets()[1].deep_csv();

            if (ak8jets.size() >= 2) {
                ak8jet2_pt = ak8jets[1]->pt();
                ak8jet2_eta = ak8jets[1]->eta();
                ak8jet2_phi = ak8jets[1]->phi();
                ak8jet2_E = ak8jets[1]->E();
                ak8jet2_mSD = ak8jets[1]->sdMom().mass();
                ak8jet2_tau1 = ak8jets[1]->tau1();
                ak8jet2_tau2 = ak8jets[1]->tau2();
                ak8jet2_sj0_pt = ak8jets[1]->subJets()[0].pt();
                ak8jet2_sj0_eta = ak8jets[1]->subJets()[0].eta();
                ak8jet2_sj0_phi = ak8jets[1]->subJets()[0].phi();
                ak8jet2_sj0_csv = ak8jets[1]->subJets()[0].deep_csv();
                ak8jet2_sj1_pt = ak8jets[1]->subJets()[1].pt();
                ak8jet2_sj1_eta = ak8jets[1]->subJets()[1].eta();
                ak8jet2_sj1_phi = ak8jets[1]->subJets()[1].phi();
                ak8jet2_sj1_csv = ak8jets[1]->subJets()[1].deep_csv();
            }
        }

        if (ak8LSjets.size()) {
            ak8LSjet1_pt = ak8LSjets[0]->pt();
            ak8LSjet1_eta = ak8LSjets[0]->eta();
            ak8LSjet1_phi = ak8LSjets[0]->phi();
            ak8LSjet1_E = ak8LSjets[0]->E();
            ak8LSjet1_mSD = ak8LSjets[0]->sdMom().mass();
            ak8LSjet1_tau1 = ak8LSjets[0]->tau1();
            ak8LSjet1_tau2 = ak8LSjets[0]->tau2();
            ak8LSjet1_sj0_pt = ak8LSjets[0]->subJets()[0].pt();
            ak8LSjet1_sj0_eta = ak8LSjets[0]->subJets()[0].eta();
            ak8LSjet1_sj0_phi = ak8LSjets[0]->subJets()[0].phi();
            ak8LSjet1_sj0_csv = ak8LSjets[0]->subJets()[0].deep_csv();
            ak8LSjet1_sj1_pt = ak8LSjets[0]->subJets()[1].pt();
            ak8LSjet1_sj1_eta = ak8LSjets[0]->subJets()[1].eta();
            ak8LSjet1_sj1_phi = ak8LSjets[0]->subJets()[1].phi();
            ak8LSjet1_sj1_csv = ak8LSjets[0]->subJets()[1].deep_csv();

            if (ak8LSjets.size() >= 2) {
                ak8LSjet2_pt = ak8LSjets[1]->pt();
                ak8LSjet2_eta = ak8LSjets[1]->eta();
                ak8LSjet2_phi = ak8LSjets[1]->phi();
                ak8LSjet2_E = ak8LSjets[1]->E();
                ak8LSjet2_mSD = ak8LSjets[1]->sdMom().mass();
                ak8LSjet2_tau1 = ak8LSjets[1]->tau1();
                ak8LSjet2_tau2 = ak8LSjets[1]->tau2();
                ak8LSjet2_sj0_pt = ak8LSjets[1]->subJets()[0].pt();
                ak8LSjet2_sj0_eta = ak8LSjets[1]->subJets()[0].eta();
                ak8LSjet2_sj0_phi = ak8LSjets[1]->subJets()[0].phi();
                ak8LSjet2_sj0_csv = ak8LSjets[1]->subJets()[0].deep_csv();
                ak8LSjet2_sj1_pt = ak8LSjets[1]->subJets()[1].pt();
                ak8LSjet2_sj1_eta = ak8LSjets[1]->subJets()[1].eta();
                ak8LSjet2_sj1_phi = ak8LSjets[1]->subJets()[1].phi();
                ak8LSjet2_sj1_csv = ak8LSjets[1]->subJets()[1].deep_csv();
            }
        }

        return true;
    }

    //Event information
    size64 event       = 0;
    size lumiSec     = 0;
    size8 run         = 1;
    size n_mu = 0;
    size n_el = 0;
    size n_ak4     = 0;
    size n_ak8    = 0;
    size n_ak8LS = 0;
    float pu_weight=0;
    float mc_weight=0;

    //Muon info
    float mu1_pt = 0;
    float mu1_eta = 0;
    float mu1_phi = 0;
    float mu1_E = 0;
    float mu1_charge = 0;
    float mu1_miniIso = 0;
    float mu1_sip = 0;
    float mu1_dz = 0;
    float mu1_dxy = 0;
    float mu1_dxyAbs = 0;
    float mu2_pt = 0;
    float mu2_eta = 0;
    float mu2_phi = 0;
    float mu2_E = 0;
    float mu2_charge = 0;
    float mu2_miniIso = 0;
    float mu2_sip = 0;
    float mu2_dz = 0;
    float mu2_dxy = 0;
    float mu2_dxyAbs = 0;

    //Electron info
    float el1_pt = 0;
    float el1_eta = 0;
    float el1_phi = 0;
    float el1_E = 0;
    float el1_charge = 0;
    float el1_miniIso = 0;
    float el1_sip = 0;
    float el1_dz = 0;
    float el1_dxy = 0;
    float el1_dxyAbs = 0;
    size8 el1_MVAID = 0;
    float el2_pt = 0;
    float el2_eta = 0;
    float el2_phi = 0;
    float el2_E = 0;
    float el2_charge = 0;
    float el2_miniIso = 0;
    float el2_sip = 0;
    float el2_dz = 0;
    float el2_dxy = 0;
    float el2_dxyAbs = 0;
    size8 el2_MVAID = 0;

    //AK4 jet info
    float ak4jet1_pt = 0;
    float ak4jet1_eta = 0;
    float ak4jet1_phi = 0;
    float ak4jet1_E = 0;
    float ak4jet1_csv = 0;
    float ak4jet2_pt = 0;
    float ak4jet2_eta = 0;
    float ak4jet2_phi = 0;
    float ak4jet2_E = 0;
    float ak4jet2_csv = 0;
    float ak4jet3_pt = 0;
    float ak4jet3_eta = 0;
    float ak4jet3_phi = 0;
    float ak4jet3_E = 0;
    float ak4jet3_csv = 0;
    float ak4jet4_pt = 0;
    float ak4jet4_eta = 0;
    float ak4jet4_phi = 0;
    float ak4jet4_E = 0;
    float ak4jet4_csv = 0;

    //AK8 jet info
    float ak8jet1_pt=0;
    float ak8jet1_eta=0;
    float ak8jet1_phi=0;
    float ak8jet1_E=0;
    float ak8jet1_mSD=0;
    float ak8jet1_tau1=0;
    float ak8jet1_tau2=0;
    float ak8jet1_sj0_pt=0;
    float ak8jet1_sj0_eta=0;
    float ak8jet1_sj0_phi=0;
    float ak8jet1_sj0_csv=0;
    float ak8jet1_sj1_pt=0;
    float ak8jet1_sj1_eta=0;
    float ak8jet1_sj1_phi=0;
    float ak8jet1_sj1_csv=0;

    float ak8jet2_pt=0;
    float ak8jet2_eta=0;
    float ak8jet2_phi=0;
    float ak8jet2_E=0;
    float ak8jet2_mSD=0;
    float ak8jet2_tau1=0;
    float ak8jet2_tau2=0;
    float ak8jet2_sj0_pt=0;
    float ak8jet2_sj0_eta=0;
    float ak8jet2_sj0_phi=0;
    float ak8jet2_sj0_csv=0;
    float ak8jet2_sj1_pt=0;
    float ak8jet2_sj1_eta=0;
    float ak8jet2_sj1_phi=0;
    float ak8jet2_sj1_csv=0;

    //AK8 LS jet info
    float ak8LSjet1_pt=0;
    float ak8LSjet1_eta=0;
    float ak8LSjet1_phi=0;
    float ak8LSjet1_E=0;
    float ak8LSjet1_mSD=0;
    float ak8LSjet1_tau1=0;
    float ak8LSjet1_tau2=0;
    float ak8LSjet1_sj0_pt=0;
    float ak8LSjet1_sj0_eta=0;
    float ak8LSjet1_sj0_phi=0;
    float ak8LSjet1_sj0_csv=0;
    float ak8LSjet1_sj1_pt=0;
    float ak8LSjet1_sj1_eta=0;
    float ak8LSjet1_sj1_phi=0;
    float ak8LSjet1_sj1_csv=0;

    float ak8LSjet2_pt=0;
    float ak8LSjet2_eta=0;
    float ak8LSjet2_phi=0;
    float ak8LSjet2_E=0;
    float ak8LSjet2_mSD=0;
    float ak8LSjet2_tau1=0;
    float ak8LSjet2_tau2=0;
    float ak8LSjet2_sj0_pt=0;
    float ak8LSjet2_sj0_eta=0;
    float ak8LSjet2_sj0_phi=0;
    float ak8LSjet2_sj0_csv=0;
    float ak8LSjet2_sj1_pt=0;
    float ak8LSjet2_sj1_eta=0;
    float ak8LSjet2_sj1_phi=0;
    float ak8LSjet2_sj1_csv=0;

    //MET
    float met_pt=0;
    float met_phi=0;

    // counting (not for trees)
    int numE_preselMu = 0;
    int numE_preselEl = 0;
    int numE_preselAk4 = 0;
    int numE_preselAk8 = 0;

};




void doOne(const std::string& fileName, const int treeInt,int randSeed,  const std::string& outFileName, float xSec = -1, float numEvent = -1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    if(xSec > 0) a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();

    printf("numE_preselMu = %d\n",a.numE_preselMu);
    printf("numE_preselEl = %d\n",a.numE_preselEl);
    printf("numE_preselAk4 = %d\n",a.numE_preselAk4);
    printf("numE_preselAk8 = %d\n",a.numE_preselAk8);

}

#endif

void makeSyncTrees(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    doOne(fileName,treeInt,randSeed,outFileName,xSec,numEvent);

}
