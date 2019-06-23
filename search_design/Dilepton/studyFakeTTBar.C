
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
#include "TreeReaders/interface/JetReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/BTagging.h"

using namespace TAna;
using namespace std;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);
    }

    const Lepton* getFakeLep(const Lepton* lep1, const Lepton* lep2) {
    	const Lepton* unmatchedLep;
    	if (smDecayEvt.topDecays.size() != 2) return 0;

        const GenParticle *genLep;
        const GenParticle *firstLep;

        for (const auto& t : smDecayEvt.topDecays) {
        	if (t.type < TopDecay::TAU_MU) continue;
        	if (ParticleInfo::isLepton(t.W_decay.dau1->pdgId())) {
        		firstLep = t.W_decay.dau1;
        		break;
        	}
        }
        if (!firstLep) return 0;

    	if (firstLep->absPdgId() == ParticleInfo::p_muminus || firstLep->absPdgId() == ParticleInfo::p_eminus) genLep = firstLep;
    	else {
    		for (unsigned int k=0; k < firstLep->numberOfDaughters(); k++) {
    			if (!ParticleInfo::isLepton(firstLep->daughter(k)->pdgId())) continue;
    			genLep = firstLep->daughter(k);
    		}
    	}
    	if (!genLep) return 0;

    	float dr1 = PhysicsUtilities::deltaR(*genLep,*lep1);
    	float dr2 = PhysicsUtilities::deltaR(*genLep,*lep2);

    	unmatchedLep = (dr1 < dr2) ? lep2 : lep1;

    	return unmatchedLep;
    }

    void plot(TString sn, const Lepton* lep1, const Lepton* lep2) {

    	if (smDecayEvt.topDecays.size() != 2) return;
    	const GenParticle *b1 = smDecayEvt.topDecays[0].b;
    	const GenParticle *b2 = smDecayEvt.topDecays[1].b;

        double drll = PhysicsUtilities::deltaR(*lep1,*lep2);
        double mll  = (lep1->p4()+lep2->p4()).mass();
        double dphi_metll = PhysicsUtilities::absDeltaPhi(reader_event->met.p4(),(lep1->p4()+lep2->p4()));

    	double dr11 = PhysicsUtilities::deltaR(*lep1,*b1);
    	double dr12 = PhysicsUtilities::deltaR(*lep1,*b2);
    	double dr21 = PhysicsUtilities::deltaR(*lep2,*b1);
    	double dr22 = PhysicsUtilities::deltaR(*lep2,*b2);

    	int nCloseLeps2 = (dr11 < 0.2) + (dr12 < 0.2) + (dr21 < 0.2) + (dr22 < 0.2);
    	int nCloseLeps1 = (dr11 < 0.1) + (dr12 < 0.1) + (dr21 < 0.1) + (dr22 < 0.1);

        plotter.getOrMake1DPre(sn,"dr",";dr",500,0,3)->Fill(dr11,weight);
        plotter.getOrMake1DPre(sn,"dr",";dr",500,0,3)->Fill(dr12,weight);
        plotter.getOrMake1DPre(sn,"dr",";dr",500,0,3)->Fill(dr21,weight);
        plotter.getOrMake1DPre(sn,"dr",";dr",500,0,3)->Fill(dr22,weight);
        plotter.getOrMake1DPre(sn,"nClose1",";num",4,-0.5,3.5)->Fill(nCloseLeps1,weight);
        plotter.getOrMake1DPre(sn,"nClose2",";num",4,-0.5,3.5)->Fill(nCloseLeps2,weight);

        if (drll > 1.6) return;
        if (mll < 12 || mll > 75) return;
        if (dphi_metll > TMath::PiOver2()) return;
        if (nMedBTags_HbbV_2l > 0) return;

        sn += "_pass2lSel";
        plotter.getOrMake1DPre(sn,"dr",";dr",500,0,3)->Fill(dr11,weight);
        plotter.getOrMake1DPre(sn,"dr",";dr",500,0,3)->Fill(dr12,weight);
        plotter.getOrMake1DPre(sn,"dr",";dr",500,0,3)->Fill(dr21,weight);
        plotter.getOrMake1DPre(sn,"dr",";dr",500,0,3)->Fill(dr22,weight);
        plotter.getOrMake1DPre(sn,"nClose1",";num",4,-0.5,3.5)->Fill(nCloseLeps1,weight);
        plotter.getOrMake1DPre(sn,"nClose2",";num",4,-0.5,3.5)->Fill(nCloseLeps2,weight);

        const GenParticle *genLep;
        for (const auto& t : smDecayEvt.topDecays) {
        	if (t.type < TopDecay::MU) continue;
        	if (ParticleInfo::isLepton(t.W_decay.dau1->pdgId())) {
        		genLep = t.W_decay.dau1;
        		break;
        	}
        }
        double dr1 = PhysicsUtilities::deltaR(*lep1,*genLep);
        double dr2 = PhysicsUtilities::deltaR(*lep2,*genLep);

        const Lepton* unmatchedLep = (dr1 < dr2) ? lep2 : lep1;
        double dr_b1 = PhysicsUtilities::deltaR(*smDecayEvt.topDecays[0].b,*unmatchedLep);
        double dr_b2 = PhysicsUtilities::deltaR(*smDecayEvt.topDecays[1].b,*unmatchedLep);

        double drOfInterest = (dr_b1 < dr_b2) ? dr_b1 : dr_b2;
        plotter.getOrMake1DPre(sn,"dr_fake_b",";dr",500,0,2)->Fill(drOfInterest,weight);

    }

    const Jet* getMatchedJet(const GenParticle* b) {

    	if (!reader_jet) return 0;
    	float dr = 9999;
    	int idx = -1;
    	for (const auto& j : reader_jet->jets) {

    		if (PhysicsUtilities::deltaR2(j,*b) < dr*dr) {
    			dr = PhysicsUtilities::deltaR(j,*b);
    			idx = j.index();
    		}
    	}

    	if (idx < 0) return 0;
    	if (dr > 0.4) return 0;
    	return &reader_jet->jets[idx];
    }

    void studyBtagging(TString sn, const Lepton* lep1, const Lepton* lep2) {

    	const Lepton *fakelep = getFakeLep(lep1,lep2);
    	const auto b1 = smDecayEvt.topDecays[0].b;
    	const auto b2 = smDecayEvt.topDecays[1].b;

    	if (b1->pt() < 25 || b2->pt() < 25) return;
    	const Jet* bjet1 = getMatchedJet(b1);
    	const Jet* bjet2 = getMatchedJet(b2);
    	if (!bjet1 || !bjet2) return;
    	if (bjet1->pt() < 25 || bjet2->pt() < 25) return;

		plotter.getOrMake1DPre(sn+"base","pt_fakelep",";pt",50,0,500)->Fill(fakelep->pt(),weight);

		if (BTagging::passJetBTagWP(parameters.jets,*bjet1)) {
		    plotter.getOrMake1DPre(sn+"base_passWP","pt_fakelep",";pt",50,0,500)->Fill(fakelep->pt(),weight);
		}
		if (BTagging::passJetBTagWP(parameters.jets,*bjet2)) {
		    plotter.getOrMake1DPre(sn+"base_passWP","pt_fakelep",";pt",50,0,500)->Fill(fakelep->pt(),weight);
		}

    	double drFakeb1 = PhysicsUtilities::deltaR(*bjet1,*fakelep);
    	double drFakeb2 = PhysicsUtilities::deltaR(*bjet2,*fakelep);

    	if (drFakeb1 < 0.4 && drFakeb2 < 0.4) {
    		plotter.getOrMake1DPre(sn+"doubleMatched","pt_fakelep",";#DeltaR",50,0,500)->Fill(fakelep->pt(),weight);
    		return;
    	}

    	if (drFakeb1 < 0.4) {
    		plotter.getOrMake1DPre(sn+"bMatched","dr_FakeB",";#DeltaR",20,0,0.4)->Fill(drFakeb1,weight);
    		plotter.getOrMake1DPre(sn+"bMatched","pt_fakelep",";pt",50,0,500)->Fill(fakelep->pt(),weight);

    		if (BTagging::passJetBTagWP(parameters.jets,*bjet1)) {
        		plotter.getOrMake1DPre(sn+"bMatched_passWP","dr_FakeB",";#DeltaR",20,0,0.4)->Fill(drFakeb1,weight);
        		plotter.getOrMake1DPre(sn+"bMatched_passWP","pt_fakelep",";pt",50,0,500)->Fill(fakelep->pt(),weight);

    		}
    	} else if (drFakeb2 < 0.4) {
    		plotter.getOrMake1DPre(sn+"bMatched","dr_FakeB",";#DeltaR",20,0,0.4)->Fill(drFakeb2,weight);
    		plotter.getOrMake1DPre(sn+"bMatched","pt_fakelep",";pt",50,0,500)->Fill(fakelep->pt(),weight);

    		if (BTagging::passJetBTagWP(parameters.jets,*bjet2)) {
        		plotter.getOrMake1DPre(sn+"bMatched_passWP","dr_FakeB",";#DeltaR",20,0,0.4)->Fill(drFakeb2,weight);
        		plotter.getOrMake1DPre(sn+"bMatched_passWP","pt_fakelep",";pt",50,0,500)->Fill(fakelep->pt(),weight);

    		}
    	}


    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        TString prefix = smpName;

        if (mcProc != FillerConstants::TTBAR) return false;
        if (nLepsTT != 1) return false;
        if (smDecayEvt.topDecays.size() != 2) return false;

        if (selectedDileptons.size() != 2) return false;
        if (!hbbCand_2l) return false;
        float ptthresh = (selectedDileptons[0]->isMuon() ? 27 : 30);
        if (selectedDileptons[0]->pt() < ptthresh) return false;

        if (hbbCSVCat_2l < BTagging::CSVSJ_MF) return false;
        if (hbbMass_2l < 30 || hbbMass_2l > 210) return false;
        if (hh_2l.mass() < 700) return false;

        plot(prefix,selectedDileptons.front(), selectedDileptons.back());
        studyBtagging(prefix+"_",selectedDileptons.front(),selectedDileptons.back());


        return true;
    }

    std::unique_ptr<FatJetProcessor>        proc     ;
    size64 triggerAccepts = 0;
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void studyFakeTTBar(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
