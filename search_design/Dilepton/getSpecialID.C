
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
        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);
    }

    bool passGEN_requirements(const GenParticle* lep, const GenParticle* q1, const GenParticle* q2) {
    	const GenParticle* b1 = diHiggsEvt.b1;
    	const GenParticle* b2 = diHiggsEvt.b2;

    	if(b1 && b2 && q1 && q2 && lep) {
        	if(!ParticleInfo::isQuark(q1->pdgId()) || !ParticleInfo::isQuark(q2->pdgId())) {
        		cout << "one of the quarks is not actually a quark" << endl;
        		cout << "q1 and q2 = " << q1->pdgId() << ", " << q2->pdgId() << endl;
        		cout << "lep and nu = " << lep->pdgId() << ", " << diHiggsEvt.w1_d2->pdgId() << endl << endl;
        		return false;
        	}
        	if(!ParticleInfo::isLepton(lep->pdgId())) {
        		cout << "lep is not a lepton" << endl;
        		return false;
        	}

        	if (PhysicsUtilities::deltaR2(*b1,*b2) > 0.8*0.8) return false;
        	if (PhysicsUtilities::deltaR2(*q1,*q2) > 0.8*0.8) return false;
        	if (PhysicsUtilities::deltaR2((q1->p4()+q2->p4()),*lep) > 1.2*1.2) return false;
        	return true;

    	} else return false;
    }

    void plot(TString sn, const Lepton *lep1, const Lepton *lep2, const FatJet *hbb) {
    	if (!hbb) return;
    	MomentumF dilepMOM = lep1->p4() + lep2->p4();
    	double dR_ll = PhysicsUtilities::deltaR(*lep1,*lep2);

    	auto getHHMass = [&] (const FatJet* hbbmom, MomentumF& llMom) {
    		double pz = reader_event->met.pt() / TMath::Tan(llMom.theta());
    		if ((pz>0) != (llMom.pz()>0)) pz = (-1.0)*pz;
    		ASTypes::CartLorentzVector pnunu(reader_event->met.px(),reader_event->met.py(),pz,
    				sqrt(pow(reader_event->met.px(),2)+pow(reader_event->met.py(),2)+pz*pz+40*40));

    		double mass = (hbbmom->p4() + llMom.p4() + pnunu).mass();
    		return mass;
    	};

    	double mhh = getHHMass(hbb,dilepMOM);

    	auto plt = [&](TString pre) {
        	plotter.getOrMake1DPre(pre,"sip",";sip",500,0,100)->Fill(lep1->sip3D(),weight);
        	plotter.getOrMake1DPre(pre,"sip",";sip",500,0,100)->Fill(lep2->sip3D(),weight);
        	plotter.getOrMake1DPre(pre,"sip_1",";sip",500,0,100)->Fill(lep1->sip3D(),weight);
        	plotter.getOrMake1DPre(pre,"sip_2",";sip",500,0,100)->Fill(lep2->sip3D(),weight);
    	};

    	plt(sn);

    	if (mhh < 700) return;
    	if (dR_ll > 1.6) return;
    	if (reader_event->met.pt() < 40) return;
    	if (dilepMOM.mass() < 12 || dilepMOM.mass() > 75) return;

        const auto jetsVeto = PhysicsUtilities::selObjsD(jets,
                [&](const Jet* j){return  PhysicsUtilities::deltaR2(*j,*hbb ) >= 1.2*1.2;});
        const auto nBtags = PhysicsUtilities::selObjsMomD(jetsVeto,
                parameters.jets.minBtagJetPT,parameters.jets.maxBTagJetETA,
                [&](const Jet* j){return BTagging::passJetBTagWP(parameters.jets,*j);} ).size();

    	if (nBtags > 0) return;

    	plt(sn+"_pass2lSel");
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        TString prefix = smpName;

        plotter.getOrMake1DPre(smpName,"nTruePUints","",100,0,100)->Fill(reader_event->nTruePUInts.val(),weight);

        vector<double> sips = {4,5,6,7,8,9,10,11,12,13,14,15,17,20,25,30,35,40,9999};
        if (!isSignal()) {
        	DileptonParameters param = parameters.dileptons;

        	for (const auto& sip : sips) {
    			param.mu_maxSip3D = sip;
    			param.mu_maxSip3D = sip;

    	    	const auto leps = DileptonProcessor::getLeptons(param,*reader_muon,*reader_electron);
    	    	if (leps.size() == 2) {
    	    		const Lepton *lep1 = leps.front();
    	    		const Lepton *lep2 = leps.back();
    	    		float pt = lep1->isMuon() ? 27 : 30;
    	    		if (lep1->pt() < pt) continue;

    	    		proc.reset(new FatJetProcessor ());
    	    		proc->loadDilepFatJet(parameters.fatJets,*reader_fatjet,lep1,lep2);
    	    		const FatJet *bbjet = fjProc->getDilepHbbCand();

        			if (bbjet) {
        				float bbmass = bbjet->sdMom().mass();
        				if (bbmass > 30 && bbmass < 210) {
        					plot(prefix+"_stdLepSel",lep1,lep2,bbjet);
        				}
        			}
    	    	}
        	}
        } else {

        }


        return true;
    }

    std::unique_ptr<FatJetProcessor>        proc     ;
    size64 triggerAccepts = 0;
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getScrap(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
