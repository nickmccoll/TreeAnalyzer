
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

    void plot(TString sn, const GenParticle* lep1, const GenParticle* lep2) {

    	TString l1S = lep1->absPdgId() == 13 ? "mu" : "e";
    	TString l2S = lep2->absPdgId() == 13 ? "mu" : "e";

        plotter.getOrMake1DPre(sn,"pt1",";pt",3000,0,3000)->Fill(lep1->pt(),weight);
        plotter.getOrMake1DPre(sn,"pt2",";pt",3000,0,3000)->Fill(lep2->pt(),weight);
        plotter.getOrMake1DPre(sn,l1S+"_pt1",";pt",3000,0,3000)->Fill(lep1->pt(),weight);
        plotter.getOrMake1DPre(sn,l2S+"_pt2",";pt",3000,0,3000)->Fill(lep2->pt(),weight);
        plotter.getOrMake1DPre(sn,"eta1",";eta",500,-5,5)->Fill(lep1->eta(),weight);
        plotter.getOrMake1DPre(sn,"eta2",";eta",500,-5,5)->Fill(lep2->eta(),weight);
        plotter.getOrMake1DPre(sn,l1S+"_eta1",";eta",500,-5,5)->Fill(lep1->eta(),weight);
        plotter.getOrMake1DPre(sn,l2S+"_eta2",";eta",500,-5,5)->Fill(lep2->eta(),weight);

        double drll = PhysicsUtilities::deltaR(*lep1,*lep2);
        double dphill = PhysicsUtilities::absDeltaPhi(*lep1,*lep2);
        double mll  = (lep1->p4()+lep2->p4()).mass();

        plotter.getOrMake1DPre(sn,"drll",";#DeltaR",100,0,5)->Fill(drll,weight);
        plotter.getOrMake1DPre(sn,"dphill",";#Delta#phi",100,0,3.14)->Fill(dphill,weight);
        plotter.getOrMake1DPre(sn,"mll",";M_{ll}",50,0,200)->Fill(mll,weight);

        MomentumF ll(lep1->p4()+lep2->p4());
        MomentumF hww(diHiggsEvt.hww->p4());
        ASTypes::CartLorentzVector invMom(hww.px()-ll.px(),hww.py()-ll.py(),hww.pz()-ll.pz(),hww.E()-ll.E());
        MomentumF inv(invMom);

        double dphibl1 = PhysicsUtilities::absDeltaPhi(*diHiggsEvt.hbb,*lep1);
        double dphibl2 = PhysicsUtilities::absDeltaPhi(*diHiggsEvt.hbb,*lep2);
        double dphibll = PhysicsUtilities::absDeltaPhi(*diHiggsEvt.hbb,ll);

        plotter.getOrMake1DPre(sn,"minv",";M_{inv}",50,0,200)->Fill(inv.mass(),weight);
        plotter.getOrMake1DPre(sn,"dphibl1",";#Delta#Phi",50,0,TMath::Pi())->Fill(dphibl1,weight);
        plotter.getOrMake1DPre(sn,"dphibl2",";#Delta#Phi",50,0,TMath::Pi())->Fill(dphibl2,weight);
        plotter.getOrMake1DPre(sn,"dphibll",";#Delta#Phi",50,0,TMath::Pi())->Fill(dphibll,weight);

        double dtheta = ll.theta() - inv.theta();
        plotter.getOrMake1DPre(sn,"dthetallinv",";#Delta#theta",500,-5,5)->Fill(dtheta,weight);

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

    bool passFullRecoSel() {
	    if(!EventSelection::passTriggerSuite2017(*reader_event)) return false;
	    if(ht_puppi < 400) return false;
	    if(selectedDileptons.size() != 2) return false;
	    if(!hbbCand_2l) return false;
	    if(hbbCSVCat_2l < BTagging::CSVSJ_MF) return false;
	    if(hbbMass_2l < 30 || hbbMass_2l > 210) return false;
	    if(dilep1->pt() < (dilep1->isMuon() ? 27 : 30)) return false;
	    if(nMedBTags_HbbV_2l != 0) return false;

	    double mll = (dilep1->p4()+dilep2->p4()).mass();
	    double drll = PhysicsUtilities::deltaR(*dilep1,*dilep2);
	    double dphillmet = PhysicsUtilities::absDeltaPhi(reader_event->met,dilep1->p4()+dilep2->p4());

	    if (mll < 6 || mll > 75) return false;
	    if (drll > 1.0) return false;
	    if (dphillmet > TMath::Pi()) return false;
	    if (reader_event->met.pt() / hh_2l.mass() < 0.1) return false;

	    return true;
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        TString prefix = smpName;

        if (!isSignal()) return false;
        if (diHiggsEvt.type != DiHiggsEvent::DILEP) return false;

        const GenParticle* glep1_0 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1;
        const GenParticle* glep2_0 = diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1;

		if (glep1_0->pdgId()<0 == glep2_0->pdgId()<0) {
			cout<<"same sign leps???"<<endl;
			return false;
		}

		const GenParticle *glep1=0;
		const GenParticle *glep2=0;

		if (glep1_0->absPdgId() == 15) glep1 = getGenLepFromTau(glep1_0);
		else glep1 = glep1_0;

		if (glep2_0->absPdgId() == 15) glep2 = getGenLepFromTau(glep2_0);
		else glep2 = glep2_0;

		if (!glep1 || !glep2) return false;

        plot(prefix,glep1,glep2);
        if (passFullRecoSel()) plot(prefix+"_fullSel",glep1,glep2);

        if (glep1_0->absPdgId()==15 || glep2_0->absPdgId()==15) return false;
        plot(prefix+"_noTaus",glep1,glep2);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getGenVars(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
