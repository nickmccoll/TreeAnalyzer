
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

using namespace TAna;
using namespace std;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

//        turnOffCorr(CORR_TRIG);
        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);
    }

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

    void testNm1Plots(TString id, const Electron *el) {
    	int idx = el->index();
    	float rho = reader_event->rho.val();
    	float scE = (reader_electron->scE)[idx];

    	float hoe                   = (reader_electron->HoE)[idx];
    	float hoe_bc                = (reader_electron->HoE_BC)[idx];
    	float abs_dEtaSeed          = (reader_electron->abs_dEtaSeed)[idx];
    	float abs_dPhiIn            = (reader_electron->abs_dPhiIn)[idx];
    	float full5x5_sigmaIetaIeta = (reader_electron->full5x5_sigmaIetaIeta)[idx];
    	float abs_1oEm1op           = (reader_electron->abs_1oEm1op)[idx];
//    	float e2x5OverE5x5          = (reader_electron->e2x5OverE5x5)[idx];
//    	float e1x5OverE5x5          = (reader_electron->e1x5OverE5x5)[idx];
    	size8 missInnerHits         = (reader_electron->missInnerHits)[idx];
    	size8 passConvVeto          = (reader_electron->passConvVeto)[idx];

    	double vars[5] = {hoe,abs_dEtaSeed,abs_dPhiIn,full5x5_sigmaIetaIeta,abs_1oEm1op};
    	vector<TString> varstrings = {"hoe","dEtaSeed","dPhiIn","sigmaIetaIeta","1oEm1op"};

    	// following 2017 WPs are taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
    	float hoe_WPB = 0.026 + 1.15/scE + 0.0324*rho/scE;
    	float hoe_WPE = 0.0188 + 2.06/scE + 0.183*rho/scE;
    	double WPT_B[5] = {hoe_WPB, 0.00255, 0.022, 0.0104, 0.159};
    	double WPT_E[5] = {hoe_WPE, 0.00501, 0.0236, .0353, 0.0197};

    	auto makeNm1_Tight = [&](TString id, unsigned int relaxN) {
    		bool pass = true;

    		for (unsigned int i=0; i<5; i++) {
    			if (i == relaxN) continue;
    			if (el->absEta() <= 1.479) {
    				if (vars[i] > WPT_B[i]) pass = false;
    			} else {
    				if (vars[i] > WPT_E[i]) pass = false;
    			}
    			if (missInnerHits > 1) pass = false;
    			if (!passConvVeto) pass = false;
    		}
    		if (pass) plotter.getOrMake1DPre(id,"relax_"+varstrings[relaxN]+"_pt",";p_{T}",500,0,500)->Fill(el->pt(),weight);
    	};

    	makeNm1_Tight(id,0);
    	makeNm1_Tight(id,1);
    	makeNm1_Tight(id,2);
    	makeNm1_Tight(id,3);
    	makeNm1_Tight(id,4);
    }

    void plotIdVars(TString id, const Electron* el) {

    	int idx = el->index();
    	float rho = reader_event->rho.val();
    	float scE = (reader_electron->scE)[idx];

    	// electron vars
    	float hoe                   = (reader_electron->HoE)[idx];
    	float hoe_bc                = (reader_electron->HoE_BC)[idx];
    	float abs_dEtaSeed          = (reader_electron->abs_dEtaSeed)[idx];
    	float abs_dPhiIn            = (reader_electron->abs_dPhiIn)[idx];
    	float full5x5_sigmaIetaIeta = (reader_electron->full5x5_sigmaIetaIeta)[idx];
    	float abs_1oEm1op           = (reader_electron->abs_1oEm1op)[idx];
    	float e2x5OverE5x5          = (reader_electron->e2x5OverE5x5)[idx];
    	float e1x5OverE5x5          = (reader_electron->e1x5OverE5x5)[idx];
    	size8 missInnerHits         = (reader_electron->missInnerHits)[idx];
    	size8 passConvVeto          = (reader_electron->passConvVeto)[idx];

    	// electron ids
    	bool passMedium = el->passMedID_noIso();
    	bool passTight = el->passTightID_noIso();
    	bool passHEEP = el->passHEEPID_noIso();
    	bool passMVA90 = el->passMVA90ID_noIso();

    	// working points for electron Tight ID 2017
    	float hoe_WPB = 0.026 + 1.15/scE + 0.0324*rho/scE;
    	float hoe_WPE = 0.0188 + 2.06/scE + 0.183*rho/scE;
    	float sigmaietaieta5x5_WPB = 0.0104;
    	float sigmaietaieta5x5_WPE = 0.0353;

    	plotter.getOrMake1DPre(id,"hoe",";hoe",200,0,1)->Fill(hoe,weight);
    	plotter.getOrMake1DPre(id,"hoe_bc",";hoe_bc",100,0,1)->Fill(hoe_bc,weight);
    	plotter.getOrMake1DPre(id,"abs_dEtaSeed",";abs_dEtaSeed",100,0,0.75)->Fill(abs_dPhiIn,weight);
    	plotter.getOrMake1DPre(id,"abs_dPhiIn",";abs_dPhiIn",100,0,0.75)->Fill(abs_dPhiIn,weight);
    	plotter.getOrMake1DPre(id,"full5x5_sigmaIetaIeta",";full5x5_sigmaIetaIeta",100,0,0.1)->Fill(full5x5_sigmaIetaIeta,weight);
    	plotter.getOrMake1DPre(id,"abs_1oEm1op",";abs_1oEm1op",100,0,1)->Fill(abs_1oEm1op,weight);
    	plotter.getOrMake1DPre(id,"e2x5OverE5x5",";e2x5OverE5x5",100,0,1)->Fill(e2x5OverE5x5,weight);
    	plotter.getOrMake1DPre(id,"e1x5OverE5x5",";e1x5OverE5x5",100,0,1)->Fill(e1x5OverE5x5,weight);

    	plotter.getOrMake1DPre(id,"missInnerHits",";missInnerHits",10,-0.5,9.5)->Fill(missInnerHits,weight);
    	plotter.getOrMake1DPre(id,"passConvVeto",";passConvVeto",10,-0.5,1.5)->Fill(passConvVeto,weight);

		plotter.getOrMake1DPre(id,"pt",";p_{T}",500,0,500)->Fill(el->pt(),weight);
		plotter.getOrMake1DPre(id,"scE",";E_{SC}",500,0,500)->Fill(scE,weight);

    	if(passMedium) plotter.getOrMake1DPre(id+"_MedID","pt",";pt",500,0,500)->Fill(el->pt(),weight);
    	if(passTight)  plotter.getOrMake1DPre(id+"_TightID","pt",";pt",500,0,500)->Fill(el->pt(),weight);
    	if(passHEEP)   plotter.getOrMake1DPre(id+"_HEEPID","pt",";pt",500,0,500)->Fill(el->pt(),weight);
    	if(passMVA90)  plotter.getOrMake1DPre(id+"_MVAID","pt",";pt",500,0,500)->Fill(el->pt(),weight);

    	if (el->absEta() <= 1.479) {
    		if (hoe < hoe_WPB) {
    			plotter.getOrMake1DPre(id+"_passHoE","pt",";p_{T}",500,0,500)->Fill(el->pt(),weight);
    			plotter.getOrMake1DPre(id+"_passHoE","scE",";E_{SC}",500,0,500)->Fill(scE,weight);
    		}
    	} else {
    		if (hoe < hoe_WPE) {
    			plotter.getOrMake1DPre(id+"_passHoE","pt",";pt",500,0,500)->Fill(el->pt(),weight);
    			plotter.getOrMake1DPre(id+"_passHoE","scE",";E_{SC}",500,0,500)->Fill(scE,weight);
    		}
    	}
		plotter.getOrMake2DPre(id,"scE_x_hoe",";E_{SC};H/E",500,0,500,500,0,1)->Fill(scE,hoe,weight);

		// hh mass with this lepton
		MomentumF nu;
		float hhmass = 0;
		if (wjjCand && hbbCand) {
			nu = HiggsSolver::getInvisible(reader_event->met, (el->p4() + wjjCand->p4()) );
			MomentumF hh_mom = el->p4() + nu.p4() + wjjCand->p4() + hbbCand->p4();
			hhmass = hh_mom.mass();
	    	plotter.getOrMake1DPre(id,"hhmass",";M_{HH}",110,0,5500)->Fill(hhmass,weight);
		}
    }

    bool passIP(const Lepton* lep) {
    	if (fabs(lep->d0()) > 0.05) return false;
    	if (fabs(lep->dz()) > 0.1 ) return false;
    	if (fabs(lep->sip3D()) > 4.0) return false;

    	return true;
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

    const Lepton * getMatchedLepton(const GenParticle& genLepton,
            const std::vector<const Muon *> muons, const std::vector<const Electron*> electrons){
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

    bool passFullSel(const Electron* el) {
    	if (!passTriggerPreselection) return false;
    	if (!wjjCand || !hbbCand) return false;
    	if (wwDM > 125) return false;
    	if (nMedBTags_HbbV > 0) return false;

    	MomentumF hww = el->p4() + wjjCand->p4() + HiggsSolver::getInvisible(reader_event->met, (el->p4() + wjjCand->p4()) ).p4();
    	MomentumF hhmom = hww.p4() + hbbCand->p4();
    	if (hww.pt() / hhmom.mass() < 0.3) return false;
    	if (hbbCSVCat < BTagging::CSVSJ_MF) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhmom.mass() < 700) return false;
    	if (wjjCand->tau2otau1() > 0.75) return false;

    	return true;
    }

    void testElWPs(TString prefix, const Electron* el) {
        auto mkEtaPlot = [&](TString prefix, const Electron* el) {
        	plotIdVars(prefix,el);
        	if (el->absEta() < 2.1) plotIdVars(prefix+"_maxEta2p1",el);
        	if (el->absEta() < 1.5) plotIdVars(prefix+"_maxEta1p5",el);
        };

        if (el->miniIso() < 0.1) mkEtaPlot(prefix+"_miniIso0p1",el);
        if (el->miniIso() < 0.2) mkEtaPlot(prefix+"_miniIso0p2",el);
        if (el->miniIsoFP() < 0.1) mkEtaPlot(prefix+"_miniIsoFP0p1",el);
        if (el->miniIsoFP() < 0.2) mkEtaPlot(prefix+"_miniIsoFP0p2",el);

        if (!passFullSel(el)) return;
        if (el->miniIso() < 0.1) mkEtaPlot(prefix+"_FULLSEL_miniIso0p1",el);
        if (el->miniIso() < 0.2) mkEtaPlot(prefix+"_FULLSEL_miniIso0p2",el);
        if (el->miniIsoFP() < 0.1) mkEtaPlot(prefix+"_FULLSEL_miniIsoFP0p1",el);
        if (el->miniIsoFP() < 0.2) mkEtaPlot(prefix+"_FULLSEL_miniIsoFP0p2",el);
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht_chs < 400) return false;
        TString prefix = smpName;

        if(isSignal() && diHiggsEvt.type == DiHiggsEvent::E){
            const GenParticle *genlep, *q1, *q2;
            if (ParticleInfo::isLepton(diHiggsEvt.w1_d1->pdgId())) {
            	genlep = diHiggsEvt.w1_d1;
            	q1  = diHiggsEvt.w2_d1;
            	q2  = diHiggsEvt.w2_d2;
            } else {
            	genlep = diHiggsEvt.w2_d1;
            	q1  = diHiggsEvt.w1_d1;
            	q2  = diHiggsEvt.w1_d2;
            }
            if (!passGEN_requirements(genlep,q1,q2)) return false;
            if (genlep->pt() < 30 || genlep->absEta() > 2.5) return false;

            const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,20,2.4);
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,20,2.4);
            const auto* recoL = getMatchedLepton(*genlep,muons,electrons);

            if (!recoL) return false;
            if (recoL->pt() < 30 || recoL->absEta() > 2.5) return false;
            if (!recoL->isElectron()) {cout << "matched lepton is not electron" << endl; return false;}

            plotIdVars(prefix, (Electron*)recoL);
            if (!passIP(recoL)) return false;
            prefix += "_passIP";
            plotIdVars(prefix, (Electron*)recoL);

            testElWPs(prefix, (Electron*)recoL);
            if (recoL->miniIso() < 0.1) testNm1Plots(prefix+"andISO", (Electron*)recoL);

        }
        if(!isSignal()){
        	LeptonParameters testparam = parameters.leptons;
        	testparam.el_getID = &Electron::passInclID;
        	testparam.el_maxISO = 9999.;

        	vector<const Lepton*> leps;
        	const Lepton* lep;
            if(reader_electron && reader_muon){
                leps = LeptonProcessor::getLeptons(testparam,*reader_muon,*reader_electron);
                lep  = leps.size() ? leps.front() : 0;
            }
            if (!leps.size() || lep==0) return false;
            if (!lep->isElectron()) return false;

            plotIdVars("bkg_passIP", (Electron*)lep);
            testElWPs("bkg_passIP", (Electron*)lep);
        	if (lep->miniIso() < 0.1) testNm1Plots("bkg_passIPandISO", (Electron*)lep);
        }

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getElectronIDvars(std::string fileName, int treeInt, int randSeed, std::string outFileName,
        float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
