
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Variables/interface/HiggsSolver.h"

using namespace TAna;
using namespace FillerConstants;
using namespace std;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
//        turnOffCorr(CORR_TRIG);
//        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
//        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_TOPPT);
        turnOffCorr(CORR_JER);
    }


    void loadVariables() override {
        reader_event       =loadReader<EventReader>   ("event",isRealData());
        reader_fatjet      =loadReader<FatJetReader>  ("ak8PuppiJet",isRealData());
        reader_fatjet_noLep=loadReader<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false,false);
        reader_jet_chs     =loadReader<JetReader>     ("ak4Jet",isRealData());
        reader_jet         =loadReader<JetReader>     ("ak4PuppiJet",isRealData(),false);
        reader_electron    =loadReader<ElectronReader>("electron");
        reader_muon        =loadReader<MuonReader>    ("muon");

        if(!isRealData()){
            reader_genpart =loadReader<GenParticleReader>   ("genParticle");
        }

        checkConfig();
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

    const Lepton* getLepton(bool isSig) {
    	const Lepton *lep = 0;
    	if (isSig) {
            const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,20,2.4);
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,20,2.5);
            const auto *recoL = getMatchedLepton(*diHiggsEvt.w1_d1,muons,electrons);

            if (!recoL) return 0;
            bool passKIN = recoL->isElectron() ? (recoL->pt() > 30 && recoL->absEta() <= 1.479) : (recoL->pt() > 27 && recoL->absEta() < 2.4);
            if (!passKIN) return 0;
            bool passIPandISO = (fabs(recoL->d0()) < 0.05 && fabs(recoL->dz()) < 0.1 && recoL->sip3D() < 4.0 && recoL->miniIso() < 0.2);
            bool passID = recoL->isMuon() ? ((Muon*)recoL)->passMedID() : ((Electron*)recoL)->passMVA90ID_noIso();

            if (!passID || !passIPandISO) return 0;

            MomentumF hbbMom(diHiggsEvt.hbb->p4());
            MomentumF wjjMom((diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4()));

        	hbb = getMatchedFJ(hbbMom,true,reader_fatjet->jets,reader_fatjet_noLep->jets);
            wjj = getMatchedFJ(wjjMom,true,reader_fatjet->jets,reader_fatjet_noLep->jets);

        	if (!hbb || !wjj) return 0;
            if (PhysicsUtilities::deltaR2(*wjj,*recoL) > 1.2*1.2) return 0;
            if (PhysicsUtilities::deltaR2(*wjj,*hbb) < 1.6*1.6) return 0;
            if (PhysicsUtilities::absDeltaPhi(*hbb,*recoL) < 2.0) return 0;

            auto goodSubStructure = [&](const FatJet* j) {
            	if (j->subJets().size() != 2) return false;
            	bool good = true;
            	for (const auto& sj : j->subJets()) {
            		if (sj.absEta() > 2.4 || sj.pt() < 20) good = false;
            	}
            	return good;
            };

            if (!goodSubStructure(hbb) || !goodSubStructure(wjj)) return 0;
            if (hbb->pt() < 200) return 0;

            lep = recoL;
            return lep;

    	} else {
    		if (selectedLeptons.size() != 1) return 0;
    		if (!selectedLepton) return 0;
    		if (!hbbCand || !wjjCand) return 0;

    		hbb = hbbCand;
    		wjj = wjjCand;
    		return selectedLepton;
    	}
    }

    void testJets(TString prefix, TString lepid, const FatJet *hbbJet) {
    	bool pass = false;
    	if (hbbJet->nSubJets() == 2) {
    		float b1 = hbbJet->subJets()[0].deep_csv();
    		float b2 = hbbJet->subJets()[1].deep_csv();
    		float looseWP = parameters.jets.DeepCSV_WP[1];
    		float medWP   = parameters.jets.DeepCSV_WP[2];

    		if ( (b1 >= medWP || b2 >= medWP) || (b1 >= looseWP && b2 >= looseWP) ) pass = true;
    	}
    	if (!pass) return;

    	vector<float> thresh_pt = {20,25,30};
    	vector<TString> btagWPs = {"bLoose","bMed","bTight"};

    	for (float pt : thresh_pt) {
    		TString ptid = TString::Format("jetPt%.0f_",pt);
    		auto jets_ = PhysicsUtilities::selObjsMom(reader_jet->jets,pt,2.4);

    		// test btag WP efficiencies with the jets in the event
    		for (const auto& j : jets_) {
                if(!j->passTightID()) continue;
                auto flvI = BTagging::jetFlavor(*j);
                TString flvS = "Fl";
                if(flvI == BTagging::FLV_B) flvS = "Fb";
                else if(flvI == BTagging::FLV_C) flvS = "Fc";

                const float ptj = j->pt();
                const float absETAj = j->absEta();

                auto fill = [&](const TString& label) {
                    plotter.getOrMake2DPre(prefix, flvS+"_"+label+"_pt_eta",";jet p_{T}[GeV];jet |#eta|",196,20,1000,12,0,2.4)->Fill(ptj,absETAj,weight);
                    plotter.getOrMake2DPre(prefix+"_"+lepid, flvS+"_"+label+"_pt_eta",";jet p_{T}[GeV];jet |#eta|",196,20,1000,12,0,2.4)->Fill(ptj,absETAj,weight);
                    plotter.getOrMake1DPre(prefix, flvS+"_"+label+"_mbb",";m_{bb}",90,30,210)->Fill(hbbJet->sdMom().mass(),weight);
                    plotter.getOrMake1DPre(prefix+"_"+lepid, flvS+"_"+label+"_mbb",";m_{bb}",90,30,210)->Fill(hbbJet->sdMom().mass(),weight);
                };

                fill(ptid+"bIncl");
            	for (unsigned int ic = 0; ic < btagWPs.size(); ic++) {
            		if (j->csv() >= parameters.jets.CSV_WP[ic+1] && j->csv()) fill(ptid+btagWPs[ic]+"_CSV");
            	}
            	for (unsigned int id = 0; id < btagWPs.size(); id++) {
            		if (j->deep_csv() >= parameters.jets.DeepCSV_WP[id+1]) fill(ptid+btagWPs[id]+"_Deep");
            	}
    		}

    		// test btag veto on jets which do not overlap with H->bb fatjet
    		// For each of L,M,T working points of DeepCSV and CSV btagging, get dist of num_jets
    		// in an event, not overlapping hbbJet, with btag values at least as great as the WP
            auto jets_HbbV_ = PhysicsUtilities::selObjsD(jets_,
                    [&](const Jet* j){return  PhysicsUtilities::deltaR2(*j,*hbbJet ) >= 1.2*1.2;});

            auto getN_bWP = [&] (bool doDeep, int idx)->int {
            	int num = 0;
            	for (const auto& j : jets_HbbV_) {
            		if (doDeep) {
            			if (j->deep_csv() >= parameters.jets.DeepCSV_WP[idx+1]) num++;
            		} else {
            			if (j->csv() >= parameters.jets.CSV_WP[idx+1]) num++;
            		}
            	}
            	return num;
            };

            auto fillVeto = [&] (TString pref, bool isDeep) {
            	TString btype = isDeep ? "DeepCSV" : "CSV";

            	for (unsigned int k=0; k<btagWPs.size(); k++) {
            		TString bWP = btagWPs[k];
            		int num = getN_bWP(isDeep,k);
                	plotter.getOrMake1DPre(pref,ptid+bWP+"_"+btype+"_nSepAK4","",6,-0.5,5.5)->Fill(num,weight);
            	}
            };

            fillVeto(prefix,true);
            fillVeto(prefix+"_"+lepid,true);
            fillVeto(prefix,false);
            fillVeto(prefix+"_"+lepid,false);
    	}
    }

    void testHbbBtagging(TString prefix, TString lepid, const FatJet* hbb) {
    	vector<TString> btagWPs = {"bLoose","bMed","bTight"};

    	int nN_csv = 0;
    	int nN_deep = 0;

        int nL_csv = 0;
        int nM_csv = 0;
        int nT_csv = 0;
        int nL_deep = 0;
        int nM_deep = 0;
        int nT_deep = 0;

        int numTrueB = 0;

        auto fill = [&](const TString pre, const TString& label, float pt, float absEta) {
            plotter.getOrMake2DPre(pre, label+"_pt_eta",";jet p_{T}[GeV];jet |#eta|",196,20,1000,12,0,2.4)->Fill(pt,absEta,weight);
            plotter.getOrMake2DPre(pre+"_"+lepid, label+"_pt_eta",";jet p_{T}[GeV];jet |#eta|",196,20,1000,12,0,2.4)->Fill(pt,absEta,weight);
            plotter.getOrMake1DPre(pre, label+"_mbb",";m_{bb}",90,30,210)->Fill(hbb->sdMom().mass(),weight);
            plotter.getOrMake1DPre(pre+"_"+lepid, label+"_mbb",";m_{bb}",90,30,210)->Fill(hbb->sdMom().mass(),weight);
        };

        for (const auto& sj : hbb->subJets()) {
        	auto flvI = BTagging::jetFlavor(sj);

        	TString flvS = "Fl_";
            if(flvI == BTagging::FLV_B) {
            	flvS = "Fb_";
            	numTrueB++;
            }
            else if(flvI == BTagging::FLV_C) flvS = "Fc_";

            const float pt = sj.pt();
            const float absEta = sj.absEta();

            if (sj.csv() >= parameters.jets.CSV_WP[3]) nT_csv++;
            else if (sj.csv() >= parameters.jets.CSV_WP[2]) nM_csv++;
            else if (sj.csv() >= parameters.jets.CSV_WP[1]) nL_csv++;
            else nN_csv++;

            if (sj.deep_csv() >= parameters.jets.DeepCSV_WP[3]) nT_deep++;
            else if (sj.deep_csv() >= parameters.jets.DeepCSV_WP[2]) nM_deep++;
            else if (sj.deep_csv() >= parameters.jets.DeepCSV_WP[1]) nL_deep++;
            else nN_deep++;

            fill(prefix,"sj_"+flvS+"bIncl",pt,absEta);
        	for (unsigned int ic = 0; ic < btagWPs.size(); ic++) {
        		if (sj.csv() >= parameters.jets.CSV_WP[ic+1]) fill(prefix,"sj_"+flvS+btagWPs[ic]+"_CSV",pt,absEta);
        	}
        	for (unsigned int id = 0; id < btagWPs.size(); id++) {
        		if (sj.deep_csv() >= parameters.jets.DeepCSV_WP[id+1]) fill(prefix,"sj_"+flvS+btagWPs[id]+"_Deep",pt,absEta);
        	}

        }

        // NN = 0; LN = 1; LL = 2; MN = 3; ML = 4; MM = 5, TN = 6, TL = 7, TM = 8, TT = 9
        auto fillCat = [&](const TString btype, int val) {
            plotter.getOrMake1DPre(prefix+"_sj",btype+"_bCats","",10,-0.5,9.5)->Fill(val,weight);
            plotter.getOrMake1DPre(prefix+"_sj_"+lepid,btype+"_bCats","",10,-0.5,9.5)->Fill(val,weight);
            plotter.getOrMake2DPre(prefix+"_sj",btype+"_bCat_nTB","",10,-0.5,9.5,3,-0.5,2.5)->Fill(val,numTrueB,weight);
            plotter.getOrMake2DPre(prefix+"_sj_"+lepid,btype+"_bCat_nTB","",10,-0.5,9.5,3,-0.5,2.5)->Fill(val,numTrueB,weight);

            plotter.getOrMake1DPre(prefix+"_sj", btype+"_bin"+TString::Format("%d_mbb",val),";m_{bb}",90,30,210)->Fill(hbb->sdMom().mass(),weight);
            plotter.getOrMake1DPre(prefix+"_sj_"+lepid, btype+"_bin"+TString::Format("%d_mbb",val),";m_{bb}",90,30,210)->Fill(hbb->sdMom().mass(),weight);
        };

        if (nN_csv == 2)                fillCat("csv",0);
        if (nN_csv == 1 && nL_csv == 1) fillCat("csv",1);
        if (nL_csv == 2)                fillCat("csv",2);
        if (nN_csv == 1 && nM_csv == 1) fillCat("csv",3);
        if (nL_csv == 1 && nM_csv == 1) fillCat("csv",4);
        if (nM_csv == 2)                fillCat("csv",5);
        if (nT_csv == 1 && nN_csv == 1) fillCat("csv",6);
        if (nT_csv == 1 && nL_csv == 1) fillCat("csv",7);
        if (nT_csv == 1 && nM_csv == 1) fillCat("csv",8);
        if (nT_csv == 2) fillCat("csv",9);

        if (nN_deep == 2)                 fillCat("deep",0);
        if (nN_deep == 1 && nL_deep == 1) fillCat("deep",1);
        if (nL_deep == 2)                 fillCat("deep",2);
        if (nN_deep == 1 && nM_deep == 1) fillCat("deep",3);
        if (nL_deep == 1 && nM_deep == 1) fillCat("deep",4);
        if (nM_deep == 2)                 fillCat("deep",5);
        if (nT_deep == 1 && nN_deep == 1) fillCat("deep",6);
        if (nT_deep == 1 && nL_deep == 1) fillCat("deep",7);
        if (nT_deep == 1 && nM_deep == 1) fillCat("deep",8);
        if (nT_deep == 2) fillCat("deep",9);

        vector<float> pts = {20,30,80};
        JetParameters jparam = parameters.jets;
        for (const auto& pt : pts) {
        	jparam.minBtagJetPT = pt;
        	int sjCat = BTagging::getCSVSJCat(jparam,hbb->subJets());
        	TString ptid = TString::Format("pt_%.0f",pt);
            plotter.getOrMake1DPre(prefix,ptid+"_CSVCat","",7,-0.5,6.5)->Fill(sjCat,weight);
        }
    }

    double getMassHH(const Lepton* lep, const FatJet *bb, const FatJet *jj) {
        HiggsSolverInfo hwwInfo;
        double chi   = higgsSolver->hSolverMinimization(lep->p4(),jj->p4(),
                reader_event->met.p4(),jj->sdMom().mass() <60,parameters.hww, &hwwInfo);

        MomentumF hww = hwwInfo.hWW;
        return (bb->p4()+hww.p4()).mass();
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;

    	if(!reader_fatjet || !reader_fatjet_noLep) return false;

        if (isSignal()) {
        	if (diHiggsEvt.type < DiHiggsEvent::MU) return false;
        	if (!passGEN_requirements(diHiggsEvt.w1_d1,diHiggsEvt.w2_d1,diHiggsEvt.w2_d2)) return false;

        	bool passGEN = diHiggsEvt.w1_d1->absPdgId() == 11 ? (diHiggsEvt.w1_d1->pt() > 30 && diHiggsEvt.w1_d1->absEta() < 2.5)
        			: (diHiggsEvt.w1_d1->pt() > 27 && diHiggsEvt.w1_d1->absEta() < 2.4);
            if (!passGEN) return false;
        }

    	const Lepton *lep = getLepton(isSignal());
    	if(!lep) return false;
//        if(wjj->sdMom().mass() < 10) return false;
    	float bbmass = isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(hbb) : hbb->sdMom().mass();
        if(bbmass < 30 || bbmass > 210) return false;

        HiggsSolverInfo hwwInfo;
        double chi   = higgsSolver->hSolverMinimization(lep->p4(),wjj->p4(),
                reader_event->met.p4(),wjj->sdMom().mass() <60,parameters.hww, &hwwInfo);
        MomentumF hww = hwwInfo.hWW;
        double hhmass = (hbb->p4()+hww.p4()).mass();

        if (hhmass < 700 || hhmass > 4000) return false;
        TString pre = smpName;
        TString lepid = lep->isMuon() ? "mu" : "el";

        testJets(pre,lepid,hbb);
        testHbbBtagging(pre,lepid,hbb);

        if (hww.pt()/hhmass < 0.3) return false;
        if (wjj->tau2otau1() > 0.75) return false;
        if (hwwInfo.chiSq > 11) return false;

        const auto jets_veto = PhysicsUtilities::selObjsD(jets,
                [&](const Jet* j){return  PhysicsUtilities::deltaR2(*j,*hbb ) >= 1.2*1.2;});
        const auto nMedBTags_veto = PhysicsUtilities::selObjsMomD(jets_veto,
                parameters.jets.minBtagJetPT,parameters.jets.maxBTagJetETA,
                [&](const Jet* j){return BTagging::passJetBTagWP(parameters.jets,*j);} ).size();

        if (nMedBTags_veto > 0) return false;

        testJets(pre+"_SR",lepid,hbb);
        testHbbBtagging(pre+"_SR",lepid,hbb);

        return true;
    }

    const FatJet *hbb = 0;
    const FatJet *wjj = 0;
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getBtagDists(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
