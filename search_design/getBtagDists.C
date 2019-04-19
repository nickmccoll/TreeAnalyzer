
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

using namespace TAna;
using namespace FillerConstants;
using namespace std;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
        turnOffCorr(CORR_TRIG);
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
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,20,2.4);
            const auto *recoL = getMatchedLepton(*diHiggsEvt.w1_d1,muons,electrons);

            if (!recoL) return 0;
            bool passKIN = recoL->isElectron() ? (recoL->pt() > 30 && recoL->absEta() < 1.5) : (recoL->pt() > 27 && recoL->absEta() < 2.4);
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

            lep = recoL;
            return lep;

    	} else {
    		LeptonParameters paramLep = parameters.leptons;
    		paramLep.el_maxETA = 1.5;
    		const auto leps = LeptonProcessor::getLeptons(paramLep,*reader_muon,*reader_electron);
    		const auto lep = leps.size() ? leps.front() : 0;
    		if (!lep) return 0;

        	std::unique_ptr<FatJetProcessor> Proc;
            Proc.reset(new FatJetProcessor ());
    		Proc->loadFatJets(parameters.fatJets,*reader_fatjet,*reader_fatjet_noLep,lep);
    		hbb = Proc->getHBBCand();
    		wjj = Proc->getWjjCand();

    		if (!hbb || !wjj) return 0;
    		else return lep;
    	}
    }

    void testJets(TString prefix, TString lepid) {
    	vector<float> thresh_pt = {20,25,30};
    	vector<TString> btagWPs = {"bLoose","bMed","bTight"};

    	for (float pt : thresh_pt) {
    		TString ptid = TString::Format("jetPt%.0f_",pt);
    		auto jets_ = PhysicsUtilities::selObjsMom(reader_jet->jets,pt,2.4);

    		for (const auto& j : jets_) {
                if(!j->passTightID()) continue;
                auto flvI = BTagging::jetFlavor(*j);
                TString flvS = "Fl";
                if(flvI == BTagging::FLV_B) flvS = "Fb";
                else if(flvI == BTagging::FLV_C) flvS = "Fc";

                const float ptj = j->pt();
                const float absETAj = j->absEta();

                auto fill = [&](const TString& label) {
                    plotter.getOrMake2DPre(prefix, flvS+"_"+label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,12,0,2.4)->Fill(ptj,absETAj,weight);
                    plotter.getOrMake2DPre(prefix+"_"+lepid, flvS+"_"+label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,12,0,2.4)->Fill(ptj,absETAj,weight);
                };

                fill(ptid+"bIncl_pt_eta");
            	for (unsigned int ic = 0; ic < btagWPs.size(); ic++) {
            		if (j->csv() >= parameters.jets.CSV_WP[ic+1] && j->csv()) fill(ptid+btagWPs[ic]+"_CSV_pt_eta");
            	}
            	for (unsigned int id = 0; id < btagWPs.size(); id++) {
            		if (j->deep_csv() >= parameters.jets.DeepCSV_WP[id+1]) fill(ptid+btagWPs[id]+"_Deep_pt_eta");
            	}
    		}
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

        auto fill = [&](const TString pre, const TString& label, float pt, float absEta) {
            plotter.getOrMake2DPre(pre, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,12,0,2.4)->Fill(pt,absEta,weight);
            plotter.getOrMake2DPre(pre+"_"+lepid, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,12,0,2.4)->Fill(pt,absEta,weight);
        };

        for (const auto& sj : hbb->subJets()) {
        	auto flvI = BTagging::jetFlavor(sj);

        	TString flvS = "Fl_";
            if(flvI == BTagging::FLV_B) flvS = "Fb_";
            else if(flvI == BTagging::FLV_C) flvS = "Fc_";

            const float pt = sj.pt();
            const float absEta = sj.absEta();

            if (sj.csv() < parameters.jets.CSV_WP[1])  nN_csv++;
            if (sj.csv() >= parameters.jets.CSV_WP[1] && sj.csv() < parameters.jets.CSV_WP[2]) nL_csv++;
            if (sj.csv() >= parameters.jets.CSV_WP[2] && sj.csv() < parameters.jets.CSV_WP[3]) nM_csv++;
            if (sj.csv() >= parameters.jets.CSV_WP[3]) nT_csv++;
            if (sj.deep_csv() < parameters.jets.DeepCSV_WP[1])  nN_deep++;
            if (sj.deep_csv() >= parameters.jets.DeepCSV_WP[1] && sj.deep_csv() < parameters.jets.DeepCSV_WP[2]) nL_deep++;
            if (sj.deep_csv() >= parameters.jets.DeepCSV_WP[2] && sj.deep_csv() < parameters.jets.DeepCSV_WP[3]) nM_deep++;
            if (sj.deep_csv() >= parameters.jets.DeepCSV_WP[3]) nT_deep++;

            fill(prefix,flvS+"bIncl_pt_eta",pt,absEta);
        	for (unsigned int ic = 0; ic < btagWPs.size(); ic++) {
        		if (sj.csv() >= parameters.jets.CSV_WP[ic+1]) fill(prefix,flvS+btagWPs[ic]+"_CSV_pt_eta",pt,absEta);
        	}
        	for (unsigned int id = 0; id < btagWPs.size(); id++) {
        		if (sj.deep_csv() >= parameters.jets.DeepCSV_WP[id+1]) fill(prefix,flvS+btagWPs[id]+"_Deep_pt_eta",pt,absEta);
        	}
        }

        // NN = 0; LN = 1; LL = 2; MN = 3; ML = 4; MM = 5, TN = 6, TL = 7, TM = 8, TT = 9
        auto fillCat = [&](const TString btype, int val) {
            plotter.getOrMake1DPre(prefix,btype+"_bCats","",10,-0.5,9.5)->Fill(val,weight);
            plotter.getOrMake1DPre(prefix+"_"+lepid,btype+"_bCats","",10,-0.5,9.5)->Fill(val,weight);
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
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht_chs < 400) return false;

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
        if(hbb->sdMom().mass() < 10 || wjj->sdMom().mass() < 10) return false;

        TString pre = smpName;
        TString lepid = lep->isMuon() ? "mu" : "el";
        testJets(pre,lepid);
        testHbbBtagging(pre,lepid,hbb);

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
