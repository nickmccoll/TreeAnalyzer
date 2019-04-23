
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
#include "Processors/Variables/interface/FatJetSelection.h"

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

    void addPlots(TString id, const Electron* el) {
    	int idx = el->index();
//    	float rho = reader_event->rho.val();
    	float scE = (reader_electron->scE)[idx];

		MomentumF nu = HiggsSolver::getInvisible(reader_event->met, (el->p4() + wjj.p4()) );
		MomentumF hh_mom = el->p4() + nu.p4() + wjj.p4() + hbb->p4();
		float hhmass =hh_mom.mass();

        // electron vars
//    	float hoe                   = (reader_electron->HoE)[idx];
//    	float hoe_bc                = (reader_electron->HoE_BC)[idx];
//    	float abs_dEtaSeed          = (reader_electron->abs_dEtaSeed)[idx];
//    	float abs_dPhiIn            = (reader_electron->abs_dPhiIn)[idx];
//    	float full5x5_sigmaIetaIeta = (reader_electron->full5x5_sigmaIetaIeta)[idx];
//    	float abs_1oEm1op           = (reader_electron->abs_1oEm1op)[idx];
//    	float e2x5OverE5x5          = (reader_electron->e2x5OverE5x5)[idx];
//    	float e1x5OverE5x5          = (reader_electron->e1x5OverE5x5)[idx];
//    	size8 missInnerHits         = (reader_electron->missInnerHits)[idx];
//    	size8 passConvVeto          = (reader_electron->passConvVeto)[idx];
//
//    	plotter.getOrMake1DPre(id,"hoe",";hoe",500,0,0.2)->Fill(hoe,weight);
//    	plotter.getOrMake1DPre(id,"hoe_bc",";hoe_bc",500,0,0.2)->Fill(hoe_bc,weight);
//    	plotter.getOrMake1DPre(id,"abs_dEtaSeed",";abs_dEtaSeed",500,0,0.3)->Fill(abs_dPhiIn,weight);
//    	plotter.getOrMake1DPre(id,"abs_dPhiIn",";abs_dPhiIn",500,0,0.3)->Fill(abs_dPhiIn,weight);
//    	plotter.getOrMake1DPre(id,"full5x5_sigmaIetaIeta",";full5x5_sigmaIetaIeta",500,0,0.1)->Fill(full5x5_sigmaIetaIeta,weight);
//    	plotter.getOrMake1DPre(id,"abs_1oEm1op",";abs_1oEm1op",500,0,0.2)->Fill(abs_1oEm1op,weight);
//    	plotter.getOrMake1DPre(id,"e2x5OverE5x5",";e2x5OverE5x5",500,0,1)->Fill(e2x5OverE5x5,weight);
//    	plotter.getOrMake1DPre(id,"e1x5OverE5x5",";e1x5OverE5x5",500,0,1)->Fill(e1x5OverE5x5,weight);
//
//    	plotter.getOrMake1DPre(id,"missInnerHits",";missInnerHits",10,-0.5,9.5)->Fill(missInnerHits,weight);
//    	plotter.getOrMake1DPre(id,"passConvVeto",";passConvVeto",10,-0.5,1.5)->Fill(passConvVeto,weight);

        plotter.getOrMake1DPre(id,"pt",";p_{T}",500,0,500)->Fill(el->pt(),weight);
        plotter.getOrMake1DPre(id,"scE",";E_{SC}",500,0,500)->Fill(scE,weight);
        plotter.getOrMake1DPre(id,"ht",";H_{T}",2000,0,2000)->Fill(ht_chs,weight);
        plotter.getOrMake1DPre(id,"hhmass",";M_{HH}",110,0,5500)->Fill(hhmass,weight);
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

    bool getFatJets(const Electron* el) {
    	wjj.setP4(0.1,0.1,0.1,0.1);
        if(reader_fatjet && reader_fatjet_noLep && el){
        	std::unique_ptr<FatJetProcessor> Proc;
            Proc.reset(new FatJetProcessor ());
            auto fjs = FatJetSelHelpers::selectFatJets(parameters.fatJets,*reader_fatjet);
            auto fjs_noLep = FatJetSelHelpers::selectFatJets(parameters.fatJets,*reader_fatjet_noLep);

            double minDR = 100000;
            int fjIDX = PhysicsUtilities::findNearestDRDeref(*el,fjs_noLep,minDR,parameters.fatJets.wjj_minPT);
            if(fjIDX >= 0 ) {
            	wjj = fjs_noLep[fjIDX]->p4();
            	if (!el->passMVA90ID_noIso() && minDR < 0.8) wjj.p4() -= el->p4();
                if(wjj.pt() < parameters.fatJets.wjj_minPT) wjj.setP4(0.1,0.1,0.1,0.1);
                minDR = PhysicsUtilities::deltaR(wjj,*el);
                if(parameters.fatJets.wjj_maxLepDR > 0 && minDR > parameters.fatJets.wjj_maxLepDR) wjj.setP4(0.1,0.1,0.1,0.1);
            }
            const FatJet wjjFJ(wjj.p4(),0);
            hbb = FatJetSelHelpers::getHbbCand(parameters.fatJets,&wjjFJ,el,fjs);

        	if (hbb) hbbCat = BTagging::getCSVSJCat(parameters.jets,hbb->subJets());
        } else {
        	hbb = 0;
        }
    	if (hbb && wjj.pt() > 10) return true;
    	else return false;
    }

    void prepPlots(TString id, const Electron* el) {

    	if (!getFatJets(el)) return; // inside this function, the hbb and wjj objects are selected

    	float massHbb = hbb->sdMom().mass();
    	addPlots(id,el);
    	if (massHbb > 30 && massHbb < 210) {
    		addPlots(id+"_Mbb30to210",el);
    		if (hbbCat >= BTagging::CSVSJ_MF) addPlots(id+"_Mbb30to210_BTag",el);
    	}
    }

    void testElWorkingPointsBKG(TString prefix) {

    	vector<double> etaWPs = {2.5,2.1,1.5};
    	vector<TString> etaStrs = {"maxEta2p5","maxEta2p1","maxEta1p5"};

    	vector<double> isoWPs = {0.1, 0.2};
    	vector<LeptonProcessor::elFunFloat> isos = {&Electron::miniIso, &Electron::miniIsoFP};
    	vector<TString> isoStrs = {"miniIso","miniIsoFP"};

    	vector<LeptonProcessor::elFunBool> ids = {&Electron::passInclID, &Electron::passMedID_noIso, &Electron::passTightID_noIso, &Electron::passMVA90ID_noIso};
    	vector<TString> idStrs = {"InclID","MedID","TightID","MVAID"};

    	for (unsigned int ieta=0; ieta<etaWPs.size(); ieta++) {

    		LeptonParameters testparam = parameters.leptons;
    		testparam.el_maxETA = etaWPs[ieta];
    		testparam.el_getISO = &Electron::inclIso;
    		testparam.el_getID  = &Electron::passInclID;

    		const auto leps = LeptonProcessor::getLeptons(testparam,*reader_muon,*reader_electron);
    		const auto lep = leps.size() ? leps.front() : 0;

    		if (lep && lep->isElectron()) prepPlots(prefix+"_passIP_"+etaStrs[ieta],(Electron*)lep);

    		for(unsigned int iD=0; iD<ids.size(); iD++) for(unsigned int iS=0; iS<isos.size(); iS++) for(unsigned int iWP=0; iWP<isoWPs.size(); iWP++) {
        		testparam.el_maxISO = isoWPs[iWP];
        		testparam.el_getISO = isos[iS];
        		testparam.el_getID  = ids[iD];

        		const auto leps_ = LeptonProcessor::getLeptons(testparam,*reader_muon,*reader_electron);
        		const auto lep_ = leps_.size() ? leps_.front() : 0;
        		if (lep_ == 0 || lep_->isMuon()) continue;

        		TString isoVal = TString::Format("%.1f",isoWPs[iWP]); isoVal.ReplaceAll(".","p");

        		prepPlots(prefix+"_passIP_"+etaStrs[ieta]+"_"+isoStrs[iS]+isoVal+"_"+idStrs[iD], (Electron*)lep_);

    		}
    	}
    }

    void testElWorkingPointsSIG(TString prefix, const Electron* el) {

        auto passIP = [&](const Lepton* lep)->bool {
        	if (fabs(lep->d0()) > 0.05) return false;
        	if (fabs(lep->dz()) > 0.1 ) return false;
        	if (fabs(lep->sip3D()) > 4.0) return false;
        	return true;
        };

        if (!passIP(el)) return;
        prefix += "_passIP_";

    	vector<double> etaWPs = {2.5,2.1,1.5};
    	vector<TString> etaStrs = {"maxEta2p5","maxEta2p1","maxEta1p5"};

    	vector<double> isoWPs = {0.1, 0.2};
    	vector<TString> isoStrs = {"miniIso","miniIsoFP"};

    	vector<TString> idStrs = {"InclID","MedID","TightID","MVAID"};

    	for (unsigned int ieta=0; ieta<etaWPs.size(); ieta++) {

    		if (el->absEta() > etaWPs[ieta]) continue;
    		prepPlots(prefix+etaStrs[ieta],el);

    		for(unsigned int iD=0; iD<idStrs.size(); iD++) for(unsigned int iS=0; iS<isoStrs.size(); iS++) for(unsigned int iWP=0; iWP<isoWPs.size(); iWP++) {

    			bool pass = true;
    			double maxISO = isoWPs[iWP];

    			if(iS==0 && el->miniIso() > maxISO)   pass = false;
    			if(iS==1 && el->miniIsoFP() > maxISO) pass = false;

    			if(iD==0 && !el->passInclID())        pass = false;
    			if(iD==1 && !el->passMedID_noIso())   pass = false;
    			if(iD==2 && !el->passTightID_noIso()) pass = false;
    			if(iD==3 && !el->passMVA90ID_noIso()) pass = false;

    			if(!pass) continue;

        		TString isoVal = TString::Format("%.1f",isoWPs[iWP]); isoVal.ReplaceAll(".","p");
        		prepPlots(prefix+etaStrs[ieta]+"_"+isoStrs[iS]+isoVal+"_"+idStrs[iD], el);
    		}
    	}
    }

    bool runEvent() override {
        islurm++;
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht_chs < 400) return false;
        TString prefix = smpName;

    	if (!reader_muon || !reader_electron) return false;

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
            const auto *recoL = getMatchedLepton(*genlep,muons,electrons);

            if (!recoL) return false;
            if (recoL->pt() < 30 || recoL->absEta() > 2.5) return false;
            if (!recoL->isElectron()) {cout << "matched lepton is not electron" << endl; return false;}

            plotter.getOrMake1DPre(prefix,"slurm",";M_{HH}",500,0,5000)->Fill(recoL->pt(),weight);
            if(recoL->passInclID()) plotter.getOrMake1DPre(prefix+"_InclID","slurm",";M_{HH}",500,0,5000)->Fill(recoL->pt(),weight);
            if(((Electron*)recoL)->passMedID_noIso())  plotter.getOrMake1DPre(prefix+"_MedID","slurm",";M_{HH}",500,0,5000)->Fill(recoL->pt(),weight);
            if(((Electron*)recoL)->passTightID_noIso())   plotter.getOrMake1DPre(prefix+"_TightID","slurm",";M_{HH}",500,0,5000)->Fill(recoL->pt(),weight);
            if(((Electron*)recoL)->passMVA90ID_noIso())  plotter.getOrMake1DPre(prefix+"_MVAID","slurm",";M_{HH}",500,0,5000)->Fill(recoL->pt(),weight);

            TH1* hi = plotter.getOrMake1DPre(smpName+"_passIP_maxEta2p5_miniIso0p1_InclID","pt",";p_{T}",500,0,500);
            TH1* hm = plotter.getOrMake1DPre(smpName+"_passIP_maxEta2p5_miniIso0p1_MVAID","pt",";p_{T}",500,0,500);


            if (hi->GetEntries() < hm->GetEntries() && itsHappened == false) {
            	cout<<islurm<<endl;
            	itsHappened = true;
            }

            testElWorkingPointsSIG(prefix,(Electron*)recoL);
//            if (recoL->miniIso() < 0.1) testNm1Plots(prefix+"andISO", (Electron*)recoL);

        }
        if(!isSignal()){

        	testElWorkingPointsBKG(prefix);
//        	if (lep->miniIso() < 0.1) testNm1Plots("bkg_passIPandISO", (Electron*)lep);
        }

        return true;
    }

    bool itsHappened = false;
    int islurm = 0;
    const FatJet *hbb = 0;
    MomentumF wjj;
    BTagging::CSVSJ_CAT        hbbCat   = BTagging::CSVSJ_INCL;
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
