
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace TAna;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){

        leptonProcNoISO .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonProcNoISO);
        leptonProcNoISO->lepSelParams.el_maxISO =-1;
        leptonProcNoISO->lepSelParams.mu_maxISO =-1;
        leptonProcNoISO->lepSelParams_dataABCDEF.el_maxISO =-1;
        leptonProcNoISO->lepSelParams_dataABCDEF.mu_maxISO =-1;
    }

    void loadVariables() override  {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        }
        reader_jet_chs =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jet_chs );;
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );
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

    void doSignal(const TString& prefix, const GenParticle* genP,const std::vector<const Muon*>& muons,const std::vector<const Electron*>& electrons,
            const std::vector<const Jet*>& jets, bool hasMedBtag) {

    		const auto* recoL = getMatchedLepton(*genP, muons, electrons);
    		// if there is no reconstructed lepton, then skip the evt (maybe not do this??)
    		if (recoL == 0) return;

    		const bool passID = leptonProcNoISO->isGoodLepton(*reader_event,recoL);
    		if(!passID) return;
		const bool passISO = leptonProc->isGoodLepton(*reader_event,recoL);

		float miniIso = recoL->miniIso();
		plotSIGeff(prefix, passISO, hasMedBtag, miniIso);
    }

    void plotSIGeff(const TString& prefix, bool passISO, bool hasMedBtag, float miniIso) {
    		TString pre = prefix;
    		plotter.getOrMake1DPre(pre+"InclB_","den",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    		if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_","den",";M_{X}",100,0,5000)->Fill(signal_mass,weight);

    		if (passISO) {
    			plotter.getOrMake1DPre(pre+"InclB_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    			if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    		}

    		if (passISO || (ht_chs > 1200 && miniIso < 0.5)) {
    			plotter.getOrMake1DPre(pre+"InclB_sel05_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    			if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel05_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    		}

    		if (passISO || (ht_chs>1200 && miniIso < 0.4)) {
    			   plotter.getOrMake1DPre(pre+"InclB_sel04_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    			   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel04_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    		}
    		if (passISO || (ht_chs>1200 && miniIso < 0.3)) {
    			   plotter.getOrMake1DPre(pre+"InclB_sel03_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    			   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel03_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    		}
    		if (passISO || (ht_chs>1200 && miniIso < 0.2)) {
    			   plotter.getOrMake1DPre(pre+"InclB_sel02_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    			   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel02_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    		}

    		if (passISO || ht_chs > 1200) {
    			plotter.getOrMake1DPre(pre+"InclB_sel0_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    			if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel0_","num",";M_{X}",100,0,5000)->Fill(signal_mass,weight);
    		}
    }

    bool goodGenLepton(const GenParticle* genP) {
        if(genP->absEta() > 2.4) return false;
        if(genP->pt() < ( genP->absPdgId() == ParticleInfo::p_eminus ? 33.0 : 29  )  ) return false;
        return true;
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht_chs < 400) return false;
        bool hasMedBtag = false;

        const std::vector<const Muon    *> muons     = PhysicsUtilities::selObjsMom(reader_muon->muons,26,2.4);
        const std::vector<const Electron*> electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,30,2.4);
        const std::vector<const Jet     *> jets      = PhysicsUtilities::selObjsMom(reader_jet_chs->jets,20,10);

        for (const auto& jet : jets) {
			if (BTagging::isMediumCSVTagged(*jet)) hasMedBtag = true;
        }

        if(reader_event->process == FillerConstants::SIGNAL){
            if(diHiggsEvt.type < DiHiggsEvent::TAU_MU) return false;
            if(!goodGenLepton(diHiggsEvt.w1_d1)) return false;

            TString sN;
            if(diHiggsEvt.type == DiHiggsEvent::TAU_MU) sN = smpName + "_mu_";
            else if(diHiggsEvt.type == DiHiggsEvent::TAU_E) sN = smpName + "_e_";
            else if(diHiggsEvt.type == DiHiggsEvent::MU) sN = smpName + "_mu_";
            else sN = smpName + "_e_";

            doSignal(sN, diHiggsEvt.w1_d1, muons,electrons,jets,hasMedBtag);

        } else if(!isRealData()){
			TString sN = smpName;
			if (muons.size() != 0) sN += "_mu_";
			if (electrons.size() != 0) sN += "_e_";
			doBKG(sN, muons, electrons, jets, hasMedBtag);
        }
        return true;
    }

   void doBKG(const TString& prefix, const std::vector<const Muon*>& muons, const std::vector<const Electron*>& electrons, const std::vector<const Jet*>& Jets,
                bool hasMedBtag) {
	
	// get set of muons and electrons that pass ID
	std::vector<const Muon*> muonsID;
	std::vector<const Electron*> electronsID;

	for (const auto& m : muons) {
	    if (leptonProcNoISO->isGoodLepton(*reader_event, m)) muonsID.push_back(m);
	}
	for (const auto& e : electrons) {
	    if (leptonProcNoISO->isGoodLepton(*reader_event, e)) electronsID.push_back(e);
	}
	if (muonsID.size() == 0 && electronsID.size() == 0) return;
	
	// check evt for leptons that pass the standard isolation criteria, and also grab the max miniIso in the evt
	bool passISO = false;
	float maxminiIso = 0;
	for (const auto& m : muonsID) {
	    if (leptonProc->isGoodLepton(*reader_event, m)) passISO = true;
	    if (m->miniIso() > maxminiIso) maxminiIso = m->miniIso();
	}
	for (const auto& e : electronsID) {
	    if (leptonProc->isGoodLepton(*reader_event, e)) passISO = true;
	    if (e->miniIso() > maxminiIso) maxminiIso = e->miniIso();
	}		
	
	plotBKGeff(prefix, passISO, hasMedBtag, maxminiIso);
   }

   void plotBKGeff(const TString& prefix, bool passISO, bool hasMedBtag, float maxminiIso) {
	   TString pre = prefix;

	   plotter.getOrMake1DPre(pre+"InclB_","den",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
	   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_","den",";H_{T}",100,0,3000)->Fill(ht_chs,weight);

	   if (passISO) {
		   plotter.getOrMake1DPre(pre+"InclB_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
		   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
	   }

	   if (passISO || (ht_chs > 1200 && maxminiIso < 0.5)) {
		   plotter.getOrMake1DPre(pre+"InclB_sel05_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
		   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel05_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
	   }

	   if (passISO || (ht_chs>1200 && maxminiIso < 0.4)) {
			   plotter.getOrMake1DPre(pre+"InclB_sel04_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
			   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel04_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
	   }
	   if (passISO || (ht_chs>1200 && maxminiIso < 0.3)) {
			   plotter.getOrMake1DPre(pre+"InclB_sel03_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
			   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel03_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
	   }
	   if (passISO || (ht_chs>1200 && maxminiIso < 0.2)) {
			   plotter.getOrMake1DPre(pre+"InclB_sel02_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
			   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel02_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
	   }

	   if (passISO || ht_chs > 1200) {
		   plotter.getOrMake1DPre(pre+"InclB_sel0_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
		   if (hasMedBtag) plotter.getOrMake1DPre(pre+"LMTB_sel0_","num",";H_{T}",100,0,3000)->Fill(ht_chs,weight);
	   }
   }

   
    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<LeptonProcessor>     leptonProcNoISO ;

};

#endif

void studyRemoveIso(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void studyRemoveIso(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000);
    a.write(outFileName);
}
