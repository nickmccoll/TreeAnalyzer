
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

    void makeHist(const TString& prefix, const Lepton* lep, bool passISO, bool passPtRel, float ht){

    		int val = getBin(lep, ht);
    		plotter.getOrMake2DPre(prefix+"ID", "quad",";M_{sig}",32,800,4000,4,0.5,4.5)->Fill(signal_mass, val, weight);
    		plotter.getOrMake1DPre(prefix+"ID", "eff", ";M_{sig}",32,800,4000)->Fill(signal_mass, weight);

    		if (passISO) {             plotter.getOrMake2DPre(prefix+"Iso", "quad", ";M_{sig}",32,800,4000,4,0.5,4.5)->Fill(signal_mass, val, weight);
    								plotter.getOrMake1DPre(prefix+"Iso", "eff", ";M_{sig}",32,800,4000)->Fill(signal_mass, weight);
    		}
    		if (passPtRel)      {
    			plotter.getOrMake2DPre(prefix+"PtRel","quad", ";M_{sig}",32,800,4000,4,0.5,4.5)->Fill(signal_mass, val, weight);
    			plotter.getOrMake1DPre(prefix+"PtRel", "eff", ";M_{sig}",32,800,4000)->Fill(signal_mass, weight);
    		}
    		if (passISO && passPtRel) plotter.getOrMake2DPre(prefix+"IsoANDPtRel","quad", ";M_{sig}",32,800,4000,4,0.5,4.5)->Fill(signal_mass, val, weight);
    }

	int getBin(const Lepton* lep, float ht) {
		float pt_thresh = 60;
		if (lep->pt() < pt_thresh) {
			if (ht < 1200) return 1;
			else return 2;
		} else {
			if (ht < 1200) return 3;
			else return 4;
		}
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

    void doLepton(const TString& prefix, const GenParticle* genP,const std::vector<const Muon*>& muons,const std::vector<const Electron*>& electrons,
            const std::vector<const Jet*>& jets, float ht){

        TString pre = prefix;
        const auto* recoL = getMatchedLepton(*genP,muons,electrons);
        if (recoL == 0) return;

        // First check if the lepton passes ID
        	const bool passID = leptonProcNoISO->isGoodLepton(*reader_event,recoL);
        if (!passID) return;

        // get appropriately lepton-subtracted and then, according to Pt cut, slimmed jets
        std::vector<MomentumF> LSjets;
        for (const auto* jet : jets) {
        		if (PhysicsUtilities::deltaR2(*jet,*recoL) < 0.4*0.4) {
        			MomentumF jml(jet->p4() - recoL->p4());
        			LSjets.push_back(jml);
        		} else {
        			MomentumF jetmom(jet->p4());
        			LSjets.push_back(jetmom);
        		}
        }
        float ptjcut = 20;
        std::vector<MomentumF> slimjets;
        for (auto jet : LSjets) {
            if (jet.pt() < ptjcut) continue;
            slimjets.push_back(jet);
        }

        bool passPtRel = false;
        float PtRel_Cut = 20;

        // check if the lepton passes PtRel: either greater than the cut or the nearest jet is > dR=0.4 away
        	double nearDRToReco = 20;
		int idx = PhysicsUtilities::findNearestDR(*recoL, slimjets,nearDRToReco,0.4);
		if (idx >=0) {
			TVector3 lep(recoL->px(), recoL->py(), recoL->pz());
			TVector3 nearestjet(slimjets[idx].px(), slimjets[idx].py(), slimjets[idx].pz());
			float ptrel = getPTrel(lep, nearestjet);
			if (ptrel > PtRel_Cut) passPtRel = true;
		} else passPtRel = true;

		// check if lepton with same flavor as genLepton passes miniISO
		const bool passISO = leptonProc->isGoodLepton(*reader_event,recoL);

		makeHist(pre, recoL, passISO, passPtRel, ht);
    }

    bool goodGenLepton(const GenParticle* genP) {
        if(genP->absEta() > 2.4) return false;
        if(genP->pt() < ( genP->absPdgId() == ParticleInfo::p_eminus ? 33.0 : 29  )  ) return false;
        return true;
    }

    float getPTrel(const TVector3 lep, const TVector3 jet) {
    		return (lep.Cross(jet)).Mag() / jet.Mag();
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht_chs < 400) return false;
        bool hasMedBtag = false;

        const std::vector<const Muon    *> muons     = PhysicsUtilities::selObjsMom(reader_muon->muons,26,2.4);
        const std::vector<const Electron*> electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,30,2.4);
        const std::vector<const Jet     *> jets      = PhysicsUtilities::selObjsMom(reader_jet_chs->jets,20,10);

        if(reader_event->process == FillerConstants::SIGNAL){
            if(diHiggsEvt.type < DiHiggsEvent::TAU_MU) return false;
            if(!goodGenLepton(diHiggsEvt.w1_d1)) return false;

            TString sN;
            if(diHiggsEvt.type == DiHiggsEvent::TAU_MU) sN = "taumu_";
            else if(diHiggsEvt.type == DiHiggsEvent::TAU_E) sN = "taue_";
            else if(diHiggsEvt.type == DiHiggsEvent::MU) sN = "mu_";
            else sN = "e_";

            doLepton(sN, diHiggsEvt.w1_d1, muons,electrons,jets, ht_chs);
        }
        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<LeptonProcessor>     leptonProcNoISO ;

};

#endif

void studyIsoPtRel(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void studyIsoPtRel(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000);
    a.write(outFileName);
}
