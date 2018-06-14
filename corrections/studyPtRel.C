
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
    const Jet * getMatchedJet(const MomentumF* mom, const std::vector<const Jet*>& jets, double& nearestDR ){
        int idx = PhysicsUtilities::findNearestDRDeref(mom->p4(),jets,nearestDR);
        if(idx < 0) return 0;
        else return jets[idx];
    }

    void makePlots(const TString& prefix, const MomentumF* leptonMom, const MomentumF* nearestJet=0){
        //PT rel
        	if (nearestJet) {
			TVector3 jml(nearestJet->px(),nearestJet->py(),nearestJet->pz());
			TVector3 lep(leptonMom->px(),leptonMom->py(),leptonMom->pz());
			auto ptrel = (lep.Cross(jml)).Mag() / jml.Mag();
			plotter.getOrMake1DPre(prefix,"ptrel",";P_{T,rel}",100,0,100)->Fill(ptrel, weight);
        	} else {
    			plotter.getOrMake1DPre(prefix,"ptrel",";P_{T,rel}",100,0,100)->Fill(999., weight);
        	}
    }

    void plotSIGwithsel(const TString& prefix, const MomentumF* leptonMom,MomentumF* nearestJet=0, bool hasMedBtag=false) {
    		// Plotting scenarios for BTag and HT cuts
    	    	TString pre;
    	    	if (ht_chs >= 400) {
	    		pre = prefix+"InclB_htgt400_";
	    		makePlots(pre,leptonMom,nearestJet);
    	    		if(hasMedBtag) {
    	    			pre = prefix+"LMTB_htgt400_";
    	    			makePlots(pre,leptonMom,nearestJet);
    	    		}
    	    	}
    	    	if (ht_chs>= 1200) {
	    		pre = prefix+"InclB_htgt1200_";
	    		makePlots(pre,leptonMom,nearestJet);
    	    		if (hasMedBtag) {
    	    			pre = prefix+"LMTB_htgt1200_";
    	    			makePlots(pre,leptonMom,nearestJet);
    	    		}
    	    	}
    }

    void plotBKGwithsel(const TString& prefix, float maxptrel, bool hasMedBtag) {
    		TString pre;
    	    	if (ht_chs >= 400){
    	    		pre = prefix + "InclB_htgt400_";
    	    		plotter.getOrMake1DPre(pre,"maxptrel",";Max P_{T,rel}",100,0,100)->Fill(maxptrel,weight);
    	    		if (hasMedBtag) {
    	    			pre = prefix + "LMTB_htgt400_";
    	    			plotter.getOrMake1DPre(pre,"maxptrel",";Max P_{T,rel}",100,0,100)->Fill(maxptrel,weight);
    	    		}
    	    	}
    	    	if (ht_chs >= 1200) {
    	    		pre = prefix + "InclB_htgt1200_";
    	    		plotter.getOrMake1DPre(pre,"maxptrel",";Max P_{T,rel}",100,0,100)->Fill(maxptrel,weight);
    	    		if (hasMedBtag) {
    	    			pre = prefix + "LMTB_htgt1200_";
    	    			plotter.getOrMake1DPre(pre,"maxptrel",";Max P_{T,rel}",100,0,100)->Fill(maxptrel,weight);
    	    		}
    	    	}
    }

    void doLepton(const TString& prefix, const GenParticle* genP,const std::vector<const Muon*>& muons,const std::vector<const Electron*>& electrons,
            const std::vector<const Jet*>& jets, bool hasMedBtag){

        const auto* recoL = getMatchedLepton(*genP,muons,electrons);

        TString pre;
        //reconstruction eff
        if(recoL == 0) return;
		const bool passID = leptonProcNoISO->isGoodLepton(*reader_event,recoL);
		if(!passID) return; // filter for events that pass ID selection

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
        std::vector<float> ptjcuts = {20,50,150};

        for (const auto& cut : ptjcuts) {

        		if(cut == 20) pre = prefix + "ptjgt20_";
        		else if (cut == 50) pre = prefix + "ptjgt50_";
        		else if (cut == 150) pre = prefix + "ptjgt150_";

        		// vector of jets that pass Pt cut
        		std::vector<MomentumF> slimjets;
        		for (auto jet : LSjets) {
        		    	if (jet.pt() < cut) continue;
        		    	slimjets.push_back(jet);
        		}
        		double nearDRToReco = 20;
			int idx = PhysicsUtilities::findNearestDR(*recoL, slimjets,nearDRToReco,0.4);

			plotSIGwithsel(pre +"ID_",recoL, idx >= 0 ? &slimjets[idx] : 0, hasMedBtag);

			//ISO eff
			const bool passISO = leptonProc->isGoodLepton(*reader_event,recoL);
			if(!passISO) { // filter for leptons that pass ID & miniISO
				plotSIGwithsel(pre+"ID_failISO_", recoL, idx >= 0 ? &slimjets[idx] : 0, hasMedBtag);
				continue;
			}
			plotSIGwithsel(pre +"IDISO_",recoL,idx >= 0 ? &slimjets[idx] : 0, hasMedBtag);
        }
    }

    bool goodGenLepton(const GenParticle* genP) {
        if(genP->absEta() > 2.4) return false;
        if(genP->pt() < ( genP->absPdgId() == ParticleInfo::p_eminus ? 33.0 : 29  )  ) return false;
        return true;
    }

    float getPTrel(const TVector3 lep, const TVector3 jet) {
    		return (lep.Cross(jet)).Mag() / jet.Mag();
    }

    float getMaxPTrel(const std::vector<MomentumF>& leps, const std::vector<MomentumF>& jets) {
    		float maxptrel = 0;
    		for (const auto& lep : leps) {
    			double closestDR = 50.;
    			int iJ = PhysicsUtilities::findNearestDR(lep,jets,closestDR,0.4);
    			if(iJ < 0) return 999;

    			TVector3 l(lep.px(), lep.py(), lep.pz());
			TVector3 j(jets[iJ].px(), jets[iJ].py(), jets[iJ].pz());
			float ptrel = getPTrel(l, j);
			if (ptrel > maxptrel) maxptrel = ptrel;
    		}
    		return maxptrel;
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
            if(diHiggsEvt.type == DiHiggsEvent::TAU_MU) sN = smpName + "_taumu_";
            else if(diHiggsEvt.type == DiHiggsEvent::TAU_E) sN = smpName + "_taue_";
            else if(diHiggsEvt.type == DiHiggsEvent::MU) sN = smpName + "_mu_";
            else sN = smpName + "_e_";

            doLepton(sN, diHiggsEvt.w1_d1, muons,electrons,jets,hasMedBtag);

        } else if(!isRealData()){
			TString sN = smpName;
			std::vector<float> ptcuts = {20,50,150};

			if (muons.size() != 0) {

				// correct the jets with lepton subtraction if necessary: for each lepton in the evt, if any jet is within dR = 0.4 of it, subtract
				// the 4-momentum of the lepton from the jet
				doMuonBKG(sN+"_mu_", muons, electrons, jets, hasMedBtag, ptcuts);
			}
			if (electrons.size() != 0) {
				doElectronBKG(sN+"_e_", muons, electrons, jets, hasMedBtag, ptcuts);
			}
        }
        return true;
    }
    void doMuonBKG(const TString& prefix, const std::vector<const Muon*>& muons, const std::vector<const Electron*>& electrons, const std::vector<const Jet*>& Jets,
    		bool hasMedBtag, std::vector<float> ptcuts) {

    	// get set of muons that pass ID
    			std::vector<const Muon*> muonsID;
    			std::vector<MomentumF> muIDmom;
    			for (const auto& m : muons) {
    				const bool passID = leptonProcNoISO->isGoodLepton(*reader_event,m);
    			    	if(passID) {
    			    		muonsID.push_back(m);
    			    		MomentumF mu(m->p4());
    			    		muIDmom.push_back(mu);
    			    	}
    			}
    			if(muonsID.size() == 0) return;

    			// get a bool that tells whether a muon also passes ISO or if no electron passes ISO
    			bool passISO = false;
    			for (const auto& mu : muonsID) {
    				const bool pass = leptonProc->isGoodLepton(*reader_event,mu);
    				if(pass) {
    					passISO = true;
    					break;
    				}
    			}
    			// get vector of appropriately lepton-subtracted jets using the ID-passing electrons
    			std::vector<MomentumF> LSjets;
    			for (const auto* jet : jets) {
    				double closestDR = 10;
    				int iM = PhysicsUtilities::findNearestDRDeref(*jet,muonsID,closestDR,0.4);
    				if(iM < 0){
    					MomentumF jetmom(jet->p4());
    					LSjets.push_back(jetmom);
    				} else {
    					MomentumF jml(jet->p4() - muonsID[iM]->p4());
    					LSjets.push_back(jml);
    				}
    			}

    	        	// Calculate max PtRel for each of a set of Jet Pt cuts we want to look at
    			TString pre;
    	        	for (auto cut : ptcuts) {
    	        		pre = prefix + TString::Format("ptjgt%.0f_",cut);

    	            	// vector of jets that pass Pt cut
    	            	std::vector<MomentumF> slimjets;
    	        		for (auto jet : LSjets) {
    	        			if (jet.pt() < cut) continue;
    	        			slimjets.push_back(jet);
    	        		}

    	        	   float maxptrelID = getMaxPTrel(muIDmom, slimjets);
    	        	   plotBKGwithsel(pre+"ID_", maxptrelID, hasMedBtag);

    	        	   if(passISO) plotBKGwithsel(pre+"IDISO_", maxptrelID, hasMedBtag);
    	        	   else        plotBKGwithsel(pre+"IDfailISO_", maxptrelID, hasMedBtag);
    	        	}
    }

    void doElectronBKG(const TString& prefix, const std::vector<const Muon*>& muons, const std::vector<const Electron*>& electrons, const std::vector<const Jet*>& jets,
        	bool hasMedBtag, std::vector<float> ptcuts) {


		// get set of electrons that pass ID
		std::vector<const Electron*> electronsID;
		std::vector<MomentumF> eIDmom;
		for (const auto& e : electrons) {
			const bool passID = leptonProcNoISO->isGoodLepton(*reader_event,e);
		    	if(passID) {
		    		electronsID.push_back(e);
		    		MomentumF el(e->p4());
		    		eIDmom.push_back(el);
		    	}
		}
		if(electronsID.size() == 0) return;

		// get a bool that tells whether an electron also passes ISO or if no electron passes ISO
		bool passISO = false;
		for (const auto& e : electronsID) {
			const bool pass = leptonProc->isGoodLepton(*reader_event,e);
			if(pass) {
				passISO = true;
				break;
			}
		}
		// get vector of appropriately lepton-subtracted jets using the ID-passing electrons
		std::vector<MomentumF> LSjets;
		for (const auto* jet : jets) {
			double closestDR = 10;
			int iE = PhysicsUtilities::findNearestDRDeref(*jet,electronsID,closestDR,0.4);
			if(iE < 0){
				MomentumF jetmom(jet->p4());
				LSjets.push_back(jetmom);
			} else {
				MomentumF jml(jet->p4() - electronsID[iE]->p4());
				LSjets.push_back(jml);
			}
		}

        	// Calculate max PtRel for each of a set of Jet Pt cuts we want to look at
		TString pre;
        	for (auto cut : ptcuts) {
        		pre = prefix + TString::Format("ptjgt%.0f_",cut);

            	// vector of jets that pass Pt cut
            	std::vector<MomentumF> slimjets;
        		for (auto jet : LSjets) {
        			if (jet.pt() < cut) continue;
        			slimjets.push_back(jet);
        		}

        	   float maxptrelID = getMaxPTrel(eIDmom, slimjets);
        	   plotBKGwithsel(pre+"ID_", maxptrelID, hasMedBtag);

        	   if(passISO) plotBKGwithsel(pre+"IDISO_", maxptrelID, hasMedBtag);
        	   else        plotBKGwithsel(pre+"IDfailISO_", maxptrelID, hasMedBtag);
        	}
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<LeptonProcessor>     leptonProcNoISO ;

};

#endif

void studyPtRel(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void studyPtRel(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000);
    a.write(outFileName);
}
