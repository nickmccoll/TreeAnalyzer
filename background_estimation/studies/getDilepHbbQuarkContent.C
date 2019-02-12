#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DileptonSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "TreeReaders/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


#include "TSystem.h"
using namespace TAna;
using namespace std;

class Analyzer : public DileptonSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DileptonSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    }
    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event);
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiJet",isRealData());       load(reader_fatjet);
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiJet",isRealData(),false); load(reader_jet);
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon);

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");          load(reader_genpart   );
        }
    }

	TString getDilepChan(const Lepton* lep1, const Lepton* lep2) {
		if (lep1->isMuon() && lep2->isMuon()) return "_mumu_";
		else if (lep1->isElectron() && lep2->isElectron()) return "_ee_";
		else return "_emu_";
	}

	void getQuarkContent() {
		double maxDR2 = 0.8*0.8;

        int topDecayType = 0; // NONE b wj wjb wjj wjjb bb wjbb wjjbb
        int maxQuarksFromTop = 0;
        int totQuarksFromTops = 0;
        int numB = 0;

        for (const auto& d : smDecayEvt.topDecays) {
        	if (d.type == TopDecay::BAD) continue;
        	if (d.type > TopDecay::HAD) {
        		if (PhysicsUtilities::deltaR2(*d.b,*hbbCand) < maxDR2) {
        			totQuarksFromTops++;
        			numB++;
        			if (maxQuarksFromTop == 0) maxQuarksFromTop = 1;
        		}
        	} else {
        		if (!d.b) continue;
        		bool passB = false;
        		if (PhysicsUtilities::deltaR2(*d.b,*hbbCand) < maxDR2) passB = true;
        		int nW = (PhysicsUtilities::deltaR2(*d.W_decay.dau1,*hbbCand) < maxDR2) + (PhysicsUtilities::deltaR2(*d.W_decay.dau2,*hbbCand) < maxDR2);
        		int nT = nW + passB;
        		totQuarksFromTops += nT;
        		numB += passB;
        		if (nT > maxQuarksFromTop) maxQuarksFromTop = nT;
        	}
        }

        if (totQuarksFromTops == 1){
            if(numB == 1) topDecayType = 1;
            else topDecayType = 2;
        } else if(totQuarksFromTops == 2){
            if (numB == 1) topDecayType = 3;
            else if (numB == 2) topDecayType = 6;
            else topDecayType = 4;
        } else if(totQuarksFromTops == 3) {
        	if (numB == 1) topDecayType = 5;
        	else topDecayType = 7;
        } else if (totQuarksFromTops == 4) topDecayType = 8;

        int WDecayType = 0; // NONE b wj wjb wjj wjjb bb wjbb wjjbb
        int maxQuarksFromW = 0;
        int totQuarksFromWs = 0;

        for (const auto& d : smDecayEvt.bosonDecays) {
            if(d.type != BosonDecay::Z_HAD && d.type != BosonDecay::W_HAD ) continue;
            int nW = (PhysicsUtilities::deltaR2(*d.dau1,*hbbCand) < maxDR2) +  (PhysicsUtilities::deltaR2(*d.dau2,*hbbCand) < maxDR2);
            totQuarksFromWs += nW;
            if (nW > maxQuarksFromW) maxQuarksFromW = nW;
        }

        if (totQuarksFromWs == 1) WDecayType = 2;
        else if (totQuarksFromWs == 2) WDecayType = 4;

        int decayType = 0;
//        int nExtraQuarks = 0;

        if(maxQuarksFromTop >= maxQuarksFromW){
            decayType = topDecayType;
//            nExtraQuarks = (totQuarksFromTops - maxQuarksFromTop) + totQuarksFromWs;
        } else {
            decayType = WDecayType;
//            nExtraQuarks = (totQuarksFromWs - maxQuarksFromW) + totQuarksFromTops;
        }

        TString ch = getDilepChan(selectedDileptons[0],selectedDileptons[1]);
        plotter.getOrMake1DPre(smpName+ch,"decayType",";decay type",9,-0.5,8.5)->Fill(decayType,weight);
//        plotter.getOrMake1DPre(smpName+ch,"extraQuarks",";N_{extra}",6,0,6)->Fill(nExtraQuarks,weight);

        plotter.getOrMake1DPre("bkg"+ch,"decayType",";decay type",9,-0.5,8.5)->Fill(decayType,weight);
//        plotter.getOrMake1DPre("bkg"+ch,"extraQuarks",";N_{extra}",6,0,6)->Fill(nExtraQuarks,weight);
	}

	bool passKinCutsLL(const Lepton *lep1, const Lepton *lep2) {
		double mll = (lep1->p4()+lep2->p4()).mass();
		double dPhi_metLL = PhysicsUtilities::deltaPhi(lep1->p4()+lep2->p4(),reader_event->met);
		double dR2_ll = PhysicsUtilities::deltaR2(*lep1,*lep2);

    	if (mll < 12 || mll > 75) return false;
    	if (dR2_ll > 1.6*1.6) return false;
    	if (fabs(dPhi_metLL) > TMath::PiOver2()) return false;

    	return true;
	}

    bool runEvent() override {
        if(!DileptonSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;

        // cuts before separating into SR and CR
        if(!hbbCand) return false;
        if(selectedDileptons.size() != 2) return false;
        if (getDilepChan(selectedDileptons[0],selectedDileptons[1]).Contains("_ee_")) {
        	if (!(((const Electron*)selectedDileptons[0])->passMedID_noISO() && ((const Electron*)selectedDileptons[1])->passMedID_noISO())) return false;
        }
        if (hbbMass < 30 || hbbMass > 210) return false;
        if (ht_puppi < 400) return false;
        if (hh.mass() < 700) return false;
        if (hbbCSVCat < BTagging::CSVSJ_MF) return false;

        if (!passKinCutsLL(selectedDileptons[0],selectedDileptons[1])) return false;
        if (!isRealData()) getQuarkContent();

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void getDilepHbbQuarkContent(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.analyze();
    a.write(outFileName);
}
void getDilepHbbQuarkContent(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
