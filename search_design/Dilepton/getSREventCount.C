#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DileptonSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "Configuration/interface/FillerConstants.h"

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

    void plotVars(TString sn, TString plotID) {
    	plotter.getOrMake1DPre(smpName+sn+plotID,"ht",";H_{T}",100,0,4000)->Fill(ht_puppi,weight);
    	plotter.getOrMake1DPre(smpName+sn+plotID,"Mbb",";M_{bb}",40,30,210)->Fill(hbbMass,weight);
    	plotter.getOrMake1DPre(smpName+sn+plotID,"Mhh",";M_{HH}",100,0,4500)->Fill(hh.mass(),weight);
    	plotter.getOrMake1DPre(smpName+sn+plotID,"Mww",";M_{WW}",100,0,400)->Fill(hww.mass(),weight);

    	plotter.getOrMake1DPre(smpName+sn+plotID,"hbbCSVCat",";hbbCSVCat",7,-0.5,6.5)->Fill(hbbCSVCat,weight);

    	if (!isSignal()) {
        	plotter.getOrMake1DPre("bkg"+sn+plotID,"ht",";H_{T}",100,0,4000)->Fill(ht_puppi,weight);
        	plotter.getOrMake1DPre("bkg"+sn+plotID,"Mbb",";M_{bb}",40,30,210)->Fill(hbbMass,weight);
        	plotter.getOrMake1DPre("bkg"+sn+plotID,"Mhh",";M_{HH}",100,0,4500)->Fill(hh.mass(),weight);
        	plotter.getOrMake1DPre("bkg"+sn+plotID,"Mww",";M_{WW}",100,0,400)->Fill(hww.mass(),weight);

        	plotter.getOrMake1DPre("bkg"+sn+plotID,"hbbCSVCat",";hbbCSVCat",7,-0.5,6.5)->Fill(hbbCSVCat,weight);
    	}

    }
    void plotByMassWindows(TString sn) {
    	plotVars(sn,"_incl");

    	bool inHbbWindow = (hbbMass > 100 && hbbMass < 150);
    	bool inHHWindow  = (hh.mass() > 2000);
    	if (inHbbWindow)               plotVars(sn,"_Mbb100to150");
    	if (inHHWindow)                plotVars(sn,"_MhhGt2000");
    	if (inHbbWindow && inHHWindow) plotVars(sn,"_Mbb100to150_MhhGt2000");
    }

    void makePlots(TString sn, const Lepton* lep1, const Lepton* lep2) {
    	double mll = (lep1->p4()+lep2->p4()).mass();
    	double dPhi_metLL = PhysicsUtilities::deltaPhi(lep1->p4()+lep2->p4(),reader_event->met);
    	double dR2_ll = PhysicsUtilities::deltaR2(*lep1,*lep2);

    	if (mll < 12 || mll > 75) return;
    	if (dR2_ll > 1.6*1.6) return;

    	plotByMassWindows(sn+"_relax_dPhiBvetoBtag");
        if (nMedBTags_HbbV == 0 && hbbCSVCat >= BTagging::CSVSJ_MF)                   plotByMassWindows(sn+"_relax_dPhi");
    	if (fabs(dPhi_metLL) < TMath::PiOver2() && hbbCSVCat >= BTagging::CSVSJ_MF)   plotByMassWindows(sn+"_relax_Bveto");

        if (nMedBTags_HbbV != 0) return;
    	if (fabs(dPhi_metLL) > TMath::PiOver2()) return;

    	plotByMassWindows(sn+"_relax_Btag");

    	for (int btag = 1; btag < 4; btag++) {
    		if (hbbCSVCat < btag) continue;
    		plotByMassWindows(sn+"_BtagCat_"+TString::Format("%i",btag));
    	}

        if (hbbCSVCat < BTagging::CSVSJ_MF) return;
        plotByMassWindows(sn);
    }

    void plotSR(const Lepton* lep1, const Lepton* lep2) {
    	makePlots("_SR_btagLMT",lep1,lep2);

    	TString chan = dilepMap[dilepChan];
    	makePlots(chan+"SR_btagLMT",lep1,lep2);

    	if (hbbCSVCat == BTagging::CSVSJ_MF) {
        	makePlots("_SR_btagL",lep1,lep2);
        	makePlots(chan+"SR_btagL",lep1,lep2);
    	} else if (hbbCSVCat == BTagging::CSVSJ_ML) {
        	makePlots("_SR_btagM",lep1,lep2);
        	makePlots(chan+"SR_btagM",lep1,lep2);
    	} else if (hbbCSVCat == BTagging::CSVSJ_MM) {
        	makePlots("_SR_btagT",lep1,lep2);
        	makePlots(chan+"SR_btagT",lep1,lep2);
    	}
    }

    bool passMMeID(const Lepton* lep1, const Lepton* lep2) {
    	bool pass = false;
    	if (((const Electron*)lep1)->passMedID_noISO() && ((const Electron*)lep2)->passMedID_noISO()) pass = true;
    	return pass;
    }

    bool runEvent() override {
        if(!DileptonSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;

        // cuts before separating into SR and CR
        if(!hbbCand) return false;
        if(selectedDileptons.size() != 2) return false;

        if(isSignal()) plotter.getOrMake1DPre(smpName,"hbbCSVCat",";hbbCSVCat",7,-0.5,6.5)->Fill(hbbCSVCat,weight);
        else plotter.getOrMake1DPre("bkg","hbbCSVCat",";hbbCSVCat",7,-0.5,6.5)->Fill(hbbCSVCat,weight);

        if (dilepChan == ee && !passMMeID(selectedDileptons[0],selectedDileptons[1])) return false;

        if (hbbMass < 30 || hbbMass > 210) return false;
        if (ht_puppi < 400) return false;
        if (hh.mass() < 700) return false;
        if (reader_event->met.pt() < 40) return false;

        if (!isRealData()) plotSR(selectedDileptons[0],selectedDileptons[1]);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void getSREventCount(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.analyze();
    a.write(outFileName);
}
void getSREventCount(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
