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
    void makeTest(TString sn, const Lepton* lep1, const Lepton* lep2) {

    	auto testPlots = [&](TString sn, double ht, double dR, double dPhi, double mll, int numB) {
    		plotter.getOrMake1DPre(sn,"dr_ll",";#DeltaR_{l,l}",50,0,4)->Fill(dR,weight);
    		plotter.getOrMake1DPre(sn,"dPhi_metLL",";#Delta#phi_{ll,met}",50,-3.14,3.14)->Fill(dPhi,weight);
    		plotter.getOrMake1DPre(sn,"mll",";M_{ll}",50,80,100)->Fill(mll,weight);
    		plotter.getOrMake1DPre(sn,"numB",";numB",5,-0.5,4.5)->Fill(numB,weight);
    	};

    	double mll = (lep1->p4()+lep2->p4()).mass();
    	double dPhi_metLL = PhysicsUtilities::deltaPhi(lep1->p4()+lep2->p4(),reader_event->met);
    	double dR_ll = PhysicsUtilities::deltaR(*lep1,*lep2);

    	if (mll > 100 || mll < 80) return;

    	bool passDphi  = fabs(dPhi_metLL) < TMath::PiOver2();
    	bool passDR    = dR_ll < 1.6;
    	bool passBveto = (nMedBTags_HbbV == 0);

    	if (passDR)                          testPlots(sn+"_passDR",ht_puppi,dR_ll,dPhi_metLL,mll,nMedBTags_HbbV);
    	if (passDphi)                        testPlots(sn+"_passDphi",ht_puppi,dR_ll,dPhi_metLL,mll,nMedBTags_HbbV);
    	if (passBveto)                       testPlots(sn+"_passB",ht_puppi,dR_ll,dPhi_metLL,mll,nMedBTags_HbbV);
    	if (passDR && passDphi)              testPlots(sn+"_passDR_Dphi",ht_puppi,dR_ll,dPhi_metLL,mll,nMedBTags_HbbV);
    	if (passDR && passBveto)             testPlots(sn+"_passDR_B",ht_puppi,dR_ll,dPhi_metLL,mll,nMedBTags_HbbV);
    	if (passDphi && passBveto)           testPlots(sn+"_passDphi_B",ht_puppi,dR_ll,dPhi_metLL,mll,nMedBTags_HbbV);
    	if (passDR && passDphi && passBveto) testPlots(sn+"_passDR_Dphi_B",ht_puppi,dR_ll,dPhi_metLL,mll,nMedBTags_HbbV);
    }
    void makePlots(TString sn, const Lepton* lep1, const Lepton* lep2) {
    	// N-1 plots for discriminating variables
    	double mll = (lep1->p4()+lep2->p4()).mass();
    	double dPhi_metLL = PhysicsUtilities::deltaPhi(lep1->p4()+lep2->p4(),reader_event->met);
    	double dR2_ll = PhysicsUtilities::deltaR2(*lep1,*lep2);

    	bool passMLL  = mll > 12 && mll < 75;
    	bool passDphi = fabs(dPhi_metLL) < TMath::PiOver2();
    	bool passDR   = dR2_ll < 1.6*1.6;

    	if (passMLL && passDphi) plotter.getOrMake1DPre(sn,"dr_ll",";#DeltaR_{l,l}",50,0,4)->Fill(sqrt(dR2_ll),weight);
    	if (passMLL && passDR)   plotter.getOrMake1DPre(sn,"dPhi_metLL",";#Delta#phi_{ll,met}",50,-3.14,3.14)->Fill(dPhi_metLL,weight);
    	if (passDR && passDphi)  plotter.getOrMake1DPre(sn,"mll",";M_{ll}",50,0,120)->Fill(mll,weight);
    	if (passDR && passDphi && passMLL) {
    		plotter.getOrMake1DPre(sn+"_fullSel","dr_ll",";#DeltaR_{l,l}",50,0,4)->Fill(sqrt(dR2_ll),weight);
    		plotter.getOrMake1DPre(sn+"_fullSel","dPhi_metLL",";#Delta#phi_{ll,met}",50,-3.14,3.14)->Fill(dPhi_metLL,weight);
    		plotter.getOrMake1DPre(sn+"_fullSel","mll",";M_{ll}",50,0,120)->Fill(mll,weight);
    		plotter.getOrMake1DPre(sn+"_fullSel","ht",";H_{T}",50,0,3000)->Fill(ht_puppi,weight);
    		plotter.getOrMake1DPre(sn+"_fullSel","Mbb",";M_{bb}",20,30,210)->Fill(hbbMass,weight);
    		plotter.getOrMake1DPre(sn+"_fullSel","Mhh",";M_{HH}",30,0,4500)->Fill(hh.mass(),weight);
    		plotter.getOrMake1DPre(sn+"_fullSel","met",";E_{T}^{miss}",50,0,3000)->Fill(reader_event->met.pt(),weight);
    		plotter.getOrMake1DPre(sn+"_fullSel","pt2",";p_{T}",50,0,3000)->Fill(lep2->pt(),weight);
    		plotter.getOrMake1DPre(sn+"_fullSel","maxLepEta",";|#eta|",20,0,3)->Fill(lep1->absEta() > lep2->absEta() ? lep1->absEta() : lep2->absEta(),weight);
    	}
    }
    void plotTopCR(TString sn, const Lepton* lep1, const Lepton* lep2) {
    	makePlots(sn+"_TopCR_btagLMT",lep1,lep2);

    	TString chan = getDilepChan(lep1,lep2);
    	makePlots(sn+chan+"TopCR_btagLMT",lep1,lep2);

    	if (hbbCSVCat == BTagging::CSVSJ_MF) {
        	makePlots(sn+"_TopCR_btagL",lep1,lep2);
        	makePlots(sn+chan+"TopCR_btagL",lep1,lep2);
    	} else if (hbbCSVCat == BTagging::CSVSJ_ML) {
        	makePlots(sn+"_TopCR_btagM",lep1,lep2);
        	makePlots(sn+chan+"TopCR_btagM",lep1,lep2);
    	} else if (hbbCSVCat == BTagging::CSVSJ_MM) {
        	makePlots(sn+"_TopCR_btagT",lep1,lep2);
        	makePlots(sn+chan+"TopCR_btagT",lep1,lep2);
    	}
    }
    void plotSR(TString sn, const Lepton* lep1, const Lepton* lep2) {
    	makePlots(sn+"_SR_btagLMT",lep1,lep2);

    	TString chan = getDilepChan(lep1,lep2);
    	makePlots(sn+chan+"SR_btagLMT",lep1,lep2);

    	if (hbbCSVCat == BTagging::CSVSJ_MF) {
        	makePlots(sn+"_SR_btagL",lep1,lep2);
        	makePlots(sn+chan+"SR_btagL",lep1,lep2);
    	} else if (hbbCSVCat == BTagging::CSVSJ_ML) {
        	makePlots(sn+"_SR_btagM",lep1,lep2);
        	makePlots(sn+chan+"SR_btagM",lep1,lep2);
    	} else if (hbbCSVCat == BTagging::CSVSJ_MM) {
        	makePlots(sn+"_SR_btagT",lep1,lep2);
        	makePlots(sn+chan+"SR_btagT",lep1,lep2);
    	}
    }
    void plotQgCR(TString sn, const Lepton* lep1, const Lepton* lep2) {
    	makePlots(sn+"_QgCR",lep1,lep2);

    	TString chan = getDilepChan(lep1,lep2);
    	makePlots(sn+chan+"QgCR",lep1,lep2);
    }

    void testQgCR(TString sn, const Lepton* lep1, const Lepton* lep2) {
    	makeTest(sn+"_QgCR",lep1,lep2);

    	TString chan = getDilepChan(lep1,lep2);
    	makeTest(sn+chan+"QgCR",lep1,lep2);
    }

	TString getDilepChan(const Lepton* lep1, const Lepton* lep2) {
		if (lep1->isMuon() && lep2->isMuon()) return "_mumu_";
		else if (lep1->isElectron() && lep2->isElectron()) return "_ee_";
		else return "_emu_";
	}

    bool runEvent() override {
//    	cout<<"jansen"<<endl;
        if(!DileptonSearchRegionAnalyzer::runEvent()) return false;
//    	cout<<"jansen01"<<endl;

        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        if(!passEventFilters) return false;
//    	cout<<"jansen02"<<endl;

        if(!passTriggerPreselection) return false;
//    	cout<<"jansen03"<<endl;

        TString sn = smpName;

        // cuts before separating into SR and CR
        if(!hbbCand) return false;
        if(selectedDileptons.size() != 2) return false;
        if (getDilepChan(selectedDileptons[0],selectedDileptons[1]).Contains("_ee_")) {
        	if (!(((const Electron*)selectedDileptons[0])->passMedID_noISO() && ((const Electron*)selectedDileptons[1])->passMedID_noISO())) return false;
        }
        if (hbbMass < 30 || hbbMass > 210) return false;
        if (ht_puppi < 400) return false;
        if (hh.mass() < 700) return false;

        // trying out straight MET cut
//        if(reader_event->met.pt() < 50) return false;
//    	cout<<"jansen2"<<endl;

        // separate into SR and CR
        if (nMedBTags_HbbV != 0 && hbbCSVCat >= BTagging::CSVSJ_MF)                  plotTopCR (sn, selectedDileptons[0],selectedDileptons[1]);
        if (nMedBTags_HbbV == 0 && hbbCSVCat == BTagging::CSVSJ_FF)                  plotQgCR  (sn, selectedDileptons[0],selectedDileptons[1]);
        if (nMedBTags_HbbV == 0 && hbbCSVCat >= BTagging::CSVSJ_MF && !isRealData()) plotSR    (sn, selectedDileptons[0],selectedDileptons[1]);

        // debugging
//        if (hbbCSVCat == BTagging::CSVSJ_FF) testQgCR (sn, selectedDileptons[0],selectedDileptons[1]);
//    	cout<<"jansen3"<<endl;

        double DR1 = PhysicsUtilities::deltaR(*hbbCand,*selectedDileptons[0]);
        double DR2 = PhysicsUtilities::deltaR(*hbbCand,*selectedDileptons[1]);
        plotter.getOrMake1DPre(sn,"minDR_lepHbb","",200,0,4)->Fill(DR1 > DR2 ? DR2 : DR1,weight);
        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void getSRandCRVars(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.analyze();
    a.write(outFileName);
}
void getSRandCRVars(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
