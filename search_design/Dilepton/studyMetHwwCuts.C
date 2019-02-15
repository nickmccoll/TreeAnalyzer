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
    double getnunuZ(const Lepton* lep1, const Lepton* lep2, double m_nunu) {
    	MomentumF llMOM = lep1->p4() + lep2->p4();
    	MomentumF met = reader_event->met.p4();

    	double A = (0.5*(125*125 - m_nunu*m_nunu - llMOM.mass()*llMOM.mass()) + met.px()*llMOM.px() + met.py()*llMOM.py()) / llMOM.E();
    	double B = met.pt()*met.pt();
    	double C = llMOM.pz() / llMOM.E();
    	double delta = 4*A*A*C*C - 4*(C*C-1)*(A*A-B);

    	double nunuZ = 0;
    	if (delta < 0) {
    		nunuZ = -A*C / (C*C-1);
    	} else {
    		double pos = (-2*A*C + sqrt(delta)) / (2*C*C-2);
    		double neg = (-2*A*C - sqrt(delta)) / (2*C*C-2);
    		if (fabs(pos) < fabs(neg)) nunuZ = pos;
    		else nunuZ = neg;
    	}
    	return nunuZ;
    }
    void plotAltHwwReco(TString sn, const Lepton* lep1, const Lepton* lep2) {
    	// make theta_nunu = theta_ll to get pz component of neutrinos with Mnunu = 0
    	double pznunu_thetaLL = reader_event->met.pt() / TMath::Tan((lep1->p4()+lep2->p4()).theta());
    	ASTypes::CartLorentzVector nunuMOM(reader_event->met.px(),reader_event->met.py(),pznunu_thetaLL,sqrt(reader_event->met.px()*reader_event->met.px()+reader_event->met.py()*reader_event->met.py()+pznunu_thetaLL*pznunu_thetaLL+40*40));
    	double hhmass_thetaLL = (lep1->p4() + lep2->p4() + nunuMOM + hbbCand->p4()).mass();

    	// choose some nu-nu invariant mass, then solve for nu-nu pz given Higgs mass constraint
    	double m_nunu = 45;
    	double pz_nunu = getnunuZ(lep1,lep2,m_nunu);
    	ASTypes::CartLorentzVector nunu(reader_event->met.px(),reader_event->met.py(),pz_nunu,sqrt(reader_event->met.px()*reader_event->met.px()+reader_event->met.py()*reader_event->met.py()+pz_nunu*pz_nunu+m_nunu*m_nunu));
    	double hhmass_mnunu = (lep1->p4() + lep2->p4() + nunu + hbbCand->p4()).mass();

    	plotter.getOrMake1DPre(sn,"hhmass_thetaLL",";GeV",100,0,4000)->Fill(hhmass_thetaLL,weight);
    	plotter.getOrMake1DPre(sn,"hhmass_mnunu",";GeV",100,0,4000)->Fill(hhmass_mnunu,weight);

    	MomentumF metLL = reader_event->met.p4() + lep1->p4() + lep2->p4();
    	if (hhmass_thetaLL > 700) {
        	plotter.getOrMake1DPre(sn,"met2",";GeV",100,0,3000)->Fill(reader_event->met.pt(),weight);
        	plotter.getOrMake1DPre(sn,"metLL2",";GeV",100,0,2000)->Fill(metLL.Et(),weight);
    	}
    	if (hhmass_mnunu > 700) {
        	plotter.getOrMake1DPre(sn,"met3",";GeV",100,0,3000)->Fill(reader_event->met.pt(),weight);
        	plotter.getOrMake1DPre(sn,"metLL3",";GeV",100,0,2000)->Fill(metLL.Et(),weight);
    	}
    }
    void plotGen(TString sn, const Lepton *lep1, const Lepton *lep2) {
    	double mnunu_gen = (diHiggsEvt.w1_d2->p4() + diHiggsEvt.w2_d2->p4()).mass();
    	double mnunu_reco = (hwwInfo.nu1 + hwwInfo.nu2).mass();
    	double mll_gen = (diHiggsEvt.w1_d1->p4() + diHiggsEvt.w2_d1->p4()).mass();
    	double mll_reco = (lep1->p4() + lep2->p4()).mass();

    	plotter.getOrMake1DPre(sn,"Mh_min_mll_min_mnunu_GEN",";mass",100,0,125)->Fill(125-mll_gen-mnunu_gen,weight);
    	plotter.getOrMake1DPre(sn,"Mh_min_mll_min_mnunu_RECO",";mass",100,0,125)->Fill(125-mll_reco-mnunu_gen,weight);
    	plotter.getOrMake1DPre(sn,"Mnunu_reco",";mass",100,0,100)->Fill(mnunu_reco,weight);
    	plotter.getOrMake2DPre(sn,"Mnunu_x_Mll",";Mll;Mnunu",50,0,100,50,0,100)->Fill(mll_gen,mnunu_gen,weight);
    }

    void plotMet(TString sn, const Lepton *lep1, const Lepton *lep2) {
    	plotter.getOrMake1DPre(sn,"hhmass",";GeV",100,0,4000)->Fill(hh.mass(),weight);

    	if (hh.mass() < 700) return;
    	MomentumF nunu = hwwInfo.nu1 + hwwInfo.nu2;
    	MomentumF metLL = reader_event->met.p4() + lep1->p4() + lep2->p4();
    	plotter.getOrMake1DPre(sn,"met",";GeV",100,0,3000)->Fill(reader_event->met.pt(),weight);
    	plotter.getOrMake1DPre(sn,"nunu_pt",";GeV",100,0,2000)->Fill(nunu.pt(),weight);
    	plotter.getOrMake1DPre(sn,"nu1_pt",";GeV",100,0,2000)->Fill(hwwInfo.nu1.Pt(),weight);
    	plotter.getOrMake1DPre(sn,"nu2_pt",";GeV",100,0,2000)->Fill(hwwInfo.nu2.Pt(),weight);

    	plotter.getOrMake1DPre(sn,"nunu_mass",";GeV",50,0,100)->Fill(nunu.mass(),weight);
    	plotter.getOrMake1DPre(sn,"metLL",";GeV",100,0,2000)->Fill(metLL.Et(),weight);
    	plotter.getOrMake1DPre(sn,"testStat","",100,0,20)->Fill(HwwTestStat,weight);
    }

	TString getDilepChan(const Lepton* lep1, const Lepton* lep2) {
		if (lep1->isMuon() && lep2->isMuon()) return "_mumu_";
		else if (lep1->isElectron() && lep2->isElectron()) return "_ee_";
		else return "_emu_";
	}

    bool runEvent() override {
        if(!DileptonSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type != DiHiggsEvent::DILEP) return false;
        if(!hbbCand) return false;
        if(selectedDileptons.size() != 2) return false;
        if(hbbCSVCat < BTagging::CSVSJ_MF) return false;
        if (hbbMass < 30 || hbbMass > 210) return false;
        if (ht_puppi < 400) return false;
        if (nMedBTags_HbbV != 0) return false;
        if (PhysicsUtilities::deltaR2(*selectedDileptons[0],*selectedDileptons[1]) > 1.6*1.6) return false;

        // assuming the input dataset is a skim using the Dilepton + Hbb selection
        TString sn = smpName + getDilepChan(selectedDileptons[0],selectedDileptons[1]);
        if (sn.Contains("_ee_")) {
        	if (!((const Electron*)selectedDileptons[0])->passMedID_noISO() && ((const Electron*)selectedDileptons[1])->passMedID_noISO()) return false;
        }

        double mll = (selectedDileptons[0]->p4()+selectedDileptons[1]->p4()).mass();
        plotter.getOrMake1DPre(sn,"mll",";M_{ll}",100,0,200)->Fill(mll,weight);
        if (mll > 75 || mll < 12) return false;

        if (isSignal()) plotGen(sn,selectedDileptons[0],selectedDileptons[1]);
        plotAltHwwReco(sn,selectedDileptons[0],selectedDileptons[1]);
        plotMet(sn,selectedDileptons[0],selectedDileptons[1]);
        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void studyMetHwwCuts(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.analyze();
    a.write(outFileName);
}
void studyMetHwwCuts(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
