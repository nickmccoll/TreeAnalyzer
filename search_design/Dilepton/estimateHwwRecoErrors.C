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
    void plot(TString sn, const GenParticle *lep1, const GenParticle *lep2, const GenParticle *nu1, const GenParticle *nu2) {

        // W bosons
        double mw1 = diHiggsEvt.w1 ? diHiggsEvt.w1->mass() : (lep1->p4() + nu1->p4()).mass();
        double mw2 = diHiggsEvt.w2 ? diHiggsEvt.w2->mass() : (lep2->p4() + nu2->p4()).mass();

        plotter.getOrMake1DPre(sn, "Mw", ";M_{W} (GeV)",100,0,120)->Fill(mw1, weight);
        plotter.getOrMake1DPre(sn, "Mw", ";M_{W} (GeV)",100,0,120)->Fill(mw2, weight);

        // Higgs mass with RECO dileptons and GEN neutrinos
        const MomentumF pH = selectedDileptons[0]->p4() + selectedDileptons[1]->p4() + nu1->p4() + nu2->p4();
        plotter.getOrMake1DPre(sn, "Mh", ";M_{H} (GeV)",100,50,200)->Fill(pH.mass(), weight);

        double metPerp, metPara, nunuPerp, nunuPara;
        setNeutrinoMetProjections(metPerp,metPara,nunuPerp,nunuPara,nu1,nu2);

        plotter.getOrMake1DPre(sn, "deltaMiss_perpNORM", ";",100,-5,5)->Fill((nunuPerp-metPerp)/pH.pt(), weight);
        plotter.getOrMake1DPre(sn, "deltaMiss_paraNORM", ";",100,-4,4)->Fill((nunuPara-metPara)/pH.pt(), weight);

        const MomentumF nunu = nu1->p4() + nu2->p4();
        const MomentumF ll = lep1->p4() + lep2->p4();
        plotter.getOrMake1DPre(sn, "dR_llvv", ";",100,0,5)->Fill(PhysicsUtilities::deltaR(nunu,ll), weight);
        plotter.getOrMake1DPre(sn, "dPhi_llvv", ";",100,-3.14,3.14)->Fill(PhysicsUtilities::deltaPhi(nunu,ll), weight);
        plotter.getOrMake1DPre(sn, "dTheta_llvv", ";",100,-3.14,3.14)->Fill(nunu.theta()-ll.theta(), weight);

        plotter.getOrMake1DPre(sn, "Delta_dR_lvlv", ";",100,0,5)->Fill(PhysicsUtilities::deltaR(*lep1,*nu1) - PhysicsUtilities::deltaR(*lep2,*nu2), weight);
        plotter.getOrMake1DPre(sn, "Delta_dPhi_lvlv", ";",100,-3.14,3.14)->Fill(PhysicsUtilities::deltaPhi(*lep1,*nu1) - PhysicsUtilities::deltaPhi(*lep2,*nu2), weight);

        double theta_lv1 = (lep1->p4()+nu1->p4()).theta();
        double theta_lv2 = (lep2->p4()+nu2->p4()).theta();
        double dTheta_lv1 = lep1->theta()-nu1->theta();
        double dTheta_lv2 = lep2->theta()-nu2->theta();
        plotter.getOrMake1DPre(sn, "Delta_thetaLV", ";",100,-3.14,3.14)->Fill(theta_lv1-theta_lv2, weight);
        plotter.getOrMake1DPre(sn, "dTheta_lv_ON", ";",100,-3.14,3.14)->Fill(mw1 > mw2 ? dTheta_lv1 : dTheta_lv2, weight);
        plotter.getOrMake1DPre(sn, "dTheta_lv_OFF", ";",100,-3.14,3.14)->Fill(mw1 > mw2 ? dTheta_lv2 : dTheta_lv1, weight);
        plotter.getOrMake1DPre(sn, "Delta_dthetaLV", ";",100,-3.14,3.14)->Fill(dTheta_lv1-dTheta_lv2, weight);
        plotter.getOrMake1DPre(sn, "Theta_hww", ";",100,0,3.14)->Fill((lep1->p4()+lep2->p4()+nu1->p4()+nu2->p4()).theta(), weight);

        plotter.getOrMake1DPre(sn, "vW_Pratio_ON", ";",100,0,1)->Fill(mw1 > mw2 ? nu1->p()/(lep1->p4()+nu1->p4()).P() : nu2->p()/(lep2->p4()+nu2->p4()).P(), weight);
        plotter.getOrMake1DPre(sn, "vW_Pratio_OFF", ";",100,0,1)->Fill(mw1 > mw2 ? nu2->p()/(lep2->p4()+nu2->p4()).P() : nu1->p()/(lep1->p4()+nu1->p4()).P(), weight);

    }

    void setNeutrinoMetProjections(double &metPerp, double &metPara, double &nunuPerp, double &nunuPara, const GenParticle *nu1, const GenParticle *nu2) {
    	double axisX = selectedDileptons[0]->px() + selectedDileptons[1]->px() + reader_event->met.px();
    	double axisY = selectedDileptons[0]->py() + selectedDileptons[1]->py() + reader_event->met.py();

    	axisX = axisX / sqrt(axisX*axisX + axisY*axisY);
    	axisY = axisY / sqrt(axisX*axisX + axisY*axisY);
    	double axisX_perp = (-1)*axisY;
    	double axisY_perp = axisX;

    	auto getPerp = [&](const double px, const double py)->double{return px*axisX_perp + py*axisY_perp;};
    	auto getPara = [&](const double px, const double py)->double{return px*axisX + py*axisY;};

    	metPerp = getPerp(reader_event->met.px(),reader_event->met.py());
    	metPara = getPara(reader_event->met.px(),reader_event->met.py());
    	nunuPerp = getPerp(nu1->px() + nu2->px(), nu1->py() + nu2->py());
    	nunuPara = getPara(nu1->px() + nu2->px(), nu1->py() + nu2->py());
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

        // kinematic cuts
        if (PhysicsUtilities::deltaR2(*selectedDileptons[0],*selectedDileptons[1]) > 1.6*1.6) return false;
        double mll = (selectedDileptons[0]->p4()+selectedDileptons[1]->p4()).mass();
        if (mll > 75 || mll < 12) return false;

        // assuming the input dataset is a skim using the Dilepton + Hbb selection
        TString sn = smpName;

        plot(sn,diHiggsEvt.w1_d1,diHiggsEvt.w2_d1,diHiggsEvt.w1_d2,diHiggsEvt.w2_d2);
        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void estimateHwwRecoErrors(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.analyze();
    a.write(outFileName);
}
void estimateHwwRecoErrors(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
