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
    void plotGen(TString sn, const GenParticle *lep1, const GenParticle *lep2, const GenParticle *nu1, const GenParticle *nu2) {

        // neutrinos
        plotter.getOrMake1DPre(sn, "nu_pt", ";P_{T} (GeV)",100,0,1000)->Fill(nu1->pt(), weight);
        plotter.getOrMake1DPre(sn, "nu_pt", ";P_{T} (GeV)",100,0,1000)->Fill(nu2->pt(), weight);
        plotter.getOrMake1DPre(sn, "nu_pz", ";P_{z} (GeV)",100,0,1000)->Fill(nu1->pz(), weight);
        plotter.getOrMake1DPre(sn, "nu_pz", ";P_{z} (GeV)",100,0,1000)->Fill(nu2->pz(), weight);
        plotter.getOrMake1DPre(sn, "nu_deltapz", ";#DeltaP_{z} (GeV)",100,0,100)->Fill(abs(nu1->pz()-nu2->pz()), weight);
        plotter.getOrMake1DPre(sn, "nu_totalpz", ";Total P_{z} (GeV)",100,0,1000)->Fill(nu1->pz()+nu2->pz(), weight);
        plotter.getOrMake1DPre(sn, "M_nunu", ";M_{#nu#nu} (GeV)",100,0,400)->Fill((nu1->p4()+nu2->p4()).mass(), weight);

        // W bosons
        double mw1 = diHiggsEvt.w1 ? diHiggsEvt.w1->mass() : (lep1->p4() + nu1->p4()).mass();
        double mw2 = diHiggsEvt.w2 ? diHiggsEvt.w2->mass() : (lep2->p4() + nu2->p4()).mass();

        plotter.getOrMake1DPre(sn, "Mw", ";M_{W} (GeV)",100,0,120)->Fill(mw1, weight);
        plotter.getOrMake1DPre(sn, "Mw", ";M_{W} (GeV)",100,0,120)->Fill(mw2, weight);

        // MET and neutrinos
        plotter.getOrMake1DPre(sn, "nu_totalpx", ";P_{x} (GeV)",100,-500,500)->Fill(nu1->px()+nu2->px(), weight);
        plotter.getOrMake1DPre(sn, "nu_totalpy", ";P_{y} (GeV)",100,-500,500)->Fill(nu1->py()+nu2->py(), weight);
        plotter.getOrMake1DPre(sn, "delta_Misspx", ";#DeltaP_{x} (GeV)",100,0,200)->Fill(abs(nu1->px()+nu2->px()-reader_event->met.px()), weight);
        plotter.getOrMake1DPre(sn, "delta_Misspy", ";#DeltaP_{y} (GeV)",100,0,200)->Fill(abs(nu1->py()+nu2->py()-reader_event->met.py()), weight);

        MomentumF nunuMOM = nu1->p4() + nu2->p4();
        MomentumF llMOM = lep1->p4() + lep2->p4();
        plotter.getOrMake1DPre(sn, "dTheta_ll_nunu", ";",100,-1,1)->Fill(nunuMOM.theta()-llMOM.theta(), weight);

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

        if (reader_event->met.pt() < 40) return false;
        if (mll > 75 || mll < 12) return false;
        if (fabs(PhysicsUtilities::deltaPhi(reader_event->met.p4(),
        		selectedDileptons[0]->p4()+selectedDileptons[1]->p4())) > TMath::PiOver2() ) return false;

        // assuming the input dataset is a skim using the Dilepton + Hbb selection
        TString sn = smpName+dilepMap[dilepChan];

        plotGen(sn,diHiggsEvt.w1_d1,diHiggsEvt.w2_d1,diHiggsEvt.w1_d2,diHiggsEvt.w2_d2);
        plotGen(smpName,diHiggsEvt.w1_d1,diHiggsEvt.w2_d1,diHiggsEvt.w1_d2,diHiggsEvt.w2_d2);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void plotGenDilepton(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.analyze();
    a.write(outFileName);
}
void plotGenDilepton(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
