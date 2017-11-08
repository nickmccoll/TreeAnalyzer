
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "TreeReaders/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/EventSelection/interface/EventSelection.h"

#include "TPRegexp.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt)
{

}




    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;


        if(reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)
            smpName = "other";

        if(!passTriggerPreselection) return false;
        if(!passEventFilters) return false;
        if(selectedLeptons.size() != 1) return false;
        if(!wjjCand) return false;
        if(!hbbCand) return false;
        if(nMedBTags_HbbV) return false;



        auto pltSet = [&](const TString& prefix){
            plotter.getOrMake1DPre(prefix,"hh_mass",";HH mass [TeV]; N. events / 5 GeV",1000,0.,5.)->Fill(hh.mass() / 1000.,weight);
            plotter.getOrMake1DPre(prefix, "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 700 && hh.mass() < 900)
                plotter.getOrMake1DPre(prefix+"_hh700to900", "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 900 && hh.mass() <1100)
                plotter.getOrMake1DPre(prefix+"_hh900to1100", "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 1400 && hh.mass() < 1800)
                plotter.getOrMake1DPre(prefix+"_hh1400to1800", "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
            if(hh.mass() >= 2500 && hh.mass() < 3500)
                plotter.getOrMake1DPre(prefix+"_hh2500to3000", "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbbMass,weight);
        };

        auto pltRecoSet = [&](const TString& prefix){
            pltSet(prefix + "_lW");
            if(passHbbSel) pltSet(prefix +"_lW_LI");
            if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MF) pltSet(prefix +"_lW_L");
            if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_ML) pltSet(prefix +"_lW_M");
            if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MM) pltSet(prefix +"_lW_T");

            if(passWjjSel && passWlnuDR && passWWDM){
                pltSet(prefix);
                if(passHbbSel) pltSet(prefix +"_LI");
                if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MF) pltSet(prefix +"_L");
                if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_ML) pltSet(prefix +"_M");
                if(passHbbSel && hbbCSVCat == BTagging::CSVSJ_MM) pltSet(prefix +"_T");
            }
        };

        auto pltLepSet = [&](const TString& prefix){
            pltRecoSet(prefix + "_emu");
            if(selectedLepton->isMuon()) pltRecoSet(prefix+"_mu");
            else pltRecoSet(prefix+"_e");
        };

        bool w_in = false;
        bool wqq_in = false;
        bool wb_in = false;
        bool wqqb_in = false;


        bool q_in = false;
        bool qL_in = false;


        const float matchR = 0.8*0.8;
        const float matchLR = 1.0;

        for(const auto& d : smDecayEvt.bosonDecays  ){
            if(d.type != BosonDecay::Z_HAD && d.type != BosonDecay::W_HAD ) continue;
            if(PhysicsUtilities::deltaR2(*d.boson,*hbbCand) < matchR) w_in = true;
            int nQIn = (PhysicsUtilities::deltaR2(*d.dau1,*hbbCand) < matchR) + (PhysicsUtilities::deltaR2(*d.dau2,*hbbCand) < matchR);
            int nQLIn = (PhysicsUtilities::deltaR2(*d.dau1,*hbbCand) < matchLR) + (PhysicsUtilities::deltaR2(*d.dau2,*hbbCand) < matchLR);
            if(nQIn == 2 ){
                wqq_in = true;
            }
            if(nQIn > 0) q_in = true;
            if(nQLIn > 0) qL_in = true;
        }




        for(const auto& d : smDecayEvt.topDecays  ){
            if(d.type != TopDecay::HAD ) continue;
            bool passB = PhysicsUtilities::deltaR2(*d.b,*hbbCand) < matchR;
            bool passBL = PhysicsUtilities::deltaR2(*d.b,*hbbCand) < matchLR;
            int nWIn = (PhysicsUtilities::deltaR2(*d.W_decay.dau1,*hbbCand) < matchR) + (PhysicsUtilities::deltaR2(*d.W_decay.dau2,*hbbCand) < matchR);

            int nWLIn = (PhysicsUtilities::deltaR2(*d.W_decay.dau1,*hbbCand) < matchLR) + (PhysicsUtilities::deltaR2(*d.W_decay.dau2,*hbbCand) < matchLR);
            if(PhysicsUtilities::deltaR2(*d.W_decay.boson,*hbbCand) < matchR ){
                w_in = true;
                if(passB) wb_in = true;
            }
            if(nWIn == 2){
                wqq_in = true;
                if(passB) wqqb_in = true;
            }

            int nQIn = nWIn + passB;
            int nQLIn = nWLIn + passBL;
            if(nQIn > 0) q_in = true;
            if(nQLIn > 0) qL_in = true;

        }

        auto pltTruthSet= [&](const TString& prefix){
            pltLepSet(prefix + "_incl");
            if(w_in) pltLepSet(prefix + "_w");
            if(!w_in) pltLepSet(prefix + "_nw");
            if(w_in && !wb_in) pltLepSet(prefix + "_wnwb");
            if(wb_in) pltLepSet(prefix + "_wb");

            if(wqq_in) pltLepSet(prefix + "_wqq");
            if(!wqq_in) pltLepSet(prefix + "_nwqq");
            if(wqq_in && !wqqb_in) pltLepSet(prefix + "_wqqnwqqb");
            if(wqqb_in) pltLepSet(prefix + "_wqqb");
            if(!wqqb_in) pltLepSet(prefix + "_nwqqb");

            if(qL_in) pltLepSet(prefix + "_wqL");
            if(!qL_in) pltLepSet(prefix + "_nwqL");

            if(q_in) pltLepSet(prefix + "_wq");
            if(!q_in) pltLepSet(prefix + "_nwq");
            if(q_in && !wqqb_in) pltLepSet(prefix + "_wqnwqqb");
        };


        pltTruthSet(smpName);
        pltTruthSet("bkg");
        if(reader_event->process != FillerConstants::QCD) pltTruthSet("bkgNoQCD");




        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;


};

#endif

void testResCats(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void testResCats(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
