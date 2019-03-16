
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Configuration/interface/FillerConstants.h"
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

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){


        DefaultFatJetSelections::setDefaultFatJetProcessor(fjProc_inclWBtag);

        fjProc_inclWBtag.param.wjj_maxCSVWP =BTagging::CSV_INCL;
        DefaultFatJetSelections::setDefaultFatJetProcessor(fjProc_inclHbbBtag);
        fjProc_inclHbbBtag.param.hbb_l_firMinCSVWP =BTagging::CSV_INCL;
        fjProc_inclHbbBtag.param.hbb_l_secMinCSVWP =BTagging::CSV_INCL;
    }

    void makeCRPlots(TString pre, const MomentumF& hWW, const MomentumF& hbb){
        auto hh = hWW.p4() + hbb.p4();
        plotter.getOrMake1DPre(pre, "hh_mass" ,";hh mass [TeV]; arbitrary units",50,0,5)->Fill(hh.mass()/1000,weight);
        if(hbb.mass() < 90)
            plotter.getOrMake1DPre(pre, "hbblt90_hh_mass" ,";hh mass [TeV]; arbitrary units",50,0,5)->Fill(hh.mass()/1000,weight);
        else if(hbb.mass() < 140)
            plotter.getOrMake1DPre(pre, "hbb90to140_hh_mass" ,";hh mass [TeV]; arbitrary units",50,0,5)->Fill(hh.mass()/1000,weight);
        else
            plotter.getOrMake1DPre(pre, "hbbgeq140_hh_mass" ,";hh mass [TeV]; arbitrary units",50,0,5)->Fill(hh.mass()/1000,weight);


        plotter.getOrMake1DPre(pre, "hbb_mass" ,";h(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbb.mass(),weight);
        if(hh.mass() >= 700 && hh.mass() < 900)
            plotter.getOrMake1DPre(pre, "hh700to900_hbb_mass" ,";h(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbb.mass(),weight);
        if(hh.mass() >= 900 && hh.mass() < 1100)
            plotter.getOrMake1DPre(pre, "hh900to1100_hbb_mass" ,";h(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbb.mass(),weight);
        if(hh.mass() >= 1400 && hh.mass() < 1800)
            plotter.getOrMake1DPre(pre, "hh1400to1800_hbb_mass" ,";h(bb) mass [GeV]; arbitrary units",50,0,250)->Fill(hbb.mass(),weight);
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;

        if(reader_event->process !=
                FillerConstants::SIGNAL){
            prefix = (reader_event->process >= FillerConstants::ZJETS &&  reader_event->process != FillerConstants::QCD)  ? "rare" :
                    FillerConstants::MCProcessNames[reader_event->process];
        }
        if(!passTriggerPreselection) return false;
        if(!passEventFilters) return false;
        if(selectedLeptons.size() != 1) return false;

        auto lepton = selectedLepton->p4();

        if(reader_event->process == FillerConstants::SIGNAL) weight *=  20.0/1000.0;

        const auto* hbbfj = fjProc->getHBBCand();
        const auto* wjjfj = fjProc->getWjjCand();
        const bool goodHBBFJ  = fjProc->passHbbSel();
        const bool goodHBBFJT = fjProc->passHbbSelTightBTag();
        const bool goodWJJFJ  = fjProc->passWjjSel();



        fjProc_inclHbbBtag.loadFatJets(*reader_fatjet,selectedLepton);
        fjProc_inclWBtag.loadFatJets(*reader_fatjet,selectedLepton);

        if(!wjjfj || !hbbfj) return false;
        auto neutrino = HiggsSolver::getInvisible(reader_event->met,(lepton + wjjfj->sdMom().p4()) );
        MomentumF Wlnu = lepton + neutrino.p4();
        MomentumF hWW = Wlnu.p4() + wjjfj->sdMom().p4();
        if(PhysicsUtilities::deltaR(Wlnu, wjjfj->sdMom().p4()) > 0.5) return false;

        if(goodWJJFJ && goodHBBFJ) makeCRPlots(prefix+"_stdWjj_stdHBB",hWW,hbbfj->sdMom());
        if(goodWJJFJ && goodHBBFJT) makeCRPlots(prefix+"_stdWjj_stdHBBT",hWW,hbbfj->sdMom());


        auto minSJCSV = [](const FatJet* fj)->float{
            if(!fj->nSubJets()) return 0;
            float minCSV = 1;
            for(const auto& sj : fj->subJets()) { minCSV = std::min(minCSV,sj.csv());}
            return std::max(minCSV,float(0.0));
        };
        auto  maxSJCSV = [](const FatJet* fj)->float{
            float maxCSV = 0;
            for(const auto& sj : fj->subJets()) { maxCSV = std::max(maxCSV,sj.csv());}
            return maxCSV;
        };




        if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && maxSJCSV(hbbfj) >= BTagging::BBT_VALS[BTagging::CSV_M] && minSJCSV(hbbfj) < BTagging::BBT_VALS[BTagging::CSV_L]  )
            makeCRPlots(prefix+"_stdWjj_oneBHBB",hWW,hbbfj->sdMom());
        if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && maxSJCSV(hbbfj) < BTagging::BBT_VALS[BTagging::CSV_L])
            makeCRPlots(prefix+"_stdWjj_noBHBB",hWW,hbbfj->sdMom());
        if(fjProc_inclWBtag.passWjjSel() && maxSJCSV(wjjfj)>= BTagging::BBT_VALS[BTagging::CSV_T] && goodHBBFJ) makeCRPlots(prefix+"_oneBWjj_stdHBB",hWW,hbbfj->sdMom());
        if(fjProc_inclWBtag.passWjjSel() && maxSJCSV(wjjfj)>= BTagging::BBT_VALS[BTagging::CSV_T]  && goodHBBFJT) makeCRPlots(prefix+"_oneBWjj_stdHBBT",hWW,hbbfj->sdMom());

        if(selectedLepton->isMuon()){
            TString tempP = prefix + "_mu";
            if(goodWJJFJ && goodHBBFJ) makeCRPlots(tempP+"_stdWjj_stdHBB",hWW,hbbfj->sdMom());
            if(goodWJJFJ && goodHBBFJT) makeCRPlots(tempP+"_stdWjj_stdHBBT",hWW,hbbfj->sdMom());
            if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && maxSJCSV(hbbfj) >= BTagging::BBT_VALS[BTagging::CSV_M] && minSJCSV(hbbfj)  < BTagging::BBT_VALS[BTagging::CSV_L]  )
                makeCRPlots(tempP+"_stdWjj_oneBHBB",hWW,hbbfj->sdMom());
            if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && maxSJCSV(hbbfj) < BTagging::BBT_VALS[BTagging::CSV_L])
                makeCRPlots(tempP+"_stdWjj_noBHBB",hWW,hbbfj->sdMom());
            if(fjProc_inclWBtag.passWjjSel() && maxSJCSV(wjjfj)>= BTagging::BBT_VALS[BTagging::CSV_T] && goodHBBFJ) makeCRPlots(tempP+"_oneBWjj_stdHBB",hWW,hbbfj->sdMom());
            if(fjProc_inclWBtag.passWjjSel() && maxSJCSV(wjjfj)>= BTagging::BBT_VALS[BTagging::CSV_T]  && goodHBBFJT) makeCRPlots(tempP+"_oneBWjj_stdHBBT",hWW,hbbfj->sdMom());
        } else {
            TString tempP = prefix + "_el";
            if(goodWJJFJ && goodHBBFJ) makeCRPlots(tempP+"_stdWjj_stdHBB",hWW,hbbfj->sdMom());
            if(goodWJJFJ && goodHBBFJT) makeCRPlots(tempP+"_stdWjj_stdHBBT",hWW,hbbfj->sdMom());
            if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && maxSJCSV(hbbfj) >= BTagging::BBT_VALS[BTagging::CSV_M] && minSJCSV(hbbfj) < BTagging::BBT_VALS[BTagging::CSV_L]  )
                makeCRPlots(tempP+"_stdWjj_oneBHBB",hWW,hbbfj->sdMom());
            if(goodWJJFJ &&  fjProc_inclHbbBtag.passHbbSel() && maxSJCSV(hbbfj) < BTagging::BBT_VALS[BTagging::CSV_L])
                makeCRPlots(tempP+"_stdWjj_noBHBB",hWW,hbbfj->sdMom());
            if(fjProc_inclWBtag.passWjjSel() && maxSJCSV(wjjfj)>= BTagging::BBT_VALS[BTagging::CSV_T] && goodHBBFJ) makeCRPlots(tempP+"_oneBWjj_stdHBB",hWW,hbbfj->sdMom());
            if(fjProc_inclWBtag.passWjjSel() && maxSJCSV(wjjfj)>= BTagging::BBT_VALS[BTagging::CSV_T]  && goodHBBFJT) makeCRPlots(tempP+"_oneBWjj_stdHBBT",hWW,hbbfj->sdMom());
        }


        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

//    float   ht         =0;
//    float   selLep_pt  =0;
//    float   selLep_eta =0;
//    float   selLep_phi =0;
//    size8   selLep_muon=0;
    HistGetter plotter;
    TString prefix;

    FatJetProcessor fjProc_inclWBtag      ;
    FatJetProcessor fjProc_inclHbbBtag ;



};

#endif

void checkControlRegions(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void checkControlRegions(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
