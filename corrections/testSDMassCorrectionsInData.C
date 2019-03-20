
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "Configuration/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "Processors/Variables/interface/FatJetSelection.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"


#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
        turnOffCorr(CORR_SDMASS);
        fromWCorr.reset(new SoftDropMassScaleFactorsFromW (dataDirectory));
    }

    void makePlots(const TString& prefix, const float hbbMass, const float weight){
        plotter.getOrMake1DPre(prefix,"njets",";number of jets",20,-0.5,19.5)->Fill(jets.size(),weight );
        plotter.getOrMake1DPre(prefix,"njetswlep",";number of jets (with leptons)",20,-0.5,19.5)->Fill(jets_chs.size(),weight );
        plotter.getOrMake1DPre(prefix,"ht_wlep",";#it{H}_{T} (with leptons) [GeV]",500,500,3000)->Fill(ht_chs,weight );
        plotter.getOrMake1DPre(prefix,"ht_wlep_plusmet",";#it{H}_{T} (with leptons) + #slash{#it{E}}_{T} [GeV]",500,500,3000)->Fill(ht_chs + reader_event->met.pt() ,weight );
        plotter.getOrMake1DPre(prefix,"hbbMass",";#it{m}_{H#rightarrowbb} [GeV]",210,0,210)->Fill(hbbMass,weight );
        plotter.getOrMake1DPre(prefix,"hhMass",";#it{m}_{HH} [GeV]",168,800,5000)->Fill(hh.mass(),weight );

        if(hh.mass() < 900){
            plotter.getOrMake1DPre(prefix + "_hh800to900","hbbMass",";#it{m}_{H#rightarrowbb} [GeV]",210,0,210)->Fill(hbbMass,weight );
        }else if(hh.mass() < 1000){
            plotter.getOrMake1DPre(prefix + "_hh900to1000","hbbMass",";#it{m}_{H#rightarrowbb} [GeV]",210,0,210)->Fill(hbbMass,weight );
        }

        if(hh.mass() < 1000){
            plotter.getOrMake1DPre(prefix + "_hh800to1000","hbbMass",";#it{m}_{H#rightarrowbb} [GeV]",210,0,210)->Fill(hbbMass,weight );
        } else if (hh.mass() < 1200){
            plotter.getOrMake1DPre(prefix + "_hh1000to1200","hbbMass",";#it{m}_{H#rightarrowbb} [GeV]",210,0,210)->Fill(hbbMass,weight );
        } else if (hh.mass() < 1500){
            plotter.getOrMake1DPre(prefix + "_hh1200to1500","hbbMass",";#it{m}_{H#rightarrowbb} [GeV]",210,0,210)->Fill(hbbMass,weight );
        } else if (hh.mass() < 2000){
            plotter.getOrMake1DPre(prefix + "_hh1500to2000","hbbMass",";#it{m}_{H#rightarrowbb} [GeV]",210,0,210)->Fill(hbbMass,weight );
        } else
            plotter.getOrMake1DPre(prefix + "_hhgt2000","hbbMass",";#it{m}_{H#rightarrowbb} [GeV]",210,0,210)->Fill(hbbMass,weight );

//        auto btagSJs = PhysicsUtilities::selObjsMom(wjjCand->subJets(),
//                fjProc->param.sj_minBTagPT, fjProc->param.sj_maxBTagETA < 0 ? 999.0 : fjProc->param.sj_maxBTagETA);
//
//        std::sort(btagSJs.begin(),btagSJs.end(), [](const SubJet* a,const SubJet* b) {return a->csv() > b->csv();} );
//
//        plotter.getOrMake1DPre(prefix,"wjjCSV",";wjj CSV [GeV]",100,0,1)->Fill(btagSJs.size() ? btagSJs.front()->csv() : 0,weight );


//        plotter.getOrMake2DPre(prefix,"njets_hbbMass"          ,";H(bb) mass [GeV] ;number of jets"                                      ,100,0,500,20,-0.5,19.5)-> Fill(hbbMass, jets.size(),weight );
//        plotter.getOrMake2DPre(prefix,"njetswlep_hbbMass"      ,";H(bb) mass [GeV] ;number of jets (with leptons)"                       ,100,0,500,20,-0.5,19.5)-> Fill(hbbMass, jets_wlep.size(),weight );
//        plotter.getOrMake2DPre(prefix,"ht_wlep_hbbMass"        ,";H(bb) mass [GeV] ;#it{H}_{T} (with leptons) [GeV]"                     ,100,0,500,500,500,3000)->Fill(hbbMass, ht_wlep,weight );
//        plotter.getOrMake2DPre(prefix,"ht_wlep_plusmet_hbbMass",";H(bb) mass [GeV] ;#it{H}_{T} (with leptons) + #slash{#it{E}}_{T} [GeV]",100,0,500,500,500,3000)->Fill(hbbMass, ht_wlep + reader_event->met.pt(),weight );
    }

    void SRPlots(TString prefix, float hbbMass, bool tight) {

            auto plt = [&](const TString& pre){
                plotter.getOrMake1DPre(pre, "hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",200,0,200)->Fill(hbbMass,weight);
                if(hh.mass() >= 700 && hh.mass() < 900)
                    plotter.getOrMake1DPre(pre, "hh700to900_hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",200,0,200)->Fill(hbbMass,weight);
                if(hh.mass() >= 900 && hh.mass() < 1100)
                    plotter.getOrMake1DPre(pre, "hh900to1100_hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",200,0,200)->Fill(hbbMass,weight);
                if(hh.mass() >= 1400 && hh.mass() < 1800)
                    plotter.getOrMake1DPre(pre, "hh1400to1800_hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",200,0,200)->Fill(hbbMass,weight);
                if(hh.mass() >= 2500 && hh.mass() < 3500)
                    plotter.getOrMake1DPre(pre, "hh2500to3500_hbb_mass" ,";H(bb) mass [GeV]; arbitrary units",200,0,200)->Fill(hbbMass,weight);
            };

            if(!tight){ plt(prefix + "_hbbL");}
            if(tight){ plt(prefix + "_hbbT");}
            plt(prefix + "_hbbI");
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passTriggerPreselection) return false;
        if(!passEventFilters) return false;

        if(reader_event->process >= FillerConstants::SINGLET && reader_event->process <= FillerConstants::TTX )
            smpName = "other";

        //select for ttbar control region
        if(selectedLeptons.size() != 1) return false;
        if(!hbbCand) return false;
        if(!wjjCand) return false;
        if(hbbNSJs < 2 || wjjNSJs < 2) return false;
        if(nMedBTags_HbbV != 1) return false;
        if(hh.mass() < 800) return false;
        if(hbbCSVCat > BTagging::CSVSJ_MF) return false;

        float wCorrMass = fromWCorr->getCorrSDMass(hbbCand);
        float corrMass = hbbFJSFProc->getCorrSDMass(hbbCand);

        float smearedHBBCorr =  corrMass;
        float scaledUpHBBCorr =  corrMass;
        float scaledDownHBBCorr =  corrMass;
        if(!isRealData()){
            smearedHBBCorr *= (1 + randGen->Gaus(0,0.12)*0.6633);
            scaledUpHBBCorr   *= (1 +0.0094);
            scaledDownHBBCorr *= (1 -0.0094);
        }

        auto doSet = [&](const TString& prefix ){
            makePlots(prefix + "_raw"  ,hbbMass,weight);
            makePlots(prefix + "_wCorr",wCorrMass,weight);
            makePlots(prefix + "_corr" ,corrMass,weight);
            makePlots(prefix + "_scorr" ,smearedHBBCorr,weight);
            makePlots(prefix + "_sUcorr" ,scaledUpHBBCorr,weight);
            makePlots(prefix + "_sDcorr" ,scaledDownHBBCorr,weight);


        };

        std::string purStr = hbbCSVCat >= BTagging::CSVSJ_MF ? "L" : "AL";

            doSet(smpName+"_1l_" + purStr);
            if(selectedLepton->isMuon()) doSet(smpName+"_1m_" + purStr);
            if(!selectedLepton->isMuon()) doSet(smpName+"_1e_" + purStr);

            if(!isRealData() && reader_event->process != FillerConstants::SIGNAL){
            doSet(std::string("bkg_1l_") + purStr);
            if(selectedLepton->isMuon()) doSet(std::string("bkg_1m_") + purStr);
            if(!selectedLepton->isMuon()) doSet(std::string("bkg_1e_") + purStr);
            }



//        if(!isRealData() && passWjjSel && passHbbSel ){
//            if(reader_event->process == FillerConstants::ZJETS) smpName = "other";
//            bool tightHBB = passHbbTSel;
//
//            auto doSRPlots = [&](const TString& prefix ){
//                SRPlots(prefix + "_raw"  ,hbbMass,tightHBB);
//                SRPlots(prefix + "_wCorr",wCorrMass,tightHBB);
//                SRPlots(prefix + "_corr" ,corrMass,tightHBB);
//            };
//
//            doSRPlots(smpName+"_1l");
//            if(selectedLepton->isMuon()) doSRPlots(smpName+"_1m");
//            if(!selectedLepton->isMuon()) doSRPlots(smpName+"_1e");
//        }

        return true;


    }

    void write(TString fileName){ plotter.write(fileName);}


    HistGetter plotter;
    std::unique_ptr<SoftDropMassScaleFactorsFromW>    fromWCorr ;

};

#endif

void testSDMassCorrectionsInData(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void testSDMassCorrectionsInData(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
