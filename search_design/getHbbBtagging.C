
#if !defined(__CINT__) || defined(__MAKECINT__)


#include "TreeReaders/interface/EventReader.h"
#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "DataFormats/interface/Lepton.h"
#include "DataFormats/interface/FatJet.h"

using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) :
        DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {


//        turnOffCorr(CORR_LEP );
//        turnOffCorr(CORR_SJBTAG);
//        turnOffCorr(CORR_AK4BTAG);
//        turnOffCorr(CORR_JER);



    }

    void setupParameters() override {
        DefaultSearchRegionAnalyzer::setupParameters();
        parameters.fatJets.hbb_minSJs=0;
        parameters.event.doTTBarStitching = false;
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(lepChan != SINGLELEP) return false;
        if(selectedLeptons.size()!=1) return false;
        if(!hbbCand || !wjjCand) return false;

        if(hh.mass() < 700) return false;
        if(hbbMass < 30 || hbbMass > 210 ) return false;
        if(nMedBTags_HbbV) return false;

//        if(wjjCand->tau2otau1() > 0.75) return false;
//        if(hwwLi>= 11) return false;

        const bool HP = (wjjCand->tau2otau1() < 0.55 && hwwLi < 2.5);
        const bool LP = !HP;



        auto mkBTAGPlts = [&](const TString& prefix){

            plotter.getOrMake1DPre(prefix, "bbt",";deep double-b",500,0,1)
                    ->Fill(hbbCand->bbt(),weight);

            plotter.getOrMake1DPre(prefix, "mdZHbb",";deep mass dec. Z/H#rightarrow b#bar{b}"
                    ,500,0,1)->Fill(hbbCand->deep_MDZHbb(),weight);

            plotter.getOrMake1DPre(prefix, "mdHbb",";deep mass dec. H#rightarrow b#bar{b}",500,0,1)
                    ->Fill(hbbCand->deep_MDHbb(),weight);

            plotter.getOrMake1DPre(prefix, "Hbb",";deep H#rightarrow b#bar{b}",500,0,1)
                                ->Fill(hbbCand->deep_Hbb(),weight);

            std::vector<float> sjBs(2,-1);
            if(hbbCand->nSubJets()) sjBs[0] = hbbCand->subJet(0).deep_csv();
            if(hbbCand->nSubJets()>1) sjBs[1] = hbbCand->subJet(1).deep_csv();
            std::sort(sjBs.begin(), sjBs.end(),[](float a, float b){return a>b;} );

            bool leadM = sjBs[0] >= parameters.jets.DeepCSV_WP [BTagging::BTAG_M];
            bool subM = sjBs[1] >= parameters.jets.DeepCSV_WP [BTagging::BTAG_M];
            bool subLNM = !subM && sjBs[1] >=
                    parameters.jets.DeepCSV_WP [BTagging::BTAG_L];

            plotter.getOrMake1DPre(prefix, "sjMaxDeepCSV",";max subjet deep CSV",500,0,1)
                    ->Fill(sjBs[0],weight);

            if(subM||subLNM)
                plotter.getOrMake1DPre(prefix, "subLM_sjMaxDeepCSV",";max subjet deep CSV",500,0,1)
                        ->Fill(sjBs[0],weight);

            plotter.getOrMake1DPre(prefix, "sjMinDeepCSV",";min subjet deep CSV",500,0,1)
                    ->Fill(sjBs[1],weight);
            if(leadM)
                plotter.getOrMake1DPre(prefix, "leadM_sjMinDeepCSV",";min subjet deep CSV",500,0,1)
                        ->Fill(sjBs[1],weight);

            double catFlt = 0;
            if(leadM && subM) catFlt = 1;
            else if(leadM && subLNM) catFlt = 0.75;
            else if(leadM) catFlt = 0.5;

            plotter.getOrMake1DPre(prefix, "sJBTagging",";subjet deep CSV",500,0,1)
                    ->Fill(catFlt,weight);


        };



        auto mkHHPlts = [&](const TString& prefix){
            mkBTAGPlts(prefix+"_mHHInlc");
            const auto hhM = hh.mass();
            auto doHH = [&](const int mC){
                if(!isSignal()){
                    if(std::abs(hh.mass()-mC) < mC*0.10)
                        mkBTAGPlts(prefix+"_m"+ ASTypes::int2Str(mC));
                } else {
                    if(mC == signal_mass)
                        mkBTAGPlts(prefix+"_m"+ ASTypes::int2Str(mC));
                }
            };
            doHH(800);
            doHH(1000);
            doHH(2000);
            doHH(2500);
            doHH(3000);
        };

        auto mkHbbPlots = [&](const TString& prefix){
            mkHHPlts(prefix+"_HbbIncl");
            int nSJs = PhysicsUtilities::selObjsMom(hbbCand->subJets(),
                    parameters.fatJets.sj_minPT, parameters.fatJets.sj_maxETA < 0 ? 999.0
                            : parameters.fatJets.sj_maxETA).size();
            if(nSJs == 2) mkHHPlts(prefix+"_Hbb2SJ");
            if(hbbMass > 100 && hbbMass < 150)
                mkHHPlts(prefix+"_HbbTMass");
            if(hbbMass > 100 && hbbMass < 150 && nSJs==2)
                mkHHPlts(prefix+"_Hbb2SJTMass");
        };


        auto mkHADPlts = [&](const TString& prefix){
            mkHbbPlots(prefix+"_LHP");
            if(LP) mkHbbPlots(prefix+"_LP");
            if(HP) mkHbbPlots(prefix+"_HP");
        };

        auto mkCutLEPPlts = [&](const TString& prefix){
            mkHADPlts(prefix+"_emu");
            if(selectedLepton->isMuon() == 1)
                mkHADPlts(prefix+"_mu");
            else
                mkHADPlts(prefix+"_e");
        };


        if(isSignal()){
            mkCutLEPPlts(smpName);
        } else {
            if(mcProc <= FillerConstants::WJETS)
                mkCutLEPPlts("bkg");
            mkCutLEPPlts("bkgWQCD");
            mkCutLEPPlts(smpName);
        }

        return true;

    }


    void write(TString fileName){ plotter.write(fileName);}


    HistGetter plotter;



};

#endif

void getHbbBtagging(std::string fileName, int treeInt, int randSeed, std::string outFileName,
        float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
