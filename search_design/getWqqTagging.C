
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


        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
        turnOffCorr(CORR_JER);



    }

    void setupParameters() override {
        DefaultSearchRegionAnalyzer::setupParameters();
        parameters.fatJets.wjj_minSJs = 0;
        parameters.event.doTTBarStitching = false;
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(lepChan != SINGLELEP) return false;
        if(selectedLeptons.size()!=1) return false;
        if(!hbbCand || !wjjCand) return false;

        if(hh.mass() < 700) return false;
        if(hbbMass < 30 || hbbMass > 210 ) return false;
        if(hwwLi>= 11) return false;
        if(nMedBTags_HbbV) return false;
        if(hbbCSVCat < 4) return false;

        const bool HP = hwwLi < 2.5;
        const bool LP = !HP;


        auto mkWqqPlts = [&](const TString& prefix){

            const double sdMass = wjjCand->sdMom().mass();
            const double sdPT = wjjCand->sdMom().pt();
            const double t2oT1 = wjjCand->tau2otau1();

            plotter.getOrMake1DPre(prefix, "wqq",";deep wqq",500,0,1)
                    ->Fill(wjjCand->deep_W(),weight);

            plotter.getOrMake1DPre(prefix, "tau2otau1",";#tau_2 / #tau_1",500,0,1)
                    ->Fill(t2oT1,weight);



            const double tau_21DDT = t2oT1 + 0.080 * std::log( sdMass*sdMass/(sdPT) );

            plotter.getOrMake1DPre(prefix, "tau21ddt",";#DDT tau_2 / #tau_1"
                    ,500,0,1)->Fill(tau_21DDT,weight);


            double oldBinning = 0;
            if (t2oT1 < 0.55)  oldBinning = 1;
            else if (t2oT1 < 0.75)  oldBinning = .5;

            plotter.getOrMake1DPre(prefix, "oldTau21Bins",";#tau_2 / #tau_1",500,0,1)
                    ->Fill(oldBinning,weight);

            double newBinning = 0;
            if (t2oT1 < 0.45)  newBinning = 1;
            else if (t2oT1 < 0.75)  newBinning = .5;

            plotter.getOrMake1DPre(prefix, "newTau21Bins",";#tau_2 / #tau_1",500,0,1)
                    ->Fill(newBinning,weight);

            double ddtBinning = 0;
            if (tau_21DDT < 0.43)  ddtBinning = 1;
            else if (tau_21DDT < 0.79)  ddtBinning = .5;

            plotter.getOrMake1DPre(prefix, "DDTTau21Bins",";#tau_2 / #tau_1",500,0,1)
                    ->Fill(ddtBinning,weight);

        };



        auto mkHHPlts = [&](const TString& prefix){
            mkWqqPlts(prefix+"_mHHInlc");
            const auto hhM = hh.mass();
            auto doHH = [&](const int mC){
                if(!isSignal()){
                    if(std::abs(hh.mass()-mC) < mC*0.10)
                        mkWqqPlts(prefix+"_m"+ ASTypes::int2Str(mC));
                } else {
                    if(mC == signal_mass)
                        mkWqqPlts(prefix+"_m"+ ASTypes::int2Str(mC));
                }
            };
            doHH(800);
            doHH(1000);
            doHH(2000);
            doHH(2500);
            doHH(3000);
        };



        auto mkHADPlts = [&](const TString& prefix){
            mkHHPlts(prefix+"_LHP");
            if(LP) mkHHPlts(prefix+"_LP");
            if(HP) mkHHPlts(prefix+"_HP");

            int nSJs = PhysicsUtilities::selObjsMom(wjjCand->subJets(),
                    parameters.fatJets.sj_minPT, parameters.fatJets.sj_maxETA < 0 ? 999.0
                            : parameters.fatJets.sj_maxETA).size();


            if(nSJs == 2){
                mkHHPlts(prefix+"_LHP_2SJ");
                if(LP) mkHHPlts(prefix+"_LP_2SJ");
                if(HP) mkHHPlts(prefix+"_HP_2SJ");
            }
        };

        auto mkHbbPlots = [&](const TString& prefix){

            mkHADPlts(prefix+"_bLMT");
            if(hbbCSVCat == 4) mkHADPlts(prefix+"_bL");
            if(hbbCSVCat == 5) mkHADPlts(prefix+"_bM");
            if(hbbCSVCat == 6) mkHADPlts(prefix+"_bT");

        };

        auto mkCutLEPPlts = [&](const TString& prefix){
            mkHbbPlots(prefix+"_emu");
            if(selectedLepton->isMuon() == 1)
                mkHbbPlots(prefix+"_mu");
            else
                mkHbbPlots(prefix+"_e");
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

void getWqqTagging(std::string fileName, int treeInt, int randSeed, std::string outFileName,
        float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
