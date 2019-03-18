
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "../framework/Processors/Corrections/interface/BTagScaleFactors.h"
#include "DataFormats/interface/FatJet.h"




using namespace TAna;
using namespace std;
using namespace TAna::CorrHelp;
using namespace TAna::BTagging;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) :
        DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    }

    void assignVals(const BTagging::FLAVOR flv,const float pt, const float eta,const float csv, const CorrHelp::CORRTYPE corrT,
            float& lE,float& hE,float& lSF,float& hSF) const {
        auto gE = [&](const BTagging::BTAGWP wp)->float{return sjbtagSFProc->getJetEff(flv,pt,eta,wp);};
        auto gS = [&](const BTagging::BTAGWP wp)->float{return sjbtagSFProc->getJetSF(flv,pt,eta,wp,corrT);};

        if(csv <  parameters.jets.sjBtagCorrWP[BTAG_L]){
            lE = 1.0;
            hE = gE(BTAG_L);
            lSF = 1.0;
            hSF = gS(BTAG_L);
        } else if(csv <  parameters.jets.sjBtagCorrWP[BTAG_M]){
            lE = gE(BTAG_L);
            hE = gE(BTAG_M);
            lSF = gS(BTAG_L);
            hSF = gS(BTAG_M);
        } else {
            lE = gE(BTAG_M);
            hE = 0;
            lSF = gS(BTAG_M);
            hSF = 0;

        }
    }

    float getJetCorr(const BaseRecoJet* jet,  CorrHelp::CORRTYPE lightT,  CorrHelp::CORRTYPE heavyT, float& lE, float& hE, float& lSF, float& hSF) const {
        const float pt = jet->pt();
        const float eta = jet->eta();
        const auto flv = jetFlavor(*jet);
        const float csv = jet->csv();

        lE = 1;
        hE = 1;
        lSF = 1;
        hSF = 1;

        const CorrHelp::CORRTYPE corrT =  flv == FLV_L ? lightT : heavyT;
        if(corrT == NONE) return 1.0;



        assignVals(flv,pt,eta,csv,corrT,lE,hE,lSF,hSF);

        const float origEff = lE - hE;
        const float newEff = lSF*lE - hSF*hE;
        if(origEff == 0 || newEff == 0 || origEff > 1.0 || newEff > 1.0){
            TString errStr = TString::Format(
                    "BTagScaleFactors::getJetCorr() -> Bad efficiencies: corr(pt,eta,flv):lowEff,highEff,lowSF,highSF,origEff,newEff -> %u(%f,%f,%u)%f,%f,%f,%f,%f,%f",
                    corrT,pt,eta,flv,lE,hE,lSF,hSF,origEff,newEff);
            throw std::invalid_argument(errStr.Data());
        }

        return newEff/origEff;
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        cout <<"---------------------------------------------------"<<endl;

        cout <<std::endl;
        cout <<"Top decays: ";
        for(const auto& d : smDecayEvt.topDecays ){
            cout << "["<< d.type; if(d.b) cout << *d.b; cout <<"] ";
        }
        cout <<std::endl;

        std::vector<const FatJet*> fjs;
        if(hbbCand) fjs.push_back(hbbCand);
        if(wjjCand) fjs.push_back(wjjCand);
        std::vector<const BaseRecoJet*> subjets;
        for(const auto* fj : fjs){
            if(fj){
                for(const auto& sj : fj->subJets()){
                    if(sj.absEta() < 2.4 && sj.pt() >= 30) subjets.push_back(static_cast<const BaseRecoJet*>(&sj) );
                }
            }
        }


        float lE = 1;
        float hE = 1;
        float lSF = 1;
        float hSF = 1;
        float sf = 1;

        for(const auto* j : subjets){
            std::cout <<*j <<" ->"<< j->hadronFlv()<<","<< j->csv() <<"-> "<< sjbtagSFProc->getJetCorr(j);

            auto dS = [&](CorrHelp::CORRTYPE lightT,  CorrHelp::CORRTYPE heavyT){
                sf = getJetCorr(j,lightT,heavyT,lE,hE,lSF,hSF);
                std::cout <<"["<<sf<<","<<lE<<","<<hE<<","<<lSF<<","<<hSF<<"] ";
            };
            dS(DOWN   ,NONE);
            dS(NOMINAL,NONE);
            dS(UP     ,NONE);
            dS(NONE,DOWN   );
            dS(NONE,NOMINAL);
            dS(NONE,UP     );
            std::cout <<endl;
        }
        cout <<"Weights: " <<sjbtagSFProc->getSF(parameters.jets,{hbbCand,wjjCand});
        return true;
    }




    void write(TString fileName){ plotter.write(fileName);}

    HistGetter plotter;

};

#endif


void testSJBTagSFs(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
