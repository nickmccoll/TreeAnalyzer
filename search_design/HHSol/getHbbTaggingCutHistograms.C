#include <AnalysisSupport/Utilities/interface/Types.h>
#include <Math/GenVector/LorentzVector.h>
#include <TH1.h>
#include <TMath.h>
#include "DataFormats/interface/FatJet.h"
#include "DataFormats/interface/GenParticle.h"
#include "DataFormats/interface/Lepton.h"
#include "DataFormats/interface/Momentum.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Configuration/interface/ReaderConstants.h"
#include <TString.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include "TreeAnalyzer/framework/Processors/Variables/interface/HiggsSolver.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#if !defined(__CINT__) || defined(__MAKECINT__)


#include "HHSolTreeAnalyzer.h"



using namespace TAna;
using namespace FillerConstants;


class Analyzer : public HHSolTreeAnalyzer {
public:

    Analyzer(std::string fileName) :
        HHSolTreeAnalyzer(fileName,"treeMaker/Events",1, 1),
        HSolver("") {

        setParameters(fileName);

        HSolver.setParamters(parameters.hww);
    }

    void setParameters(const std::string& fileName){
        if(fileName.find("2016") !=  std::string::npos){
            parameters = ReaderConstants::set2016Parameters();
            hpTauCut = 0.55;
        } else if(fileName.find("2017") !=  std::string::npos){
            parameters = ReaderConstants::set2017Parameters();
        } else if(fileName.find("2018") !=  std::string::npos){
            parameters = ReaderConstants::set2018Parameters();
        } else {
            parameters.hww.liFileName = "ERROR";
        }
    }

    bool runEvent() override {
        if(!HHSolTreeAnalyzer::runEvent()) return false;

        if(!passPreHHMassCuts()) return false;
        fillEventVariables();
        if(!passPostHHMassCuts()) return false;

        makePlots();

        return true;
    }

    bool passPreHHMassCuts(){
        if((*bbJet_SDmass < 100 || *bbJet_SDmass>150)) return false;
        if(*nAK4Btags != 0) return false;
        if(*qqJet_t2ot1 >= 0.75) return false;
        return true;
    }

    void fillEventVariables() {
        HSolverLiInfo hwwInfo;
        hWWLikeli = HSolver.minimize(lepton,met,qqJet,*qqJet_SDmass,hwwInfo);

        hhMass = (hwwInfo.hWW + bbJet.p4()).mass();
        ptom = hwwInfo.hWW.pt() / hhMass;

        passHPCut         = *qqJet_t2ot1 < hpTauCut && hWWLikeli < 2.5;
        passLPCut         = !passHPCut;
    }

    bool passPostHHMassCuts(){
        if(ptom < 0.3) return false;
        if(hhMass< 700) return false;

        if(hWWLikeli >= 11) return false;
        return true;

    }

    void makePlots(){
        if(isSignal()) {
            addPurityCuts(smpName);

        } else {
            if(*process <= FillerConstants::WJETS)
                addPurityCuts("bkg");
            addPurityCuts("bkgWQCD");
            addPurityCuts(smpName);
        }
    }

    void addPurityCuts(const TString& prefix){
        addLeptonCuts(prefix+"_LHP");
        if(passHPCut)
            addLeptonCuts(prefix+"_HP");
        if(passLPCut)
            addLeptonCuts(prefix+"_LP");
    }

    void addLeptonCuts(const TString& prefix){
        addMassCuts(prefix+"_emu");
        if(*isMuon == 1)
            addMassCuts(prefix+"_mu");
        else
            addMassCuts(prefix+"_e");
    }

    void addMassCuts(const TString& prefix) {
        makeTaggingPlots(prefix+"_mIncl");

        for(auto targetMass : targetMasses) {
            if(isInMassWindow(targetMass)) {
                TString targetPrefix = prefix + "_m" + ASTypes::int2Str(targetMass);
                makeTaggingPlots(targetPrefix);
            }
        }
    }

    void makeTaggingPlots(const TString& prefix){
        plotter.getOrMake1DPre(prefix, "ak8Tag","; Deep AK8 #it{D}",501,0,1.002)
                ->Fill(sanitizedDeepAK8Tag,weight);
        plotter.getOrMake1DPre(prefix, "dCSVCat","; Deep CSV cat.",11,-0.5,10.5)
                ->Fill(*hbbCat,weight);
    }

    bool isInMassWindow(double targetResonanceMass){
        if(isSignal()) return true;
        return std::abs(hhMass - targetResonanceMass) < targetResonanceMass*0.10;
    }


    void write(TString fileName) {
        plotter.write(fileName);
    }


    HistGetter plotter;
    double hpTauCut = 0.45;
    std::vector<int> targetMasses = {800,1000,1200,1400,2000,2500,3000};

    bool passHPCut = false;
    bool passLPCut = false;
    double hhMass = 0;
    double ptom = 0;
    double hWWLikeli = 0;

    HSolverLi HSolver;
    ParameterSet parameters;
};

#endif

void getHbbTaggingCutHistograms(std::string fileName, std::string outFileName,
        float xSec=-1, float numEvent=-1){

    Analyzer a(fileName);
    std::string outN = "hSolTrees_hbbTaggingCutHistograms";


    outN += "_"+ outFileName+".root";

    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outN);

}
