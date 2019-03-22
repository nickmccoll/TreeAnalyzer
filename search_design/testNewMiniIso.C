
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"

using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
    : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

        turnOffCorr(CORR_TRIG);
        turnOffCorr(CORR_PU  );
        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
        turnOffCorr(CORR_SDMASS);
        turnOffCorr(CORR_JER);
    }



    void mkPlots(const std::string& prefix, const Electron* electron){
        bool highEta = electron->absEta() > 1.5;
        const std::string etaString = electron->absEta() > 1.5 ? "eta_gt1p5" : "eta_lt1p5";
        const std::string htString = ht_chs  > 2000 ? "ht_gt2000" :
                ( ht_chs > 1000 ? "ht_1000to2000" : "ht_500to1000" );
        const std::string ptString = electron->pt()  > 200 ? "pt_gt200" :
                ( electron->pt() > 100 ? "pt_100to200" :
                        ( electron->pt() > 50 ? "pt_50to100" : "pt_2to50" ) );

        plotter.getOrMake1DPre((prefix + "_"+ "eta_incl").c_str(),"ht",";#it{H}_{T} [GeV]",100,500,5500)->Fill(ht_chs,weight);
        plotter.getOrMake1DPre((prefix + "_"+ etaString).c_str(),"ht",";#it{H}_{T} [GeV]",100,500,5500)->Fill(ht_chs,weight);

        plotter.getOrMake1DPre((prefix + "_"+ "eta_incl").c_str(),"pt",";lepton #it{p}_{T} [GeV]",100,0,500)->Fill(electron->pt(),weight);
        plotter.getOrMake1DPre((prefix + "_"+ etaString).c_str(),"pt",";lepton #it{p}_{T} [GeV]",100,0,500)->Fill(electron->pt(),weight);

        plotter.getOrMake1DPre((prefix + "_"+ "ht_incl").c_str(),"eta",";lepton |#eta|",50,0,2.5)->Fill(electron->absEta(),weight);
        plotter.getOrMake1DPre((prefix + "_"+ htString).c_str(),"eta",";lepton |#eta|",50,0,2.5)->Fill(electron->absEta(),weight);
        plotter.getOrMake1DPre((prefix + "_"+ "pt_incl").c_str(),"eta",";lepton |#eta|",50,0,2.5)->Fill(electron->absEta(),weight);
        plotter.getOrMake1DPre((prefix + "_"+ ptString).c_str(),"eta",";lepton |#eta|",50,0,2.5)->Fill(electron->absEta(),weight);

        plotter.getOrMake1DPre((prefix).c_str(),"count",";count",1,-0.5,0.5)->Fill(0.0,weight);

    }

    void testISOVal(const std::string& prefix, const LeptonParameters& noValParameters,  const Electron* electron=0){
        LeptonParameters val0p1 = noValParameters;
        val0p1.el_maxISO =0.1;
        LeptonParameters val0p2 = noValParameters;
        val0p2.el_maxISO =0.2;

        if(electron == 0){
            auto val0p1Iso_electrons = LeptonProcessor::getElectrons(val0p1,*reader_electron);
            auto val0p2Iso_electrons = LeptonProcessor::getElectrons(val0p2,*reader_electron);

            if(val0p1Iso_electrons.size())mkPlots(prefix+"_0p1",val0p1Iso_electrons.front());
            if(val0p2Iso_electrons.size())mkPlots(prefix+"_0p2",val0p2Iso_electrons.front());
        } else {
            if(LeptonProcessor::isGoodElectron(val0p1,electron))
                mkPlots(prefix+"_0p1",electron);
            if(LeptonProcessor::isGoodElectron(val0p2,electron))
                mkPlots(prefix+"_0p2",electron);
        }
    }

    void testISOType(const std::string& prefix, const LeptonParameters& noISOParameters,  const Electron* electron=0){
        LeptonParameters oldIso = noISOParameters;
        oldIso.el_getISO =&Electron::miniIso;
        LeptonParameters newIso = noISOParameters;
        newIso.el_getISO =&Electron::miniIsoFP;

        auto noIso_electrons = LeptonProcessor::getElectrons(noISOParameters,*reader_electron);
        mkPlots(prefix+"_noIso",noIso_electrons.front());
        testISOVal(prefix+"_oldISO",oldIso,electron);
        testISOVal(prefix+"_newISO",newIso,electron);
    }

    void testID(const std::string& prefix, const LeptonParameters& inclParameters,  const Electron* electron=0){
        LeptonParameters medIDNoIso = inclParameters;
        medIDNoIso.el_getID =&Electron::passMedID_noIso;
        LeptonParameters tightIDNoIso = inclParameters;
        tightIDNoIso.el_getID =&Electron::passTightID_noIso;


        if(electron == 0){
            auto medNoIso_electrons = LeptonProcessor::getElectrons(medIDNoIso,*reader_electron);
            auto tightNoIso_electrons = LeptonProcessor::getElectrons(tightIDNoIso,*reader_electron);

            if(medNoIso_electrons.size())testISOType(prefix+"_medID",medIDNoIso);
            if(tightNoIso_electrons.size())testISOType(prefix+"_tightID",tightIDNoIso);
        } else {
            if(LeptonProcessor::isGoodElectron(medIDNoIso,electron))
                testISOType(prefix+"_medID",medIDNoIso,electron);
            if(LeptonProcessor::isGoodElectron(tightIDNoIso,electron))
                testISOType(prefix+"_tightID",tightIDNoIso,electron);
        }
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(ht_chs < 500) return false;
        std::string prefix = smpName.Data();

        const Electron * recoElectron = 0;

        if(isSignal()){
            if(diHiggsEvt.type != DiHiggsEvent::E) return false;
            const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,20,2.4);
            double nearestDR =10;
            int idx = PhysicsUtilities::findNearestDRDeref(*diHiggsEvt.w1_d1,electrons,nearestDR,0.2);
            if(idx < 0) return false;
            recoElectron = electrons[idx];
            mkPlots(prefix+"_reco",recoElectron);
        }

        LeptonParameters paramsNoIDNoISO = parameters.leptons;
        paramsNoIDNoISO.el_maxISO = -1;
        paramsNoIDNoISO.el_getID = &Electron::passInclID;
        paramsNoIDNoISO.el_getISO = &Electron::inclIso;

        testID(prefix,paramsNoIDNoISO,recoElectron);
        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void testNewMiniIso(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
