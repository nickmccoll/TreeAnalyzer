
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeCopier.h"
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
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Corrections/interface/TriggerScaleFactors.h"
#include "Processors/Corrections/interface/BTagScaleFactors.h"
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "Processors/Corrections/interface/JetAndMETCorrections.h"

#include "TPRegexp.h"
using namespace TAna;
using namespace CorrHelp;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {
        addUncVariables = (treeType == TREE_OTHER);
    }

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed, CORRTYPE jerUNC, CORRTYPE jesUNC,CORRTYPE metUNC) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
        JERAK4PuppiProc ->setCorrType(jerUNC);
        JERAK4CHSProc   ->setCorrType(jerUNC);
        JERAK8PuppiProc ->setCorrType(jerUNC);
        if(jesUNC == UP || jesUNC == DOWN ){
            JESUncProc ->setCorrType(jesUNC);
            turnOnCorr(CORR_JES);
        }
        if(metUNC == UP || metUNC == DOWN){
            METUncProc ->setCorrType(metUNC);
            turnOnCorr(CORR_MET);
        }

    }

    virtual BaseEventAnalyzer * setupEventAnalyzer() override {return new CopierEventAnalyzer();}

    virtual void bookOutputVariables() override {

        if(isRealData()){
            i_dataset             =  outTree->add<size8>  ("","dataset"  ,"b",0);
            i_dataRun             =  outTree->add<size8>  ("","dataRun"  ,"b",0);
        } else {
            i_process     =outTree->add<size8>  ("","process","b",0);
            i_dhType      =outTree->add<size8>  ("","dhType" ,"b",0);
            i_xsec        =outTree->add<float>  ("","xsec"   ,"F",0);
            i_trig_N      =outTree->add<float>  ("","trig_N" ,"F",0);
            i_pu_N        =outTree->add<float>  ("","pu_N"   ,"F",0);
            i_lep_N       =outTree->add<float>  ("","lep_N"  ,"F",0);
            i_btag_N      =outTree->add<float>  ("","btag_N" ,"F",0);
        }
        i_ht          =outTree->add<float>  ("","ht"        ,"F",0);
        i_met         =outTree->add<float>  ("","met"       ,"F",0);
        i_isMuon      =outTree->add<size8>  ("","isMuon"    ,"b",0);
        i_lepPT       =outTree->add<float>  ("","lepPT"     ,"F",0);
        i_lepETA      =outTree->add<float>  ("","lepETA"    ,"F",0);

        i_hbbMass     =outTree->add<float>  ("","hbbMass"   ,"F",0);
        i_hbbPT       =outTree->add<float>  ("","hbbPT"     ,"F",0);
        i_hbbNSJs     =outTree->add<size8>  ("","hbbNSJs"   ,"b",0);
        i_hbbCSVCat   =outTree->add<size8>  ("","hbbCSVCat" ,"b",0);
        i_hbbBBT      =outTree->add<float>  ("","hbbBBT"    ,"F",0);
        i_hbbTau2o1   =outTree->add<float>  ("","hbbTau2o1" ,"F",0);

        i_hhMass      =outTree->add<float>  ("","hhMass"    ,"F",0);
        i_wwDM        =outTree->add<float>  ("","wwDM"      ,"F",0);
        i_hwwPT       =outTree->add<float>  ("","hwwPT"     ,"F",0);
        i_hwwETA      =outTree->add<float>  ("","hwwETA"    ,"F",0);
        i_wjjCSVCat   =outTree->add<size8>  ("","wjjCSVCat" ,"b",0);
        i_wjjTau2o1   =outTree->add<float>  ("","wjjTau2o1" ,"F",0);
        i_wjjMass     =outTree->add<float>  ("","wjjMass"   ,"F",0);
        i_wjjPT       =outTree->add<float>  ("","wjjPT"     ,"F",0);
        i_wjjNSJs     =outTree->add<size8>  ("","wjjNSJs"   ,"b",0);
        i_wlnuMass    =outTree->add<float>  ("","wlnuMass"  ,"F",0);
        i_wlnuPT      =outTree->add<float>  ("","wlnuPT"    ,"F",0);
        i_nAK4Btags   =outTree->add<size8>  ("","nAK4Btags" ,"b",0);

        if(!isRealData()){
            i_hbbWQuark   =outTree->add<size8>  ("","hbbWQuark"  ,"b",0);
            i_hbbWEQuark  =outTree->add<size8>  ("","hbbWEQuark"  ,"b",0);
        }

        if(addUncVariables){
            i_w_muIDUp        = outTree->add<float>("","w_muIDUp","F",0);
            i_w_muISOUp       = outTree->add<float>("","w_muISOUp","F",0);
            i_w_elRecoUp      = outTree->add<float>("","w_elRecoUp","F",0);
            i_w_elIDUp        = outTree->add<float>("","w_elIDUp","F",0);
            i_w_elISOUp       = outTree->add<float>("","w_elISOUp","F",0);
            i_w_b_realUp     = outTree->add<float>("","w_b_realUp","F",0);
            i_w_b_fakeUp      = outTree->add<float>("","w_b_fakeUp","F",0);
            i_w_muIDDown      = outTree->add<float>("","w_muIDDown","F",0);
            i_w_muISODown     = outTree->add<float>("","w_muISODown","F",0);
            i_w_elRecoDown    = outTree->add<float>("","w_elRecoDown","F",0);
            i_w_elIDDown      = outTree->add<float>("","w_elIDDown","F",0);
            i_w_elISODown     = outTree->add<float>("","w_elISODown","F",0);
            i_w_b_realDown   = outTree->add<float>("","w_b_realDown","F",0);
            i_w_b_fakeDown    = outTree->add<float>("","w_b_fakeDown","F",0);
            i_w_scale        = outTree->addMulti<float>("","w_scale",0);
            i_w_pdf          = outTree->addMulti<float>("","w_pdf",0);
        }

    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passTriggerPreselection) return false;
        if(!passEventFilters) return false;
        if(selectedLeptons.size() != 1) return false;
        if(!hbbCand) return false;
        if(!wjjCand) return false;


        if(isRealData()){
            outTree->fill(i_dataset     ,reader_event->dataset);
            outTree->fill(i_dataRun     ,reader_event->dataRun);
        } else {
            outTree->fill(i_process     ,reader_event->process);
            outTree->fill(i_dhType      ,size8(diHiggsEvt.type));
            outTree->fill(i_xsec        ,float( EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi())));
            outTree->fill(i_trig_N      ,float(smDecayEvt.promptElectrons.size() + smDecayEvt.promptMuons.size() ? trigSFProc->getLeptonTriggerSF(ht_chs, (selectedLepton && selectedLepton->isMuon())) : 1.0 ));
            outTree->fill(i_pu_N        ,float(puSFProc->getCorrection(reader_event->nTruePUInts,CorrHelp::NOMINAL)));
            outTree->fill(i_lep_N       ,float(leptonSFProc->getSF()));
            outTree->fill(i_btag_N      ,float(sjbtagSFProc->getSF({hbbCand})*ak4btagSFProc->getSF(jets_HbbV)));
        }

        float ljchef = -1;
        if(reader_jet_chs->jets.size()){
            ljchef = reader_jet_chs->chef->at(reader_jet_chs->jets[0].index());
        }

        outTree->fill(i_ht     ,float(ht_chs));
        outTree->fill(i_met    ,float(reader_event->met.pt()));

        outTree->fill(i_isMuon      ,size8(selectedLepton->isMuon()));
        outTree->fill(i_lepPT       ,float(selectedLepton->pt()));
        outTree->fill(i_lepETA      ,float(selectedLepton->eta()));
        outTree->fill(i_hbbMass     ,float(hbbMass));
        outTree->fill(i_hbbPT       ,float(hbbCand->pt()));
        outTree->fill(i_hbbNSJs     ,size8(hbbNSJs));
        outTree->fill(i_hbbCSVCat   ,size8(hbbCSVCat));
        outTree->fill(i_hbbBBT       ,float(hbbCand->bbt()));
        outTree->fill(i_hbbTau2o1   ,float(hbbCand->tau2otau1()));

        outTree->fill(i_hhMass      ,float(hh.mass()));
        outTree->fill(i_wlnuMass    ,float(wlnu.mass()));
        outTree->fill(i_wwDM        ,float(wwDM));
        outTree->fill(i_hwwPT       ,float(hWW.pt()));
        outTree->fill(i_hwwETA      ,float(hWW.eta()));
        outTree->fill(i_wjjCSVCat   ,size8(wjjCSVCat));
        outTree->fill(i_wjjTau2o1   ,float(wjjCand->tau2otau1()));
        outTree->fill(i_wjjMass     ,float(wjjCand->sdMom().mass()));
        outTree->fill(i_wjjPT       ,float(wjjCand->pt()));
        outTree->fill(i_wjjNSJs     ,size8(wjjNSJs));
        outTree->fill(i_wlnuPT      ,float(wlnu.pt()));
        outTree->fill(i_nAK4Btags   ,size8(std::min(nMedBTags_HbbV,250)));

        if(!isRealData()){
            const float matchR = 0.8*0.8;

            int topDecayType = 0; //NONE b wj wjb wjj wjjb
            int maxNTopQuarks = 0;
            int totNTopQuarks = 0;

            for(const auto& d : smDecayEvt.topDecays  ){
                if(d.type != TopDecay::HAD){
                    if(d.type != TopDecay::BAD){
                        if(PhysicsUtilities::deltaR2(*d.b,*hbbCand) < matchR) totNTopQuarks++;
                    }
                    continue;
                }
                bool passB = (PhysicsUtilities::deltaR2(*d.b,*hbbCand) < matchR);
                int nW =  (PhysicsUtilities::deltaR2(*d.W_decay.dau1,*hbbCand) < matchR) + (PhysicsUtilities::deltaR2(*d.W_decay.dau2,*hbbCand) < matchR);
                int nT = nW + passB;
                totNTopQuarks += nT;
                if(nT <= maxNTopQuarks) continue;
                maxNTopQuarks = nT;
                if(nT == 1){
                    if(passB) topDecayType = 1;
                    else topDecayType = 2;
                } else if(nT == 2){
                    if(passB) topDecayType = 3;
                    else topDecayType = 4;
                } else if(nT == 3) topDecayType = 5;
            }

            int maxNWQuarks = 0;
            int totNWQuarks = 0;
            int WDecayType  = 0;//NONE b wj wjb wjj wjjb
            for(const auto& d : smDecayEvt.bosonDecays  ){
                if(d.type != BosonDecay::Z_HAD && d.type != BosonDecay::W_HAD ) continue;
                int nW = (PhysicsUtilities::deltaR2(*d.dau1,*hbbCand) < matchR) +  (PhysicsUtilities::deltaR2(*d.dau2,*hbbCand) < matchR);
                totNWQuarks += nW;
                if(nW <= maxNWQuarks) continue;
                maxNWQuarks = nW;
                if(nW == 1) WDecayType = 2;
                else if(nW == 2) WDecayType = 4;
            }
            size8 decayType = 0;
            size8 nExtraQuarks = 0;
            if(maxNTopQuarks >= maxNWQuarks){
                decayType = topDecayType;
                nExtraQuarks = (totNTopQuarks - maxNTopQuarks) + totNWQuarks;
            } else {
                decayType = WDecayType;
                nExtraQuarks = (totNWQuarks - maxNWQuarks) + totNTopQuarks;
            }

            outTree->fill(i_hbbWQuark   ,decayType);
            outTree->fill(i_hbbWEQuark   ,nExtraQuarks);
        }

        if(addUncVariables){
            const float nomMu = leptonSFProc->getMuonSF();
            const float nomEl = leptonSFProc->getElectronSF();
            outTree->fill(i_w_muIDUp     ,float(leptonSFProc->getMuonSF(NONE,UP,NOMINAL)*nomEl));
            outTree->fill(i_w_muISOUp    ,float(leptonSFProc->getMuonSF(NONE,NOMINAL,UP)*nomEl));
            outTree->fill(i_w_elRecoUp   ,float(leptonSFProc->getElectronSF(UP,NOMINAL,NOMINAL)*nomMu));
            outTree->fill(i_w_elIDUp     ,float(leptonSFProc->getElectronSF(NOMINAL,UP,NOMINAL)*nomMu));
            outTree->fill(i_w_elISOUp    ,float(leptonSFProc->getElectronSF(NOMINAL,NOMINAL,UP)*nomMu));
            outTree->fill(i_w_b_realUp   ,float(ak4btagSFProc->getSF(jets_HbbV,NOMINAL,UP)* sjbtagSFProc->getSF({hbbCand},NOMINAL,UP)));
            outTree->fill(i_w_b_fakeUp   ,float(ak4btagSFProc->getSF(jets_HbbV,UP,NOMINAL)* sjbtagSFProc->getSF({hbbCand},UP,NOMINAL)));

            outTree->fill(i_w_muIDDown   ,float(leptonSFProc->getMuonSF(NONE,DOWN,NOMINAL)*nomEl));
            outTree->fill(i_w_muISODown  ,float(leptonSFProc->getMuonSF(NONE,NOMINAL,DOWN)*nomEl));
            outTree->fill(i_w_elRecoDown ,float(leptonSFProc->getElectronSF(DOWN,NOMINAL,NOMINAL)*nomMu));
            outTree->fill(i_w_elIDDown   ,float(leptonSFProc->getElectronSF(NOMINAL,DOWN,NOMINAL)*nomMu));
            outTree->fill(i_w_elISODown  ,float(leptonSFProc->getElectronSF(NOMINAL,NOMINAL,DOWN)*nomMu));
            outTree->fill(i_w_b_realDown ,float(ak4btagSFProc->getSF(jets_HbbV,NOMINAL,DOWN)* sjbtagSFProc->getSF({hbbCand},NOMINAL,DOWN)));
            outTree->fill(i_w_b_fakeDown ,float(ak4btagSFProc->getSF(jets_HbbV,DOWN,NOMINAL)* sjbtagSFProc->getSF({hbbCand},DOWN,NOMINAL)));


            for(unsigned int i = 1; i < 9; ++i){
                if(i == 5 || i ==7)continue; //told to ignore
                outTree->fillMulti(i_w_scale,(*reader_event->genWeights)[i]);
            }
            for(unsigned int i = 111; i < 211; ++i) //Number 110 is the nominal
                outTree->fillMulti(i_w_pdf,(*reader_event->genWeights)[i]);


        }




        return true;
    }



    //Event information and weights
    size i_process    = 0;
    size i_dhType     = 0;
    size i_dataset    = 0;
    size i_dataRun    = 0;
    size i_xsec       = 0;
    size i_trig_N = 0;
    size i_pu_N   = 0;
    size i_lep_N  = 0;
    size i_btag_N = 0;

    //SR variables
    size i_ht        = 0;
    size i_met       = 0;
    size i_isMuon    = 0;
    size i_lepPT     = 0;
    size i_lepETA    = 0;
    size i_hbbMass   = 0;
    size i_hbbPT     = 0;
    size i_hbbNSJs   =0;
    size i_hbbCSVCat = 0;
    size i_hbbBBT    = 0;
    size i_hbbTau2o1 = 0;

    size i_hhMass    = 0;
    size i_wwDM      = 0;
    size i_hwwPT     = 0;
    size i_hwwETA    = 0;
    size i_wjjCSVCat = 0;
    size i_wjjTau2o1 = 0;
    size i_wjjMass   = 0;
    size i_wjjPT     = 0;
    size i_wjjNSJs   =0;
    size i_wlnuMass  = 0;
    size i_wlnuPT    = 0;
    size i_nAK4Btags = 0;

    //BE extra variables
    size i_hbbWQuark   =0;
    size i_hbbWEQuark   =0;

    //systematic variables
    size i_w_muIDUp      = 0;
    size i_w_muISOUp     = 0;
    size i_w_elRecoUp    = 0;
    size i_w_elIDUp      = 0;
    size i_w_elISOUp     = 0;
    size i_w_b_realUp   = 0;
    size i_w_b_fakeUp    = 0;
    size i_w_muIDDown    = 0;
    size i_w_muISODown   = 0;
    size i_w_elRecoDown  = 0;
    size i_w_elIDDown    = 0;
    size i_w_elISODown   = 0;
    size i_w_b_realDown = 0;
    size i_w_b_fakeDown  = 0;

    size i_w_scale       = 0;
    size i_w_pdf         = 0;

    bool addUncVariables = false;



};




void doOne(const std::string& fileName, const int treeInt,int randSeed,  const std::string& outFileName, float xSec = -1, float numEvent = -1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    if(xSec > 0) a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
}

void doOneVar(const std::string& fileName, const int treeInt,int randSeed, const std::string& outFileName, const CORRTYPE jerType, CORRTYPE jesUNC, CORRTYPE metUNC, float xSec = -1, float numEvent = -1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed,jerType,jesUNC,metUNC);
    if(xSec > 0) a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
}

#endif



void makeBETrees(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    doOne(fileName,treeInt,randSeed,outFileName);
    if(treeInt==2){
        size_t lastindex = outFileName.find_last_of(".");
        std::string extLessName = outFileName.substr(0, lastindex);
        doOneVar(fileName,treeInt,randSeed+1,extLessName+"_JERUp.root"  ,UP     ,NONE,NONE);
        doOneVar(fileName,treeInt,randSeed+2,extLessName+"_JERDown.root",DOWN   ,NONE,NONE);
        doOneVar(fileName,treeInt,randSeed+3,extLessName+"_JESDOWN.root",NOMINAL,DOWN,NONE);
        doOneVar(fileName,treeInt,randSeed+4,extLessName+"_JESUp.root"  ,NOMINAL,UP  ,NONE);
        doOneVar(fileName,treeInt,randSeed+5,extLessName+"_METDOWN.root",NOMINAL,NONE,DOWN);
        doOneVar(fileName,treeInt,randSeed+6,extLessName+"_METUp.root"  ,NOMINAL,NONE,UP  );

    }
}
void makeBETrees(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    doOne(fileName,treeInt,randSeed,outFileName,xSec,numEvent);
    if(treeInt==2){
        size_t lastindex = outFileName.find_last_of(".");
        std::string extLessName = outFileName.substr(0, lastindex);
        doOneVar(fileName,treeInt,randSeed+1,extLessName+"_JERUp.root"  ,UP     ,NONE,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+2,extLessName+"_JERDown.root",DOWN   ,NONE,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+3,extLessName+"_JESDOWN.root",NOMINAL,DOWN,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+4,extLessName+"_JESUp.root"  ,NOMINAL,UP  ,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+5,extLessName+"_METDOWN.root",NOMINAL,NONE,DOWN,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+6,extLessName+"_METUp.root"  ,NOMINAL,NONE,UP  ,xSec,numEvent);
    }

}
