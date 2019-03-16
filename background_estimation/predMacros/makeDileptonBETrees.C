
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "TreeAnalyzer/interface/DileptonSearchRegionAnalyzer.h"
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
#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Variables/interface/DileptonSelection.h"
#include "Processors/Corrections/interface/TriggerScaleFactors.h"
#include "Processors/Corrections/interface/BTagScaleFactors.h"
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "Processors/Corrections/interface/JetAndMETCorrections.h"

#include "TPRegexp.h"
using namespace TAna;
using namespace CorrHelp;
class Analyzer : public DileptonSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DileptonSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {
        addUncVariables = (treeType == TREE_OTHER);
        filename = fileName;
    }

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed, CORRTYPE jerUNC, CORRTYPE jesUNC,CORRTYPE metUNC) : DileptonSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
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
        filename = fileName;

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
        i_passPre     =outTree->add<size8>  ("","passPre"   ,"b",0);
        i_ht          =outTree->add<float>  ("","ht"        ,"F",0);
        i_met         =outTree->add<float>  ("","met"       ,"F",0);
        i_isMuon1      =outTree->add<size8>  ("","isMuon1"    ,"b",0);
        i_isMuon2      =outTree->add<size8>  ("","isMuon2"    ,"b",0);
        i_lep1PT       =outTree->add<float>  ("","lep1PT"     ,"F",0);
        i_lep2PT       =outTree->add<float>  ("","lep2PT"     ,"F",0);
        i_dilepPT     =outTree->add<float>  ("","dilepPT"    ,"F",0);
        i_dilepMass   =outTree->add<float>  ("","dilepMass"  ,"F",0);
        i_dR_ll       =outTree->add<float>  ("","dilepDR"    ,"F",0);
        i_dPhi_metll  =outTree->add<float>  ("","dPhi_met_dilep"    ,"F",0);

        i_hbbMass     =outTree->add<float>  ("","hbbMass"   ,"F",0);
        i_hbbPT       =outTree->add<float>  ("","hbbPT"     ,"F",0);
        i_hbbCSVCat   =outTree->add<size8>  ("","hbbCSVCat" ,"b",0);

        i_hhMass      =outTree->add<float>  ("","hhMass"    ,"F",0);
        i_hwwPT       =outTree->add<float>  ("","hwwPT"     ,"F",0);

        i_nAK4Btags   =outTree->add<size8>  ("","nAK4Btags" ,"b",0);

        if(!isRealData()){
            i_hbbDecayTypeMC   =outTree->add<size8>  ("","hbbDecayTypeMC"  ,"b",0);
            i_nLepsTT  =outTree->add<size8>  ("","nLepsTT" ,"b",0);
        }

        if(addUncVariables){
            i_w_muIDUp        = outTree->add<float>("","w_muIDUp","F",0);
            i_w_muISOUp       = outTree->add<float>("","w_muISOUp","F",0);
            i_w_elRecoUp      = outTree->add<float>("","w_elRecoUp","F",0);
            i_w_elIDUp        = outTree->add<float>("","w_elIDUp","F",0);
            i_w_elISOUp       = outTree->add<float>("","w_elISOUp","F",0);
            i_w_b_realUp     = outTree->add<float>("","w_b_realUp","F",0);
            i_w_b_fakeUp      = outTree->add<float>("","w_b_fakeUp","F",0);
            i_w_puUp          = outTree->add<float>("","w_puUp","F",0);
            i_w_muIDDown      = outTree->add<float>("","w_muIDDown","F",0);
            i_w_muISODown     = outTree->add<float>("","w_muISODown","F",0);
            i_w_elRecoDown    = outTree->add<float>("","w_elRecoDown","F",0);
            i_w_elIDDown      = outTree->add<float>("","w_elIDDown","F",0);
            i_w_elISODown     = outTree->add<float>("","w_elISODown","F",0);
            i_w_b_realDown   = outTree->add<float>("","w_b_realDown","F",0);
            i_w_b_fakeDown    = outTree->add<float>("","w_b_fakeDown","F",0);
            i_w_puDown        = outTree->add<float>("","w_puDown","F",0);
            i_w_scale        = outTree->addMulti<float>("","w_scale",0);
            i_w_pdf          = outTree->addMulti<float>("","w_pdf",0);
        }

    }

    bool runEvent() override {
        bool passPre = true;
        if(!DileptonSearchRegionAnalyzer::runEvent()) passPre = false;
        if(!passTriggerPreselection) passPre = false;
        if(!passEventFilters) passPre = false;
        if(selectedDileptons.size() != 2) passPre = false;
        if(!hbbCand) passPre = false;
        if (ht_puppi < 400) passPre = false;

        if(dilepChan == ee && passPre) {
        	if (!(((const Electron*)selectedDileptons[0])->passMedID_noISO() && ((const Electron*)selectedDileptons[1])->passMedID_noISO())) passPre = false;
        }

        if(!addUncVariables && !passPre) return false;
        outTree->fill(i_passPre     ,size8(passPre));


        if(isRealData()){
            outTree->fill(i_dataset     ,reader_event->dataset);
            outTree->fill(i_dataRun     ,reader_event->dataRun);
        } else {
            outTree->fill(i_process     ,reader_event->process);
            outTree->fill(i_dhType      ,size8(diHiggsEvt.type));
            outTree->fill(i_xsec        ,float( EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi())));
            outTree->fill(i_trig_N      ,float(smDecayEvt.promptElectrons.size() + smDecayEvt.promptMuons.size() ? trigSFProc->getLeptonTriggerSF(ht_puppi, (selectedDileptons[0] && selectedDileptons[0]->isMuon())) : 1.0 ));
            outTree->fill(i_pu_N        ,float(puSFProc->getCorrection(reader_event->nTruePUInts,CorrHelp::NOMINAL)));
            outTree->fill(i_lep_N       ,float(leptonSFProc->getSF()));
            outTree->fill(i_btag_N      ,float(sjbtagSFProc->getSF({hbbCand})*ak4btagSFProc->getSF(jets_HbbV)));
        }

        outTree->fill(i_ht     ,float(ht_puppi));
        outTree->fill(i_met    ,float(reader_event->met.pt()));
        if(selectedDileptons.front() && selectedDileptons.back()){
            outTree->fill(i_isMuon1      ,size8(selectedDileptons.front()->isMuon()));
            outTree->fill(i_isMuon2      ,size8(selectedDileptons.back()->isMuon()));
            outTree->fill(i_lep1PT      ,float(selectedDileptons.front()->pt()));
            outTree->fill(i_lep2PT      ,float(selectedDileptons.back()->pt()));

            outTree->fill(i_dilepPT, float((selectedDileptons[0]->p4()+selectedDileptons[1]->p4()).pt()));
            outTree->fill(i_dilepMass, float((selectedDileptons[0]->p4()+selectedDileptons[1]->p4()).mass()));
            outTree->fill(i_dR_ll, float(PhysicsUtilities::deltaR(*selectedDileptons[0],*selectedDileptons[1])));
            outTree->fill(i_dPhi_metll, float(PhysicsUtilities::deltaPhi(reader_event->met.p4(),(selectedDileptons[0]->p4()+selectedDileptons[1]->p4()))));

        }

        outTree->fill(i_hbbMass     ,float(hbbMass));
        if(hbbCand) {
            outTree->fill(i_hbbPT       ,float(hbbCand->pt()));
            outTree->fill(i_hbbCSVCat   ,size8(hbbCSVCat));

            outTree->fill(i_hhMass      ,float(hh.mass()));
            outTree->fill(i_hwwPT       ,float(hww.pt()));

            outTree->fill(i_nAK4Btags   ,size8(std::min(nMedBTags_HbbV,250)));
        }

        if(!isRealData() && hbbCand){
    		double maxDR2 = 0.8*0.8;

            int topDecayType = 0; // NONE b wj wjb wjj wjjb bb wjbb wjjbb
            int maxQuarksFromTop = 0;
            int totQuarksFromTops = 0;
            int numB = 0;

            for (const auto& d : smDecayEvt.topDecays) {
            	if (d.type == TopDecay::BAD) continue;
            	if (d.type > TopDecay::HAD) {
            		if (PhysicsUtilities::deltaR2(*d.b,*hbbCand) < maxDR2) {
            			totQuarksFromTops++;
            			numB++;
            			if (maxQuarksFromTop == 0) maxQuarksFromTop = 1;
            		}
            	} else {
            		if (!d.b) continue;
            		bool passB = false;
            		if (PhysicsUtilities::deltaR2(*d.b,*hbbCand) < maxDR2) passB = true;
            		int nW = (PhysicsUtilities::deltaR2(*d.W_decay.dau1,*hbbCand) < maxDR2) + (PhysicsUtilities::deltaR2(*d.W_decay.dau2,*hbbCand) < maxDR2);
            		int nT = nW + passB;
            		totQuarksFromTops += nT;
            		numB += passB;
            		if (nT > maxQuarksFromTop) maxQuarksFromTop = nT;
            	}
            }

            if (totQuarksFromTops == 1){
                if(numB == 1) topDecayType = 1;
                else topDecayType = 2;
            } else if(totQuarksFromTops == 2){
                if (numB == 1) topDecayType = 3;
                else if (numB == 2) topDecayType = 6;
                else topDecayType = 4;
            } else if(totQuarksFromTops == 3) {
            	if (numB == 1) topDecayType = 5;
            	else topDecayType = 7;
            } else if (totQuarksFromTops == 4) topDecayType = 8;

            int WDecayType = 0; // NONE b wj wjb wjj wjjb bb wjbb wjjbb
            int maxQuarksFromW = 0;
            int totQuarksFromWs = 0;

            for (const auto& d : smDecayEvt.bosonDecays) {
                if(d.type != BosonDecay::Z_HAD && d.type != BosonDecay::W_HAD ) continue;
                int nW = (PhysicsUtilities::deltaR2(*d.dau1,*hbbCand) < maxDR2) +  (PhysicsUtilities::deltaR2(*d.dau2,*hbbCand) < maxDR2);
                totQuarksFromWs += nW;
                if (nW > maxQuarksFromW) maxQuarksFromW = nW;
            }

            if (totQuarksFromWs == 1) WDecayType = 2;
            else if (totQuarksFromWs == 2) WDecayType = 4;

            int decayType = 0;

            if(maxQuarksFromTop >= maxQuarksFromW){
                decayType = topDecayType;
            } else {
                decayType = WDecayType;
            }
            outTree->fill(i_hbbDecayTypeMC, size8(decayType));

            if (mcProc == FillerConstants::TTBAR) {
                TPRegexp r1("(.*)(tbar)(.*)(\\d)(l)(.*)");
                auto b = r1.MatchS(filename);
                int nLepsTT = (((TObjString *)b->At(4))->GetString()).Atoi();
            	outTree->fill(i_nLepsTT, size8(nLepsTT));
            }
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
            outTree->fill(i_w_puUp       ,float(puSFProc->getCorrection(reader_event->nTruePUInts,CorrHelp::UP)));

            outTree->fill(i_w_muIDDown   ,float(leptonSFProc->getMuonSF(NONE,DOWN,NOMINAL)*nomEl));
            outTree->fill(i_w_muISODown  ,float(leptonSFProc->getMuonSF(NONE,NOMINAL,DOWN)*nomEl));
            outTree->fill(i_w_elRecoDown ,float(leptonSFProc->getElectronSF(DOWN,NOMINAL,NOMINAL)*nomMu));
            outTree->fill(i_w_elIDDown   ,float(leptonSFProc->getElectronSF(NOMINAL,DOWN,NOMINAL)*nomMu));
            outTree->fill(i_w_elISODown  ,float(leptonSFProc->getElectronSF(NOMINAL,NOMINAL,DOWN)*nomMu));
            outTree->fill(i_w_b_realDown ,float(ak4btagSFProc->getSF(jets_HbbV,NOMINAL,DOWN)* sjbtagSFProc->getSF({hbbCand},NOMINAL,DOWN)));
            outTree->fill(i_w_b_fakeDown ,float(ak4btagSFProc->getSF(jets_HbbV,DOWN,NOMINAL)* sjbtagSFProc->getSF({hbbCand},DOWN,NOMINAL)));
            outTree->fill(i_w_puDown     ,float(puSFProc->getCorrection(reader_event->nTruePUInts,CorrHelp::DOWN)));
        }
        return true;
    }



    //Event information and weights
    size i_process    = 0;
    size i_dhType     = 0;
    size i_dataset    = 0;
    size i_dataRun    = 0;
    size i_xsec       = 0;
    size i_trig_N     = 0;
    size i_pu_N       = 0;
    size i_lep_N      = 0;
    size i_btag_N     = 0;
    size i_passPre    = 0;

    //SR variables
    size i_ht         = 0;
    size i_met        = 0;
    size i_isMuon1    = 0;
    size i_isMuon2    = 0;
    size i_lep1PT     = 0;
    size i_lep2PT     = 0;
    size i_dilepPT    = 0;
    size i_dilepMass  = 0;
    size i_dR_ll      = 0;
    size i_dPhi_metll = 0;

    size i_hbbMass   = 0;
    size i_hbbPT     = 0;
    size i_hbbCSVCat = 0;

    size i_hhMass    = 0;
    size i_hwwPT     = 0;

    size i_nAK4Btags = 0;

    //BE extra variables
    size i_hbbDecayTypeMC   =0;
    size i_nLepsTT  =0;

    //systematic variables
    size i_w_muIDUp      = 0;
    size i_w_muISOUp     = 0;
    size i_w_elRecoUp    = 0;
    size i_w_elIDUp      = 0;
    size i_w_elISOUp     = 0;
    size i_w_b_realUp    = 0;
    size i_w_b_fakeUp    = 0;
    size i_w_puUp        = 0;
    size i_w_muIDDown    = 0;
    size i_w_muISODown   = 0;
    size i_w_elRecoDown  = 0;
    size i_w_elIDDown    = 0;
    size i_w_elISODown   = 0;
    size i_w_b_realDown  = 0;
    size i_w_b_fakeDown  = 0;
    size i_w_puDown      = 0;

    size i_w_scale       = 0;
    size i_w_pdf         = 0;

    bool addUncVariables = false;
    TString filename = "";



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



void makeDileptonBETrees(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    doOne(fileName,treeInt,randSeed,outFileName);
//    if(treeInt==2){
//        size_t lastindex = outFileName.find_last_of(".");
//        std::string extLessName = outFileName.substr(0, lastindex);
//        doOneVar(fileName,treeInt,randSeed+1,extLessName+"_JERUp.root"  ,UP     ,NONE,NONE);
//        doOneVar(fileName,treeInt,randSeed+2,extLessName+"_JERDown.root",DOWN   ,NONE,NONE);
//        doOneVar(fileName,treeInt,randSeed+3,extLessName+"_JESDOWN.root",NOMINAL,DOWN,NONE);
//        doOneVar(fileName,treeInt,randSeed+4,extLessName+"_JESUp.root"  ,NOMINAL,UP  ,NONE);
//        doOneVar(fileName,treeInt,randSeed+5,extLessName+"_METDOWN.root",NOMINAL,NONE,DOWN);
//        doOneVar(fileName,treeInt,randSeed+6,extLessName+"_METUp.root"  ,NOMINAL,NONE,UP  );
//    }
}
void makeDileptonBETrees(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    doOne(fileName,treeInt,randSeed,outFileName,xSec,numEvent);
//    if(treeInt==2){
//        size_t lastindex = outFileName.find_last_of(".");
//        std::string extLessName = outFileName.substr(0, lastindex);
//        doOneVar(fileName,treeInt,randSeed+1,extLessName+"_JERUp.root"  ,UP     ,NONE,NONE,xSec,numEvent);
//        doOneVar(fileName,treeInt,randSeed+2,extLessName+"_JERDown.root",DOWN   ,NONE,NONE,xSec,numEvent);
//        doOneVar(fileName,treeInt,randSeed+3,extLessName+"_JESDOWN.root",NOMINAL,DOWN,NONE,xSec,numEvent);
//        doOneVar(fileName,treeInt,randSeed+4,extLessName+"_JESUp.root"  ,NOMINAL,UP  ,NONE,xSec,numEvent);
//        doOneVar(fileName,treeInt,randSeed+5,extLessName+"_METDOWN.root",NOMINAL,NONE,DOWN,xSec,numEvent);
//        doOneVar(fileName,treeInt,randSeed+6,extLessName+"_METUp.root"  ,NOMINAL,NONE,UP  ,xSec,numEvent);
//    }

}
