
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeCopier.h"
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
#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Variables/interface/DileptonSelection.h"
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
        filename = fileName;
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
        	outTree->addSingle(dataset_,  "",  "dataset");
        	outTree->addSingle(dataRun_,  "",  "dataRun");
        } else {
        	outTree->addSingle(process_,  "",  "process");
        	outTree->addSingle(dhType_,  "",  "dhType");
        	outTree->addSingle(xsec_,  "",  "xsec");
        	outTree->addSingle(trig_N_,  "",  "trig_N");
        	outTree->addSingle(pu_N_,  "",  "pu_N");
        	outTree->addSingle(lep_N_,  "",  "lep_N");
        	outTree->addSingle(btag_N_,  "",  "btag_N");
        }
        outTree->addSingle(passPre_, "", "passPre");
        outTree->addSingle(ht_, "", "ht");
        outTree->addSingle(met_, "", "met");
        outTree->addSingle(isMuon1_, "", "isMuon1");
        outTree->addSingle(isMuon2_, "", "isMuon2");
        outTree->addSingle(lep1PT_, "", "lep1PT");
        outTree->addSingle(lep2PT_, "", "lep2PT");
        outTree->addSingle(dilepPT_, "", "dilepPT");
        outTree->addSingle(dilepMass_, "", "dilepMass");
        outTree->addSingle(dilepDR_, "", "dilepDR");
        outTree->addSingle(dPhi_metll_, "", "dPhi_metll");
        outTree->addSingle(hbbMass_, "", "hbbMass");
        outTree->addSingle(hbbPT_, "", "hbbPT");
        outTree->addSingle(hbbCSVCat_, "", "hbbCSVCat");
        outTree->addSingle(hhMass_, "", "hhMass");
        outTree->addSingle(hwwPT_, "", "hwwPT");
        outTree->addSingle(nAK4Btags_, "", "nAK4Btags");

        if(!isRealData()){
            outTree->addSingle(hbbDecayTypeMC_, "", "hbbDecayTypeMC");
            outTree->addSingle(numBinHbb_, "", "numBinHbb");
            outTree->addSingle(nLepsTT_, "", "nLepsTT");
        }

        if(addUncVariables){
        	outTree->addSingle(w_muIDUp_,  "",  "w_muIDUp");
        	outTree->addSingle(w_muISOUp_,  "",  "w_muISOUp");
        	outTree->addSingle(w_elRecoUp_,  "",  "w_elRecoUp");
        	outTree->addSingle(w_elIDUp_,  "",  "w_elIDUp");
        	outTree->addSingle(w_elISOUp_,  "",  "w_elISOUp");
        	outTree->addSingle(w_b_realUp_,  "",  "w_b_realUp");
        	outTree->addSingle(w_b_fakeUp_,  "",  "w_b_fakeUp");
        	outTree->addSingle(w_puUp_,  "",  "w_puUp");
        	outTree->addSingle(w_puDown_,  "",  "w_puDown");
        	outTree->addSingle(w_muIDDown_,  "",  "w_muIDDown");
        	outTree->addSingle(w_muISODown_,  "",  "w_muISODown");
        	outTree->addSingle(w_elRecoDown_,  "",  "w_elRecoDown");
        	outTree->addSingle(w_elIDDown_,  "",  "w_elIDDown");
        	outTree->addSingle(w_elISODown_,  "",  "w_elISODown");
        }

    }

    bool runEvent() override {
        bool passPre = true;
        if(!DefaultSearchRegionAnalyzer::runEvent()) passPre = false;
        if(!passTriggerPreselection) passPre = false;
        if(!passEventFilters) passPre = false;
        if(selectedDileptons.size() != 2) passPre = false;
        if(!hbbCand_2l) passPre = false;

        // muon ID preselection
        if(selectedDileptons.size() == 2) {
        	if (selectedDileptons[0]->isMuon() && selectedDileptons[1]->isMuon()) {
        		bool pass1 = ((Muon*)selectedDileptons[0])->passMedID() && selectedDileptons[0]->pt() > 27;
        		bool pass2 = ((Muon*)selectedDileptons[1])->passMedID() && selectedDileptons[1]->pt() > 27;
        		if (!pass1 && !pass2) passPre = false;
        	}
        }
        if (selectedDileptons.size() == 2) {
            float ptthresh = (selectedDileptons[0]->isMuon() ? 27 : 30);
            if (selectedDileptons[0]->pt() < ptthresh) passPre = false;
        }

        if(!addUncVariables && !passPre) return false;
        passPre_ = size8(passPre);
        ht_ = float(ht_puppi);

        if(isRealData()){
        	dataset_ = size8(*reader_event->dataset);
        	dataRun_ = size8(*reader_event->dataRun);

        } else {
        	process_ = size8(*reader_event->process);
        	dhType_  = size8(diHiggsEvt.type);
        	xsec_    = float( EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),parameters.event,smDecayEvt.genMtt,smDecayEvt.nLepsTT));
        	trig_N_  = float(smDecayEvt.promptElectrons.size() + smDecayEvt.promptMuons.size() ? trigSFProc->getLeptonTriggerSF(ht_, (selectedDileptons.size() && selectedDileptons[0]->isMuon())) : 1.0 );
        	pu_N_    = float(puSFProc->getCorrection(*reader_event->nTruePUInts,CorrHelp::NOMINAL));
        	lep_N_   = 1.0 /*float(leptonSFProc->getSF())*/;
        	btag_N_  = 1.0 /*float(sjbtagSFProc->getSF(parameters.jets,{hbbCand})*ak4btagSFProc->getSF(jets_HbbV))*/;

        }

        met_ = reader_event->met.pt();

        if(selectedDileptons.size() == 2 && selectedDileptons.front() && selectedDileptons.back()){
        	isMuon1_ = size8(selectedDileptons.front()->isMuon());
        	isMuon2_ = size8(selectedDileptons.back()->isMuon());
        	lep1PT_  = float(selectedDileptons.front()->pt());
        	lep2PT_  = float(selectedDileptons.back()->pt());

        	dilepPT_   = float((selectedDileptons[0]->p4()+selectedDileptons[1]->p4()).pt());
        	dilepMass_ = float((selectedDileptons[0]->p4()+selectedDileptons[1]->p4()).mass());
        	dilepDR_   = float(PhysicsUtilities::deltaR(*selectedDileptons[0],*selectedDileptons[1]));
        	dPhi_metll_= float(PhysicsUtilities::deltaPhi(reader_event->met.p4(),(selectedDileptons[0]->p4()+selectedDileptons[1]->p4())));
        }

        hbbMass_ = float(hbbMass_2l);
        if(hbbCand_2l) {
        	hbbPT_ = float(hbbCand_2l->pt());
        	hbbCSVCat_ = size8(hbbCSVCat_2l);
        	hhMass_ = float(hh_2l.mass());
        	hwwPT_  = float(hWW_2l.mass());
        	nAK4Btags_ = size8(std::min(nMedBTags_HbbV_2l,250));

        }

        if(!isRealData() && hbbCand_2l){
    		double maxDR2 = 0.8*0.8;

            int topDecayType = 0; // NONE b wj wjb wjj wjjb bb wjbb wjjbb
            int maxQuarksFromTop = 0;
            int totQuarksFromTops = 0;
            int numB = 0;
            int numB_ll = 0;

            for (const auto& d : smDecayEvt.topDecays) {
            	if (d.type == TopDecay::BAD) continue;
            	if (d.type > TopDecay::HAD) {
            		if (PhysicsUtilities::deltaR2(*d.b,*hbbCand_2l) < maxDR2) {
            			totQuarksFromTops++;
            			numB++; numB_ll++;
            			if (maxQuarksFromTop == 0) maxQuarksFromTop = 1;
            		}
            	} else {
            		if (!d.b) continue;
            		bool passB = false;
            		if (PhysicsUtilities::deltaR2(*d.b,*hbbCand_2l) < maxDR2) passB = true;
            		int nW = (PhysicsUtilities::deltaR2(*d.W_decay.dau1,*hbbCand_2l) < maxDR2) + (PhysicsUtilities::deltaR2(*d.W_decay.dau2,*hbbCand_2l) < maxDR2);
            		int nT = nW + passB;
            		totQuarksFromTops += nT;
            		numB += passB;
            		numB_ll += passB;
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
                int nW = (PhysicsUtilities::deltaR2(*d.dau1,*hbbCand_2l) < maxDR2) +  (PhysicsUtilities::deltaR2(*d.dau2,*hbbCand_2l) < maxDR2);
                if (d.dau1->absPdgId() == ParticleInfo::p_b && PhysicsUtilities::deltaR2(*d.dau1,*hbbCand_2l) < maxDR2) numB_ll++;
                if (d.dau2->absPdgId() == ParticleInfo::p_b && PhysicsUtilities::deltaR2(*d.dau2,*hbbCand_2l) < maxDR2) numB_ll++;

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
            hbbDecayTypeMC_ = size8(decayType);
            numBinHbb_ = size8(numB_ll);
            if (smDecayEvt.nLepsTT == -1) nLepsTT_ = 255;
            else               nLepsTT_ = smDecayEvt.nLepsTT;

        }

        if(addUncVariables){
            const float nomMu = leptonSFProc->getMuonSF();
            const float nomEl = leptonSFProc->getElectronSF();

            w_muIDUp_     = float(leptonSFProc->getMuonSF(NONE,UP,NOMINAL)*nomEl);
            w_muISOUp_    = float(leptonSFProc->getMuonSF(NONE,NOMINAL,UP)*nomEl);
            w_elRecoUp_   = float(leptonSFProc->getElectronSF(UP,NOMINAL,NOMINAL)*nomMu);
            w_elIDUp_     = float(leptonSFProc->getElectronSF(NOMINAL,UP,NOMINAL)*nomMu);
            w_elISOUp_    = float(leptonSFProc->getElectronSF(NOMINAL,NOMINAL,UP)*nomMu);
            w_b_realUp_   = float(ak4btagSFProc->getSF(jets_HbbV,NOMINAL,UP)* sjbtagSFProc->getSF(parameters.jets,{hbbCand},NOMINAL,UP));
            w_b_fakeUp_   = float(ak4btagSFProc->getSF(jets_HbbV,UP,NOMINAL)* sjbtagSFProc->getSF(parameters.jets,{hbbCand},UP,NOMINAL));
            w_puUp_       = float(puSFProc->getCorrection(*reader_event->nTruePUInts,CorrHelp::UP));

            w_muIDDown_   = float(leptonSFProc->getMuonSF(NONE,DOWN,NOMINAL)*nomEl);
            w_muISODown_  = float(leptonSFProc->getMuonSF(NONE,NOMINAL,DOWN)*nomEl);
            w_elRecoDown_ = float(leptonSFProc->getElectronSF(DOWN,NOMINAL,NOMINAL)*nomMu);
            w_elIDDown_   = float(leptonSFProc->getElectronSF(NOMINAL,DOWN,NOMINAL)*nomMu);
            w_elISODown_  = float(leptonSFProc->getElectronSF(NOMINAL,NOMINAL,DOWN)*nomMu);
            w_b_realDown_ = float(ak4btagSFProc->getSF(jets_HbbV,NOMINAL,DOWN)* sjbtagSFProc->getSF(parameters.jets,{hbbCand},NOMINAL,DOWN));
            w_b_fakeDown_ = float(ak4btagSFProc->getSF(jets_HbbV,DOWN,NOMINAL)* sjbtagSFProc->getSF(parameters.jets,{hbbCand},DOWN,NOMINAL));
            w_puDown_     = float(puSFProc->getCorrection(*reader_event->nTruePUInts,CorrHelp::DOWN));
        }

        return true;
    }



    //Event information and weights
    size8 process_    = 0;
    size8 dhType_     = 0;
    size8 dataset_    = 0;
    size8 dataRun_    = 0;
    float xsec_       = 0;
    float trig_N_     = 0;
    float pu_N_       = 0;
    float lep_N_      = 0;
    float btag_N_     = 0;
    size8 passPre_    = 0;

    //SR variables
    float ht_         = 0;
    float met_        = 0;
    size8 isMuon1_    = 0;
    size8 isMuon2_    = 0;
    float lep1PT_     = 0;
    float lep2PT_     = 0;
    float dilepPT_    = 0;
    float dilepMass_  = 0;
    float dilepDR_    = 0;
    float dPhi_metll_ = 0;

    float hbbMass_   = 0;
    float hbbPT_     = 0;
    size8 hbbCSVCat_ = 0;

    float hhMass_    = 0;
    float hwwPT_     = 0;

    size8 nAK4Btags_ = 0;

    //BE extra variables
    size8 hbbDecayTypeMC_   =0;
    size8 numBinHbb_        =0;
    size8 nLepsTT_ = 0;

    //systematic variables
    float w_muIDUp_      = 0;
    float w_muISOUp_     = 0;
    float w_elRecoUp_    = 0;
    float w_elIDUp_      = 0;
    float w_elISOUp_     = 0;
    float w_b_realUp_    = 0;
    float w_b_fakeUp_    = 0;
    float w_puUp_        = 0;
    float w_muIDDown_    = 0;
    float w_muISODown_   = 0;
    float w_elRecoDown_  = 0;
    float w_elIDDown_    = 0;
    float w_elISODown_   = 0;
    float w_b_realDown_  = 0;
    float w_b_fakeDown_  = 0;
    float w_puDown_      = 0;

//    spv_float w_scale_   = make_spv_float();
//    spv_float w_pdf_     = make_spv_float();

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


void makeDileptonBETrees(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
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
