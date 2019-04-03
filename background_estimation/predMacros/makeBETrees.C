
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

    	outTree->addSingle(passPre_,  "",  "passPre");
    	outTree->addSingle(ht_,  "",  "ht");
    	outTree->addSingle(met_,  "",  "met");
    	outTree->addSingle(isMuon_,  "",  "isMuon");
    	outTree->addSingle(lepPT_,  "",  "lepPT");
    	outTree->addSingle(lepETA_,  "",  "lepETA");

    	outTree->addSingle(hbbMass_,  "",  "hbbMass");
    	outTree->addSingle(hbbPT_,  "",  "hbbPT");
    	outTree->addSingle(hbbCSVCat_,  "",  "hbbCSVCat");

    	outTree->addSingle(hhMass_,  "",  "hhMass");
    	outTree->addSingle(wwDM_,  "",  "wwDM");
    	outTree->addSingle(hwwPT_,  "",  "hwwPT");

    	outTree->addSingle(wjjTau2o1_,  "",  "wjjTau2o1");
    	outTree->addSingle(wjjMass_,  "",  "wjjMass");
    	outTree->addSingle(wjjPT_,  "",  "wjjPT");
    	outTree->addSingle(wlnuMass_,  "",  "wlnuMass");
    	outTree->addSingle(wlnuPT_,  "",  "wlnuPT");
    	outTree->addSingle(nAK4Btags_,  "",  "nAK4Btags");

        if(!isRealData()){
        	outTree->addSingle(hbbWQuark_,  "",  "hbbWQuark");
        	outTree->addSingle(hbbWEQuark_,  "",  "hbbWEQuark");
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

//        	outTree->addVector(w_scale_,  "", "w_scale_N", "w_scale",10);
//        	outTree->addVector(w_pdf_,  "", "w_pdf_N", "w_pdf",10);

//            for(unsigned int i = 1; i < 9; ++i){
//                if(i == 5 || i ==7)continue; //told to ignore
//                outTree->addVector((*reader_event->genWeights)[i],  "",   "w_scale");
//            }
//            for(unsigned int i = 111; i < 211; ++i) //Number 110 is the nominal
//                outTree->addVector((*reader_event->genWeights)[i],  "",   "w_pdf");

        }

    }

    bool runEvent() override {
        bool passPre = true;
        if(!DefaultSearchRegionAnalyzer::runEvent()) passPre = false;
        if(!passTriggerPreselection) passPre = false;
        if(!passEventFilters) passPre = false;
        if(selectedLeptons.size() != 1) passPre = false;
        if(!hbbCand)  passPre = false;
        if(!wjjCand)  passPre = false;

        if(!addUncVariables && !passPre) return false;
        passPre_ = size8(passPre);


        if(isRealData()){
        	dataset_ = size8(*reader_event->dataset);
        	dataRun_ = size8(*reader_event->dataRun);

        } else {
        	process_ = size8(*reader_event->process);
        	dhType_  = size8(diHiggsEvt.type);
        	xsec_    = float( EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),parameters.event.lumi));
        	trig_N_  = float(smDecayEvt.promptElectrons.size() + smDecayEvt.promptMuons.size() ? trigSFProc->getLeptonTriggerSF(ht_chs, (selectedLepton && selectedLepton->isMuon())) : 1.0 );
        	pu_N_    = float(puSFProc->getCorrection(*reader_event->nTruePUInts,CorrHelp::NOMINAL));
        	lep_N_   = float(leptonSFProc->getSF());
        	btag_N_  = float(sjbtagSFProc->getSF(parameters.jets,{hbbCand})*ak4btagSFProc->getSF(jets_HbbV));
        }

        ht_  = float(ht_chs);
        met_ = float(reader_event->met.pt());

        if(selectedLepton){
        	isMuon_  = size8(selectedLepton->isMuon());
        	lepPT_   = float(selectedLepton->pt());
        	lepETA_  = float(selectedLepton->eta());
        }

        hbbMass_ = float(hbbMass);

        if(hbbCand) hbbPT_ = float(hbbCand->pt());
        hbbCSVCat_ = size8(hbbCSVCat);

        hhMass_   = float(hh.mass());
        wlnuMass_ = float(wlnu.mass());
        wwDM_     = float(wwDM);
        hwwPT_    = float(hWW.pt());

        if(wjjCand){
        	wjjTau2o1_ = float(wjjCand->tau2otau1());
        	wjjMass_   = float(wjjCand->sdMom().mass());
        	wjjPT_     = float(wjjCand->pt());
        }

        wlnuPT_    = float(wlnu.pt());
        nAK4Btags_ = size8(std::min(nMedBTags_HbbV,250));

        if(!isRealData() && hbbCand){
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

            hbbWQuark_   = size8(decayType);
            hbbWEQuark_  = size8(nExtraQuarks);
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

//            for(unsigned int i = 1; i < 9; ++i){
//                if(i == 5 || i ==7)continue; //told to ignore
//                outTree->fillMulti(i_w_scale,(*reader_event->genWeights)[i]);
//            }
//            for(unsigned int i = 111; i < 211; ++i) //Number 110 is the nominal
//                outTree->fillMulti(i_w_pdf,(*reader_event->genWeights)[i]);


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
    float ht_        = 0;
    float met_       = 0;
    size8 isMuon_    = 0;
    float lepPT_     = 0;
    float lepETA_    = 0;
    float hbbMass_   = 0;
    float hbbPT_     = 0;
    size8 hbbCSVCat_ = 0;

    float hhMass_    = 0;
    float wwDM_      = 0;
    float hwwPT_     = 0;

    float wjjTau2o1_ = 0;
    float wjjMass_   = 0;
    float wjjPT_     = 0;
    float wlnuMass_  = 0;
    float wlnuPT_    = 0;
    size8 nAK4Btags_ = 0;

    //BE extra variables
    size8 hbbWQuark_   =0;
    size8 hbbWEQuark_   =0;

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

void makeBETrees(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
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
