
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

#include "TPRegexp.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
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
        i_isMuon      =outTree->add<size8>  ("","isMuon"    ,"b",0);
        i_hbbMass     =outTree->add<float>  ("","hbbMass"   ,"F",0);
        i_hbbPT       =outTree->add<float>  ("","hbbPT"     ,"F",0);
        i_hbbNSJs     =outTree->add<size8>  ("","hbbNSJs"   ,"b",0);
        i_hbbCSVCat   =outTree->add<size8>  ("","hbbCSVCat" ,"b",0);
        i_hbbTau2o1   =outTree->add<float>  ("","hbbTau2o1" ,"F",0);

        i_hhMass      =outTree->add<float>  ("","hhMass"    ,"F",0);
        i_wlnuDR      =outTree->add<float>  ("","wlnuDR"    ,"F",0);
        i_wwDM        =outTree->add<float>  ("","wwDM"      ,"F",0);
        i_hwwPT       =outTree->add<float>  ("","hwwPT"   ,"F",0);
        i_wjjCSVCat   =outTree->add<size8>  ("","wjjCSVCat" ,"b",0);
        i_wjjTau2o1   =outTree->add<float>  ("","wjjTau2o1" ,"F",0);
        i_wjjMass     =outTree->add<float>  ("","wjjMass"   ,"F",0);
        i_wjjPT       =outTree->add<float>  ("","wjjPT"     ,"F",0);
        i_wjjNSJs     =outTree->add<size8>  ("","wjjNSJs"   ,"b",0);
        i_wlnuPT      =outTree->add<float>  ("","wlnuPT"    ,"F",0);
        i_nAK4Btags   =outTree->add<size8>  ("","nAK4Btags" ,"b",0);
        i_minBtagMT   =outTree->add<float>  ("","minBtagMT" ,"F",0);

        if(!isRealData()){
            i_hbbGenPT    =outTree->add<float>  ("","hbbGenPT"   ,"F",0);
            i_hbbGenMass  =outTree->add<float>  ("","hbbGenMass" ,"F",0);
            i_hbbWQuark   =outTree->add<size8>  ("","hbbWQuark"  ,"b",0);
            i_hbbWEQuark  =outTree->add<size8>  ("","hbbWEQuark"  ,"b",0);
//            i_hhHT        =outTree->add<float>  ("","hhHT"       ,"F",0);
//            i_wjjlepGenPT =outTree->add<float>  ("","wjjlepGenPT","F",0);
//            i_genMET      =outTree->add<float>  ("","genMET"     ,"F",0);
            i_genhhMass   =outTree->add<float>  ("","genhhMass"  ,"F",0);
            i_genhhMass2   =outTree->add<float>  ("","genhhMass2"  ,"F",0);

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

        outTree->fill(i_ht     ,float(ht_chs));


        outTree->fill(i_isMuon      ,size8(selectedLepton->isMuon()));
        outTree->fill(i_hbbMass     ,float(hbbMass));
        outTree->fill(i_hbbPT       ,float(hbbCand->pt()));
        outTree->fill(i_hbbNSJs     ,size8(hbbNSJs));
        outTree->fill(i_hbbCSVCat   ,size8(hbbCSVCat));
        outTree->fill(i_hbbTau2o1   ,float(hbbCand->tau2otau1()));

        outTree->fill(i_hhMass      ,float(hh.mass()));
        outTree->fill(i_wlnuDR      ,float(wlnuDR));
        outTree->fill(i_wwDM        ,float(wwDM));
        outTree->fill(i_hwwPT       ,float(hWW.pt()));
        outTree->fill(i_wjjCSVCat   ,size8(wjjCSVCat));
        outTree->fill(i_wjjTau2o1   ,float(wjjCand->tau2otau1()));
        outTree->fill(i_wjjMass     ,float(wjjCand->sdMom().mass()));
        outTree->fill(i_wjjPT       ,float(wjjCand->pt()));
        outTree->fill(i_wjjNSJs     ,size8(wjjNSJs));
        outTree->fill(i_wlnuPT      ,float(wlnu.pt()));
        outTree->fill(i_nAK4Btags   ,size8(std::min(nMedBTags_HbbV,250)));
        double minMT = -1;
        for(const auto* j: jets_HbbV ){
            if(!BTagging::isMediumCSVTagged(*j)) continue;
            const double mt = JetKinematics::massiveTransverseMass(j->p4()+selectedLepton->p4(),reader_event->met);
            if(minMT < 0 || mt < minMT) minMT = mt;
        }
        outTree->fill(i_minBtagMT  ,float(minMT));


        if(!isRealData()){
            double nearestDR = 999.;
            int hbbGenIDX = PhysicsUtilities::findNearestDR(*hbbCand,reader_fatjet->genJets,nearestDR,0.8);

            outTree->fill(i_hbbGenPT    ,float(hbbGenIDX < 0 ? 0.0 : reader_fatjet->genJets[hbbGenIDX].pt()));
            outTree->fill(i_hbbGenMass  ,float(hbbGenIDX < 0 ? 0.0 : reader_fatjet->genJets[hbbGenIDX].mass()));

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

//            MomentumF genVisWW;
//            MomentumF genMET;
//            const float matchGJR2 = 1.2*1.2;
//            for(const auto& j : reader_jet_chs->genJets){
//                if(j.pt() < 20) continue;
//                if(j.absEta() > 5.0) continue;
//                genMET.p4() -= j.p4();
//                if(PhysicsUtilities::deltaR2(j,*selectedLepton) < matchGJR2)  genVisWW.p4() += j.p4();
//            }
//
//            const MomentumF genNeutrino    = HiggsSolver::getInvisible(genMET,genVisWW );
//            const MomentumF genhWW   = genNeutrino.p4() + genVisWW.p4();
//            MomentumF genHH = genhWW.p4();
//            if(hbbGenIDX >= 0 ) genHH.p4() += reader_fatjet->genJets[hbbGenIDX].p4();

            MomentumF genMET = reader_event->met.p4() + hbbCand->p4();
            if(hbbGenIDX >= 0) genMET.p4() -=  reader_fatjet->genJets[hbbGenIDX].p4();

            const MomentumF genNeutrino    = HiggsSolver::getInvisible(genMET,(selectedLepton->p4() + wjjCand->p4()) );
            const MomentumF genhWW   = genNeutrino.p4() + selectedLepton->p4() + wjjCand->p4();
            MomentumF genHH = genhWW.p4();
            if(hbbGenIDX >= 0 ) genHH.p4() += reader_fatjet->genJets[hbbGenIDX].p4();
            MomentumF genHH2 = hWW.p4();
            if(hbbGenIDX >= 0 ) genHH2.p4() += reader_fatjet->genJets[hbbGenIDX].p4();

            outTree->fill(i_hbbWQuark   ,decayType);
            outTree->fill(i_hbbWEQuark   ,nExtraQuarks);
//            outTree->fill(i_hhHT        ,float( genhWW.pt()  + (hbbGenIDX >= 0 ? reader_fatjet->genJets[hbbGenIDX].pt() : 0.0)  ));
//            outTree->fill(i_wjjlepGenPT ,float(genVisWW.pt()));
//            outTree->fill(i_genMET      ,float(genMET.pt()));
            outTree->fill(i_genhhMass   ,float(genHH.mass()));
            outTree->fill(i_genhhMass2   ,float(genHH2.mass()));
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
    size i_isMuon    = 0;
    size i_hbbMass   = 0;
    size i_hbbPT     = 0;
    size i_hbbNSJs   =0;
    size i_hbbCSVCat = 0;
    size i_hbbTau2o1 = 0;

    size i_hhMass    = 0;
    size i_wlnuDR    = 0;
    size i_wwDM      = 0;
    size i_hwwPT     = 0;
    size i_wjjCSVCat = 0;
    size i_wjjTau2o1 = 0;
    size i_wjjMass   = 0;
    size i_wjjPT     = 0;
    size i_wjjNSJs   =0;
    size i_wlnuPT    = 0;
    size i_nAK4Btags = 0;
    size i_minBtagMT = 0;

    //BE extra variables
    size i_hbbGenPT    =0;
    size i_hbbGenMass  =0;
    size i_hbbWQuark   =0;
    size i_hbbWEQuark   =0;
//    size i_hhHT        =0;
//    size i_wjjlepGenPT =0;
//    size i_genMET      =0;
    size i_genhhMass   =0;
    size i_genhhMass2   =0;


};

#endif

void makeBETrees(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
}
void makeBETrees(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
}
