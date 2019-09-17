
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

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
    : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {
        turnOffCorr(CORR_TRIG);
        turnOffCorr(CORR_PU);
    }


    virtual void setupParameters() override {
        DefaultSearchRegionAnalyzer::setupParameters();
        parameters.event.doTTBarStitching = false;
    }


    virtual BaseEventAnalyzer * setupEventAnalyzer() override {return new CopierEventAnalyzer();}

    virtual void bookOutputVariables() override {

        outTree->addSingle(process        ,  "",  "process");
        outTree->addSingle(sampParam      ,  "",  "sampParam");
        outTree->addSingle(dhType         ,  "",  "dhType");
        outTree->addSingle(hbbCat         ,  "",  "hbbCat");
        outTree->addSingle(isMuon         ,  "",  "isMuon");
        outTree->addSingle(weight         ,  "",  "weight");
        outTree->addSingle(nAK4Btags      ,  "",  "nAK4Btags");
        outTree->addSingle(hh_orig        ,  "",  "hh_orig");
        outTree->addSingle(hh_chi2        ,  "",  "hh_chi2");
        outTree->addSingle(md             ,  "",  "md");
        outTree->addSingle(chi2           ,  "",  "chi2");
        outTree->addSingle(wqqDR          ,  "",  "wqqDR");
        outTree->addSingle(qqJet_pt       ,  "",  "qqJet_pt");
        outTree->addSingle(qqJet_eta      ,  "",  "qqJet_eta");
        outTree->addSingle(qqJet_phi      ,  "",  "qqJet_phi");
        outTree->addSingle(qqJet_mass     ,  "",  "qqJet_mass");
        outTree->addSingle(qqJet_SDmass   ,  "",  "qqJet_SDmass");
        outTree->addSingle(qqJet_t2ot1    ,  "",  "qqJet_t2ot1");
        outTree->addSingle(lep_pt         ,  "",  "lep_pt");
        outTree->addSingle(lep_eta        ,  "",  "lep_eta");
        outTree->addSingle(lep_phi        ,  "",  "lep_phi");
        outTree->addSingle(met_pt         ,  "",  "met_pt");
        outTree->addSingle(met_phi        ,  "",  "met_phi");
        outTree->addSingle(bbJet_pt       ,  "",  "bbJet_pt");
        outTree->addSingle(bbJet_eta      ,  "",  "bbJet_eta");
        outTree->addSingle(bbJet_phi      ,  "",  "bbJet_phi");
        outTree->addSingle(bbJet_mass     ,  "",  "bbJet_mass");
        outTree->addSingle(bbJet_SDmass   ,  "",  "bbJet_SDmass");
        outTree->addSingle(true_neut_pt   ,  "",  "true_neut_pt");
        outTree->addSingle(true_neut_eta  ,  "",  "true_neut_eta");
        outTree->addSingle(true_neut_phi  ,  "",  "true_neut_phi");
        outTree->addSingle(true_neut_mass ,  "",  "true_neut_mass");
        outTree->addSingle(true_lep_pt    ,  "",  "true_lep_pt");
        outTree->addSingle(true_lep_eta   ,  "",  "true_lep_eta");
        outTree->addSingle(true_lep_phi   ,  "",  "true_lep_phi");
        outTree->addSingle(true_jet_pt    ,  "",  "true_jet_pt");
        outTree->addSingle(true_jet_eta   ,  "",  "true_jet_eta");
        outTree->addSingle(true_jet_phi   ,  "",  "true_jet_phi");
        outTree->addSingle(true_jet_mass  ,  "",  "true_jet_mass");

    }

    void fillTau(const GenParticle* tau, ASTypes::CylLorentzVectorF& np,
            ASTypes::CylLorentzVectorF& lp){
        auto finalTau = ParticleInfo::getFinal(tau);
        for (int i = 0; i<finalTau->numberOfDaughters(); i++) {
            int pp = finalTau->daughter(i)->absPdgId();
            if(ParticleInfo::isANeutrino(finalTau->daughter(i)->absPdgId()))
                np += finalTau->daughter(i)->p4();
            else
                lp += finalTau->daughter(i)->p4();
    }
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passTriggerPreselection) return false;
        if(selectedLeptons.size()!=1) return false;
        if(!hbbCand || !wjjCand) return false;
        if(hbbCSVCat <4) return false;

        //Start with gen info

        ASTypes::CylLorentzVectorF tLepton;
        ASTypes::CylLorentzVectorF tNeutrino;
        ASTypes::CylLorentzVectorF tJet;
        if(isSignal()){
            if( diHiggsEvt.type < DiHiggsEvent::TAU_MU) return false;
            dhType = size8(diHiggsEvt.type);
            if( diHiggsEvt.type < DiHiggsEvent::MU){
                tNeutrino += diHiggsEvt.w1_d2->p4();
                fillTau(diHiggsEvt.w1_d1,tNeutrino,tLepton);
            } else{
                tLepton   += diHiggsEvt.w1_d1->p4();
                tNeutrino += diHiggsEvt.w1_d2->p4();
            }
            tJet = diHiggsEvt.w2_d1->p4()+diHiggsEvt.w2_d2->p4();
        } else if(mcProc == FillerConstants::TTBAR){
            int nLep = 0;
            const TopDecay * lepDecay = 0;
            for(const auto& d : smDecayEvt.topDecays){
                if(d.type==TopDecay::BAD) return false;
                if(d.type==TopDecay::TAU_HAD) return false;
                if(d.type==TopDecay::HAD) continue;
                nLep++;
                lepDecay = &d;
            }
            if(nLep!= 1) return false;
            dhType = size8(lepDecay->type);
            tJet = lepDecay->b->p4();
            if(lepDecay->type<TopDecay::MU){
                tNeutrino += lepDecay->W_decay.dau2->p4();
                fillTau(lepDecay->W_decay.dau1,tNeutrino,tLepton);
            } else{
                tLepton   += lepDecay->W_decay.dau1->p4();
                tNeutrino += lepDecay->W_decay.dau2->p4();
            }
        } else if(mcProc == FillerConstants::WJETS){
            if(smDecayEvt.bosonDecays.size() != 1) return false;
            const auto& wDecay = smDecayEvt.bosonDecays[0];
            if(wDecay.type < BosonDecay::W_TAU_MU ) return false;
            dhType = size8(wDecay.type);
            if(wDecay.type < BosonDecay::W_MU){
                tNeutrino += wDecay.dau2->p4();
                fillTau(wDecay.dau1,tNeutrino,tLepton);
            } else{
                tLepton   += wDecay.dau1->p4();
                tNeutrino += wDecay.dau2->p4();
            }
        } else{
            dhType =0;
        }

        //Now the regular stuff

        process = *reader_event->process;
        sampParam = *reader_event->sampParam;
        hbbCat = hbbCSVCat;
        isMuon = selectedLepton->isMuon();
        weight  =  EventWeights::getNormalizedEventWeight(
                *reader_event,xsec(),nSampEvt(),parameters.event,smDecayEvt.genMtt,smDecayEvt.nLepsTT);

        nAK4Btags = size8(std::min(nMedBTags_HbbV,250));
        hh_orig =hh_basic.mass();
        hh_chi2 =hh_chi.mass();
        md      = wwDM;
        chi2    =hwwChi;
        wqqDR = isSignal() ? std::max(PhysicsUtilities::deltaR(*wjjCand,*diHiggsEvt.w2_d1),
                                      PhysicsUtilities::deltaR(*wjjCand,*diHiggsEvt.w2_d2))
                : -1;

        qqJet_pt     = wjjCand->pt();
        qqJet_eta    = wjjCand->eta();
        qqJet_phi    = wjjCand->phi();
        qqJet_mass   = wjjCand->mass();
        qqJet_SDmass = wjjCand->sdMom().mass();
        qqJet_t2ot1  = wjjCand->tau2otau1();

        lep_pt      = selectedLepton->pt();
        lep_eta     = selectedLepton->eta();
        lep_phi     = selectedLepton->phi();

        met_pt      = reader_event->met.pt();
        met_phi     = reader_event->met.phi();

        bbJet_pt      = hbbCand->pt();
        bbJet_eta     = hbbCand->eta();
        bbJet_phi     = hbbCand->phi();
        bbJet_mass    = hbbCand->mass();
        bbJet_SDmass  = hbbMass;


        true_neut_pt      = tNeutrino.Pt();
        true_neut_eta     = tNeutrino.Eta();
        true_neut_phi     = tNeutrino.Phi();
        true_neut_mass    = tNeutrino.M();

        true_lep_pt      = tLepton.Pt();
        true_lep_eta     = tLepton.Eta();
        true_lep_phi     = tLepton.Phi();

        true_jet_pt    = tJet.Pt();
        true_jet_eta   = tJet.Eta();
        true_jet_phi   = tJet.Phi();
        true_jet_mass  = tJet.M();

//        if(qqJet_t2ot1 < 0.75 && (bbJet_SDmass >= 100 && bbJet_SDmass<=150)  ){
//            double true_sf =  tJet.Pt()/wjjCand->pt();
//            auto tW = tNeutrino + tLepton;
//            bool isVirtualWqq = tJet.mass() < tW.mass();
//
//            ASTypes::CylLorentzVectorF approxJ (true_sf*wjjCand->pt(),wjjCand->eta(),wjjCand->phi(),
//                    isVirtualWqq ? 31 : 80.0);
//
//            auto approxW = selectedLepton->p4() + tNeutrino;
//            auto approxH = approxW + approxJ;
//
//
//
//            if(approxH.mass()>135){
//                        std:: cout << sampParam<<" "<<int(dhType)<<" "<< approxH.pt() <<" "<< approxH.mass()<<" "<< approxW.pt() <<" "<< approxW.mass()
//                                <<" "<< true_sf <<" "<< tJet.Pt() <<" "<<wjjCand->pt() <<" "<< selectedLepton->pt()<< " "<< tLepton.pt()
//                                <<" "<< PhysicsUtilities::deltaR(*wjjCand,tJet) <<" "<<
//
//                                std::max(PhysicsUtilities::deltaR(*wjjCand,*diHiggsEvt.w2_d1),PhysicsUtilities::deltaR(*wjjCand,*diHiggsEvt.w2_d2)) <<std::endl;
//
//                        ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
//
//                        std::cout << "JETS!"<<std::endl;
//                        for(unsigned int iJ = 0; iJ <reader_fatjet_noLep->jets.size(); ++iJ){
//                            std::cout << reader_fatjet_noLep->jets[iJ]<<std::endl;
//                        }
//
//        }
//        }



        return true;
    }


    size8 process        = 0;
    int   sampParam      = 0;
    size8 dhType         = 0;
    size8 hbbCat         = 0;
    size8 isMuon         = 0;
    float weight         = 0;

    size8 nAK4Btags      = 0;
    float hh_orig        = 0;
    float hh_chi2        = 0;
    float md             = 0;
    float chi2           = 0;
    float wqqDR          = 0;

    float qqJet_pt       = 0;
    float qqJet_eta      = 0;
    float qqJet_phi      = 0;
    float qqJet_mass     = 0;
    float qqJet_SDmass   = 0;
    float qqJet_t2ot1    = 0;

    float lep_pt         = 0;
    float lep_eta        = 0;
    float lep_phi        = 0;

    float met_pt         = 0;
    float met_phi        = 0;

    float bbJet_pt       = 0;
    float bbJet_eta      = 0;
    float bbJet_phi      = 0;
    float bbJet_mass     = 0;
    float bbJet_SDmass   = 0;


    float true_neut_pt   = 0;
    float true_neut_eta  = 0;
    float true_neut_phi  = 0;
    float true_neut_mass = 0;


    float true_lep_pt    = 0;
    float true_lep_eta   = 0;
    float true_lep_phi   = 0;

    float true_jet_pt    = 0;
    float true_jet_eta   = 0;
    float true_jet_phi   = 0;
    float true_jet_mass  = 0;
};


#endif

void makeHHSolTree(std::string fileName, int treeInt, int randSeed, std::string outFileName,
        float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    if(xSec > 0) a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
}
