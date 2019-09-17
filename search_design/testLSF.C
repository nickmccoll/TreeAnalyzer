
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
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/JetReader.h"

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

    virtual void setupParameters() override {
        DefaultSearchRegionAnalyzer::setupParameters();
        parameters.leptons.mu_getISO = &Muon::inclIso;
        parameters.leptons.el_getISO = &Electron::inclIso;
        parameters.event.doTTBarStitching = false;
    }

    virtual void loadVariables() override {
        reader_event       =loadReader<EventReader>   ("event",isRealData());
        reader_fatjet      =loadReader<FatJetReader>  ("ak8PuppiJet",isRealData(),false,true,true);
        reader_fatjet_noLep      =loadReader<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false,false,false);
        reader_jet_chs     =loadReader<JetReader>     ("ak4Jet",isRealData());
        reader_electron    =loadReader<ElectronReader>("electron");
        reader_muon        =loadReader<MuonReader>    ("muon");

        if(!isRealData()){
            reader_genpart =loadReader<GenParticleReader>   ("genParticle");
        }

        checkConfig();
    }


    const Lepton * getMatchedLepton(const GenParticle& genLepton){
        bool isGenMuon = genLepton.absPdgId() == ParticleInfo::p_muminus;
        auto recoLeptons = PhysicsUtilities::selObjsD(selectedLeptons,
                [&](const Lepton* o){return o->isMuon() == isGenMuon;} );
        double nearestDR =10;
        int idx = PhysicsUtilities::findNearestDRDeref(genLepton,recoLeptons,nearestDR,0.2);
        if(idx < 0) return 0;
        return recoLeptons[idx];
    }

    const FatJet * getMatchedJet(const Lepton* l){
        double nearestDR =10;
        int idx = PhysicsUtilities::findNearestDR(*l,reader_fatjet->jets,nearestDR,1.5);
        if(idx < 0) return 0;
        return &(reader_fatjet->jets[idx]);
    }

    const FatJet * getMatchedNoLepJet(const Lepton* l){
        double nearestDR =10;
        int idx = PhysicsUtilities::findNearestDR(*l,reader_fatjet_noLep->jets,nearestDR,1.5);
        if(idx < 0) return 0;
        return &(reader_fatjet_noLep->jets[idx]);
    }

    double getLSF(const Lepton* l){
        auto fj = getMatchedJet(l);
        return fj ? reader_fatjet->lsf3[fj->index()] : -1;
    }


    void plotVars(TString prefix, const Lepton* l ){
        plotter.getOrMake1DPre(prefix,"miniIso","miniIso", 808,-0.01,1)->Fill(l->miniIso(),weight);
        plotter.getOrMake1DPre(prefix,"pfIso","pfIso", 808,-0.01,1)->Fill(l->pfIso(),weight);

        auto fj = getMatchedJet(l);
        auto fjNL = getMatchedNoLepJet(l);

        std::cout<< l->isMuon() <<" "<< *l << " "<<l->miniIso() <<std::endl;

        if(fj)
        std::cout<< "FJ "<< PhysicsUtilities::deltaR(*fj,*l) <<" " << *fj <<" "<<fj->nSubJets()<<" "<< fj->sdMom().mass() << " "<<  reader_fatjet->lsf3[fj->index()]
                 <<" "<<  reader_fatjet->dRLep[fj->index()]<<" "<<  reader_fatjet->tau3[fj->index()]
                 <<" "<<  reader_fatjet->tau2[fj->index()]<<" "<<  reader_fatjet->tau1[fj->index()]
                 <<" "<<  reader_fatjet->tau3[fj->index()]/reader_fatjet->tau2[fj->index()]
                 <<" "<<  reader_fatjet->tau2[fj->index()]/reader_fatjet->tau1[fj->index()]
                                                                               << std::endl;
        if(fjNL)
        std::cout<< "FJNL "<< PhysicsUtilities::deltaR(*fjNL,*l) <<" "<< *fjNL <<" "<<fjNL->nSubJets()<<" "<< fjNL->sdMom().mass()
                 <<" "<<  reader_fatjet_noLep->tau2[fj->index()]<<" "<<  reader_fatjet_noLep->tau1[fj->index()]
                 <<" "<<  reader_fatjet_noLep->tau2[fj->index()]/reader_fatjet_noLep->tau1[fj->index()]
                                                                               << std::endl;

        for(unsigned int iJ = 0; iJ < reader_fatjet->jets.size(); ++iJ){
            std:: cout << reader_fatjet->jets[iJ] <<" " <<  reader_fatjet->lsf3[reader_fatjet->jets[iJ].index()]  <<" "<<  reader_fatjet->dRLep[reader_fatjet->jets[iJ].index()]<<std::endl;
        }

        plotter.getOrMake1DPre(prefix,"lsf","lsf", 808,-0.01,1)
                ->Fill(fj ? reader_fatjet->lsf3[fj->index()] : -1,weight);

        plotter.getOrMake1DPre(prefix,"tau32","tau32", 808,-0.01,1)
                ->Fill(fj ? reader_fatjet->lsf3[fj->index()] : -1,weight);


        if(fj){
            plotter.getOrMake1DPre(prefix,"wfj_miniIso","miniIso", 808,-0.01,1)->Fill(l->miniIso(),weight);
            plotter.getOrMake1DPre(prefix,"wfj_pfIso","pfIso", 808,-0.01,1)->Fill(l->pfIso(),weight);
            plotter.getOrMake1DPre(prefix,"wfj_lsf","lsf", 808,-0.01,1)
                    ->Fill(fj ? reader_fatjet->lsf3[fj->index()] : -1,weight);
            plotter.getOrMake1DPre(prefix,"wfj_tau32","tau32", 808,-0.01,1)
                    ->Fill(fj ? reader_fatjet->lsf3[fj->index()] : -1,weight);

        }
    }

    void plotLepton(TString prefix, const Lepton* l ){
        TString lepType = l->isMuon()? "mu" : "e";
        plotVars(prefix+"_emu_inclPT",l);
//        plotVars(prefix+"_"+ lepType +"_inclPT",l);
//        TString ptStr = "ptgeq100";
//        if(l->pt() < 50) ptStr = "ptlt50";
//        else if(l->pt() < 100) ptStr = "pt50to100";
//        plotVars(prefix+"_emu_" + ptStr,l);
//        plotVars(prefix+"_"+ lepType +"_"+ptStr,l);
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht_chs < 500) return false;
        const bool passTight = ht_chs >= 1500;
        TString prefix = smpName;

        if(isSignal() && diHiggsEvt.type >= DiHiggsEvent::TAU_MU){
            const auto* recoL = getMatchedLepton(*diHiggsEvt.w1_d1);
            if(recoL) {
//                std::cout << recoL <<" "<< recoL->isMuon()<<" "<< recoL->pt() <<std::endl;
//                std::cout << recoL->miniIso() <<std::endl;
//                std::cout << recoL->pfIso() <<std::endl;
//
//                auto fj = getMatchedJet(recoL);
//                std::cout << fj <<std::endl;
//                if(fj) std::cout << fj->index() <<std::endl;
//                std::cout << reader_fatjet->lsf3.size()<<std::endl;
//                std::cout << "pass" <<std::endl;

//                std::cout << getLSF(recoL) <<std::endl;
                plotLepton(smpName,recoL);
//                if(passTight)
//                    plotLepton(smpName+"_htgeq1500",recoL);
            }
        }

        if(!isSignal()){
            for(const auto *l : selectedLeptons){
                plotLepton(smpName,l);
                if(passTight)
                    plotLepton(smpName+"_htgeq1500",l);
            }

        }

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void testLSF(std::string fileName, int treeInt, int randSeed, std::string outFileName,
        float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
