
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/FatJetReader.h"

using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {
    }

//    bool runEvent() override {
//        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
//        if(selectedLeptons.size() != 1) return false;
//        if(diHiggsEvt.type < DiHiggsEvent::MU) return false;
//
//        const float minHT = 400;
//        const float minElePT = 30;
//        const float minMuPT = 26;
//
//        if(ht_chs < minHT) return false;
//
//        float maxElePT = 0;
//        float maxMuPT = 0;
//        for(const auto * l : selectedLeptons ) {
//            if(l->isMuon()) maxMuPT = std::max(maxMuPT, l->pt());
//            else maxElePT = std::max(maxElePT, l->pt());
//        }
//        if(maxElePT < minElePT && maxMuPT < minMuPT ) return false;
//
//
//        std::string prefix = std::string(diHiggsEvt.type == DiHiggsEvent::MU ? "mu" : "e") + "_"+ smpName.Data();
//
//        plotter.getOrMake1DPre(prefix.c_str(),"incl",";HT",36,400,4000)->Fill(ht_chs,weight);
//
//        if(!EventSelection::passTriggerSuite(*reader_event)) return false;
//
//        plotter.getOrMake1DPre(prefix.c_str(),"pass",";HT",36,400,4000)->Fill(ht_chs,weight);
//
//        return true;
//    }



    const Lepton * getMatchedLepton(const GenParticle& genLepton,const std::vector<const Muon *> muons, const std::vector<const Electron*> electrons){
       if(genLepton.absPdgId() == ParticleInfo::p_muminus){
           double nearestDR =10;
           int idx = PhysicsUtilities::findNearestDRDeref(genLepton,muons,nearestDR,0.2);
           if(idx < 0) return 0;
           else return muons[idx];
       } else {
           double nearestDR =10;
           int idx = PhysicsUtilities::findNearestDRDeref(genLepton,electrons,nearestDR,0.2);
           if(idx < 0) return 0;
           else return electrons[idx];
       }
    }

    bool runEvent() override {
        bool passPre = true;
        if(!DefaultSearchRegionAnalyzer::runEvent()) passPre = false;
        const GenParticle * genLep = 0;
        std::string typeStr = "";
        if(diHiggsEvt.type >= DiHiggsEvent::MU ){
            if(diHiggsEvt.type == DiHiggsEvent::MU) typeStr = "HHMu";
            else typeStr = "HHEl";
            genLep = diHiggsEvt.w1_d1;
        } else {
            int nWLeps = 0;
            int nHad = 0;
            int nBad = 0;
            bool isE = false;
            const GenParticle * genLepTemp = 0;
            for(const auto& t : smDecayEvt.topDecays){
                if(t.type >= TopDecay::MU){
                    nWLeps++;
                    genLepTemp = t.W_decay.dau1;
                }
                else if(t.type == TopDecay::HAD) nHad++;
                else nBad++;
                if(t.type == TopDecay::E) isE = true;

            }
            if(nWLeps == 1 && nHad==1 && nBad==0){
                typeStr = isE ? "ZpEl": "ZpMu";
                genLep = genLepTemp;
            } else {
                typeStr = "Bad";
            }
        }

        if(genLep==0) return false;
        const GenParticle* nearQ = 0;
        double minDR2 = 6;
        for(const auto& g : reader_genpart->genParticles){
            if(!ParticleInfo::isDoc(g.status())) continue;
            if(!ParticleInfo::isLastInChain(&g)) continue;
            const int pdg = std::abs(g.pdgId());
            if(!(ParticleInfo::p_d <= pdg && pdg < ParticleInfo::p_t)) continue;
            const float dr2 = PhysicsUtilities::deltaR2(g,*genLep);
            if(dr2 > minDR2) continue;
            minDR2 = dr2;
            nearQ = &g;
        }


        const auto muons = PhysicsUtilities::selObjsMom(reader_muon->muons,20,2.4);
        const auto electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,20,2.4);
        const auto* recoL = getMatchedLepton(*genLep,muons,electrons);

        auto mkPlots = [&](const std::string& pref){
            plotter.getOrMake1DPre(pref.c_str(),"genLepPT",";generated lepton #it{p}_{T}",30,0,600)->Fill(genLep->pt(),weight);
            plotter.getOrMake1DPre(pref.c_str(),"genLepQDR",";#DeltaR (lepton,Q)",30,0,5)->Fill(nearQ ? std::sqrt(minDR2) : 6,weight);
            plotter.getOrMake1DPre(pref.c_str(),"genLepQRelPT",";Q pT / lep PT",30,0,3)->Fill(nearQ ? nearQ->pt()/genLep->pt() : 6,weight);



            if(nearQ){
                plotter.getOrMake2DPre(pref.c_str(),"genLepQDR_genLepQRelPT",";#DeltaR (lepton,nearest quark);quark #it{p}_{T} / lepton #it{p}_{T}",20,0,1,5,0,2)->Fill(std::sqrt(minDR2),nearQ->pt()/genLep->pt(),weight);

            }
            if(nearQ && nearQ->pt()/genLep->pt() > 0.5){
                plotter.getOrMake1DPre(pref.c_str(),"nearQ_genLepQDR",";#DeltaR (lepton,nearest quark)",20,0,1)->Fill(nearQ ? std::sqrt(minDR2) : 6,weight);
            }

        };
        mkPlots(typeStr);
        if(recoL) mkPlots(typeStr+"_reco");

        if(recoL){
            plotter.getOrMake1DPre(typeStr.c_str(),"genLepRecoLepDR",";#DeltaR (gen. lepton,reco. lepton)",40,0,.2)->Fill( PhysicsUtilities::deltaR(*recoL,*genLep),weight);
        }

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getLeptonTest(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void getLeptonTest(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
