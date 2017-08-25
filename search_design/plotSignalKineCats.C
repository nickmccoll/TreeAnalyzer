
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeAnalyzer.h"
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

using namespace TAna;

class Analyzer : public BaseTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
    }
    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );
        reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jetwlep );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData());  load(reader_fatjet  );
    }

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
        const float weight = EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());
        DiHiggsEvent diHiggsEvt; diHiggsEvt.setDecayInfo(reader_genpart->genParticles);

        plotter.getOrMake1D("type",";type; arbitrary units",20,-0.5,19.5 )->Fill(diHiggsEvt.type,weight);

        if(diHiggsEvt.type < DiHiggsEvent::MU) return false;

        bool isMuon = diHiggsEvt.w1_d1->absPdgId() == ParticleInfo::p_muminus;

        double ht = 0;
        double maxAK8 = 0;
        double maxMassiveAK8 = 0;

        for(const auto& j : reader_jetwlep->jets){
            if(j.pt() < 30) continue;
            ht += j.pt();
        }
        double maxAK8PT        = 0;
        double maxTrimmedAK8PT = 0;
        for(const auto& j : reader_fatjet->jets){
            if(j.absEta() <5.0)
                if(j.pt() >maxAK8PT ) maxAK8PT = j.pt();
            if(j.rawSdMom().absEta() < 5.0)
                if(j.rawSdMom().mass() > 60 && j.rawSdMom().pt() > maxTrimmedAK8PT ) maxTrimmedAK8PT = j.rawSdMom().pt();
        }


        auto makePlots =[&](TString lepP, TString prefix, int ecIDX) {
            plotter.getOrMake1D(TString::Format("%s_event_count",lepP.Data()),";event counts; arbitrary units",20,-0.5,19.5 )->Fill(ecIDX,weight);
            plotter.getOrMake1D(TString::Format("%s_%s_genlep_pt",lepP.Data(),prefix.Data()),";gen. lepton #it{p}_{T} [GeV]; arbitrary units",200,0,1000 )->Fill(diHiggsEvt.w1_d1->pt(),weight);

            plotter.getOrMake1D(TString::Format("%s_%s_genlep_pt",lepP.Data(),prefix.Data()),";gen. lepton #it{p}_{T} [GeV]; arbitrary units",200,0,1000 )->Fill(diHiggsEvt.w1_d1->pt(),weight);
            plotter.getOrMake1D(TString::Format("%s_%s_genW_pt"  ,lepP.Data(),prefix.Data()),";gen. W #it{p}_{T} [GeV]; arbitrary units",200,0,2000 )->Fill((diHiggsEvt.w1_d1->p4()+diHiggsEvt.w1_d2->p4()).pt(),weight);
            plotter.getOrMake1D(TString::Format("%s_%s_genH_pt"  ,lepP.Data(),prefix.Data()),";gen. H #it{p}_{T} [GeV]; arbitrary units",300,0,3000 )->Fill(diHiggsEvt.hww->pt(),weight);

            double deltaR = PhysicsUtilities::deltaR(*diHiggsEvt.w1_d1,diHiggsEvt.w2_d1->p4() + diHiggsEvt.w1_d2->p4());
            plotter.getOrMake1D(TString::Format("%s_%s_genLepWDR"        ,lepP.Data(),prefix.Data()),";gen. #DeltaR(l,had W); arbitrary units",600,0,3 )->Fill(deltaR,weight);
            plotter.getOrMake2D(TString::Format("%s_%s_pt_vs_genLepWDR"  ,lepP.Data(),prefix.Data()),"gen. lep #it{p}_{T} [GeV]; gen. #DeltaR(l,had W)",100,0,1000,300,0,3 )->Fill(diHiggsEvt.w1_d1->pt(),deltaR,weight);

            plotter.getOrMake1D(TString::Format("%s_%s_ht"       ,lepP.Data(),prefix.Data()),";#it{H}_{T} [GeV]; arbitrary units",500,0,5000 )->Fill(ht,weight);
            plotter.getOrMake1D(TString::Format("%s_%s_ak8pt"    ,lepP.Data(),prefix.Data()),";max. AK8 #it{p}_{T} [GeV]; arbitrary units",300,0,3000 )->Fill(maxAK8PT,weight);
            plotter.getOrMake1D(TString::Format("%s_%s_sdak8pt"  ,lepP.Data(),prefix.Data()),";max. sd AK8 #it{p}_{T} [GeV]; arbitrary units",300,0,3000 )->Fill(maxTrimmedAK8PT,weight);
        };
        auto make3Plots =[&](TString prefix, int ecIDX) {
            makePlots("all",prefix,ecIDX);
            if(isMuon) makePlots("mu",prefix,ecIDX);
            else  makePlots("el",prefix,ecIDX);
        };



        std::vector<const Electron*> recoElectrons;
        std::vector<const Electron*> IDElectrons;
        std::vector<const Electron*> isoElectrons;
        for(const auto& l : reader_electron->electrons){
            if(l.pt() < 20 ) continue;
            if(l.absEta() > 2.4) continue;
            recoElectrons.push_back(&l);
            if(std::fabs(l.dz()) > .1) continue;
            if(std::fabs(l.d0()) > .05) continue;
            if(!l.passMedID_noISO()) continue;
            IDElectrons.push_back(&l);
            if(l.miniIso() > 0.1) continue;
            isoElectrons.push_back(&l);
        }
        std::vector<const Muon*> recoMuons;
        std::vector<const Muon*> IDMuons;
        std::vector<const Muon*> isoMuons;
        for(const auto& l : reader_muon->muons){
            if(l.pt() < 20 ) continue;
            if(l.absEta() > 2.4) continue;
            recoMuons.push_back(&l);
            if(std::fabs(l.dz()) > .1) continue;
            if(std::fabs(l.d0()) > .05) continue;
            if(!l.passMedID()) continue;
            IDMuons.push_back(&l);
            if(l.miniIso() > 0.2) continue;
            isoMuons.push_back(&l);
        }

        const auto* recoL = getMatchedLepton(*diHiggsEvt.w1_d1,recoMuons,recoElectrons);
        const auto* idL   = getMatchedLepton(*diHiggsEvt.w1_d1,IDMuons,IDElectrons);
        const auto* isoL  = getMatchedLepton(*diHiggsEvt.w1_d1,isoMuons,isoElectrons);


//        if(diHiggsEvt.w1_d1->pt() > 20 && diHiggsEvt.w1_d1->absEta() < 2.4 && recoL == 0){
//            ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
//            std::cout <<"\nlepton: " << *diHiggsEvt.w1_d1 <<" ("<< diHiggsEvt.w1_d1->pdgId() <<")\n";
//            std::cout <<"Jets: ";
//            for(const auto& j : reader_jet->jets){
//                std::cout << j <<" ";
//            }
//            std::cout <<"\nFatJets: " ;
//            for(const auto& j : reader_fatjet->jets){
//                std::cout << j <<" ";
//            }
//            std::cout <<"\nElectrons: " ;
//            for(const auto& j : reader_electron->electrons){
//                std::cout << j <<" ";
//            }
//            std::cout <<"\nMuons: " ;
//            for(const auto& j : reader_muon->muons){
//                std::cout << j <<" ";
//            }
//            std::cout <<"\n---------------------------------------\n\n";
//        }
        make3Plots("incl",0);
        if(diHiggsEvt.w1_d1->pt() > 20 && diHiggsEvt.w1_d1->absEta() < 2.4) make3Plots("genAcc",1);
        if(recoL)  make3Plots("reco",2);
        if(idL)    make3Plots("id",3);
        if(isoL)    make3Plots("iso",4);

        if( (isoL && ht < 1200) || (ht >1200))           make3Plots("iso_incHT",5);
        if( (isoL && ht < 1200) || (idL && ht >1200))    make3Plots("iso_genHT",6);

        if( (isoL && maxAK8PT < 675) || (maxAK8PT >675))    make3Plots("iso_incAK8",7);
        if( (isoL && maxAK8PT < 675) || (idL && maxAK8PT >675))    make3Plots("iso_idAK8",8);
        if( (isoL && maxTrimmedAK8PT < 500) || (maxTrimmedAK8PT >500))    make3Plots("iso_incAK8SD",9);
        if( (isoL && maxTrimmedAK8PT < 500) || (idL && maxTrimmedAK8PT >500))    make3Plots("iso_idAK8SD",10);

        if(isoL && ht > 400)     make3Plots("isoHT400",11);
        if(isoL) if((ht > 400) || (isoL->pt() > 30 && ht < 400))
            make3Plots("isoHT400LowPT",12);

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}

    std::shared_ptr<EventReader      > reader_event    ;
    std::shared_ptr<GenParticleReader> reader_genpart  ;
    std::shared_ptr<ElectronReader   > reader_electron ;
    std::shared_ptr<MuonReader       > reader_muon     ;
    std::shared_ptr<JetReader        > reader_jetwlep      ;
    std::shared_ptr<FatJetReader     > reader_fatjet   ;
    HistGetter plotter;

};

#endif

void plotSignalKineCats(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void plotSignalKineCats(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
