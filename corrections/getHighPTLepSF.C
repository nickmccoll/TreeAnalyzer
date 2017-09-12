
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/FillerConstants.h"
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
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Processors/Variables/interface/BTagging.h"

using namespace TAna;
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){

        leptonProcNoISO .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonProcNoISO);
        leptonProcNoISO->lepSelParams.el_maxISO =-1;
        leptonProcNoISO->lepSelParams.mu_maxISO =-1;
        leptonProcNoISO->lepSelParams_dataABCDEF.el_maxISO =-1;
        leptonProcNoISO->lepSelParams_dataABCDEF.mu_maxISO =-1;
    }


    void loadVariables() override  {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData());  load(reader_fatjet  );
        reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jetwlep );
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );

//        if(!isRealData()){
//            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
//        }
    }
    void doFJ(TString prefix, const FatJet *fj,const std::vector<const Muon*>& muons, const std::vector<const Jet     *>& jets){
        std::vector<const Jet*> btags;
        for(const auto* j : jets){
            if(j->absEta() > 2.4) continue;
            if(j->csv() < BTagging::CSVWP_VALS[BTagging::CSV_M]) continue;
            if(PhysicsUtilities::deltaPhi(*fj,*j) < TMath::PiOver2() ) continue;
            btags.push_back(j);
        }
        if(btags.size() != 1) return;
        std::vector<const Muon*> gmuons;
        for(const auto* m : muons){
            if(PhysicsUtilities::deltaPhi(*fj,*m) < TMath::PiOver2() ) continue;
            if(JetKinematics::transverseMass(*m,reader_event->met) >= 100) continue;
            double nearDR = 10;
            int iJ = PhysicsUtilities::findNearestDRDeref(*m,jets,nearDR);
            if(nearDR > 0.4) continue;
            if(jets[iJ]->index() != btags.front()->index()) continue;
            gmuons.push_back(m);
        }
        plotter.getOrMake1DPre(prefix,"nGoodMuons",";nGoodMuons",10,-0.5,9.5)->Fill(gmuons.size(),weight);

        if(gmuons.size() == 0) return;

        const auto * j = btags.front();
        const auto *l = gmuons.front();



        plotter.getOrMake1DPre(prefix,"mt",";mt",125,0,125)->Fill(JetKinematics::transverseMass(*l,reader_event->met),weight);
        const MomentumF jetMLep = j->p4() - l->p4();
        const float ratio= jetMLep.pt()/l->pt();
        const float rC  = std::max(0.05,std::min(0.2, 10.0/l->pt()));
        const float dR = PhysicsUtilities::deltaR(l->p4(),jetMLep.p4());
        plotter.getOrMake1DPre(prefix,"jml_drN_l",";#DeltaR(jet - lep,lep)/(iso. cone size)",200,0,4)->Fill(dR/rC,weight );
        plotter.getOrMake1DPre(prefix,"jml_dr_l",";#DeltaR(jet - lep,lep)",200,0,0.4)->Fill(dR,weight );

        if(jetMLep.pt() >= 0.5 * l->pt()){
            plotter.getOrMake1DPre(prefix,"pt0p5_jml_drN_l",";#DeltaR(jet - lep,lep)/(iso. cone size)",200,0,4)->Fill(dR/rC,weight );
            plotter.getOrMake1DPre(prefix,"pt0p5_jml_dr_l",";#DeltaR(jet - lep,lep)",200,0,0.4)->Fill(dR,weight );
        }
        if(jetMLep.pt() >=  l->pt()){
            plotter.getOrMake1DPre(prefix,"pt1_jml_drN_l",";#DeltaR(jet - lep,lep)/(iso. cone size)",200,0,4)->Fill(dR/rC,weight );
            plotter.getOrMake1DPre(prefix,"pt1_jml_dr_l",";#DeltaR(jet - lep,lep)",200,0,0.4)->Fill(dR,weight );
        }

        const bool passISO = leptonProc->isGoodLepton(*reader_event,l);
        if(!passISO) return;

        prefix += "_pass";
        plotter.getOrMake1DPre(prefix,"jml_drN_l",";#DeltaR(jet - lep,lep)/(iso. cone size)",200,0,4)->Fill(dR/rC,weight );
        plotter.getOrMake1DPre(prefix,"jml_dr_l",";#DeltaR(jet - lep,lep)",200,0,0.4)->Fill(dR,weight );
        if(jetMLep.pt() >= 0.5 * l->pt()){
            plotter.getOrMake1DPre(prefix,"pt0p5_jml_drN_l",";#DeltaR(jet - lep,lep)/(iso. cone size)",200,0,4)->Fill(dR/rC,weight );
            plotter.getOrMake1DPre(prefix,"pt0p5_jml_dr_l",";#DeltaR(jet - lep,lep)",200,0,0.4)->Fill(dR,weight );
        }
        if(jetMLep.pt() >=  l->pt()){
            plotter.getOrMake1DPre(prefix,"pt1_jml_drN_l",";#DeltaR(jet - lep,lep)/(iso. cone size)",200,0,4)->Fill(dR/rC,weight );
            plotter.getOrMake1DPre(prefix,"pt1_jml_dr_l",";#DeltaR(jet - lep,lep)",200,0,0.4)->Fill(dR,weight );
        }

    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        const bool passJetHT = FillerConstants::doesPass(reader_event->triggerAccepts,FillerConstants::HLT_PFHT800) || FillerConstants::doesPass(reader_event->triggerAccepts,FillerConstants::HLT_PFHT900) ;
        if(!passJetHT) return false;
        if(ht_wlep < 1200) return false;

        if(!isRealData() && reader_event->process >= MCProcess::ZJETS && reader_event->process <= MCProcess::TTX  ) smpName = "other";

        const std::vector<const Jet     *> jets      = JetKinematics::selectObjectsConst(reader_jetwlep->jets,20,10);
        const std::vector<const FatJet  *> fatJets   = JetKinematics::selectObjectsConst(reader_fatjet->jets,400,2.4);
        const std::vector<const Muon*> iDMuons = leptonProcNoISO->getMuons(*reader_event,*reader_muon);

        const FatJet *topJet = 0;
        const FatJet *topJetT = 0;
        for(const auto* fj: fatJets){
            if(fj->sdMom().mass() < 105) continue;
            if(fj->sdMom().mass() > 210) continue;
            if(fj->tau3otau2() > 0.8) continue;
            if(fj->maxSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_L] && topJet == 0 ) {
                topJet = fj;
            }
            if(fj->maxSJCSV() >= BTagging::CSVWP_VALS[BTagging::CSV_M] && topJetT == 0) {
                topJetT = fj;
            }

        }

        if(topJet) doFJ(smpName + "_ltj",topJet,iDMuons,jets);
        if(topJetT) doFJ(smpName + "_mtj",topJet,iDMuons,jets);
        if(topJet){
            plotter.getOrMake1DPre(smpName,"topmass",";topmass",165,75,240)->Fill(topJet->sdMom().mass(),weight);
        plotter.getOrMake1DPre(smpName,"tau32",";tau32",200,0,1)->Fill(topJet->tau3otau2(),weight);
        plotter.getOrMake1DPre(smpName,"csv",";csv",200,0,1)->Fill(topJet->maxSJCSV(),weight);
        plotter.getOrMake1DPre(smpName,"toppt","toppt",200,0,2000)->Fill(topJet->pt(),weight);
        }


        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<LeptonProcessor>     leptonProcNoISO ;
};

#endif

void getHighPTLepSF(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getHighPTLepSF(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000);
    a.write(outFileName);
}
