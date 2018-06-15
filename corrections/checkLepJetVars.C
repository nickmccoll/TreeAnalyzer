
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
#include "Processors/GenTools/interface/SMDecayEvent.h"
#include "Processors/Variables/interface/BTagging.h"
#include "TVector3.h"

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

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
        }
        reader_jet_chs =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jet_chs );;
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );
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
    const Jet * getMatchedJet(const MomentumF* mom, const std::vector<const Jet*>& jets, double& nearestDR ){
        int idx = PhysicsUtilities::findNearestDRDeref(mom->p4(),jets,nearestDR);
        if(idx < 0) return 0;
        else return jets[idx];
    }

    void makePlots(const TString& prefix, const MomentumF* leptonMom, const Jet* nearestJet, float ISO = 0){

        const float projLep = nearestJet ? nearestJet->pt()*TMath::Cos( PhysicsUtilities::deltaPhi(*leptonMom,*nearestJet)) : 0.0;
        const MomentumF jetMLep(nearestJet ? nearestJet->p4() - leptonMom->p4() : ASTypes::CylLorentzVectorF() );
        const float projJetMLep = nearestJet ? jetMLep.pt()*TMath::Cos( PhysicsUtilities::deltaPhi(*leptonMom,jetMLep)) : 0.0;
        double rV  = std::max(0.05,std::min(0.2, 10.0/leptonMom->pt()));

        plotter.getOrMake1DPre(prefix,"projlep",";P_{T}",200,0,1000)->Fill(nearestJet ? projLep : 0, weight);
        plotter.getOrMake1DPre(prefix,"projJetMLep",";P_{T}",200,0,1000)->Fill(nearestJet ? projJetMLep : 0, weight);
        plotter.getOrMake1DPre(prefix,"jml_drN_l",";#DeltaR(jet - lep,lep)",200,0,10)->Fill(projLep > leptonMom->pt()  ? PhysicsUtilities::deltaR(jetMLep, *leptonMom)/rV : 10,weight );
        plotter.getOrMake1DPre(prefix,"jml_dr_l",";#DeltaR(jet - lep,lep)",100,0,1)->Fill( projLep > leptonMom->pt()  ? PhysicsUtilities::deltaR(jetMLep, *leptonMom) : 10,weight );


        if(jetMLep.pt() > 0.5*leptonMom->pt()){
            plotter.getOrMake1DPre(prefix,"pt0p5_jml_drN_l",";#DeltaR(jet - lep,lep)",200,0,10)->Fill( projLep > leptonMom->pt()  ? PhysicsUtilities::deltaR(jetMLep, *leptonMom)/rV : 10,weight );
            plotter.getOrMake1DPre(prefix,"pt0p5_jml_dr_l",";#DeltaR(jet - lep,lep)",100,0,1)->Fill( projLep > leptonMom->pt()  ? PhysicsUtilities::deltaR(jetMLep, *leptonMom) : 10,weight );
        }

        // Added by Brent
        plotter.getOrMake1DPre(prefix,"j_dr_l",";#DeltaR(jet,lep)",100,0,0.8)->Fill(nearestJet ? PhysicsUtilities::deltaR(nearestJet->p4(), leptonMom->p4()) : 0,weight);
        plotter.getOrMake1DPre(prefix,"ptj",";p_{T}(jet)",100,0,800)->Fill(nearestJet ? nearestJet->pt() : 0,weight);
        plotter.getOrMake1DPre(prefix,"eta_j",";#eta(jet)",100,-5.,5.)->Fill(nearestJet ? nearestJet->eta() : 0,weight);
        plotter.getOrMake1DPre(prefix,"ptjml",";p_{T}(jet-lep)",100,0,800)->Fill(nearestJet ? jetMLep.pt() : 0,weight);

        //PT rel
        if(nearestJet) {
			TVector3 j(nearestJet->px(),nearestJet->py(),nearestJet->pz());
			TVector3 lep(leptonMom->px(),leptonMom->py(),leptonMom->pz());
			auto jcrossl = j.Cross(lep);
			auto ptrel = jcrossl.Mag() / j.Mag();
			plotter.getOrMake1DPre(prefix,"ptrel",";P_{T,rel}",100,0,150)->Fill(ptrel, weight);
			plotter.getOrMake1DPre(prefix,"jldphi",";#Delta#phi(jet,lep)",100,-3.14, 3.14)->Fill(j.DeltaPhi(lep), weight);
	        plotter.getOrMake1DPre(prefix,"jml_dphi_l",";#Delta#phi(jet-lep,lep)",100,-3.14,3.14)->Fill((j-lep).DeltaPhi(lep), weight);

        }
    }

    void doLepton(const TString& prefix, const GenParticle* genP,const std::vector<const Muon*>& muons,const std::vector<const Electron*>& electrons,
            const std::vector<const Jet*>& jets){

        const auto* recoL = getMatchedLepton(*genP,muons,electrons);
        double nearDRToGen = 20;
        const auto* jetNearGen = getMatchedJet(genP,jets,nearDRToGen);
        double nearDRToReco = 20;
        const auto* jetNearReco = recoL ? getMatchedJet(recoL,jets,nearDRToReco) : 0;

        //reconstruction eff
        makePlots(prefix +"_gen_dr0p4_incl",genP,nearDRToGen < 0.4 ? jetNearGen : 0);
        if(recoL == 0) return;
        makePlots(prefix +"_gen_dr0p4_reco",genP,nearDRToGen < 0.4 ? jetNearGen : 0);

        //ID eff
        makePlots(prefix +"_reco_dr0p4_incl",recoL,nearDRToReco < 0.4 ? jetNearReco : 0);
        plotter.getOrMake1DPre(prefix+"_reco_dr0p4_incl","miniIso",";miniIso",100,0,1)->Fill(recoL->miniIso(),weight);

        const bool passID = leptonProcNoISO->isGoodLepton(*reader_event,recoL);
        if(!passID) return;
        makePlots(prefix +"_reco_dr0p4_ID",recoL,nearDRToReco < 0.4 ? jetNearReco : 0,recoL->miniIso());
        plotter.getOrMake1DPre(prefix+"_reco_dr0p4_ID","miniIso",";miniIso",100,0,1)->Fill(recoL->miniIso(),weight);

        //ISO eff
        const bool passISO = leptonProc->isGoodLepton(*reader_event,recoL);
        if(!passISO) return;
        makePlots(prefix +"_reco_dr0p4_ISOID",recoL,nearDRToReco < 0.4 ? jetNearReco : 0);
        plotter.getOrMake1DPre(prefix+"_reco_dr0p4_ISOID","miniIso",";miniIso",100,0,1)->Fill(recoL->miniIso(),weight);
        
    };

    bool goodGenLepton(const GenParticle* genP) {
        if(genP->absEta() > 2.4) return false;
        if(genP->pt() < ( genP->absPdgId() == ParticleInfo::p_eminus ? 33.0 : 29  )  ) return false;
        return true;
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;

        const std::vector<const Muon    *> muons     = PhysicsUtilities::selObjsMom(reader_muon->muons,26,2.4);
        const std::vector<const Electron*> electrons = PhysicsUtilities::selObjsMom(reader_electron->electrons,30,2.4);
        const std::vector<const Jet     *> jets      = PhysicsUtilities::selObjsMom(reader_jet_chs->jets,20,10);

        if(reader_event->process == FillerConstants::SIGNAL){
            if(diHiggsEvt.type < DiHiggsEvent::TAU_MU) return false;
            if(!goodGenLepton(diHiggsEvt.w1_d1)) return false;

            TString sN;
            if(diHiggsEvt.type == DiHiggsEvent::TAU_MU) {sN = smpName + "_taumu"; ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);}
            else if(diHiggsEvt.type == DiHiggsEvent::TAU_E) sN = smpName + "_taue";
            else if(diHiggsEvt.type == DiHiggsEvent::MU) sN = smpName + "_mu";
            else sN = smpName + "_e";

            plotter.getOrMake1DPre(smpName,"lepPT",";pt",100,0,1000)->Fill(diHiggsEvt.w1_d1->pt(),weight );
            pTS = diHiggsEvt.hww->pt();
            if(PhysicsUtilities::deltaR2(*diHiggsEvt.w1_d1,*diHiggsEvt.w2_d1) < PhysicsUtilities::deltaR2(*diHiggsEvt.w1_d1,*diHiggsEvt.w2_d2)){
                pQ = diHiggsEvt.w2_d1->p4();
                pQ2 = diHiggsEvt.w2_d2->p4();
            } else {
                pQ = diHiggsEvt.w2_d2->p4();
                pQ2 = diHiggsEvt.w2_d1->p4();
            }
            doLepton(sN, diHiggsEvt.w1_d1, muons,electrons,jets);
            if(diHiggsEvt.w1_d1->pt() < 50)doLepton(sN+"_lt50", diHiggsEvt.w1_d1, muons,electrons,jets);
            else if(diHiggsEvt.w1_d1->pt() < 100)doLepton(sN+"_50to100", diHiggsEvt.w1_d1, muons,electrons,jets);
            else if(diHiggsEvt.w1_d1->pt() < 200)doLepton(sN+"_100to200", diHiggsEvt.w1_d1, muons,electrons,jets);
            else if(diHiggsEvt.w1_d1->pt() < 500)doLepton(sN+"_200to500", diHiggsEvt.w1_d1, muons,electrons,jets);
            else doLepton(sN+"_gt500", diHiggsEvt.w1_d1, muons,electrons,jets);
        } else if(!isRealData()){
            SMDecayEvent smEvent; smEvent.setDecayInfo(reader_genpart->genParticles);
            for(const auto& td : smEvent.topDecays){
                if(td.type >= TopDecay::MU && goodGenLepton(td.W_decay.dau1)  ){
                    TString sN = smpName +  ( td.W_decay.dau1->absPdgId() == ParticleInfo::p_eminus ? "_e" : "_mu"  );
                    pTS = td.top->pt();
                    pQ = MomentumF(td.b->p4());
                    pQ2 = MomentumF();
                    doLepton(sN, td.W_decay.dau1, muons,electrons,jets);

                    if(td.W_decay.dau1->pt() < 50)      doLepton(sN+"_lt50", td.W_decay.dau1, muons,electrons,jets);
                    else if(td.W_decay.dau1->pt() < 100)doLepton(sN+"_50to100", td.W_decay.dau1, muons,electrons,jets);
                    else if(td.W_decay.dau1->pt() < 200)doLepton(sN+"_100to200", td.W_decay.dau1, muons,electrons,jets);
                    else if(td.W_decay.dau1->pt() < 500)doLepton(sN+"_200to500", td.W_decay.dau1, muons,electrons,jets);
                    else doLepton(sN+"_gt500", td.W_decay.dau1, muons,electrons,jets);
                }
            }
            plotter.getOrMake1DPre(smpName,"nTops",";# of tops",10,-0.5,9.5)->Fill( smEvent.topDecays.size(),weight );
            plotter.getOrMake1DPre(smpName,"nBosons",";# of bosons",10,-0.5,9.5)->Fill( smEvent.bosonDecays.size(),weight );
            if(smEvent.topDecays.size() == 2){
                int type = 0; //BAD //HAD //DiLEP// W_TAU_HAD,W_TAU_MU,W_TAU_E,W_MU,W_E

                if(smEvent.topDecays[0].type == TopDecay::BAD || smEvent.topDecays[1].type == TopDecay::BAD) type = 0;
                else if(smEvent.topDecays[0].type == TopDecay::HAD && smEvent.topDecays[1].type == TopDecay::HAD) type = 1;
                else if(smEvent.topDecays[0].type != TopDecay::HAD && smEvent.topDecays[1].type != TopDecay::HAD) type = 2;
                else if(smEvent.topDecays[0].type == TopDecay::TAU_HAD || smEvent.topDecays[1].type == TopDecay::TAU_HAD) type = 3;
                else if(smEvent.topDecays[0].type == TopDecay::TAU_MU || smEvent.topDecays[1].type == TopDecay::TAU_MU) type = 4;
                else if(smEvent.topDecays[0].type == TopDecay::TAU_E || smEvent.topDecays[1].type == TopDecay::TAU_E) type = 5;
                else if(smEvent.topDecays[0].type == TopDecay::MU || smEvent.topDecays[1].type == TopDecay::MU) type = 6;
                else if(smEvent.topDecays[0].type == TopDecay::E || smEvent.topDecays[1].type == TopDecay::E) type = 7;

                plotter.getOrMake1DPre(smpName,"diTop_type",";BAD//HAD//DILEP//TAU_HAD//TAU_MU//TAU_E//MU//E",10,-0.5,9.5)->Fill(type,weight );

            }


        }
        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<LeptonProcessor>     leptonProcNoISO ;
    float pTS = 0;
    MomentumF pQ;
    MomentumF pQ2;

};

#endif

void checkLepJetVars(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void checkLepJetVars(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000);
    a.write(outFileName);
}
