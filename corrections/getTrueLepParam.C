
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
        reader_jetwlep =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jetwlep );;
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

        plotter.getOrMake1DPre(prefix,"ptj_o_ptl",";p_{T}(jet)/p_{T}(lep)",100,0,10)->Fill( nearestJet ? nearestJet->pt() / leptonMom->pt() : 0,weight  );
        const float projLep = nearestJet ? nearestJet->pt()*TMath::Cos( PhysicsUtilities::deltaPhi(*leptonMom,*nearestJet)) : 0.0;
        plotter.getOrMake1DPre(prefix,"ptjpl_o_ptl",";p_{T}(jet proj. on lep)/p_{T}(lep)",100,0,10)->Fill( projLep > 0 ? projLep / leptonMom->pt() : 0,weight );

        const MomentumF jetMLep(nearestJet ? nearestJet->p4() - leptonMom->p4() : ASTypes::CylLorentzVectorF() );
        plotter.getOrMake1DPre(prefix,"ptjml_o_ptl",";p_{T}(jet - lep)/p_{T}(lep)",100,0,10)->Fill( projLep > leptonMom->pt()  ? jetMLep.pt()   / leptonMom->pt() : 0,weight );

        plotter.getOrMake1DPre(prefix,"jml_dr_l",";#DeltaR(jet - lep,lep)",100,0,1)->Fill( projLep > leptonMom->pt()  ? PhysicsUtilities::deltaR(jetMLep, *leptonMom) : 10,weight );


        const float projJetMLep = nearestJet ? jetMLep.pt()*TMath::Cos( PhysicsUtilities::deltaPhi(*leptonMom,jetMLep)) : 0.0;
        plotter.getOrMake1DPre(prefix,"ptjmlpl_o_ptl",";p_{T}(jet - lep)/p_{T}(lep)",100,0,10)->Fill( projLep > leptonMom->pt()  ? projJetMLep / leptonMom->pt() : 0,weight );


        plotter.getOrMake1DPre(prefix,"pts",";p_{T}(s)",300,0,3000)->Fill(pTS,weight );
        plotter.getOrMake1DPre(prefix,"ptl_o_pts",";p_{T}(l)/p_{T}(s)",100,0,2)->Fill(pTS > 0 ? leptonMom->pt() / pTS : 0,weight );

        plotter.getOrMake2DPre(prefix,"ptl_o_pt1q_v_dr",";p_{T}(l)/p_{T}(l+q) ; dr", 100,0,1, 80,0,0.4)->Fill(leptonMom->pt()/(leptonMom->p4() + pQ.p4()).pt(),PhysicsUtilities::deltaR(leptonMom->p4(),pQ.p4()) ,weight );
        plotter.getOrMake2DPre(prefix,"ptq_o_ptl_v_dr",";p_{T}(q)/p_{T}(l) ; dr", 100,0,10, 80,0,0.4)->Fill(pQ.pt()/(leptonMom->pt()),PhysicsUtilities::deltaR(leptonMom->p4(),pQ.p4()) ,weight );
        double rV  = std::max(0.05,std::min(0.2, 10.0/leptonMom->pt()));
        plotter.getOrMake2DPre(prefix,"ptq_o_ptl_v_drN",";p_{T}(q)/p_{T}(l) ; drN", 100,0,10, 100,0,2)->Fill(pQ.pt()/(leptonMom->pt()),PhysicsUtilities::deltaR(leptonMom->p4(),pQ.p4())/rV ,weight );

        MomentumF totQ = pQ.p4() + pQ2.p4();
        plotter.getOrMake2DPre(prefix,"ptqt_o_ptl_v_dr",";p_{T}(qt)/p_{T}(l) ; dr", 100,0,10, 80,0,0.4)->Fill(totQ.pt()/(leptonMom->pt()),PhysicsUtilities::deltaR(leptonMom->p4(),totQ.p4()) ,weight );
        plotter.getOrMake2DPre(prefix,"ptqt_o_ptl_v_drN",";p_{T}(qt)/p_{T}(l) ; drN", 100,0,10, 100,0,2)->Fill(totQ.pt()/(leptonMom->pt()),PhysicsUtilities::deltaR(leptonMom->p4(),totQ.p4())/rV ,weight );


            plotter.getOrMake1DPre(prefix,"jml_drN_l",";#DeltaR(jet - lep,lep)",200,0,4)->Fill( projLep > leptonMom->pt()  ? PhysicsUtilities::deltaR(jetMLep, *leptonMom)/rV : 10,weight );

            if(jetMLep.pt() > 0.5*leptonMom->pt()){
                plotter.getOrMake1DPre(prefix,"pt0p5_jml_drN_l",";#DeltaR(jet - lep,lep)",200,0,4)->Fill( projLep > leptonMom->pt()  ? PhysicsUtilities::deltaR(jetMLep, *leptonMom)/rV : 10,weight );
                plotter.getOrMake1DPre(prefix,"pt0p5_jml_dr_l",";#DeltaR(jet - lep,lep)",100,0,1)->Fill( projLep > leptonMom->pt()  ? PhysicsUtilities::deltaR(jetMLep, *leptonMom) : 10,weight );
            }


            if(projLep > leptonMom->pt())
            plotter.getOrMake2DPre(prefix,"jml_o_ptl_v_dr",";p_{T}(jml)/p_{T}(l) ; dr", 100,0,10, 80,0,0.4)->Fill( jetMLep.pt()   / leptonMom->pt() ,PhysicsUtilities::deltaR(leptonMom->p4(),jetMLep.p4()) ,weight );
            if(projLep > leptonMom->pt())
            plotter.getOrMake2DPre(prefix,"jml_o_ptl_v_drN",";p_{T}(jml)/p_{T}(l) ; dr", 100,0,10, 100,0,2)->Fill( jetMLep.pt()   / leptonMom->pt(),PhysicsUtilities::deltaR(leptonMom->p4(),jetMLep.p4())/rV ,weight );


//        plotter.getOrMake1DPre(prefix,"eiso_004_o_ptl",";p_{T}(jet - lep)/p_{T}(lep)",100,0,10)->Fill( projLep > leptonMom->pt()  ? jetMLep.pt()*TMath::Exp(-PhysicsUtilities::deltaR2(jetMLep, *leptonMom)/0.004) / leptonMom->pt() : 0,weight );
//        plotter.getOrMake1DPre(prefix,"eiso_04_o_ptl",";p_{T}(jet - lep)/p_{T}(lep)",100,0,10)->Fill( projLep > leptonMom->pt()  ? jetMLep.pt()*TMath::Exp(-PhysicsUtilities::deltaR2(jetMLep, *leptonMom)/0.04) / leptonMom->pt() : 0,weight );
//        plotter.getOrMake1DPre(prefix,"eiso_4_o_ptl",";p_{T}(jet - lep)/p_{T}(lep)",100,0,10)->Fill( projLep > leptonMom->pt()  ? jetMLep.pt()*TMath::Exp(-PhysicsUtilities::deltaR2(jetMLep, *leptonMom)/0.4) / leptonMom->pt() : 0,weight );


//                    plotter.getOrMake2DPre(prefix,"ptjmlpl_o_ptl_v_ptl",";ptl ;ptjmlpl_o_ptl",50,0,1000,50,0,10)->Fill( leptonMom->pt(), projLep > leptonMom->pt()  ? projJetMLep / leptonMom->pt() : 0,weight );
//                    plotter.getOrMake2DPre(prefix,"ptjmlpl_v_ptl",";ptl ;ptjmlpl",50,0,1000,50,0,1000)->Fill( leptonMom->pt(), projLep > leptonMom->pt()  ? projJetMLep : 0 ,weight );
//                    plotter.getOrMake2DPre(prefix,"ptj_v_ptl",";ptl ;ptj",50,0,1000,50,0,1000)->Fill( leptonMom->pt(), nearestJet ? nearestJet->pt()   : 0 ,weight );


//        if(projLep > leptonMom->pt() && ISO > 0 ){
//            float val = -2*TMath::Log(leptonMom->pt()*ISO/jetMLep.pt());
//            plotter.getOrMake2DPre(prefix,"log_v_r2",";r2 ;log",200,0,.2,100,-20,20)->Fill( PhysicsUtilities::deltaR2(jetMLep, *leptonMom),val ,weight );
//        }
    }

    void doLepton(const TString& prefix, const GenParticle* genP,const std::vector<const Muon*>& muons,const std::vector<const Electron*>& electrons,
            const std::vector<const Jet*>& jets){

        const auto* recoL = getMatchedLepton(*genP,muons,electrons);
        double nearDRToGen = 20;
        const auto* jetNearGen = getMatchedJet(genP,jets,nearDRToGen);
        double nearDRToReco = 20;
        const auto* jetNearReco = recoL ? getMatchedJet(recoL,jets,nearDRToReco) : 0;

        //reconstruction eff
//        makePlots(prefix +"_gen_dr0p6_incl",genP,nearDRToGen < 0.6 ? jetNearGen : 0);
        makePlots(prefix +"_gen_dr0p4_incl",genP,nearDRToGen < 0.4 ? jetNearGen : 0);
//        makePlots(prefix +"_gen_dr0p2_incl",genP,nearDRToGen < 0.2 ? jetNearGen : 0);
        if(recoL == 0) return;
//        makePlots(prefix +"_gen_dr0p6_reco",genP,nearDRToGen < 0.6 ? jetNearGen : 0);
        makePlots(prefix +"_gen_dr0p4_reco",genP,nearDRToGen < 0.4 ? jetNearGen : 0);
//        makePlots(prefix +"_gen_dr0p2_reco",genP,nearDRToGen < 0.2 ? jetNearGen : 0);
        //ID eff
//        makePlots(prefix +"_reco_dr0p6_incl",recoL,nearDRToReco < 0.6 ? jetNearReco : 0);
        makePlots(prefix +"_reco_dr0p4_incl",recoL,nearDRToReco < 0.4 ? jetNearReco : 0);
//        makePlots(prefix +"_reco_dr0p2_incl",recoL,nearDRToReco < 0.2 ? jetNearReco : 0);
        const bool passID = leptonProcNoISO->isGoodLepton(*reader_event,recoL);
        if(!passID) return;
//        makePlots(prefix +"_reco_dr0p6_ID",recoL,nearDRToReco < 0.6 ? jetNearReco : 0,recoL->miniIso());
        makePlots(prefix +"_reco_dr0p4_ID",recoL,nearDRToReco < 0.4 ? jetNearReco : 0,recoL->miniIso());
//        makePlots(prefix +"_reco_dr0p2_ID",recoL,nearDRToReco < 0.2 ? jetNearReco : 0,recoL->miniIso());
        //ISO eff
        const bool passISO = leptonProc->isGoodLepton(*reader_event,recoL);
        if(!passISO) return;
//        makePlots(prefix +"_reco_dr0p6_ISO",recoL,nearDRToReco < 0.6 ? jetNearReco : 0);
        makePlots(prefix +"_reco_dr0p4_ISO",recoL,nearDRToReco < 0.4 ? jetNearReco : 0);
//        makePlots(prefix +"_reco_dr0p2_ISO",recoL,nearDRToReco < 0.2 ? jetNearReco : 0);
    };

    bool goodGenLepton(const GenParticle* genP) {
        if(genP->absEta() > 2.4) return false;
        if(genP->pt() < ( genP->absPdgId() == ParticleInfo::p_eminus ? 33.0 : 29  )  ) return false;
        return true;
    }


    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;

        const std::vector<const Muon    *> muons     = JetKinematics::selectObjectsConst(reader_muon->muons,26,2.4);
        const std::vector<const Electron*> electrons = JetKinematics::selectObjectsConst(reader_electron->electrons,30,2.4);
        const std::vector<const Jet     *> jets      = JetKinematics::selectObjectsConst(reader_jetwlep->jets,20,10);

        if(reader_event->process == FillerConstants::SIGNAL){
            if(diHiggsEvt.type < DiHiggsEvent::MU) return false;
            if(!goodGenLepton(diHiggsEvt.w1_d1)) return false;
            TString sN = smpName +  ( diHiggsEvt.w1_d1->absPdgId() == ParticleInfo::p_eminus ? "_e" : "_mu"  );
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

void getTrueLepParam(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getTrueLepParam(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze(10000);
    a.write(outFileName);
}
