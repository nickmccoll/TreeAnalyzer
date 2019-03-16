
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "Configuration/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "Processors/Variables/interface/FatJetSelection.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/HiggsSolver.h"


#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
    }

    MomentumF getInvisible2(const MomentumF& met, const MomentumF& vis, const double hMass, const MomentumF& lep){
        const double a = hMass*hMass - vis.mass()*vis.mass() +2*vis.x()*met.x() +2*vis.y()*met.y();
        const double A = 4*(vis.E()*vis.E() - vis.z()*vis.z());
        const double B = -4*a* vis.z();
        const double C = 4*vis.E()*vis.E()*(met.x()*met.x() + met.y()*met.y()) - a*a;
        const double delta = B*B -4*A*C;

        double metZ = 0;
        if(delta < 0) {
            metZ= -B/(2*A);
        } else {
            const double pos = (-B + std::sqrt(delta))/(2*A);
            const double neg = (-B - std::sqrt(delta))/(2*A);

            ASTypes::CartLorentzVector neutP(met.x(),met.y(),pos,std::sqrt(met.x()*met.x()+met.y()*met.y()+pos*pos));
            ASTypes::CartLorentzVector neutN(met.x(),met.y(),neg,std::sqrt(met.x()*met.x()+met.y()*met.y()+neg*neg));

            if((neutP + lep.p4()).mass() < (neutN + lep.p4()).mass())metZ = pos;
            else metZ = neg;
        }
        ASTypes::CartLorentzVector neutrino(met.x(),met.y(),metZ,std::sqrt(met.x()*met.x()+met.y()*met.y()+metZ*metZ));
        return MomentumF(neutrino);
    }

    void SRPlots(TString prefix) {
        auto wlnu = neutrino.p4() + selectedLepton->p4();
        auto hWW = wlnu + wjjCand->p4();
        const double hWWDR = PhysicsUtilities::deltaR(wlnu,*wjjCand);

        auto visHWW = selectedLepton->p4() + wjjCand->p4();
        auto visHWWSD = selectedLepton->p4() + wjjCand->sdMom().p4();

        auto semVtransverseMass = [](const ASTypes::CylLorentzVectorF& visible, const ASTypes::CylLorentzVectorF& invisible)->float
        {
          const double    cosDPhi   = TMath::Cos( PhysicsUtilities::deltaPhi(visible.phi(), invisible.phi()) );
          return TMath::Sqrt( visible.mass()*visible.mass() +  2 * invisible.pt() *(TMath::Sqrt(visible.mass()*visible.mass() + visible.pt()*visible.pt()) - visible.pt()*cosDPhi) );
        };

        float hWWTM = JetKinematics::transverseMass(visHWW,reader_event->met);
        float hWWTM2 = semVtransverseMass(visHWW,reader_event->met.p4());
        float WTM = JetKinematics::transverseMass(*wjjCand,reader_event->met);


        const double wlnuNuDR = PhysicsUtilities::deltaR(neutrino,wlnu);
        const double wlnuNuDPhi = PhysicsUtilities::absDeltaPhi(neutrino,wlnu);
        const double hWWNuDR = PhysicsUtilities::deltaR(neutrino,hWW);
        const double hWWNuDPhi = PhysicsUtilities::absDeltaPhi(neutrino,hWW);

        const double lepNuDR = PhysicsUtilities::deltaR(neutrino,*selectedLepton);
        const double lepNuDPhi = PhysicsUtilities::absDeltaPhi(neutrino,*selectedLepton);
        const double nuHBBDPhi = PhysicsUtilities::absDeltaPhi(reader_event->met,*hbbCand);
        const double nuHBBDR = PhysicsUtilities::deltaR(neutrino,*hbbCand);

        ASTypes::CartLorentzVector vM(reader_event->met.x(),reader_event->met.y(),0,std::sqrt(reader_event->met.x()*reader_event->met.x()+reader_event->met.y()*reader_event->met.y()));
        auto plt = [&](const TString& pre){
            plotter.getOrMake1DPre(pre, "WTM_mass"   ,";W transverse mass [GeV]; arbitrary units",50,0,500)->Fill(WTM,weight);
            plotter.getOrMake1DPre(pre, "hWWTM2_mass"   ,";hWW transverse mass [GeV]; arbitrary units",50,0,500)->Fill(hWWTM2,weight);
            plotter.getOrMake1DPre(pre, "hWWTM_mass"   ,";hWW transverse mass [GeV]; arbitrary units",50,0,500)->Fill(hWWTM,weight);
            if(hWWDR < 0.5)
                plotter.getOrMake1DPre(pre, "drm_hWWv_mass"   ,";HWW vis mass [GeV]; arbitrary units",50,0,500)->Fill(visHWW.mass(),weight);
            plotter.getOrMake1DPre(pre, "hWWv_mass"   ,";HWW vis mass [GeV]; arbitrary units",50,0,500)->Fill(visHWW.mass(),weight);
            plotter.getOrMake1DPre(pre, "hWWvsd_mass"   ,";HWW vis mass [GeV]; arbitrary units",50,0,500)->Fill(visHWWSD.mass(),weight);
            plotter.getOrMake1DPre(pre, "wjj_wlnu_dR" ,";#DeltaR(Wjj,Wl#nu); arbitrary units",128,0,6.4)->Fill(hWWDR,weight);
            plotter.getOrMake1DPre(pre, "wjj_wlnu_dRM" ,";#DeltaR(Wjj,Wl#nu)*p_{T}(HWW)/125; arbitrary units",100,0,10)->Fill(hWWDR*hWW.pt()/125.0,weight);

            plotter.getOrMake1DPre(pre, "lep_nu_dRM" ,";#DeltaR(#nu,lep)*p_{T}(W)/80; arbitrary units",100,0,10)->Fill(lepNuDR*wlnu.pt()/80.0,weight);
            plotter.getOrMake1DPre(pre, "wlnu_pt_o_hWW_pt" ,";p_{T} W / p_{T} HWW; arbitrary units",100,0,2)->Fill(wlnu.pt()/hWW.pt(),weight);
            plotter.getOrMake1DPre(pre, "wlnu_pt_o_wjj_pt" ,";p_{T} Wl#nu / p_{T} Wjj; arbitrary units",500,0,5)->Fill(wlnu.pt()/wjjCand->pt(),weight);


            plotter.getOrMake1DPre(pre, "wlnu_nu_dR" ,";#DeltaR(#nu,Wlnu); arbitrary units",128,0,6.4)->Fill(wlnuNuDR,weight);
            plotter.getOrMake1DPre(pre, "wlnu_nu_dPhi" ,";#DeltaPhi(#nu,Wlnu); arbitrary units",64,0,3.2)->Fill(wlnuNuDPhi,weight);
            plotter.getOrMake1DPre(pre, "hWW_nu_dR" ,";#DeltaR(#nu,HWW); arbitrary units",128,0,6.4)->Fill(hWWNuDR,weight);
            plotter.getOrMake1DPre(pre, "hWW_nu_dPhi" ,";#DeltaPhi(#nu,HWW); arbitrary units",64,0,3.2)->Fill(hWWNuDPhi,weight);

            plotter.getOrMake1DPre(pre, "lep_nu_dR" ,";#DeltaR(#nu,lep); arbitrary units",128,0,6.4)->Fill(lepNuDR,weight);
            plotter.getOrMake1DPre(pre, "lep_nu_dPhi" ,";|#Delta#phi(#nu,lep)|; arbitrary units",64,0,3.2)->Fill(lepNuDPhi,weight);
            plotter.getOrMake1DPre(pre, "hbb_nu_dPhi" ,";|#Delta#phi(#nu,Hbb)|; arbitrary units",64,0,3.2)->Fill(nuHBBDPhi,weight);
            plotter.getOrMake1DPre(pre, "hbb_nu_dR" ,";#DeltaR(#nu,Hbb); arbitrary units",128,0,6.4)->Fill(nuHBBDR,weight);
            plotter.getOrMake1DPre(pre, "wlnu_pt" ,";p_{T} W; arbitrary units",500,0,1000)->Fill(wlnu.pt(),weight);
        };

        plt(prefix+"_hhInc");
        if(hh.mass() >= 700 && hh.mass() < 900) plt(prefix+"_hh700to900");
        if(hh.mass() >= 900 && hh.mass() < 1100) plt(prefix+"_hh900to1100");
        if(hh.mass() >= 1400 && hh.mass() < 1800) plt(prefix+"_hh1400to1800");
        if(hh.mass() >= 2500 && hh.mass() < 3500) plt(prefix+"_hh2500to3500");
    }

    void radiusPlots(TString prefix) {

        auto doOneDecay=[&](const TString& pre,const ASTypes::CylLorentzVectorF& a,const ASTypes::CylLorentzVectorF& b, const ASTypes::CylLorentzVectorF& tot, const double deltaR ){
            plotter.getOrMake1DPre(pre, "pt" ,";p_{T}(a+b); arbitrary units",400,0,2000)->Fill(tot.pt(),weight);
            plotter.getOrMake1DPre(pre, "dR" ,";#DeltaR(a,b); arbitrary units",64,0,3.2)->Fill(deltaR,weight);
            plotter.getOrMake1DPre(pre, "dPhi" ,";#Delta#phi(a,b); arbitrary units",64,0,3.2)->Fill(PhysicsUtilities::absDeltaPhi(a,b),weight);
            plotter.getOrMake1DPre(pre, "minPToTotPT" ,";min[p_{T}]/a+b p_{T}; arbitrary units",120,0,1.2)->Fill(std::min(a.pt(),b.pt())/tot.pt(),weight);
            plotter.getOrMake1DPre(pre, "deltaRPToMH" ,";#DeltaR(a,b)*(a+b p_{T})/M(H); arbitrary units",100,0,10)->Fill(deltaR*tot.pt()/125.0,weight);
            plotter.getOrMake1DPre(pre, "deltaRPToMT" ,";#DeltaR(a,b)*(a+b p_{T})/M(a+b); arbitrary units",100,0,10)->Fill(deltaR*tot.pt()/tot.mass(),weight);

        };


        const MomentumF wlnu = neutrino.p4() + selectedLepton->p4();
        const double WlnuDR = PhysicsUtilities::deltaR(neutrino,*selectedLepton);
        const MomentumF wjjSD = wjjCand->sdMom();
        const MomentumF hWW   =wjjCand->p4() + wlnu.p4();
        const MomentumF hWWSD =wjjSD.p4() + wlnu.p4();
        const double hWWDR = PhysicsUtilities::deltaR(wlnu,*wjjCand);
        const double hWWDRSD = PhysicsUtilities::deltaR(wlnu,wjjSD);


        auto plt = [&](const TString& pre){
            doOneDecay(pre+"_wlnu",neutrino.p4(),selectedLepton->p4(),wlnu.p4(),WlnuDR);
            doOneDecay(pre+"_hWW",wjjCand->p4(),wlnu.p4(),hWW.p4(),hWWDR);
            doOneDecay(pre+"_hWWSD",wjjSD.p4(),wlnu.p4(),hWWSD.p4(),hWWDRSD);
        };

        plt(prefix+"_hWWInc");
        if(hWW.pt() < 500)  plt(prefix+"_hWW_pt_lt500");
        if(hWW.pt() >= 500  && hWW.pt() < 750)  plt(prefix+"_hWW_pt_500to750");
        if(hWW.pt() >= 750 && hWW.pt() < 1000)  plt(prefix+"_hWW_pt_750to1000");
        if(hWW.pt() >= 1000 && hWW.pt() < 1500) plt(prefix+"_hWW_pt_1000to1500");
        if(hWW.pt() >= 1500)                    plt(prefix+"_hWW_pt_geq1500");
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process >= FillerConstants::ZJETS && reader_event->process <= FillerConstants::TTX )
            smpName = "other";

        if(!passEventFilters) return false;
        if(!passTriggerPreselection) return false;
        if(selectedLeptons.size() != 1) return false;
        if(!passWjjSel || !hbbCand) return false;
        if(hbbMass >= 105 && hbbMass < 145){
            SRPlots(smpName+"_hbbI");
            if(passHbbSel)SRPlots(smpName+"_hbbLT");
            if(passHbbSel&&!passHbbTSel)SRPlots(smpName+"_hbbL");
            if(passHbbTSel)SRPlots(smpName+"_hbbT");
        }



        radiusPlots(smpName+"_hbbI");
        if(isSignal())radiusPlots("signal_hbbI");
        else radiusPlots("bkg_hbbI");
        return true;


    }

    void write(TString fileName){ plotter.write(fileName);}

    HistGetter plotter;

};

#endif

void getRecoHHCuts(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getRecoHHCuts(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
