#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "TreeReaders/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


#include "TSystem.h"
using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    }

    void plotSpectra(TString sn, const GenParticle *lep1, const GenParticle *lep2, const GenParticle *nu1, const GenParticle *nu2) {

        // some variables for plotting
        double lepminpt = lep2->pt();
        double lepmaxpt = lep1->pt();

        double numinpt = nu2->pt();
        double numaxpt = nu1->pt();

        double dR_ll = PhysicsUtilities::deltaR(*lep1, *lep2);
        double dPhi_ll = PhysicsUtilities::deltaPhi(*lep1, *lep2);
        double dEta_ll = PhysicsUtilities::deltaEta(*lep1, *lep2);

        const MomentumF dilepton = lep1->p4() + lep2->p4();
        double dR_bbll = PhysicsUtilities::deltaR(*diHiggsEvt.hbb,dilepton);
        double dphi_bbll = PhysicsUtilities::deltaPhi(*diHiggsEvt.hbb,dilepton);

        const MomentumF dineutrino = nu1->p4() + nu2->p4();
        double met = dineutrino.pt();

        // fill plots
        plotter.getOrMake1DPre(sn, "Min_lep_pt", ";P_{T} (GeV)",100,0,1000)->Fill(lepminpt, weight);
        plotter.getOrMake1DPre(sn, "Max_lep_pt", ";P_{T} (GeV)",100,0,1000)->Fill(lepmaxpt, weight);
        plotter.getOrMake1DPre(sn, "Min_nu_pt", ";P_{T} (GeV)",100,0,1000)->Fill(numinpt, weight);
        plotter.getOrMake1DPre(sn, "Max_nu_pt", ";P_{T} (GeV)",100,0,1000)->Fill(numaxpt, weight);
        plotter.getOrMake1DPre(sn, "dR_ll", ";#DeltaR",50,0,5)->Fill(dR_ll, weight);
        plotter.getOrMake1DPre(sn, "dPhi_ll", ";|#Delta#phi|",50,0,3.14)->Fill(abs(dPhi_ll), weight);
        plotter.getOrMake1DPre(sn, "dEta_ll", ";|#Delta#eta|",50,0,5)->Fill(abs(dEta_ll), weight);
        plotter.getOrMake1DPre(sn, "genmet",";E_{T}^{miss} (GeV)",120,0,1200)->Fill(met, weight);

        plotter.getOrMake1DPre(sn, "Lep_pt", ";P_{T} (GeV)",100,0,1000)->Fill(lepminpt, weight);
        plotter.getOrMake1DPre(sn, "Lep_pt", ";P_{T} (GeV)",100,0,1000)->Fill(lepmaxpt, weight);
        plotter.getOrMake1DPre(sn, "dR_bbll", ";#DeltaR",100,0,7)->Fill(dR_bbll, weight);
        plotter.getOrMake1DPre(sn, "dphi_bbll", ";|#Delta#phi|",50,0,3.14)->Fill(abs(dphi_bbll), weight);
    }

    bool goodGenLepton(const GenParticle* genP) {
        if(genP->absEta() > (genP->absPdgId() == ParticleInfo::p_eminus ? 2.5 : 2.4)) return false;
        if(genP->pt() < (genP->absPdgId() == ParticleInfo::p_eminus ? 30.0 : 26.0)) return false;
        return true;
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process >= FillerConstants::ZJETS && reader_event->process <= FillerConstants::TTX )
            smpName = "other";
        TString sn = smpName + "_";

        // take only dilepton events
        if(diHiggsEvt.type != DiHiggsEvent::DILEP) return false;

        // order the variable names by pt
        const GenParticle *lep1 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w1_d1 : diHiggsEvt.w2_d1);
        const GenParticle *lep2 = (diHiggsEvt.w1_d1->pt() > diHiggsEvt.w2_d1->pt() ? diHiggsEvt.w2_d1 : diHiggsEvt.w1_d1);
        const GenParticle *nu1 = (diHiggsEvt.w1_d2->pt() > diHiggsEvt.w2_d2->pt() ? diHiggsEvt.w1_d2 : diHiggsEvt.w2_d2);
        const GenParticle *nu2 = (diHiggsEvt.w1_d2->pt() > diHiggsEvt.w2_d2->pt() ? diHiggsEvt.w2_d2 : diHiggsEvt.w1_d2);

        // keep track of all combos of dilepton events for e, mu, tau
        if (lep1->absPdgId() == 11 && lep2->absPdgId() == 11) sn += "ee";
        else if (lep1->absPdgId() == 13 && lep2->absPdgId() == 13) sn += "mumu";
        else if ((lep1->absPdgId() == 11 && lep2->absPdgId() == 13) || (lep1->absPdgId() == 13 && lep2->absPdgId() == 11)) sn += "emu";
        else return false;

        plotSpectra(smpName+"_incl",lep1,lep2,nu1,nu2);
        plotSpectra(sn,lep1,lep2,nu1,nu2);

        // event filter and trigger preselection for all dilepton events
        if (passEventFilters && passTriggerPreselection) plotSpectra(smpName+"_incl_trig_evtfltr",lep1,lep2,nu1,nu2);

        // ht selection
        if (ht_chs < 400) return false;
        plotSpectra(smpName+"_incl_ht400",lep1,lep2,nu1,nu2);
        plotSpectra(sn+"_ht400",lep1,lep2,nu1,nu2);

        // ensure leading lepton is good GenLepton
        if (!goodGenLepton(lep1)) return false;
        plotSpectra(smpName+"_incl_ht400_goodlep",lep1,lep2,nu1,nu2);
        plotSpectra(sn+"_ht400_goodlep",lep1,lep2,nu1,nu2);

        // for orthogonality with single-lepton channel, require second lepton to have pt>20
        if (lep2->pt() < 10) return false;
        plotSpectra(smpName+"_incl_baseline",lep1,lep2,nu1,nu2);
        plotSpectra(sn+"_baseline",lep1,lep2,nu1,nu2);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void plotGenDilepton(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.analyze();
    a.write(outFileName);
}
void plotGenDilepton(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt, randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
