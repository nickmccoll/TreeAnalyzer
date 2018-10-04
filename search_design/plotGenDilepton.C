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

#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/EventWeights.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"


#include "TSystem.h"
using namespace TAna;

// Macro Description:

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
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

        const MomentumF dineutrino = nu1->p4() - nu2->p4();
        double met = dineutrino.pt();

        // fill plots
        plotter.getOrMake1DPre(sn, "Min_lep_pt", ";P_{T} (GeV)",50,0,600)->Fill(lepminpt, weight);
        plotter.getOrMake1DPre(sn, "Max_lep_pt", ";P_{T} (GeV)",50,0,600)->Fill(lepmaxpt, weight);
        plotter.getOrMake1DPre(sn, "Min_nu_pt", ";P_{T} (GeV)",50,0,600)->Fill(numinpt, weight);
        plotter.getOrMake1DPre(sn, "Max_nu_pt", ";P_{T} (GeV)",50,0,600)->Fill(numaxpt, weight);
        plotter.getOrMake1DPre(sn, "dR_ll", ";#DeltaR",50,0,5)->Fill(dR_ll, weight);
        plotter.getOrMake1DPre(sn, "dPhi_ll", ";#Delta#phi",50,-3.2,3.2)->Fill(dPhi_ll, weight);
        plotter.getOrMake1DPre(sn, "dEta_ll", ";#Delta#eta",50,-5,5)->Fill(dEta_ll, weight);
        plotter.getOrMake1DPre(sn, "genmet",";E_{T}^{miss} (GeV)",50,0,600)->Fill(met, weight);
    }

    bool goodGenLepton(const GenParticle* genP) {
        if(genP->absEta() > 2.4) return false;
        if(genP->pt() < ( genP->absPdgId() == ParticleInfo::p_eminus ? 33.0 : 29  )  ) return false;
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
        else if (lep1->absPdgId() == 15 && lep2->absPdgId() == 15) sn += "tautau";
        else if ((lep1->absPdgId() == 11 && lep2->absPdgId() == 13) || (lep1->absPdgId() == 13 && lep2->absPdgId() == 11)) sn += "emu";
        else if ((lep1->absPdgId() == 11 && lep2->absPdgId() == 15) || (lep1->absPdgId() == 15 && lep2->absPdgId() == 11)) sn += "etau";
        else if ((lep1->absPdgId() == 13 && lep2->absPdgId() == 15) || (lep1->absPdgId() == 15 && lep2->absPdgId() == 13)) sn += "mutau";
        else {
        	std::cout << "Error encountered: dilep not found" << std::endl;
        	std::cout << "decaytype = " << diHiggsEvt.type << std::endl;
        	std::cout << "lep1 = " << lep1->absPdgId() << std::endl;
        	std::cout << "lep2 = " << lep2->absPdgId() << std::endl;
        	std::cout << "diHiggsEvt.w1_d1 = " << diHiggsEvt.w1_d1->absPdgId() << std::endl;
        	std::cout << "diHiggsEvt.w1_d2 = " << diHiggsEvt.w1_d2->absPdgId() << std::endl;
        	std::cout << "diHiggsEvt.w2_d1 = " << diHiggsEvt.w2_d1->absPdgId() << std::endl;
        	std::cout << "diHiggsEvt.w2_d2 = " << diHiggsEvt.w2_d2->absPdgId() << std::endl;
        	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
        }
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
        if (lep2->pt() < 20) return false;
        plotSpectra(smpName+"_incl_baseline",lep1,lep2,nu1,nu2);
        plotSpectra(sn+"_baseline",lep1,lep2,nu1,nu2);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void plotGenDilepton(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void plotGenDilepton(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);

}
