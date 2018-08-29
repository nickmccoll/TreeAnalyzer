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
        double lepminpt = lep1->pt() > lep2->pt() ? lep2->pt() : lep1->pt();
        double lepmaxpt = lep1->pt() > lep2->pt() ? lep1->pt() : lep2->pt();

        double numinpt = nu1->pt() > nu2->pt() ? nu2->pt() : nu1->pt();
        double numaxpt = nu1->pt() > nu2->pt() ? nu1->pt() : nu2->pt();

        double dR_ll = PhysicsUtilities::deltaR(lep1, lep2);
        double dPhi_ll = PhysicsUtilities::deltaPhi(lep1, lep2);
        double dEta_ll = PhysicsUtilities::deltaEta(lep1, lep2);

        const MomentumF dineutrino = nu1->p4() - nu2->p4();
        double met = dineutrino.pt();

        // fill plots
        plotter.getOrMake1DPre(sn, "Min_lep_pt", ";P_{T}",50,0,500)->Fill(lepminpt, weight);
        plotter.getOrMake1DPre(sn, "Max_lep_pt", ";P_{T}",50,0,500)->Fill(lepmaxpt, weight);
        plotter.getOrMake1DPre(sn, "Min_nu_pt", ";P_{T}",50,0,500)->Fill(numinpt, weight);
        plotter.getOrMake1DPre(sn, "Min_nu_pt", ";P_{T}",50,0,500)->Fill(numinpt, weight);
        plotter.getOrMake1DPre(sn, "dR_ll", ";#deltaR",50,0,5)->Fill(dR_ll, weight);
        plotter.getOrMake1DPre(sn, "dPhi_ll", ";#delta#phi",50,-3.2,3.2)->Fill(dPhi_ll, weight);
        plotter.getOrMake1DPre(sn, "dEta_ll", ";#delta#eta",50,-5,5)->Fill(dEta_ll, weight);
        plotter.getOrMake1DPre(sn, "genmet",";E_T^{miss}",50,0,500)->Fill(met, weight);
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process >= FillerConstants::ZJETS && reader_event->process <= FillerConstants::TTX )
            smpName = "other";
        TString sn = smpName + "_";

        // take dilepton events
        if(diHiggsEvt.type != DiHiggsEvent::DILEP) return false;

        // keep track of ee, mumu, and emu events
        const GenParticle *lep1 = diHiggsEvt.w1_d1;
        if (lep1->absPdgId() == 11) sn += "e";
        else if (lep1->absPdgId() == 13) sn += "mu";
        else std::cout << "lep1 neither e nor mu" << std::endl;

        const GenParticle *lep2 = diHiggsEvt.w2_d1;
        if (lep2->absPdgId() == 11) sn += "e";
        else if (lep2->absPdgId() == 13) sn += "mu";
        else std::cout << "lep2 neither e nor mu" << std::endl;
        sn += "_";

        const GenParticle *nu1 = diHiggsEvt.w1_d2;
        const GenParticle *nu2 = diHiggsEvt.w2_d2;

        plotSpectra(smpName,lep1,lep2,nu1,nu2);

        // cut flow

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
