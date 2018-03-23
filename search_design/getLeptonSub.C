
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
#include "TreeReaders/interface/FatJetReader.h"

using namespace TAna;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt){
        leptonNoISOProc  .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonNoISOProc);
        leptonNoISOProc->lepSelParams.mu_maxISO  = -1;
        leptonNoISOProc->lepSelParams.el_maxISO  = -1;
        leptonNoISOProc->lepSelParams_dataABCDEF.mu_maxISO  = -1;
        leptonNoISOProc->lepSelParams_dataABCDEF.el_maxISO  = -1;
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type < DiHiggsEvent::TAU_MU) return false;
        if(selectedLeptons.size() != 1) return false;
        if(!wjjCand) return false;
        if(!hbbCand) return false;
        if(nMedBTags_HbbV) return false;

        auto leps = leptonNoISOProc->getLeptons(*reader_event,*reader_muon,*reader_electron);
        TString prefix = smpName;

        int nTot = 0;
        int nDPhiPio2 = 0;
        int nDR8 = 0;
        int nTotPT2 = 0;
        int nDPhiPio2PT2 = 0;
        int nDR8PT2 = 0;
        double maxOPT=0;
        for(const auto* l : leps){
            if(l->isMuon() == selectedLepton->isMuon() && l->index() == selectedLepton->index()) continue;
            nTot++;
            double dphi = PhysicsUtilities::absDeltaPhi(*l,*wjjCand);
            double dr2 = PhysicsUtilities::deltaR2(*l,*wjjCand);
            bool passPT = l->pt() >= (l->isMuon() ? 26.0 : 30.0);
            if(passPT) nTotPT2++;
            if(dphi < TMath::PiOver2() ) nDPhiPio2++;
            if(dr2 < 0.8*0.8){
                nDR8++;
                if(l->pt() > maxOPT) maxOPT = l->pt();
            }
            if(dphi < TMath::PiOver2() && passPT ) nDPhiPio2PT2++;
            if(dr2 < 0.8*0.8 && passPT){
                nDR8PT2++;
            }
        }

        plotter.getOrMake1DPre(smpName,"all",";nLep",10,-0.5,9.5)->Fill(nTot, weight);
        plotter.getOrMake1DPre(smpName,"dPhiLtPio2",";nLep",10,-0.5,9.5)->Fill(nDPhiPio2, weight);
        plotter.getOrMake1DPre(smpName,"dPhiGtPio2",";nLep",10,-0.5,9.5)->Fill(nTot-nDPhiPio2, weight);
        plotter.getOrMake1DPre(smpName,"drLt0p8",";nLep",10,-0.5,9.5)->Fill(nDR8, weight);

        plotter.getOrMake1DPre(smpName,"pt_all",";nLep",10,-0.5,9.5)->Fill(nTotPT2, weight);
        plotter.getOrMake1DPre(smpName,"pt_dPhiLtPio2",";nLep",10,-0.5,9.5)->Fill(nDPhiPio2PT2, weight);
        plotter.getOrMake1DPre(smpName,"pt_dPhiGtPio2",";nLep",10,-0.5,9.5)->Fill(nTotPT2-nDPhiPio2PT2, weight);
        plotter.getOrMake1DPre(smpName,"pt_drLt0p8",";nLep",10,-0.5,9.5)->Fill(nDR8PT2, weight);

        plotter.getOrMake1DPre(smpName,"drWjjHbb",";JJDR",500,0,5)->Fill(PhysicsUtilities::deltaR(*wjjCand,*hbbCand), weight);

        plotter.getOrMake1DPre(smpName,"drLt0p8",";nLep",10,-0.5,9.5)->Fill(nDR8, weight);
        plotter.getOrMake1DPre(smpName,"drLt0p8_pt",";oLep_pt",300,-1,2)->Fill(maxOPT == 0 ? -0.1 : maxOPT/wjjCand->pt(), weight);
        if(maxOPT > 0)
            plotter.getOrMake1DPre(smpName,"wjjMass",";wjjMass",200,0,200)->Fill(wjjCand->sdMom().mass(), weight);

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    std::unique_ptr<LeptonProcessor>     leptonNoISOProc ;

};

#endif

void getLeptonSub(std::string fileName, int treeInt, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.analyze();
    a.write(outFileName);
}
void getLeptonSub(std::string fileName, int treeInt, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
